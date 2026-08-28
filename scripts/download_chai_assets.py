#!/usr/bin/env python3
"""Download the complete Chai-1 0.6.1 runtime asset set safely."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Sequence


BASE_URL = "https://chaiassets.com/chai1-inference-depencencies"
ASSETS = (
    ("conformers_v1.apkl", f"{BASE_URL}/conformers_v1.apkl", 124_726_115),
    ("models_v2/feature_embedding.pt", f"{BASE_URL}/models_v2/feature_embedding.pt", 4_587_962),
    ("models_v2/bond_loss_input_proj.pt", f"{BASE_URL}/models_v2/bond_loss_input_proj.pt", 5_502),
    ("models_v2/token_embedder.pt", f"{BASE_URL}/models_v2/token_embedder.pt", 6_286_772),
    ("models_v2/trunk.pt", f"{BASE_URL}/models_v2/trunk.pt", 633_042_135),
    ("models_v2/diffusion_module.pt", f"{BASE_URL}/models_v2/diffusion_module.pt", 476_546_304),
    ("models_v2/confidence_head.pt", f"{BASE_URL}/models_v2/confidence_head.pt", 55_259_496),
    (
        "esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt",
        f"{BASE_URL}/esm2/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt",
        5_678_940_466,
    ),
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _curl_base_args() -> list[str]:
    return [
        "curl",
        "--location",
        "--fail",
        "--show-error",
        "--retry",
        "5",
        "--retry-all-errors",
        "--retry-delay",
        "5",
        "--connect-timeout",
        "15",
        "--speed-time",
        "60",
        "--speed-limit",
        "10240",
    ]


def _download_parallel_ranges(
    *, url: str, partial: Path, expected_bytes: int, parts: int
) -> None:
    segment_size = (expected_bytes + parts - 1) // parts
    segments: list[tuple[int, int, Path]] = []
    for index in range(parts):
        start = index * segment_size
        if start >= expected_bytes:
            break
        end = min(expected_bytes - 1, start + segment_size - 1)
        segments.append((start, end, partial.with_name(f"{partial.name}.segment-{index:02d}")))

    def fetch(segment: tuple[int, int, Path]) -> None:
        start, end, path = segment
        expected = end - start + 1
        observed = path.stat().st_size if path.exists() else 0
        if observed > expected:
            raise RuntimeError(f"oversized range segment: {path}")
        if observed == expected:
            return
        request_start = start + observed
        with path.open("ab") as handle:
            subprocess.run(
                [*_curl_base_args(), "--range", f"{request_start}-{end}", url],
                stdout=handle,
                check=True,
            )
        if path.stat().st_size != expected:
            raise RuntimeError(
                f"range size mismatch for {path}: expected {expected}, "
                f"got {path.stat().st_size}"
            )

    with ThreadPoolExecutor(max_workers=len(segments)) as pool:
        list(pool.map(fetch, segments))
    with partial.open("wb") as destination:
        for _, _, segment_path in segments:
            with segment_path.open("rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
    if partial.stat().st_size != expected_bytes:
        raise RuntimeError(
            f"combined range size mismatch: expected {expected_bytes}, "
            f"got {partial.stat().st_size}"
        )
    for _, _, segment_path in segments:
        segment_path.unlink()


def download_asset(
    *,
    relative: str,
    url: str,
    expected_bytes: int,
    staging: Path,
    final: Path,
    parallel_parts: int,
) -> dict[str, object]:
    destination = final / relative
    partial = staging / f"{relative}.part"
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and destination.stat().st_size == expected_bytes:
        state = "reused_verified_size"
    else:
        if expected_bytes > 400_000_000 and parallel_parts > 1:
            _download_parallel_ranges(
                url=url,
                partial=partial,
                expected_bytes=expected_bytes,
                parts=parallel_parts,
            )
        else:
            subprocess.run(
                [
                    *_curl_base_args(),
                    "--continue-at",
                    "-",
                    "--output",
                    str(partial),
                    url,
                ],
                check=True,
            )
        observed = partial.stat().st_size
        if observed != expected_bytes:
            raise RuntimeError(
                f"size mismatch for {relative}: expected {expected_bytes}, got {observed}"
            )
        os.replace(partial, destination)
        state = "downloaded_verified_and_promoted"
    return {
        "relative_path": relative,
        "path": str(destination),
        "url": url,
        "bytes": destination.stat().st_size,
        "sha256": file_sha256(destination),
        "state": state,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staging-dir", required=True, type=Path)
    parser.add_argument("--final-dir", required=True, type=Path)
    parser.add_argument("--parallel-parts", type=int, default=8)
    parser.add_argument(
        "--terms-accepted",
        action="store_true",
        help="Required explicit acknowledgement before downloading Chai assets",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.terms_accepted:
        print("ERROR: --terms-accepted is required")
        return 2
    if args.parallel_parts <= 0:
        print("ERROR: --parallel-parts must be positive")
        return 2
    staging = args.staging_dir.expanduser().resolve()
    final = args.final_dir.expanduser().resolve()
    staging.mkdir(parents=True, exist_ok=True)
    final.mkdir(parents=True, exist_ok=True)
    rows = []
    for relative, url, expected_bytes in ASSETS:
        print(f"asset={relative} expected_bytes={expected_bytes}", flush=True)
        rows.append(
            download_asset(
                relative=relative,
                url=url,
                expected_bytes=expected_bytes,
                staging=staging,
                final=final,
                parallel_parts=args.parallel_parts,
            )
        )
    receipt = {
        "schema": "nfl_ab_design.chai_asset_download_receipt.v1",
        "chai_lab_version": "0.6.1",
        "terms_accepted_by_user": True,
        "assets": rows,
    }
    receipt_path = final / "asset_download_receipt.json"
    temporary = final / f".{receipt_path.name}.tmp-{os.getpid()}"
    temporary.write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, receipt_path)
    print(f"receipt={receipt_path}")
    print(f"receipt_sha256={file_sha256(receipt_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
