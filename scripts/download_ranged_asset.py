#!/usr/bin/env python3
"""Download one public asset with strict parallel HTTP Range validation."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import os
from pathlib import Path
import shutil
import time
from typing import Sequence
from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit
from urllib.request import Request, urlopen


def with_query(url: str, **values: str) -> str:
    parts = urlsplit(url)
    query = dict(parse_qsl(parts.query, keep_blank_values=True))
    query.update(values)
    return urlunsplit((parts.scheme, parts.netloc, parts.path, urlencode(query), parts.fragment))


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_segment(
    *, url: str, start: int, end: int, path: Path, index: int, retries: int
) -> None:
    expected = end - start + 1
    if path.is_file() and path.stat().st_size == expected:
        return
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    for attempt in range(1, retries + 1):
        temporary.unlink(missing_ok=True)
        request_url = with_query(
            url,
            range_part=str(index),
            range_attempt=str(attempt),
            range_nonce="plan_fixA0828",
        )
        request = Request(
            request_url,
            headers={
                "Range": f"bytes={start}-{end}",
                "Accept-Encoding": "identity",
                "User-Agent": "NFL_AB_design-runtime-asset-downloader/1",
            },
        )
        try:
            with urlopen(request, timeout=90) as response:
                status = getattr(response, "status", None)
                content_range = response.headers.get("Content-Range")
                expected_range = f"bytes {start}-{end}/"
                if status != 206 or not content_range or not content_range.startswith(
                    expected_range
                ):
                    raise RuntimeError(
                        f"range {index} expected HTTP 206 and {expected_range!r}; "
                        f"got status={status}, Content-Range={content_range!r}"
                    )
                with temporary.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            observed = temporary.stat().st_size
            if observed != expected:
                raise RuntimeError(
                    f"range {index} size mismatch: expected {expected}, got {observed}"
                )
            os.replace(temporary, path)
            return
        except Exception:
            if attempt == retries:
                raise
            time.sleep(5)
    raise AssertionError("unreachable")


def download(
    *, url: str, expected_bytes: int, staging_dir: Path, output: Path, parts: int
) -> tuple[str, int]:
    staging_dir.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    segment_size = (expected_bytes + parts - 1) // parts
    segments: list[tuple[int, int, Path, int]] = []
    for index in range(parts):
        start = index * segment_size
        if start >= expected_bytes:
            break
        end = min(expected_bytes - 1, start + segment_size - 1)
        segments.append(
            (start, end, staging_dir / f"{output.name}.segment-{index:03d}", index)
        )
    with ThreadPoolExecutor(max_workers=len(segments)) as pool:
        list(
            pool.map(
                lambda row: fetch_segment(
                    url=url,
                    start=row[0],
                    end=row[1],
                    path=row[2],
                    index=row[3],
                    retries=5,
                ),
                segments,
            )
        )
    partial = output.with_name(f".{output.name}.combining-{os.getpid()}")
    with partial.open("wb") as destination:
        for _, _, segment, _ in segments:
            with segment.open("rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
    observed = partial.stat().st_size
    if observed != expected_bytes:
        raise RuntimeError(
            f"combined size mismatch: expected {expected_bytes}, got {observed}"
        )
    os.replace(partial, output)
    digest = file_sha256(output)
    for _, _, segment, _ in segments:
        segment.unlink()
    return digest, observed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--expected-bytes", type=int, required=True)
    parser.add_argument("--staging-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--parts", type=int, default=16)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.expected_bytes <= 0 or args.parts <= 0:
        print("ERROR: --expected-bytes and --parts must be positive")
        return 2
    digest, observed = download(
        url=args.url,
        expected_bytes=args.expected_bytes,
        staging_dir=args.staging_dir.expanduser().resolve(),
        output=args.output.expanduser().resolve(),
        parts=args.parts,
    )
    print(f"bytes={observed}")
    print(f"sha256={digest}")
    print(f"output={args.output.expanduser().resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
