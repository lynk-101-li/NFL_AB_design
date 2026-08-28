#!/usr/bin/env python3
"""Build a deterministic, path-bound runtime asset manifest.

The manifest is intentionally environment-specific and belongs beside the
runtime assets, not in Git.  It records exact file paths, byte sizes, and
SHA-256 digests so a runtime attestation can bind an explainable manifest hash
instead of an opaque logical bundle name.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Sequence


SCHEMA = "nfl_ab_design.runtime_asset_manifest.v1"


class AssetManifestError(ValueError):
    """Raised when a runtime asset manifest cannot be built safely."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_manifest(
    *,
    family: str,
    variant: str,
    source_revision: str,
    loader: str,
    license_state: str,
    runtime_selection_flags: Sequence[str],
    assets: Sequence[str | Path],
) -> dict[str, object]:
    if not family.strip() or not variant.strip() or not source_revision.strip():
        raise AssetManifestError("family, variant, and source_revision are required")
    if not loader.strip() or not license_state.strip():
        raise AssetManifestError("loader and license_state are required")
    resolved_assets = sorted({Path(value).expanduser().resolve() for value in assets})
    if not resolved_assets:
        raise AssetManifestError("at least one asset file is required")
    rows: list[dict[str, object]] = []
    for path in resolved_assets:
        if not path.is_file():
            raise AssetManifestError(f"runtime asset is not a file: {path}")
        rows.append(
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }
        )
    return {
        "schema": SCHEMA,
        "family": family.strip(),
        "variant": variant.strip(),
        "source_revision": source_revision.strip(),
        "runtime_selection_flags": list(runtime_selection_flags),
        "loader": loader.strip(),
        "license_state": license_state.strip(),
        "assets": rows,
    }


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    destination = path.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n"
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    try:
        temporary.write_text(payload, encoding="utf-8")
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", required=True)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--loader", required=True)
    parser.add_argument("--license-state", required=True)
    parser.add_argument("--runtime-selection-flag", action="append", default=[])
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("assets", nargs="+")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        manifest = build_manifest(
            family=args.family,
            variant=args.variant,
            source_revision=args.source_revision,
            loader=args.loader,
            license_state=args.license_state,
            runtime_selection_flags=args.runtime_selection_flag,
            assets=args.assets,
        )
        write_manifest(args.output, manifest)
    except (AssetManifestError, OSError) as exc:
        print(f"ERROR: {exc}")
        return 2
    print(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
