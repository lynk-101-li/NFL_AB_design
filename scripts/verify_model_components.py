#!/usr/bin/env python3
"""Fail-closed verification of pinned third-party model source submodules."""

from __future__ import annotations

import argparse
import configparser
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = ROOT / "config" / "model_components.json"


class VerificationError(RuntimeError):
    pass


def _git(*args: str, cwd: Path) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=cwd, text=True, capture_output=True, check=False
    )
    if completed.returncode:
        raise VerificationError(completed.stderr.strip() or "git command failed")
    return completed.stdout.strip()


def _submodule_urls(root: Path) -> dict[str, str]:
    parser = configparser.ConfigParser()
    parser.read(root / ".gitmodules", encoding="utf-8")
    result: dict[str, str] = {}
    for section in parser.sections():
        path = parser.get(section, "path", fallback="")
        url = parser.get(section, "url", fallback="")
        if path:
            result[path] = url
    return result


def verify(registry_path: Path, *, allow_uninitialized: bool = False) -> dict:
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    if registry.get("schema") != "nfl_ab_design.model_components.v1":
        raise VerificationError("unsupported model component registry schema")
    urls = _submodule_urls(ROOT)
    rows = []
    for item in registry.get("components", []):
        path_text = item["submodule_path"]
        expected_url = item["repository"]
        expected_revision = item["revision"]
        if urls.get(path_text) != expected_url:
            raise VerificationError(f".gitmodules URL mismatch: {path_text}")
        path = ROOT / path_text
        git_marker = path / ".git"
        if not git_marker.exists():
            if allow_uninitialized:
                rows.append({"id": item["id"], "state": "uninitialized", "path": path_text})
                continue
            raise VerificationError(
                f"submodule is not initialized: {path_text}; run "
                "git submodule update --init --recursive"
            )
        actual_revision = _git("rev-parse", "HEAD", cwd=path)
        if actual_revision != expected_revision:
            raise VerificationError(
                f"revision mismatch for {path_text}: {actual_revision} != {expected_revision}"
            )
        tracked_changes = _git("status", "--porcelain", "--untracked-files=no", cwd=path)
        if tracked_changes:
            raise VerificationError(f"tracked files differ from the pinned commit: {path_text}")
        rows.append(
            {
                "id": item["id"],
                "state": "source_verified",
                "path": path_text,
                "revision": actual_revision,
            }
        )
    return {
        "schema": "nfl_ab_design.model_component_verification.v1",
        "registry": str(registry_path),
        "components": rows,
        "all_sources_initialized": all(row["state"] == "source_verified" for row in rows),
        "model_execution_performed": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument("--allow-uninitialized", action="store_true")
    parser.add_argument("--json", action="store_true", dest="as_json")
    args = parser.parse_args()
    try:
        result = verify(args.registry.resolve(), allow_uninitialized=args.allow_uninitialized)
    except (OSError, KeyError, ValueError, VerificationError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if args.as_json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        for row in result["components"]:
            print(f"{row['state']}\t{row['id']}\t{row['path']}")
    return 0 if result["all_sources_initialized"] or args.allow_uninitialized else 1


if __name__ == "__main__":
    raise SystemExit(main())
