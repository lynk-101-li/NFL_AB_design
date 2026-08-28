#!/usr/bin/env python3
"""Print the authoritative real-model execution roster without running it."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Show only engines[].execution_jobs from a unified handoff. "
            "This command never executes a model."
        )
    )
    parser.add_argument("--manifest", required=True)
    parser.add_argument(
        "--engine",
        choices=("RFantibody", "IgGM", "Germinal"),
        action="append",
        help="Repeat to show multiple engines; default shows all selected engines.",
    )
    parser.add_argument("--json", action="store_true", help="Print selected jobs as JSON")
    args = parser.parse_args()

    path = Path(args.manifest).expanduser().resolve()
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema") != "nfl_ab_design.real_model_handoff.v1":
        raise SystemExit("not a unified real-model handoff")
    requested = set(args.engine or ())
    selected = []
    for engine in value.get("engines", []):
        name = engine.get("engine")
        if requested and name not in requested:
            continue
        for job in engine.get("execution_jobs", []):
            if job.get("selected_for_execution") is not True:
                raise SystemExit(f"unsafe manifest: execution job is not selected: {job.get('job_id')}")
            selected.append(job)
    if not selected:
        raise SystemExit("no authoritative execution_jobs matched the requested engine")
    if args.json:
        print(json.dumps(selected, indent=2, ensure_ascii=False))
        return 0
    for job in selected:
        print(f"[{job['engine']}] {job['job_id']}")
        print(f"  template={job['template_id']} epitope={job['epitope_id']}")
        for index, command in enumerate(job["commands"], start=1):
            cwd = command.get("working_directory") or "<inherit>"
            print(f"  {index:02d}. {command['stage']}")
            print(f"      cwd: {cwd}")
            print(f"      argv: {shlex.join(command['argv'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
