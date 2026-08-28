#!/usr/bin/env python3
"""Run Germinal's initial filter stack on one existing trajectory.

This recovery utility never designs a new sequence. It is intended to expose
post-hallucination integration failures before an expensive multi-trajectory
retry. Run it with the same pinned/patched Germinal checkout and model assets
that will be attested for the corrected attempt.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pyrosetta
import yaml

from germinal.filters import filter_utils
from germinal.utils import utils
from germinal.utils.io import IO, RunLayout, Trajectory


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--final-config", required=True)
    parser.add_argument("--trajectory-pdb", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--external-target-hotspots", required=True)
    args = parser.parse_args()

    config_path = Path(args.final_config).expanduser().resolve()
    trajectory_pdb = Path(args.trajectory_pdb).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    if not config_path.is_file() or not trajectory_pdb.is_file():
        raise SystemExit("final config and trajectory PDB must both exist")
    if output_dir.exists():
        raise SystemExit(f"refusing to overwrite filter-only output: {output_dir}")

    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    run_settings = dict(config["run_settings"])
    target_settings = dict(config["target_settings"])
    target_settings["external_target_hotspots"] = args.external_target_hotspots
    initial_filter_path = (
        Path(__import__("germinal").__file__).resolve().parents[1]
        / "configs"
        / "filter"
        / "initial"
        / "scfv.yaml"
    )
    initial_filters = yaml.safe_load(initial_filter_path.read_text(encoding="utf-8"))

    sequences = utils.get_sequence_from_pdb(str(trajectory_pdb))
    target_chain = target_settings["target_chain"]
    binder_chain = target_settings["binder_chain"]
    if "," in target_chain:
        raise SystemExit("filter-only validator currently requires one target chain")
    target_len = len(sequences[target_chain])
    trajectory_sequence = sequences[binder_chain]
    design_name = trajectory_pdb.stem

    layout = RunLayout.create(output_dir)
    io = IO(layout)
    trajectory = Trajectory(
        design_name,
        "attempt07_filter_only",
        ",".join(map(str, run_settings["cdr_lengths"])),
        target_settings["target_hotspots"],
    )
    trajectory.set_save_location("trajectories")

    pyrosetta.init(
        f"-ignore_unrecognized_res -ignore_zero_occupancy -mute all "
        f"-holes:dalphaball {run_settings['dalphaball_path']} "
        f"-corrections::beta_nov16 true -relax:default_repeats 1"
    )
    metrics, results, accepted, final_structure, extra = filter_utils.run_filters(
        trajectory,
        run_settings,
        target_settings,
        initial_filters,
        io,
        trajectory_sequence,
        str(trajectory_pdb),
        target_len,
        multi_relax=False,
        select_mode=run_settings.get("af3_structure_select_mode", "best"),
        af3_seed_size=3,
    )
    report = {
        "schema": "nfl_ab_design.germinal_filter_only_validation.v1",
        "source_final_config": str(config_path),
        "source_trajectory_pdb": str(trajectory_pdb),
        "external_target_hotspots": args.external_target_hotspots,
        "accepted_initial_filters": bool(accepted),
        "metrics": _jsonable(metrics),
        "filter_results": _jsonable(results),
        "final_structure": str(final_structure),
        "extra": _jsonable(extra),
    }
    report_path = output_dir / "filter_only_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
