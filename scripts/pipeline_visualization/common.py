"""Shared plotting utilities for the student-facing NfL pipeline guide."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Callable

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="nfl-mpl-"))
os.environ.setdefault("XDG_CACHE_HOME", tempfile.mkdtemp(prefix="nfl-cache-"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COLORS = {
    "blue": "#2563EB",
    "cyan": "#0891B2",
    "green": "#16A34A",
    "amber": "#D97706",
    "red": "#DC2626",
    "gray": "#64748B",
    "light": "#E2E8F0",
    "purple": "#7C3AED",
}


def parser(stage: str) -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=f"Visualize pipeline stage: {stage}")
    value.add_argument("--config", default="config/pipeline_visualization.example.json")
    value.add_argument("--output-dir", default="visualizations/pipeline")
    value.add_argument("--strict", action="store_true")
    return value


def load_context(args: argparse.Namespace) -> tuple[Path, dict[str, Any], Path]:
    config_path = Path(args.config).expanduser().resolve()
    if not config_path.is_file():
        raise SystemExit(f"Visualization config does not exist: {config_path}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repo = Path(config.get("repository_root", config_path.parents[1])).expanduser().resolve()
    output = Path(args.output_dir).expanduser()
    if not output.is_absolute():
        output = repo / output
    output.mkdir(parents=True, exist_ok=True)
    return repo, config, output.resolve()


def resolve(repo: Path, config: dict[str, Any], key: str) -> Path | None:
    raw = config.get("paths", {}).get(key)
    if not raw:
        return None
    path = Path(str(raw)).expanduser()
    return (repo / path).resolve() if not path.is_absolute() else path.resolve()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def missing(stage: str, output: Path, required: list[Path | None], strict: bool) -> Path:
    absent = [str(path) if path else "<not configured>" for path in required if not path or not path.exists()]
    if not absent:
        raise RuntimeError("missing() called without a missing path")
    if strict:
        raise SystemExit(f"{stage}: missing required inputs: {absent}")
    fig, ax = plt.subplots(figsize=(10, 4.8))
    ax.axis("off")
    ax.text(0.5, 0.68, stage, ha="center", fontsize=20, weight="bold")
    ax.text(0.5, 0.50, "NOT RUN / DATA NOT CONFIGURED", ha="center", fontsize=15, color=COLORS["amber"])
    ax.text(0.5, 0.26, "\n".join(absent), ha="center", va="center", fontsize=9, color=COLORS["gray"])
    return save(fig, output)


def save(fig: plt.Figure, output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(output)
    return output


def atoms(path: Path) -> list[dict[str, Any]]:
    records = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        try:
            records.append(
                {
                    "atom": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21],
                    "residue": int(line[22:26]),
                    "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                }
            )
        except ValueError:
            continue
    return records


def stage_antigen(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    fragments = resolve(repo, cfg, "antigen_fragments_csv")
    epitopes = resolve(repo, cfg, "epitope_windows_csv")
    target = out / "01_antigen_and_epitopes.png"
    if not fragments or not epitopes or not fragments.exists() or not epitopes.exists():
        return missing("01 Antigen and epitope selection", target, [fragments, epitopes], strict)
    frag_rows, epi_rows = read_csv(fragments), read_csv(epitopes)
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)
    top = sorted(frag_rows, key=lambda row: float(row["antigen_rank"]))[:8]
    axes[0].barh([r["fragment"] for r in top][::-1], [float(r["antigen_confidence_score"]) for r in top][::-1], color=COLORS["blue"])
    axes[0].set(xlabel="Antigen confidence score (proxy)", title="Prioritized NfL fragments")
    axes[0].grid(axis="x", alpha=.2)
    for row in epi_rows:
        start, end = int(row["start"]), int(row["end"])
        selected = row.get("configured_design_target", "").lower() == "true"
        axes[1].plot([start, end], [float(row["epitope_priority_score"])] * 2, lw=7, color=COLORS["green"] if selected else COLORS["gray"], solid_capstyle="round")
        axes[1].text((start + end) / 2, float(row["epitope_priority_score"]) + .7, row["epitope_id"], ha="center", fontsize=8)
    axes[1].set(xlabel="NEFL residue number", ylabel="Epitope priority (proxy)", title="Epitope windows; green = configured target")
    axes[1].grid(alpha=.2)
    fig.suptitle("Stage 01 — Antigen context and epitope selection", fontsize=16, weight="bold")
    return save(fig, target)


def stage_templates(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    path = resolve(repo, cfg, "template_frameworks_csv")
    target = out / "02_antibody_templates.png"
    if not path or not path.exists():
        return missing("02 Antibody framework templates", target, [path], strict)
    rows = read_csv(path)
    cdrs = ["H1", "H2", "H3", "L1", "L2", "L3"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    width = .35
    for idx, row in enumerate(rows):
        regions = json.loads(row["region_coordinates_json"])
        lengths = [int(regions[name]["length"]) for name in cdrs]
        axes[0].bar(np.arange(6) + (idx - .5) * width, lengths, width, label=row["template_id"])
        axes[1].bar(row["template_id"].replace("template_", ""), row["vh_framework_masked"].count("X") + row["vl_framework_masked"].count("X"), color=[COLORS["blue"], COLORS["purple"]][idx % 2])
    axes[0].set_xticks(range(6), cdrs)
    axes[0].set(ylabel="CDR length (aa)", title="Six designed CDR lengths")
    axes[0].legend(fontsize=8)
    axes[1].set(ylabel="Masked residues", title="Framework-only request check")
    fig.suptitle("Stage 02 — Antibody templates", fontsize=16, weight="bold")
    return save(fig, target)


def plot_structure(path: Path, hotspots: list[int], title: str, output: Path) -> Path:
    records = atoms(path)
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection="3d")
    chain_colors = [COLORS["blue"], COLORS["purple"], COLORS["green"], COLORS["amber"]]
    chains = sorted({row["chain"] for row in records})
    for idx, chain in enumerate(chains):
        ca = [row for row in records if row["chain"] == chain and row["atom"] == "CA"]
        if ca:
            xyz = np.stack([row["xyz"] for row in ca])
            ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], lw=1.7, color=chain_colors[idx % len(chain_colors)], label=f"chain {chain}")
    hot = [row for row in records if row["atom"] == "CA" and row["residue"] in hotspots]
    if hot:
        xyz = np.stack([row["xyz"] for row in hot])
        ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=70, color=COLORS["red"], label="hotspots")
        for row in hot:
            ax.text(*row["xyz"], f"{row['chain']}{row['residue']}", fontsize=8)
    ax.set_title(title, weight="bold")
    ax.legend(fontsize=8)
    ax.set_axis_off()
    return save(fig, output)


def stage_target(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    path = resolve(repo, cfg, "target_pdb")
    target = out / "03_target_structure.png"
    if not path or not path.exists():
        return missing("03 Target structure", target, [path], strict)
    return plot_structure(path, list(map(int, cfg.get("hotspots", [368, 372, 375]))), "Stage 03 — NfL target structure and hotspots", target)


def stage_proxy_funnel(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    path = resolve(repo, cfg, "screening_funnel_csv")
    target = out / "04_proxy_screening_funnel.png"
    if not path or not path.exists():
        return missing("04 Proxy screening funnel", target, [path], strict)
    rows = read_csv(path)
    labels = [r["stage"] for r in rows]
    passed, removed = [int(r["pass_count"]) for r in rows], [int(r["removed_count"]) for r in rows]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(labels, passed, label="pass", color=COLORS["green"])
    ax.bar(labels, removed, bottom=passed, label="removed", color=COLORS["red"])
    ax.set(ylabel="Candidate count", title="Stage 04 — Simulated proxy funnel (not real model confidence)")
    ax.legend()
    ax.tick_params(axis="x", rotation=15)
    return save(fig, target)


def stage_handoff(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    manifest = resolve(repo, cfg, "handoff_manifest")
    attestation = resolve(repo, cfg, "runtime_attestation")
    target = out / "05_handoff_and_attestation.png"
    if not manifest or not attestation or not manifest.exists() or not attestation.exists():
        return missing("05 Handoff and runtime attestation", target, [manifest, attestation], strict)
    hand, att = read_json(manifest), read_json(attestation)
    att_names = {row["engine"] for row in att.get("engines", [])}
    rows = []
    for engine in hand.get("engines", []):
        name = engine["engine"]
        rows.append((name, int(engine.get("planned_job_count", 0)), int(engine.get("selected_job_count", 0)), name in att_names))
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.axis("off")
    for idx, (name, planned, selected, attested) in enumerate(rows):
        x = .15 + idx * .34
        color = COLORS["green"] if attested else COLORS["amber"]
        ax.add_patch(plt.Rectangle((x - .11, .34), .22, .34, facecolor=color, alpha=.18, edgecolor=color, lw=2))
        ax.text(x, .58, name, ha="center", weight="bold", fontsize=13)
        ax.text(x, .48, f"planned {planned}\nselected {selected}\nattested {attested}", ha="center", va="center")
        if idx < len(rows) - 1:
            ax.annotate("", xy=(x + .23, .51), xytext=(x + .12, .51), arrowprops={"arrowstyle": "->", "color": COLORS["gray"]})
    ax.text(.5, .84, "Stage 05 — Immutable handoff + independent runtime attestation", ha="center", fontsize=16, weight="bold")
    ax.text(.5, .16, f"handoff: {hand.get('handoff_id', 'unknown')}\nstate: {att.get('attestation_state', 'unknown')}", ha="center", color=COLORS["gray"])
    return save(fig, target)


def inventory_plot(root: Path, title: str, output: Path) -> Path:
    files = [p for p in root.rglob("*") if p.is_file()]
    groups = Counter((p.suffix.lower() or "no suffix") for p in files)
    sizes = Counter()
    for p in files:
        sizes[p.suffix.lower() or "no suffix"] += p.stat().st_size
    labels = sorted(groups, key=lambda key: sizes[key], reverse=True)[:10]
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    axes[0].bar(labels, [groups[k] for k in labels], color=COLORS["blue"])
    axes[0].set(ylabel="Files", title="Artifact counts")
    axes[1].bar(labels, [sizes[k] / 1024**2 for k in labels], color=COLORS["cyan"])
    axes[1].set(ylabel="MiB", title="Artifact volume")
    for ax in axes:
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle(title, fontsize=16, weight="bold")
    return save(fig, output)


def stage_rfantibody(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    root = resolve(repo, cfg, "rfantibody_results")
    target = out / "06_rfantibody_results.png"
    if not root or not root.is_dir():
        return missing("06 RFantibody results", target, [root], strict)
    return inventory_plot(root, "Stage 06 — RFantibody: diffusion → MPNN → RF2", target)


def stage_iggm(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    root = resolve(repo, cfg, "iggm_results")
    target = out / "07_iggm_results.png"
    if not root or not root.is_dir():
        return missing("07 IgGM results", target, [root], strict)
    pdbs, fastas = list(root.rglob("*.pdb")), list(root.rglob("*.fasta"))
    if not pdbs:
        return missing("07 IgGM results", target, [root / "*.pdb"], strict)
    chain_counts = Counter()
    for pdb in pdbs:
        chain_counts.update({chain: len({row["residue"] for row in atoms(pdb) if row["chain"] == chain}) for chain in {row["chain"] for row in atoms(pdb)}})
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    axes[0].bar(chain_counts.keys(), chain_counts.values(), color=COLORS["purple"])
    axes[0].set(xlabel="Chain", ylabel="Residues across PDBs", title="IgGM chain content")
    axes[1].bar(["PDB", "FASTA"], [len(pdbs), len(fastas)], color=[COLORS["blue"], COLORS["green"]])
    axes[1].set(ylabel="Files", title="Canonical companions")
    fig.suptitle("Stage 07 — IgGM paired-chain designs", fontsize=16, weight="bold")
    return save(fig, target)


def stage_seed(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    path = resolve(repo, cfg, "germinal_seed_pdb")
    target = out / "08_germinal_seed_geometry.png"
    if not path or not path.exists():
        return missing("08 Germinal seed geometry", target, [path], strict)
    external = list(map(int, cfg.get("external_hotspots", [89, 93, 96])))
    return plot_structure(path, external + list(map(int, cfg.get("hotspots", []))), "Stage 08 — Hotspot-aware Germinal seed", target)


def stage_optimization(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    root = resolve(repo, cfg, "germinal_results")
    target = out / "09_germinal_optimization.png"
    files = list(root.rglob("*_metrics.csv")) if root and root.is_dir() else []
    if not files:
        return missing("09 Germinal optimization", target, [root], strict)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    requested = [("plddt", "pLDDT"), ("i_ptm", "iPTM"), ("i_pae", "iPAE"), ("loss", "Loss")]
    for file in files:
        rows = read_csv(file)
        for ax, (key, label) in zip(axes.flat, requested):
            values = [float(row[key]) for row in rows if row.get(key, "") not in {"", "None", "nan"}]
            if values:
                ax.plot(values, lw=1.2, alpha=.8, label=file.stem.replace("_metrics", ""))
                ax.set(ylabel=label, xlabel="Optimization record")
    for ax in axes.flat:
        ax.grid(alpha=.2)
    if len(files) <= 6:
        axes[0, 0].legend(fontsize=6)
    fig.suptitle("Stage 09 — Germinal optimization trajectories", fontsize=16, weight="bold")
    return save(fig, target)


def stage_filters(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    report = resolve(repo, cfg, "filter_only_report")
    target = out / "10_filters_and_abmpnn.png"
    if not report or not report.exists():
        return missing("10 Chai/PyRosetta filters and AbMPNN", target, [report], strict)
    data = read_json(report)
    results = data.get("filter_results", {})
    metrics = data.get("metrics", {})
    names = [key.removesuffix("_filter") for key in results]
    values = [1 if results[key] else 0 for key in results]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    axes[0].barh(names, values, color=[COLORS["green"] if value else COLORS["red"] for value in values])
    axes[0].set(xlim=(0, 1.1), xticks=[0, 1], xticklabels=["fail", "pass"], title="Filter decisions")
    chosen = [key for key in ("external_plddt", "external_iptm", "pdockq2", "percent_interface_cdr", "interface_shape_comp") if key in metrics]
    axes[1].barh(chosen, [float(metrics[key]) for key in chosen], color=COLORS["cyan"])
    axes[1].set(title="Selected real-model metrics")
    fig.suptitle("Stage 10 — Chai/PyRosetta filters; AbMPNN follows initial pass", fontsize=15, weight="bold")
    return save(fig, target)


def stage_summary(repo: Path, cfg: dict[str, Any], out: Path, strict: bool) -> Path:
    report = resolve(repo, cfg, "execution_report")
    status = resolve(repo, cfg, "status_json")
    target = out / "11_pipeline_summary.png"
    if not report or not status or not report.exists() or not status.exists():
        return missing("11 Pipeline summary", target, [report, status], strict)
    execution, state = read_json(report), read_json(status)
    jobs = execution.get("jobs", [])
    stages = ["Inputs", "Handoff", "Design", "Cofold", "Filters", "Candidate"]
    colors = [COLORS["green"]] * 5 + ([COLORS["green"]] if state.get("completed_outputs", 0) else [COLORS["gray"]])
    fig, ax = plt.subplots(figsize=(13, 4.5))
    ax.axis("off")
    for idx, (label, color) in enumerate(zip(stages, colors)):
        x = .08 + idx * .17
        ax.add_patch(plt.Rectangle((x - .06, .38), .12, .24, facecolor=color, alpha=.18, edgecolor=color, lw=2))
        ax.text(x, .5, label, ha="center", va="center", weight="bold")
        if idx < len(stages) - 1:
            ax.annotate("", xy=(x + .105, .5), xytext=(x + .065, .5), arrowprops={"arrowstyle": "->", "color": COLORS["gray"]})
    ax.text(.5, .82, "Stage 11 — Real pipeline execution summary", ha="center", fontsize=17, weight="bold")
    ax.text(.5, .18, f"execution={execution.get('status')} | jobs={len(jobs)} | run_state={state.get('state')} | outputs={state.get('completed_outputs')}/{state.get('planned_outputs')} | automatic promotion={state.get('automatic_candidate_promotion')}", ha="center", fontsize=10)
    return save(fig, target)


STAGES: dict[str, Callable[[Path, dict[str, Any], Path, bool], Path]] = {
    "01_antigen": stage_antigen,
    "02_templates": stage_templates,
    "03_target": stage_target,
    "04_proxy_funnel": stage_proxy_funnel,
    "05_handoff": stage_handoff,
    "06_rfantibody": stage_rfantibody,
    "07_iggm": stage_iggm,
    "08_germinal_seed": stage_seed,
    "09_germinal_optimization": stage_optimization,
    "10_filters": stage_filters,
    "11_summary": stage_summary,
}


def run(stage: str) -> int:
    args = parser(stage).parse_args()
    repo, cfg, out = load_context(args)
    STAGES[stage](repo, cfg, out, args.strict)
    return 0


def run_all() -> int:
    args = parser("all").parse_args()
    repo, cfg, out = load_context(args)
    for name, function in STAGES.items():
        function(repo, cfg, out, args.strict)
    return 0
