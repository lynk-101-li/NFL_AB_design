"""Reproducible NfL epitope-conditioned antibody design workflow.

This module converts the repository project context, antigen-inference
resources, design constraints, and validation antibodies into an executable
de novo design simulation and real-model handoff. It intentionally uses the Python standard
library so it can run in a clean workspace without structure-prediction
dependencies.

The scoring is a deterministic proxy. It is meant to organize the work, create
auditable intermediate tables, and export input templates for external structure
tools such as IgFold/ABodyBuilder3, AlphaFold3, Chai-1, Boltz, and Rosetta.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable

from . import __version__ as PACKAGE_VERSION

SCRIPT_PATH = Path(__file__).resolve()
PACKAGE_ROOT = SCRIPT_PATH.parents[2]
RESOURCE_DIR = PACKAGE_ROOT / "resources"
PROJECT_CONTEXT_DIR = RESOURCE_DIR / "project_context"
ANTIGEN_DIR = RESOURCE_DIR / "antigen_inference"

STORY_PATH = PROJECT_CONTEXT_DIR / "storyline.txt"
RESEARCH_PLAN_PATH = PROJECT_CONTEXT_DIR / "research_plan.txt"
ANTIGEN_REPORT_PATH = ANTIGEN_DIR / "NFL_22kDa_disulfide_dimer_cathepsin_truncation.md"
FRAGMENT_CANDIDATES_PATH = ANTIGEN_DIR / "nfl_22kda_disulfide_dimer_fragment_candidates.csv"
CLEAVAGE_SITES_PATH = ANTIGEN_DIR / "nfl_medium_high_cathepsin_candidate_sites.csv"
GENPEPT_PATH = ANTIGEN_DIR / "nfl_cathepsin_annotated_for_snapgene.gp"
ANTIBODY_FASTA_PATH = PACKAGE_ROOT / "validation" / "experimentally_validated_antibodies.fasta"
TRUNCATION_CONSTRAINTS_PATH = PACKAGE_ROOT / "input" / "antigen_truncation" / "truncation_constraints.json"
ANTIBODY_TEMPLATE_FASTA_PATH = PACKAGE_ROOT / "input" / "antibody_templates" / "template_fv_backgrounds.fasta"
CONFIG_DIR = PACKAGE_ROOT / "config"
EXTERNAL_PIPELINE_CONFIG_PATH = CONFIG_DIR / "external_pipelines.example.json"
DESIGN_CAMPAIGN_CONFIG_PATH = CONFIG_DIR / "design_campaign.json"
OUTPUT_DIR = PACKAGE_ROOT / "outputs"
EXPORT_DIR = OUTPUT_DIR / "exports"
REAL_RUNS_DIR = PACKAGE_ROOT / "real_runs"


AA_HYDROPATHY = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

AA_MASS_AVG = {
    "A": 89.09,
    "R": 174.20,
    "N": 132.12,
    "D": 133.10,
    "C": 121.16,
    "Q": 146.15,
    "E": 147.13,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "L": 131.17,
    "K": 146.19,
    "M": 149.21,
    "F": 165.19,
    "P": 115.13,
    "S": 105.09,
    "T": 119.12,
    "W": 204.23,
    "Y": 181.19,
    "V": 117.15,
}

HYDROPHOBIC = set("AVLIMFWYC")
ALIPHATIC_HYDROPHOBIC = set("AVLIM")
AROMATIC_OR_LEU = set("FWYL")
POLAR = set("NQSTYC")
AROMATIC = set("FWY")
BASIC = set("KRH")
ACIDIC = set("DE")
CANONICAL_AA = "ACDEFGHIKLMNPQRSTVWY"
MUTATION_AA = "ADEGIKLNQRSTVY"
CYS322_POSITION = 322
SUPPORTED_CDR_ANNOTATION_METHOD = "ANARCI 2020.04.23 Chothia"


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str


@dataclass
class Antibody:
    antibody_id: str
    vh_id: str
    vl_id: str
    vh: str
    vl: str
    experimental_status: str = "validated"
    generation_method: str = "experimentally_validated_parent"
    parent_id: str = ""
    mutation_count: int = 0


@dataclass(frozen=True)
class Cdr:
    name: str
    chain: str
    start: int
    end: int
    sequence: str


def clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def ensure_clean_output_dir() -> None:
    """Remove only reproducible proxy artifacts and preserve real-run data.

    Real model handoffs/results must live under ``real_runs/``.  This routine
    nevertheless uses an allowlist instead of deleting ``outputs/`` wholesale,
    so a misplaced external artifact is not silently destroyed.
    """

    reproducible_files = [
        *OUTPUT_DIR.glob("00_antigen_truncation_*"),
        *OUTPUT_DIR.glob("01_antigen_fragment_prioritization.csv"),
        *OUTPUT_DIR.glob("02_epitope_windows.csv"),
        *OUTPUT_DIR.glob("03_template_frameworks.csv"),
        *OUTPUT_DIR.glob("03_antibody_developability.csv"),
        *OUTPUT_DIR.glob("04_backbone_generation.csv"),
        *OUTPUT_DIR.glob("04_candidate_library.csv"),
        *OUTPUT_DIR.glob("05_sequence_candidates.csv"),
        *OUTPUT_DIR.glob("05_candidate_ranking.csv"),
        *OUTPUT_DIR.glob("06_structure_interface_screen.csv"),
        *OUTPUT_DIR.glob("06_sandwich_pair_report.md"),
        *OUTPUT_DIR.glob("07_developability_screen.csv"),
        *OUTPUT_DIR.glob("08_screening_funnel.csv"),
        *OUTPUT_DIR.glob("09_prospective_candidates.csv"),
        *OUTPUT_DIR.glob("10_retrospective_demo_candidates.csv"),
        *OUTPUT_DIR.glob("11_sandwich_pair_*"),
        *OUTPUT_DIR.glob("workflow_report.md"),
    ]
    for path in reproducible_files:
        if path.is_file() or path.is_symlink():
            path.unlink()
    reproducible_directories = [
        OUTPUT_DIR / "exports" / "fasta",
        OUTPUT_DIR / "exports" / "af3_json",
        OUTPUT_DIR / "exports" / "design_requests",
        OUTPUT_DIR / "exports" / "external_jobs",
    ]
    for directory in reproducible_directories:
        if directory.exists():
            shutil.rmtree(directory)
    for path in (
        OUTPUT_DIR / "exports" / "external_tool_manifest.json",
        OUTPUT_DIR / "intermediate" / "source_manifest.json",
        OUTPUT_DIR / "intermediate" / "run_manifest.json",
    ):
        if path.is_file() or path.is_symlink():
            path.unlink()
    (OUTPUT_DIR / "exports" / "fasta").mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "exports" / "af3_json").mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "intermediate").mkdir(parents=True, exist_ok=True)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def parse_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    name: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append(FastaRecord(name=name, sequence="".join(chunks).upper()))
            name = line[1:].strip()
            chunks = []
        else:
            chunks.append(re.sub(r"[^A-Za-z]", "", line))
    if name is not None:
        records.append(FastaRecord(name=name, sequence="".join(chunks).upper()))
    return records


def load_antibodies(path: Path) -> list[Antibody]:
    grouped: dict[str, dict[str, str]] = {}
    ids: dict[str, dict[str, str]] = {}
    order: list[str] = []
    for record in parse_fasta(path):
        parts = record.name.split("|")
        if len(parts) != 3:
            raise ValueError(f"Expected FASTA header clone|chain|sequence_id, got: {record.name}")
        antibody_id, chain, sequence_id = parts
        chain = chain.upper()
        if chain not in {"VH", "VL"}:
            raise ValueError(f"Expected VH or VL chain in header, got: {record.name}")
        if antibody_id not in grouped:
            order.append(antibody_id)
        grouped.setdefault(antibody_id, {})[chain] = record.sequence
        ids.setdefault(antibody_id, {})[chain] = sequence_id

    antibodies: list[Antibody] = []
    for antibody_id in order:
        chains = grouped[antibody_id]
        if "VH" not in chains or "VL" not in chains:
            raise ValueError(f"Antibody {antibody_id} must have both VH and VL sequences")
        antibodies.append(
            Antibody(
                antibody_id=antibody_id,
                vh_id=ids[antibody_id]["VH"],
                vl_id=ids[antibody_id]["VL"],
                vh=chains["VH"],
                vl=chains["VL"],
                parent_id=antibody_id,
            )
        )
    return antibodies


def parse_genpept_sequence(path: Path) -> str:
    in_origin = False
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip()
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin and line.startswith("//"):
            break
        if in_origin:
            chunks.append(re.sub(r"[^A-Za-z]", "", line))
    seq = "".join(chunks).upper()
    if not seq:
        raise ValueError(f"No ORIGIN protein sequence found in {path}")
    return seq


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_truncation_constraints(path: Path) -> dict[str, Any]:
    default = {
        "experimental_observation": {"non_reducing_band_kDa": 22.0},
        "fragment_constraints": {
            "required_cysteine_position": CYS322_POSITION,
            "monomer_mass_min_kDa": 10.0,
            "monomer_mass_max_kDa": 12.6,
            "preferred_region_start": 280,
            "preferred_region_end": 377,
            "core_candidate_fragments": ["280-375", "281-376", "282-377"],
        },
    }
    if not path.exists():
        return default
    loaded = json.loads(path.read_text(encoding="utf-8"))
    default.update({key: value for key, value in loaded.items() if isinstance(value, dict)})
    return default


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_json(data: Any) -> str:
    payload = json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def parse_range(fragment: str) -> tuple[int, int]:
    match = re.fullmatch(r"(\d+)-(\d+)", fragment.strip())
    if not match:
        raise ValueError(f"Cannot parse fragment range: {fragment}")
    start, end = int(match.group(1)), int(match.group(2))
    if start > end:
        raise ValueError(f"Invalid fragment range: {fragment}")
    return start, end


def subseq(full_sequence: str, start: int, end: int) -> str:
    return full_sequence[start - 1 : end]


def sequence_features(sequence: str) -> dict[str, Any]:
    seq = sequence.upper()
    length = len(seq)
    if length == 0:
        raise ValueError("Cannot score an empty sequence")
    counts = Counter(seq)
    invalid = set(counts) - set(CANONICAL_AA)
    if invalid:
        raise ValueError(f"Unsupported amino-acid letters: {sorted(invalid)}")

    net_charge = counts["K"] + counts["R"] + 0.1 * counts["H"] - counts["D"] - counts["E"]
    hydrophobic_fraction = sum(counts[aa] for aa in HYDROPHOBIC) / length
    polar_fraction = sum(counts[aa] for aa in POLAR) / length
    aromatic_fraction = sum(counts[aa] for aa in AROMATIC) / length
    acidic_fraction = sum(counts[aa] for aa in ACIDIC) / length
    basic_fraction = sum(counts[aa] for aa in BASIC) / length
    charged_fraction = sum(counts[aa] for aa in BASIC | ACIDIC) / length
    low_complexity_fraction = max(counts.values()) / length
    gravy = sum(AA_HYDROPATHY[aa] for aa in seq) / length
    mass_kda = (sum(AA_MASS_AVG[aa] for aa in seq) - 18.015 * (length - 1)) / 1000.0
    glyco = find_motifs(seq, r"N[^P][ST]")
    deamidation = find_motifs(seq, r"N[GST]|DG")
    oxidation_count = counts["M"] + counts["W"]

    return {
        "length": length,
        "net_charge_pH7_proxy": round(net_charge, 3),
        "hydrophobic_fraction": round(hydrophobic_fraction, 4),
        "polar_fraction": round(polar_fraction, 4),
        "aromatic_fraction": round(aromatic_fraction, 4),
        "acidic_fraction": round(acidic_fraction, 4),
        "basic_fraction": round(basic_fraction, 4),
        "charged_fraction": round(charged_fraction, 4),
        "low_complexity_fraction": round(low_complexity_fraction, 4),
        "gravy": round(gravy, 4),
        "avg_mass_kDa_proxy": round(mass_kda, 3),
        "cysteine_count": counts["C"],
        "n_glyco_motif_count": len(glyco),
        "n_glyco_motifs": ";".join(f"{pos}:{motif}" for pos, motif in glyco),
        "deamidation_motif_count": len(deamidation),
        "deamidation_motifs": ";".join(f"{pos}:{motif}" for pos, motif in deamidation),
        "oxidation_MW_count": oxidation_count,
    }


def find_motifs(sequence: str, pattern: str) -> list[tuple[int, str]]:
    hits: list[tuple[int, str]] = []
    for match in re.finditer(pattern, sequence):
        hits.append((match.start() + 1, match.group(0)))
    return hits


def peptide_mass_kda(sequence: str) -> float:
    seq = sequence.upper()
    if not seq:
        raise ValueError("Cannot calculate mass for an empty sequence")
    return (sum(AA_MASS_AVG[aa] for aa in seq) - 18.015 * (len(seq) - 1)) / 1000.0


def domain_at_cut(cut_after_position: int) -> str:
    pos = cut_after_position
    if pos == 1:
        return "N-terminus/transition"
    if 2 <= pos <= 92:
        return "head"
    if 93 <= pos <= 124:
        return "rod:coil1A"
    if 125 <= pos <= 137:
        return "rod:linker1"
    if 138 <= pos <= 234:
        return "rod:coil1B"
    if 235 <= pos <= 252:
        return "rod:linker12"
    if 253 <= pos <= 271:
        return "rod:coil2A"
    if 272 <= pos <= 280:
        return "rod:linker2"
    if 281 <= pos <= 396:
        return "rod:coil2B"
    if 397 <= pos <= 400:
        return "rod"
    if 401 <= pos <= 543:
        return "tail"
    return "outside_annotated_range"


def padded_aa(sequence: str, position_1_based: int) -> str:
    if position_1_based < 1 or position_1_based > len(sequence):
        return "-"
    return sequence[position_1_based - 1]


def cleavage_context(sequence: str, cut_after_position: int) -> tuple[str, dict[str, str]]:
    residues = {
        "P4": padded_aa(sequence, cut_after_position - 3),
        "P3": padded_aa(sequence, cut_after_position - 2),
        "P2": padded_aa(sequence, cut_after_position - 1),
        "P1": padded_aa(sequence, cut_after_position),
        "P1prime": padded_aa(sequence, cut_after_position + 1),
        "P2prime": padded_aa(sequence, cut_after_position + 2),
        "P3prime": padded_aa(sequence, cut_after_position + 3),
        "P4prime": padded_aa(sequence, cut_after_position + 4),
    }
    context = (
        f"{residues['P4']}{residues['P3']}{residues['P2']}{residues['P1']}"
        f"|{residues['P1prime']}{residues['P2prime']}{residues['P3prime']}{residues['P4prime']}"
    )
    return context, residues


def score_cysteine_cathepsin(residues: dict[str, str]) -> tuple[float, str, str]:
    score = 0.0
    notes: list[str] = []
    bias: list[str] = []
    p4 = residues["P4"]
    p3 = residues["P3"]
    p2 = residues["P2"]
    p1 = residues["P1"]
    p1p = residues["P1prime"]

    if p2 in AROMATIC:
        score += 2.4
        notes.append("P2 aromatic")
        bias.append("CatL/V favored")
    elif p2 in ALIPHATIC_HYDROPHOBIC:
        score += 1.6
        notes.append("P2 aliphatic hydrophobic")
        bias.append("CatB/S/K possible")
    elif p2 in BASIC:
        score += 1.0
        notes.append("P2 basic/dibasic-context")

    if p3 in HYDROPHOBIC:
        score += 0.8
        notes.append("P3 hydrophobic")
    if p4 in HYDROPHOBIC:
        score += 0.5
        notes.append("P4 hydrophobic")

    if p1 in BASIC:
        score += 1.4
        notes.append("P1 basic")
    elif p1 in POLAR or p1 in ACIDIC:
        score += 0.8
        notes.append("P1 polar/acidic")
    elif p1 in HYDROPHOBIC:
        score += 0.7
        notes.append("P1 hydrophobic")

    if p1p == "P":
        score -= 1.2
        notes.append("P1prime Pro unfavorable")
    elif p1p != "-":
        score += 0.4
        notes.append("P1prime tolerated")

    if p1 in BASIC and p1p in BASIC:
        score += 0.6
        notes.append("dibasic bond")
        bias.append("CatL/V dibasic possible")

    return round(min(max(score, 0.0), 4.8), 1), "; ".join(dict.fromkeys(notes)), "; ".join(dict.fromkeys(bias))


def score_cathepsin_d_e(residues: dict[str, str]) -> tuple[float, str]:
    score = 0.0
    notes: list[str] = []
    p2 = residues["P2"]
    p1 = residues["P1"]
    p1p = residues["P1prime"]
    p2p = residues["P2prime"]

    if p1 in AROMATIC_OR_LEU:
        score += 2.6
        notes.append("P1 aromatic/Leu")
    elif p1 in HYDROPHOBIC:
        score += 1.7
        notes.append("P1 hydrophobic")
    elif p1 != "-" and p1 != "P":
        score += 0.7
        notes.append("P1 tolerated")

    if p1p in AROMATIC_OR_LEU:
        score += 2.5
        notes.append("P1prime aromatic/Leu")
    elif p1p in HYDROPHOBIC:
        score += 1.8
        notes.append("P1prime hydrophobic")
    elif p1p == "P":
        score -= 1.4
        notes.append("P1prime Pro unfavorable")
    elif p1p != "-":
        score += 0.4
        notes.append("P1prime tolerated")

    if p2 in HYDROPHOBIC:
        score += 0.4
        notes.append("P2 hydrophobic")
    if p2p in HYDROPHOBIC:
        score += 0.4
        notes.append("P2prime hydrophobic")

    return round(max(score, 0.0), 1), "; ".join(dict.fromkeys(notes))


def cleavage_priority(cysteine_score: float, catd_score: float) -> str:
    best = max(cysteine_score, catd_score)
    if best >= 5.2:
        return "high"
    if best >= 4.0:
        return "medium"
    if best >= 3.0:
        return "weak"
    return "low"


def score_all_cathepsin_sites(sequence: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for cut_after_position in range(1, len(sequence)):
        context, residues = cleavage_context(sequence, cut_after_position)
        cys_score, cys_notes, cys_bias = score_cysteine_cathepsin(residues)
        catd_score, catd_notes = score_cathepsin_d_e(residues)
        p1 = padded_aa(sequence, cut_after_position)
        p1p = padded_aa(sequence, cut_after_position + 1)
        rows.append(
            {
                "cut_after_position": cut_after_position,
                "bond": f"{p1}{cut_after_position}|{p1p}{cut_after_position + 1}",
                "context_P4_to_P4prime": context,
                "domain_at_cut": domain_at_cut(cut_after_position),
                "contains_or_near_C322": "yes" if 281 <= cut_after_position <= 396 else "no",
                "cysteine_cathepsin_like_score": cys_score,
                "cysteine_cathepsin_notes": cys_notes,
                "cysteine_cathepsin_possible_bias": cys_bias,
                "cathepsin_D_E_like_score": catd_score,
                "cathepsin_D_E_notes": catd_notes,
                "best_cathepsin_like_score": round(max(cys_score, catd_score), 1),
                "overall_priority": cleavage_priority(cys_score, catd_score),
            }
        )
    return rows


def infer_truncation_fragments(
    sequence: str,
    cleavage_rows: list[dict[str, Any]],
    constraints: dict[str, Any],
) -> list[dict[str, Any]]:
    fragment_constraints = constraints.get("fragment_constraints", {})
    observation = constraints.get("experimental_observation", {})
    required_cys = int(fragment_constraints.get("required_cysteine_position", CYS322_POSITION))
    monomer_min = float(fragment_constraints.get("monomer_mass_min_kDa", 10.0))
    monomer_max = float(fragment_constraints.get("monomer_mass_max_kDa", 12.6))
    preferred_start = int(fragment_constraints.get("preferred_region_start", 280))
    preferred_end = int(fragment_constraints.get("preferred_region_end", 377))
    core_fragments = set(fragment_constraints.get("core_candidate_fragments", ["280-375", "281-376", "282-377"]))
    target_dimer_mass = float(observation.get("non_reducing_band_kDa", 22.0))
    boundary_sites = [
        row
        for row in cleavage_rows
        if 250 <= int(row["cut_after_position"]) <= 400 and float(row["best_cathepsin_like_score"]) >= 3.0
    ]
    candidates: list[dict[str, Any]] = []

    for n_site in boundary_sites:
        start = int(n_site["cut_after_position"]) + 1
        if not 250 <= start <= required_cys:
            continue
        for c_site in boundary_sites:
            end = int(c_site["cut_after_position"])
            if not required_cys <= end <= 405:
                continue
            if start > end:
                continue
            fragment_sequence = subseq(sequence, start, end)
            length = len(fragment_sequence)
            monomer_mass = peptide_mass_kda(fragment_sequence)
            dimer_mass = 2 * monomer_mass - 0.002
            if not monomer_min <= monomer_mass <= monomer_max:
                continue
            if required_cys < start or required_cys > end:
                continue

            n_score = float(n_site["best_cathepsin_like_score"])
            c_score = float(c_site["best_cathepsin_like_score"])
            combined_score = round(n_score + c_score, 1)
            mass_error = abs(dimer_mass - target_dimer_mass)
            fragment_label = f"{start}-{end}"
            core_bonus = 1.0 if fragment_label in core_fragments else 0.0
            preferred_region_bonus = 1.0 if preferred_start <= start <= required_cys <= end <= preferred_end else 0.0
            mass_score = clamp(1.0 - mass_error / (target_dimer_mass * 0.10))
            boundary_score = clamp(combined_score / 10.2)
            length_score = clamp(1.0 - abs(length - 97) / 32.0)
            truncation_score = 100.0 * (
                0.38 * mass_score
                + 0.30 * boundary_score
                + 0.16 * length_score
                + 0.10 * core_bonus
                + 0.06 * preferred_region_bonus
            )
            candidates.append(
                {
                    "fragment": fragment_label,
                    "length_aa": length,
                    "monomer_avg_mass_kDa": round(monomer_mass, 3),
                    "disulfide_homodimer_avg_mass_kDa": round(dimer_mass, 3),
                    "sequence": fragment_sequence,
                    "N_terminal_cut": n_site["bond"],
                    "N_context": n_site["context_P4_to_P4prime"],
                    "N_best_cathepsin_like_score": n_score,
                    "N_cys_cat_score": n_site["cysteine_cathepsin_like_score"],
                    "N_catD_E_score": n_site["cathepsin_D_E_like_score"],
                    "N_notes": format_boundary_notes(n_site),
                    "C_terminal_cut": c_site["bond"],
                    "C_context": c_site["context_P4_to_P4prime"],
                    "C_best_cathepsin_like_score": c_score,
                    "C_cys_cat_score": c_site["cysteine_cathepsin_like_score"],
                    "C_catD_E_score": c_site["cathepsin_D_E_like_score"],
                    "C_notes": format_boundary_notes(c_site),
                    "combined_boundary_score": combined_score,
                    "mass_error_from_22kDa": round(mass_error, 3),
                    "truncation_candidate_score": round(truncation_score, 2),
                    "contains_C322": "yes",
                }
            )

    candidates.sort(
        key=lambda row: (
            -float(row["truncation_candidate_score"]),
            float(row["mass_error_from_22kDa"]),
            abs(int(row["length_aa"]) - 97),
            row["fragment"],
        )
    )
    return candidates


def format_boundary_notes(site: dict[str, Any]) -> str:
    cys_notes = site.get("cysteine_cathepsin_notes", "")
    catd_notes = site.get("cathepsin_D_E_notes", "")
    notes: list[str] = []
    if cys_notes:
        notes.append(f"Cys-cat: {cys_notes}")
    if catd_notes:
        notes.append(f"D/E: {catd_notes}")
    return " | ".join(notes)


def write_truncation_report(
    cleavage_rows: list[dict[str, Any]],
    medium_high_sites: list[dict[str, Any]],
    fragment_rows: list[dict[str, Any]],
    path: Path,
) -> None:
    top_sites = sorted(medium_high_sites, key=lambda row: float(row["best_cathepsin_like_score"]), reverse=True)[:12]
    top_fragments = fragment_rows[:10]
    body = f"""# NfL Antigen Truncation Inference

## Objective

Infer NfL rod/coil-2B antigen fragments compatible with a non-reducing 22 kDa band interpreted as a Cys322-linked disulfide homodimer.

## Constraints

- Protein sequence: canonical human NEFL/NfL, 543 aa.
- Cys322 is the only cysteine in canonical human NfL.
- A 22 kDa disulfide-linked dimer implies an approximately 11 kDa monomeric fragment.
- Candidate monomer fragments must include Cys322.
- Boundary support is estimated from cathepsin-like cleavage preferences.

## Cleavage Model

Cysteine cathepsin-like scoring emphasizes hydrophobic/aromatic P2 preference, basic P1 support, dibasic context, and P1' Pro penalty. Cathepsin D/E-like scoring emphasizes hydrophobic or aromatic P1/P1' boundaries. Scores are deterministic substrate-preference proxies and require experimental validation.

## Summary

- Total peptide bonds scored: {len(cleavage_rows)}
- Medium/high cleavage candidates: {len(medium_high_sites)}
- 22 kDa/Cys322-compatible fragment candidates: {len(fragment_rows)}

## Top Cleavage Sites

{markdown_table(top_sites, ['cut_after_position', 'bond', 'context_P4_to_P4prime', 'domain_at_cut', 'cysteine_cathepsin_like_score', 'cathepsin_D_E_like_score', 'overall_priority'])}

## Top Fragment Candidates

{markdown_table(top_fragments, ['fragment', 'length_aa', 'monomer_avg_mass_kDa', 'disulfide_homodimer_avg_mass_kDa', 'N_terminal_cut', 'C_terminal_cut', 'combined_boundary_score', 'mass_error_from_22kDa'])}

## Interpretation

The highest-priority fragments cluster around NfL aa 280-377, consistent with a rod/coil-2B fragment containing Cys322. These inferred fragments are passed into the antibody epitope and ranking modules.
"""
    path.write_text(body, encoding="utf-8")


def infer_antigen_truncation(
    sequence: str,
    constraints: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    cleavage_rows = score_all_cathepsin_sites(sequence)
    medium_high_sites = [row for row in cleavage_rows if row["overall_priority"] in {"medium", "high"}]
    fragment_rows = infer_truncation_fragments(sequence, cleavage_rows, constraints)
    return cleavage_rows, medium_high_sites, fragment_rows


def prioritize_antigen_fragments(rows: list[dict[str, str]], constraints: dict[str, Any]) -> list[dict[str, Any]]:
    max_boundary = max(float(row["combined_boundary_score"]) for row in rows)
    fragment_constraints = constraints.get("fragment_constraints", {})
    observation = constraints.get("experimental_observation", {})
    report_core = set(fragment_constraints.get("core_candidate_fragments", ["280-375", "281-376", "282-377"]))
    target_dimer_mass = float(observation.get("non_reducing_band_kDa", 22.0))
    prioritized: list[dict[str, Any]] = []

    for row in rows:
        fragment = row["fragment"]
        length = int(row["length_aa"])
        dimer_mass = float(row["disulfide_homodimer_avg_mass_kDa"])
        boundary = float(row["combined_boundary_score"])
        contains_c322 = 1.0 if row["contains_C322"].lower() == "yes" else 0.0
        mass_score = clamp(1.0 - abs(dimer_mass - target_dimer_mass) / (target_dimer_mass * 0.10))
        boundary_score = boundary / max_boundary
        length_score = clamp(1.0 - abs(length - 97) / 32.0)
        report_bonus = 1.0 if fragment in report_core else 0.0
        confidence = 100.0 * (
            0.34 * mass_score
            + 0.30 * boundary_score
            + 0.16 * length_score
            + 0.12 * contains_c322
            + 0.08 * report_bonus
        )
        prioritized.append(
            {
                "fragment": fragment,
                "length_aa": length,
                "sequence": row["sequence"],
                "monomer_avg_mass_kDa": float(row["monomer_avg_mass_kDa"]),
                "disulfide_homodimer_avg_mass_kDa": dimer_mass,
                "mass_closeness_to_22kDa_score": round(mass_score * 100, 2),
                "combined_boundary_score": boundary,
                "contains_C322": row["contains_C322"],
                "N_terminal_cut": row["N_terminal_cut"],
                "C_terminal_cut": row["C_terminal_cut"],
                "report_core_fragment": "yes" if fragment in report_core else "no",
                "antigen_confidence_score": round(confidence, 2),
            }
        )
    prioritized.sort(key=lambda item: item["antigen_confidence_score"], reverse=True)
    for rank, item in enumerate(prioritized, start=1):
        item["antigen_rank"] = rank
    return prioritized


def epitope_score(sequence: str, start: int, end: int, label: str) -> dict[str, Any]:
    features = sequence_features(sequence)
    length = features["length"]
    hydrophobic = float(features["hydrophobic_fraction"])
    polar = float(features["polar_fraction"])
    charged = float(features["charged_fraction"])
    aromatic = float(features["aromatic_fraction"])
    low_complexity = float(features["low_complexity_fraction"])
    net_charge = float(features["net_charge_pH7_proxy"])
    glyco_count = int(features["n_glyco_motif_count"])
    deamidation_count = int(features["deamidation_motif_count"])

    exposure_proxy = clamp(0.40 + 0.55 * (polar + charged) - 0.20 * max(0.0, hydrophobic - 0.45))
    specificity_proxy = clamp(0.45 + 0.35 * charged + 0.20 * aromatic - 0.10 * max(0.0, low_complexity - 0.25))
    rod_stability_proxy = 0.92 if 90 <= start <= end <= 400 else 0.55
    ptm_low_risk_proxy = clamp(1.0 - 0.20 * glyco_count - 0.08 * deamidation_count)
    anchor_boost = 0.0
    if start <= 322 <= end:
        anchor_boost += 0.35
    if any(abs(start - cut) <= 2 or abs(end - cut) <= 2 for cut in (279, 280, 281, 282, 375, 376, 377)):
        anchor_boost += 0.25
    if "boundary" in label or "helix_surface" in label:
        anchor_boost += 0.15
    anchor_score = clamp(0.35 + anchor_boost)

    score = 100.0 * (
        0.25 * exposure_proxy
        + 0.20 * specificity_proxy
        + 0.20 * rod_stability_proxy
        + 0.20 * anchor_score
        + 0.15 * ptm_low_risk_proxy
    )

    return {
        "epitope_id": label,
        "start": start,
        "end": end,
        "length_aa": length,
        "sequence": sequence,
        "net_charge_pH7_proxy": features["net_charge_pH7_proxy"],
        "hydrophobic_fraction": features["hydrophobic_fraction"],
        "polar_fraction": features["polar_fraction"],
        "aromatic_fraction": features["aromatic_fraction"],
        "exposure_proxy": round(exposure_proxy, 3),
        "nfl_specificity_proxy": round(specificity_proxy, 3),
        "rod_stability_proxy": round(rod_stability_proxy, 3),
        "ptm_low_risk_proxy": round(ptm_low_risk_proxy, 3),
        "anchor_score": round(anchor_score, 3),
        "epitope_priority_score": round(score, 2),
        "notes": epitope_notes(start, end, net_charge, hydrophobic, label),
    }


def epitope_notes(start: int, end: int, net_charge: float, hydrophobic: float, label: str) -> str:
    notes: list[str] = []
    if label.startswith("helix_surface"):
        notes.append("structure-reviewed monomer alpha-helical surface")
    if start <= 282 and end >= 279:
        notes.append("near inferred N-terminal cathepsin boundary")
    if start <= 377 and end >= 368:
        notes.append("near inferred C-terminal cathepsin boundary")
    if net_charge < -1:
        notes.append("acidic surface proxy favors basic antibody contacts")
    if hydrophobic > 0.45:
        notes.append("hydrophobic/coiled-coil contact proxy")
    if label.startswith("sliding"):
        notes.append("sliding-window candidate")
    return "; ".join(notes)


def build_epitope_windows(full_sequence: str, primary_fragment: str) -> list[dict[str, Any]]:
    start, end = parse_range(primary_fragment)
    windows: dict[tuple[int, int, str], dict[str, Any]] = {}
    anchors = [
        ("N_boundary_279_290", 279, 290),
        ("N_boundary_280_291", 280, 291),
        ("helix_surface_323_331", 323, 331),
        ("Cys322_core_319_327", 319, 327),
        ("C_boundary_368_377", 368, 377),
        ("C_boundary_369_380", 369, 380),
    ]
    for label, ep_start, ep_end in anchors:
        if 1 <= ep_start <= ep_end <= len(full_sequence):
            seq = subseq(full_sequence, ep_start, ep_end)
            windows[(ep_start, ep_end, label)] = epitope_score(seq, ep_start, ep_end, label)

    for ep_start in range(start, end - 11 + 1, 4):
        ep_end = ep_start + 11
        seq = subseq(full_sequence, ep_start, ep_end)
        label = f"sliding_{ep_start}_{ep_end}"
        windows[(ep_start, ep_end, label)] = epitope_score(seq, ep_start, ep_end, label)

    epitope_rows = list(windows.values())
    epitope_rows.sort(key=lambda row: row["epitope_priority_score"], reverse=True)
    for rank, row in enumerate(epitope_rows, start=1):
        row["epitope_rank"] = rank
    return epitope_rows


def find_heavy_cdr3(sequence: str) -> tuple[int, int] | None:
    region_start = max(0, len(sequence) - 55)
    tail = sequence[region_start:]
    match = re.search(r"C([A-Z]{3,26}?)(W[A-Z]QG|WGQG)", tail)
    if not match:
        return None
    start = region_start + match.start(1) + 1
    end = region_start + match.end(1)
    return start, end


def find_light_cdr3(sequence: str) -> tuple[int, int] | None:
    region_start = max(0, len(sequence) - 45)
    tail = sequence[region_start:]
    match = re.search(r"C([A-Z]{3,18}?)(FGGG|WGGG)", tail)
    if not match:
        return None
    start = region_start + match.start(1) + 1
    end = region_start + match.end(1)
    return start, end


def safe_slice_cdr(sequence: str, name: str, chain: str, start: int, end: int) -> Cdr:
    start = max(1, min(start, len(sequence)))
    end = max(start, min(end, len(sequence)))
    return Cdr(name=name, chain=chain, start=start, end=end, sequence=sequence[start - 1 : end])


def annotate_cdrs(vh: str, vl: str) -> list[Cdr]:
    cdrs: list[Cdr] = [
        safe_slice_cdr(vh, "HCDR1_proxy", "VH", 26, 35),
        safe_slice_cdr(vh, "HCDR2_proxy", "VH", 50, 65),
        safe_slice_cdr(vl, "LCDR1_proxy", "VL", 24, 34),
        safe_slice_cdr(vl, "LCDR2_proxy", "VL", 50, 56),
    ]
    h3 = find_heavy_cdr3(vh)
    if h3 is None:
        cdrs.append(safe_slice_cdr(vh, "HCDR3_proxy", "VH", max(1, len(vh) - 20), max(1, len(vh) - 8)))
    else:
        cdrs.append(safe_slice_cdr(vh, "HCDR3_proxy", "VH", h3[0], h3[1]))
    l3 = find_light_cdr3(vl)
    if l3 is None:
        cdrs.append(safe_slice_cdr(vl, "LCDR3_proxy", "VL", max(1, len(vl) - 18), max(1, len(vl) - 8)))
    else:
        cdrs.append(safe_slice_cdr(vl, "LCDR3_proxy", "VL", l3[0], l3[1]))
    return cdrs


def cdr_summary(cdrs: list[Cdr]) -> dict[str, str]:
    return {cdr.name: cdr.sequence for cdr in cdrs}


def antibody_developability(antibody: Antibody) -> dict[str, Any]:
    cdrs = annotate_cdrs(antibody.vh, antibody.vl)
    cdr_seq = "".join(cdr.sequence for cdr in cdrs)
    vh_features = sequence_features(antibody.vh)
    vl_features = sequence_features(antibody.vl)
    cdr_features = sequence_features(cdr_seq)

    vh_glyco = find_motifs(antibody.vh, r"N[^P][ST]")
    vl_glyco = find_motifs(antibody.vl, r"N[^P][ST]")
    cdr_ranges = {(cdr.chain, pos) for cdr in cdrs for pos in range(cdr.start, cdr.end + 1)}
    glyco_notes: list[str] = []
    for pos, motif in vh_glyco:
        site = "CDR" if any(("VH", p) in cdr_ranges for p in range(pos, pos + len(motif))) else "FR"
        glyco_notes.append(f"VH:{pos}:{motif}:{site}")
    for pos, motif in vl_glyco:
        site = "CDR" if any(("VL", p) in cdr_ranges for p in range(pos, pos + len(motif))) else "FR"
        glyco_notes.append(f"VL:{pos}:{motif}:{site}")

    unpaired_cys_penalty = 0
    if int(vh_features["cysteine_count"]) % 2 != 0:
        unpaired_cys_penalty += 10
    if int(vl_features["cysteine_count"]) % 2 != 0:
        unpaired_cys_penalty += 10

    hydrophobic_patch_proxy = clamp(
        0.45 * float(cdr_features["hydrophobic_fraction"])
        + 0.25 * max(0.0, float(cdr_features["gravy"]))
        + 0.30 * float(cdr_features["low_complexity_fraction"])
    )
    cdr_oxidation = sum(cdr.sequence.count("M") + cdr.sequence.count("W") for cdr in cdrs)
    glyco_penalty = 7 * len(vh_glyco) + 7 * len(vl_glyco)
    developability_score = 100.0
    developability_score -= glyco_penalty
    developability_score -= 3.0 * int(cdr_features["deamidation_motif_count"])
    developability_score -= 2.0 * cdr_oxidation
    developability_score -= 18.0 * hydrophobic_patch_proxy
    developability_score -= unpaired_cys_penalty
    developability_score = clamp(developability_score / 100.0) * 100.0

    summary = cdr_summary(cdrs)
    return {
        "antibody_id": antibody.antibody_id,
        "vh_id": antibody.vh_id,
        "vl_id": antibody.vl_id,
        "vh_length": len(antibody.vh),
        "vl_length": len(antibody.vl),
        "HCDR1_proxy": summary.get("HCDR1_proxy", ""),
        "HCDR2_proxy": summary.get("HCDR2_proxy", ""),
        "HCDR3_proxy": summary.get("HCDR3_proxy", ""),
        "LCDR1_proxy": summary.get("LCDR1_proxy", ""),
        "LCDR2_proxy": summary.get("LCDR2_proxy", ""),
        "LCDR3_proxy": summary.get("LCDR3_proxy", ""),
        "cdr_net_charge_pH7_proxy": cdr_features["net_charge_pH7_proxy"],
        "cdr_aromatic_fraction": cdr_features["aromatic_fraction"],
        "cdr_hydrophobic_fraction": cdr_features["hydrophobic_fraction"],
        "vh_n_glyco_motifs": ";".join(f"{pos}:{motif}" for pos, motif in vh_glyco),
        "vl_n_glyco_motifs": ";".join(f"{pos}:{motif}" for pos, motif in vl_glyco),
        "n_glyco_risk_notes": ";".join(glyco_notes) if glyco_notes else "none",
        "cdr_oxidation_MW_count": cdr_oxidation,
        "hydrophobic_patch_proxy": round(hydrophobic_patch_proxy, 3),
        "developability_score": round(developability_score, 2),
    }


def candidate_to_row(candidate: Antibody) -> dict[str, Any]:
    return {
        "candidate_id": candidate.antibody_id,
        "parent_id": candidate.parent_id,
        "vh_id": candidate.vh_id,
        "vl_id": candidate.vl_id,
        "experimental_status": candidate.experimental_status,
        "generation_method": candidate.generation_method,
        "mutation_count": candidate.mutation_count,
        "vh_sequence": candidate.vh,
        "vl_sequence": candidate.vl,
    }


def all_cdr_positions(vh: str, vl: str) -> list[tuple[str, int]]:
    positions: list[tuple[str, int]] = []
    for cdr in annotate_cdrs(vh, vl):
        for pos in range(cdr.start, cdr.end + 1):
            positions.append((cdr.chain, pos))
    return positions


def mutate_sequence(sequence: str, positions_1_based: Iterable[int], rng: random.Random, n_mutations: int) -> tuple[str, int]:
    seq = list(sequence)
    positions = list(positions_1_based)
    rng.shuffle(positions)
    changed = 0
    for pos in positions:
        if changed >= n_mutations:
            break
        old = seq[pos - 1]
        if old == "C":
            continue
        choices = [aa for aa in MUTATION_AA if aa != old]
        seq[pos - 1] = rng.choice(choices)
        changed += 1
    return "".join(seq), changed


def make_candidate_library(validated: list[Antibody], variants_per_parent: int = 18) -> list[Antibody]:
    rng = random.Random(20260708)
    candidates: list[Antibody] = list(validated)

    for parent in validated:
        cdr_positions = all_cdr_positions(parent.vh, parent.vl)
        vh_positions = [pos for chain, pos in cdr_positions if chain == "VH"]
        vl_positions = [pos for chain, pos in cdr_positions if chain == "VL"]
        for idx in range(1, variants_per_parent + 1):
            requested = 1 + (idx % 7)
            vh_mutations = max(1, requested // 2)
            vl_mutations = max(0, requested - vh_mutations)
            vh, vh_changed = mutate_sequence(parent.vh, vh_positions, rng, vh_mutations)
            vl, vl_changed = mutate_sequence(parent.vl, vl_positions, rng, vl_mutations)
            total_changed = vh_changed + vl_changed
            candidates.append(
                Antibody(
                    antibody_id=f"{parent.antibody_id}-cdr-perturb-{idx:02d}",
                    vh_id=f"{parent.vh_id}_var{idx:02d}",
                    vl_id=f"{parent.vl_id}_var{idx:02d}",
                    vh=vh,
                    vl=vl,
                    experimental_status="in_silico_decoy",
                    generation_method="deterministic_CDR_perturbation",
                    parent_id=parent.antibody_id,
                    mutation_count=total_changed,
                )
            )

    synthetic_templates = [
        ("negative_low_basic", "SSGSSGSS", "QQSSSTP"),
        ("negative_hydrophobic_patch", "LLVLLVLL", "LLYYLLT"),
        ("negative_short_loop", "GSG", "QQS"),
        ("negative_acidic_cdr", "DEDEGDE", "QDEDDST"),
    ]
    for name, h3, l3 in synthetic_templates:
        parent = validated[0]
        vh = replace_cdr3(parent.vh, "VH", h3)
        vl = replace_cdr3(parent.vl, "VL", l3)
        candidates.append(
            Antibody(
                antibody_id=f"synthetic-{name}",
                vh_id=f"{name}_VH",
                vl_id=f"{name}_VL",
                vh=vh,
                vl=vl,
                experimental_status="in_silico_negative_control",
                generation_method="synthetic_CDR3_negative_control",
                parent_id=parent.antibody_id,
                mutation_count=99,
            )
        )
    return candidates


def replace_cdr3(sequence: str, chain: str, new_cdr3: str) -> str:
    if chain == "VH":
        cdr3 = find_heavy_cdr3(sequence)
    else:
        cdr3 = find_light_cdr3(sequence)
    if cdr3 is None:
        return sequence
    start, end = cdr3
    return sequence[: start - 1] + new_cdr3 + sequence[end:]


def cdr_contact_features(antibody: Antibody) -> dict[str, Any]:
    cdrs = annotate_cdrs(antibody.vh, antibody.vl)
    h3 = next(cdr.sequence for cdr in cdrs if cdr.name == "HCDR3_proxy")
    l3 = next(cdr.sequence for cdr in cdrs if cdr.name == "LCDR3_proxy")
    cdr_seq = "".join(cdr.sequence for cdr in cdrs)
    features = sequence_features(cdr_seq)
    h3_features = sequence_features(h3)
    l3_features = sequence_features(l3)
    return {
        "cdrs": cdrs,
        "h3": h3,
        "l3": l3,
        "cdr_features": features,
        "h3_features": h3_features,
        "l3_features": l3_features,
    }


def score_candidate_against_epitope(candidate: Antibody, epitope: dict[str, Any], dev_score: float) -> dict[str, Any]:
    contact = cdr_contact_features(candidate)
    cdr_features = contact["cdr_features"]
    h3 = contact["h3"]
    ep_features = sequence_features(epitope["sequence"])

    cdr_charge = float(cdr_features["net_charge_pH7_proxy"])
    ep_charge = float(ep_features["net_charge_pH7_proxy"])
    charge_score = clamp(0.50 + 0.065 * cdr_charge * (-ep_charge))

    cdr_aromatic = float(cdr_features["aromatic_fraction"])
    ep_aromatic = float(ep_features["aromatic_fraction"])
    ep_hydrophobic = float(ep_features["hydrophobic_fraction"])
    aromatic_contact = clamp(0.35 + 1.6 * cdr_aromatic + 0.45 * ep_aromatic + 0.20 * ep_hydrophobic)

    h3_len = len(h3)
    ep_len = int(epitope["length_aa"])
    if ep_len <= 11:
        ideal_h3 = 7.0
    elif ep_len <= 13:
        ideal_h3 = 9.0
    else:
        ideal_h3 = 12.0
    geometry_fit = clamp(1.0 - abs(h3_len - ideal_h3) / 12.0)

    epitope_quality = float(epitope["epitope_priority_score"]) / 100.0
    mutation_penalty = min(0.42, 0.045 * candidate.mutation_count)
    consensus = clamp(0.62 + 0.22 * epitope_quality + 0.10 * geometry_fit - mutation_penalty)
    developability = dev_score / 100.0

    binding = 100.0 * (
        0.28 * charge_score
        + 0.22 * aromatic_contact
        + 0.18 * geometry_fit
        + 0.16 * epitope_quality
        + 0.10 * consensus
        + 0.06 * developability
    )
    binding -= 100.0 * mutation_penalty
    binding = clamp(binding / 100.0) * 100.0

    off_target_penalty = 100.0 * clamp(
        0.10 * max(0.0, float(cdr_features["hydrophobic_fraction"]) - 0.45)
        + 0.10 * max(0.0, float(cdr_features["low_complexity_fraction"]) - 0.22)
        + 0.03 * int(cdr_features["n_glyco_motif_count"])
        + 0.006 * candidate.mutation_count
    )
    total = (
        0.48 * binding
        + 0.18 * float(epitope["epitope_priority_score"])
        + 0.14 * (100.0 * consensus)
        + 0.12 * dev_score
        + 0.08 * (100.0 - off_target_penalty)
    )

    return {
        "candidate_id": candidate.antibody_id,
        "parent_id": candidate.parent_id,
        "experimental_status": candidate.experimental_status,
        "generation_method": candidate.generation_method,
        "mutation_count": candidate.mutation_count,
        "best_epitope_id": epitope["epitope_id"],
        "best_epitope_start": epitope["start"],
        "best_epitope_end": epitope["end"],
        "best_epitope_sequence": epitope["sequence"],
        "HCDR3_proxy": h3,
        "HCDR3_length": h3_len,
        "CDR_net_charge_pH7_proxy": cdr_features["net_charge_pH7_proxy"],
        "CDR_aromatic_fraction": cdr_features["aromatic_fraction"],
        "charge_complementarity_score": round(charge_score * 100.0, 2),
        "aromatic_contact_score": round(aromatic_contact * 100.0, 2),
        "geometry_fit_score": round(geometry_fit * 100.0, 2),
        "model_consensus_proxy": round(consensus * 100.0, 2),
        "developability_score": round(dev_score, 2),
        "off_target_penalty_proxy": round(off_target_penalty, 2),
        "binding_confidence_score": round(binding, 2),
        "total_rank_score": round(total, 2),
    }


def rank_candidates(candidates: list[Antibody], epitope_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    candidate_rows: list[dict[str, Any]] = []
    scorable_epitopes = [row for row in epitope_rows if int(row["length_aa"]) <= 18]
    for candidate in candidates:
        dev = antibody_developability(candidate)["developability_score"]
        all_scores = [score_candidate_against_epitope(candidate, epitope, float(dev)) for epitope in scorable_epitopes]
        all_scores.sort(key=lambda row: row["total_rank_score"], reverse=True)
        candidate_rows.append(all_scores[0])
    candidate_rows.sort(key=lambda row: row["total_rank_score"], reverse=True)
    for rank, row in enumerate(candidate_rows, start=1):
        row["rank"] = rank
        row["top_tier"] = "yes" if rank <= max(5, math.ceil(len(candidate_rows) * 0.10)) else "no"
    return candidate_rows


def pair_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> tuple[int, float, int]:
    intersection = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    shorter = min(a_end - a_start + 1, b_end - b_start + 1)
    overlap_ratio = intersection / shorter if shorter else 0.0
    if intersection > 0:
        distance = 0
    else:
        distance = max(a_start, b_start) - min(a_end, b_end) - 1
    return intersection, overlap_ratio, distance


def sandwich_pair_analysis(validated: list[Antibody], epitope_rows: list[dict[str, Any]]) -> dict[str, Any]:
    if len(validated) < 2:
        return {"status": "not_enough_validated_antibodies"}

    ab1, ab2 = validated[:2]
    dev1 = antibody_developability(ab1)["developability_score"]
    dev2 = antibody_developability(ab2)["developability_score"]
    scorable_epitopes = [row for row in epitope_rows if int(row["length_aa"]) <= 18]
    scores1 = [score_candidate_against_epitope(ab1, ep, float(dev1)) for ep in scorable_epitopes]
    scores2 = [score_candidate_against_epitope(ab2, ep, float(dev2)) for ep in scorable_epitopes]

    best_pair: dict[str, Any] | None = None
    for s1 in scores1:
        for s2 in scores2:
            intersection, overlap_ratio, distance = pair_overlap(
                int(s1["best_epitope_start"]),
                int(s1["best_epitope_end"]),
                int(s2["best_epitope_start"]),
                int(s2["best_epitope_end"]),
            )
            non_overlap = 1.0 - overlap_ratio
            distance_score = clamp(distance / 35.0)
            clash_score = 100.0 * clamp(overlap_ratio * 0.85 + max(0.0, 8 - distance) / 40.0)
            pair_score = (
                0.30 * float(s1["binding_confidence_score"])
                + 0.30 * float(s2["binding_confidence_score"])
                + 0.20 * (100.0 * non_overlap)
                + 0.10 * (100.0 * distance_score)
                + 0.10 * ((float(dev1) + float(dev2)) / 2.0)
            )
            item = {
                "antibody_1": ab1.antibody_id,
                "antibody_1_epitope_id": s1["best_epitope_id"],
                "antibody_1_epitope_start": s1["best_epitope_start"],
                "antibody_1_epitope_end": s1["best_epitope_end"],
                "antibody_1_binding_score": s1["binding_confidence_score"],
                "antibody_2": ab2.antibody_id,
                "antibody_2_epitope_id": s2["best_epitope_id"],
                "antibody_2_epitope_start": s2["best_epitope_start"],
                "antibody_2_epitope_end": s2["best_epitope_end"],
                "antibody_2_binding_score": s2["binding_confidence_score"],
                "epitope_overlap_aa": intersection,
                "epitope_overlap_ratio": round(overlap_ratio, 3),
                "linear_epitope_gap_aa": distance,
                "clash_score_proxy": round(clash_score, 2),
                "sandwich_compatibility_score": round(pair_score, 2),
            }
            if best_pair is None or item["sandwich_compatibility_score"] > best_pair["sandwich_compatibility_score"]:
                best_pair = item

    assert best_pair is not None
    if float(dev1) >= float(dev2):
        best_pair["recommended_capture"] = ab1.antibody_id
        best_pair["recommended_detection"] = ab2.antibody_id
        best_pair["orientation_reason"] = "capture uses the antibody with the higher sequence developability proxy; detection antibody should be checked for label-site interference"
    else:
        best_pair["recommended_capture"] = ab2.antibody_id
        best_pair["recommended_detection"] = ab1.antibody_id
        best_pair["orientation_reason"] = "capture uses the antibody with the higher sequence developability proxy; detection antibody should be checked for label-site interference"
    return best_pair


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for name, sequence in records:
        lines.append(f">{name}")
        for index in range(0, len(sequence), 80):
            lines.append(sequence[index : index + 80])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def mask_antibody_cdrs(antibody: Antibody) -> dict[str, Any]:
    """Return framework-only VH/VL templates with all six proxy CDRs masked."""

    aliases = {
        "HCDR1_proxy": "H1",
        "HCDR2_proxy": "H2",
        "HCDR3_proxy": "H3",
        "LCDR1_proxy": "L1",
        "LCDR2_proxy": "L2",
        "LCDR3_proxy": "L3",
    }
    masked = {"VH": list(antibody.vh), "VL": list(antibody.vl)}
    regions: list[dict[str, Any]] = []
    for cdr in annotate_cdrs(antibody.vh, antibody.vl):
        for index in range(cdr.start - 1, cdr.end):
            masked[cdr.chain][index] = "X"
        regions.append(
            {
                "region": aliases[cdr.name],
                "chain": cdr.chain,
                "start_1_based": cdr.start,
                "end_1_based_inclusive": cdr.end,
                "length_aa": cdr.end - cdr.start + 1,
            }
        )
    regions.sort(key=lambda row: ("HL".index(str(row["chain"])[1]), int(row["start_1_based"])))
    source_digest = hashlib.sha256(f"{antibody.vh}|{antibody.vl}".encode("ascii")).hexdigest()
    return {
        "template_id": f"template_{antibody.antibody_id}",
        "framework_source_id": antibody.antibody_id,
        "template_role": "framework_source_only",
        "masked_vh": "".join(masked["VH"]),
        "masked_vl": "".join(masked["VL"]),
        "design_regions": regions,
        "designed_regions": "H1;H2;H3;L1;L2;L3",
        "source_sequence_sha256": source_digest,
    }


def prepared_template_request_rows(
    prepared_rows: list[dict[str, Any]],
    antibodies: list[Antibody],
) -> list[dict[str, Any]]:
    """Translate the exact pipeline masks into normalized adapter templates."""

    source_by_id = {antibody.antibody_id: antibody for antibody in antibodies}
    result: list[dict[str, Any]] = []
    expected_regions = ("H1", "H2", "H3", "L1", "L2", "L3")
    for row in prepared_rows:
        source_id = str(row["framework_source_antibody_id"])
        if source_id not in source_by_id:
            raise ValueError(f"Prepared template refers to unknown framework source: {source_id}")
        region_map = row.get("region_coordinates_json")
        if isinstance(region_map, str):
            region_map = json.loads(region_map)
        if not isinstance(region_map, dict) or set(region_map) != set(expected_regions):
            raise ValueError(f"Prepared template {row.get('template_id')} must contain all six CDR ranges")
        regions: list[dict[str, Any]] = []
        for name in expected_regions:
            spec = region_map[name]
            chain = str(spec["chain"])
            start = int(spec["start"])
            end = int(spec["end"])
            regions.append(
                {
                    "region": name,
                    "chain": chain,
                    "start_1_based": start,
                    "end_1_based_inclusive": end,
                    "length_aa": end - start + 1,
                }
            )
        source = source_by_id[source_id]
        result.append(
            {
                "template_id": str(row["template_id"]),
                "framework_source_id": source_id,
                "template_role": "framework_source_only",
                "masked_vh": str(row["vh_framework_masked"]),
                "masked_vl": str(row["vl_framework_masked"]),
                "design_regions": regions,
                "designed_regions": ";".join(expected_regions),
                "source_sequence_sha256": hashlib.sha256(
                    f"{source.vh}|{source.vl}".encode("ascii")
                ).hexdigest(),
            }
        )
    if len(result) < 2:
        raise ValueError("Real-model requests require at least two prepared framework templates")
    return result


def _configured_epitope_ids(config: dict[str, Any]) -> list[str]:
    raw = config.get("target_epitopes")
    if raw is None:
        raw = config.get("epitopes", [])
    if isinstance(raw, dict):
        raw = raw.get("targets", raw.get("ids", []))
    ids: list[str] = []
    if isinstance(raw, list):
        for item in raw:
            if isinstance(item, str):
                ids.append(item)
            elif isinstance(item, dict):
                value = item.get("epitope_id", item.get("id", ""))
                if value:
                    ids.append(str(value))
    if not ids:
        raise ValueError("Design campaign must configure at least one target epitope")
    if len(ids) != len(set(ids)):
        raise ValueError("Design campaign target epitope IDs must be unique")
    return ids


def export_de_novo_model_requests(
    templates: list[Antibody],
    epitope_rows: list[dict[str, Any]],
    full_sequence: str,
    design_config: dict[str, Any],
    *,
    prepared_template_rows: list[dict[str, Any]] | None = None,
    run_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Export normalized, non-executed RFantibody/IgGM/Germinal requests.

    The repository does not contain an antigen structure or model checkpoints.
    These manifests therefore stop at an explicit blocked/not-run handoff.
    """

    request_dir = EXPORT_DIR / "design_requests"
    fasta_dir = EXPORT_DIR / "fasta"
    request_dir.mkdir(parents=True, exist_ok=True)
    template_rows = (
        prepared_template_request_rows(prepared_template_rows, templates)
        if prepared_template_rows is not None
        else [mask_antibody_cdrs(antibody) for antibody in templates]
    )

    configured_ids = _configured_epitope_ids(design_config)
    by_id = {str(row["epitope_id"]): row for row in epitope_rows}
    selected = [by_id[item] for item in configured_ids if item in by_id]
    if not selected:
        selected = epitope_rows[:2]

    template_fasta = fasta_dir / "design_templates_six_cdr_masked.fasta"
    fasta_records: list[tuple[str, str]] = []
    for row in template_rows:
        fasta_records.extend(
            [
                (f"{row['template_id']}|VH|six_CDR_masked", str(row["masked_vh"])),
                (f"{row['template_id']}|VL|six_CDR_masked", str(row["masked_vl"])),
            ]
        )
    write_fasta(template_fasta, fasta_records)

    epitope_requests = []
    for row in selected:
        start = int(row["start"])
        end = int(row["end"])
        epitope_requests.append(
            {
                "epitope_id": row["epitope_id"],
                "sequence": row["sequence"],
                "start_1_based": start,
                "end_1_based_inclusive": end,
                "candidate_hotspot_residue_indices": list(range(start, end + 1)),
                "selection_basis": "configured_target_epitope",
            }
        )

    common = {
        "schema": "nfl_ab_design.normalized_de_novo_request.v1",
        "campaign_mode": "paired_Fv_six_CDR_de_novo_design",
        "execution_state": "not_run",
        "result_provenance": "adapter_request_only",
        "antigen": {
            "protein": "NEFL",
            "full_sequence": full_sequence,
            "antigen_pdb_path": "",
            "structure_input_state": "blocked_missing_antigen_pdb",
        },
        "epitopes": epitope_requests,
        "templates": template_rows,
        "cdr_annotation": _configured_cdr_annotation_metadata(design_config),
        "important_note": (
            "Only framework residues and CDR masks are exported. Known CDR amino-acid identities are not "
            "included in this generation request. Translate this normalized request to the checked-out model version."
        ),
        "run_metadata": dict(run_metadata or {}),
    }
    engine_specs = [
        {
            "engine": "RFantibody",
            "role": "primary_structure_conditioned_generator",
            "required_external_inputs": [
                "antigen_pdb_path",
                "per_template_Chothia_HLT_PDB",
                "full_antigen_coordinate_to_PDB_residue_map",
                "curated_PDB_hotspots",
                "RFantibody_environment_and_checkpoints",
            ],
            "adapter_state": "blocked_missing_antigen_pdb_and_runtime",
        },
        {
            "engine": "IgGM",
            "role": "paired_VH_VL_template_conditioned_generator",
            "required_external_inputs": [
                "antigen_pdb_path",
                "PDB_antigen_chain_sequence",
                "full_antigen_coordinate_to_PDB_local_position_map",
                "IgGM_environment_and_checkpoints",
            ],
            "adapter_state": "blocked_missing_antigen_pdb_and_runtime",
        },
        {
            "engine": "Germinal",
            "role": "parallel_epitope_conditioned_scFv_generator",
            "required_external_inputs": [
                "antigen_pdb_path",
                "per_template_scFv_coordinate_pdb",
                "full_antigen_coordinate_to_PDB_residue_map",
                "curated_PDB_hotspots",
                "Germinal_environment_and_structure_backend",
            ],
            "adapter_state": "blocked_missing_antigen_pdb_scFv_templates_and_runtime",
            "geometry_note": (
                "Germinal designs a single-chain scFv (VH-linker-VL), not a native two-chain paired Fv. "
                "Its candidates must remain a separate geometry track until converted and revalidated."
            ),
        },
    ]
    request_files: list[str] = []
    request_hashes: dict[str, str] = {}
    for spec in engine_specs:
        request_path = request_dir / f"{spec['engine'].lower()}_design_request.json"
        write_json(request_path, {**common, **spec})
        request_files.append(relative_to_package(request_path))
        request_hashes[spec["engine"]] = sha256_file(request_path)

    index_path = request_dir / "design_request_index.json"
    write_json(
        index_path,
        {
            "schema": "nfl_ab_design.design_request_index.v1",
            "execution_state": "not_run",
            "run_metadata": dict(run_metadata or {}),
            "masked_template_fasta": relative_to_package(template_fasta),
            "request_files": request_files,
            "template_count": len(template_rows),
            "epitope_count": len(epitope_requests),
            "designed_regions": ["H1", "H2", "H3", "L1", "L2", "L3"],
            "engines": [spec["engine"] for spec in engine_specs],
            "request_sha256_by_engine": request_hashes,
            "shared_three_engine_handoff_supported_state": "single_chain_monomer",
            "future_engine_specific_target_question": (
                "single_chain_monomer_only; any oligomeric target requires a separate campaign"
            ),
        },
    )
    return {
        "masked_template_fasta": relative_to_package(template_fasta),
        "design_request_files": request_files,
        "design_request_index": relative_to_package(index_path),
        "design_request_index_sha256": sha256_file(index_path),
        "design_request_sha256_by_engine": request_hashes,
        "design_adapter_execution_state": "not_run",
    }


def af3_json(name: str, chains: list[tuple[str, str]]) -> dict[str, Any]:
    return {
        "name": name,
        "modelSeeds": [1, 2, 3, 4, 5],
        "sequences": [{"protein": {"id": chain_id, "sequence": sequence}} for chain_id, sequence in chains],
        "dialect": "alphafold3",
        "version": 1,
        "notes": "AF3-style template generated by NFL_AB_design; verify against the active AF3 runner schema before submission.",
    }


def relative_to_package(path: Path) -> str:
    try:
        return str(path.relative_to(PACKAGE_ROOT))
    except ValueError:
        return str(path)


def load_external_pipeline_config(config_path: Path | None) -> dict[str, Any]:
    if config_path is None:
        return {"schema": "nfl_ab_design.external_pipelines.v1", "pipelines": []}
    if not config_path.is_file():
        raise FileNotFoundError(f"External pipeline config does not exist: {config_path}")
    loaded = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise ValueError("External pipeline config root must be an object")
    if loaded.get("schema") != "nfl_ab_design.external_pipelines.v1":
        raise ValueError("Unsupported or missing external pipeline config schema")
    pipelines = loaded.get("pipelines")
    if not isinstance(pipelines, list):
        raise ValueError("External pipeline config pipelines must be a list")
    supported_selectors = {
        "design_request_index",
        "design_request_files",
        "masked_template_fasta",
        "candidate_fv_chains_fasta",
        "complex_fastas",
        "af3_json_files",
        "sandwich_fasta",
    }
    for index, item in enumerate(pipelines):
        if not isinstance(item, dict):
            raise ValueError(f"External pipeline entry {index} must be an object")
        if not isinstance(item.get("enabled"), bool):
            raise ValueError(f"External pipeline {item.get('name', index)!r} enabled must be boolean")
        selector = item.get("input_selector")
        if selector not in supported_selectors:
            raise ValueError(
                f"External pipeline {item.get('name', index)!r} uses unsupported input_selector {selector!r}"
            )
        if not str(item.get("name", "")).strip():
            raise ValueError(f"External pipeline entry {index} lacks a name")
        if not str(item.get("command_template", "")).strip():
            raise ValueError(f"External pipeline {item['name']!r} lacks a command_template")
    return loaded


def load_design_campaign_config(config_path: Path | None) -> dict[str, Any]:
    default = {
        "schema_version": "1.0",
        "seed": 20260812,
        "simulation": {
            "mode": "deterministic_proxy_simulation",
            "designs_per_template_epitope": 24,
            "real_model_execution": False,
        },
        "target_epitopes": ["helix_surface_323_331", "C_boundary_368_377"],
    }
    if config_path is None:
        return default
    if not config_path.is_file():
        raise FileNotFoundError(f"Design campaign config does not exist: {config_path}")
    loaded = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict):
        raise ValueError(f"Design campaign config must contain a JSON object: {config_path}")
    return {**default, **loaded}


def _configured_template_specs(config: dict[str, Any]) -> list[dict[str, Any]] | None:
    raw = config.get("templates")
    if raw is None:
        return None
    if not isinstance(raw, list) or not raw:
        raise ValueError("design_campaign.templates must be a non-empty list")
    if not all(isinstance(item, dict) for item in raw):
        raise ValueError("Every design_campaign.templates entry must be an object")
    return [dict(item) for item in raw]


def _configured_cdr_ranges(config: dict[str, Any]) -> dict[str, Any] | None:
    raw = config.get("cdr_ranges")
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError("design_campaign.cdr_ranges must be an object")
    coordinate_system = raw.get("coordinate_system")
    if coordinate_system != "1_based_inclusive_chain_positions":
        raise ValueError(
            "Only cdr_ranges.coordinate_system='1_based_inclusive_chain_positions' is supported"
        )
    _configured_cdr_annotation_metadata(config)
    by_template = raw.get("by_template")
    if not isinstance(by_template, dict) or not by_template:
        raise ValueError("design_campaign.cdr_ranges.by_template must be a non-empty object")
    return {str(key): value for key, value in by_template.items()}


def _configured_cdr_annotation_metadata(config: dict[str, Any]) -> dict[str, Any]:
    """Validate and summarize the pinned ANARCI/Chothia numbering evidence."""

    raw = config.get("cdr_ranges")
    if raw is None:
        return {
            "annotation_method": "legacy_heuristic_fallback",
            "coordinate_system": "1_based_inclusive_chain_positions",
            "annotation_evidence_path": "",
            "annotation_evidence_sha256": "",
        }
    if not isinstance(raw, dict):
        raise ValueError("design_campaign.cdr_ranges must be an object")
    method = raw.get("annotation_method")
    if method != SUPPORTED_CDR_ANNOTATION_METHOD:
        raise ValueError(
            "cdr_ranges.annotation_method must be "
            f"{SUPPORTED_CDR_ANNOTATION_METHOD!r}; got {method!r}"
        )
    evidence_value = raw.get("annotation_evidence_path")
    if not isinstance(evidence_value, str) or not evidence_value.strip():
        raise ValueError("cdr_ranges.annotation_evidence_path must be a non-empty relative path")
    evidence_rel = Path(evidence_value)
    if evidence_rel.is_absolute():
        raise ValueError("cdr_ranges.annotation_evidence_path must be relative to the package root")
    evidence_path = (PACKAGE_ROOT / evidence_rel).resolve()
    try:
        evidence_path.relative_to(PACKAGE_ROOT.resolve())
    except ValueError as exc:
        raise ValueError("cdr_ranges.annotation_evidence_path escapes the package root") from exc
    if not evidence_path.is_file():
        raise FileNotFoundError(f"CDR annotation evidence does not exist: {evidence_path}")
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if not isinstance(evidence, dict):
        raise ValueError("CDR annotation evidence root must be an object")
    if evidence.get("annotation_method") != method:
        raise ValueError("CDR annotation method differs between campaign and evidence")
    tool = evidence.get("tool")
    if not isinstance(tool, dict) or tool.get("name") != "ANARCI" or tool.get("version") != "2020.04.23":
        raise ValueError("CDR annotation evidence must identify ANARCI version 2020.04.23")
    evidence_input = evidence.get("input")
    if not isinstance(evidence_input, dict):
        raise ValueError("CDR annotation evidence lacks input provenance")
    source_value = evidence_input.get("path")
    if source_value != relative_to_package(ANTIBODY_FASTA_PATH):
        raise ValueError("CDR annotation evidence input path does not match the validation FASTA")
    if evidence_input.get("sha256") != sha256_file(ANTIBODY_FASTA_PATH):
        raise ValueError("CDR annotation evidence input hash does not match the validation FASTA")

    by_template = raw.get("by_template")
    evidence_templates = evidence.get("templates")
    if not isinstance(by_template, dict) or not isinstance(evidence_templates, dict):
        raise ValueError("Campaign and evidence must both contain per-template CDR definitions")
    if set(by_template) != set(evidence_templates):
        raise ValueError("CDR annotation evidence template IDs do not match the campaign")
    for template_id, regions in by_template.items():
        if not isinstance(regions, dict):
            raise ValueError(f"Campaign CDR definitions for {template_id} must be an object")
        template_evidence = evidence_templates.get(template_id)
        if not isinstance(template_evidence, dict):
            raise ValueError(f"CDR annotation evidence is missing template {template_id}")
        chains = template_evidence.get("chains")
        if not isinstance(chains, dict):
            raise ValueError(f"CDR annotation evidence is missing chains for {template_id}")
        for name, spec in regions.items():
            if not isinstance(spec, dict):
                raise ValueError(f"Campaign CDR definition {template_id}/{name} must be an object")
            chain = spec.get("chain")
            chain_evidence = chains.get(chain) if isinstance(chain, str) else None
            cdrs = chain_evidence.get("cdrs") if isinstance(chain_evidence, dict) else None
            cdr_evidence = cdrs.get(name) if isinstance(cdrs, dict) else None
            if not isinstance(cdr_evidence, dict):
                raise ValueError(f"CDR annotation evidence is missing {template_id}/{name}")
            if (
                cdr_evidence.get("raw_start_1_based") != spec.get("start")
                or cdr_evidence.get("raw_end_1_based_inclusive") != spec.get("end")
            ):
                raise ValueError(f"Campaign range differs from CDR evidence for {template_id}/{name}")

    return {
        "annotation_method": method,
        "coordinate_system": raw.get("coordinate_system"),
        "annotation_evidence_path": evidence_value,
        "annotation_evidence_sha256": sha256_file(evidence_path),
    }


def resolve_campaign_epitopes(
    config: dict[str, Any],
    computed_rows: list[dict[str, Any]],
    full_sequence: str,
) -> list[dict[str, Any]]:
    """Resolve configured targets in configuration order and verify coordinates."""

    raw = config.get("target_epitopes")
    if raw is None:
        raw = config.get("epitopes", [])
    if not isinstance(raw, list):
        raise ValueError("design_campaign.target_epitopes must be a list")
    by_id = {str(row["epitope_id"]): row for row in computed_rows}
    resolved: list[dict[str, Any]] = []
    for index, item in enumerate(raw, start=1):
        if isinstance(item, str):
            if item not in by_id:
                raise ValueError(f"Configured target epitope is not present in computed windows: {item}")
            resolved.append(dict(by_id[item]))
            continue
        if not isinstance(item, dict):
            raise ValueError("Every target_epitopes entry must be a string or object")
        epitope_id = str(item.get("epitope_id", item.get("id", ""))).strip()
        if not epitope_id:
            raise ValueError(f"target_epitopes[{index}] lacks id/epitope_id")
        try:
            start = int(item["start"])
            end = int(item["end"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Target epitope {epitope_id} requires integer start/end") from exc
        if start < 1 or end < start or end > len(full_sequence):
            raise ValueError(f"Target epitope {epitope_id} coordinates are outside NEFL")
        expected_sequence = subseq(full_sequence, start, end)
        configured_sequence = str(item.get("sequence", expected_sequence)).strip().upper()
        if configured_sequence != expected_sequence:
            raise ValueError(
                f"Target epitope {epitope_id} sequence does not match NEFL {start}-{end}"
            )
        if epitope_id in by_id:
            base = dict(by_id[epitope_id])
            if (
                int(base["start"]) != start
                or int(base["end"]) != end
                or str(base["sequence"]) != expected_sequence
            ):
                raise ValueError(
                    f"Target epitope {epitope_id} conflicts with the computed epitope definition"
                )
        else:
            base = {
                "epitope_rank": index,
                "epitope_id": epitope_id,
                "start": start,
                "end": end,
                "sequence": expected_sequence,
                "epitope_priority_score": float(item.get("epitope_priority_score", 0.0)),
                "notes": str(item.get("notes", "explicit campaign target")),
            }
        resolved.append(base)
    configured_ids = _configured_epitope_ids(config)
    if [str(row["epitope_id"]) for row in resolved] != configured_ids:
        raise ValueError("Resolved epitope order does not match configured target epitope IDs")
    return resolved


def design_context_fragment(
    primary_fragment: str,
    config: dict[str, Any],
    full_sequence: str,
) -> str:
    """Expand the biochemical lead fragment to contain every configured epitope."""

    start, end = parse_range(primary_fragment)
    raw = config.get("target_epitopes")
    if isinstance(raw, list):
        for item in raw:
            if isinstance(item, dict):
                try:
                    target_start = int(item["start"])
                    target_end = int(item["end"])
                except (KeyError, TypeError, ValueError) as exc:
                    raise ValueError("Configured target epitope requires integer start/end") from exc
                start = min(start, target_start)
                end = max(end, target_end)
    if start < 1 or end > len(full_sequence):
        raise ValueError("Modeling antigen context falls outside the NEFL sequence")
    return f"{start}-{end}"


def validate_design_campaign_contract(
    config: dict[str, Any],
    target_epitopes: list[dict[str, Any]],
    validated_antibodies: list[Antibody],
) -> None:
    """Fail closed when a declared campaign field is unsupported or inconsistent."""

    if str(config.get("schema_version", "1.0")) != "1.0":
        raise ValueError("Unsupported design campaign schema_version")
    regions = config.get("design_regions")
    expected_regions = ["H1", "H2", "H3", "L1", "L2", "L3"]
    if regions is not None and list(regions) != expected_regions:
        raise ValueError("This campaign requires design_regions H1,H2,H3,L1,L2,L3 in order")
    preprocessing = config.get("template_preprocessing", {})
    if isinstance(preprocessing, dict):
        if preprocessing.get("mask_source_cdrs_before_generation") is False:
            raise ValueError("Source CDR masking cannot be disabled")
        if preprocessing.get("allow_known_cdr_feature_leakage") is True:
            raise ValueError("Known CDR feature leakage cannot be enabled")
        if preprocessing.get("allow_full_known_sequence_feature_leakage") is True:
            raise ValueError("Known full-sequence feature leakage cannot be enabled")
    simulation = config.get("simulation", {})
    if isinstance(simulation, dict) and simulation.get("real_model_execution") is True:
        raise ValueError("The local design campaign cannot relabel proxy simulation as real-model execution")
    retrospective = config.get("retrospective_controls", {})
    if isinstance(retrospective, dict):
        if retrospective.get("enabled") is False:
            raise ValueError("This workflow contract requires the separate retrospective control track")
        if retrospective.get("eligible_for_prospective_selection") is True:
            raise ValueError("Retrospective controls cannot be eligible for prospective selection")
        if retrospective.get("injection_stage") not in (None, "after_prospective_candidate_ranking"):
            raise ValueError("Retrospective controls must be injected after prospective ranking")
        if retrospective.get("status_label") not in (None, "retrospective_positive_control"):
            raise ValueError("Unsupported retrospective control status_label")
        if retrospective.get("blind_discovery_claim_allowed") is True:
            raise ValueError("A blind-discovery claim cannot be enabled for retrospective controls")
        configured_controls = retrospective.get("known_positive_ids")
        loaded_controls = [antibody.antibody_id for antibody in validated_antibodies]
        if configured_controls is not None and list(configured_controls) != loaded_controls:
            raise ValueError(
                "retrospective_controls.known_positive_ids must match the loaded validation FASTA in order"
            )
        source_fasta = retrospective.get("source_fasta")
        if source_fasta is not None and str(source_fasta) != relative_to_package(ANTIBODY_FASTA_PATH):
            raise ValueError("retrospective_controls.source_fasta does not match the loaded validation FASTA")
    template_specs = _configured_template_specs(config)
    if template_specs is not None:
        for spec in template_specs:
            source_fasta = spec.get("source_fasta")
            if source_fasta is not None and str(source_fasta) != relative_to_package(ANTIBODY_FASTA_PATH):
                raise ValueError(
                    f"Template {spec.get('template_id')} source_fasta does not match the loaded framework source"
                )
    template_count = len(template_specs) if template_specs is not None else 2
    designs_per = _simulation_setting(config, "designs_per_template_epitope", 24)
    if isinstance(simulation, dict):
        expected_combinations = template_count * len(target_epitopes)
        declared_combinations = simulation.get("template_epitope_combinations")
        if declared_combinations is not None and int(declared_combinations) != expected_combinations:
            raise ValueError("simulation.template_epitope_combinations is inconsistent with templates × epitopes")
        declared_designs = simulation.get("planned_prospective_designs")
        if declared_designs is not None and int(declared_designs) != expected_combinations * designs_per:
            raise ValueError("simulation.planned_prospective_designs is inconsistent with campaign scale")


def _simulation_setting(config: dict[str, Any], key: str, default: int) -> int:
    for section_name in ("simulation", "generation", "campaign"):
        section = config.get(section_name, {})
        if isinstance(section, dict) and key in section:
            return int(section[key])
    return int(config.get(key, default))


def select_external_inputs(exported: dict[str, Any], selector: str) -> list[str]:
    value = exported.get(selector)
    if isinstance(value, str):
        return [value] if value else []
    if isinstance(value, list):
        return [str(item) for item in value]
    return []


def render_external_command(template: str, input_path: str, output_dir: str, tool_name: str) -> str:
    return template.format(
        input=input_path,
        input_path=input_path,
        output_dir=output_dir,
        stem=Path(input_path).stem,
        tool=tool_name,
        # Generated command sheets are repository artifacts and must remain
        # portable across student checkouts. They are intended to run from the
        # repository root, so never embed the maintainer's absolute path.
        package_root=".",
    )


def prepare_external_pipeline_handoff(
    exported: dict[str, Any],
    config_path: Path | None,
    *,
    run_id: str,
) -> dict[str, Any]:
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", run_id):
        raise ValueError(f"Unsafe external handoff run_id: {run_id!r}")
    config = load_external_pipeline_config(config_path)
    job_dir = EXPORT_DIR / "external_jobs"
    job_dir.mkdir(parents=True, exist_ok=True)
    jobs: list[dict[str, Any]] = []

    for item in config.get("pipelines", []):
        name = str(item.get("name", "")).strip()
        if not name:
            continue
        selector = str(item.get("input_selector", "")).strip()
        command_template = str(item.get("command_template", "")).strip()
        stage = str(item.get("stage", "external_structure")).strip()
        enabled = item["enabled"]
        selected_inputs = select_external_inputs(exported, selector)
        if not selected_inputs:
            raise ValueError(f"External pipeline {name!r} selector {selector!r} resolved no inputs")
        for input_path in selected_inputs:
            if name == "compile_real_generation_handoff":
                result_dir = REAL_RUNS_DIR / "handoffs" / run_id
            else:
                result_dir = REAL_RUNS_DIR / "results" / run_id / name / Path(input_path).stem
            output_dir = relative_to_package(result_dir)
            command = render_external_command(command_template, input_path, output_dir, name) if command_template else ""
            jobs.append(
                {
                    "run_id": run_id,
                    "stage": stage,
                    "tool": name,
                    "enabled": "yes" if enabled else "no",
                    "input_selector": selector,
                    "input_path": input_path,
                    "output_dir": output_dir,
                    "command": command,
                }
            )

    job_table = job_dir / "pipeline_jobs.tsv"
    with job_table.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["run_id", "stage", "tool", "enabled", "input_selector", "input_path", "output_dir", "command"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(jobs)

    runner = job_dir / "run_external_pipelines.sh"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        "# Generated command sheet for external structure and docking tools.",
        "# Review config/external_pipelines.example.json before enabling jobs.",
        "",
    ]
    for job in jobs:
        lines.append(f"mkdir -p {json.dumps(job['output_dir'])}")
        prefix = "" if job["enabled"] == "yes" else "# "
        command = job["command"] or f"# no command template configured for {job['tool']}"
        lines.append(f"{prefix}{command}")
        lines.append("")
    runner.write_text("\n".join(lines), encoding="utf-8")

    return {
        "config_path": relative_to_package(config_path) if config_path and config_path.exists() else "",
        "job_table": relative_to_package(job_table),
        "runner_script": relative_to_package(runner),
        "job_count": len(jobs),
        "enabled_job_count": sum(1 for job in jobs if job["enabled"] == "yes"),
        "jobs": jobs,
    }


def export_structure_inputs(
    candidates: list[Antibody],
    prioritized_fragments: list[dict[str, Any]],
    full_sequence: str,
    sandwich: dict[str, Any],
    external_config_path: Path | None,
    *,
    template_antibodies: list[Antibody] | None = None,
    prepared_template_rows: list[dict[str, Any]] | None = None,
    sandwich_antibodies: list[Antibody] | None = None,
    epitope_rows: list[dict[str, Any]] | None = None,
    design_config: dict[str, Any] | None = None,
    run_metadata: dict[str, Any] | None = None,
    modeling_fragment: str | None = None,
) -> dict[str, Any]:
    fasta_dir = EXPORT_DIR / "fasta"
    af3_dir = EXPORT_DIR / "af3_json"
    top_fragments = [dict(row) for row in prioritized_fragments[:3]]
    if modeling_fragment is not None and modeling_fragment not in {
        str(row["fragment"]) for row in top_fragments
    }:
        context_start, context_end = parse_range(modeling_fragment)
        top_fragments.insert(
            0,
            {
                "fragment": modeling_fragment,
                "modeling_context_only": True,
                "length_aa": context_end - context_start + 1,
            },
        )
    antigen_records: list[tuple[str, str]] = []
    for row in top_fragments:
        start, end = parse_range(row["fragment"])
        antigen_records.append((f"NEFL_{row['fragment']}_P07196", subseq(full_sequence, start, end)))
    antigen_fasta_path = fasta_dir / "antigen_fragments.fasta"
    write_fasta(antigen_fasta_path, antigen_records)

    antibody_records: list[tuple[str, str]] = []
    for antibody in candidates:
        antibody_records.append((f"{antibody.antibody_id}|VH|{antibody.vh_id}", antibody.vh))
        antibody_records.append((f"{antibody.antibody_id}|VL|{antibody.vl_id}", antibody.vl))
    candidate_fv_path = fasta_dir / "selected_candidate_fv_chains.fasta"
    write_fasta(candidate_fv_path, antibody_records)

    primary_fragment = modeling_fragment or str(top_fragments[0]["fragment"])
    primary_start, primary_end = parse_range(primary_fragment)
    primary_antigen = subseq(full_sequence, primary_start, primary_end)
    exported: dict[str, Any] = {
        "primary_antigen_fragment": primary_fragment,
        "primary_antigen_fragment_role": "modeling_context_covering_all_configured_epitopes",
        "antigen_fragments_fasta": relative_to_package(antigen_fasta_path),
        "candidate_fv_chains_fasta": relative_to_package(candidate_fv_path),
        "complex_fastas": [],
        "sandwich_fasta": "",
        "sandwich_export_scope": "",
        "sandwich_candidate_ids": [],
        "fasta_files": [],
        "af3_json_files": [],
    }

    for antibody in candidates:
        safe_id = antibody.antibody_id.replace("/", "_")
        chains = [
            ("A", primary_antigen),
            ("H", antibody.vh),
            ("L", antibody.vl),
        ]
        fasta_path = fasta_dir / f"complex_{safe_id}_NEFL_{primary_fragment}.fasta"
        write_fasta(
            fasta_path,
            [
                (f"A|NEFL_{primary_fragment}", primary_antigen),
                (f"H|{antibody.antibody_id}|VH", antibody.vh),
                (f"L|{antibody.antibody_id}|VL", antibody.vl),
            ],
        )
        json_path = af3_dir / f"af3_complex_{safe_id}_NEFL_{primary_fragment}.json"
        write_json(json_path, af3_json(f"{safe_id}_NEFL_{primary_fragment}", chains))
        exported["complex_fastas"].append(relative_to_package(fasta_path))
        exported["fasta_files"].append(relative_to_package(fasta_path))
        exported["af3_json_files"].append(relative_to_package(json_path))

    sandwich_pool = sandwich_antibodies if sandwich_antibodies is not None else candidates
    if len(sandwich_pool) >= 2:
        candidate_by_id = {candidate.antibody_id: candidate for candidate in sandwich_pool}
        ab1 = candidate_by_id.get(str(sandwich.get("antibody_1", "")))
        ab2 = candidate_by_id.get(str(sandwich.get("antibody_2", "")))
    else:
        ab1 = None
        ab2 = None
    if ab1 is not None and ab2 is not None and ab1.antibody_id != ab2.antibody_id:
        chains = [
            ("A", primary_antigen),
            ("H", ab1.vh),
            ("L", ab1.vl),
            ("I", ab2.vh),
            ("M", ab2.vl),
        ]
        fasta_path = fasta_dir / f"sandwich_{ab1.antibody_id}_{ab2.antibody_id}_NEFL_{primary_fragment}.fasta"
        write_fasta(
            fasta_path,
            [
                (f"A|NEFL_{primary_fragment}", primary_antigen),
                (f"H|{ab1.antibody_id}|VH", ab1.vh),
                (f"L|{ab1.antibody_id}|VL", ab1.vl),
                (f"I|{ab2.antibody_id}|VH", ab2.vh),
                (f"M|{ab2.antibody_id}|VL", ab2.vl),
            ],
        )
        json_path = af3_dir / f"af3_sandwich_{ab1.antibody_id}_{ab2.antibody_id}_NEFL_{primary_fragment}.json"
        write_json(json_path, af3_json(f"sandwich_{ab1.antibody_id}_{ab2.antibody_id}_NEFL_{primary_fragment}", chains))
        exported["sandwich_fasta"] = relative_to_package(fasta_path)
        exported["sandwich_export_scope"] = str(sandwich.get("claim_scope", "unlabeled"))
        exported["sandwich_candidate_ids"] = [ab1.antibody_id, ab2.antibody_id]
        exported["fasta_files"].append(relative_to_package(fasta_path))
        exported["af3_json_files"].append(relative_to_package(json_path))

    if template_antibodies and epitope_rows is not None:
        exported.update(
            export_de_novo_model_requests(
                templates=template_antibodies,
                epitope_rows=epitope_rows,
                full_sequence=full_sequence,
                design_config=design_config or {},
                prepared_template_rows=prepared_template_rows,
                run_metadata=run_metadata,
            )
        )

    external_handoff = prepare_external_pipeline_handoff(
        exported,
        external_config_path,
        run_id=str((run_metadata or {}).get("run_id", "unversioned_proxy_handoff")),
    )

    manifest = {
        "schema": "nfl_ab_design.external_tool_manifest.v2",
        "run_metadata": dict(run_metadata or {}),
        "purpose": "Structure-tool input handoff for the NfL antibody design workflow.",
        "limitations": [
            "VH/VL Fv chains are exported without constant regions.",
            "RFantibody, IgGM, and Germinal requests were not executed; all remain blocked until validated coordinate inputs and verified model runtimes are supplied.",
            "The selected candidate Fv export contains prospective designs only; the separately labeled sandwich export may contain retrospective positive controls.",
            "AF3 JSON files are schema templates and should be checked against the active runner.",
            "Proxy ranking metrics should be replaced by measured or modeled ipTM, pTM, interface PAE, pDockQ, buried surface area, Rosetta interface dG, and clash metrics when structures are available.",
        ],
        "recommended_tool_order": [
            "RFantibody as the primary native paired-Fv structure-conditioned generator and IgGM as a paired-Fv template-conditioned generator",
            "Germinal as an independent scFv-only design track whose candidates require paired-Fv rebuilding and revalidation",
            "IgFold or ABodyBuilder3 for Fv/Fab sanity checks",
            "AF3, Chai-1, or Boltz co-folding for antibody-antigen complexes",
            "Rosetta relax/interface analyzer for post-prediction interface metrics",
            "Pair-aware trimer prediction for sandwich compatibility",
        ],
        "sandwich_pair_proxy": sandwich,
        "exports": exported,
        "external_pipeline_handoff": external_handoff,
        "external_pipeline_config": (
            {
                "path": relative_to_package(external_config_path),
                "sha256": sha256_file(external_config_path),
            }
            if external_config_path is not None and external_config_path.is_file()
            else None
        ),
    }
    write_json(EXPORT_DIR / "external_tool_manifest.json", manifest)
    return manifest


def markdown_table(rows: list[dict[str, Any]], columns: list[str], limit: int | None = None) -> str:
    visible = rows if limit is None else rows[:limit]
    if not visible:
        return "_No rows._"

    def cell(value: Any) -> str:
        text = str(value)
        return text.replace("|", "\\|").replace("\n", "<br>")

    lines = ["|" + "|".join(columns) + "|", "|" + "|".join("---" for _ in columns) + "|"]
    for row in visible:
        lines.append("|" + "|".join(cell(row.get(column, "")) for column in columns) + "|")
    return "\n".join(lines)


def write_sandwich_report(sandwich: dict[str, Any], path: Path) -> None:
    if sandwich.get("status") == "not_enough_validated_antibodies":
        path.write_text("# Sandwich pair report\n\nNot enough validated antibodies.\n", encoding="utf-8")
        return
    body = f"""# NfL sandwich pair compatibility report

This report uses pair-aware epitope assignment across the two experimentally validated antibodies. Scores are deterministic proxy metrics intended for ranking and handoff to structure-modeling tools.

| Metric | Value |
|---|---:|
| Antibody 1 | {sandwich['antibody_1']} |
| Antibody 1 epitope | {sandwich['antibody_1_epitope_id']} ({sandwich['antibody_1_epitope_start']}-{sandwich['antibody_1_epitope_end']}) |
| Antibody 1 binding proxy | {sandwich['antibody_1_binding_score']} |
| Antibody 2 | {sandwich['antibody_2']} |
| Antibody 2 epitope | {sandwich['antibody_2_epitope_id']} ({sandwich['antibody_2_epitope_start']}-{sandwich['antibody_2_epitope_end']}) |
| Antibody 2 binding proxy | {sandwich['antibody_2_binding_score']} |
| Epitope overlap ratio | {sandwich['epitope_overlap_ratio']} |
| Linear epitope gap aa | {sandwich['linear_epitope_gap_aa']} |
| Clash proxy score | {sandwich['clash_score_proxy']} |
| Sandwich compatibility proxy | {sandwich['sandwich_compatibility_score']} |

Recommended orientation:

- Capture antibody: `{sandwich['recommended_capture']}`
- Detection antibody: `{sandwich['recommended_detection']}`

Reason: {sandwich['orientation_reason']}.

Interpretation: low epitope overlap and a low clash proxy support taking this pair forward into trimer co-folding or docking. Final pair selection should use actual Fab1:NfL:Fab2 structures and interface metrics when available.
"""
    path.write_text(body, encoding="utf-8")


def _row_to_antibody(row: dict[str, Any]) -> Antibody:
    candidate_id = str(row["candidate_id"])
    return Antibody(
        antibody_id=candidate_id,
        vh_id=str(row.get("vh_id", f"{candidate_id}_VH")),
        vl_id=str(row.get("vl_id", f"{candidate_id}_VL")),
        vh=str(row["vh_sequence"]),
        vl=str(row["vl_sequence"]),
        experimental_status=str(row.get("control_status", "prospective_design")),
        generation_method=str(row.get("generation_method", "de_novo_epitope_conditioned")),
        parent_id=str(row.get("framework_source_antibody_id", row.get("template_id", ""))),
    )


def sandwich_pair_ranking(candidate_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rank candidate pairs from candidate-level epitope and simulation evidence."""

    survivors = [row for row in candidate_rows if str(row.get("funnel_status", "pass")) != "fail"]
    if len(survivors) < 2:
        survivors = candidate_rows
    pair_rows: list[dict[str, Any]] = []
    for first, second in combinations(survivors[:12], 2):
        first_start = int(first.get("best_epitope_start", first.get("epitope_start", 0)))
        first_end = int(first.get("best_epitope_end", first.get("epitope_end", first_start)))
        second_start = int(second.get("best_epitope_start", second.get("epitope_start", 0)))
        second_end = int(second.get("best_epitope_end", second.get("epitope_end", second_start)))
        intersection, overlap_ratio, distance = pair_overlap(first_start, first_end, second_start, second_end)
        first_score = float(first.get("total_rank_score", first.get("composite_score", 0.0)))
        second_score = float(second.get("total_rank_score", second.get("composite_score", 0.0)))
        first_dev = float(first.get("developability_score", 0.0))
        second_dev = float(second.get("developability_score", 0.0))
        independent_bonus = min(
            5.0,
            0.025
            * (
                float(first.get("independent_evidence_score", 0.0))
                + float(second.get("independent_evidence_score", 0.0))
            ),
        )
        clash_score = 100.0 * clamp(overlap_ratio * 0.80 + max(0.0, 8 - distance) / 42.0)
        pair_score = clamp(
            (
                0.32 * first_score
                + 0.32 * second_score
                + 0.17 * (100.0 * (1.0 - overlap_ratio))
                + 0.08 * (100.0 * clamp(distance / 35.0))
                + 0.11 * ((first_dev + second_dev) / 2.0)
                + independent_bonus
            )
            / 100.0
        ) * 100.0
        # Capture orientation is a deterministic heuristic, not an experimental
        # fact: first prefer stronger retrospective/functional evidence, then
        # developability, then a content-stable sequence digest.
        def capture_key(row: dict[str, Any]) -> tuple[float, float, int]:
            sequence_digest = int(
                hashlib.sha256(f"{row.get('vh_sequence', '')}|{row.get('vl_sequence', '')}".encode("ascii")).hexdigest()[:12],
                16,
            )
            return (
                float(row.get("independent_evidence_score", 0.0)),
                float(row.get("developability_score", 0.0)),
                -sequence_digest,
            )

        if capture_key(first) >= capture_key(second):
            capture, detection = first, second
        else:
            capture, detection = second, first
        pair_rows.append(
            {
                "antibody_1": first["candidate_id"],
                "antibody_1_control_status": first.get("control_status", ""),
                "antibody_1_epitope_id": first.get("best_epitope_id", ""),
                "antibody_1_epitope_start": first_start,
                "antibody_1_epitope_end": first_end,
                "antibody_1_binding_score": first.get("binding_confidence_score", first_score),
                "antibody_2": second["candidate_id"],
                "antibody_2_control_status": second.get("control_status", ""),
                "antibody_2_epitope_id": second.get("best_epitope_id", ""),
                "antibody_2_epitope_start": second_start,
                "antibody_2_epitope_end": second_end,
                "antibody_2_binding_score": second.get("binding_confidence_score", second_score),
                "epitope_overlap_aa": intersection,
                "epitope_overlap_ratio": round(overlap_ratio, 3),
                "linear_epitope_gap_aa": distance,
                "clash_score_proxy": round(clash_score, 2),
                "sandwich_compatibility_score": round(pair_score, 2),
                "recommended_capture": capture["candidate_id"],
                "recommended_detection": detection["candidate_id"],
                "orientation_reason": (
                    "capture orientation prioritizes explicit retrospective evidence when present, then simulated developability; "
                    "experimental label-site and orientation tests remain required"
                ),
                "data_status": "simulated",
                "claim_scope": "retrospective_demo" if (
                    first.get("control_status") == "retrospective_positive_control"
                    or second.get("control_status") == "retrospective_positive_control"
                ) else "prospective_simulation",
            }
        )
    pair_rows.sort(
        key=lambda row: (
            -float(row["sandwich_compatibility_score"]),
            str(row["antibody_1"]),
            str(row["antibody_2"]),
        )
    )
    for rank, row in enumerate(pair_rows, start=1):
        row["pair_rank"] = rank
    return pair_rows


def _write_legacy_workflow_report(
    story: str,
    research_plan: str,
    cleavage_rows: list[dict[str, Any]],
    medium_high_sites: list[dict[str, Any]],
    inferred_fragments: list[dict[str, Any]],
    prioritized_fragments: list[dict[str, Any]],
    epitope_rows: list[dict[str, Any]],
    dev_rows: list[dict[str, Any]],
    ranking_rows: list[dict[str, Any]],
    sandwich: dict[str, Any],
    manifest: dict[str, Any],
) -> None:
    validated_ranks = [row for row in ranking_rows if row["experimental_status"] == "validated"]
    validated_summary = markdown_table(
        validated_ranks,
        ["rank", "candidate_id", "best_epitope_id", "binding_confidence_score", "developability_score", "total_rank_score", "top_tier"],
    )
    report = f"""# NfL 抗体设计计算流程报告

生成日期：2026-07-08

## 1. 输入材料

- 故事线：`{relative_to_package(STORY_PATH)}`
- 研究计划：`{relative_to_package(RESEARCH_PLAN_PATH)}`
- 抗原推断目录：`{relative_to_package(ANTIGEN_DIR)}`
- 已验证抗体 FASTA：`{relative_to_package(ANTIBODY_FASTA_PATH)}`

本流程先执行 NfL 抗原截断推断，再进入表位、抗体和 sandwich pair 计算。当前工作假设是 NfL rod/coil-2B 的 aa 280-377 附近片段包含 Cys322，并形成二硫键同源二聚体。

## 2. 计算模块

1. 抗原截断推断：对 NfL 全长逐肽键进行 cathepsin-like 打分，并按 Cys322 + 22 kDa 二硫键二聚体约束枚举片段。
2. 抗原结构可靠性分层：用 rod/coil-2B 区域和截断推断结果定义 NfL aa 280-377 的优先抗原上下文。
3. NfL 特异性表位图谱：在候选片段内生成边界/Cys322/滑窗表位，并用暴露度、带电性、rod 稳定性和 PTM 风险代理指标打分。
4. 抗体序列建模准备：解析 VH/VL，做启发式 CDR 注释和 developability 体检。
5. 计算复现筛选：把两株真实抗体和 CDR 扰动阴性候选放在同一评分体系下排序。
6. Sandwich pair 兼容性：对两株真实抗体做 pair-aware 表位分配，计算表位重叠、线性距离和 clash 代理指标。
7. 外部结构工具衔接：导出 Fv/抗原复合物 FASTA、AF3-style JSON、job table 和 shell 任务模板。

## 3. 抗原截断推断

- 全长 NfL 肽键打分数：`{len(cleavage_rows)}`。
- 中高优先级 cathepsin-like 切点数：`{len(medium_high_sites)}`。
- 满足 Cys322 + 22 kDa 二硫键二聚体约束的候选片段数：`{len(inferred_fragments)}`。
- 详细报告：`00_antigen_truncation_report.md`。

{markdown_table(inferred_fragments, ['fragment', 'length_aa', 'monomer_avg_mass_kDa', 'disulfide_homodimer_avg_mass_kDa', 'N_terminal_cut', 'C_terminal_cut', 'combined_boundary_score', 'mass_error_from_22kDa'], limit=8)}

## 4. 抗原片段优先级

{markdown_table(prioritized_fragments, ['antigen_rank', 'fragment', 'length_aa', 'disulfide_homodimer_avg_mass_kDa', 'combined_boundary_score', 'N_terminal_cut', 'C_terminal_cut', 'antigen_confidence_score'], limit=8)}

## 5. 候选表位窗口

{markdown_table(epitope_rows, ['epitope_rank', 'epitope_id', 'start', 'end', 'sequence', 'epitope_priority_score', 'notes'], limit=12)}

## 6. 已验证抗体序列体检

{markdown_table(dev_rows, ['antibody_id', 'HCDR3_proxy', 'LCDR3_proxy', 'n_glyco_risk_notes', 'hydrophobic_patch_proxy', 'developability_score'])}

## 7. 计算复现排序

真实抗体在包含 CDR 扰动 decoy 的候选库中的排名：

{validated_summary}

Top 10 候选：

{markdown_table(ranking_rows, ['rank', 'candidate_id', 'experimental_status', 'best_epitope_id', 'binding_confidence_score', 'developability_score', 'off_target_penalty_proxy', 'total_rank_score'], limit=10)}

## 8. Sandwich Pair 结论

- 抗体 1：`{sandwich.get('antibody_1', 'NA')}`，表位 `{sandwich.get('antibody_1_epitope_id', 'NA')}`。
- 抗体 2：`{sandwich.get('antibody_2', 'NA')}`，表位 `{sandwich.get('antibody_2_epitope_id', 'NA')}`。
- 表位重叠比例：`{sandwich.get('epitope_overlap_ratio', 'NA')}`。
- 线性表位间隔：`{sandwich.get('linear_epitope_gap_aa', 'NA')}` aa。
- sandwich 兼容性代理评分：`{sandwich.get('sandwich_compatibility_score', 'NA')}`。
- 建议 capture：`{sandwich.get('recommended_capture', 'NA')}`。
- 建议 detection：`{sandwich.get('recommended_detection', 'NA')}`。

## 9. 外部结构工具 Handoff

外部工具输入已经写入：

- `{Path(manifest['exports']['fasta_files'][0]).parent}` for FASTA templates
- `{Path(manifest['exports']['af3_json_files'][0]).parent}` for AF3-style JSON templates
- `{manifest['external_pipeline_handoff']['job_table']}` for external pipeline jobs
- `{manifest['external_pipeline_handoff']['runner_script']}` for editable command sheet

当前流程中的以下指标是代理指标，建议由外部结构或实验结果替换：

- Fv/Fab 模型质量：应替换为 IgFold/ABodyBuilder3 的模型质量、CDR loop 收敛性和 VH/VL packing 结果。
- 复合物可信度：应替换为 AF3/Chai-1/Boltz 的 ipTM、pTM、interface PAE、interface pLDDT、pDockQ/DockQ。
- 界面物理量：应替换为 buried surface area、Rosetta interface ΔG、shape complementarity、氢键/盐桥和 buried unsatisfied polar atoms。
- sandwich 空间兼容性：应替换为 Fab1:NfL:Fab2 三元复合物结构的实际 clash、Fab-Fab 最小距离和标记端可及性。

## 10. 方法学边界

本流程用于计算复现、候选排序和外部结构预测输入准备。候选库中的扰动序列、表位分配和多目标评分均为确定性代理计算；正式结论应结合真实结构建模、亲和力测定、交叉反应性实验和 sandwich assay 数据。
"""
    (OUTPUT_DIR / "workflow_report.md").write_text(report, encoding="utf-8")


def _run_legacy_workflow(external_config_path: Path | None = EXTERNAL_PIPELINE_CONFIG_PATH) -> dict[str, Any]:
    ensure_clean_output_dir()

    story = read_text(STORY_PATH)
    research_plan = read_text(RESEARCH_PLAN_PATH)
    antigen_report = read_text(ANTIGEN_REPORT_PATH)
    full_sequence = parse_genpept_sequence(GENPEPT_PATH)
    cleavage_sites = read_csv_dicts(CLEAVAGE_SITES_PATH)
    validated = load_antibodies(ANTIBODY_FASTA_PATH)
    truncation_constraints = load_truncation_constraints(TRUNCATION_CONSTRAINTS_PATH)

    cleavage_rows, medium_high_sites, inferred_fragments = infer_antigen_truncation(full_sequence, truncation_constraints)
    prioritized_fragments = prioritize_antigen_fragments(inferred_fragments, truncation_constraints)
    primary_fragment = prioritized_fragments[0]["fragment"]
    epitope_rows = build_epitope_windows(full_sequence, primary_fragment)
    dev_rows = [antibody_developability(ab) for ab in validated]
    candidates = make_candidate_library(validated)
    candidate_rows = [candidate_to_row(candidate) for candidate in candidates]
    ranking_rows = rank_candidates(candidates, epitope_rows)
    sandwich = sandwich_pair_analysis(validated, epitope_rows)
    manifest = export_structure_inputs(validated, prioritized_fragments, full_sequence, sandwich, external_config_path)

    write_csv(OUTPUT_DIR / "00_antigen_truncation_all_peptide_bonds.csv", cleavage_rows)
    write_csv(OUTPUT_DIR / "00_antigen_truncation_medium_high_sites.csv", medium_high_sites)
    write_csv(OUTPUT_DIR / "00_antigen_truncation_fragment_candidates.csv", inferred_fragments)
    write_truncation_report(
        cleavage_rows=cleavage_rows,
        medium_high_sites=medium_high_sites,
        fragment_rows=inferred_fragments,
        path=OUTPUT_DIR / "00_antigen_truncation_report.md",
    )
    write_csv(OUTPUT_DIR / "01_antigen_fragment_prioritization.csv", prioritized_fragments)
    write_csv(OUTPUT_DIR / "02_epitope_windows.csv", epitope_rows)
    write_csv(OUTPUT_DIR / "03_antibody_developability.csv", dev_rows)
    write_csv(OUTPUT_DIR / "04_candidate_library.csv", candidate_rows)
    write_csv(OUTPUT_DIR / "05_candidate_ranking.csv", ranking_rows)
    write_sandwich_report(sandwich, OUTPUT_DIR / "06_sandwich_pair_report.md")

    write_json(
        OUTPUT_DIR / "intermediate" / "source_manifest.json",
        {
            "story_path": relative_to_package(STORY_PATH),
            "research_plan_path": relative_to_package(RESEARCH_PLAN_PATH),
            "antigen_report_path": relative_to_package(ANTIGEN_REPORT_PATH),
            "fragment_candidates_path": relative_to_package(FRAGMENT_CANDIDATES_PATH),
            "cleavage_sites_path": relative_to_package(CLEAVAGE_SITES_PATH),
            "genpept_path": relative_to_package(GENPEPT_PATH),
            "truncation_constraints_path": relative_to_package(TRUNCATION_CONSTRAINTS_PATH),
            "antibody_template_fasta_path": relative_to_package(ANTIBODY_TEMPLATE_FASTA_PATH),
            "validation_antibody_fasta_path": relative_to_package(ANTIBODY_FASTA_PATH),
            "nfl_sequence_length": len(full_sequence),
            "upstream_cleavage_site_rows": len(cleavage_sites),
            "computed_peptide_bond_rows": len(cleavage_rows),
            "computed_medium_high_site_rows": len(medium_high_sites),
            "computed_fragment_candidate_rows": len(inferred_fragments),
            "antigen_report_characters": len(antigen_report),
        },
    )
    write_workflow_report(
        story=story,
        research_plan=research_plan,
        cleavage_rows=cleavage_rows,
        medium_high_sites=medium_high_sites,
        inferred_fragments=inferred_fragments,
        prioritized_fragments=prioritized_fragments,
        epitope_rows=epitope_rows,
        dev_rows=dev_rows,
        ranking_rows=ranking_rows,
        sandwich=sandwich,
        manifest=manifest,
    )

    return {
        "output_dir": OUTPUT_DIR,
        "primary_fragment": primary_fragment,
        "ranking_rows": ranking_rows,
        "sandwich": sandwich,
        "manifest": manifest,
    }


def write_workflow_report(
    *,
    cleavage_rows: list[dict[str, Any]],
    medium_high_sites: list[dict[str, Any]],
    inferred_fragments: list[dict[str, Any]],
    epitope_rows: list[dict[str, Any]],
    design_result: dict[str, Any],
    pair_rows: list[dict[str, Any]],
    sandwich: dict[str, Any],
    manifest: dict[str, Any],
    run_timestamp: str,
    biochemical_lead_fragment: str,
    modeling_context_fragment: str,
) -> None:
    prospective = design_result["prospective_ranking_rows"]
    selected_prospective = [row for row in prospective if row["selected_for_export"]]
    retrospective = design_result["retrospective_ranking_rows"]
    controls = [row for row in retrospective if row["control_status"] == "retrospective_positive_control"]
    report = f"""# NfL 抗体从头设计与回顾性对照演示报告

运行时间：`{run_timestamp}`

> **证据边界：** 本次生成、结构/界面、可开发性和 sandwich 数值均为确定性 `simulated proxy`。
> RFantibody、IgGM、Germinal、tFold 及后续结构工具均未在本次运行中执行。
> 两株已知抗体在前瞻排名完成后才以 `retrospective_positive_control` 注入；其 Top 2 不是盲法从头发现。

## 1. 输入与设计边界

- 抗原推断资源：`{relative_to_package(ANTIGEN_DIR)}`
- 设计 campaign：`{relative_to_package(DESIGN_CAMPAIGN_CONFIG_PATH)}`
- 回顾性阳性对照：`{relative_to_package(ANTIBODY_FASTA_PATH)}`

工作抗原采用 NfL rod/coil-2B 的 aa280–377 单链单体上下文；设计热点不包含 Cys322，也不要求抗体接触半胱氨酸。两株已知抗体只在生成轨道中提供两个不同的配对 VH/VL framework；H1/H2/H3/L1/L2/L3 全部遮罩并重新设计。已知 CDR 氨基酸和完整已知 VH/VL 不作为 prospective generation feature。

CDR 坐标由 `ANARCI 2020.04.23 Chothia` 编号后映射到链内 1-based inclusive raw 坐标；模拟生成和真实模型请求共用同一组精确遮罩。编号 labels、工具版本与输入哈希见 `input/antibody_templates/chothia_numbering_evidence.json`。

## 2. 抗原截断与表位

- 全长 NfL 肽键 proxy：`{len(cleavage_rows)}`
- 中高优先级 cathepsin-like 切点：`{len(medium_high_sites)}`
- 约束内候选截断片段：`{len(inferred_fragments)}`
- 生化截断排序第一名：`NEFL {biochemical_lead_fragment}`
- 覆盖全部配置表位的建模上下文：`NEFL {modeling_context_fragment}`

{markdown_table(epitope_rows, ['epitope_rank', 'epitope_id', 'start', 'end', 'sequence', 'epitope_priority_score', 'notes'], limit=8)}

## 3. 双模板、六 CDR 与模拟生成

{markdown_table(design_result['template_rows'], ['template_id', 'framework_source_antibody_id', 'template_role', 'design_regions', 'known_cdr_sequences_used_for_generation', 'data_status'])}

本地 proxy 生成 `{len(design_result['generation_rows'])}` 个 prospective candidates。这不是 RFantibody、IgGM 或 Germinal 的实际生成量。

## 4. 分步筛选漏斗

{markdown_table(design_result['funnel_rows'], ['stage_order', 'stage', 'metric', 'threshold', 'input_count', 'pass_count', 'removed_count', 'data_status'])}

`06` 和 `07` 表中的分数均保留 `*_is_simulated=True` 与 `metric_provenance`。

## 5. Prospective 分层短名单

`09_prospective_candidates.csv` 不含两株已知阳性全序列。
最终导出采用 `template × epitope` 分层配额；表中的 `rank` 仍是全局模拟分数排名，不能把入选状态解释为纯全局 Top12。

{markdown_table(selected_prospective, ['rank', 'candidate_id', 'template_id', 'best_epitope_id', 'selection_stratum_rank', 'binding_confidence_score', 'developability_score', 'total_rank_score', 'selection_reason'], limit=12)}

## 6. Retrospective 阳性对照 Top 2

已知阳性仅在 prospective ranking 完成后注入，并用明确的回顾性独立证据字段给分。

{markdown_table(controls, ['rank', 'candidate_id', 'control_status', 'best_epitope_id', 'independent_evidence_score', 'independent_evidence_provenance', 'total_rank_score'])}

## 7. Sandwich pair 模拟优先级

- Top pair：`{sandwich.get('antibody_1', 'NA')}` + `{sandwich.get('antibody_2', 'NA')}`
- 表位重叠：`{sandwich.get('epitope_overlap_ratio', 'NA')}`
- 线性间隔：`{sandwich.get('linear_epitope_gap_aa', 'NA')}` aa
- 兼容性 proxy：`{sandwich.get('sandwich_compatibility_score', 'NA')}`
- 建议 capture/detection：`{sandwich.get('recommended_capture', 'NA')}` / `{sandwich.get('recommended_detection', 'NA')}`

{markdown_table(pair_rows, ['pair_rank', 'antibody_1', 'antibody_2', 'epitope_overlap_ratio', 'linear_epitope_gap_aa', 'sandwich_compatibility_score', 'data_status', 'claim_scope'], limit=8)}

## 8. 真实模型 Handoff

已产生六 CDR 遮罩模板与 RFantibody/IgGM/Germinal 规范化请求，但当前缺少经验证的抗原 PDB、坐标映射、模型 runtime 和 checkpoint，所以保持 `not_run/blocked` 状态。Germinal 是独立 scFv 轨道，不视为 native paired-Fv 结果。

- 遮罩模板：`{manifest['exports'].get('masked_template_fasta', 'NA')}`
- 请求索引：`{manifest['exports'].get('design_request_index', 'NA')}`
- job table：`{manifest['external_pipeline_handoff']['job_table']}`
- command sheet：`{manifest['external_pipeline_handoff']['runner_script']}`

## 9. 下一步才能取代 proxy 的证据

- RFantibody/IgGM/Germinal 真实生成结果、日志、版本和 checkpoint。
- Fv/Fab 结构质量、复合物 PAE/ipTM/pTM/DockQ、埋藏表面积和界面能量。
- 亲和力、特异性、交叉反应、可开发性和 sandwich assay 实验。
"""
    (OUTPUT_DIR / "workflow_report.md").write_text(report, encoding="utf-8")


def run_workflow(
    external_config_path: Path | None = EXTERNAL_PIPELINE_CONFIG_PATH,
    design_config_path: Path | None = DESIGN_CAMPAIGN_CONFIG_PATH,
) -> dict[str, Any]:
    from .design_pipeline import run_design_pipeline

    # Validate explicit config paths before touching any existing artifacts.
    design_config = load_design_campaign_config(design_config_path)
    load_external_pipeline_config(external_config_path)
    story = read_text(STORY_PATH)
    research_plan = read_text(RESEARCH_PLAN_PATH)
    antigen_report = read_text(ANTIGEN_REPORT_PATH)
    full_sequence = parse_genpept_sequence(GENPEPT_PATH)
    cleavage_sites = read_csv_dicts(CLEAVAGE_SITES_PATH)
    validated = load_antibodies(ANTIBODY_FASTA_PATH)
    truncation_constraints = load_truncation_constraints(TRUNCATION_CONSTRAINTS_PATH)
    run_timestamp = datetime.now().astimezone().isoformat(timespec="microseconds")
    run_id = "nfl_design_" + re.sub(r"[^0-9A-Za-z]+", "", run_timestamp)
    design_config_sha256 = (
        sha256_file(design_config_path)
        if design_config_path is not None and design_config_path.is_file()
        else sha256_json(design_config)
    )
    run_metadata = {
        "run_id": run_id,
        "generated_at": run_timestamp,
        "nfl_ab_design_version": PACKAGE_VERSION,
        "design_campaign_sha256": design_config_sha256,
        "workflow_source_sha256": sha256_file(SCRIPT_PATH),
        "design_pipeline_source_sha256": sha256_file(SCRIPT_PATH.with_name("design_pipeline.py")),
    }

    cleavage_rows, medium_high_sites, inferred_fragments = infer_antigen_truncation(full_sequence, truncation_constraints)
    prioritized_fragments = prioritize_antigen_fragments(inferred_fragments, truncation_constraints)
    primary_fragment = prioritized_fragments[0]["fragment"]
    modeling_fragment = design_context_fragment(primary_fragment, design_config, full_sequence)
    epitope_rows = build_epitope_windows(full_sequence, modeling_fragment)
    target_epitope_rows = resolve_campaign_epitopes(design_config, epitope_rows, full_sequence)
    validate_design_campaign_contract(design_config, target_epitope_rows, validated)
    target_epitope_ids = _configured_epitope_ids(design_config)
    target_id_set = set(target_epitope_ids)
    for row in epitope_rows:
        row["configured_design_target"] = str(row["epitope_id"]) in target_id_set
    for row in target_epitope_rows:
        if str(row["epitope_id"]) not in {str(item["epitope_id"]) for item in epitope_rows}:
            epitope_rows.append({**row, "configured_design_target": True})
    thresholds = design_config.get("stage_thresholds", {})
    design_run = run_design_pipeline(
        validated,
        target_epitope_rows,
        seed=_simulation_setting(design_config, "seed", 20260812),
        designs_per_template_epitope=_simulation_setting(design_config, "designs_per_template_epitope", 24),
        selection_count=int(thresholds.get("selection_count", 12)),
        thresholds=thresholds,
        template_specs=_configured_template_specs(design_config),
        cdr_ranges_by_template=_configured_cdr_ranges(design_config),
        epitope_ids=target_epitope_ids,
    )
    design_result = design_run.as_dict()
    prospective_rows = design_result["prospective_ranking_rows"]
    retrospective_rows = design_result["retrospective_ranking_rows"]
    pair_rows = sandwich_pair_ranking(retrospective_rows)
    sandwich = pair_rows[0] if pair_rows else {"status": "not_enough_candidates"}
    export_count = int(thresholds.get("selection_count", 12))
    selected_prospective_rows = [row for row in prospective_rows if bool(row.get("selected_for_export"))]
    exported_candidates = [_row_to_antibody(row) for row in selected_prospective_rows[:export_count]]
    ensure_clean_output_dir()
    manifest = export_structure_inputs(
        exported_candidates,
        prioritized_fragments,
        full_sequence,
        sandwich,
        external_config_path,
        template_antibodies=validated,
        prepared_template_rows=design_result["template_rows"],
        sandwich_antibodies=validated,
        epitope_rows=target_epitope_rows,
        design_config=design_config,
        run_metadata=run_metadata,
        modeling_fragment=modeling_fragment,
    )

    write_csv(OUTPUT_DIR / "00_antigen_truncation_all_peptide_bonds.csv", cleavage_rows)
    write_csv(OUTPUT_DIR / "00_antigen_truncation_medium_high_sites.csv", medium_high_sites)
    write_csv(OUTPUT_DIR / "00_antigen_truncation_fragment_candidates.csv", inferred_fragments)
    write_truncation_report(cleavage_rows, medium_high_sites, inferred_fragments, OUTPUT_DIR / "00_antigen_truncation_report.md")
    write_csv(OUTPUT_DIR / "01_antigen_fragment_prioritization.csv", prioritized_fragments)
    write_csv(OUTPUT_DIR / "02_epitope_windows.csv", epitope_rows)
    write_csv(OUTPUT_DIR / "03_template_frameworks.csv", design_result["template_rows"])
    backbone_rows = [
        {
            "candidate_id": row["candidate_id"],
            "template_id": row["template_id"],
            "framework_source_antibody_id": row["framework_source_antibody_id"],
            "target_epitope_id": row["target_epitope_id"],
            "design_regions": row["design_regions"],
            "generation_stage": "simulated_backbone_and_six_CDR_proposal",
            "requested_real_engines": "RFantibody;IgGM;Germinal",
            "real_engine_execution_state": "not_run",
            "data_status": "simulated",
            "metric_provenance": row["metric_provenance"],
        }
        for row in design_result["generation_rows"]
    ]
    write_csv(OUTPUT_DIR / "04_backbone_generation.csv", backbone_rows)
    write_csv(OUTPUT_DIR / "05_sequence_candidates.csv", design_result["generation_rows"])
    write_csv(OUTPUT_DIR / "06_structure_interface_screen.csv", design_result["structure_rows"])
    write_csv(OUTPUT_DIR / "07_developability_screen.csv", design_result["developability_rows"])
    write_csv(OUTPUT_DIR / "08_screening_funnel.csv", design_result["funnel_rows"])
    write_csv(OUTPUT_DIR / "09_prospective_candidates.csv", design_result["prospective_ranking_rows"])
    write_csv(OUTPUT_DIR / "10_retrospective_demo_candidates.csv", retrospective_rows)
    write_csv(OUTPUT_DIR / "11_sandwich_pair_ranking.csv", pair_rows)
    write_sandwich_report(sandwich, OUTPUT_DIR / "11_sandwich_pair_report.md")

    source_manifest = {
        "schema": "nfl_ab_design.run_manifest.v2",
        **run_metadata,
        "run_timestamp": run_timestamp,
        "run_mode": "deterministic_proxy_simulation_with_retrospective_positive_controls",
        "real_model_execution": False,
        "story_path": relative_to_package(STORY_PATH),
        "research_plan_path": relative_to_package(RESEARCH_PLAN_PATH),
        "antigen_report_path": relative_to_package(ANTIGEN_REPORT_PATH),
        "cleavage_sites_path": relative_to_package(CLEAVAGE_SITES_PATH),
        "genpept_path": relative_to_package(GENPEPT_PATH),
        "truncation_constraints_path": relative_to_package(TRUNCATION_CONSTRAINTS_PATH),
        "design_campaign_config_path": relative_to_package(design_config_path) if design_config_path else "",
        "design_campaign_config": design_config,
        "design_campaign_config_sha256": design_config_sha256,
        "validation_antibody_fasta_path": relative_to_package(ANTIBODY_FASTA_PATH),
        "validation_sequence_usage": "framework_only_during_generation;full_sequence_only_after_prospective_ranking",
        "nfl_sequence_length": len(full_sequence),
        "upstream_cleavage_site_rows": len(cleavage_sites),
        "computed_peptide_bond_rows": len(cleavage_rows),
        "computed_medium_high_site_rows": len(medium_high_sites),
        "computed_fragment_candidate_rows": len(inferred_fragments),
        "generation_candidate_count": len(design_result["generation_rows"]),
        "biochemical_lead_fragment": primary_fragment,
        "modeling_context_fragment": modeling_fragment,
        "prospective_survivor_count": len(design_result["prospective_ranking_rows"]),
        "prospective_export_count": len(exported_candidates),
        "retrospective_control_count": sum(
            row["control_status"] == "retrospective_positive_control" for row in retrospective_rows
        ),
        "story_characters": len(story),
        "research_plan_characters": len(research_plan),
        "antigen_report_characters": len(antigen_report),
        "design_request_index_sha256": manifest["exports"].get("design_request_index_sha256", ""),
        "external_tool_manifest_sha256": sha256_file(EXPORT_DIR / "external_tool_manifest.json"),
    }
    write_json(OUTPUT_DIR / "intermediate" / "source_manifest.json", source_manifest)
    write_json(OUTPUT_DIR / "intermediate" / "run_manifest.json", source_manifest)
    write_workflow_report(
        cleavage_rows=cleavage_rows,
        medium_high_sites=medium_high_sites,
        inferred_fragments=inferred_fragments,
        epitope_rows=epitope_rows,
        design_result=design_result,
        pair_rows=pair_rows,
        sandwich=sandwich,
        manifest=manifest,
        run_timestamp=run_timestamp,
        biochemical_lead_fragment=primary_fragment,
        modeling_context_fragment=modeling_fragment,
    )
    return {
        "output_dir": OUTPUT_DIR,
        "primary_fragment": primary_fragment,
        "modeling_fragment": modeling_fragment,
        "ranking_rows": prospective_rows,
        "ranking_rows_scope": "prospective_simulation",
        "prospective_ranking_rows": design_result["prospective_ranking_rows"],
        "retrospective_ranking_rows": retrospective_rows,
        "design_result": design_result,
        "target_epitope_rows": target_epitope_rows,
        "pair_rows": pair_rows,
        "sandwich": sandwich,
        "manifest": manifest,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nfl-ab-design",
        description="Run the NfL epitope-conditioned de novo design workflow and prepare real-model handoff inputs.",
    )
    parser.add_argument(
        "--external-config",
        type=Path,
        default=EXTERNAL_PIPELINE_CONFIG_PATH,
        help="JSON file describing external structure/docking pipeline adapters.",
    )
    parser.add_argument(
        "--design-config",
        type=Path,
        default=DESIGN_CAMPAIGN_CONFIG_PATH,
        help="JSON campaign defining templates, all-six-CDR coordinates, target epitopes, simulation scale, and thresholds.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    result = run_workflow(
        external_config_path=args.external_config,
        design_config_path=args.design_config,
    )
    prospective_rows = result["prospective_ranking_rows"]
    retrospective_controls = [
        row
        for row in result["retrospective_ranking_rows"]
        if row.get("control_status") == "retrospective_positive_control"
    ]

    print(f"NfL antibody workflow complete. Outputs written to: {result['output_dir']}")
    print(f"Primary antigen fragment: NEFL {result['primary_fragment']}")
    print(f"Modeling context covering all target epitopes: NEFL {result['modeling_fragment']}")
    print("Prospective simulated candidates (not real-model results):")
    for row in prospective_rows[:5]:
        print(
            f"Rank {row['rank']:>2}: {row['candidate_id']} | "
            f"{row['best_epitope_id']} | total={row['total_rank_score']}"
        )
    print("Retrospective positive-control demonstration (not blind discovery):")
    for row in retrospective_controls:
        print(
            f"Demo rank {row['rank']:>2}: {row['candidate_id']} | "
            f"status={row['control_status']} | total={row['total_rank_score']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
