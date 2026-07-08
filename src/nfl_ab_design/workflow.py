"""Reproducible NfL antibody design workflow.

This module converts the repository project context, antigen-inference
resources, design constraints, and validation antibodies into an executable
computational replay pipeline. It intentionally uses only the Python standard
library so it can run in a clean workspace without structure-prediction
dependencies.

The scoring is a deterministic proxy. It is meant to organize the work, create
auditable intermediate tables, and export input templates for external structure
tools such as IgFold/ABodyBuilder3, AlphaFold3, Chai-1, Boltz, and Rosetta.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import shutil
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


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
OUTPUT_DIR = PACKAGE_ROOT / "outputs"
EXPORT_DIR = OUTPUT_DIR / "exports"


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
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
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
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


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
    if "boundary" in label or "Cys322" in label:
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
    if start <= 322 <= end:
        notes.append("contains Cys322 disulfide-anchor region")
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
        ("Cys322_anchor_316_331", 316, 331),
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
    if config_path is None or not config_path.exists():
        return {"pipelines": []}
    return json.loads(config_path.read_text(encoding="utf-8"))


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
        package_root=str(PACKAGE_ROOT),
    )


def prepare_external_pipeline_handoff(exported: dict[str, Any], config_path: Path | None) -> dict[str, Any]:
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
        enabled = bool(item.get("enabled", False))
        for input_path in select_external_inputs(exported, selector):
            result_dir = OUTPUT_DIR / "external_results" / name / Path(input_path).stem
            output_dir = relative_to_package(result_dir)
            command = render_external_command(command_template, input_path, output_dir, name) if command_template else ""
            jobs.append(
                {
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
            fieldnames=["stage", "tool", "enabled", "input_selector", "input_path", "output_dir", "command"],
            delimiter="\t",
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
    validated: list[Antibody],
    prioritized_fragments: list[dict[str, Any]],
    full_sequence: str,
    sandwich: dict[str, Any],
    external_config_path: Path | None,
) -> dict[str, Any]:
    fasta_dir = EXPORT_DIR / "fasta"
    af3_dir = EXPORT_DIR / "af3_json"
    top_fragments = prioritized_fragments[:3]
    antigen_records: list[tuple[str, str]] = []
    for row in top_fragments:
        start, end = parse_range(row["fragment"])
        antigen_records.append((f"NEFL_{row['fragment']}_P07196", subseq(full_sequence, start, end)))
    antigen_fasta_path = fasta_dir / "antigen_fragments.fasta"
    write_fasta(antigen_fasta_path, antigen_records)

    antibody_records: list[tuple[str, str]] = []
    for antibody in validated:
        antibody_records.append((f"{antibody.antibody_id}|VH|{antibody.vh_id}", antibody.vh))
        antibody_records.append((f"{antibody.antibody_id}|VL|{antibody.vl_id}", antibody.vl))
    validated_fv_path = fasta_dir / "validated_fv_chains.fasta"
    write_fasta(validated_fv_path, antibody_records)

    primary_fragment = top_fragments[0]["fragment"]
    primary_start, primary_end = parse_range(primary_fragment)
    primary_antigen = subseq(full_sequence, primary_start, primary_end)
    exported: dict[str, Any] = {
        "primary_antigen_fragment": primary_fragment,
        "antigen_fragments_fasta": relative_to_package(antigen_fasta_path),
        "validated_fv_chains_fasta": relative_to_package(validated_fv_path),
        "complex_fastas": [],
        "sandwich_fasta": "",
        "fasta_files": [],
        "af3_json_files": [],
    }

    for antibody in validated:
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

    if len(validated) >= 2:
        ab1, ab2 = validated[:2]
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
        exported["fasta_files"].append(relative_to_package(fasta_path))
        exported["af3_json_files"].append(relative_to_package(json_path))

    external_handoff = prepare_external_pipeline_handoff(exported, external_config_path)

    manifest = {
        "purpose": "Structure-tool input handoff for the NfL antibody design workflow.",
        "limitations": [
            "VH/VL Fv chains are exported without constant regions.",
            "AF3 JSON files are schema templates and should be checked against the active runner.",
            "Proxy ranking metrics should be replaced by measured or modeled ipTM, pTM, interface PAE, pDockQ, buried surface area, Rosetta interface dG, and clash metrics when structures are available.",
        ],
        "recommended_tool_order": [
            "IgFold or ABodyBuilder3 for Fv/Fab sanity checks",
            "AF3, Chai-1, or Boltz co-folding for antibody-antigen complexes",
            "Rosetta relax/interface analyzer for post-prediction interface metrics",
            "Pair-aware trimer prediction for sandwich compatibility",
        ],
        "sandwich_pair_proxy": sandwich,
        "exports": exported,
        "external_pipeline_handoff": external_handoff,
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


def write_workflow_report(
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


def run_workflow(external_config_path: Path | None = EXTERNAL_PIPELINE_CONFIG_PATH) -> dict[str, Any]:
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


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nfl-ab-design",
        description="Run the NfL antibody design replay workflow and prepare external structure-pipeline inputs.",
    )
    parser.add_argument(
        "--external-config",
        type=Path,
        default=EXTERNAL_PIPELINE_CONFIG_PATH,
        help="JSON file describing external structure/docking pipeline adapters.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    result = run_workflow(external_config_path=args.external_config)
    ranking_rows = result["ranking_rows"]

    print(f"NfL antibody workflow complete. Outputs written to: {result['output_dir']}")
    print(f"Primary antigen fragment: NEFL {result['primary_fragment']}")
    for row in ranking_rows[:5]:
        print(
            f"Rank {row['rank']:>2}: {row['candidate_id']} | "
            f"{row['best_epitope_id']} | total={row['total_rank_score']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
