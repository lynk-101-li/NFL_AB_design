#!/usr/bin/env python3
"""Build audited RFantibody HLT and Germinal scFv template coordinates.

The production path is deliberately local and deterministic: it consumes
already-generated coordinate files, never invokes a model, subprocess, or
network service, and leaves every generated manifest blocked for human review.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import tempfile
from typing import Any, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_SCHEMA = "nfl_ab_design.antibody_template_evidence.v1"
MANIFEST_FRAGMENT_SCHEMA = "nfl_ab_design.antibody_template_manifest_fragment.v1"
UPSTREAM_PROVENANCE_SCHEMA = "nfl_ab_design.template_coordinate_generation_provenance.v1"
AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
AA20 = frozenset(AA3_TO_AA1.values())
BACKBONE = frozenset({"N", "CA", "C", "O"})
EXPECTED_HEAVY_ATOMS = {
    "A": BACKBONE | {"CB"},
    "C": BACKBONE | {"CB", "SG"},
    "D": BACKBONE | {"CB", "CG", "OD1", "OD2"},
    "E": BACKBONE | {"CB", "CG", "CD", "OE1", "OE2"},
    "F": BACKBONE | {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "G": BACKBONE,
    "H": BACKBONE | {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "I": BACKBONE | {"CB", "CG1", "CG2", "CD1"},
    "K": BACKBONE | {"CB", "CG", "CD", "CE", "NZ"},
    "L": BACKBONE | {"CB", "CG", "CD1", "CD2"},
    "M": BACKBONE | {"CB", "CG", "SD", "CE"},
    "N": BACKBONE | {"CB", "CG", "OD1", "ND2"},
    "P": BACKBONE | {"CB", "CG", "CD"},
    "Q": BACKBONE | {"CB", "CG", "CD", "OE1", "NE2"},
    "R": BACKBONE | {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "S": BACKBONE | {"CB", "OG"},
    "T": BACKBONE | {"CB", "OG1", "CG2"},
    "V": BACKBONE | {"CB", "CG1", "CG2"},
    "W": BACKBONE | {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "Y": BACKBONE | {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
}
STANDARD_GGGGS3 = "GGGGSGGGGSGGGGS"
# Germinal's pinned official ``pdbs/scfv.pdb`` yields this 15-residue interval
# under the adapter's explicit [first len(VH)]/[last len(VL)] slicing contract.
LINKER = "AGGGGSGGGGSGGGS"
GENERIC_LINKER_OBSERVED = LINKER
CDR_ORDER = ("H1", "H2", "H3", "L1", "L2", "L3")
TEMPLATE_IDS = ("template_7-H11-D3-2-C7", "template_15-C12-H6")

# These are neutral, deterministic coordinate seeds, not candidate antibodies.
# Only the CDR positions proved by the pinned ANARCI tables are replaced.
NEUTRAL_CDR_SEEDS: dict[str, dict[str, str]] = {
    "template_7-H11-D3-2-C7": {
        "H1": "ASGTQSS", "H2": "GSTASQ", "H3": "GST",
        "L1": "QASSSGTSTLA", "L2": "GASTSAS", "L3": "QQSGTSPRT",
    },
    "template_15-C12-H6": {
        "H1": "GASSTQSA", "H2": "STAGQ", "H3": "SGTQADSGTAY",
        "L1": "RASQSGTSTLA", "L2": "SASTGAS", "L3": "QQAGTSPAT",
    },
}


class AntibodyTemplatePreparationError(ValueError):
    """Raised when an input cannot satisfy the audited coordinate contract."""


@dataclass(frozen=True)
class Atom:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    number: int
    icode: str
    x: float
    y: float
    z: float
    occupancy: float
    bfactor: float
    element: str

    @property
    def xyz(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z


@dataclass(frozen=True)
class Residue:
    chain: str
    number: int
    icode: str
    amino_acid: str
    atoms: tuple[Atom, ...]

    @property
    def label(self) -> str:
        return f"{self.number}{self.icode}"

    @property
    def by_name(self) -> dict[str, Atom]:
        return {atom.name: atom for atom in self.atoms}


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _portable(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _json_bytes(value: Mapping[str, Any]) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode()


def _read_fasta(path: Path) -> dict[tuple[str, str], str]:
    records: dict[tuple[str, str], str] = {}
    header: str | None = None
    sequence: list[str] = []

    def commit() -> None:
        if header is None:
            return
        fields = header.split("|")
        if len(fields) < 2 or fields[1] not in {"VH", "VL"}:
            raise AntibodyTemplatePreparationError(f"Invalid VH/VL FASTA header: {header}")
        value = "".join(sequence).replace(" ", "").upper()
        invalid = sorted(set(value) - AA20)
        if not value or invalid:
            raise AntibodyTemplatePreparationError(
                f"Invalid sequence for {header}; unsupported residues={invalid}"
            )
        key = (fields[0], fields[1])
        if key in records:
            raise AntibodyTemplatePreparationError(f"Duplicate FASTA record: {key}")
        records[key] = value

    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            commit()
            header, sequence = line[1:], []
        else:
            if header is None:
                raise AntibodyTemplatePreparationError("FASTA sequence precedes header")
            sequence.append(line)
    commit()
    return records


def _is_chothia_label(value: str) -> bool:
    return bool(re.fullmatch(r"\d+[A-Za-z]?", value))


def _read_anarci_tables(
    heavy_csv: Path, light_csv: Path
) -> dict[tuple[str, str], list[tuple[str, str]]]:
    result: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for path in (heavy_csv, light_csv):
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                fields = str(row.get("Id", "")).split("|")
                if len(fields) < 2 or fields[1] not in {"VH", "VL"}:
                    raise AntibodyTemplatePreparationError(f"Invalid ANARCI Id in {path}")
                mapped = [
                    (label, value.upper())
                    for label, value in row.items()
                    if _is_chothia_label(str(label)) and value not in {None, "", "-"}
                ]
                key = (fields[0], fields[1])
                if key in result or not mapped:
                    raise AntibodyTemplatePreparationError(f"Duplicate/empty ANARCI row: {key}")
                result[key] = mapped
    return result


def _read_neutral_anarci_tables(
    heavy_csv: Path, light_csv: Path
) -> dict[tuple[str, str], list[tuple[str, str]]]:
    expected_ids = {
        f"{template_id}_H": (template_id, "VH") for template_id in TEMPLATE_IDS
    } | {
        f"{template_id}_L": (template_id, "VL") for template_id in TEMPLATE_IDS
    }
    result: dict[tuple[str, str], list[tuple[str, str]]] = {}
    observed_ids: set[str] = set()
    for path, expected_chain in ((heavy_csv, "VH"), (light_csv, "VL")):
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                identifier = str(row.get("Id", ""))
                if identifier not in expected_ids:
                    raise AntibodyTemplatePreparationError(
                        f"Invalid neutral ANARCI Id {identifier!r}; expected one of {sorted(expected_ids)}"
                    )
                key = expected_ids[identifier]
                if key[1] != expected_chain:
                    raise AntibodyTemplatePreparationError(
                        f"Neutral ANARCI Id {identifier!r} is in the wrong chain table"
                    )
                mapped = [
                    (label, value.upper())
                    for label, value in row.items()
                    if _is_chothia_label(str(label)) and value not in {None, "", "-"}
                ]
                if identifier in observed_ids or not mapped:
                    raise AntibodyTemplatePreparationError(
                        f"Duplicate/empty neutral ANARCI row: {identifier}"
                    )
                observed_ids.add(identifier)
                result[key] = mapped
    if observed_ids != set(expected_ids):
        raise AntibodyTemplatePreparationError(
            f"Neutral ANARCI tables must contain exactly {sorted(expected_ids)}; "
            f"found {sorted(observed_ids)}"
        )
    return result


def _read_numbering_evidence(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AntibodyTemplatePreparationError(f"Cannot read numbering evidence: {exc}") from exc
    if value.get("annotation_method") != "ANARCI 2020.04.23 Chothia":
        raise AntibodyTemplatePreparationError("Numbering evidence is not pinned ANARCI Chothia")
    return value


def _read_and_bind_upstream_provenance(
    path: Path,
    *, paired_paths: Mapping[str, Path], prepared_scfv_paths: Mapping[str, Path],
    generic_path: Path, broad_paths: Mapping[str, Path],
) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AntibodyTemplatePreparationError(f"Cannot read upstream provenance: {exc}") from exc
    if not isinstance(value, dict) or value.get("schema") != UPSTREAM_PROVENANCE_SCHEMA:
        raise AntibodyTemplatePreparationError(
            f"Upstream provenance must be a JSON object with schema {UPSTREAM_PROVENANCE_SCHEMA!r}"
        )
    hashes = value.get("input_artifact_sha256")
    if not isinstance(hashes, Mapping):
        raise AntibodyTemplatePreparationError("Upstream provenance lacks input_artifact_sha256")
    expected: dict[str, Any] = {
        "paired_fv_pdbs": {template_id: _sha256(paired_paths[template_id]) for template_id in TEMPLATE_IDS},
        "prepared_scfv_pdbs": {
            template_id: _sha256(prepared_scfv_paths[template_id]) for template_id in TEMPLATE_IDS
        } if prepared_scfv_paths else {},
        "generic_scfv_pdb": _sha256(generic_path),
        "official_broad_hlt_pdbs": {
            template_id: _sha256(broad_paths[template_id]) for template_id in TEMPLATE_IDS
        } if broad_paths else {},
    }
    normalized = {key: hashes.get(key) for key in expected}
    if normalized != expected or set(hashes) != set(expected):
        raise AntibodyTemplatePreparationError(
            "Upstream provenance input_artifact_sha256 does not exactly bind CLI inputs"
        )
    return {
        "path": _portable(path),
        "sha256": _sha256(path),
        "schema": UPSTREAM_PROVENANCE_SCHEMA,
        "input_hashes_verified": True,
        "model_claims_trusted_by_preparation_cli": False,
    }


def _parse_atom_line(line: str, line_number: int) -> Atom:
    if len(line) < 54:
        raise AntibodyTemplatePreparationError(f"Truncated PDB atom line {line_number}")
    try:
        return Atom(
            record=line[:6].strip(), serial=int(line[6:11]), name=line[12:16].strip(),
            altloc=line[16:17].strip(), resname=line[17:20].strip().upper(),
            chain=line[21:22].strip(), number=int(line[22:26]),
            icode=line[26:27].strip().upper(), x=float(line[30:38]),
            y=float(line[38:46]), z=float(line[46:54]),
            occupancy=float(line[54:60] or 0), bfactor=float(line[60:66] or 0),
            element=(line[76:78].strip().upper() if len(line) >= 78 else line[12:14].strip()),
        )
    except ValueError as exc:
        raise AntibodyTemplatePreparationError(f"Malformed PDB atom line {line_number}") from exc


def _parse_pdb(path: Path, *, allowed_chains: set[str] | None = None) -> dict[str, list[Residue]]:
    if not path.is_file():
        raise AntibodyTemplatePreparationError(f"PDB does not exist: {path}")
    atom_groups: dict[tuple[str, int, str], list[Atom]] = {}
    order: list[tuple[str, int, str]] = []
    model_count = 0
    for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if line.startswith("MODEL "):
            model_count += 1
            if model_count > 1:
                raise AntibodyTemplatePreparationError(f"Multiple MODEL records in {path}")
        if not line.startswith("ATOM  "):
            continue
        atom = _parse_atom_line(line, line_number)
        if not atom.chain:
            raise AntibodyTemplatePreparationError(f"Blank chain ID in {path}:{line_number}")
        if atom.altloc not in {"", "A"}:
            raise AntibodyTemplatePreparationError(f"Unsupported altloc in {path}:{line_number}")
        if atom.resname not in AA3_TO_AA1:
            raise AntibodyTemplatePreparationError(f"Noncanonical residue {atom.resname} in {path}")
        key = (atom.chain, atom.number, atom.icode)
        if key not in atom_groups:
            atom_groups[key] = []
            order.append(key)
        if any(existing.name == atom.name for existing in atom_groups[key]):
            raise AntibodyTemplatePreparationError(f"Duplicate atom {atom.name} at {key} in {path}")
        atom_groups[key].append(atom)
    if not order:
        raise AntibodyTemplatePreparationError(f"No ATOM records in {path}")
    chains: dict[str, list[Residue]] = {}
    seen: set[tuple[str, int, str]] = set()
    for key in order:
        if key in seen:
            raise AntibodyTemplatePreparationError(f"Non-contiguous duplicate residue {key} in {path}")
        seen.add(key)
        atoms = tuple(atom_groups[key])
        residue_names = {a.resname for a in atoms}
        if len(residue_names) != 1:
            raise AntibodyTemplatePreparationError(f"Mixed residue names at {key}")
        amino_acid = AA3_TO_AA1[atoms[0].resname]
        names = {a.name for a in atoms if not a.element.startswith("H")}
        missing_atoms = EXPECTED_HEAVY_ATOMS[amino_acid] - names
        if missing_atoms:
            raise AntibodyTemplatePreparationError(
                f"Residue {key} lacks canonical heavy atom(s): {sorted(missing_atoms)}"
            )
        if any(
            not all(math.isfinite(value) for value in (*atom.xyz, atom.occupancy, atom.bfactor))
            for atom in atoms
        ):
            raise AntibodyTemplatePreparationError(f"Non-finite atom value at {key}")
        chains.setdefault(key[0], []).append(
            Residue(key[0], key[1], key[2], amino_acid, atoms)
        )
    if allowed_chains is not None and set(chains) != allowed_chains:
        raise AntibodyTemplatePreparationError(
            f"{path} must contain exactly chains {sorted(allowed_chains)}, found {sorted(chains)}"
        )
    for chain, residues in chains.items():
        labels = [r.label for r in residues]
        if len(labels) != len(set(labels)):
            raise AntibodyTemplatePreparationError(f"Duplicate PDB residue labels on chain {chain}")
    return chains


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))


def _peptide_geometry(residues: Sequence[Residue]) -> dict[str, Any]:
    values: list[float] = []
    outliers: list[dict[str, Any]] = []
    for left, right in zip(residues, residues[1:]):
        value = _distance(left.by_name["C"].xyz, right.by_name["N"].xyz)
        if not 1.15 <= value <= 1.55:
            raise AntibodyTemplatePreparationError(
                f"Peptide C-N distance outside strict 1.15-1.55 A window: {value:.3f} A at "
                f"{left.chain}{left.label}->{right.chain}{right.label}"
            )
        values.append(value)
    return {
        "bond_count": len(values),
        "minimum_cn_angstrom": round(min(values), 6) if values else None,
        "maximum_cn_angstrom": round(max(values), 6) if values else None,
        "ideal_window_angstrom": [1.15, 1.55],
        "hard_reject_window_angstrom": [1.15, 1.55],
        "outside_ideal_count": len(outliers),
        "outside_ideal": outliers,
        "requires_human_review": bool(outliers),
    }


def _assert_no_severe_clashes(
    chains: Mapping[str, Sequence[Residue]], *, threshold: float = 1.8
) -> dict[str, Any]:
    """Reject nonbonded heavy atoms closer than a conservative hard cutoff."""
    cell_size = threshold
    cells: dict[tuple[int, int, int], list[tuple[str, int, str, str, tuple[float, float, float]]]] = {}
    atom_count = 0
    minimum = math.inf
    for chain, residues in chains.items():
        for sequence_index, residue in enumerate(residues):
            for atom in residue.atoms:
                if atom.element.startswith("H"):
                    continue
                atom_count += 1
                cell = tuple(math.floor(value / cell_size) for value in atom.xyz)
                for dx in (-1, 0, 1):
                    for dy in (-1, 0, 1):
                        for dz in (-1, 0, 1):
                            for other in cells.get((cell[0] + dx, cell[1] + dy, cell[2] + dz), []):
                                other_chain, other_index, other_label, other_atom, other_xyz = other
                                if chain == other_chain and abs(sequence_index - other_index) <= 1:
                                    continue
                                distance = _distance(atom.xyz, other_xyz)
                                minimum = min(minimum, distance)
                                if distance < threshold:
                                    raise AntibodyTemplatePreparationError(
                                        f"Severe nonbonded heavy-atom clash {distance:.3f} A: "
                                        f"{other_chain}{other_label}/{other_atom} vs "
                                        f"{chain}{residue.label}/{atom.name}"
                                    )
                cells.setdefault(cell, []).append(
                    (chain, sequence_index, residue.label, atom.name, atom.xyz)
                )
    return {
        "heavy_atom_count": atom_count,
        "hard_clash_cutoff_angstrom": threshold,
        "severe_clash_count": 0,
        "minimum_checked_nonbonded_distance_angstrom": (
            round(minimum, 6) if math.isfinite(minimum) else None
        ),
    }


def _replace_cdrs(
    source: str, cdrs: Mapping[str, Mapping[str, Any]], seeds: Mapping[str, str]
) -> str:
    output = list(source)
    for name, record in cdrs.items():
        start = int(record["raw_start_1_based"])
        end = int(record["raw_end_1_based_inclusive"])
        seed = seeds[name]
        if len(seed) != end - start + 1:
            raise AntibodyTemplatePreparationError(f"Neutral seed length mismatch for {name}")
        if source[start - 1:end] != record["raw_sequence"]:
            raise AntibodyTemplatePreparationError(f"Known CDR evidence mismatch for {name}")
        if seed == record["raw_sequence"]:
            raise AntibodyTemplatePreparationError(f"Neutral seed equals known CDR for {name}")
        output[start - 1:end] = seed
    return "".join(output)


def _validate_contracts(
    fasta: dict[tuple[str, str], str],
    tables: dict[tuple[str, str], list[tuple[str, str]]],
    evidence: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    templates: dict[str, dict[str, Any]] = {}
    evidence_templates = evidence.get("templates")
    if not isinstance(evidence_templates, Mapping):
        raise AntibodyTemplatePreparationError("Numbering evidence lacks templates")
    for template_id in TEMPLATE_IDS:
        if template_id not in evidence_templates:
            raise AntibodyTemplatePreparationError(f"Missing evidence for {template_id}")
        template_evidence = evidence_templates[template_id]
        neutral: dict[str, str] = {}
        labels: dict[str, list[str]] = {}
        raw_ranges: dict[str, list[int]] = {}
        for chain_name in ("VH", "VL"):
            key = (template_id.removeprefix("template_"), chain_name)
            source = fasta.get(key)
            table = tables.get(key)
            if source is None or table is None:
                raise AntibodyTemplatePreparationError(f"Missing FASTA/ANARCI record {key}")
            table_sequence = "".join(aa for _, aa in table)
            if table_sequence != source:
                raise AntibodyTemplatePreparationError(f"ANARCI sequence mismatch for {key}")
            chain_evidence = template_evidence["chains"][chain_name]
            expected_hash = chain_evidence["sequence_sha256"]
            if sha256(source.encode()).hexdigest() != expected_hash:
                raise AntibodyTemplatePreparationError(f"Sequence hash mismatch for {key}")
            chain_cdrs = chain_evidence["cdrs"]
            neutral[chain_name] = _replace_cdrs(source, chain_cdrs, NEUTRAL_CDR_SEEDS[template_id])
            for name, cdr in chain_cdrs.items():
                start, end = int(cdr["raw_start_1_based"]), int(cdr["raw_end_1_based_inclusive"])
                actual_labels = [label for label, _ in table[start - 1:end]]
                if actual_labels != list(cdr["chothia_labels"]):
                    raise AntibodyTemplatePreparationError(f"ANARCI CDR labels mismatch for {template_id}/{name}")
                labels[name] = actual_labels
                raw_ranges[name] = [start, end]
        templates[template_id] = {
            "neutral_vh": neutral["VH"], "neutral_vl": neutral["VL"],
            "cdr_chothia_labels": labels, "cdr_raw_ranges": raw_ranges,
            "anarci_by_chain": {
                chain: tables[(template_id.removeprefix("template_"), chain)]
                for chain in ("VH", "VL")
            },
        }
    return templates


def _pose_indices(chains: Mapping[str, Sequence[Residue]]) -> dict[tuple[str, str], int]:
    result: dict[tuple[str, str], int] = {}
    index = 0
    for chain in ("H", "L"):
        for residue in chains[chain]:
            index += 1
            result[(chain, residue.label)] = index
    return result


def _format_atom(atom: Atom, *, serial: int, chain: str, number: int, icode: str,
                 xyz: tuple[float, float, float] | None = None,
                 resname: str | None = None) -> str:
    x, y, z = xyz or atom.xyz
    element = (atom.element or atom.name[0]).rjust(2)
    return (
        f"ATOM  {serial:5d} {atom.name:>4s}{atom.altloc:1s}{(resname or atom.resname):>3s} "
        f"{chain:1s}{number:4d}{icode:1s}   {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {element:>2s}  "
    )


def _write_hlt(
    path: Path, chains: Mapping[str, Sequence[Residue]], cdr_labels: Mapping[str, Sequence[str]]
) -> tuple[bytes, dict[str, list[int]]]:
    pose = _pose_indices(chains)
    absolute: dict[str, list[int]] = {}
    remarks: list[str] = []
    for cdr in CDR_ORDER:
        chain = cdr[0]
        values: list[int] = []
        for label in cdr_labels[cdr]:
            key = (chain, label)
            if key not in pose:
                raise AntibodyTemplatePreparationError(f"PDB lacks Chothia residue {chain}{label} for {cdr}")
            values.append(pose[key])
            remarks.append(f"REMARK PDBinfo-LABEL: {pose[key]:4d} {cdr}")
        absolute[cdr] = values
    lines = remarks[:]
    serial = 1
    for chain in ("H", "L"):
        for residue in chains[chain]:
            for atom in residue.atoms:
                lines.append(_format_atom(atom, serial=serial, chain=chain,
                                          number=residue.number, icode=residue.icode))
                serial += 1
        lines.append("TER")
    lines.append("END")
    return ("\n".join(lines) + "\n").encode("ascii"), absolute


def _kabsch(mobile: Sequence[tuple[float, float, float]], reference: Sequence[tuple[float, float, float]]):
    try:
        import numpy as np
    except ImportError as exc:
        raise AntibodyTemplatePreparationError("NumPy is required for rigid-body alignment") from exc
    if len(mobile) != len(reference) or len(mobile) < 3:
        raise AntibodyTemplatePreparationError("Alignment requires at least three paired CA atoms")
    mob, ref = np.asarray(mobile, float), np.asarray(reference, float)
    mob_center, ref_center = mob.mean(axis=0), ref.mean(axis=0)
    u, _, vt = np.linalg.svd((mob - mob_center).T @ (ref - ref_center))
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    translation = ref_center - mob_center @ rotation
    fitted = mob @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((fitted - ref) ** 2, axis=1))))
    return rotation, translation, rmsd


def _apply_transform(xyz, rotation, translation) -> tuple[float, float, float]:
    import numpy as np
    value = np.asarray(xyz, float) @ rotation + translation
    return tuple(float(x) for x in value)


def _framework_indices(length: int, cdr_raw_ranges: Mapping[str, Sequence[int]], prefix: str) -> list[int]:
    excluded: set[int] = set()
    for name, (start, end) in cdr_raw_ranges.items():
        if name.startswith(prefix):
            excluded.update(range(int(start), int(end) + 1))
    return [index for index in range(1, length + 1) if index not in excluded]


def _domain_boundaries(
    generic: Sequence[Residue], _vh_len: int, _vl_len: int
) -> tuple[list[Residue], list[Residue], list[Residue]]:
    sequence = "".join(residue.amino_acid for residue in generic)
    occurrences = [m.start() for m in re.finditer(f"(?={LINKER})", sequence)]
    if len(occurrences) != 1:
        raise AntibodyTemplatePreparationError(
            "Generic scFv must contain exactly one locked upstream domain-junction segment; "
            f"found {len(occurrences)}"
        )
    start = occurrences[0]
    end = start + len(LINKER)
    if start < 70 or len(generic) - end < 70:
        raise AntibodyTemplatePreparationError("Generic scFv VH/VL domains are implausibly short")
    return list(generic[:start]), list(generic[start:end]), list(generic[end:])


def _global_alignment_pairs(left: str, right: str) -> list[tuple[int, int]]:
    """Return deterministic Needleman-Wunsch raw-index correspondences."""
    rows, columns = len(left) + 1, len(right) + 1
    scores = [[0] * columns for _ in range(rows)]
    trace = [[""] * columns for _ in range(rows)]
    for i in range(1, rows):
        scores[i][0], trace[i][0] = -2 * i, "U"
    for j in range(1, columns):
        scores[0][j], trace[0][j] = -2 * j, "L"
    for i in range(1, rows):
        for j in range(1, columns):
            options = (
                (scores[i - 1][j - 1] + (2 if left[i - 1] == right[j - 1] else -1), "D"),
                (scores[i - 1][j] - 2, "U"),
                (scores[i][j - 1] - 2, "L"),
            )
            scores[i][j], trace[i][j] = max(options, key=lambda item: (item[0], "DLU".index(item[1]) * -1))
    pairs: list[tuple[int, int]] = []
    i, j = len(left), len(right)
    while i or j:
        direction = trace[i][j]
        if direction == "D":
            pairs.append((i, j))
            i, j = i - 1, j - 1
        elif direction == "U":
            i -= 1
        elif direction == "L":
            j -= 1
        else:
            raise AntibodyTemplatePreparationError("Internal sequence-alignment failure")
    return list(reversed(pairs))


def _build_scfv(
    paired: Mapping[str, Sequence[Residue]], generic: Sequence[Residue], contract: Mapping[str, Any]
) -> tuple[bytes, dict[str, Any]]:
    vh, vl = list(paired["H"]), list(paired["L"])
    generic_vh, generic_linker, generic_vl = _domain_boundaries(generic, len(vh), len(vl))
    transforms: dict[str, tuple[Any, Any, float]] = {}
    anchored_framework_rmsd: dict[str, float] = {}
    junction_anchor_translation: dict[str, float] = {}
    alignment_counts: dict[str, int] = {}
    for chain_name, source, reference, prefix in (
        ("VH", vh, generic_vh, "H"), ("VL", vl, generic_vl, "L")
    ):
        framework_indices = set(
            _framework_indices(len(source), contract["cdr_raw_ranges"], prefix)
        )
        pairs = [
            pair
            for pair in _global_alignment_pairs(
                "".join(r.amino_acid for r in source),
                "".join(r.amino_acid for r in reference),
            )
            if pair[0] in framework_indices
        ]
        if len(pairs) < 40:
            raise AntibodyTemplatePreparationError(
                f"Only {len(pairs)} framework CA correspondences for {chain_name}"
            )
        mobile = [source[i - 1].by_name["CA"].xyz for i, _ in pairs]
        target = [reference[j - 1].by_name["CA"].xyz for _, j in pairs]
        rotation, translation, best_fit_rmsd = _kabsch(mobile, target)
        anchor_source = source[-1].by_name["C"].xyz if chain_name == "VH" else source[0].by_name["N"].xyz
        anchor_target = reference[-1].by_name["C"].xyz if chain_name == "VH" else reference[0].by_name["N"].xyz
        transformed_anchor = _apply_transform(anchor_source, rotation, translation)
        correction = tuple(wanted - observed for wanted, observed in zip(anchor_target, transformed_anchor, strict=True))
        import numpy as np
        anchored_translation = translation + np.asarray(correction, float)
        fitted = np.asarray(mobile, float) @ rotation + anchored_translation
        ref = np.asarray(target, float)
        anchored_framework_rmsd[chain_name] = float(
            np.sqrt(np.mean(np.sum((fitted - ref) ** 2, axis=1)))
        )
        junction_anchor_translation[chain_name] = math.sqrt(sum(value * value for value in correction))
        transforms[chain_name] = (rotation, anchored_translation, best_fit_rmsd)
        alignment_counts[chain_name] = len(pairs)
    lines: list[str] = []
    serial, residue_number = 1, 1
    output_sequence: list[str] = []
    output_residues: list[Residue] = []
    for segment_name, residues in (("VH", vh), ("linker", generic_linker), ("VL", vl)):
        for residue in residues:
            output_sequence.append(residue.amino_acid)
            transformed_atoms: list[Atom] = []
            for atom in residue.atoms:
                xyz = atom.xyz
                if segment_name in transforms:
                    rotation, translation, _ = transforms[segment_name]
                    xyz = _apply_transform(xyz, rotation, translation)
                lines.append(_format_atom(atom, serial=serial, chain="A", number=residue_number,
                                          icode="", xyz=xyz))
                transformed_atoms.append(
                    replace(
                        atom, serial=serial, chain="A", number=residue_number,
                        icode="", x=xyz[0], y=xyz[1], z=xyz[2],
                    )
                )
                serial += 1
            output_residues.append(
                Residue("A", residue_number, "", residue.amino_acid, tuple(transformed_atoms))
            )
            residue_number += 1
    lines.extend(("TER", "END"))
    expected = contract["neutral_vh"] + LINKER + contract["neutral_vl"]
    actual = "".join(output_sequence)
    if actual != expected:
        raise AntibodyTemplatePreparationError("Constructed scFv sequence is not neutral VH-linker-VL")
    # The generic linker is preserved byte-for-coordinate (only chain/residue numbering changes).
    n_interface = _distance(vh[-1].by_name["C"].xyz, generic_linker[0].by_name["N"].xyz)
    c_interface = _distance(generic_linker[-1].by_name["C"].xyz, vl[0].by_name["N"].xyz)
    # Interface distances above must use transformed domains.
    n_c = _apply_transform(vh[-1].by_name["C"].xyz, *transforms["VH"][:2])
    c_n = _apply_transform(vl[0].by_name["N"].xyz, *transforms["VL"][:2])
    n_interface = _distance(n_c, generic_linker[0].by_name["N"].xyz)
    c_interface = _distance(generic_linker[-1].by_name["C"].xyz, c_n)
    if not 1.15 <= n_interface <= 1.55 or not 1.15 <= c_interface <= 1.55:
        raise AntibodyTemplatePreparationError(
            f"scFv peptide interfaces fail: VH-linker={n_interface:.3f}, linker-VL={c_interface:.3f} A"
        )
    output_peptide_geometry = _peptide_geometry(output_residues)
    output_clash_geometry = _assert_no_severe_clashes({"A": output_residues})
    return ("\n".join(lines) + "\n").encode("ascii"), {
        "sequence": actual, "linker_sequence": LINKER,
        "generic_linker_observed": GENERIC_LINKER_OBSERVED,
        "linker_is_standard_GGGGS3": LINKER == STANDARD_GGGGS3,
        "generic_linker_coordinate_or_identity_changes": [],
        "generic_linker_geometry_preserved_verbatim": True,
        "vh_framework_ca_alignment_pairs": alignment_counts["VH"],
        "vl_framework_ca_alignment_pairs": alignment_counts["VL"],
        "vh_framework_ca_best_fit_rmsd_to_generic_angstrom": round(transforms["VH"][2], 6),
        "vl_framework_ca_best_fit_rmsd_to_generic_angstrom": round(transforms["VL"][2], 6),
        "vh_framework_ca_junction_anchored_rmsd_angstrom": round(anchored_framework_rmsd["VH"], 6),
        "vl_framework_ca_junction_anchored_rmsd_angstrom": round(anchored_framework_rmsd["VL"], 6),
        "vh_junction_anchor_translation_angstrom": round(junction_anchor_translation["VH"], 6),
        "vl_junction_anchor_translation_angstrom": round(junction_anchor_translation["VL"], 6),
        "vh_linker_peptide_cn_angstrom": round(n_interface, 6),
        "linker_vl_peptide_cn_angstrom": round(c_interface, 6),
        "full_chain_peptide_geometry": output_peptide_geometry,
        "full_chain_clash_geometry": output_clash_geometry,
        "generic_linker_geometry_preserved": True,
    }


def _write_atomic(path: Path, payload: bytes) -> None:
    with tempfile.NamedTemporaryFile(mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        handle.write(payload)
        temporary = Path(handle.name)
    temporary.replace(path)


def _count_source_remarks(path: Path) -> int:
    return sum(line.startswith("REMARK PDBinfo-LABEL:") for line in path.read_text(errors="replace").splitlines())


def _remark_counts(path: Path) -> dict[str, int]:
    pattern = re.compile(r"^REMARK\s+PDBinfo-LABEL:\s+\d+\s+(H1|H2|H3|L1|L2|L3)\s*$")
    counts = {name: 0 for name in CDR_ORDER}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("REMARK PDBinfo-LABEL:"):
            continue
        match = pattern.match(line)
        if match is None:
            raise AntibodyTemplatePreparationError(f"Malformed HLT REMARK in {path}: {line}")
        counts[match.group(1)] += 1
    return counts


def _validate_coordinate_identity(
    reference: Mapping[str, Sequence[Residue]], candidate: Mapping[str, Sequence[Residue]],
    *, label: str, tolerance: float = 0.002,
) -> dict[str, Any]:
    maximum = 0.0
    atom_count = 0
    for chain in ("H", "L"):
        if len(reference[chain]) != len(candidate[chain]):
            raise AntibodyTemplatePreparationError(f"{label} changes residue count on chain {chain}")
        for left, right in zip(reference[chain], candidate[chain], strict=True):
            if left.amino_acid != right.amino_acid:
                raise AntibodyTemplatePreparationError(f"{label} changes sequence at chain {chain}")
            left_atoms, right_atoms = left.by_name, right.by_name
            if set(left_atoms) != set(right_atoms):
                raise AntibodyTemplatePreparationError(f"{label} changes atom set at {chain}{left.label}")
            for name in left_atoms:
                delta = _distance(left_atoms[name].xyz, right_atoms[name].xyz)
                maximum = max(maximum, delta)
                atom_count += 1
                if delta > tolerance:
                    raise AntibodyTemplatePreparationError(
                        f"{label} coordinate drift {delta:.6f} A exceeds {tolerance:.3f} A"
                    )
    return {
        "compared_atom_count": atom_count,
        "maximum_coordinate_delta_angstrom": round(maximum, 6),
        "tolerance_angstrom": tolerance,
        "comparison_key": "chain_and_sequence_order_and_atom_name",
        "numbering_may_be_changed_to_absolute_pose_indices": True,
    }


def prepare_antibody_template_inputs(
    *, paired_fv_pdbs: Mapping[str, str | Path], generic_scfv_pdb: str | Path,
    official_broad_hlt_pdbs: Mapping[str, str | Path] | None = None,
    prepared_scfv_pdbs: Mapping[str, str | Path] | None = None,
    output_dir: str | Path = REPO_ROOT / "input" / "template_structures",
    fasta_path: str | Path = REPO_ROOT / "validation" / "experimentally_validated_antibodies.fasta",
    numbering_evidence_path: str | Path = REPO_ROOT / "input" / "antibody_templates" / "chothia_numbering_evidence.json",
    heavy_anarci_csv: str | Path = REPO_ROOT / "input" / "antibody_templates" / "nfl_H.csv",
    light_anarci_csv: str | Path = REPO_ROOT / "input" / "antibody_templates" / "nfl_KL.csv",
    neutral_heavy_anarci_csv: str | Path | None = None,
    neutral_light_anarci_csv: str | Path | None = None,
    upstream_provenance_path: str | Path | None = None,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Prepare reviewed-input candidates without executing external software."""
    supplied = set(paired_fv_pdbs)
    if supplied != set(TEMPLATE_IDS):
        raise AntibodyTemplatePreparationError(
            f"paired_fv_pdbs must define exactly {TEMPLATE_IDS}; got {sorted(supplied)}"
        )
    paths = {key: Path(value).expanduser().resolve() for key, value in paired_fv_pdbs.items()}
    broad_paths: dict[str, Path] = {}
    if official_broad_hlt_pdbs is not None:
        if set(official_broad_hlt_pdbs) != set(TEMPLATE_IDS):
            raise AntibodyTemplatePreparationError(
                f"official_broad_hlt_pdbs must define exactly {TEMPLATE_IDS} when supplied"
            )
        broad_paths = {
            key: Path(value).expanduser().resolve()
            for key, value in official_broad_hlt_pdbs.items()
        }
    prepared_scfv_paths: dict[str, Path] = {}
    if prepared_scfv_pdbs is not None:
        if set(prepared_scfv_pdbs) != set(TEMPLATE_IDS):
            raise AntibodyTemplatePreparationError(
                f"prepared_scfv_pdbs must define exactly {TEMPLATE_IDS}; partial input is forbidden"
            )
        prepared_scfv_paths = {
            key: Path(value).expanduser().resolve()
            for key, value in prepared_scfv_pdbs.items()
        }
    generic_path = Path(generic_scfv_pdb).expanduser().resolve()
    fasta_path, numbering_evidence_path = Path(fasta_path), Path(numbering_evidence_path)
    heavy_anarci_csv, light_anarci_csv = Path(heavy_anarci_csv), Path(light_anarci_csv)
    neutral_h_path = Path(neutral_heavy_anarci_csv) if neutral_heavy_anarci_csv else None
    neutral_l_path = Path(neutral_light_anarci_csv) if neutral_light_anarci_csv else None
    if (neutral_h_path is None) != (neutral_l_path is None):
        raise AntibodyTemplatePreparationError(
            "neutral heavy/light ANARCI CSVs must be supplied together"
        )
    optional_paths = tuple(path for path in (neutral_h_path, neutral_l_path) if path is not None)
    upstream_path = Path(upstream_provenance_path).expanduser().resolve() if upstream_provenance_path else None
    for path in (*paths.values(), *broad_paths.values(), *prepared_scfv_paths.values(),
                 *optional_paths, *((upstream_path,) if upstream_path else ()),
                 generic_path, fasta_path, numbering_evidence_path,
                 heavy_anarci_csv, light_anarci_csv):
        if not path.is_file():
            raise AntibodyTemplatePreparationError(f"Required input is not a file: {path}")
    evidence_source = _read_numbering_evidence(numbering_evidence_path)
    declared = {x["logical_name"]: x["sha256"] for x in evidence_source["raw_annotation_artifacts"]}
    for path in (heavy_anarci_csv, light_anarci_csv):
        if declared.get(path.name) != _sha256(path):
            raise AntibodyTemplatePreparationError(f"Pinned ANARCI table hash mismatch: {path}")
    fasta = _read_fasta(fasta_path)
    contracts = _validate_contracts(
        fasta, _read_anarci_tables(heavy_anarci_csv, light_anarci_csv), evidence_source
    )
    neutral_tables = (
        _read_neutral_anarci_tables(neutral_h_path, neutral_l_path)
        if neutral_h_path is not None and neutral_l_path is not None else None
    )
    upstream_provenance = (
        _read_and_bind_upstream_provenance(
            upstream_path, paired_paths=paths, prepared_scfv_paths=prepared_scfv_paths,
            generic_path=generic_path, broad_paths=broad_paths,
        ) if upstream_path is not None else None
    )
    paired_structures: dict[str, dict[str, list[Residue]]] = {}
    paired_geometry: dict[str, dict[str, Any]] = {}
    broad_evidence: dict[str, dict[str, Any]] = {}
    for template_id, path in paths.items():
        chains = _parse_pdb(path, allowed_chains={"H", "L"})
        expected = {"H": contracts[template_id]["neutral_vh"], "L": contracts[template_id]["neutral_vl"]}
        for chain in ("H", "L"):
            sequence = "".join(r.amino_acid for r in chains[chain])
            if sequence != expected[chain]:
                raise AntibodyTemplatePreparationError(
                    f"{template_id} chain {chain} sequence is not the locked neutral seed"
                )
            expected_labels = [x[0] for x in contracts[template_id]["anarci_by_chain"]["VH" if chain == "H" else "VL"]]
            if neutral_tables is not None:
                neutral_key = (template_id, "VH" if chain == "H" else "VL")
                neutral_record = neutral_tables.get(neutral_key)
                if neutral_record is None:
                    raise AntibodyTemplatePreparationError(
                        f"Neutral ANARCI tables lack {neutral_key}"
                    )
                neutral_sequence = "".join(amino_acid for _, amino_acid in neutral_record)
                if neutral_sequence != expected[chain]:
                    raise AntibodyTemplatePreparationError(
                        f"Neutral ANARCI sequence mismatch for {template_id} chain {chain}"
                    )
                neutral_labels = [label for label, _ in neutral_record]
                if neutral_labels != expected_labels:
                    raise AntibodyTemplatePreparationError(
                        f"Neutral mutation changed Chothia layout for {template_id} chain {chain}"
                    )
            actual_labels = [r.label for r in chains[chain]]
            if actual_labels != expected_labels:
                raise AntibodyTemplatePreparationError(
                    f"{template_id} chain {chain} PDB numbering/insertions do not match pinned ANARCI"
                )
        paired_geometry[template_id] = {
            "peptide_geometry": {
                chain: _peptide_geometry(chains[chain]) for chain in ("H", "L")
            },
            "clash_geometry": _assert_no_severe_clashes(chains),
        }
        paired_structures[template_id] = chains
        if template_id in broad_paths:
            broad_chains = _parse_pdb(broad_paths[template_id], allowed_chains={"H", "L"})
            counts = _remark_counts(broad_paths[template_id])
            if any(count == 0 for count in counts.values()):
                raise AntibodyTemplatePreparationError(
                    f"Official broad HLT lacks one or more CDR REMARK classes: {counts}"
                )
            broad_evidence[template_id] = {
                "path": _portable(broad_paths[template_id]),
                "sha256": _sha256(broad_paths[template_id]),
                "remark_counts": counts,
                "coordinate_identity_to_paired_fv": _validate_coordinate_identity(
                    chains, broad_chains, label=f"{template_id} official broad HLT"
                ),
            }
    generic_chains = _parse_pdb(generic_path)
    if len(generic_chains) != 1:
        raise AntibodyTemplatePreparationError("Generic scFv PDB must contain exactly one chain")
    generic = next(iter(generic_chains.values()))
    generic_geometry = {
        "peptide_geometry": _peptide_geometry(generic),
        "clash_geometry": _assert_no_severe_clashes(generic_chains),
    }

    output_root = Path(output_dir).expanduser().resolve()
    output_paths = {
        template_id: {
            "hlt": output_root / f"{template_id}.exact_chothia.hlt.pdb",
            "scfv": output_root / f"{template_id}.neutral_seed.scfv.pdb",
        } for template_id in TEMPLATE_IDS
    }
    evidence_path = output_root / "antibody_template_evidence.json"
    fragment_path = output_root / "antibody_template_manifest_fragment.blocked.json"
    all_outputs = [evidence_path, fragment_path] + [p for group in output_paths.values() for p in group.values()]
    conflicts = [p for p in all_outputs if p.exists()]
    if conflicts and not overwrite:
        raise AntibodyTemplatePreparationError(
            "Refusing to overwrite existing template artifacts: " + ", ".join(map(str, conflicts))
        )
    output_root.mkdir(parents=True, exist_ok=True)
    records: dict[str, Any] = {}
    payloads: dict[Path, bytes] = {}
    for template_id in TEMPLATE_IDS:
        hlt_bytes, absolute = _write_hlt(
            output_paths[template_id]["hlt"], paired_structures[template_id],
            contracts[template_id]["cdr_chothia_labels"],
        )
        scfv_coordinate_source = "local_framework_CA_domain_graft"
        prepared_scfv_input: dict[str, Any] | None = None
        if template_id in prepared_scfv_paths:
            prepared_path = prepared_scfv_paths[template_id]
            prepared_chains = _parse_pdb(prepared_path, allowed_chains={"A"})
            prepared_residues = prepared_chains["A"]
            if [residue.number for residue in prepared_residues] != list(
                range(1, len(prepared_residues) + 1)
            ) or any(residue.icode for residue in prepared_residues):
                raise AntibodyTemplatePreparationError(
                    f"Prepared scFv {template_id} must use continuous chain-A numbering from 1"
                )
            expected_scfv_sequence = (
                contracts[template_id]["neutral_vh"] + LINKER + contracts[template_id]["neutral_vl"]
            )
            if "".join(residue.amino_acid for residue in prepared_residues) != expected_scfv_sequence:
                raise AntibodyTemplatePreparationError(
                    f"Prepared scFv {template_id} sequence is not neutral VH+observed linker+neutral VL"
                )
            prepared_peptide = _peptide_geometry(prepared_residues)
            prepared_clash = _assert_no_severe_clashes(prepared_chains)
            scfv_bytes = prepared_path.read_bytes()
            scfv_metrics = {
                "sequence": expected_scfv_sequence,
                "linker_sequence": LINKER,
                "generic_linker_observed": GENERIC_LINKER_OBSERVED,
                "linker_is_standard_GGGGS3": False,
                "generic_linker_coordinate_or_identity_changes": [],
                "source_prepared_scfv_peptide_geometry": prepared_peptide,
                "source_prepared_scfv_clash_geometry": prepared_clash,
                "external_refinement_method_claim_trusted": False,
                "external_refinement_input_hash_bound_by_provenance": (
                    upstream_provenance is not None
                ),
                "geometry_independently_revalidated_by_preparation_cli": True,
            }
            scfv_coordinate_source = "externally_refined_domain_graft"
            prepared_scfv_input = {
                "path": _portable(prepared_path), "sha256": _sha256(prepared_path)
            }
        else:
            scfv_bytes, scfv_metrics = _build_scfv(
                paired_structures[template_id], generic, contracts[template_id]
            )
        payloads[output_paths[template_id]["hlt"]] = hlt_bytes
        payloads[output_paths[template_id]["scfv"]] = scfv_bytes
        records[template_id] = {
            "source_paired_fv_pdb": _portable(paths[template_id]),
            "source_paired_fv_sha256": _sha256(paths[template_id]),
            "source_paired_fv_hlt_remark_count": _count_source_remarks(paths[template_id]),
            "official_converter_broad_hlt": broad_evidence.get(template_id),
            "official_converter_broad_labels_discarded_and_normalized": template_id in broad_evidence,
            "output_hlt_remark_counts": {
                name: len(contracts[template_id]["cdr_chothia_labels"][name])
                for name in CDR_ORDER
            },
            "neutral_cdr_seeds": NEUTRAL_CDR_SEEDS[template_id],
            "neutral_vh": contracts[template_id]["neutral_vh"],
            "neutral_vl": contracts[template_id]["neutral_vl"],
            "source_paired_fv_geometry": paired_geometry[template_id],
            "cdr_chothia_labels": contracts[template_id]["cdr_chothia_labels"],
            "hlt_absolute_pose_indices": absolute,
            "hlt_path": _portable(output_paths[template_id]["hlt"]),
            "hlt_sha256": sha256(hlt_bytes).hexdigest(),
            "scfv_path": _portable(output_paths[template_id]["scfv"]),
            "scfv_sha256": sha256(scfv_bytes).hexdigest(),
            "scfv_coordinate_source": scfv_coordinate_source,
            "prepared_scfv_input": prepared_scfv_input,
            "scfv_metrics": scfv_metrics,
        }
    for path, payload in payloads.items():
        _write_atomic(path, payload)
    evidence = {
        "schema": OUTPUT_SCHEMA,
        "execution_state": "blocked_pending_human_review",
        "preparation_cli_model_or_network_execution_performed": False,
        "annotation_contract": "ANARCI 2020.04.23 Chothia exact campaign CDRs",
        "known_positive_cdr_sequences_used": False,
        "fasta": {"path": _portable(fasta_path), "sha256": _sha256(fasta_path)},
        "numbering_evidence": {"path": _portable(numbering_evidence_path), "sha256": _sha256(numbering_evidence_path)},
        "upstream_coordinate_generation_provenance": upstream_provenance,
        "anarci_tables": [
            {"path": _portable(heavy_anarci_csv), "sha256": _sha256(heavy_anarci_csv)},
            {"path": _portable(light_anarci_csv), "sha256": _sha256(light_anarci_csv)},
        ],
        "neutral_anarci_tables": (
            [
                {"path": _portable(neutral_h_path), "sha256": _sha256(neutral_h_path)},
                {"path": _portable(neutral_l_path), "sha256": _sha256(neutral_l_path)},
            ] if neutral_h_path is not None and neutral_l_path is not None else None
        ),
        "generic_scfv": {
            "path": _portable(generic_path), "sha256": _sha256(generic_path),
            "observed_linker": LINKER, "is_standard_GGGGS3": False,
            "geometry": generic_geometry,
        },
        "templates": records,
        "review": {"status": "blocked_pending_human_review", "reviewed_by": "", "reviewed_at": "", "contracts_acknowledged": False},
        "real_model_handoff_authorized": False,
    }
    evidence_bytes = _json_bytes(evidence)
    _write_atomic(evidence_path, evidence_bytes)
    fragment = {
        "schema": MANIFEST_FRAGMENT_SCHEMA,
        "execution_state": "blocked_pending_human_review",
        "evidence_path": _portable(evidence_path),
        "evidence_sha256": sha256(evidence_bytes).hexdigest(),
        "framework_coordinate_inputs": {
            "rfantibody_hlt_pdbs": {t: _portable(output_paths[t]["hlt"]) for t in TEMPLATE_IDS},
            "germinal_scfv_pdbs": {t: _portable(output_paths[t]["scfv"]) for t in TEMPLATE_IDS},
        },
        "review": {"reviewed_by": "", "reviewed_at": "", "contracts_acknowledged": False},
        "real_model_handoff_authorized": False,
        "note": "Fragment only; do not treat as a formal target manifest.",
    }
    fragment_bytes = _json_bytes(fragment)
    _write_atomic(fragment_path, fragment_bytes)
    return {
        "execution_state": "blocked_pending_human_review",
        "output_dir": _portable(output_root),
        "evidence": _portable(evidence_path),
        "manifest_fragment": _portable(fragment_path),
        "templates": records,
        "real_model_handoff_authorized": False,
    }


def _parse_template_path(values: Sequence[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise AntibodyTemplatePreparationError("--paired-fv must be TEMPLATE_ID=PATH")
        key, path = value.split("=", 1)
        if key in result or not path:
            raise AntibodyTemplatePreparationError(f"Invalid/duplicate --paired-fv: {value}")
        result[key] = path
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare audited antibody coordinate templates")
    parser.add_argument("--paired-fv", action="append", required=True, metavar="TEMPLATE_ID=PATH")
    parser.add_argument(
        "--official-broad-hlt", action="append", default=[], metavar="TEMPLATE_ID=PATH",
        help="Optional official chothia2HLT.py output retained only as normalization evidence",
    )
    parser.add_argument(
        "--prepared-scfv", action="append", default=[], metavar="TEMPLATE_ID=PATH",
        help="Two externally refined neutral scFv PDBs; partial input is forbidden",
    )
    parser.add_argument("--generic-scfv", required=True)
    parser.add_argument("--output-dir", default=str(REPO_ROOT / "input" / "template_structures"))
    parser.add_argument("--fasta", default=str(REPO_ROOT / "validation" / "experimentally_validated_antibodies.fasta"))
    parser.add_argument("--numbering-evidence", default=str(REPO_ROOT / "input" / "antibody_templates" / "chothia_numbering_evidence.json"))
    parser.add_argument("--heavy-anarci-csv", default=str(REPO_ROOT / "input" / "antibody_templates" / "nfl_H.csv"))
    parser.add_argument("--light-anarci-csv", default=str(REPO_ROOT / "input" / "antibody_templates" / "nfl_KL.csv"))
    parser.add_argument("--neutral-heavy-anarci-csv")
    parser.add_argument("--neutral-light-anarci-csv")
    parser.add_argument("--upstream-provenance")
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = prepare_antibody_template_inputs(
            paired_fv_pdbs=_parse_template_path(args.paired_fv), generic_scfv_pdb=args.generic_scfv,
            official_broad_hlt_pdbs=(
                _parse_template_path(args.official_broad_hlt)
                if args.official_broad_hlt else None
            ),
            prepared_scfv_pdbs=(
                _parse_template_path(args.prepared_scfv)
                if args.prepared_scfv else None
            ),
            output_dir=args.output_dir, fasta_path=args.fasta,
            numbering_evidence_path=args.numbering_evidence,
            heavy_anarci_csv=args.heavy_anarci_csv, light_anarci_csv=args.light_anarci_csv,
            neutral_heavy_anarci_csv=args.neutral_heavy_anarci_csv,
            neutral_light_anarci_csv=args.neutral_light_anarci_csv,
            upstream_provenance_path=args.upstream_provenance,
            overwrite=args.overwrite,
        )
    except AntibodyTemplatePreparationError as exc:
        print(f"ERROR: {exc}", file=__import__("sys").stderr)
        return 2
    print(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
