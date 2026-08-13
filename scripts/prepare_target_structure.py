#!/usr/bin/env python3
"""Prepare a review-blocked NEFL 280--377 AlphaFold DB structure bundle.

This script validates the pinned AlphaFold DB P07196 v6 PDB before extracting
chain A residues 280--377.  It writes a cropped PDB, FASTA, machine-readable
evidence, and a *candidate* target manifest.  The candidate manifest is
deliberately blocked and cannot be consumed by ``prepare_real_model_jobs.py``
until a human reviewer creates the formal reviewed manifest.

Only the Python standard library is used.  No network request or model is run.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import statistics
import sys
import tempfile
from typing import Any, Iterable, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILENAME = "AF-P07196-F1-model_v6.pdb"
SOURCE_URL = "https://alphafold.ebi.ac.uk/files/AF-P07196-F1-model_v6.pdb"
SOURCE_SHA256 = "37912aac5cefd85b177e754a7c55c10c0f50166baf7f30b012151492eae300b1"
SOURCE_ACCESSION = "P07196"
SOURCE_CHAIN = "A"
SOURCE_LENGTH = 543
TARGET_START = 280
TARGET_END = 377
TARGET_SEQUENCE = (
    "FKSRFTVLTESAAKNTDAVRAAKDEVSESRRLLKAKTLEIEACRGMNEALEK"
    "QLQELEDKQNADISAMQDTINKLENELRTTKSEMARYLKEYQDLLN"
)
TARGET_PDB_FILENAME = "NEFL_P07196_AFDB_v6_280-377_chainA.pdb"
TARGET_FASTA_FILENAME = "NEFL_P07196_AFDB_v6_280-377.fasta"
EVIDENCE_FILENAME = "NEFL_P07196_AFDB_v6_280-377_evidence.json"
CANDIDATE_MANIFEST_FILENAME = "target_structure_manifest.candidate.blocked.json"
EVIDENCE_SCHEMA = "nfl_ab_design.target_structure_evidence.v1"
TARGET_MANIFEST_SCHEMA = "nfl_ab_design.target_structure_manifest.v1"
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})
BACKBONE_AND_CB = frozenset({"N", "CA", "C", "O", "CB"})

# Canonical heavy-atom topology used to reject incomplete AlphaFold coordinates.
# Extra terminal atoms such as OXT are permitted, but every atom below is required.
EXPECTED_HEAVY_ATOMS_BY_AA = {
    "A": BACKBONE_ATOMS | {"CB"},
    "C": BACKBONE_ATOMS | {"CB", "SG"},
    "D": BACKBONE_ATOMS | {"CB", "CG", "OD1", "OD2"},
    "E": BACKBONE_ATOMS | {"CB", "CG", "CD", "OE1", "OE2"},
    "F": BACKBONE_ATOMS | {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "G": BACKBONE_ATOMS,
    "H": BACKBONE_ATOMS | {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "I": BACKBONE_ATOMS | {"CB", "CG1", "CG2", "CD1"},
    "K": BACKBONE_ATOMS | {"CB", "CG", "CD", "CE", "NZ"},
    "L": BACKBONE_ATOMS | {"CB", "CG", "CD1", "CD2"},
    "M": BACKBONE_ATOMS | {"CB", "CG", "SD", "CE"},
    "N": BACKBONE_ATOMS | {"CB", "CG", "OD1", "ND2"},
    "P": BACKBONE_ATOMS | {"CB", "CG", "CD"},
    "Q": BACKBONE_ATOMS | {"CB", "CG", "CD", "OE1", "NE2"},
    "R": BACKBONE_ATOMS | {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "S": BACKBONE_ATOMS | {"CB", "OG"},
    "T": BACKBONE_ATOMS | {"CB", "OG1", "CG2"},
    "V": BACKBONE_ATOMS | {"CB", "CG1", "CG2"},
    "W": BACKBONE_ATOMS
    | {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "Y": BACKBONE_ATOMS
    | {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
}

AA3_TO_AA1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

EPITOPES: tuple[dict[str, Any], ...] = (
    {
        "epitope_id": "helix_surface_323_331",
        "start_1_based": 323,
        "end_1_based_inclusive": 331,
        "sequence": "RGMNEALEK",
        "proposed_hotspots": (325, 329),
        "proposal_basis": (
            "Met325 and Leu329 form an i/i+4 pair on the same predicted "
            "solvent-accessible alpha-helical face; no cysteine contact is required."
        ),
    },
    {
        "epitope_id": "C_boundary_368_377",
        "start_1_based": 368,
        "end_1_based_inclusive": 377,
        "sequence": "YLKEYQDLLN",
        "proposed_hotspots": (368, 372, 375),
        "proposal_basis": (
            "Tyr368, Tyr372, and Leu375 provide separated side-chain anchors "
            "across the configured C-terminal boundary epitope."
        ),
    },
)


class TargetStructurePreparationError(ValueError):
    """Raised when the pinned source or extracted structure violates a contract."""


@dataclass(frozen=True)
class Atom:
    """One parsed PDB ATOM record."""

    serial: int
    name: str
    altloc: str
    residue_name: str
    chain_id: str
    residue_number: int
    insertion_code: str
    x: float
    y: float
    z: float
    occupancy: float
    b_factor: float
    element: str
    source_line: str

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)


@dataclass(frozen=True)
class Residue:
    """One canonical residue and its coordinate records."""

    number: int
    residue_name: str
    amino_acid: str
    atoms: tuple[Atom, ...]

    @property
    def by_name(self) -> dict[str, Atom]:
        return {atom.name: atom for atom in self.atoms}

    @property
    def plddt(self) -> float:
        values = [atom.b_factor for atom in self.atoms]
        if max(values) - min(values) > 0.011:
            raise TargetStructurePreparationError(
                f"Residue A{self.number} has inconsistent per-atom pLDDT/B-factor values"
            )
        return values[0]


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _bytes_sha256(value: bytes) -> str:
    return sha256(value).hexdigest()


def _text_sha256(value: str) -> str:
    return sha256(value.encode("utf-8")).hexdigest()


def _portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _parse_atom_line(line: str, *, line_number: int) -> Atom:
    if len(line) < 66:
        raise TargetStructurePreparationError(
            f"Truncated ATOM record at source line {line_number}"
        )
    try:
        atom = Atom(
            serial=int(line[6:11]),
            name=line[12:16].strip(),
            altloc=line[16:17].strip(),
            residue_name=line[17:20].strip().upper(),
            chain_id=line[21:22],
            residue_number=int(line[22:26]),
            insertion_code=line[26:27].strip(),
            x=float(line[30:38]),
            y=float(line[38:46]),
            z=float(line[46:54]),
            occupancy=float(line[54:60]),
            b_factor=float(line[60:66]),
            element=(line[76:78].strip().upper() if len(line) >= 78 else ""),
            source_line=line.rstrip("\r\n"),
        )
    except ValueError as exc:
        raise TargetStructurePreparationError(
            f"Malformed ATOM record at source line {line_number}"
        ) from exc
    if atom.altloc:
        raise TargetStructurePreparationError(
            f"Alternate location {atom.altloc!r} at {atom.chain_id}{atom.residue_number} "
            "is unsupported by this unambiguous coordinate contract"
        )
    if atom.residue_name not in AA3_TO_AA1:
        raise TargetStructurePreparationError(
            f"Non-canonical residue {atom.residue_name!r} in ATOM records"
        )
    if not atom.name:
        raise TargetStructurePreparationError(
            f"Blank atom name at source line {line_number}"
        )
    if not all(math.isfinite(value) for value in (*atom.xyz, atom.occupancy, atom.b_factor)):
        raise TargetStructurePreparationError(
            f"Non-finite coordinate/metadata at source line {line_number}"
        )
    return atom


def _parse_seqres(lines: Sequence[str], *, chain_id: str) -> tuple[str, int]:
    residues: list[str] = []
    declared_lengths: set[int] = set()
    serials: list[int] = []
    for line_number, line in enumerate(lines, start=1):
        if not line.startswith("SEQRES") or line[11:12] != chain_id:
            continue
        try:
            serials.append(int(line[7:10]))
            declared_lengths.add(int(line[13:17]))
        except ValueError as exc:
            raise TargetStructurePreparationError(
                f"Malformed SEQRES header at source line {line_number}"
            ) from exc
        for residue_name in line[19:70].split():
            try:
                residues.append(AA3_TO_AA1[residue_name.upper()])
            except KeyError as exc:
                raise TargetStructurePreparationError(
                    f"Non-canonical SEQRES residue {residue_name!r}"
                ) from exc
    if not residues:
        raise TargetStructurePreparationError(
            f"No SEQRES records found for chain {chain_id!r}"
        )
    if len(declared_lengths) != 1:
        raise TargetStructurePreparationError(
            f"SEQRES chain {chain_id!r} has inconsistent declared lengths"
        )
    if serials != list(range(1, len(serials) + 1)):
        raise TargetStructurePreparationError(
            f"SEQRES chain {chain_id!r} record serials are not contiguous"
        )
    declared_length = next(iter(declared_lengths))
    if len(residues) != declared_length:
        raise TargetStructurePreparationError(
            f"SEQRES chain {chain_id!r} declares {declared_length} residues but "
            f"contains {len(residues)}"
        )
    return "".join(residues), declared_length


def _parse_and_validate_source(
    source_path: Path,
    *,
    expected_sha256: str,
) -> tuple[bytes, str, dict[int, Residue]]:
    if not source_path.is_file():
        raise TargetStructurePreparationError(
            f"Pinned AlphaFold DB source PDB is missing: {source_path}"
        )
    source_bytes = source_path.read_bytes()
    actual_sha256 = _bytes_sha256(source_bytes)
    if actual_sha256 != expected_sha256:
        raise TargetStructurePreparationError(
            "AlphaFold DB source SHA-256 mismatch: expected "
            f"{expected_sha256}, observed {actual_sha256}. Refusing unpinned coordinates."
        )
    try:
        source_text = source_bytes.decode("ascii")
    except UnicodeDecodeError as exc:
        raise TargetStructurePreparationError("Source PDB is not ASCII") from exc
    lines = source_text.splitlines()
    if not any(
        line.startswith("DBREF") and SOURCE_ACCESSION in line for line in lines
    ):
        raise TargetStructurePreparationError(
            f"Source PDB lacks a DBREF for {SOURCE_ACCESSION}"
        )
    model_records = [line for line in lines if line.startswith("MODEL")]
    end_model_records = [line for line in lines if line.startswith("ENDMDL")]
    if len(model_records) != 1 or len(end_model_records) != 1:
        raise TargetStructurePreparationError(
            "Source PDB must contain exactly one MODEL/ENDMDL coordinate model"
        )
    try:
        model_number = int(model_records[0][10:14])
    except ValueError as exc:
        raise TargetStructurePreparationError("Malformed MODEL record") from exc
    if model_number != 1:
        raise TargetStructurePreparationError(
            f"Expected MODEL 1, observed MODEL {model_number}"
        )

    seqres_sequence, declared_length = _parse_seqres(lines, chain_id=SOURCE_CHAIN)
    if declared_length != SOURCE_LENGTH:
        raise TargetStructurePreparationError(
            f"Expected {SOURCE_LENGTH} SEQRES residues for P07196, observed {declared_length}"
        )

    atom_groups: dict[tuple[str, int, str], list[Atom]] = {}
    seen_atom_keys: set[tuple[str, int, str, str]] = set()
    other_atom_chains: set[str] = set()
    for line_number, line in enumerate(lines, start=1):
        if line.startswith("HETATM"):
            chain = line[21:22] if len(line) > 21 else ""
            if chain == SOURCE_CHAIN:
                raise TargetStructurePreparationError(
                    f"Unexpected HETATM in source chain {SOURCE_CHAIN} at line {line_number}"
                )
            continue
        if not line.startswith("ATOM"):
            continue
        atom = _parse_atom_line(line, line_number=line_number)
        if atom.chain_id != SOURCE_CHAIN:
            other_atom_chains.add(atom.chain_id)
            continue
        key = (atom.chain_id, atom.residue_number, atom.insertion_code, atom.name)
        if key in seen_atom_keys:
            raise TargetStructurePreparationError(
                f"Duplicate ATOM name {atom.name!r} at A{atom.residue_number}"
            )
        seen_atom_keys.add(key)
        atom_groups.setdefault(
            (atom.chain_id, atom.residue_number, atom.insertion_code), []
        ).append(atom)
    if other_atom_chains:
        raise TargetStructurePreparationError(
            f"Pinned source unexpectedly contains ATOM chains besides A: {sorted(other_atom_chains)}"
        )

    observed_addresses = sorted((number, insertion) for _, number, insertion in atom_groups)
    expected_addresses = [(position, "") for position in range(1, SOURCE_LENGTH + 1)]
    if observed_addresses != expected_addresses:
        missing = [address for address in expected_addresses if address not in observed_addresses]
        extra = [address for address in observed_addresses if address not in expected_addresses]
        raise TargetStructurePreparationError(
            "Source coordinate residues are not the complete contiguous P07196 1--543 "
            f"chain A set; missing={missing[:10]}, extra={extra[:10]}"
        )

    residues: dict[int, Residue] = {}
    atom_sequence: list[str] = []
    for position in range(1, SOURCE_LENGTH + 1):
        atoms = tuple(atom_groups[(SOURCE_CHAIN, position, "")])
        residue_names = {atom.residue_name for atom in atoms}
        if len(residue_names) != 1:
            raise TargetStructurePreparationError(
                f"ATOM records at A{position} contain multiple residue identities"
            )
        residue_name = next(iter(residue_names))
        amino_acid = AA3_TO_AA1[residue_name]
        atom_names = {atom.name for atom in atoms}
        missing_backbone = sorted(BACKBONE_ATOMS - atom_names)
        if missing_backbone:
            raise TargetStructurePreparationError(
                f"Incomplete backbone at A{position}; missing {missing_backbone}"
            )
        missing_heavy_atoms = sorted(
            EXPECTED_HEAVY_ATOMS_BY_AA[amino_acid] - atom_names
        )
        if missing_heavy_atoms:
            raise TargetStructurePreparationError(
                f"Incomplete canonical heavy-atom coordinates at A{position}; "
                f"missing {missing_heavy_atoms}"
            )
        residue = Residue(position, residue_name, amino_acid, atoms)
        if not 0.0 <= residue.plddt <= 100.0:
            raise TargetStructurePreparationError(
                f"pLDDT/B-factor outside 0--100 at A{position}: {residue.plddt}"
            )
        residues[position] = residue
        atom_sequence.append(amino_acid)

    atom_sequence_text = "".join(atom_sequence)
    if atom_sequence_text != seqres_sequence:
        mismatch = next(
            (
                position
                for position, (left, right) in enumerate(
                    zip(atom_sequence_text, seqres_sequence, strict=True), start=1
                )
                if left != right
            ),
            None,
        )
        raise TargetStructurePreparationError(
            f"ATOM sequence does not equal SEQRES sequence; first mismatch at {mismatch}"
        )
    target_sequence = atom_sequence_text[TARGET_START - 1 : TARGET_END]
    if target_sequence != TARGET_SEQUENCE:
        raise TargetStructurePreparationError(
            "Extracted P07196 280--377 sequence does not match the pinned campaign "
            f"contract: expected {TARGET_SEQUENCE}, observed {target_sequence}"
        )
    for epitope in EPITOPES:
        start = int(epitope["start_1_based"])
        end = int(epitope["end_1_based_inclusive"])
        observed = atom_sequence_text[start - 1 : end]
        if observed != epitope["sequence"]:
            raise TargetStructurePreparationError(
                f"Epitope {epitope['epitope_id']} sequence mismatch: {observed}"
            )
        proposed = tuple(int(value) for value in epitope["proposed_hotspots"])
        if not proposed or len(set(proposed)) != len(proposed):
            raise TargetStructurePreparationError(
                f"Epitope {epitope['epitope_id']} has invalid hotspot proposals"
            )
        if any(position < start or position > end for position in proposed):
            raise TargetStructurePreparationError(
                f"Epitope {epitope['epitope_id']} hotspot proposal falls outside its window"
            )
    return source_bytes, atom_sequence_text, residues


def _distance(
    left: tuple[float, float, float], right: tuple[float, float, float]
) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right, strict=True)))


def _centroid(points: Iterable[tuple[float, float, float]]) -> tuple[float, float, float]:
    values = list(points)
    if not values:
        raise TargetStructurePreparationError("Cannot calculate a centroid of no points")
    return tuple(sum(point[axis] for point in values) / len(values) for axis in range(3))  # type: ignore[return-value]


def _sidechain_probe_atom(residue: Residue) -> Atom:
    by_name = residue.by_name
    if residue.amino_acid == "C" and "SG" in by_name:
        return by_name["SG"]
    ca = by_name["CA"]
    candidates = [
        atom
        for atom in residue.atoms
        if atom.name not in BACKBONE_AND_CB and not atom.element.startswith("H")
    ]
    if not candidates:
        candidates = [by_name.get("CB", ca)]
    return max(candidates, key=lambda atom: (_distance(atom.xyz, ca.xyz), atom.name))


def _hotspot_evidence(
    *,
    residues: Mapping[int, Residue],
    epitope: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    start = int(epitope["start_1_based"])
    end = int(epitope["end_1_based_inclusive"])
    proposed = [int(value) for value in epitope["proposed_hotspots"]]
    epitope_centroid = _centroid(
        residues[position].by_name["CA"].xyz for position in range(start, end + 1)
    )
    rows: list[dict[str, Any]] = []
    all_heavy_atoms = [
        (position, atom)
        for position, residue in residues.items()
        for atom in residue.atoms
        if not atom.element.startswith("H")
    ]
    for position in proposed:
        residue = residues[position]
        by_name = residue.by_name
        ca = by_name["CA"]
        probe = _sidechain_probe_atom(residue)
        nonlocal_ca_neighbors = sorted(
            other_position
            for other_position, other_residue in residues.items()
            if abs(other_position - position) > 2
            and _distance(ca.xyz, other_residue.by_name["CA"].xyz) <= 10.0
        )
        contacting_residues = sorted(
            {
                other_position
                for other_position, atom in all_heavy_atoms
                if abs(other_position - position) > 1
                and _distance(probe.xyz, atom.xyz) <= 6.0
            }
        )
        nonlocal_distances = [
            _distance(probe.xyz, atom.xyz)
            for other_position, atom in all_heavy_atoms
            if abs(other_position - position) > 1
        ]
        atom_names = sorted(by_name)
        rows.append(
            {
                "full_coordinate_1_based": position,
                "pdb_address": {
                    "chain_id": SOURCE_CHAIN,
                    "residue_number": position,
                    "insertion_code": "",
                },
                "local_coordinate_1_based": position - TARGET_START + 1,
                "residue_name_3_letter": residue.residue_name,
                "amino_acid_1_letter": residue.amino_acid,
                "atom_presence": {
                    "atom_names": atom_names,
                    "backbone_N_CA_C_O_complete": BACKBONE_ATOMS.issubset(by_name),
                    "sidechain_probe_atom": probe.name,
                    "sidechain_probe_present": probe.name in by_name,
                    "cys_anchor_SG_present": (
                        "SG" in by_name if residue.amino_acid == "C" else None
                    ),
                },
                "alphafold_plddt": round(residue.plddt, 2),
                "geometry": {
                    "ca_distance_to_epitope_ca_centroid_A": round(
                        _distance(ca.xyz, epitope_centroid), 3
                    ),
                    "sidechain_probe_nearest_nonlocal_heavy_atom_distance_A": round(
                        min(nonlocal_distances), 3
                    ),
                },
                "exposure_proxies_full_length_model_context": {
                    "nonlocal_ca_neighbor_count_within_10A": len(
                        nonlocal_ca_neighbors
                    ),
                    "nonlocal_ca_neighbor_residues_within_10A": nonlocal_ca_neighbors,
                    "sidechain_probe_contacting_residue_count_within_6A": len(
                        contacting_residues
                    ),
                    "sidechain_probe_contacting_residues_within_6A": contacting_residues,
                },
                "proposal_status": "proposed_pending_human_review",
            }
        )

    pairwise: list[dict[str, Any]] = []
    for left_index, left in enumerate(proposed):
        for right in proposed[left_index + 1 :]:
            pairwise.append(
                {
                    "residue_1": left,
                    "residue_2": right,
                    "ca_distance_A": round(
                        _distance(
                            residues[left].by_name["CA"].xyz,
                            residues[right].by_name["CA"].xyz,
                        ),
                        3,
                    ),
                    "sidechain_probe_distance_A": round(
                        _distance(
                            _sidechain_probe_atom(residues[left]).xyz,
                            _sidechain_probe_atom(residues[right]).xyz,
                        ),
                        3,
                    ),
                }
            )
    return rows, pairwise


def _seqres_lines(sequence: str) -> list[str]:
    aa1_to_aa3 = {value: key for key, value in AA3_TO_AA1.items()}
    residues = [aa1_to_aa3[amino_acid] for amino_acid in sequence]
    lines: list[str] = []
    for index in range(0, len(residues), 13):
        serial = index // 13 + 1
        chunk = " ".join(residues[index : index + 13])
        lines.append(f"SEQRES {serial:3d} A {len(residues):4d}  {chunk}")
    return lines


def _cropped_pdb_text(
    residues: Mapping[int, Residue], *, source_sha256: str
) -> str:
    target_residues = [residues[position] for position in range(TARGET_START, TARGET_END + 1)]
    atom_lines = [atom.source_line for residue in target_residues for atom in residue.atoms]
    last = target_residues[-1]
    ter_serial = max(atom.serial for residue in target_residues for atom in residue.atoms) + 1
    lines = [
        "HEADER    THEORETICAL MODEL / ALPHAFOLD DB",
        "TITLE     NEFL P07196 ALPHAFOLD DB V6, CHAIN A RESIDUES 280-377",
        "COMPND    MOL_ID: 1; MOLECULE: NEUROFILAMENT LIGHT POLYPEPTIDE; CHAIN: A",
        "SOURCE    MOL_ID: 1; ORGANISM_SCIENTIFIC: HOMO SAPIENS; ORGANISM_TAXID: 9606",
        "REMARK 900 SOURCE ALPHAFOLD DB AF-P07196-F1 MODEL V6",
        f"REMARK 900 SOURCE URL {SOURCE_URL}",
        f"REMARK 900 SOURCE SHA256 {source_sha256}",
        "REMARK 900 THEORETICAL PREDICTION; NOT AN EXPERIMENTALLY DETERMINED STRUCTURE",
        "REMARK 900 B-FACTOR FIELD STORES ALPHAFOLD pLDDT (0-100)",
        "REMARK 900 FULL-LENGTH UNIPROT RESIDUE NUMBERS 280-377 ARE RETAINED",
        "DBREF  XXXX A  280   377  UNP    P07196   NFL_HUMAN      280    377",
        *_seqres_lines(TARGET_SEQUENCE),
        "MODEL        1",
        *atom_lines,
        f"TER   {ter_serial:5d}      {last.residue_name:>3s} A{last.number:4d}",
        "ENDMDL",
        "END",
    ]
    return "\n".join(lines) + "\n"


def _round_stats(values: Sequence[float]) -> dict[str, Any]:
    return {
        "count": len(values),
        "mean": round(statistics.fmean(values), 2),
        "median": round(statistics.median(values), 2),
        "minimum": round(min(values), 2),
        "maximum": round(max(values), 2),
        "counts_by_confidence_band": {
            "very_high_ge_90": sum(value >= 90.0 for value in values),
            "confident_70_to_lt_90": sum(70.0 <= value < 90.0 for value in values),
            "low_50_to_lt_70": sum(50.0 <= value < 70.0 for value in values),
            "very_low_lt_50": sum(value < 50.0 for value in values),
        },
    }


def _build_evidence(
    *,
    source_copy_path: Path,
    target_pdb_path: Path,
    target_pdb_text: str,
    target_fasta_path: Path,
    target_fasta_text: str,
    full_sequence: str,
    residues: Mapping[int, Residue],
    expected_source_sha256: str,
) -> dict[str, Any]:
    target_positions = list(range(TARGET_START, TARGET_END + 1))
    plddt_values = [residues[position].plddt for position in target_positions]
    per_residue_plddt = {
        str(position): round(residues[position].plddt, 2)
        for position in target_positions
    }
    epitope_rows: list[dict[str, Any]] = []
    for epitope in EPITOPES:
        hotspot_rows, pairwise = _hotspot_evidence(
            residues=residues, epitope=epitope
        )
        start = int(epitope["start_1_based"])
        end = int(epitope["end_1_based_inclusive"])
        epitope_rows.append(
            {
                "epitope_id": epitope["epitope_id"],
                "full_coordinate_range_1_based_inclusive": [start, end],
                "sequence": epitope["sequence"],
                "structure_basis": "AlphaFold DB P07196 v6 theoretical monomer prediction",
                "proposal_basis": epitope["proposal_basis"],
                "proposed_hotspots_pending_human_review": list(
                    epitope["proposed_hotspots"]
                ),
                "alphafold_plddt_statistics": _round_stats(
                    [residues[position].plddt for position in range(start, end + 1)]
                ),
                "hotspot_evidence": hotspot_rows,
                "pairwise_proposed_hotspot_geometry": pairwise,
            }
        )

    full_to_pdb = {
        str(position): {
            "chain_id": SOURCE_CHAIN,
            "residue_number": position,
            "insertion_code": "",
        }
        for position in target_positions
    }
    full_to_local = {
        str(position): position - TARGET_START + 1 for position in target_positions
    }
    return {
        "schema": EVIDENCE_SCHEMA,
        "execution_state": "blocked_pending_human_review",
        "scientific_status": "theoretical_prediction_not_experimentally_determined",
        "source": {
            "database": "AlphaFold Protein Structure Database",
            "accession": SOURCE_ACCESSION,
            "model_id": "AF-P07196-F1-model_v6",
            "model_version": "v6",
            "download_url": SOURCE_URL,
            "source_pdb_path": _portable_path(source_copy_path),
            "sha256": expected_source_sha256,
            "license_noted_in_source_pdb": "CC-BY-4.0",
            "source_disclaimer": (
                "AlphaFold DB coordinates are theoretical predictions and are not "
                "a substitute for experimentally determined structure evidence."
            ),
        },
        "extraction": {
            "target_pdb_path": _portable_path(target_pdb_path),
            "target_pdb_sha256": _text_sha256(target_pdb_text),
            "target_fasta_path": _portable_path(target_fasta_path),
            "target_fasta_sha256": _text_sha256(target_fasta_text),
            "source_chain": SOURCE_CHAIN,
            "full_coordinate_range_1_based_inclusive": [TARGET_START, TARGET_END],
            "local_coordinate_range_1_based_inclusive": [1, len(target_positions)],
            "residue_count": len(target_positions),
            "full_residue_numbering_retained_in_pdb": True,
            "antigen_sequence_in_pdb_order": TARGET_SEQUENCE,
        },
        "validation": {
            "pinned_source_sha256_match": True,
            "source_seqres_length": len(full_sequence),
            "source_atom_residue_count": len(residues),
            "source_seqres_equals_atom_sequence": True,
            "source_coordinates_contiguous_1_through_543": True,
            "target_coordinates_contiguous_280_through_377": True,
            "target_sequence_equals_campaign_contract": True,
            "target_backbone_N_CA_C_O_complete_for_all_residues": True,
            "source_canonical_heavy_atom_topology_complete_for_all_residues": True,
            "target_canonical_heavy_atom_topology_complete_for_all_residues": True,
            "epitope_sequences_match_campaign_contract": True,
        },
        "coordinate_mappings": {
            "coordinate_system": "UniProt_P07196_1_based_inclusive",
            "full_coordinate_to_pdb": full_to_pdb,
            "full_coordinate_to_local_1_based": full_to_local,
        },
        "alphafold_plddt": {
            "field_source": "PDB B-factor column",
            "interpretation": "AlphaFold confidence, not experimental B-factor",
            "target_statistics": _round_stats(plddt_values),
            "per_full_coordinate": per_residue_plddt,
        },
        "exposure_proxy_definition": {
            "scope": "all 543 residues of the full-length AlphaFold DB source model",
            "ca_neighbor_proxy": (
                "Count residues with CA within 10 A after excluding sequence-near "
                "positions |delta|<=2; fewer neighbors can suggest less geometric "
                "occlusion but do not establish solvent exposure."
            ),
            "sidechain_contact_proxy": (
                "Use Cys SG or the farthest resolved side-chain heavy atom from CA; "
                "count nonlocal residues with any heavy atom within 6 A after excluding "
                "|delta|<=1. Fewer contacts can suggest exposure."
            ),
            "not_sasa": True,
            "limitations": (
                "These are deterministic geometry proxies, not solvent-accessible "
                "surface area, dynamics, epitope accessibility, or binding evidence."
            ),
        },
        "hotspot_proposals": epitope_rows,
        "review_gate": {
            "status": "human_review_required",
            "required_action": (
                "Inspect structure, confidence, local chemistry, exposure (preferably "
                "with an independent SASA calculation), and model limitations before "
                "copying any proposal into selected_hotspots_by_epitope."
            ),
            "whole_epitope_window_must_not_be_selected_as_hotspots": True,
        },
    }


def _build_candidate_manifest(
    *,
    target_pdb_path: Path,
    evidence_path: Path,
    evidence_sha256: str,
    expected_source_sha256: str,
) -> dict[str, Any]:
    target_positions = list(range(TARGET_START, TARGET_END + 1))
    proposed = {
        str(epitope["epitope_id"]): list(epitope["proposed_hotspots"])
        for epitope in EPITOPES
    }
    return {
        "schema": TARGET_MANIFEST_SCHEMA,
        "execution_state": "blocked_pending_human_review",
        "review": {
            "reviewed_by": "",
            "reviewed_at": "",
            "contracts_acknowledged": False,
        },
        "conformation_id": "AFDB-P07196-F1-v6-NEFL-280-377-monomer-candidate",
        "oligomeric_state": "single_chain_monomer",
        "target_pdb_path": _portable_path(target_pdb_path),
        "target_chain": SOURCE_CHAIN,
        "antigen_sequence_in_pdb_order": TARGET_SEQUENCE,
        "full_coordinate_to_pdb": {
            str(position): {
                "chain_id": SOURCE_CHAIN,
                "residue_number": position,
                "insertion_code": "",
            }
            for position in target_positions
        },
        "full_coordinate_to_local_1_based": {
            str(position): position - TARGET_START + 1 for position in target_positions
        },
        "selected_hotspots_by_epitope": {
            str(epitope["epitope_id"]): [] for epitope in EPITOPES
        },
        "proposed_hotspots_by_epitope_pending_human_review": proposed,
        "hotspot_review_state": "proposed_not_selected",
        "hotspot_evidence": {
            "path": _portable_path(evidence_path),
            "sha256": evidence_sha256,
        },
        "framework_coordinate_inputs": {
            "rfantibody_hlt_pdbs": {
                "template_7-H11-D3-2-C7": "",
                "template_15-C12-H6": "",
            },
            "germinal_scfv_pdbs": {
                "template_7-H11-D3-2-C7": "",
                "template_15-C12-H6": "",
            },
        },
        "candidate_provenance": {
            "source_database": "AlphaFold Protein Structure Database",
            "source_accession": SOURCE_ACCESSION,
            "source_model_version": "v6",
            "source_sha256": expected_source_sha256,
            "target_pdb_sha256": _file_sha256(target_pdb_path)
            if target_pdb_path.is_file()
            else None,
            "scientific_status": "theoretical_prediction_not_experimentally_determined",
        },
        "review_contract": [
            "The AlphaFold DB structure is a theoretical monomer prediction, not an experimentally determined structure.",
            "PDB sequence and residue numbering were extracted from coordinates and independently checked against SEQRES.",
            "Every selected hotspot must have an explicit full-coordinate to PDB-residue mapping and be present in ATOM records.",
            "Proposed hotspots must be reviewed for solvent accessibility and local chemistry; geometry proxies are not SASA.",
            "An entire epitope window must not be passed blindly as selected hotspots.",
            "RFantibody HLT files must contain H/L chains and all six CDR REMARK labels.",
            "Germinal templates must be independent single-chain VH-linker-VL scFv PDBs.",
            "Schema v1 represents one single-chain monomer; oligomeric targets require a separate multichain campaign.",
        ],
        "blocking_reasons": [
            "Human structure and hotspot review has not been recorded.",
            "selected_hotspots_by_epitope is intentionally empty.",
            "Two distinct RFantibody HLT coordinate templates are still required.",
            "Two distinct Germinal scFv coordinate templates are still required.",
        ],
    }


def _json_text(value: Any) -> str:
    return json.dumps(value, indent=2, ensure_ascii=False, sort_keys=False) + "\n"


def prepare_target_structure(
    source_pdb: str | Path,
    output_dir: str | Path,
    *,
    expected_source_sha256: str = SOURCE_SHA256,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Validate the pinned source and materialize a blocked review bundle."""

    source_path = Path(source_pdb).expanduser().resolve()
    output_root = Path(output_dir).expanduser().resolve()
    if not isinstance(expected_source_sha256, str) or len(expected_source_sha256) != 64:
        raise TargetStructurePreparationError("expected_source_sha256 must be 64 hex characters")
    try:
        int(expected_source_sha256, 16)
    except ValueError as exc:
        raise TargetStructurePreparationError(
            "expected_source_sha256 must contain only hexadecimal characters"
        ) from exc

    source_bytes, full_sequence, residues = _parse_and_validate_source(
        source_path, expected_sha256=expected_source_sha256.lower()
    )
    source_copy_path = output_root / SOURCE_FILENAME
    target_pdb_path = output_root / TARGET_PDB_FILENAME
    target_fasta_path = output_root / TARGET_FASTA_FILENAME
    evidence_path = output_root / EVIDENCE_FILENAME
    candidate_manifest_path = output_root / CANDIDATE_MANIFEST_FILENAME

    target_pdb_text = _cropped_pdb_text(
        residues, source_sha256=expected_source_sha256.lower()
    )
    target_fasta_text = (
        ">NEFL_P07196_280-377|source=AlphaFoldDB_v6|"
        "status=theoretical_prediction_candidate\n"
        f"{TARGET_SEQUENCE}\n"
    )
    evidence = _build_evidence(
        source_copy_path=source_copy_path,
        target_pdb_path=target_pdb_path,
        target_pdb_text=target_pdb_text,
        target_fasta_path=target_fasta_path,
        target_fasta_text=target_fasta_text,
        full_sequence=full_sequence,
        residues=residues,
        expected_source_sha256=expected_source_sha256.lower(),
    )
    evidence_text = _json_text(evidence)

    payloads: dict[Path, bytes] = {
        source_copy_path: source_bytes,
        target_pdb_path: target_pdb_text.encode("ascii"),
        target_fasta_path: target_fasta_text.encode("ascii"),
        evidence_path: evidence_text.encode("utf-8"),
    }
    conflicts = sorted(path for path in payloads if path.exists())
    if candidate_manifest_path.exists():
        conflicts.append(candidate_manifest_path)
    if conflicts and not overwrite:
        rendered = "\n  - ".join(str(path) for path in conflicts)
        raise TargetStructurePreparationError(
            "Refusing to overwrite existing structure artifacts; rerun with "
            f"--overwrite after review:\n  - {rendered}"
        )

    output_root.mkdir(parents=True, exist_ok=True)
    for path, payload in payloads.items():
        if path.resolve() == source_path and payload == source_bytes:
            continue
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=output_root, prefix=f".{path.name}.", delete=False
        ) as handle:
            handle.write(payload)
            temporary_path = Path(handle.name)
        temporary_path.replace(path)

    # The candidate manifest records the materialized cropped-PDB hash, so it is
    # constructed only after that immutable payload has been written.
    candidate_manifest = _build_candidate_manifest(
        target_pdb_path=target_pdb_path,
        evidence_path=evidence_path,
        evidence_sha256=_file_sha256(evidence_path),
        expected_source_sha256=expected_source_sha256.lower(),
    )
    candidate_manifest_text = _json_text(candidate_manifest)
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=output_root, prefix=f".{candidate_manifest_path.name}.", delete=False
    ) as handle:
        handle.write(candidate_manifest_text.encode("utf-8"))
        temporary_manifest = Path(handle.name)
    temporary_manifest.replace(candidate_manifest_path)

    return {
        "execution_state": "prepared_blocked_pending_human_review",
        "source_pdb": _portable_path(source_copy_path),
        "source_sha256": _file_sha256(source_copy_path),
        "target_pdb": _portable_path(target_pdb_path),
        "target_pdb_sha256": _file_sha256(target_pdb_path),
        "target_fasta": _portable_path(target_fasta_path),
        "evidence": _portable_path(evidence_path),
        "candidate_manifest": _portable_path(candidate_manifest_path),
        "residue_count": TARGET_END - TARGET_START + 1,
        "selected_hotspots": None,
        "proposed_hotspots_pending_human_review": {
            str(epitope["epitope_id"]): list(epitope["proposed_hotspots"])
            for epitope in EPITOPES
        },
        "real_model_handoff_authorized": False,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate pinned AlphaFold DB P07196 v6 coordinates and prepare a "
            "human-review-blocked NEFL 280-377 single-chain monomer bundle."
        )
    )
    parser.add_argument(
        "--source-pdb",
        default=str(REPO_ROOT / "input" / "structures" / SOURCE_FILENAME),
        help="Pinned full-length AF-P07196-F1-model_v6 PDB (no download is performed)",
    )
    parser.add_argument(
        "--output-dir",
        default=str(REPO_ROOT / "input" / "structures"),
        help="Directory for source copy, cropped PDB, evidence, FASTA, and candidate manifest",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing generated candidate bundle after review",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = prepare_target_structure(
            args.source_pdb,
            args.output_dir,
            overwrite=args.overwrite,
        )
    except TargetStructurePreparationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(_json_text(summary), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
