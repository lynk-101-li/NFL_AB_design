"""Generate a non-executing Germinal handoff from a normalized design request.

This adapter deliberately stops before running Germinal.  It validates local PDB
inputs, converts the repository's paired VH/VL template description into the
single-chain scFv geometry expected by Germinal, and emits target-YAML text plus
Hydra argument vectors.  Calling code remains responsible for reviewing and
staging the generated jobs in a separately installed Germinal checkout.

Upstream documentation: https://github.com/SantiagoMille/germinal
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence


NORMALIZED_REQUEST_SCHEMA = "nfl_ab_design.normalized_de_novo_request.v1"
GERMINAL_HANDOFF_SCHEMA = "nfl_ab_design.germinal_handoff.v1"
EXPECTED_CDR_ORDER = ("H1", "H2", "H3", "L1", "L2", "L3")
BACKEND_PATH_KEYS = frozenset(
    {
        "af3_repo_path",
        "af3_sif_path",
        "af3_model_dir",
        "af3_db_dir",
        "msa_db_dir",
    }
)
BACKEND_IDENTIFIER_KEYS = frozenset(
    {
        "protenix_conda_env",
        "protenix_model_name",
    }
)

# Germinal has no releases or tags at the time this integration was researched.
# Pinning a full commit makes the adapter contract auditable despite that fact.
UPSTREAM_PROVENANCE: dict[str, Any] = {
    "project": "Germinal",
    "public_repository": "https://github.com/SantiagoMille/germinal",
    "repository_owner": "SantiagoMille",
    "repository_note": (
        "The public paper implementation is SantiagoMille/germinal; "
        "RosettaCommons/Germinal is not a public GitHub repository."
    ),
    "pinned_commit": "1e1c1a5b79884ae45abae030c9df90d9423a990a",
    "pinned_commit_url": (
        "https://github.com/SantiagoMille/germinal/commit/"
        "1e1c1a5b79884ae45abae030c9df90d9423a990a"
    ),
    "last_user_validated_commit_reported_by_upstream": (
        "2c0a13b76833b6463cb59c571cfeadf17fd710c1"
    ),
    "package_version_in_setup_py": "0.0.1",
    "release_state": "no_tags_or_releases_observed",
    "research_date": "2026-08-12",
    "supported_binder_geometries": ["VHH", "single-chain scFv"],
    "native_paired_chain_fv_supported": False,
    "code_license": "Apache-2.0",
    "license_caveats": [
        "PyRosetta is a separately obtained dependency with non-commercial/non-profit terms; commercial use needs a separate license.",
        "IgLM is distributed under a separate non-commercial academic license.",
        "AlphaFold 3 source and model parameters have separate license/terms and are not required by the smoke profile.",
        "Every model, weight set, and dependency must be reviewed under its own current terms.",
    ],
    "official_source_urls": [
        "https://github.com/SantiagoMille/germinal",
        "https://github.com/SantiagoMille/germinal/blob/main/README.md",
        "https://github.com/SantiagoMille/germinal/blob/main/environment_setup.md",
        "https://github.com/SantiagoMille/germinal/blob/main/configs/run/scfv.yaml",
        "https://github.com/SantiagoMille/germinal/blob/main/configs/target/pdl1.yaml",
        "https://github.com/SantiagoMille/germinal/blob/main/LICENSE",
    ],
}


PROFILE_SETTINGS: dict[str, dict[str, Any]] = {
    "smoke": {
        "purpose": "Minimal integration smoke test; not a scientific design campaign.",
        "structure_backend": "chai",
        "scientific_use_allowed": False,
        "hydra_overrides": {
            "max_trajectories": 1,
            "max_hallucinated_trajectories": 1,
            "max_passing_designs": 1,
            "structure_model": "chai",
            "logits_steps": 1,
            "softmax_steps": 0,
            "search_steps": 0,
            "num_seqs": 1,
            "max_mpnn_sequences": 1,
            "multi_relax": False,
            "save_design_animations": False,
            "save_design_trajectory_plots": False,
        },
        "notes": [
            "Chai is selected so the smoke test does not require AlphaFold 3 weights or license acceptance.",
            "Germinal still requires its AlphaFold-Multimer hallucination parameters and imports PyRosetta.",
            "One optimization step is intended to expose installation, input, and shape errors; it is not a quality-producing run.",
        ],
    },
    "full": {
        "purpose": "Upstream-scale job specification requiring resource and license review before execution.",
        "structure_backend": "af3",
        "scientific_use_allowed": True,
        "hydra_overrides": {
            "max_trajectories": 10000,
            "max_hallucinated_trajectories": 1000,
            "max_passing_designs": 100,
            "structure_model": "af3",
            "af3_structure_select_mode": "worst",
            "multi_relax": False,
            "save_design_animations": False,
            "save_design_trajectory_plots": True,
        },
        "notes": [
            "Counts mirror Germinal's upstream top-level defaults and are expensive; review them per epitope before launch.",
            "AF3 is the upstream-recommended calibrated filter backend, but requires a separate installation, parameters, databases, container, and applicable license/terms.",
            "Germinal documents scFv configs and filters as preliminary/experimental.",
        ],
    },
}


AA3_TO_AA1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
AA20 = frozenset(AA3_TO_AA1.values())


class GerminalAdapterError(ValueError):
    """Raised when a request cannot be translated without ambiguity."""


@dataclass(frozen=True)
class PdbResidue:
    chain: str
    number: int
    insertion_code: str
    amino_acid: str

    @property
    def germinal_hotspot(self) -> str:
        if self.insertion_code:
            raise GerminalAdapterError(
                "Germinal hotspot syntax is not safely defined for PDB insertion "
                f"codes; cannot emit {self.chain}{self.number}{self.insertion_code}."
            )
        return f"{self.chain}{self.number}"

    @property
    def chai_hotspot(self) -> str:
        if self.insertion_code:
            raise GerminalAdapterError(
                "Chai restraint hotspot syntax is not safely defined for PDB "
                f"insertion codes; cannot emit {self.chain}{self.number}{self.insertion_code}."
            )
        return f"{self.amino_acid}{self.number}"


def load_normalized_request(path: str | Path) -> dict[str, Any]:
    """Load and minimally validate a normalized JSON design request."""

    request_path = Path(path).expanduser()
    if not request_path.is_file():
        raise GerminalAdapterError(f"Normalized request does not exist: {request_path}")
    try:
        value = json.loads(request_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GerminalAdapterError(
            f"Could not read normalized request JSON {request_path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise GerminalAdapterError("Normalized request root must be a JSON object.")
    _validate_request(value)
    return value


def _request_dict(request: Mapping[str, Any] | str | Path) -> dict[str, Any]:
    if isinstance(request, (str, Path)):
        return load_normalized_request(request)
    if not isinstance(request, Mapping):
        raise GerminalAdapterError("request must be a mapping or JSON path")
    value = dict(request)
    _validate_request(value)
    return value


def _require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise GerminalAdapterError(f"{label} must be an object/mapping.")
    return value


def _require_nonempty_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise GerminalAdapterError(f"{label} must be a non-empty string.")
    return value.strip()


def _validate_amino_acid_sequence(sequence: Any, label: str, *, masked: bool = False) -> str:
    text = _require_nonempty_string(sequence, label).upper()
    allowed = AA20 | ({"X"} if masked else set())
    invalid = sorted(set(text) - allowed)
    if invalid:
        raise GerminalAdapterError(
            f"{label} contains unsupported residues: {', '.join(invalid)}"
        )
    return text


def _validate_request(request: Mapping[str, Any]) -> None:
    schema = request.get("schema")
    if schema != NORMALIZED_REQUEST_SCHEMA:
        raise GerminalAdapterError(
            f"Unsupported normalized request schema {schema!r}; expected "
            f"{NORMALIZED_REQUEST_SCHEMA!r}."
        )

    antigen = _require_mapping(request.get("antigen"), "request.antigen")
    full_sequence = _validate_amino_acid_sequence(
        antigen.get("full_sequence"), "request.antigen.full_sequence"
    )
    _require_nonempty_string(antigen.get("protein"), "request.antigen.protein")

    epitopes = request.get("epitopes")
    if not isinstance(epitopes, list) or not epitopes:
        raise GerminalAdapterError("request.epitopes must be a non-empty list.")
    seen_epitopes: set[str] = set()
    for index, raw_epitope in enumerate(epitopes):
        epitope = _require_mapping(raw_epitope, f"request.epitopes[{index}]")
        epitope_id = _require_nonempty_string(
            epitope.get("epitope_id"), f"request.epitopes[{index}].epitope_id"
        )
        if epitope_id in seen_epitopes:
            raise GerminalAdapterError(f"Duplicate epitope_id: {epitope_id}")
        seen_epitopes.add(epitope_id)
        try:
            start = int(epitope["start_1_based"])
            end = int(epitope["end_1_based_inclusive"])
        except (KeyError, TypeError, ValueError) as exc:
            raise GerminalAdapterError(
                f"Epitope {epitope_id} requires integer 1-based start/end coordinates."
            ) from exc
        if start < 1 or end < start or end > len(full_sequence):
            raise GerminalAdapterError(
                f"Epitope {epitope_id} coordinates {start}-{end} are outside "
                f"the antigen sequence (length {len(full_sequence)})."
            )
        epitope_sequence = _validate_amino_acid_sequence(
            epitope.get("sequence"), f"epitope {epitope_id} sequence"
        )
        expected = full_sequence[start - 1 : end]
        if epitope_sequence != expected:
            raise GerminalAdapterError(
                f"Epitope {epitope_id} sequence does not match antigen "
                f"coordinates {start}-{end}: expected {expected}, got {epitope_sequence}."
            )
        raw_hotspots = epitope.get("candidate_hotspot_residue_indices")
        if not isinstance(raw_hotspots, list) or not raw_hotspots:
            raise GerminalAdapterError(
                f"Epitope {epitope_id} requires candidate_hotspot_residue_indices."
            )
        try:
            hotspots = [int(item) for item in raw_hotspots]
        except (TypeError, ValueError) as exc:
            raise GerminalAdapterError(
                f"Epitope {epitope_id} hotspot indices must be integers."
            ) from exc
        if len(set(hotspots)) != len(hotspots):
            raise GerminalAdapterError(f"Epitope {epitope_id} has duplicate hotspots.")
        outside = [position for position in hotspots if not start <= position <= end]
        if outside:
            raise GerminalAdapterError(
                f"Epitope {epitope_id} hotspots fall outside its coordinate range: {outside}"
            )

    templates = request.get("templates")
    if not isinstance(templates, list) or not templates:
        raise GerminalAdapterError("request.templates must be a non-empty list.")
    seen_templates: set[str] = set()
    for index, raw_template in enumerate(templates):
        template = _require_mapping(raw_template, f"request.templates[{index}]")
        template_id = _require_nonempty_string(
            template.get("template_id"), f"request.templates[{index}].template_id"
        )
        if template_id in seen_templates:
            raise GerminalAdapterError(f"Duplicate template_id: {template_id}")
        seen_templates.add(template_id)
        role = template.get("template_role")
        if role != "framework_source_only":
            raise GerminalAdapterError(
                f"Template {template_id} role must be 'framework_source_only', got {role!r}."
            )
        _validate_amino_acid_sequence(
            template.get("masked_vh"), f"template {template_id} masked_vh", masked=True
        )
        _validate_amino_acid_sequence(
            template.get("masked_vl"), f"template {template_id} masked_vl", masked=True
        )
        _validated_regions(template)


def _validated_regions(template: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    template_id = str(template.get("template_id", "<unknown>"))
    raw_regions = template.get("design_regions")
    if not isinstance(raw_regions, list):
        raise GerminalAdapterError(
            f"Template {template_id} design_regions must be a list."
        )
    regions: dict[str, dict[str, Any]] = {}
    for index, raw_region in enumerate(raw_regions):
        region = _require_mapping(
            raw_region, f"template {template_id} design_regions[{index}]"
        )
        name = _require_nonempty_string(
            region.get("region"), f"template {template_id} region name"
        )
        if name in regions:
            raise GerminalAdapterError(f"Template {template_id} repeats region {name}.")
        if name not in EXPECTED_CDR_ORDER:
            raise GerminalAdapterError(
                f"Template {template_id} has unsupported design region {name!r}."
            )
        expected_chain = "VH" if name.startswith("H") else "VL"
        if region.get("chain") != expected_chain:
            raise GerminalAdapterError(
                f"Template {template_id} region {name} must use chain {expected_chain}."
            )
        try:
            start = int(region["start_1_based"])
            end = int(region["end_1_based_inclusive"])
            recorded_length = int(region["length_aa"])
        except (KeyError, TypeError, ValueError) as exc:
            raise GerminalAdapterError(
                f"Template {template_id} region {name} requires integer coordinates and length."
            ) from exc
        chain_sequence = str(template["masked_vh" if expected_chain == "VH" else "masked_vl"])
        if start < 1 or end < start or end > len(chain_sequence):
            raise GerminalAdapterError(
                f"Template {template_id} region {name} coordinates {start}-{end} "
                f"are invalid for {expected_chain} length {len(chain_sequence)}."
            )
        if recorded_length != end - start + 1:
            raise GerminalAdapterError(
                f"Template {template_id} region {name} length does not match its coordinates."
            )
        if set(chain_sequence[start - 1 : end]) != {"X"}:
            raise GerminalAdapterError(
                f"Template {template_id} region {name} is not fully masked with X."
            )
        regions[name] = {
            "region": name,
            "chain": expected_chain,
            "start": start,
            "end": end,
            "length": recorded_length,
        }
    missing = [name for name in EXPECTED_CDR_ORDER if name not in regions]
    if missing:
        raise GerminalAdapterError(
            f"Template {template_id} is missing CDR design regions: {', '.join(missing)}"
        )
    for names in (("H1", "H2", "H3"), ("L1", "L2", "L3")):
        previous_end = 0
        for name in names:
            if regions[name]["start"] <= previous_end:
                raise GerminalAdapterError(
                    f"Template {template_id} CDR regions overlap or are out of order."
                )
            previous_end = int(regions[name]["end"])
    return regions


def _parse_pdb(path: str | Path, label: str) -> dict[str, list[PdbResidue]]:
    pdb_path = Path(path).expanduser().resolve()
    if not pdb_path.is_file():
        raise GerminalAdapterError(f"{label} PDB does not exist: {pdb_path}")
    if pdb_path.suffix.lower() not in {".pdb", ".ent"}:
        raise GerminalAdapterError(
            f"{label} must be a PDB file (.pdb/.ent), not {pdb_path.suffix or '<no suffix>'}."
        )
    chains: dict[str, list[PdbResidue]] = {}
    seen: dict[tuple[str, int, str], str] = {}
    model_count = 0
    active_model = True
    try:
        lines = pdb_path.read_text(encoding="utf-8", errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise GerminalAdapterError(f"Could not read {label} PDB {pdb_path}: {exc}") from exc
    for line_number, line in enumerate(lines, start=1):
        record = line[:6].strip().upper()
        if record == "MODEL":
            model_count += 1
            active_model = model_count == 1
            continue
        if record == "ENDMDL":
            active_model = False
            continue
        if record != "ATOM" or not active_model or len(line) < 27:
            continue
        chain = line[21].strip()
        if not chain:
            raise GerminalAdapterError(
                f"{label} PDB has a blank chain ID at line {line_number}."
            )
        try:
            number = int(line[22:26].strip())
        except ValueError as exc:
            raise GerminalAdapterError(
                f"{label} PDB has an invalid residue number at line {line_number}."
            ) from exc
        insertion_code = line[26].strip()
        residue_name = line[17:20].strip().upper()
        if residue_name not in AA3_TO_AA1:
            raise GerminalAdapterError(
                f"{label} PDB contains unsupported ATOM residue {residue_name!r} "
                f"at {chain}{number}{insertion_code}."
            )
        key = (chain, number, insertion_code)
        amino_acid = AA3_TO_AA1[residue_name]
        if key in seen:
            if seen[key] != amino_acid:
                raise GerminalAdapterError(
                    f"{label} PDB assigns conflicting residue types to {key}."
                )
            continue
        seen[key] = amino_acid
        chains.setdefault(chain, []).append(
            PdbResidue(chain, number, insertion_code, amino_acid)
        )
    if model_count > 1:
        raise GerminalAdapterError(
            f"{label} PDB contains {model_count} MODEL records; provide one model only."
        )
    if not chains:
        raise GerminalAdapterError(f"{label} PDB contains no supported ATOM residues.")
    return chains


def _normalize_chain_id(value: str, label: str) -> str:
    chain = _require_nonempty_string(value, label)
    if len(chain) != 1 or not chain.isalnum():
        raise GerminalAdapterError(
            f"{label} must be one alphanumeric PDB chain character, got {chain!r}."
        )
    return chain


def _normalize_residue_ref(
    value: Any, *, default_chain: str, source_position: int
) -> tuple[str, int, str]:
    chain = default_chain
    insertion_code = ""
    number: Any = None
    if isinstance(value, bool):
        number = None
    elif isinstance(value, int):
        number = value
    elif isinstance(value, str):
        match = re.fullmatch(r"(?:(?P<chain>[A-Za-z0-9]):?)?(?P<number>-?\d+)(?P<icode>[A-Za-z]?)", value.strip())
        if match:
            chain = match.group("chain") or default_chain
            number = int(match.group("number"))
            insertion_code = match.group("icode") or ""
    elif isinstance(value, Mapping):
        chain = str(value.get("chain", default_chain))
        number = value.get("number", value.get("resseq", value.get("residue_number")))
        insertion_code = str(value.get("insertion_code", value.get("icode", "")))
    elif isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        if len(value) in (2, 3):
            chain = str(value[0])
            number = value[1]
            insertion_code = str(value[2]) if len(value) == 3 else ""
    try:
        number = int(number)
    except (TypeError, ValueError) as exc:
        raise GerminalAdapterError(
            f"Invalid PDB residue mapping for antigen position {source_position}: {value!r}"
        ) from exc
    chain = _normalize_chain_id(chain, f"residue-map chain for {source_position}")
    if len(insertion_code) > 1 or (insertion_code and not insertion_code.isalpha()):
        raise GerminalAdapterError(
            f"Invalid insertion code for antigen position {source_position}: {insertion_code!r}"
        )
    return chain, number, insertion_code


def _map_epitope_hotspots(
    epitope: Mapping[str, Any],
    *,
    full_sequence: str,
    target_residues: Mapping[tuple[str, int, str], PdbResidue],
    target_chain: str,
    target_residue_map: Mapping[Any, Any] | None,
) -> list[tuple[int, PdbResidue]]:
    epitope_id = str(epitope["epitope_id"])
    source_positions = _selected_germinal_hotspots(epitope)
    normalized_map: dict[int, Any] | None = None
    if target_residue_map is not None:
        normalized_map = {}
        for raw_key, mapped_value in target_residue_map.items():
            try:
                key = int(raw_key)
            except (TypeError, ValueError) as exc:
                raise GerminalAdapterError(
                    f"Target residue-map key must be an antigen integer position: {raw_key!r}"
                ) from exc
            if key in normalized_map:
                raise GerminalAdapterError(f"Target residue map repeats source position {key}.")
            normalized_map[key] = mapped_value
        missing = [position for position in source_positions if position not in normalized_map]
        if missing:
            raise GerminalAdapterError(
                f"Explicit target residue map is missing {epitope_id} hotspot positions: {missing}"
            )

    mapped: list[tuple[int, PdbResidue]] = []
    seen_pdb_refs: set[tuple[str, int, str]] = set()
    for source_position in source_positions:
        raw_ref: Any = source_position
        if normalized_map is not None:
            raw_ref = normalized_map[source_position]
        ref = _normalize_residue_ref(
            raw_ref, default_chain=target_chain, source_position=source_position
        )
        if ref[0] != target_chain:
            raise GerminalAdapterError(
                f"Epitope {epitope_id} maps antigen position {source_position} to "
                f"chain {ref[0]}, but this single-chain adapter target is {target_chain}."
            )
        residue = target_residues.get(ref)
        if residue is None:
            suffix = ref[2]
            raise GerminalAdapterError(
                f"Epitope {epitope_id} maps antigen position {source_position} to "
                f"missing PDB residue {ref[0]}{ref[1]}{suffix}. Provide a correct "
                "target_residue_map when target PDB numbering differs from NEFL numbering."
            )
        expected_aa = full_sequence[source_position - 1]
        if residue.amino_acid != expected_aa:
            raise GerminalAdapterError(
                f"Target mapping mismatch for {epitope_id}: antigen {expected_aa}{source_position} "
                f"maps to PDB {residue.amino_acid}{residue.number}{residue.insertion_code}."
            )
        pdb_key = (residue.chain, residue.number, residue.insertion_code)
        if pdb_key in seen_pdb_refs:
            raise GerminalAdapterError(
                f"Epitope {epitope_id} maps multiple antigen positions to PDB residue {pdb_key}."
            )
        seen_pdb_refs.add(pdb_key)
        # Accessing this property rejects unsupported insertion-code hotspots now.
        residue.germinal_hotspot
        mapped.append((source_position, residue))
    return mapped


def _selected_germinal_hotspots(epitope: Mapping[str, Any]) -> list[int]:
    """Return curated hotspots, never the full candidate window by default."""

    epitope_id = str(epitope.get("epitope_id", "<unknown>"))
    selected_key = next(
        (
            key
            for key in (
                "germinal_hotspot_residue_indices",
                "selected_hotspot_residue_indices",
            )
            if key in epitope
        ),
        None,
    )
    if selected_key is None:
        raise GerminalAdapterError(
            f"Epitope {epitope_id} has only a candidate hotspot window. A real "
            "Germinal job requires a curated non-empty "
            "germinal_hotspot_residue_indices (preferred) or "
            "selected_hotspot_residue_indices list; refusing to treat every "
            "window residue as an experimentally justified hotspot."
        )
    raw_selected = epitope.get(selected_key)
    if not isinstance(raw_selected, list) or not raw_selected:
        raise GerminalAdapterError(
            f"Epitope {epitope_id} {selected_key} must be a non-empty list."
        )
    if any(isinstance(item, bool) or not isinstance(item, int) for item in raw_selected):
        raise GerminalAdapterError(
            f"Epitope {epitope_id} {selected_key} must contain JSON integers only."
        )
    selected = list(raw_selected)
    candidates = set(epitope["candidate_hotspot_residue_indices"])
    if len(set(selected)) != len(selected):
        raise GerminalAdapterError(
            f"Epitope {epitope_id} {selected_key} contains duplicate positions."
        )
    outside = [position for position in selected if position not in candidates]
    if outside:
        raise GerminalAdapterError(
            f"Epitope {epitope_id} curated Germinal hotspots are not contained "
            f"in candidate_hotspot_residue_indices: {outside}"
        )
    return selected


def _resolve_template_pdbs(
    templates: Sequence[Mapping[str, Any]], supplied: Mapping[str, str | Path]
) -> dict[str, Path]:
    if not isinstance(supplied, Mapping) or not supplied:
        raise GerminalAdapterError(
            "template_scfv_pdbs must map every template_id (or framework_source_id) "
            "to a per-template scFv PDB."
        )
    resolved: dict[str, Path] = {}
    used_keys: set[str] = set()
    for template in templates:
        template_id = str(template["template_id"])
        source_id = str(template.get("framework_source_id", ""))
        matching_keys = [key for key in (template_id, source_id) if key and key in supplied]
        if not matching_keys:
            raise GerminalAdapterError(
                f"Missing per-template scFv PDB for {template_id}; accepted mapping keys "
                f"are {template_id!r} or {source_id!r}."
            )
        unique_values = {str(Path(supplied[key]).expanduser()) for key in matching_keys}
        if len(unique_values) != 1:
            raise GerminalAdapterError(
                f"Conflicting scFv PDB paths supplied for template {template_id}."
            )
        used_keys.update(matching_keys)
        path = Path(next(iter(unique_values))).resolve()
        if not path.is_file():
            raise GerminalAdapterError(
                f"Per-template scFv PDB for {template_id} does not exist: {path}"
            )
        resolved[template_id] = path
    unexpected = sorted(str(key) for key in supplied if str(key) not in used_keys)
    if unexpected:
        raise GerminalAdapterError(
            f"Unknown template_scfv_pdbs keys (possible typo): {', '.join(unexpected)}"
        )
    return resolved


def _scfv_layout(
    template: Mapping[str, Any], pdb_path: Path, *, scfv_chain: str
) -> dict[str, Any]:
    template_id = str(template["template_id"])
    chains = _parse_pdb(pdb_path, f"template {template_id} scFv")
    if scfv_chain != "A":
        raise GerminalAdapterError(
            "Germinal create_starting_structure reads the template binder from chain A; "
            f"this adapter therefore requires scfv_chain='A', got {scfv_chain!r}."
        )
    if scfv_chain not in chains:
        raise GerminalAdapterError(
            f"Template {template_id} scFv PDB lacks required chain {scfv_chain}; "
            f"found {sorted(chains)}."
        )
    if len(chains) != 1:
        raise GerminalAdapterError(
            f"Template {template_id} scFv PDB must contain exactly one chain; "
            f"found {sorted(chains)}."
        )
    pdb_sequence = "".join(residue.amino_acid for residue in chains[scfv_chain])
    masked_vh = str(template["masked_vh"])
    masked_vl = str(template["masked_vl"])
    linker_length = len(pdb_sequence) - len(masked_vh) - len(masked_vl)
    if linker_length < 1:
        raise GerminalAdapterError(
            f"Template {template_id} PDB chain length {len(pdb_sequence)} cannot encode "
            f"VH ({len(masked_vh)}) + linker + VL ({len(masked_vl)})."
        )
    pdb_vh = pdb_sequence[: len(masked_vh)]
    pdb_linker = pdb_sequence[len(masked_vh) : len(masked_vh) + linker_length]
    pdb_vl = pdb_sequence[-len(masked_vl) :]
    for chain_name, masked, observed in (
        ("VH", masked_vh, pdb_vh),
        ("VL", masked_vl, pdb_vl),
    ):
        mismatch = next(
            (
                index
                for index, (expected, actual) in enumerate(zip(masked, observed), start=1)
                if expected != "X" and expected != actual
            ),
            None,
        )
        if mismatch is not None:
            raise GerminalAdapterError(
                f"Template {template_id} scFv PDB does not preserve masked-request "
                f"framework sequence: {chain_name}{mismatch} expected {masked[mismatch - 1]}, "
                f"observed {observed[mismatch - 1]}."
            )

    regions = _validated_regions(template)
    cdr_lengths = [int(regions[name]["length"]) for name in EXPECTED_CDR_ORDER]
    framework_lengths = [
        int(regions["H1"]["start"]) - 1,
        int(regions["H2"]["start"]) - int(regions["H1"]["end"]) - 1,
        int(regions["H3"]["start"]) - int(regions["H2"]["end"]) - 1,
        len(masked_vh)
        - int(regions["H3"]["end"])
        + linker_length
        + int(regions["L1"]["start"])
        - 1,
        int(regions["L2"]["start"]) - int(regions["L1"]["end"]) - 1,
        int(regions["L3"]["start"]) - int(regions["L2"]["end"]) - 1,
        len(masked_vl) - int(regions["L3"]["end"]),
    ]
    if any(length < 0 for length in framework_lengths):
        raise GerminalAdapterError(
            f"Template {template_id} produces a negative Germinal framework segment: "
            f"{framework_lengths}"
        )
    if sum(cdr_lengths) + sum(framework_lengths) != len(pdb_sequence):
        raise GerminalAdapterError(
            f"Template {template_id} Germinal length decomposition does not match its PDB."
        )
    return {
        "scfv_pdb_path": str(pdb_path),
        "scfv_pdb_chain": scfv_chain,
        "scfv_pdb_length": len(pdb_sequence),
        "vh_first": True,
        "vh_length": len(masked_vh),
        "linker_length": linker_length,
        "linker_sequence": pdb_linker,
        "vl_length": len(masked_vl),
        "cdr_order": list(EXPECTED_CDR_ORDER),
        "cdr_lengths": cdr_lengths,
        "framework_lengths": framework_lengths,
        "framework_sequence_validated": True,
        "cdr_seed_identity_validated_as_non_control": False,
    }


def _slug(value: str, *, max_length: int = 80) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()
    if not slug:
        raise GerminalAdapterError(f"Cannot form a safe identifier from {value!r}.")
    return slug[:max_length].rstrip("_")


def _yaml_scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    if isinstance(value, (int, float)):
        return str(value)
    return json.dumps(str(value), ensure_ascii=False)


def _target_yaml(values: Mapping[str, Any]) -> str:
    ordered_keys = (
        "target_name",
        "target_pdb_path",
        "target_chain",
        "binder_chain",
        "target_hotspots",
        "hotspot_residue",
        "dimer",
    )
    lines = [
        "# Generated by nfl_ab_design.adapters.germinal; review before staging.",
        "# Source antigen coordinates were validated against the supplied target PDB.",
    ]
    for key in ordered_keys:
        if key in values:
            lines.append(f"{key}: {_yaml_scalar(values[key])}")
    return "\n".join(lines) + "\n"


def _hydra_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(_hydra_value(item) for item in value) + "]"
    return str(value)


def _validated_optional_directory(path: str | Path | None, label: str) -> str | None:
    if path is None:
        return None
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_dir():
        raise GerminalAdapterError(f"{label} directory does not exist: {resolved}")
    return str(resolved)


def build_germinal_jobs(
    request: Mapping[str, Any] | str | Path,
    *,
    target_pdb_path: str | Path,
    template_scfv_pdbs: Mapping[str, str | Path],
    target_residue_map: Mapping[Any, Any] | None = None,
    handoff_root: str | Path = "germinal_handoff",
    profile: str = "smoke",
    target_chain: str = "A",
    binder_chain: str = "B",
    scfv_chain: str = "A",
    germinal_repo_dir: str | Path | None = None,
    af_params_dir: str | Path | None = None,
    backend_paths: Mapping[str, str | Path] | None = None,
    profile_overrides: Mapping[str, Any] | None = None,
    python_executable: str = "python",
) -> dict[str, Any]:
    """Build, but never execute, Germinal jobs for every template/epitope pair.

    ``target_residue_map`` maps normalized antigen positions to target-PDB
    residues. Values may be integers (identity chain), strings such as ``A322``,
    two/three-item sequences, or mappings with ``chain``, ``number``/``resseq``,
    and optional ``insertion_code``. If omitted, identity numbering on
    ``target_chain`` is required and validated residue-by-residue.
    """

    normalized = _request_dict(request)
    if profile not in PROFILE_SETTINGS:
        raise GerminalAdapterError(
            f"Unknown Germinal profile {profile!r}; choose one of {sorted(PROFILE_SETTINGS)}."
        )
    target_chain = _normalize_chain_id(target_chain, "target_chain")
    binder_chain = _normalize_chain_id(binder_chain, "binder_chain")
    scfv_chain = _normalize_chain_id(scfv_chain, "scfv_chain")
    if target_chain == binder_chain:
        raise GerminalAdapterError("target_chain and binder_chain must differ.")

    target_pdb = Path(target_pdb_path).expanduser().resolve()
    target_chains = _parse_pdb(target_pdb, "target")
    if len(target_chains) != 1:
        raise GerminalAdapterError(
            "This NFL Germinal adapter supports one monomeric target protein "
            "chain only and always emits target dimer: false; target PDB contains "
            f"chains {sorted(target_chains)}. Extract the reviewed antigen chain "
            "into a separate PDB instead of silently discarding other chains."
        )
    if target_chain not in target_chains:
        raise GerminalAdapterError(
            f"Target PDB lacks chain {target_chain}; found {sorted(target_chains)}."
        )
    target_lookup = {
        (residue.chain, residue.number, residue.insertion_code): residue
        for residue in target_chains[target_chain]
    }

    templates = [dict(item) for item in normalized["templates"]]
    epitopes = [dict(item) for item in normalized["epitopes"]]
    resolved_template_pdbs = _resolve_template_pdbs(templates, template_scfv_pdbs)
    scfv_layouts = {
        str(template["template_id"]): _scfv_layout(
            template,
            resolved_template_pdbs[str(template["template_id"])],
            scfv_chain=scfv_chain,
        )
        for template in templates
    }

    full_sequence = str(normalized["antigen"]["full_sequence"])
    mapped_epitopes: dict[str, list[tuple[int, PdbResidue]]] = {}
    for epitope in epitopes:
        epitope_id = str(epitope["epitope_id"])
        mapped_epitopes[epitope_id] = _map_epitope_hotspots(
            epitope,
            full_sequence=full_sequence,
            target_residues=target_lookup,
            target_chain=target_chain,
            target_residue_map=target_residue_map,
        )

    handoff_path = Path(handoff_root).expanduser().resolve()
    germinal_root = _validated_optional_directory(germinal_repo_dir, "Germinal repository")
    if germinal_root is not None and not (Path(germinal_root) / "run_germinal.py").is_file():
        raise GerminalAdapterError(
            f"Germinal repository lacks run_germinal.py: {germinal_root}"
        )
    af_parameter_path = _validated_optional_directory(
        af_params_dir, "AlphaFold-Multimer parameters"
    )
    backend_paths = dict(backend_paths or {})
    normalized_backend_paths: dict[str, str] = {}
    unknown_backend_keys = sorted(
        str(key)
        for key in backend_paths
        if str(key) not in BACKEND_PATH_KEYS | BACKEND_IDENTIFIER_KEYS
    )
    if unknown_backend_keys:
        raise GerminalAdapterError(
            "Unknown backend_paths keys (possible typo): "
            + ", ".join(unknown_backend_keys)
        )
    for raw_key, raw_value in backend_paths.items():
        key = str(raw_key)
        if key in BACKEND_PATH_KEYS:
            resolved = Path(raw_value).expanduser().resolve()
            if not resolved.exists():
                raise GerminalAdapterError(
                    f"backend_paths[{key!r}] does not exist: {resolved}"
                )
            normalized_backend_paths[key] = str(resolved)
        else:
            normalized_backend_paths[key] = _require_nonempty_string(
                raw_value, f"backend_paths[{key!r}]"
            )

    profile_spec = json.loads(json.dumps(PROFILE_SETTINGS[profile]))
    settings = dict(profile_spec["hydra_overrides"])
    if profile_overrides:
        forbidden = {
            "type",
            "cdr_lengths",
            "fw_lengths",
            "vh_first",
            "vh_len",
            "vl_len",
            "pdb_dir",
            "target",
            "run",
        }
        conflicts = sorted(forbidden.intersection(profile_overrides))
        if conflicts:
            raise GerminalAdapterError(
                "profile_overrides may not replace adapter-critical fields: "
                + ", ".join(conflicts)
            )
        settings.update(profile_overrides)
    backend = str(settings.get("structure_model", profile_spec["structure_backend"]))
    if backend not in {"chai", "af3", "protenix"}:
        raise GerminalAdapterError(
            f"Unsupported Germinal structure backend {backend!r}; use chai, af3, or protenix."
        )
    if profile == "smoke" and backend != "chai":
        raise GerminalAdapterError(
            "The smoke profile is intentionally pinned to Chai to avoid AF3 "
            "weights/license requirements; use profile='full' for another backend."
        )

    jobs: list[dict[str, Any]] = []
    for template in templates:
        template_id = str(template["template_id"])
        layout = scfv_layouts[template_id]
        for epitope in epitopes:
            epitope_id = str(epitope["epitope_id"])
            mapped = mapped_epitopes[epitope_id]
            job_id = _slug(f"germinal_{template_id}_{epitope_id}", max_length=120)
            target_config_name = _slug(f"nfl_{template_id}_{epitope_id}")
            job_workspace = handoff_path / "jobs" / job_id
            staged_pdb_dir = job_workspace / "pdbs"
            staged_scfv_pdb = staged_pdb_dir / "scfv.pdb"
            target_yaml_file = handoff_path / "target_configs" / f"{target_config_name}.yaml"
            installed_target_yaml = (
                Path(germinal_root) / "configs" / "target" / f"{target_config_name}.yaml"
                if germinal_root
                else Path("<GERMINAL_REPO_DIR>")
                / "configs"
                / "target"
                / f"{target_config_name}.yaml"
            )
            median_residue = mapped[len(mapped) // 2][1]
            target_values = {
                "target_name": target_config_name,
                "target_pdb_path": str(target_pdb),
                "target_chain": target_chain,
                "binder_chain": binder_chain,
                "target_hotspots": ",".join(
                    residue.germinal_hotspot for _, residue in mapped
                ),
                "hotspot_residue": median_residue.chai_hotspot,
                "dimer": False,
            }

            common_overrides: dict[str, Any] = {
                "project_dir": str(handoff_path),
                "results_dir": "results",
                "experiment_name": "nfl_germinal",
                "run_config": job_id,
                "pdb_dir": str(staged_pdb_dir),
                "af_params_dir": af_parameter_path or "<GERMINAL_AF_MULTIMER_PARAMS_DIR>",
                "type": "scfv",
                "cdr_lengths": layout["cdr_lengths"],
                "fw_lengths": layout["framework_lengths"],
                "vh_first": True,
                "vh_len": layout["vh_length"],
                "vl_len": layout["vl_length"],
                **settings,
            }
            if backend == "af3":
                common_overrides.update(
                    {
                        "af3_repo_path": normalized_backend_paths.get(
                            "af3_repo_path", "<AF3_REPO_PATH>"
                        ),
                        "af3_sif_path": normalized_backend_paths.get(
                            "af3_sif_path", "<AF3_SIF_PATH>"
                        ),
                        "af3_model_dir": normalized_backend_paths.get(
                            "af3_model_dir", "<AF3_MODEL_DIR>"
                        ),
                        "af3_db_dir": normalized_backend_paths.get(
                            "af3_db_dir", "<AF3_DB_DIR>"
                        ),
                        "msa_db_dir": normalized_backend_paths.get(
                            "msa_db_dir", "<COLABFOLD_MSA_DB_DIR>"
                        ),
                    }
                )
            elif backend == "protenix":
                common_overrides.update(
                    {
                        "protenix_conda_env": normalized_backend_paths.get(
                            "protenix_conda_env", "<PROTENIX_CONDA_ENV>"
                        ),
                        "protenix_model_name": normalized_backend_paths.get(
                            "protenix_model_name", "<PROTENIX_MODEL_NAME>"
                        ),
                    }
                )

            argv = [
                python_executable,
                "run_germinal.py",
                "run=scfv",
                "filter/initial=scfv",
                "filter/final=scfv",
                f"target={target_config_name}",
            ]
            argv.extend(
                f"{key}={_hydra_value(value)}" for key, value in common_overrides.items()
            )
            placeholders = sorted(
                {
                    match.group(0)
                    for argument in argv
                    for match in re.finditer(r"<[A-Z0-9_]+>", argument)
                }
            )
            expected_output_dir = (
                handoff_path / "results" / "nfl_germinal" / job_id
            )
            jobs.append(
                {
                    "job_id": job_id,
                    "execution_state": "not_run",
                    "ready_for_execution": not placeholders and germinal_root is not None,
                    "template_id": template_id,
                    "framework_source_id": template.get("framework_source_id", ""),
                    "epitope_id": epitope_id,
                    "profile": profile,
                    "structure_backend": backend,
                    "geometry_conversion": {
                        "source_geometry": "paired_VH_VL_Fv",
                        "germinal_geometry": "single_chain_VH_linker_VL_scFv",
                        "native_paired_fv_geometry_preserved": False,
                        "linker_sequence": layout["linker_sequence"],
                        "linker_length": layout["linker_length"],
                        "limitation": (
                            "Germinal does not natively design a two-chain paired Fv. "
                            "The supplied VH and VL frameworks are represented as one "
                            "VH-linker-VL chain, which changes inter-domain/linker geometry. "
                            "Selected sequences must be rebuilt and revalidated as a true "
                            "paired Fv/Fab before scientific interpretation."
                        ),
                        "cdr_seed_warning": (
                            "A coordinate-complete scFv PDB necessarily contains residues "
                            "inside all six CDR loops. This adapter validates framework "
                            "identity and CDR coordinates but cannot prove those seed residues "
                            "are independent of known positive-control CDRs. Use a neutral or "
                            "independently generated loop seed and record its provenance."
                        ),
                    },
                    "inputs": {
                        "normalized_request_schema": normalized["schema"],
                        "normalized_request_engine": normalized.get("engine", ""),
                        "target_pdb_path": str(target_pdb),
                        "template_scfv_pdb_path": layout["scfv_pdb_path"],
                        "template_scfv_source_chain": scfv_chain,
                        "combined_complex_binder_chain": binder_chain,
                        "scfv_layout": layout,
                    },
                    "target_mapping": {
                        "coordinate_source": "normalized antigen 1-based positions",
                        "pdb_chain": target_chain,
                        "target_geometry": "single_chain_monomer",
                        "target_pdb_chain_count": 1,
                        "germinal_dimer": False,
                        "identity_mapping_used": target_residue_map is None,
                        "curated_hotspot_source": next(
                            key
                            for key in (
                                "germinal_hotspot_residue_indices",
                                "selected_hotspot_residue_indices",
                            )
                            if key in epitope
                        ),
                        "mapped_hotspots": [
                            {
                                "antigen_position": source_position,
                                "antigen_amino_acid": full_sequence[source_position - 1],
                                "pdb_chain": residue.chain,
                                "pdb_residue_number": residue.number,
                                "pdb_insertion_code": residue.insertion_code,
                                "pdb_amino_acid": residue.amino_acid,
                                "germinal_hotspot": residue.germinal_hotspot,
                            }
                            for source_position, residue in mapped
                        ],
                    },
                    "target_yaml": {
                        "config_name": target_config_name,
                        "handoff_path": str(target_yaml_file),
                        "germinal_install_path": str(installed_target_yaml),
                        "values": target_values,
                        "content": _target_yaml(target_values),
                    },
                    "staging_plan": [
                        {
                            "action": "copy_file",
                            "source": layout["scfv_pdb_path"],
                            "destination": str(staged_scfv_pdb),
                            "reason": (
                                "Germinal hard-codes the scFv template filename as "
                                "<pdb_dir>/scfv.pdb; each template therefore needs an "
                                "isolated pdb_dir."
                            ),
                        },
                        {
                            "action": "copy_generated_target_yaml",
                            "source": str(target_yaml_file),
                            "destination": str(installed_target_yaml),
                        },
                    ],
                    "hydra_job": {
                        "working_directory": germinal_root or "<GERMINAL_REPO_DIR>",
                        "config_groups": {
                            "run": "scfv",
                            "filter/initial": "scfv",
                            "filter/final": "scfv",
                            "target": target_config_name,
                        },
                        "overrides": common_overrides,
                        "argv": argv,
                        "command_preview": shlex.join(argv),
                        "unresolved_placeholders": placeholders,
                    },
                    "expected_outputs": {
                        "run_directory": str(expected_output_dir),
                        "final_config": str(expected_output_dir / "final_config.yaml"),
                        "all_trajectories": str(expected_output_dir / "all_trajectories.csv"),
                        "failure_counts": str(expected_output_dir / "failure_counts.csv"),
                        "accepted_designs": str(
                            expected_output_dir / "accepted" / "designs.csv"
                        ),
                        "accepted_structures_glob": str(
                            expected_output_dir / "accepted" / "structures" / "*.pdb"
                        ),
                    },
                }
            )

    return {
        "schema": GERMINAL_HANDOFF_SCHEMA,
        "adapter": "nfl_ab_design.adapters.germinal",
        "execution_state": "not_run",
        "does_not_execute_external_code": True,
        "profile": profile,
        "profile_provenance": profile_spec,
        "upstream": json.loads(json.dumps(UPSTREAM_PROVENANCE)),
        "input_request": {
            "schema": normalized["schema"],
            "campaign_mode": normalized.get("campaign_mode", ""),
            "engine": normalized.get("engine", ""),
            "result_provenance": normalized.get("result_provenance", ""),
        },
        "resource_requirements": {
            "gpu": "NVIDIA CUDA GPU; upstream recommends CUDA 12+",
            "vram": "40 GB minimum documented; 60 GB+ recommended for larger runs",
            "storage": "50 GB+ recommended",
            "required_for_all_profiles": [
                "Pinned Germinal checkout",
                "JAX with GPU support",
                "AlphaFold-Multimer parameters for ColabDesign hallucination",
                "PyRosetta installation and applicable license",
                "AbMPNN weight bundled by upstream",
            ],
            "smoke_backend": [
                "chai-lab 0.6.1 as pinned by upstream environment instructions",
                "Chai model cache/download availability",
            ],
            "full_af3_backend": [
                "AlphaFold 3 repository and Singularity image",
                "AlphaFold 3 model parameters and databases",
                "ColabFold MSA databases when configured",
                "Acceptance of all applicable AlphaFold 3 terms",
            ],
        },
        "global_limitations": [
            "No Germinal, Chai, AF3, Protenix, AbMPNN, ColabDesign, or PyRosetta process was run by this adapter.",
            "Germinal's scFv path is preliminary/experimental according to upstream documentation.",
            "The smoke profile tests integration only and must not be interpreted as a scientific result.",
            "A generated scFv is not geometrically equivalent to the native two-chain Fv/Fab required by the NfL campaign.",
            "This adapter accepts exactly one target PDB chain and emits dimer: false; target dimers, oligomers, and multi-chain antigens are out of scope.",
        ],
        "job_count": len(jobs),
        "template_count": len(templates),
        "epitope_count": len(epitopes),
        "jobs": jobs,
    }


def build_germinal_handoff(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility alias for :func:`build_germinal_jobs`."""

    return build_germinal_jobs(*args, **kwargs)


def write_germinal_handoff(
    handoff: Mapping[str, Any],
    output_dir: str | Path,
    *,
    stage_template_pdbs: bool = False,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Materialize YAML/JSON handoff artifacts without executing Germinal.

    When ``stage_template_pdbs`` is true, the validated source scFv files are
    copied to the isolated ``<job>/pdbs/scfv.pdb`` locations described by each
    job.  Generated YAML is *not* installed into a Germinal checkout; that
    deliberate final staging step remains visible in each job specification.
    """

    if handoff.get("schema") != GERMINAL_HANDOFF_SCHEMA:
        raise GerminalAdapterError(
            f"Unsupported handoff schema: {handoff.get('schema')!r}"
        )
    jobs = handoff.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        raise GerminalAdapterError("Handoff contains no jobs to write.")
    root = Path(output_dir).expanduser().resolve()
    manifest_path = root / "germinal_jobs.json"
    yaml_dir = root / "target_configs"
    if not overwrite:
        existing = [manifest_path]
        existing.extend(
            yaml_dir / f"{job['target_yaml']['config_name']}.yaml" for job in jobs
        )
        conflicts = [path for path in existing if path.exists()]
        if conflicts:
            raise GerminalAdapterError(
                "Refusing to overwrite existing handoff artifacts: "
                + ", ".join(str(path) for path in conflicts)
            )
    yaml_dir.mkdir(parents=True, exist_ok=True)
    written_yamls: list[str] = []
    staged_pdbs: list[str] = []
    for job in jobs:
        yaml_path = yaml_dir / f"{job['target_yaml']['config_name']}.yaml"
        yaml_path.write_text(str(job["target_yaml"]["content"]), encoding="utf-8")
        written_yamls.append(str(yaml_path))
        if stage_template_pdbs:
            destination = root / "jobs" / str(job["job_id"]) / "pdbs" / "scfv.pdb"
            if destination.exists() and not overwrite:
                raise GerminalAdapterError(
                    f"Refusing to overwrite staged template PDB: {destination}"
                )
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(job["inputs"]["template_scfv_pdb_path"], destination)
            staged_pdbs.append(str(destination))
    root.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(handoff, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "execution_state": "not_run",
        "manifest": str(manifest_path),
        "target_yamls": written_yamls,
        "staged_template_pdbs": staged_pdbs,
    }


def _parse_key_value(items: Sequence[str], label: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise GerminalAdapterError(f"{label} must use KEY=VALUE syntax: {item!r}")
        key, value = item.split("=", 1)
        if not key or not value:
            raise GerminalAdapterError(f"{label} must use non-empty KEY=VALUE: {item!r}")
        if key in parsed:
            raise GerminalAdapterError(f"{label} repeats key {key!r}.")
        parsed[key] = value
    return parsed


def main(argv: Sequence[str] | None = None) -> int:
    """CLI that writes reviewed handoff artifacts and never launches Germinal."""

    parser = argparse.ArgumentParser(
        description="Build a non-executing Germinal scFv handoff from a normalized request."
    )
    parser.add_argument("request", help="Normalized design request JSON")
    parser.add_argument("--target-pdb", required=True)
    parser.add_argument(
        "--template-scfv-pdb",
        action="append",
        default=[],
        metavar="TEMPLATE_ID=PATH",
        help="Repeat once per normalized template.",
    )
    parser.add_argument(
        "--target-residue-map",
        help="Optional JSON object mapping antigen positions to target-PDB residues.",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--profile", choices=sorted(PROFILE_SETTINGS), default="smoke")
    parser.add_argument("--target-chain", default="A")
    parser.add_argument("--binder-chain", default="B")
    parser.add_argument("--germinal-repo-dir")
    parser.add_argument("--af-params-dir")
    parser.add_argument("--stage-template-pdbs", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)

    template_pdbs = _parse_key_value(
        args.template_scfv_pdb, "--template-scfv-pdb"
    )
    residue_map: Mapping[Any, Any] | None = None
    if args.target_residue_map:
        residue_map_path = Path(args.target_residue_map).expanduser()
        try:
            loaded_map = json.loads(residue_map_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise GerminalAdapterError(
                f"Could not load target residue map {residue_map_path}: {exc}"
            ) from exc
        if not isinstance(loaded_map, Mapping):
            raise GerminalAdapterError("Target residue-map JSON root must be an object.")
        residue_map = loaded_map

    handoff = build_germinal_jobs(
        args.request,
        target_pdb_path=args.target_pdb,
        template_scfv_pdbs=template_pdbs,
        target_residue_map=residue_map,
        handoff_root=args.output_dir,
        profile=args.profile,
        target_chain=args.target_chain,
        binder_chain=args.binder_chain,
        germinal_repo_dir=args.germinal_repo_dir,
        af_params_dir=args.af_params_dir,
    )
    result = write_germinal_handoff(
        handoff,
        args.output_dir,
        stage_template_pdbs=args.stage_template_pdbs,
        overwrite=args.overwrite,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
