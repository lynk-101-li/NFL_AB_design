#!/usr/bin/env python3
"""Assemble a hash-bound, deliberately blocked target-manifest review candidate.

The command joins the independently prepared antigen, SASA-review, and antibody
template bundles.  It is intentionally *not* a promotion command: reviewer
identity fields remain empty, authorization remains false, and the output state
remains ``blocked_pending_human_review``.  A human must create the formal
``config/target_structure_manifest.json`` separately after completing review.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CANDIDATE = Path(
    "input/structures/target_structure_manifest.candidate.blocked.json"
)
DEFAULT_TEMPLATE_FRAGMENT = Path(
    "input/template_structures/antibody_template_manifest_fragment.blocked.json"
)
DEFAULT_SASA_REVIEW = Path(
    "input/structures/NEFL_P07196_AFDB_v6_280-377_sasa_review.json"
)
DEFAULT_OUTPUT = Path(
    "input/structures/target_structure_manifest.review_candidate.blocked.json"
)
FORMAL_MANIFEST = (REPO_ROOT / "config/target_structure_manifest.json").resolve()

TARGET_MANIFEST_SCHEMA = "nfl_ab_design.target_structure_manifest.v1"
TARGET_EVIDENCE_SCHEMA = "nfl_ab_design.target_structure_evidence.v1"
TEMPLATE_FRAGMENT_SCHEMA = (
    "nfl_ab_design.antibody_template_manifest_fragment.v1"
)
TEMPLATE_EVIDENCE_SCHEMA = "nfl_ab_design.antibody_template_evidence.v1"
SASA_REVIEW_SCHEMA = "nfl_ab_design.target_structure_sasa_review.v1"
BLOCKED_STATE = "blocked_pending_human_review"
SASA_BLOCKED_STATE = "computational_proposal_pending_human_review"

TEMPLATE_IDS = (
    "template_7-H11-D3-2-C7",
    "template_15-C12-H6",
)
EXPECTED_WINDOWS = {
    "helix_surface_323_331": (323, 331),
    "C_boundary_368_377": (368, 377),
}
EXPECTED_HOTSPOTS = {
    "helix_surface_323_331": [325, 329],
    "C_boundary_368_377": [368, 372, 375],
}
EXPECTED_TARGET_POSITIONS = tuple(range(280, 378))
SHA256_RE = re.compile(r"[0-9a-f]{64}")


class ReviewCandidateAssemblyError(ValueError):
    """Raised when a review candidate cannot be assembled safely."""


def _resolve(path: str | Path) -> Path:
    value = Path(path).expanduser()
    if not value.is_absolute():
        value = REPO_ROOT / value
    return value.resolve()


def _portable(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ReviewCandidateAssemblyError(f"{label} must be a JSON object")
    return value


def _list(value: Any, *, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ReviewCandidateAssemblyError(f"{label} must be a JSON list")
    return value


def _string(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ReviewCandidateAssemblyError(f"{label} must be a non-empty string")
    return value.strip()


def _load_json(path: Path, *, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise ReviewCandidateAssemblyError(f"{label} does not exist: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReviewCandidateAssemblyError(
            f"Cannot read valid JSON from {label} {path}: {exc}"
        ) from exc
    if not isinstance(value, dict):
        raise ReviewCandidateAssemblyError(f"{label} root must be a JSON object")
    return value


def _file_sha256(path: Path) -> str:
    digest = sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise ReviewCandidateAssemblyError(f"Cannot hash file {path}: {exc}") from exc
    return digest.hexdigest()


def _expected_sha256(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise ReviewCandidateAssemblyError(
            f"{label} must be a lowercase 64-character SHA-256 digest"
        )
    return value


def _require_file_hash(
    raw_path: Any,
    raw_hash: Any,
    *,
    label: str,
    require_repository_file: bool = True,
) -> tuple[Path, str]:
    if isinstance(raw_path, Path):
        path = _resolve(raw_path)
    else:
        path = _resolve(_string(raw_path, label=f"{label}.path"))
    if require_repository_file:
        try:
            path.relative_to(REPO_ROOT)
        except ValueError as exc:
            raise ReviewCandidateAssemblyError(
                f"{label}.path must resolve inside the repository: {path}"
            ) from exc
    if not path.is_file():
        raise ReviewCandidateAssemblyError(f"{label}.path does not exist: {path}")
    expected = _expected_sha256(raw_hash, label=f"{label}.sha256")
    observed = _file_sha256(path)
    if observed != expected:
        raise ReviewCandidateAssemblyError(
            f"{label} hash mismatch for {path}: expected {expected}, observed {observed}"
        )
    return path, observed


def _require_pristine_review(value: Any, *, label: str, status: bool = False) -> None:
    review = _mapping(value, label=label)
    if review.get("reviewed_by") != "" or review.get("reviewed_at") != "":
        raise ReviewCandidateAssemblyError(
            f"{label} must retain empty reviewer identity and timestamp"
        )
    if review.get("contracts_acknowledged") is not False:
        raise ReviewCandidateAssemblyError(
            f"{label}.contracts_acknowledged must remain JSON false"
        )
    if status and review.get("status") != BLOCKED_STATE:
        raise ReviewCandidateAssemblyError(
            f"{label}.status must remain {BLOCKED_STATE!r}"
        )


def _require_schema_state(
    value: Mapping[str, Any],
    *,
    schema: str,
    label: str,
) -> None:
    if value.get("schema") != schema:
        raise ReviewCandidateAssemblyError(
            f"{label}.schema must be exactly {schema!r}, got {value.get('schema')!r}"
        )
    if value.get("execution_state") != BLOCKED_STATE:
        raise ReviewCandidateAssemblyError(
            f"{label}.execution_state must remain {BLOCKED_STATE!r}"
        )


def _validate_path_hash_objects(value: Any, *, label: str) -> None:
    """Validate every conventional ``{path, sha256}`` evidence object recursively.

    Persisted repository references are checked against their files.  Absolute
    scratch-space paths describe upstream provenance, not handoff inputs; their
    digests are syntax-checked but are not required to survive on another host.
    The four final template coordinate files are validated separately below.
    """

    if isinstance(value, Mapping):
        if "path" in value and "sha256" in value:
            path = _resolve(_string(value["path"], label=f"{label}.path"))
            expected = _expected_sha256(value["sha256"], label=f"{label}.sha256")
            try:
                path.relative_to(REPO_ROOT)
            except ValueError:
                pass
            else:
                _require_file_hash(path, expected, label=label)
        for key, child in value.items():
            _validate_path_hash_objects(child, label=f"{label}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _validate_path_hash_objects(child, label=f"{label}[{index}]")


def _validate_hotspot_mapping(value: Any, *, label: str) -> None:
    mapping = _mapping(value, label=label)
    if set(mapping) != set(EXPECTED_HOTSPOTS):
        raise ReviewCandidateAssemblyError(
            f"{label} must contain exactly {sorted(EXPECTED_HOTSPOTS)}"
        )
    for epitope_id, expected in EXPECTED_HOTSPOTS.items():
        observed = mapping.get(epitope_id)
        if observed != expected:
            raise ReviewCandidateAssemblyError(
                f"{label}[{epitope_id!r}] must be exactly {expected}, got {observed!r}"
            )
        start, end = EXPECTED_WINDOWS[epitope_id]
        if any(
            isinstance(position, bool)
            or not isinstance(position, int)
            or not start <= position <= end
            for position in observed
        ):
            raise ReviewCandidateAssemblyError(
                f"{label}[{epitope_id!r}] falls outside request window {start}-{end}"
            )


def _validate_candidate(
    candidate: Mapping[str, Any], candidate_path: Path
) -> dict[str, Any]:
    _require_schema_state(
        candidate, schema=TARGET_MANIFEST_SCHEMA, label="candidate manifest"
    )
    _require_pristine_review(candidate.get("review"), label="candidate manifest.review")
    selected = _mapping(
        candidate.get("selected_hotspots_by_epitope"),
        label="candidate manifest.selected_hotspots_by_epitope",
    )
    if set(selected) != set(EXPECTED_HOTSPOTS) or any(selected.values()):
        raise ReviewCandidateAssemblyError(
            "candidate manifest selected_hotspots_by_epitope must remain the pristine "
            "empty blocked mapping"
        )
    _validate_hotspot_mapping(
        candidate.get("proposed_hotspots_by_epitope_pending_human_review"),
        label="candidate manifest proposed hotspots",
    )
    if candidate.get("oligomeric_state") != "single_chain_monomer":
        raise ReviewCandidateAssemblyError(
            "candidate manifest must remain the reviewed-scope single_chain_monomer"
        )
    target_chain = _string(
        candidate.get("target_chain"), label="candidate manifest.target_chain"
    )
    if target_chain != "A":
        raise ReviewCandidateAssemblyError("candidate target chain must be exactly A")
    full_to_pdb = _mapping(
        candidate.get("full_coordinate_to_pdb"),
        label="candidate manifest.full_coordinate_to_pdb",
    )
    full_to_local = _mapping(
        candidate.get("full_coordinate_to_local_1_based"),
        label="candidate manifest.full_coordinate_to_local_1_based",
    )
    expected_position_keys = {str(position) for position in EXPECTED_TARGET_POSITIONS}
    if set(full_to_pdb) != expected_position_keys or set(full_to_local) != expected_position_keys:
        raise ReviewCandidateAssemblyError(
            "candidate coordinate mappings must contain exactly full positions 280-377"
        )
    for position in EXPECTED_TARGET_POSITIONS:
        address = _mapping(
            full_to_pdb.get(str(position)),
            label=f"candidate full_coordinate_to_pdb[{position}]",
        )
        if address != {
            "chain_id": "A",
            "residue_number": position,
            "insertion_code": "",
        }:
            raise ReviewCandidateAssemblyError(
                f"candidate PDB-coordinate mapping is invalid at residue {position}"
            )
        if full_to_local.get(str(position)) != position - 279:
            raise ReviewCandidateAssemblyError(
                f"candidate local-coordinate mapping is invalid at residue {position}"
            )

    provenance = _mapping(
        candidate.get("candidate_provenance"),
        label="candidate manifest.candidate_provenance",
    )
    target_path, target_hash = _require_file_hash(
        candidate.get("target_pdb_path"),
        provenance.get("target_pdb_sha256"),
        label="candidate target PDB",
    )
    hotspot_reference = _mapping(
        candidate.get("hotspot_evidence"), label="candidate manifest.hotspot_evidence"
    )
    evidence_path, evidence_hash = _require_file_hash(
        hotspot_reference.get("path"),
        hotspot_reference.get("sha256"),
        label="target structure evidence",
    )
    evidence = _load_json(evidence_path, label="target structure evidence")
    _require_schema_state(
        evidence, schema=TARGET_EVIDENCE_SCHEMA, label="target structure evidence"
    )
    _validate_path_hash_objects(evidence, label="target structure evidence")

    source = _mapping(evidence.get("source"), label="target structure evidence.source")
    extraction = _mapping(
        evidence.get("extraction"), label="target structure evidence.extraction"
    )
    source_path, source_hash = _require_file_hash(
        source.get("source_pdb_path"),
        source.get("sha256"),
        label="target structure full-length source",
    )
    evidence_target_path, evidence_target_hash = _require_file_hash(
        extraction.get("target_pdb_path"),
        extraction.get("target_pdb_sha256"),
        label="target structure cropped target",
    )
    _require_file_hash(
        extraction.get("target_fasta_path"),
        extraction.get("target_fasta_sha256"),
        label="target structure cropped FASTA",
    )
    validation = _mapping(
        evidence.get("validation"), label="target structure evidence.validation"
    )
    required_true_validation = (
        "pinned_source_sha256_match",
        "source_seqres_equals_atom_sequence",
        "source_coordinates_contiguous_1_through_543",
        "target_coordinates_contiguous_280_through_377",
        "target_sequence_equals_campaign_contract",
        "target_backbone_N_CA_C_O_complete_for_all_residues",
        "source_canonical_heavy_atom_topology_complete_for_all_residues",
        "target_canonical_heavy_atom_topology_complete_for_all_residues",
        "epitope_sequences_match_campaign_contract",
    )
    failed_validation = sorted(
        key for key in required_true_validation if validation.get(key) is not True
    )
    if failed_validation:
        raise ReviewCandidateAssemblyError(
            "target structure evidence contains failed validation flags: "
            + ", ".join(failed_validation)
        )
    if validation.get("source_seqres_length") != 543:
        raise ReviewCandidateAssemblyError(
            "target structure evidence source_seqres_length must be exactly 543"
        )
    if validation.get("source_atom_residue_count") != 543:
        raise ReviewCandidateAssemblyError(
            "target structure evidence source_atom_residue_count must be exactly 543"
        )
    if extraction.get("full_coordinate_range_1_based_inclusive") != [280, 377]:
        raise ReviewCandidateAssemblyError(
            "target structure extraction full-coordinate range must be exactly 280-377"
        )
    if extraction.get("local_coordinate_range_1_based_inclusive") != [1, 98]:
        raise ReviewCandidateAssemblyError(
            "target structure extraction local-coordinate range must be exactly 1-98"
        )
    if extraction.get("residue_count") != 98:
        raise ReviewCandidateAssemblyError(
            "target structure extraction residue_count must be exactly 98"
        )
    if evidence_target_path != target_path or evidence_target_hash != target_hash:
        raise ReviewCandidateAssemblyError(
            "candidate target PDB is not hash-bound to target structure evidence"
        )
    if provenance.get("source_sha256") != source_hash:
        raise ReviewCandidateAssemblyError(
            "candidate source hash is not bound to target structure evidence"
        )
    if candidate.get("antigen_sequence_in_pdb_order") != extraction.get(
        "antigen_sequence_in_pdb_order"
    ):
        raise ReviewCandidateAssemblyError(
            "candidate antigen sequence is not bound to target structure evidence"
        )
    proposals = _list(
        evidence.get("hotspot_proposals"),
        label="target structure evidence.hotspot_proposals",
    )
    if len(proposals) != len(EXPECTED_HOTSPOTS):
        raise ReviewCandidateAssemblyError(
            "target structure evidence must contain exactly two epitope proposals"
        )
    proposal_mapping: dict[str, Any] = {}
    for row in proposals:
        row_map = _mapping(row, label="target structure evidence hotspot proposal")
        epitope_id = _string(
            row_map.get("epitope_id"), label="target structure evidence epitope_id"
        )
        if epitope_id in proposal_mapping:
            raise ReviewCandidateAssemblyError(
                f"Duplicate target structure epitope proposal {epitope_id!r}"
            )
        proposal_mapping[epitope_id] = row_map.get(
            "proposed_hotspots_pending_human_review"
        )
        if row_map.get("full_coordinate_range_1_based_inclusive") != list(
            EXPECTED_WINDOWS.get(epitope_id, ())
        ):
            raise ReviewCandidateAssemblyError(
                f"Target structure request window mismatch for {epitope_id!r}"
            )
    _validate_hotspot_mapping(
        proposal_mapping, label="target structure evidence hotspot proposals"
    )
    review_gate = _mapping(
        evidence.get("review_gate"), label="target structure evidence.review_gate"
    )
    if review_gate.get("status") != "human_review_required":
        raise ReviewCandidateAssemblyError(
            "target structure evidence review gate must require human review"
        )
    return {
        "candidate_path": candidate_path,
        "candidate_hash": _file_sha256(candidate_path),
        "target_path": target_path,
        "target_hash": target_hash,
        "source_path": source_path,
        "source_hash": source_hash,
        "evidence_path": evidence_path,
        "evidence_hash": evidence_hash,
    }


def _validate_sasa(
    sasa: Mapping[str, Any], sasa_path: Path, candidate_info: Mapping[str, Any]
) -> dict[str, Any]:
    if sasa.get("schema") != SASA_REVIEW_SCHEMA:
        raise ReviewCandidateAssemblyError(
            f"SASA review schema must be exactly {SASA_REVIEW_SCHEMA!r}"
        )
    if sasa.get("analysis_state") != SASA_BLOCKED_STATE:
        raise ReviewCandidateAssemblyError(
            f"SASA review analysis_state must remain {SASA_BLOCKED_STATE!r}"
        )
    conclusion = _mapping(
        sasa.get("review_conclusion"), label="SASA review.review_conclusion"
    )
    if conclusion.get("human_review_still_required") is not True:
        raise ReviewCandidateAssemblyError(
            "SASA review must retain human_review_still_required=true"
        )
    if conclusion.get("current_proposals_supported") is not True:
        raise ReviewCandidateAssemblyError(
            "SASA review must explicitly support the staged recommendations"
        )
    inputs = _mapping(sasa.get("inputs"), label="SASA review.inputs")
    full = _mapping(inputs.get("full_length_context"), label="SASA full-length input")
    crop = _mapping(inputs.get("cropped_design_target"), label="SASA cropped input")
    full_path, full_hash = _require_file_hash(
        full.get("path"), full.get("sha256"), label="SASA full-length input"
    )
    crop_path, crop_hash = _require_file_hash(
        crop.get("path"), crop.get("sha256"), label="SASA cropped input"
    )
    if full.get("chain") != "A" or full.get(
        "residue_range_1_based_inclusive"
    ) != [1, 543]:
        raise ReviewCandidateAssemblyError(
            "SASA full-length input must be chain A with residue range 1-543"
        )
    if crop.get("chain") != "A" or crop.get(
        "residue_range_1_based_inclusive"
    ) != [280, 377]:
        raise ReviewCandidateAssemblyError(
            "SASA cropped input must be chain A with residue range 280-377"
        )
    if (
        full_path != candidate_info["source_path"]
        or full_hash != candidate_info["source_hash"]
    ):
        raise ReviewCandidateAssemblyError(
            "SASA full-length input is not hash-bound to target structure evidence"
        )
    if (
        crop_path != candidate_info["target_path"]
        or crop_hash != candidate_info["target_hash"]
    ):
        raise ReviewCandidateAssemblyError(
            "SASA cropped input is not hash-bound to the candidate target PDB"
        )
    reviews = _list(sasa.get("epitope_reviews"), label="SASA epitope_reviews")
    if len(reviews) != len(EXPECTED_HOTSPOTS):
        raise ReviewCandidateAssemblyError(
            "SASA review must contain exactly two epitope reviews"
        )
    recommendations: dict[str, Any] = {}
    for row in reviews:
        row_map = _mapping(row, label="SASA epitope review")
        epitope_id = _string(row_map.get("epitope_id"), label="SASA epitope_id")
        if epitope_id in recommendations:
            raise ReviewCandidateAssemblyError(
                f"Duplicate SASA epitope review {epitope_id!r}"
            )
        if row_map.get("full_coordinate_range_1_based_inclusive") != list(
            EXPECTED_WINDOWS.get(epitope_id, ())
        ):
            raise ReviewCandidateAssemblyError(
                f"SASA request window mismatch for {epitope_id!r}"
            )
        recommendations[epitope_id] = row_map.get(
            "recommended_hotspots_pending_human_review"
        )
    _validate_hotspot_mapping(recommendations, label="SASA recommended hotspots")
    return {
        "path": sasa_path,
        "hash": _file_sha256(sasa_path),
        "full_path": full_path,
        "full_hash": full_hash,
        "crop_path": crop_path,
        "crop_hash": crop_hash,
    }


def _validate_templates(
    fragment: Mapping[str, Any], fragment_path: Path
) -> dict[str, Any]:
    _require_schema_state(
        fragment, schema=TEMPLATE_FRAGMENT_SCHEMA, label="template fragment"
    )
    _require_pristine_review(fragment.get("review"), label="template fragment.review")
    if fragment.get("real_model_handoff_authorized") is not False:
        raise ReviewCandidateAssemblyError(
            "template fragment real_model_handoff_authorized must remain false"
        )
    evidence_path, evidence_hash = _require_file_hash(
        fragment.get("evidence_path"),
        fragment.get("evidence_sha256"),
        label="antibody template evidence",
    )
    evidence = _load_json(evidence_path, label="antibody template evidence")
    _require_schema_state(
        evidence, schema=TEMPLATE_EVIDENCE_SCHEMA, label="antibody template evidence"
    )
    _require_pristine_review(
        evidence.get("review"), label="antibody template evidence.review", status=True
    )
    if evidence.get("real_model_handoff_authorized") is not False:
        raise ReviewCandidateAssemblyError(
            "antibody template evidence authorization must remain false"
        )
    if evidence.get("known_positive_cdr_sequences_used") is not False:
        raise ReviewCandidateAssemblyError(
            "antibody template evidence must confirm known-positive CDRs were not used"
        )
    _validate_path_hash_objects(evidence, label="antibody template evidence")

    framework = _mapping(
        fragment.get("framework_coordinate_inputs"),
        label="template fragment.framework_coordinate_inputs",
    )
    evidence_templates = _mapping(
        evidence.get("templates"), label="antibody template evidence.templates"
    )
    expected_ids = set(TEMPLATE_IDS)
    if set(evidence_templates) != expected_ids:
        raise ReviewCandidateAssemblyError(
            "antibody template evidence must contain exactly the two campaign templates"
        )
    validated_files: dict[str, dict[str, dict[str, str]]] = {}
    all_paths: list[Path] = []
    for engine_key, evidence_path_key, evidence_hash_key in (
        ("rfantibody_hlt_pdbs", "hlt_path", "hlt_sha256"),
        ("germinal_scfv_pdbs", "scfv_path", "scfv_sha256"),
    ):
        paths = _mapping(framework.get(engine_key), label=f"framework.{engine_key}")
        if set(paths) != expected_ids:
            raise ReviewCandidateAssemblyError(
                f"framework.{engine_key} must contain exactly {sorted(expected_ids)}"
            )
        validated_files[engine_key] = {}
        engine_paths: list[Path] = []
        for template_id in TEMPLATE_IDS:
            record = _mapping(
                evidence_templates[template_id],
                label=f"antibody template evidence.templates.{template_id}",
            )
            path, digest = _require_file_hash(
                paths[template_id],
                record.get(evidence_hash_key),
                label=f"{engine_key}.{template_id}",
            )
            evidence_record_path = _resolve(
                _string(
                    record.get(evidence_path_key),
                    label=f"templates.{template_id}.{evidence_path_key}",
                )
            )
            if path != evidence_record_path:
                raise ReviewCandidateAssemblyError(
                    f"{engine_key}.{template_id} path does not match template evidence"
                )
            engine_paths.append(path)
            all_paths.append(path)
            validated_files[engine_key][template_id] = {
                "path": _portable(path),
                "sha256": digest,
            }
        if len(set(engine_paths)) != len(TEMPLATE_IDS):
            raise ReviewCandidateAssemblyError(
                f"framework.{engine_key} reuses one file for both templates"
            )
    if len(set(all_paths)) != 2 * len(TEMPLATE_IDS):
        raise ReviewCandidateAssemblyError(
            "RFantibody and Germinal inputs must be four distinct coordinate files"
        )
    return {
        "fragment_path": fragment_path,
        "fragment_hash": _file_sha256(fragment_path),
        "evidence_path": evidence_path,
        "evidence_hash": evidence_hash,
        "framework": deepcopy(dict(framework)),
        "validated_files": validated_files,
    }


def _write_atomic(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def assemble_review_candidate_manifest(
    *,
    candidate_manifest: str | Path = DEFAULT_CANDIDATE,
    template_fragment: str | Path = DEFAULT_TEMPLATE_FRAGMENT,
    sasa_review: str | Path = DEFAULT_SASA_REVIEW,
    output: str | Path = DEFAULT_OUTPUT,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Validate and combine the three blocked bundles into one blocked candidate."""

    candidate_path = _resolve(candidate_manifest)
    fragment_path = _resolve(template_fragment)
    sasa_path = _resolve(sasa_review)
    output_path = _resolve(output)
    if output_path == FORMAL_MANIFEST:
        raise ReviewCandidateAssemblyError(
            "This assembler can never create or overwrite "
            "config/target_structure_manifest.json"
        )
    if output_path.exists() and not overwrite:
        raise ReviewCandidateAssemblyError(
            f"Refusing to overwrite existing output {output_path}; use --overwrite "
            "only for the blocked review candidate"
        )

    candidate = _load_json(candidate_path, label="candidate manifest")
    fragment = _load_json(fragment_path, label="template fragment")
    sasa = _load_json(sasa_path, label="SASA review")
    candidate_info = _validate_candidate(candidate, candidate_path)
    template_info = _validate_templates(fragment, fragment_path)
    sasa_info = _validate_sasa(sasa, sasa_path, candidate_info)

    assembled = deepcopy(candidate)
    assembled["execution_state"] = BLOCKED_STATE
    assembled["review"] = {
        "reviewed_by": "",
        "reviewed_at": "",
        "contracts_acknowledged": False,
    }
    assembled["real_model_handoff_authorized"] = False
    assembled["selected_hotspots_by_epitope"] = deepcopy(EXPECTED_HOTSPOTS)
    assembled["hotspot_review_state"] = (
        "computational_recommendations_staged_pending_human_review"
    )
    assembled["framework_coordinate_inputs"] = template_info["framework"]
    assembled["blocking_reasons"] = [
        "Human review of the theoretical target structure and SASA evidence "
        "has not been recorded.",
        "Human review of both neutral antibody coordinate templates has not been recorded.",
        "Reviewer identity, timezone-qualified review time, and contract "
        "acknowledgement remain empty/false.",
        "This review candidate is not the formal config/target_structure_manifest.json.",
    ]
    assembled["review_evidence"] = {
        "candidate_manifest": {
            "path": _portable(candidate_info["candidate_path"]),
            "sha256": candidate_info["candidate_hash"],
        },
        "target_structure_evidence": {
            "path": _portable(candidate_info["evidence_path"]),
            "sha256": candidate_info["evidence_hash"],
        },
        "sasa_review": {
            "path": _portable(sasa_info["path"]),
            "sha256": sasa_info["hash"],
        },
        "antibody_template_manifest_fragment": {
            "path": _portable(template_info["fragment_path"]),
            "sha256": template_info["fragment_hash"],
        },
        "antibody_template_evidence": {
            "path": _portable(template_info["evidence_path"]),
            "sha256": template_info["evidence_hash"],
        },
        "bound_structure_inputs": {
            "full_length_context": {
                "path": _portable(sasa_info["full_path"]),
                "sha256": sasa_info["full_hash"],
            },
            "cropped_design_target": {
                "path": _portable(sasa_info["crop_path"]),
                "sha256": sasa_info["crop_hash"],
            },
        },
        "bound_template_coordinate_inputs": template_info["validated_files"],
    }
    assembled["promotion_instructions"] = [
        "Inspect the AlphaFold DB monomer, full-length-context SASA evidence, "
        "crop-boundary limitation, monomer surface geometry, and both neutral template geometries.",
        "Confirm or revise each staged hotspot list with a written rationale; "
        "if any input or recommendation changes, regenerate this candidate so "
        "every SHA-256 binding is current.",
        "After review, create config/target_structure_manifest.json separately; "
        "do not rename this blocked candidate.",
        "Only in that separately reviewed file set execution_state to "
        "reviewed_ready_for_handoff, fill reviewed_by and a timezone-qualified "
        "reviewed_at, set contracts_acknowledged true, and explicitly authorize "
        "the real-model handoff.",
    ]

    payload = (
        json.dumps(assembled, indent=2, ensure_ascii=False, sort_keys=False) + "\n"
    ).encode("utf-8")
    _write_atomic(output_path, payload)
    return {
        "schema": TARGET_MANIFEST_SCHEMA,
        "execution_state": BLOCKED_STATE,
        "output_path": _portable(output_path),
        "output_sha256": sha256(payload).hexdigest(),
        "selected_hotspots_by_epitope": deepcopy(EXPECTED_HOTSPOTS),
        "framework_coordinate_inputs": deepcopy(template_info["framework"]),
        "real_model_handoff_authorized": False,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-manifest", default=str(DEFAULT_CANDIDATE))
    parser.add_argument("--template-fragment", default=str(DEFAULT_TEMPLATE_FRAGMENT))
    parser.add_argument("--sasa-review", default=str(DEFAULT_SASA_REVIEW))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace only the blocked review-candidate output after full revalidation.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        summary = assemble_review_candidate_manifest(
            candidate_manifest=args.candidate_manifest,
            template_fragment=args.template_fragment,
            sasa_review=args.sasa_review,
            output=args.output,
            overwrite=args.overwrite,
        )
    except ReviewCandidateAssemblyError as exc:
        print(f"ERROR: {exc}")
        return 2
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
