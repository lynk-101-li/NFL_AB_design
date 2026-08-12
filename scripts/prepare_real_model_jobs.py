#!/usr/bin/env python3
"""Compile reviewed NfL structure inputs into non-executing model handoffs.

This command is intentionally a compiler, not a model runner.  It reads the
three normalized requests exported by :mod:`nfl_ab_design.workflow`, injects
manually curated epitope hotspots into the RFantibody and Germinal requests,
validates every coordinate-bearing input through the model-specific adapters,
and writes reviewable manifests and input files.  It never invokes a shell,
subprocess, model executable, package installer, or network client.

All relative paths are interpreted from the repository root, independent of
the caller's current working directory.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
import tempfile
from typing import Any, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from nfl_ab_design import __version__ as PACKAGE_VERSION  # noqa: E402
from nfl_ab_design.adapters.germinal import (  # noqa: E402
    GerminalAdapterError,
    build_germinal_jobs,
    write_germinal_handoff,
)
from nfl_ab_design.adapters.iggm import (  # noqa: E402
    IgGMAdapterError,
    build_iggm_jobs,
)
from nfl_ab_design.adapters.rfantibody import (  # noqa: E402
    RFANTIBODY_OFFICIAL_MAIN_SHA,
    RFantibodyAdapterError,
    build_rfantibody_plan,
)


TARGET_MANIFEST_SCHEMA = "nfl_ab_design.target_structure_manifest.v1"
UNIFIED_HANDOFF_SCHEMA = "nfl_ab_design.real_model_handoff.v1"
DESIGN_REQUEST_INDEX_SCHEMA = "nfl_ab_design.design_request_index.v1"
REQUEST_FILENAMES = {
    "RFantibody": "rfantibody_design_request.json",
    "IgGM": "iggm_design_request.json",
    "Germinal": "germinal_design_request.json",
}
NORMALIZED_REQUEST_SCHEMA = "nfl_ab_design.normalized_de_novo_request.v1"
PLACEHOLDER_PATTERN = re.compile(r"<[^<>]+>")
PLACEHOLDER_PREFIXES = ("REPLACE_WITH", "PLACEHOLDER", "TODO")
AA20 = frozenset("ACDEFGHIKLMNPQRSTVWY")
COMPILER_SOURCE_PATH = Path(__file__).resolve()
ADAPTER_SOURCE_PATHS = {
    "RFantibody": SRC_ROOT / "nfl_ab_design" / "adapters" / "rfantibody.py",
    "IgGM": SRC_ROOT / "nfl_ab_design" / "adapters" / "iggm.py",
    "Germinal": SRC_ROOT / "nfl_ab_design" / "adapters" / "germinal.py",
}
SELF_TEST_AA1_TO_AA3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


class RealModelHandoffError(ValueError):
    """Raised when the unified handoff cannot be compiled safely."""


def _repo_path(value: str | Path, *, label: str) -> Path:
    """Resolve a filesystem path relative to :data:`REPO_ROOT`."""

    text = str(value).strip()
    if not text:
        raise RealModelHandoffError(f"{label} must not be empty")
    _reject_placeholder(text, label)
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def _safe_path_component(value: Any, *, label: str) -> str:
    """Validate a provenance identifier before using it in a default path."""

    text = _require_nonempty_string(value, label=label)
    if text in {".", ".."} or not re.fullmatch(r"[A-Za-z0-9._-]+", text):
        raise RealModelHandoffError(
            f"{label} must contain only ASCII letters, digits, '.', '_', or '-' "
            f"when used in the automatic handoff path, got {text!r}. Supply an "
            "explicit --output-dir if the reviewed identifier requires other characters."
        )
    return text


def _reject_placeholder(value: str, label: str) -> None:
    stripped = value.strip()
    if PLACEHOLDER_PATTERN.search(stripped) or stripped.upper().startswith(
        PLACEHOLDER_PREFIXES
    ):
        raise RealModelHandoffError(
            f"{label} still contains a placeholder: {value!r}. "
            "Replace it with a reviewed value before compiling real-model jobs."
        )


def _load_json_object(path: Path, *, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise RealModelHandoffError(f"{label} does not exist or is not a file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise RealModelHandoffError(f"Cannot read {label}: {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise RealModelHandoffError(f"Invalid JSON in {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RealModelHandoffError(f"{label} root must be a JSON object: {path}")
    return value


def _require_mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise RealModelHandoffError(f"{label} must be a JSON object")
    return value


def _require_nonempty_string(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise RealModelHandoffError(f"{label} must be a non-empty string")
    text = value.strip()
    _reject_placeholder(text, label)
    return text


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _compiler_provenance() -> dict[str, Any]:
    missing = [str(path) for path in ADAPTER_SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise RealModelHandoffError(
            "Cannot hash adapter source provenance because files are missing: "
            + ", ".join(missing)
        )
    return {
        "package_version": PACKAGE_VERSION,
        "compiler_source": {
            "path": str(COMPILER_SOURCE_PATH),
            "sha256": _file_sha256(COMPILER_SOURCE_PATH),
        },
        "adapter_source_by_engine": {
            engine: {"path": str(path), "sha256": _file_sha256(path)}
            for engine, path in ADAPTER_SOURCE_PATHS.items()
        },
    }


def _text_sha256(value: str) -> str:
    return sha256(value.encode("utf-8")).hexdigest()


def _resolve_path_mapping(value: Any, *, label: str) -> dict[str, Path]:
    mapping = _require_mapping(value, label=label)
    if not mapping:
        raise RealModelHandoffError(f"{label} must not be empty")
    resolved: dict[str, Path] = {}
    for raw_key, raw_value in mapping.items():
        key = _require_nonempty_string(str(raw_key), label=f"{label} key")
        if not isinstance(raw_value, (str, Path)):
            raise RealModelHandoffError(f"{label}[{key!r}] must be a path string")
        path = _repo_path(raw_value, label=f"{label}[{key!r}]")
        if not path.is_file():
            raise RealModelHandoffError(
                f"{label}[{key!r}] does not exist or is not a file: {path}"
            )
        resolved[key] = path
    return resolved


def _normalize_hotspot_mapping(
    value: Any,
    *,
    epitopes: Sequence[Mapping[str, Any]],
) -> dict[str, list[int]]:
    mapping = _require_mapping(value, label="selected_hotspots_by_epitope")
    expected_ids = {
        _require_nonempty_string(
            epitope.get("epitope_id"), label="request epitope_id"
        )
        for epitope in epitopes
    }
    actual_ids = {str(key) for key in mapping}
    missing = sorted(expected_ids - actual_ids)
    extra = sorted(actual_ids - expected_ids)
    if missing or extra:
        raise RealModelHandoffError(
            "selected_hotspots_by_epitope must match the normalized request exactly; "
            f"missing={missing}, extra={extra}"
        )

    normalized: dict[str, list[int]] = {}
    for epitope in epitopes:
        epitope_id = str(epitope["epitope_id"])
        raw_positions = mapping[epitope_id]
        if (
            not isinstance(raw_positions, Sequence)
            or isinstance(raw_positions, (str, bytes, bytearray))
            or not raw_positions
        ):
            raise RealModelHandoffError(
                f"selected_hotspots_by_epitope[{epitope_id!r}] must be a non-empty "
                "list curated after structure review"
            )
        positions: list[int] = []
        for raw_position in raw_positions:
            if isinstance(raw_position, bool):
                raise RealModelHandoffError(
                    f"Hotspot for epitope {epitope_id} is boolean, not an integer: "
                    f"{raw_position!r}"
                )
            try:
                position = int(raw_position)
            except (TypeError, ValueError) as exc:
                raise RealModelHandoffError(
                    f"Hotspot for epitope {epitope_id} is not an integer: "
                    f"{raw_position!r}"
                ) from exc
            positions.append(position)
        if len(set(positions)) != len(positions):
            raise RealModelHandoffError(
                f"Epitope {epitope_id} contains duplicate selected hotspots"
            )
        positions = sorted(positions)
        start = int(epitope["start_1_based"])
        end = int(epitope["end_1_based_inclusive"])
        outside = [position for position in positions if not start <= position <= end]
        if outside:
            raise RealModelHandoffError(
                f"Selected hotspots for {epitope_id} fall outside {start}-{end}: "
                f"{outside}"
            )
        candidate_window = epitope.get("candidate_hotspot_residue_indices", [])
        try:
            candidate_positions = sorted({int(item) for item in candidate_window})
        except (TypeError, ValueError) as exc:
            raise RealModelHandoffError(
                f"Normalized candidate hotspot window is invalid for {epitope_id}"
            ) from exc
        if candidate_positions and positions == candidate_positions:
            raise RealModelHandoffError(
                f"Selected hotspots for {epitope_id} equal the entire candidate epitope "
                "window. Curate a smaller solvent-accessible/chemically justified set; "
                "the adapters intentionally reject blind whole-window conditioning."
            )
        normalized[epitope_id] = positions
    return normalized


def _inject_hotspots(
    request: Mapping[str, Any],
    hotspots: Mapping[str, Sequence[int]],
    *,
    field: str,
) -> dict[str, Any]:
    injected = deepcopy(dict(request))
    raw_epitopes = injected.get("epitopes")
    if not isinstance(raw_epitopes, list):
        raise RealModelHandoffError("Normalized request epitopes must be a list")
    for epitope in raw_epitopes:
        if not isinstance(epitope, dict):
            raise RealModelHandoffError("Normalized request epitope rows must be objects")
        epitope_id = str(epitope.get("epitope_id", ""))
        if epitope_id not in hotspots:
            raise RealModelHandoffError(
                f"No curated hotspot list was provided for epitope {epitope_id!r}"
            )
        epitope[field] = list(hotspots[epitope_id])
    return injected


def _validate_request_set(
    request_dir: Path,
) -> tuple[
    dict[str, dict[str, Any]],
    dict[str, Path],
    dict[str, Any],
    dict[str, Any],
]:
    index_path = request_dir / "design_request_index.json"
    index: dict[str, Any] | None = None
    candidate_paths: list[Path]
    if index_path.is_file():
        index = _load_json_object(index_path, label="design request index")
        if index.get("schema") != DESIGN_REQUEST_INDEX_SCHEMA:
            raise RealModelHandoffError(
                f"Unsupported design request index schema {index.get('schema')!r}; "
                f"expected {DESIGN_REQUEST_INDEX_SCHEMA!r}: {index_path}"
            )
        raw_engines = index.get("engines")
        expected_engines = list(REQUEST_FILENAMES)
        if raw_engines != expected_engines:
            raise RealModelHandoffError(
                "design request index engines must be exactly "
                f"{expected_engines}, got {raw_engines!r}: {index_path}"
            )
        raw_hashes = _require_mapping(
            index.get("request_sha256_by_engine"),
            label="design request index request_sha256_by_engine",
        )
        if set(raw_hashes) != set(expected_engines):
            raise RealModelHandoffError(
                "design request index request_sha256_by_engine must contain exactly "
                f"{expected_engines}, got {sorted(str(key) for key in raw_hashes)}"
            )
        raw_files = index.get("request_files")
        if (
            not isinstance(raw_files, Sequence)
            or isinstance(raw_files, (str, bytes, bytearray))
            or len(raw_files) != len(expected_engines)
        ):
            raise RealModelHandoffError(
                "design request index request_files must contain exactly three "
                f"paths: {index_path}"
            )
        candidate_paths = [
            _repo_path(raw_path, label=f"design request index request_files[{index}]")
            for index, raw_path in enumerate(raw_files)
        ]
    else:
        candidate_paths = sorted(request_dir.glob("*_design_request.json"))
        if not candidate_paths:
            candidate_paths = [
                request_dir / filename for filename in REQUEST_FILENAMES.values()
            ]

    requests: dict[str, dict[str, Any]] = {}
    paths: dict[str, Path] = {}
    for path in candidate_paths:
        request = _load_json_object(path, label="normalized design request")
        actual_engine = request.get("engine")
        if actual_engine not in REQUEST_FILENAMES:
            raise RealModelHandoffError(
                f"{path} declares unsupported engine {actual_engine!r}; expected one of "
                f"{list(REQUEST_FILENAMES)}"
            )
        engine = str(actual_engine)
        if engine in requests:
            raise RealModelHandoffError(
                f"More than one normalized request declares engine {engine!r}: "
                f"{paths[engine]} and {path}"
            )
        if request.get("schema") != NORMALIZED_REQUEST_SCHEMA:
            raise RealModelHandoffError(
                f"{path} has unsupported schema {request.get('schema')!r}; expected "
                f"{NORMALIZED_REQUEST_SCHEMA!r}"
            )
        requests[engine] = request
        paths[engine] = path.resolve()

    missing_engines = sorted(set(REQUEST_FILENAMES) - set(requests))
    if missing_engines:
        raise RealModelHandoffError(
            f"Request set is missing normalized engine request(s): {missing_engines}"
        )

    run_metadata_by_engine: dict[str, dict[str, Any]] = {}
    for engine, request in requests.items():
        metadata = _require_mapping(
            request.get("run_metadata"),
            label=f"{engine} normalized request run_metadata",
        )
        run_id = _require_nonempty_string(
            metadata.get("run_id"),
            label=f"{engine} normalized request run_metadata.run_id",
        )
        run_metadata_by_engine[engine] = dict(metadata)
        if run_id != str(metadata["run_id"]):
            raise RealModelHandoffError(
                f"{engine} run_metadata.run_id contains unsupported whitespace"
            )

    source_run_metadata = run_metadata_by_engine["RFantibody"]
    for engine in ("IgGM", "Germinal"):
        if run_metadata_by_engine[engine] != source_run_metadata:
            raise RealModelHandoffError(
                "Normalized requests contain different run_metadata; RFantibody "
                f"does not match {engine}. Refusing to mix request generations."
            )

    reference = requests["RFantibody"]
    for engine in ("IgGM", "Germinal"):
        for field in ("schema", "campaign_mode", "antigen", "epitopes", "templates"):
            if requests[engine].get(field) != reference.get(field):
                raise RealModelHandoffError(
                    f"Normalized requests disagree in {field!r}: RFantibody versus "
                    f"{engine}. Regenerate all requests in one workflow run."
                )

    if index is not None:
        expected_hashes = _require_mapping(
            index["request_sha256_by_engine"],
            label="design request index request_sha256_by_engine",
        )
        for engine, path in paths.items():
            recorded_hash = expected_hashes.get(engine)
            if not isinstance(recorded_hash, str) or not re.fullmatch(
                r"[0-9a-f]{64}", recorded_hash
            ):
                raise RealModelHandoffError(
                    f"Invalid SHA-256 recorded for {engine} in {index_path}: "
                    f"{recorded_hash!r}"
                )
            actual_hash = _file_sha256(path)
            if actual_hash != recorded_hash:
                raise RealModelHandoffError(
                    f"Normalized request hash mismatch for {engine}: index records "
                    f"{recorded_hash}, file is {actual_hash}: {path}"
                )
        index_run_metadata = _require_mapping(
            index.get("run_metadata"), label="design request index run_metadata"
        )
        if dict(index_run_metadata) != source_run_metadata:
            raise RealModelHandoffError(
                "design request index run_metadata does not exactly match the three "
                "normalized requests; refusing a stale or mixed request set"
            )
        index_provenance = {
            "present": True,
            "path": str(index_path.resolve()),
            "schema": DESIGN_REQUEST_INDEX_SCHEMA,
            "sha256": _file_sha256(index_path),
        }
    else:
        index_provenance = {
            "present": False,
            "path": None,
            "schema": None,
            "sha256": None,
        }
    return requests, paths, index_provenance, source_run_metadata


def _validate_target_manifest(
    manifest: Mapping[str, Any],
    *,
    request: Mapping[str, Any],
) -> dict[str, Any]:
    if manifest.get("schema") != TARGET_MANIFEST_SCHEMA:
        raise RealModelHandoffError(
            f"Unsupported target manifest schema {manifest.get('schema')!r}; "
            f"expected {TARGET_MANIFEST_SCHEMA!r}"
        )
    manifest_state = _require_nonempty_string(
        manifest.get("execution_state"), label="target manifest execution_state"
    )
    if manifest_state != "reviewed_ready_for_handoff":
        raise RealModelHandoffError(
            "target manifest execution_state must be exactly "
            "'reviewed_ready_for_handoff', got "
            f"{manifest_state!r}. Complete structure, mapping, hotspot, and review "
            "contract checks before compiling real-model jobs."
        )
    review = _require_mapping(
        manifest.get("review"), label="target manifest review"
    )
    reviewed_by = _require_nonempty_string(
        review.get("reviewed_by"), label="target manifest review.reviewed_by"
    )
    reviewed_at = _require_nonempty_string(
        review.get("reviewed_at"), label="target manifest review.reviewed_at"
    )
    try:
        reviewed_datetime = datetime.fromisoformat(reviewed_at.replace("Z", "+00:00"))
    except ValueError as exc:
        raise RealModelHandoffError(
            "target manifest review.reviewed_at must be an ISO-8601 datetime, got "
            f"{reviewed_at!r}"
        ) from exc
    if reviewed_datetime.tzinfo is None or reviewed_datetime.utcoffset() is None:
        raise RealModelHandoffError(
            "target manifest review.reviewed_at must include an explicit timezone"
        )
    if review.get("contracts_acknowledged") is not True:
        raise RealModelHandoffError(
            "target manifest review.contracts_acknowledged must be JSON true after "
            "the reviewer checks every review_contract item"
        )
    normalized_review = {
        "reviewed_by": reviewed_by,
        "reviewed_at": reviewed_at,
        "contracts_acknowledged": True,
    }
    conformation_id = _require_nonempty_string(
        manifest.get("conformation_id"), label="target manifest conformation_id"
    )
    oligomeric_state = _require_nonempty_string(
        manifest.get("oligomeric_state"),
        label="target manifest oligomeric_state",
    )
    if oligomeric_state != "single_chain_monomer":
        raise RealModelHandoffError(
            "target manifest v1 supports only "
            "oligomeric_state='single_chain_monomer', got "
            f"{oligomeric_state!r}. Do not collapse a dimer or oligomer into the "
            "single-chain Germinal track; use a separately reviewed manifest/adapter."
        )
    target_pdb = _repo_path(
        _require_nonempty_string(
            manifest.get("target_pdb_path"), label="target manifest target_pdb_path"
        ),
        label="target manifest target_pdb_path",
    )
    if not target_pdb.is_file():
        raise RealModelHandoffError(
            f"Target PDB does not exist or is not a file: {target_pdb}"
        )
    target_chain = _require_nonempty_string(
        manifest.get("target_chain"), label="target manifest target_chain"
    )
    if len(target_chain) != 1 or not target_chain.isalnum():
        raise RealModelHandoffError(
            f"target_chain must be one alphanumeric PDB chain ID, got {target_chain!r}"
        )
    antigen_sequence = _require_nonempty_string(
        manifest.get("antigen_sequence_in_pdb_order"),
        label="target manifest antigen_sequence_in_pdb_order",
    ).upper()
    invalid_residues = sorted(set(antigen_sequence) - AA20)
    if invalid_residues:
        raise RealModelHandoffError(
            "antigen_sequence_in_pdb_order contains non-canonical residues: "
            + ", ".join(invalid_residues)
        )

    full_to_pdb = _require_mapping(
        manifest.get("full_coordinate_to_pdb"),
        label="target manifest full_coordinate_to_pdb",
    )
    if not full_to_pdb:
        raise RealModelHandoffError(
            "target manifest full_coordinate_to_pdb is empty; provide explicit "
            "full-antigen to PDB residue addresses"
        )
    full_to_local = _require_mapping(
        manifest.get("full_coordinate_to_local_1_based"),
        label="target manifest full_coordinate_to_local_1_based",
    )
    if not full_to_local:
        raise RealModelHandoffError(
            "target manifest full_coordinate_to_local_1_based is empty; provide the "
            "alignment-derived IgGM coordinate mapping"
        )

    raw_epitopes = request.get("epitopes")
    if not isinstance(raw_epitopes, list) or not all(
        isinstance(row, Mapping) for row in raw_epitopes
    ):
        raise RealModelHandoffError("Normalized request epitopes are malformed")
    epitopes = [dict(row) for row in raw_epitopes]
    hotspots = _normalize_hotspot_mapping(
        manifest.get("selected_hotspots_by_epitope"), epitopes=epitopes
    )

    framework_inputs = _require_mapping(
        manifest.get("framework_coordinate_inputs"),
        label="target manifest framework_coordinate_inputs",
    )
    rfantibody_hlt = _resolve_path_mapping(
        framework_inputs.get("rfantibody_hlt_pdbs"),
        label="framework_coordinate_inputs.rfantibody_hlt_pdbs",
    )
    germinal_scfv = _resolve_path_mapping(
        framework_inputs.get("germinal_scfv_pdbs"),
        label="framework_coordinate_inputs.germinal_scfv_pdbs",
    )
    template_ids = {
        str(row.get("template_id"))
        for row in request.get("templates", [])
        if isinstance(row, Mapping)
    }
    for label, values in (
        ("rfantibody_hlt_pdbs", rfantibody_hlt),
        ("germinal_scfv_pdbs", germinal_scfv),
    ):
        missing = sorted(template_ids - set(values))
        extra = sorted(set(values) - template_ids)
        if missing or extra:
            raise RealModelHandoffError(
                f"{label} must contain exactly one distinct file per template; "
                f"missing={missing}, extra={extra}"
            )
        resolved_values = [str(values[key]) for key in sorted(values)]
        if len(set(resolved_values)) != len(resolved_values):
            raise RealModelHandoffError(
                f"{label} reuses one coordinate file for multiple templates; "
                "provide the two distinct template geometries explicitly"
            )

    return {
        "manifest_state": manifest_state,
        "review": normalized_review,
        "conformation_id": conformation_id,
        "oligomeric_state": oligomeric_state,
        "target_pdb": target_pdb,
        "target_chain": target_chain,
        "antigen_sequence": antigen_sequence,
        "full_to_pdb": dict(full_to_pdb),
        "full_to_local": dict(full_to_local),
        "hotspots": hotspots,
        "rfantibody_hlt": rfantibody_hlt,
        "germinal_scfv": germinal_scfv,
    }


def _required_runtime_directory(
    value: str | Path | None,
    *,
    label: str,
    sentinel_file: str | None = None,
) -> Path:
    if value is None:
        raise RealModelHandoffError(
            f"Missing {label}. Supply its CLI option so emitted jobs contain no "
            "runtime path placeholders."
        )
    path = _repo_path(value, label=label)
    if not path.is_dir():
        raise RealModelHandoffError(f"{label} is not an existing directory: {path}")
    if sentinel_file and not (path / sentinel_file).is_file():
        raise RealModelHandoffError(
            f"{label} lacks expected {sentinel_file}: {path}"
        )
    return path


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def _ensure_writable_targets(paths: Sequence[Path], *, overwrite: bool) -> None:
    conflicts = sorted({path.resolve() for path in paths if path.exists()})
    if conflicts and not overwrite:
        rendered = "\n  - ".join(str(path) for path in conflicts)
        raise RealModelHandoffError(
            "Refusing to overwrite existing handoff artifacts. Use --overwrite only "
            f"after reviewing them:\n  - {rendered}"
        )


def _under(path: Path, root: Path, *, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise RealModelHandoffError(
            f"{label} escaped the selected output directory: {resolved}"
        ) from exc
    return resolved


def _standardized_jobs(
    *,
    rfantibody_native: Mapping[str, Any],
    iggm_native: Mapping[str, Any],
    germinal_native: Mapping[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    """Normalize three native job schemas without discarding their manifests."""

    rf_jobs: list[dict[str, Any]] = []
    rf_target_hash = str(rfantibody_native["target_pdb_sha256"])
    rf_framework_hashes = rfantibody_native["framework_pdb_sha256"]
    for index, native_job in enumerate(rfantibody_native["jobs"]):
        commands = [
            {
                "stage": command["stage"],
                "argv": list(command["argv"]),
                "working_directory": None,
                "expected_outputs": [
                    {"path": output} for output in command["expected_outputs"]
                ],
            }
            for command in native_job["commands"]
        ]
        rf_jobs.append(
            {
                "job_id": native_job["job_id"],
                "engine": "RFantibody",
                "profile": native_job["mode"],
                "geometry": "native_two_chain_paired_Fv",
                "execution_state": "planned_not_executed",
                "template_id": native_job["template_id"],
                "epitope_id": native_job["epitope_id"],
                "commands": commands,
                "expected_outputs": [
                    output
                    for command in commands
                    for output in command["expected_outputs"]
                ],
                "input_artifacts": [
                    {
                        "role": "target_pdb",
                        "path": native_job["target_pdb"],
                        "sha256": rf_target_hash,
                    },
                    {
                        "role": "framework_hlt_pdb",
                        "path": native_job["framework_hlt_pdb"],
                        "sha256": rf_framework_hashes[native_job["template_id"]],
                    },
                ],
                "unresolved_blockers": [],
                "native_job_reference": f"jobs[{index}]",
            }
        )

    iggm_jobs: list[dict[str, Any]] = []
    iggm_target_hash = str(iggm_native["target"]["pdb_sha256"])
    for index, native_job in enumerate(iggm_native["jobs"]):
        iggm_jobs.append(
            {
                "job_id": native_job["job_id"],
                "engine": "IgGM",
                "profile": native_job["profile"],
                "geometry": "native_two_chain_paired_VH_VL_plus_antigen",
                "execution_state": "planned_not_executed",
                "template_id": native_job["template_id"],
                "epitope_id": native_job["epitope_id"],
                "commands": [
                    {
                        "stage": "iggm_design",
                        "argv": list(native_job["command_argv"]),
                        "working_directory": native_job["working_directory"],
                        "expected_outputs": list(native_job["expected_outputs"]),
                    }
                ],
                "expected_outputs": list(native_job["expected_outputs"]),
                "input_artifacts": [
                    {
                        "role": "target_pdb",
                        "path": native_job["target_pdb_path"],
                        "sha256": iggm_target_hash,
                    },
                    {
                        "role": "masked_hla_fasta",
                        "path": native_job["input_fasta_path"],
                        "sha256": _text_sha256(native_job["input_fasta_content"]),
                    },
                ],
                "unresolved_blockers": [],
                "native_job_reference": f"jobs[{index}]",
            }
        )

    germinal_jobs: list[dict[str, Any]] = []
    for index, native_job in enumerate(germinal_native["jobs"]):
        target_pdb = Path(native_job["inputs"]["target_pdb_path"])
        scfv_pdb = Path(native_job["inputs"]["template_scfv_pdb_path"])
        normalized_outputs = [
            {"name": name, "path": value}
            for name, value in native_job["expected_outputs"].items()
        ]
        germinal_jobs.append(
            {
                "job_id": native_job["job_id"],
                "engine": "Germinal",
                "profile": native_job["profile"],
                "geometry": "single_chain_VH_linker_VL_scFv_separate_track",
                "execution_state": "planned_not_executed",
                "template_id": native_job["template_id"],
                "epitope_id": native_job["epitope_id"],
                "commands": [
                    {
                        "stage": "germinal_scfv_design_and_filter",
                        "argv": list(native_job["hydra_job"]["argv"]),
                        "working_directory": native_job["hydra_job"][
                            "working_directory"
                        ],
                        "expected_outputs": normalized_outputs,
                    }
                ],
                "expected_outputs": normalized_outputs,
                "input_artifacts": [
                    {
                        "role": "target_pdb",
                        "path": str(target_pdb),
                        "sha256": _file_sha256(target_pdb),
                    },
                    {
                        "role": "template_scfv_pdb",
                        "path": str(scfv_pdb),
                        "sha256": _file_sha256(scfv_pdb),
                    },
                    {
                        "role": "generated_target_yaml",
                        "path": native_job["target_yaml"]["handoff_path"],
                        "sha256": _text_sha256(native_job["target_yaml"]["content"]),
                    },
                ],
                "unresolved_blockers": [],
                "native_job_reference": f"jobs[{index}]",
            }
        )
    return {
        "RFantibody": rf_jobs,
        "IgGM": iggm_jobs,
        "Germinal": germinal_jobs,
    }


def _apply_execution_selection(
    standardized_jobs: Mapping[str, list[dict[str, Any]]],
    *,
    profile: str,
    requested_scope: str | None,
    canary_template_id: str | None,
    canary_epitope_id: str | None,
) -> dict[str, Any]:
    """Mark an auditable execution subset while retaining every planned job."""

    job_scope = requested_scope or ("canary" if profile == "smoke" else "all")
    if job_scope not in {"canary", "all"}:
        raise RealModelHandoffError(
            f"Unsupported job scope {job_scope!r}; expected 'canary' or 'all'"
        )
    if (canary_template_id is None) != (canary_epitope_id is None):
        raise RealModelHandoffError(
            "--canary-template-id and --canary-epitope-id must be supplied together"
        )
    if job_scope == "all" and canary_template_id is not None:
        raise RealModelHandoffError(
            "Canary identity options cannot be used with --job-scope all"
        )

    combinations_by_engine: dict[str, set[tuple[str, str]]] = {}
    for engine in REQUEST_FILENAMES:
        jobs = standardized_jobs.get(engine)
        if not isinstance(jobs, list) or not jobs:
            raise RealModelHandoffError(
                f"Standardized wrapper has no jobs for engine {engine}"
            )
        combinations: list[tuple[str, str]] = [
            (str(job["template_id"]), str(job["epitope_id"])) for job in jobs
        ]
        if len(set(combinations)) != len(combinations):
            raise RealModelHandoffError(
                f"Engine {engine} repeats a template-by-epitope job identity"
            )
        combinations_by_engine[engine] = set(combinations)

    reference_engine = next(iter(REQUEST_FILENAMES))
    reference_combinations = combinations_by_engine[reference_engine]
    for engine, combinations in combinations_by_engine.items():
        if combinations != reference_combinations:
            missing = sorted(reference_combinations - combinations)
            extra = sorted(combinations - reference_combinations)
            raise RealModelHandoffError(
                "All engines must expose the same template-by-epitope combinations "
                f"before execution selection; {engine} missing={missing}, extra={extra}"
            )

    canary_identity: tuple[str, str] | None = None
    if job_scope == "canary":
        if canary_template_id is not None and canary_epitope_id is not None:
            template_id = _require_nonempty_string(
                canary_template_id, label="canary template ID"
            )
            epitope_id = _require_nonempty_string(
                canary_epitope_id, label="canary epitope ID"
            )
            canary_identity = (template_id, epitope_id)
        else:
            canary_identity = sorted(reference_combinations)[0]
        missing_from = sorted(
            engine
            for engine, combinations in combinations_by_engine.items()
            if canary_identity not in combinations
        )
        if missing_from:
            raise RealModelHandoffError(
                "Requested canary template-by-epitope identity "
                f"{canary_identity!r} is missing from engines: {missing_from}"
            )

    selected_job_ids: dict[str, list[str]] = {}
    excluded_job_ids: dict[str, list[str]] = {}
    selected_job_counts: dict[str, int] = {}
    for engine, jobs in standardized_jobs.items():
        selected_ids: list[str] = []
        excluded_ids: list[str] = []
        for job in jobs:
            identity = (str(job["template_id"]), str(job["epitope_id"]))
            selected = job_scope == "all" or identity == canary_identity
            job["selected_for_execution"] = selected
            job["execution_disposition"] = (
                "selected_for_execution"
                if selected
                else "excluded_by_job_scope_do_not_execute"
            )
            if selected:
                selected_ids.append(str(job["job_id"]))
            else:
                excluded_ids.append(str(job["job_id"]))
        expected_selected_count = len(jobs) if job_scope == "all" else 1
        if len(selected_ids) != expected_selected_count:
            raise RealModelHandoffError(
                f"Internal execution-selection error for {engine}: expected "
                f"{expected_selected_count} selected jobs, observed {len(selected_ids)}"
            )
        selected_job_ids[engine] = selected_ids
        excluded_job_ids[engine] = excluded_ids
        selected_job_counts[engine] = len(selected_ids)

    return {
        "requested_job_scope": requested_scope or "auto",
        "resolved_job_scope": job_scope,
        "default_rule": "smoke->canary; full->all",
        "canary_identity": (
            {
                "template_id": canary_identity[0],
                "epitope_id": canary_identity[1],
                "selection_method": (
                    "explicit_cli"
                    if canary_template_id is not None
                    else "lexicographically_first_common_identity"
                ),
            }
            if canary_identity is not None
            else None
        ),
        "selected_job_count_per_engine": selected_job_counts,
        "selected_job_ids_by_engine": selected_job_ids,
        "excluded_job_ids_by_engine": excluded_job_ids,
        "executor_contract": {
            "authoritative_jobs": "engines[].execution_jobs",
            "native_plan_role": (
                "validation_and_provenance_only; native plans contain all 2x2 jobs "
                "and MUST NOT be submitted wholesale when resolved_job_scope=canary"
            ),
            "selection_rule": (
                "Execute only jobs with selected_for_execution=true. Never execute "
                "jobs marked excluded_by_job_scope_do_not_execute."
            ),
        },
    }


def compile_handoff(args: argparse.Namespace) -> dict[str, Any]:
    """Validate all inputs, build all plans, then materialize handoff artifacts."""

    request_dir = _repo_path(args.request_dir, label="request directory")
    if not request_dir.is_dir():
        raise RealModelHandoffError(
            f"Request directory does not exist or is not a directory: {request_dir}"
        )
    target_manifest_path = _repo_path(
        args.target_manifest, label="target structure manifest"
    )
    (
        requests,
        request_paths,
        design_request_index,
        source_run_metadata,
    ) = _validate_request_set(request_dir)
    target_manifest = _load_json_object(
        target_manifest_path, label="target structure manifest"
    )
    target = _validate_target_manifest(
        target_manifest, request=requests["RFantibody"]
    )
    if args.output_dir is None:
        safe_run_id = _safe_path_component(
            source_run_metadata["run_id"], label="source run_id"
        )
        safe_conformation_id = _safe_path_component(
            target["conformation_id"], label="target conformation_id"
        )
        output_root = (
            REPO_ROOT
            / "real_runs"
            / "handoffs"
            / safe_run_id
            / safe_conformation_id
            / args.profile
        ).resolve()
        output_path_source = "automatic_source_run_conformation_profile"
    else:
        output_root = _repo_path(args.output_dir, label="output directory")
        output_path_source = "explicit_cli"
    disposable_outputs_root = (REPO_ROOT / "outputs").resolve()
    if output_root == disposable_outputs_root or disposable_outputs_root in output_root.parents:
        raise RealModelHandoffError(
            "Real-model handoffs must not be written under the repository outputs/ "
            f"tree because reproducible workflow cleanup may remove it: {output_root}. "
            "Use the automatic real_runs/handoffs/<run-id>/<conformation>/<profile> "
            "path or another durable reviewed location."
        )

    iggm_repo = _required_runtime_directory(
        args.iggm_repo_dir, label="IgGM repository", sentinel_file="design.py"
    )
    germinal_repo = _required_runtime_directory(
        args.germinal_repo_dir,
        label="Germinal repository",
        sentinel_file="run_germinal.py",
    )
    germinal_af_params = _required_runtime_directory(
        args.germinal_af_params_dir,
        label="Germinal AlphaFold-Multimer parameter directory",
    )

    backend_paths: dict[str, str] = {}
    if args.profile == "full":
        for key in (
            "af3_repo_path",
            "af3_sif_path",
            "af3_model_dir",
            "af3_db_dir",
            "msa_db_dir",
        ):
            raw_value = getattr(args, key)
            if raw_value is None:
                raise RealModelHandoffError(
                    f"--{key.replace('_', '-')} is required for the Germinal full/AF3 profile"
                )
            path = _repo_path(raw_value, label=f"Germinal backend {key}")
            if not path.exists():
                raise RealModelHandoffError(
                    f"Germinal backend {key} does not exist: {path}"
                )
            backend_paths[key] = str(path)

    rfantibody_request = _inject_hotspots(
        requests["RFantibody"],
        target["hotspots"],
        field="rfantibody_hotspot_residue_indices",
    )
    germinal_request = _inject_hotspots(
        requests["Germinal"],
        target["hotspots"],
        field="germinal_hotspot_residue_indices",
    )

    rfantibody_root = output_root / "rfantibody"
    iggm_root = output_root / "iggm"
    germinal_root = output_root / "germinal"

    # Build every adapter output before creating a directory or writing a file.
    # This prevents a late coordinate error from leaving a misleading partial
    # handoff that looks ready to run.
    rfantibody_plan = build_rfantibody_plan(
        rfantibody_request,
        target_pdb=target["target_pdb"],
        framework_hlt_pdbs=target["rfantibody_hlt"],
        full_coordinate_to_pdb=target["full_to_pdb"],
        output_root=rfantibody_root / "results",
        mode=args.profile,
        seed=args.seed,
        command_prefix=tuple(args.rfantibody_command_prefix),
        runtime_ref=args.rfantibody_runtime_ref,
    )
    iggm_handoff = build_iggm_jobs(
        requests["IgGM"],
        target_pdb_path=target["target_pdb"],
        pdb_antigen_sequence=target["antigen_sequence"],
        pdb_antigen_chain=target["target_chain"],
        full_to_local_residue_map=target["full_to_local"],
        profile=args.profile,
        iggm_repo_dir=iggm_repo,
        python_executable=args.iggm_python,
        input_dir=iggm_root / "inputs",
        output_dir=iggm_root / "results",
    )
    germinal_handoff = build_germinal_jobs(
        germinal_request,
        target_pdb_path=target["target_pdb"],
        template_scfv_pdbs=target["germinal_scfv"],
        target_residue_map=target["full_to_pdb"],
        handoff_root=germinal_root,
        profile=args.profile,
        target_chain=target["target_chain"],
        binder_chain=args.germinal_binder_chain,
        scfv_chain=args.germinal_scfv_chain,
        germinal_repo_dir=germinal_repo,
        af_params_dir=germinal_af_params,
        backend_paths=backend_paths,
        python_executable=args.germinal_python,
    )

    unresolved = sorted(
        {
            placeholder
            for job in germinal_handoff["jobs"]
            for placeholder in job["hydra_job"]["unresolved_placeholders"]
        }
    )
    if unresolved:
        raise RealModelHandoffError(
            "Germinal plan still contains unresolved runtime placeholders after "
            f"validation: {unresolved}"
        )

    rfantibody_native = rfantibody_plan.as_dict()
    standardized_jobs = _standardized_jobs(
        rfantibody_native=rfantibody_native,
        iggm_native=iggm_handoff,
        germinal_native=germinal_handoff,
    )
    execution_selection = _apply_execution_selection(
        standardized_jobs,
        profile=args.profile,
        requested_scope=args.job_scope,
        canary_template_id=args.canary_template_id,
        canary_epitope_id=args.canary_epitope_id,
    )

    rfantibody_plan_path = rfantibody_root / "rfantibody_plan.json"
    iggm_plan_path = iggm_root / "iggm_jobs.json"
    germinal_plan_path = germinal_root / "germinal_jobs.json"
    unified_path = output_root / "unified_handoff_manifest.json"
    normalized_request_paths = {
        "RFantibody": rfantibody_root / "normalized_request.json",
        "IgGM": iggm_root / "normalized_request.json",
        "Germinal": germinal_root / "normalized_request.json",
    }
    fasta_paths = [
        _under(
            Path(job["input_fasta_path"]),
            output_root,
            label=f"IgGM FASTA for {job['job_id']}",
        )
        for job in iggm_handoff["jobs"]
    ]
    germinal_yaml_paths = [
        germinal_root
        / "target_configs"
        / f"{job['target_yaml']['config_name']}.yaml"
        for job in germinal_handoff["jobs"]
    ]
    germinal_staged_pdbs = [
        germinal_root / "jobs" / str(job["job_id"]) / "pdbs" / "scfv.pdb"
        for job in germinal_handoff["jobs"]
    ]
    planned_writes = [
        rfantibody_plan_path,
        iggm_plan_path,
        germinal_plan_path,
        unified_path,
        *normalized_request_paths.values(),
        *fasta_paths,
        *germinal_yaml_paths,
        *germinal_staged_pdbs,
    ]
    _ensure_writable_targets(planned_writes, overwrite=args.overwrite)

    # Materialization starts only after all three adapters and every destination
    # have passed validation.  The emitted argv arrays remain data; they are not
    # passed to subprocess or a shell here.
    _write_json(rfantibody_plan_path, rfantibody_native)
    _write_json(iggm_plan_path, iggm_handoff)
    _write_json(normalized_request_paths["RFantibody"], rfantibody_request)
    _write_json(normalized_request_paths["IgGM"], requests["IgGM"])
    _write_json(normalized_request_paths["Germinal"], germinal_request)
    for job, fasta_path in zip(iggm_handoff["jobs"], fasta_paths, strict=True):
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        fasta_path.write_text(str(job["input_fasta_content"]), encoding="utf-8")
    germinal_written = write_germinal_handoff(
        germinal_handoff,
        germinal_root,
        stage_template_pdbs=True,
        overwrite=args.overwrite,
    )

    source_hashes = {
        engine: {
            "path": str(path),
            "sha256": _file_sha256(path),
        }
        for engine, path in request_paths.items()
    }
    source_integrity = {
        "state": (
            "verified_by_design_request_index"
            if design_request_index["present"]
            else "unverified_missing_design_request_index"
        ),
        "ready_for_execution": bool(design_request_index["present"]),
        "blockers": (
            []
            if design_request_index["present"]
            else [
                "No design_request_index.json was present; request file hashes were "
                "not authenticated against a workflow-generated index. Regenerate the "
                "request set before external execution."
            ]
        ),
    }
    compiler_provenance = _compiler_provenance()
    generated_at = datetime.now().astimezone().isoformat()
    invocation_config = {
        "profile": args.profile,
        "requested_job_scope": args.job_scope,
        "resolved_job_scope": execution_selection["resolved_job_scope"],
        "canary_template_id": args.canary_template_id,
        "canary_epitope_id": args.canary_epitope_id,
        "resolved_canary_identity": execution_selection["canary_identity"],
        "seed": args.seed,
        "rfantibody_command_prefix": list(args.rfantibody_command_prefix),
        "rfantibody_runtime_ref": args.rfantibody_runtime_ref,
        "iggm_python": args.iggm_python,
        "germinal_python": args.germinal_python,
        "germinal_binder_chain": args.germinal_binder_chain,
        "germinal_scfv_chain": args.germinal_scfv_chain,
        "iggm_repo_dir": str(iggm_repo),
        "germinal_repo_dir": str(germinal_repo),
        "germinal_af_params_dir": str(germinal_af_params),
        "germinal_backend_paths": dict(sorted(backend_paths.items())),
        "output_root": str(output_root),
        "output_path_source": output_path_source,
    }
    native_plan_semantic_sha256_by_engine = {
        "RFantibody": _json_sha256(rfantibody_native),
        "IgGM": _json_sha256(iggm_handoff),
        "Germinal": _json_sha256(germinal_handoff),
    }
    standardized_jobs_semantic_sha256_by_engine = {
        engine: _json_sha256(jobs)
        for engine, jobs in standardized_jobs.items()
    }
    handoff_identity = {
        "source_run_id": source_run_metadata["run_id"],
        "request_sha256_by_engine": {
            engine: value["sha256"] for engine, value in source_hashes.items()
        },
        "target_manifest_sha256": _file_sha256(target_manifest_path),
        "target_pdb_sha256": _file_sha256(target["target_pdb"]),
        "invocation_config": invocation_config,
        "native_plan_semantic_sha256_by_engine": (
            native_plan_semantic_sha256_by_engine
        ),
        "standardized_jobs_semantic_sha256_by_engine": (
            standardized_jobs_semantic_sha256_by_engine
        ),
        "compiler_source_sha256": compiler_provenance["compiler_source"]["sha256"],
        "adapter_source_sha256_by_engine": {
            engine: value["sha256"]
            for engine, value in compiler_provenance[
                "adapter_source_by_engine"
            ].items()
        },
    }
    identity_sha256 = _json_sha256(handoff_identity)
    handoff_id = f"nfl_handoff_{identity_sha256[:24]}"
    unified_manifest = {
        "schema": UNIFIED_HANDOFF_SCHEMA,
        "handoff_id": handoff_id,
        "handoff_generated_at": generated_at,
        "identity_sha256": identity_sha256,
        "handoff_identity": handoff_identity,
        "invocation_config": invocation_config,
        "native_plan_semantic_sha256_by_engine": (
            native_plan_semantic_sha256_by_engine
        ),
        "standardized_jobs_semantic_sha256_by_engine": (
            standardized_jobs_semantic_sha256_by_engine
        ),
        "compiler": compiler_provenance,
        "execution_state": "planned_not_executed",
        "does_not_execute_external_models": True,
        "profile": args.profile,
        "job_scope": execution_selection["resolved_job_scope"],
        "seed": args.seed,
        "run_id": source_run_metadata["run_id"],
        "source_run": {
            "run_id": source_run_metadata["run_id"],
            "run_metadata": source_run_metadata,
            "design_request_index": design_request_index,
            "source_integrity": source_integrity,
        },
        "handoff_location": {
            "path": str(output_root),
            "path_source": output_path_source,
            "automatic_layout": (
                "real_runs/handoffs/<source-run-id>/<conformation-id>/<profile>"
            ),
        },
        "repository_root": str(REPO_ROOT),
        "target_structure": {
            "manifest_path": str(target_manifest_path),
            "manifest_schema": TARGET_MANIFEST_SCHEMA,
            "manifest_execution_state": target["manifest_state"],
            "review": target["review"],
            "manifest_sha256": _file_sha256(target_manifest_path),
            "conformation_id": target["conformation_id"],
            "oligomeric_state": target["oligomeric_state"],
            "target_pdb": str(target["target_pdb"]),
            "target_pdb_sha256": _file_sha256(target["target_pdb"]),
            "target_chain": target["target_chain"],
            "curated_hotspots_by_epitope": target["hotspots"],
        },
        "normalized_request_sources": source_hashes,
        "injected_request_sha256": {
            "RFantibody": _json_sha256(rfantibody_request),
            "IgGM": _json_sha256(requests["IgGM"]),
            "Germinal": _json_sha256(germinal_request),
        },
        "execution_selection": execution_selection,
        "engines": [
            {
                "engine": "RFantibody",
                "profile": args.profile,
                "geometry": "native_two_chain_paired_Fv",
                "execution_state": "planned_not_executed",
                "adapter_schema": rfantibody_native["schema"],
                "required_upstream_revision": args.rfantibody_runtime_ref,
                "runtime_verification_state": "unverified_not_executed",
                "ready_for_execution": False,
                "request_source": source_hashes["RFantibody"],
                "native_plan": {
                    "path": str(rfantibody_plan_path),
                    "schema": rfantibody_native["schema"],
                    "execution_state": rfantibody_native["execution_state"],
                    "submission_allowed": False,
                    "reason": "Contains all validation jobs; use execution_jobs only.",
                },
                "planned_job_count": len(rfantibody_plan.jobs),
                "selected_job_count": execution_selection[
                    "selected_job_count_per_engine"
                ]["RFantibody"],
                "selected_job_ids": execution_selection[
                    "selected_job_ids_by_engine"
                ]["RFantibody"],
                "excluded_job_ids": execution_selection[
                    "excluded_job_ids_by_engine"
                ]["RFantibody"],
                "execution_jobs": [
                    job
                    for job in standardized_jobs["RFantibody"]
                    if job["selected_for_execution"]
                ],
                "jobs": standardized_jobs["RFantibody"],
                "unresolved_blockers": [
                    "Verify installed RFantibody revision and checkpoint hashes before execution."
                ],
            },
            {
                "engine": "IgGM",
                "profile": args.profile,
                "geometry": "native_two_chain_paired_VH_VL_plus_antigen",
                "execution_state": "planned_not_executed",
                "adapter_schema": iggm_handoff["schema"],
                "required_upstream_revision": iggm_handoff["upstream"]["commit"],
                "runtime_verification_state": "unverified_not_executed",
                "ready_for_execution": False,
                "request_source": source_hashes["IgGM"],
                "native_plan": {
                    "path": str(iggm_plan_path),
                    "schema": iggm_handoff["schema"],
                    "execution_state": iggm_handoff["execution_state"],
                    "submission_allowed": False,
                    "reason": "Contains all validation jobs; use execution_jobs only.",
                },
                "planned_job_count": int(iggm_handoff["job_count"]),
                "selected_job_count": execution_selection[
                    "selected_job_count_per_engine"
                ]["IgGM"],
                "selected_job_ids": execution_selection[
                    "selected_job_ids_by_engine"
                ]["IgGM"],
                "excluded_job_ids": execution_selection[
                    "excluded_job_ids_by_engine"
                ]["IgGM"],
                "execution_jobs": [
                    job
                    for job in standardized_jobs["IgGM"]
                    if job["selected_for_execution"]
                ],
                "jobs": standardized_jobs["IgGM"],
                "unresolved_blockers": [
                    "Verify installed IgGM git revision and checkpoint hashes before execution."
                ],
            },
            {
                "engine": "Germinal",
                "profile": args.profile,
                "geometry": "single_chain_VH_linker_VL_scFv_separate_track",
                "execution_state": "planned_not_executed",
                "adapter_schema": germinal_handoff["schema"],
                "required_upstream_revision": germinal_handoff["upstream"]["pinned_commit"],
                "runtime_verification_state": "unverified_not_executed",
                "ready_for_execution": False,
                "request_source": source_hashes["Germinal"],
                "native_plan": {
                    "path": str(germinal_plan_path),
                    "schema": germinal_handoff["schema"],
                    "execution_state": germinal_handoff["execution_state"],
                    "submission_allowed": False,
                    "reason": "Contains all validation jobs; use execution_jobs only.",
                },
                "planned_job_count": int(germinal_handoff["job_count"]),
                "selected_job_count": execution_selection[
                    "selected_job_count_per_engine"
                ]["Germinal"],
                "selected_job_ids": execution_selection[
                    "selected_job_ids_by_engine"
                ]["Germinal"],
                "excluded_job_ids": execution_selection[
                    "excluded_job_ids_by_engine"
                ]["Germinal"],
                "execution_jobs": [
                    job
                    for job in standardized_jobs["Germinal"]
                    if job["selected_for_execution"]
                ],
                "jobs": standardized_jobs["Germinal"],
                "unresolved_blockers": [
                    "Verify installed Germinal revision and all backend/model parameter hashes before execution."
                ],
            },
        ],
        "outputs": {
            "RFantibody": {
                "plan": str(rfantibody_plan_path),
                "normalized_request": str(normalized_request_paths["RFantibody"]),
                "job_count": len(rfantibody_plan.jobs),
                "selected_job_count": execution_selection[
                    "selected_job_count_per_engine"
                ]["RFantibody"],
                "command_count": len(rfantibody_plan.commands),
                "expected_candidate_count": rfantibody_plan.expected_candidate_count,
            },
            "IgGM": {
                "plan": str(iggm_plan_path),
                "normalized_request": str(normalized_request_paths["IgGM"]),
                "input_fastas": [str(path) for path in fasta_paths],
                "job_count": int(iggm_handoff["job_count"]),
                "selected_job_count": execution_selection[
                    "selected_job_count_per_engine"
                ]["IgGM"],
            },
            "Germinal": {
                "plan": str(germinal_plan_path),
                "normalized_request": str(normalized_request_paths["Germinal"]),
                "target_yamls": germinal_written["target_yamls"],
                "staged_template_pdbs": germinal_written["staged_template_pdbs"],
                "job_count": int(germinal_handoff["job_count"]),
                "selected_job_count": execution_selection[
                    "selected_job_count_per_engine"
                ]["Germinal"],
            },
        },
        "execution_contract": [
            "Only engines[].execution_jobs is an executable selection; engines[].jobs and native plans include validation-only unselected jobs.",
            "Never execute a job whose selected_for_execution flag is false or whose execution_disposition is excluded_by_job_scope_do_not_execute.",
            "Before execution, resolve every engine unresolved_blocker and independently record installed revision/checkpoint hashes; this planner does not verify runtimes.",
            "Every command is stored as an argv array for later review.",
            "This compiler did not invoke a shell, subprocess, model, installer, or network client.",
            "Smoke plans validate integration only and are not scientific results.",
            "Generated files remain planned inputs until separately executed and provenance-captured.",
        ],
    }
    _write_json(unified_path, unified_manifest)
    return unified_manifest


def _self_test_write_pdb(
    path: Path, sequence: str, *, chain: str = "A", start: int = 1
) -> None:
    lines: list[str] = []
    for serial, (offset, amino_acid) in enumerate(enumerate(sequence), start=1):
        residue_number = start + offset
        lines.append(
            f"ATOM  {serial:5d}  CA  {SELF_TEST_AA1_TO_AA3[amino_acid]:>3s} "
            f"{chain}{residue_number:4d}    {float(serial):8.3f}"
            f"{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        )
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def _self_test_write_hlt(path: Path) -> None:
    cdr_names = ("H1", "H2", "H3", "L1", "L2", "L3")
    remarks = [
        f"REMARK PDBinfo-LABEL: {index} {name}"
        for index, name in enumerate(cdr_names, start=1)
    ]
    atoms = [
        "ATOM      1  CA  ALA H   1       1.000   0.000   0.000  1.00 20.00           C",
        "ATOM      2  CA  GLY L   1       2.000   0.000   0.000  1.00 20.00           C",
    ]
    path.write_text("\n".join([*remarks, *atoms, "END", ""]), encoding="utf-8")


def _self_test_request(engine: str, run_metadata: Mapping[str, Any]) -> dict[str, Any]:
    cdr_names = ("H1", "H2", "H3", "L1", "L2", "L3")
    regions = [
        {
            "region": name,
            "chain": "VH" if name.startswith("H") else "VL",
            "start_1_based": position,
            "end_1_based_inclusive": position,
            "length_aa": 1,
        }
        for name, position in zip(cdr_names, (2, 5, 8, 2, 5, 8), strict=True)
    ]
    return {
        "schema": NORMALIZED_REQUEST_SCHEMA,
        "campaign_mode": "paired_Fv_six_CDR_de_novo_design",
        "execution_state": "not_run",
        "result_provenance": "adapter_request_only",
        "engine": engine,
        "run_metadata": dict(run_metadata),
        "antigen": {
            "protein": "self_test_antigen",
            "full_sequence": "ACDEFGHIKLMN",
            "antigen_pdb_path": "",
        },
        "epitopes": [
            {
                "epitope_id": "epitope_N",
                "sequence": "CD",
                "start_1_based": 2,
                "end_1_based_inclusive": 3,
                "candidate_hotspot_residue_indices": [2, 3],
            },
            {
                "epitope_id": "epitope_C",
                "sequence": "KL",
                "start_1_based": 9,
                "end_1_based_inclusive": 10,
                "candidate_hotspot_residue_indices": [9, 10],
            },
        ],
        "templates": [
            {
                "template_id": "template_A",
                "framework_source_id": "source_A",
                "template_role": "framework_source_only",
                "masked_vh": "AXAAXAAXA",
                "masked_vl": "GXGGXGGXG",
                "design_regions": regions,
                "designed_regions": ";".join(cdr_names),
            },
            {
                "template_id": "template_B",
                "framework_source_id": "source_B",
                "template_role": "framework_source_only",
                "masked_vh": "CXCCXCCXC",
                "masked_vl": "EXEEXEEXE",
                "design_regions": regions,
                "designed_regions": ";".join(cdr_names),
            },
        ],
    }


def run_self_test() -> dict[str, Any]:
    """Run a fresh-directory compiler regression without executing any model."""

    with tempfile.TemporaryDirectory(prefix="nfl_real_handoff_self_test_") as raw_root:
        root = Path(raw_root)
        request_dir = root / "requests"
        request_dir.mkdir()
        run_metadata = {
            "run_id": "self_test_run_20260813T0000000000000800",
            "generated_at": "2026-08-13T00:00:00+08:00",
            "nfl_ab_design_version": PACKAGE_VERSION,
            "design_campaign_sha256": "a" * 64,
            "workflow_source_sha256": "b" * 64,
            "design_pipeline_source_sha256": "c" * 64,
        }
        request_paths: dict[str, Path] = {}
        for engine in REQUEST_FILENAMES:
            path = request_dir / REQUEST_FILENAMES[engine]
            _write_json(path, _self_test_request(engine, run_metadata))
            request_paths[engine] = path
        index_path = request_dir / "design_request_index.json"
        _write_json(
            index_path,
            {
                "schema": DESIGN_REQUEST_INDEX_SCHEMA,
                "execution_state": "not_run",
                "run_metadata": dict(run_metadata),
                "engines": list(REQUEST_FILENAMES),
                "request_files": [
                    str(request_paths[engine]) for engine in REQUEST_FILENAMES
                ],
                "request_sha256_by_engine": {
                    engine: _file_sha256(path)
                    for engine, path in request_paths.items()
                },
            },
        )

        antigen = "ACDEFGHIKLMN"
        target_pdb = root / "target.pdb"
        _self_test_write_pdb(target_pdb, antigen)
        hlt_pdbs: dict[str, str] = {}
        scfv_pdbs: dict[str, str] = {}
        scfv_sequences = {
            "template_A": "ASAASAASA" + "GGGGS" + "GSGGSGGSG",
            "template_B": "CSCCSCCSC" + "GGGGS" + "ESEESEESE",
        }
        for template_id, scfv_sequence in scfv_sequences.items():
            hlt_path = root / f"{template_id}.hlt.pdb"
            scfv_path = root / f"{template_id}.scfv.pdb"
            _self_test_write_hlt(hlt_path)
            _self_test_write_pdb(scfv_path, scfv_sequence)
            hlt_pdbs[template_id] = str(hlt_path)
            scfv_pdbs[template_id] = str(scfv_path)

        iggm_repo = root / "IgGM"
        germinal_repo = root / "germinal"
        af_params = root / "af_params"
        iggm_repo.mkdir()
        germinal_repo.mkdir()
        af_params.mkdir()
        (iggm_repo / "design.py").write_text("", encoding="utf-8")
        (germinal_repo / "run_germinal.py").write_text("", encoding="utf-8")
        target_manifest_path = root / "target_manifest.json"
        _write_json(
            target_manifest_path,
            {
                "schema": TARGET_MANIFEST_SCHEMA,
                "execution_state": "reviewed_ready_for_handoff",
                "review": {
                    "reviewed_by": "built_in_self_test",
                    "reviewed_at": "2026-08-13T00:00:00+08:00",
                    "contracts_acknowledged": True,
                },
                "conformation_id": "synthetic_monomer_test",
                "oligomeric_state": "single_chain_monomer",
                "target_pdb_path": str(target_pdb),
                "target_chain": "A",
                "antigen_sequence_in_pdb_order": antigen,
                "full_coordinate_to_pdb": {
                    str(position): {
                        "chain_id": "A",
                        "residue_number": position,
                    }
                    for position in range(1, len(antigen) + 1)
                },
                "full_coordinate_to_local_1_based": {
                    str(position): position
                    for position in range(1, len(antigen) + 1)
                },
                "selected_hotspots_by_epitope": {
                    "epitope_N": [2],
                    "epitope_C": [9],
                },
                "framework_coordinate_inputs": {
                    "rfantibody_hlt_pdbs": hlt_pdbs,
                    "germinal_scfv_pdbs": scfv_pdbs,
                },
            },
        )
        output_dir = root / "fresh_handoff"
        arguments = build_parser().parse_args(
            [
                "--request-dir",
                str(request_dir),
                "--target-manifest",
                str(target_manifest_path),
                "--output-dir",
                str(output_dir),
                "--profile",
                "smoke",
                "--iggm-repo-dir",
                str(iggm_repo),
                "--germinal-repo-dir",
                str(germinal_repo),
                "--germinal-af-params-dir",
                str(af_params),
            ]
        )
        result = compile_handoff(arguments)
        if result["job_scope"] != "canary":
            raise RealModelHandoffError("Self-test smoke profile did not select canary scope")
        if any(
            engine["planned_job_count"] != 4
            or engine["selected_job_count"] != 1
            for engine in result["engines"]
        ):
            raise RealModelHandoffError(
                "Self-test wrapper did not retain four plans/select one canary per engine"
            )
        if not (output_dir / "germinal" / "germinal_jobs.json").is_file():
            raise RealModelHandoffError(
                "Self-test did not materialize Germinal manifest on a fresh output"
            )
        different_seed_output = root / "different_seed_handoff"
        different_seed_arguments = build_parser().parse_args(
            [
                "--request-dir",
                str(request_dir),
                "--target-manifest",
                str(target_manifest_path),
                "--output-dir",
                str(different_seed_output),
                "--profile",
                "smoke",
                "--seed",
                "20260814",
                "--iggm-repo-dir",
                str(iggm_repo),
                "--germinal-repo-dir",
                str(germinal_repo),
                "--germinal-af-params-dir",
                str(af_params),
            ]
        )
        different_seed_result = compile_handoff(different_seed_arguments)
        if different_seed_result["handoff_id"] == result["handoff_id"]:
            raise RealModelHandoffError(
                "Self-test changing the campaign seed did not change handoff_id"
            )
        if different_seed_result["identity_sha256"] == result["identity_sha256"]:
            raise RealModelHandoffError(
                "Self-test changing the campaign seed did not change identity_sha256"
            )
        try:
            compile_handoff(arguments)
        except RealModelHandoffError as exc:
            if "Refusing to overwrite" not in str(exc):
                raise
        else:
            raise RealModelHandoffError(
                "Self-test second no-overwrite compile unexpectedly succeeded"
            )
        request_paths["IgGM"].write_text(
            request_paths["IgGM"].read_text(encoding="utf-8") + "\n",
            encoding="utf-8",
        )
        try:
            _validate_request_set(request_dir)
        except RealModelHandoffError as exc:
            if "hash mismatch for IgGM" not in str(exc):
                raise
        else:
            raise RealModelHandoffError(
                "Self-test tampered request passed indexed SHA-256 validation"
            )
        return {
            "execution_state": "self_test_passed",
            "fresh_compile_without_overwrite": True,
            "second_compile_refused_overwrite": True,
            "request_hash_tamper_rejected": True,
            "different_seed_changes_handoff_id": True,
            "smoke_canary_selected_job_count_per_engine": 1,
            "does_not_execute_external_models": True,
        }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compile reviewed target/template structures into RFantibody, IgGM, "
            "and Germinal handoffs without executing any external model. Relative "
            "paths are resolved from the repository root."
        )
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help=(
            "Run a temporary fresh-output/no-overwrite/index-tamper regression; "
            "does not execute external models"
        ),
    )
    parser.add_argument(
        "--request-dir",
        default="outputs/exports/design_requests",
        help="Directory containing the three normalized request JSON files",
    )
    parser.add_argument(
        "--target-manifest",
        default="config/target_structure_manifest.json",
        help="Reviewed target-structure manifest (copy and fill the .example.json)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Durable destination for plans and staged inputs. Dynamic default: "
            "real_runs/handoffs/<source-run-id>/<conformation-id>/<profile>. "
            "Paths under repo outputs/ are rejected."
        ),
    )
    parser.add_argument(
        "--profile",
        choices=("smoke", "full"),
        default="smoke",
        help="Adapter profile to compile; neither choice executes a model",
    )
    parser.add_argument(
        "--job-scope",
        choices=("canary", "all"),
        default=None,
        help=(
            "Execution selection recorded in the wrapper. Default is canary for "
            "--profile smoke and all for --profile full; all 2x2 jobs are still "
            "validated and retained in native plans."
        ),
    )
    parser.add_argument(
        "--canary-template-id",
        help=(
            "Explicit canary template ID; requires --canary-epitope-id and a "
            "canary job scope"
        ),
    )
    parser.add_argument(
        "--canary-epitope-id",
        help=(
            "Explicit canary epitope ID; requires --canary-template-id and a "
            "canary job scope"
        ),
    )
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument(
        "--iggm-repo-dir",
        help="Existing pinned IgGM checkout containing design.py (required)",
    )
    parser.add_argument(
        "--germinal-repo-dir",
        help="Existing pinned Germinal checkout containing run_germinal.py (required)",
    )
    parser.add_argument(
        "--germinal-af-params-dir",
        help="Existing AlphaFold-Multimer parameter directory required by Germinal",
    )
    parser.add_argument(
        "--germinal-binder-chain",
        default="B",
        help="Binder chain ID in Germinal's planned complex (default: B)",
    )
    parser.add_argument(
        "--germinal-scfv-chain",
        default="A",
        help="Chain ID used by each input scFv PDB (default: A)",
    )
    parser.add_argument(
        "--rfantibody-command-prefix",
        action="append",
        default=[],
        metavar="ARGV_TOKEN",
        help=(
            "Repeat for each optional argv prefix token, e.g. "
            "--rfantibody-command-prefix uv --rfantibody-command-prefix run"
        ),
    )
    parser.add_argument(
        "--rfantibody-runtime-ref",
        default=RFANTIBODY_OFFICIAL_MAIN_SHA,
        help="RFantibody git SHA/tag recorded in plan provenance",
    )
    parser.add_argument(
        "--iggm-python",
        default="python",
        help="Executable token placed in IgGM argv arrays",
    )
    parser.add_argument(
        "--germinal-python",
        default="python",
        help="Executable token placed in Germinal argv arrays",
    )
    parser.add_argument("--af3-repo-path", help="Required for --profile full")
    parser.add_argument("--af3-sif-path", help="Required for --profile full")
    parser.add_argument("--af3-model-dir", help="Required for --profile full")
    parser.add_argument("--af3-db-dir", help="Required for --profile full")
    parser.add_argument("--msa-db-dir", help="Required for --profile full")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite only the known generated artifacts; never delete directories",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.self_test:
        try:
            result = run_self_test()
        except (
            RealModelHandoffError,
            RFantibodyAdapterError,
            IgGMAdapterError,
            GerminalAdapterError,
            OSError,
        ) as exc:
            print(f"SELF-TEST ERROR: {exc}", file=sys.stderr)
            return 2
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    try:
        result = compile_handoff(args)
    except (
        RealModelHandoffError,
        RFantibodyAdapterError,
        IgGMAdapterError,
        GerminalAdapterError,
        OSError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    output = Path(result["outputs"]["RFantibody"]["plan"]).parents[1]
    print(
        f"Prepared {result['profile']}/{result['job_scope']} handoff for "
        "RFantibody, IgGM, and Germinal: "
        f"{output / 'unified_handoff_manifest.json'}"
    )
    print("No external model was executed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
