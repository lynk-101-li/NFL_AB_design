"""Compile normalized NFL antibody-design requests into pinned IgGM jobs.

This module is deliberately an adapter, not an IgGM runner.  It validates the
biological coordinate boundary, emits the exact three-chain FASTA text expected
by IgGM, and produces shell-free CLI argument vectors.  It never downloads a
checkpoint, imports IgGM, starts a subprocess, or claims that a model ran.

The adapter targets the official TencentAI4S/IgGM ``master`` revision recorded
in :data:`UPSTREAM`.  IgGM uses ``X`` residues as its design mask and expects the
antigen to be the final FASTA chain.  Its ``--epitope`` values are local,
one-based residue indices in that antigen chain, not UniProt coordinates.
"""

from __future__ import annotations

import hashlib
import json
import re
import shlex
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


NORMALIZED_REQUEST_SCHEMA = "nfl_ab_design.normalized_de_novo_request.v1"
HANDOFF_SCHEMA = "nfl_ab_design.iggm_handoff.v1"
DESIGN_REGIONS = ("H1", "H2", "H3", "L1", "L2", "L3")
CANONICAL_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Pinned after inspecting the official repository on 2026-08-12.  A caller
# should check out this exact commit before executing any emitted command.
UPSTREAM: dict[str, Any] = {
    "name": "TencentAI4S/IgGM",
    "repository": "https://github.com/TencentAI4S/IgGM",
    "branch": "master",
    "commit": "06abc563b3fc8c7ea020543add16b69b6f8a1c8d",
    "commit_url": (
        "https://github.com/TencentAI4S/IgGM/commit/"
        "06abc563b3fc8c7ea020543add16b69b6f8a1c8d"
    ),
    "license": "MIT",
    "license_url": "https://github.com/TencentAI4S/IgGM/blob/master/LICENSE",
    "official_readme": "https://github.com/TencentAI4S/IgGM/blob/master/README.md",
    "official_cli": "https://github.com/TencentAI4S/IgGM/blob/master/design.py",
    "official_all_cdr_example": (
        "https://github.com/TencentAI4S/IgGM/blob/master/examples/"
        "fasta.files.design/8hpu_M_N_A/8hpu_M_N_A_CDR_All.fasta"
    ),
    "weights_record": "https://zenodo.org/records/16909543",
    "environment": {
        "python": "3.10.14",
        "pytorch": "2.0.1",
        "cuda_runtime": "11.7",
    },
    "required_design_weights": [
        {
            "filename": "esm_ppi_650m_ab.pth",
            "size": "2.6 GB",
            "md5": "9f332b21296d8182c6159ba7833d3a74",
        },
        {
            "filename": "antibody_design_trunk.pth",
            "size": "101.7 MB",
            "md5": "975baa1f0f5d9ae5cb7afdd4ed179da7",
        },
        {
            "filename": "igso3_buffer.pth",
            "size": "371.4 MB",
            "md5": "8963fa425002a5a65c0b13ddaa443e9e",
        },
    ],
}

# ``smoke`` minimizes sampling while exercising model loading, sequence design,
# and PDB/FASTA export.  ``full`` matches this repository's planned campaign of
# 24 designs per template-epitope pair.  Both retain the upstream chunk and
# temperature defaults and omit optional PyRosetta relaxation.
PROFILE_PARAMETERS: dict[str, dict[str, Any]] = {
    "smoke": {
        "run_task": "design",
        "steps": 1,
        "num_samples": 1,
        "chunk_size": 64,
        "temperature": 1.0,
        "max_antigen_size": 384,
        "relax": False,
        "purpose": "minimal end-to-end checkpoint/input/output validation",
    },
    "full": {
        "run_task": "design",
        "steps": 10,
        "num_samples": 24,
        "chunk_size": 64,
        "temperature": 1.0,
        "max_antigen_size": 384,
        "relax": False,
        "purpose": "NFL campaign generation (24 samples per template-epitope pair)",
    },
}


class IgGMAdapterError(ValueError):
    """Raised when a normalized request cannot be mapped safely to IgGM."""


def load_normalized_request(path: str | Path) -> dict[str, Any]:
    """Load a normalized request JSON object without changing it."""

    request_path = Path(path)
    try:
        data = json.loads(request_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise IgGMAdapterError(f"Cannot read normalized request: {request_path}") from exc
    except json.JSONDecodeError as exc:
        raise IgGMAdapterError(f"Invalid normalized request JSON: {request_path}: {exc}") from exc
    if not isinstance(data, dict):
        raise IgGMAdapterError("Normalized request JSON must contain one object")
    return data


def build_iggm_jobs(
    request: Mapping[str, Any],
    *,
    target_pdb_path: str | Path,
    pdb_antigen_sequence: str,
    pdb_antigen_chain: str,
    full_to_local_residue_map: Mapping[int | str, int | str],
    profile: str = "smoke",
    iggm_repo_dir: str | Path = ".",
    python_executable: str = "python",
    input_dir: str | Path = "real_runs/adapter_plans/iggm/inputs",
    output_dir: str | Path = "real_runs/adapter_plans/iggm/results",
) -> dict[str, Any]:
    """Build deterministic, non-executed IgGM CLI job specifications.

    Args:
        request: ``nfl_ab_design.normalized_de_novo_request.v1`` object.
        target_pdb_path: Existing PDB file containing the antigen chain.
        pdb_antigen_sequence: Explicit sequence represented by that PDB chain.
        pdb_antigen_chain: Explicit one-character PDB chain ID.  It is also used
            as the final FASTA record ID; ``H`` and ``L`` are reserved.
        full_to_local_residue_map: Mapping from full-length antigen coordinates
            to local one-based indices in ``pdb_antigen_sequence``.  Every local
            residue must be represented exactly once.
        profile: ``smoke`` or ``full``.
        iggm_repo_dir: Working directory of the pinned IgGM checkout.  The
            adapter records but does not inspect or mutate this checkout.
        python_executable: Executable token placed at argv[0].
        input_dir: Planned location of generated FASTA inputs.
        output_dir: Planned root of IgGM result directories.

    Returns:
        A JSON-serializable handoff containing one job for every sorted
        template-by-epitope combination, including FASTA content, argv, expected
        output paths, validation contracts, and pinned provenance.

    Raises:
        IgGMAdapterError: On any schema, sequence, PDB-chain, CDR-mask, epitope,
            alignment, or coordinate-mapping ambiguity.
    """

    if profile not in PROFILE_PARAMETERS:
        raise IgGMAdapterError(
            f"Unsupported profile {profile!r}; expected one of {sorted(PROFILE_PARAMETERS)}"
        )
    if not isinstance(request, Mapping):
        raise IgGMAdapterError("Normalized request must be a mapping")

    _validate_request_header(request)
    full_sequence = _validate_request_antigen(request)
    templates = _validate_templates(request.get("templates"))
    epitopes = _validate_epitopes(request.get("epitopes"), full_sequence)

    pdb_path = Path(target_pdb_path).expanduser()
    if not pdb_path.is_file():
        raise IgGMAdapterError(f"Target antigen PDB does not exist or is not a file: {pdb_path}")
    pdb_path = pdb_path.resolve()

    chain_id = _validate_antigen_chain_id(pdb_antigen_chain)
    pdb_sequence = _normalize_sequence(
        pdb_antigen_sequence,
        label="pdb_antigen_sequence",
        allow_mask=False,
    )
    sequence_from_file = _read_pdb_chain_sequence(pdb_path, chain_id)
    if sequence_from_file != pdb_sequence:
        raise IgGMAdapterError(
            "Explicit PDB antigen sequence does not match the residues parsed "
            f"from chain {chain_id!r}: explicit length={len(pdb_sequence)}, "
            f"parsed length={len(sequence_from_file)}"
        )

    coordinate_map = _normalize_coordinate_map(
        full_to_local_residue_map,
        full_sequence=full_sequence,
        pdb_sequence=pdb_sequence,
    )
    mapped_epitopes = [
        _map_epitope(epitope, coordinate_map, pdb_sequence) for epitope in epitopes
    ]

    params = dict(PROFILE_PARAMETERS[profile])
    input_root = Path(input_dir)
    output_root = Path(output_dir)
    working_directory = str(Path(iggm_repo_dir))
    python_token = str(python_executable).strip()
    if not python_token:
        raise IgGMAdapterError("python_executable must not be empty")

    request_digest = _json_digest(request)
    pdb_digest = _file_sha256(pdb_path)
    map_digest = _json_digest({str(key): value for key, value in coordinate_map.items()})
    jobs: list[dict[str, Any]] = []
    used_job_ids: set[str] = set()

    for template in sorted(templates, key=lambda row: row["template_id"]):
        for epitope in sorted(mapped_epitopes, key=lambda row: row["epitope_id"]):
            identity = {
                "template_id": template["template_id"],
                "epitope_id": epitope["epitope_id"],
                "masked_vh": template["masked_vh"],
                "masked_vl": template["masked_vl"],
                "pdb_sequence_sha256": hashlib.sha256(pdb_sequence.encode("ascii")).hexdigest(),
                "local_epitope": epitope["local_positions_1_based"],
                "profile": profile,
            }
            suffix = _json_digest(identity)[:10]
            job_id = (
                f"iggm_{_slug(template['template_id'])}__"
                f"{_slug(epitope['epitope_id'])}__{suffix}"
            )
            if job_id in used_job_ids:
                raise IgGMAdapterError(f"Deterministic IgGM job ID collision: {job_id}")
            used_job_ids.add(job_id)

            fasta_path = input_root / f"{job_id}.fasta"
            result_path = output_root / job_id
            fasta_content = render_iggm_fasta(
                template["masked_vh"],
                template["masked_vl"],
                pdb_sequence,
                antigen_chain=chain_id,
            )
            argv = _build_cli_argv(
                python_executable=python_token,
                fasta_path=fasta_path,
                pdb_path=pdb_path,
                epitope_positions=epitope["local_positions_1_based"],
                result_dir=result_path,
                parameters=params,
            )
            expected_outputs = _expected_outputs(
                result_path=result_path,
                fasta_stem=fasta_path.stem,
                num_samples=int(params["num_samples"]),
            )
            jobs.append(
                {
                    "job_id": job_id,
                    "engine": "IgGM",
                    "execution_state": "planned_not_run",
                    "data_status": "real_model_job_spec_not_executed",
                    "template_id": template["template_id"],
                    "framework_source_id": template.get("framework_source_id", ""),
                    "epitope_id": epitope["epitope_id"],
                    "epitope_full_positions_1_based": epitope["full_positions_1_based"],
                    "epitope_local_positions_1_based": epitope["local_positions_1_based"],
                    "profile": profile,
                    "parameters": dict(params),
                    "input_fasta_path": str(fasta_path),
                    "input_fasta_content": fasta_content,
                    "target_pdb_path": str(pdb_path),
                    "target_pdb_chain": chain_id,
                    "working_directory": working_directory,
                    "command_argv": argv,
                    "command_display": shlex.join(argv),
                    "environment": {},
                    "expected_outputs": expected_outputs,
                    "acceptance_checks": [
                        "all expected PDB and FASTA files exist and are non-empty",
                        "output FASTA contains H, L, and the configured antigen chain",
                        "output H and L contain no X and retain the masked-input chain lengths",
                        "output antigen sequence is identical to pdb_antigen_sequence",
                        "output PDB contains H, L, and the configured antigen chain",
                        "retain stdout/stderr, checkpoint hashes, and GPU/runtime metadata",
                    ],
                }
            )

    return {
        "schema": HANDOFF_SCHEMA,
        "execution_state": "planned_not_run",
        "result_provenance": "adapter_compilation_only_no_iggm_output",
        "request_schema": NORMALIZED_REQUEST_SCHEMA,
        "request_sha256": request_digest,
        "profile": profile,
        "selected_profile_parameters": params,
        "available_profiles": {key: dict(value) for key, value in PROFILE_PARAMETERS.items()},
        "upstream": json.loads(json.dumps(UPSTREAM)),
        "target": {
            "pdb_path": str(pdb_path),
            "pdb_sha256": pdb_digest,
            "antigen_chain": chain_id,
            "antigen_sequence": pdb_sequence,
            "antigen_sequence_sha256": hashlib.sha256(pdb_sequence.encode("ascii")).hexdigest(),
            "full_to_local_residue_map": {
                str(key): value for key, value in coordinate_map.items()
            },
            "full_to_local_residue_map_sha256": map_digest,
        },
        "job_count": len(jobs),
        "jobs": jobs,
        "execution_requirements": {
            "checkout_exact_upstream_commit": UPSTREAM["commit"],
            "checkpoint_directory_behavior": (
                "official code reads/downloads required .pth files under "
                "<IgGM working directory>/checkpoints"
            ),
            "gpu_note": (
                "official environment is CUDA 11.7; upstream publishes no minimum VRAM. "
                "design.py falls back to CPU but CPU practicality is not established"
            ),
            "seed_note": (
                "upstream design.py exposes no seed CLI option and time-shuffles batches; "
                "job specifications are deterministic but model-output byte identity is not promised"
            ),
            "device_note": (
                "upstream parses --device but selects CUDA whenever torch.cuda.is_available(); "
                "select a GPU externally (for example with CUDA_VISIBLE_DEVICES)"
            ),
            "relaxation": "disabled; PyRosetta is optional and outside this smoke/full contract",
        },
    }


def render_iggm_fasta(
    masked_vh: str,
    masked_vl: str,
    antigen_sequence: str,
    *,
    antigen_chain: str,
) -> str:
    """Render the official IgGM H/L/antigen FASTA representation."""

    vh = _normalize_sequence(masked_vh, label="masked_vh", allow_mask=True)
    vl = _normalize_sequence(masked_vl, label="masked_vl", allow_mask=True)
    antigen = _normalize_sequence(
        antigen_sequence,
        label="antigen_sequence",
        allow_mask=False,
    )
    chain_id = _validate_antigen_chain_id(antigen_chain)
    if "X" not in vh or "X" not in vl:
        raise IgGMAdapterError("Both H and L must contain X design masks")
    return f">H\n{vh}\n>L\n{vl}\n>{chain_id}\n{antigen}\n"


def _validate_request_header(request: Mapping[str, Any]) -> None:
    schema = request.get("schema")
    if schema != NORMALIZED_REQUEST_SCHEMA:
        raise IgGMAdapterError(
            f"Unsupported request schema {schema!r}; expected {NORMALIZED_REQUEST_SCHEMA!r}"
        )
    engine = request.get("engine")
    if engine not in (None, "IgGM"):
        raise IgGMAdapterError(f"Request engine must be 'IgGM', not {engine!r}")
    mode = request.get("campaign_mode")
    if mode != "paired_Fv_six_CDR_de_novo_design":
        raise IgGMAdapterError(
            "IgGM adapter requires campaign_mode='paired_Fv_six_CDR_de_novo_design'"
        )


def _validate_request_antigen(request: Mapping[str, Any]) -> str:
    antigen = request.get("antigen")
    if not isinstance(antigen, Mapping):
        raise IgGMAdapterError("request.antigen must be an object")
    return _normalize_sequence(
        antigen.get("full_sequence"),
        label="request.antigen.full_sequence",
        allow_mask=False,
    )


def _validate_templates(value: Any) -> list[dict[str, Any]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise IgGMAdapterError("request.templates must be a non-empty list")
    rows: list[dict[str, Any]] = []
    template_ids: set[str] = set()
    for index, raw in enumerate(value):
        if not isinstance(raw, Mapping):
            raise IgGMAdapterError(f"templates[{index}] must be an object")
        template_id = _required_identifier(raw.get("template_id"), f"templates[{index}].template_id")
        if template_id in template_ids:
            raise IgGMAdapterError(f"Duplicate template_id: {template_id}")
        template_ids.add(template_id)
        if raw.get("template_role") != "framework_source_only":
            raise IgGMAdapterError(
                f"Template {template_id} is not marked framework_source_only"
            )
        vh = _normalize_sequence(
            raw.get("masked_vh"),
            label=f"{template_id}.masked_vh",
            allow_mask=True,
        )
        vl = _normalize_sequence(
            raw.get("masked_vl"),
            label=f"{template_id}.masked_vl",
            allow_mask=True,
        )
        regions = _validate_six_cdr_regions(template_id, raw.get("design_regions"), vh, vl)
        declared = raw.get("designed_regions")
        if declared is not None:
            declared_regions = tuple(item for item in str(declared).split(";") if item)
            if set(declared_regions) != set(DESIGN_REGIONS) or len(declared_regions) != 6:
                raise IgGMAdapterError(
                    f"Template {template_id} designed_regions must name all six CDRs exactly once"
                )
        rows.append(
            {
                "template_id": template_id,
                "framework_source_id": str(raw.get("framework_source_id", "")),
                "masked_vh": vh,
                "masked_vl": vl,
                "design_regions": regions,
            }
        )
    return rows


def _validate_six_cdr_regions(
    template_id: str,
    value: Any,
    vh: str,
    vl: str,
) -> list[dict[str, Any]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise IgGMAdapterError(f"Template {template_id} design_regions must be a list")
    if len(value) != 6:
        raise IgGMAdapterError(f"Template {template_id} must declare exactly six CDR regions")

    seen_names: set[str] = set()
    masked_positions: dict[str, set[int]] = {"VH": set(), "VL": set()}
    normalized: list[dict[str, Any]] = []
    for index, raw in enumerate(value):
        if not isinstance(raw, Mapping):
            raise IgGMAdapterError(
                f"Template {template_id} design_regions[{index}] must be an object"
            )
        name = str(raw.get("region", ""))
        if name not in DESIGN_REGIONS or name in seen_names:
            raise IgGMAdapterError(
                f"Template {template_id} has invalid or duplicate design region {name!r}"
            )
        seen_names.add(name)
        expected_chain = "VH" if name.startswith("H") else "VL"
        chain = str(raw.get("chain", ""))
        if chain != expected_chain:
            raise IgGMAdapterError(
                f"Template {template_id} region {name} must use chain {expected_chain}"
            )
        sequence = vh if chain == "VH" else vl
        start = _positive_int(raw.get("start_1_based"), f"{template_id}.{name}.start")
        end = _positive_int(
            raw.get("end_1_based_inclusive"),
            f"{template_id}.{name}.end",
        )
        if start > end or end > len(sequence):
            raise IgGMAdapterError(
                f"Template {template_id} region {name} is outside its {chain} sequence"
            )
        positions = set(range(start, end + 1))
        if positions & masked_positions[chain]:
            raise IgGMAdapterError(
                f"Template {template_id} has overlapping CDR masks on {chain}"
            )
        masked_positions[chain].update(positions)
        if any(sequence[position - 1] != "X" for position in positions):
            raise IgGMAdapterError(
                f"Template {template_id} region {name} is not completely X-masked"
            )
        declared_length = raw.get("length_aa")
        if declared_length is not None and _positive_int(
            declared_length, f"{template_id}.{name}.length_aa"
        ) != len(positions):
            raise IgGMAdapterError(
                f"Template {template_id} region {name} length_aa disagrees with its coordinates"
            )
        normalized.append(
            {
                "region": name,
                "chain": chain,
                "start_1_based": start,
                "end_1_based_inclusive": end,
                "length_aa": len(positions),
            }
        )

    if seen_names != set(DESIGN_REGIONS):
        missing = sorted(set(DESIGN_REGIONS) - seen_names)
        raise IgGMAdapterError(f"Template {template_id} lacks CDR masks: {missing}")
    for chain, sequence in (("VH", vh), ("VL", vl)):
        actual_x = {index for index, residue in enumerate(sequence, start=1) if residue == "X"}
        if actual_x != masked_positions[chain]:
            raise IgGMAdapterError(
                f"Template {template_id} {chain} contains X outside declared CDRs "
                "or declares unmasked CDR positions"
            )
    normalized.sort(key=lambda row: DESIGN_REGIONS.index(row["region"]))
    return normalized


def _validate_epitopes(value: Any, full_sequence: str) -> list[dict[str, Any]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise IgGMAdapterError("request.epitopes must be a non-empty list")
    rows: list[dict[str, Any]] = []
    epitope_ids: set[str] = set()
    for index, raw in enumerate(value):
        if not isinstance(raw, Mapping):
            raise IgGMAdapterError(f"epitopes[{index}] must be an object")
        epitope_id = _required_identifier(raw.get("epitope_id"), f"epitopes[{index}].epitope_id")
        if epitope_id in epitope_ids:
            raise IgGMAdapterError(f"Duplicate epitope_id: {epitope_id}")
        epitope_ids.add(epitope_id)
        start = _positive_int(raw.get("start_1_based"), f"{epitope_id}.start_1_based")
        end = _positive_int(
            raw.get("end_1_based_inclusive"),
            f"{epitope_id}.end_1_based_inclusive",
        )
        if start > end or end > len(full_sequence):
            raise IgGMAdapterError(f"Epitope {epitope_id} is outside the full antigen sequence")
        sequence = _normalize_sequence(
            raw.get("sequence"),
            label=f"{epitope_id}.sequence",
            allow_mask=False,
        )
        expected_sequence = full_sequence[start - 1 : end]
        if sequence != expected_sequence:
            raise IgGMAdapterError(
                f"Epitope {epitope_id} sequence does not match full antigen coordinates {start}-{end}"
            )
        hotspots = raw.get("candidate_hotspot_residue_indices")
        if not isinstance(hotspots, Sequence) or isinstance(hotspots, (str, bytes)):
            raise IgGMAdapterError(
                f"Epitope {epitope_id} must provide candidate_hotspot_residue_indices"
            )
        positions = [
            _positive_int(position, f"{epitope_id}.candidate_hotspot_residue_indices")
            for position in hotspots
        ]
        expected_positions = list(range(start, end + 1))
        if positions != expected_positions:
            raise IgGMAdapterError(
                f"Epitope {epitope_id} hotspots must exactly and monotonically cover {start}-{end}"
            )
        rows.append(
            {
                "epitope_id": epitope_id,
                "sequence": sequence,
                "full_positions_1_based": positions,
            }
        )
    return rows


def _normalize_coordinate_map(
    value: Mapping[int | str, int | str],
    *,
    full_sequence: str,
    pdb_sequence: str,
) -> dict[int, int]:
    if not isinstance(value, Mapping) or not value:
        raise IgGMAdapterError("full_to_local_residue_map must be a non-empty mapping")
    normalized: dict[int, int] = {}
    for raw_full, raw_local in value.items():
        full_position = _positive_int(raw_full, "full-coordinate map key")
        local_position = _positive_int(raw_local, "full-coordinate map value")
        if full_position in normalized:
            raise IgGMAdapterError(
                f"Duplicate full coordinate after integer normalization: {full_position}"
            )
        normalized[full_position] = local_position

    local_values = list(normalized.values())
    expected_local_values = list(range(1, len(pdb_sequence) + 1))
    if sorted(local_values) != expected_local_values:
        raise IgGMAdapterError(
            "Coordinate map must cover every local PDB sequence position exactly once "
            f"(expected 1..{len(pdb_sequence)})"
        )
    by_local = sorted(normalized.items(), key=lambda item: item[1])
    full_positions_in_pdb_order = [item[0] for item in by_local]
    if full_positions_in_pdb_order != sorted(full_positions_in_pdb_order):
        raise IgGMAdapterError("Coordinate map reverses or reorders the antigen sequence")

    for full_position, local_position in normalized.items():
        if full_position > len(full_sequence):
            raise IgGMAdapterError(
                f"Full coordinate {full_position} exceeds antigen length {len(full_sequence)}"
            )
        full_residue = full_sequence[full_position - 1]
        local_residue = pdb_sequence[local_position - 1]
        if full_residue != local_residue:
            raise IgGMAdapterError(
                "Coordinate-map sequence mismatch: "
                f"full[{full_position}]={full_residue}, local[{local_position}]={local_residue}"
            )
    return dict(sorted(normalized.items()))


def _map_epitope(
    epitope: Mapping[str, Any],
    coordinate_map: Mapping[int, int],
    pdb_sequence: str,
) -> dict[str, Any]:
    full_positions = list(epitope["full_positions_1_based"])
    missing = [position for position in full_positions if position not in coordinate_map]
    if missing:
        raise IgGMAdapterError(
            f"Epitope {epitope['epitope_id']} is absent from the PDB coordinate map at {missing}"
        )
    local_positions = [coordinate_map[position] for position in full_positions]
    if local_positions != sorted(local_positions) or len(set(local_positions)) != len(local_positions):
        raise IgGMAdapterError(
            f"Epitope {epitope['epitope_id']} maps non-monotonically or ambiguously into the PDB"
        )
    local_sequence = "".join(pdb_sequence[position - 1] for position in local_positions)
    if local_sequence != epitope["sequence"]:
        raise IgGMAdapterError(
            f"Epitope {epitope['epitope_id']} does not align to the explicit PDB sequence"
        )
    return {
        "epitope_id": epitope["epitope_id"],
        "sequence": epitope["sequence"],
        "full_positions_1_based": full_positions,
        "local_positions_1_based": local_positions,
    }


def _build_cli_argv(
    *,
    python_executable: str,
    fasta_path: Path,
    pdb_path: Path,
    epitope_positions: Sequence[int],
    result_dir: Path,
    parameters: Mapping[str, Any],
) -> list[str]:
    argv = [
        python_executable,
        "design.py",
        "--fasta",
        str(fasta_path),
        "--antigen",
        str(pdb_path),
        "--epitope",
        *(str(position) for position in epitope_positions),
        "--run_task",
        str(parameters["run_task"]),
        "--num_samples",
        str(parameters["num_samples"]),
        "--steps",
        str(parameters["steps"]),
        "--chunk_size",
        str(parameters["chunk_size"]),
        "--temperature",
        str(parameters["temperature"]),
        "--max_antigen_size",
        str(parameters["max_antigen_size"]),
        "--output",
        str(result_dir),
    ]
    if parameters.get("relax"):
        argv.append("--relax")
    return argv


def _expected_outputs(
    *,
    result_path: Path,
    fasta_stem: str,
    num_samples: int,
) -> list[dict[str, str | int]]:
    outputs: list[dict[str, str | int]] = []
    for sample_index in range(num_samples):
        for file_type, suffix in (("structure", ".pdb"), ("sequence", ".fasta")):
            outputs.append(
                {
                    "sample_index": sample_index,
                    "type": file_type,
                    "path": str(result_path / f"{fasta_stem}_{sample_index}{suffix}"),
                }
            )
    return outputs


_PDB_RESIDUE_TO_ONE = {
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
    "MSE": "M",
}


def _read_pdb_chain_sequence(path: Path, chain_id: str) -> str:
    residues: list[str] = []
    seen_residue_keys: set[tuple[str, str]] = set()
    saw_model = False
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise IgGMAdapterError(f"Cannot read target PDB: {path}") from exc
    for line in lines:
        record = line[:6].strip().upper()
        if record == "MODEL":
            if saw_model:
                break
            saw_model = True
            continue
        if record == "ENDMDL" and saw_model:
            break
        if record not in {"ATOM", "HETATM"} or len(line) < 27:
            continue
        if line[21:22] != chain_id:
            continue
        residue_name = line[17:20].strip().upper()
        if record == "HETATM" and residue_name not in _PDB_RESIDUE_TO_ONE:
            continue
        residue_key = (line[22:26].strip(), line[26:27].strip())
        if residue_key in seen_residue_keys:
            continue
        if residue_name not in _PDB_RESIDUE_TO_ONE:
            raise IgGMAdapterError(
                f"Unsupported residue {residue_name!r} in PDB chain {chain_id} at {residue_key}"
            )
        seen_residue_keys.add(residue_key)
        residues.append(_PDB_RESIDUE_TO_ONE[residue_name])
    if not residues:
        raise IgGMAdapterError(f"PDB contains no parseable residues for antigen chain {chain_id!r}")
    return "".join(residues)


def _validate_antigen_chain_id(value: Any) -> str:
    chain_id = str(value)
    if len(chain_id) != 1 or chain_id.isspace() or not chain_id.isalnum():
        raise IgGMAdapterError("pdb_antigen_chain must be one alphanumeric PDB chain ID")
    if chain_id in {"H", "L"}:
        raise IgGMAdapterError("pdb_antigen_chain cannot be H or L (reserved for antibody chains)")
    return chain_id


def _normalize_sequence(value: Any, *, label: str, allow_mask: bool) -> str:
    if not isinstance(value, str):
        raise IgGMAdapterError(f"{label} must be a string")
    sequence = re.sub(r"\s+", "", value).upper()
    if not sequence:
        raise IgGMAdapterError(f"{label} must not be empty")
    allowed = set(CANONICAL_AMINO_ACIDS)
    if allow_mask:
        allowed.add("X")
    invalid = sorted(set(sequence) - allowed)
    if invalid:
        raise IgGMAdapterError(f"{label} contains unsupported residues: {invalid}")
    return sequence


def _positive_int(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise IgGMAdapterError(f"{label} must be a positive integer")
    if isinstance(value, int):
        result = value
    elif isinstance(value, str) and re.fullmatch(r"[1-9]\d*", value.strip()):
        result = int(value.strip())
    else:
        raise IgGMAdapterError(f"{label} must be a positive integer")
    if result <= 0:
        raise IgGMAdapterError(f"{label} must be a positive integer")
    return result


def _required_identifier(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise IgGMAdapterError(f"{label} must be a non-empty string")
    return value.strip()


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "-", value).strip("-").lower()
    return slug[:64] or "unnamed"


def _json_digest(value: Any) -> str:
    try:
        encoded = json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise IgGMAdapterError("Request contains values that are not JSON-serializable") from exc
    return hashlib.sha256(encoded).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as exc:
        raise IgGMAdapterError(f"Cannot hash target PDB: {path}") from exc
    return digest.hexdigest()


__all__ = [
    "DESIGN_REGIONS",
    "HANDOFF_SCHEMA",
    "IgGMAdapterError",
    "NORMALIZED_REQUEST_SCHEMA",
    "PROFILE_PARAMETERS",
    "UPSTREAM",
    "build_iggm_jobs",
    "load_normalized_request",
    "render_iggm_fasta",
]
