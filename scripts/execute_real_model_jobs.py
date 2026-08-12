#!/usr/bin/env python3
"""Fail-closed executor for a reviewed unified real-model handoff.

The compiler in :mod:`prepare_real_model_jobs` deliberately emits more native
jobs than a smoke run is allowed to submit.  This executor therefore consumes
*only* ``engines[].execution_jobs`` from ``unified_handoff_manifest.json``.  It
never follows ``native_plan`` and never falls back to ``engines[].jobs``.

Execution is opt-in.  Without ``--execute`` the command validates the complete
handoff/runtime-attestation contract and prints a dry-run preview without
launching a process or staging a file.  With ``--execute`` commands are run one
at a time with ``shell=False`` while an exclusive handoff lock is held.  A
hash-bound execution report makes resuming already-successful jobs safe.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from time import monotonic
from typing import Any, Iterator, Mapping, Sequence


UNIFIED_HANDOFF_SCHEMA = "nfl_ab_design.real_model_handoff.v1"
RUNTIME_ATTESTATION_SCHEMA = "nfl_ab_design.runtime_attestation.v1"
EXECUTION_REPORT_SCHEMA = "nfl_ab_design.real_model_execution_report.v1"
DRY_RUN_SCHEMA = "nfl_ab_design.real_model_execution_preview.v1"
EXPECTED_ENGINES = ("RFantibody", "IgGM", "Germinal")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
SAFE_ID_PATTERN = re.compile(r"[A-Za-z0-9._-]+")
PLACEHOLDER_PATTERN = re.compile(r"<[^<>]+>")


class RealModelExecutorError(ValueError):
    """Raised before or during execution when a safety contract is violated."""


class RealModelExecutionFailed(RealModelExecutorError):
    """Raised after recording a failed command and stopping the serial run."""


def _now() -> str:
    return datetime.now().astimezone().isoformat()


def _json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_object(path: Path, *, label: str) -> dict[str, Any]:
    if not path.is_file():
        raise RealModelExecutorError(f"{label} is not a file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise RealModelExecutorError(f"Cannot read {label} {path}: {exc}") from exc
    except json.JSONDecodeError as exc:
        raise RealModelExecutorError(f"Invalid JSON in {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RealModelExecutorError(f"{label} root must be a JSON object: {path}")
    return value


def _mapping(value: Any, *, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise RealModelExecutorError(f"{label} must be a JSON object")
    return value


def _string(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise RealModelExecutorError(f"{label} must be a non-empty string")
    text = value.strip()
    if PLACEHOLDER_PATTERN.search(text):
        raise RealModelExecutorError(f"{label} contains an unresolved placeholder: {text!r}")
    return text


def _sha256(value: Any, *, label: str) -> str:
    text = _string(value, label=label)
    if not SHA256_PATTERN.fullmatch(text):
        raise RealModelExecutorError(f"{label} must be a lowercase SHA-256 hex digest")
    return text


def _safe_id(value: Any, *, label: str) -> str:
    text = _string(value, label=label)
    if text in {".", ".."} or not SAFE_ID_PATTERN.fullmatch(text):
        raise RealModelExecutorError(
            f"{label} must contain only ASCII letters, digits, '.', '_', or '-': {text!r}"
        )
    return text


def _path(value: Any, *, label: str) -> Path:
    text = _string(value, label=label)
    return Path(text).expanduser().resolve()


def _under(path: Path, root: Path, *, label: str) -> Path:
    resolved = path.expanduser().resolve()
    root_resolved = root.expanduser().resolve()
    try:
        resolved.relative_to(root_resolved)
    except ValueError as exc:
        raise RealModelExecutorError(
            f"{label} escapes its authorized root {root_resolved}: {resolved}"
        ) from exc
    return resolved


def _atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, raw_temp_path = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temp_path = Path(raw_temp_path)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temp_path.replace(path)
    except BaseException:
        try:
            temp_path.unlink()
        except FileNotFoundError:
            pass
        raise


@contextmanager
def _exclusive_execution_lock(lock_path: Path) -> Iterator[None]:
    """Prevent two executors from sharing the same single-GPU handoff."""

    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RealModelExecutorError(
                f"Another executor already holds the single-GPU lock: {lock_path}"
            ) from exc
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _validate_input_artifacts(job: Mapping[str, Any], *, label: str) -> list[dict[str, str]]:
    raw_artifacts = job.get("input_artifacts")
    if not isinstance(raw_artifacts, list) or not raw_artifacts:
        raise RealModelExecutorError(f"{label}.input_artifacts must be a non-empty list")
    artifacts: list[dict[str, str]] = []
    roles: set[str] = set()
    for index, raw_artifact in enumerate(raw_artifacts):
        artifact = _mapping(raw_artifact, label=f"{label}.input_artifacts[{index}]")
        role = _string(artifact.get("role"), label=f"{label}.input_artifacts[{index}].role")
        if role in roles:
            raise RealModelExecutorError(f"{label} repeats input artifact role {role!r}")
        roles.add(role)
        path = _path(artifact.get("path"), label=f"{label}.{role}.path")
        expected_hash = _sha256(artifact.get("sha256"), label=f"{label}.{role}.sha256")
        if not path.is_file():
            raise RealModelExecutorError(f"{label} input artifact is not a file: {path}")
        actual_hash = _file_sha256(path)
        if actual_hash != expected_hash:
            raise RealModelExecutorError(
                f"{label} input artifact hash mismatch for {role}: "
                f"expected {expected_hash}, observed {actual_hash}: {path}"
            )
        artifacts.append({"role": role, "path": str(path), "sha256": actual_hash})
    return artifacts


def _validate_command(
    raw_command: Any,
    *,
    label: str,
    handoff_root: Path,
) -> dict[str, Any]:
    command = _mapping(raw_command, label=label)
    stage = _safe_id(command.get("stage"), label=f"{label}.stage")
    raw_argv = command.get("argv")
    if not isinstance(raw_argv, list) or not raw_argv:
        raise RealModelExecutorError(f"{label}.argv must be a non-empty JSON array")
    argv: list[str] = []
    for index, raw_token in enumerate(raw_argv):
        token = _string(raw_token, label=f"{label}.argv[{index}]")
        if "\x00" in token:
            raise RealModelExecutorError(f"{label}.argv[{index}] contains a NUL byte")
        argv.append(token)

    raw_working_directory = command.get("working_directory")
    working_directory: str | None
    if raw_working_directory is None:
        working_directory = None
    else:
        directory = _path(raw_working_directory, label=f"{label}.working_directory")
        if not directory.is_dir():
            raise RealModelExecutorError(
                f"{label}.working_directory is not an existing directory: {directory}"
            )
        working_directory = str(directory)

    raw_outputs = command.get("expected_outputs")
    if not isinstance(raw_outputs, list) or not raw_outputs:
        raise RealModelExecutorError(f"{label}.expected_outputs must be a non-empty list")
    outputs: list[dict[str, str]] = []
    for index, raw_output in enumerate(raw_outputs):
        output = _mapping(raw_output, label=f"{label}.expected_outputs[{index}]")
        path_text = _string(
            output.get("path"), label=f"{label}.expected_outputs[{index}].path"
        )
        # Expected output globs are data, not expanded.  Their non-glob prefix
        # must still remain inside this handoff.
        prefix = re.split(r"[*?[]", path_text, maxsplit=1)[0]
        prefix_path = Path(prefix or path_text).expanduser()
        if not prefix_path.is_absolute():
            prefix_path = handoff_root / prefix_path
        _under(prefix_path, handoff_root, label=f"{label}.expected_outputs[{index}]")
        normalized = {"path": path_text}
        if "name" in output:
            normalized["name"] = _safe_id(
                output.get("name"), label=f"{label}.expected_outputs[{index}].name"
            )
        outputs.append(normalized)
    return {
        "stage": stage,
        "argv": argv,
        "working_directory": working_directory,
        "expected_outputs": outputs,
    }


def _validate_execution_job(
    raw_job: Any,
    *,
    engine: str,
    index: int,
    handoff_root: Path,
) -> dict[str, Any]:
    label = f"{engine}.execution_jobs[{index}]"
    job = _mapping(raw_job, label=label)
    if job.get("engine") != engine:
        raise RealModelExecutorError(
            f"{label}.engine must be {engine!r}, got {job.get('engine')!r}"
        )
    if job.get("selected_for_execution") is not True:
        raise RealModelExecutorError(f"{label} is not explicitly selected_for_execution=true")
    if job.get("execution_disposition") != "selected_for_execution":
        raise RealModelExecutorError(
            f"{label} has unauthorized execution_disposition "
            f"{job.get('execution_disposition')!r}"
        )
    if job.get("execution_state") != "planned_not_executed":
        raise RealModelExecutorError(
            f"{label}.execution_state must be 'planned_not_executed'"
        )
    blockers = job.get("unresolved_blockers")
    if blockers != []:
        raise RealModelExecutorError(f"{label}.unresolved_blockers must be an empty list")
    job_id = _safe_id(job.get("job_id"), label=f"{label}.job_id")
    template_id = _string(job.get("template_id"), label=f"{label}.template_id")
    epitope_id = _string(job.get("epitope_id"), label=f"{label}.epitope_id")
    artifacts = _validate_input_artifacts(job, label=label)
    raw_commands = job.get("commands")
    if not isinstance(raw_commands, list) or not raw_commands:
        raise RealModelExecutorError(f"{label}.commands must be a non-empty list")
    commands = [
        _validate_command(
            command,
            label=f"{label}.commands[{command_index}]",
            handoff_root=handoff_root,
        )
        for command_index, command in enumerate(raw_commands)
    ]
    normalized = {
        "job_id": job_id,
        "engine": engine,
        "profile": _string(job.get("profile"), label=f"{label}.profile"),
        "geometry": _string(job.get("geometry"), label=f"{label}.geometry"),
        "template_id": template_id,
        "epitope_id": epitope_id,
        "selected_for_execution": True,
        "execution_disposition": "selected_for_execution",
        "input_artifacts": artifacts,
        "commands": commands,
    }
    normalized["job_sha256"] = _json_sha256(normalized)
    return normalized


def _validate_handoff(
    manifest: Mapping[str, Any],
    *,
    manifest_path: Path,
) -> dict[str, Any]:
    if manifest.get("schema") != UNIFIED_HANDOFF_SCHEMA:
        raise RealModelExecutorError(
            "Executor accepts only unified_handoff_manifest.json with schema "
            f"{UNIFIED_HANDOFF_SCHEMA!r}; native plans/jobs are not executable inputs"
        )
    handoff_identity = _mapping(manifest.get("handoff_identity"), label="handoff_identity")
    identity_sha256 = _sha256(manifest.get("identity_sha256"), label="identity_sha256")
    observed_identity_sha256 = _json_sha256(handoff_identity)
    if observed_identity_sha256 != identity_sha256:
        raise RealModelExecutorError(
            "handoff identity hash mismatch: the manifest or identity was modified"
        )
    handoff_id = _string(manifest.get("handoff_id"), label="handoff_id")
    expected_handoff_id = f"nfl_handoff_{identity_sha256[:24]}"
    if handoff_id != expected_handoff_id:
        raise RealModelExecutorError(
            f"handoff_id does not match identity_sha256: expected {expected_handoff_id!r}"
        )
    if manifest.get("does_not_execute_external_models") is not True:
        raise RealModelExecutorError(
            "handoff must be an unexecuted compiler artifact before executor submission"
        )
    if manifest.get("execution_state") != "planned_not_executed":
        raise RealModelExecutorError(
            "unified handoff execution_state must be 'planned_not_executed'"
        )
    source_integrity = _mapping(
        _mapping(manifest.get("source_run"), label="source_run").get("source_integrity"),
        label="source_run.source_integrity",
    )
    if source_integrity.get("ready_for_execution") is not True:
        raise RealModelExecutorError(
            "source request integrity is not ready_for_execution; runtime attestation "
            "cannot override an unauthenticated request set"
        )

    location = _mapping(manifest.get("handoff_location"), label="handoff_location")
    handoff_root = _path(location.get("path"), label="handoff_location.path")
    if handoff_root != manifest_path.parent.resolve():
        raise RealModelExecutorError(
            "handoff_location.path must exactly match the directory containing the "
            f"unified manifest: {handoff_root} != {manifest_path.parent.resolve()}"
        )

    raw_engines = manifest.get("engines")
    if not isinstance(raw_engines, list) or len(raw_engines) != len(EXPECTED_ENGINES):
        raise RealModelExecutorError("handoff engines must contain exactly three entries")
    engines: list[dict[str, Any]] = []
    observed_engine_names: list[str] = []
    all_job_ids: set[str] = set()
    for engine_index, raw_engine in enumerate(raw_engines):
        engine_manifest = _mapping(raw_engine, label=f"engines[{engine_index}]")
        engine = _string(engine_manifest.get("engine"), label=f"engines[{engine_index}].engine")
        observed_engine_names.append(engine)
        if engine not in EXPECTED_ENGINES:
            raise RealModelExecutorError(f"Unsupported engine in handoff: {engine!r}")
        if engine_manifest.get("execution_state") != "planned_not_executed":
            raise RealModelExecutorError(
                f"{engine}.execution_state must be 'planned_not_executed'"
            )
        # Deliberately no fallback to engine_manifest['jobs'] or native_plan.
        raw_execution_jobs = engine_manifest.get("execution_jobs")
        if not isinstance(raw_execution_jobs, list) or not raw_execution_jobs:
            raise RealModelExecutorError(
                f"{engine} has no non-empty execution_jobs selection; native_plan/jobs "
                "are rejected as submission sources"
            )
        jobs = [
            _validate_execution_job(
                raw_job,
                engine=engine,
                index=job_index,
                handoff_root=handoff_root,
            )
            for job_index, raw_job in enumerate(raw_execution_jobs)
        ]
        selected_ids = engine_manifest.get("selected_job_ids")
        if selected_ids != [job["job_id"] for job in jobs]:
            raise RealModelExecutorError(
                f"{engine}.selected_job_ids must exactly match execution_jobs order"
            )
        if engine_manifest.get("selected_job_count") != len(jobs):
            raise RealModelExecutorError(
                f"{engine}.selected_job_count does not match execution_jobs"
            )
        for job in jobs:
            if job["job_id"] in all_job_ids:
                raise RealModelExecutorError(
                    f"job_id must be globally unique across engines: {job['job_id']}"
                )
            all_job_ids.add(job["job_id"])
        engines.append(
            {
                "engine": engine,
                "required_upstream_revision": _string(
                    engine_manifest.get("required_upstream_revision"),
                    label=f"{engine}.required_upstream_revision",
                ),
                "manifest_ready_for_execution": engine_manifest.get(
                    "ready_for_execution"
                )
                is True,
                "execution_jobs": jobs,
                "execution_jobs_sha256": _json_sha256(raw_execution_jobs),
            }
        )
    if tuple(observed_engine_names) != EXPECTED_ENGINES:
        raise RealModelExecutorError(
            f"engines must be ordered exactly as {list(EXPECTED_ENGINES)}, "
            f"got {observed_engine_names}"
        )
    return {
        "handoff_id": handoff_id,
        "identity_sha256": identity_sha256,
        "manifest_sha256": _file_sha256(manifest_path),
        "handoff_root": handoff_root,
        "engines": engines,
        "execution_plan_sha256": _json_sha256(
            [
                {
                    "engine": engine["engine"],
                    "execution_jobs_sha256": engine["execution_jobs_sha256"],
                    "jobs": engine["execution_jobs"],
                }
                for engine in engines
            ]
        ),
    }


def _validate_attestation(
    attestation: Mapping[str, Any],
    *,
    attestation_path: Path,
    handoff: Mapping[str, Any],
) -> dict[str, Any]:
    if attestation.get("schema") != RUNTIME_ATTESTATION_SCHEMA:
        raise RealModelExecutorError(
            f"runtime attestation schema must be {RUNTIME_ATTESTATION_SCHEMA!r}"
        )
    if attestation.get("handoff_id") != handoff["handoff_id"]:
        raise RealModelExecutorError("runtime attestation handoff_id mismatch")
    if attestation.get("identity_sha256") != handoff["identity_sha256"]:
        raise RealModelExecutorError("runtime attestation identity_sha256 mismatch")
    if attestation.get("handoff_manifest_sha256") != handoff["manifest_sha256"]:
        raise RealModelExecutorError("runtime attestation handoff_manifest_sha256 mismatch")
    raw_engines = attestation.get("engines")
    if not isinstance(raw_engines, list) or len(raw_engines) != len(EXPECTED_ENGINES):
        raise RealModelExecutorError(
            "runtime attestation must contain exactly one entry for each of the three engines"
        )
    attested_by_name: dict[str, Mapping[str, Any]] = {}
    for index, raw_engine in enumerate(raw_engines):
        value = _mapping(raw_engine, label=f"attestation.engines[{index}]")
        engine = _string(value.get("engine"), label=f"attestation.engines[{index}].engine")
        if engine in attested_by_name:
            raise RealModelExecutorError(f"runtime attestation repeats engine {engine!r}")
        attested_by_name[engine] = value
    if set(attested_by_name) != set(EXPECTED_ENGINES):
        raise RealModelExecutorError(
            f"runtime attestation engines must be exactly {list(EXPECTED_ENGINES)}"
        )

    normalized_engines: list[dict[str, Any]] = []
    for engine_handoff in handoff["engines"]:
        engine = engine_handoff["engine"]
        value = attested_by_name[engine]
        revision = _string(value.get("revision"), label=f"{engine} attested revision")
        if revision != engine_handoff["required_upstream_revision"]:
            raise RealModelExecutorError(
                f"{engine} attested revision {revision!r} does not match handoff "
                f"requirement {engine_handoff['required_upstream_revision']!r}"
            )
        jobs_hash = _sha256(
            value.get("execution_jobs_sha256"),
            label=f"{engine} attested execution_jobs_sha256",
        )
        if jobs_hash != engine_handoff["execution_jobs_sha256"]:
            raise RealModelExecutorError(
                f"{engine} attested execution_jobs_sha256 does not match handoff"
            )
        checkpoints = _mapping(
            value.get("checkpoint_sha256"), label=f"{engine} checkpoint_sha256"
        )
        if not checkpoints:
            raise RealModelExecutorError(
                f"{engine} checkpoint_sha256 must contain at least one verified checkpoint"
            )
        normalized_checkpoints: dict[str, str] = {}
        for raw_name, raw_hash in checkpoints.items():
            name = _safe_id(raw_name, label=f"{engine} checkpoint logical name")
            normalized_checkpoints[name] = _sha256(
                raw_hash, label=f"{engine} checkpoint_sha256[{name!r}]"
            )
        if value.get("ready") is not True:
            raise RealModelExecutorError(f"{engine} runtime attestation ready must be true")
        if (
            not engine_handoff["manifest_ready_for_execution"]
            and value.get("overrides_manifest_ready_for_execution") is not True
        ):
            raise RealModelExecutorError(
                f"{engine} handoff ready_for_execution=false; a reviewed attestation "
                "must explicitly set overrides_manifest_ready_for_execution=true"
            )
        normalized_engines.append(
            {
                "engine": engine,
                "revision": revision,
                "execution_jobs_sha256": jobs_hash,
                "checkpoint_sha256": dict(sorted(normalized_checkpoints.items())),
                "ready": True,
                "overrides_manifest_ready_for_execution": bool(
                    value.get("overrides_manifest_ready_for_execution")
                ),
            }
        )
    return {
        "path": str(attestation_path),
        "sha256": _file_sha256(attestation_path),
        "schema": RUNTIME_ATTESTATION_SCHEMA,
        "engines": normalized_engines,
    }


def _artifact_by_role(job: Mapping[str, Any], role: str) -> Mapping[str, str]:
    matching = [item for item in job["input_artifacts"] if item["role"] == role]
    if len(matching) != 1:
        raise RealModelExecutorError(
            f"{job['job_id']} must contain exactly one input artifact with role {role!r}"
        )
    return matching[0]


def _argv_assignment(argv: Sequence[str], key: str, *, job_id: str) -> str:
    prefix = f"{key}="
    matching = [token[len(prefix) :] for token in argv if token.startswith(prefix)]
    if len(matching) != 1 or not matching[0]:
        raise RealModelExecutorError(
            f"Germinal job {job_id} must contain exactly one non-empty {key}= argv token"
        )
    return matching[0]


def _germinal_staging_plan(job: Mapping[str, Any], *, handoff_root: Path) -> list[dict[str, str]]:
    if len(job["commands"]) != 1:
        raise RealModelExecutorError(
            f"Germinal job {job['job_id']} must contain exactly one authorized command"
        )
    command = job["commands"][0]
    if command["working_directory"] is None:
        raise RealModelExecutorError(
            f"Germinal job {job['job_id']} requires an explicit working_directory"
        )
    germinal_repo = Path(command["working_directory"]).resolve()
    if not (germinal_repo / "run_germinal.py").is_file():
        raise RealModelExecutorError(
            f"Germinal working_directory lacks run_germinal.py: {germinal_repo}"
        )
    argv = command["argv"]
    target_name = _safe_id(
        _argv_assignment(argv, "target", job_id=job["job_id"]),
        label=f"Germinal {job['job_id']} target config name",
    )
    pdb_dir = _path(
        _argv_assignment(argv, "pdb_dir", job_id=job["job_id"]),
        label=f"Germinal {job['job_id']} pdb_dir",
    )
    expected_pdb_dir = (
        handoff_root / "germinal" / "jobs" / job["job_id"] / "pdbs"
    ).resolve()
    if pdb_dir != expected_pdb_dir:
        raise RealModelExecutorError(
            f"Germinal {job['job_id']} pdb_dir is not its isolated manifest workspace: "
            f"{pdb_dir} != {expected_pdb_dir}"
        )
    _under(pdb_dir, handoff_root, label=f"Germinal {job['job_id']} pdb_dir")

    scfv_source = Path(_artifact_by_role(job, "template_scfv_pdb")["path"]).resolve()
    yaml_source = Path(_artifact_by_role(job, "generated_target_yaml")["path"]).resolve()
    yaml_source_root = (handoff_root / "germinal" / "target_configs").resolve()
    _under(
        yaml_source,
        yaml_source_root,
        label=f"Germinal {job['job_id']} generated target YAML",
    )
    if yaml_source.name != f"{target_name}.yaml":
        raise RealModelExecutorError(
            f"Germinal {job['job_id']} target= does not match generated YAML basename"
        )

    installed_target_root = (germinal_repo / "configs" / "target").resolve()
    if not installed_target_root.is_dir():
        raise RealModelExecutorError(
            f"Germinal checkout lacks configs/target directory: {installed_target_root}"
        )
    installed_yaml = (installed_target_root / f"{target_name}.yaml").resolve()
    _under(
        installed_yaml,
        installed_target_root,
        label=f"Germinal {job['job_id']} installed target YAML",
    )
    staged_scfv = (pdb_dir / "scfv.pdb").resolve()
    _under(staged_scfv, pdb_dir, label=f"Germinal {job['job_id']} staged scFv")
    return [
        {
            "role": "template_scfv_pdb",
            "source": str(scfv_source),
            "destination": str(staged_scfv),
            "sha256": _artifact_by_role(job, "template_scfv_pdb")["sha256"],
        },
        {
            "role": "generated_target_yaml",
            "source": str(yaml_source),
            "destination": str(installed_yaml),
            "sha256": _artifact_by_role(job, "generated_target_yaml")["sha256"],
        },
    ]


def _stage_file(action: Mapping[str, str]) -> dict[str, str]:
    source = Path(action["source"])
    destination = Path(action["destination"])
    expected_hash = action["sha256"]
    if _file_sha256(source) != expected_hash:
        raise RealModelExecutorError(
            f"Staging source changed after validation: {source}"
        )
    if destination.exists():
        if not destination.is_file():
            raise RealModelExecutorError(
                f"Staging destination exists but is not a file: {destination}"
            )
        if _file_sha256(destination) != expected_hash:
            raise RealModelExecutorError(
                "Refusing to overwrite a different Germinal staged artifact: "
                f"{destination}"
            )
        status = "reused_verified_existing_file"
    else:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        if _file_sha256(destination) != expected_hash:
            raise RealModelExecutorError(
                f"Staged Germinal artifact failed post-copy hash verification: {destination}"
            )
        status = "copied_and_hash_verified"
    return {**dict(action), "status": status}


def _revalidate_job_inputs(job: Mapping[str, Any]) -> list[dict[str, str]]:
    """Close the gap between initial validation and a later serial job start."""

    verified: list[dict[str, str]] = []
    for artifact in job["input_artifacts"]:
        path = Path(artifact["path"])
        if not path.is_file():
            raise RealModelExecutorError(
                f"{job['job_id']} input disappeared before execution: {path}"
            )
        actual_hash = _file_sha256(path)
        if actual_hash != artifact["sha256"]:
            raise RealModelExecutorError(
                f"{job['job_id']} input changed before execution: {path}"
            )
        verified.append(
            {
                "role": artifact["role"],
                "path": str(path),
                "sha256": actual_hash,
                "status": "reverified_immediately_before_job",
            }
        )
    return verified


def _preview(handoff: Mapping[str, Any], attestation: Mapping[str, Any]) -> dict[str, Any]:
    jobs: list[dict[str, Any]] = []
    for engine in handoff["engines"]:
        for job in engine["execution_jobs"]:
            staging = (
                _germinal_staging_plan(job, handoff_root=handoff["handoff_root"])
                if engine["engine"] == "Germinal"
                else []
            )
            jobs.append(
                {
                    "engine": engine["engine"],
                    "job_id": job["job_id"],
                    "job_sha256": job["job_sha256"],
                    "staging": [{**item, "status": "planned_not_written"} for item in staging],
                    "commands": job["commands"],
                    "status": "planned_not_executed",
                }
            )
    return {
        "schema": DRY_RUN_SCHEMA,
        "handoff_id": handoff["handoff_id"],
        "identity_sha256": handoff["identity_sha256"],
        "execution_plan_sha256": handoff["execution_plan_sha256"],
        "runtime_attestation_sha256": attestation["sha256"],
        "status": "dry_run_validated_no_process_or_staging",
        "single_gpu_policy": "exclusive_lock_and_strictly_serial_commands",
        "jobs": jobs,
    }


def _new_report(handoff: Mapping[str, Any], attestation: Mapping[str, Any]) -> dict[str, Any]:
    timestamp = _now()
    return {
        "schema": EXECUTION_REPORT_SCHEMA,
        "handoff_id": handoff["handoff_id"],
        "identity_sha256": handoff["identity_sha256"],
        "handoff_manifest_sha256": handoff["manifest_sha256"],
        "runtime_attestation": {
            "path": attestation["path"],
            "sha256": attestation["sha256"],
            "schema": attestation["schema"],
            "engines": attestation["engines"],
        },
        "execution_plan_sha256": handoff["execution_plan_sha256"],
        "single_gpu_policy": "exclusive_lock_and_strictly_serial_commands",
        "created_at": timestamp,
        "updated_at": timestamp,
        "status": "running",
        "jobs": [],
    }


def _validate_recorded_success(
    record: Mapping[str, Any],
    *,
    current_job: Mapping[str, Any],
) -> None:
    """Require a complete zero-exit attempt before a resume may skip a job."""

    attempts = record.get("attempts")
    if not isinstance(attempts, list) or not attempts:
        raise RealModelExecutorError(
            f"Cannot resume: succeeded job {current_job['job_id']} has no attempts"
        )
    attempt = _mapping(
        attempts[-1], label=f"succeeded job {current_job['job_id']} final attempt"
    )
    if attempt.get("status") != "succeeded":
        raise RealModelExecutorError(
            f"Cannot resume: succeeded job {current_job['job_id']} lacks a "
            "succeeded final attempt"
        )
    recorded_commands = attempt.get("commands")
    current_commands = current_job["commands"]
    if not isinstance(recorded_commands, list) or len(recorded_commands) != len(
        current_commands
    ):
        raise RealModelExecutorError(
            f"Cannot resume: succeeded job {current_job['job_id']} command count mismatch"
        )
    for index, (raw_recorded, expected) in enumerate(
        zip(recorded_commands, current_commands, strict=True)
    ):
        recorded = _mapping(
            raw_recorded,
            label=f"succeeded job {current_job['job_id']} command[{index}]",
        )
        expected_identity = {
            "stage": expected["stage"],
            "argv": expected["argv"],
            "working_directory": expected["working_directory"],
            "shell": False,
        }
        if any(recorded.get(key) != value for key, value in expected_identity.items()):
            raise RealModelExecutorError(
                f"Cannot resume: succeeded job {current_job['job_id']} command[{index}] "
                "identity differs from the authorized command"
            )
        if recorded.get("status") != "succeeded" or recorded.get("exit_code") != 0:
            raise RealModelExecutorError(
                f"Cannot resume: succeeded job {current_job['job_id']} command[{index}] "
                "does not contain a zero-exit success record"
            )
        if not isinstance(recorded.get("stdout"), str) or not isinstance(
            recorded.get("stderr"), str
        ):
            raise RealModelExecutorError(
                f"Cannot resume: succeeded job {current_job['job_id']} command[{index}] "
                "does not preserve stdout/stderr"
            )


def _resume_report(
    state_path: Path,
    *,
    handoff: Mapping[str, Any],
    attestation: Mapping[str, Any],
) -> dict[str, Any]:
    report = _load_object(state_path, label="execution report")
    expected = {
        "schema": EXECUTION_REPORT_SCHEMA,
        "handoff_id": handoff["handoff_id"],
        "identity_sha256": handoff["identity_sha256"],
        "handoff_manifest_sha256": handoff["manifest_sha256"],
        "execution_plan_sha256": handoff["execution_plan_sha256"],
    }
    for key, value in expected.items():
        if report.get(key) != value:
            raise RealModelExecutorError(
                f"Cannot resume: execution report {key} does not match current handoff"
            )
    recorded_attestation = _mapping(
        report.get("runtime_attestation"), label="execution report runtime_attestation"
    )
    if recorded_attestation.get("sha256") != attestation["sha256"]:
        raise RealModelExecutorError(
            "Cannot resume: runtime attestation file hash differs from the original run"
        )
    raw_jobs = report.get("jobs")
    if not isinstance(raw_jobs, list):
        raise RealModelExecutorError("Cannot resume: execution report jobs must be a list")
    current_jobs = {
        job["job_id"]: job
        for engine in handoff["engines"]
        for job in engine["execution_jobs"]
    }
    seen: set[str] = set()
    for index, raw_record in enumerate(raw_jobs):
        record = _mapping(raw_record, label=f"execution report jobs[{index}]")
        job_id = _string(record.get("job_id"), label=f"execution report jobs[{index}].job_id")
        if job_id in seen or job_id not in current_jobs:
            raise RealModelExecutorError(
                f"Cannot resume: unknown or duplicate execution report job_id {job_id!r}"
            )
        seen.add(job_id)
        if record.get("job_sha256") != current_jobs[job_id]["job_sha256"]:
            raise RealModelExecutorError(
                f"Cannot resume: job hash changed for {job_id}"
            )
        if record.get("status") not in {"running", "failed", "succeeded"}:
            raise RealModelExecutorError(
                f"Cannot resume: invalid prior status for {job_id}: {record.get('status')!r}"
            )
        if not isinstance(record.get("attempts"), list):
            raise RealModelExecutorError(
                f"Cannot resume: attempts for {job_id} must be a list"
            )
        if record.get("status") == "succeeded":
            _validate_recorded_success(record, current_job=current_jobs[job_id])
    report["status"] = "running"
    report["updated_at"] = _now()
    return report


def _run_serial(
    *,
    handoff: Mapping[str, Any],
    attestation: Mapping[str, Any],
    state_path: Path,
    resume: bool,
) -> dict[str, Any]:
    if state_path.exists():
        if not resume:
            raise RealModelExecutorError(
                f"Execution report already exists; use --resume after review: {state_path}"
            )
        report = _resume_report(
            state_path, handoff=handoff, attestation=attestation
        )
    else:
        report = _new_report(handoff, attestation)
    records = {record["job_id"]: record for record in report["jobs"]}
    _atomic_write_json(state_path, report)

    for engine in handoff["engines"]:
        for job in engine["execution_jobs"]:
            existing = records.get(job["job_id"])
            if existing is not None and existing["status"] == "succeeded":
                existing["resume_disposition"] = "skipped_previously_succeeded"
                report["updated_at"] = _now()
                _atomic_write_json(state_path, report)
                continue
            if existing is None:
                record: dict[str, Any] = {
                    "engine": engine["engine"],
                    "job_id": job["job_id"],
                    "job_sha256": job["job_sha256"],
                    "status": "running",
                    "attempts": [],
                }
                report["jobs"].append(record)
                records[job["job_id"]] = record
            else:
                record = existing
                record["status"] = "running"
                record.pop("resume_disposition", None)
            attempt: dict[str, Any] = {
                "attempt": len(record["attempts"]) + 1,
                "started_at": _now(),
                "status": "running",
                "input_artifacts": [],
                "staging": [],
                "commands": [],
            }
            record["attempts"].append(attempt)
            report["updated_at"] = _now()
            _atomic_write_json(state_path, report)
            try:
                attempt["input_artifacts"] = _revalidate_job_inputs(job)
                report["updated_at"] = _now()
                _atomic_write_json(state_path, report)
                if engine["engine"] == "Germinal":
                    actions = _germinal_staging_plan(
                        job, handoff_root=handoff["handoff_root"]
                    )
                    attempt["staging"] = [_stage_file(action) for action in actions]
                    report["updated_at"] = _now()
                    _atomic_write_json(state_path, report)
                for command in job["commands"]:
                    started_at = _now()
                    timer = monotonic()
                    command_result: dict[str, Any] = {
                        "stage": command["stage"],
                        "argv": list(command["argv"]),
                        "working_directory": command["working_directory"],
                        "shell": False,
                        "started_at": started_at,
                        "status": "running",
                    }
                    attempt["commands"].append(command_result)
                    _atomic_write_json(state_path, report)
                    try:
                        completed = subprocess.run(
                            list(command["argv"]),
                            cwd=command["working_directory"],
                            shell=False,
                            capture_output=True,
                            text=True,
                            check=False,
                        )
                    except OSError as exc:
                        command_result.update(
                            {
                                "stdout": "",
                                "stderr": str(exc),
                                "exit_code": None,
                                "status": "failed_to_launch",
                                "finished_at": _now(),
                                "duration_seconds": round(monotonic() - timer, 6),
                            }
                        )
                        raise RealModelExecutionFailed(
                            f"{engine['engine']} job {job['job_id']} failed to launch "
                            f"stage {command['stage']}: {exc}"
                        ) from exc
                    command_result.update(
                        {
                            "stdout": completed.stdout if isinstance(completed.stdout, str) else "",
                            "stderr": completed.stderr if isinstance(completed.stderr, str) else "",
                            "exit_code": int(completed.returncode),
                            "status": "succeeded" if completed.returncode == 0 else "failed",
                            "finished_at": _now(),
                            "duration_seconds": round(monotonic() - timer, 6),
                        }
                    )
                    report["updated_at"] = _now()
                    _atomic_write_json(state_path, report)
                    if completed.returncode != 0:
                        raise RealModelExecutionFailed(
                            f"{engine['engine']} job {job['job_id']} stage "
                            f"{command['stage']} exited with code {completed.returncode}"
                        )
            except (RealModelExecutorError, OSError) as exc:
                attempt["status"] = "failed"
                attempt["finished_at"] = _now()
                attempt["error"] = str(exc)
                record["status"] = "failed"
                report["status"] = "failed_stopped_on_first_error"
                report["updated_at"] = _now()
                _atomic_write_json(state_path, report)
                if isinstance(exc, RealModelExecutionFailed):
                    raise
                raise RealModelExecutionFailed(str(exc)) from exc
            attempt["status"] = "succeeded"
            attempt["finished_at"] = _now()
            record["status"] = "succeeded"
            report["updated_at"] = _now()
            _atomic_write_json(state_path, report)

    report["status"] = "succeeded"
    report["finished_at"] = _now()
    report["updated_at"] = report["finished_at"]
    _atomic_write_json(state_path, report)
    return report


def execute_handoff(
    *,
    handoff_manifest_path: str | Path,
    runtime_attestation_path: str | Path,
    execute: bool = False,
    resume: bool = False,
    state_file: str | Path | None = None,
) -> dict[str, Any]:
    """Validate and optionally execute one unified handoff serially."""

    manifest_path = Path(handoff_manifest_path).expanduser().resolve()
    attestation_path = Path(runtime_attestation_path).expanduser().resolve()
    manifest = _load_object(manifest_path, label="unified handoff manifest")
    handoff = _validate_handoff(manifest, manifest_path=manifest_path)
    attestation_raw = _load_object(attestation_path, label="runtime attestation")
    attestation = _validate_attestation(
        attestation_raw,
        attestation_path=attestation_path,
        handoff=handoff,
    )
    # Dry-run validation includes the Germinal staging derivation but performs no
    # filesystem mutation and never acquires or launches a process.
    if not execute:
        if resume:
            raise RealModelExecutorError("--resume is meaningful only with --execute")
        return _preview(handoff, attestation)

    if state_file is None:
        state_path = handoff["handoff_root"] / "execution" / "execution_report.json"
    else:
        state_path = Path(state_file).expanduser().resolve()
    _under(state_path, handoff["handoff_root"], label="execution state file")
    # One fixed lock per handoff prevents a custom report filename from
    # accidentally creating a second concurrent single-GPU submission lane.
    lock_path = handoff["handoff_root"] / "execution" / ".single_gpu.lock"
    with _exclusive_execution_lock(lock_path):
        return _run_serial(
            handoff=handoff,
            attestation=attestation,
            state_path=state_path,
            resume=resume,
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate and optionally execute only engines[].execution_jobs from a "
            "unified RFantibody/IgGM/Germinal handoff. Default is dry-run."
        )
    )
    parser.add_argument(
        "--handoff-manifest",
        required=True,
        help="Path to unified_handoff_manifest.json (native plans are rejected)",
    )
    parser.add_argument(
        "--runtime-attestation",
        required=True,
        help="Independent hash-bound runtime attestation JSON",
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually run selected commands serially; without this flag, dry-run only",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume a matching report and skip jobs already recorded as succeeded",
    )
    parser.add_argument(
        "--state-file",
        help=(
            "Execution report path; defaults to "
            "<handoff>/execution/execution_report.json and must remain in the handoff"
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = execute_handoff(
            handoff_manifest_path=args.handoff_manifest,
            runtime_attestation_path=args.runtime_attestation,
            execute=args.execute,
            resume=args.resume,
            state_file=args.state_file,
        )
    except RealModelExecutionFailed as exc:
        print(f"EXECUTION FAILED: {exc}", file=sys.stderr)
        return 1
    except (RealModelExecutorError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
