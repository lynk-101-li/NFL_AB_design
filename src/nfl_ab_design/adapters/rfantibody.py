"""Plan RFantibody v1 jobs from the repository's normalized design request.

This module is deliberately a *planner*, not a runner.  It validates all
structure-dependent inputs and returns argv-safe command specifications for the
official three-stage RFantibody workflow::

    rfdiffusion -> proteinmpnn -> rf2

RFantibody cannot consume the normalized masked-FASTA request directly.  It
requires a target PDB, one HLT-annotated coordinate framework per template,
and an explicit mapping from full-antigen sequence positions to PDB chain and
residue numbers.  Missing or inconsistent inputs therefore raise
``RFantibodyAdapterError`` instead of producing a plausible-looking but invalid
command.

Only Python's standard library is used.  No command is executed and no output
directory is created by this module.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
import shlex
from typing import Any, Literal


RFANTIBODY_REPOSITORY = "https://github.com/RosettaCommons/RFantibody"
RFANTIBODY_OFFICIAL_MAIN_SHA = "8fe311415754e0276d1a39c87c57e69c88927a2d"
RFANTIBODY_PACKAGE_VERSION = "1.0.0"
ADAPTER_SCHEMA = "nfl_ab_design.rfantibody_plan.v1"
SUPPORTED_REQUEST_SCHEMA = "nfl_ab_design.normalized_de_novo_request.v1"
CDR_NAMES = ("H1", "H2", "H3", "L1", "L2", "L3")

RunMode = Literal["smoke", "full"]


class RFantibodyAdapterError(ValueError):
    """Raised when a normalized request cannot be mapped safely to RFantibody."""


@dataclass(frozen=True)
class RFantibodyProfile:
    """Execution parameters for a smoke or full RFantibody campaign."""

    num_backbones: int
    sequences_per_backbone: int
    diffuser_timesteps: int
    rf2_recycles: int
    proteinmpnn_temperature: float = 0.2
    omit_amino_acids: str = "CX"
    rf2_hotspot_show_fraction: float = 0.0

    def __post_init__(self) -> None:
        if self.num_backbones <= 0:
            raise RFantibodyAdapterError("num_backbones must be positive")
        if self.sequences_per_backbone <= 0:
            raise RFantibodyAdapterError("sequences_per_backbone must be positive")
        if self.diffuser_timesteps <= 0:
            raise RFantibodyAdapterError("diffuser_timesteps must be positive")
        if self.rf2_recycles <= 0:
            raise RFantibodyAdapterError("rf2_recycles must be positive")
        if self.proteinmpnn_temperature <= 0:
            raise RFantibodyAdapterError("proteinmpnn_temperature must be positive")
        if not 0.0 <= self.rf2_hotspot_show_fraction <= 1.0:
            raise RFantibodyAdapterError(
                "rf2_hotspot_show_fraction must be between 0 and 1"
            )

    @property
    def expected_sequences(self) -> int:
        """Return the expected sequence count after ProteinMPNN."""

        return self.num_backbones * self.sequences_per_backbone

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable profile."""

        return {
            "num_backbones": self.num_backbones,
            "sequences_per_backbone": self.sequences_per_backbone,
            "expected_sequences": self.expected_sequences,
            "diffuser_timesteps": self.diffuser_timesteps,
            "rf2_recycles": self.rf2_recycles,
            "proteinmpnn_temperature": self.proteinmpnn_temperature,
            "omit_amino_acids": self.omit_amino_acids,
            "rf2_hotspot_show_fraction": self.rf2_hotspot_show_fraction,
        }


DEFAULT_PROFILES: dict[RunMode, RFantibodyProfile] = {
    # The smoke profile exercises every model while keeping GPU work small.  It
    # is an integration test, not a scientifically useful design campaign.
    "smoke": RFantibodyProfile(
        num_backbones=2,
        sequences_per_backbone=1,
        diffuser_timesteps=10,
        rf2_recycles=2,
    ),
    # These values follow the scale of the official full-pipeline example.  The
    # RFantibody README notes that many real campaigns may need ~10k designs.
    "full": RFantibodyProfile(
        num_backbones=1_000,
        sequences_per_backbone=4,
        diffuser_timesteps=50,
        rf2_recycles=10,
    ),
}


@dataclass(frozen=True)
class PDBResidue:
    """A residue address in the target PDB."""

    chain_id: str
    residue_number: int
    insertion_code: str = ""

    def __post_init__(self) -> None:
        if len(self.chain_id) != 1 or not self.chain_id.isalpha():
            raise RFantibodyAdapterError(
                f"PDB chain_id must be one alphabetic character, got {self.chain_id!r}"
            )
        if self.insertion_code and (
            len(self.insertion_code) != 1 or not self.insertion_code.isalpha()
        ):
            raise RFantibodyAdapterError(
                f"PDB insertion_code must be empty or one letter, got {self.insertion_code!r}"
            )

    @property
    def hotspot_token(self) -> str:
        """Return RFantibody's ``A100`` hotspot syntax.

        RFantibody v1 parses everything after the first character with
        ``int()`` and therefore cannot represent insertion codes.  Rejecting
        them here avoids addressing the wrong residue silently.
        """

        if self.insertion_code:
            raise RFantibodyAdapterError(
                "RFantibody v1 hotspot syntax does not support PDB insertion "
                f"codes ({self.chain_id}{self.residue_number}{self.insertion_code})"
            )
        return f"{self.chain_id}{self.residue_number}"

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable residue address."""

        return {
            "chain_id": self.chain_id,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "rfantibody_token": self.hotspot_token,
        }


@dataclass(frozen=True)
class CommandSpec:
    """One non-executed command and the artifacts it is expected to create."""

    stage: str
    argv: tuple[str, ...]
    expected_outputs: tuple[str, ...]
    expected_record_count: int | None = None
    notes: str = ""

    @property
    def rendered(self) -> str:
        """Return a shell-display form; callers should execute ``argv`` directly."""

        return shlex.join(self.argv)

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable command description."""

        return {
            "stage": self.stage,
            "argv": list(self.argv),
            "rendered_for_review": self.rendered,
            "expected_outputs": list(self.expected_outputs),
            "expected_record_count": self.expected_record_count,
            "notes": self.notes,
            "execution_state": "planned_not_executed",
        }


@dataclass(frozen=True)
class RFantibodyJobSpec:
    """A single template-by-epitope RFantibody job."""

    job_id: str
    template_id: str
    framework_source_id: str
    epitope_id: str
    mode: RunMode
    target_pdb: str
    framework_hlt_pdb: str
    output_dir: str
    full_coordinate_hotspots: tuple[int, ...]
    pdb_hotspots: tuple[str, ...]
    hotspot_mapping: tuple[dict[str, Any], ...]
    hotspot_selection_source: str
    design_loops: str
    loop_lengths: dict[str, str]
    rfdiffusion_design_startnum: int
    seed: int
    profile: RFantibodyProfile
    commands: tuple[CommandSpec, ...]

    @property
    def expected_candidate_count(self) -> int:
        """Return expected RF2 best-model records for this job."""

        return self.profile.expected_sequences

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable job."""

        return {
            "job_id": self.job_id,
            "template_id": self.template_id,
            "framework_source_id": self.framework_source_id,
            "epitope_id": self.epitope_id,
            "mode": self.mode,
            "target_pdb": self.target_pdb,
            "framework_hlt_pdb": self.framework_hlt_pdb,
            "output_dir": self.output_dir,
            "full_coordinate_hotspots": list(self.full_coordinate_hotspots),
            "pdb_hotspots": list(self.pdb_hotspots),
            "hotspot_mapping": list(self.hotspot_mapping),
            "hotspot_selection_source": self.hotspot_selection_source,
            "designed_regions": list(CDR_NAMES),
            "design_loops": self.design_loops,
            "loop_lengths": dict(self.loop_lengths),
            "rfdiffusion_design_startnum": self.rfdiffusion_design_startnum,
            "seed": self.seed,
            "profile": self.profile.as_dict(),
            "expected_candidate_count": self.expected_candidate_count,
            "commands": [command.as_dict() for command in self.commands],
            "execution_state": "planned_not_executed",
        }


@dataclass(frozen=True)
class RFantibodyPlan:
    """Four-job dual-template/dual-epitope RFantibody execution plan."""

    request_schema: str
    request_sha256: str
    mode: RunMode
    seed: int
    target_pdb: str
    target_pdb_sha256: str
    framework_pdb_sha256: dict[str, str]
    profile: RFantibodyProfile
    jobs: tuple[RFantibodyJobSpec, ...]
    command_prefix: tuple[str, ...]
    runtime_ref: str

    @property
    def expected_candidate_count(self) -> int:
        """Return the total number of expected sequence/model records."""

        return sum(job.expected_candidate_count for job in self.jobs)

    @property
    def commands(self) -> tuple[CommandSpec, ...]:
        """Return all job commands in deterministic execution order."""

        return tuple(command for job in self.jobs for command in job.commands)

    def as_dict(self) -> dict[str, Any]:
        """Return a complete JSON-serializable plan manifest."""

        return {
            "schema": ADAPTER_SCHEMA,
            "execution_state": "planned_not_executed",
            "result_provenance": "adapter_plan_only_no_model_results",
            "normalized_request_schema": self.request_schema,
            "normalized_request_sha256": self.request_sha256,
            "mode": self.mode,
            "seed": self.seed,
            "target_pdb": self.target_pdb,
            "target_pdb_sha256": self.target_pdb_sha256,
            "framework_pdb_sha256": dict(self.framework_pdb_sha256),
            "profile": self.profile.as_dict(),
            "job_count": len(self.jobs),
            "expected_candidate_count": self.expected_candidate_count,
            "designed_regions": list(CDR_NAMES),
            "command_prefix": list(self.command_prefix),
            "jobs": [job.as_dict() for job in self.jobs],
            "provenance": {
                "adapter": "nfl_ab_design.adapters.rfantibody",
                "adapter_schema": ADAPTER_SCHEMA,
                "official_repository": RFANTIBODY_REPOSITORY,
                "official_package_version": RFANTIBODY_PACKAGE_VERSION,
                "official_main_sha_verified_during_adapter_development": (
                    RFANTIBODY_OFFICIAL_MAIN_SHA
                ),
                "requested_runtime_ref": self.runtime_ref,
                "license": "MIT",
                "model_pipeline": ["rfdiffusion", "proteinmpnn", "rf2"],
                "determinism_notes": [
                    "RFantibody v1 RFdiffusion has no numeric --seed option; "
                    "the adapter uses --deterministic plus inference.design_startnum.",
                    "RFantibody v1 ProteinMPNN --deterministic uses the runtime's "
                    "fixed deterministic seed.",
                    "RF2 receives the campaign seed through --seed.",
                ],
            },
            "required_runtime": {
                "operating_system": "Linux (Ubuntu 22.04 recommended upstream)",
                "gpu": "NVIDIA GPU with CUDA 11.8+",
                "python": "3.10",
                "preflight_commands": [
                    "nvidia-smi",
                    shlex.join((*self.command_prefix, "rfdiffusion", "--help")),
                    shlex.join((*self.command_prefix, "proteinmpnn", "--help")),
                    shlex.join((*self.command_prefix, "rf2", "--help")),
                ],
            },
        }

    def to_json(self, *, indent: int = 2) -> str:
        """Serialize the plan without writing it to disk."""

        return json.dumps(self.as_dict(), indent=indent, ensure_ascii=False)


_HLT_REMARK = re.compile(
    r"^REMARK\s+PDBinfo-LABEL:\s+(?P<absolute_index>\d+)\s+"
    r"(?P<cdr>H1|H2|H3|L1|L2|L3)\s*$"
)
_LOOP_LENGTH = re.compile(r"^(?P<minimum>\d+)(?:-(?P<maximum>\d+))?$")


def _canonical_request_hash(request: Mapping[str, Any]) -> str:
    payload = json.dumps(
        request,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_file(value: str | Path | None, *, label: str) -> Path:
    if value is None or not str(value).strip():
        raise RFantibodyAdapterError(f"Missing required {label}")
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise RFantibodyAdapterError(f"Required {label} is not a file: {path}")
    return path


def _pdb_residues(path: Path) -> set[tuple[str, int, str]]:
    residues: set[tuple[str, int, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or len(line) < 27:
                continue
            chain = line[21].strip()
            number_text = line[22:26].strip()
            insertion_code = line[26].strip()
            if not chain or not number_text:
                continue
            try:
                number = int(number_text)
            except ValueError:
                continue
            residues.add((chain, number, insertion_code))
    if not residues:
        raise RFantibodyAdapterError(f"No ATOM residues found in PDB: {path}")
    return residues


def _parse_hlt_framework(path: Path) -> dict[str, str]:
    """Validate an HLT PDB and return fixed loop lengths from its REMARKs."""

    chains: set[str] = set()
    loop_indices: dict[str, set[int]] = {name: set() for name in CDR_NAMES}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM  ") and len(line) >= 22:
                chain = line[21].strip()
                if chain:
                    chains.add(chain)
            match = _HLT_REMARK.match(line.rstrip("\n"))
            if match:
                loop_indices[match.group("cdr")].add(
                    int(match.group("absolute_index"))
                )
    missing_chains = {"H", "L"} - chains
    if missing_chains:
        raise RFantibodyAdapterError(
            f"HLT framework {path} is missing chain(s): {sorted(missing_chains)}"
        )
    extra = chains - {"H", "L"}
    if extra:
        raise RFantibodyAdapterError(
            f"Framework-only HLT PDB {path} contains unexpected chain(s) "
            f"{sorted(extra)}; pass the target separately with --target"
        )
    missing_loops = [name for name, values in loop_indices.items() if not values]
    if missing_loops:
        raise RFantibodyAdapterError(
            f"HLT framework {path} lacks PDBinfo-LABEL REMARKs for "
            f"{', '.join(missing_loops)}"
        )
    return {name: str(len(loop_indices[name])) for name in CDR_NAMES}


def _validate_loop_lengths(values: Mapping[str, str | int]) -> dict[str, str]:
    missing = [name for name in CDR_NAMES if name not in values]
    extra = sorted(set(values) - set(CDR_NAMES))
    if missing or extra:
        raise RFantibodyAdapterError(
            f"Loop lengths must define exactly {CDR_NAMES}; missing={missing}, extra={extra}"
        )
    normalized: dict[str, str] = {}
    for name in CDR_NAMES:
        value = str(values[name]).strip()
        match = _LOOP_LENGTH.fullmatch(value)
        if not match:
            raise RFantibodyAdapterError(
                f"Invalid RFantibody loop length/range for {name}: {value!r}"
            )
        minimum = int(match.group("minimum"))
        maximum = int(match.group("maximum") or minimum)
        if minimum <= 0 or maximum < minimum:
            raise RFantibodyAdapterError(
                f"Invalid RFantibody loop length/range for {name}: {value!r}"
            )
        normalized[name] = value
    return normalized


def _normalize_mapping_entry(value: Mapping[str, Any], *, position: int) -> PDBResidue:
    chain = value.get("chain_id", value.get("chain"))
    residue_number = value.get(
        "residue_number", value.get("pdb_residue_number", value.get("resseq"))
    )
    insertion_code = value.get("insertion_code", value.get("icode", ""))
    if chain is None or residue_number is None:
        raise RFantibodyAdapterError(
            f"Mapping for full position {position} must contain chain_id and residue_number"
        )
    try:
        number = int(residue_number)
    except (TypeError, ValueError) as exc:
        raise RFantibodyAdapterError(
            f"Invalid PDB residue number for full position {position}: {residue_number!r}"
        ) from exc
    return PDBResidue(str(chain).strip(), number, str(insertion_code).strip())


def _normalize_coordinate_mapping(
    raw: Mapping[Any, Any] | Sequence[Mapping[str, Any]] | None,
) -> dict[int, PDBResidue]:
    if raw is None:
        raise RFantibodyAdapterError(
            "Missing explicit full-coordinate-to-PDB residue mapping"
        )
    normalized: dict[int, PDBResidue] = {}
    if isinstance(raw, Mapping):
        items = raw.items()
        for raw_position, value in items:
            try:
                position = int(raw_position)
            except (TypeError, ValueError) as exc:
                raise RFantibodyAdapterError(
                    f"Invalid full-sequence coordinate: {raw_position!r}"
                ) from exc
            if not isinstance(value, Mapping):
                raise RFantibodyAdapterError(
                    f"Mapping value for full position {position} must be an object"
                )
            normalized[position] = _normalize_mapping_entry(value, position=position)
    elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        for value in raw:
            if not isinstance(value, Mapping):
                raise RFantibodyAdapterError("Coordinate mapping rows must be objects")
            raw_position = value.get(
                "full_position", value.get("full_coordinate", value.get("position"))
            )
            try:
                position = int(raw_position)
            except (TypeError, ValueError) as exc:
                raise RFantibodyAdapterError(
                    f"Invalid full-sequence coordinate: {raw_position!r}"
                ) from exc
            if position in normalized:
                raise RFantibodyAdapterError(
                    f"Duplicate mapping for full-sequence coordinate {position}"
                )
            normalized[position] = _normalize_mapping_entry(value, position=position)
    else:
        raise RFantibodyAdapterError(
            "full-coordinate-to-PDB mapping must be an object or a list of objects"
        )
    if not normalized:
        raise RFantibodyAdapterError("Full-coordinate-to-PDB mapping is empty")
    return normalized


def _epitope_hotspots(epitope: Mapping[str, Any]) -> tuple[tuple[int, ...], str]:
    for field in (
        "rfantibody_hotspot_residue_indices",
        "selected_hotspot_residue_indices",
    ):
        raw = epitope.get(field)
        if raw is None:
            continue
        if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray)):
            raise RFantibodyAdapterError(
                f"Epitope {epitope.get('epitope_id', epitope.get('id'))!r} field "
                f"{field} must be a list of full-sequence positions"
            )
        try:
            positions = tuple(sorted({int(value) for value in raw}))
        except (TypeError, ValueError) as exc:
            raise RFantibodyAdapterError(
                f"Epitope hotspot field {field} contains a non-integer position"
            ) from exc
        if not positions or positions[0] <= 0:
            raise RFantibodyAdapterError(
                f"Epitope hotspot field {field} must contain positive positions"
            )
        return positions, field
    raise RFantibodyAdapterError(
        f"Epitope {epitope.get('epitope_id', epitope.get('id'))!r} has no explicit "
        "RFantibody hotspot positions. candidate_hotspot_residue_indices is only "
        "an epitope window and is intentionally not accepted as a selected hotspot list"
    )


def _slug(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return cleaned or "unnamed"


def _string_id(row: Mapping[str, Any], *keys: str, label: str) -> str:
    for key in keys:
        value = row.get(key)
        if value is not None and str(value).strip():
            return str(value).strip()
    raise RFantibodyAdapterError(f"Missing {label}; expected one of {keys}")


def _command(
    prefix: tuple[str, ...],
    executable: str,
    *arguments: str,
) -> tuple[str, ...]:
    return (*prefix, executable, *arguments)


def build_rfantibody_plan(
    request: Mapping[str, Any],
    *,
    target_pdb: str | Path | None = None,
    framework_hlt_pdbs: Mapping[str, str | Path] | None = None,
    full_coordinate_to_pdb: (
        Mapping[Any, Any] | Sequence[Mapping[str, Any]] | None
    ) = None,
    output_root: str | Path = "real_runs/adapter_plans/rfantibody",
    mode: RunMode = "smoke",
    seed: int = 20260812,
    command_prefix: Sequence[str] = (),
    profile: RFantibodyProfile | None = None,
    loop_length_overrides: (
        Mapping[str, Mapping[str, str | int]] | None
    ) = None,
    runtime_ref: str = RFANTIBODY_OFFICIAL_MAIN_SHA,
) -> RFantibodyPlan:
    """Build a validated, non-executed RFantibody plan.

    Parameters
    ----------
    request:
        A ``nfl_ab_design.normalized_de_novo_request.v1`` request.
    target_pdb:
        Coordinate PDB used by ``rfdiffusion --target``.  If omitted, the
        adapter reads ``request['antigen']['antigen_pdb_path']``.
    framework_hlt_pdbs:
        Mapping from each template ID (or framework source ID) to a distinct
        framework-only HLT PDB.  A template may alternatively contain
        ``framework_hlt_pdb_path`` itself.
    full_coordinate_to_pdb:
        Explicit mapping from 1-based full-antigen coordinates to
        ``chain_id``, ``residue_number``, and optional ``insertion_code``.
        It may be a dictionary keyed by full coordinate or a list whose rows
        contain ``full_position``.  If omitted, the mapping is read from the
        antigen object's ``full_coordinate_to_pdb`` or
        ``full_coordinate_to_pdb_mapping`` field.
    output_root:
        Planned output root.  It is not created.
    mode:
        ``smoke`` runs 2 backbones x 1 sequence for each of the four jobs;
        ``full`` defaults to 1000 x 4.  Pass ``profile`` to override.
    seed:
        Campaign seed.  RFdiffusion v1 uses it as deterministic design-start
        numbering, while RF2 receives it via ``--seed``.
    command_prefix:
        Optional argv prefix such as ``('uv', 'run')`` or an Apptainer wrapper.
    loop_length_overrides:
        Optional mapping keyed by template ID/source ID.  Each value must
        define all six CDRs using fixed lengths (``7``) or ranges (``5-13``).
        Without overrides, fixed lengths are derived from HLT REMARK labels.
    runtime_ref:
        Installed RFantibody git SHA/tag recorded in plan provenance.

    Returns
    -------
    RFantibodyPlan
        Four deterministic template-by-epitope jobs and their exact argv.

    Raises
    ------
    RFantibodyAdapterError
        If any required structure, HLT annotation, residue mapping, hotspot,
        or campaign-shape constraint is missing or inconsistent.
    """

    if mode not in DEFAULT_PROFILES:
        raise RFantibodyAdapterError(f"Unsupported RFantibody mode: {mode!r}")
    selected_profile = profile or DEFAULT_PROFILES[mode]
    if not isinstance(selected_profile, RFantibodyProfile):
        raise RFantibodyAdapterError("profile must be an RFantibodyProfile")
    try:
        seed = int(seed)
    except (TypeError, ValueError) as exc:
        raise RFantibodyAdapterError(f"Invalid seed: {seed!r}") from exc
    if seed < 0:
        raise RFantibodyAdapterError("seed must be non-negative")

    request_schema = str(request.get("schema", "")).strip()
    if request_schema != SUPPORTED_REQUEST_SCHEMA:
        raise RFantibodyAdapterError(
            f"Unsupported normalized request schema {request_schema!r}; expected "
            f"{SUPPORTED_REQUEST_SCHEMA!r}"
        )
    engine = str(request.get("engine", "RFantibody")).strip().casefold()
    if engine != "rfantibody":
        raise RFantibodyAdapterError(
            f"Request engine is {request.get('engine')!r}, not RFantibody"
        )

    antigen = request.get("antigen")
    if not isinstance(antigen, Mapping):
        raise RFantibodyAdapterError("Normalized request lacks an antigen object")
    target_value = target_pdb or antigen.get("antigen_pdb_path")
    target_path = _require_file(target_value, label="target antigen PDB")
    target_residues = _pdb_residues(target_path)

    raw_mapping = full_coordinate_to_pdb
    if raw_mapping is None:
        raw_mapping = antigen.get(
            "full_coordinate_to_pdb",
            antigen.get("full_coordinate_to_pdb_mapping"),
        )
    coordinate_mapping = _normalize_coordinate_mapping(raw_mapping)

    raw_templates = request.get("templates")
    raw_epitopes = request.get("epitopes")
    if not isinstance(raw_templates, Sequence) or isinstance(
        raw_templates, (str, bytes, bytearray)
    ):
        raise RFantibodyAdapterError("Normalized request templates must be a list")
    if not isinstance(raw_epitopes, Sequence) or isinstance(
        raw_epitopes, (str, bytes, bytearray)
    ):
        raise RFantibodyAdapterError("Normalized request epitopes must be a list")
    templates = tuple(raw_templates)
    epitopes = tuple(raw_epitopes)
    if len(templates) != 2 or len(epitopes) != 2:
        raise RFantibodyAdapterError(
            "This campaign adapter requires exactly 2 templates and 2 epitopes "
            f"(received {len(templates)} and {len(epitopes)})"
        )
    if not all(isinstance(row, Mapping) for row in (*templates, *epitopes)):
        raise RFantibodyAdapterError("Template and epitope list entries must be objects")

    prefix = tuple(str(part) for part in command_prefix)
    if any(not part for part in prefix):
        raise RFantibodyAdapterError("command_prefix cannot contain empty argv entries")
    output_path = Path(output_root).expanduser().resolve()
    framework_values = framework_hlt_pdbs or {}
    loop_overrides = loop_length_overrides or {}

    prepared_templates: list[tuple[str, str, Path, dict[str, str]]] = []
    framework_hashes: dict[str, str] = {}
    seen_template_ids: set[str] = set()
    for raw_template in templates:
        template = raw_template  # narrowed by the all(isinstance(...)) check above
        template_id = _string_id(template, "template_id", label="template ID")
        source_id = _string_id(
            template,
            "framework_source_id",
            "source_antibody_id",
            label=f"framework source ID for {template_id}",
        )
        if template_id in seen_template_ids:
            raise RFantibodyAdapterError(f"Duplicate template ID: {template_id}")
        seen_template_ids.add(template_id)
        framework_value = (
            framework_values.get(template_id)
            or framework_values.get(source_id)
            or template.get("framework_hlt_pdb_path")
        )
        framework_path = _require_file(
            framework_value, label=f"HLT framework PDB for template {template_id}"
        )
        derived_lengths = _parse_hlt_framework(framework_path)
        override = loop_overrides.get(template_id) or loop_overrides.get(source_id)
        loop_lengths = _validate_loop_lengths(override or derived_lengths)
        prepared_templates.append(
            (template_id, source_id, framework_path, loop_lengths)
        )
        framework_hashes[template_id] = _file_sha256(framework_path)

    prepared_epitopes: list[
        tuple[str, tuple[int, ...], tuple[PDBResidue, ...], str]
    ] = []
    seen_epitope_ids: set[str] = set()
    for raw_epitope in epitopes:
        epitope = raw_epitope
        epitope_id = _string_id(
            epitope, "epitope_id", "id", label="epitope ID"
        )
        if epitope_id in seen_epitope_ids:
            raise RFantibodyAdapterError(f"Duplicate epitope ID: {epitope_id}")
        seen_epitope_ids.add(epitope_id)
        positions, selection_source = _epitope_hotspots(epitope)
        missing_positions = [
            position for position in positions if position not in coordinate_mapping
        ]
        if missing_positions:
            raise RFantibodyAdapterError(
                f"Epitope {epitope_id} hotspot positions lack explicit PDB mapping: "
                f"{missing_positions}"
            )
        mapped = tuple(coordinate_mapping[position] for position in positions)
        for position, pdb_residue in zip(positions, mapped, strict=True):
            # Calling hotspot_token also rejects insertion codes unsupported by v1.
            pdb_residue.hotspot_token
            address = (
                pdb_residue.chain_id,
                pdb_residue.residue_number,
                pdb_residue.insertion_code,
            )
            if address not in target_residues:
                raise RFantibodyAdapterError(
                    f"Full position {position} maps to absent target PDB residue "
                    f"{pdb_residue.chain_id}{pdb_residue.residue_number}"
                    f"{pdb_residue.insertion_code}: {target_path}"
                )
        prepared_epitopes.append((epitope_id, positions, mapped, selection_source))

    jobs: list[RFantibodyJobSpec] = []
    job_index = 0
    for template_id, source_id, framework_path, loop_lengths in prepared_templates:
        for epitope_id, positions, mapped, selection_source in prepared_epitopes:
            digest = sha256(f"{template_id}\0{epitope_id}".encode("utf-8")).hexdigest()[:8]
            job_id = f"rfab_{_slug(template_id)}__{_slug(epitope_id)}__{digest}"
            job_dir = output_path / job_id
            diffusion_qv = job_dir / "1_rfdiffusion.qv"
            diffusion_sc = diffusion_qv.with_suffix(".sc")
            proteinmpnn_qv = job_dir / "2_proteinmpnn.qv"
            proteinmpnn_sc = proteinmpnn_qv.with_suffix(".sc")
            rf2_qv = job_dir / "3_rf2.qv"
            rf2_sc = rf2_qv.with_suffix(".sc")
            final_pdb_dir = job_dir / "final_pdbs"
            design_loops = ",".join(
                f"{name}:{loop_lengths[name]}" for name in CDR_NAMES
            )
            pdb_hotspots = tuple(residue.hotspot_token for residue in mapped)
            design_startnum = seed + job_index * selected_profile.num_backbones
            rf2_seed = seed + job_index

            commands = (
                CommandSpec(
                    stage="prepare_output_directory",
                    argv=("mkdir", "-p", str(job_dir)),
                    expected_outputs=(str(job_dir),),
                    notes=(
                        "RFantibody's Quiver output path does not create its parent "
                        "directory; execute this preparation step first."
                    ),
                ),
                CommandSpec(
                    stage="rfdiffusion_backbone_design",
                    argv=_command(
                        prefix,
                        "rfdiffusion",
                        "--target",
                        str(target_path),
                        "--framework",
                        str(framework_path),
                        "--output-quiver",
                        str(diffusion_qv),
                        "--num-designs",
                        str(selected_profile.num_backbones),
                        "--design-loops",
                        design_loops,
                        "--hotspots",
                        ",".join(pdb_hotspots),
                        "--diffuser-t",
                        str(selected_profile.diffuser_timesteps),
                        "--deterministic",
                        "--no-trajectory",
                        "--extra",
                        f"inference.design_startnum={design_startnum}",
                    ),
                    expected_outputs=(str(diffusion_qv),),
                    expected_record_count=selected_profile.num_backbones,
                    notes="All six CDR backbones and their target-relative docks are designed.",
                ),
                CommandSpec(
                    stage="rfdiffusion_score_export",
                    argv=_command(prefix, "qvscorefile", str(diffusion_qv)),
                    expected_outputs=(str(diffusion_sc),),
                    expected_record_count=selected_profile.num_backbones,
                ),
                CommandSpec(
                    stage="proteinmpnn_six_cdr_sequence_design",
                    argv=_command(
                        prefix,
                        "proteinmpnn",
                        "--input-quiver",
                        str(diffusion_qv),
                        "--output-quiver",
                        str(proteinmpnn_qv),
                        "--loops",
                        ",".join(CDR_NAMES),
                        "--seqs-per-struct",
                        str(selected_profile.sequences_per_backbone),
                        "--temperature",
                        str(selected_profile.proteinmpnn_temperature),
                        "--omit-aas",
                        selected_profile.omit_amino_acids,
                        "--deterministic",
                    ),
                    expected_outputs=(str(proteinmpnn_qv),),
                    expected_record_count=selected_profile.expected_sequences,
                    notes="ProteinMPNN redesigns H1,H2,H3,L1,L2,L3 explicitly.",
                ),
                CommandSpec(
                    stage="proteinmpnn_score_export",
                    argv=_command(prefix, "qvscorefile", str(proteinmpnn_qv)),
                    expected_outputs=(str(proteinmpnn_sc),),
                    expected_record_count=selected_profile.expected_sequences,
                ),
                CommandSpec(
                    stage="rf2_structure_prediction",
                    argv=_command(
                        prefix,
                        "rf2",
                        "--input-quiver",
                        str(proteinmpnn_qv),
                        "--output-quiver",
                        str(rf2_qv),
                        "--num-recycles",
                        str(selected_profile.rf2_recycles),
                        "--seed",
                        str(rf2_seed),
                        "--hotspot-show-prop",
                        str(selected_profile.rf2_hotspot_show_fraction),
                    ),
                    expected_outputs=(str(rf2_qv),),
                    expected_record_count=selected_profile.expected_sequences,
                    notes=(
                        "RF2 Quiver scores should include interaction_pae, pred_lddt, "
                        "target_aligned_antibody_rmsd, and target_aligned_cdr_rmsd."
                    ),
                ),
                CommandSpec(
                    stage="rf2_score_export",
                    argv=_command(prefix, "qvscorefile", str(rf2_qv)),
                    expected_outputs=(str(rf2_sc),),
                    expected_record_count=selected_profile.expected_sequences,
                    notes=(
                        "Upstream minimal filters are RF2 pAE < 10 and design-versus-"
                        "prediction RMSD < 2 Angstrom; retain raw metrics for audit."
                    ),
                ),
                CommandSpec(
                    stage="rf2_pdb_extraction",
                    argv=_command(
                        prefix,
                        "qvextract",
                        str(rf2_qv),
                        "--output-dir",
                        str(final_pdb_dir),
                    ),
                    expected_outputs=(str(final_pdb_dir),),
                    expected_record_count=selected_profile.expected_sequences,
                ),
            )
            hotspot_mapping = tuple(
                {
                    "full_position_1_based": position,
                    **residue.as_dict(),
                }
                for position, residue in zip(positions, mapped, strict=True)
            )
            jobs.append(
                RFantibodyJobSpec(
                    job_id=job_id,
                    template_id=template_id,
                    framework_source_id=source_id,
                    epitope_id=epitope_id,
                    mode=mode,
                    target_pdb=str(target_path),
                    framework_hlt_pdb=str(framework_path),
                    output_dir=str(job_dir),
                    full_coordinate_hotspots=positions,
                    pdb_hotspots=pdb_hotspots,
                    hotspot_mapping=hotspot_mapping,
                    hotspot_selection_source=selection_source,
                    design_loops=design_loops,
                    loop_lengths=dict(loop_lengths),
                    rfdiffusion_design_startnum=design_startnum,
                    seed=seed,
                    profile=selected_profile,
                    commands=commands,
                )
            )
            job_index += 1

    if len(jobs) != 4:
        raise RFantibodyAdapterError(
            f"Internal campaign-shape error: expected 4 jobs, built {len(jobs)}"
        )
    return RFantibodyPlan(
        request_schema=request_schema,
        request_sha256=_canonical_request_hash(request),
        mode=mode,
        seed=seed,
        target_pdb=str(target_path),
        target_pdb_sha256=_file_sha256(target_path),
        framework_pdb_sha256=framework_hashes,
        profile=selected_profile,
        jobs=tuple(jobs),
        command_prefix=prefix,
        runtime_ref=str(runtime_ref),
    )


def plan_rfantibody_jobs(
    request: Mapping[str, Any], **kwargs: Any
) -> RFantibodyPlan:
    """Compatibility alias for :func:`build_rfantibody_plan`."""

    return build_rfantibody_plan(request, **kwargs)


__all__ = [
    "ADAPTER_SCHEMA",
    "CDR_NAMES",
    "CommandSpec",
    "DEFAULT_PROFILES",
    "PDBResidue",
    "RFANTIBODY_OFFICIAL_MAIN_SHA",
    "RFANTIBODY_PACKAGE_VERSION",
    "RFANTIBODY_REPOSITORY",
    "RFantibodyAdapterError",
    "RFantibodyJobSpec",
    "RFantibodyPlan",
    "RFantibodyProfile",
    "SUPPORTED_REQUEST_SCHEMA",
    "build_rfantibody_plan",
    "plan_rfantibody_jobs",
]
