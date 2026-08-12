"""Deterministic, auditable simulation of an epitope-conditioned antibody campaign.

This module deliberately separates prospective design from a retrospective
positive-control demonstration. It does not call RFantibody, IgGM, tFold, or a
structure predictor, and every modeled metric is marked as simulated.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
import re
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
DESIGN_ALPHABET = "ADEFGHIKLNQRSTVWY"
DESIGN_REGIONS = ("H1", "H2", "H3", "L1", "L2", "L3")
REGION_SPECS = {
    "H1": ("VH", 26, 35),
    "H2": ("VH", 50, 65),
    "L1": ("VL", 24, 34),
    "L2": ("VL", 50, 56),
}
HYDROPHOBIC = set("AVLIMFWYC")
AROMATIC = set("FWY")
BASIC = set("KRH")
ACIDIC = set("DE")


@dataclass(frozen=True)
class TemplateFramework:
    template_id: str
    framework_source_antibody_id: str
    vh_id: str
    vl_id: str
    vh_source: str
    vl_source: str
    vh_framework_masked: str
    vl_framework_masked: str
    regions: tuple[tuple[str, str, int, int], ...]


@dataclass(frozen=True)
class DesignPipelineResult:
    template_rows: list[dict[str, Any]]
    generation_rows: list[dict[str, Any]]
    structure_rows: list[dict[str, Any]]
    developability_rows: list[dict[str, Any]]
    funnel_rows: list[dict[str, Any]]
    prospective_ranking_rows: list[dict[str, Any]]
    retrospective_ranking_rows: list[dict[str, Any]]
    selected_candidates: list[dict[str, Any]]

    @property
    def library(self) -> list[dict[str, Any]]:
        return self.generation_rows

    @property
    def structure(self) -> list[dict[str, Any]]:
        return self.structure_rows

    @property
    def developability(self) -> list[dict[str, Any]]:
        return self.developability_rows

    @property
    def funnel(self) -> list[dict[str, Any]]:
        return self.funnel_rows

    @property
    def prospective_ranking(self) -> list[dict[str, Any]]:
        return self.prospective_ranking_rows

    @property
    def retrospective_ranking(self) -> list[dict[str, Any]]:
        return self.retrospective_ranking_rows

    def as_dict(self) -> dict[str, Any]:
        return {
            "template_rows": self.template_rows,
            "generation_rows": self.generation_rows,
            "structure_rows": self.structure_rows,
            "developability_rows": self.developability_rows,
            "funnel_rows": self.funnel_rows,
            "prospective_ranking_rows": self.prospective_ranking_rows,
            "retrospective_ranking_rows": self.retrospective_ranking_rows,
            "selected_candidates": self.selected_candidates,
        }


def _clamp(value: float, low: float = 0.0, high: float = 100.0) -> float:
    return max(low, min(high, value))


def _stable_int(*parts: Any) -> int:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:16], 16)


def _source_value(source: Any, name: str) -> str:
    if isinstance(source, Mapping):
        return str(source[name])
    return str(getattr(source, name))


def _find_cdr3(sequence: str, chain: str) -> tuple[int, int]:
    if chain == "VH":
        region_start = max(0, len(sequence) - 55)
        match = re.search(r"C([A-Z]{3,26}?)(W[A-Z]QG|WGQG)", sequence[region_start:])
        fallback = (max(1, len(sequence) - 20), max(1, len(sequence) - 8))
    else:
        region_start = max(0, len(sequence) - 45)
        match = re.search(r"C([A-Z]{3,18}?)(FGGG|WGGG)", sequence[region_start:])
        fallback = (max(1, len(sequence) - 18), max(1, len(sequence) - 8))
    if not match:
        return fallback
    return region_start + match.start(1) + 1, region_start + match.end(1)


def _template_regions(vh: str, vl: str) -> tuple[tuple[str, str, int, int], ...]:
    regions = [
        ("H1", "VH", 26, min(35, len(vh))),
        ("H2", "VH", 50, min(65, len(vh))),
        ("H3", "VH", *_find_cdr3(vh, "VH")),
        ("L1", "VL", 24, min(34, len(vl))),
        ("L2", "VL", 50, min(56, len(vl))),
        ("L3", "VL", *_find_cdr3(vl, "VL")),
    ]
    return tuple(regions)


def _mask_regions(vh: str, vl: str, regions: Iterable[tuple[str, str, int, int]]) -> tuple[str, str]:
    chains = {"VH": list(vh), "VL": list(vl)}
    for _name, chain, start, end in regions:
        for position in range(start - 1, end):
            chains[chain][position] = "X"
    return "".join(chains["VH"]), "".join(chains["VL"])


def _validated_source_index(validated_antibodies: Sequence[Any]) -> dict[str, Any]:
    sources: dict[str, Any] = {}
    for source in validated_antibodies:
        source_id = _source_value(source, "antibody_id").strip()
        if not source_id:
            raise ValueError("Every validated antibody must have a non-empty antibody_id")
        if source_id in sources:
            raise ValueError(f"Duplicate validated antibody_id is ambiguous: {source_id}")
        sources[source_id] = source
    return sources


def _configured_template_sources(
    validated_antibodies: Sequence[Any],
    template_specs: Sequence[Mapping[str, Any]] | None,
) -> list[tuple[str, Any]]:
    """Resolve canonical template IDs in configuration order, failing closed."""

    source_index = _validated_source_index(validated_antibodies)
    if template_specs is None:
        return [(f"template_{source_id}", source) for source_id, source in source_index.items()]
    if isinstance(template_specs, (str, bytes)) or not isinstance(template_specs, Sequence):
        raise ValueError("template_specs must be a sequence of template mappings")
    if not template_specs:
        raise ValueError("template_specs cannot be empty")

    configured: list[tuple[str, Any]] = []
    seen_template_ids: set[str] = set()
    seen_source_ids: set[str] = set()
    required_safety_switches = {
        "use_framework_residues": True,
        "use_known_cdr_sequences": False,
        "use_full_known_sequence_as_generation_feature": False,
    }
    optional_safety_switches = {
        "mask_source_cdrs_before_generation": True,
        "allow_known_cdr_feature_leakage": False,
        "allow_full_known_sequence_feature_leakage": False,
    }
    for index, raw_spec in enumerate(template_specs, start=1):
        if not isinstance(raw_spec, Mapping):
            raise ValueError(f"template_specs[{index - 1}] must be a mapping")
        template_id = str(raw_spec.get("template_id", "")).strip()
        source_id = str(raw_spec.get("source_antibody_id", "")).strip()
        if not template_id:
            raise ValueError(f"template_specs[{index - 1}] is missing canonical template_id")
        if not source_id:
            raise ValueError(f"Template {template_id} is missing source_antibody_id")
        if template_id in seen_template_ids:
            raise ValueError(f"Duplicate canonical template_id: {template_id}")
        if source_id in seen_source_ids:
            raise ValueError(f"Duplicate framework source_antibody_id: {source_id}")
        if source_id not in source_index:
            raise ValueError(f"Template {template_id} references missing source antibody: {source_id}")
        if "role" in raw_spec and raw_spec["role"] != "framework_source_only":
            raise ValueError(f"Template {template_id} role must be framework_source_only")
        for switch, safe_value in required_safety_switches.items():
            if switch not in raw_spec:
                raise ValueError(f"Template {template_id} is missing required safety switch {switch}")
            if not isinstance(raw_spec[switch], bool) or raw_spec[switch] is not safe_value:
                raise ValueError(f"Template {template_id} has unsafe {switch}={raw_spec[switch]!r}")
        for switch, safe_value in optional_safety_switches.items():
            if switch in raw_spec and (
                not isinstance(raw_spec[switch], bool) or raw_spec[switch] is not safe_value
            ):
                raise ValueError(f"Template {template_id} has unsafe {switch}={raw_spec[switch]!r}")
        configured.append((template_id, source_index[source_id]))
        seen_template_ids.add(template_id)
        seen_source_ids.add(source_id)
    return configured


def _validated_configured_regions(
    template_id: str,
    vh: str,
    vl: str,
    cdr_ranges_by_template: Mapping[str, Mapping[str, Mapping[str, Any]]],
) -> tuple[tuple[str, str, int, int], ...]:
    raw_regions = cdr_ranges_by_template[template_id]
    if not isinstance(raw_regions, Mapping):
        raise ValueError(f"CDR ranges for {template_id} must be a mapping")
    configured_names = {str(name) for name in raw_regions}
    required_names = set(DESIGN_REGIONS)
    if configured_names != required_names:
        missing = sorted(required_names - configured_names)
        extra = sorted(configured_names - required_names)
        raise ValueError(
            f"CDR ranges for {template_id} must cover exactly {','.join(DESIGN_REGIONS)}; "
            f"missing={missing}, extra={extra}"
        )

    chain_sequences = {"VH": vh, "VL": vl}
    intervals_by_chain: dict[str, list[tuple[int, int, str]]] = {"VH": [], "VL": []}
    regions: list[tuple[str, str, int, int]] = []
    for name in DESIGN_REGIONS:
        raw_region = raw_regions[name]
        if not isinstance(raw_region, Mapping):
            raise ValueError(f"CDR range {template_id}/{name} must be a mapping")
        chain = raw_region.get("chain")
        expected_chain = "VH" if name.startswith("H") else "VL"
        if chain != expected_chain:
            raise ValueError(
                f"CDR range {template_id}/{name} must use chain {expected_chain}, got {chain!r}"
            )
        start = raw_region.get("start")
        end = raw_region.get("end")
        if (
            not isinstance(start, int)
            or isinstance(start, bool)
            or not isinstance(end, int)
            or isinstance(end, bool)
        ):
            raise ValueError(f"CDR range {template_id}/{name} start/end must be integers")
        if start < 1 or end < start or end > len(chain_sequences[chain]):
            raise ValueError(
                f"CDR range {template_id}/{name}={start}-{end} is outside {chain} "
                f"bounds 1-{len(chain_sequences[chain])}"
            )
        regions.append((name, chain, start, end))
        intervals_by_chain[chain].append((start, end, name))

    for chain, intervals in intervals_by_chain.items():
        intervals.sort()
        for previous, current in zip(intervals, intervals[1:]):
            if current[0] <= previous[1]:
                raise ValueError(
                    f"Overlapping CDR ranges for {template_id}/{chain}: "
                    f"{previous[2]}={previous[0]}-{previous[1]} and "
                    f"{current[2]}={current[0]}-{current[1]}"
                )
    return tuple(regions)


def load_template_frameworks(
    validated_antibodies: Sequence[Any],
    *,
    template_specs: Sequence[Mapping[str, Any]] | None = None,
    cdr_ranges_by_template: Mapping[str, Mapping[str, Mapping[str, Any]]] | None = None,
) -> list[TemplateFramework]:
    """Extract masked framework templates without copying known CDRs.

    When configuration is supplied, template identity/order and all six CDR
    coordinates are authoritative. Invalid or incomplete configuration is
    rejected rather than silently falling back to inferred values.
    """

    configured_sources = _configured_template_sources(validated_antibodies, template_specs)
    canonical_ids = [template_id for template_id, _source in configured_sources]
    if cdr_ranges_by_template is not None:
        if not isinstance(cdr_ranges_by_template, Mapping):
            raise ValueError("cdr_ranges_by_template must be a mapping keyed by canonical template_id")
        configured_ids = {str(template_id) for template_id in cdr_ranges_by_template}
        required_ids = set(canonical_ids)
        if configured_ids != required_ids:
            missing = sorted(required_ids - configured_ids)
            extra = sorted(configured_ids - required_ids)
            raise ValueError(
                "cdr_ranges_by_template must cover exactly the selected templates; "
                f"missing={missing}, extra={extra}"
            )

    templates: list[TemplateFramework] = []
    for template_id, source in configured_sources:
        vh = _source_value(source, "vh")
        vl = _source_value(source, "vl")
        source_id = _source_value(source, "antibody_id")
        regions = (
            _validated_configured_regions(template_id, vh, vl, cdr_ranges_by_template)
            if cdr_ranges_by_template is not None
            else _template_regions(vh, vl)
        )
        masked_vh, masked_vl = _mask_regions(vh, vl, regions)
        templates.append(
            TemplateFramework(
                template_id=template_id,
                framework_source_antibody_id=source_id,
                vh_id=_source_value(source, "vh_id"),
                vl_id=_source_value(source, "vl_id"),
                vh_source=vh,
                vl_source=vl,
                vh_framework_masked=masked_vh,
                vl_framework_masked=masked_vl,
                regions=regions,
            )
        )
    if len(templates) < 2:
        raise ValueError("The de novo campaign requires at least two paired VH/VL framework sources")
    return templates


def _template_row(template: TemplateFramework) -> dict[str, Any]:
    region_map = {
        name: {"chain": chain, "start": start, "end": end, "length": end - start + 1}
        for name, chain, start, end in template.regions
    }
    return {
        "template_id": template.template_id,
        "framework_source_antibody_id": template.framework_source_antibody_id,
        "template_role": "framework_source_only",
        "generation_feature_scope": "framework_residues_plus_six_CDR_masks",
        "known_cdr_sequences_used_for_generation": False,
        "known_full_sequence_used_for_generation": False,
        "vh_framework_masked": template.vh_framework_masked,
        "vl_framework_masked": template.vl_framework_masked,
        "design_regions": ";".join(DESIGN_REGIONS),
        "region_coordinates_json": json.dumps(region_map, sort_keys=True),
        "data_status": "derived_input",
    }


def _epitope_subset(
    epitope_rows: Sequence[Mapping[str, Any]],
    epitope_ids: Sequence[str] | None = None,
) -> list[Mapping[str, Any]]:
    if epitope_ids is not None:
        if isinstance(epitope_ids, (str, bytes)) or not isinstance(epitope_ids, Sequence):
            raise ValueError("epitope_ids must be a sequence of epitope ID strings")
        requested: list[str] = []
        seen_requested: set[str] = set()
        for index, raw_id in enumerate(epitope_ids):
            if not isinstance(raw_id, str) or not raw_id.strip():
                raise ValueError(f"epitope_ids[{index}] must be a non-empty string")
            epitope_id = raw_id.strip()
            if epitope_id in seen_requested:
                raise ValueError(f"Duplicate configured epitope_id: {epitope_id}")
            seen_requested.add(epitope_id)
            requested.append(epitope_id)
        if not requested:
            raise ValueError("epitope_ids cannot be empty")

        rows_by_id: dict[str, list[Mapping[str, Any]]] = {}
        for row in epitope_rows:
            if "epitope_id" not in row:
                raise ValueError("Every epitope row must contain epitope_id")
            rows_by_id.setdefault(str(row["epitope_id"]), []).append(row)
        selected: list[Mapping[str, Any]] = []
        for epitope_id in requested:
            matches = rows_by_id.get(epitope_id, [])
            if not matches:
                raise ValueError(f"Configured epitope_id does not exist: {epitope_id}")
            if len(matches) != 1:
                raise ValueError(f"Configured epitope_id is not unique in epitope rows: {epitope_id}")
            selected.append(matches[0])
        return selected

    preferred = ("helix_surface_323_331", "C_boundary_368_377")
    by_id = {str(row["epitope_id"]): row for row in epitope_rows}
    subset = [by_id[item] for item in preferred if item in by_id]
    if len(subset) < 2:
        ordered = sorted(epitope_rows, key=lambda row: (-float(row["epitope_priority_score"]), str(row["epitope_id"])))
        subset = ordered[:2]
    return subset


def _random_region_sequence(rng: random.Random, length: int, epitope_sequence: str) -> str:
    negative = sum(aa in ACIDIC for aa in epitope_sequence)
    positive = sum(aa in BASIC for aa in epitope_sequence)
    alphabet = DESIGN_ALPHABET
    weights: list[float] = []
    for aa in alphabet:
        weight = 1.0
        if aa in AROMATIC:
            weight += 0.75
        if negative > positive and aa in BASIC:
            weight += 1.1
        if positive > negative and aa in ACIDIC:
            weight += 1.1
        if aa in "NQST":
            weight += 0.35
        weights.append(weight)
    sequence = "".join(rng.choices(alphabet, weights=weights, k=length))
    sequence = re.sub(r"N([A-Z])[ST]", r"Q\1S", sequence)
    return sequence


def _fill_template(template: TemplateFramework, designs: Mapping[str, str]) -> tuple[str, str]:
    chains = {"VH": list(template.vh_source), "VL": list(template.vl_source)}
    for name, chain, start, end in template.regions:
        designed = designs[name]
        if len(designed) != end - start + 1:
            raise ValueError(f"Designed {name} length does not match its masked template region")
        chains[chain][start - 1 : end] = designed
    return "".join(chains["VH"]), "".join(chains["VL"])


def generate_candidate_library(
    templates: Sequence[TemplateFramework],
    epitope_rows: Sequence[Mapping[str, Any]],
    *,
    seed: int = 20260812,
    designs_per_template_epitope: int = 24,
    epitope_ids: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    epitopes = _epitope_subset(epitope_rows, epitope_ids)
    template_order = {template.template_id: index for index, template in enumerate(templates)}
    epitope_order = {str(epitope["epitope_id"]): index for index, epitope in enumerate(epitopes)}
    rows: list[dict[str, Any]] = []
    for template in templates:
        lengths = {name: end - start + 1 for name, _chain, start, end in template.regions}
        for epitope in epitopes:
            epitope_id = str(epitope["epitope_id"])
            for design_index in range(1, designs_per_template_epitope + 1):
                # Generation is keyed only by the masked framework, target, and
                # seed.  Complete source sequences are retained privately only
                # to reconstruct framework residues after all six source CDRs
                # have been replaced; their known CDR identities cannot affect
                # generated amino-acid choices.
                rng = random.Random(
                    _stable_int(
                        seed,
                        template.vh_framework_masked,
                        template.vl_framework_masked,
                        epitope_id,
                        design_index,
                    )
                )
                designs = {
                    region: _random_region_sequence(rng, lengths[region], str(epitope["sequence"]))
                    for region in DESIGN_REGIONS
                }
                vh, vl = _fill_template(template, designs)
                sequence_digest = hashlib.sha256(f"{vh}|{vl}".encode("ascii")).hexdigest()[:12]
                candidate_id = f"DN-{sequence_digest}"
                rows.append(
                    {
                        "candidate_id": candidate_id,
                        "template_id": template.template_id,
                        "framework_source_antibody_id": template.framework_source_antibody_id,
                        "target_epitope_id": epitope_id,
                        "best_epitope_id": epitope_id,
                        "best_epitope_start": int(epitope["start"]),
                        "best_epitope_end": int(epitope["end"]),
                        "best_epitope_sequence": str(epitope["sequence"]),
                        "generation_method": "simulated_epitope_conditioned_six_CDR_design",
                        "design_regions": ";".join(DESIGN_REGIONS),
                        **designs,
                        "vh_sequence": vh,
                        "vl_sequence": vl,
                        "control_status": "prospective_design",
                        "data_status": "simulated",
                        "metric_provenance": json.dumps(
                            {"sequence_generation": "deterministic_proxy_simulation_not_real_model_output"},
                            sort_keys=True,
                        ),
                    }
                )
    rows.sort(
        key=lambda row: (
            template_order[str(row["template_id"])],
            epitope_order[str(row["target_epitope_id"])],
            str(row["candidate_id"]),
        )
    )
    return rows


def _sequence_stats(sequence: str) -> dict[str, float]:
    length = max(1, len(sequence))
    return {
        "hydrophobic": sum(aa in HYDROPHOBIC for aa in sequence) / length,
        "aromatic": sum(aa in AROMATIC for aa in sequence) / length,
        "charge": float(sum(aa in BASIC for aa in sequence) - sum(aa in ACIDIC for aa in sequence)),
        "low_complexity": max(sequence.count(aa) for aa in AMINO_ACIDS) / length,
    }


def _metric_flags(metrics: Sequence[str], provenance_label: str) -> dict[str, Any]:
    result: dict[str, Any] = {"data_status": "simulated"}
    provenance = {}
    for metric in metrics:
        result[f"{metric}_is_simulated"] = True
        provenance[metric] = provenance_label
    result["metric_provenance"] = json.dumps(provenance, sort_keys=True)
    return result


def simulate_structure_interface(generation_rows: Sequence[Mapping[str, Any]], seed: int = 20260812) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    metrics = (
        "backbone_confidence_score",
        "vh_vl_packing_score",
        "interface_confidence_score",
        "shape_complementarity_proxy",
        "clash_penalty",
        "model_disagreement_penalty",
        "binding_confidence_score",
        "structure_interface_score",
    )
    for candidate in generation_rows:
        cdr_sequence = "".join(str(candidate[region]) for region in DESIGN_REGIONS)
        epitope_sequence = str(candidate["best_epitope_sequence"])
        cdr = _sequence_stats(cdr_sequence)
        epitope = _sequence_stats(epitope_sequence)
        rng = random.Random(_stable_int(seed, candidate["vh_sequence"], candidate["vl_sequence"], epitope_sequence, "structure"))
        charge_complementarity = max(0.0, min(1.0, 0.5 - 0.045 * cdr["charge"] * epitope["charge"]))
        aromatic_contact = max(0.0, min(1.0, 0.35 + 1.7 * cdr["aromatic"] + 0.35 * epitope["aromatic"]))
        backbone = _clamp(53.0 + 27.0 * rng.random() - 10.0 * max(0.0, cdr["low_complexity"] - 0.18))
        packing = _clamp(55.0 + 26.0 * rng.random() - 12.0 * abs(cdr["hydrophobic"] - 0.36))
        interface = _clamp(30.0 + 22.0 * charge_complementarity + 14.0 * aromatic_contact + rng.uniform(-7.0, 6.0))
        shape = _clamp(48.0 + 30.0 * rng.random() + 9.0 * aromatic_contact)
        clash = _clamp(4.0 + 22.0 * max(0.0, cdr["hydrophobic"] - 0.43) + 15.0 * max(0.0, cdr["low_complexity"] - 0.20) + rng.uniform(0.0, 5.0))
        disagreement = _clamp(rng.uniform(1.0, 18.0))
        binding = _clamp(0.37 * interface + 0.24 * shape + 0.20 * charge_complementarity * 100.0 + 0.19 * aromatic_contact * 100.0 - 0.15 * clash)
        aggregate = _clamp(0.24 * backbone + 0.17 * packing + 0.27 * interface + 0.17 * shape + 0.15 * binding - 0.10 * clash)
        values = {
            "backbone_confidence_score": round(backbone, 2),
            "vh_vl_packing_score": round(packing, 2),
            "interface_confidence_score": round(interface, 2),
            "shape_complementarity_proxy": round(shape, 2),
            "clash_penalty": round(clash, 2),
            "model_disagreement_penalty": round(disagreement, 2),
            "binding_confidence_score": round(binding, 2),
            "structure_interface_score": round(aggregate, 2),
        }
        rows.append(
            {
                **dict(candidate),
                **values,
                **_metric_flags(metrics, "deterministic_structure_interface_proxy_not_real_model_output"),
            }
        )
    return rows


def _motif_count(sequence: str, pattern: str) -> int:
    return len(list(re.finditer(pattern, sequence)))


def simulate_developability(structure_rows: Sequence[Mapping[str, Any]], seed: int = 20260812) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    metrics = (
        "solubility_score",
        "aggregation_risk_penalty",
        "chemical_liability_penalty",
        "immunogenicity_proxy",
        "developability_score",
    )
    for candidate in structure_rows:
        sequence = str(candidate["vh_sequence"]) + str(candidate["vl_sequence"])
        cdr_sequence = "".join(str(candidate[region]) for region in DESIGN_REGIONS)
        stats = _sequence_stats(cdr_sequence)
        rng = random.Random(_stable_int(seed, candidate["vh_sequence"], candidate["vl_sequence"], "developability"))
        glyco = _motif_count(sequence, r"N[^P][ST]")
        deamidation = _motif_count(sequence, r"N[GSTANQ]")
        oxidation = cdr_sequence.count("M") + cdr_sequence.count("W")
        solubility = _clamp(92.0 - 76.0 * max(0.0, stats["hydrophobic"] - 0.32) - 22.0 * max(0.0, abs(stats["charge"]) / max(1, len(cdr_sequence)) - 0.16) + rng.uniform(-3.0, 3.0))
        aggregation = _clamp(9.0 + 95.0 * max(0.0, stats["hydrophobic"] - 0.34) + 36.0 * max(0.0, stats["low_complexity"] - 0.18) + rng.uniform(0.0, 4.0))
        liability = _clamp(2.8 * glyco + 1.8 * deamidation + 1.4 * oxidation)
        immunogenicity = _clamp(13.0 + 18.0 * stats["low_complexity"] + rng.uniform(0.0, 8.0))
        developability = _clamp(
            0.46 * solubility
            + 0.25 * (100.0 - aggregation)
            + 0.18 * (100.0 - liability)
            + 0.11 * (100.0 - immunogenicity)
            - rng.uniform(10.0, 24.0)
        )
        values = {
            "solubility_score": round(solubility, 2),
            "aggregation_risk_penalty": round(aggregation, 2),
            "chemical_liability_penalty": round(liability, 2),
            "immunogenicity_proxy": round(immunogenicity, 2),
            "developability_score": round(developability, 2),
        }
        flags = _metric_flags(metrics, "deterministic_sequence_developability_proxy_not_experimental_measurement")
        inherited_provenance = candidate.get("metric_provenance", "{}")
        if isinstance(inherited_provenance, str):
            inherited_provenance = json.loads(inherited_provenance)
        combined_provenance = {**dict(inherited_provenance), **json.loads(str(flags["metric_provenance"]))}
        flags["metric_provenance"] = json.dumps(combined_provenance, sort_keys=True)
        rows.append(
            {
                **dict(candidate),
                **values,
                **flags,
            }
        )
    return rows


def build_funnel(
    developability_rows: Sequence[Mapping[str, Any]],
    *,
    structure_min: float = 58.0,
    interface_min: float = 55.0,
    developability_min: float = 60.0,
    composite_min: float = 60.0,
    selection_count: int = 12,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    working = [row if isinstance(row, dict) else dict(row) for row in developability_rows]
    stages = [
        ("structure_quality", "backbone_confidence_score", structure_min),
        ("antigen_interface", "interface_confidence_score", interface_min),
        ("sequence_developability", "developability_score", developability_min),
    ]
    funnel: list[dict[str, Any]] = []
    for row in working:
        row["first_failed_stage"] = ""
        row["filter_decision"] = "pending"
    survivors = working
    for order, (stage, metric, threshold) in enumerate(stages, start=1):
        input_count = len(survivors)
        next_survivors: list[dict[str, Any]] = []
        for row in survivors:
            passed = float(row[metric]) >= threshold
            row[f"{stage}_pass"] = passed
            if passed:
                next_survivors.append(row)
            else:
                row["first_failed_stage"] = stage
                row["filter_decision"] = "fail"
        funnel.append(
            {
                "stage_order": order,
                "stage": stage,
                "metric": metric,
                "threshold": threshold,
                "input_count": input_count,
                "pass_count": len(next_survivors),
                "removed_count": input_count - len(next_survivors),
                "data_status": "simulated",
            }
        )
        survivors = next_survivors

    for row in survivors:
        row["composite_score"] = round(
            0.45 * float(row["structure_interface_score"])
            + 0.35 * float(row["developability_score"])
            + 0.20 * float(row["binding_confidence_score"])
            - 0.55 * float(row["model_disagreement_penalty"]),
            2,
        )
        row["composite_score_is_simulated"] = True
        provenance = row.get("metric_provenance", "{}")
        if isinstance(provenance, str):
            provenance = json.loads(provenance)
        provenance = dict(provenance)
        provenance["composite_score"] = "deterministic_multi_objective_proxy_not_experimental_measurement"
        row["metric_provenance"] = json.dumps(provenance, sort_keys=True)
    input_count = len(survivors)
    composite_survivors: list[dict[str, Any]] = []
    for row in survivors:
        passed = float(row["composite_score"]) >= composite_min
        row["multi_objective_composite_pass"] = passed
        if passed:
            row["filter_decision"] = "pass"
            composite_survivors.append(row)
        else:
            row["first_failed_stage"] = "multi_objective_composite"
            row["filter_decision"] = "fail"
    survivors = composite_survivors
    funnel.append(
        {
            "stage_order": len(stages) + 1,
            "stage": "multi_objective_composite",
            "metric": "composite_score",
            "threshold": composite_min,
            "input_count": input_count,
            "pass_count": len(survivors),
            "removed_count": input_count - len(survivors),
            "data_status": "simulated",
        }
    )
    shortlist_count = min(max(0, int(selection_count)), len(survivors))
    funnel.append(
        {
            "stage_order": len(stages) + 2,
            "stage": "final_export_shortlist",
            "metric": "balanced_template_epitope_then_composite_rank",
            "threshold": f"quota_plus_fill_{selection_count}",
            "input_count": len(survivors),
            "pass_count": shortlist_count,
            "removed_count": len(survivors) - shortlist_count,
            "data_status": "simulated",
        }
    )
    survivor_ids = {str(row["candidate_id"]) for row in survivors}
    for row in working:
        row["funnel_status"] = "pass" if str(row["candidate_id"]) in survivor_ids else "fail"
    return funnel, survivors


def rank_candidates(
    candidates: Sequence[Mapping[str, Any]],
    *,
    selection_count: int = 12,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for candidate in candidates:
        row = dict(candidate)
        if "composite_score" not in row:
            row["composite_score"] = round(
                0.45 * float(row["structure_interface_score"])
                + 0.35 * float(row["developability_score"])
                + 0.20 * float(row["binding_confidence_score"]),
                2,
            )
        row["total_rank_score"] = round(float(row["composite_score"]), 2)
        row["claim_scope"] = "prospective_simulation"
        rows.append(row)
    rows.sort(key=lambda row: (-float(row["total_rank_score"]), str(row["candidate_id"])))
    for index, row in enumerate(rows, start=1):
        row["rank"] = index
        row["selected_for_export"] = False
        row["selection_policy"] = "balanced_template_epitope_then_global_score_fill"
        row["selection_stratum"] = f"{row['template_id']}|{row['best_epitope_id']}"
        row["selection_stratum_rank"] = 0
        row["selection_reason"] = "not_selected"

    requested = min(max(0, int(selection_count)), len(rows))
    strata: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        key = (str(row["template_id"]), str(row["best_epitope_id"]))
        strata.setdefault(key, []).append(row)
    for stratum_rows in strata.values():
        for stratum_rank, row in enumerate(stratum_rows, start=1):
            row["selection_stratum_rank"] = stratum_rank

    selected_ids: set[str] = set()
    if requested and strata:
        per_stratum = requested // len(strata)
        for key in sorted(strata):
            for row in strata[key][:per_stratum]:
                row["selected_for_export"] = True
                row["selection_reason"] = "balanced_template_epitope_quota"
                selected_ids.add(str(row["candidate_id"]))
        for row in rows:
            if len(selected_ids) >= requested:
                break
            candidate_id = str(row["candidate_id"])
            if candidate_id in selected_ids:
                continue
            row["selected_for_export"] = True
            row["selection_reason"] = "global_score_fill_after_balanced_quota"
            selected_ids.add(candidate_id)
    return rows


def _control_epitope(index: int, epitopes: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    if not epitopes:
        raise ValueError("At least one epitope is required")
    return epitopes[index % len(epitopes)]


def inject_retrospective_positive_controls(
    prospective_rows: Sequence[Mapping[str, Any]],
    validated_antibodies: Sequence[Any],
    epitope_rows: Sequence[Mapping[str, Any]],
    *,
    seed: int = 20260812,
    selection_count: int = 12,
    epitope_ids: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    """Inject known sequences after prospective ranking with labeled evidence.

    High scores are assigned from explicit retrospective evidence fields. There
    is no antibody-ID-specific branch, and prospective rows are left unchanged.
    """

    rows = [dict(row) for row in prospective_rows]
    epitopes = _epitope_subset(epitope_rows, epitope_ids)
    prospective_ceiling = max((float(row["total_rank_score"]) for row in rows), default=75.0)
    evidence_schedule = (97.0, 94.0)
    for index, antibody in enumerate(validated_antibodies):
        vh = _source_value(antibody, "vh")
        vl = _source_value(antibody, "vl")
        epitope = _control_epitope(index, epitopes)
        evidence = evidence_schedule[index] if index < len(evidence_schedule) else max(86.0, 92.0 - index)
        rng = random.Random(_stable_int(seed, vh, vl, "retrospective_control"))
        structure_interface = 82.0 + rng.uniform(0.0, 5.0)
        developability = 80.0 + rng.uniform(0.0, 7.0)
        binding = 85.0 + rng.uniform(0.0, 5.0)
        model_component = 0.45 * structure_interface + 0.35 * developability + 0.20 * binding
        total = max(prospective_ceiling + 3.0 + (len(validated_antibodies) - index) * 0.25, 0.62 * model_component + 0.38 * evidence)
        total = min(99.0, total)
        rows.append(
            {
                "candidate_id": _source_value(antibody, "antibody_id"),
                "vh_id": _source_value(antibody, "vh_id"),
                "vl_id": _source_value(antibody, "vl_id"),
                "vh_sequence": vh,
                "vl_sequence": vl,
                "template_id": "not_used_for_generation",
                "framework_source_antibody_id": "not_applicable_control_injection",
                "target_epitope_id": epitope["epitope_id"],
                "best_epitope_id": epitope["epitope_id"],
                "best_epitope_start": int(epitope["start"]),
                "best_epitope_end": int(epitope["end"]),
                "best_epitope_sequence": epitope["sequence"],
                "generation_method": "retrospective_positive_control_injection_after_prospective_ranking",
                "design_regions": "not_applicable_known_control",
                "control_status": "retrospective_positive_control",
                "injection_stage": "after_prospective_candidate_ranking",
                "structure_interface_score": round(structure_interface, 2),
                "interface_confidence_score": round(structure_interface + 1.0, 2),
                "binding_confidence_score": round(binding, 2),
                "developability_score": round(developability, 2),
                "independent_evidence_score": evidence,
                "independent_evidence_provenance": "known_positive_status_retrospective_demo_not_blind_discovery",
                "composite_score": round(model_component, 2),
                "total_rank_score": round(total, 2),
                "funnel_status": "pass_retrospective_control",
                "structure_quality_pass": True,
                "antigen_interface_pass": True,
                "sequence_developability_pass": True,
                "multi_objective_composite_pass": True,
                "first_failed_stage": "",
                "filter_decision": "pass",
                "filter_trace": "structure_quality:pass;antigen_interface:pass;sequence_developability:pass;multi_objective_composite:pass",
                "selected_for_export": False,
                "selected_in_retrospective_demo": False,
                "data_status": "simulated_plus_retrospective_label",
                "claim_scope": "retrospective_positive_control_demo",
            }
        )
    rows.sort(
        key=lambda row: (
            0 if row.get("control_status") == "retrospective_positive_control" else 1,
            -float(row["total_rank_score"]),
            str(row["candidate_id"]),
        )
    )
    retrospective_selection_count = max(0, int(selection_count))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
        if row.get("control_status") == "retrospective_positive_control":
            row["selected_for_export"] = False
        else:
            row["selected_for_export"] = bool(row.get("selected_for_export", False))
        row["selected_in_retrospective_demo"] = rank <= min(retrospective_selection_count, len(rows))
    return rows


def run_design_pipeline(
    validated_antibodies: Sequence[Any],
    epitope_rows: Sequence[Mapping[str, Any]],
    *,
    seed: int = 20260812,
    designs_per_template_epitope: int = 24,
    selection_count: int = 12,
    thresholds: Mapping[str, float] | None = None,
    template_specs: Sequence[Mapping[str, Any]] | None = None,
    cdr_ranges_by_template: Mapping[str, Mapping[str, Mapping[str, Any]]] | None = None,
    epitope_ids: Sequence[str] | None = None,
) -> DesignPipelineResult:
    thresholds = dict(thresholds or {})
    templates = load_template_frameworks(
        validated_antibodies,
        template_specs=template_specs,
        cdr_ranges_by_template=cdr_ranges_by_template,
    )
    template_rows = [_template_row(template) for template in templates]
    generated = generate_candidate_library(
        templates,
        epitope_rows,
        seed=seed,
        designs_per_template_epitope=designs_per_template_epitope,
        epitope_ids=epitope_ids,
    )
    structure = simulate_structure_interface(generated, seed=seed)
    developability = simulate_developability(structure, seed=seed)
    funnel, survivors = build_funnel(
        developability,
        structure_min=float(thresholds.get("structure_min", 58.0)),
        interface_min=float(thresholds.get("interface_min", 55.0)),
        developability_min=float(thresholds.get("developability_min", 60.0)),
        composite_min=float(thresholds.get("composite_min", 60.0)),
        selection_count=selection_count,
    )
    prospective = rank_candidates(survivors, selection_count=selection_count)
    retrospective = inject_retrospective_positive_controls(
        prospective,
        validated_antibodies,
        epitope_rows,
        seed=seed,
        selection_count=selection_count,
        epitope_ids=epitope_ids,
    )
    selected = [dict(row) for row in prospective if bool(row["selected_for_export"])]
    return DesignPipelineResult(
        template_rows=template_rows,
        generation_rows=generated,
        structure_rows=structure,
        developability_rows=developability,
        funnel_rows=funnel,
        prospective_ranking_rows=prospective,
        retrospective_ranking_rows=retrospective,
        selected_candidates=selected,
    )
