import json
import sys
import unittest
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from nfl_ab_design.design_pipeline import run_design_pipeline
from nfl_ab_design.workflow import (
    ANTIBODY_FASTA_PATH,
    GENPEPT_PATH,
    build_epitope_windows,
    load_antibodies,
    parse_genpept_sequence,
)


CDR_REGIONS = {"H1", "H2", "H3", "L1", "L2", "L3"}
EXPECTED_CHOTHIA_RAW_RANGES = {
    "template_7-H11-D3-2-C7": {
        "H1": ("VH", 26, 32),
        "H2": ("VH", 52, 57),
        "H3": ("VH", 99, 101),
        "L1": ("VL", 24, 34),
        "L2": ("VL", 50, 56),
        "L3": ("VL", 89, 97),
    },
    "template_15-C12-H6": {
        "H1": ("VH", 26, 33),
        "H2": ("VH", 53, 57),
        "H3": ("VH", 99, 109),
        "L1": ("VL", 24, 34),
        "L2": ("VL", 50, 56),
        "L3": ("VL", 89, 97),
    },
}


def _regions(value: Any) -> set[str]:
    if isinstance(value, str):
        return {item.strip() for item in value.replace(",", ";").split(";") if item.strip()}
    if isinstance(value, Iterable):
        return {str(item) for item in value}
    return set()


def _candidate_ids(rows: list[Any]) -> list[str]:
    ids: list[str] = []
    for row in rows:
        if isinstance(row, str):
            ids.append(row)
        else:
            ids.append(str(row["candidate_id"]))
    return ids


def _score_fields(row: dict[str, Any]) -> dict[str, float]:
    """Return score-bearing fields while intentionally excluding rank and identity."""

    suffixes = ("_score", "_proxy", "_penalty")
    return {
        key: float(value)
        for key, value in row.items()
        if not isinstance(value, bool)
        and isinstance(value, (int, float))
        and key.endswith(suffixes)
    }


class DesignPipelineContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.validated = load_antibodies(ANTIBODY_FASTA_PATH)
        full_sequence = parse_genpept_sequence(GENPEPT_PATH)
        cls.epitope_rows = build_epitope_windows(full_sequence, "280-375")
        cls.seed = 20260812
        cls.designs_per_template_epitope = 3
        cls.result = run_design_pipeline(
            cls.validated,
            cls.epitope_rows,
            seed=cls.seed,
            designs_per_template_epitope=cls.designs_per_template_epitope,
        ).as_dict()

    def test_two_framework_templates_have_distinct_sources_and_all_six_cdrs_are_design_regions(self) -> None:
        templates = self.result["template_rows"]

        self.assertGreaterEqual(len(templates), 2)
        self.assertGreaterEqual(
            len({row["framework_source_antibody_id"] for row in templates}),
            2,
            "The two framework templates must not silently come from the same source antibody.",
        )
        for row in templates:
            self.assertEqual(row["template_role"], "framework_source_only")
            self.assertEqual(_regions(row["design_regions"]), CDR_REGIONS)
            self.assertIn("X", row["vh_framework_masked"])
            self.assertIn("X", row["vl_framework_masked"])

        for row in self.result["generation_rows"]:
            self.assertEqual(_regions(row["design_regions"]), CDR_REGIONS)
            for region in CDR_REGIONS:
                self.assertTrue(
                    row.get(region),
                    f"{row.get('candidate_id', '<unknown>')} has no independently recorded {region} design.",
                )

    def test_fixed_seed_reproduces_every_pipeline_table(self) -> None:
        replay = run_design_pipeline(
            self.validated,
            self.epitope_rows,
            seed=self.seed,
            designs_per_template_epitope=self.designs_per_template_epitope,
        ).as_dict()

        self.assertEqual(replay, self.result)

    def test_configured_anarci_chothia_masks_and_generation_change_only_cdrs(self) -> None:
        config = json.loads((PROJECT_ROOT / "config/design_campaign.json").read_text(encoding="utf-8"))
        configured_ranges = config["cdr_ranges"]["by_template"]
        self.assertEqual(config["cdr_ranges"]["annotation_method"], "ANARCI 2020.04.23 Chothia")
        self.assertEqual(
            {
                template_id: {
                    name: (spec["chain"], spec["start"], spec["end"])
                    for name, spec in regions.items()
                }
                for template_id, regions in configured_ranges.items()
            },
            EXPECTED_CHOTHIA_RAW_RANGES,
        )

        result = run_design_pipeline(
            self.validated,
            self.epitope_rows,
            seed=self.seed,
            designs_per_template_epitope=2,
            template_specs=config["templates"],
            cdr_ranges_by_template=configured_ranges,
            epitope_ids=[row["id"] for row in config["target_epitopes"]],
        ).as_dict()
        sources = {antibody.antibody_id: antibody for antibody in self.validated}
        template_rows = {row["template_id"]: row for row in result["template_rows"]}

        for template_id, expected_regions in EXPECTED_CHOTHIA_RAW_RANGES.items():
            template_row = template_rows[template_id]
            source = sources[template_row["framework_source_antibody_id"]]
            coordinates = json.loads(template_row["region_coordinates_json"])
            masked_by_chain = {
                "VH": template_row["vh_framework_masked"],
                "VL": template_row["vl_framework_masked"],
            }
            source_by_chain = {"VH": source.vh, "VL": source.vl}
            masked_positions = {"VH": set(), "VL": set()}
            for name, (chain, start, end) in expected_regions.items():
                self.assertEqual(
                    coordinates[name],
                    {"chain": chain, "start": start, "end": end, "length": end - start + 1},
                )
                masked_positions[chain].update(range(start, end + 1))
            for chain in ("VH", "VL"):
                for position, (masked_aa, source_aa) in enumerate(
                    zip(masked_by_chain[chain], source_by_chain[chain], strict=True), start=1
                ):
                    self.assertEqual(masked_aa, "X" if position in masked_positions[chain] else source_aa)

        for candidate in result["generation_rows"]:
            source = sources[candidate["framework_source_antibody_id"]]
            source_by_chain = {"VH": source.vh, "VL": source.vl}
            candidate_by_chain = {"VH": candidate["vh_sequence"], "VL": candidate["vl_sequence"]}
            regions = EXPECTED_CHOTHIA_RAW_RANGES[candidate["template_id"]]
            designed_positions = {"VH": set(), "VL": set()}
            for name, (chain, start, end) in regions.items():
                designed_positions[chain].update(range(start, end + 1))
                self.assertEqual(candidate_by_chain[chain][start - 1 : end], candidate[name])
            for chain in ("VH", "VL"):
                for position, (candidate_aa, source_aa) in enumerate(
                    zip(candidate_by_chain[chain], source_by_chain[chain], strict=True), start=1
                ):
                    if position not in designed_positions[chain]:
                        self.assertEqual(candidate_aa, source_aa)

    def test_prospective_tables_contain_no_known_positive(self) -> None:
        validated_ids = {antibody.antibody_id for antibody in self.validated}
        validated_sequences = {(antibody.vh, antibody.vl) for antibody in self.validated}

        self.assertTrue(self.result["generation_rows"])
        self.assertTrue(self.result["prospective_ranking_rows"])
        for table_name in ("generation_rows", "prospective_ranking_rows"):
            for row in self.result[table_name]:
                self.assertEqual(row["control_status"], "prospective_design")
                self.assertNotIn(row["candidate_id"], validated_ids)
                self.assertNotIn((row["vh_sequence"], row["vl_sequence"]), validated_sequences)

    def test_retrospective_table_has_exactly_two_positive_controls_ranked_top_two(self) -> None:
        ranking = sorted(self.result["retrospective_ranking_rows"], key=lambda row: int(row["rank"]))
        controls = [row for row in ranking if row["control_status"] == "retrospective_positive_control"]

        self.assertEqual(len(controls), 2)
        self.assertEqual(
            {row["candidate_id"] for row in controls},
            {antibody.antibody_id for antibody in self.validated},
        )
        self.assertEqual(
            [row["control_status"] for row in ranking[:2]],
            ["retrospective_positive_control", "retrospective_positive_control"],
        )
        for row in controls:
            self.assertGreater(float(row["independent_evidence_score"]), 0.0)
            self.assertTrue(row["independent_evidence_provenance"])
            self.assertIs(row["selected_for_export"], False)
            self.assertTrue(row["selected_in_retrospective_demo"])

    def test_every_simulated_metric_has_machine_readable_provenance(self) -> None:
        for table_name in ("structure_rows", "developability_rows"):
            rows = self.result[table_name]
            self.assertTrue(rows, f"{table_name} unexpectedly empty")
            for row in rows:
                self.assertEqual(row["data_status"], "simulated")
                provenance = row["metric_provenance"]
                if isinstance(provenance, str):
                    provenance = json.loads(provenance)
                self.assertIsInstance(provenance, dict)

                simulated_metrics = {
                    key[: -len("_is_simulated")]
                    for key, value in row.items()
                    if key.endswith("_is_simulated") and value is True
                }
                self.assertTrue(simulated_metrics, f"{table_name} row has no explicitly simulated metrics")
                for metric in simulated_metrics:
                    self.assertIn(metric, row)
                    self.assertIn(metric, provenance)
                    self.assertTrue(str(provenance[metric]).strip())

                for metric in _score_fields(row):
                    self.assertIs(
                        row.get(f"{metric}_is_simulated"),
                        True,
                        f"{metric} is a simulated score but is not marked as such.",
                    )
                    self.assertIn(metric, provenance)

    def test_candidate_id_and_input_order_do_not_change_content_scores(self) -> None:
        renamed_reversed = []
        for index, antibody in enumerate(reversed(self.validated), start=1):
            renamed_reversed.append(
                replace(
                    antibody,
                    antibody_id=f"anonymous-source-{index}",
                    vh_id=f"anonymous-vh-{index}",
                    vl_id=f"anonymous-vl-{index}",
                    parent_id=f"anonymous-source-{index}",
                )
            )

        perturbed_result = run_design_pipeline(
            renamed_reversed,
            list(reversed(self.epitope_rows)),
            seed=self.seed,
            designs_per_template_epitope=self.designs_per_template_epitope,
        ).as_dict()

        def scores_by_sequence(result: dict[str, Any], table_name: str) -> dict[tuple[str, str, str], dict[str, float]]:
            return {
                (row["vh_sequence"], row["vl_sequence"], row["best_epitope_id"]): _score_fields(row)
                for row in result[table_name]
            }

        self.assertEqual(
            scores_by_sequence(self.result, "prospective_ranking_rows"),
            scores_by_sequence(perturbed_result, "prospective_ranking_rows"),
        )

    def test_funnel_counts_are_internally_consistent_and_monotone(self) -> None:
        funnel = sorted(self.result["funnel_rows"], key=lambda row: int(row["stage_order"]))
        self.assertGreaterEqual(len(funnel), 2)

        pass_counts: list[int] = []
        for row in funnel:
            input_count = int(row["input_count"])
            pass_count = int(row["pass_count"])
            removed_count = int(row["removed_count"])
            self.assertGreaterEqual(input_count, pass_count)
            self.assertEqual(removed_count, input_count - pass_count)
            pass_counts.append(pass_count)
        self.assertEqual(pass_counts, sorted(pass_counts, reverse=True))
        self.assertTrue(
            any(int(row["removed_count"]) > 0 for row in funnel[:-1]),
            "At least one scientific screening stage should remove simulated candidates.",
        )
        self.assertEqual(funnel[-1]["stage"], "final_export_shortlist")
        self.assertLessEqual(int(funnel[-1]["pass_count"]), 12)

    def test_export_selection_is_balanced_across_template_epitope_strata(self) -> None:
        ranking = sorted(self.result["prospective_ranking_rows"], key=lambda row: int(row["rank"]))
        selected_ids = _candidate_ids(self.result["selected_candidates"])
        flagged_ids = [row["candidate_id"] for row in ranking if row["selected_for_export"]]

        self.assertTrue(selected_ids)
        self.assertEqual(selected_ids, flagged_ids)
        selected_rows = [row for row in ranking if row["selected_for_export"]]
        available_strata = {
            (row["template_id"], row["best_epitope_id"])
            for row in ranking
        }
        selected_strata = {
            (row["template_id"], row["best_epitope_id"])
            for row in selected_rows
        }
        self.assertEqual(selected_strata, available_strata)
        self.assertTrue(
            all(
                row["selection_policy"]
                == "balanced_template_epitope_then_global_score_fill"
                for row in selected_rows
            )
        )


if __name__ == "__main__":
    unittest.main()
