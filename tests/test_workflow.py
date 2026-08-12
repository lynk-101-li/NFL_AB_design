import unittest
from pathlib import Path
import sys
import json
import tempfile


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from nfl_ab_design.workflow import (
    ANTIBODY_FASTA_PATH,
    OUTPUT_DIR,
    load_antibodies,
    load_external_pipeline_config,
    run_workflow,
)


class WorkflowCampaignTest(unittest.TestCase):
    def test_validation_antibodies_rank_top_two(self) -> None:
        result = run_workflow()
        ranking_rows = result["retrospective_ranking_rows"]

        self.assertEqual(result["primary_fragment"], "280-375")
        self.assertEqual(result["modeling_fragment"], "280-377")
        self.assertEqual(result["ranking_rows_scope"], "prospective_simulation")
        self.assertTrue(
            all(row["control_status"] == "prospective_design" for row in result["ranking_rows"])
        )
        self.assertGreaterEqual(len(ranking_rows), 2)
        self.assertEqual(ranking_rows[0]["candidate_id"], "7-H11-D3-2-C7")
        self.assertEqual(ranking_rows[1]["candidate_id"], "15-C12-H6")

    def test_sandwich_pair_is_recommended(self) -> None:
        result = run_workflow()
        sandwich = result["sandwich"]

        self.assertEqual(sandwich["antibody_1"], "7-H11-D3-2-C7")
        self.assertEqual(sandwich["antibody_2"], "15-C12-H6")
        self.assertEqual(sandwich["recommended_capture"], "7-H11-D3-2-C7")
        self.assertEqual(sandwich["recommended_detection"], "15-C12-H6")
        self.assertLessEqual(float(sandwich["epitope_overlap_ratio"]), 0.25)
        self.assertGreaterEqual(float(sandwich["sandwich_compatibility_score"]), 80.0)

    def test_structure_candidate_export_is_prospective_only(self) -> None:
        result = run_workflow()
        prospective = result["prospective_ranking_rows"]
        expected_ids = {
            row["candidate_id"]
            for row in prospective
            if row["selected_for_export"]
        }
        candidate_fasta = PROJECT_ROOT / result["manifest"]["exports"]["candidate_fv_chains_fasta"]
        exported_ids = {
            line[1:].split("|", 1)[0]
            for line in candidate_fasta.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        }

        self.assertEqual(exported_ids, expected_ids)
        self.assertNotIn("7-H11-D3-2-C7", exported_ids)
        self.assertNotIn("15-C12-H6", exported_ids)
        self.assertEqual(
            result["manifest"]["exports"]["sandwich_export_scope"],
            "retrospective_demo",
        )
        self.assertEqual(
            result["manifest"]["exports"]["primary_antigen_fragment"],
            "280-377",
        )
        for fasta_path in result["manifest"]["exports"]["complex_fastas"]:
            fasta_text = (PROJECT_ROOT / fasta_path).read_text(encoding="utf-8")
            self.assertIn(">A|NEFL_280-377", fasta_text)

    def test_campaign_config_drives_canonical_template_ids_and_epitopes(self) -> None:
        result = run_workflow()
        configured = json.loads((PROJECT_ROOT / "config/design_campaign.json").read_text())
        expected_template_ids = [row["template_id"] for row in configured["templates"]]
        expected_epitope_ids = [row["id"] for row in configured["target_epitopes"]]

        self.assertEqual(
            [row["template_id"] for row in result["design_result"]["template_rows"]],
            expected_template_ids,
        )
        self.assertEqual(
            sorted({row["target_epitope_id"] for row in result["design_result"]["generation_rows"]}),
            sorted(expected_epitope_ids),
        )
        for request_path in result["manifest"]["exports"]["design_request_files"]:
            request = json.loads((PROJECT_ROOT / request_path).read_text())
            self.assertEqual(
                [row["template_id"] for row in request["templates"]],
                expected_template_ids,
            )
            self.assertEqual(
                request["cdr_annotation"]["annotation_method"],
                "ANARCI 2020.04.23 Chothia",
            )

        iggm_request = json.loads(
            (PROJECT_ROOT / "outputs/exports/design_requests/iggm_design_request.json").read_text(
                encoding="utf-8"
            )
        )
        source_by_id = {antibody.antibody_id: antibody for antibody in load_antibodies(ANTIBODY_FASTA_PATH)}
        configured_ranges = configured["cdr_ranges"]["by_template"]
        for template in iggm_request["templates"]:
            source = source_by_id[template["framework_source_id"]]
            source_by_chain = {"VH": source.vh, "VL": source.vl}
            masked_by_chain = {"VH": template["masked_vh"], "VL": template["masked_vl"]}
            expected_regions = configured_ranges[template["template_id"]]
            self.assertEqual(
                {
                    row["region"]: (
                        row["chain"],
                        row["start_1_based"],
                        row["end_1_based_inclusive"],
                    )
                    for row in template["design_regions"]
                },
                {
                    name: (spec["chain"], spec["start"], spec["end"])
                    for name, spec in expected_regions.items()
                },
            )
            masked_positions = {"VH": set(), "VL": set()}
            for spec in expected_regions.values():
                masked_positions[spec["chain"]].update(range(spec["start"], spec["end"] + 1))
            for chain in ("VH", "VL"):
                for position, (masked_aa, source_aa) in enumerate(
                    zip(masked_by_chain[chain], source_by_chain[chain], strict=True), start=1
                ):
                    self.assertEqual(masked_aa, "X" if position in masked_positions[chain] else source_aa)

    def test_proxy_rerun_preserves_misplaced_external_results(self) -> None:
        sentinel = OUTPUT_DIR / "external_results" / "preserve-contract.txt"
        export_sentinel = OUTPUT_DIR / "exports" / "reviewed-human-note.txt"
        sentinel.parent.mkdir(parents=True, exist_ok=True)
        export_sentinel.parent.mkdir(parents=True, exist_ok=True)
        sentinel.write_text("must survive proxy regeneration\n", encoding="utf-8")
        export_sentinel.write_text("unknown export artifact must survive\n", encoding="utf-8")
        self.addCleanup(lambda: sentinel.unlink(missing_ok=True))
        self.addCleanup(lambda: export_sentinel.unlink(missing_ok=True))

        run_workflow()

        self.assertEqual(sentinel.read_text(encoding="utf-8"), "must survive proxy regeneration\n")
        self.assertEqual(
            export_sentinel.read_text(encoding="utf-8"),
            "unknown export artifact must survive\n",
        )

    def test_explicit_missing_configs_fail_closed(self) -> None:
        with self.assertRaises(FileNotFoundError):
            run_workflow(design_config_path=PROJECT_ROOT / "config/missing-campaign.json")
        with self.assertRaises(FileNotFoundError):
            run_workflow(external_config_path=PROJECT_ROOT / "config/missing-external.json")

    def test_external_enabled_must_be_boolean(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "external.json"
            path.write_text(
                json.dumps(
                    {
                        "schema": "nfl_ab_design.external_pipelines.v1",
                        "pipelines": [
                            {
                                "name": "unsafe-string-false",
                                "enabled": "false",
                                "input_selector": "candidate_fv_chains_fasta",
                                "command_template": "tool {input}",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "enabled must be boolean"):
                load_external_pipeline_config(path)


if __name__ == "__main__":
    unittest.main()
