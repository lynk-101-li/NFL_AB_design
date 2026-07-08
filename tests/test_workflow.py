import unittest
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from nfl_ab_design.workflow import run_workflow


class WorkflowReplayTest(unittest.TestCase):
    def test_validation_antibodies_rank_top_two(self) -> None:
        result = run_workflow()
        ranking_rows = result["ranking_rows"]

        self.assertEqual(result["primary_fragment"], "280-375")
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


if __name__ == "__main__":
    unittest.main()
