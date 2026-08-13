import json
from pathlib import Path
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
import assemble_review_candidate_manifest as assemble  # noqa: E402
import prepare_real_model_jobs as real_jobs  # noqa: E402


class ReviewCandidateManifestTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.candidate = (
            PROJECT_ROOT
            / "input/structures/target_structure_manifest.candidate.blocked.json"
        )
        cls.fragment = (
            PROJECT_ROOT
            / "input/template_structures/antibody_template_manifest_fragment.blocked.json"
        )
        cls.sasa = (
            PROJECT_ROOT
            / "input/structures/NEFL_P07196_AFDB_v6_280-377_sasa_review.json"
        )

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    @staticmethod
    def _copy_json(source: Path, destination: Path) -> dict:
        value = json.loads(source.read_text(encoding="utf-8"))
        destination.write_text(
            json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        return value

    @staticmethod
    def _write_json(path: Path, value: dict) -> None:
        path.write_text(
            json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )

    def _assemble(self, **overrides):
        arguments = {
            "candidate_manifest": self.candidate,
            "template_fragment": self.fragment,
            "sasa_review": self.sasa,
            "output": self.root / "review_candidate.blocked.json",
        }
        arguments.update(overrides)
        return assemble.assemble_review_candidate_manifest(**arguments)

    def test_assembles_hash_bound_candidate_but_keeps_all_review_gates_closed(self):
        summary = self._assemble()
        output = self.root / "review_candidate.blocked.json"
        manifest = json.loads(output.read_text(encoding="utf-8"))

        self.assertEqual(summary["execution_state"], assemble.BLOCKED_STATE)
        self.assertEqual(manifest["execution_state"], assemble.BLOCKED_STATE)
        self.assertFalse(manifest["real_model_handoff_authorized"])
        self.assertEqual(
            manifest["review"],
            {
                "reviewed_by": "",
                "reviewed_at": "",
                "contracts_acknowledged": False,
            },
        )
        self.assertEqual(
            manifest["selected_hotspots_by_epitope"], assemble.EXPECTED_HOTSPOTS
        )
        self.assertEqual(
            set(manifest["framework_coordinate_inputs"]["rfantibody_hlt_pdbs"]),
            set(assemble.TEMPLATE_IDS),
        )
        self.assertEqual(
            set(manifest["framework_coordinate_inputs"]["germinal_scfv_pdbs"]),
            set(assemble.TEMPLATE_IDS),
        )
        self.assertIn("sasa_review", manifest["review_evidence"])
        self.assertEqual(len(manifest["promotion_instructions"]), 4)

        # Exercise the exact gate used by prepare_real_model_jobs: populated
        # recommendations and coordinates still cannot bypass human review.
        with self.assertRaisesRegex(
            real_jobs.RealModelHandoffError, "reviewed_ready_for_handoff"
        ):
            real_jobs._validate_target_manifest(manifest, request={"epitopes": []})

    def test_overwrite_is_explicit_and_deterministic(self):
        first = self._assemble()
        output = self.root / "review_candidate.blocked.json"
        first_bytes = output.read_bytes()
        with self.assertRaisesRegex(
            assemble.ReviewCandidateAssemblyError, "Refusing to overwrite"
        ):
            self._assemble()
        second = self._assemble(overwrite=True)
        self.assertEqual(first_bytes, output.read_bytes())
        self.assertEqual(first["output_sha256"], second["output_sha256"])

    def test_rejects_candidate_hash_tampering(self):
        candidate_copy = self.root / "candidate.json"
        value = self._copy_json(self.candidate, candidate_copy)
        value["candidate_provenance"]["target_pdb_sha256"] = "0" * 64
        self._write_json(candidate_copy, value)
        with self.assertRaisesRegex(
            assemble.ReviewCandidateAssemblyError, "hash mismatch"
        ):
            self._assemble(candidate_manifest=candidate_copy)

    def test_rejects_sasa_binding_or_recommendation_tampering(self):
        sasa_copy = self.root / "sasa.json"
        value = self._copy_json(self.sasa, sasa_copy)
        value["inputs"]["cropped_design_target"]["sha256"] = "0" * 64
        self._write_json(sasa_copy, value)
        with self.assertRaisesRegex(
            assemble.ReviewCandidateAssemblyError, "hash mismatch"
        ):
            self._assemble(sasa_review=sasa_copy)

        value = self._copy_json(self.sasa, sasa_copy)
        value["epitope_reviews"][0][
            "recommended_hotspots_pending_human_review"
        ] = [325, 330]
        self._write_json(sasa_copy, value)
        with self.assertRaisesRegex(
            assemble.ReviewCandidateAssemblyError, "must be exactly"
        ):
            self._assemble(sasa_review=sasa_copy)

    def test_rejects_missing_or_nonexact_template_set(self):
        fragment_copy = self.root / "fragment.json"
        value = self._copy_json(self.fragment, fragment_copy)
        value["framework_coordinate_inputs"]["rfantibody_hlt_pdbs"][
            "template_7-H11-D3-2-C7"
        ] = "input/template_structures/does-not-exist.pdb"
        self._write_json(fragment_copy, value)
        with self.assertRaisesRegex(
            assemble.ReviewCandidateAssemblyError, "does not exist"
        ):
            self._assemble(template_fragment=fragment_copy)

        value = self._copy_json(self.fragment, fragment_copy)
        value["framework_coordinate_inputs"]["germinal_scfv_pdbs"][
            "unexpected-template"
        ] = value["framework_coordinate_inputs"]["germinal_scfv_pdbs"][
            "template_7-H11-D3-2-C7"
        ]
        self._write_json(fragment_copy, value)
        with self.assertRaisesRegex(
            assemble.ReviewCandidateAssemblyError, "must contain exactly"
        ):
            self._assemble(template_fragment=fragment_copy)

    def test_never_writes_formal_manifest_even_with_overwrite(self):
        with self.assertRaisesRegex(
            assemble.ReviewCandidateAssemblyError, "can never create or overwrite"
        ):
            self._assemble(output=assemble.FORMAL_MANIFEST, overwrite=True)


if __name__ == "__main__":
    unittest.main()
