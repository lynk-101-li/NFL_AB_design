import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import prepare_real_model_jobs  # noqa: E402
import prepare_target_structure as target_prep  # noqa: E402


class TargetStructurePreparationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source_pdb = (
            PROJECT_ROOT / "input" / "structures" / target_prep.SOURCE_FILENAME
        )
        if not cls.source_pdb.is_file():
            raise RuntimeError(f"Pinned test source is missing: {cls.source_pdb}")

    def _temporary_output(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        return temporary, Path(temporary.name) / "bundle"

    def test_pinned_source_builds_complete_review_blocked_bundle(self) -> None:
        _, output_dir = self._temporary_output()
        summary = target_prep.prepare_target_structure(self.source_pdb, output_dir)

        self.assertEqual(
            summary["execution_state"],
            "prepared_blocked_pending_human_review",
        )
        self.assertIs(summary["real_model_handoff_authorized"], False)
        self.assertEqual(summary["source_sha256"], target_prep.SOURCE_SHA256)
        self.assertEqual(summary["residue_count"], 98)

        source_copy = output_dir / target_prep.SOURCE_FILENAME
        target_pdb = output_dir / target_prep.TARGET_PDB_FILENAME
        target_fasta = output_dir / target_prep.TARGET_FASTA_FILENAME
        evidence_path = output_dir / target_prep.EVIDENCE_FILENAME
        manifest_path = output_dir / target_prep.CANDIDATE_MANIFEST_FILENAME
        for path in (source_copy, target_pdb, target_fasta, evidence_path, manifest_path):
            self.assertTrue(path.is_file(), path)

        self.assertEqual(
            hashlib.sha256(source_copy.read_bytes()).hexdigest(),
            target_prep.SOURCE_SHA256,
        )
        fasta_lines = target_fasta.read_text(encoding="ascii").splitlines()
        self.assertIn("theoretical_prediction_candidate", fasta_lines[0])
        self.assertEqual("".join(fasta_lines[1:]), target_prep.TARGET_SEQUENCE)

        pdb_lines = target_pdb.read_text(encoding="ascii").splitlines()
        atom_lines = [line for line in pdb_lines if line.startswith("ATOM")]
        observed_residue_numbers = sorted({int(line[22:26]) for line in atom_lines})
        self.assertEqual(observed_residue_numbers, list(range(280, 378)))
        self.assertIn(
            "THEORETICAL PREDICTION; NOT AN EXPERIMENTALLY DETERMINED STRUCTURE",
            "\n".join(pdb_lines),
        )
        self.assertTrue(
            any(
                line[12:16].strip() == "SG"
                and int(line[22:26]) == 322
                and line[21:22] == "A"
                for line in atom_lines
            )
        )

        evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
        self.assertEqual(evidence["schema"], target_prep.EVIDENCE_SCHEMA)
        self.assertEqual(
            evidence["scientific_status"],
            "theoretical_prediction_not_experimentally_determined",
        )
        self.assertTrue(all(evidence["validation"].values()))
        self.assertEqual(
            evidence["extraction"]["antigen_sequence_in_pdb_order"],
            target_prep.TARGET_SEQUENCE,
        )
        self.assertEqual(
            evidence["coordinate_mappings"]["full_coordinate_to_local_1_based"]["280"],
            1,
        )
        self.assertEqual(
            evidence["coordinate_mappings"]["full_coordinate_to_local_1_based"]["377"],
            98,
        )
        self.assertEqual(
            len(evidence["coordinate_mappings"]["full_coordinate_to_pdb"]), 98
        )
        self.assertEqual(evidence["alphafold_plddt"]["target_statistics"]["count"], 98)
        self.assertEqual(evidence["alphafold_plddt"]["target_statistics"]["mean"], 95.15)
        self.assertIs(evidence["exposure_proxy_definition"]["not_sasa"], True)

        hotspot_rows = {
            row["full_coordinate_1_based"]: row
            for epitope in evidence["hotspot_proposals"]
            for row in epitope["hotspot_evidence"]
        }
        self.assertEqual(set(hotspot_rows), {325, 329, 368, 372, 375})
        for row in hotspot_rows.values():
            self.assertTrue(row["atom_presence"]["backbone_N_CA_C_O_complete"])
            self.assertTrue(row["atom_presence"]["sidechain_probe_present"])
            self.assertGreater(
                row["geometry"][
                    "sidechain_probe_nearest_nonlocal_heavy_atom_distance_A"
                ],
                0.0,
            )
            self.assertIn(
                "nonlocal_ca_neighbor_count_within_10A",
                row["exposure_proxies_full_length_model_context"],
            )

        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema"], target_prep.TARGET_MANIFEST_SCHEMA)
        self.assertEqual(manifest["execution_state"], "blocked_pending_human_review")
        self.assertEqual(
            manifest["review"],
            {"reviewed_by": "", "reviewed_at": "", "contracts_acknowledged": False},
        )
        self.assertEqual(
            manifest["selected_hotspots_by_epitope"],
            {"helix_surface_323_331": [], "C_boundary_368_377": []},
        )
        self.assertEqual(
            manifest["proposed_hotspots_by_epitope_pending_human_review"],
            {
                "helix_surface_323_331": [325, 329],
                "C_boundary_368_377": [368, 372, 375],
            },
        )
        self.assertEqual(len(manifest["full_coordinate_to_pdb"]), 98)
        self.assertEqual(len(manifest["full_coordinate_to_local_1_based"]), 98)
        with self.assertRaisesRegex(
            prepare_real_model_jobs.RealModelHandoffError,
            "reviewed_ready_for_handoff",
        ):
            prepare_real_model_jobs._validate_target_manifest(manifest, request={})

    def test_regeneration_is_deterministic_and_overwrite_is_explicit(self) -> None:
        _, output_dir = self._temporary_output()
        target_prep.prepare_target_structure(self.source_pdb, output_dir)
        artifact_names = (
            target_prep.SOURCE_FILENAME,
            target_prep.TARGET_PDB_FILENAME,
            target_prep.TARGET_FASTA_FILENAME,
            target_prep.EVIDENCE_FILENAME,
            target_prep.CANDIDATE_MANIFEST_FILENAME,
        )
        before = {
            name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
            for name in artifact_names
        }
        with self.assertRaisesRegex(
            target_prep.TargetStructurePreparationError, "Refusing to overwrite"
        ):
            target_prep.prepare_target_structure(self.source_pdb, output_dir)

        target_prep.prepare_target_structure(
            self.source_pdb, output_dir, overwrite=True
        )
        after = {
            name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
            for name in artifact_names
        }
        self.assertEqual(after, before)

    def test_source_hash_mismatch_fails_before_creating_output(self) -> None:
        temporary, output_dir = self._temporary_output()
        corrupted = Path(temporary.name) / "corrupted_source.pdb"
        corrupted.write_bytes(self.source_pdb.read_bytes() + b"\n")
        with self.assertRaisesRegex(
            target_prep.TargetStructurePreparationError, "SHA-256 mismatch"
        ):
            target_prep.prepare_target_structure(corrupted, output_dir)
        self.assertFalse(output_dir.exists())

    def test_missing_backbone_atom_fails_closed_after_hash_validation(self) -> None:
        temporary, output_dir = self._temporary_output()
        incomplete = Path(temporary.name) / "incomplete_source.pdb"
        lines = self.source_pdb.read_text(encoding="ascii").splitlines()
        filtered = [
            line
            for line in lines
            if not (
                line.startswith("ATOM")
                and line[21:22] == "A"
                and int(line[22:26]) == 322
                and line[12:16].strip() == "CA"
            )
        ]
        incomplete.write_text("\n".join(filtered) + "\n", encoding="ascii")
        digest = hashlib.sha256(incomplete.read_bytes()).hexdigest()
        with self.assertRaisesRegex(
            target_prep.TargetStructurePreparationError,
            "Incomplete backbone at A322.*CA",
        ):
            target_prep.prepare_target_structure(
                incomplete,
                output_dir,
                expected_source_sha256=digest,
            )
        self.assertFalse(output_dir.exists())

    def test_atom_seqres_identity_mismatch_fails_closed(self) -> None:
        temporary, output_dir = self._temporary_output()
        mismatch = Path(temporary.name) / "sequence_mismatch_source.pdb"
        changed: list[str] = []
        for line in self.source_pdb.read_text(encoding="ascii").splitlines():
            if (
                line.startswith("ATOM")
                and line[21:22] == "A"
                and int(line[22:26]) == 322
            ):
                # ALA's required heavy-atom set is a subset of the original CYS
                # atoms, allowing this fixture to reach the independent
                # SEQRES-versus-ATOM identity check.
                line = line[:17] + "ALA" + line[20:]
            changed.append(line)
        mismatch.write_text("\n".join(changed) + "\n", encoding="ascii")
        digest = hashlib.sha256(mismatch.read_bytes()).hexdigest()
        with self.assertRaisesRegex(
            target_prep.TargetStructurePreparationError,
            "ATOM sequence does not equal SEQRES sequence.*322",
        ):
            target_prep.prepare_target_structure(
                mismatch,
                output_dir,
                expected_source_sha256=digest,
            )
        self.assertFalse(output_dir.exists())


if __name__ == "__main__":
    unittest.main()
