import json
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from nfl_ab_design import workflow
from nfl_ab_design.adapters.germinal import (
    PROFILE_SETTINGS as GERMINAL_PROFILES,
    UPSTREAM_PROVENANCE as GERMINAL_UPSTREAM,
    GerminalAdapterError,
    build_germinal_jobs,
    write_germinal_handoff,
)
from nfl_ab_design.adapters.iggm import (
    UPSTREAM as IGGM_UPSTREAM,
    IgGMAdapterError,
    build_iggm_jobs,
)
from nfl_ab_design.adapters.rfantibody import (
    RFANTIBODY_OFFICIAL_MAIN_SHA,
    RFantibodyAdapterError,
    build_rfantibody_plan,
)


AA1_TO_AA3 = {
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
CDR_NAMES = ("H1", "H2", "H3", "L1", "L2", "L3")
EXPECTED_RFANTIBODY_SHA = "8fe311415754e0276d1a39c87c57e69c88927a2d"
EXPECTED_IGGM_SHA = "06abc563b3fc8c7ea020543add16b69b6f8a1c8d"
EXPECTED_GERMINAL_SHA = "1e1c1a5b79884ae45abae030c9df90d9423a990a"


def _write_pdb(path: Path, sequence: str, *, chain: str = "A", start: int = 1) -> None:
    lines = []
    for serial, (offset, amino_acid) in enumerate(enumerate(sequence), start=1):
        residue_number = start + offset
        lines.append(
            f"ATOM  {serial:5d}  CA  {AA1_TO_AA3[amino_acid]:>3s} "
            f"{chain}{residue_number:4d}    {float(serial):8.3f}"
            f"{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        )
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def _write_hlt_pdb(path: Path) -> None:
    remarks = [
        f"REMARK PDBinfo-LABEL: {index} {name}"
        for index, name in enumerate(CDR_NAMES, start=1)
    ]
    atoms = []
    serial = 1
    for residue_number, chain in enumerate(("H", "H", "H", "L", "L", "L"), start=1):
        base = float(residue_number * 3)
        for atom, coordinates in (
            ("N", (base, 0.0, 0.0)),
            ("CA", (base + 1.0, 1.0, 0.0)),
            ("C", (base + 2.0, 0.0, 0.0)),
        ):
            atoms.append(
                f"ATOM  {serial:5d} {atom:^4s} ALA {chain}{residue_number:4d}    "
                f"{coordinates[0]:8.3f}{coordinates[1]:8.3f}{coordinates[2]:8.3f}"
                "  1.00 20.00           C"
            )
            serial += 1
    path.write_text("\n".join([*remarks, *atoms, "END", ""]), encoding="utf-8")


def _design_regions() -> list[dict[str, object]]:
    return [
        {
            "region": name,
            "chain": "VH" if name.startswith("H") else "VL",
            "start_1_based": position,
            "end_1_based_inclusive": position,
            "length_aa": 1,
        }
        for name, position in zip(CDR_NAMES, (2, 5, 8, 2, 5, 8), strict=True)
    ]


def _normalized_request(engine: str, *, curated_hotspots: bool = True) -> dict[str, object]:
    epitopes: list[dict[str, object]] = [
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
    ]
    if curated_hotspots:
        epitopes[0]["selected_hotspot_residue_indices"] = [2]
        epitopes[1]["selected_hotspot_residue_indices"] = [9]

    return {
        "schema": "nfl_ab_design.normalized_de_novo_request.v1",
        "campaign_mode": "paired_Fv_six_CDR_de_novo_design",
        "execution_state": "not_run",
        "result_provenance": "adapter_request_only",
        "engine": engine,
        "antigen": {
            "protein": "fixture_antigen",
            "full_sequence": "ACDEFGHIKLMN",
            "antigen_pdb_path": "",
        },
        "epitopes": epitopes,
        "templates": [
            {
                "template_id": "template_A",
                "framework_source_id": "source_A",
                "template_role": "framework_source_only",
                "masked_vh": "AXAAXAAXA",
                "masked_vl": "GXGGXGGXG",
                "design_regions": _design_regions(),
                "designed_regions": ";".join(CDR_NAMES),
            },
            {
                "template_id": "template_B",
                "framework_source_id": "source_B",
                "template_role": "framework_source_only",
                "masked_vh": "CXCCXCCXC",
                "masked_vl": "EXEEXEEXE",
                "design_regions": _design_regions(),
                "designed_regions": ";".join(CDR_NAMES),
            },
        ],
    }


class RealModelAdapterContractTest(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

        self.antigen_sequence = "ACDEFGHIKLMN"
        self.target_pdb = self.root / "target.pdb"
        self.offset_target_pdb = self.root / "target_offset_numbering.pdb"
        _write_pdb(self.target_pdb, self.antigen_sequence)
        _write_pdb(self.offset_target_pdb, self.antigen_sequence, start=101)

        self.hlt_pdbs = {
            "template_A": self.root / "template_A_hlt.pdb",
            "template_B": self.root / "template_B_hlt.pdb",
        }
        for path in self.hlt_pdbs.values():
            _write_hlt_pdb(path)

        scfv_sequences = {
            "template_A": "ASAASAASA" + "GGGGS" + "GSGGSGGSG",
            "template_B": "CSCCSCCSC" + "GGGGS" + "ESEESEESE",
        }
        self.scfv_pdbs = {
            template_id: self.root / f"{template_id}_scfv.pdb"
            for template_id in scfv_sequences
        }
        for template_id, sequence in scfv_sequences.items():
            _write_pdb(self.scfv_pdbs[template_id], sequence)

        self.rf_mapping = {
            position: {"chain_id": "A", "residue_number": position}
            for position in range(1, len(self.antigen_sequence) + 1)
        }
        self.iggm_mapping = {
            position: position for position in range(1, len(self.antigen_sequence) + 1)
        }

    def _rfantibody_plan(self):
        return build_rfantibody_plan(
            _normalized_request("RFantibody"),
            target_pdb=self.target_pdb,
            framework_hlt_pdbs=self.hlt_pdbs,
            full_coordinate_to_pdb=self.rf_mapping,
            output_root=self.root / "rfantibody_results",
            mode="smoke",
        )

    def _iggm_plan(self):
        return build_iggm_jobs(
            _normalized_request("IgGM"),
            target_pdb_path=self.target_pdb,
            pdb_antigen_sequence=self.antigen_sequence,
            pdb_antigen_chain="A",
            full_to_local_residue_map=self.iggm_mapping,
            input_dir=self.root / "iggm_inputs",
            output_dir=self.root / "iggm_results",
            profile="smoke",
        )

    def _germinal_plan(self):
        return build_germinal_jobs(
            _normalized_request("Germinal"),
            target_pdb_path=self.target_pdb,
            template_scfv_pdbs=self.scfv_pdbs,
            handoff_root=self.root / "germinal_handoff",
            profile="smoke",
        )

    def test_upstream_pins_and_four_job_campaign_shape_are_stable(self) -> None:
        rf_plan = self._rfantibody_plan()
        iggm_plan = self._iggm_plan()
        germinal_plan = self._germinal_plan()

        self.assertEqual(RFANTIBODY_OFFICIAL_MAIN_SHA, EXPECTED_RFANTIBODY_SHA)
        self.assertEqual(rf_plan.runtime_ref, EXPECTED_RFANTIBODY_SHA)
        self.assertEqual(IGGM_UPSTREAM["commit"], EXPECTED_IGGM_SHA)
        self.assertEqual(iggm_plan["upstream"]["commit"], EXPECTED_IGGM_SHA)
        self.assertEqual(GERMINAL_UPSTREAM["pinned_commit"], EXPECTED_GERMINAL_SHA)
        self.assertEqual(germinal_plan["upstream"]["pinned_commit"], EXPECTED_GERMINAL_SHA)

        expected_pairs = {
            (template_id, epitope_id)
            for template_id in ("template_A", "template_B")
            for epitope_id in ("epitope_N", "epitope_C")
        }
        self.assertEqual(len(rf_plan.jobs), 4)
        self.assertEqual(iggm_plan["job_count"], 4)
        self.assertEqual(germinal_plan["job_count"], 4)
        for jobs in (rf_plan.jobs, iggm_plan["jobs"], germinal_plan["jobs"]):
            observed = {
                (job.template_id, job.epitope_id)
                if hasattr(job, "template_id")
                else (job["template_id"], job["epitope_id"])
                for job in jobs
            }
            self.assertEqual(observed, expected_pairs)

        for job in rf_plan.jobs:
            stages = [command.stage for command in job.commands]
            self.assertNotIn("proteinmpnn_score_export", stages)
            self.assertEqual(
                stages[stages.index("proteinmpnn_six_cdr_sequence_design") + 1],
                "rf2_structure_prediction",
            )
            score_exports = [
                command.argv[-1]
                for command in job.commands
                if command.argv and "qvscorefile" in command.argv
            ]
            self.assertEqual(len(score_exports), 2)
            self.assertTrue(score_exports[0].endswith("1_rfdiffusion.qv"))
            self.assertTrue(score_exports[1].endswith("3_rf2.qv"))

    def test_germinal_pilot_uses_full_optimization_and_materializes_hashed_seed(self) -> None:
        cdr_positions = {2, 5, 8, 16, 19, 22}
        for path in self.scfv_pdbs.values():
            rewritten = []
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.startswith("ATOM"):
                    residue = int(line[22:26])
                    y = 10.0 if residue in cdr_positions else 0.0
                    line = f"{line[:38]}{y:8.3f}{line[46:]}"
                rewritten.append(line)
            path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
        handoff = build_germinal_jobs(
            _normalized_request("Germinal"),
            target_pdb_path=self.target_pdb,
            template_scfv_pdbs=self.scfv_pdbs,
            handoff_root=self.root / "pilot_handoff",
            profile="pilot",
        )
        self.assertEqual(GERMINAL_PROFILES["pilot"]["hydra_overrides"]["logits_steps"], 60)
        self.assertEqual(GERMINAL_PROFILES["pilot"]["hydra_overrides"]["softmax_steps"], 35)
        self.assertEqual(GERMINAL_PROFILES["pilot"]["hydra_overrides"]["search_steps"], 10)
        for job in handoff["jobs"]:
            seed = job["inputs"]["prepositioned_starting_complex"]
            self.assertGreaterEqual(seed["minimum_heavy_atom_distance_angstrom"], 2.5)
            self.assertLessEqual(seed["minimum_hotspot_to_cdr_distance_angstrom"], 12.0)
            expected_external = "A2" if job["epitope_id"] == "epitope_N" else "A9"
            expected_chai_hotspot = (
                "C2" if job["epitope_id"] == "epitope_N" else "K9"
            )
            self.assertEqual(
                job["target_yaml"]["values"]["external_target_hotspots"],
                expected_external,
            )
            self.assertIn(
                f'external_target_hotspots: "{expected_external}"',
                job["target_yaml"]["content"],
            )
            self.assertEqual(
                job["target_yaml"]["values"]["hotspot_residue"],
                expected_chai_hotspot,
            )
        written = write_germinal_handoff(
            handoff,
            self.root / "pilot_handoff",
            stage_template_pdbs=True,
        )
        self.assertEqual(len(written["generated_starting_complexes"]), 4)
        for path in written["generated_starting_complexes"]:
            self.assertTrue(Path(path).is_file())

    def test_planners_only_return_commands_and_never_execute_or_stage_them(self) -> None:
        forbidden_calls = (
            mock.patch.object(
                subprocess,
                "Popen",
                side_effect=AssertionError("planner attempted to launch a process"),
            ),
            mock.patch.object(
                subprocess,
                "run",
                side_effect=AssertionError("planner attempted to run a process"),
            ),
            mock.patch.object(
                os,
                "system",
                side_effect=AssertionError("planner attempted to invoke a shell"),
            ),
            mock.patch(
                "nfl_ab_design.adapters.germinal.shutil.copy2",
                side_effect=AssertionError("planner attempted to stage a file"),
            ),
        )
        with ExitStack() as stack:
            for patcher in forbidden_calls:
                stack.enter_context(patcher)
            rf_plan = self._rfantibody_plan()
            iggm_plan = self._iggm_plan()
            germinal_plan = self._germinal_plan()

        rf_manifest = rf_plan.as_dict()
        self.assertEqual(rf_manifest["execution_state"], "planned_not_executed")
        self.assertTrue(rf_plan.commands)
        self.assertTrue(
            all(
                command["execution_state"] == "planned_not_executed"
                for job in rf_manifest["jobs"]
                for command in job["commands"]
            )
        )
        self.assertEqual(iggm_plan["execution_state"], "planned_not_run")
        self.assertTrue(
            all(
                job["execution_state"] == "planned_not_run"
                for job in iggm_plan["jobs"]
            )
        )
        self.assertEqual(germinal_plan["execution_state"], "not_run")
        self.assertIs(germinal_plan["does_not_execute_external_code"], True)
        self.assertTrue(all(job["execution_state"] == "not_run" for job in germinal_plan["jobs"]))

        for path in (
            self.root / "rfantibody_results",
            self.root / "iggm_inputs",
            self.root / "iggm_results",
            self.root / "germinal_handoff",
        ):
            self.assertFalse(path.exists(), f"planner unexpectedly materialized {path}")

    def test_rfantibody_fails_loudly_without_target_mapping_or_curated_hotspots(self) -> None:
        valid_request = _normalized_request("RFantibody")
        valid_arguments = {
            "framework_hlt_pdbs": self.hlt_pdbs,
            "full_coordinate_to_pdb": self.rf_mapping,
        }
        with self.subTest("actual target PDB is required"):
            with self.assertRaisesRegex(RFantibodyAdapterError, "target antigen PDB.*not a file"):
                build_rfantibody_plan(
                    valid_request,
                    target_pdb=self.root / "missing_target.pdb",
                    **valid_arguments,
                )

        with self.subTest("explicit full-coordinate mapping is required"):
            with self.assertRaisesRegex(
                RFantibodyAdapterError,
                "Missing explicit full-coordinate-to-PDB residue mapping",
            ):
                build_rfantibody_plan(
                    valid_request,
                    target_pdb=self.target_pdb,
                    framework_hlt_pdbs=self.hlt_pdbs,
                    full_coordinate_to_pdb=None,
                )

        with self.subTest("candidate window cannot silently become selected hotspots"):
            candidate_only = _normalized_request("RFantibody", curated_hotspots=False)
            self.assertTrue(candidate_only["epitopes"][0]["candidate_hotspot_residue_indices"])
            with self.assertRaisesRegex(
                RFantibodyAdapterError,
                "candidate_hotspot_residue_indices is only an epitope window",
            ):
                build_rfantibody_plan(
                    candidate_only,
                    target_pdb=self.target_pdb,
                    **valid_arguments,
                )

    def test_rfantibody_rejects_parser_incompatible_hlt_frameworks(self) -> None:
        request = _normalized_request("RFantibody")

        def assert_rejected(path: Path, pattern: str) -> None:
            frameworks = dict(self.hlt_pdbs)
            frameworks["template_A"] = path
            with self.assertRaisesRegex(RFantibodyAdapterError, pattern):
                build_rfantibody_plan(
                    request,
                    target_pdb=self.target_pdb,
                    framework_hlt_pdbs=frameworks,
                    full_coordinate_to_pdb=self.rf_mapping,
                    output_root=self.root / "rejected_hlt_results",
                    mode="smoke",
                )

        source_lines = self.hlt_pdbs["template_A"].read_text(
            encoding="utf-8"
        ).splitlines()

        with self.subTest("Chothia insertion codes would collapse parser keys"):
            insertion_path = self.root / "framework_with_H1A.pdb"
            insertion_lines = []
            for line in source_lines:
                if (
                    line.startswith("ATOM  ")
                    and line[21] == "H"
                    and line[22:26].strip() == "2"
                ):
                    line = f"{line[:22]}{1:4d}A{line[27:]}"
                insertion_lines.append(line)
            insertion_path.write_text(
                "\n".join(insertion_lines) + "\n", encoding="utf-8"
            )
            assert_rejected(insertion_path, "retains insertion code H1A")

        with self.subTest("every residue needs N CA C"):
            missing_path = self.root / "framework_missing_N.pdb"
            missing_path.write_text(
                "\n".join(
                    line
                    for line in source_lines
                    if not (
                        line.startswith("ATOM  ")
                        and line[21] == "H"
                        and line[22:26].strip() == "1"
                        and line[12:16].strip() == "N"
                    )
                )
                + "\n",
                encoding="utf-8",
            )
            assert_rejected(missing_path, "H1 is missing backbone atom")

        with self.subTest("null backbone frames are rejected"):
            degenerate_path = self.root / "framework_degenerate_frame.pdb"
            degenerate_lines = []
            for line in source_lines:
                if (
                    line.startswith("ATOM  ")
                    and line[21] == "H"
                    and line[22:26].strip() == "1"
                    and line[12:16].strip() in {"N", "CA", "C"}
                ):
                    line = (
                        f"{line[:30]}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{line[54:]}"
                    )
                degenerate_lines.append(line)
            degenerate_path.write_text(
                "\n".join(degenerate_lines) + "\n", encoding="utf-8"
            )
            assert_rejected(degenerate_path, "null or degenerate N-CA-C")

        with self.subTest("CDR labels must address real absolute residues"):
            out_of_range_path = self.root / "framework_bad_remark.pdb"
            out_of_range_lines = list(source_lines)
            out_of_range_lines[0] = "REMARK PDBinfo-LABEL: 99 H1"
            out_of_range_path.write_text(
                "\n".join(out_of_range_lines) + "\n", encoding="utf-8"
            )
            assert_rejected(out_of_range_path, "indices outside 1..6")

    def test_iggm_fails_loudly_without_target_or_complete_coordinate_mapping(self) -> None:
        request = _normalized_request("IgGM")
        with self.subTest("actual target PDB is required"):
            with self.assertRaisesRegex(IgGMAdapterError, "Target antigen PDB does not exist"):
                build_iggm_jobs(
                    request,
                    target_pdb_path=self.root / "missing_target.pdb",
                    pdb_antigen_sequence=self.antigen_sequence,
                    pdb_antigen_chain="A",
                    full_to_local_residue_map=self.iggm_mapping,
                )

        with self.subTest("mapping must be explicit and non-empty"):
            with self.assertRaisesRegex(IgGMAdapterError, "must be a non-empty mapping"):
                build_iggm_jobs(
                    request,
                    target_pdb_path=self.target_pdb,
                    pdb_antigen_sequence=self.antigen_sequence,
                    pdb_antigen_chain="A",
                    full_to_local_residue_map={},
                )

        with self.subTest("mapping must cover the complete PDB chain"):
            incomplete_mapping = dict(self.iggm_mapping)
            incomplete_mapping.pop(len(self.antigen_sequence))
            with self.assertRaisesRegex(
                IgGMAdapterError,
                "cover every local PDB sequence position",
            ):
                build_iggm_jobs(
                    request,
                    target_pdb_path=self.target_pdb,
                    pdb_antigen_sequence=self.antigen_sequence,
                    pdb_antigen_chain="A",
                    full_to_local_residue_map=incomplete_mapping,
                )

    def test_germinal_fails_loudly_without_required_structure_inputs(self) -> None:
        valid_request = _normalized_request("Germinal")
        with self.subTest("actual target PDB is required"):
            with self.assertRaisesRegex(GerminalAdapterError, "target PDB does not exist"):
                build_germinal_jobs(
                    valid_request,
                    target_pdb_path=self.root / "missing_target.pdb",
                    template_scfv_pdbs=self.scfv_pdbs,
                )

        with self.subTest("each template needs a real coordinate PDB"):
            missing_template = dict(self.scfv_pdbs)
            missing_template["template_A"] = self.root / "missing_scfv.pdb"
            with self.assertRaisesRegex(GerminalAdapterError, "scFv PDB.*does not exist"):
                build_germinal_jobs(
                    valid_request,
                    target_pdb_path=self.target_pdb,
                    template_scfv_pdbs=missing_template,
                )

        with self.subTest("non-identity target numbering requires a mapping"):
            with self.assertRaisesRegex(
                GerminalAdapterError,
                "Provide a correct target_residue_map",
            ):
                build_germinal_jobs(
                    valid_request,
                    target_pdb_path=self.offset_target_pdb,
                    template_scfv_pdbs=self.scfv_pdbs,
                    target_residue_map=None,
                )

        with self.subTest("candidate window cannot silently become selected hotspots"):
            candidate_only = _normalized_request("Germinal", curated_hotspots=False)
            self.assertTrue(candidate_only["epitopes"][0]["candidate_hotspot_residue_indices"])
            with self.assertRaisesRegex(
                GerminalAdapterError,
                "refusing to treat every window residue",
            ):
                build_germinal_jobs(
                    candidate_only,
                    target_pdb_path=self.target_pdb,
                    template_scfv_pdbs=self.scfv_pdbs,
                )

    def test_request_export_emits_three_blocked_normalized_engine_requests(self) -> None:
        full_sequence = workflow.parse_genpept_sequence(workflow.GENPEPT_PATH)
        templates = workflow.load_antibodies(workflow.ANTIBODY_FASTA_PATH)
        epitope_rows = workflow.build_epitope_windows(full_sequence, "280-375")
        export_dir = self.root / "exports"
        design_config = {
            "epitopes": ["helix_surface_323_331", "C_boundary_368_377"]
        }

        with mock.patch.object(workflow, "EXPORT_DIR", export_dir):
            exported = workflow.export_de_novo_model_requests(
                templates,
                epitope_rows,
                full_sequence,
                design_config,
            )

        request_dir = export_dir / "design_requests"
        index = json.loads(
            (request_dir / "design_request_index.json").read_text(encoding="utf-8")
        )
        self.assertEqual(index["execution_state"], "not_run")
        self.assertEqual(index["engines"], ["RFantibody", "IgGM", "Germinal"])
        self.assertEqual(index["template_count"], 2)
        self.assertEqual(index["epitope_count"], 2)
        self.assertEqual(len(index["request_files"]), 3)
        self.assertEqual(len(exported["design_request_files"]), 3)

        for engine in index["engines"]:
            request = json.loads(
                (request_dir / f"{engine.lower()}_design_request.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(request["schema"], "nfl_ab_design.normalized_de_novo_request.v1")
            self.assertEqual(request["engine"], engine)
            self.assertEqual(request["execution_state"], "not_run")
            self.assertTrue(request["adapter_state"].startswith("blocked_missing_"))
            self.assertEqual(len(request["templates"]), 2)
            self.assertEqual(len(request["epitopes"]), 2)
            self.assertEqual(
                {region["region"] for region in request["templates"][0]["design_regions"]},
                set(CDR_NAMES),
            )
            self.assertNotIn("selected_hotspot_residue_indices", request["epitopes"][0])
            self.assertTrue(request["epitopes"][0]["candidate_hotspot_residue_indices"])


if __name__ == "__main__":
    unittest.main()
