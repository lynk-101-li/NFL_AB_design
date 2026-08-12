import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import execute_real_model_jobs as executor


class ExecutorFixture:
    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self.manifest_path = self.root / "unified_handoff_manifest.json"
        self.attestation_path = self.root / "runtime_attestation.json"
        self.runtime_root = self.root / "runtime"
        self.runtime_root.mkdir(parents=True)
        self.iggm_repo = self.runtime_root / "iggm"
        self.germinal_repo = self.runtime_root / "germinal"
        self.iggm_repo.mkdir()
        self.germinal_repo.mkdir()
        (self.iggm_repo / "design.py").write_text("", encoding="utf-8")
        (self.germinal_repo / "run_germinal.py").write_text("", encoding="utf-8")
        (self.germinal_repo / "configs" / "target").mkdir(parents=True)

        inputs = self.root / "inputs"
        inputs.mkdir()
        self.target = inputs / "target.pdb"
        self.rf_framework = inputs / "framework_hlt.pdb"
        self.iggm_fasta = inputs / "masked.fasta"
        self.scfv = inputs / "template_scfv.pdb"
        self.target_yaml = self.root / "germinal" / "target_configs" / "nfl_canary.yaml"
        self.target.write_text("TARGET\n", encoding="utf-8")
        self.rf_framework.write_text("FRAMEWORK\n", encoding="utf-8")
        self.iggm_fasta.write_text(">H\nAXA\n", encoding="utf-8")
        self.scfv.write_text("SCFV\n", encoding="utf-8")
        self.target_yaml.parent.mkdir(parents=True)
        self.target_yaml.write_text("target_name: nfl_canary\n", encoding="utf-8")

        germinal_pdb_dir = (
            self.root / "germinal" / "jobs" / "germinal_job" / "pdbs"
        )
        self.jobs = {
            "RFantibody": [
                self._job(
                    "RFantibody",
                    "rf_job",
                    [
                        self._command(
                            "rf_backbone",
                            ["rfantibody-backbone", "--out", str(self.root / "results" / "rf")],
                        ),
                        self._command(
                            "rf_sequence",
                            ["rfantibody-sequence", "--in", str(self.root / "results" / "rf")],
                        ),
                    ],
                    [
                        self._artifact("target_pdb", self.target),
                        self._artifact("framework_hlt_pdb", self.rf_framework),
                    ],
                )
            ],
            "IgGM": [
                self._job(
                    "IgGM",
                    "iggm_job",
                    [
                        self._command(
                            "iggm_design",
                            ["python", "design.py", "--output", str(self.root / "results" / "iggm")],
                            cwd=self.iggm_repo,
                        )
                    ],
                    [
                        self._artifact("target_pdb", self.target),
                        self._artifact("masked_hla_fasta", self.iggm_fasta),
                    ],
                )
            ],
            "Germinal": [
                self._job(
                    "Germinal",
                    "germinal_job",
                    [
                        self._command(
                            "germinal_scfv_design_and_filter",
                            [
                                "python",
                                "run_germinal.py",
                                "run=scfv",
                                "target=nfl_canary",
                                f"pdb_dir={germinal_pdb_dir}",
                            ],
                            cwd=self.germinal_repo,
                        )
                    ],
                    [
                        self._artifact("target_pdb", self.target),
                        self._artifact("template_scfv_pdb", self.scfv),
                        self._artifact("generated_target_yaml", self.target_yaml),
                    ],
                )
            ],
        }
        identity = {"fixture": "executor", "seed": 7}
        identity_hash = executor._json_sha256(identity)
        engines = []
        for engine in executor.EXPECTED_ENGINES:
            execution_jobs = self.jobs[engine]
            engines.append(
                {
                    "engine": engine,
                    "execution_state": "planned_not_executed",
                    "ready_for_execution": False,
                    "required_upstream_revision": f"{engine.lower()}-revision",
                    "selected_job_count": 1,
                    "selected_job_ids": [execution_jobs[0]["job_id"]],
                    "execution_jobs": execution_jobs,
                    # Poison data proves the executor neither follows nor falls
                    # back to native/full validation plans.
                    "native_plan": {
                        "path": str(self.root / f"{engine}.native.json"),
                        "submission_allowed": False,
                    },
                    "jobs": [
                        {
                            "job_id": "unauthorized_job",
                            "commands": [{"argv": ["MUST_NOT_RUN"]}],
                        }
                    ],
                }
            )
        self.manifest = {
            "schema": executor.UNIFIED_HANDOFF_SCHEMA,
            "handoff_id": f"nfl_handoff_{identity_hash[:24]}",
            "identity_sha256": identity_hash,
            "handoff_identity": identity,
            "execution_state": "planned_not_executed",
            "does_not_execute_external_models": True,
            "source_run": {
                "source_integrity": {
                    "ready_for_execution": True,
                    "state": "verified_by_design_request_index",
                }
            },
            "handoff_location": {"path": str(self.root)},
            "engines": engines,
        }
        self.attestation = {
            "schema": executor.RUNTIME_ATTESTATION_SCHEMA,
            "handoff_id": self.manifest["handoff_id"],
            "identity_sha256": identity_hash,
            "handoff_manifest_sha256": "",
            "engines": [],
        }
        self.rewrite()

    def _artifact(self, role: str, path: Path) -> dict[str, str]:
        return {
            "role": role,
            "path": str(path),
            "sha256": executor._file_sha256(path),
        }

    def _command(
        self,
        stage: str,
        argv: list[str],
        *,
        cwd: Path | None = None,
    ) -> dict[str, object]:
        return {
            "stage": stage,
            "argv": argv,
            "working_directory": str(cwd) if cwd is not None else None,
            "expected_outputs": [
                {"path": str(self.root / "results" / stage / "result.json")}
            ],
        }

    def _job(
        self,
        engine: str,
        job_id: str,
        commands: list[dict[str, object]],
        artifacts: list[dict[str, str]],
    ) -> dict[str, object]:
        return {
            "job_id": job_id,
            "engine": engine,
            "profile": "smoke",
            "geometry": f"{engine}_fixture_geometry",
            "execution_state": "planned_not_executed",
            "template_id": "template_7-H11-D3-2-C7",
            "epitope_id": "helix_surface_323_331",
            "selected_for_execution": True,
            "execution_disposition": "selected_for_execution",
            "commands": commands,
            "input_artifacts": artifacts,
            "unresolved_blockers": [],
        }

    def rewrite(self) -> None:
        self.manifest_path.write_text(
            json.dumps(self.manifest, indent=2) + "\n", encoding="utf-8"
        )
        self.attestation["handoff_id"] = self.manifest["handoff_id"]
        self.attestation["identity_sha256"] = self.manifest["identity_sha256"]
        self.attestation["handoff_manifest_sha256"] = executor._file_sha256(
            self.manifest_path
        )
        self.attestation["engines"] = [
            {
                "engine": engine["engine"],
                "revision": engine["required_upstream_revision"],
                "execution_jobs_sha256": executor._json_sha256(
                    engine["execution_jobs"]
                ),
                "checkpoint_sha256": {"fixture_checkpoint": "a" * 64},
                "ready": True,
                "overrides_manifest_ready_for_execution": True,
            }
            for engine in self.manifest["engines"]
        ]
        self.attestation_path.write_text(
            json.dumps(self.attestation, indent=2) + "\n", encoding="utf-8"
        )


class RealModelExecutorTest(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.fixture = ExecutorFixture(Path(temporary.name))

    def _run(self, **kwargs):
        return executor.execute_handoff(
            handoff_manifest_path=self.fixture.manifest_path,
            runtime_attestation_path=self.fixture.attestation_path,
            **kwargs,
        )

    def test_default_is_strict_dry_run_and_ignores_native_and_full_jobs(self) -> None:
        installed_yaml = (
            self.fixture.germinal_repo / "configs" / "target" / "nfl_canary.yaml"
        )
        staged_scfv = (
            self.fixture.root
            / "germinal"
            / "jobs"
            / "germinal_job"
            / "pdbs"
            / "scfv.pdb"
        )
        with mock.patch.object(
            executor.subprocess,
            "run",
            side_effect=AssertionError("dry-run launched a process"),
        ), mock.patch.object(
            executor.shutil,
            "copy2",
            side_effect=AssertionError("dry-run staged a file"),
        ):
            preview = self._run()
        self.assertEqual(preview["status"], "dry_run_validated_no_process_or_staging")
        self.assertEqual([job["job_id"] for job in preview["jobs"]], [
            "rf_job",
            "iggm_job",
            "germinal_job",
        ])
        self.assertNotIn("unauthorized_job", json.dumps(preview))
        self.assertNotIn("MUST_NOT_RUN", json.dumps(preview))
        self.assertFalse(installed_yaml.exists())
        self.assertFalse(staged_scfv.exists())

    def test_execute_is_serial_shell_false_stages_germinal_and_records_streams(self) -> None:
        completed = [
            subprocess.CompletedProcess([], 0, stdout=f"out-{index}", stderr=f"err-{index}")
            for index in range(4)
        ]
        with mock.patch.object(executor.subprocess, "run", side_effect=completed) as run:
            report = self._run(execute=True)
        self.assertEqual(report["status"], "succeeded")
        self.assertEqual(run.call_count, 4)
        observed_first_tokens = [call.args[0][0] for call in run.call_args_list]
        self.assertEqual(
            observed_first_tokens,
            ["rfantibody-backbone", "rfantibody-sequence", "python", "python"],
        )
        for call in run.call_args_list:
            self.assertIs(call.kwargs["shell"], False)
            self.assertIs(call.kwargs["check"], False)
            self.assertIs(call.kwargs["capture_output"], True)
            self.assertIs(call.kwargs["text"], True)
        self.assertIsNone(run.call_args_list[0].kwargs["cwd"])
        self.assertEqual(run.call_args_list[2].kwargs["cwd"], str(self.fixture.iggm_repo))
        self.assertEqual(run.call_args_list[3].kwargs["cwd"], str(self.fixture.germinal_repo))

        state_path = self.fixture.root / "execution" / "execution_report.json"
        saved = json.loads(state_path.read_text(encoding="utf-8"))
        self.assertEqual(saved["status"], "succeeded")
        first_command = saved["jobs"][0]["attempts"][0]["commands"][0]
        self.assertEqual(first_command["argv"][0], "rfantibody-backbone")
        self.assertEqual(first_command["stdout"], "out-0")
        self.assertEqual(first_command["stderr"], "err-0")
        self.assertEqual(first_command["exit_code"], 0)
        germinal_attempt = saved["jobs"][2]["attempts"][0]
        self.assertEqual(len(germinal_attempt["staging"]), 2)
        self.assertTrue(
            (self.fixture.germinal_repo / "configs" / "target" / "nfl_canary.yaml").is_file()
        )
        self.assertTrue(
            (
                self.fixture.root
                / "germinal"
                / "jobs"
                / "germinal_job"
                / "pdbs"
                / "scfv.pdb"
            ).is_file()
        )

    def test_failure_is_recorded_and_stops_before_next_job(self) -> None:
        completed = [
            subprocess.CompletedProcess([], 0, stdout="ok", stderr=""),
            subprocess.CompletedProcess([], 17, stdout="partial", stderr="boom"),
        ]
        with mock.patch.object(executor.subprocess, "run", side_effect=completed) as run:
            with self.assertRaisesRegex(executor.RealModelExecutionFailed, "code 17"):
                self._run(execute=True)
        self.assertEqual(run.call_count, 2)
        saved = json.loads(
            (self.fixture.root / "execution" / "execution_report.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(saved["status"], "failed_stopped_on_first_error")
        self.assertEqual(len(saved["jobs"]), 1)
        self.assertEqual(saved["jobs"][0]["status"], "failed")
        failed = saved["jobs"][0]["attempts"][0]["commands"][1]
        self.assertEqual(failed["stdout"], "partial")
        self.assertEqual(failed["stderr"], "boom")
        self.assertEqual(failed["exit_code"], 17)

    def test_resume_skips_hash_matching_successful_jobs(self) -> None:
        successes = [
            subprocess.CompletedProcess([], 0, stdout="ok", stderr="")
            for _ in range(4)
        ]
        with mock.patch.object(executor.subprocess, "run", side_effect=successes):
            self._run(execute=True)
        with mock.patch.object(
            executor.subprocess,
            "run",
            side_effect=AssertionError("resume reran a successful job"),
        ) as run:
            resumed = self._run(execute=True, resume=True)
        run.assert_not_called()
        self.assertEqual(resumed["status"], "succeeded")
        self.assertTrue(
            all(
                job["resume_disposition"] == "skipped_previously_succeeded"
                for job in resumed["jobs"]
            )
        )

    def test_resume_refuses_a_forged_success_record(self) -> None:
        successes = [
            subprocess.CompletedProcess([], 0, stdout="ok", stderr="")
            for _ in range(4)
        ]
        with mock.patch.object(executor.subprocess, "run", side_effect=successes):
            self._run(execute=True)
        state_path = self.fixture.root / "execution" / "execution_report.json"
        report = json.loads(state_path.read_text(encoding="utf-8"))
        report["jobs"][0]["attempts"][-1]["commands"][0]["exit_code"] = 9
        state_path.write_text(json.dumps(report), encoding="utf-8")
        with mock.patch.object(executor.subprocess, "run") as run:
            with self.assertRaisesRegex(
                executor.RealModelExecutorError, "zero-exit success record"
            ):
                self._run(execute=True, resume=True)
        run.assert_not_called()

    def test_never_falls_back_to_jobs_or_native_plan(self) -> None:
        del self.fixture.manifest["engines"][0]["execution_jobs"]
        self.fixture.manifest_path.write_text(
            json.dumps(self.fixture.manifest, indent=2) + "\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(
            executor.RealModelExecutorError,
            "native_plan/jobs are rejected",
        ):
            self._run()

        native_path = self.fixture.root / "rfantibody_plan.json"
        native_path.write_text(
            json.dumps({"schema": "rfantibody.native.v1", "jobs": []}),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(
            executor.RealModelExecutorError, "native plans/jobs are not executable",
        ):
            executor.execute_handoff(
                handoff_manifest_path=native_path,
                runtime_attestation_path=self.fixture.attestation_path,
            )

    def test_attestation_must_bind_all_hashes_revision_checkpoints_and_override(self) -> None:
        mutations = {
            "wrong handoff": lambda value: value.__setitem__("handoff_id", "wrong"),
            "wrong manifest hash": lambda value: value.__setitem__(
                "handoff_manifest_sha256", "b" * 64
            ),
            "wrong revision": lambda value: value["engines"][0].__setitem__(
                "revision", "wrong"
            ),
            "wrong jobs hash": lambda value: value["engines"][1].__setitem__(
                "execution_jobs_sha256", "c" * 64
            ),
            "empty checkpoints": lambda value: value["engines"][2].__setitem__(
                "checkpoint_sha256", {}
            ),
            "not ready": lambda value: value["engines"][0].__setitem__(
                "ready", False
            ),
            "no explicit override": lambda value: value["engines"][0].__setitem__(
                "overrides_manifest_ready_for_execution", False
            ),
        }
        original = json.loads(json.dumps(self.fixture.attestation))
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                value = json.loads(json.dumps(original))
                mutate(value)
                self.fixture.attestation_path.write_text(
                    json.dumps(value), encoding="utf-8"
                )
                with self.assertRaises(executor.RealModelExecutorError):
                    self._run()

    def test_germinal_staging_cannot_escape_isolated_manifest_workspace(self) -> None:
        command = self.fixture.manifest["engines"][2]["execution_jobs"][0]["commands"][0]
        command["argv"][-1] = f"pdb_dir={self.fixture.root.parent / 'escaped_pdb_dir'}"
        self.fixture.rewrite()
        with self.assertRaisesRegex(
            executor.RealModelExecutorError,
            "pdb_dir is not its isolated manifest workspace",
        ):
            self._run()
        self.assertFalse((self.fixture.root.parent / "escaped_pdb_dir").exists())


if __name__ == "__main__":
    unittest.main()
