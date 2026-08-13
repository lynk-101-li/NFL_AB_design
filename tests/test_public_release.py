import json
from pathlib import Path
import re
import subprocess
import unittest


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class PublicReleaseContractTest(unittest.TestCase):
    def test_deployment_scripts_are_valid_shell(self) -> None:
        for path in sorted((PROJECT_ROOT / "deploy/autodl").glob("*.sh")):
            completed = subprocess.run(
                ["bash", "-n", str(path)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_model_registry_matches_submodule_gitlinks(self) -> None:
        config = json.loads(
            (PROJECT_ROOT / "config/model_components.json").read_text(encoding="utf-8")
        )
        script = (PROJECT_ROOT / "deploy/autodl/bootstrap_models.sh").read_text(
            encoding="utf-8"
        )
        gitmodules = (PROJECT_ROOT / ".gitmodules").read_text(encoding="utf-8")
        self.assertEqual(config["schema"], "nfl_ab_design.model_components.v1")
        self.assertEqual(len(config["components"]), 9)
        self.assertEqual(
            {item["category"] for item in config["components"]},
            {"antibody_design", "structure_prediction"},
        )
        for component in config["components"]:
            with self.subTest(component=component["id"]):
                self.assertIn(component["submodule_path"], gitmodules)
                self.assertIn(component["repository"], gitmodules)
                staged = subprocess.run(
                    ["git", "ls-files", "--stage", component["submodule_path"]],
                    cwd=PROJECT_ROOT,
                    check=True,
                    capture_output=True,
                    text=True,
                ).stdout.split()
                self.assertEqual(staged[0], "160000")
                self.assertEqual(staged[1], component["revision"])
                self.assertIn(component["id"], script)

    def test_public_text_has_no_private_machine_or_assistant_trace(self) -> None:
        forbidden = re.compile(
            r"(?i)(codex|chatgpt|co-authored-by|connect\.westd|seetacloud|"
            r"/Volumes/IXUNICS|/Users/rinck|autodl-pro-[0-9a-f]+)"
        )
        roots = (
            PROJECT_ROOT / "README.md",
            PROJECT_ROOT / "THIRD_PARTY_NOTICES.md",
            PROJECT_ROOT / "config",
            PROJECT_ROOT / "deploy",
            PROJECT_ROOT / "docs",
            PROJECT_ROOT / "input",
            PROJECT_ROOT / "scripts",
            PROJECT_ROOT / "src",
        )
        failures: list[str] = []
        for root in roots:
            paths = [root] if root.is_file() else root.rglob("*")
            for path in paths:
                if not path.is_file() or path.suffix.lower() in {".pdb", ".pyc"}:
                    continue
                try:
                    text = path.read_text(encoding="utf-8")
                except UnicodeDecodeError:
                    continue
                match = forbidden.search(text)
                if match:
                    failures.append(f"{path.relative_to(PROJECT_ROOT)}: {match.group(0)}")
        self.assertEqual(failures, [])

    def test_private_execution_files_are_not_present(self) -> None:
        self.assertFalse((PROJECT_ROOT / "config/target_structure_manifest.json").exists())
        self.assertFalse((PROJECT_ROOT / "runtime_attestation.json").exists())
        for directory in ("handoffs", "results", "final"):
            path = PROJECT_ROOT / "real_runs" / directory
            self.assertFalse(path.exists() and any(path.rglob("*")))


if __name__ == "__main__":
    unittest.main()
