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

    def test_bootstrap_revisions_match_backend_manifest(self) -> None:
        config = json.loads(
            (PROJECT_ROOT / "config/real_model_backends.json").read_text(encoding="utf-8")
        )
        script = (PROJECT_ROOT / "deploy/autodl/bootstrap_models.sh").read_text(
            encoding="utf-8"
        )
        for backend in config["backends"]:
            with self.subTest(backend=backend["name"]):
                self.assertIn(backend["repository"], script)
                self.assertIn(backend["revision"], script)

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
