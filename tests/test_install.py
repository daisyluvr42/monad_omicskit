import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class InstallTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="omics 安装 ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.repo = self.root / "repo"
        self.repo.mkdir()
        shutil.copy2(ROOT / "install.py", self.repo / "install.py")
        for name in ("src", "workbuddy-connector"):
            shutil.copytree(ROOT / name, self.repo / name, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info"))
        self.workbuddy = self.root / "workbuddy"
        self.skill = self.workbuddy / "skills/monadomics-analysis"

    def invoke(self, action, *arguments):
        return subprocess.run([sys.executable, str(self.repo / "install.py"), "--workbuddy-dir", str(self.workbuddy), action, *arguments], cwd=self.root, capture_output=True, text=True)

    def git(self, *arguments, cwd=None):
        return subprocess.run(["git", *arguments], cwd=cwd or self.repo, check=True, capture_output=True, text=True).stdout

    def test_install_bundles_exact_runtime_and_runs_outside_checkout(self):
        result = self.invoke("install")
        self.assertEqual(result.returncode, 0, result.stderr)
        source = self.repo / "src/monadomics"
        for path in source.rglob("*"):
            if path.is_file():
                self.assertEqual(path.read_bytes(), (self.skill / "scripts/monadomics" / path.relative_to(source)).read_bytes())
        for arguments in [("--version",), ("capabilities",), ("schema", "deg")]:
            result = self.invoke("run", *arguments)
            self.assertEqual(result.returncode, 0, result.stderr)
        command = json.loads((self.skill / "_monadomics-install.json").read_text())["command"]
        shutil.rmtree(self.repo)
        result = subprocess.run([*command, "capabilities"], cwd=self.root, env={**os.environ, "PATH": "", "PYTHONPATH": "/nonexistent"}, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["count"], 21)

    def test_migration_preserves_other_connectors_outputs_and_backups(self):
        self.workbuddy.mkdir()
        config = {"setting": "keep", "mcpServers": {"scholar": {"command": "keep"}, "omics": {"command": "python3", "args": ["/old repo/mcp/omics_mcp.py"], "env": {"X": "keep"}}}}
        (self.workbuddy / "mcp.json").write_text(json.dumps(config))
        old = self.workbuddy / "skills/omics-analysis"
        old.mkdir(parents=True)
        (old / "SKILL.md").write_text("old skill")
        output = self.workbuddy / "workspace/omics/result.csv"
        output.parent.mkdir(parents=True)
        output.write_text("keep data")
        result = self.invoke("install")
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = json.loads(json.dumps(config))
        expected["mcpServers"]["omics"]["disabled"] = True
        self.assertEqual(json.loads((self.workbuddy / "mcp.json").read_text()), expected)
        backup = next((self.workbuddy / "monadomics-backups").glob("install-*"))
        self.assertEqual(json.loads((backup / "mcp.json").read_text()), config)
        self.assertEqual((backup / "omics-analysis/SKILL.md").read_text(), "old skill")
        self.assertFalse(old.exists())
        self.assertEqual(self.invoke("uninstall").returncode, 0)
        self.assertFalse(self.skill.exists())
        self.assertEqual(output.read_text(), "keep data")
        self.assertEqual(json.loads((self.workbuddy / "mcp.json").read_text()), expected)

    def test_invalid_config_and_failed_stage_keep_existing_install(self):
        self.assertEqual(self.invoke("install").returncode, 0)
        original = (self.skill / "SKILL.md").read_bytes()
        (self.workbuddy / "mcp.json").write_text("invalid json")
        self.assertNotEqual(self.invoke("install").returncode, 0)
        self.assertEqual((self.skill / "SKILL.md").read_bytes(), original)
        (self.workbuddy / "mcp.json").write_text("{}")
        (self.repo / "src/monadomics/cli.py").write_text("raise RuntimeError('failed stage')")
        self.assertNotEqual(self.invoke("install").returncode, 0)
        self.assertEqual((self.skill / "SKILL.md").read_bytes(), original)
        self.assertEqual(self.invoke("run", "--version").returncode, 0)

    def test_cli_failure_propagates(self):
        self.assertEqual(self.invoke("install").returncode, 0)
        result = self.invoke("run", "deg", "--params", str(self.root / "missing.json"))
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(json.loads(result.stdout)["ok"])

    def test_uninstall_refuses_unmanaged_skill(self):
        self.skill.mkdir(parents=True)
        (self.skill / "SKILL.md").write_text("personal skill")
        self.assertNotEqual(self.invoke("uninstall").returncode, 0)
        self.assertEqual((self.skill / "SKILL.md").read_text(), "personal skill")

    def test_update_fast_forwards_and_uses_new_installer(self):
        self.git("init", "-b", "main")
        self.git("add", ".")
        self.git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-m", "initial")
        remote = self.root / "remote.git"
        self.git("clone", "--bare", str(self.repo), str(remote))
        self.git("remote", "add", "origin", str(remote))
        self.git("fetch", "origin")
        self.git("branch", "--set-upstream-to=origin/main", "main")
        author = self.root / "author"
        self.git("clone", str(remote), str(author))
        installer = author / "install.py"
        installer.write_text(installer.read_text().replace('Installed: {destination}', 'Fresh installer: {destination}'))
        reference = author / "workbuddy-connector/skills/monadomics-analysis/references/data-preparation.md"
        reference.write_text(reference.read_text() + "\nUpdated reference\n")
        self.git("add", ".", cwd=author)
        self.git("-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-m", "update", cwd=author)
        self.git("push", "origin", "main", cwd=author)
        result = self.invoke("update")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Fresh installer:", result.stdout)
        self.assertTrue((self.skill / "references/data-preparation.md").read_text().endswith("Updated reference\n"))
        (self.repo / "install.py").write_text("# local change\n" + (self.repo / "install.py").read_text())
        result = self.invoke("update")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("local changes", result.stderr)


if __name__ == "__main__":
    unittest.main()
