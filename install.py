#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent
REPOSITORY = "https://github.com/daisyluvr42/monad_omicskit"
SKILL_NAME = "monadomics-analysis"
MARKER = "_monadomics-install.json"
LAUNCHER = "import sys\nfrom pathlib import Path\nsys.path.insert(0, str(Path(__file__).resolve().parent))\nfrom monadomics.cli import main\n\nif __name__ == '__main__':\n    raise SystemExit(main())\n"


def read_json(path):
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def write_json(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def install(workbuddy):
    source = ROOT / "workbuddy-connector/skills" / SKILL_NAME
    package = ROOT / "src/monadomics"
    if not (source / "SKILL.md").is_file() or not (package / "cli.py").is_file():
        raise ValueError("Run install.py from a complete MonadOmics checkout or source archive.")
    config_path = workbuddy / "mcp.json"
    config = read_json(config_path) if config_path.exists() else {}
    servers = config.get("mcpServers", {})
    if not isinstance(servers, dict):
        raise ValueError(f"mcpServers must be an object: {config_path}")
    legacy = servers.get("omics", {})
    legacy_args = legacy.get("args", []) if isinstance(legacy, dict) else []
    migrate_mcp = isinstance(legacy_args, list) and any(
        isinstance(arg, str) and arg.replace("\\", "/").endswith("/mcp/omics_mcp.py")
        for arg in legacy_args
    ) and not legacy.get("disabled", False)
    destination = workbuddy / "skills" / SKILL_NAME
    workbuddy.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".monadomics-stage-", dir=workbuddy) as temporary:
        stage = Path(temporary) / SKILL_NAME
        shutil.copytree(source, stage, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", ".DS_Store"))
        shutil.copy2(ROOT / "LICENSE", stage / "LICENSE")
        scripts = stage / "scripts"
        scripts.mkdir(exist_ok=True)
        shutil.copytree(package, scripts / "monadomics", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        (scripts / "run.py").write_text(LAUNCHER, encoding="utf-8")
        command = [sys.executable, "-I", str(destination / "scripts/run.py")]
        entrypoint = (
            "\n## 本机命令入口（GitHub 安装）\n\n"
            "本 Skill 及 references 中的 `monadomics` 均使用以下命令参数数组作为入口，随后追加子命令和参数。"
            "使用实际路径并按宿主 shell 正确引用，不依赖 PATH，也不切换分析任务的工作目录。\n\n"
            "```json\n" + json.dumps(command, ensure_ascii=False) + "\n```\n"
        )
        skill_path = stage / "SKILL.md"
        content = skill_path.read_text(encoding="utf-8")
        frontmatter_end = content.index("\n---\n", 4) + len("\n---\n")
        skill_path.write_text(content[:frontmatter_end] + entrypoint + content[frontmatter_end:], encoding="utf-8")
        write_json(stage / MARKER, {"repository": REPOSITORY, "command": command})
        write_json(stage / "_user_meta.json", {"name": "MonadOmics 生信分析", "source": "userImport", "installedAt": int(time.time() * 1000)})
        # Verify the staged runtime without importing a different installed copy.
        probe = subprocess.run([sys.executable, "-I", str(scripts / "run.py"), "capabilities"], capture_output=True, text=True)
        if probe.returncode:
            raise RuntimeError(f"Staged CLI check failed: {probe.stderr or probe.stdout}")
        if not json.loads(probe.stdout).get("ok"):
            raise RuntimeError("Staged CLI did not report success.")
        backups = workbuddy / "monadomics-backups"
        backups.mkdir(exist_ok=True)
        backup = Path(tempfile.mkdtemp(prefix="install-", dir=backups))
        if destination.exists() or destination.is_symlink():
            shutil.move(str(destination), backup / SKILL_NAME)
        destination.parent.mkdir(exist_ok=True)
        try:
            shutil.move(str(stage), destination)
            if migrate_mcp:
                shutil.copy2(config_path, backup / "mcp.json")
                legacy["disabled"] = True
                write_json(config_path, config)
            old_skill = workbuddy / "skills/omics-analysis"
            if (old_skill / "SKILL.md").is_file():
                shutil.move(str(old_skill), backup / "omics-analysis")
        except OSError:
            if destination.exists():
                shutil.rmtree(destination)
            if (backup / SKILL_NAME).exists():
                shutil.move(str(backup / SKILL_NAME), destination)
            if (backup / "mcp.json").exists():
                shutil.copy2(backup / "mcp.json", config_path)
            raise
        if any(backup.iterdir()):
            print(f"Previous Skill/config backup: {backup}")
        else:
            backup.rmdir()
    print(f"Installed: {destination}")
    print("CLI check passed. R and R packages were not installed or checked.")
    print("Next: python3 install.py run doctor --group deg")
    print("If needed: python3 install.py run setup-r deg")
    print("Refresh WorkBuddy Skills/restart WorkBuddy to load the new Skill and stop the old MCP.")


def installed_command(workbuddy):
    path = workbuddy / "skills" / SKILL_NAME / MARKER
    if not path.is_file():
        raise ValueError("No GitHub installation found. Run install.py install first.")
    record = read_json(path)
    if record.get("repository") != REPOSITORY:
        raise ValueError(f"Installation belongs to a different repository: {path}")
    return record["command"]


def main(argv=None):
    parser = argparse.ArgumentParser(description="Install the MonadOmics CLI and Skill from GitHub without pip or marketplace approval.")
    parser.add_argument("--workbuddy-dir", type=Path, default=Path.home() / ".workbuddy")
    parser.add_argument("action", choices=("install", "update", "uninstall", "run"))
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="Arguments passed to the installed CLI by run.")
    args = parser.parse_args(argv)
    if args.action != "run" and args.arguments:
        parser.error(f"Unexpected arguments for {args.action}: {' '.join(args.arguments)}")
    workbuddy = args.workbuddy_dir.expanduser().resolve()
    try:
        if sys.version_info < (3, 11):
            raise ValueError("Python 3.11 or newer is required.")
        if args.action == "install":
            install(workbuddy)
        elif args.action == "update":
            if not (ROOT / ".git").exists():
                raise ValueError("Update needs a git clone. For a downloaded archive, download the new version and run install.")
            if subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip():
                raise ValueError("Checkout has local changes. Commit or preserve them before updating.")
            subprocess.run(["git", "pull", "--ff-only"], cwd=ROOT, check=True)
            os.execv(sys.executable, [sys.executable, str(ROOT / "install.py"), "--workbuddy-dir", str(workbuddy), "install"])
        elif args.action == "uninstall":
            installed_command(workbuddy)
            backups = workbuddy / "monadomics-backups"
            backups.mkdir(exist_ok=True)
            backup = Path(tempfile.mkdtemp(prefix="uninstall-", dir=backups))
            shutil.move(str(workbuddy / "skills" / SKILL_NAME), backup / SKILL_NAME)
            print(f"Uninstalled Skill and bundled CLI; backup: {backup}")
            print("R packages, analysis outputs and other connectors were kept. Old MCP remains disabled.")
        else:
            return subprocess.call([*installed_command(workbuddy), *args.arguments])
        return 0
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
