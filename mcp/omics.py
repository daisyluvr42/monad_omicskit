#!/usr/bin/env python3
"""Install, update, or uninstall Omics Kit for WorkBuddy."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MCP = ROOT / "mcp" / "omics_mcp.py"
WORKBUDDY_CONFIG = Path.home() / ".workbuddy" / "mcp.json"
WORKBUDDY_SKILLS_DIR = Path.home() / ".workbuddy" / "skills"
OUTPUT_DIR = Path.home() / ".workbuddy" / "workspace" / "omics"

SKILLS = {
    "omics-analysis": "Omics 生信分析",
}


def read_config(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit(f"WorkBuddy MCP config is not valid JSON: {path}\n{exc}") from exc
    if not isinstance(value, dict):
        raise SystemExit(f"WorkBuddy MCP config must contain a JSON object: {path}")
    return value


def write_config(path: Path, config: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def find_python3() -> str:
    configured = os.getenv("OMICS_BOOTSTRAP_PYTHON", "").strip()
    candidates = [configured] if configured else []
    candidates.extend(["python3.12", "python3.13", "python3.11", "python3"])
    for name in candidates:
        path = shutil.which(name)
        if path:
            return path
    raise SystemExit("Python 3.11 or newer is required.")


def find_rscript() -> str | None:
    found = shutil.which("Rscript")
    if found:
        return found
    for candidate in [
        Path("/Library/Frameworks/R.framework/Resources/bin/Rscript"),
        Path("/opt/homebrew/bin/Rscript"),
        Path("/usr/local/bin/Rscript"),
    ]:
        if candidate.exists():
            return str(candidate)
    return None


def ensure_runtime() -> None:
    if not MCP.exists():
        raise SystemExit(f"Missing MCP server: {MCP}")


def server_config() -> dict:
    environment = {"OMICS_OUTPUT_DIR": str(OUTPUT_DIR)}
    rscript = find_rscript()
    # Pin the interpreter: GUI-launched clients often have a minimal PATH.
    if rscript:
        environment["OMICS_RSCRIPT"] = rscript
    return {
        "command": find_python3(),
        "args": [str(MCP)],
        "env": environment,
        "disabled": False,
    }


def install_mcp_config() -> Path:
    config = read_config(WORKBUDDY_CONFIG)
    servers = config.setdefault("mcpServers", {})
    if not isinstance(servers, dict):
        raise SystemExit(f"mcpServers must be an object in {WORKBUDDY_CONFIG}")
    servers["omics"] = server_config()
    write_config(WORKBUDDY_CONFIG, config)
    return WORKBUDDY_CONFIG


def uninstall_mcp_config() -> bool:
    config = read_config(WORKBUDDY_CONFIG)
    servers = config.get("mcpServers")
    if not isinstance(servers, dict) or "omics" not in servers:
        return False
    del servers["omics"]
    write_config(WORKBUDDY_CONFIG, config)
    return True


def install_skills() -> list[Path]:
    WORKBUDDY_SKILLS_DIR.mkdir(parents=True, exist_ok=True)
    installed = []
    installed_at = int(time.time() * 1000)
    for skill_name, display_name in SKILLS.items():
        source = ROOT / "skills" / skill_name
        if not (source / "SKILL.md").exists():
            raise SystemExit(f"Missing Skill: {source / 'SKILL.md'}")
        destination = WORKBUDDY_SKILLS_DIR / skill_name
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source, destination)
        metadata = {"name": display_name, "installedAt": installed_at, "source": "userImport"}
        (destination / "_user_meta.json").write_text(
            json.dumps(metadata, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        installed.append(destination)
    return installed


def uninstall_skills() -> list[Path]:
    removed = []
    for skill_name in SKILLS:
        destination = WORKBUDDY_SKILLS_DIR / skill_name
        if destination.exists():
            shutil.rmtree(destination)
            removed.append(destination)
    return removed


def stop_existing_mcp() -> int:
    pgrep = shutil.which("pgrep")
    if not pgrep:
        return 0
    result = subprocess.run([pgrep, "-f", str(MCP)], text=True, capture_output=True, check=False)
    stopped = 0
    for line in result.stdout.splitlines():
        try:
            pid = int(line.strip())
        except ValueError:
            continue
        if pid == os.getpid():
            continue
        try:
            os.kill(pid, signal.SIGTERM)
            stopped += 1
        except ProcessLookupError:
            pass
    return stopped


def report_r_status() -> None:
    rscript = find_rscript()
    if not rscript:
        print("\nR was not found. The MCP server installs fine, but analyses need R:")
        print("  brew install r")
        print(f"  cd {ROOT} && Rscript r/bootstrap.R")
        return
    print(f"\nFound R: {rscript}")
    probe = "pkgs <- commandArgs(trailingOnly=TRUE); cat(sum(sapply(pkgs, requireNamespace, quietly=TRUE)), length(pkgs))"
    packages = ["jsonlite", "ggplot2", "DESeq2", "limma", "clusterProfiler", "org.Hs.eg.db", "pheatmap", "glmnet", "survival"]
    result = subprocess.run(
        [rscript, "--vanilla", "-e", probe, "--args", *packages],
        capture_output=True, text=True, timeout=180, check=False,
    )
    counts = (result.stdout or "").strip().split()
    if len(counts) == 2 and counts[0] == counts[1]:
        print(f"R packages: {counts[0]}/{counts[1]} core packages present.")
    else:
        have = counts[0] if counts else "?"
        total = counts[1] if len(counts) > 1 else str(len(packages))
        print(f"R packages: {have}/{total} core packages present.")
        print(f"  Install the rest with: cd {ROOT} && Rscript r/bootstrap.R")


def update_repo() -> None:
    if (ROOT / ".git").exists():
        subprocess.run(["git", "pull", "--ff-only"], cwd=str(ROOT), check=True)


def install_workbuddy(_: argparse.Namespace) -> None:
    ensure_runtime()
    config_path = install_mcp_config()
    skill_paths = install_skills()
    stopped = stop_existing_mcp()
    print(f"Installed Omics MCP for WorkBuddy: {config_path}")
    print(f"Installed {len(skill_paths)} Omics skill(s) into: {WORKBUDDY_SKILLS_DIR}")
    if stopped:
        print(f"Stopped {stopped} existing Omics MCP process(es).")
    report_r_status()
    print("\nOpen WorkBuddy -> Connectors -> Custom Connector -> MCP management, then trust/enable omics if prompted.")


def update_workbuddy(_: argparse.Namespace) -> None:
    ensure_runtime()
    update_repo()
    config_path = install_mcp_config()
    skill_paths = install_skills()
    stopped = stop_existing_mcp()
    print(f"Updated Omics MCP for WorkBuddy: {config_path}")
    print(f"Updated {len(skill_paths)} Omics skill(s) into: {WORKBUDDY_SKILLS_DIR}")
    if stopped:
        print(f"Stopped {stopped} existing Omics MCP process(es).")
    report_r_status()
    print(f"\nOmics outputs were kept at: {OUTPUT_DIR}")


def uninstall_workbuddy(_: argparse.Namespace) -> None:
    removed_config = uninstall_mcp_config()
    removed_skills = uninstall_skills()
    stopped = stop_existing_mcp()
    if removed_config:
        print(f"Removed Omics MCP from WorkBuddy: {WORKBUDDY_CONFIG}")
    if removed_skills:
        print(f"Removed {len(removed_skills)} Omics skill(s) from: {WORKBUDDY_SKILLS_DIR}")
    if stopped:
        print(f"Stopped {stopped} Omics MCP process(es).")
    if not removed_config and not removed_skills:
        print("Omics Kit was not installed for WorkBuddy.")
    print(f"Omics outputs and installed R packages were kept. Outputs: {OUTPUT_DIR}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Install Omics Kit Skill + MCP for WorkBuddy.")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command, handler in [
        ("install", install_workbuddy),
        ("update", update_workbuddy),
        ("uninstall", uninstall_workbuddy),
    ]:
        command_parser = subparsers.add_parser(command)
        targets = command_parser.add_subparsers(dest="target", required=True)
        workbuddy = targets.add_parser("workbuddy")
        workbuddy.set_defaults(func=handler)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
