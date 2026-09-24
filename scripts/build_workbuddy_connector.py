import json
import re
import tomllib
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONNECTOR = ROOT / "workbuddy-connector"


def build():
    meta = json.loads((CONNECTOR / "connector-meta.json").read_text(encoding="utf-8"))
    required = {"name", "name_en", "description", "description_zh", "description_en", "source", "type", "version", "minWorkbuddyVersion", "examples_zh", "examples_en"}
    if required - meta.keys():
        raise ValueError(f"Missing connector metadata: {sorted(required - meta.keys())}")
    for key in required - {"examples_zh", "examples_en"}:
        if not isinstance(meta[key], str) or not meta[key].strip():
            raise ValueError(f"Metadata {key} must be a non-empty string")
    if meta["source"] != "monadomics" or meta["type"] != "cli":
        raise ValueError("Expected the monadomics CLI connector")
    version = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["version"]
    if meta["version"] != version:
        raise ValueError("Python and connector versions differ")
    if tuple(map(int, meta["minWorkbuddyVersion"].split("."))) < (5, 0, 0):
        raise ValueError("Python runtime requires WorkBuddy 5.0.0 or newer")
    for key in ("examples_zh", "examples_en"):
        if not isinstance(meta[key], list) or not 2 <= len(meta[key]) <= 5 or not all(isinstance(item, str) and item.strip() for item in meta[key]):
            raise ValueError(f"{key} needs 2-5 non-empty examples")
    config = json.loads((CONNECTOR / "cli.json").read_text(encoding="utf-8"))
    if config.get("runtime") != {"type": "python", "version": "3.12"}:
        raise ValueError("Expected the declared Python 3.12 runtime")
    for platform in ("darwin", "linux", "win32"):
        if config.get("init", {}).get(platform) != f"python -m pip install --upgrade monadomics=={version}":
            raise ValueError(f"Unexpected install command for {platform}")
    if (CONNECTOR / "mcp.json").exists():
        raise ValueError("Do not mix CLI and MCP in this connector")
    icon = ET.parse(CONNECTOR / "icon.svg").getroot()
    if icon.tag != "{http://www.w3.org/2000/svg}svg":
        raise ValueError("Invalid SVG icon")
    skills = list((CONNECTOR / "skills").glob("*/SKILL.md"))
    if len(skills) != 1 or skills[0].parent.name != "monadomics-analysis":
        raise ValueError("Expected exactly the monadomics-analysis main Skill")
    skill = skills[0]
    content = skill.read_text(encoding="utf-8")
    if not content.startswith("---\n"):
        raise ValueError("Skill needs YAML frontmatter")
    fields = dict(line.split(":", 1) for line in content.split("---\n", 2)[1].splitlines() if ":" in line)
    for key in ("name", "description", "description_zh", "description_en", "version", "author"):
        if not fields.get(key, "").strip():
            raise ValueError(f"Missing Skill field: {key}")
    if fields["version"].strip() != version or fields["name"].strip() != skill.parent.name:
        raise ValueError("Skill name/version mismatch")
    for path in CONNECTOR.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlinks cannot be packaged: {path}")
        if path.suffix not in {".json", ".md", ".svg"} or not path.is_file():
            continue
        content = path.read_text(encoding="utf-8")
        if re.search(r"/Users/[^/\s]+/|[A-Z]:\\\\Users\\\\", content):
            raise ValueError(f"Developer-specific absolute path: {path}")
        for target in re.findall(r"@(references/[^\s)]+\.md)", content):
            if not (skill.parent / target).is_file():
                raise ValueError(f"Broken Skill reference: {target}")
    output = ROOT / "dist" / f"monadomics-workbuddy-{version}.zip"
    output.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.write(ROOT / "LICENSE", "LICENSE")
        for path in sorted(CONNECTOR.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(CONNECTOR))
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:
            raise ValueError("Corrupt connector ZIP")
    return output


if __name__ == "__main__":
    print(build())
