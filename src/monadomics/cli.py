import argparse
import json
import subprocess
import sys
from pathlib import Path

from . import __version__, backend


ANALYSES = ("deg", "enrich", "plot", "survival")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="monadomics", description="Local R-backed bioinformatics analysis.")
    parser.add_argument("--version", action="version", version=f"monadomics {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    doctor_parser = commands.add_parser("doctor", help="Check R and required packages without installing anything.")
    doctor_parser.add_argument("--group", choices=["all", *backend.R_PACKAGE_GROUPS], default="all")
    commands.add_parser("capabilities", help="List the 21 supported capabilities.")
    schema_parser = commands.add_parser("schema", help="Show an analysis command's JSON input schema.")
    schema_parser.add_argument("analysis", choices=ANALYSES)
    setup_parser = commands.add_parser("setup-r", help="Install R packages for a selected analysis group.")
    setup_parser.add_argument("group", choices=["all", *backend.R_PACKAGE_GROUPS])
    for name in ANALYSES:
        analysis = commands.add_parser(name, help=next(item["description"] for item in backend.COMMANDS if item["name"] == name))
        analysis.add_argument("--params", required=True, help="UTF-8 JSON file, or - for standard input. Relative data paths use the current working directory.")
        analysis.add_argument("--output-dir", type=Path, help="Artifact directory; otherwise OMICS_OUTPUT_DIR or ~/.workbuddy/workspace/omics.")
        analysis.add_argument("--timeout", type=int, default=600 if name == "plot" else 900, help="Maximum R execution time in seconds.")
    args = parser.parse_args(argv)

    try:
        if args.command == "doctor":
            result = backend.doctor({"group": args.group})
            ok = result.get("ready", False)
        elif args.command == "capabilities":
            result = backend.feature_menu({})
            ok = True
        elif args.command == "schema":
            result = next(item for item in backend.COMMANDS if item["name"] == args.analysis)
            ok = True
        elif args.command == "setup-r":
            packages = sorted({package for group in backend.R_PACKAGE_GROUPS.values() for package in group}) if args.group == "all" else sorted(set(backend.R_PACKAGE_GROUPS["core"] + backend.R_PACKAGE_GROUPS[args.group]))
            process = subprocess.run(
                [backend._find_rscript(), "--vanilla", str(backend.R_DIR / "bootstrap.R"), *packages],
                stdout=sys.stderr, stderr=sys.stderr, check=False,
            )
            if process.returncode:
                raise RuntimeError(f"R package installation failed (exit {process.returncode}); see stderr.")
            result = backend.doctor({"group": args.group})
            ok = result.get("ready", False)
        else:
            if args.timeout <= 0:
                raise ValueError("--timeout must be a positive number of seconds.")
            params = json.loads(sys.stdin.read() if args.params == "-" else Path(args.params).expanduser().read_text(encoding="utf-8-sig"))
            if not isinstance(params, dict):
                raise ValueError("--params must contain a JSON object.")
            schema = next(item["inputSchema"] for item in backend.COMMANDS if item["name"] == args.command)
            missing = [key for key in schema.get("required", []) if key not in params]
            if missing:
                raise ValueError(f"Missing required parameters: {', '.join(missing)}")
            unknown = sorted(set(params) - set(schema["properties"]))
            if unknown:
                raise ValueError(f"Unknown parameters: {', '.join(unknown)}. See monadomics schema {args.command}.")
            for key, value in params.items():
                allowed = schema["properties"][key].get("enum")
                if allowed and value not in allowed:
                    raise ValueError(f"{key} must be one of: {', '.join(map(str, allowed))}")
                if key.endswith("_path"):
                    if not isinstance(value, str):
                        raise ValueError(f"{key} must be a file path string.")
                    params[key] = str(Path(value).expanduser().resolve())
            params["timeout"] = args.timeout
            backend.OUTPUT_DIR = (args.output_dir or backend.OUTPUT_DIR).expanduser().resolve()
            result = backend.HANDLERS[args.command](params)
            ok = True
        print(json.dumps({"ok": bool(ok), **result}, ensure_ascii=False, allow_nan=False))
        return 0 if ok else 1
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(json.dumps({"ok": False, "error": type(exc).__name__, "message": str(exc)}, ensure_ascii=False))
        return 1
