import argparse
import csv
import json
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Exercise an installed MonadOmics package with synthetic data.")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    inputs = output / "inputs"
    responses = output / "responses"
    inputs.mkdir(parents=True, exist_ok=True)
    responses.mkdir(exist_ok=True)
    for fixture in (Path(__file__).parent / "fixtures").glob("*.csv"):
        shutil.copy2(fixture, inputs / fixture.name)
    report = {"ok": False, "data": "Synthetic fixed-seed fixtures; not biological or clinical evidence.", "platform": platform.platform(), "python": sys.version, "steps": []}
    report_path = output / "acceptance.json"

    def run(name, command, params=None):
        invocation = [sys.executable, "-m", "monadomics", *command]
        if params is not None:
            params_file = inputs / f"{name}.json"
            params_file.write_text(json.dumps(params, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            invocation += ["--params", str(params_file), "--output-dir", str(output / "artifacts")]
        started = time.monotonic()
        result = subprocess.run(invocation, cwd=output, text=True, capture_output=True, timeout=960, check=False)
        (responses / f"{name}.json").write_text(result.stdout, encoding="utf-8")
        (responses / f"{name}.stderr.txt").write_text(result.stderr, encoding="utf-8")
        payload = json.loads(result.stdout)
        report["steps"].append({"name": name, "exit_code": result.returncode, "ok": payload.get("ok", False), "seconds": round(time.monotonic() - started, 3)})
        report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        if result.returncode or not payload.get("ok"):
            raise RuntimeError(f"{name} failed: {payload}")
        print(f"PASS {name}", flush=True)
        return payload

    with (inputs / "counts.csv").open(encoding="utf-8", newline="") as handle:
        source_genes = {row["gene"] for row in csv.DictReader(handle)}
    doctor = run("doctor", ["doctor", "--group", "core"])
    report["r_version"] = doctor["r_version"]
    matrix = {"matrix_path": str(inputs / "counts.csv"), "matrix_type": "counts", "coldata_path": str(inputs / "coldata.csv")}
    run("pca", ["plot"], {**matrix, "type": "pca", "colour_column": "group", "output_name": "sample_pca"})
    deg = run("deg", ["deg"], {**matrix, "method": "deseq2", "group_column": "group", "treat": "disease", "control": "control", "output_name": "disease_vs_control"})
    with Path(deg["result_table"]).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["gene"] for row in rows} == source_genes
    assert deg["significant"]["up"] > 0 and deg["significant"]["down"] > 0
    assert deg["comparison"] == "disease vs control"
    selected = [row["gene"] for row in deg["top_genes"]][:15]
    assert set(selected) <= source_genes
    run("volcano", ["plot"], {"type": "volcano", "deg_path": deg["result_table"], "output_name": "volcano"})
    run("heatmap", ["plot"], {**matrix, "type": "heatmap", "genes": selected, "annotation_columns": ["group"], "output_name": "heatmap"})
    up = [row["gene"] for row in rows if row["direction"] == "up"]
    enrichment = run("go", ["enrich"], {"method": "go", "species": "human", "id_type": "SYMBOL", "genes": up, "universe": sorted(source_genes), "output_name": "go_up"})
    assert enrichment["genes_mapped"] <= len(up)
    survival = run("lasso_cox", ["survival"], {"method": "lasso_cox", "data_path": str(inputs / "survival.csv"), "time": "time", "event": "event", "predictors": [f"gene{i}" for i in range(1, 9)], "output_name": "lasso_cox", "seed": 42})
    assert survival["n"] == 200 and survival["selected_variables"]
    assert any("optimistic" in warning for warning in survival["warnings"])
    artifacts = [path for path in (output / "artifacts").rglob("*") if path.is_file()]
    assert artifacts and all(path.stat().st_size for path in artifacts)
    for name in ("pca", "volcano", "heatmap", "go", "lasso_cox"):
        payload = json.loads((responses / f"{name}.json").read_text())
        for key in ("figure", "cv_figure", "km_figure"):
            if payload.get(key):
                for field in ("png", "svg"):
                    path = Path(payload[key][field])
                    assert path.is_file() and path.stat().st_size > 100
    report.update(ok=True, gene_count=len(rows), significant=deg["significant"], survival_n=survival["n"], artifact_count=len(artifacts), artifacts=[str(path.relative_to(output)) for path in sorted(artifacts)], not_exercised=["Fresh R/OS dependency installation", "KEGG and other external resource retrieval", "WorkBuddy marketplace installation", "Linux and Windows runtime execution"])
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(report_path)


if __name__ == "__main__":
    main()
