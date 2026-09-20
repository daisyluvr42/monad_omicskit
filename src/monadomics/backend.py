#!/usr/bin/env python3
"""Deterministic R-backed omics analyses shared by CLI commands."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parent
R_DIR = ROOT / "r"
OUTPUT_DIR = Path(os.getenv("OMICS_OUTPUT_DIR", "~/.workbuddy/workspace/omics")).expanduser()

R_PACKAGE_GROUPS: dict[str, list[str]] = {
    "core": ["jsonlite", "ggplot2", "svglite"],
    "deg": ["DESeq2", "edgeR", "limma"],
    "enrich": ["clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "ReactomePA", "GSVA"],
    "plot": ["edgeR", "pheatmap", "ggrepel", "ggvenn", "svglite"],
    "survival": ["survival", "glmnet", "timeROC", "rms"],
}


COMMANDS: list[dict[str, Any]] = [
    {
        "name": "doctor",
        "description": "Report the R runtime and which omics R packages are installed, with the exact command to install anything missing. Call this before the first analysis in a session.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "group": {"type": "string", "enum": ["all", "core", "deg", "enrich", "plot", "survival"], "default": "all"},
            },
        },
    },
    {
        "name": "deg",
        "description": "Two-group differential expression with DESeq2, edgeR, or limma; other groups are excluded before fitting. DESeq2/edgeR require raw integer counts. Returns full and significant tables, actual sample selection and optional sourced gene annotation.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "method": {"type": "string", "enum": ["deseq2", "edger", "limma"], "default": "deseq2"},
                "matrix_type": {"type": "string", "enum": ["counts", "normalized"]},
                "matrix_path": {"type": "string", "description": "Gene-by-sample matrix; gene IDs in column 1."},
                "matrix": {"type": "array", "items": {"type": "object"}},
                "coldata_path": {"type": "string", "description": "Sample metadata; sample IDs in column 1."},
                "coldata": {"type": "array", "items": {"type": "object"}},
                "group_column": {"type": "string"},
                "treat": {"type": "string"},
                "control": {"type": "string"},
                "covariates": {"type": "array", "items": {"type": "string"}},
                "log2fc": {"type": "number", "default": 1},
                "padj": {"type": "number", "default": 0.05},
                "padj_method": {"type": "string", "default": "BH"},
                "voom": {"type": "boolean", "default": False},
                "species": {"type": "string", "enum": ["human", "mouse"], "description": "Optional gene annotation; supply together with id_type."},
                "id_type": {"type": "string", "enum": ["SYMBOL", "ENSEMBL", "ENTREZID"], "description": "Original matrix ID type; required with species for annotation."},
                "output_name": {"type": "string"},
            },
            "required": ["matrix_type", "group_column", "treat", "control"],
        },
    },
    {
        "name": "enrich",
        "description": "Functional enrichment: GO/KEGG/Reactome over-representation, GSEA on a ranked list, or GSVA/ssGSEA per-sample pathway scores. Species is required because gene identifiers differ between organisms.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "method": {"type": "string", "enum": ["go", "kegg", "reactome", "gsea", "gsva", "ssgsea"], "default": "go"},
                "species": {"type": "string", "enum": ["human", "mouse"]},
                "id_type": {"type": "string", "enum": ["SYMBOL", "ENSEMBL", "ENTREZID"]},
                "genes": {"type": "array", "items": {"type": "string"}},
                "universe": {"type": "array", "items": {"type": "string"}},
                "ontology": {"type": "string", "enum": ["BP", "CC", "MF", "ALL"], "default": "BP"},
                "ranked_path": {"type": "string"},
                "ranked": {"type": "array", "items": {"type": "object"}},
                "gene_column": {"type": "string"},
                "metric_column": {"type": "string", "default": "log2FoldChange"},
                "matrix_path": {"type": "string"},
                "matrix": {"type": "array", "items": {"type": "object"}},
                "matrix_type": {"type": "string", "enum": ["counts", "normalized"]},
                "gene_sets": {"type": "object"},
                "pvalue": {"type": "number", "default": 0.05},
                "padj": {"type": "number", "default": 0.05, "description": "GSEA only: BH adjusted-P cutoff, separate from raw pvalue."},
                "qvalue": {"type": "number", "description": "ORA: default 0.2. GSEA: optional additional q-value cutoff."},
                "seed": {"type": "integer", "default": 42, "description": "GSEA random seed."},
                "top_n": {"type": "integer", "default": 10},
                "output_name": {"type": "string"},
            },
            "required": ["species", "id_type"],
        },
    },
    {
        "name": "plot",
        "description": "Publication figures for expression analyses: volcano, heatmap, Venn, or PCA. Produces 300 dpi PNG plus editable SVG, and for Venn the per-region gene membership table.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "type": {"type": "string", "enum": ["volcano", "heatmap", "venn", "pca"]},
                "deg_path": {"type": "string"},
                "deg": {"type": "array", "items": {"type": "object"}},
                "matrix_path": {"type": "string"},
                "matrix": {"type": "array", "items": {"type": "object"}},
                "matrix_type": {"type": "string", "enum": ["counts", "normalized"]},
                "coldata_path": {"type": "string"},
                "coldata": {"type": "array", "items": {"type": "object"}},
                "sets": {"type": "object", "description": "venn: 2-4 named gene vectors."},
                "genes": {"type": "array", "items": {"type": "string"}},
                "gene_column": {"type": "string", "default": "gene"},
                "annotation_columns": {"type": "array", "items": {"type": "string"}},
                "colour_column": {"type": "string"},
                "log2fc": {"type": "number", "default": 1},
                "padj": {"type": "number", "default": 0.05},
                "label_top": {"type": "integer", "default": 10},
                "label_genes": {"type": "array", "items": {"type": "string"}},
                "label_samples": {"type": "boolean", "default": False},
                "scale_rows": {"type": "boolean", "default": True},
                "show_rownames": {"type": "boolean"},
                "cluster_columns": {"type": "boolean", "default": True},
                "top_variable": {"type": "integer", "default": 2000},
                "title": {"type": "string"},
                "width": {"type": "number"},
                "height": {"type": "number"},
                "output_name": {"type": "string"},
            },
            "required": ["type"],
        },
    },
    {
        "name": "survival",
        "description": "Prognostic modelling: LASSO-Cox variable selection with risk score, time-dependent ROC, nomogram, bootstrap calibration, and decision curve analysis. Reports events-per-variable and flags optimistic training-set performance.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "method": {"type": "string", "enum": ["lasso_cox", "timeroc", "nomogram", "calibration", "dca"], "default": "lasso_cox"},
                "data_path": {"type": "string"},
                "data": {"type": "array", "items": {"type": "object"}},
                "time": {"type": "string", "default": "time"},
                "event": {"type": "string", "default": "event", "description": "0/1 coded; 1 = event occurred."},
                "predictors": {"type": "array", "items": {"type": "string"}},
                "id_column": {"type": "string", "description": "Sample identifier column; defaults to the first column."},
                "risk_column": {"type": "string", "default": "risk_score"},
                "alpha": {"type": "number", "default": 1},
                "nfolds": {"type": "integer", "default": 10},
                "lambda": {"type": "string", "enum": ["1se", "min"], "default": "1se"},
                "seed": {"type": "integer", "default": 42},
                "times": {"type": "array", "items": {"type": "number"}},
                "thresholds": {"type": "array", "items": {"type": "number"}},
                "bootstrap": {"type": "integer", "default": 200},
                "groups": {"type": "integer", "description": "Requested number of calibration groups."},
                "width": {"type": "number"},
                "height": {"type": "number"},
                "output_name": {"type": "string"},
            },
            "required": ["method"],
        },
    },
    {
        "name": "capabilities",
        "description": "Return the omics capability menu and routing hints.",
        "inputSchema": {
            "type": "object",
            "properties": {"context": {"type": "string", "default": "general"}},
        },
    },
]



def _find_rscript() -> str:
    configured = os.getenv("OMICS_RSCRIPT", "").strip()
    if configured:
        found = shutil.which(str(Path(configured).expanduser()))
        if not found:
            raise RuntimeError("OMICS_RSCRIPT does not point to an executable Rscript. Correct the path and retry.")
        return found
    found = shutil.which("Rscript")
    if found:
        return found
    # The macOS CRAN build is not on PATH for GUI-launched processes.
    for candidate in [
        Path("/Library/Frameworks/R.framework/Resources/bin/Rscript"),
        Path("/opt/homebrew/bin/Rscript"),
        Path("/usr/local/bin/Rscript"),
    ]:
        if candidate.exists():
            return str(candidate)
    raise RuntimeError(
        "Rscript was not found. Install R (macOS: brew install r), then reopen the client. "
        "If R is installed somewhere unusual, set OMICS_RSCRIPT to its Rscript path."
    )


def _run_r(script: str, payload: dict[str, Any], timeout: int = 900) -> dict[str, Any]:
    rscript = _find_rscript()
    script_path = R_DIR / script
    if not script_path.exists():
        raise RuntimeError(f"Missing R script: {script_path}")
    environment = dict(os.environ)
    environment["OMICS_OUTPUT_DIR"] = str(OUTPUT_DIR)
    environment["OMICS_R_DIR"] = str(R_DIR)
    with tempfile.TemporaryDirectory() as workdir:
        args_path = Path(workdir) / "args.json"
        out_path = Path(workdir) / "result.json"
        args_path.write_text(json.dumps(payload, ensure_ascii=False, default=str), encoding="utf-8")
        try:
            process = subprocess.run(
                [rscript, "--vanilla", str(script_path), str(args_path), str(out_path)],
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=workdir,
                env=environment,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                f"R script '{script}' exceeded {timeout}s. Reduce the input size, or run this step directly in R."
            ) from exc
        # The R wrapper writes a JSON payload even for handled errors; prefer it over stderr noise.
        if out_path.exists():
            try:
                result = json.loads(out_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                result = None
            if isinstance(result, dict):
                if result.get("error"):
                    raise RuntimeError(str(result["error"]))
                if process.returncode == 0:
                    return result
        stderr = (process.stderr or "").strip()
        raise RuntimeError(f"R script '{script}' failed (exit {process.returncode}): {stderr[-1500:] or 'no output'}")


def _r_vector(values: list[str]) -> str:
    """Build an R character vector literal.

    Package names are passed inside the -e expression rather than after --args,
    because `commandArgs(trailingOnly=TRUE)` includes the literal "--args" when R
    is invoked with -e, which silently poisons any check over that vector.
    """
    safe = [name for name in values if re.fullmatch(r"[A-Za-z0-9._]+", name)]
    return "c(" + ",".join(f'"{name}"' for name in safe) + ")"


def doctor(args: dict[str, Any]) -> dict[str, Any]:
    group = str(args.get("group", "all") or "all")
    if group == "all":
        packages = sorted({name for names in R_PACKAGE_GROUPS.values() for name in names})
    elif group in R_PACKAGE_GROUPS:
        packages = sorted(set(R_PACKAGE_GROUPS["core"]) | set(R_PACKAGE_GROUPS[group]))
    else:
        raise ValueError(f"Unknown group: {group}. Use all, core, deg, enrich, plot, or survival.")

    try:
        rscript = _find_rscript()
    except RuntimeError as exc:
        return {
            "r_available": False,
            "message": str(exc),
            "install_r": "https://cran.r-project.org/",
            "groups": R_PACKAGE_GROUPS,
        }

    probe = (
        "if (getRversion() < '4.2.0') stop('R 4.2 or newer is required'); "
        f"pkgs <- {_r_vector(packages)}; "
        "cat(R.version.string, '\\n'); "
        "for (p in pkgs) cat(p, as.integer(requireNamespace(p, quietly=TRUE)), '\\n')"
    )
    process = subprocess.run(
        [rscript, "--vanilla", "-e", probe],
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )
    if process.returncode:
        raise RuntimeError(f"R package check failed (exit {process.returncode}): {process.stderr.strip()[-1500:]}")
    lines = [line.strip() for line in (process.stdout or "").splitlines() if line.strip()]
    version = lines[0] if lines else "unknown"
    installed: list[str] = []
    missing: list[str] = []
    for line in lines[1:]:
        parts = line.rsplit(" ", 1)
        if len(parts) != 2:
            continue
        (installed if parts[1] == "1" else missing).append(parts[0])

    if set(installed + missing) != set(packages):
        raise RuntimeError("R package check returned incomplete results; dependency readiness is unknown.")

    return {
        "r_available": True,
        "rscript": rscript,
        "r_version": version,
        "group": group,
        "installed": installed,
        "missing": missing,
        "ready": not missing,
        "install_command": None if not missing else f"monadomics setup-r {group}",
        "output_dir": str(OUTPUT_DIR),
    }


def deg(args: dict[str, Any]) -> dict[str, Any]:
    return _run_r("deg.R", args, timeout=int(args.get("timeout", 900)))


def enrich(args: dict[str, Any]) -> dict[str, Any]:
    return _run_r("enrich.R", args, timeout=int(args.get("timeout", 900)))


def plot(args: dict[str, Any]) -> dict[str, Any]:
    return _run_r("plots.R", args, timeout=int(args.get("timeout", 600)))


def survival(args: dict[str, Any]) -> dict[str, Any]:
    return _run_r("survival.R", args, timeout=int(args.get("timeout", 900)))


def feature_menu(args: dict[str, Any]) -> dict[str, Any]:
    groups = [
        {
            "module": "数据准备与质控",
            "items": [
                "1. 表达矩阵与样本表核对",
                "2. 数据类型判断（count / TPM / 芯片 / 已 log2）",
                "3. 样本 PCA 与批次结构检查",
            ],
        },
        {
            "module": "差异表达",
            "items": [
                "4. DESeq2 差异分析（原始 count）",
                "5. edgeR 差异分析（原始 count）",
                "6. limma / limma-voom 差异分析（芯片或已标准化数据）",
                "7. 协变量校正与多重检验",
                "8. 火山图",
                "9. 表达热图",
                "10. 基因集 Venn 与区域归属表",
            ],
        },
        {
            "module": "功能富集",
            "items": [
                "11. GO 过表达分析（BP/CC/MF）",
                "12. KEGG 通路富集",
                "13. Reactome 通路富集",
                "14. GSEA（完整排序列表）",
                "15. GSVA / ssGSEA 单样本通路打分",
            ],
        },
        {
            "module": "预后模型",
            "items": [
                "16. LASSO-Cox 变量筛选与风险评分",
                "17. 高低危分组生存曲线",
                "18. 时间依赖 ROC",
                "19. 列线图",
                "20. Bootstrap 校准曲线",
                "21. 决策曲线分析（DCA）",
            ],
        },
    ]
    return {
        "context": args.get("context", "general"),
        "title": "Omics 分析工作台",
        "groups": groups,
        "count": sum(len(group["items"]) for group in groups),
        "routing": {
            "species": "物种必须显式确认，人鼠基因符号不通用。",
            "data_type": "count 走 DESeq2/edgeR；芯片、TPM 或已 log2 数据走 limma。判断错会导致整条结果链失效。",
            "single_cell": "Seurat/CellChat/Monocle 属于小时级重流程，不在工具层运行；生成脚本交用户执行，再把下游结果拿回来分析。",
            "clinical_stats": "普通临床表格统计、基础 KM 与单变量 Cox 应交给独立统计工具；本 kit 不依赖其他连接器。",
            "integrity": "基因符号、通路条目、模型系数只能来自工具输出，不得由模型凭记忆产生。",
        },
    }


HANDLERS: dict[str, Callable[[dict[str, Any]], Any]] = {
    "doctor": doctor,
    "deg": deg,
    "enrich": enrich,
    "plot": plot,
    "survival": survival,
    "capabilities": feature_menu,
}
