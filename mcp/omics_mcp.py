#!/usr/bin/env python3
"""Deterministic R-backed omics analysis tools exposed over MCP."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
R_DIR = ROOT / "r"
OUTPUT_DIR = Path(os.getenv("OMICS_OUTPUT_DIR", "~/.workbuddy/workspace/omics")).expanduser()
MESSAGE_MODE = "headers"
PROTOCOL_VERSION = "2024-11-05"
SERVER_VERSION = "0.1.0"

R_PACKAGE_GROUPS: dict[str, list[str]] = {
    "core": ["jsonlite", "ggplot2"],
    "deg": ["DESeq2", "edgeR", "limma"],
    "enrich": ["clusterProfiler", "enrichplot", "org.Hs.eg.db", "org.Mm.eg.db", "ReactomePA", "GSVA", "msigdbr"],
    "plot": ["pheatmap", "ggrepel", "ggvenn", "RColorBrewer"],
    "survival": ["survival", "survminer", "glmnet", "timeROC", "rms"],
}


TOOLS: list[dict[str, Any]] = [
    {
        "name": "omics_env",
        "description": "Report the R runtime and which omics R packages are installed, with the exact command to install anything missing. Call this before the first analysis in a session.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "group": {"type": "string", "enum": ["all", "deg", "enrich", "plot", "survival"], "default": "all"},
            },
        },
    },
    {
        "name": "omics_deg",
        "description": "Differential expression with DESeq2, edgeR, or limma. DESeq2/edgeR require raw integer counts; use limma for microarray or already-normalised data. Returns the full result table plus the significant subset.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "method": {"type": "string", "enum": ["deseq2", "edger", "limma"], "default": "deseq2"},
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
                "output_name": {"type": "string"},
            },
            "required": ["group_column", "treat", "control"],
        },
    },
    {
        "name": "omics_enrich",
        "description": "Functional enrichment: GO/KEGG/Reactome over-representation, GSEA on a ranked list, or GSVA/ssGSEA per-sample pathway scores. Species is required because gene identifiers differ between organisms.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "method": {"type": "string", "enum": ["go", "kegg", "reactome", "gsea", "gsva", "ssgsea"], "default": "go"},
                "species": {"type": "string", "enum": ["human", "mouse"]},
                "genes": {"type": "array", "items": {"type": "string"}},
                "universe": {"type": "array", "items": {"type": "string"}},
                "ontology": {"type": "string", "enum": ["BP", "CC", "MF", "ALL"], "default": "BP"},
                "ranked_path": {"type": "string"},
                "ranked": {"type": "array", "items": {"type": "object"}},
                "gene_column": {"type": "string"},
                "metric_column": {"type": "string", "default": "log2FoldChange"},
                "matrix_path": {"type": "string"},
                "matrix": {"type": "array", "items": {"type": "object"}},
                "gene_sets": {"type": "object"},
                "pvalue": {"type": "number", "default": 0.05},
                "qvalue": {"type": "number", "default": 0.2},
                "top_n": {"type": "integer", "default": 10},
                "output_name": {"type": "string"},
            },
            "required": ["species"],
        },
    },
    {
        "name": "omics_plot",
        "description": "Publication figures for expression analyses: volcano, heatmap, Venn, or PCA. Produces 300 dpi PNG plus editable SVG, and for Venn the per-region gene membership table.",
        "inputSchema": {
            "type": "object",
            "properties": {
                "type": {"type": "string", "enum": ["volcano", "heatmap", "venn", "pca"]},
                "deg_path": {"type": "string"},
                "deg": {"type": "array", "items": {"type": "object"}},
                "matrix_path": {"type": "string"},
                "matrix": {"type": "array", "items": {"type": "object"}},
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
        "name": "omics_survival",
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
                "groups": {"type": "integer"},
                "width": {"type": "number"},
                "height": {"type": "number"},
                "output_name": {"type": "string"},
            },
            "required": ["method"],
        },
    },
    {
        "name": "omics_feature_menu",
        "description": "Return the omics capability menu and routing hints.",
        "inputSchema": {
            "type": "object",
            "properties": {"context": {"type": "string", "default": "general"}},
        },
    },
]


class MessageParseError(Exception):
    pass


def _find_rscript() -> str:
    configured = os.getenv("OMICS_RSCRIPT", "").strip()
    if configured:
        return configured
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
                cwd=str(ROOT),
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
                return result
        stderr = (process.stderr or "").strip()
        raise RuntimeError(f"R script '{script}' failed (exit {process.returncode}): {stderr[-1500:] or 'no output'}")


def env_tool(args: dict[str, Any]) -> dict[str, Any]:
    group = str(args.get("group", "all") or "all")
    if group == "all":
        packages = sorted({name for names in R_PACKAGE_GROUPS.values() for name in names})
    elif group in R_PACKAGE_GROUPS:
        packages = sorted(set(R_PACKAGE_GROUPS["core"]) | set(R_PACKAGE_GROUPS[group]))
    else:
        raise ValueError(f"Unknown group: {group}. Use all, deg, enrich, plot, or survival.")

    try:
        rscript = _find_rscript()
    except RuntimeError as exc:
        return {
            "r_available": False,
            "message": str(exc),
            "install_r": "brew install r",
            "groups": R_PACKAGE_GROUPS,
        }

    probe = (
        "pkgs <- commandArgs(trailingOnly=TRUE); "
        "cat(R.version.string, '\\n'); "
        "for (p in pkgs) cat(p, as.integer(requireNamespace(p, quietly=TRUE)), '\\n')"
    )
    process = subprocess.run(
        [rscript, "--vanilla", "-e", probe, "--args", *packages],
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )
    lines = [line.strip() for line in (process.stdout or "").splitlines() if line.strip()]
    version = lines[0] if lines else "unknown"
    installed: list[str] = []
    missing: list[str] = []
    for line in lines[1:]:
        parts = line.rsplit(" ", 1)
        if len(parts) != 2:
            continue
        (installed if parts[1] == "1" else missing).append(parts[0])

    return {
        "r_available": True,
        "rscript": rscript,
        "r_version": version,
        "group": group,
        "installed": installed,
        "missing": missing,
        "ready": not missing,
        "install_command": None if not missing else f"cd {ROOT} && Rscript r/bootstrap.R {' '.join(sorted(set(missing)))}",
        "output_dir": str(OUTPUT_DIR),
    }


def deg_tool(args: dict[str, Any]) -> dict[str, Any]:
    return _run_r("deg.R", args, timeout=int(args.get("timeout", 900)))


def enrich_tool(args: dict[str, Any]) -> dict[str, Any]:
    return _run_r("enrich.R", args, timeout=int(args.get("timeout", 900)))


def plot_tool(args: dict[str, Any]) -> dict[str, Any]:
    return _run_r("plots.R", args, timeout=int(args.get("timeout", 600)))


def survival_tool(args: dict[str, Any]) -> dict[str, Any]:
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
            "clinical_stats": "普通临床表格统计、基础 KM 与单变量 Cox 属于 Scholar 的 sci-stats，不在这里重复。",
            "integrity": "基因符号、通路条目、模型系数只能来自工具输出，不得由模型凭记忆产生。",
        },
    }


CALLS: dict[str, Callable[[dict[str, Any]], Any]] = {
    "omics_env": env_tool,
    "omics_deg": deg_tool,
    "omics_enrich": enrich_tool,
    "omics_plot": plot_tool,
    "omics_survival": survival_tool,
    "omics_feature_menu": feature_menu,
}


def read_message() -> dict[str, Any] | None:
    global MESSAGE_MODE
    stream = getattr(sys.stdin, "buffer", sys.stdin)
    line = stream.readline()
    if not line:
        return None
    if isinstance(line, bytes):
        line_text = line.decode("utf-8")
    else:
        line_text = line
    stripped = line_text.strip()
    if not stripped:
        return read_message()
    if stripped.lower().startswith("content-length:"):
        MESSAGE_MODE = "headers"
        try:
            length = int(stripped.split(":", 1)[1].strip())
        except ValueError as exc:
            raise MessageParseError(f"Invalid Content-Length: {stripped}") from exc
        while True:
            header = stream.readline()
            if isinstance(header, bytes):
                header_text = header.decode("utf-8")
            else:
                header_text = header
            if not header or header_text.strip() == "":
                break
        body = stream.read(length)
        if not body:
            return None
        if isinstance(body, bytes):
            body = body.decode("utf-8")
        try:
            return json.loads(body)
        except json.JSONDecodeError as exc:
            raise MessageParseError(str(exc)) from exc
    MESSAGE_MODE = "lines"
    try:
        return json.loads(stripped)
    except json.JSONDecodeError as exc:
        raise MessageParseError(str(exc)) from exc


def write_message(message: dict[str, Any]) -> None:
    body = json.dumps(message, ensure_ascii=False)
    if MESSAGE_MODE == "headers":
        sys.stdout.write(f"Content-Length: {len(body.encode('utf-8'))}\r\n\r\n{body}")
    else:
        sys.stdout.write(body + "\n")
    sys.stdout.flush()


def result_response(request_id: Any, result: Any) -> dict[str, Any]:
    return {"jsonrpc": "2.0", "id": request_id, "result": result}


def error_response(request_id: Any, code: int, message: str) -> dict[str, Any]:
    return {"jsonrpc": "2.0", "id": request_id, "error": {"code": code, "message": message}}


def tool_content(data: Any, is_error: bool = False) -> dict[str, Any]:
    text = data if isinstance(data, str) else json.dumps(data, ensure_ascii=False, indent=2)
    return {"content": [{"type": "text", "text": text}], "isError": is_error}


def tool_error(error_type: str, message: str, retryable: bool = False) -> dict[str, Any]:
    return tool_content({"error": error_type, "message": message, "retryable": retryable}, is_error=True)


def handle(message: dict[str, Any]) -> dict[str, Any] | None:
    method = message.get("method")
    request_id = message.get("id")
    params = message.get("params") or {}

    if method == "initialize":
        return result_response(
            request_id,
            {
                "protocolVersion": PROTOCOL_VERSION,
                "capabilities": {"tools": {}},
                "serverInfo": {"name": "omics", "version": SERVER_VERSION},
            },
        )
    if method in {"notifications/initialized", "initialized"}:
        return None
    if method == "tools/list":
        return result_response(request_id, {"tools": TOOLS})
    if method == "tools/call":
        name = params.get("name")
        arguments = params.get("arguments") or {}
        handler = CALLS.get(name)
        if handler is None:
            return result_response(request_id, tool_error("unknown_tool", f"Unknown tool: {name}"))
        try:
            return result_response(request_id, tool_content(handler(arguments)))
        except ValueError as exc:
            return result_response(request_id, tool_error("invalid_input", str(exc)))
        except RuntimeError as exc:
            return result_response(request_id, tool_error("analysis_failed", str(exc)))
        except Exception as exc:  # noqa: BLE001 - surface unexpected failures to the client
            return result_response(request_id, tool_error("internal_error", f"{type(exc).__name__}: {exc}"))
    if request_id is None:
        return None
    return error_response(request_id, -32601, f"Method not found: {method}")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    while True:
        try:
            message = read_message()
        except MessageParseError as exc:
            write_message(error_response(None, -32700, f"Parse error: {exc}"))
            continue
        if message is None:
            return 0
        response = handle(message)
        if response is not None:
            write_message(response)


if __name__ == "__main__":
    raise SystemExit(main())
