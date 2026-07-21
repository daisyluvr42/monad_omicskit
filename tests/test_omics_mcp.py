#!/usr/bin/env python3
"""Tests for the Omics MCP layer.

Protocol and validation tests always run. Tests that need R and its
Bioconductor packages skip cleanly when those are absent, so the suite is
useful before the 20-minute dependency install completes.
"""

from __future__ import annotations

import json
import io
import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
sys.path.insert(0, str(ROOT / "mcp"))

import omics_mcp  # noqa: E402


def r_packages_present(*packages: str) -> bool:
    try:
        rscript = omics_mcp._find_rscript()
    except RuntimeError:
        return False
    probe = f"cat(all(sapply({omics_mcp._r_vector(list(packages))}, requireNamespace, quietly=TRUE)))"
    result = subprocess.run(
        [rscript, "--vanilla", "-e", probe],
        capture_output=True, text=True, timeout=180, check=False,
    )
    return "TRUE" in (result.stdout or "")


HAS_R = r_packages_present("jsonlite")
HAS_DEG = HAS_R and r_packages_present("DESeq2", "limma")
HAS_PLOT = HAS_R and r_packages_present("ggplot2", "pheatmap")
HAS_SURVIVAL = HAS_R and r_packages_present("glmnet", "survival")


class ProtocolTests(unittest.TestCase):
    def test_every_tool_has_a_handler(self) -> None:
        declared = {tool["name"] for tool in omics_mcp.TOOLS}
        self.assertEqual(declared, set(omics_mcp.CALLS))

    def test_tool_schemas_are_objects(self) -> None:
        for tool in omics_mcp.TOOLS:
            self.assertEqual(tool["inputSchema"]["type"], "object", tool["name"])
            self.assertTrue(tool["description"].strip(), tool["name"])

    def test_feature_menu_counts_match(self) -> None:
        menu = omics_mcp.feature_menu({})
        total = sum(len(group["items"]) for group in menu["groups"])
        self.assertEqual(menu["count"], total)
        self.assertGreaterEqual(total, 20)

    def test_unknown_tool_is_reported(self) -> None:
        response = omics_mcp.handle(
            {"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": "nope", "arguments": {}}}
        )
        self.assertTrue(response["result"]["isError"])

    def test_initialize_reports_server_name(self) -> None:
        response = omics_mcp.handle({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {}})
        self.assertEqual(response["result"]["serverInfo"]["name"], "omics")

    def test_content_length_uses_utf8_bytes(self) -> None:
        body = json.dumps(
            {
                "jsonrpc": "2.0",
                "id": 7,
                "method": "tools/call",
                "params": {
                    "name": "omics_feature_menu",
                    "arguments": {"context": "\u7ec4\u5b66"},
                },
            },
            ensure_ascii=False,
        )
        payload = (
            f"Content-Length: {len(body.encode('utf-8'))}\r\n\r\n".encode("utf-8")
            + body.encode("utf-8")
        )
        stdin = io.TextIOWrapper(io.BytesIO(payload), encoding="utf-8")
        with mock.patch.object(sys, "stdin", stdin):
            omics_mcp.MESSAGE_MODE = "headers"
            message = omics_mcp.read_message()
        self.assertEqual(message["params"]["arguments"]["context"], "\u7ec4\u5b66")


class EnvTests(unittest.TestCase):
    def test_env_reports_availability_either_way(self) -> None:
        result = omics_mcp.env_tool({"group": "deg"})
        self.assertIn("r_available", result)
        if result["r_available"]:
            self.assertIn("missing", result)
            self.assertIn("installed", result)
        else:
            self.assertIn("install_r", result)

    def test_unknown_group_rejected(self) -> None:
        with self.assertRaises(ValueError):
            omics_mcp.env_tool({"group": "bogus"})


@unittest.skipUnless(HAS_DEG, "R with DESeq2/limma not available")
class DegTests(unittest.TestCase):
    def base_args(self, **overrides) -> dict:
        args = {
            "matrix_path": str(FIXTURES / "counts.csv"),
            "coldata_path": str(FIXTURES / "coldata.csv"),
            "group_column": "group",
            "treat": "disease",
            "control": "control",
            "output_name": "test_deg",
        }
        args.update(overrides)
        return args

    def test_deseq2_finds_the_planted_signal(self) -> None:
        result = omics_mcp.deg_tool(self.base_args(method="deseq2"))
        self.assertEqual(result["comparison"], "disease vs control")
        self.assertGreater(result["significant"]["up"], 0)
        self.assertGreater(result["significant"]["down"], 0)
        self.assertTrue(Path(result["result_table"]).exists())
        genes = {row["gene"] for row in result["top_genes"]}
        self.assertTrue(genes & {"IL6", "TNF", "IL1B"}, "planted up-regulated genes should surface")

    def test_covariates_enter_the_design(self) -> None:
        result = omics_mcp.deg_tool(self.base_args(method="deseq2", covariates=["batch"]))
        self.assertIn("batch", result["design"])

    def test_log_transformed_input_is_rejected(self) -> None:
        # The guard that stops the most common fatal bioinformatics mistake.
        with self.assertRaises(RuntimeError) as ctx:
            omics_mcp.deg_tool(self.base_args(method="deseq2", matrix_path=str(FIXTURES / "logged.csv")))
        self.assertIn("count", str(ctx.exception).lower())

    def test_unknown_group_level_is_reported(self) -> None:
        with self.assertRaises(RuntimeError) as ctx:
            omics_mcp.deg_tool(self.base_args(treat="nonexistent"))
        self.assertIn("nonexistent", str(ctx.exception))


@unittest.skipUnless(HAS_PLOT, "R with ggplot2/pheatmap not available")
class PlotTests(unittest.TestCase):
    def test_pca_reports_variance_explained(self) -> None:
        result = omics_mcp.plot_tool({
            "type": "pca",
            "matrix_path": str(FIXTURES / "counts.csv"),
            "coldata_path": str(FIXTURES / "coldata.csv"),
            "colour_column": "group",
            "output_name": "test_pca",
        })
        self.assertTrue(Path(result["figure"]["png"]).exists())
        self.assertTrue(Path(result["figure"]["svg"]).exists())
        self.assertGreater(result["variance_explained"]["PC1"], 0)

    def test_venn_returns_region_membership(self) -> None:
        result = omics_mcp.plot_tool({
            "type": "venn",
            "sets": {"A": ["TP53", "EGFR", "MYC"], "B": ["MYC", "PTEN"]},
            "output_name": "test_venn",
        })
        self.assertEqual(result["intersection_all"], ["MYC"])
        self.assertTrue(Path(result["membership_table"]).exists())

    def test_venn_rejects_too_many_sets(self) -> None:
        with self.assertRaises(RuntimeError):
            omics_mcp.plot_tool({
                "type": "venn",
                "sets": {name: ["G1"] for name in "ABCDE"},
            })


@unittest.skipUnless(HAS_SURVIVAL, "R with glmnet/survival not available")
class SurvivalTests(unittest.TestCase):
    def test_lasso_cox_selects_variables_and_warns(self) -> None:
        result = omics_mcp.survival_tool({
            "method": "lasso_cox",
            "data_path": str(FIXTURES / "survival.csv"),
            "time": "time",
            "event": "event",
            "predictors": [f"gene{i + 1}" for i in range(8)],
            "output_name": "test_lasso",
        })
        self.assertGreater(len(result["selected_variables"]), 0)
        self.assertGreater(result["concordance_index"], 0.5)
        self.assertTrue(any("optimistic" in w for w in result["warnings"]))
        self.assertTrue(Path(result["coefficient_table"]).exists())

    def test_too_few_events_is_rejected(self) -> None:
        rows = [{"time": 10 + i, "event": 1 if i < 2 else 0, "g1": i, "g2": -i} for i in range(20)]
        with self.assertRaises(RuntimeError) as ctx:
            omics_mcp.survival_tool({
                "method": "lasso_cox", "data": rows,
                "time": "time", "event": "event", "predictors": ["g1", "g2"],
            })
        self.assertIn("event", str(ctx.exception).lower())


if __name__ == "__main__":
    unittest.main()
