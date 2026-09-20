#!/usr/bin/env python3
"""Tests for the Omics analysis backend.

Schema and validation tests always run. Tests that need R and its
Bioconductor packages skip cleanly when those are absent, so the suite is
useful before the 20-minute dependency install completes.
"""

from __future__ import annotations

import json
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).resolve().parent / "fixtures"
sys.path.insert(0, str(ROOT / "src"))

from monadomics import backend


def r_packages_present(*packages: str) -> bool:
    try:
        rscript = backend._find_rscript()
    except RuntimeError:
        return False
    probe = f"cat(all(sapply({backend._r_vector(list(packages))}, requireNamespace, quietly=TRUE)))"
    result = subprocess.run(
        [rscript, "--vanilla", "-e", probe],
        capture_output=True, text=True, timeout=180, check=False,
    )
    return "TRUE" in (result.stdout or "")


HAS_R = r_packages_present("jsonlite")
HAS_DEG = HAS_R and r_packages_present("DESeq2", "limma")
HAS_PLOT = HAS_R and r_packages_present("ggplot2", "pheatmap", "edgeR")
HAS_SURVIVAL = HAS_R and r_packages_present("glmnet", "survival")
HAS_TIMEROC = HAS_R and r_packages_present("timeROC", "survival")
HAS_RMS = HAS_R and r_packages_present("rms", "survival")
HAS_ENRICH = HAS_R and r_packages_present("clusterProfiler", "org.Hs.eg.db", "GSVA")


class SchemaTests(unittest.TestCase):
    def test_every_tool_has_a_handler(self) -> None:
        declared = {tool["name"] for tool in backend.COMMANDS}
        self.assertEqual(declared, set(backend.HANDLERS))

    def test_tool_schemas_are_objects(self) -> None:
        for tool in backend.COMMANDS:
            self.assertEqual(tool["inputSchema"]["type"], "object", tool["name"])
            self.assertTrue(tool["description"].strip(), tool["name"])

    def test_scientific_input_semantics_are_explicit(self) -> None:
        schemas = {tool["name"]: tool["inputSchema"] for tool in backend.COMMANDS}
        self.assertIn("matrix_type", schemas["deg"]["required"])
        self.assertIn("id_type", schemas["enrich"]["required"])
        self.assertIn("matrix_type", schemas["plot"]["properties"])

    def test_feature_menu_counts_match(self) -> None:
        menu = backend.feature_menu({})
        total = sum(len(group["items"]) for group in menu["groups"])
        self.assertEqual(menu["count"], total)
        self.assertGreaterEqual(total, 20)



class EnvTests(unittest.TestCase):
    def test_env_reports_availability_either_way(self) -> None:
        result = backend.doctor({"group": "deg"})
        self.assertIn("r_available", result)
        if result["r_available"]:
            self.assertIn("missing", result)
            self.assertIn("installed", result)
        else:
            self.assertIn("install_r", result)

    def test_unknown_group_rejected(self) -> None:
        with self.assertRaises(ValueError):
            backend.doctor({"group": "bogus"})

    def test_install_command_uses_the_packaged_cli(self) -> None:
        packages = sorted(set(backend.R_PACKAGE_GROUPS["core"] + backend.R_PACKAGE_GROUPS["deg"]))
        probe = mock.Mock(stdout="R version 4.5.0\n" + "".join(f"{package} 0\n" for package in packages), stderr="", returncode=0)
        with (
            mock.patch.object(backend, "_find_rscript", return_value="/path with spaces/Rscript"),
            mock.patch.object(backend.subprocess, "run", return_value=probe),
        ):
            result = backend.doctor({"group": "deg"})
        self.assertEqual(result["install_command"], "monadomics setup-r deg")

    def test_failed_dependency_probe_cannot_report_ready(self) -> None:
        probe = mock.Mock(stdout="", stderr="R failed to start", returncode=1)
        with mock.patch.object(backend, "_find_rscript", return_value="Rscript"), mock.patch.object(backend.subprocess, "run", return_value=probe):
            with self.assertRaisesRegex(RuntimeError, "package check failed"):
                backend.doctor({"group": "deg"})

    def test_incomplete_dependency_probe_cannot_report_ready(self) -> None:
        probe = mock.Mock(stdout="R version 4.6.1\njsonlite 1\n", stderr="", returncode=0)
        with mock.patch.object(backend, "_find_rscript", return_value="Rscript"), mock.patch.object(backend.subprocess, "run", return_value=probe):
            with self.assertRaisesRegex(RuntimeError, "incomplete"):
                backend.doctor({"group": "deg"})


@unittest.skipUnless(HAS_R, "R/jsonlite not available")
class RHelperTests(unittest.TestCase):
    def run_r(self, expression: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                backend._find_rscript(),
                "--vanilla",
                "-e",
                f"source('{backend.R_DIR / 'lib/common.R'}'); {expression}",
            ],
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )

    def test_low_integer_counts_are_valid(self) -> None:
        result = self.run_r("omics_assert_counts(matrix(c(0, 2, 8, 30), nrow=2), 'test'); cat('OK')")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("OK", result.stdout)

    def test_duplicate_count_rows_require_explicit_resolution(self) -> None:
        expression = (
            "params <- list(matrix=data.frame(gene=c('G1','G1'), S1=c(1,2), S2=c(3,4), "
            "check.names=FALSE)); omics_read_matrix(params, matrix_type='counts')"
        )
        result = self.run_r(expression)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unique", result.stderr)


@unittest.skipUnless(HAS_DEG, "R with DESeq2/limma not available")
class DegTests(unittest.TestCase):
    def base_args(self, **overrides) -> dict:
        args = {
            "matrix_path": str(FIXTURES / "counts.csv"),
            "matrix_type": "counts",
            "coldata_path": str(FIXTURES / "coldata.csv"),
            "group_column": "group",
            "treat": "disease",
            "control": "control",
            "output_name": "test_deg",
        }
        args.update(overrides)
        return args

    def test_deseq2_finds_the_planted_signal(self) -> None:
        result = backend.deg(self.base_args(method="deseq2"))
        self.assertEqual(result["comparison"], "disease vs control")
        self.assertGreater(result["significant"]["up"], 0)
        self.assertGreater(result["significant"]["down"], 0)
        self.assertTrue(Path(result["result_table"]).exists())
        genes = {row["gene"] for row in result["top_genes"]}
        self.assertTrue(genes & {"IL6", "TNF", "IL1B"}, "planted up-regulated genes should surface")

    def test_covariates_enter_the_design(self) -> None:
        result = backend.deg(self.base_args(method="deseq2", covariates=["batch"]))
        self.assertIn("batch", result["design"])

    def test_log_transformed_input_is_rejected(self) -> None:
        # The guard that stops the most common fatal bioinformatics mistake.
        with self.assertRaises(RuntimeError) as ctx:
            backend.deg(self.base_args(method="deseq2", matrix_path=str(FIXTURES / "logged.csv")))
        self.assertIn("count", str(ctx.exception).lower())

    def test_unknown_group_level_is_reported(self) -> None:
        with self.assertRaises(RuntimeError) as ctx:
            backend.deg(self.base_args(treat="nonexistent"))
        self.assertIn("nonexistent", str(ctx.exception))

    def test_formula_like_group_column_is_not_executed(self) -> None:
        marker = Path(tempfile.gettempdir()) / "omics-deg-formula-injection-marker"
        marker.unlink(missing_ok=True)
        malicious = f"system('touch {marker}')"
        with (FIXTURES / "coldata.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            row[malicious] = row.pop("group")
        try:
            backend.deg({
                "method": "deseq2",
                "matrix_type": "counts",
                "matrix_path": str(FIXTURES / "counts.csv"),
                "coldata": rows,
                "group_column": malicious,
                "treat": "disease",
                "control": "control",
                "output_name": "test_safe_deg_formula",
            })
        finally:
            self.assertFalse(marker.exists())


@unittest.skipUnless(HAS_PLOT, "R with ggplot2/pheatmap not available")
class PlotTests(unittest.TestCase):
    def test_pca_reports_variance_explained(self) -> None:
        result = backend.plot({
            "type": "pca",
            "matrix_path": str(FIXTURES / "counts.csv"),
            "matrix_type": "counts",
            "coldata_path": str(FIXTURES / "coldata.csv"),
            "colour_column": "group",
            "output_name": "test_pca",
        })
        self.assertTrue(Path(result["figure"]["png"]).exists())
        self.assertTrue(Path(result["figure"]["svg"]).exists())
        self.assertGreater(result["variance_explained"]["PC1"], 0)
        self.assertEqual(result["transformation"], "log2 CPM (TMM, prior.count=2)")

    def test_non_numeric_matrix_cell_is_rejected(self) -> None:
        with self.assertRaises(RuntimeError) as ctx:
            backend.plot({
                "type": "pca",
                "matrix_type": "normalized",
                "matrix": [
                    {"gene": "G1", "S1": 1, "S2": "n/a", "S3": 3},
                    {"gene": "G2", "S1": 2, "S2": 3, "S3": 4},
                ],
            })
        self.assertIn("non-numeric", str(ctx.exception))

    def test_heatmap_rejects_partial_sample_annotation(self) -> None:
        with self.assertRaises(RuntimeError) as ctx:
            backend.plot({
                "type": "heatmap",
                "matrix_type": "normalized",
                "matrix": [
                    {"gene": "G1", "S1": 1, "S2": 2, "S3": 3},
                    {"gene": "G2", "S1": 3, "S2": 2, "S3": 1},
                ],
                "coldata": [
                    {"sample": "S1", "group": "A"},
                    {"sample": "S2", "group": "A"},
                ],
                "annotation_columns": ["group"],
            })
        self.assertIn("S3", str(ctx.exception))

    def test_venn_returns_region_membership(self) -> None:
        result = backend.plot({
            "type": "venn",
            "sets": {"A": ["TP53", "EGFR", "MYC"], "B": ["MYC", "PTEN"]},
            "output_name": "test_venn",
        })
        self.assertEqual(result["intersection_all"], ["MYC"])
        self.assertTrue(Path(result["membership_table"]).exists())

    def test_venn_rejects_too_many_sets(self) -> None:
        with self.assertRaises(RuntimeError):
            backend.plot({
                "type": "venn",
                "sets": {name: ["G1"] for name in "ABCDE"},
            })


@unittest.skipUnless(HAS_SURVIVAL, "R with glmnet/survival not available")
class SurvivalTests(unittest.TestCase):
    def test_lasso_cox_selects_variables_and_warns(self) -> None:
        (backend.ROOT / "Rplots.pdf").unlink(missing_ok=True)
        result = backend.survival({
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
        self.assertFalse((backend.ROOT / "Rplots.pdf").exists())

    def test_formula_like_predictor_name_is_not_executed(self) -> None:
        marker = Path(tempfile.gettempdir()) / "omics-formula-injection-marker"
        marker.unlink(missing_ok=True)
        with (FIXTURES / "survival.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        malicious = f"system('touch {marker}')"
        for row in rows:
            row[malicious] = row["gene1"]
        try:
            try:
                backend.survival({
                    "method": "lasso_cox",
                    "data": rows,
                    "time": "time",
                    "event": "event",
                    "predictors": [malicious, "gene2"],
                    "lambda": "min",
                    "output_name": "test_safe_formula",
                })
            except RuntimeError:
                pass
        finally:
            self.assertFalse(marker.exists())

    def test_too_few_events_is_rejected(self) -> None:
        rows = [{"time": 10 + i, "event": 1 if i < 2 else 0, "g1": i, "g2": -i} for i in range(20)]
        with self.assertRaises(RuntimeError) as ctx:
            backend.survival({
                "method": "lasso_cox", "data": rows,
                "time": "time", "event": "event", "predictors": ["g1", "g2"],
            })
        self.assertIn("event", str(ctx.exception).lower())


@unittest.skipUnless(HAS_TIMEROC, "R with timeROC/survival not available")
class TimeRocTests(unittest.TestCase):
    def test_timeroc_loads_survival_dependency(self) -> None:
        with (FIXTURES / "survival.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            row["risk_score"] = row["gene1"]
        result = backend.survival({
            "method": "timeroc",
            "data": rows,
            "time": "time",
            "event": "event",
            "risk_column": "risk_score",
            "times": [5, 10],
            "output_name": "test_timeroc",
        })
        self.assertEqual(result["method"], "timeroc")
        self.assertEqual(len(result["auc"]), 2)


@unittest.skipUnless(HAS_RMS, "R with rms/survival not available")
class CalibrationTests(unittest.TestCase):
    def test_groups_means_requested_number_of_groups(self) -> None:
        result = backend.survival({
            "method": "calibration",
            "data_path": str(FIXTURES / "survival.csv"),
            "time": "time",
            "event": "event",
            "predictors": ["gene1", "gene2"],
            "groups": 5,
            "bootstrap": 2,
            "times": [5],
            "output_name": "test_calibration",
        })
        self.assertEqual(result["groups"], 5)
        self.assertEqual(result["per_group"], 40)


@unittest.skipUnless(HAS_ENRICH, "R enrichment packages not available")
class EnrichmentTests(unittest.TestCase):
    def test_ensembl_identifiers_are_supported(self) -> None:
        result = backend.enrich({
            "method": "go",
            "species": "human",
            "id_type": "ENSEMBL",
            "genes": [
                "ENSG00000141510",
                "ENSG00000146648",
                "ENSG00000136997",
                "ENSG00000171862",
                "ENSG00000157764",
            ],
            "output_name": "test_ensembl",
        })
        self.assertEqual(result["id_type"], "ENSEMBL")
        self.assertGreater(result["genes_mapped"], 0)

    def test_gsva_uses_gaussian_for_normalized_matrix(self) -> None:
        matrix = []
        for index, gene in enumerate(["G1", "G2", "G3", "G4", "G5", "G6"], start=1):
            matrix.append({"gene": gene, "S1": 50.5 + index, "S2": 60.5 + index, "S3": 70.5 + index})
        result = backend.enrich({
            "method": "gsva",
            "species": "human",
            "id_type": "SYMBOL",
            "matrix_type": "normalized",
            "matrix": matrix,
            "gene_sets": {"set_a": ["G1", "G2", "G3", "G4", "G5", "G6"]},
            "output_name": "test_gsva",
        })
        self.assertEqual(result["kcdf"], "Gaussian")


if __name__ == "__main__":
    unittest.main()
