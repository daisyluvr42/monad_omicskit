import contextlib
import importlib.util
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from monadomics import backend, cli


class CliTests(unittest.TestCase):
    def invoke(self, arguments, input_text=""):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), mock.patch.object(sys, "stdin", io.StringIO(input_text)):
            code = cli.main(arguments)
        return code, json.loads(output.getvalue())

    def test_cli_subprocess_outputs_valid_json(self):
        result = subprocess.run([sys.executable, "-m", "monadomics", "capabilities"], env={**os.environ, "PYTHONPATH": str(ROOT / "src")}, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        payload = json.loads(result.stdout)
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["count"], 21)

    def test_schema_exposes_scientific_constraints(self):
        code, result = self.invoke(["schema", "enrich"])
        self.assertEqual(code, 0)
        self.assertIn("species", result["inputSchema"]["required"])
        self.assertIn("id_type", result["inputSchema"]["required"])

    def test_invalid_json_is_structured_error(self):
        code, result = self.invoke(["plot", "--params", "-"], "{")
        self.assertEqual(code, 1)
        self.assertFalse(result["ok"])

    def test_non_object_params_are_rejected(self):
        code, result = self.invoke(["plot", "--params", "-"], "[]")
        self.assertEqual(code, 1)
        self.assertIn("JSON object", result["message"])

    def test_missing_required_parameter_does_not_launch_r(self):
        with mock.patch.object(backend, "_run_r") as run:
            code, result = self.invoke(["deg", "--params", "-"], "{}")
        self.assertEqual(code, 1)
        self.assertIn("matrix_type", result["message"])
        run.assert_not_called()

    def test_unknown_parameter_is_rejected(self):
        code, result = self.invoke(["plot", "--params", "-"], '{"type":"pca","matrix_typ":"counts"}')
        self.assertEqual(code, 1)
        self.assertIn("matrix_typ", result["message"])

    def test_invalid_method_is_rejected(self):
        code, result = self.invoke(["plot", "--params", "-"], '{"type":"invented"}')
        self.assertEqual(code, 1)
        self.assertIn("type must be one of", result["message"])

    def test_relative_and_unicode_paths_resolve_before_r(self):
        with tempfile.TemporaryDirectory(prefix="omics 空格 ") as directory:
            params_path = Path(directory) / "参数.json"
            params_path.write_text(json.dumps({"type": "pca", "matrix_type": "counts", "matrix_path": "数据/counts.csv"}), encoding="utf-8-sig")
            output_dir = Path(directory) / "结果"
            with mock.patch.dict(backend.HANDLERS, {"plot": mock.Mock(return_value={"figure": {}})}), mock.patch.object(backend, "OUTPUT_DIR", backend.OUTPUT_DIR):
                handler = backend.HANDLERS["plot"]
                code, _ = self.invoke(["plot", "--params", str(params_path), "--output-dir", str(output_dir), "--timeout", "1200"])
                self.assertEqual(code, 0)
                self.assertEqual(handler.call_args.args[0]["matrix_path"], str(Path("数据/counts.csv").resolve()))
                self.assertEqual(handler.call_args.args[0]["timeout"], 1200)
                self.assertEqual(backend.OUTPUT_DIR, output_dir.resolve())

    def test_doctor_missing_r_returns_nonzero_without_install(self):
        with mock.patch.dict(os.environ, {"OMICS_RSCRIPT": "/nonexistent/monadomics/Rscript"}):
            code, result = self.invoke(["doctor", "--group", "deg"])
        self.assertEqual(code, 1)
        self.assertFalse(result["r_available"])

    def test_analysis_failure_is_not_a_success_response(self):
        with mock.patch.dict(backend.HANDLERS, {"plot": mock.Mock(side_effect=RuntimeError("R rejected matrix"))}):
            code, result = self.invoke(["plot", "--params", "-"], '{"type":"pca"}')
        self.assertEqual(code, 1)
        self.assertFalse(result["ok"])
        self.assertEqual(result["message"], "R rejected matrix")

    def test_nonzero_r_exit_is_not_hidden_by_result_file(self):
        def fail_with_result(command, **kwargs):
            Path(command[-1]).write_text('{"result_table":"partial.csv"}')
            return subprocess.CompletedProcess(command, 1, "", "R failed after writing output")
        with mock.patch.object(backend, "_find_rscript", return_value="Rscript"), mock.patch.object(backend.subprocess, "run", side_effect=fail_with_result):
            with self.assertRaisesRegex(RuntimeError, "exit 1"):
                backend.plot({"type": "pca"})

    def test_r_timeout_is_reported(self):
        with mock.patch.object(backend, "_find_rscript", return_value="Rscript"), mock.patch.object(backend.subprocess, "run", side_effect=subprocess.TimeoutExpired("Rscript", 1)):
            with self.assertRaisesRegex(RuntimeError, "exceeded 1s"):
                backend.plot({"type": "pca", "timeout": 1})


class PackagingTests(unittest.TestCase):
    def test_connector_zip_has_required_root_and_resolved_references(self):
        spec = importlib.util.spec_from_file_location("connector_build", ROOT / "scripts/build_workbuddy_connector.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with zipfile.ZipFile(module.build()) as archive:
            required = {"connector-meta.json", "cli.json", "icon.svg", "skills/monadomics-analysis/SKILL.md"}
            self.assertTrue(required.issubset(archive.namelist()))
            self.assertNotIn("mcp.json", archive.namelist())
            self.assertIsNone(archive.testzip())


if __name__ == "__main__":
    unittest.main()
