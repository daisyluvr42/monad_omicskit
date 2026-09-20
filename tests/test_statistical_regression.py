import subprocess
import unittest

from test_backend import ROOT, backend, r_packages_present


@unittest.skipUnless(r_packages_present("DESeq2", "edgeR", "ggplot2", "pheatmap", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db"), "R statistical dependencies unavailable")
class StatisticalRegressionTests(unittest.TestCase):
    def test_lightweight_checkpoint_contracts(self):
        result = subprocess.run(
            [backend._find_rscript(), "--vanilla", str(ROOT / "tests/checkpoint_regression.R"), str(ROOT)],
            capture_output=True, text=True, timeout=600, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_general_statistical_contracts(self):
        result = subprocess.run(
            [backend._find_rscript(), "--vanilla", str(ROOT / "tests/statistical_regression.R"), str(ROOT)],
            capture_output=True, text=True, timeout=600, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
