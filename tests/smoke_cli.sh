#!/usr/bin/env bash
set -euo pipefail
OMICS_TEST_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OMICS_TEST_PYTHON="${OMICS_PYTHON:-python3}"
export PYTHONPATH="$OMICS_TEST_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
"$OMICS_TEST_PYTHON" -m monadomics --version
"$OMICS_TEST_PYTHON" -m unittest discover -s "$OMICS_TEST_ROOT/tests" -p test_cli.py
