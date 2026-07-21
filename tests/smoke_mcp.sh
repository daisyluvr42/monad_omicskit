#!/usr/bin/env bash
# Verify the Omics MCP server initializes, lists tools, and answers a call.
# Runs without R: omics_env reports R availability instead of failing.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON="${OMICS_PYTHON:-python3}"

request() {
  printf '%s\n' "$1"
}

output="$(
  {
    request '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{}}}'
    request '{"jsonrpc":"2.0","method":"notifications/initialized"}'
    request '{"jsonrpc":"2.0","id":2,"method":"tools/list"}'
    request '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"omics_feature_menu","arguments":{}}}'
    request '{"jsonrpc":"2.0","id":4,"method":"tools/call","params":{"name":"omics_env","arguments":{"group":"deg"}}}'
    request '{"jsonrpc":"2.0","id":5,"method":"tools/call","params":{"name":"omics_deg","arguments":{}}}'
  } | "$PYTHON" "$ROOT/mcp/omics_mcp.py"
)"

check() {
  if ! grep -q "$1" <<<"$output"; then
    echo "FAILED: expected $2" >&2
    echo "--- server output ---" >&2
    echo "$output" >&2
    exit 1
  fi
}

check '"name": "omics"' "server identity in initialize"
check 'omics_deg' "tool list to include omics_deg"
check 'omics_survival' "tool list to include omics_survival"
check 'Omics 分析工作台' "feature menu title"
check 'r_available' "omics_env to report R availability"
check '"isError": true' "a bad omics_deg call to return a tool error"

echo "Omics MCP smoke test passed"
