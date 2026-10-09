#!/usr/bin/env bash
# ==============================================================================
# HyPhy 3: Adaptive Branch-Site Random Effects Likelihood (aBSREL) Example
# Tests for episodic diversifying selection on specific lineages and branches.
# ==============================================================================
set -euo pipefail

HYPHY3_BIN="${HYPHY3_BIN:-../../build/hyphy3}"

if [[ ! -x "${HYPHY3_BIN}" ]]; then
    echo "Error: hyphy3 binary not found at ${HYPHY3_BIN}. Please build the project first."
    exit 1
fi

DATA_DIR="../../benchmarks/data"

echo "=== 1. Standard aBSREL on Beta-Globin (Testing All Branches) ==="
"${HYPHY3_BIN}" absrel \
    --alignment "${DATA_DIR}/bglobin.nex" \
    --max-rates 3 \
    --pvalue 0.05 \
    --threads 4 \
    --output "bglobin.ABSREL.json"

echo -e "\n=== 2. aBSREL on CD2 (Testing Internal Branches Only) ==="
"${HYPHY3_BIN}" absrel \
    --alignment "${DATA_DIR}/cd2.fna" \
    --tree "${DATA_DIR}/cd2.nwk" \
    --branches Internal \
    --max-rates 3 \
    --pvalue 0.05 \
    --threads 4 \
    --output "cd2_internal.ABSREL.json"

echo -e "\n=== aBSREL CLI Example Completed Successfully ==="
