#!/usr/bin/env bash
# ==============================================================================
# HyPhy 3: Branch-site Unrestricted Statistical Test (BUSTED) Example
# Tests for gene-wide episodic diversifying selection.
# ==============================================================================
set -euo pipefail

HYPHY3_BIN="${HYPHY3_BIN:-../../build/hyphy3}"

if [[ ! -x "${HYPHY3_BIN}" ]]; then
    echo "Error: hyphy3 binary not found at ${HYPHY3_BIN}. Please build the project first."
    exit 1
fi

DATA_DIR="../../benchmarks/data"

echo "=== 1. Standard BUSTED (K=3 rate categories) on CD2 ==="
"${HYPHY3_BIN}" busted \
    --alignment "${DATA_DIR}/cd2.fna" \
    --tree "${DATA_DIR}/cd2.nwk" \
    --rates 3 \
    --threads 4 \
    --output "cd2.BUSTED.json"

echo -e "\n=== 2. BUSTED with Automatic Model Selection (Auto-K via AICc) on ADH ==="
"${HYPHY3_BIN}" busted \
    --alignment "${DATA_DIR}/adh.fna" \
    --tree "${DATA_DIR}/adh.nwk" \
    --auto-k \
    --threads 4 \
    --output "adh.BUSTED_autoK.json"

echo -e "\n=== 3. BUSTED-S: Synonymous Rate Variation (SRV) on ADH ==="
"${HYPHY3_BIN}" busted \
    --alignment "${DATA_DIR}/adh.fna" \
    --tree "${DATA_DIR}/adh.nwk" \
    --srv \
    --syn-rates 3 \
    --rates 3 \
    --threads 4 \
    --output "adh.BUSTED_S.json"

echo -e "\nBUSTED analyses completed successfully!"
