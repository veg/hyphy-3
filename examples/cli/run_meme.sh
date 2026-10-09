#!/usr/bin/env bash
# ==============================================================================
# HyPhy 3: Mixed Effects Model of Evolution (MEME) Example
# Tests for site-by-site episodic diversifying positive selection.
# ==============================================================================
set -euo pipefail

HYPHY3_BIN="${HYPHY3_BIN:-../../build/hyphy3}"

if [[ ! -x "${HYPHY3_BIN}" ]]; then
    echo "Error: hyphy3 binary not found at ${HYPHY3_BIN}. Please build the project first."
    exit 1
fi

DATA_DIR="../../benchmarks/data"

echo "=== Running MEME on CD2 Dataset ==="
"${HYPHY3_BIN}" meme \
    --alignment "${DATA_DIR}/cd2.fna" \
    --tree "${DATA_DIR}/cd2.nwk" \
    --pvalue 0.1 \
    --threads 4 \
    --output "cd2.MEME.json"

echo -e "\n=== Running MEME on ADH Dataset ==="
"${HYPHY3_BIN}" meme \
    --alignment "${DATA_DIR}/adh.fna" \
    --tree "${DATA_DIR}/adh.nwk" \
    --pvalue 0.1 \
    --threads 4 \
    --output "adh.MEME.json"

echo -e "\nMEME analysis completed successfully! Results saved to cd2.MEME.json and adh.MEME.json."
