#!/usr/bin/env bash
# ==============================================================================
# HyPhy 3: Fixed Effects Likelihood (FEL) Example
# Tests for site-by-site pervasive diversifying and purifying selection.
# ==============================================================================
set -euo pipefail

# Ensure hyphy3 executable is built
HYPHY3_BIN="${HYPHY3_BIN:-../../build/hyphy3}"

if [[ ! -x "${HYPHY3_BIN}" ]]; then
    echo "Error: hyphy3 binary not found at ${HYPHY3_BIN}. Please build the project first."
    exit 1
fi

DATA_DIR="../../benchmarks/data"

echo "=== Running FEL on CD2 Dataset (10 taxa, 187 codons) ==="
"${HYPHY3_BIN}" fel \
    --alignment "${DATA_DIR}/cd2.fna" \
    --tree "${DATA_DIR}/cd2.nwk" \
    --pvalue 0.1 \
    --threads 4 \
    --output "cd2.FEL.json"

echo -e "\n=== Running FEL on ADH Dataset (23 taxa, 254 codons) ==="
"${HYPHY3_BIN}" fel \
    --alignment "${DATA_DIR}/adh.fna" \
    --tree "${DATA_DIR}/adh.nwk" \
    --pvalue 0.05 \
    --threads 4 \
    --output "adh.FEL.json"

echo -e "\nFEL analysis completed successfully! Results saved to cd2.FEL.json and adh.FEL.json."
