#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${REPO_ROOT}/build-wasm"
WEB_DIR="${REPO_ROOT}/web"

# Auto-detect Homebrew Emscripten toolchain paths if not already set in the environment
if [[ -z "${EMSDK_PYTHON:-}" ]] && [[ -x "/opt/homebrew/bin/python3" ]]; then
    export EMSDK_PYTHON="/opt/homebrew/bin/python3"
fi

if [[ -z "${EM_LLVM_ROOT:-}" ]]; then
    BREW_EM_LLVM="$(find /opt/homebrew/Cellar/emscripten -maxdepth 4 -type d -name bin -path "*/libexec/llvm/bin" 2>/dev/null | sort -V | tail -n 1 || true)"
    if [[ -n "${BREW_EM_LLVM}" ]]; then
        export EM_LLVM_ROOT="${BREW_EM_LLVM}"
    fi
fi

if [[ -z "${EM_BINARYEN_ROOT:-}" ]]; then
    BREW_EM_BIN="$(find /opt/homebrew/Cellar/emscripten -maxdepth 3 -type d -name binaryen -path "*/libexec/binaryen" 2>/dev/null | sort -V | tail -n 1 || true)"
    if [[ -n "${BREW_EM_BIN}" ]]; then
        export EM_BINARYEN_ROOT="${BREW_EM_BIN}"
    fi
fi

echo "==> Configuring HyPhy-3 WASM target in ${BUILD_DIR}..."
emcmake cmake -S "${REPO_ROOT}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release

echo "==> Building hyphy3_wasm (hyphy3_walrus.js + hyphy3_walrus.wasm)..."
cmake --build "${BUILD_DIR}" --target hyphy3_wasm -j8

mkdir -p "${WEB_DIR}"
cp -f "${BUILD_DIR}/hyphy3_walrus.js" "${WEB_DIR}/hyphy3_walrus.js"
cp -f "${BUILD_DIR}/hyphy3_walrus.wasm" "${WEB_DIR}/hyphy3_walrus.wasm"

echo "==> WASM artifacts copied to ${WEB_DIR}:"
ls -lh "${WEB_DIR}/hyphy3_walrus.js" "${WEB_DIR}/hyphy3_walrus.wasm"
