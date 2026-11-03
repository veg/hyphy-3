// Owner translation unit for the shared precompiled header (src/pch.hpp).
// Compiling this file produces the .gch artifact every executable target
// reuses via target_precompile_headers(... REUSE_FROM hyphy3_pch).
#include "pch.hpp"

namespace hyphy::detail {
int pch_anchor_symbol = 0;
} // namespace hyphy::detail

