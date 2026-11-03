// Aggregated precompiled header.
//
// This lists the stable, heavily-included portion of the header stack. Every
// executable target reuses the compiled artifact of this file (see the
// HYPHY3_PCH block at the bottom of CMakeLists.txt), so these headers are
// parsed once per build instead of once per translation unit.
//
// Deliberately NOT included here: the per-analysis emitters
// (hyphy/analyses/*.hpp). They change frequently and are not shared by all
// targets, so precompiling them would add work rather than remove it.
#pragma once

#include <Eigen/Core>

#include "hyphy/core/types.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/alignment.hpp"
#include "nlohmann/json.hpp"
