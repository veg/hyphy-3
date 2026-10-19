# Phase 1: Gradient Computation Infrastructure for β Parameter

## Overview
Phase 1 implements the gradient computation infrastructure for the stringency parameter (β) in the ExpCM model. This phase establishes the foundation for parameter optimization.

## Components Implemented

### Phase 1A: Beta Gradient Computation Infrastructure
**Commits:**
- `56366da` - Phase 1A: Implement beta gradient computation infrastructure
- `ae8b83a` - Phase 1B: Add spectral Fréchet derivative infrastructure for beta gradient
- `7417482` - Make matrix_exp_scaling_squaring public for testing

**Key Files Modified:**
- `include/hyphy/core/likelihood.hpp`
- `src/core/likelihood.cpp`
- `include/hyphy/core/matrix_exp.hpp`
- `src/core/matrix_exp.cpp`

**Functionality:**
- Spectral decomposition-based matrix exponential computation
- Fréchet derivative computation for matrix exponentials
- Beta gradient computation for site-specific likelihood functions
- Numerical stability through matrix scaling-and-squaring method

**Technical Details:**
- Uses eigenvalue decomposition for efficient matrix exponential
- Implements Fréchet linear operator for differential computation
- Handles numerical precision through careful scaling

**Testing:**
- Matrix exponential verification tests
- Stationary distribution validation
- F matrix computation tests
- Per-site likelihood tests

## Mathematical Foundation

The ExpCM model defines rate matrix Q with:
```
Q_{xy} = β · μ_r · ω_r · R_xy · π_y
```

Where:
- β: stringency parameter (fitness effect scale)
- μ_r: site-specific synonymous rate
- ω_r: site-specific dN/dS ratio
- R_xy: reversible substitution rate
- π_y: stationary frequency

Phase 1 enables computing ∂Q/∂β for gradient-based optimization.

## Dependencies
- Eigen library for linear algebra
- GSL for special function approximations
- C++17 standard library

## Status
✅ **COMPLETE** - All components tested and integrated into optimization pipeline

## Notes for Future Development
- Fréchet derivative computation is foundational for all parameter gradients
- Matrix exponential stability is critical for numerical accuracy
- Spectral method chosen for efficiency in small-dimension problems (61×61 for codons)
