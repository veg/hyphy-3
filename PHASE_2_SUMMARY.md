# Phase 2: Gradient Computation for κ and ω Parameters

## Overview
Phase 2 extends the gradient computation framework to the transition/transversion rate ratio (κ) and non-synonymous/synonymous ratio (ω) parameters, enabling joint parameter optimization.

## Components Implemented

### Phase 2A: Omega (ω) Gradient Computation
**Commit:** `3ec639d` - Phase 2: Implement omega gradient computation

**Key Files Modified:**
- `include/hyphy/core/likelihood.hpp`
- `src/core/likelihood.cpp`

**Functionality:**
- Gradient of log-likelihood with respect to ω parameter
- Handles per-site and global ω computation
- Supports both fixed and variable ω across sites

**Technical Details:**
- Uses chain rule with previously computed β gradient infrastructure
- Efficient computation through matrix derivative reuse
- Numerical stability through vectorization

### Phase 2B: Kappa (κ) Gradient Computation  
**Commit:** `895b456` - Phase 2: Implement kappa gradient computation

**Key Files Modified:**
- `include/hyphy/core/likelihood.hpp`
- `src/core/likelihood.cpp`

**Functionality:**
- Gradient of log-likelihood with respect to κ parameter
- Transition/transversion rate ratio differentiation
- Supports DNA mutation model gradients

**Technical Details:**
- Computes ∂Q/∂κ for substitution rate matrix
- Handles constraint that κ > 0
- Integrated into unified gradient computation framework

## Mathematical Foundation

Gradients enable optimization of:
```
log L = Σ_sites log L_site(β, ω, κ, ...)

∂ log L / ∂θ = Σ_sites ∂ log L_site / ∂θ
```

Where θ ∈ {β, ω, κ, branch_lengths, ...}

## Dependencies
- Phase 1 infrastructure (Fréchet derivatives, matrix exponentials)
- Eigen library
- C++17 standard library

## Testing Status
✅ **COMPLETE**
- Gradient verification tests implemented
- Numerical gradient checking with finite differences
- Integration tests with optimization routines

## Integration Notes
- Works seamlessly with Phase 1 Fréchet derivative infrastructure
- Provides foundation for Phase 3 L-BFGS optimization
- Enables both global and per-site parameter fitting

## Notes for Future Development
- Consider implementing automatic differentiation (AD) for gradient verification
- May benefit from algorithmic differentiation for new parameters
- Gradient computation can be parallelized further for multiple sites
