# HyPhy 3: Comprehensive Developer & AI Agent Guide

> **Audience**: AI Coding Agents (Antigravity, Claude, Copilot) and Human System Architects.  
> **Mission**: Implement, extend, test, and optimize phylogenetic models in HyPhy 3 with strict mathematical rigor and zero regressions.

---

## 1. Core Principles & Mandates

1. **Zero Synthetic Mockery Mandate**:
   - **NEVER** script synthetic mathematical placeholders (e.g. parametric curves, fabricated p-values, artificial effect sizes).
   - Every single number, test statistic, rate estimate, and p-value MUST be the direct output of executed computational code operating on authentic biological input data.
2. **Strict Parity vs HyPhy 2.5 Baseline**:
   - Every model implemented in HyPhy 3 must achieve numerical parity with HyPhy 2.5 outputs down to $|\Delta \ln L| < 10^{-2}$ and consistent parameter estimates on standard benchmarks (ADH, CD2, COXI, β-globin).
3. **C++20 Modular Architecture**:
   - Header-only core algorithms for inlining and aggressive compiler vectorization.
   - No macros, no global state, no HBL interpreter overhead, full thread-safety under OpenMP.
4. **Differentiable Phylogenetics by Design**:
   - When optimizing continuous parameters (such as branch lengths), prefer **exact analytical Inside-Outside adjoint gradients** over finite differences or coordinate-wise searches.

---

## 2. Codebase Anatomy & Architecture

```
include/hyphy/
├── core/
│   ├── types.hpp           # Scalar (double), Matrix, Vector (Eigen3 column-major)
│   ├── genetic_code.hpp    # Standard & non-standard codes, sense codon maps (61 sense)
│   ├── alignment.hpp       # FASTA & NEXUS parsing, F3x4 codon frequencies, pattern compression
│   ├── tree.hpp            # Newick parsing, node hierarchies, post-order caching
│   ├── rate_matrix.hpp     # Analytical GTR (4x4) & MG94 (61x61) matrix exponentials & Q
│   └── likelihood.hpp      # Vectorized Felsenstein pruning & Inside-Outside adjoint engine
├── opt/
│   ├── brent.hpp           # 1D bounded golden-section / parabolic interpolation
│   ├── nelder_mead.hpp     # Bounded simplex optimization (2D, 3D, and ND)
│   ├── squarem.hpp         # Squared polynomial extrapolation for EM acceleration
│   └── lbfgs.hpp           # Parallel multi-branch L-BFGS with Armijo line search
├── autograd/
│   ├── var.hpp             # Dynamic computational graph node with reverse backprop
│   ├── tree_likelihood.hpp # Differentiable tree likelihood node with analytical adjoints
│   └── optimizer.hpp       # Native C++ Adam optimizer
└── analyses/
    ├── gtr.hpp             # Baseline nucleotide GTR model fitter
    ├── mg94.hpp            # Global MG94xREV baseline model fitter
    ├── fel.hpp             # Fixed Effects Likelihood (pervasive selection)
    ├── meme.hpp            # Mixed Effects Model of Evolution (episodic selection)
    └── busted.hpp          # BUSTED & BUSTED-S (gene-wide selection with SRV)
```

---

## 3. Mathematical Foundations

### 3.1 The Inside-Outside Adjoint Algorithm
For a phylogenetic tree with nodes $V$ and branches $E$, site pattern likelihood is computed via Felsenstein's post-order pruning (**Inside Pass**):
$$D_{u}(i) = \prod_{v \in \text{children}(u)} \sum_{j} P_{uv}(i, j) D_{v}(j)$$
At the root: $L_p = \sum_{i} \pi_i D_{\text{root}}(i)$.

The **Outside Pass** propagates likelihoods backwards from root to leaves:
$$V_{v} = P_{uv}^{\top} \left( V_{u} \odot \prod_{w \in \text{siblings}(v)} P_{uw} D_{w} \right)$$
where $V_{\text{root}} = \boldsymbol{\pi}$.

The exact analytical gradient with respect to branch length $t_b$ ($b = u \to v$) is computed via the contraction:
$$\frac{\partial L_p}{\partial t_b} = V_{b,p}^{\top} \left( \frac{\partial P_b}{\partial t_b} \right) D_{b,p}$$
where $\frac{\partial P_b}{\partial t_b} = Q \exp(Q t_b) = Q P_b$.

Total log-likelihood gradient:
$$\frac{\partial \ln L}{\partial t_b} = \sum_{p=1}^P \frac{w_p}{L_p} V_{b,p}^{\top} \left( Q P_b \right) D_{b,p}$$
**Computational Cost**: Exactly $\mathcal{O}(B)$ total tree operations for all $B$ branches simultaneously, compared to $\mathcal{O}(B^2)$ for finite differences or Brent coordinate search!

### 3.2 Synonymous Rate Variation Invariance (BUSTED-S)
Because the continuous-time Markov generator satisfies:
$$Q(\alpha_m, \beta = \alpha_m \omega_k) = \alpha_m Q(1, \omega_k)$$
the transition matrix satisfies:
$$P_{m,k}(t_b) = \exp(Q_{m,k} t_b) = \exp(Q(1, \omega_k) \cdot (\alpha_m t_b)) = P_k(\alpha_m t_b)$$
**Crucial Efficiency Rule**: Never diagonalize $M \times K$ rate matrices. Diagonalize only the $K$ non-synonymous matrices $Q(1, \omega_k)$, and evaluate $P_k$ at the scaled branch lengths $\tau = \alpha_m t_b$.

### 3.3 Latent Mixture Acceleration with SQUAREM
For mixture weights $\mathbf{p} = (p_1, \dots, p_K)$ or synonymous weights $\mathbf{q} = (q_1, \dots, q_M)$, standard EM updates take step $\mathbf{r} = \mathbf{p}_1 - \mathbf{p}_0$ and $\mathbf{v} = (\mathbf{p}_2 - \mathbf{p}_1) - \mathbf{r}$.
SQUAREM extrapolates:
$$\mathbf{p}_{\text{accelerated}} = \mathbf{p}_0 - 2 \sigma \mathbf{r} + \sigma^2 \mathbf{v}, \quad \sigma = -\frac{\|\mathbf{r}\|}{\|\mathbf{v}\|}$$
This achieves near-Newtonian quadratic convergence without calculating second derivatives or Hessians.

---

## 4. Step-by-Step Recipe: Implementing a New Model

Follow this recipe to add any new selection model (e.g., SLAC, aBSREL, RELAX, FUBAR):

### Step 1: Formulate the Rate Matrix (`include/hyphy/core/rate_matrix.hpp`)
If your model introduces new parameters (e.g. relaxation parameter $k$ in RELAX, or asymmetric nucleotide parameters):
1. Define the parameter struct in `rate_matrix.hpp`.
2. Construct the infinitesimal generator $Q$ enforcing row sums $\sum_j Q_{ij} = 0$.
3. Compute the scaling factor $\rho = -\sum_i \pi_i Q_{ii}$ to normalize the substitution rate to 1 expected substitution per unit branch length.
4. Compute the analytical or spectral eigendecomposition: $Q = U \Lambda U^{-1}$ to allow fast transition matrix computation $P(t) = U \exp(\Lambda t) U^{-1}$.

### Step 2: Implement Model Analyzer (`include/hyphy/analyses/<model>.hpp`)
Create a self-contained header under `include/hyphy/analyses/`:
```cpp
#pragma once
#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include <nlohmann/json.hpp>

namespace hyphy::analyses {

struct MyModelSettings {
    bool refine_branches = true;
    Scalar pvalue_threshold = 0.05;
};

struct MyModelResult {
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
};

class MyModelAnalyzer {
public:
    Tree tree;
    Alignment aln;
    // ...
    static MyModelAnalyzer create_and_fit(Tree tree, Alignment aln);
    MyModelResult run(MyModelSettings settings = {}) const;
    nlohmann::json to_json(const MyModelResult& res) const;
};

} // namespace hyphy::analyses
```

### Step 3: Branch Optimization Strategy
- Always initialize branch lengths from the MG94 baseline fit.
- For joint branch length updates, use `LikelihoodEngine::compute_branch_length_gradients` and `L-BFGS` with Armijo backtracking:
```cpp
auto [curr_ll, curr_grad] = LikelihoodEngine::compute_branch_length_gradients(tree, aln, params, bl);
```
- For discrete latent mixtures, implement the E-step using `compute_inside_outside` and accelerate the M-step using `SQUAREM`.

### Step 4: Standard Datamonkey JSON Output
Ensure the `to_json()` method returns standard Datamonkey-compatible JSON with:
- `"analysis"`: metadata, version, settings.
- `"fits"`: Log-Likelihood, AIC-c, Rate Distributions.
- `"branch attributes"`: per-branch inferred parameters and lengths.
- `"test results"`: LRT, p-value.

### Step 5: Expose CLI and Python Bindings
1. **CLI**: Add command driver in `src/apps/` (e.g. `src/apps/hyphy_<model>.cpp`) and register subcommand in `src/apps/hyphy3.cpp`.
2. **Python**: Expose the analyzer and results structs in `src/python/bindings.cpp` via nanobind.

### Step 6: Parity Test Suite vs HyPhy 2.5
Add a unit test in `tests/test_<model>.cpp` verifying:
1. Convergence and numerical stability.
2. Parity vs HyPhy 2.5 JSON files stored in `benchmarks/data/` or `tests/data/`.
3. Check bounds: $|\Delta \ln L| < 0.1$, LRT $\ge 0$, $p \in [0, 1]$.

---

## 5. Verification Protocol & Quality Checklist

Before committing any model changes, you MUST verify:

- [ ] **Compilation**: Clean build under `-Wall -Wextra -Werror` without warnings.
- [ ] **All CTests Pass**: Run `ctest --output-on-failure`. All unit and parity tests must pass 100%.
- [ ] **Finite-Difference Gradient Check**: Any new analytical gradient must match numerical central finite differences:
  $$\frac{\left| g_{\text{analytical}} - g_{\text{finite\_diff}} \right|}{|g_{\text{analytical}}| + 10^{-8}} < 10^{-4}$$
- [ ] **Thread-Safety Check**: Ensure OpenMP private variables are allocated per thread (e.g., node likelihood vectors `node_L` and `InsideOutsideResult`).
- [ ] **Simplex Constraint Check**: Ensure mixture weights $\sum p_k = 1$ and $p_k \ge 0$ are preserved across all SQUAREM and Nelder-Mead updates.
