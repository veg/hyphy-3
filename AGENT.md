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
│   ├── types.hpp                 # Scalar (double), Matrix, Vector (Eigen3 column-major)
│   ├── genetic_code.hpp          # Standard & non-standard codes, sense codon maps (61 sense)
│   ├── alignment.hpp             # FASTA & NEXUS parsing, F3x4 codon frequencies, pattern compression
│   ├── tree.hpp                  # Newick parsing, WheatstoneBridge hubs, eval_internal_node
│   ├── neighbor_joining.hpp      # Saitou-Nei NJ from TN93 / Jukes-Cantor distances
│   ├── ohmic_network.hpp         # Foster effective resistance, ScalableWoodbury, HyphAeon 384D bridge
│   ├── wheatstone_likelihood.hpp # Linear-time 4-way Wheatstone pruning & Inside-Outside Fréchet adjoints
│   ├── rate_matrix.hpp           # Analytical GTR (4x4) & MG94 (61x61) matrix exponentials & Q
│   └── likelihood.hpp            # Vectorized Felsenstein pruning & Inside-Outside adjoint engine
├── opt/
│   ├── brent.hpp                 # 1D bounded golden-section / parabolic interpolation
│   ├── nelder_mead.hpp           # Bounded simplex optimization (2D, 3D, and ND)
│   ├── squarem.hpp               # Squared polynomial extrapolation for EM acceleration
│   └── lbfgs.hpp                 # Parallel multi-branch L-BFGS with Armijo line search
├── autograd/
│   ├── var.hpp                   # Dynamic computational graph node with reverse backprop
│   ├── tree_likelihood.hpp       # Differentiable tree likelihood node with analytical adjoints
│   └── optimizer.hpp             # Native C++ Adam optimizer
└── analyses/
    ├── gtr.hpp                   # Baseline nucleotide GTR model fitter
    ├── mg94.hpp                  # Global MG94xREV baseline model fitter
    ├── fel.hpp                   # Fixed Effects Likelihood (pervasive selection)
    ├── contrast_fel.hpp          # Contrast-FEL (multi-group selection with Empirical Bayes Ridge)
    ├── directional_fel.hpp       # Directional-FEL (non-stationary Softmax regimes & ML auto-rooting)
    ├── fubar.hpp                 # FUBAR (bivariate rate grid with VB0 / Collapsed Gibbs / MH)
    ├── prime.hpp                 # PRIME (biophysical property selection with L2 Ridge shrinkage)
    ├── meme.hpp                  # Mixed Effects Model of Evolution (episodic selection)
    ├── busted.hpp                # BUSTED, BUSTED-S, BUSTED-MH, BUSTED-E & 4-Gate BUSTED-Auto
    ├── absrel.hpp                # aBSREL (adaptive branch-site random effects likelihood)
    ├── relax.hpp                 # RELAX (selection relaxation / intensification test)
    ├── slac.hpp                  # SLAC (Joint ML ancestral state reconstruction & binomial counting)
    ├── gard.hpp                  # GARD (Genetic Algorithm for Recombination Detection)
    └── ohmic.hpp                 # OHMIC (Ohmic resistor network & Wheatstone bridge ML refinement)
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

### 3.4 Local Likelihood Projection for Branch-Site Models (aBSREL)
When optimizing parameters on an individual branch $b = (u \to v)$ while keeping the rest of the tree fixed (e.g. during aBSREL branch complexity step-up or constrained null hypothesis testing):
The conditional subtree likelihoods $D_{v,p}$ and ancestral complement likelihoods $V_{v,p}$ are **strictly invariant** with respect to branch $b$'s parameters.
The pattern likelihood simplifies to an immediate dot-product projection:
$$L_p = V_{v,p}^{\top} \bar{P}_b D_{v,p}$$
where $\bar{P}_b = \sum_{k=1}^K p_{b,k} P(t_b, \omega_{b,k})$.
Evaluating candidate mixtures or optimizing rates on branch $b$ requires **zero tree traversals**, reducing evaluation time to $\approx 10\,\mu\text{s}$.
Furthermore, the asymptotic null distribution for aBSREL hypothesis testing ($\omega_{b,\max} \le 1$) is a mixture:
$$\text{Null} \sim \frac{1}{2} \chi^2_0 + \frac{1}{2}\left(0.4 \chi^2_1 + 0.6 \chi^2_2\right)$$
with closed-form $p$-value:
$$p = \frac{1}{2}\left[ 0.4 \cdot \text{erfc}\left(\sqrt{\frac{\text{LRT}}{2}}\right) + 0.6 \cdot \exp\left(-\frac{\text{LRT}}{2}\right) \right]$$

### 3.5 Ohmic Networks, Linear-Time Wheatstone Pruning & HyphAeon 384D
To prevent homoplasy and incomplete lineage sorting (ILS) on short internal branches from inflating terminal branch lengths or generating false-positive selection calls:
1. **Resistor Network & Foster's Theorem**: Taxa form an electrical conductance network $\mathbf{L} = \mathbf{D} - \mathbf{C}$ (augmented by the **HyphAeon 384D Multi-Scale Transformer** embeddings $\mathbf{E} \in \mathbb{R}^{N \times 384}$ split into $192\text{D}$ $dS$ and $192\text{D}$ $dN$ channels). Effective resistance $\Omega_{ij} = L^+_{ii} + L^+_{jj} - 2L^+_{ij}$ satisfies $\sum_{i<j} c_{ij}\Omega_{ij} = N - 1$ (`ExactDense`) or is evaluated in $\mathcal{O}(N \log N)$ time via equipotential busbars and the rank-$p$ **Woodbury Matrix Identity** (`ScalableWoodbury`).
2. **Series-Parallel Wheatstone Bridges**: High-stress internal splits ($\tau_e^{\text{Ohm}} \ge \tau_{\min}$) are promoted into 3-child `WheatstoneBridge` hubs embedded directly on `hyphy::core::Tree`. During Felsenstein pruning (`Tree::eval_internal_node`), each hub evaluates the convex combination over the 3 quartet resolutions:
   $$\mathbf{D}_k^{\text{Ohm}} = \sum_{r=1}^3 w_{k,r} \left[ \left(\mathbf{P}(t_{k,r})\,(\mathbf{M}_{a_r} \odot \mathbf{M}_{b_r})\right) \odot \mathbf{M}_{c_r} \right]$$
   preserving $\mathcal{O}(N)$ complexity and enabling every standard HyPhy 3 analysis to run with `--ohmic` / `--omhic` or on tree-free FASTA alignments via `OhmicNetwork::resolve_tree`.

---

## 4. Step-by-Step Recipe: Implementing a New Model

Follow this recipe to add any new selection model (e.g., RELAX, FUBAR, SLAC, BUSTED-PH):

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

---

## 6. Shared Infrastructure: CI, Tests & Commits

Other people maintain parts of this repository (CI, test harness, report viewers). Changes to them land on `main` between your sessions.

1. **Start from the current `main`, never from an older copy.** Before editing, `git fetch origin && git rebase origin/main` (or start a fresh branch from `origin/main`). Never copy files into the repo from another checkout, worktree or backup: that silently reverts everyone else's changes to those files.
2. **Check the diff before every commit.** Run `git diff --stat origin/main` and make sure every file listed is one the task meant to change. If a file you did not intend to touch appears (especially `.github/`, `CMakeLists.txt`, `tests/`), restore it with `git checkout origin/main -- <file>` before committing.
3. **Do not modify `.github/workflows/` or test labels unless the task is about CI.** CI has two tiers, and both are deliberate:
   - **fast**: every push to `main` and every PR. Cached build, `ctest -LE deep`, no Python bindings. A few minutes.
   - **deep**: the full suite, including tests labelled `deep`. Runs nightly, on manual dispatch, and on PRs labelled `full-ci`.
4. **Label slow tests.** Any test that takes a minute or more goes in the `set_tests_properties(... PROPERTIES LABELS deep)` list in `CMakeLists.txt`. Before pushing, run the fast tier locally: `ctest --test-dir build -LE deep --output-on-failure`.
5. **Test data lives in this repo.** Load test inputs from `HYPHY3_SOURCE_DIR "/tests/data/..."` (or `benchmarks/data/`). Never reference a sibling checkout such as `../hyphy/`: CI does not have one. Copy the file into `tests/data/` instead.
6. **`assert()` is compiled out in Release builds** (CI builds Release). Checks a test relies on must throw or return a nonzero exit code, not only `assert()`.
7. **Keep the CI build cache working.** CI uses ccache, keyed on compiler flags and preprocessed source, and a cached build takes about a minute instead of about 30.
   - Never put a value that changes per commit or per run (commit hash, date, build number) into a compile definition or compiler flag. If it reaches every target, every commit becomes a full rebuild.
   - The commit hash reported in JSON output (`HYPHY3_GIT_COMMIT`) is the one exception, handled by the `if(NOT DEFINED HYPHY3_GIT_COMMIT)` guard in `CMakeLists.txt`. CI passes `-DHYPHY3_GIT_COMMIT=ci` so the definition never changes. Keep that guard. Local builds still get the real hash.
   - The `ci` pin never needs bumping: it is a fixed label, not a version. Leave it alone when you release or change versions.
   - To start CI from an empty cache (suspected stale or corrupt objects, a compiler change that ccache misses), bump the cache key prefix in `.github/workflows/test.yml`. Change both `key:` and `restore-keys:` for the tier, e.g. `ccache-fast-` → `ccache-fast-v2-`. Expect one slow build after that.
   - To check the cache, look for `Hits:` in the Build step's `ccache --show-stats` output. After the first build on a new key, a commit that only touches a few files should hit nearly every object. Hits near zero on every run mean something per-commit has reached the compile flags.
