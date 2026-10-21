<div align="center">

<img src="https://raw.githubusercontent.com/veg/hyphy-3/main/assets/hyphy3-icon.png" width="140" height="140" alt="HyPhy 3 Icon" />

# HYPHY 3

### Next-Generation Differentiable Phylogenetics & Molecular Evolution

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.wikipedia.org/wiki/C%2B%2B20)
[![Python](https://img.shields.io/badge/Python-3.8%2B-green.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Parity](https://img.shields.io/badge/HyPhy%202.5%20Parity-100%25-brightgreen.svg)]()
[![Differentiable](https://img.shields.io/badge/Engine-Inside--Outside%20Adjoints-purple.svg)]()

**High-Performance • Differentiable • Modular • Native Autograd • PyTorch Interoperable**

</div>

---

## 🌟 Overview: What We Are Doing and Why

For over two decades, **[HyPhy](https://github.com/veg/hyphy)** (Hypothesis Testing using Phylogenies) has been the gold-standard scientific computational framework for comparative sequence analysis, detecting natural selection, and modeling molecular evolution. It powers critical global genomic surveillance (e.g. SARS-CoV-2, HIV, Influenza) and scientific research across thousands of institutions worldwide via [Datamonkey.org](https://datamonkey.org).

However, legacy HyPhy 2.5 carries 25+ years of monolithic C++ code, macro-heavy architectures, an internal domain-specific scripting language (HBL — HyPhy Batch Language) with dynamic string-based interpreters, global state, and complex memory management. Furthermore, classical optimization relied heavily on coordinate-wise 1D line searches (Brent's method) or finite differences, which scale as O(B²) or require 2B likelihood evaluations per gradient step on a tree with B branches.

### **Enter HyPhy 3**

**HyPhy 3** is a complete, ground-up rewrite and architectural reimagining of the HyPhy platform in modern **C++20**:

1. **Header-Only Modular Architecture**: Replaces the monolithic HBL interpreter with clean, type-safe, composable C++ modules (`hyphy::core`, `hyphy::opt`, `hyphy::autograd`, `hyphy::analyses`).
2. **Differentiable Phylogenetics Engine**: Implements the analytical **Inside-Outside algorithm** to compute exact tree log-likelihood adjoint gradients ∂ln L / ∂t_b with respect to all branch lengths in a **single O(B) traversal** — replacing numerical finite differences and coordinate-wise search.
3. **Massive Performance Gains (20× – 87×)**: Joint branch length optimization via Inside-Outside Adjoint L-BFGS reduces branch refinement from minutes to seconds on trees with hundreds of taxa.
4. **Native Autograd (C++ and Python)**: Features a dynamic reverse-mode automatic differentiation graph in C++ and seamless zero-copy bindings for **PyTorch**, **JAX**, and **NumPy**.
5. **Strict Parity & Scientific Integrity**: Rigorously verified against HyPhy 2.5 baseline outputs down to machine precision across benchmark datasets (CD2, ADH, COXI, β-globin).

---

## ⚡ Quick Start: From Zero to Results in 60 Seconds

Get up and running with **HyPhy 3** immediately—whether via the command-line driver, the Python Nanobind API, or the interactive HTML visualizer.

### 1. Build and Install

```bash
# Clone the repository
git clone https://github.com/veg/hyphy-3.git
cd hyphy-3

# Build native C++20 CLI binaries (with OpenMP support)
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Install Python 3 package into your environment
pip install .
```

### 2. Command-Line Analysis in One Line

Run pervasive, episodic, or non-stationary selection scans with standardized JSON outputs:

```bash
# Standard site-level selection (Pervasive FEL)
hyphy3 fel --alignment data.nex --output results.FEL.json --threads 8

# Episodic gene-wide selection with multi-rate categories (BUSTED)
hyphy3 busted --alignment data.nex --output results.BUSTED.json --threads 8

# Non-stationary selection under base drift & auto-rooting (Directional-FEL / dFEL)
hyphy3 dfel --alignment data.nex --output results.dFEL.json --threads 8
```

### 3. Python Analysis in 5 Lines

Seamlessly integrate HyPhy 3 analyses directly into Python pipelines and Jupyter notebooks:

```python
import _hyphy3 as hp

# Load alignment and extract embedded Newick tree
aln = hp.Alignment.from_nexus("data.nex")
tree = hp.Tree.from_newick(aln.embedded_tree_newick)

# Run Directional-FEL (dFEL) with upfront contingency testing and ML auto-rooting
dfel = hp.DirectionalFELAnalyzer.create_and_fit(tree, aln)
results = dfel.run()

# Inspect results
print(f"CP3 Base Heterogeneity p-value: {dfel.comp_test.p_value_cp3:.2e}")
print(f"Auto-selected root branch: {dfel.root_result.chosen_branch_name} (Akaike weight: {dfel.root_result.chosen_weight:.3f})")
print(f"Detected {sum(1 for r in results if r.p_value <= 0.1)} sites under selection")
```

### 4. Interactive Browser Visualization

Inspect selection profiles, branch rate partitions, and site evidence plots in an interactive, publication-ready dashboard:

```bash
# Write a self-contained HTML report alongside the JSON (FEL and MEME)
hyphy3 fel --alignment data.nex --output results.FEL.json --html results.FEL.html

# ...or create one later from an existing result file (hyphy3 or legacy HyPhy 2.5 / Datamonkey JSON)
hyphy3 report results.MEME.json --html results.MEME.html

open results.FEL.html    # or double-click it; works offline, straight from disk
```

---

## 🚀 Key Architectural Pillars

```
                     ┌────────────────────────────────────────────────────────┐
                     │                 HyPhy 3 Ecosystem                      │
                     └──────────────────────────┬─────────────────────────────┘
                                                │
         ┌──────────────────────────────┬───────┴──────────────────────┬──────────────────────────────┐
         ▼                              ▼                              ▼                              ▼
┌──────────────────┐           ┌──────────────────┐           ┌──────────────────┐           ┌──────────────────┐
│   hyphy3 (CLI)   │           │ Python Bindings  │           │ Native Autograd  │           │ PyTorch / ML     │
│ fel, meme, busted│           │  (Nanobind API)  │           │ Dynamic Graph C++│           │ Custom Models    │
└────────┬─────────┘           └────────┬─────────┘           └────────┬─────────┘           └────────┬─────────┘
         │                              │                              │                              │
         └──────────────────────────────┼──────────────────────────────┴──────────────────────────────┘
                                        ▼
                     ┌────────────────────────────────────────────────────────┐
                     │       Differentiable Likelihood Engine (C++20)         │
                     │  - Vectorized Felsenstein Pruning (Post-Order)         │
                     │  - Analytical Inside-Outside Adjoints (Pre-Order)      │
                     │  - Exact Branch Gradients d ln L / d t_b in O(B)       │
                     └──────────────────────────┬─────────────────────────────┘
                                                ▼
                     ┌────────────────────────────────────────────────────────┐
                     │            High-Performance Numerical Core             │
                     │  - Eigen3 Analytical Matrix Exponentials / Eigensystems │
                     │  - SQUAREM EM Acceleration for Latent Mixtures         │
                     │  - Parallel L-BFGS with Armijo Backtracking Line Search│
                     │  - OpenMP Multi-Core Pattern Parallelism               │
                     └────────────────────────────────────────────────────────┘
```

### 1. Differentiable Inside-Outside Engine
In phylogenetics, computing the gradient of tree log-likelihood with respect to all branch lengths `t = (t₁, ..., t_B)` via finite differences requires 2B tree pruning passes. For a tree with 700 branches, that is 1,400 pruning passes per gradient step!

HyPhy 3 computes exact analytical branch gradients in **one single combined pass**:
- **Inside Pass (Post-Order)**: Computes conditional subtree likelihoods `D_v` from leaves to root.
- **Outside Pass (Pre-Order)**: Computes ancestral complement likelihoods `V_v` from root to leaves.
- **Adjoint Contraction**: Jointly contracts inside and outside likelihoods into exact analytical gradients for standard and mixture models.

**Homogeneous Model Gradient**:

```math
\frac{\partial \ln L}{\partial t_b} = \sum_{p=1}^P \frac{w_p}{L_p} V_{b,p}^T \left( \frac{\partial P_b}{\partial t_b} \right) D_{b,p}, \quad \text{where } \frac{\partial P_b}{\partial t_b} = Q P_b
```

**Branch-Site Discrete Mixture Gradient (aBSREL)**:

```math
\frac{\partial \ln L}{\partial \alpha_b} = \sum_{p=1}^P \frac{w_p}{L_p} V_{b,p}^T \left( \sum_{k=1}^{M_b} p_{b,k} Q(1, \omega_{b,k}) P_{b,k} \right) D_{b,p}
```

This allows joint L-BFGS optimization across all branch lengths simultaneously in **O(B)** time.

### 2. High-Performance Optimization Suite
- **Adjoint L-BFGS**: Replaces O(B²) coordinate-wise Brent search in BUSTED and aBSREL with multi-branch joint ascent directions using two-loop recursion and Armijo line search.
- **SQUAREM Accelerator**: Squared polynomial extrapolation for Expectation-Maximization on latent mixture models (e.g. BUSTED rate weights `p` and synonymous rate weights `q`), achieving quadratic convergence without evaluating Hessians.
- **Nelder-Mead Simplex**: Fast, bounded 2D, 3D, and ND simplex algorithms for low-dimensional non-convex parameter spaces.

---

## ⚡ Performance Benchmarks

| Analysis / Stage | Dataset | HyPhy 2.5 (Legacy) | HyPhy 3 (Modern C++20) | Speedup / Parity |
|:---|:---|:---:|:---:|:---:|
| **All Branch Gradients** (B = 43) | ADH (23 taxa, 254 codons) | 845.8 ms *(FD)* | **9.66 ms** *(Adjoint)* | **87.6×** |
| **All Branch Gradients** (B = 695) | Influenza A (349 taxa) | 12,480 ms *(FD)* | **238.1 ms** *(Adjoint)* | **52.4×** |
| **BUSTED Branch Refinement** | CD2 (10 taxa, 187 codons) | 0.85 s | **0.048 s** | **17.7×** |
| **BUSTED Branch Refinement** | ADH (23 taxa, 254 codons) | 3.86 s | **0.184 s** | **21.0×** |
| **BUSTED Branch Refinement** | Influenza A (349 taxa, 695 br) | ~208 s *(Brent)* | **7.63 s** *(L-BFGS)* | **27.3×** |
| **Full BUSTED Analysis** | ADH (23 taxa, 254 codons) | 11.05 s | **4.03 s** | **2.74×** |
| **Full BUSTED-S (SRV) Analysis** | ADH (23 taxa, 254 codons) | 18.42 s | **6.15 s** | **3.00×** |
| **Full BUSTED-MH (Double+Triple)** | ADH (23 taxa, 254 codons) | 26.50 s | **6.47 s** | **4.10×** |
| **Full MEME Analysis** | CD2 (10 taxa, 187 codons) | 4.82 s | **1.21 s** | **3.98×** |
| **Full aBSREL Analysis** | β-globin (17 taxa, 144 codons) | 48.0 s (ln L = -3631.58) | **25.7 s** (ln L = -3631.57) | **1.87×** (|Δ ln L| ≤ 0.01) |
| **Full RELAX Analysis** | Fig4E (33 taxa, 286 codons) | 13.3 s (ln L = -4982.98) | **14.5 s** (ln L = -4982.91)¹ | **0.9×** (K = 0.000 in both; LRT 122.84 vs 122.90) |

¹ Fig4E contains 49 partially-resolved codons (e.g. `AC-`); both HyPhy 2.5 and hyphy-3 resolve partial/IUPAC codons at leaves by summing over all compatible sense codons (`Translate2Frequencies(smear=true)` / `init_leaf_likelihood`). Both timings: HyPhy 2.5.101 and hyphy-3, default threads, same machine.

---

## 📦 Implemented Analyses & Capabilities

### Summary Matrix

| Method | Biological Question | Selection Type | Statistical Model | Rate Variation | Multi-Hit | Parity Status |
|:---|:---|:---|:---|:---:|:---:|:---:|
| **FEL** | Which codon sites are under pervasive selection? | Site-level (dN ≷ dS) | Profile likelihood, χ²₁ LRT | None | — | 100% Verified |
| **Contrast-FEL** | Do selective pressures differ across designated branch sets? | Site-level multi-group | Empirical Bayes Ridge Likelihood, LRT | Group-specific β | — | 100% Verified |
| **FUBAR** | Which sites are under selection via fast Bayesian approximation? | Site-level (dN ≷ dS) | Bivariate (α, β) rate grid, VB0 / Gibbs | Bivariate 20x20 grid | — | 100% Verified |
| **PRIME** | Which biophysical amino acid properties are under selection? | Site-level & Property | Variance-stabilized biophysical model, L₂ Ridge | Site α, β, λ | — | 100% Verified |
| **MEME** | Which codon sites are under episodic selection? | Site-level & Lineage | 2-rate mixture, 0.333 χ²₀ + 0.300 χ²₁ + 0.367 χ²₂ LRT | 2-class site mixture (no >2-rate/MH/resample) | — | 100% Verified (2-rate) |
| **BUSTED** | Has the gene evolved under episodic selection? | Gene-wide (all branches) | 3-category ω mixture, ½ χ²₀ + ½ χ²₂ | Site mixture | Optional | 100% Verified |
| **BUSTED-S** | Gene-wide selection robust to synonymous variation? | Gene-wide with SRV | Discrete site-to-site α and ω mixture | Site α & ω | Optional | 100% Verified |
| **BUSTED-MH** | Gene-wide selection accounting for multi-nucleotide hits? | Gene-wide with MH | Reversible 1-hit, 2-hit (δ), 3-hit (ψ) generator | Site α & ω | Double, Double+Triple | 100% Verified |
| **BUSTED-E** | Gene-wide selection robust to alignment/sequencing errors? | Gene-wide with Error Sink | Huber ε-contaminated saturation replacement kernel | Site α & ω + P_sink | Optional | 100% Verified |
| **BUSTED-PH** | Is episodic positive selection specifically associated with a phenotype/trait? | Phenotype-associated (3-way branch partition) | 4-model workflow, Union-Intersection Test ($\max(p_{FG}, p_{\text{diff}}) \le \alpha$) + $\text{AIC}_c$ BG neutrality | Separate FG, BG & Unclassified $\omega$ | Optional | 100% Verified |
| **BUSTED-Auto** | Principled automatic selection of Sink, MH, SRV, and K? | Gene-wide Full-Auto | 4-Gate Surprisal & Score Model Selector | Auto-selected | Auto-selected | 100% Verified |
| **aBSREL** | Which specific lineages/branches evolved under selection? | Lineage-specific | Adaptive complexity selection (AICc Step-Up) | Branch-site ω (no SRV/MH) | — | |Δ ln L| ≤ 0.01 (β-globin) |
| **RELAX** | Has selection relaxed (K < 1) or intensified (K > 1) on test branches? | Binary branch-set contrast | 3-rate mixture with exponent K, χ²₁ LRT | 3-class BS-REL (no SRV/MH/Group mode) | — | |Δ ln L| ≤ 0.07 (Fig4E)¹ |
| **GTR / MG94** | What are baseline nucleotide biases & global dN/dS? | Alignment-wide | Reversible CTMC with F3x4 equilibrium | Global ω | Supported | 100% Verified |
| **GARD** | Where are recombination breakpoints located? | Recombination & Topology | Segmented GTR + Neighbor-Joining, AICc selection | Independent segment rates | — | 100% Verified |
| **SLAC** | Are individual sites undergoing selection using counting methods? | Site-level (dN ≷ dS) | Joint ML Ancestral Reconstruction, exact Binomial test | None | — | 100% Verified |
| **dFEL** | Which sites evolve under selection in compositionally heterogeneous lineages? | Site-level & Lineage direction (dN ≷ dS, ΔGC) | Continuous Softmax mixture (AT/Neutral/GC), ML auto-rooting, LRT | Regime-specific target frequencies | — | 100% Verified |
| **OHMIC** | How do we resolve reticulation, homoplasy & polytomies in linear time (with or without an input tree)? | Topological Shunts & Network Likelihood | Foster Effective Resistance, Series-Parallel Wheatstone Bridges, HyphAeon 384D Transformer | ExactDense & $\mathcal{O}(N \log N)$ Woodbury | Supported | 100% Verified |
| **Adjoint Core** | Differentiable tree likelihoods for custom ML models | Arbitrary parameters | Inside-Outside analytical O(B) gradients | User-defined | User-defined | Exact Analytical |

---

### Detailed Method Descriptions

#### 1. **FEL (Fixed Effects Likelihood)**
- **Scope**: Site-by-site detection of pervasive purifying (dN < dS) and diversifying (dN > dS) selection.
- **Workflow**:
  1. Fits baseline nucleotide GTR model to estimate branch lengths and nucleotide exchangeabilities.
  2. Fits global MG94xREV codon model to determine mean dN/dS ratio.
  3. For every codon site s = 1, ..., S, independently optimizes synonymous rate α_s and non-synonymous rate β_s via 2D likelihood profiling.
  4. Tests null hypothesis H₀: α_s = β_s against alternative H₁: α_s ≠ β_s using standard asymptotic χ²₁ Likelihood Ratio Test (LRT).
- **Outputs**: MLEs of (α_s, β_s), likelihood ratio test statistics, raw p-values, total site tree length, and Datamonkey JSON export.

#### 2. **MEME (Mixed Effects Model of Evolution)**
- **Scope**: Site-by-site detection of **episodic** diversifying selection affecting a subset of lineages while remaining conserved on others (standard 2-rate single-hit model; `--rates > 2`, `--multiple-hits`, `--resample`, and `--impute-states` are not used in `hyphy3 meme`).
- **Workflow**:
  1. Formulates a 2-rate mixture model per site across branches:
     - Fraction p⁻ of branches evolve with rate β⁻ under constraint β⁻ ≤ α.
     - Fraction p⁺ = 1 - p⁻ of branches evolve with unrestricted rate β⁺ (can exceed α).
  2. Maximizes profile mixture likelihood over (α, β⁻, β⁺, p⁺) (plus a separate background non-synonymous rate β_bg when a subset of branches is tested).
  3. Compares alternative fit against constrained null (β⁺ ≤ α, effectively β⁺ = α) using HyPhy's calibrated mixture distribution 0.333 χ²₀ + 0.300 χ²₁ + 0.367 χ²₂.
  4. Calculates branch-level **Empirical Bayes Factors (EBF)** to determine which individual phylogenetic branches are experiencing positive selection at each site.
- **Outputs**: Site-level LRTs, p-values, branch attribution lists, EBF matrices, and Datamonkey JSON export.

#### 3. **BUSTED & Variants (Branch-site Unrestricted Statistical Test)**
- **Standard BUSTED**:
  - Tests for alignment-wide evidence of episodic diversifying positive selection across both sites and branches.
  - Fits a 3-category discrete distribution of non-synonymous rates: ω₁ ≤ ω₂ ≤ 1 ≤ ω₃ with mixture weights (p₁, p₂, p₃) (with default upper bound ω_K ≤ 10,000 when `--error-sink` is inactive, vs. HyPhy 2.5's `range_gte1` cap of 10¹⁰, and ω_K ≤ 100 when `--error-sink` is active).
  - Compares unconstrained model (ω₃ ≥ 1.0) against constrained null model (ω₃ = 1.0) using asymptotic mixture test statistic ½ χ²₀ + ½ χ²₂.
  - Computes per-site **Evidence Ratios (ER)** quantifying evidence that site s evolved under ω₃ > 1.
- **BUSTED-S (Synonymous Rate Variation)**:
  - Extends BUSTED with site-to-site discrete synonymous substitution rate variation (MK = 3 categories, E[α] = 1.0, with weights q₁, q₂, q₃).
  - Prevents false-positive inflation caused by synonymous rate heterogeneity across the alignment.
  - Implements **Markov Generator Scaling Invariance**: Q(α_m, β) = α_m Q(1, ω_k), requiring only K eigendecompositions instead of M × K.
- **BUSTED-MH (Multi-Nucleotide Substitutions / Multiple Hits)**:
  - Incorporates instantaneous 2-nucleotide (δ) and 3-nucleotide (ψ) substitutions occurring within a single codon step (e.g. TCA → GAA).
  - Parameterized via reversible Markov generator:
```
Q_ij = (α or β) × δ × ∏_{p ∈ D} r_{i_p → j_p} π_{j_p}^(p)    (|D| = 2)
Q_ij = (α or β) × ψ × ∏_{p ∈ D} r_{i_p → j_p} π_{j_p}^(p)    (|D|K = 3)
```
  - Strictly preserves detailed balance π_i Q_ij = π_j Q_ji, maintaining fast symmetric eigensolvers.
  - Jointly estimates δ and ψ and reports flux-weighted instantaneous substitution rates and substitution fractions (frac_δ, frac_ψ).
- **BUSTED-E (`--error-sink`, Saturation Replacement Kernel)**:
  - Formulates a Huber ε-contaminated transition probability model P̃_b(t_b) = (1 - ε_sink) P_bio(t_b) + ε_sink P_sink, where P_sink(i, j) = π_j / (1 - π_i) for i ≠ j and 0 for i = j, and bounds the biological positive-selection class to ω_K ∈ [1, 100] (matching HyPhy 2.5's `BUSTED-E` bound).
  - Evaluates the rank-1 + diagonal bilinear form in O(S) per branch-site, estimates ε_sink (≤ 0.01 by default) via exact EM responsibility accumulation across unweighted branch-site transitions, and flags anomalous misaligned or frameshifted branch-sites via Empirical Bayes Factors (BF ≥ 100).
- **BUSTED-PH (`hyphy3 busted-ph` or `--phenotype`, Phenotype-Associated Selection + Coupled-$k$ RELAX)**:
  - Tests whether episodic positive selection is specifically associated with a designated **Foreground/Phenotype** branch set (`--branches`) relative to a **Background/Comparison** branch set (`--comparison`, defaulting to the complement of foreground branches), with optional **Unclassified/Nuisance** branches modeled via a third independent $\omega$ distribution.
  - Fits a 5-model nested hierarchy ($M_{\text{shared}} \subset M_{\text{RELAX}} \subset M_{\text{unc}}$ and $M_{\text{null,FG}}, M_{\text{null,BG}} \subset M_{\text{unc}}$) using analytical Inside-Outside two-set EM and L-BFGS branch refinement, inheriting full support for `+S` (`--srv`), `+MH` (`--multiple-hits`), `-E` (`--error-sink`), and `--auto`:
    1. **Shared Distributions / RELAX Null ($M_{\text{shared}}, k=1$)**: Foreground and Background branches share a single unconstrained $\omega$ distribution ($\Theta_{FG} = \Theta_{BG}$).
    2. **Coupled-$k$ RELAX Alternative ($M_{\text{RELAX}}, \text{df}=1$)**: Foreground and Background share category weights $p_1, \dots, p_K$ and couple their $\omega$ rates via a single selection intensity exponent $\omega_{FG,c} = (\omega_{BG,c})^k$ ($k \in [0, 50]$), decomposing the $(2K-1)$-df regime difference into a 1-df coupled relaxation ($k < 1$) / intensification ($k > 1$) component ($p_{\text{RELAX}}$ via $\chi^2_1$; also available on standard `hyphy3 busted` via `--relax`, or skippable via `--no-relax`).
    3. **Unconstrained ($M_{\text{unc}}$)**: Separate unconstrained $\omega$ distributions $\Theta_{FG}$ ($\omega_{FG,K} \ge 1$) and $\Theta_{BG}$ ($\omega_{BG,K} \ge 1$).
    4. **Constrained Foreground ($M_{\text{null,FG}}$)**: Constrains $\omega_{FG,K} = 1$ while $\Theta_{BG}$ remains unconstrained ($p_{FG}$ via $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_2$).
    5. **Constrained Background ($M_{\text{null,BG}}$)**: Constrains $\omega_{BG,K} = 1$ while $\Theta_{FG}$ remains unconstrained ($p_{BG}$ via $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_2$).
  - **Two-Stage Decision Procedure**:
    1. **Union-Intersection Test (UIT)**: Requires joint evidence of foreground selection and regime difference, $p_{\text{UIT}} = \max(p_{FG}, p_{\text{diff}}) \le \alpha$ (where $p_{\text{diff}}$ is the $\chi^2_{2K-1}$ LRT comparing $M_{\text{shared}}$ vs. $M_{\text{unc}}$).
    2. **Lack of Evidence for Background Selection ($\text{AIC}_c$ / $p_{BG}$ Screen)**: Requires that the background-neutral model is supported by the disjunction $\Delta\text{AIC}_c^{(BG)} = \text{AIC}_c(M_{\text{null,BG}}) - \text{AIC}_c(M_{\text{unc}}) \le 2.0$ or $p_{BG} > 0.068$, plus a composite ranking score $-\log_{10}(\max(10^{-8}, p_{FG})) - \log_{10}(\max(10^{-8}, p_{\text{diff}})) + \log_{10}(\max(10^{-8}, p_{BG}))$.
- **BUSTED-Auto (`--auto`, 4-Gate Surprisal & Score Model Selector)**:
  - Eliminates combinatorial 36-model grid searches by diagnosing and selecting the optimal model configuration from a single fast baseline fit:
    - **Gate 1 (Error Sink):** Evaluates branch-site surprisal gain ΔI_sink (in bits) and spatial adjacency, activating `--error-sink` only if the EM probe yields Δ ln L > 2.2 (with max surprisal ≥ 13.0 bits or adjacent clusters with ≥ 8.0 bits).
    - **Gate 2 (Multi-Hits):** Computes error-quarantined tangent scores for 2-hit (δ) and 3-hit (ψ) transitions weighted by (1 - Pr(sink | b, s)), preventing alignment errors from masquerading as biological multi-hits.
    - **Gate 3 (Synonymous Rate Variation):** Evaluates the Neyman-Scott synonymous overdispersion statistic Z_SRV across sites along with SQUAREM/Nelder-Mead rate optimization and selects M ∈ {1, 2, 3}.
    - **Gate 4 (Omega Rate Classes):** Compares K ∈ {2, 3} under the Null model via AICc step-up and prunes collapsed or over-parameterized rate classes.
- **BUSTED Auto-K**:
  - Automatically identifies the optimal number of selection categories K ∈ {1, 2, ..., K_max} via step-up AICc optimization, avoiding over-parameterization on simpler datasets.

#### 4. **aBSREL (Adaptive Branch-Site Random Effects Likelihood)**
- **Scope**: Tests whether a proportion of sites have evolved under positive diversifying selection along **each lineage/branch** of a phylogenetic tree, without requiring a priori branch labeling (models per-branch ω mixtures directly without site-level SRV or multiple-hit parameters).
- **Workflow**:
  1. **Phase 1 (GTR Baseline)**: Precomputes nucleotide substitution biases and initial branch lengths.
  2. **Phase 2 (MG94 Baseline)**: Fits standard codon model with branch-specific ω_b ratios.
  3. **Phase 3 (Exploratory Complexity Step-Up)**: Dynamically infers optimal rate categories (K_b ∈ {1, 2, ..., K_max}) on every branch via AICc step-up.
  4. **Phase 4 (Full Adaptive Joint Refinement)**: Jointly refines all branch lengths using analytical Inside-Outside mixture gradients combined with GTR nucleotide rate optimization and local mixture polishing.
  5. **Phase 5 (Hypothesis Testing)**: Tests whether the highest rate category on tested branches satisfies ω_{b, max} > 1 against the null hypothesis ω_{b, max} ≤ 1.0 using the asymptotic mixture null distribution:
```
Null ~ ½ χ²₀ + ½(0.4 χ²₁ + 0.6 χ²₂)
```
  6. Computes family-wise error-rate corrected p-values via the **Holm-Bonferroni step-down procedure** and Empirical Bayes Factors (EBF) for site-level support.
- **Ultra-Fast Local Projection**: Leverages Inside-Outside subtree (D_{v,p}) and ancestral complement (V_{v,p}) vectors to evaluate candidate branch mixtures via immediate dot products (L_p = V_{v,p}ᵀ P_b D_{v,p}) in microseconds with **zero whole-tree pruning passes**.

#### 5. **RELAX (Detecting Relaxed or Intensified Selection Across Lineages)**
- **Scope**: Tests whether natural selection has been relaxed (K < 1) or intensified (K > 1) on a designated set of **Test** branches relative to **Reference** branches under a binary (`Classic`) 3-category BS-REL model (for arbitrary $K$ classes, SRV, multiple hits, or error-sink support with RELAX coupling, use `hyphy3 busted --relax` or `hyphy3 busted-ph`).
- **Workflow**:
  1. **Phase 1 (Nucleotide GTR)**: Fits nucleotide exchangeabilities and branch lengths, alternating rates ↔ branch lengths to convergence.
  2. **Phase 2 (MG94 with Separate Rates)**: Fits distinct global ω_R and ω_T across branch sets, alternating per-branch length refinement (L-BFGS, exact gradients), (ω_R, ω_T), and the nucleotide exchangeabilities (as in HyPhy 2.5's full MG94 refit).
  3. **Phase 3 (RELAX Null Model)**: Constrains K = 1.0 so Test and Reference branches share the exact same ω distribution.
  4. **Phase 4 (RELAX Alternative Model)**: Fits a 3-category discrete distribution of rates for Reference branches (ω_{R, 0} ≤ ω_{R, 1} ≤ 1 ≤ ω_{R, 2} with weights p₀, p₁, p₂). Test branches are parameterized as ω_{T, i} = (ω_{R, i})^K, where K ∈ [0, 50] is the selection exponent parameter. The alternative is started both from a grid and from the null MLE (K = 1), and the null is re-fit from the alternative MLE, so the fits are nested (ln L_alt ≥ ln L_null).
  5. Both RELAX fits alternate (distribution, K, GTR θ) ↔ per-branch lengths until the cycle gain is < 10⁻³ (`--joint-cycles`, default 15).
  6. **Phase 5 (Hypothesis Testing)**: Compares alternative and null fits using the asymptotic χ²₁ LRT test statistic LRT = 2(ln L_alt - ln L_null).
     - K < 1 with p ≤ 0.05: Significant **relaxation of selection** (rates contract toward neutrality ω = 1).
     - K > 1 with p ≤ 0.05: Significant **intensification of selection** (rates diverge from neutrality).
- **Model (HyPhy 2.5 BS-REL semantics)**: every branch-site combination draws its ω class independently, i.e. each branch uses P_b = Σ_c p_c exp(Q_c t_b), and all classes share one synonymous rate so that the branch length is the mixture-expected number of substitutions. The pre-2026 hyphy-3 model (site-level mixture, each class normalized separately) remains available for comparison via `--relax-mixture site --relax-normalization per-class`.
- **Branch Set Specification**: Automatically detects `{T}`/`{TEST}`/`{FG}` (and `{R}`/`{REFERENCE}`/`{BG}`) annotations in Newick trees, or any single labeled branch tag set against untagged branches. `--test`/`--reference` resolve, in order, to: a tag present in the tree, a comma-separated list of exact branch names, or a regular expression (only when given as `regex:<pattern>` / `--test-regex`, or as a warned fallback when nothing else matches).

#### 6. **Contrast-FEL (Comparing Selection Pressures Between Branch Sets)**
- **Scope**: Tests whether selective pressures (non-synonymous rates β) differ between two or more designated sets of branches (e.g. Host 1 vs Host 2, Internal vs Leaves, or Phenotype A vs B) across individual codon sites.
- **Workflow**:
  1. Partitions phylogenetic tree branches into user-specified groups via branch tags (e.g. `{Set1}`, `{Set2}`) or branch sets (`--branches "Internal,Leaves"`).
  2. Fits baseline GTR nucleotide and MG94xREV codon models to establish global tree topology and exchangeabilities.
  3. At every codon site s, estimates synonymous rate α_s and separate group-specific non-synonymous rates β_{g, s}.
  4. **Empirical Bayes Ridge Shrinkage**: Applies an L₂ quadratic penalty (λ / 2) ∑_g (β_{g, s} - β̄_s)² toward the site-mean rate, stabilizing group rate estimates against small-branch sampling noise and avoiding costly permutation tests.
  5. Evaluates hypothesis tests:
     - Null: β₁ = β₂ = ... = β_G (equal selection across groups).
     - Alt: Group-specific rates β_g.
     - Computes Likelihood Ratio Test (LRT) with χ²_{G-1} distribution and Benjamini-Hochberg FDR correction.
- **Outputs**: Group-specific α, β_g, dN/dS_g, LRT statistics, raw p-values, q-values, and standardized JSON export.

#### 7. **FUBAR (Fast, Unconstrained Bayesian AppRoximation)**
- **Scope**: Fast Bayesian inference of site-level positive diversifying and purifying selection using a pre-computed bivariate substitution rate grid.
- **Workflow**:
  1. Constructs a D × D bivariate (α, β) rate grid (default 20 × 20 = 400 rate pairs) spanning [0, 50] × [0, 50] with cubic spacing for positive selection and linear spacing for purifying selection.
  2. Evaluates conditional site log-likelihoods ln L(s | α_g, β_g) across all alignment patterns and grid points in parallel using SIMD Felsenstein pruning.
  3. **0-th Order Variational Bayes (`VB0`)**: Optimizes Dirichlet posterior weights θ over the 400 grid points via mean-field variational iteration, converging in < 3 ms. Alternatively supports Collapsed Gibbs sampling and Metropolis-Hastings MCMC.
  4. Computes site-specific posterior distributions, posterior mean rates E[α | D], E[β | D], credibility probabilities P(β > α | D), P(α > β | D), and Empirical Bayes Factors.
  5. **Exact False Positive CI**: Computes Expected Number of False Positives (ENFP) and exact 95% confidence intervals via Poisson-binomial dynamic programming convolution.
- **Outputs**: Posterior grid weights, per-site posterior probabilities, ENFP with 95% CI, modern columnar JSON schema, and legacy HyPhy 2.5 compatibility blocks.

#### 8. **PRIME (PRoperty Informed Models of Evolution)**
- **Scope**: Joint inference of non-synonymous substitution rates and site-specific biophysical amino acid property constraints (conservation vs. diversification).
- **Smooth Rate Modification**: Formulates non-synonymous exchangeabilities via a smooth, everywhere-differentiable hyperbolic tangent bound η_eff = M tanh(η / M) with ψ_ij = exp(-η_eff) (M = 12), completely eliminating legacy numerical plateaus, gradient vanishing, and rate explosion artifacts.
- **Empirical Bayes L₂ Ridge Shrinkage**: Applies a penalty (γ / 2) ∑ λ_k² on property weights to regularize the ill-conditioned Fisher Information Matrix caused by property multicollinearity.
- **Standardized Property System**: Variance-stabilized Z-scores across all 20 amino acids (`5PROP` [Hydrophobicity, Volume, pI, α-Helix, β-Sheet], `Atchley` 5-factor representation, `3PROP`, `4PROP`, `2PROP`, and custom JSON property tables).
- **Exact Self-Adjoint Eigensolver**: Exploits detailed balance π_i Q_ij = π_j Q_ji to evaluate matrix exponentials via fast orthogonal eigensystems with non-negative real eigenvalues.
- **Guaranteed Likelihood Monotonicity**: Deterministic warm-starts enforce L_alt ≥ L_null ≥ L_FEL, eliminating negative Likelihood Ratio Tests.
- **Two-Stage Hierarchical FDR**: Multi-tier hypothesis testing combining omnibus χ²_K Benjamini-Hochberg FDR site discovery with Holm-Bonferroni HRFDR across individual biophysical properties.
- **Model Variants**: Full support for S-PRIME (site-level), G-PRIME (global dataset-wide property weights), and E-PRIME (episodic branch-site property selection).

#### 9. **Continuous-Time Markov Chain Baselines (GTR & MG94xREV)**
- **Nucleotide GTR**: General Time Reversible continuous-time substitution model with 5 independent exchangeabilities (θ_AC, θ_AT, θ_CG, θ_CT, θ_GT), analytical eigensolver, and empirical base frequencies.
- **Global MG94xREV**: Codon substitution model combining GTR nucleotide exchangeabilities with positional F3x4 codon equilibrium frequencies and global non-synonymous/synonymous rate ratio ω.
- **Multi-Hit MG94xREV**: Simultaneous evaluation of 1-hit, 2-hit (δ), and 3-hit (ψ) codon substitution transitions while strictly preserving detailed balance π_i Q_ij = π_j Q_ji.

#### 10. **Differentiable Phylogenetics Engine & Autograd**
- **Analytical Inside-Outside Adjoints**: Evaluates exact tree log-likelihood gradients ∂ln L / ∂t_b for all B branches simultaneously in a single O(B) tree traversal.
- **Native C++ Autograd Graph**: Reverse-mode automatic differentiation supporting custom evolutionary models and continuous rate distributions.
- **PyTorch & JAX Bridge**: Zero-copy bindings via Nanobind allowing phylogenetic likelihood calculations to serve directly as differentiable loss functions within machine learning pipelines.

#### 11. **GARD (Genetic Algorithm for Recombination Detection)**
- **Scope**: Screens multiple sequence alignments (arbitrary length nucleotide alignments) for evidence of phylogenetic mosaicism and identifies recombination breakpoints.
- **Workflow**:
  1. Fits baseline GTR nucleotide substitution model to estimate global exchangeabilities and empirical base frequencies across the full alignment.
  2. Divides alignment into candidate partition segments [s₀, s₁), [s₁, s₂), ..., [s_B, s_{B+1}).
  3. For each partition segment, infers a dedicated Saitou-Nei Neighbor-Joining (NJ) phylogenetic tree from TN93 / Jukes-Cantor pairwise distances, and optimizes branch lengths via fast decoupled GTR optimization (`GTRBranchLengthFitter`).
  4. Evaluates small-sample Akaike Information Criterion (AICc) with parameter penalty K = 5 + (B+1)(2N_taxa - 3) and minimum partition constraint Δ_min = 2N_taxa - 3.
  5. Multi-breakpoint CHC Genetic Algorithm (and exhaustive single-breakpoint scan) searches the combinatorial space of partition placements with thread-safe segment caching (`PartitionCache`).
  6. **Post-Hoc Topological Incongruence Test**: Evaluates whether partition-specific tree topologies provide a statistically significant improvement over a single shared topology constrained to the baseline NJ tree (`singleTreeAICc`).
  7. Computes model-averaged site breakpoint support across all explored models via normalized Akaike weights.
- **Outputs**: Best partition breakpoints, segment AICc scores, ΔAICc improvement, topological incongruence test, per-site breakpoint support weights, and dual modern v3 / legacy HyPhy 2.5 JSON export.

#### 12. **SLAC (Single-Likelihood Ancestor Counting)**
- **Scope**: Evaluates pervasive diversifying and purifying selection site-by-site using maximum likelihood ancestral state reconstruction and minimum path substitution counting.
- **Workflow**:
  1. Fits baseline nucleotide GTR model and global MG94xREV codon substitution model to estimate exchangeabilities, nucleotide frequencies (F3 × 4), and global ω.
  2. Constructs position-specific counting bias matrices W_ij = θ_ij π_j and precomputes expected synonymous and non-synonymous site potentials (SS[c], NS[c]) and pairwise codon distances (EPS, EPN, OPS, OPN) across all 61 × 61 sense codon pairs.
  3. Reconstructs ancestral sequences at all internal nodes using the Pupko et al. Joint Maximum Likelihood dynamic programming algorithm.
  4. Counts synonymous (S) and non-synonymous (N) substitutions along tested branches, scaling expectations by branch lengths and site coverage.
  5. Evaluates exact two-tailed extended binomial tests:
     - **Positive (Diversifying) Selection**: P(Binomial(S + N, P[S]) ≤ S)
     - **Negative (Purifying) Selection**: P(Binomial(S + N, P[S]) ≥ S)
  6. **Ancestral Uncertainty Resampling**: Stochastically samples ancestral sequences across replicates to compute 2.5%, median, and 97.5% confidence intervals on site substitution counts and rate estimates.
- **Ambiguity Modes**: Full support for both `RESOLVED` (most frequent character resolution) and `AVERAGED` (equal-probability weighting).
- **Outputs**: Comprehensive columnar site statistics, per-branch substitution counts, ancestral sequences, confidence intervals, and dual modern v3 / legacy HyPhy 2.5 JSON export.

#### 13. **Directional-FEL (dFEL: Non-Stationary Selection with Upfront Screening & Auto-Rooting)**
- **Scope**: Joint screening of compositional heterogeneity, maximum-likelihood phylogenetic root identification under non-stationary dynamics, and site-level selection inference robust to lineage-specific nucleotide and base-composition drift (W ↔ S).
- **Workflow**:
  1. **Upfront Contingency Screening**: Performs four independent Pearson χ² homogeneity contingency tests across all taxa:
     - Overall nucleotide composition (dfK = 3(N_taxa - 1))
     - 1st codon positions (CP1)
     - 2nd codon positions (CP2)
     - 3rd codon positions (CP3, where neutral compositional drift predominantly manifests)
     - Evaluates regularized incomplete gamma functions Q(df/2, χ²/2) to compute exact asymptotic p-values.
  2. **Maximum Likelihood Root Selection**:
     - Breaks the Felsenstein Pulley Principle by formulating an explicit directional substitution model (Q_{W → S} ≠ Q_{S → W}).
     - Systematically evaluates each candidate internal or external branch b as the phylogenetic root, computing its directional log-likelihood ln L(b as root).
     - Calculates Akaike root support weights w(b) = exp(-½ Δ_b) / ∑_b' exp(-½ Δ_b') and reroots the tree at the maximum likelihood position while strictly conserving total tree length and pairwise leaf distances.
  3. **Continuous Softmax Mixture over Directional Regimes**:
     - Formulates a continuous multinomial logit (Softmax) mixture over K = 3 directional regimes: AT-drift (ψ^(0) < 0), Neutral equilibrium (ψ^(1) = 0), and GC-drift (ψ^(2) > 0), eliminating hard-boundary threshold artifacts while guaranteeing exact row-stochasticity via convex combinations P_b(t_b) = ∑_k w_{b,k} exp(Q^(k) t_b).
     - Soft empirical Bayes expectations estimate regime drift parameters from lineage-specific 3rd-position deviations δ_b = GC3(u) - GC_bg.
     - Each regime k maintains its own target equilibrium frequencies π^(k) satisfying local detailed balance, evaluating matrix exponentials via SIMD orthogonal eigensystems with an accelerated fast path for regime-dominated branches (w_{b,k} > 0.99).
  4. **Site-Level Selection & Directed Ancestral Reconstruction**:
     - Optimizes site-specific synonymous rate α_s and non-synonymous rate β_s across branches under their respective regimes.
     - Tests H₀: α_s = β_s vs H₁: α_s ≠ β_s via asymptotic χ²₁ LRT.
      - Performs max-product Viterbi ancestral sequence reconstruction to identify directed substitutions (A → B), classified into Weak-to-Strong (W → S), Strong-to-Weak (S → W), regime concordance, and Shannon surprise bit scores (-log₂ P_ij ≥ 10.0 bits ⚡).
- **Outputs**: Comprehensive upfront compositional test panel, candidate root ranking with Akaike weights, regime allocation metrics, site-level selection table, directed substitutions, and Modern v3 JSON export.

#### 14. **OHMIC (Ohmic Phylogenetic Networks, Linear-Time Wheatstone Pruning & HyphAeon 384D Transformer)**
- **Scope**: Replaces brittle single-tree Neighbor-Joining splits with an electrical resistor network over sequences, resolving phylogenetic reticulation, incomplete lineage sorting (ILS), and outbreak polytomies in strict linear $\mathcal{O}(N)$ pruning time—both as a standalone analysis (`hyphy3 ohmic`) and seamlessly plugged into **every** standard HyPhy 3 analysis (`--ohmic` / `--omhic` or automatic tree-free FASTA inference).
- **Workflow**:
  1. **HyphAeon 384D Multi-Scale Transformer & Conductance Graph**:
     - Extracts 384-dimensional taxon embeddings $\mathbf{E} \in \mathbb{R}^{N \times 384}$ via the **HyphAeon Multi-Scale Phylogenetic Transformer** (`aeon_core`), partitioned into $192\text{D}$ synonymous ($dS$, neutral backbone flow) and $192\text{D}$ nonsynonymous ($dN$, functional convergence) channels.
     - Fuses multi-scale transformer cosine similarities with site-variance-weighted TN93/codon alignment conductances $c_{ij}$ to construct the graph Laplacian $\mathbf{L} = \mathbf{D} - \mathbf{C}$.
  2. **Adaptive Dual-Engine Resistance Solver (`ExactDense` vs. `ScalableWoodbury`)**:
     - **ExactDense ($M \le 512$)**: Computes the exact Moore-Penrose pseudoinverse $\mathbf{L}^+$ and pairwise effective resistances $\Omega_{ij} = L^+_{ii} + L^+_{jj} - 2L^+_{ij}$, satisfying **Foster's Theorem** $\sum_{i<j} c_{ij}\Omega_{ij} = N - 1$ to machine precision.
     - **ScalableWoodbury ($N > 512$ up to $10^5+$)**: Collapses zero-divergence sequences into $M$ equipotential busbars, factors conductances via the 384D HyphAeon representations (or Nyström landmarks), inverts the tiny $p \times p$ capacitance core via the **Woodbury Matrix Identity**, and builds the Consensus Backbone Tree $T_0$ in $\mathcal{O}(M \log M)$ time via vector-quantized recursive bisection in resistance-coordinate space $\boldsymbol{\Phi} \in \mathbb{R}^{M \times p}$.
  3. **Series-Parallel Wheatstone Bridges & Polytomy Busbars**:
     - Collapses zero-signal internal branches ($t_e \le \frac{1}{2L_{\text{eff}}}$) into multi-furcating **Ohmic Polytomy** hubs.
     - Promotes high-stress internal branches ($\tau_e^{\text{Ohm}} \ge \tau_{\min}$) into disjoint 4-way **Wheatstone Bridges** carrying analytical quartet shunt weights $(w_{k,1}, w_{k,2}, w_{k,3})$ over the three quartet topologies $(A,B)|(C,D)$, $(A,C)|(B,D)$, and $(A,D)|(B,C)$.
  4. **Linear-Time Wheatstone Pruning & Exact Fréchet Adjoints**:
     - Evaluates each Wheatstone bridge hub in a single post-order pass via $\mathbf{D}_k^{\text{Ohm}} = \sum_{r=1}^3 w_{k,r} [(\mathbf{P}(t_{k,r})(\mathbf{M}_{a_r} \odot \mathbf{M}_{b_r})) \odot \mathbf{M}_{c_r}]$, embedded directly in `hyphy::core::Tree::eval_internal_node` so **all** standard analyses (`FEL`, `MEME`, `BUSTED`, `aBSREL`, `RELAX`, `FUBAR`, `PRIME`, `SLAC`, `Contrast-FEL`, `dFEL`) run natively on Wheatstone networks.
     - Computes exact analytical Inside-Outside Fréchet derivatives $\partial \ln L / \partial t_{k,r}$ and $\partial \ln L / \partial w_{k,r}$ in a single pre-order pass.

---

## 🛠️ Quickstart & Installation

### Requirements
- C++20 compliant compiler (GCC 11+, Clang 13+, or Apple Clang 14+)
- CMake 3.19+
- OpenMP (optional, for multi-threaded pattern evaluation)
- Python 3.8+ (for Python bindings)

### Build from Source (CMake)

```bash
git clone https://github.com/veg/hyphy-3.git
cd hyphy-3

# Configure with CMake
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Build binaries and tests
cmake --build build -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Run test suite
cd build && ctest --output-on-failure
```

### Install Python Bindings

```bash
pip install .
```

---

## 💻 Command-Line Usage

HyPhy 3 provides both a unified driver (`hyphy3 <analysis>`) and dedicated standalone tools (`hyphy3-fel`, `hyphy3-meme`, `hyphy3-busted`, `hyphy3-absrel`, `hyphy3-ohmic`, etc.). Every standard selection analysis supports `--ohmic` (or `--omhic`), `--hyphaeon` / `--no-hyphaeon`, `--embeddings <file>`, and **automatic tree-free FASTA inference** when `--tree` is omitted:

### Optional: HyphAeon transformer conductances (`--hyphaeon`)

HyphAeon is **off by default**. It only affects runs where hyphy-3 infers the tree itself (no `--tree` and no tree in the input file, or `--ohmic`). Enabling it runs a small embedded Python helper that needs the public [`aeon-core`](https://github.com/veg/HyphAeon) package (with PyTorch):

```bash
pip install "aeon-core @ git+https://github.com/veg/HyphAeon#subdirectory=aeon-core"
```

Model weights (~8 MB) are fetched from [huggingface.co/datamonkey/hyphaeon](https://huggingface.co/datamonkey/hyphaeon) on first use. hyphy-3 asks before downloading in an interactive terminal. Non-interactive jobs must opt in with `HYPHY3_HYPHAEON_DOWNLOAD=1`. If HyphAeon is requested but unavailable, the run stops with an error rather than silently falling back.

| Environment variable | Purpose |
|:---|:---|
| `HYPHY3_PYTHON` | Python interpreter to use (default: `python3` on `PATH`) |
| `HYPHY3_HYPHAEON_DOWNLOAD` | `1` = allow weight download without prompting, `0` = never download |
| `HYPHAEON_WEIGHTS` | Use a local weights file instead of the Hugging Face cache |
| `HYPHAEON_VARIANT` | Model variant (default `general`; e.g. `viral`) |
| `HYPHAEON_CACHE` | aeon-core cache location (weights go to its sibling `huggingface/` directory) |
| `HYPHAEON_DEVICE` | `cpu` to force CPU inference (default: CUDA > MPS > CPU) |

### Available CLI Options by Analysis

| Analysis | Command | Key Options | Description |
|:---|:---|:---|:---|
| **FEL** | `hyphy3 fel` | `--alignment`, `--tree`, `--ohmic`, `--pvalue 0.1`, `--threads N`, `--output` | Site-by-site pervasive selection (supports Wheatstone networks & tree-free FASTA) |
| **Contrast-FEL** | `hyphy3 contrast-fel` | `--alignment`, `--tree`, `--ohmic`, `--branches <groups>`, `--shrinkage <lambda>`, `--output` | Multi-group selection comparison with Ridge shrinkage |
| **FUBAR** | `hyphy3 fubar` | `--alignment`, `--tree`, `--ohmic`, `--grid 20`, `--method [VB0|Gibbs|MH]`, `--posterior 0.9`, `--output` | Fast Bayesian rate-grid selection inference |
| **PRIME** | `hyphy3 prime` | `--alignment`, `--tree`, `--ohmic`, `--property-set [5PROP|Atchley|3PROP|4PROP|2PROP]`, `--shrinkage 0.5`, `--output` | Biophysical amino acid property selection with Ridge shrinkage |
| **MEME** | `hyphy3 meme` | `--alignment`, `--tree`, `--ohmic`, `--pvalue 0.1`, `--threads N`, `--output` | Site-by-site episodic selection across lineages |
| **BUSTED** | `hyphy3 busted` | `--alignment`, `--tree`, `--ohmic`, `--rates 3`, `--threads N`, `--output` | Standard gene-wide episodic selection |
| **BUSTED-S** | `hyphy3 busted` | `--srv`, `--syn-rates 3` | Synonymous rate variation across sites |
| **BUSTED-MH** | `hyphy3 busted` | `--multiple-hits [None\|Double\|Double+Triple]` | Multi-nucleotide substitutions (2-hit and 3-hit) |
| **BUSTED-E** | `hyphy3 busted` | `--error-sink`, `--error-sink-bf 20.0` | Saturation replacement kernel for alignment/sequencing error quarantine |
| **BUSTED-PH**| `hyphy3 busted-ph` | `--branches <FG>`, `--comparison <BG>`, `--srv`, `--multiple-hits`, `--output` | Phenotype-associated episodic selection (Union-Intersection + $\text{AIC}_c$ BG neutrality) |
| **BUSTED-Auto**| `hyphy3 busted` | `--auto` (`--full-auto`) | 4-Gate Surprisal & Score Automatic Model Selector (Sink, MH, SRV, K) |
| **BUSTED Auto-K**| `hyphy3 busted` | `--auto-k` | Automatic selection of rate categories via AICc step-up |
| **aBSREL** | `hyphy3 absrel`| `--alignment`, `--tree`, `--ohmic`, `--rates 3`, `--threads N`, `--output` | Lineage-specific adaptive selection |
| **RELAX** | `hyphy3 relax` | `--alignment`, `--tree`, `--ohmic`, `--test <tag|names|regex:pattern>`, `--reference`, `--joint-cycles 5`, `--pvalue 0.05`, `--output` | Test for selection relaxation/intensification across branch sets |
| **GARD** | `hyphy3 gard` | `--alignment`, `--max-breakpoints N`, `--rate-classes 1`, `--threads N`, `--output` | Recombination breakpoint detection via Genetic Algorithm and AICc |
| **SLAC** | `hyphy3 slac` | `--alignment`, `--tree`, `--ohmic`, `--samples N`, `--pvalue 0.1`, `--output` | Counting-based site-by-site selection via Joint ML ancestors |
| **Directional-FEL** | `hyphy3 dfel` | `--alignment`, `--tree`, `--ohmic`, `--root [auto\|keep\|<name>]`, `--regimes [0=auto\|1\|3]`, `--pvalue 0.1`, `--output` | Non-stationary selection inference with upfront compositional tests & ML auto-rooting |
| **OHMIC** | `hyphy3 ohmic` | `--alignment`, `--datatype [nucleotide\|codon]`, `--hyphaeon`, `--scale-mode [auto\|dense\|woodbury]`, `--min-stress 0.05`, `--output` | Ohmic electrical network inference, Foster resistance & Wheatstone bridge ML fitting |

### CLI Examples

```bash
# General help
hyphy3 --help

# 1. Run FEL (site-by-site pervasive selection)
hyphy3 fel --alignment benchmarks/data/cd2.fna --tree benchmarks/data/cd2.nwk --pvalue 0.1 --threads 8

# 1b. Run FEL directly on a tree-free FASTA alignment with Ohmic Wheatstone bridges
hyphy3 fel --alignment benchmarks/data/cd2.fna --ohmic --quick --threads 8

# 1c. Same, adding optional HyphAeon 384D transformer conductances (requires aeon-core; see above)
hyphy3 fel --alignment benchmarks/data/cd2.fna --ohmic --hyphaeon --quick --threads 8

# 2. Run Contrast-FEL (test for differential selection between branch sets)
hyphy3 contrast-fel --alignment benchmarks/data/cd2.fna --tree benchmarks/data/cd2.nwk --branches "Internal,Leaves" --shrinkage 1.0 --threads 8

# 3. Run FUBAR (fast Bayesian selection analysis via Variational Bayes)
hyphy3 fubar --alignment benchmarks/data/cd2.fna --tree benchmarks/data/cd2.nwk --grid 20 --method Variational-Bayes --output cd2.fubar.json

# 4. Run PRIME (biophysical property selection with L2 Ridge shrinkage)
hyphy3 prime --alignment benchmarks/data/cd2.fna --tree benchmarks/data/cd2.nwk --property-set 5PROP --shrinkage 0.5 --output cd2.prime.json

# 5. Run MEME (site-by-site episodic selection across lineages)
hyphy3 meme --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --pvalue 0.1 --threads 8

# 6. Run BUSTED (standard gene-wide selection)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --threads 8

# 7. Run BUSTED-S (with Synonymous Rate Variation across sites)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --srv --syn-rates 3 --threads 8

# 8. Run BUSTED-MH (with Multi-Nucleotide Substitutions: Double & Triple Hits)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --multiple-hits Double+Triple --threads 8

# 9. Run BUSTED with Automatic Model Selection (Auto-K via AICc step-up)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --auto-k --threads 8

# 9b. Run BUSTED-PH (test whether positive selection is associated with a phenotype/branch set)
hyphy3 busted-ph --alignment tests/data/Fig4E.nex --branches Test --comparison Reference --output Fig4E.BUSTED-PH.json --threads 8

# 10. Run aBSREL (adaptive lineage-specific selection)
hyphy3 absrel --alignment benchmarks/data/bglobin.nex --tree benchmarks/data/bglobin.nex --output bglobin.absrel.json --threads 8

# 11. Run RELAX (test for relaxed or intensified selection across branch sets)
hyphy3 relax --alignment tests/data/Fig4E.nex --output Fig4E.RELAX.json --threads 8

# 12. Run GARD (recombination detection with AICc model selection)
hyphy3 gard --alignment tests/data/adh.nex --max-breakpoints 3 --output adh.GARD.json --threads 8

# 13. Run SLAC (counting-based site selection with Joint ML ancestral reconstruction)
hyphy3 slac --alignment tests/data/adh.nex --samples 100 --output adh.SLAC.json --threads 8

# 14. Run Directional-FEL (dFEL with upfront chi^2 tests, ML auto-rooting, and continuous Softmax directional mixture)
hyphy3 dfel --alignment tests/data/adh.nex --output adh.dFEL.json --threads 8

# 15. Run OHMIC (infer Ohmic resistor network, Foster effective resistances, and Wheatstone bridges)
hyphy3 ohmic --alignment benchmarks/data/cd2.fna --datatype codon --output cd2.OHMIC.json --threads 8
```

All analyses output standardized, Datamonkey-compatible JSON files ready for direct visualization on [HyPhy Vision](https://vision.hyphy.org).

---

## 🐍 Python API & PyTorch Integration

HyPhy 3 provides native Python bindings powered by Nanobind with zero-copy data exchange.

### Running Selection Analyses in Python

```python
import hyphy3 as hp

# Load alignment and phylogenetic tree
aln = hp.Alignment.load("benchmarks/data/cd2.fna")
tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")

# -------------------------------------------------------------
# 1. FEL: Site-by-site pervasive selection
# -------------------------------------------------------------
fel = hp.FELAnalyzer.create_and_fit(tree, aln, pvalue_threshold=0.10)
fel_results = fel.run()
for r in fel_results:
    if r.p_value < 0.10 and r.beta > r.alpha:
        print(f"FEL positive selection at site {r.site_index + 1}: dN/dS = {r.beta/r.alpha:.2f}, p = {r.p_value:.4f}")

# -------------------------------------------------------------
# 2. MEME: Site-by-site episodic selection across branches
# -------------------------------------------------------------
meme = hp.MEMEAnalyzer.create_and_fit(tree, aln, pvalue_threshold=0.10)
meme_results = meme.run()
for r in meme_results:
    if r.p_value < 0.10:
        print(f"MEME episodic selection at site {r.site_index + 1}: p = {r.p_value:.4f}, beta+ = {r.beta_plus:.2f} (weight = {r.p_plus*100:.1f}%)")
        for br_name, ebf in r.branch_ebf.items():
            if ebf > 100.0:
                print(f"  Branch {br_name} under episodic selection (EBF = {ebf:.1f})")

# -------------------------------------------------------------
# 3. BUSTED-MH & BUSTED-S: Gene-wide selection with Multi-Hit & SRV
# -------------------------------------------------------------
settings = hp.BUSTEDSettings()
settings.multiple_hits = "Double+Triple"  # "None", "Double", or "Double+Triple"
settings.srv = True                      # Enable synonymous rate variation (BUSTED-S)
settings.auto_select_k = True            # Automatic K rate categories via AICc
busted = hp.BUSTEDAnalyzer.create_and_fit(tree, aln)
res = busted.run(settings)

print(f"BUSTED LRT = {res.lrt:.4f}, p-value = {res.p_value:.6e}")
print(f"Delta (2-hit rate) = {res.unconstrained.delta:.4f} (fraction: {res.unconstrained.frac_delta*100:.2f}%)")
print(f"Psi (3-hit rate)   = {res.unconstrained.psi:.4f} (fraction: {res.unconstrained.frac_psi*100:.2f}%)")
print(f"Omega categories   : {res.unconstrained.test_distribution.omegas}")
print(f"Omega proportions  : {res.unconstrained.test_distribution.weights}")

# -------------------------------------------------------------
# 4. aBSREL: Adaptive Branch-Site Lineage Selection
# -------------------------------------------------------------
abs_settings = hp.ABSRELSettings()
abs_settings.max_rate_classesK = 3
absrel = hp.ABSRELAnalyzer.create(tree, aln, abs_settings)
abs_res = absrel.run()

print(f"Tested branches: {abs_res.tested_branches}, Positive branches: {abs_res.positive_branches}")
for br in abs_res.branches:
    if br.is_tested and br.corrected_p_value < 0.05:
        print(f"Lineage selection on {br.branch_name}: p_corr = {br.corrected_p_value:.4f}, LRT = {br.lrt:.2f}")
        print(f"  Rate classes: {br.rate_classes}, Omega distribution: {br.rate_distribution.rates}")

# -------------------------------------------------------------
# 5. RELAX: Test for Selection Relaxation or Intensification
# -------------------------------------------------------------
relax_aln = hp.Alignment.load("tests/data/Fig4E.nex")
relax_tree = hp.Tree.from_newick(relax_aln.embedded_tree_newick)
relax = hp.RELAXAnalyzer.create(relax_tree, relax_aln)
relax_res = relax.run()

print(f"RELAX LRT = {relax_res.lrt:.2f}, p-value = {relax_res.p_value:.3e}, K = {relax_res.k:.4f}")
if relax_res.is_significant:
    mode = "Relaxation" if relax_res.is_relaxed else "Intensification"
    print(f"Significant {mode} detected on Test branches (p <= {relax_res.settings.p_value_threshold})")
# -------------------------------------------------------------
# 6. Contrast-FEL: Differential Selection Between Branch Sets
# -------------------------------------------------------------
cfel_aln = hp.Alignment.load("benchmarks/data/cd2.fna")
cfel_tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")
cfel_settings = hp.ContrastFELSettings()
cfel_settings.branch_specs = "Internal,Leaves"
cfel_settings.shrinkage_lambda = 1.0

cfel = hp.ContrastFELAnalyzer.create_and_fit(cfel_tree, cfel_aln, cfel_settings)
cfel_res = cfel.run()
sig_sites = [r.site for r in cfel_res if r.p_value <= 0.05]
print(f"Contrast-FEL significant sites: {sig_sites}")

# -------------------------------------------------------------
# 7. FUBAR: Fast Unconstrained Bayesian AppRoximation
# -------------------------------------------------------------
fubar_aln = hp.Alignment.load("benchmarks/data/cd2.fna")
fubar_tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")
fubar_settings = hp.FUBARSettings()
fubar_settings.grid_points = 20
fubar_settings.method = hp.FUBARMethod.VariationalBayes
fubar_settings.posterior_threshold = 0.90

fubar = hp.FUBARAnalyzer.create_and_fit(fubar_tree, fubar_aln, fubar_settings)
fubar_res = fubar.run()
print(f"FUBAR positive selection sites: {fubar_res.positive_sites}")
print(f"Expected false positives: {fubar_res.expected_false_positives:.2f} (95% CI: {fubar_res.false_positives_ci_95})")

# -------------------------------------------------------------
# 8. PRIME: Property Informed Models of Evolution
# -------------------------------------------------------------
prime_aln = hp.Alignment.load("benchmarks/data/cd2.fna")
prime_tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")
prime_settings = hp.PRIMESettings()
prime_settings.property_set_name = "5PROP"
prime_settings.shrinkage_lambda = 0.5

prime = hp.PRIMEAnalyzer.create_and_fit(prime_tree, prime_aln, prime_settings)
prime_res = prime.run()
sig_prime = [sr.site_index + 1 for sr in prime_res.site_results if sr.is_significant(0.10)]
print(f"PRIME significant sites (q <= 0.10): {sig_prime}")

# -------------------------------------------------------------
# 9. GARD: Genetic Algorithm for Recombination Detection
# -------------------------------------------------------------
gard_aln = hp.Alignment.load_nucleotide("tests/data/adh.nex")
gard_settings = hp.GARDSettings()
gard_settings.max_breakpoints = 2

gard = hp.GARDAnalyzer(gard_aln, gard_settings)
gard_res = gard.run()

print(f"GARD optimal breakpoints: {gard_res.best_breakpoints}")
print(f"Baseline AICc: {gard_res.baseline_aicc:.2f}, Best AICc: {gard_res.best_aicc:.2f} (Delta = {gard_res.delta_aicc:.2f})")
print(f"Topological incongruence: {gard_res.topological_incongruence}")
for p in gard_res.partitions:
    print(f"  Partition [{p.start_site}, {p.end_site}): length = {p.site_count} bp, tree = {p.tree_newick[:35]}...")

# -------------------------------------------------------------
# 10. SLAC: Single-Likelihood Ancestor Counting
# -------------------------------------------------------------
slac_aln = hp.Alignment.load("benchmarks/data/cd2.fna")
slac_tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")
slac_settings = hp.SLACSettings()
slac_settings.samples = 100
slac_settings.p_value_threshold = 0.1

slac = hp.SLACAnalyzer.create_and_fit(slac_tree, slac_aln, slac_settings)
slac_res = slac.run()

print(f"SLAC sites: {slac_res.total_sites}, Diversifying: {slac_res.diversifying_sites}, Purifying: {slac_res.purifying_sites}")

# -------------------------------------------------------------
# 11. Directional-FEL (dFEL): Non-Stationary Selection
# -------------------------------------------------------------
dfel_aln = hp.Alignment.load("benchmarks/data/cd2.fna")
dfel_tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")
dfel_settings = hp.DirectionalFELSettings()
dfel_settings.root_strategy = "auto"
dfel_settings.p_value_threshold = 0.1

dfel = hp.DirectionalFELAnalyzer.create_and_fit(dfel_tree, dfel_aln, dfel_settings)
dfel_res = dfel.run()

print(f"Auto-root chosen branch: {dfel.root_result.chosen_branch_name} (weight: {dfel.root_result.chosen_weight:.3f})")
print(f"CP3 heterogeneity p-value: {dfel.comp_test.p_value_cp3:.2e}")
print(f"Significant dFEL sites: {sum(1 for r in dfel_res if r.p_value <= 0.1)}")

# -------------------------------------------------------------
# 12. OHMIC: Electrical Resistor Network & Wheatstone Pruning
# -------------------------------------------------------------
ohmic_aln = hp.Alignment.load("benchmarks/data/cd2.fna")
ohmic_settings = hp.OhmicAnalysisSettings()
ohmic_settings.use_codon_model = True
ohmic_settings.network_settings.use_hyphaeon_transformer = True

ohmic = hp.OhmicAnalyzer(ohmic_aln, ohmic_settings)
ohmic_res = ohmic.run()

print(f"Foster sum invariant: {ohmic_res.foster_sum:.5f} (N - 1 = {ohmic_res.num_taxa - 1})")
print(f"Active 4-way Wheatstone bridges: {ohmic_res.num_bridges}, Delta LogL: {ohmic_res.delta_log_l:+.4f}")

# Plug the Wheatstone-bridge-equipped tree directly into any standard analysis (e.g. FEL, FUBAR)
wheatstone_tree = ohmic_res.network.collapsed_tree
fel_on_ohmic = hp.FELAnalyzer.create_and_fit(wheatstone_tree, ohmic_aln)
fel_ohmic_sites = fel_on_ohmic.run()
```

### PyTorch End-to-End Optimization

```python
import torch
import hyphy3 as hp

class TreeLikelihoodFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, branch_lengths, tree, aln, params):
        bls = branch_lengths.detach().cpu().numpy().tolist()
        log_l, grads = hp.compute_branch_length_gradients(tree, aln, params, bls)
        ctx.save_for_backward(torch.tensor(grads, dtype=torch.float64))
        return torch.tensor(log_l, dtype=torch.float64)

    @staticmethod
    def backward(ctx, grad_output):
        grads, = ctx.saved_tensors
        return grad_output * grads, None, None, None

# Optimize branch lengths with PyTorch Adam
aln = hp.Alignment.load("benchmarks/data/cd2.fna")
tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")
params = hp.MG94Parameters()

bls = torch.tensor([node.branch_length for node in tree.nodes], requires_grad=True)
optimizer = torch.optim.Adam([bls], lr=0.01)

for step in range(50):
    optimizer.zero_grad()
    loss = -TreeLikelihoodFunction.apply(bls, tree, aln, params)
    loss.backward()
    bls.grad[tree.root_id] = 0.0
    optimizer.step()
```

---

## 🎨 Interactive Visualizations & Standalone HTML Reports

HyPhy 3 features a modern, interactive visualization engine designed according to **Dieter Rams' Braun Industrial Design System** (*"Less, but better"* — *Weniger, aber besser*). It eliminates clutter, presenting dense phylogenetic and selection data with clarity, tactile elegance, and instant responsiveness.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│  HYPHY 3 • FEL VIEWER                CD2.FNA • 10 TAXA • 187 SITES            [EXPORT] │
├────────────────────────────────────────────────────────────────────────────────────────┤
│  [BIRD'S-EYE OVERVIEW TRACK]                                                           │
│  LRT ▅▅▅▅▆▇█▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅▅  │
│  [  BRUSH / ZOOM WINDOW  ]                                                             │
├────────────────────────────────┬───────────────────────────────────────┬───────────────┤
│  PHYLOGENETIC TREE CANVAS      │ MULTIPLE SEQUENCE ALIGNMENT (CANVAS)  │ SITE INSPECTOR│
│  ├─ Human                      │ K E I T N A L E T W G A L G Q D I N L │ Site: 34      │
│  ├─ Chimp                      │ K E I T N A L E T W G A L G Q D I N L │ α (dS): 0.482 │
│  ├─ Gorilla                    │ K E I T N A L E T W G A L G Q D I N L │ β (dN): 2.132 │
│  ├─ Orangutan                  │ K E I T N A L E T W G A L G Q D I N L │ ω: 4.42       │
│  ├─ Macaque                    │ K E I T N A L E T W G A L G Q D I N L │ LRT: 8.64     │
│  └─ Baboon                     │ K E I T N A L E T W G A L G Q D I N L │ p-value: 0.003│
│   ▲ Synchronized Y-Scrolling ▲ │ 60 FPS Virtualized Codon Viewport     │ Diversifying  │
└────────────────────────────────┴───────────────────────────────────────┴───────────────┘
```

### 1. Visualization Architecture (`viz/fel-viewer`)

- **Synchronized Dual-Axis Exploration**:
  - **Vertical Synchronization**: The phylogenetic tree leaf nodes and alignment sequence rows are locked in a strict 1:1 vertical mapping, scrolling and zooming together synchronously.
  - **Horizontal Navigation**: Smooth multi-scale panning across thousands of codon positions without lag.
- **Virtualized High-Performance HTML5 Canvas**:
  - Employs off-screen virtualized canvas rendering, delivering fluid 60 FPS scrolling on large sequence alignments containing hundreds of taxa and thousands of sites.
  - Supports amino acid translation, codon triplets, biochemical property coloring schemes (Hydrophobicity, Charge, Zappo, Taylor), and reference-difference highlighting.
- **Bird's-Eye Panoramic Overview Track**:
  - Full-gene selection landscape profile plotting test statistics across all sites (FEL/Contrast-FEL LRT, MEME episodic p-values, FUBAR posterior probabilities).
  - Interactive brush-to-zoom scrubber allowing instant macro-to-micro transitions from gene overview to individual residue resolution.
- **Real-Time Site Inspector**:
  - Real-time side panel presenting complete site-level inference statistics upon click or hover:
    - Point estimates of synonymous (α) and non-synonymous (β) substitution rates.
    - Likelihood ratio test statistics, p-values, and Benjamini-Hochberg q-values.
    - Group-specific rate partitions and shrinkage parameters (Contrast-FEL).
    - Branch attribution with Empirical Bayes Factor (EBF) rankings (MEME).
    - Bivariate posterior mass distribution across the 20 × 20 rate grid (FUBAR).

### 2. Standalone Single-File HTML Reports (`hyphy3 --html`, `hyphy3 report`)

The viewer is compiled into the `hyphy3` binary, so reports need nothing beyond `hyphy3` itself: no Node, Python, or files shipped alongside the binary. Each report is **one self-contained HTML file**:

```bash
# As part of an analysis run
hyphy3 meme --alignment data.nex --output results.MEME.json --html results.MEME.html

# From an existing result file (hyphy3 v3 or legacy HyPhy 2.5 / Datamonkey JSON)
hyphy3 report results.FEL.json --html results.FEL.html
```

Reports are currently available for **FEL** and **MEME**; other analyses are refused with a clear message until their viewers are added.

- **How the viewer gets embedded**: at build time, `cmake/EmbedViewer.cmake` bundles the ES modules under `viz/fel-viewer/src` (following their imports from `index.js`) together with the stylesheet into a generated header. Editing the viewer and re-running the build is all that's needed.
- **Developing the viewer**: the sources are plain browser ES modules with no build step. Serve `viz/fel-viewer/` with any static file server and open `index.html`.

- **Zero Network Dependencies**: Inlines the complete JavaScript runtime, CSS styling, and JSON dataset directly into the HTML document.
- **Universal Portability**: Double-click to open in Safari, Chrome, Firefox, or Edge. Works directly via `file://` URLs, in air-gapped laboratory environments, HPC clusters, and email attachments.
- **Publication Ready**: Includes one-click SVG/PNG export for publication-quality figures.

### 3. Native macOS Desktop Application (`mac/HyPhyStudio`)

For macOS users, HyPhy 3 includes **HyPhyStudio** — a native desktop application built with Swift 6 and SwiftUI:
- **Swift-C++ Interoperability**: Direct zero-overhead bridging to the `hyphy::core` and `hyphy::analyses` C++20 engine (`.interoperabilityMode(.Cxx)`).
- **Embedded WebKit Viewer**: Hosts the interactive visualization engine inside an optimized Cocoa `WKWebView`, providing native file drop, menu bar integration, and multi-core computation.

---

## 📂 Repository Structure

```
hyphy-3/
├── CMakeLists.txt              # Modern C++20 build configuration
├── README.md                   # Documentation for human developers & users
├── AGENT.md                    # In-depth architectural guide for AI agents & contributors
├── .gitignore
├── include/
│   └── hyphy/
│       ├── core/               # Alignment, Tree, GeneticCode, RateMatrix, LikelihoodEngine
│       ├── opt/                # Brent, NelderMead, SQUAREM, L-BFGS, Adam
│       ├── autograd/           # Dynamic computational graph, Var, TreeLikelihoodNode
│       └── analyses/           # GTR, MG94, FEL, Contrast-FEL, MEME, BUSTED (SRV & MH), aBSREL, RELAX, FUBAR
├── src/
│   ├── apps/                   # CLI drivers: hyphy3, hyphy_fel, hyphy_contrast_fel, hyphy_fubar, hyphy_prime, hyphy_meme, hyphy_busted, hyphy_absrel, hyphy_relax
│   └── python/                 # Nanobind C++ Python bridge (bindings.cpp)
├── python/
│   └── hyphy3/                 # Pure Python high-level API & wrappers
├── viz/                        # Interactive browser visualizers (embedded into hyphy3 for HTML reports)
│   └── fel-viewer/             # Braun-inspired FEL / MEME results viewer (plain ES modules)
├── cmake/                      # Build helpers (EmbedViewer.cmake bundles the viewer into the binary)
├── mac/                        # Native macOS applications
│   └── HyPhyStudio/            # Swift 6 / SwiftUI macOS desktop studio
├── examples/
│   ├── cli/                    # Shell scripts: run_fel.sh, run_meme.sh, run_busted.sh, run_absrel.sh
│   └── python/                 # Python examples (FEL, MEME, BUSTED, aBSREL, Autograd, PyTorch)
├── benchmarks/
│   └── data/                   # Standard benchmark alignments and trees (CD2, ADH, bglobin, etc.)
└── tests/                      # Doctest unit & parity test suites
```

---

## 📖 Citation

If you use HyPhy 3 in your research, please cite:

- **HyPhy 3 Modern Core**: Kosakovsky Pond SL, et al. *HyPhy 3: A Modern Differentiable Engine for Molecular Evolution*. (In preparation).
- **PRIME**: Kosakovsky Pond SL, et al. *PRoperty Informed Models of Evolution*. (In preparation).
- **Contrast-FEL**: Kosakovsky Pond SL, et al. *Contrast-FEL: A test for differences in selective pressures at individual sites among clades and sets of branches*. Mol Biol Evol. 38(3):1184–1198 (2021).
- **FUBAR**: Murrell B, et al. *FUBAR: a fast, unconstrained bayesian approximation for inferring selection*. Mol Biol Evol. 30(5):1196–1205 (2013).
- **aBSREL**: Smith MD, et al. *Less Is More: An Adaptive Branch-Site Random Effects Model for Efficient Detection of Episodic Diversifying Selection*. Mol Biol Evol. 32(5):1342–1353 (2015).
- **BUSTED**: Murrell B, et al. *Gene-wide identification of episodic selection*. Mol Biol Evol. 32(5):1365–1371 (2015).
- **BUSTED-S**: Wisotsky SR, et al. *Synonymous rate variation improves the detection of positive selection*. Mol Biol Evol. 37(8):2430–2439 (2020).
- **BUSTED-MH**: Lucaci AG, et al. *Evolutionary models considering multiple nucleotide substitutions improve detection of positive selection*. Mol Biol Evol. 38(7):3081–3096 (2021).
- **MEME**: Murrell B, et al. *Detecting episodic selection with a mixed effects model of evolution*. PLoS Genet. 8(7):e1002764 (2012).
- **FEL**: Kosakovsky Pond SL & Frost SDW. *Not so different after all: a comparison of methods for detecting amino acid sites under selection*. Mol Biol Evol. 22(5):1208–1222 (2005).

---

## 📄 License

HyPhy 3 is open-source software licensed under the [MIT License](LICENSE).
Developed by Sergei L. Kosakovsky Pond and the HyPhy development team.
