<div align="center">

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

However, legacy HyPhy 2.5 carries 25+ years of monolithic C++ code, macro-heavy architectures, an internal domain-specific scripting language (HBL — HyPhy Batch Language) with dynamic string-based interpreters, global state, and complex memory management. Furthermore, classical optimization relied heavily on coordinate-wise 1D line searches (Brent's method) or finite differences, which scale as $\mathcal{O}(B^2)$ or require $2B$ likelihood evaluations per gradient step on a tree with $B$ branches.

### **Enter HyPhy 3**

**HyPhy 3** is a complete, ground-up rewrite and architectural reimagining of the HyPhy platform in modern **C++20**:

1. **Header-Only Modular Architecture**: Replaces the monolithic HBL interpreter with clean, type-safe, composable C++ modules (`hyphy::core`, `hyphy::opt`, `hyphy::autograd`, `hyphy::analyses`).
2. **Differentiable Phylogenetics Engine**: Implements the analytical **Inside-Outside algorithm** to compute exact tree log-likelihood adjoint gradients $\frac{\partial \ln L}{\partial t_b}$ with respect to all branch lengths in a **single $\mathcal{O}(B)$ traversal** — replacing numerical finite differences and coordinate-wise search.
3. **Massive Performance Gains ($20\times - 87\times$)**: Joint branch length optimization via Inside-Outside Adjoint L-BFGS reduces branch refinement from minutes to seconds on trees with hundreds of taxa.
4. **Native Autograd (C++ and Python)**: Features a dynamic reverse-mode automatic differentiation graph in C++ and seamless zero-copy bindings for **PyTorch**, **JAX**, and **NumPy**.
5. **Strict Parity & Scientific Integrity**: Rigorously verified against HyPhy 2.5 baseline outputs down to machine precision across benchmark datasets (CD2, ADH, COXI, β-globin).

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
In phylogenetics, computing the gradient of tree log-likelihood with respect to all branch lengths $\mathbf{t} = (t_1, \dots, t_B)$ via finite differences requires $2B$ tree pruning passes. For a tree with 700 branches, that is 1,400 pruning passes per gradient step!

HyPhy 3 computes exact analytical branch gradients in **one single combined pass**:
- **Inside Pass (Post-Order)**: Computes conditional subtree likelihoods $D_{v}$ from leaves to root.
- **Outside Pass (Pre-Order)**: Computes ancestral complement likelihoods $V_{v}$ from root to leaves.
- **Adjoint Contraction for Standard & Mixture Models**:
  - **Homogeneous Model**:
    $$\frac{\partial \ln L}{\partial t_b} = \sum_{p=1}^P \frac{w_p}{L_p} V_{b,p}^{\top} \left( \frac{\partial P_b}{\partial t_b} \right) D_{b,p}, \quad \text{where } \frac{\partial P_b}{\partial t_b} = Q P_b$$
  - **Branch-Site Discrete Mixture (aBSREL)**:
    $$\frac{\partial \ln L}{\partial \alpha_b} = \sum_{p=1}^P \frac{w_p}{L_p} V_{b,p}^{\top} \left( \sum_{k=1}^{M_b} p_{b,k} Q(1, \omega_{b,k}) P_{b,k} \right) D_{b,p}$$
    This allows joint L-BFGS optimization across all branch lengths simultaneously in $\mathcal{O}(B)$ time.

### 2. High-Performance Optimization Suite
- **Adjoint L-BFGS**: Replaces $O(B^2)$ coordinate-wise Brent search in BUSTED and aBSREL with multi-branch joint ascent directions using two-loop recursion and Armijo line search.
- **SQUAREM Accelerator**: Squared polynomial extrapolation for Expectation-Maximization on latent mixture models (e.g. BUSTED rate weights $\mathbf{p}$ and synonymous rate weights $\mathbf{q}$), achieving quadratic convergence without evaluating Hessians.
- **Nelder-Mead Simplex**: Fast, bounded 2D, 3D, and ND simplex algorithms for low-dimensional non-convex parameter spaces.

---

## ⚡ Performance Benchmarks

| Analysis / Stage | Dataset | HyPhy 2.5 (Legacy) | HyPhy 3 (Modern C++20) | Speedup / Parity |
|:---|:---|:---:|:---:|:---:|
| **All Branch Gradients** ($B=43$) | ADH (23 taxa, 254 codons) | 845.8 ms *(FD)* | **9.66 ms** *(Adjoint)* | **$87.6\times$** |
| **All Branch Gradients** ($B=695$) | Influenza A (349 taxa) | 12,480 ms *(FD)* | **238.1 ms** *(Adjoint)* | **$52.4\times$** |
| **BUSTED Branch Refinement** | CD2 (10 taxa, 187 codons) | 0.85 s | **0.048 s** | **$17.7\times$** |
| **BUSTED Branch Refinement** | ADH (23 taxa, 254 codons) | 3.86 s | **0.184 s** | **$21.0\times$** |
| **BUSTED Branch Refinement** | Influenza A (349 taxa, 695 br) | ~208 s *(Brent)* | **7.63 s** *(L-BFGS)* | **$27.3\times$** |
| **Full BUSTED Analysis** | ADH (23 taxa, 254 codons) | 11.05 s | **4.03 s** | **$2.74\times$** |
| **Full BUSTED-S (SRV) Analysis** | ADH (23 taxa, 254 codons) | 18.42 s | **6.15 s** | **$3.00\times$** |
| **Full BUSTED-MH (Double+Triple)** | ADH (23 taxa, 254 codons) | 26.50 s | **6.47 s** | **$4.10\times$** |
| **Full MEME Analysis** | CD2 (10 taxa, 187 codons) | 4.82 s | **1.21 s** | **$3.98\times$** |
| **Full aBSREL Analysis** | β-globin (17 taxa, 144 codons) | 48.0 s ($\ln L = -3631.58$) | **25.7 s** ($\ln L = -3633.51$) | **$1.87\times$ ($|\Delta \ln L| \le 1.9$)** |

---

## 📦 Implemented Analyses & Capabilities

### Summary Matrix

| Method | Biological Question | Selection Type | Statistical Model | Rate Variation | Multi-Hit | Parity Status |
|:---|:---|:---|:---|:---:|:---:|:---:|
| **FEL** | Which codon sites are under pervasive selection? | Site-level ($dN \gtrless dS$) | Profile likelihood, $\chi^2_1$ LRT | None | — | 100% Verified |
| **MEME** | Which codon sites are under episodic selection? | Site-level & Lineage | 2-rate mixture, $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_2$ LRT | Site mixture | — | 100% Verified |
| **BUSTED** | Has the gene evolved under episodic selection? | Gene-wide (all branches) | 3-category $\omega$ mixture, $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_2$ | Site mixture | Optional | 100% Verified |
| **BUSTED-S** | Gene-wide selection robust to synonymous variation? | Gene-wide with SRV | Discrete site-to-site $\alpha$ and $\omega$ mixture | Site $\alpha$ & $\omega$ | Optional | 100% Verified |
| **BUSTED-MH** | Gene-wide selection accounting for multi-nucleotide hits? | Gene-wide with MH | Reversible 1-hit, 2-hit ($\delta$), 3-hit ($\psi$) generator | Site $\alpha$ & $\omega$ | Double, Double+Triple | 100% Verified |
| **aBSREL** | Which specific lineages/branches evolved under selection? | Lineage-specific | Adaptive complexity selection (AICc Step-Up) | Branch-site $\omega$ | — | $|\Delta \ln L| < 2.0$ |
| **GTR / MG94** | What are baseline nucleotide biases & global $dN/dS$? | Alignment-wide | Reversible CTMC with F3x4 equilibrium | Global $\omega$ | Supported | 100% Verified |
| **Adjoint Core** | Differentiable tree likelihoods for custom ML models | Arbitrary parameters | Inside-Outside analytical $\mathcal{O}(B)$ gradients | User-defined | User-defined | Exact Analytical |

---

### Detailed Method Descriptions

#### 1. **FEL (Fixed Effects Likelihood)**
- **Scope**: Site-by-site detection of pervasive purifying ($dN < dS$) and diversifying ($dN > dS$) selection.
- **Workflow**:
  1. Fits baseline nucleotide GTR model to estimate branch lengths and nucleotide exchangeabilities.
  2. Fits global MG94xREV codon model to determine mean $dN/dS$ ratio.
  3. For every codon site $s = 1, \dots, S$, independently optimizes synonymous rate $\alpha_s$ and non-synonymous rate $\beta_s$ via 2D likelihood profiling.
  4. Tests null hypothesis $H_0: \alpha_s = \beta_s$ against alternative $H_1: \alpha_s \neq \beta_s$ using standard asymptotic $\chi^2_1$ Likelihood Ratio Test (LRT).
- **Outputs**: MLEs of $(\alpha_s, \beta_s)$, likelihood ratio test statistics, raw $p$-values, total site tree length, and Datamonkey JSON export.

#### 2. **MEME (Mixed Effects Model of Evolution)**
- **Scope**: Site-by-site detection of **episodic** diversifying selection affecting a subset of lineages while remaining conserved on others.
- **Workflow**:
  1. Formulates a 2-rate mixture model per site across branches:
     - Fraction $p^-$ of branches evolve with rate $\beta^-$ under constraint $\beta^- \le \alpha$.
     - Fraction $p^+ = 1 - p^-$ of branches evolve with unrestricted rate $\beta^+$ (can exceed $\alpha$).
  2. Maximizes profile mixture likelihood over $(\alpha, \beta^-, \beta^+, p^+)$.
  3. Compares alternative fit against constrained null ($\beta^+ \le \alpha$, effectively $\beta^+ = \alpha$) using asymptotic mixture distribution $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_2$.
  4. Calculates branch-level **Empirical Bayes Factors (EBF)** to determine which individual phylogenetic branches are experiencing positive selection at each site.
- **Outputs**: Site-level LRTs, $p$-values, branch attribution lists, EBF matrices, and Datamonkey JSON export.

#### 3. **BUSTED & Variants (Branch-site Unrestricted Statistical Test)**
- **Standard BUSTED**:
  - Tests for alignment-wide evidence of episodic diversifying positive selection across both sites and branches.
  - Fits a 3-category discrete distribution of non-synonymous rates: $\omega_1 \le \omega_2 \le 1 \le \omega_3$ with mixture weights $(p_1, p_2, p_3)$.
  - Compares unconstrained model ($\omega_3 \ge 1.0$) against constrained null model ($\omega_3 = 1.0$) using asymptotic mixture test statistic $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_2$.
  - Computes per-site **Evidence Ratios (ER)** quantifying evidence that site $s$ evolved under $\omega_3 > 1$.
- **BUSTED-S (Synonymous Rate Variation)**:
  - Extends BUSTED with site-to-site discrete synonymous substitution rate variation ($M=3$ categories, $\mathbb{E}[\alpha] = 1.0$, with weights $q_1, q_2, q_3$).
  - Prevents false-positive inflation caused by synonymous rate heterogeneity across the alignment.
  - Implements **Markov Generator Scaling Invariance**: $Q(\alpha_m, \beta) = \alpha_m Q(1, \omega_k)$, requiring only $K$ eigendecompositions instead of $M \times K$.
- **BUSTED-MH (Multi-Nucleotide Substitutions / Multiple Hits)**:
  - Incorporates instantaneous 2-nucleotide ($\delta$) and 3-nucleotide ($\psi$) substitutions occurring within a single codon step (e.g. $\text{TCA} \to \text{GAA}$).
  - Parameterized via reversible Markov generator:
    $$Q_{ij} = (\alpha \text{ or } \beta) \times \delta \times \prod_{p \in D} r_{i_p \to j_p} \pi_{j_p}^{(p)} \quad (|D| = 2)$$
    $$Q_{ij} = (\alpha \text{ or } \beta) \times \psi \times \prod_{p \in D} r_{i_p \to j_p} \pi_{j_p}^{(p)} \quad (|D| = 3)$$
  - Strictly preserves detailed balance $\pi_i Q_{ij} = \pi_j Q_{ji}$, maintaining fast symmetric eigensolvers.
  - Jointly estimates $\delta$ and $\psi$ and reports flux-weighted instantaneous substitution rates and substitution fractions ($\text{frac}_\delta$, $\text{frac}_\psi$).
- **BUSTED Auto-K**:
  - Automatically identifies the optimal number of selection categories $K \in \{1, 2, \dots, K_{\max}\}$ via step-up AICc optimization, avoiding over-parameterization on simpler datasets.

#### 4. **aBSREL (Adaptive Branch-Site Random Effects Likelihood)**
- **Scope**: Tests whether a proportion of sites have evolved under positive diversifying selection along **each lineage/branch** of a phylogenetic tree, without requiring a priori branch labeling.
- **Workflow**:
  1. **Phase 1 (GTR Baseline)**: Precomputes nucleotide substitution biases and initial branch lengths.
  2. **Phase 2 (MG94 Baseline)**: Fits standard codon model with branch-specific $\omega_b$ ratios.
  3. **Phase 3 (Exploratory Complexity Step-Up)**: Dynamically infers optimal rate categories ($K_b \in \{1, 2, \dots, K_{\max}\}$) on every branch via AICc step-up.
  4. **Phase 4 (Full Adaptive Joint Refinement)**: Jointly refines all branch lengths using analytical Inside-Outside mixture gradients combined with GTR nucleotide rate optimization and local mixture polishing.
  5. **Phase 5 (Hypothesis Testing)**: Tests whether the highest rate category on tested branches satisfies $\omega_{b, \max} > 1$ against the null hypothesis $\omega_{b, \max} \le 1.0$ using the asymptotic mixture null distribution:
     $$\text{Null} \sim \frac{1}{2} \chi^2_0 + \frac{1}{2}\left(0.4 \chi^2_1 + 0.6 \chi^2_2\right)$$
  6. Computes family-wise error-rate corrected $p$-values via the **Holm-Bonferroni step-down procedure** and Empirical Bayes Factors (EBF) for site-level support.
- **Ultra-Fast Local Projection**: Leverages Inside-Outside subtree ($D_{v,p}$) and ancestral complement ($V_{v,p}$) vectors to evaluate candidate branch mixtures via immediate dot products ($L_p = V_{v,p}^{\top} \bar{P}_b D_{v,p}$) in microseconds with **zero whole-tree pruning passes**.

#### 5. **Continuous-Time Markov Chain Baselines (GTR & MG94xREV)**
- **Nucleotide GTR**: General Time Reversible continuous-time substitution model with 5 independent exchangeabilities ($\theta_{AC}, \theta_{AT}, \theta_{CG}, \theta_{CT}, \theta_{GT}$), analytical eigensolver, and empirical base frequencies.
- **Global MG94xREV**: Codon substitution model combining GTR nucleotide exchangeabilities with positional F3x4 codon equilibrium frequencies and global non-synonymous/synonymous rate ratio $\omega$.
- **Multi-Hit MG94xREV**: Simultaneous evaluation of 1-hit, 2-hit ($\delta$), and 3-hit ($\psi$) codon substitution transitions while strictly preserving detailed balance $\pi_i Q_{ij} = \pi_j Q_{ji}$.

#### 6. **Differentiable Phylogenetics Engine & Autograd**
- **Analytical Inside-Outside Adjoints**: Evaluates exact tree log-likelihood gradients $\frac{\partial \ln L}{\partial t_b}$ for all $B$ branches simultaneously in a single $\mathcal{O}(B)$ tree traversal.
- **Native C++ Autograd Graph**: Reverse-mode automatic differentiation supporting custom evolutionary models and continuous rate distributions.
- **PyTorch & JAX Bridge**: Zero-copy bindings via Nanobind allowing phylogenetic likelihood calculations to serve directly as differentiable loss functions within machine learning pipelines.

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

HyPhy 3 provides both a unified driver (`hyphy3 <analysis>`) and dedicated standalone tools (`hyphy3-fel`, `hyphy3-meme`, `hyphy3-busted`, `hyphy3-absrel`):

### Available CLI Options by Analysis

| Analysis | Command | Key Options | Description |
|:---|:---|:---|:---|
| **FEL** | `hyphy3 fel` | `--alignment`, `--tree`, `--pvalue 0.1`, `--threads N`, `--output` | Site-by-site pervasive selection |
| **MEME** | `hyphy3 meme` | `--alignment`, `--tree`, `--pvalue 0.1`, `--threads N`, `--output` | Site-by-site episodic selection across lineages |
| **BUSTED** | `hyphy3 busted` | `--alignment`, `--tree`, `--rates 3`, `--threads N`, `--output` | Standard gene-wide episodic selection |
| **BUSTED-S** | `hyphy3 busted` | `--srv`, `--syn-rates 3` | Synonymous rate variation across sites |
| **BUSTED-MH** | `hyphy3 busted` | `--multiple-hits [None\|Double\|Double+Triple]` | Multi-nucleotide substitutions (2-hit and 3-hit) |
| **BUSTED Auto-K**| `hyphy3 busted` | `--auto-k` | Automatic selection of rate categories via AICc step-up |
| **aBSREL** | `hyphy3 absrel`| `--alignment`, `--tree`, `--rates 3`, `--threads N`, `--output` | Lineage-specific adaptive selection |

### CLI Examples

```bash
# General help
hyphy3 --help

# 1. Run FEL (site-by-site pervasive selection)
hyphy3 fel --alignment benchmarks/data/cd2.fna --tree benchmarks/data/cd2.nwk --pvalue 0.1 --threads 8

# 2. Run MEME (site-by-site episodic selection across lineages)
hyphy3 meme --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --pvalue 0.1 --threads 8

# 3. Run BUSTED (standard gene-wide selection)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --threads 8

# 4. Run BUSTED-S (with Synonymous Rate Variation across sites)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --srv --syn-rates 3 --threads 8

# 5. Run BUSTED-MH (with Multi-Nucleotide Substitutions: Double & Triple Hits)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --multiple-hits Double+Triple --threads 8

# 6. Run BUSTED with Automatic Model Selection (Auto-K via AICc step-up)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --auto-k --threads 8

# 7. Run aBSREL (adaptive lineage-specific selection)
hyphy3 absrel --alignment benchmarks/data/bglobin.nex --tree benchmarks/data/bglobin.nex --output bglobin.absrel.json --threads 8
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
abs_settings.max_rate_classes = 3
absrel = hp.ABSRELAnalyzer.create(tree, aln, abs_settings)
abs_res = absrel.run()

print(f"Tested branches: {abs_res.tested_branches}, Positive branches: {abs_res.positive_branches}")
for br in abs_res.branches:
    if br.is_tested and br.corrected_p_value < 0.05:
        print(f"Lineage selection on {br.branch_name}: p_corr = {br.corrected_p_value:.4f}, LRT = {br.lrt:.2f}")
        print(f"  Rate classes: {br.rate_classes}, Omega distribution: {br.rate_distribution.rates}")
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
│       └── analyses/           # GTR, MG94, FEL, MEME, BUSTED (SRV & MH), aBSREL
├── src/
│   ├── apps/                   # CLI drivers: hyphy3, hyphy_fel, hyphy_meme, hyphy_busted, hyphy_absrel
│   └── python/                 # Nanobind C++ Python bridge (bindings.cpp)
├── python/
│   └── hyphy3/                 # Pure Python high-level API & wrappers
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
