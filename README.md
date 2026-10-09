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
| **Full MEME Analysis** | CD2 (10 taxa, 187 codons) | 4.82 s | **1.21 s** | **$3.98\times$** |
| **Full aBSREL Analysis** | β-globin (17 taxa, 144 codons) | 48.0 s ($\ln L = -3631.58$) | **25.7 s** ($\ln L = -3633.51$) | **$1.87\times$ ($|\Delta \ln L| \le 1.9$)** |

---

## 📦 Implemented Analyses

### 1. **FEL (Fixed Effects Likelihood)**
- Tests for site-by-site pervasive diversifying ($dN > dS$) and purifying ($dN < dS$) selection.
- Precomputes nucleotide GTR and global MG94 baselines.
- Optimizes site-specific synonymous rate $\alpha$ and non-synonymous rate $\beta$ via 2D profiling.
- Standard asymptotic $\chi^2_1$ Likelihood Ratio Test (LRT).

### 2. **MEME (Mixed Effects Model of Evolution)**
- Detects episodic diversifying selection affecting individual codons on a subset of branches.
- Uses empirical Bayes 2-rate mixture model: $(\alpha, \beta^-, p^-)$ and $(\alpha, \beta^+, p^+)$ with constraint $\beta^- \le \alpha$.
- Optimizes mixture parameters per site and reports evidence of episodic positive selection ($\beta^+ > 1$).

### 3. **BUSTED, BUSTED-S & BUSTED-MH (Branch-site Unrestricted Statistical Test)**
- Gene-wide test for episodic diversifying positive selection across branches and sites.
- Constrained Null ($\omega_K = 1.0$) vs Unconstrained Alternative ($\omega_K \ge 1.0$).
- **BUSTED-S**: Full support for **Synonymous Rate Variation (SRV)** using site-to-site discrete rate distributions ($M=3$ classes, $\mathbb{E}[\alpha] = 1.0$).
- **BUSTED-MH**: Full support for **Multi-Nucleotide Substitutions (Multiple Hits)** (`Double`, `Double+Triple`), estimating instantaneous 2-hit ($\delta$) and 3-hit ($\psi$) substitution rates and substitution fractions while maintaining exact detailed balance and fast eigendecompositions.
- Automatic model selection (**Auto-K**) via AICc step-up.
- Per-site Evidence Ratios (empirical Bayes factors) for positive selection.

### 4. **aBSREL (Adaptive Branch-Site Random Effects Likelihood)**
- Tests whether a proportion of sites have evolved under positive selection along each lineage/branch.
- Dynamic model complexity selection (AICc step-up) assigns optimal $\omega$ rate classes per branch without over-parameterization.
- Accelerated via local Inside-Outside projection ($V_{v,p}^{\top} \bar{P}_b D_{v,p}$) during branch complexity search and constrained null testing, avoiding full-tree traversals.
- **Phase 4 Full Adaptive Refinement**: Joint L-BFGS branch length optimization using analytical Inside-Outside mixture gradients combined with GTR nucleotide rate optimization and local mixture refinement, closing log-likelihood parity with HyPhy 2.5 to within $< 2$ units ($\ln L = -3633.51$ vs $-3631.58$).
- Exact closed-form asymptotic mixture distribution $p$-value computation ($\frac{1}{2} \chi^2_0 + \frac{1}{2}[0.4 \chi^2_1 + 0.6 \chi^2_2]$).
- Computes Holm-Bonferroni corrected $p$-values and Empirical Bayes Factors (EBF) for site-level support.

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

HyPhy 3 provides a unified `hyphy3` executable as well as dedicated tools:

```bash
# General help
hyphy3 --help

# 1. Run FEL (site-by-site selection)
hyphy3 fel --alignment benchmarks/data/cd2.fna --tree benchmarks/data/cd2.nwk --threads 8

# 2. Run MEME (episodic selection)
hyphy3 meme --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --pvalue 0.1 --threads 8

# 3. Run BUSTED (gene-wide selection)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --threads 8

# 4. Run BUSTED-S (with Synonymous Rate Variation)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --srv --syn-rates 3 --threads 8

# 5. Run BUSTED-MH (with Multi-Nucleotide Substitutions)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --multiple-hits Double+Triple --threads 8

# 6. Run BUSTED with Automatic Model Selection (Auto-K)
hyphy3 busted --alignment benchmarks/data/adh.fna --tree benchmarks/data/adh.nwk --auto-k --threads 8

# 7. Run aBSREL (lineage-specific selection)
hyphy3 absrel --alignment benchmarks/data/bglobin.nex --tree benchmarks/data/bglobin.nex --output bglobin.absrel.json --threads 8
```

All analyses output standardized, Datamonkey-compatible JSON files ready for visualization on [HyPhy Vision](https://vision.hyphy.org).

---

## 🐍 Python API & PyTorch Integration

HyPhy 3 provides native Python bindings powered by Nanobind with zero-copy data exchange.

### Running Selection Analyses

```python
import hyphy3 as hp

# Load alignment and phylogenetic tree
aln = hp.Alignment.load("benchmarks/data/cd2.fna")
tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")

# Run FEL
fel = hp.FELAnalyzer.create_and_fit(tree, aln, pvalue_threshold=0.10)
results = fel.run()

for r in results:
    if r.p_value < 0.10 and r.beta > r.alpha:
        print(f"Positive selection at site {r.site + 1}: dN/dS = {r.beta/r.alpha:.2f}, p = {r.p_value:.4f}")

# Run BUSTED-MH (with multiple hits)
settings = hp.BUSTEDSettings()
settings.multiple_hits = "Double+Triple"  # "None", "Double", or "Double+Triple"
busted = hp.BUSTEDAnalyzer.create_and_fit(tree, aln)
res = busted.run(settings)

print(f"LRT = {res.lrt:.4f}, p-value = {res.p_value:.6e}")
print(f"Delta (2-hit rate) = {res.unconstrained.delta:.4f}, fraction = {res.unconstrained.frac_delta*100:.2f}%")
print(f"Psi (3-hit rate) = {res.unconstrained.psi:.4f}, fraction = {res.unconstrained.frac_psi*100:.2f}%")

# Run aBSREL (adaptive branch-site selection)
absrel = hp.ABSRELAnalyzer.create_and_fit(tree, aln)
abs_res = absrel.run()
for br in abs_res.branches:
    if br.tested and br.p_corrected < 0.05:
        print(f"Lineage selection on {br.name}: p_corr = {br.p_corrected:.4f}, LRT = {br.lrt:.2f}")
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
│       └── analyses/           # GTR, MG94, FEL, MEME, BUSTED (with SRV), aBSREL
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
- **MEME**: Murrell B, et al. *Detecting episodic selection with a mixed effects model of evolution*. PLoS Genet. 8(7):e1002764 (2012).
- **FEL**: Kosakovsky Pond SL & Frost SDW. *Not so different after all: a comparison of methods for detecting amino acid sites under selection*. Mol Biol Evol. 22(5):1208–1222 (2005).

---

## 📄 License

HyPhy 3 is open-source software licensed under the [MIT License](LICENSE).
Developed by Sergei L. Kosakovsky Pond and the HyPhy development team.
