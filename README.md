# 🧬 HyPhy 3 & WALRUS: Differentiable Molecular Evolution & Unified Selection Workbench

[![WALRUS Web App](https://img.shields.io/badge/Launch_Web_App-WALRUS_Workbench-0284c7?style=for-the-badge&logo=webassembly&logoColor=white)](https://veg.github.io/hyphy-3/)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg?style=flat-square&logo=c%2B%2B)](https://en.cppreference.com/w/cpp/20)
[![WebAssembly](https://img.shields.io/badge/WASM-SIMD128-654FF0.svg?style=flat-square&logo=webassembly)](https://veg.github.io/hyphy-3/)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-3776AB.svg?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg?style=flat-square)](LICENSE)

> ### 🌐 **[Launch the Live WALRUS WebAssembly Workbench → `https://veg.github.io/hyphy-3/`](https://veg.github.io/hyphy-3/)**
> Run end-to-end codon alignment, frameshift repair, phylogenetic rooting, recombination screening, and unified selection forensics (**BUSTED**, **FEL**, **MEME**, **Contrast-FEL**, **RELAX**, **aBSREL**, **PRIME**, ancestral state reconstruction, and epistasis) **100% locally in your browser** via SIMD128 WebAssembly. Zero installation, zero server queues, and complete data privacy.

---

## ⚡ TL;DR

**HyPhy 3** is a complete ground-up C++20 re-architecture of the [HyPhy](https://hyphy.org) molecular evolution platform, built around an analytical **Inside-Outside adjoint engine** and the flagship **WALRUS** (*Wheatstone Alignment and Likelihood Rate Unification for Selection*) pipeline:

1. **🦭 WALRUS End-to-End Pipeline ([Web App](https://veg.github.io/hyphy-3/) & `hyphy3 walrus`)**: Takes raw, unaligned, or frameshifted FASTA/NEXUS sequences and executes a 4-stage checkpointed workflow—**6-frame QA & homology admission** $\rightarrow$ **Codon Pair-HMM + POA alignment with column homology posteriors ($P_{\text{hom}}$)** $\rightarrow$ **MG94×REV phylogeny, rooting & $\Delta\text{BIC}$ recombination routing** $\rightarrow$ **unified gene/branch/site/property selection & epistatic forensics**.
2. **🚀 $\mathcal{O}(B)$ Differentiable Phylogenetic Core**: Computes exact analytical log-likelihood derivatives across all $B$ branches in a single two-pass traversal (**20×–87× faster** than numerical finite differences or 1D Brent loops), vectorized with Apple Accelerate AMX, AVX2/NEON, and WASM SIMD128 $61 \times P$ batched GEMM kernels.
3. **🧰 18+ Gold-Standard Evolutionary Methods in One Binary**: Unified CLI (`hyphy3 <method>`), zero-copy Python 3 bindings (`import _hyphy3`), self-contained interactive HTML reports (`--html`), and browser WebAssembly module—all sharing the same deterministic C++20 engine.

---

## 🦭 Flagship Component: WALRUS

**WALRUS** (*Wheatstone Alignment and Likelihood Rate Unification for Selection*) replaces fragile multi-tool bioinformatics shell pipelines with a single likelihood-consistent engine where alignment uncertainty ($P_{\text{hom}}$), frameshift repairs, tree topology, and recombination breakpoints propagate directly into downstream selection tests.

### The 4 Interactive Stages

| Stage | Name | What It Does | Key Outputs & Safeguards |
|:---|:---|:---|:---|
| **Stage 0** | **Sequence QA & Homology Admission** | Screens raw sequences across all 6 reading frames against a 3-mer amino acid seed profile; auto-detects genetic code, rescues reverse-complemented reads, trims 5'/3' UTR overhangs, and quarantines non-homologous outliers. | Per-sequence homology Z-scores, frame/strand rescue badges, internal stop-codon quarantine, and interactive admission checkboxes. |
| **Stage 1** | **Codon Pair-HMM + POA Alignment** | Aligns admitted sequences via frame-preserving Codon Partial Order Alignment (POA), computes Forward–Backward per-column homology posteriors $P_{\text{hom}}$, repairs 1-bp/2-bp frameshift indels via codon DP, and evaluates competing gap placements under MG94 likelihood. | Synchronized Codon/AA canvas, per-site $P_{\text{hom}}$ confidence track (down-weighting ambiguous columns in Stage 3), frameshift repair markers, and alternative indel path drawer. |
| **Stage 2** | **Phylogeny, Rooting & Recombination** | Fits full $61 \times 61$ MG94×REV codon branch lengths via Inside-Outside L-BFGS + coupled SPR search, audits GC3 compositional drift, and screens for recombination via sliding-window topological $\Delta\text{BIC}$. | Interactive radial/rectangular phylogram, click-to-tag `{Foreground}` clade selector, outgroup/midpoint rooting, and automatic partition routing when $\Delta\text{BIC} > 10$. |
| **Stage 3** | **Unified Selection & Evolutionary Forensics** | Jointly executes gene-wide, lineage-specific, site-level, biophysical property, and co-evolutionary tests on the curated alignment and tree(s), weighted by Stage 1 $P_{\text{hom}}$. | Unified topline verdict, synchronized site selection tracks (**FEL**, **MEME**, **Contrast-FEL**, **PRIME**), **BUSTED** / **RELAX** / **aBSREL** branch cards, Empirical Bayes ancestral codon trees, and hypergeometric epistasis network. |

### Running WALRUS

- **🌐 In the Browser (Zero Install):** Open **[https://veg.github.io/hyphy-3/](https://veg.github.io/hyphy-3/)**, drop any FASTA/NEXUS file (or click a built-in showcase: *CD2*, *$\beta$-globin*, *ADH*, *IFITM3*), and step through Stages 0–3 interactively or click **Run All Stages**. Includes a built-in pre-flight memory/runtime estimator for large alignments and one-click standalone HTML dossier export.
- **💻 From the Command Line:**
  ```bash
  # End-to-end automated run from raw FASTA to JSON + self-contained HTML dossier
  hyphy3 walrus --input unaligned_seqs.fasta --output result.walrus.json --html result.walrus.html --threads 8

  # With foreground branches tagged for Contrast-FEL, RELAX, and targeted BUSTED/aBSREL
  hyphy3 walrus --input data.nex --branches "Primate_Clade" --html primate_dossier.html
  ```

---

## 🏗️ Main Architectural Components

1. **`hyphy::core` — Differentiable Inside-Outside Codon Engine** (`include/hyphy/core/`)
   - **Analytical Adjoints**: Exact post-order (Inside $\mathcal{L}_u$) and pre-order (Outside $\mathcal{U}_u$) message passing computes all branch derivatives $\partial \ln L / \partial t_b = \sum_s w_s (\mathbf{u}_u^{(s)})^T \mathbf{Q}\mathbf{P}(t_b)\mathbf{l}_u^{(s)} / L^{(s)}$ in $\mathcal{O}(B)$.
   - **Batched Matrix Kernels**: Packs all $P$ unique site patterns into contiguous $61 \times P$ matrices (`Accelerate cblas_dgemm` on macOS AMX, AVX2/NEON on Linux, SIMD128 in WebAssembly).
   - **Tree-Free Ohmic & Wheatstone Networks**: Infers exact maximum-likelihood graphical resistor networks (`--ohmic`) and 4-way reticulate Wheatstone bridges directly from FASTA alignments without requiring an input tree.
2. **`hyphy::opt` & `hyphy::autograd` — High-Order Optimization** (`include/hyphy/opt/`, `include/hyphy/autograd/`)
   - L-BFGS with analytical tree gradients, SQUAREM-accelerated EM for mixture weights, 2D Newton-Raphson with exact analytical Hessians for site models, and reverse-mode `Var` computational graphs exportable to PyTorch.
3. **`hyphy::analyses` — Unified Selection & Alignment Suite** (`include/hyphy/analyses/`)
   - All 18+ evolutionary analyses implemented as clean, composable C++20 classes emitting standardized Datamonkey-compatible JSON and interactive standalone HTML reports.
4. **Multi-Target Frontends** (`web/`, `src/apps/`, `src/python/`, `viz/`)
   - **WebAssembly Workbench** (`web/`): Pure client-side ES module + Web Worker driving `hyphy3_walrus.js`.
   - **Native CLI** (`hyphy3` and standalone binaries `hyphy3-walrus`, `hyphy3-fel`, `hyphy3-meme`, `hyphy3-busted`, etc.).
   - **Python 3 Nanobind Extension** (`_hyphy3`): Zero-copy C++/Python objects for scripting and deep-learning integration.

### Engine Performance (Apple Silicon M-Series, 8 Threads, Full Convergence)

| Analysis | Dataset ($N$ Taxa × $S$ Codons) | Legacy HyPhy 2.5 | **HyPhy 3** | Speedup | LogL Agreement |
|:---|:---|---:|---:|---:|:---|
| **FEL** | `CD2.nex` (10 × 187) | 1.85 s | **0.09 s** | **20.5×** | $\Delta \ln L < 10^{-4}$ |
| **MEME** | `CD2.nex` (10 × 187) | 3.40 s | **0.14 s** | **24.3×** | $\Delta \ln L < 10^{-4}$ |
| **BUSTED-MH + SRV** | `adh.nex` (23 × 254) | 18.20 s | **0.52 s** | **35.0×** | $\Delta \ln L < 10^{-3}$ |
| **aBSREL** | `bglobin.nex` (17 × 144) | 9.60 s | **0.31 s** | **31.0×** | $\Delta \ln L < 10^{-3}$ |
| **FUBAR (VB0)** | `CD2.nex` (10 × 187) | 1.92 s | **0.022 s** | **87.3×** | $\Delta \ln L < 10^{-4}$ |

---

## 📊 Implemented Methods & CLI Reference

Every method below can be invoked via `hyphy3 <subcommand> [options]`, from Python (`_hyphy3`), or—for core selection and alignment workflows—directly inside the **[WALRUS Web App](https://veg.github.io/hyphy-3/)**. All selection methods accept `--alignment <file>`, optional `--tree <file>` (or infer a tree/network automatically), `--threads <N>`, `--output <json>`, and `--html <report.html>`.

| Method / CLI Subcommand | Scope | Biological Question | Statistical Model & Key Features |
|:---|:---|:---|:---|
| **`walrus`** (`characterize`) | End-to-End Pipeline | *How do we go from raw, uncurated sequences to alignment-aware phylogeny and unified selection forensics?* | **4-Stage Pipeline**: 6-frame QA $\rightarrow$ Codon Pair-HMM + POA alignment ($P_{\text{hom}}$, frameshift DP repair) $\rightarrow$ MG94×REV + SPR + $\Delta\text{BIC}$ recombination $\rightarrow$ joint gene/branch/site/property selection, ancestral codons & epistasis. **[Live Web App](https://veg.github.io/hyphy-3/)** |
| **`fel`** | Site (Pervasive) | *Which individual sites experience pervasive purifying ($\beta < \alpha$) or diversifying ($\beta > \alpha$) selection?* | Site-by-site fixed effects ML $(\alpha_s, \beta_s)$ via 2D analytical Newton-Raphson + $\chi^2_1$ LRT; supports column reliability weights $P_{\text{hom}}$ and `--srv`. |
| **`contrast-fel`** | Site × Branch Group | *Which sites evolve under different selective pressures ($\beta_1 \neq \beta_2 \dots$) between pre-defined clades or phenotypes?* | Multi-group site ML with $L_2$ Ridge shrinkage (`--shrinkage`) to prevent infinite $\hat{\beta}$ on sparse branch sets; omnibus + pairwise permutation/LRT q-values. |
| **`meme`** | Site (Episodic) | *Which sites undergo episodic positive selection ($\beta^+ > \alpha$) on a subset of lineages?* | Mixed-effects branch-site mixture at each codon ($\omega^- \le 1$ vs. unrestricted $\beta^+$); reports $\chi^2$ mixture LRT and per-branch Empirical Bayes Factors (EBF). |
| **`busted`** | Gene-Wide / Foreground | *Has any site on any branch experienced episodic diversifying selection ($\omega_3 \ge 1$)?* | $K$-category random-effects $\omega$ mixture (`--srv` synonymous rate variation, `--multiple-hits Double+Triple`, `--error-sink` alignment artifact filter, and `--auto` 4-gate model selector). |
| **`busted-ph`** | Gene × Phenotype | *Is episodic positive selection specifically associated with a foreground phenotype/clade rather than background evolution?* | 3-tree-fit Union-Intersection LRT (Unconstrained vs. FG-Null $\omega_3^{\text{FG}}=1$ vs. BG-Null $\omega_3^{\text{BG}}=1$) + $\text{AIC}_c$ background neutrality check and optional `--relax` link. |
| **`absrel`** | Branch (Lineage) | *Which specific branches in the phylogeny experienced episodic diversifying selection?* | Adaptive branch-site random effects; steps up optimal $\omega$ rate classes ($K_b \in \{1,2,3\}$) per branch via $\text{AIC}_c$, followed by Holm-Bonferroni corrected LRTs and LoF purifying-core erosion Bayes Factors. |
| **`relax`** | Branch Set (Regime) | *Has natural selection been relaxed ($K < 1$, $\omega \to 1$) or intensified ($K > 1$) on `{Test}` vs. `{Reference}` branches?* | Power-law transformation $\omega_{\text{Test}} = \omega_{\text{Ref}}^K$ across a 3-class $\omega$ distribution; optimized via joint coordinate-ascent + L-BFGS cycles. |
| **`fubar`** | Site (Bayesian) | *What is the posterior probability of pervasive positive ($\beta > \alpha$) or purifying ($\beta < \alpha$) selection at each site?* | $G \times G$ ($20 \times 20$) pre-computed codon likelihood grid with Dirichlet-prior Variational Bayes (`VB0`), collapsed Gibbs, or Metropolis-Hastings sampling; reports $\mathbb{E}[\text{FP}]$ and 95% CIs. |
| **`prime`** | Site × Biophysical Property | *Which amino-acid physicochemical properties (charge, volume, polarity, hydrophobicity, $\alpha/\beta$ propensity) are conserved ($\lambda > 0$) or radically shifting ($\lambda < 0$)?* | Property-informed codon exchangeability $\omega_{ij} = \exp(-\sum_p \lambda_p d_p(i,j))$ across `5PROP` (Conchillos-Solé) or `Atchley` scales with $L_2$ Ridge stabilization and omnibus + per-property BH FDR. |
| **`dfel`** | Site (Non-Stationary) | *Are specific residues undergoing directional selection toward target amino acids along rooted lineages?* | Upfront $\chi^2$ compositional heterogeneity screen (CP12 vs. CP3), non-reversible ML root-finding, and continuous Softmax amino-acid propensity mixtures across temporal tree strata. |
| **`slac`** | Site (Counting) | *What are the site-by-site $dN$ and $dS$ counts reconstructed directly on the ancestral phylogeny?* | Joint fractional/stochastic ML ancestral codon state reconstruction (`--samples N`) with exact binomial tests for selection. |
| **`gard`** | Alignment Partitions | *Where are recombination breakpoints causing topological incongruence across the alignment?* | CHC Genetic Algorithm + $\text{AIC}_c$ over candidate breakpoint models with per-partition NJ/ML tree refitting and Kishino-Hasegawa topological incongruence tests. |
| **`mss`** | Synonymous Codon Fitness | *Do synonymous codons differ in translational efficiency or mRNA stability beyond nucleotide mutation bias?* | MG94×REV with codon-specific relative synonymous fitness parameters ($\theta_c$ or `SynREVM` exchangeabilities) and neutrality LRTs. |
| **`expcm`** | Site × DMS Preferences | *How does natural evolutionary stringency ($\beta$) scale experimentally measured Deep Mutational Scanning (DMS) amino-acid preferences?* | Experimentally informed codon models (Bloom ExpCM) fitting stringency $\beta$, $\omega$, and site-specific differential preference shifts ($\Delta \pi_{s,a}$). |
| **`ohmic`** | Reticulate Network | *Can non-tree-like evolutionary signals (recombination, hybridization, viral reassortment) be captured without a fixed tree?* | Graphical Lasso conductances + Foster's Theorem effective resistance ($R_{\text{eff}}$) + 4-way **Wheatstone bridge** cyclic pruning with exact $2 \times 2$ Schur-complement elimination (optional `--hyphaeon` 384D transformer embeddings). |
| **`joint-align`** | Codon MSA + Tree | *How can we jointly infer a frame-aware multiple codon alignment and phylogenetic tree from unaligned sequences?* | Iterative Codon Pair-HMM + POA alignment coupled to MG94×REV tree optimization and frameshift indel repair. |
| **`tree`** / **`report`** | Utilities | *Fast ML codon/nucleotide tree inference (`hyphy3 tree`) and standalone interactive HTML generation (`hyphy3 report`).* | Converts any HyPhy 3 or legacy HyPhy 2.5 JSON into a zero-dependency interactive HTML5 canvas viewer (`--html`). |

---

## 🚀 Quick Start

### 1. Zero-Install Web App (Recommended for Interactive Exploration)
👉 **[https://veg.github.io/hyphy-3/](https://veg.github.io/hyphy-3/)**
Runs the full C++20 WALRUS engine compiled to SIMD128 WebAssembly inside a Web Worker. Upload a FASTA/NEXUS file or pick a pre-loaded benchmark dataset, inspect alignment quality and frameshift repairs, tag `{Foreground}` branches on the tree, and explore interactive selection tracks, 3D-ready property heatmaps, ancestral codon trees, and epistatic pairs.

### 2. Build Native C++20 CLI from Source

```bash
git clone https://github.com/veg/hyphy-3.git
cd hyphy-3
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j 8

# Run unit & parity tests
ctest --output-on-failure

# Run WALRUS or any individual analysis
./hyphy3 walrus --input ../benchmarks/data/CD2.nex --html cd2_walrus.html --threads 8
./hyphy3 busted --alignment ../benchmarks/data/adh.nex --auto --threads 8
```

### 3. Python 3 API (`_hyphy3`)

```bash
pip install .
```

```python
import hyphy3 as hp

aln  = hp.Alignment.load("benchmarks/data/cd2.fna")
tree = hp.Tree.from_newick_file("benchmarks/data/cd2.nwk")

# Run FEL, MEME, or BUSTED directly from Python with zero-copy C++ objects
fel = hp.FELAnalyzer.create_and_fit(tree, aln, pvalue_threshold=0.10)
for r in fel.run():
    if r.p_value <= 0.10 and r.beta > r.alpha:
        print(f"Site {r.site_index + 1}: dN/dS = {r.beta / r.alpha:.2f}, p = {r.p_value:.4g}")
```

### 4. Rebuild the WebAssembly Bundle (`web/hyphy3_walrus.js`)

Requires the [Emscripten SDK (`em++`)](https://emscripten.org):

```bash
./scripts/build_wasm.sh
python3 -m http.server 8787 --directory web
# Open http://localhost:8787/
```

---

## 📂 Repository Layout

```text
hyphy-3/
├── web/                        # 🌐 Live WALRUS WebAssembly Workbench (https://veg.github.io/hyphy-3/)
├── include/hyphy/
│   ├── core/                   # Alignment, Tree, GeneticCode, MG94/GTR RateMatrix, Inside-Outside LikelihoodEngine
│   ├── opt/                    # L-BFGS, SQUAREM, 2D Newton-Raphson, Brent, Nelder-Mead
│   ├── autograd/               # Reverse-mode computational graph & TreeLikelihoodNode
│   └── analyses/               # WALRUS, BUSTED, MEME, FEL, Contrast-FEL, aBSREL, RELAX, FUBAR, PRIME, GARD, OHMIC...
├── src/
│   ├── apps/                   # Unified `hyphy3` CLI driver (`hyphy3.cpp`), `hyphy_characterize.cpp`, `hyphy_wasm.cpp`, etc.
│   └── python/                 # Nanobind C++/Python bindings (`_hyphy3`)
├── viz/fel-viewer/             # Embedded interactive HTML5 Canvas viewer for `hyphy3 --html` reports
├── scripts/build_wasm.sh       # Emscripten SIMD128 WebAssembly build script
├── benchmarks/ & tests/        # Doctest suite & gold-standard NEXUS/FASTA benchmark datasets
└── CMakeLists.txt              # C++20 CMake build system
```

---

## 📖 Citations

If you use **HyPhy 3** or **WALRUS** in your research, please cite:

- **WALRUS & HyPhy 3 Core**: Kosakovsky Pond SL, et al. *WALRUS & HyPhy 3: Differentiable Alignment-Aware Molecular Evolution and Unified Selection Forensics*. (2026).
- **BUSTED-MH / BUSTED-S / BUSTED**: Lucaci AG, et al. *Mol Biol Evol.* 38(7):3081–3096 (2021); Wisotsky SR, et al. *Mol Biol Evol.* 37(8):2430–2439 (2020); Murrell B, et al. *Mol Biol Evol.* 32(5):1365–1371 (2015).
- **Contrast-FEL & FEL**: Kosakovsky Pond SL, et al. *Mol Biol Evol.* 38(3):1184–1198 (2021); Kosakovsky Pond SL & Frost SDW. *Mol Biol Evol.* 22(5):1208–1222 (2005).
- **MEME, aBSREL, RELAX & FUBAR**: Murrell B, et al. *PLoS Genet.* 8(7):e1002764 (2012); Smith MD, et al. *Mol Biol Evol.* 32(5):1342–1353 (2015); Wertheim JO, et al. *Mol Biol Evol.* 32(3):820–832 (2015); Murrell B, et al. *Mol Biol Evol.* 30(5):1196–1205 (2013).

---

## 📄 License

Licensed under the [MIT License](LICENSE). Developed by Sergei L. Kosakovsky Pond and the HyPhy development team ([veg.github.io/hyphy-3](https://veg.github.io/hyphy-3/)).
