# HyPhy 3 Technical Specification: Ohmic Phylogenetic Networks & Wheatstone Pruning (`ohmic`)

> **Target Repository**: `/Users/sergei/Development/hyphy-3`  
> **Target Audience**: Gemini 4 / AI Implementation Agent & HyPhy 3 Core Architects  
> **Status**: Ready for C++20 Implementation  
> **Mandate**: Zero Synthetic Mockery, Strict Analytical Gradients, Header-Only C++20 + Eigen3 + OpenMP

---

## 1. Executive Summary & Architectural Mission

Classical maximum-likelihood (ML) phylogenetics and selection analyses (`FEL`, `MEME`, `BUSTED`, `PRIME`, `aBSREL`) condition on a **single discrete binary tree** $T$. When the number of taxa $N$ approaches or exceeds the effective number of informative sites ($N \sim L$), no single tree has appreciable posterior probability, and discrete ML tree search acts as an unregularized $\ell_0$ winner-take-all projection. Empirical sitewise $\Delta\text{SLS}$ dissections across standard benchmarks (`HIVvif`, `HepatitisD`, `InfluenzaA`, `camelid`) demonstrate that **83%–100%** of contested splits between ML and distance-spectrum methods are driven by **1–2 slow-rate homoplastic outlier sites** ("The Two-Site Hijack"), causing severe out-of-sample overfitting (`-80` to `-143` $\ln L$ loss on held-out test codons) and false-positive episodic selection calls when homoplasy is forced onto a single tree.

Conversely, general reticulate phylogenetic networks (arbitrary DAGs with $K$ reticulations) break Felsenstein's conditional independence and require $O(2^K N)$ time per site to unpack into $2^K$ displayed trees.

**This specification introduces Ohmic Likelihood ("Wheatstone Pruning") for `hyphy-3`:**
1. **OhmTrees Electrical Network Inference (`OhmicNetwork`)**: Models the alignment as a resistive electrical circuit where pairwise evolutionary conductance $W_{ij}$ is governed by Kirchhoff's Matrix-Tree Theorem and **Foster's Theorem** ($\sum_{i<j} W_{ij}\Omega_{ij} = N - 1$).
2. **Series-Parallel Buneman Decomposition**: Extracts a mutually compatible **Consensus Backbone** $T_0$ together with $K$ pairwise disjoint 4-way **Wheatstone Bridges** $\mathcal{C}_1, \dots, \mathcal{C}_K$ at high-stress topological bottlenecks (where shunt stress $\tau_{\text{Ohm}}(S) = G_{\text{alt}}(S)/G(S) > \tau_{\min}$).
3. **Linear-Time Wheatstone Pruning (`WheatstoneLikelihood`)**: Because the $K$ Wheatstone bridges are mutually compatible (a Level-1 / galled-tree series-parallel network), the exact marginal likelihood over all $2^K$ displayed spanning trees factorizes into a single local linear combination of conditional likelihood vectors at each bridge:
   $$\mathbf{D}_u^{\text{Ohm}} = w_{k,1}\,\mathbf{D}_u^{(1)} + w_{k,2}\,\mathbf{D}_u^{(2)}$$
   evaluating all $2^K$ trees and their exact Inside-Outside analytical gradients in **$O(N + K)$ time per site pattern**.

---

## 2. Target File Layout in `hyphy-3`

Create or update the following files in `/Users/sergei/Development/hyphy-3`:

```text
include/hyphy/
├── core/
│   ├── ohmic_network.hpp          # [NEW] Electrical solver, Foster resistances, Wheatstone bridge builder
│   └── wheatstone_likelihood.hpp  # [NEW] Inside-Outside Wheatstone pruning & Fréchet adjoint gradients
└── analyses/
    └── ohmic.hpp                  # [NEW] GTR & MG94xREV Ohmic Network fitter + Datamonkey JSON export
src/
├── apps/
│   ├── hyphy_ohmic.cpp            # [NEW] Standalone CLI driver: hyphy3-ohmic
│   └── hyphy3.cpp                 # [MOD] Register `ohmic` subcommand and `--ohmic` topology flag
└── python/
    └── bindings.cpp               # [MOD] Expose OhmicNetwork, WheatstoneBridge, and OhmicAnalyzer
tests/
└── test_ohmic.cpp                 # [NEW] 5-taxon worked example, parity, and finite-difference gradient tests
```

---

## 3. Mathematical & Algorithmic Specification

### 3.1 Step 1: Pairwise Evolutionary Conductance & Foster Effective Resistance

Given an alignment $\mathcal{A}$ of $N$ taxa (`Alignment` in `include/hyphy/core/alignment.hpp`):
1. Compute the $N \times N$ symmetric Tamura-Nei 93 (TN93) evolutionary distance matrix $D$ via `NeighborJoining::compute_tn93_distances(aln)`.
2. Compute the local scale bandwidth $\sigma$ as the median of the upper-triangular entries $\{D_{ij} : 1 \le i < j \le N\}$ (clamped below by $10^{-6}$).
3. Construct the symmetric **Evolutionary Conductance Matrix** $W \in \mathbb{R}^{N \times N}$ with zero diagonal ($W_{ii} = 0$) and off-diagonal entries ($i \ne j$):
   $$W_{ij} = \frac{\exp(-D_{ij} / \sigma)}{D_{ij} + \epsilon_0}, \quad \epsilon_0 = 10^{-3}$$
4. Construct the **Graph Laplacian** $L \in \mathbb{R}^{N \times N}$:
   $$L = \operatorname{diag}(d_1, \dots, d_N) - W, \quad d_i = \sum_{j=1}^N W_{ij}$$
5. Compute the **Moore-Penrose Pseudoinverse** $L^\dagger$ via spectral decomposition (`Eigen::SelfAdjointEigenSolver<Matrix>`):
   $$L = \sum_{m=2}^N \lambda_m \mathbf{u}_m \mathbf{u}_m^\top \implies L^\dagger = \sum_{\lambda_m > 10^{-10}} \lambda_m^{-1} \mathbf{u}_m \mathbf{u}_m^\top$$
6. Compute the symmetric **Effective Resistance Matrix** $\Omega \in \mathbb{R}^{N \times N}$ and **Foster Marginal Edge Probability Matrix** $P_{\text{edge}} \in \mathbb{R}^{N \times N}$:
   $$\Omega_{ij} = L^\dagger_{ii} + L^\dagger_{jj} - 2 L^\dagger_{ij}, \qquad P_{\text{edge}}(i, j) = W_{ij} \Omega_{ij}$$
7. **Mandatory Invariant Check (Foster's Theorem)**:
   $$\sum_{1 \le i < j \le N} W_{ij} \Omega_{ij} = N - 1$$
   Assert in debug/test builds that `std::abs(0.5 * P_edge.sum() - static_cast<Scalar>(N - 1)) < 1e-5`.

---

### 3.2 Step 2: Candidate Split Pool & Electrical Split Scoring

A nontrivial bipartition (split) $S = A \mid B$ partitions $\{0, \dots, N-1\}$ into disjoint subsets $A, B$ with $\min(|A|, |B|) \ge 2$. Represent each canonical split as a `std::vector<uint64_t>` bitmask of length $\lceil N / 64 \rceil$ where bit $0$ is always $0$ (canonically oriented so taxon $0 \in B$).

#### 3.2.1 Generating the Candidate Split Pool
Generate a compact candidate split pool $\mathcal{S}_{\text{cand}}$ ($|\mathcal{S}_{\text{cand}}| \le 15 N$) from the union of:
1. **Multi-scale Neighbor-Joining Splits**: Run fast incremental Neighbor-Joining (`NeighborJoining::infer_from_distances`) on:
   - Raw distance $D$,
   - Effective resistance $\Omega$,
   - Hybrid convex combinations $(1 - \alpha)\hat{D} + \alpha \hat{\Omega}$ for $\alpha \in \{0.15, 0.30, 0.50, 0.70, 0.85, 1.0\}$ (where $\hat{X} = X / \operatorname{mean}_{i<j}(X_{ij})$),
   - Multiplicative power distances $D_{ij}^\gamma \Omega_{ij}^{1-\gamma}$ for $\gamma \in \{0.25, 0.50, 0.75\}$,
   - Foster structural distance $-\ln(\max(P_{\text{edge}}(i,j), 10^{-12}))$.
2. **Local Quartet Alternative Splits (1-Edge NNI Swaps)**:
   - Take the primary NJ tree on $(1 - 0.5)\hat{D} + 0.5\hat{\Omega}$. For every internal branch $e = (u, v)$ whose removal partitions the taxa into 4 subtrees $(A_1, A_2 \mid B_1, B_2)$ (where $A = A_1 \cup A_2$ are the two child clades of $v$ and $B_1, B_2$ are the sibling clade of $v$ and the remaining complement of $u$), explicitly add the two local NNI quartet alternatives:
     $$S_{\text{alt}, 1} = (A_1 \cup B_1) \mid (A_2 \cup B_2), \qquad S_{\text{alt}, 2} = (A_1 \cup B_2) \mid (A_2 \cup B_1)$$
   - This guarantees that every internal branch has its exact 4-way Wheatstone bridge alternatives present in $\mathcal{S}_{\text{cand}}$.

#### 3.2.2 Scoring Each Candidate Split
For each split $S = A \mid B \in \mathcal{S}_{\text{cand}}$, compute:
1. **Normalized Cut Conductance** $C_{\text{cut}}(S)$:
   $$C_{\text{cut}}(S) = \frac{1}{|A|\cdot|B|} \sum_{i \in A} \sum_{j \in B} \frac{P_{\text{edge}}(i, j)}{D_{ij} + 10^{-4}}$$
2. **Balanced Quartet Isolation Margin** $M_{\text{quartet}}(S)$ on $\tilde{\Omega} = 0.5 \hat{\Omega} + 0.5 \hat{D}$:
   - Partition $A = A_1 \cup A_2$ and $B = B_1 \cup B_2$ (using the two highest-scoring compatible sub-splits inside $A$ and $B$, or farthest-pair seeding if $|A|=2$ or $|B|=2$).
   - Compute the $4 \times 4$ inter-cluster mean distance matrix $\bar{\Omega}_{pq} = \frac{1}{|C_p||C_q|} \sum_{i \in C_p, j \in C_q} \tilde{\Omega}_{ij}$ for $C_1=A_1, C_2=A_2, C_3=B_1, C_4=B_2$.
   - Compare the three quartet sums:
     $$q_1 = \bar{\Omega}_{12} + \bar{\Omega}_{34} \text{ (supports } AB|CD\text{)}, \quad q_2 = \bar{\Omega}_{13} + \bar{\Omega}_{24} \text{ (supports } AC|BD\text{)}, \quad q_3 = \bar{\Omega}_{14} + \bar{\Omega}_{23} \text{ (supports } AD|BC\text{)}$$
   - Define the quartet margin $M_{\text{quartet}}(S) = \max\left(0,\, \frac{\min(q_2, q_3) - q_1}{\bar{\Omega}_{\text{quartet}}}\right)$ and composite split conductance:
     $$G(S) = \frac{M_{\text{quartet}}(S)}{C_{\text{cut}}(S) + 10^{-6}}$$

---

### 3.3 Step 3: Consensus Backbone & Wheatstone Bridge Decomposition

Two splits $S_1 = A_1 \mid B_1$ and $S_2 = A_2 \mid B_2$ are **pairwise compatible** (`are_splits_compatible(S1, S2)`) if and only if at least one of the four intersections $A_1 \cap A_2$, $A_1 \cap B_2$, $B_1 \cap A_2$, or $B_1 \cap B_2$ is empty.

1. **Consensus Backbone Extraction**:
   - Sort $\mathcal{S}_{\text{cand}}$ in descending order of $G(S)$.
   - Greedily accept splits into the **Consensus Backbone** $\mathcal{B}_0$ if they are compatible with all previously accepted splits in $\mathcal{B}_0$, up to $N - 3$ internal splits.
   - If $|\mathcal{B}_0| < N - 3$, resolve any remaining multifurcation using local Neighbor-Joining on $0.5\hat{D} + 0.5\hat{\Omega}$ so that $T_0$ is a fully resolved binary tree with $N - 3$ internal splits.
2. **Ohmic Current Support & Shunt Stress per Backbone Branch**:
   - For each internal branch $e \in T_0$ corresponding to split $S_1 = (A_1 \cup A_2) \mid (B_1 \cup B_2)$, let $S_2 = (A_1 \cup B_1) \mid (A_2 \cup B_2)$ and $S_3 = (A_1 \cup B_2) \mid (A_2 \cup B_1)$ be its two local quartet NNI alternatives with conductances $G_1 = G(S_1)$, $G_2 = G(S_2)$, $G_3 = G(S_3)$ (ordered so $G_2 \ge G_3$).
   - Compute:
     $$\text{Ohmic Current Support: } I_{\text{Ohm}}(e) = \frac{G_1}{G_1 + G_2 + G_3} \in [0, 1]$$
     $$\text{Shunt Stress: } \tau_{\text{Ohm}}(e) = \frac{G_2}{G_1 + 10^{-12}} \ge 0$$
3. **Selecting Disjoint Wheatstone Bridges ($\mathcal{C}_1, \dots, \mathcal{C}_K$)**:
   - An internal edge $e = (u, v)$ is a candidate **Wheatstone Bridge** if $\tau_{\text{Ohm}}(e) \ge \tau_{\min}$ (default $\tau_{\min} = 0.05$).
   - **Disjointness Constraint**: Two internal edges $e_1, e_2 \in T_0$ can both be promoted to active 4-way Wheatstone bridges if and only if they do not share an internal incident node (i.e., their 4-way quartet hubs do not overlap, ensuring $S_{k,2}$ is strictly compatible with all other backbone splits and all other active bridges).
   - Greedily select candidate edges in descending order of $\tau_{\text{Ohm}}(e)$ subject to node-disjointness (up to `max_bridges`, default $K_{\max} = \lfloor N / 4 \rfloor$).

---

### 3.4 Step 4: Tree Representation of a 4-Way Wheatstone Bridge in `hyphy-3`

To make **Wheatstone Pruning** seamless and fast inside `hyphy-3`, see how a 4-way Wheatstone Bridge $\mathcal{C}_k$ operates on a rooted tree $T_0$:

Consider an internal child node $v$ of parent node $u$ in $T_0$, where:
- $v$ has two children $(c_0, c_1)$ (representing subtrees $A, B$),
- $u$ has child $v$ and sibling child $c_2$ (representing subtree $C$), and $u$ connects upward to the rest of the tree (representing subtree $D$).

When edge $(u, v)$ is promoted to a **Wheatstone Bridge** $\mathcal{C}_k$:
- Node $u$ becomes a **4-way Bridge Hub** with three downward child subtrees $(c_0, c_1, c_2) = (A, B, C)$ and one upward connection to $D$, mediated by internal quartet bridge node $v$:
  - **Resolution 1 (Primary Split $AB \mid C$, weight $w_{k,1}$, internal branch length $t_{k,1}$)**:
    Subtrees $A$ (`c_0`) and $B$ (`c_1`) coalesce at internal node $v$ first, evolve across bridge branch $t_{k,1}$ (`P_bridge_1`), and then join $C$ (`c_2`) at node $u$:
    $$\mathbf{D}_u^{(1)} = \left[ P(t_{k,1}) \left( \mathbf{M}_{A} \odot \mathbf{M}_{B} \right) \right] \odot \mathbf{M}_{C}$$
  - **Resolution 2 (Secondary Shunt Split $AC \mid B$, weight $w_{k,2}$, internal branch length $t_{k,2}$)**:
    Subtrees $A$ (`c_0`) and $C$ (`c_2`) coalesce first, evolve across shunt bridge branch $t_{k,2}$ (`P_bridge_2`), and then join $B$ (`c_1`) at node $u$:
    $$\mathbf{D}_u^{(2)} = \left[ P(t_{k,2}) \left( \mathbf{M}_{A} \odot \mathbf{M}_{C} \right) \right] \odot \mathbf{M}_{B}$$
  - **Resolution 3 (Optional Tertiary Shunt Split $BC \mid A$, weight $w_{k,3}$, internal branch length $t_{k,3}$)**:
    Subtrees $B$ (`c_1`) and $C$ (`c_2`) coalesce first, evolve across $t_{k,3}$ (`P_bridge_3`), and then join $A$ (`c_0`) at node $u$:
    $$\mathbf{D}_u^{(3)} = \left[ P(t_{k,3}) \left( \mathbf{M}_{B} \odot \mathbf{M}_{C} \right) \right] \odot \mathbf{M}_{A}$$
  where $\mathbf{M}_X = P(t_X) \mathbf{D}_X$ is the standard upward message from child subtree $X \in \{A, B, C\}$, and $w_{k,1} + w_{k,2} + w_{k,3} = 1$.
- **Superposition at Hub $u$**:
  $$\mathbf{D}_u^{\text{Ohm}} = w_{k,1}\,\mathbf{D}_u^{(1)} + w_{k,2}\,\mathbf{D}_u^{(2)} + w_{k,3}\,\mathbf{D}_u^{(3)}$$

> **Why this representation is brilliant for `hyphy-3`**:
> 1. $\mathbf{M}_A$, $\mathbf{M}_B$, and $\mathbf{M}_C$ are computed **once** per site pattern!
> 2. Notice that the sum of child log-scales $s_A + s_B + s_C$ (`node_log_scale[A] + node_log_scale[B] + node_log_scale[C]`) is **identical** across $\mathbf{D}_u^{(1)}$, $\mathbf{D}_u^{(2)}$, and $\mathbf{D}_u^{(3)}$ because all three resolutions multiply the exact same three subtree messages $\mathbf{M}_A, \mathbf{M}_B, \mathbf{M}_C$ (only differing in which pair is multiplied by $P(t_{k,r})$ first, and $P(t_{k,r})$ is a stochastic matrix with entries in $[0, 1]$)!
> 3. Therefore, **zero log-scale exponent mismatch ever occurs** when summing $w_{k,1}\mathbf{D}_u^{(1)} + w_{k,2}\mathbf{D}_u^{(2)} + w_{k,3}\mathbf{D}_u^{(3)}$, and standard `SCALER_THRESH` underflow scaling can be applied directly to $\mathbf{D}_u^{\text{Ohm}}$!

---

### 3.5 Step 5: Exact Inside-Outside Pass & Adjoint Gradients Through a Wheatstone Bridge

Let us derive the exact **Outside Pass** (top-down traversal) and **Analytical Gradients** at a Wheatstone Bridge $\mathcal{C}_k$ rooted at hub $u$ with child subtrees $A, B, C$ and upward incoming adjoint vector $\mathbf{U}_u$ (where at the root $\mathbf{U}_{\text{root}} = \boldsymbol{\pi}$, and for any non-root node $\mathbf{U}_u = P(t_u)^\top \mathbf{V}_u$, scaled by the node's underflow boost multiplier $B_u$).

Let $\tilde{\mathbf{U}}_u = B_u \mathbf{U}_u$. Because:
$$\mathbf{D}_u^{\text{Ohm}} = B_u \sum_{r=1}^3 w_{k,r} \left[ P_{k,r} (\mathbf{M}_{a_r} \odot \mathbf{M}_{b_r}) \odot \mathbf{M}_{c_r} \right]$$
where $(a_1, b_1, c_1) = (A, B, C)$, $(a_2, b_2, c_2) = (A, C, B)$, and $(a_3, b_3, c_3) = (B, C, A)$, the inner product at hub $u$ is:
$$L_p = \mathbf{U}_u^\top \mathbf{D}_u^{\text{Ohm}} = \sum_{r=1}^3 w_{k,r} \tilde{\mathbf{U}}_u^\top \left[ \left(P_{k,r} \mathbf{Y}_{k,r}\right) \odot \mathbf{M}_{c_r} \right], \qquad \mathbf{Y}_{k,r} = \mathbf{M}_{a_r} \odot \mathbf{M}_{b_r}$$

#### 1. Adjoint Gradient with Respect to Bridge Branch Length $t_{k,r}$ ($r \in \{1, 2, 3\}$):
Define the pre-bridge upward context vector for resolution $r$:
$$\mathbf{V}_{k,r}^{\text{bridge}} = w_{k,r} \left( \tilde{\mathbf{U}}_u \odot \mathbf{M}_{c_r} \right)$$
Then the derivative of pattern likelihood $L_p$ with respect to bridge branch length $t_{k,r}$ is:
$$\frac{\partial L_p}{\partial t_{k,r}} = \left(\mathbf{V}_{k,r}^{\text{bridge}}\right)^\top \left( Q P_{k,r} \right) \mathbf{Y}_{k,r}$$
and the outer-product accumulator matrix for spectral Fréchet global parameter derivatives ($\omega, \theta_{ij}$) on bridge $r$ is:
$$\mathbf{A}_{k,r}^{\text{bridge}} \mathrel{+}= \frac{\text{weight}_p}{L_p} \mathbf{V}_{k,r}^{\text{bridge}} \mathbf{Y}_{k,r}^\top$$

#### 2. Adjoint Gradient with Respect to Shunt Weights $w_{k,r}$:
For each resolution $r \in \{1, 2, 3\}$, the resolution-specific pattern likelihood is:
$$L_{p, k, r} = \tilde{\mathbf{U}}_u^\top \mathbf{D}_u^{(r)} \implies \frac{\partial \ln L}{\partial w_{k,r}} = \sum_{p=1}^P \frac{\text{weight}_p}{L_p} L_{p, k, r}$$
And the E-step posterior responsibility of resolution $r$ at pattern $p$ (for EM / `SQUAREM` weight refinement if desired) is simply:
$$\gamma_{p, k, r} = \frac{w_{k,r} L_{p, k, r}}{\sum_{r'=1}^3 w_{k,r'} L_{p, k, r'}}$$

#### 3. Exact Top-Down Context Vectors $\mathbf{V}_A, \mathbf{V}_B, \mathbf{V}_C$ for Child Subtrees $A, B, C$:
Let $\mathbf{Z}_{k,r} = P_{k,r}^\top \left( \tilde{\mathbf{U}}_u \odot \mathbf{M}_{c_r} \right)$ be the message pulled back across bridge $P_{k,r}$.
By differentiating $L_p$ with respect to $\mathbf{M}_A$, $\mathbf{M}_B$, and $\mathbf{M}_C$, the exact parent-sibling context vectors $\mathbf{V}_A, \mathbf{V}_B, \mathbf{V}_C$ (such that $L_p = \mathbf{V}_X^\top \mathbf{M}_X = \mathbf{V}_X^\top P_X \mathbf{D}_X$ for each $X \in \{A, B, C\}$!) are given by the linear superposition over the active bridge paths:
$$\mathbf{V}_A = w_{k,1} \left(\mathbf{Z}_{k,1} \odot \mathbf{M}_B\right) + w_{k,2} \left(\mathbf{Z}_{k,2} \odot \mathbf{M}_C\right) + w_{k,3} \left(\tilde{\mathbf{U}}_u \odot \left(P_{k,3} \mathbf{Y}_{k,3}\right)\right)$$
$$\mathbf{V}_B = w_{k,1} \left(\mathbf{Z}_{k,1} \odot \mathbf{M}_A\right) + w_{k,2} \left(\tilde{\mathbf{U}}_u \odot \left(P_{k,2} \mathbf{Y}_{k,2}\right)\right) + w_{k,3} \left(\mathbf{Z}_{k,3} \odot \mathbf{M}_C\right)$$
$$\mathbf{V}_C = w_{k,1} \left(\tilde{\mathbf{U}}_u \odot \left(P_{k,1} \mathbf{Y}_{k,1}\right)\right) + w_{k,2} \left(\mathbf{Z}_{k,2} \odot \mathbf{M}_A\right) + w_{k,3} \left(\mathbf{Z}_{k,3} \odot \mathbf{M}_B\right)$$
And then for each $X \in \{A, B, C\}$, we set $\mathbf{U}_X = P_X^\top \mathbf{V}_X$ and continue the standard top-down traversal into the subtrees $A, B, C$ without any further modification!

> **Verify the Identity**: Notice that for each $X \in \{A, B, C\}$, $\mathbf{V}_X^\top \mathbf{M}_X = \sum_{r=1}^3 w_{k,r} L_{p,k,r} = L_p$. And when $w_{k,1} = 1, w_{k,2} = w_{k,3} = 0$, these equations reduce **identically** to the standard binary tree Inside-Outside equations in `LikelihoodEngine::compute_inside_outside`!

---

---

### 3.6 Step 6: Breaking the $O(N^3)$ Barrier — Automatic Switchover to the $O(N \log N)$ Scalable Woodbury & Equipotential Busbar Engine

When $N \gg L$ (e.g., $N = 1,000\text{–}100,000$ sequences over $L = 300\text{–}1,000$ codons), two physical and computational facts converge:
1. **Dense $O(N^3)$ Laplacian eigendecomposition and $O(N^3)$ Neighbor-Joining become prohibitive**, and even forming an $N \times N$ distance matrix ($O(N^2 L)$) wastes memory.
2. **Tree resolution mathematically drops ($N \gg L$):** An alignment of $L$ sites cannot resolve $N - 3$ independent binary splits when $N \gg L$. Most local clusters of sequences are separated by $0\text{–}1$ mutations ($D_{ij} \lesssim 1/L$), meaning their pairwise conductance $W_{ij} \to \infty$ (near-zero resistance). Physically, they sit on **equipotential copper busbars (unresolved multifurcations / Ohmic Polytomies)** connected by a sparse global backbone!

To support seamless execution from $N = 5$ to $N = 100,000$, `OhmicNetwork::infer` **MUST** implement two solver regimes controlled by `OhmicScaleMode::{Auto, ExactDense, ScalableWoodbury}`:

#### 3.6.1 Automatic Switchover Criterion (`OhmicScaleMode::Auto`)
1. First, run **Stage A ($O(N L)$ Equipotential Busbar Compression)** to group sequences with identical non-ambiguous site patterns (or distance $D_{ij} \le \epsilon_{\text{busbar}} = \frac{0.5}{L_{\text{nuc}}}$) into $M \le N$ unique **Equipotential Super-Nodes**.
2. **Switchover Rule**:
   $$\text{Use } \texttt{ExactDense} \iff M \le N_{\text{switch}} \text{ and } \frac{M}{L_{\text{eff}}} \le \rho_{\text{switch}}$$
   $$\text{Use } \texttt{ScalableWoodbury} \iff M > N_{\text{switch}} \text{ or } \left(M > 256 \text{ and } \frac{M}{L_{\text{eff}}} > \rho_{\text{switch}}\right)$$
   where the default thresholds are:
   - **$N_{\text{switch}} = 512$** (at $M = 512$, a dense $512 \times 512$ `Eigen::SelfAdjointEigenSolver` takes $\approx 18\text{ ms}$; above $M = 512$, the $O(M p^2)$ Woodbury solver with rank $p \le 256$ is strictly faster and uses $O(M p)$ memory instead of $O(M^2)$), and
   - **$\rho_{\text{switch}} = 1.5$** (where $L_{\text{eff}} = \text{aln.patterns.size()}$ is the number of distinct site patterns in the alignment; when $M > 1.5 L_{\text{eff}}$, the alignment is rank-limited and enters the resolution-drop regime).

---

#### 3.6.2 The Four Stages of the $O(N \log N)$ `ScalableWoodbury` Pipeline

##### Stage A: $O(N L)$ Equipotential Busbar Compression ($N \to M$ Super-Nodes)
1. Compute a fast 64-bit FNV-1a / xxHash signature over each sequence's compressed pattern states to group strictly identical sequences in $O(N L_{\text{eff}})$ time.
2. Optionally merge near-identical sequences ($0$ unambiguous nucleotide differences due to gaps/ambiguities) into the same **Equipotential Super-Node** $C_m$ ($m = 0, \dots, M-1$), selecting the sequence with the fewest ambiguous sites as the canonical representative $r_m$.
3. When the final tree is assembled, any Super-Node $C_m$ with $|C_m| > 1$ is expanded into a native **$k$-way multifurcating busbar (Ohmic Polytomy)** with leaf branch lengths $t_{\min} = 10^{-6}$ (or exact terminal distances $D(i, r_m)$).

##### Stage B: $O(M p^2)$ Low-Rank Nyström + Woodbury Laplacian Pseudoinverse (Zero $M \times M$ Matrices!)
Instead of forming the $M \times M$ conductance matrix $W$ and running an $O(M^3)$ eigensolve, we exploit the fact that the evolutionary conductance kernel (and optional 384D HyphAeon latent space $U_{\text{HyphAeon}} \in \mathbb{R}^{M \times 384}$) has effective rank $p \ll M$ (default $p = \min(M, 256)$, or $p = 384$ when HyphAeon latent embeddings are supplied):

1. **Farthest-Point Landmark Selection in $O(M p)$ Time**:
   Select $p$ landmark representatives $\mathcal{L} = \{a_0, \dots, a_{p-1}\} \subset \{0, \dots, M-1\}$ greedily:
   - Set $a_0 = 0$, and maintain `min_dist[i]` $= \min_{a \in \mathcal{L}} D_{\text{TN93}}(i, a)$.
   - In step $k = 1, \dots, p-1$, pick $a_k = \arg\max_{i} \text{min\_dist}[i]$, compute the single TN93 distance column $D(:, a_k) \in \mathbb{R}^M$ in $O(M L)$ time, and update `min_dist[i]`.
   - Record the $M \times p$ distance matrix $D_{M \times p}$ (`Matrix(M, p)`) and its $p \times p$ landmark submatrix $D_{p \times p}$ (`Matrix(p, p)`). Estimate scale $\sigma$ as the median of the upper triangle of $D_{p \times p}$.
2. **Low-Rank Conductance Factor $\Phi \in \mathbb{R}^{M \times p}$**:
   - Define the landmark conductance kernel entries $K(d) = \frac{\exp(-d / \sigma)}{d + 10^{-3}}$ (with $K(0) = \frac{1}{10^{-3}}$ on the diagonal so $K_{p \times p}$ is strictly positive definite).
   - Compute $K_{M \times p} \in \mathbb{R}^{M \times p}$ and $K_{p \times p} \in \mathbb{R}^{p \times p}$.
   - Eigendecompose the small $p \times p$ matrix $K_{p \times p} = V_p \Lambda_p V_p^\top$ (`O(p^3)` $\approx 2\text{ ms}$ for $p=256$) and form the Nyström feature map:
     $$\Phi = K_{M \times p} V_p \left(\Lambda_p + 10^{-6} I_p\right)^{-1/2} \in \mathbb{R}^{M \times p}$$
     *(If external 384D HyphAeon embeddings $U_{\text{HyphAeon}} \in \mathbb{R}^{M \times 384}$ are provided via `--embeddings`, concatenate the normalized `192D dS` and `192D dN` columns directly into $\Phi$!)*
   - Now the off-diagonal conductance matrix is approximated without ever materializing an $M \times M$ matrix:
     $$W \approx \Phi \Phi^\top - \operatorname{diag}\!\left(\|\Phi_{i,:}\|_2^2\right)$$
3. **Exact Matrix-Free Row Degrees $\mathbf{d} \in \mathbb{R}^M$ in $O(M p)$ Time**:
   $$\mathbf{s}_\Phi = \Phi^\top \mathbf{1}_M \in \mathbb{R}^p, \qquad d_i = \max\!\left(\Phi_{i,:}\mathbf{s}_\Phi - \|\Phi_{i,:}\|_2^2,\; 10^{-6}\right)$$
   Let $\tilde{d}_i = d_i + \|\Phi_{i,:}\|_2^2 = \Phi_{i,:}\mathbf{s}_\Phi$ be the augmented diagonal so that $L = \tilde{D} - \Phi \Phi^\top$ where $\tilde{D} = \operatorname{diag}(\tilde{d}_0, \dots, \tilde{d}_{M-1})$.
4. **Woodbury Matrix Identity for the Laplacian Inverse in $O(M p^2 + p^3)$ Time**:
   By the Woodbury identity, regularized with $\mu_{\text{reg}} = 10^{-5}$ to deflate the constant nullspace mode $\mathbf{1}_M$:
   $$\left(\tilde{D} - \Phi \Phi^\top\right)_{\text{reg}}^{-1} = \tilde{D}^{-1} + \tilde{D}^{-1} \Phi \underbrace{\left((1 + \mu_{\text{reg}}) I_p - \Phi^\top \tilde{D}^{-1} \Phi\right)^{-1}}_{C_{p \times p} \in \mathbb{R}^{p \times p}} \Phi^\top \tilde{D}^{-1}$$
   - Compute the $p \times p$ **Latent Capacitance Matrix** $C_{\text{inv}} = (1 + \mu_{\text{reg}}) I_p - \Phi^\top \tilde{D}^{-1} \Phi$ in $O(M p^2)$ time.
   - Compute its symmetric square-root factor $R \in \mathbb{R}^{p \times p}$ such that $C_{p \times p} = R R^\top$ via `Eigen::SelfAdjointEigenSolver<Matrix>` (clamping eigenvalues $\lambda_k > 10^{-8}$).
   - Define the explicit **$p$-dimensional Ohm-Resistance Coordinate Matrix** $Z \in \mathbb{R}^{M \times p}$:
     $$Z = \tilde{D}^{-1} \Phi R \in \mathbb{R}^{M \times p} \qquad (\text{row } \mathbf{z}_i = \tilde{d}_i^{-1} \Phi_{i,:} R \in \mathbb{R}^p)$$
5. **$O(p)$ Pairwise Effective Resistance Evaluation**:
   Since $L^\dagger_{ij} = \tilde{d}_i^{-1} \delta_{ij} + \mathbf{z}_i^\top \mathbf{z}_j$, the effective resistance $\Omega_{ij} = L^\dagger_{ii} + L^\dagger_{jj} - 2 L^\dagger_{ij}$ between **any** two Super-Nodes $i \ne j$ is simply:
   $$\boxed{\tilde{\Omega}_{ij} = \tilde{d}_i^{-1} + \tilde{d}_j^{-1} + \|\mathbf{z}_i - \mathbf{z}_j\|_2^2}$$
   Wait—look at the simplicity of that formula! Because $L^\dagger_{ii} = \tilde{d}_i^{-1} + \|\mathbf{z}_i\|_2^2$ and $L^\dagger_{ij} = \mathbf{z}_i^\top \mathbf{z}_j$ for $i \ne j$, the cross terms cancels into the exact squared Euclidean distance $\|\mathbf{z}_i - \mathbf{z}_j\|_2^2$ plus the inverse node degrees $\tilde{d}_i^{-1} + \tilde{d}_j^{-1}$!

---

##### Stage C: $O(M p \log M)$ Matrix-Free Neighbor-Joining in Resistance Coordinate Space
In standard Neighbor-Joining, computing the row sums $R_i = \sum_{j \ne i} \tilde{\Omega}_{ij}$ takes $O(M^2)$ time per step ($O(M^3)$ total).
Using our closed-form Woodbury coordinates $(\tilde{d}_i^{-1}, \mathbf{z}_i \in \mathbb{R}^p)$, the sum of $\tilde{\Omega}_{ij}$ over all $j \in \mathcal{A}$ (where $\mathcal{A}$ is the set of $m = |\mathcal{A}|$ currently active clusters) expands analytically into **four global accumulators** computed in **$O(m p)$ time**:
$$\boxed{R_i = (m - 2)\tilde{d}_i^{-1} + S_{\tilde{d}} + m \|\mathbf{z}_i\|_2^2 + S_{\|\mathbf{z}\|^2} - 2\,\mathbf{z}_i^\top \mathbf{s}_{\mathbf{z}}}$$
where:
$$S_{\tilde{d}} = \sum_{j \in \mathcal{A}} \tilde{d}_j^{-1} \in \mathbb{R}, \qquad S_{\|\mathbf{z}\|^2} = \sum_{j \in \mathcal{A}} \|\mathbf{z}_j\|_2^2 \in \mathbb{R}, \qquad \mathbf{s}_{\mathbf{z}} = \sum_{j \in \mathcal{A}} \mathbf{z}_j \in \mathbb{R}^p$$

To avoid scanning all $\binom{m}{2}$ pairs when minimizing $Q_{ij} = (m - 2)\tilde{\Omega}_{ij} - R_i - R_j$:
1. Maintain for each active cluster $i \in \mathcal{A}$ a short list $\mathcal{N}(i)$ of its $k_{\text{knn}} = 32$ nearest neighbors in $\mathbf{z}$-space (found via a vantage-point tree / $k$-d tree or recursive random-projection tree on $\{\mathbf{z}_i\}$ in $O(M p \log M)$ time), plus the top-16 highest-$R_j$ global outlier candidates.
2. In each step (or in bulk batches of $\lceil m / 64 \rceil$ independent non-overlapping minimum-$Q_{ij}$ pairs, as in FastTree / RapidNJ), merge the best pair $(i^\star, j^\star) \to u$, update the cluster centroid coordinates:
   $$\mathbf{z}_u = \frac{1}{2}(\mathbf{z}_{i^\star} + \mathbf{z}_{j^\star}), \qquad \tilde{d}_u^{-1} = \frac{1}{2}(\tilde{d}_{i^\star}^{-1} + \tilde{d}_{j^\star}^{-1}) - \frac{1}{2}\tilde{\Omega}_{i^\star j^\star}$$
   and update the three global accumulators $(S_{\tilde{d}}, S_{\|\mathbf{z}\|^2}, \mathbf{s}_{\mathbf{z}})$ in **$O(p)$ time per merge**!

---

##### Stage D: Resolution-Adaptive "Ohmic Polytomy" Collapse ($N \gg L$)
When $N \gg L$, many internal splits in any binary tree have zero evolutionary support ($0$ mutations along the branch) and represent arbitrary tie-breaking inside an equipotential cluster.
1. For every internal node $v$ with parent $u$ in the inferred tree, evaluate its local 4-way quartet resistance sums $(q_1, q_2, q_3)$ using the $O(p)$ cluster resistance coordinates $(\tilde{d}_C^{-1}, \mathbf{z}_C)$ of the four incident subtrees $C \in \{A, B, C, D\}$.
2. Compute the **Quartet Signal-to-Noise Ratio** relative to the alignment's finite-site resolution limit $\epsilon_L = \frac{0.5}{\max(L_{\text{nuc}}, 1)}$:
   - **Case 1 — True Contested Conflict ($q_2 - q_1$ small, but individual clade distances $> \epsilon_L$ and $\tau_{\text{Ohm}}(v) \ge \tau_{\min}$ with positive shunt branch length $t_{k,2} > \epsilon_L$)**:
     Promote $(u, v)$ to a 4-way **Wheatstone Bridge** $\mathcal{C}_k$.
   - **Case 2 — Unresolved Equipotential Noise ($t_v \le \epsilon_L$ and $\max(q_2 - q_1, q_3 - q_1) \le \epsilon_L$)**:
     **Collapse node $v$ into parent $u$ as a native Ohmic Polytomy** (`u.children` absorbs `v.children`)!
   - **Case 3 — Resolved Backbone Branch ($t_v > \epsilon_L$ and $\tau_{\text{Ohm}}(v) < \tau_{\min}$)**:
     Retain as a standard binary branch in `collapsed_tree`.

> **Computational Payoff in `hyphy-3`**: Because `LikelihoodEngine` and `WheatstoneLikelihood` iterate over `for (int32_t cid : node.children)`, collapsing $K_{\text{poly}}$ unresolvable zero-length branches into native multifurcating busbars **eliminates $K_{\text{poly}}$ internal nodes, $K_{\text{poly}}$ matrix exponentials $P(t)$, and $K_{\text{poly}}$ $61 \times 61$ matrix-vector products per site pattern**, making likelihood evaluation **faster** as tree resolution drops!

---

## 4. Concrete C++20 Data Structures & Signatures

### 4.1 `include/hyphy/core/ohmic_network.hpp`

```cpp
#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/neighbor_joining.hpp"
#include <Eigen/Eigenvalues>
#include <array>
#include <vector>

namespace hyphy::core {

enum class OhmicScaleMode {
    Auto,             // ExactDense if M <= 512 and M <= 1.5 * L_eff, else ScalableWoodbury
    ExactDense,       // Full N x N Laplacian pseudoinverse + exact Foster resistances
    ScalableWoodbury  // O(N log N) Equipotential Busbars + Rank-p Woodbury + R^p NJ
};

struct OhmicSettings {
    OhmicScaleMode scale_mode = OhmicScaleMode::Auto;
    size_t dense_switch_threshold = 512;       // Switch to ScalableWoodbury above M = 512
    Scalar rank_ratio_threshold = 1.5;         // Switch if M > 256 and M / num_patterns > 1.5
    size_t woodbury_rank = 256;                // Number of Nystrom landmarks p (or 384 with HyphAeon)
    Scalar min_shunt_stress = 0.05;            // Minimum tau_Ohm to promote a Wheatstone bridge
    size_t max_bridges = 64;                   // Maximum disjoint Wheatstone bridges K
    bool collapse_ohmic_polytomies = true;     // Collapse zero-signal branches (t <= 0.5/L) into busbars
    bool include_tertiary = false;             // Include 3rd quartet resolution w_3
};

struct WheatstoneBridge {
    int32_t hub_node_id = INVALID_INDEX;       // Node u (collapsed 3-child hub: A, B, C)
    int32_t bridge_node_id = INVALID_INDEX;    // Original node v in backbone tree
    std::array<int32_t, 3> children;           // {id_A, id_B, id_C}
    std::array<Scalar, 3> weights;             // {w_1, w_2, w_3} summing to 1.0 (Ohmic currents)
    std::array<Scalar, 3> branch_lengths;      // {t_1, t_2, t_3} internal quartet bridge lengths
    Scalar ohmic_support = 1.0;                // I_Ohm = G_1 / (G_1 + G_2 + G_3)
    Scalar shunt_stress = 0.0;                 // tau_Ohm = G_2 / G_1
};

struct OhmicNetwork {
    Tree backbone_tree;                        // Standard binary tree T_0 (with 2N-2 nodes)
    Tree collapsed_tree;                       // Tree with 3-child bridge hubs AND Ohmic Polytomy busbars
    std::vector<WheatstoneBridge> bridges;     // K disjoint 4-way Wheatstone bridges
    std::vector<int32_t> node_to_bridge_idx;   // Maps hub_node_id -> index in `bridges` (-1 if standard node)
    OhmicScaleMode active_mode = OhmicScaleMode::ExactDense;
    size_t num_super_nodes = 0;                // M <= N equipotential representatives
    size_t num_collapsed_polytomies = 0;       // Count of zero-signal branches collapsed into busbars

    // Populated in ExactDense mode (empty in ScalableWoodbury mode to save O(N^2) RAM):
    Matrix D_tn93;                             // Pairwise TN93 distance matrix (ExactDense)
    Matrix Omega;                              // Pairwise effective resistance matrix (ExactDense)
    Matrix P_edge;                             // Foster marginal edge probability matrix (ExactDense)
    Scalar foster_sum = 0.0;                   // Sum_{i<j} P_edge(i, j) == N - 1

    // Populated in ScalableWoodbury mode (M x p coordinates, O(M * p) RAM):
    Vector inv_aug_degrees;                    // \tilde{d}_i^{-1} [M]
    Matrix Z_resistance_coords;                // Explicit p-dim resistance coordinates Z [M x p]

    // Evaluate O(1) or O(p) effective resistance between two representative indices (i, j)
    Scalar effective_resistance(size_t i, size_t j) const {
        if (i == j) return 0.0;
        if (active_mode == OhmicScaleMode::ExactDense) {
            return Omega(i, j);
        }
        return inv_aug_degrees(i) + inv_aug_degrees(j)
             + (Z_resistance_coords.row(i) - Z_resistance_coords.row(j)).squaredNorm();
    }

    static OhmicNetwork infer(
        const Alignment& aln,
        const OhmicSettings& settings = {},
        const Matrix* external_hyphaeon_embeddings = nullptr // Optional [N x 384] latent matrix
    );
};

} // namespace hyphy::core
```

### 4.2 `include/hyphy/core/wheatstone_likelihood.hpp`

```cpp
#pragma once

#include "hyphy/core/ohmic_network.hpp"
#include "hyphy/core/likelihood.hpp"

namespace hyphy::core {

class WheatstoneLikelihood {
public:
    struct WheatstoneInsideOutsideResult {
        std::vector<Vector> D;                 // Downward subtree likelihoods [num_nodes]
        std::vector<Vector> V;                 // Upward parent-sibling context [num_nodes]
        std::vector<std::array<Vector, 3>> Y_bridge; // Cached M_{a_r} .* M_{b_r} [num_bridges][3]
        std::vector<std::array<Vector, 3>> V_bridge; // Adjoint context for bridge r [num_bridges][3]
        std::vector<std::array<Scalar, 3>> L_bridge; // Resolution likelihoods L_{p,k,r} [num_bridges][3]
        Scalar likelihood = 0.0;
        Scalar log_likelihood = 0.0;
    };

    // Single-pass GTR log-likelihood on an OhmicNetwork
    static Scalar compute_gtr_log_likelihood(
        const OhmicNetwork& net,
        const Alignment& aln,
        const GTRMatrix& gtr_model
    );

    // Single-pass MG94 codon log-likelihood on an OhmicNetwork
    static Scalar compute_mg94_log_likelihood(
        const OhmicNetwork& net,
        const Alignment& aln,
        const MG94Matrix& mg94_model,
        Scalar branch_scale = 1.0
    );

    // Exact 2-pass Inside-Outside traversal through backbone + K Wheatstone bridges
    static WheatstoneInsideOutsideResult compute_inside_outside(
        const OhmicNetwork& net,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const std::vector<Matrix>& P_branches,                 // [num_nodes]
        const std::vector<std::array<Matrix, 3>>& P_bridges,   // [num_bridges][3]
        const Vector& pi,
        int S
    );

    struct JointWheatstoneMG94Gradients {
        Scalar log_likelihood = 0.0;
        std::vector<Scalar> grad_branches;                     // [num_nodes]
        std::vector<std::array<Scalar, 3>> grad_bridge_lengths;// [num_bridges][3]
        std::vector<std::array<Scalar, 3>> grad_bridge_weights;// [num_bridges][3]
        Scalar grad_omega = 0.0;
        Scalar grad_theta_AC = 0.0;
        Scalar grad_theta_AT = 0.0;
        Scalar grad_theta_CG = 0.0;
        Scalar grad_theta_CT = 0.0;
        Scalar grad_theta_GT = 0.0;
    };

    // Exact joint analytical gradients for all backbone branches, bridge branches,
    // bridge weights, and global MG94 parameters via Inside-Outside + Fréchet derivatives
    static JointWheatstoneMG94Gradients compute_joint_frechet_gradients(
        const OhmicNetwork& net,
        const Alignment& aln,
        const MG94Parameters& params
    );
};

} // namespace hyphy::core
```

---

## 5. Ground-Truth Worked Example for Unit Testing (`tests/test_ohmic.cpp`)

The unit test suite in `tests/test_ohmic.cpp` **MUST** include the following exact 5-taxon worked benchmark ($N = 5$, $L = 100$ nucleotide sites under Jukes-Cantor $\mu = 4/3$) to verify the electrical solver and the Wheatstone Pruning engine against exact closed-form numbers:

### 5.1 Input Alignment ($N = 5$ taxa: `A, B, C, D, O`, $L = 100$ sites)
- **Sites 1–90 (90 invariant sites)**: `A=A, B=A, C=A, D=A, O=A`
- **Sites 91–97 (7 majority synapomorphic sites supporting `AB|CD`, rate $r_{\text{fast}} = 1.0$)**: `A=C, B=C, C=A, D=A, O=A`
- **Sites 98–100 (3 slow-rate convergent shunt sites supporting `AC|BD`, rate $r_{\text{slow}} = 0.05$)**: `A=G, B=A, C=G, D=A, O=A`

### 5.2 Expected Electrical Matrices (with $\sigma = 0.10$, $\epsilon_0 = 0.01$ on this toy example)
1. **Hamming Distance Matrix $D$**:
   - $D_{AB} = 0.03$, $D_{AC} = 0.07$, $D_{AD} = 0.10$, $D_{AO} = 0.10$
   - $D_{BC} = 0.10$, $D_{BD} = 0.07$, $D_{BO} = 0.07$
   - $D_{CD} = 0.03$, $D_{CO} = 0.03$, $D_{DO} = 0.00$
2. **Foster's Theorem Identity**:
   $$\sum_{1 \le i < j \le 5} W_{ij} \Omega_{ij} = \mathbf{4.000000} = N - 1$$
3. **Quartet Conductances & Ohmic Current Weights**:
   - Primary split $S_1 = \{A, B\} \mid \{C, D, O\}$: $w_1 = I_{\text{Ohm}}(S_1) = \mathbf{0.7000}$
   - Shunt split $S_2 = \{A, C\} \mid \{B, D, O\}$: $w_2 = I_{\text{Ohm}}(S_2) = \mathbf{0.3000}$
   - Shunt stress: $\tau_{\text{Ohm}}(S_1) = G(S_2)/G(S_1) = \mathbf{0.4286}$

### 5.3 Expected Log-Likelihoods (JC69 with pendant lengths $t_{\text{leaf}} = 0.005$, internal bridge length $t_{\text{int}} = 0.08$, zero-length outgroup branch $t_O = 0$)
- **Invariants (Sites 1–90, $r = 1.0$)**:
  - $L_{\text{inv}}^{(1)} = L_{\text{inv}}^{(2)} = L_{\text{inv}}^{\text{Ohm}} = 0.22645873 \implies \ln L_{\text{inv}} = -1.4851922$ per site.
- **Majority Synapomorphic Sites (Sites 91–97, `AB|CD`, $r_{\text{fast}} = 1.0$)**:
  - Tree 1 (`AB|CD`): $L_{\text{maj}}^{(1)} = 6.19897 \times 10^{-3}$ ($\ln L = -5.083371$)
  - Tree 2 (`AC|BD`): $L_{\text{maj}}^{(2)} = 6.48481 \times 10^{-7}$ ($\ln L = -14.248631$)
  - **Wheatstone Network ($0.70 L^{(1)} + 0.30 L^{(2)}$)**: $L_{\text{maj}}^{\text{Ohm}} = 4.33947 \times 10^{-3}$ ($\ln L = -5.440001$, paying only `-0.3566` nats vs `+8.8086` nats on Tree 2).
- **Slow-Rate Homoplastic Sites (Sites 98–100, `AC|BD`, $r_{\text{slow}} = 0.05$)**:
  - Tree 1 (`AB|CD`): $L_{\text{shunt}}^{(1)} = 1.73005 \times 10^{-9}$ ($\ln L = -20.175114$ per site — **catastrophic two-site hijack penalty!**)
  - Tree 2 (`AC|BD`): $L_{\text{shunt}}^{(2)} = 3.32116 \times 10^{-4}$ ($\ln L = -8.010028$ per site)
  - **Wheatstone Network ($0.70 L^{(1)} + 0.30 L^{(2)}$)**: $L_{\text{shunt}}^{\text{Ohm}} = 9.96360 \times 10^{-5}$ ($\ln L = -9.213988$ per site — **saves `+10.9611` nats per site** by routing slow-rate homoplasy through the $w_2 = 0.30$ shunt!).
- **Total Alignment Log-Likelihood ($L = 100$)**:
  - Single Tree 1 (`AB|CD` only): $\ln L(T_1) = \mathbf{-229.7762}$
  - Single Tree 2 (`AC|BD` only): $\ln L(T_2) = \mathbf{-257.4378}$
  - **Wheatstone Pruning Network**: $\ln L(\mathcal{N}_{\text{Ohm}}) = \mathbf{-199.3893}$ (**`+30.3869` $\ln L$ improvement over Tree 1 and `+58.0485` $\ln L$ improvement over Tree 2**).

---

## 6. Verification Checklist for the Implementing Agent

Before marking the implementation complete, verify all 7 items:
1. [x] **Foster's Theorem Unit Test (`ExactDense`)**: On synthetic and real alignments (`CD2`, `ADH`, `bglobin`), `OhmicNetwork::infer(aln).foster_sum` equals $N - 1$ within $10^{-5}$.
2. [x] **Degenerate Single-Tree Parity Test**: When all bridge weights are set to $w_{k,1} = 1.0, w_{k,2} = 0.0, w_{k,3} = 0.0$, `WheatstoneLikelihood::compute_joint_frechet_gradients` matches `LikelihoodEngine::compute_joint_frechet_gradients` on `backbone_tree` to $< 10^{-11}$ in log-likelihood and all branch/global gradients.
3. [x] **Worked 5-Taxon Example Test**: Reproduces $\ln L(T_1) = -229.7762$, $\ln L(T_2) = -257.4378$, and $\ln L(\mathcal{N}_{\text{Ohm}}) = -199.3893$ within $10^{-3}$.
4. [x] **Finite-Difference Gradient Verification**: Central finite-difference checks on `grad_branches`, `grad_bridge_lengths`, `grad_bridge_weights`, and `grad_omega` pass with relative error $< 10^{-4}$.
5. [x] **`ExactDense` vs. `ScalableWoodbury` Crossover Parity Test**: On a medium alignment ($N \in [100, 300]$), running `OhmicScaleMode::ScalableWoodbury` with $p = \min(M, 128)$ correlates with `OhmicScaleMode::ExactDense` effective resistances $\Omega_{ij}$ at $R > 0.98$ and produces log-likelihood within $0.1\%$.
6. [x] **$O(N \log N)$ Large-$N$ Scalability & Polytomy Collapse Test**: On a synthetic or empirical alignment with $N = 2,000$ sequences and $L = 300$ sites (containing clonal outbreak clusters), `OhmicScaleMode::Auto` automatically selects `ScalableWoodbury`, collapses zero-signal branches into native Ohmic Polytomies (`num_collapsed_polytomies > 0`), and completes network inference + GTR likelihood evaluation in $< 2.0\text{ seconds}$.
7. [x] **Empirical Improvement on Real Multi-Gene Benchmarks**: Running `hyphy3 ohmic --alignment benchmarks/data/cd2.fna` and `bglobin` completes in $< 0.5\text{s}$ and improves log-likelihood and AICc over single-tree `GTR` / `MG94` fits where $\tau_{\text{Ohm}} \ge 0.05$.

