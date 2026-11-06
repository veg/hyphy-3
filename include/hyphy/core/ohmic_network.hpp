#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/neighbor_joining.hpp"
#include "hyphy/core/hyphaeon.hpp"
#include <Eigen/Eigenvalues>
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>

namespace hyphy::core {

enum class OhmicScaleMode {
    Auto,             // ExactDense if M <= 512 and M <= 1.5 * L_eff, else ScalableWoodbury
    ExactDense,       // Full N x N Laplacian pseudoinverse + exact Foster resistances
    ScalableWoodbury  // O(N log N) Equipotential Busbars + Rank-p Woodbury + R^p NJ
};

enum class ScalableBackboneMethod {
    SpectralFosterHybrid, // Top-down matrix-free spectral bisection + local Foster BioNJ (Default)
    WoodburyWard          // Nearest-Neighbor chain Ward variance clustering in Woodbury space
};

struct OhmicSettings {
    OhmicScaleMode scale_mode = OhmicScaleMode::Auto;
    ScalableBackboneMethod scalable_method = ScalableBackboneMethod::WoodburyWard;
    size_t dense_switch_threshold = 1024;      // Switch to ScalableWoodbury above M = 1024
    Scalar rank_ratio_threshold = 1.5;         // Switch if M > 1024 and M / num_patterns > 1.5
    size_t woodbury_rank = 256;                // Number of Nystrom landmarks p (or 384 with HyphAeon)
    size_t foster_cluster_threshold = 128;     // Cluster size threshold T to transition from Spectral Cut to Foster BioNJ
    Scalar min_shunt_stress = 0.05;            // Minimum tau_Ohm to promote a Wheatstone bridge
    size_t max_bridges = 64;                   // Maximum disjoint Wheatstone bridges K
    bool collapse_ohmic_polytomies = true;     // Collapse zero-signal branches (t <= 0.5/L) into busbars
    bool include_tertiary = false;             // Include 3rd quartet resolution w_3
    bool use_hamming_distance = false;         // Optional raw p-distance (for closed-form benchmarks)
    bool use_site_rate_weights = false;        // Empirical site-rate variability weighting to dampen homoplasy
    bool perform_nni_refinement = false;       // Deprecated: discrete combinatorial search purged in favor of pure algebraic geometry
    Scalar fixed_sigma = -1.0;                 // Optional fixed scale bandwidth sigma (> 0 overrides median)
    Scalar epsilon_0 = 1e-3;                   // Conductance kernel regularization epsilon_0
    bool use_hyphaeon_transformer = false;     // Extract 384D HyphAeon (192D dS + 192D dN) via aeon_core (enabled in CLI tools)
    std::string embeddings_filepath = "";      // Optional path to precomputed [N x 384] embeddings
};

struct OhmicNetwork {
    Tree backbone_tree;                        // Standard binary tree T_0 (with 2N-2 nodes)
    Tree collapsed_tree;                       // Tree with 3-child bridge hubs AND Ohmic Polytomy busbars
    std::vector<WheatstoneBridge> bridges;     // K disjoint 4-way Wheatstone bridges
    std::vector<int32_t> node_to_bridge_idx;   // Maps hub_node_id -> index in `bridges` (-1 if standard node)
    OhmicScaleMode active_mode = OhmicScaleMode::ExactDense;
    size_t num_super_nodes = 0;                // M <= N equipotential representatives
    size_t num_collapsed_polytomies = 0;       // Count of zero-signal branches collapsed into busbars
    bool used_hyphaeon_embeddings = false;     // True if 384D HyphAeon (192D dS + 192D dN) was fused
    std::string embeddings_row_mapping;        // "by_taxon_name" or "positional_alignment_order" (precomputed file only)
    Matrix hyphaeon_embeddings;                // [N x 384] latent representations when extracted/loaded

    // Populated in ExactDense mode (empty in ScalableWoodbury mode to save O(N^2) RAM):
    Matrix D_tn93;                             // Pairwise TN93 distance matrix (ExactDense)
    Matrix D_backbone;                         // Synonymous-stabilized backbone distance matrix (ExactDense)
    Matrix Omega;                              // Pairwise effective resistance matrix (ExactDense)
    Matrix P_edge;                             // Foster marginal edge probability matrix (ExactDense)
    Scalar foster_sum = 0.0;                   // Sum_{i<j} P_edge(i, j) == N - 1

    // Populated in ScalableWoodbury mode (M x p coordinates, O(M * p) RAM):
    std::vector<size_t> taxon_to_super_node;   // Maps alignment taxon index t -> super-node index m in [0, M)
    Vector inv_aug_degrees;                    // \tilde{d}_i^{-1} [M]
    Matrix Z_resistance_coords;                // Explicit p-dim resistance coordinates Z [M x p]

    // Evaluate O(1) or O(p) effective resistance between two representative indices (i, j)
    Scalar effective_resistance(size_t i, size_t j) const {
        if (i == j) return 0.0;
        if (active_mode == OhmicScaleMode::ExactDense) {
            if (static_cast<Eigen::Index>(i) < Omega.rows() && static_cast<Eigen::Index>(j) < Omega.cols()) {
                return Omega(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
            }
            return 1.0;
        }
        if (static_cast<Eigen::Index>(i) >= inv_aug_degrees.size() ||
            static_cast<Eigen::Index>(j) >= inv_aug_degrees.size() ||
            static_cast<Eigen::Index>(i) >= Z_resistance_coords.rows() ||
            static_cast<Eigen::Index>(j) >= Z_resistance_coords.rows()) {
            return 1.0;
        }
        return inv_aug_degrees(static_cast<Eigen::Index>(i)) +
               inv_aug_degrees(static_cast<Eigen::Index>(j)) +
               (Z_resistance_coords.row(static_cast<Eigen::Index>(i)) -
                Z_resistance_coords.row(static_cast<Eigen::Index>(j))).squaredNorm();
    }

    // Evaluate effective resistance between two alignment taxon indices (t_i, t_j) in [0, N)
    // across both ExactDense (N x N Omega) and ScalableWoodbury (M x p Z_resistance_coords) modes (#36)
    Scalar taxon_effective_resistance(size_t t_i, size_t t_j) const {
        if (t_i == t_j) return 0.0;
        if (active_mode == OhmicScaleMode::ExactDense &&
            static_cast<Eigen::Index>(t_i) < Omega.rows() &&
            static_cast<Eigen::Index>(t_j) < Omega.cols()) {
            return Omega(static_cast<Eigen::Index>(t_i), static_cast<Eigen::Index>(t_j));
        }
        size_t m_i = (t_i < taxon_to_super_node.size()) ? taxon_to_super_node[t_i] : t_i;
        size_t m_j = (t_j < taxon_to_super_node.size()) ? taxon_to_super_node[t_j] : t_j;
        if (m_i == m_j) return 0.0;
        return effective_resistance(m_i, m_j);
    }

    // Compute empirical site variability weights w_s = 1 / sqrt(r_s + 0.05)
    // where r_s is Gini impurity / heterozygosity across valid nucleotides at site s.
    static std::vector<Scalar> compute_site_weights(const Alignment& aln) {
        if (aln.num_taxa == 0) return {};
        size_t L = aln.sequences[0].raw_seq.size();
        for (const auto& s : aln.sequences) {
            L = std::min(L, s.raw_seq.size());
        }
        if (L == 0) return {};

        std::vector<Scalar> weights(L, 1.0);
        Scalar sum_w = 0.0;

        for (size_t s = 0; s < L; ++s) {
            std::array<size_t, 4> counts = {0, 0, 0, 0};
            size_t total = 0;
            for (size_t i = 0; i < aln.num_taxa; ++i) {
                int b = char_to_nuc(aln.sequences[i].raw_seq[s]);
                if (b >= 0 && b < 4) {
                    counts[b]++;
                    total++;
                }
            }
            if (total < 2) {
                weights[s] = 1.0;
            } else {
                Scalar r_s = 1.0;
                Scalar inv_tot = 1.0 / static_cast<Scalar>(total);
                for (size_t b = 0; b < 4; ++b) {
                    Scalar freq = static_cast<Scalar>(counts[b]) * inv_tot;
                    r_s -= freq * freq;
                }
                weights[s] = 1.0 / std::sqrt(std::max<Scalar>(0.0, r_s) + 0.05);
            }
            sum_w += weights[s];
        }

        if (sum_w > 0.0) {
            Scalar scale = static_cast<Scalar>(L) / sum_w;
            for (size_t s = 0; s < L; ++s) {
                weights[s] *= scale;
            }
        }
        return weights;
    }

    // Helper: compute pairwise Hamming or TN93 distance between two taxa (i, j) in O(L) time
    static Scalar compute_pair_distance(
        const Alignment& aln,
        size_t i,
        size_t j,
        const Vector4& pi,
        Scalar k1_denom,
        Scalar k2_denom,
        Scalar k3_denom,
        bool use_hamming,
        const std::vector<Scalar>& site_weights = {}
    ) {
        if (i == j) return 0.0;
        const std::string& seq_i = aln.sequences[i].raw_seq;
        const std::string& seq_j = aln.sequences[j].raw_seq;
        size_t L = std::min(seq_i.size(), seq_j.size());

        Scalar valid_weight = 0.0;
        Scalar ag_transitions = 0.0;
        Scalar ct_transitions = 0.0;
        Scalar transversions = 0.0;
        bool has_weights = (site_weights.size() >= L);

        for (size_t s = 0; s < L; ++s) {
            int b1 = char_to_nuc(seq_i[s]);
            int b2 = char_to_nuc(seq_j[s]);
            if (b1 < 0 || b2 < 0) continue;
            Scalar w = has_weights ? site_weights[s] : 1.0;
            valid_weight += w;
            if (b1 != b2) {
                bool p1 = (b1 == 0 || b1 == 2);
                bool p2 = (b2 == 0 || b2 == 2);
                if (p1 && p2) {
                    ag_transitions += w;
                } else if (!p1 && !p2) {
                    ct_transitions += w;
                } else {
                    transversions += w;
                }
            }
        }

        if (valid_weight <= 0.0) return 0.0;
        Scalar P1 = ag_transitions / valid_weight;
        Scalar P2 = ct_transitions / valid_weight;
        Scalar Q = transversions / valid_weight;
        Scalar raw_p = P1 + P2 + Q;

        if (use_hamming) {
            return raw_p;
        }

        Scalar gR = pi(0) + pi(2);
        Scalar gY = pi(1) + pi(3);
        Scalar arg1 = std::max<Scalar>(0.02, 1.0 - P1 / k1_denom - Q / (2.0 * gR));
        Scalar arg2 = std::max<Scalar>(0.02, 1.0 - P2 / k2_denom - Q / (2.0 * gY));
        Scalar arg3 = std::max<Scalar>(0.02, 1.0 - Q / k3_denom);

        Scalar term1 = (2.0 * pi(0) * pi(2) / gR) * std::log(arg1);
        Scalar term2 = (2.0 * pi(1) * pi(3) / gY) * std::log(arg2);
        Scalar term3 = 2.0 * (gR * gY - (pi(0) * pi(2) * gY / gR) - (pi(1) * pi(3) * gR / gY)) * std::log(arg3);
        Scalar d = -(term1 + term2 + term3);
        return std::max<Scalar>(0.0, d);
    }

    static Matrix compute_hamming_distances(const Alignment& aln) {
        size_t N = aln.num_taxa;
        Matrix D = Matrix::Zero(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
        Vector4 dummy_pi = Vector4::Constant(0.25);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Scalar d = compute_pair_distance(aln, i, j, dummy_pi, 0.25, 0.25, 0.5, true);
                D(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = d;
                D(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = d;
            }
        }
        return D;
    }

    // For coding alignments (L >= 30 and L % 3 == 0), downweight multi-nucleotide non-synonymous
    // codon jumps (the hallmark of functional convergence) when inferring the vertical backbone tree,
    // using the exact Tamura-Nei 1993 (TN93) purine/pyrimidine transition and transversion formula
    // with depth-adaptive blending to prevent synonymous saturation rollover on deep vertebrate clades.
    static Matrix compute_synonymous_backbone_distances(const Alignment& aln, const Matrix& D_tn93) {
        size_t N = aln.num_taxa;
        if (N == 0 || aln.sequences.empty() || aln.num_codons == 0) return D_tn93;
        size_t L = aln.sequences[0].raw_seq.size();
        if (L < 30 || (L % 3) != 0) return D_tn93;

        Vector4 pi = aln.nuc_frequencies;
        if (pi.sum() < 0.99) {
            pi = Vector4::Constant(0.25);
        }
        for (int b = 0; b < 4; ++b) {
            if (pi(b) < 1e-4) pi(b) = 1e-4;
        }
        pi /= pi.sum();
        Scalar gR = pi(0) + pi(2);
        Scalar gY = pi(1) + pi(3);
        Scalar k1_denom = 2.0 * pi(0) * pi(2) / gR;
        Scalar k2_denom = 2.0 * pi(1) * pi(3) / gY;
        Scalar k3_denom = 2.0 * gR * gY;

        auto gcode = aln.code ? aln.code : GeneticCode::universal();
        Matrix D_syn = Matrix::Zero(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
        Scalar sum_tn93 = 0.0, sum_syn = 0.0;

        for (size_t i = 0; i < N; ++i) {
            const std::string& seq_i = aln.sequences[i].raw_seq;
            for (size_t j = i + 1; j < N; ++j) {
                const std::string& seq_j = aln.sequences[j].raw_seq;
                size_t L_pair = std::min(seq_i.size(), seq_j.size());
                Scalar d0 = D_tn93(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                Scalar w_multi = 0.45 * (d0 * d0) / (d0 * d0 + 0.0144); // 0.12^2 = 0.0144
                Scalar valid_nt = 0.0;
                Scalar ag_ts = 0.0, ct_ts = 0.0, tv = 0.0;

                for (size_t s = 0; s + 2 < L_pair; s += 3) {
                    int8_t c_i = gcode->triplet_to_sense(seq_i[s], seq_i[s + 1], seq_i[s + 2]);
                    int8_t c_j = gcode->triplet_to_sense(seq_j[s], seq_j[s + 1], seq_j[s + 2]);
                    if (c_i < 0 || c_j < 0) continue;
                    valid_nt += 3.0;
                    if (c_i == c_j) continue;

                    int ndiff = gcode->count_nucleotide_differences(c_i, c_j);
                    Scalar w = 0.0;
                    if (gcode->is_synonymous(c_i, c_j)) {
                        w = 1.0;
                    } else if (ndiff == 1) {
                        w = 0.45;
                    } else {
                        w = w_multi;
                    }
                    if (w > 0.0) {
                        for (size_t pos = s; pos < s + 3; ++pos) {
                            int b1 = char_to_nuc(seq_i[pos]);
                            int b2 = char_to_nuc(seq_j[pos]);
                            if (b1 >= 0 && b2 >= 0 && b1 != b2) {
                                bool p1 = (b1 == 0 || b1 == 2);
                                bool p2 = (b2 == 0 || b2 == 2);
                                if (p1 && p2) ag_ts += w;
                                else if (!p1 && !p2) ct_ts += w;
                                else tv += w;
                            }
                        }
                    }
                }

                Scalar d_val = d0;
                if (valid_nt >= 30.0) {
                    Scalar P1 = ag_ts / valid_nt;
                    Scalar P2 = ct_ts / valid_nt;
                    Scalar Q = tv / valid_nt;
                    Scalar arg1 = std::max<Scalar>(0.02, 1.0 - P1 / k1_denom - Q / (2.0 * gR));
                    Scalar arg2 = std::max<Scalar>(0.02, 1.0 - P2 / k2_denom - Q / (2.0 * gY));
                    Scalar arg3 = std::max<Scalar>(0.02, 1.0 - Q / k3_denom);
                    Scalar term1 = k1_denom * std::log(arg1);
                    Scalar term2 = k2_denom * std::log(arg2);
                    Scalar term3 = 2.0 * (gR * gY - (pi(0) * pi(2) * gY / gR) - (pi(1) * pi(3) * gR / gY)) * std::log(arg3);
                    d_val = std::max<Scalar>(1e-6, -(term1 + term2 + term3));
                }
                D_syn(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = d_val;
                D_syn(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = d_val;
                if (d0 < 0.30) {
                    sum_tn93 += d0;
                    sum_syn += d_val;
                }
            }
        }
        if (sum_syn > 1e-12 && sum_tn93 > 1e-12) {
            D_syn *= (sum_tn93 / sum_syn);
        }
        Matrix D_out = D_tn93;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Scalar d0 = D_tn93(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                Scalar alpha_syn = 0.85 * std::exp(-(d0 / 0.35) * (d0 / 0.35));
                Scalar db = alpha_syn * D_syn(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) + (1.0 - alpha_syn) * d0;
                D_out(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = db;
                D_out(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = db;
            }
        }
        return D_out;
    }

    // Split bitmask representation
    using SplitMask = std::vector<uint64_t>;

    static SplitMask make_canonical_split(const std::vector<int>& members, size_t N) {
        size_t num_words = (N + 63) / 64;
        SplitMask mask(num_words, 0ULL);
        for (int idx : members) {
            if (idx >= 0 && static_cast<size_t>(idx) < N) {
                mask[static_cast<size_t>(idx) / 64] |= (1ULL << (static_cast<size_t>(idx) % 64));
            }
        }
        // Canonical orientation: taxon 0 must have bit 0 == 0 (taxon 0 in B)
        if (mask[0] & 1ULL) {
            for (size_t w = 0; w < num_words; ++w) {
                mask[w] = ~mask[w];
            }
            size_t rem = N % 64;
            if (rem != 0) {
                mask[num_words - 1] &= ((1ULL << rem) - 1ULL);
            }
        }
        return mask;
    }

    static size_t split_popcount(const SplitMask& mask) {
        size_t cnt = 0;
        for (uint64_t w : mask) {
            cnt += static_cast<size_t>(__builtin_popcountll(w));
        }
        return cnt;
    }

    static bool is_nontrivial_split(const SplitMask& mask, size_t N) {
        size_t a = split_popcount(mask);
        return (a >= 2 && (N - a) >= 2);
    }

    static bool are_splits_compatible(const SplitMask& s1, const SplitMask& s2, size_t N) {
        size_t num_words = s1.size();
        bool a1_a2 = false, a1_b2 = false, b1_a2 = false, b1_b2 = false;
        for (size_t w = 0; w < num_words; ++w) {
            uint64_t valid_bits = ~0ULL;
            if (w + 1 == num_words && (N % 64) != 0) {
                valid_bits = (1ULL << (N % 64)) - 1ULL;
            }
            uint64_t m1 = s1[w] & valid_bits;
            uint64_t m2 = s2[w] & valid_bits;
            uint64_t n1 = (~s1[w]) & valid_bits;
            uint64_t n2 = (~s2[w]) & valid_bits;
            if (m1 & m2) a1_a2 = true;
            if (m1 & n2) a1_b2 = true;
            if (n1 & m2) b1_a2 = true;
            if (n1 & n2) b1_b2 = true;
        }
        return !(a1_a2 && a1_b2 && b1_a2 && b1_b2);
    }

    // Collect descendant leaf taxon indices (0..N-1) for every node in a tree
    static std::vector<std::vector<int>> collect_subtree_taxa(
        const Tree& tree,
        const Alignment& aln
    ) {
        size_t num_nodes = tree.num_nodes();
        std::vector<std::vector<int>> clade_taxa(num_nodes);
        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            if (node.is_leaf) {
                auto it = aln.taxon_to_index.find(node.name);
                if (it != aln.taxon_to_index.end()) {
                    clade_taxa[nid].push_back(static_cast<int>(it->second));
                }
            } else {
                for (int32_t cid : node.children) {
                    clade_taxa[nid].insert(
                        clade_taxa[nid].end(),
                        clade_taxa[cid].begin(),
                        clade_taxa[cid].end()
                    );
                }
            }
        }
        return clade_taxa;
    }

    // Convert an unrooted/trifurcating NJ tree into a rooted tree where every internal edge
    // is represented as a parent-child pair (u, v) with u having sibling C and upward/root context D
    static Tree normalize_root_orientation(const Tree& raw_tree) {
        if (raw_tree.root_id == INVALID_INDEX || raw_tree.nodes.empty()) return raw_tree;
        const auto& root = raw_tree.nodes[raw_tree.root_id];
        if (root.children.size() <= 2) return raw_tree;

        // If root has 3 children {c0, c1, c2} and at least one is an internal node,
        // ensure that internal children are visited cleanly; if root has 3 leaf children (e.g. N=3), keep as-is.
        return raw_tree;
    }

    // Load precomputed [N x D] (default D=384) HyphAeon embeddings from .bin, .npy, or text (.tsv/.csv)
    static Matrix load_embeddings_file(const std::string& filepath, size_t expected_N,
                                       const std::vector<std::string>* taxon_names = nullptr,
                                       std::string* mapping = nullptr) {
        if (mapping) *mapping = "positional_alignment_order";  // binary formats carry no names
        std::ifstream f(filepath, std::ios::binary);
        if (!f.is_open()) {
            throw std::runtime_error("OhmicNetwork::load_embeddings_file: Cannot open " + filepath);
        }

        f.seekg(0, std::ios::end);
        size_t file_bytes = static_cast<size_t>(f.tellg());
        f.seekg(0, std::ios::beg);

        // Case 1: Binary file with 8-byte header [uint32_t rows, uint32_t cols] followed by float64
        if (file_bytes >= 8) {
            uint32_t hdr[2] = {0, 0};
            f.read(reinterpret_cast<char*>(hdr), 8);
            if (hdr[0] == expected_N && hdr[1] > 0 &&
                file_bytes == 8 + static_cast<size_t>(hdr[0]) * static_cast<size_t>(hdr[1]) * sizeof(double)) {
                Matrix M(hdr[0], hdr[1]);
                std::vector<double> buf(static_cast<size_t>(hdr[0]) * hdr[1]);
                f.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(double));
                for (size_t r = 0; r < hdr[0]; ++r) {
                    for (size_t c = 0; c < hdr[1]; ++c) {
                        M(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = buf[r * hdr[1] + c];
                    }
                }
                return M;
            }
        }

        // Case 2: Raw float64 [expected_N x 384] or float32 [expected_N x 384]
        if (expected_N > 0 && (file_bytes % (expected_N * sizeof(double))) == 0) {
            size_t cols = file_bytes / (expected_N * sizeof(double));
            if (cols == 384 || cols == 192 || cols == 256) {
                f.seekg(0, std::ios::beg);
                std::vector<double> buf(expected_N * cols);
                f.read(reinterpret_cast<char*>(buf.data()), file_bytes);
                Matrix M(static_cast<Eigen::Index>(expected_N), static_cast<Eigen::Index>(cols));
                for (size_t r = 0; r < expected_N; ++r) {
                    for (size_t c = 0; c < cols; ++c) {
                        M(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = buf[r * cols + c];
                    }
                }
                return M;
            }
        }
        if (expected_N > 0 && (file_bytes % (expected_N * sizeof(float))) == 0) {
            size_t cols = file_bytes / (expected_N * sizeof(float));
            if (cols == 384 || cols == 192 || cols == 256) {
                f.seekg(0, std::ios::beg);
                std::vector<float> buf(expected_N * cols);
                f.read(reinterpret_cast<char*>(buf.data()), file_bytes);
                Matrix M(static_cast<Eigen::Index>(expected_N), static_cast<Eigen::Index>(cols));
                for (size_t r = 0; r < expected_N; ++r) {
                    for (size_t c = 0; c < cols; ++c) {
                        M(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = static_cast<Scalar>(buf[r * cols + c]);
                    }
                }
                return M;
            }
        }

        // Case 3: Whitespace or comma-separated text file.
        // Rows may carry a leading taxon-name column (first field non-numeric); such rows are joined to the
        // alignment by name (any order; missing, unknown or duplicate names are errors). Purely numeric rows
        // are taken positionally in alignment order. An optional header line (non-numeric fields) is skipped.
        f.clear();
        f.seekg(0, std::ios::beg);
        auto parse_num = [](const std::string& tok, Scalar& out) {
            char* end = nullptr;
            out = std::strtod(tok.c_str(), &end);
            return end != tok.c_str() && *end == '\0';
        };
        std::vector<std::vector<Scalar>> rows;
        std::vector<std::string> row_names;  // empty string for unnamed rows
        std::vector<std::string> first_toks; // raw first field of each row (taxon names may look numeric)
        size_t n_named = 0;
        bool first_line = true;
        std::string line;
        while (std::getline(f, line)) {
            for (char& ch : line) {
                if (ch == ',' || ch == '\t' || ch == '\r') ch = ' ';
            }
            std::istringstream iss(line);
            std::vector<std::string> toks;
            std::string tok;
            while (iss >> tok) toks.push_back(tok);
            if (toks.empty()) continue;
            std::vector<Scalar> vals;
            std::string name;
            Scalar v;
            size_t start = 0;
            if (!parse_num(toks[0], v)) { name = toks[0]; start = 1; }
            bool ok = true;
            for (size_t k = start; k < toks.size(); ++k) {
                if (!parse_num(toks[k], v)) { ok = false; break; }
                vals.push_back(v);
            }
            if (!ok || (start == 1 && vals.empty())) {
                bool any_numeric = false;
                for (const auto& t : toks) any_numeric = any_numeric || parse_num(t, v);
                if (first_line && !any_numeric) { first_line = false; continue; }  // header line
                throw std::runtime_error("OhmicNetwork::load_embeddings_file: non-numeric embedding value in row '" +
                                         toks[0] + "' of " + filepath);
            }
            first_line = false;
            if (start == 1) ++n_named;
            row_names.push_back(name);
            first_toks.push_back(toks[0]);
            rows.push_back(std::move(vals));
        }
        // Numeric-looking taxon names: if the first fields of an all-numeric file are exactly the alignment's
        // taxon names, it is a named file, not an extra embedding column.
        if (n_named == 0 && taxon_names != nullptr && rows.size() == expected_N && rows[0].size() > 1) {
            std::unordered_set<std::string> names(taxon_names->begin(), taxon_names->end());
            std::unordered_set<std::string> firsts(first_toks.begin(), first_toks.end());
            if (names.size() == expected_N && firsts == names) {
                for (size_t r = 0; r < rows.size(); ++r) {
                    row_names[r] = first_toks[r];
                    rows[r].erase(rows[r].begin());
                }
                n_named = rows.size();
            }
        }
        if (!rows.empty() && n_named != 0) {
            if (n_named != rows.size())
                throw std::runtime_error("OhmicNetwork::load_embeddings_file: " + filepath +
                                         " mixes named and unnamed rows");
            if (taxon_names == nullptr || taxon_names->size() != expected_N)
                throw std::runtime_error("OhmicNetwork::load_embeddings_file: named rows in " + filepath +
                                         " but no alignment taxon names to join on");
            std::unordered_map<std::string, size_t> idx;
            for (size_t i = 0; i < taxon_names->size(); ++i) idx.emplace((*taxon_names)[i], i);
            const size_t cols = rows[0].size();
            Matrix M = Matrix::Zero(static_cast<Eigen::Index>(expected_N), static_cast<Eigen::Index>(cols));
            std::vector<bool> seen(expected_N, false);
            for (size_t r = 0; r < rows.size(); ++r) {
                auto it = idx.find(row_names[r]);
                if (it == idx.end())
                    throw std::runtime_error("OhmicNetwork::load_embeddings_file: taxon '" + row_names[r] +
                                             "' in " + filepath + " is not in the alignment");
                if (seen[it->second])
                    throw std::runtime_error("OhmicNetwork::load_embeddings_file: duplicate taxon '" + row_names[r] +
                                             "' in " + filepath);
                if (rows[r].size() != cols)
                    throw std::runtime_error("OhmicNetwork::load_embeddings_file: row '" + row_names[r] +
                                             "' has " + std::to_string(rows[r].size()) + " values, expected " +
                                             std::to_string(cols) + " in " + filepath);
                seen[it->second] = true;
                for (size_t c = 0; c < cols; ++c)
                    M(static_cast<Eigen::Index>(it->second), static_cast<Eigen::Index>(c)) = rows[r][c];
            }
            for (size_t i = 0; i < expected_N; ++i) {
                if (!seen[i])
                    throw std::runtime_error("OhmicNetwork::load_embeddings_file: alignment taxon '" +
                                             (*taxon_names)[i] + "' has no row in " + filepath);
            }
            if (mapping) *mapping = "by_taxon_name";
            return M;
        }
        if (rows.size() == expected_N && !rows[0].empty()) {
            size_t cols = rows[0].size();
            for (size_t r = 1; r < rows.size(); ++r) {
                if (rows[r].size() != cols)
                    throw std::runtime_error("OhmicNetwork::load_embeddings_file: row " + std::to_string(r + 1) +
                                             " has " + std::to_string(rows[r].size()) + " values, expected " +
                                             std::to_string(cols) + " in " + filepath);
            }
            Matrix M = Matrix::Zero(static_cast<Eigen::Index>(expected_N), static_cast<Eigen::Index>(cols));
            for (size_t r = 0; r < expected_N; ++r) {
                for (size_t c = 0; c < std::min(cols, rows[r].size()); ++c) {
                    M(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = rows[r][c];
                }
            }
            if (mapping) *mapping = "positional_alignment_order";
            return M;
        }
        throw std::runtime_error("OhmicNetwork::load_embeddings_file: Unrecognized format or row count mismatch in " + filepath);
    }

    // Extract [N x 384] HyphAeon Multi-Scale Transformer representations (192D dS + 192D dN) via aeon_core.
    // See hyphy/core/hyphaeon.hpp: no hard-coded paths, prompted weight download, throws on failure.
    static Matrix extract_hyphaeon_embeddings(const Alignment& aln) {
        return hyphaeon::extract_embeddings(aln);
    }

    // Split [N x 2d] HyphAeon embeddings into L2-normalized 192D dS and 192D dN tracks
    static void split_and_normalize_hyphaeon(const Matrix& H, Matrix& norm_ds, Matrix& norm_dn) {
        Eigen::Index rows = H.rows();
        Eigen::Index half = H.cols() / 2;
        norm_ds = H.leftCols(half);
        norm_dn = H.rightCols(H.cols() - half);
        for (Eigen::Index i = 0; i < rows; ++i) {
            Scalar n1 = norm_ds.row(i).norm();
            if (n1 > 1e-12) norm_ds.row(i) /= n1;
            Scalar n2 = norm_dn.row(i).norm();
            if (n2 > 1e-12) norm_dn.row(i) /= n2;
        }
    }

    // Convenience helper for all CLI tools to resolve --tree / --ohmic / --hyphaeon / --embeddings or tree-free FASTA
    static Tree resolve_tree(
        const Alignment& aln,
        const std::string& tree_file,
        bool force_ohmic,
        bool use_hyphaeon = false,
        const std::string& embeddings_file = "",
        OhmicNetwork* out_net = nullptr,
        std::string* out_source_desc = nullptr
    ) {
        bool need_ohmic = force_ohmic || !embeddings_file.empty() ||
                          (tree_file.empty() && aln.embedded_tree_newick.empty());
        if (!need_ohmic) {
            if (!tree_file.empty()) {
                if (out_source_desc) *out_source_desc = tree_file;
                return Tree::from_newick_file(tree_file);
            }
            if (out_source_desc) *out_source_desc = "Embedded in alignment";
            return Tree::from_newick(aln.embedded_tree_newick);
        }

        OhmicSettings settings;
        settings.use_hyphaeon_transformer = use_hyphaeon;
        settings.embeddings_filepath = embeddings_file;
        OhmicNetwork net = infer(aln, settings);
        if (out_source_desc) {
            std::string mode_str = (net.active_mode == OhmicScaleMode::ExactDense) ? "ExactDense" : "ScalableWoodbury";
            std::string tr_str = net.used_hyphaeon_embeddings ? "+HyphAeon384D" : "";
            *out_source_desc = "OhmicNetwork(" + mode_str + tr_str + ", bridges=" +
                               std::to_string(net.bridges.size()) + ", polytomies=" +
                               std::to_string(net.num_collapsed_polytomies) + ")";
        }
        Tree res = net.collapsed_tree;
        if (out_net) *out_net = std::move(net);
        return res;
    }

    static Tree resolve_tree(
        const Alignment& aln,
        const std::string& tree_file,
        bool force_ohmic,
        bool use_hyphaeon,
        const std::string& embeddings_file,
        std::string* out_source_desc
    ) {
        return resolve_tree(aln, tree_file, force_ohmic, use_hyphaeon, embeddings_file, nullptr, out_source_desc);
    }

    static OhmicNetwork infer(
        const Alignment& aln,
        const OhmicSettings& settings = {},
        const Matrix* external_hyphaeon_embeddings = nullptr
    ) {
        OhmicNetwork net;
        size_t N = aln.num_taxa;
        if (N < 4) {
            net.backbone_tree = NeighborJoining::infer(aln);
            net.collapsed_tree = net.backbone_tree;
            net.node_to_bridge_idx.assign(net.collapsed_tree.num_nodes(), -1);
            net.collapsed_tree.node_to_bridge_idx = net.node_to_bridge_idx;
            net.num_super_nodes = N;
            return net;
        }

        // Resolve optional 384D HyphAeon embeddings (from pointer, file, or live aeon_core transformer)
        if (external_hyphaeon_embeddings != nullptr &&
            static_cast<size_t>(external_hyphaeon_embeddings->rows()) == N &&
            external_hyphaeon_embeddings->cols() >= 2) {
            net.hyphaeon_embeddings = *external_hyphaeon_embeddings;
            net.used_hyphaeon_embeddings = true;
        } else if (!settings.embeddings_filepath.empty()) {
            net.hyphaeon_embeddings = load_embeddings_file(settings.embeddings_filepath, N, &aln.taxon_names,
                                                          &net.embeddings_row_mapping);
            net.used_hyphaeon_embeddings = (net.hyphaeon_embeddings.rows() == static_cast<Eigen::Index>(N));
        } else if (settings.use_hyphaeon_transformer && !settings.use_hamming_distance) {
            net.hyphaeon_embeddings = extract_hyphaeon_embeddings(aln);
            net.used_hyphaeon_embeddings = (net.hyphaeon_embeddings.rows() == static_cast<Eigen::Index>(N) &&
                                            net.hyphaeon_embeddings.cols() >= 2);
        }

        Matrix norm_ds, norm_dn;
        if (net.used_hyphaeon_embeddings) {
            split_and_normalize_hyphaeon(net.hyphaeon_embeddings, norm_ds, norm_dn);
        }

        size_t L_eff = std::max<size_t>(1, aln.patterns.size());
        size_t L_nuc = aln.num_nucleotides > 0 ? aln.num_nucleotides : (aln.num_codons * 3);
        if (L_nuc == 0 && !aln.sequences.empty()) {
            L_nuc = aln.sequences[0].raw_seq.size();
        }
        Scalar eps_L = 0.5 / static_cast<Scalar>(std::max<size_t>(1, L_nuc));

        // =========================================================================
        // Stage A: O(N * L_eff) Equipotential Busbar Compression (N -> M Super-Nodes)
        // =========================================================================
        std::vector<size_t> super_node_reps;
        std::vector<std::vector<size_t>> super_node_members;
        std::vector<size_t> taxon_to_super_node(N, 0);

        // Compute 64-bit FNV-1a hash over pattern states for each taxon
        std::unordered_map<uint64_t, std::vector<size_t>> hash_buckets;
        hash_buckets.reserve(N);
        for (size_t i = 0; i < N; ++i) {
            uint64_t h = 14695981039346656037ULL;
            if (!aln.patterns.empty()) {
                for (const auto& pat : aln.patterns) {
                    uint8_t st = static_cast<uint8_t>(pat.states[i] + 2);
                    h ^= st;
                    h *= 1099511628211ULL;
                }
            } else {
                for (char c : aln.sequences[i].raw_seq) {
                    h ^= static_cast<uint8_t>(c);
                    h *= 1099511628211ULL;
                }
            }
            hash_buckets[h].push_back(i);
        }

        // Verify exact pattern equality within each hash bucket to form M super-nodes
        for (size_t i = 0; i < N; ++i) {
            taxon_to_super_node[i] = static_cast<size_t>(-1);
        }
        for (size_t i = 0; i < N; ++i) {
            if (taxon_to_super_node[i] != static_cast<size_t>(-1)) continue;
            uint64_t h = 14695981039346656037ULL;
            if (!aln.patterns.empty()) {
                for (const auto& pat : aln.patterns) {
                    uint8_t st = static_cast<uint8_t>(pat.states[i] + 2);
                    h ^= st;
                    h *= 1099511628211ULL;
                }
            } else {
                for (char c : aln.sequences[i].raw_seq) {
                    h ^= static_cast<uint8_t>(c);
                    h *= 1099511628211ULL;
                }
            }
            size_t m_idx = super_node_reps.size();
            super_node_reps.push_back(i);
            super_node_members.push_back({i});
            taxon_to_super_node[i] = m_idx;

            for (size_t j : hash_buckets[h]) {
                if (j <= i || taxon_to_super_node[j] != static_cast<size_t>(-1)) continue;
                bool identical = true;
                if (!aln.patterns.empty()) {
                    for (const auto& pat : aln.patterns) {
                        if (pat.states[i] != pat.states[j]) {
                            identical = false;
                            break;
                        }
                    }
                } else {
                    identical = (aln.sequences[i].raw_seq == aln.sequences[j].raw_seq);
                }
                if (identical) {
                    super_node_members[m_idx].push_back(j);
                    taxon_to_super_node[j] = m_idx;
                }
            }
        }

        size_t M = super_node_reps.size();
        net.num_super_nodes = M;
        net.taxon_to_super_node = taxon_to_super_node;

        std::vector<Scalar> site_weights;
        if (settings.use_site_rate_weights) {
            site_weights = compute_site_weights(aln);
        }

        // Automatic Switchover Rule (Section 3.6.1)
        OhmicScaleMode mode = settings.scale_mode;
        if (mode == OhmicScaleMode::Auto) {
            Scalar ratio = static_cast<Scalar>(M) / static_cast<Scalar>(L_eff);
            if (M <= settings.dense_switch_threshold &&
                !(M > 1024 && ratio > settings.rank_ratio_threshold)) {
                mode = OhmicScaleMode::ExactDense;
            } else {
                mode = OhmicScaleMode::ScalableWoodbury;
            }
        }
        net.active_mode = mode;

        // =========================================================================
        // Regime 1: ExactDense Solver (Full N x N Laplacian Pseudoinverse & Foster)
        // =========================================================================
        if (mode == OhmicScaleMode::ExactDense) {
            net.D_tn93 = settings.use_hamming_distance
                ? compute_hamming_distances(aln)
                : NeighborJoining::compute_tn93_distances(aln);

            // If 384D HyphAeon embeddings (192D dS + 192D dN) are active, fuse cosine distances
            // matching OhmTrees (ohmtrees.py lines 68-82)
            Matrix sim_ds_mat, sim_dn_mat;
            if (net.used_hyphaeon_embeddings) {
                sim_ds_mat = ((norm_ds * norm_ds.transpose()).array() + 1.0) * 0.5;
                sim_ds_mat = sim_ds_mat.cwiseMax(0.0).cwiseMin(1.0);
                sim_dn_mat = ((norm_dn * norm_dn.transpose()).array() + 1.0) * 0.5;
                sim_dn_mat = sim_dn_mat.cwiseMax(0.0).cwiseMin(1.0);
            }

            // Compute median of upper-triangular distances
            std::vector<Scalar> upper_d;
            upper_d.reserve(N * (N - 1) / 2);
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = i + 1; j < N; ++j) {
                    Scalar d = net.D_tn93(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                    if (d > 0.0) upper_d.push_back(d);
                }
            }
            Scalar sigma = 0.1;
            if (settings.fixed_sigma > 0.0) {
                sigma = settings.fixed_sigma;
            } else if (!upper_d.empty()) {
                std::sort(upper_d.begin(), upper_d.end());
                size_t mid = upper_d.size() / 2;
                sigma = (upper_d.size() % 2 == 1)
                    ? upper_d[mid]
                    : 0.5 * (upper_d[mid - 1] + upper_d[mid]);
                sigma = std::max<Scalar>(sigma, 1e-6);
            }

            Scalar eps0 = settings.epsilon_0;
            Matrix W = Matrix::Zero(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = i + 1; j < N; ++j) {
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    Scalar d = net.D_tn93(ii, jj);
                    Scalar w = std::exp(-d / sigma) / (d + eps0);
                    if (net.used_hyphaeon_embeddings) {
                        w *= (0.50 + 0.35 * sim_ds_mat(ii, jj) + 0.15 * sim_dn_mat(ii, jj));
                    }
                    W(ii, jj) = w;
                    W(jj, ii) = w;
                }
            }

            Vector deg = W.rowwise().sum();
            Matrix L = Matrix(deg.asDiagonal()) - W;
            L = 0.5 * (L + L.transpose());

            Eigen::SelfAdjointEigenSolver<Matrix> eig(L);
            if (eig.info() != Eigen::Success) {
                throw std::runtime_error("OhmicNetwork::infer: Graph Laplacian eigendecomposition failed");
            }

            const Vector& evals = eig.eigenvalues();
            const Matrix& evecs = eig.eigenvectors();
            Matrix L_pinv = Matrix::Zero(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
            for (Eigen::Index m = 0; m < evals.size(); ++m) {
                if (evals(m) > 1e-10) {
                    L_pinv.noalias() += (1.0 / evals(m)) * (evecs.col(m) * evecs.col(m).transpose());
                }
            }

            net.Omega = Matrix::Zero(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
            net.P_edge = Matrix::Zero(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
            Scalar foster_acc = 0.0;
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = i + 1; j < N; ++j) {
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    Scalar om = std::max<Scalar>(0.0, L_pinv(ii, ii) + L_pinv(jj, jj) - 2.0 * L_pinv(ii, jj));
                    net.Omega(ii, jj) = om;
                    net.Omega(jj, ii) = om;
                    Scalar pedge = W(ii, jj) * om;
                    net.P_edge(ii, jj) = pedge;
                    net.P_edge(jj, ii) = pedge;
                    foster_acc += pedge;
                }
            }
            net.foster_sum = foster_acc;
            assert(std::abs(net.foster_sum - static_cast<Scalar>(N - 1)) < 1e-4);

            // Primary backbone tree from NJ on synonymous-stabilized distances (with branch lengths scaled to D_tn93)
            net.D_backbone = compute_synonymous_backbone_distances(aln, net.D_tn93);
            net.backbone_tree = NeighborJoining::infer_from_distances(net.D_backbone, aln.taxon_names);

            // Build collapsed_tree and extract node-disjoint Wheatstone bridges
            build_bridges_and_collapsed_tree(net, aln, settings, eps_L, true);
            return net;
        }

        // =========================================================================
        // Regime 2: ScalableWoodbury O(N log N) Solver
        // =========================================================================
        // Estimate global nucleotide frequencies for TN93
        Vector4 pi = aln.nuc_frequencies;
        if (pi.sum() < 0.99) {
            pi = Vector4::Constant(0.25);
        }
        for (int b = 0; b < 4; ++b) {
            if (pi(b) < 1e-4) pi(b) = 1e-4;
        }
        pi /= pi.sum();
        Scalar gR = pi(0) + pi(2);
        Scalar gY = pi(1) + pi(3);
        Scalar k1_denom = 2.0 * pi(0) * pi(2) / gR;
        Scalar k2_denom = 2.0 * pi(1) * pi(3) / gY;
        Scalar k3_denom = 2.0 * gR * gY;

        // 1. Farthest-Point Landmark Selection in O(M * p) time
        size_t p = std::min<size_t>(M, settings.woodbury_rank);
        std::vector<size_t> landmarks;
        landmarks.reserve(p);
        std::vector<Scalar> min_dist(M, 1e30);
        std::vector<bool> is_landmark(M, false);

        Matrix D_Mp = Matrix::Zero(static_cast<Eigen::Index>(M), static_cast<Eigen::Index>(p));
        size_t cur_lm = 0;
        for (size_t k = 0; k < p; ++k) {
            landmarks.push_back(cur_lm);
            is_landmark[cur_lm] = true;
            size_t rep_lm = super_node_reps[cur_lm];

            Scalar max_d = -1.0;
            size_t next_lm = 0;
            for (size_t i = 0; i < M; ++i) {
                Scalar d = (i == cur_lm) ? 0.0 : compute_pair_distance(
                    aln, super_node_reps[i], rep_lm, pi, k1_denom, k2_denom, k3_denom, settings.use_hamming_distance, site_weights
                );
                D_Mp(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(k)) = d;
                if (d < min_dist[i]) {
                    min_dist[i] = d;
                }
                if (!is_landmark[i] && min_dist[i] > max_d) {
                    max_d = min_dist[i];
                    next_lm = i;
                }
            }
            cur_lm = next_lm;
        }

        Matrix D_pp = Matrix::Zero(static_cast<Eigen::Index>(p), static_cast<Eigen::Index>(p));
        std::vector<Scalar> lm_upper;
        lm_upper.reserve(p * (p - 1) / 2);
        for (size_t a = 0; a < p; ++a) {
            for (size_t b = 0; b < p; ++b) {
                Scalar d = D_Mp(static_cast<Eigen::Index>(landmarks[a]), static_cast<Eigen::Index>(b));
                D_pp(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) = d;
                if (a < b && d > 0.0) lm_upper.push_back(d);
            }
        }
        D_pp = 0.5 * (D_pp + D_pp.transpose());

        Scalar sigma = 0.1;
        if (settings.fixed_sigma > 0.0) {
            sigma = settings.fixed_sigma;
        } else if (!lm_upper.empty()) {
            std::sort(lm_upper.begin(), lm_upper.end());
            size_t mid = lm_upper.size() / 2;
            sigma = (lm_upper.size() % 2 == 1)
                ? lm_upper[mid]
                : 0.5 * (lm_upper[mid - 1] + lm_upper[mid]);
            sigma = std::max<Scalar>(sigma, 1e-6);
        }

        // 2. Low-Rank Nystrom Conductance Feature Map Phi [M x p_total]
        Scalar eps0 = settings.epsilon_0;
        Matrix K_Mp(static_cast<Eigen::Index>(M), static_cast<Eigen::Index>(p));
        for (size_t i = 0; i < M; ++i) {
            for (size_t k = 0; k < p; ++k) {
                Scalar d = D_Mp(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(k));
                K_Mp(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(k)) = std::exp(-d / sigma) / (d + eps0);
            }
        }
        Matrix K_pp(static_cast<Eigen::Index>(p), static_cast<Eigen::Index>(p));
        for (size_t a = 0; a < p; ++a) {
            K_pp.row(static_cast<Eigen::Index>(a)) = K_Mp.row(static_cast<Eigen::Index>(landmarks[a]));
        }
        K_pp = 0.5 * (K_pp + K_pp.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix> eig_K(K_pp);
        Vector evals_K = eig_K.eigenvalues().cwiseMax(1e-6);
        Vector inv_sqrt_K = evals_K.array().rsqrt();
        Matrix Phi = K_Mp * eig_K.eigenvectors() * inv_sqrt_K.asDiagonal();

        // Concatenate normalized 384D HyphAeon (192D dS + 192D dN) columns into Phi (Section 3.6.2 Stage B)
        if (net.used_hyphaeon_embeddings) {
            Eigen::Index ds_cols = norm_ds.cols();
            Eigen::Index dn_cols = norm_dn.cols();
            Eigen::Index ext_cols = ds_cols + dn_cols;
            Scalar phi_scale = std::sqrt(std::max<Scalar>(1.0, Phi.rowwise().squaredNorm().mean()) * 0.25);
            Matrix Phi_aug(static_cast<Eigen::Index>(M), Phi.cols() + ext_cols);
            Phi_aug.leftCols(Phi.cols()) = Phi;
            for (size_t i = 0; i < M; ++i) {
                Eigen::Index rep_i = static_cast<Eigen::Index>(super_node_reps[i]);
                Phi_aug.block(static_cast<Eigen::Index>(i), Phi.cols(), 1, ds_cols) =
                    phi_scale * norm_ds.row(rep_i);
                Phi_aug.block(static_cast<Eigen::Index>(i), Phi.cols() + ds_cols, 1, dn_cols) =
                    phi_scale * norm_dn.row(rep_i);
            }
            Phi = std::move(Phi_aug);
        }

        Eigen::Index p_eff = Phi.cols();

        // 3. Exact Matrix-Free Row Degrees d_tilde [M] in O(M * p) time
        Vector s_Phi = Phi.colwise().sum().transpose();
        Vector row_sq_norm = Phi.rowwise().squaredNorm();
        Vector d_tilde = (Phi * s_Phi).cwiseMax((row_sq_norm.array() + 1e-6).matrix());
        net.inv_aug_degrees = d_tilde.cwiseInverse();

        // 4. Woodbury Capacitance Matrix C_inv = (1 + mu_reg) I_p - Phi^T * D_tilde^{-1} * Phi
        const Scalar mu_reg = 1e-4;
        Matrix scaled_Phi = net.inv_aug_degrees.asDiagonal() * Phi;
        Matrix C_inv = (1.0 + mu_reg) * Matrix::Identity(p_eff, p_eff) - (Phi.transpose() * scaled_Phi);
        C_inv = 0.5 * (C_inv + C_inv.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix> eig_C(C_inv);
        Vector evals_C = eig_C.eigenvalues();
        const Matrix& evecs_C = eig_C.eigenvectors();

        // Deflate the smallest nullspace mode (corresponding to 1_M) and clamp remaining modes
        Matrix R = Matrix::Zero(p_eff, p_eff);
        for (Eigen::Index k = 0; k < p_eff; ++k) {
            if (k == 0 && evals_C(k) < 5.0 * mu_reg) {
                continue;
            }
            Scalar lam = std::max<Scalar>(evals_C(k), 1e-6);
            R.col(k) = evecs_C.col(k) * (1.0 / std::sqrt(lam));
        }

        net.Z_resistance_coords = scaled_Phi * R;
        // Center Z rows so mean(z) = 0 (preserves all pairwise ||z_i - z_j||^2 exactly)
        Eigen::RowVectorXd z_mean = net.Z_resistance_coords.colwise().mean();
        net.Z_resistance_coords.rowwise() -= z_mean;

        // =========================================================================
        // Stage C: Multi-Scale Spectral-Foster Hybrid or Woodbury Ward Clustering
        // =========================================================================
        Tree super_tree;
        if (settings.scalable_method == ScalableBackboneMethod::SpectralFosterHybrid) {
            super_tree = infer_spectral_foster_hybrid(
                Phi,
                aln,
                super_node_reps,
                pi, k1_denom, k2_denom, k3_denom,
                settings.use_hamming_distance,
                settings.foster_cluster_threshold,
                settings.epsilon_0,
                site_weights
            );
        } else {
            super_tree = infer_woodbury_ward(
                net.inv_aug_degrees,
                net.Z_resistance_coords,
                aln,
                super_node_reps,
                pi, k1_denom, k2_denom, k3_denom,
                settings.use_hamming_distance,
                site_weights
            );
        }

        // Expand any Equipotential Super-Node C_m with |C_m| > 1 into a multifurcating busbar
        // and also produce a resolved binary backbone_tree for reference
        expand_super_nodes_into_trees(
            net, super_tree, aln, super_node_reps, super_node_members
        );


        // Stage D: Quartet SNR evaluation, Wheatstone bridge promotion, and Ohmic Polytomy collapse
        build_bridges_and_collapsed_tree(net, aln, settings, eps_L, false);
        net.foster_sum = static_cast<Scalar>(M - 1);
        return net;
    }

    static void rebuild_bridges(
        OhmicNetwork& net,
        const Alignment& aln,
        const OhmicSettings& settings
    ) {
        size_t L_sites = aln.num_nucleotides > 0 ? aln.num_nucleotides : std::max<size_t>(1, aln.num_codons * 3);
        if (L_sites == 0 && !aln.sequences.empty()) {
            L_sites = aln.sequences[0].raw_seq.size();
        }
        Scalar eps_L = 0.5 / static_cast<Scalar>(std::max<size_t>(1, L_sites));
        net.bridges.clear();
        net.num_collapsed_polytomies = 0;
        build_bridges_and_collapsed_tree(
            net, aln, settings, eps_L, net.active_mode == OhmicScaleMode::ExactDense
        );
    }

private:
    // Stage C (Primary): Multi-Scale Spectral-Foster Engine
    // Top-down matrix-free spectral bisection via normalized Laplacian Fiedler vector and Fisher 1D 2-means split,
    // combined with bottom-up local exact Foster BioNJ on terminal micro-clades (size <= threshold_foster).
    static Tree infer_spectral_foster_hybrid(
        const Matrix& Phi,
        const Alignment& aln,
        const std::vector<size_t>& super_node_reps,
        const Vector4& pi,
        Scalar k1_denom,
        Scalar k2_denom,
        Scalar k3_denom,
        bool use_hamming,
        size_t threshold_foster = 128,
        Scalar eps0 = 1e-3,
        const std::vector<Scalar>& site_weights = {}
    ) {
        (void)Phi;
        size_t M = super_node_reps.size();

        Tree tree;
        if (M == 1) {
            TreeNode root;
            root.id = 0;
            root.name = aln.taxon_names[super_node_reps[0]];
            root.is_leaf = true;
            root.branch_length = 0.0;
            tree.nodes.push_back(root);
            tree.root_id = 0;
            tree.leaf_name_to_id[root.name] = 0;
            tree.compute_post_order();
            return tree;
        }
        if (M == 2) {
            Scalar d01 = compute_pair_distance(
                aln, super_node_reps[0], super_node_reps[1], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights
            );
            for (size_t i = 0; i < 2; ++i) {
                TreeNode leaf;
                leaf.id = static_cast<int32_t>(i);
                leaf.name = aln.taxon_names[super_node_reps[i]];
                leaf.is_leaf = true;
                leaf.branch_length = std::max<Scalar>(1e-6, 0.5 * d01);
                leaf.parent_id = 2;
                tree.nodes.push_back(leaf);
                tree.leaf_name_to_id[leaf.name] = leaf.id;
            }
            TreeNode root;
            root.id = 2;
            root.name = "Root";
            root.is_leaf = false;
            root.children = {0, 1};
            tree.nodes.push_back(root);
            tree.root_id = 2;
            tree.compute_post_order();
            return tree;
        }

        // Initialize leaves 0 .. M-1
        tree.nodes.reserve(2 * M + 16);
        for (size_t i = 0; i < M; ++i) {
            TreeNode leaf;
            leaf.id = static_cast<int32_t>(i);
            leaf.name = aln.taxon_names[super_node_reps[i]];
            leaf.is_leaf = true;
            leaf.branch_length = 0.0;
            tree.nodes.push_back(leaf);
            tree.leaf_name_to_id[leaf.name] = leaf.id;
        }

        int internal_counter = 1;

        // Local exact Foster BioNJ on terminal micro-clades (m <= threshold_foster)
        auto solve_local_foster = [&](const std::vector<size_t>& cluster_indices) -> int32_t {
            size_t m = cluster_indices.size();
            if (m == 1) {
                return static_cast<int32_t>(cluster_indices[0]);
            }
            if (m == 2) {
                int32_t u = static_cast<int32_t>(cluster_indices[0]);
                int32_t v = static_cast<int32_t>(cluster_indices[1]);
                Scalar d01 = compute_pair_distance(
                    aln, super_node_reps[u], super_node_reps[v], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights
                );
                tree.nodes[u].branch_length = std::max<Scalar>(1e-6, 0.5 * d01);
                tree.nodes[v].branch_length = std::max<Scalar>(1e-6, 0.5 * d01);
                int32_t w_id = static_cast<int32_t>(tree.nodes.size());
                TreeNode w_node;
                w_node.id = w_id;
                w_node.name = "FosterNode_" + std::to_string(internal_counter++);
                w_node.is_leaf = false;
                w_node.children = {u, v};
                tree.nodes.push_back(w_node);
                tree.nodes[u].parent_id = w_id;
                tree.nodes[v].parent_id = w_id;
                return w_id;
            }

            // 1. Exact pairwise distances D [m x m]
            Matrix D = Matrix::Zero(m, m);
            std::vector<Scalar> upper_d;
            upper_d.reserve(m * (m - 1) / 2);
            for (size_t i = 0; i < m; ++i) {
                size_t rep_i = super_node_reps[cluster_indices[i]];
                for (size_t j = i + 1; j < m; ++j) {
                    size_t rep_j = super_node_reps[cluster_indices[j]];
                    Scalar d = compute_pair_distance(aln, rep_i, rep_j, pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights);
                    D(i, j) = d;
                    D(j, i) = d;
                    if (d > 0.0) upper_d.push_back(d);
                }
            }

            // Local scale sigma
            Scalar loc_sigma = 0.1;
            if (!upper_d.empty()) {
                std::sort(upper_d.begin(), upper_d.end());
                size_t mid = upper_d.size() / 2;
                loc_sigma = (upper_d.size() % 2 == 1) ? upper_d[mid] : 0.5 * (upper_d[mid - 1] + upper_d[mid]);
                loc_sigma = std::max<Scalar>(loc_sigma, 1e-6);
            }

            // Conductance matrix C [m x m] and degree
            Matrix C = Matrix::Zero(m, m);
            Vector deg = Vector::Zero(m);
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    if (i == j) continue;
                    Scalar d = D(i, j);
                    Scalar c = std::exp(-d / loc_sigma) / (d + eps0);
                    C(i, j) = c;
                    deg(i) += c;
                }
            }

            // Laplacian pseudoinverse
            Matrix L = -C;
            L.diagonal() += deg;
            Eigen::SelfAdjointEigenSolver<Matrix> eig_L(L);
            const Vector& evals_L = eig_L.eigenvalues();
            const Matrix& evecs_L = eig_L.eigenvectors();
            Matrix L_pinv = Matrix::Zero(m, m);
            for (Eigen::Index k = 1; k < static_cast<Eigen::Index>(m); ++k) {
                if (evals_L(k) > 1e-8) {
                    L_pinv.noalias() += (1.0 / evals_L(k)) * (evecs_L.col(k) * evecs_L.col(k).transpose());
                }
            }

            // Foster resistance Omega
            Matrix Omega = Matrix::Zero(m, m);
            Scalar sum_om = 0.0;
            Scalar sum_d = 0.0;
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    if (i != j) {
                        Scalar om = std::max<Scalar>(0.0, L_pinv(i, i) + L_pinv(j, j) - 2.0 * L_pinv(i, j));
                        Omega(i, j) = om;
                        if (i < j) {
                            sum_om += om;
                            sum_d += D(i, j);
                        }
                    }
                }
            }
            if (sum_om > 1e-12 && sum_d > 1e-12) {
                Omega *= (sum_d / sum_om);
            }

            // BioNJ agglomeration on Omega with evolutionary distances D
            std::vector<int32_t> cur_nodes(m);
            for (size_t i = 0; i < m; ++i) cur_nodes[i] = static_cast<int32_t>(cluster_indices[i]);
            Matrix cur_dist = Omega;
            Matrix cur_D = D;

            while (cur_nodes.size() > 2) {
                size_t K = cur_nodes.size();
                Vector r = cur_dist.rowwise().sum();

                Scalar best_q = std::numeric_limits<Scalar>::infinity();
                size_t best_i = 0, best_j = 1;
                Scalar K_minus_2 = static_cast<Scalar>(K - 2);

                for (size_t i = 0; i < K; ++i) {
                    for (size_t j = i + 1; j < K; ++j) {
                        Scalar q = K_minus_2 * cur_dist(i, j) - r(i) - r(j);
                        if (q < best_q) {
                            best_q = q;
                            best_i = i;
                            best_j = j;
                        }
                    }
                }

                int32_t u = cur_nodes[best_i];
                int32_t v = cur_nodes[best_j];
                Scalar d_uv = cur_D(best_i, best_j);
                Scalar om_uv = cur_dist(best_i, best_j);
                Scalar r_diff = (K > 2) ? (r(best_i) - r(best_j)) / K_minus_2 : 0.0;
                Scalar len_u = std::clamp<Scalar>(0.5 * (d_uv + r_diff), 0.01 * d_uv, 0.99 * d_uv);
                Scalar len_v = std::max<Scalar>(1e-6, d_uv - len_u);
                len_u = std::max<Scalar>(1e-6, len_u);

                tree.nodes[u].branch_length = len_u;
                tree.nodes[v].branch_length = len_v;

                int32_t w_id = static_cast<int32_t>(tree.nodes.size());
                TreeNode w_node;
                w_node.id = w_id;
                w_node.name = "FosterNode_" + std::to_string(internal_counter++);
                w_node.is_leaf = false;
                w_node.children = {u, v};
                tree.nodes.push_back(w_node);
                tree.nodes[u].parent_id = w_id;
                tree.nodes[v].parent_id = w_id;

                size_t new_K = K - 1;
                Matrix new_dist(new_K, new_K);
                Matrix new_D(new_K, new_K);
                std::vector<int32_t> next_nodes(new_K);

                size_t dest = 0;
                std::vector<size_t> old_to_new(K, 999999);
                for (size_t k = 0; k < K; ++k) {
                    if (k != best_i && k != best_j) {
                        old_to_new[k] = dest;
                        next_nodes[dest] = cur_nodes[k];
                        dest++;
                    }
                }
                next_nodes[new_K - 1] = w_id;

                for (size_t a = 0; a < K; ++a) {
                    if (a == best_i || a == best_j) continue;
                    size_t na = old_to_new[a];
                    for (size_t b = a + 1; b < K; ++b) {
                        if (b == best_i || b == best_j) continue;
                        size_t nb = old_to_new[b];
                        Scalar om_ab = cur_dist(a, b);
                        new_dist(na, nb) = om_ab;
                        new_dist(nb, na) = om_ab;
                        Scalar d_ab = cur_D(a, b);
                        new_D(na, nb) = d_ab;
                        new_D(nb, na) = d_ab;
                    }
                    Scalar om_wa = std::max<Scalar>(1e-6, 0.5 * (cur_dist(best_i, a) + cur_dist(best_j, a) - om_uv));
                    new_dist(na, new_K - 1) = om_wa;
                    new_dist(new_K - 1, na) = om_wa;

                    Scalar d_wa = std::max<Scalar>(1e-6, 0.5 * (cur_D(best_i, a) + cur_D(best_j, a) - d_uv));
                    new_D(na, new_K - 1) = d_wa;
                    new_D(new_K - 1, na) = d_wa;
                }
                new_dist(new_K - 1, new_K - 1) = 0.0;
                new_D(new_K - 1, new_K - 1) = 0.0;

                cur_nodes = std::move(next_nodes);
                cur_dist = std::move(new_dist);
                cur_D = std::move(new_D);
            }

            int32_t u = cur_nodes[0];
            int32_t v = cur_nodes[1];
            Scalar d_uv = cur_D(0, 1);
            tree.nodes[u].branch_length = std::max<Scalar>(1e-6, 0.5 * d_uv);
            tree.nodes[v].branch_length = std::max<Scalar>(1e-6, 0.5 * d_uv);

            int32_t sub_root = static_cast<int32_t>(tree.nodes.size());
            TreeNode sr;
            sr.id = sub_root;
            sr.name = "FosterSubRoot_" + std::to_string(internal_counter++);
            sr.is_leaf = false;
            sr.children = {u, v};
            tree.nodes.push_back(sr);
            tree.nodes[u].parent_id = sub_root;
            tree.nodes[v].parent_id = sub_root;
            return sub_root;
        };

        // Recursive top-down spectral bisection
        std::function<int32_t(const std::vector<size_t>&)> recurse_hybrid;
        recurse_hybrid = [&](const std::vector<size_t>& indices) -> int32_t {
            size_t n = indices.size();
            if (n <= threshold_foster) {
                return solve_local_foster(indices);
            }

            size_t p_loc = std::min<size_t>(n, 32);
            std::vector<size_t> lm_local;
            lm_local.reserve(p_loc);
            lm_local.push_back(0);

            std::vector<Scalar> min_d(n, std::numeric_limits<Scalar>::infinity());
            Matrix D_loc(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(p_loc));

            for (size_t k = 0; k < p_loc; ++k) {
                size_t lm_idx = lm_local[k];
                size_t rep_lm = super_node_reps[indices[lm_idx]];
                size_t best_next = 0;
                Scalar max_d = -1.0;

                for (size_t i = 0; i < n; ++i) {
                    size_t rep_i = super_node_reps[indices[i]];
                    Scalar d = (i == lm_idx) ? 0.0 : compute_pair_distance(
                        aln, rep_i, rep_lm, pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights
                    );
                    D_loc(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(k)) = d;
                    if (d < min_d[i]) min_d[i] = d;
                    if (k + 1 < p_loc && min_d[i] > max_d) {
                        max_d = min_d[i];
                        best_next = i;
                    }
                }
                if (k + 1 < p_loc) {
                    lm_local.push_back(best_next);
                }
            }

            // Lineage-specific substitution rate correction (Saitou-Nei / Outgroup Centering)
            Vector r_i = D_loc.rowwise().mean();
            Vector rho_k = D_loc.colwise().mean().transpose();
            Scalar r_bar = r_i.mean();

            std::vector<Scalar> d_vals;
            d_vals.reserve(n * p_loc);
            for (size_t i = 0; i < n; ++i) {
                for (size_t k = 0; k < p_loc; ++k) {
                    Scalar d_corr = std::max<Scalar>(0.0, D_loc(i, k) - 0.5 * (r_i(i) + rho_k(k)) + 0.5 * r_bar);
                    D_loc(i, k) = d_corr;
                    if (d_corr > 1e-6) d_vals.push_back(d_corr);
                }
            }
            Scalar sigma_loc = 0.05;
            if (!d_vals.empty()) {
                std::sort(d_vals.begin(), d_vals.end());
                sigma_loc = d_vals[d_vals.size() / 2];
                sigma_loc = std::max<Scalar>(sigma_loc, 1e-4);
            }

            Matrix K_loc(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(p_loc));
            for (size_t i = 0; i < n; ++i) {
                for (size_t k = 0; k < p_loc; ++k) {
                    Scalar d = D_loc(i, k);
                    K_loc(i, k) = std::exp(-d / sigma_loc) / (d + eps0);
                }
            }

            Matrix K_pp(static_cast<Eigen::Index>(p_loc), static_cast<Eigen::Index>(p_loc));
            for (size_t a = 0; a < p_loc; ++a) {
                K_pp.row(static_cast<Eigen::Index>(a)) = K_loc.row(static_cast<Eigen::Index>(lm_local[a]));
            }
            K_pp = 0.5 * (K_pp + K_pp.transpose());

            Eigen::SelfAdjointEigenSolver<Matrix> eig_K(K_pp);
            Vector evals_K = eig_K.eigenvalues().cwiseMax(1e-6);
            Vector inv_sqrt_K = evals_K.array().rsqrt();
            Matrix Phi_S = K_loc * eig_K.eigenvectors() * inv_sqrt_K.asDiagonal();
            Eigen::Index p_cur = static_cast<Eigen::Index>(p_loc);

            // Matrix-free sub-degrees d = Phi_S * (Phi_S^T * 1_n)
            Vector s_Phi = Phi_S.colwise().sum().transpose();
            Vector d = (Phi_S * s_Phi).cwiseMax(1e-12);
            Vector inv_sqrt_d = d.array().rsqrt();

            // Normalized feature matrix A = D^{-1/2} * Phi_S [n x p_loc]
            Matrix A = inv_sqrt_d.asDiagonal() * Phi_S;

            // Gram matrix G = A^T * A [p_loc x p_loc]
            Matrix G = A.transpose() * A;
            G = 0.5 * (G + G.transpose());

            Eigen::SelfAdjointEigenSolver<Matrix> eig_G(G);
            const Matrix& evecs_G = eig_G.eigenvectors();

            // The second largest eigenvalue corresponds to column p_cur - 2 (eigenvalues sorted ascending)
            Eigen::Index sec_col = std::max<Eigen::Index>(0, p_cur - 2);
            Vector g_sec = evecs_G.col(sec_col);

            // Fiedler vector v2 = D^{-1} * Phi_S * g_sec
            Vector v2 = d.cwiseInverse().asDiagonal() * (Phi_S * g_sec);

            // Sort indices by v2
            std::vector<size_t> order(n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return v2(static_cast<Eigen::Index>(a)) < v2(static_cast<Eigen::Index>(b));
            });

            // Cheeger Conductance Sweep Cut in O(n * p_loc) time
            // For kernel graph C = Phi_S * Phi_S^T:
            // vol(A) = sum_{i in A} d_i
            // cut(A, \bar{A}) = vol(A) - ||sum_{i in A} Phi_S[i]||^2
            // Conductance phi(A) = cut(A, \bar{A}) / min(vol(A), vol_tot - vol(A))
            Scalar total_vol = d.sum();
            Scalar best_conductance = std::numeric_limits<Scalar>::infinity();
            size_t best_k = n / 2;
            size_t min_sz = std::max<size_t>(2, static_cast<size_t>(0.05 * static_cast<Scalar>(n)));

            Vector s_Phi_A = Vector::Zero(p_cur);
            Scalar cur_vol = 0.0;

            for (size_t k = 1; k < n; ++k) {
                size_t u = order[k - 1];
                s_Phi_A += Phi_S.row(static_cast<Eigen::Index>(u)).transpose();
                cur_vol += d(static_cast<Eigen::Index>(u));

                if (k >= min_sz && k + min_sz <= n) {
                    Scalar min_vol = std::min(cur_vol, total_vol - cur_vol);
                    if (min_vol > 1e-12) {
                        Scalar cut = std::max<Scalar>(0.0, cur_vol - s_Phi_A.squaredNorm());
                        Scalar cond = cut / min_vol;
                        if (cond < best_conductance) {
                            best_conductance = cond;
                            best_k = k;
                        }
                    }
                }
            }

            std::vector<size_t> indices_A, indices_B;
            indices_A.reserve(best_k);
            indices_B.reserve(n - best_k);
            for (size_t i = 0; i < best_k; ++i) indices_A.push_back(indices[order[i]]);
            for (size_t i = best_k; i < n; ++i) indices_B.push_back(indices[order[i]]);

            // Inter-clade representative distance via robust stratified sampling
            size_t s_a = std::min<size_t>(best_k, 8);
            size_t s_b = std::min<size_t>(n - best_k, 8);
            Scalar sum_d = 0.0;
            for (size_t ia = 0; ia < s_a; ++ia) {
                size_t rep_a = super_node_reps[indices_A[(ia * best_k) / s_a]];
                for (size_t ib = 0; ib < s_b; ++ib) {
                    size_t rep_b = super_node_reps[indices_B[(ib * (n - best_k)) / s_b]];
                    sum_d += compute_pair_distance(aln, rep_a, rep_b, pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights);
                }
            }
            Scalar inter_d = sum_d / static_cast<Scalar>(s_a * s_b);

            int32_t cl_A = recurse_hybrid(indices_A);
            int32_t cl_B = recurse_hybrid(indices_B);

            auto compute_mean_subtree_height = [&](int32_t nid) -> Scalar {
                if (nid < 0 || tree.nodes[nid].is_leaf) return 0.0;
                Scalar sum_h = 0.0;
                size_t count = 0;
                std::vector<std::pair<int32_t, Scalar>> st = {{nid, 0.0}};
                while (!st.empty()) {
                    auto [cur, d_cur] = st.back();
                    st.pop_back();
                    if (tree.nodes[cur].is_leaf) {
                        sum_h += d_cur;
                        count++;
                    } else {
                        for (int32_t ch : tree.nodes[cur].children) {
                            st.push_back({ch, d_cur + tree.nodes[ch].branch_length});
                        }
                    }
                }
                return (count > 0) ? (sum_h / static_cast<Scalar>(count)) : 0.0;
            };

            Scalar h_A = compute_mean_subtree_height(cl_A);
            Scalar h_B = compute_mean_subtree_height(cl_B);
            Scalar bridge_len = std::max<Scalar>(2e-6, inter_d - (h_A + h_B));

            tree.nodes[cl_A].branch_length = std::max<Scalar>(1e-6, 0.5 * bridge_len);
            tree.nodes[cl_B].branch_length = std::max<Scalar>(1e-6, 0.5 * bridge_len);

            int32_t spec_id = static_cast<int32_t>(tree.nodes.size());
            TreeNode spec_node;
            spec_node.id = spec_id;
            spec_node.name = "SpectralCut_" + std::to_string(internal_counter++);
            spec_node.is_leaf = false;
            spec_node.children = {cl_A, cl_B};
            tree.nodes.push_back(spec_node);
            tree.nodes[cl_A].parent_id = spec_id;
            tree.nodes[cl_B].parent_id = spec_id;
            return spec_id;
        };

        std::vector<size_t> root_indices(M);
        std::iota(root_indices.begin(), root_indices.end(), 0);
        int32_t root_id = recurse_hybrid(root_indices);

        tree.root_id = root_id;
        tree.nodes[root_id].parent_id = -1;
        tree.compute_post_order();
        return tree;
    }
    // Stage C: O(M * p) Ward's Minimum Variance Hierarchical Clustering via Nearest-Neighbor Chain
    static Tree infer_woodbury_ward(
        const Vector& inv_aug_deg,
        const Matrix& Z_coords,
        const Alignment& aln,
        const std::vector<size_t>& super_node_reps,
        const Vector4& pi,
        Scalar k1_denom,
        Scalar k2_denom,
        Scalar k3_denom,
        bool use_hamming,
        const std::vector<Scalar>& site_weights = {}
    ) {
        size_t M = super_node_reps.size();
        Eigen::Index p = Z_coords.cols();

        Tree tree;
        if (M == 1) {
            TreeNode root;
            root.id = 0;
            root.name = aln.taxon_names[super_node_reps[0]];
            root.is_leaf = true;
            root.branch_length = 0.0;
            tree.nodes.push_back(root);
            tree.root_id = 0;
            tree.leaf_name_to_id[root.name] = 0;
            tree.compute_post_order();
            return tree;
        }
        if (M == 2) {
            Scalar d01 = compute_pair_distance(
                aln, super_node_reps[0], super_node_reps[1], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights
            );
            for (size_t i = 0; i < 2; ++i) {
                TreeNode leaf;
                leaf.id = static_cast<int32_t>(i);
                leaf.name = aln.taxon_names[super_node_reps[i]];
                leaf.is_leaf = true;
                leaf.branch_length = std::max<Scalar>(1e-6, 0.5 * d01);
                leaf.parent_id = 2;
                tree.nodes.push_back(leaf);
                tree.leaf_name_to_id[leaf.name] = leaf.id;
            }
            TreeNode root;
            root.id = 2;
            root.name = "Root";
            root.is_leaf = false;
            root.children = {0, 1};
            tree.nodes.push_back(root);
            tree.root_id = 2;
            tree.compute_post_order();
            return tree;
        }

        tree.nodes.reserve(2 * M - 1);
        for (size_t i = 0; i < M; ++i) {
            TreeNode leaf;
            leaf.id = static_cast<int32_t>(i);
            leaf.name = aln.taxon_names[super_node_reps[i]];
            leaf.is_leaf = true;
            tree.nodes.push_back(leaf);
            tree.leaf_name_to_id[leaf.name] = leaf.id;
        }

        size_t max_clusters = 2 * M;
        Matrix c_Z = Matrix::Zero(static_cast<Eigen::Index>(max_clusters), p);
        std::vector<Scalar> c_weight(max_clusters, 1.0);
        std::vector<Scalar> c_z_sq(max_clusters, 0.0);
        std::vector<size_t> c_rep_taxon(max_clusters, 0);
        std::vector<Scalar> c_internal_depth(max_clusters, 0.0);
        std::vector<bool> is_active(max_clusters, false);
        std::vector<int32_t> active_list;
        active_list.reserve(M);
        std::vector<size_t> pos_in_active(max_clusters, 0);

        std::vector<Scalar> c_psi(max_clusters, 0.0);
        Vector s_z = Vector::Zero(p);
        for (size_t i = 0; i < M; ++i) {
            c_Z.row(static_cast<Eigen::Index>(i)) = Z_coords.row(static_cast<Eigen::Index>(i));
            c_z_sq[i] = c_Z.row(static_cast<Eigen::Index>(i)).squaredNorm();
            c_weight[i] = 1.0;
            c_rep_taxon[i] = super_node_reps[i];
            c_internal_depth[i] = 0.0;
            Scalar inv_d = (static_cast<Eigen::Index>(i) < inv_aug_deg.size()) ? inv_aug_deg(static_cast<Eigen::Index>(i)) : 0.0;
            c_psi[i] = inv_d + c_z_sq[i];
            s_z += c_Z.row(static_cast<Eigen::Index>(i)).transpose();
            is_active[i] = true;
            pos_in_active[i] = active_list.size();
            active_list.push_back(static_cast<int32_t>(i));
        }

        auto remove_active = [&](int32_t id) {
            size_t idx = pos_in_active[id];
            size_t last_idx = active_list.size() - 1;
            int32_t last_id = active_list[last_idx];
            active_list[idx] = last_id;
            pos_in_active[last_id] = idx;
            active_list.pop_back();
            is_active[id] = false;
        };

        auto add_active = [&](int32_t id) {
            pos_in_active[id] = active_list.size();
            active_list.push_back(id);
            is_active[id] = true;
        };

        auto find_nn = [&](int32_t u) -> int32_t {
            Scalar best_score = std::numeric_limits<Scalar>::infinity();
            int32_t best_v = -1;
            size_t m = active_list.size();
            Scalar m_minus_2 = static_cast<Scalar>(m - 2);
            Vector w_u = s_z - m_minus_2 * c_Z.row(u).transpose();

            if (m > 256) {
                #pragma omp parallel
                {
                    Scalar local_best_score = std::numeric_limits<Scalar>::infinity();
                    int32_t local_best_v = -1;

                    #pragma omp for nowait
                    for (size_t i = 0; i < m; ++i) {
                        int32_t v = active_list[i];
                        if (v == u) continue;
                        Scalar score = c_Z.row(v).dot(w_u.transpose()) - c_psi[v];
                        if (score < local_best_score || (score == local_best_score && v < local_best_v)) {
                            local_best_score = score;
                            local_best_v = v;
                        }
                    }

                    #pragma omp critical
                    {
                        if (local_best_score < best_score || (local_best_score == best_score && local_best_v < best_v)) {
                            best_score = local_best_score;
                            best_v = local_best_v;
                        }
                    }
                }
            } else {
                for (size_t i = 0; i < m; ++i) {
                    int32_t v = active_list[i];
                    if (v == u) continue;
                    Scalar score = c_Z.row(v).dot(w_u.transpose()) - c_psi[v];
                    if (score < best_score || (score == best_score && v < best_v)) {
                        best_score = score;
                        best_v = v;
                    }
                }
            }
            return best_v;
        };

        int32_t next_cluster_id = static_cast<int32_t>(M);
        int internal_counter = 1;
        std::vector<int32_t> chain;
        chain.reserve(M + 1);

        auto do_merge = [&](int32_t a, int32_t b) {
            int32_t w_id = next_cluster_id++;
            Scalar wu = c_weight[a];
            Scalar wv = c_weight[b];
            Scalar total_w = wu + wv;
            c_weight[w_id] = total_w;

            // Exact Saitou-Nei reduction in Woodbury space:
            // z_w = 0.5 * (z_a + z_b)
            c_Z.row(w_id) = 0.5 * (c_Z.row(a) + c_Z.row(b));
            // psi_w = z_a . z_b
            c_psi[w_id] = c_Z.row(a).dot(c_Z.row(b));

            // Update global sum: s_z <- s_z - 0.5 * (z_a + z_b)
            s_z.noalias() -= 0.5 * (c_Z.row(a) + c_Z.row(b)).transpose();

            // Compute evolutionary branch lengths between cluster representatives
            Scalar d_uv = compute_pair_distance(
                aln, c_rep_taxon[a], c_rep_taxon[b], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights
            );
            size_t m = active_list.size();
            Scalar D_ab = c_psi[a] + c_psi[b] - 2.0 * c_psi[w_id];
            Vector s_z_orig = s_z + 0.5 * (c_Z.row(a) + c_Z.row(b)).transpose();
            Scalar r_diff = static_cast<Scalar>(m) * (c_psi[a] - c_psi[b]) - 2.0 * (c_Z.row(a) - c_Z.row(b)).dot(s_z_orig.transpose());
            Scalar ratio_a = 0.5;
            if (m > 2 && D_ab > 1e-12) {
                ratio_a = std::clamp<Scalar>(0.5 + r_diff / (2.0 * static_cast<Scalar>(m - 2) * D_ab), 0.05, 0.95);
            }

            Scalar net_d = std::max<Scalar>(0.0, d_uv - c_internal_depth[a] - c_internal_depth[b]);
            Scalar len_a = std::max<Scalar>(1e-6, ratio_a * net_d);
            Scalar len_b = std::max<Scalar>(1e-6, (1.0 - ratio_a) * net_d);

            tree.nodes[a].branch_length = len_a;
            tree.nodes[b].branch_length = len_b;

            TreeNode w_node;
            w_node.id = w_id;
            w_node.name = "Node" + std::to_string(internal_counter++);
            w_node.is_leaf = false;
            w_node.children = {a, b};
            tree.nodes.push_back(w_node);
            tree.nodes[a].parent_id = w_id;
            tree.nodes[b].parent_id = w_id;

            c_rep_taxon[w_id] = (wu >= wv) ? c_rep_taxon[a] : c_rep_taxon[b];
            c_internal_depth[w_id] = 0.5 * (c_internal_depth[a] + len_a + c_internal_depth[b] + len_b);

            remove_active(a);
            remove_active(b);
            add_active(w_id);
        };

        while (active_list.size() > 3) {
            if (chain.empty()) {
                chain.push_back(active_list[0]);
            }

            int32_t u = chain.back();
            int32_t v = find_nn(u);
            if (v < 0) break;

            if (chain.size() >= 2 && chain[chain.size() - 2] == v) {
                // Mutual Nearest Neighbors! Pop both and merge
                chain.pop_back();
                chain.pop_back();
                do_merge(u, v);
            } else {
                auto it = std::find(chain.begin(), chain.end(), v);
                if (it != chain.end()) {
                    // Cycle detected: pop back to it and merge (u, v)
                    while (!chain.empty() && chain.back() != v) {
                        chain.pop_back();
                    }
                    if (!chain.empty()) chain.pop_back();
                    do_merge(u, v);
                } else {
                    chain.push_back(v);
                }
            }
        }

        // Join final 3 clusters into trifurcating root
        if (active_list.size() == 3) {
            int32_t n0 = active_list[0], n1 = active_list[1], n2 = active_list[2];
            Scalar d01 = std::max<Scalar>(0.0, compute_pair_distance(aln, c_rep_taxon[n0], c_rep_taxon[n1], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights) - c_internal_depth[n0] - c_internal_depth[n1]);
            Scalar d02 = std::max<Scalar>(0.0, compute_pair_distance(aln, c_rep_taxon[n0], c_rep_taxon[n2], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights) - c_internal_depth[n0] - c_internal_depth[n2]);
            Scalar d12 = std::max<Scalar>(0.0, compute_pair_distance(aln, c_rep_taxon[n1], c_rep_taxon[n2], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights) - c_internal_depth[n1] - c_internal_depth[n2]);

            tree.nodes[n0].branch_length = std::max<Scalar>(1e-6, 0.5 * (d01 + d02 - d12));
            tree.nodes[n1].branch_length = std::max<Scalar>(1e-6, 0.5 * (d01 + d12 - d02));
            tree.nodes[n2].branch_length = std::max<Scalar>(1e-6, 0.5 * (d02 + d12 - d01));

            int32_t root_id = static_cast<int32_t>(tree.nodes.size());
            TreeNode root;
            root.id = root_id;
            root.name = "Root";
            root.is_leaf = false;
            root.children = {n0, n1, n2};
            tree.nodes.push_back(root);
            tree.nodes[n0].parent_id = root_id;
            tree.nodes[n1].parent_id = root_id;
            tree.nodes[n2].parent_id = root_id;
            tree.root_id = root_id;
        } else if (active_list.size() == 2) {
            int32_t n0 = active_list[0], n1 = active_list[1];
            Scalar d01 = std::max<Scalar>(0.0, compute_pair_distance(aln, c_rep_taxon[n0], c_rep_taxon[n1], pi, k1_denom, k2_denom, k3_denom, use_hamming, site_weights) - c_internal_depth[n0] - c_internal_depth[n1]);
            tree.nodes[n0].branch_length = std::max<Scalar>(1e-6, 0.5 * d01);
            tree.nodes[n1].branch_length = std::max<Scalar>(1e-6, 0.5 * d01);
            int32_t root_id = static_cast<int32_t>(tree.nodes.size());
            TreeNode root;
            root.id = root_id;
            root.name = "Root";
            root.is_leaf = false;
            root.children = {n0, n1};
            tree.nodes.push_back(root);
            tree.nodes[n0].parent_id = root_id;
            tree.nodes[n1].parent_id = root_id;
            tree.root_id = root_id;
        } else if (active_list.size() == 1) {
            tree.root_id = active_list[0];
        }

        tree.compute_post_order();
        return tree;
    }

    // Build initial k-NN lists in O(M * p * log M) using recursive random-projection tree leaves
    static void build_rp_knn(
        const Matrix& Z,
        const std::vector<int32_t>& active,
        size_t k_knn,
        std::vector<std::vector<int32_t>>& knn_list
    ) {
        size_t M = active.size();
        if (M <= 96) {
            // Direct leaf bucket for small M
            for (size_t i = 0; i < M; ++i) {
                int32_t ca = active[i];
                std::vector<std::pair<Scalar, int32_t>> dists;
                dists.reserve(M - 1);
                for (size_t j = 0; j < M; ++j) {
                    if (i == j) continue;
                    int32_t cb = active[j];
                    Scalar d2 = (Z.row(ca) - Z.row(cb)).squaredNorm();
                    dists.push_back({d2, cb});
                }
                size_t k = std::min(k_knn, dists.size());
                std::partial_sort(dists.begin(), dists.begin() + static_cast<std::ptrdiff_t>(k), dists.end());
                knn_list[ca].resize(k);
                for (size_t r = 0; r < k; ++r) knn_list[ca][r] = dists[r].second;
            }
            return;
        }

        // Recursive coordinate-variance split into buckets of size <= 64
        auto recurse = [&](auto& self, const std::vector<int32_t>& subset, int depth) -> void {
            if (subset.size() <= 64 || depth >= 12) {
                for (size_t i = 0; i < subset.size(); ++i) {
                    int32_t ca = subset[i];
                    for (size_t j = 0; j < subset.size(); ++j) {
                        if (i != j) knn_list[ca].push_back(subset[j]);
                    }
                }
                return;
            }
            // Pick coordinate with largest spread among first 16 dims
            Eigen::Index best_col = depth % Z.cols();
            std::vector<std::pair<Scalar, int32_t>> proj(subset.size());
            for (size_t i = 0; i < subset.size(); ++i) {
                proj[i] = {Z(subset[i], best_col), subset[i]};
            }
            size_t mid = proj.size() / 2;
            std::nth_element(proj.begin(), proj.begin() + static_cast<std::ptrdiff_t>(mid), proj.end());
            // Overlapping split (15% overlap) to preserve boundary neighbors
            size_t overlap = std::max<size_t>(4, subset.size() / 8);
            size_t left_end = std::min(proj.size(), mid + overlap);
            size_t right_start = (mid > overlap) ? (mid - overlap) : 0;

            std::vector<int32_t> left_sub, right_sub;
            left_sub.reserve(left_end);
            right_sub.reserve(proj.size() - right_start);
            for (size_t i = 0; i < left_end; ++i) left_sub.push_back(proj[i].second);
            for (size_t i = right_start; i < proj.size(); ++i) right_sub.push_back(proj[i].second);

            self(self, left_sub, depth + 1);
            self(self, right_sub, depth + 1);
        };

        recurse(recurse, active, 0);

        // Deduplicate and keep top k_knn closest in Z-space per node
        for (int32_t ca : active) {
            auto& vec = knn_list[ca];
            std::sort(vec.begin(), vec.end());
            vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
            if (vec.size() > k_knn) {
                std::partial_sort(
                    vec.begin(),
                    vec.begin() + static_cast<std::ptrdiff_t>(k_knn),
                    vec.end(),
                    [&](int32_t a, int32_t b) {
                        return (Z.row(ca) - Z.row(a)).squaredNorm() <
                               (Z.row(ca) - Z.row(b)).squaredNorm();
                    }
                );
                vec.resize(k_knn);
            }
        }
    }

    // Expand super-nodes with |C_m| > 1 into terminal busbars (for collapsed_tree)
    // and resolved zero-length binary ladders (for backbone_tree)
    static void expand_super_nodes_into_trees(
        OhmicNetwork& net,
        const Tree& super_tree,
        const Alignment& aln,
        const std::vector<size_t>& super_node_reps,
        const std::vector<std::vector<size_t>>& super_node_members
    ) {
        size_t M = super_node_reps.size();
        size_t N = aln.num_taxa;
        if (M == N) {
            net.backbone_tree = super_tree;
            net.collapsed_tree = super_tree;
            return;
        }

        std::unordered_map<std::string, size_t> rep_name_to_m;
        for (size_t m = 0; m < M; ++m) {
            rep_name_to_m[aln.taxon_names[super_node_reps[m]]] = m;
        }

        // 1. Build collapsed_tree with native k-way multifurcating busbars for |C_m| > 1
        Tree col = super_tree;
        for (size_t nid = 0; nid < super_tree.nodes.size(); ++nid) {
            if (!col.nodes[nid].is_leaf) continue;
            auto it = rep_name_to_m.find(col.nodes[nid].name);
            if (it == rep_name_to_m.end()) continue;
            size_t m = it->second;
            const auto& members = super_node_members[m];
            if (members.size() <= 1) continue;

            // Convert leaf node `nid` into an internal busbar hub with |members| leaf children
            col.leaf_name_to_id.erase(col.nodes[nid].name);
            col.nodes[nid].is_leaf = false;
            col.nodes[nid].name = "Busbar_" + std::to_string(m + 1);
            col.nodes[nid].children.clear();

            for (size_t tax_idx : members) {
                int32_t child_id = static_cast<int32_t>(col.nodes.size());
                TreeNode leaf;
                leaf.id = child_id;
                leaf.name = aln.taxon_names[tax_idx];
                leaf.is_leaf = true;
                leaf.parent_id = static_cast<int32_t>(nid);
                leaf.branch_length = 1e-6;
                col.nodes.push_back(leaf);
                col.nodes[nid].children.push_back(child_id);
                col.leaf_name_to_id[leaf.name] = child_id;
            }
            net.num_collapsed_polytomies += (members.size() - 2);
        }
        col.compute_post_order();
        net.collapsed_tree = col;

        // 2. Build resolved binary backbone_tree by binary-resolving each super-node cluster
        Tree bin = super_tree;
        for (size_t nid = 0; nid < super_tree.nodes.size(); ++nid) {
            if (!bin.nodes[nid].is_leaf) continue;
            auto it = rep_name_to_m.find(bin.nodes[nid].name);
            if (it == rep_name_to_m.end()) continue;
            size_t m = it->second;
            const auto& members = super_node_members[m];
            if (members.size() <= 1) continue;

            bin.leaf_name_to_id.erase(bin.nodes[nid].name);
            int32_t cur_hub = static_cast<int32_t>(nid);
            for (size_t k = 0; k + 1 < members.size(); ++k) {
                bin.nodes[cur_hub].is_leaf = false;
                bin.nodes[cur_hub].name = "BusBin_" + std::to_string(m + 1) + "_" + std::to_string(k + 1);
                bin.nodes[cur_hub].children.clear();

                // Left child: leaf members[k]
                int32_t leaf_id = static_cast<int32_t>(bin.nodes.size());
                TreeNode leaf;
                leaf.id = leaf_id;
                leaf.name = aln.taxon_names[members[k]];
                leaf.is_leaf = true;
                leaf.parent_id = cur_hub;
                leaf.branch_length = 1e-6;
                bin.nodes.push_back(leaf);
                bin.leaf_name_to_id[leaf.name] = leaf_id;

                // Right child: either final leaf (if k + 2 == members.size()) or next internal node
                int32_t right_id = static_cast<int32_t>(bin.nodes.size());
                TreeNode right_node;
                right_node.id = right_id;
                right_node.parent_id = cur_hub;
                right_node.branch_length = 1e-6;
                if (k + 2 == members.size()) {
                    right_node.is_leaf = true;
                    right_node.name = aln.taxon_names[members[k + 1]];
                    bin.nodes.push_back(right_node);
                    bin.leaf_name_to_id[right_node.name] = right_id;
                } else {
                    right_node.is_leaf = false;
                    right_node.name = "BusBin_" + std::to_string(m + 1) + "_" + std::to_string(k + 2);
                    bin.nodes.push_back(right_node);
                }

                bin.nodes[cur_hub].children = {leaf_id, right_id};
                cur_hub = right_id;
            }
        }
        bin.compute_post_order();
        net.backbone_tree = bin;
    }



    // Score quartet conductances and build node-disjoint 4-way Wheatstone bridges + Ohmic Polytomies
    static void build_bridges_and_collapsed_tree(
        OhmicNetwork& net,
        const Alignment& aln,
        const OhmicSettings& settings,
        Scalar eps_L,
        bool is_dense_mode
    ) {
        if (is_dense_mode) {
            net.collapsed_tree = net.backbone_tree;
        }
        Tree& tree = net.collapsed_tree;
        size_t N = aln.num_taxa;

        auto clade_taxa = collect_subtree_taxa(tree, aln);
        std::vector<int> all_taxa(N);
        std::iota(all_taxa.begin(), all_taxa.end(), 0);

        Scalar omega_scale = 1.0;
        if (is_dense_mode && net.Omega.rows() == static_cast<Eigen::Index>(N)) {
            Scalar sum_d = net.D_tn93.sum();
            Scalar sum_om = net.Omega.sum();
            if (sum_om > 1e-12 && sum_d > 1e-12) {
                omega_scale = sum_d / sum_om;
            }
        }

        // Helper to compute mean pairwise distance/resistance between two taxon subsets C1, C2
        auto mean_inter_dist = [&](const std::vector<int>& C1, const std::vector<int>& C2) -> Scalar {
            if (C1.empty() || C2.empty()) return 0.0;
            if (is_dense_mode) {
                Scalar sum = 0.0;
                for (int i : C1) {
                    for (int j : C2) {
                        if (net.used_hyphaeon_embeddings) {
                            sum += 0.85 * net.D_tn93(i, j) + 0.15 * omega_scale * net.Omega(i, j);
                        } else {
                            sum += net.D_tn93(i, j);
                        }
                    }
                }
                return sum / static_cast<Scalar>(C1.size() * C2.size());
            } else {
                // Sample up to 8 representatives per clade in ScalableWoodbury mode for O(1) evaluation
                size_t s1 = std::min<size_t>(C1.size(), 8);
                size_t s2 = std::min<size_t>(C2.size(), 8);
                Scalar sum = 0.0;
                Vector4 dummy_pi = Vector4::Constant(0.25);
                for (size_t a = 0; a < s1; ++a) {
                    size_t idx1 = C1[(a * C1.size()) / s1];
                    for (size_t b = 0; b < s2; ++b) {
                        size_t idx2 = C2[(b * C2.size()) / s2];
                        sum += compute_pair_distance(aln, idx1, idx2, dummy_pi, 0.25, 0.25, 0.5, true);
                    }
                }
                return sum / static_cast<Scalar>(s1 * s2);
            }
        };

        struct CandidateBridge {
            int32_t u_id = INVALID_INDEX;
            int32_t v_id = INVALID_INDEX;
            int32_t id_A = INVALID_INDEX;
            int32_t id_B = INVALID_INDEX;
            int32_t id_C = INVALID_INDEX;
            std::array<Scalar, 3> weights{1.0, 0.0, 0.0};
            std::array<Scalar, 3> lengths{0.01, 0.01, 0.01};
            Scalar ohmic_support = 1.0;
            Scalar shunt_stress = 0.0;
            bool is_polytomy_noise = false;
        };

        // Precompute complement representative taxa for all internal nodes in O(N) time
        size_t num_nodes = tree.nodes.size();
        std::vector<std::vector<int>> comp_taxa(num_nodes);
        std::vector<int32_t> bfs_queue;
        bfs_queue.push_back(tree.root_id);
        size_t q_head = 0;
        while (q_head < bfs_queue.size()) {
            int32_t pid = bfs_queue[q_head++];
            const auto& p_node = tree.nodes[pid];
            size_t num_ch = p_node.children.size();
            for (size_t i = 0; i < num_ch; ++i) {
                int32_t cid = p_node.children[i];
                if (!tree.nodes[cid].is_leaf) {
                    bfs_queue.push_back(cid);
                }
                auto& c_comp = comp_taxa[cid];
                // Collect representatives from siblings
                for (size_t j = 0; j < num_ch; ++j) {
                    if (i == j) continue;
                    int32_t sib = p_node.children[j];
                    const auto& s_taxa = clade_taxa[sib];
                    size_t take = is_dense_mode ? s_taxa.size() : std::min<size_t>(s_taxa.size(), 8);
                    for (size_t k = 0; k < take; ++k) {
                        c_comp.push_back(s_taxa[(k * s_taxa.size()) / take]);
                    }
                }
                // Collect representatives from parent's complement
                const auto& p_comp = comp_taxa[pid];
                size_t take_p = is_dense_mode ? p_comp.size() : std::min<size_t>(p_comp.size(), 8);
                for (size_t k = 0; k < take_p; ++k) {
                    c_comp.push_back(p_comp[(k * p_comp.size()) / take_p]);
                }
            }
        }

        std::vector<CandidateBridge> candidates;

        for (const auto& v_node : tree.nodes) {
            if (v_node.is_leaf || v_node.id == tree.root_id) continue;
            if (v_node.children.size() != 2) continue;
            // Short-circuit: large internal branches cannot qualify as unresolved Wheatstone bridges or Ohmic polytomies
            if (!is_dense_mode && v_node.branch_length > 0.05) continue;

            int32_t v_id = v_node.id;
            int32_t u_id = v_node.parent_id;
            if (u_id == INVALID_INDEX) continue;
            const auto& u_node = tree.nodes[u_id];

            int32_t id_A = v_node.children[0];
            int32_t id_B = v_node.children[1];

            std::vector<int32_t> u_other_children;
            for (int32_t cid : u_node.children) {
                if (cid != v_id) u_other_children.push_back(cid);
            }
            if (u_other_children.empty()) continue;

            int32_t best_C = INVALID_INDEX;
            CandidateBridge best_cand;
            Scalar best_stress = -1.0;

            for (int32_t cand_C : u_other_children) {
                const auto& taxa_A = clade_taxa[id_A];
                const auto& taxa_B = clade_taxa[id_B];
                const auto& taxa_C = clade_taxa[cand_C];
                const std::vector<int>* taxa_D_ptr = nullptr;

                if (u_id == tree.root_id) {
                    for (int32_t r_cid : u_node.children) {
                        if (r_cid != v_id && r_cid != cand_C) {
                            taxa_D_ptr = &clade_taxa[r_cid];
                            break;
                        }
                    }
                } else {
                    taxa_D_ptr = &comp_taxa[u_id];
                }

                if (!taxa_D_ptr || taxa_D_ptr->empty() || taxa_A.empty() || taxa_B.empty() || taxa_C.empty()) continue;
                const auto& taxa_D = *taxa_D_ptr;

                // Evaluate 4x4 inter-cluster distances
                Scalar d_AB = mean_inter_dist(taxa_A, taxa_B);
                Scalar d_CD = mean_inter_dist(taxa_C, taxa_D);
                Scalar d_AC = mean_inter_dist(taxa_A, taxa_C);
                Scalar d_BD = mean_inter_dist(taxa_B, taxa_D);
                Scalar d_AD = mean_inter_dist(taxa_A, taxa_D);
                Scalar d_BC = mean_inter_dist(taxa_B, taxa_C);

                Scalar q1 = d_AB + d_CD; // Primary: AB | CD
                Scalar q2 = d_AC + d_BD; // Secondary: AC | BD
                Scalar q3 = d_AD + d_BC; // Tertiary: AD | BC (or BC | AD)

                // Order A and B so that AC | BD has at least as much support as BC | AD (q2 <= q3)
                int32_t cur_A = id_A;
                int32_t cur_B = id_B;
                if (q3 < q2) {
                    std::swap(q2, q3);
                    std::swap(cur_A, cur_B);
                }

                Scalar q_max = std::max({q1, q2, q3});
                Scalar raw_d1 = std::max<Scalar>(0.0, q_max - q1);
                Scalar raw_d2 = std::max<Scalar>(0.0, q3 - q2);

                // Stage D check for unresolved equipotential noise (star tree: all 3 quartets within 2 * eps_L)
                if (v_node.branch_length <= eps_L && raw_d1 <= 2.0 * eps_L && raw_d2 <= 2.0 * eps_L) {
                    CandidateBridge poly_cand;
                    poly_cand.u_id = u_id;
                    poly_cand.v_id = v_id;
                    poly_cand.is_polytomy_noise = true;
                    best_cand = poly_cand;
                    best_C = cand_C;
                    break;
                }

                Scalar base_len = std::max<Scalar>(v_node.branch_length, eps_L);
                Scalar scale_eps = std::max<Scalar>(eps_L, 1e-5);
                Scalar delta1 = std::max<Scalar>(2.0 * base_len, raw_d1);
                Scalar delta2 = raw_d2 + base_len * std::exp(-std::max<Scalar>(0.0, q2 - q1) / scale_eps);
                Scalar delta3 = settings.include_tertiary
                    ? (base_len * std::exp(-std::max<Scalar>(0.0, q3 - q1) / scale_eps))
                    : 0.0;

                Scalar G1 = delta1;
                Scalar G2 = delta2;
                Scalar G3 = delta3;
                Scalar G_sum = G1 + G2 + G3;
                if (G_sum <= 1e-12) continue;

                Scalar tau = G2 / G1;
                if (tau > best_stress) {
                    best_stress = tau;
                    best_C = cand_C;
                    best_cand.u_id = u_id;
                    best_cand.v_id = v_id;
                    best_cand.id_A = cur_A;
                    best_cand.id_B = cur_B;
                    best_cand.id_C = cand_C;
                    best_cand.weights = {G1 / G_sum, G2 / G_sum, G3 / G_sum};
                    Scalar t1 = v_node.branch_length;
                    Scalar t2 = std::max<Scalar>(1e-5, 0.5 * delta2);
                    Scalar t3 = std::max<Scalar>(1e-5, 0.5 * delta3);
                    best_cand.lengths = {t1, t2, t3};
                    best_cand.ohmic_support = G1 / G_sum;
                    best_cand.shunt_stress = tau;
                    best_cand.is_polytomy_noise = false;
                }
            }

            if (best_C != INVALID_INDEX) {
                candidates.push_back(best_cand);
            }
        }

        // Sort bridge candidates by descending shunt_stress
        std::sort(candidates.begin(), candidates.end(), [](const CandidateBridge& a, const CandidateBridge& b) {
            return a.shunt_stress > b.shunt_stress;
        });

        std::unordered_set<int32_t> used_nodes;
        std::vector<CandidateBridge> selected_bridges;
        std::vector<CandidateBridge> selected_polytomies;

        for (const auto& cand : candidates) {
            if (cand.is_polytomy_noise) {
                if (settings.collapse_ohmic_polytomies &&
                    used_nodes.find(cand.v_id) == used_nodes.end()) {
                    selected_polytomies.push_back(cand);
                    used_nodes.insert(cand.v_id);
                }
                continue;
            }
            if (cand.shunt_stress < settings.min_shunt_stress) continue;
            if (selected_bridges.size() >= settings.max_bridges) continue;

            // Node-disjointness constraint: u_id, v_id, id_A, id_B, id_C must not be a hub or bridge of another bridge
            // Specifically, u_id must have exactly 2 children (or if u_id is 3-child root, we insert an intermediate 2-child hub)
            if (used_nodes.count(cand.u_id) || used_nodes.count(cand.v_id) ||
                used_nodes.count(cand.id_A) || used_nodes.count(cand.id_B) || used_nodes.count(cand.id_C)) {
                continue;
            }
            // If u_id has > 2 children (e.g. trifurcating root), only allow if u_id has 2 children or we create a dedicated hub
            selected_bridges.push_back(cand);
            used_nodes.insert(cand.u_id);
            used_nodes.insert(cand.v_id);
        }

        // Apply selected Wheatstone bridges to `collapsed_tree`
        for (const auto& sb : selected_bridges) {
            int32_t hub_id = sb.u_id;
            int32_t v_id = sb.v_id;

            // If hub_id has > 2 children (e.g. trifurcating root with {v_id, id_C, id_D}),
            // reuse v_id as the 3-child Wheatstone hub under root so root keeps {v_id, id_D}!
            if (tree.nodes[hub_id].children.size() > 2) {
                // Make v_id the 3-child hub containing {id_A, id_B, id_C}, and remove id_C from hub_id.children
                auto& root_ch = tree.nodes[hub_id].children;
                root_ch.erase(std::remove(root_ch.begin(), root_ch.end(), sb.id_C), root_ch.end());
                tree.nodes[v_id].children = {sb.id_A, sb.id_B, sb.id_C};
                tree.nodes[sb.id_A].parent_id = v_id;
                tree.nodes[sb.id_B].parent_id = v_id;
                tree.nodes[sb.id_C].parent_id = v_id;
                // The branch from root (hub_id) to v_id has zero additional length beyond D's own branch
                tree.nodes[v_id].branch_length = 0.0;

                WheatstoneBridge wb;
                wb.hub_node_id = v_id;
                wb.bridge_node_id = v_id;
                wb.children = {sb.id_A, sb.id_B, sb.id_C};
                wb.weights = sb.weights;
                wb.branch_lengths = sb.lengths;
                wb.ohmic_support = sb.ohmic_support;
                wb.shunt_stress = sb.shunt_stress;
                net.bridges.push_back(wb);
            } else {
                // Standard binary parent u_id with children {v_id, id_C}:
                // Collapse v_id into u_id so u_id becomes the 3-child hub {id_A, id_B, id_C}
                tree.nodes[hub_id].children = {sb.id_A, sb.id_B, sb.id_C};
                tree.nodes[sb.id_A].parent_id = hub_id;
                tree.nodes[sb.id_B].parent_id = hub_id;
                tree.nodes[sb.id_C].parent_id = hub_id;
                tree.nodes[v_id].children.clear();
                tree.nodes[v_id].parent_id = INVALID_INDEX;

                WheatstoneBridge wb;
                wb.hub_node_id = hub_id;
                wb.bridge_node_id = v_id;
                wb.children = {sb.id_A, sb.id_B, sb.id_C};
                wb.weights = sb.weights;
                wb.branch_lengths = sb.lengths;
                wb.ohmic_support = sb.ohmic_support;
                wb.shunt_stress = sb.shunt_stress;
                net.bridges.push_back(wb);
            }
        }

        // Apply selected Ohmic Polytomy collapses (Stage D Case 2)
        for (const auto& poly : selected_polytomies) {
            int32_t v_id = poly.v_id;
            int32_t u_id = tree.nodes[v_id].parent_id;
            if (u_id < 0 || tree.nodes[v_id].children.empty()) continue;
            auto& u_ch = tree.nodes[u_id].children;
            auto it = std::find(u_ch.begin(), u_ch.end(), v_id);
            if (it == u_ch.end()) continue;
            u_ch.erase(it);
            for (int32_t cid : tree.nodes[v_id].children) {
                u_ch.push_back(cid);
                tree.nodes[cid].parent_id = u_id;
            }
            tree.nodes[v_id].children.clear();
            tree.nodes[v_id].parent_id = INVALID_INDEX;
            net.num_collapsed_polytomies++;
        }

        tree.compute_post_order();
        net.node_to_bridge_idx.assign(tree.num_nodes(), -1);
        for (size_t k = 0; k < net.bridges.size(); ++k) {
            int32_t hid = net.bridges[k].hub_node_id;
            if (hid >= 0 && static_cast<size_t>(hid) < net.node_to_bridge_idx.size()) {
                net.node_to_bridge_idx[hid] = static_cast<int32_t>(k);
            }
        }
        tree.bridges = net.bridges;
        tree.node_to_bridge_idx = net.node_to_bridge_idx;
    }
};

} // namespace hyphy::core
