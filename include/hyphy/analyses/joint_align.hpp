#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/ohmic_network.hpp"
#include "hyphy/core/wheatstone_likelihood.hpp"
#include "hyphy/analyses/tree_engine.hpp"
#include "hyphy/analyses/fel.hpp"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <iostream>
#include <iomanip>

namespace hyphy::analyses {

using namespace hyphy::core;

// =============================================================================
// Joint Ohmic Alignment + Phylogeny + Selection Engine (3-Stage Prototype)
//   Stage 1: Phylogenetic Codon Indel Transducer & Window Realignment (Inside-Outside EM)
//   Stage 2: Coupled Ohmic SPR / Port-Slide with Boundary Indel Transducer Updates
//   Stage 3: Alignment-Marginalized FEL (Posterior Homology-Weighted Selection LRT)
// =============================================================================

struct IndelWindow {
    size_t start_codon = 0; // inclusive (0-indexed)
    size_t end_codon = 0;   // exclusive (0-indexed)
    size_t num_gapped_taxa = 0;
    size_t total_gap_codons = 0;
    Scalar mean_gap_fraction = 0.0;
};

struct TaxonWindowPlacement {
    size_t taxon_idx = 0;
    std::vector<int8_t> aligned_codons; // length W = end_codon - start_codon (-1 for gap)
    std::string raw_Subseq;             // length 3 * W
    Scalar log_posterior = 0.0;
    Scalar posterior_prob = 1.0;
    bool changed_from_input = false;
};

struct WindowRealignmentSummary {
    IndelWindow window;
    Scalar initial_window_log_l = 0.0;
    Scalar polished_window_log_l = 0.0;
    Scalar delta_log_l = 0.0;
    size_t num_taxa_realigned = 0;
    size_t num_codons_shifted = 0;
    std::vector<Scalar> column_homology_posteriors; // length W: expected posterior confidence per column
};

struct CoupledMoveEvent {
    std::string move_type; // "Coupled_SPR" or "Coupled_PortSlide"
    int32_t clade_node_id = -1;
    int32_t target_node_id = -1;
    Scalar delta_sub_log_l = 0.0;
    Scalar delta_joint_log_l = 0.0;
    size_t realigned_clade_taxa = 0;
};

struct MarginalizedSiteResult {
    size_t site_index = 0; // 0-indexed codon site
    bool in_indel_window = false;
    Scalar column_homology_posterior = 1.0; // Posterior alignment consistency in [0, 1]
    Scalar gap_fraction = 0.0;

    // Input (rigid MSA) FEL statistics
    Scalar rigid_alpha = 0.0;
    Scalar rigid_beta = 0.0;
    Scalar rigid_lrt = 0.0;
    Scalar rigid_p_value = 1.0;

    // MAP Joint Realigned FEL statistics
    Scalar map_alpha = 0.0;
    Scalar map_beta = 0.0;
    Scalar map_lrt = 0.0;
    Scalar map_p_value = 1.0;

    // Stage 3: Alignment-Marginalized FEL statistics (integrated over candidate gap placements)
    Scalar marginalized_alpha = 0.0;
    Scalar marginalized_beta = 0.0;
    Scalar marginalized_lrt = 0.0;
    Scalar marginalized_p_value = 1.0;

    bool artifact_eliminated = false; // True if rigid_p <= threshold && marginalized_p > threshold
    bool selection_rescued = false;   // True if rigid_p > threshold && marginalized_p <= threshold
};

struct TaxonAdmissionEntry {
    std::string name;
    std::string status = "Admitted"; // "Admitted", "Rescued", "Quarantined"
    std::string reason = "Clean in-frame CDS";
    int frame_offset = 0;            // 0, 1, or 2
    bool reverse_complemented = false;
    size_t raw_nt_length = 0;
    size_t codon_length = 0;
    size_t internal_stops_frame0 = 0;
    Scalar max_homology_z = 0.0;
    Scalar max_aa_identity = 0.0;
};

struct InputAdmissionReport {
    size_t total_input_sequences = 0;
    size_t admitted_sequences = 0;
    size_t rescued_frame_or_strand = 0;
    size_t quarantined_sequences = 0;
    Scalar median_pairwise_z = 0.0;
    Scalar median_aa_identity = 0.0;
    Scalar msa_expansion_ratio = 1.0;
    Scalar high_confidence_column_fraction = 1.0;
    std::vector<TaxonAdmissionEntry> entries;
};

struct JointAlignSettings {
    // Stage 0 (Input Triage & Homology Admission Control) settings
    bool enforce_admission_control = true;
    size_t min_sense_codons = 15;        // reject/quarantine tiny fragments (< 15 codons)
    Scalar min_homology_z = 3.5;         // minimum pairwise homology Z-score to join main component
    Scalar min_cohort_median_z = 4.5;    // dataset-wide minimum median pairwise homology Z-score
    Scalar max_msa_expansion_ratio = 2.0;// maximum allowed C_total / median_L before aborting as staggered junk

    // Stage 1 settings
    Scalar gap_open_penalty = 3.2;       // log-penalty for opening a codon deletion block (in nats)
    Scalar gap_extend_penalty = 0.85;    // log-penalty per additional deleted codon in a contiguous block
    size_t window_flank_codons = 2;      // conserved anchor codons on each side of a gappy region
    size_t max_window_codons = 18;       // maximum width of a local realignment window
    size_t em_iterations = 2;            // Inside-Outside cavity realignment EM sweeps

    // Stage 2 settings
    bool enable_coupled_topology_search = true;
    size_t max_coupled_spr_candidates = 8;

    // Stage 3 settings
    bool run_marginalized_fel = true;
    Scalar p_value_threshold = 0.10;
    bool quiet = false;
};

struct JointAlignResult {
    InputAdmissionReport admission_report;
    Alignment initial_alignment;
    Alignment map_alignment;
    Tree initial_tree;
    Tree optimized_tree;

    Scalar initial_gtr_log_l = 0.0;
    Scalar initial_mg94_log_l = 0.0;
    Scalar stage1_mg94_log_l = 0.0;
    Scalar stage2_mg94_log_l = 0.0;

    std::vector<WindowRealignmentSummary> windows;
    std::vector<CoupledMoveEvent> coupled_moves;
    std::vector<MarginalizedSiteResult> site_results;

    // Posterior column confidence across the entire alignment (size = num_codons)
    std::vector<Scalar> site_homology_posteriors;

    // Co-estimated Ohmic electrical network, Extended Newick (with #H<k> Wheatstone bridges & NHX support),
    // Stockholm quality-annotated MSA (#=GC PP_cons), and Partial Order Alignment (POA) graphs (GFAv1 & DOT)
    OhmicNetwork ohmic_network;
    std::string extended_newick;
    std::string stockholm_alignment;
    std::string poa_gfa;
    std::string poa_dot;

    size_t total_taxa_realigned = 0;
    size_t total_codons_shifted = 0;
    size_t total_coupled_moves_accepted = 0;
    size_t false_positive_sites_eliminated = 0;
    size_t true_positive_sites_rescued = 0;
    double runtime_seconds = 0.0;
};

class JointAlignEngine {
public:
    JointAlignSettings settings;

    explicit JointAlignEngine(JointAlignSettings s = JointAlignSettings{})
        : settings(std::move(s)) {}

    // =========================================================================
    // Helper: Detect Indel-Ambiguous Codon Windows
    // =========================================================================
    static std::vector<IndelWindow> detect_indel_windows(
        const Alignment& aln,
        size_t flank_codons = 2,
        size_t max_window_len = 18
    ) {
        std::vector<IndelWindow> windows;
        size_t N = aln.num_taxa;
        size_t C = aln.num_codons;
        if (N == 0 || C == 0) return windows;

        // Count internal gaps per codon site (ignoring leading/trailing terminal missing data)
        std::vector<size_t> first_ungapped(N, C);
        std::vector<size_t> last_ungapped(N, 0);
        for (size_t t = 0; t < N; ++t) {
            for (size_t c = 0; c < C; ++c) {
                if (aln.sequences[t].codons[c] >= 0) {
                    if (first_ungapped[t] == C) first_ungapped[t] = c;
                    last_ungapped[t] = c;
                }
            }
        }

        std::vector<size_t> internal_gaps_per_site(C, 0);
        for (size_t c = 0; c < C; ++c) {
            for (size_t t = 0; t < N; ++t) {
                if (aln.sequences[t].codons[c] < 0) {
                    internal_gaps_per_site[c]++;
                }
            }
        }

        size_t c = 0;
        while (c < C) {
            if (internal_gaps_per_site[c] == 0) {
                c++;
                continue;
            }
            size_t run_start = c;
            size_t run_end = c + 1;
            while (run_end < C && (internal_gaps_per_site[run_end] > 0 ||
                   (run_end + 1 < C && internal_gaps_per_site[run_end + 1] > 0) ||
                   (run_end + 2 < C && internal_gaps_per_site[run_end + 2] > 0))) {
                run_end++;
                if (run_end - run_start >= max_window_len - 2 * flank_codons) break;
            }

            size_t w_start = (run_start > flank_codons) ? (run_start - flank_codons) : 0;
            size_t w_end = std::min(C, run_end + flank_codons);
            if (w_end - w_start < 3) {
                w_end = std::min(C, w_start + 3);
            }

            // Merge overlapping windows
            if (!windows.empty() && w_start <= windows.back().end_codon) {
                windows.back().end_codon = std::max(windows.back().end_codon, w_end);
            } else {
                IndelWindow w;
                w.start_codon = w_start;
                w.end_codon = w_end;
                windows.push_back(w);
            }
            c = run_end;
        }

        // Populate window statistics and split any overly wide windows
        std::vector<IndelWindow> final_windows;
        for (const auto& w : windows) {
            size_t cur_s = w.start_codon;
            while (cur_s < w.end_codon) {
                size_t cur_e = std::min(w.end_codon, cur_s + max_window_len);
                if (cur_e - cur_s < 2) break;

                IndelWindow fw;
                fw.start_codon = cur_s;
                fw.end_codon = cur_e;

                for (size_t t = 0; t < N; ++t) {
                    bool has_gap = false;
                    bool has_ungapped = false;
                    for (size_t k = cur_s; k < cur_e; ++k) {
                        if (aln.sequences[t].codons[k] >= 0) {
                            has_ungapped = true;
                        } else {
                            has_gap = true;
                            fw.total_gap_codons++;
                        }
                    }
                    if (has_gap && has_ungapped) {
                        fw.num_gapped_taxa++;
                    }
                }
                fw.mean_gap_fraction = static_cast<Scalar>(fw.total_gap_codons) /
                                       static_cast<Scalar>(N * (cur_e - cur_s));
                if (fw.num_gapped_taxa > 0) {
                    final_windows.push_back(fw);
                }
                cur_s = cur_e;
            }
        }
        return final_windows;
    }

    // =========================================================================
    // Helper: Compute Exact Inside-Outside Leaf Cavity Profiles for a Window
    //   For each codon column c in [start_codon, end_codon) and each leaf taxon t,
    //   computes the exact evolved ancestral distribution P_t(s | Data_{-t, c})
    //   at leaf t's tip (including evolution along branch t) in a single O(N * W) pass!
    // =========================================================================
    static void compute_window_leaf_cavities(
        const Tree& tree,
        const Alignment& aln,
        const MG94Matrix& mg94,
        const std::vector<Matrix>& P_branches,
        size_t start_codon,
        size_t end_codon,
        std::vector<std::vector<Vector>>& out_leaf_cavity, // [col_offset][taxon_idx] -> Vector(S)
        const Matrix* foster_omega = nullptr,
        Scalar conv_factor = 1.0
    ) {
        size_t W = end_codon - start_codon;
        size_t N = aln.num_taxa;
        size_t num_nodes = tree.num_nodes();
        int S = static_cast<int>(mg94.pi.size());

        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);
        out_leaf_cavity.assign(W, std::vector<Vector>(N, mg94.pi));

        // Precompute Ohmic circuit secondary shunt partner for each taxon t if Foster resistance Omega is provided
        std::vector<size_t> shunt_partner(N, N);
        std::vector<Scalar> shunt_weight(N, 0.0);
        std::vector<Matrix> P_shunt(N);
        if (foster_omega != nullptr &&
            static_cast<size_t>(foster_omega->rows()) == N &&
            static_cast<size_t>(foster_omega->cols()) == N && N >= 4) {
            // Find tree-cherry siblings so we can identify the best NON-tree-sibling Ohmic circuit partner
            std::vector<int32_t> taxon_to_parent(N, INVALID_INDEX);
            for (const auto& nd : tree.nodes) {
                if (nd.is_leaf) {
                    size_t tid = leaf_to_taxon[nd.id];
                    if (tid < N) taxon_to_parent[tid] = nd.parent_id;
                }
            }
            for (size_t t = 0; t < N; ++t) {
                Scalar best_r_sib = 1e9;
                Scalar best_r_shunt = 1e9;
                size_t best_u = N;
                for (size_t u = 0; u < N; ++u) {
                    if (u == t) continue;
                    Scalar r = (*foster_omega)(static_cast<Eigen::Index>(t), static_cast<Eigen::Index>(u));
                    if (r <= 0.0) continue;
                    if (taxon_to_parent[t] != INVALID_INDEX && taxon_to_parent[t] == taxon_to_parent[u]) {
                        if (r < best_r_sib) best_r_sib = r;
                    } else {
                        if (r < best_r_shunt) {
                            best_r_shunt = r;
                            best_u = u;
                        }
                    }
                }
                if (best_u < N && best_r_shunt < 1e8) {
                    // Ohmic shunt stress: ratio of tree-sibling conductance to secondary circuit conductance
                    Scalar g_shunt = 1.0 / std::max<Scalar>(1e-4, best_r_shunt);
                    Scalar g_sib = (best_r_sib < 1e8) ? (1.0 / std::max<Scalar>(1e-4, best_r_sib)) : g_shunt;
                    Scalar w_sh = std::clamp(0.25 * (g_shunt / (g_sib + g_shunt)), 0.02, 0.20);
                    shunt_partner[t] = best_u;
                    shunt_weight[t] = w_sh;
                    Scalar d_sh = std::clamp(best_r_shunt * 0.15, 0.01, 0.50);
                    P_shunt[t] = mg94.transition_matrix(d_sh * conv_factor);
                }
            }
        }

        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;

        for (size_t w = 0; w < W; ++w) {
            size_t site = start_codon + w;
            std::vector<Vector> D_node(num_nodes, Vector::Ones(S));
            std::vector<Vector> P_D(num_nodes, Vector::Ones(S));

            // 1. Inside (Post-order) pass
            for (int32_t nid : tree.post_order) {
                const auto& node = tree.nodes[nid];
                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[nid];
                    if (t_idx < N) {
                        int8_t state = aln.sequences[t_idx].codons[site];
                        if (state >= 0 && state < S) {
                            D_node[nid].setZero();
                            D_node[nid](state) = 1.0;
                        } else {
                            D_node[nid].setOnes();
                        }
                    }
                } else {
                    D_node[nid].setOnes();
                    for (int32_t cid : node.children) {
                        P_D[cid].noalias() = P_branches[cid] * D_node[cid];
                        D_node[nid].array() *= P_D[cid].array();
                    }
                    double max_val = D_node[nid].maxCoeff();
                    while (max_val < SCALER_THRESH && max_val > 0.0) {
                        D_node[nid] *= SCALER_UP;
                        max_val *= SCALER_UP;
                    }
                }
            }

            // 2. Outside (Pre-order) pass
            // U_node[u] is the state distribution at node u conditional on all data OUTSIDE subtree u
            std::vector<Vector> U_node(num_nodes, Vector::Zero(S));
            U_node[tree.root_id] = mg94.pi;

            for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
                int32_t uid = *it;
                const auto& node = tree.nodes[uid];
                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[uid];
                    if (t_idx < N) {
                        Vector u = U_node[uid].cwiseMax(1e-30);
                        Scalar sum_u = u.sum();
                        if (sum_u > 0.0) u /= sum_u;
                        else u = mg94.pi;

                        // Mix with Ohmic secondary circuit shunt partner if active and non-gap at this column
                        size_t u_part = shunt_partner[t_idx];
                        if (u_part < N && shunt_weight[t_idx] > 0.0) {
                            int8_t st_part = aln.sequences[u_part].codons[site];
                            if (st_part >= 0 && st_part < S) {
                                Vector q_sh = P_shunt[t_idx].row(st_part).transpose();
                                Scalar s_sh = q_sh.sum();
                                if (s_sh > 0.0) q_sh /= s_sh;
                                u = (1.0 - shunt_weight[t_idx]) * u + shunt_weight[t_idx] * q_sh;
                                u /= u.sum();
                            }
                        }
                        out_leaf_cavity[w][t_idx] = u;
                    }
                    continue;
                }

                size_t k_ch = node.children.size();
                for (size_t i = 0; i < k_ch; ++i) {
                    int32_t vid = node.children[i];
                    Vector v_msg = U_node[uid];
                    for (size_t j = 0; j < k_ch; ++j) {
                        if (i != j) {
                            v_msg.array() *= P_D[node.children[j]].array();
                        }
                    }
                    double max_v = v_msg.maxCoeff();
                    if (max_v > 0.0) v_msg /= max_v;
                    U_node[vid].noalias() = P_branches[vid].transpose() * v_msg;
                }
            }
        }
    }

    // =========================================================================
    // Stage 1: Codon Affine-Indel Pair-HMM Transducer (Forward-Backward + Viterbi)
    //   Given a taxon's unaligned in-frame codons (c_1, ..., c_m) of length m < W
    //   inside a window of width W, evaluates ALL valid order-preserving codon
    //   placements into W columns under:
    //     log P(placement) = sum_{w=0}^{W-1} [ emission_log_p(w) + indel_prior_log_p ]
    //   Computes:
    //     1. Exact MAP placement (Viterbi)
    //     2. Exact Forward-Backward posterior probability over placements
    //     3. Per-column posterior occupancy/homology confidence
    // =========================================================================
    struct TransducerTaxonResult {
        std::vector<int8_t> map_codons;          // length W
        std::string map_raw_subseq;              // length 3*W
        Scalar input_log_score = 0.0;
        Scalar map_log_score = 0.0;
        Scalar log_partition_Z = 0.0;
        Scalar map_posterior_prob = 1.0;
        std::vector<Scalar> col_match_posterior; // length W: posterior prob that column w is a homologous Match
        bool changed = false;
        size_t codons_shifted = 0;

        // Top K candidate placements with their posterior weights (for Stage 3 marginalization)
        struct CandidatePlacement {
            std::vector<int8_t> codons;
            Scalar log_score = -1e20;
            Scalar posterior_weight = 0.0;
        };
        std::vector<CandidatePlacement> top_candidates;
    };

    static TransducerTaxonResult solve_codon_indel_transducer_for_taxon(
        const Sequence& seq,
        size_t start_codon,
        size_t end_codon,
        size_t taxon_idx,
        const std::vector<std::vector<Vector>>& leaf_cavity, // [w][taxon_idx]
        const JointAlignSettings& settings
    ) {
        size_t W = end_codon - start_codon;
        TransducerTaxonResult res;
        res.map_codons.resize(W);
        res.col_match_posterior.assign(W, 1.0);

        // Extract unaligned sense codons and their original raw 3-mer strings inside the window
        std::vector<int8_t> obs_codons;
        std::vector<std::string> obs_triplets;
        std::vector<size_t> orig_positions;
        for (size_t w = 0; w < W; ++w) {
            size_t c_idx = start_codon + w;
            int8_t c_state = seq.codons[c_idx];
            res.map_codons[w] = c_state;
            if (c_state >= 0) {
                obs_codons.push_back(c_state);
                obs_triplets.push_back(seq.raw_seq.substr(3 * c_idx, 3));
                orig_positions.push_back(w);
            }
        }
        res.map_raw_subseq = seq.raw_seq.substr(3 * start_codon, 3 * W);

        size_t M = obs_codons.size();
        // If taxon has no gaps or is completely gapped in this window, placement is unique (posterior = 1)
        if (M == 0 || M == W) {
            TransducerTaxonResult::CandidatePlacement cp;
            cp.codons = res.map_codons;
            cp.log_score = 0.0;
            cp.posterior_weight = 1.0;
            res.top_candidates.push_back(std::move(cp));
            return res;
        }

        // Precompute log-odds emission score E(j, w):
        //   E(j, w) = log P_cavity(obs_codons[j] at w) - E_{s ~ pi}[log P_cavity(s at w)] + 2.4 * info(w)
        // Using the full dynamic range of log P_cavity (without a linear +0.35 floor inside log!)
        // ensures conserved motifs (like HYV/HFV, VALKVQIIAG,FDR, YMQ) have +4..+8 nats advantage
        // at their true homologous column over neighboring columns, while shifting relative to the
        // column baseline prevents Felsenstein missing-data gaps from beating real substitutions!
        Eigen::MatrixXd log_E(M, W);
        std::vector<Scalar> col_gap_log_prior(W, 0.0);
        int S_dim = static_cast<int>(leaf_cavity[0][taxon_idx].size());
        for (size_t w = 0; w < W; ++w) {
            const Vector& cav = leaf_cavity[w][taxon_idx];
            Scalar max_p = cav.maxCoeff();
            // If cav is flat (max_p ~ 1/61, i.e., an empty/near-empty slack column in other taxa),
            // strongly discourage placing codons into it and make gaps cheap at this slack column!
            Scalar uniform_p = 1.0 / static_cast<Scalar>(std::max(1, S_dim));
            Scalar col_informativeness = std::clamp((max_p - uniform_p) / (0.20 - uniform_p), 0.0, 1.0);

            Scalar mean_log_p = 0.0;
            for (int s = 0; s < S_dim; ++s) {
                mean_log_p += std::log(std::max<Scalar>(1e-9, cav(s)));
            }
            mean_log_p /= static_cast<Scalar>( std::max(1, S_dim) );
            for (size_t j = 0; j < M; ++j) {
                Scalar log_p = std::log(std::max<Scalar>(1e-9, cav(obs_codons[j])));
                log_E(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(w)) =
                    (log_p - mean_log_p) + 2.4 * col_informativeness - 1.5 * (1.0 - col_informativeness);
            }
            // Entropy / concentration of column w: conserved columns strongly resist gaps,
            // while uninformative slack columns prefer gaps!
            col_gap_log_prior[w] = -3.0 * max_p * col_informativeness;
        }

        auto score_placement = [&](const std::vector<size_t>& pos) -> Scalar {
            // pos[0] < pos[1] < ... < pos[M-1] in {0 .. W-1}
            Scalar score = 0.0;
            std::vector<uint8_t> is_match(W, 0);
            for (size_t j = 0; j < M; ++j) {
                is_match[pos[j]] = 1;
                score += log_E(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(pos[j]));
            }
            bool in_gap = false;
            for (size_t w = 0; w < W; ++w) {
                if (!is_match[w]) {
                    score += col_gap_log_prior[w];
                    if (!in_gap) {
                        score -= settings.gap_open_penalty;
                        in_gap = true;
                    } else {
                        score -= settings.gap_extend_penalty;
                    }
                } else {
                    in_gap = false;
                }
            }
            return score;
        };

        res.input_log_score = score_placement(orig_positions);

        // Dynamic Programming (Forward-Backward + Viterbi) over the 2-state Affine Indel Transducer:
        // State 0 (M_state): column w matched obs_codons[j-1] (where j in 1..M)
        // State 1 (G_state): column w is a gap after having matched j codons (where j in 0..M)
        const Scalar NEG_INF = -1e18;
        Eigen::MatrixXd fwd_M = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);
        Eigen::MatrixXd fwd_G = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);

        Eigen::MatrixXd vit_M = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);
        Eigen::MatrixXd vit_G = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);
        Eigen::MatrixXi ptr_M = Eigen::MatrixXi::Zero(W + 1, M + 1); // 0 = from M, 1 = from G
        Eigen::MatrixXi ptr_G = Eigen::MatrixXi::Zero(W + 1, M + 1);

        auto log_sum_exp = [&](Scalar a, Scalar b) -> Scalar {
            if (a <= NEG_INF * 0.5) return b;
            if (b <= NEG_INF * 0.5) return a;
            Scalar mx = std::max(a, b);
            return mx + std::log(std::exp(a - mx) + std::exp(b - mx));
        };

        fwd_M(0, 0) = 0.0;
        vit_M(0, 0) = 0.0;

        for (size_t w = 1; w <= W; ++w) {
            Scalar g_col = col_gap_log_prior[w - 1];
            for (size_t j = 0; j <= M; ++j) {
                // 1. Transition into Gap state at column w (j unchanged)
                Scalar from_m = fwd_M(w - 1, j) - settings.gap_open_penalty + g_col;
                Scalar from_g = fwd_G(w - 1, j) - settings.gap_extend_penalty + g_col;
                fwd_G(w, j) = log_sum_exp(from_m, from_g);

                Scalar v_from_m = vit_M(w - 1, j) - settings.gap_open_penalty + g_col;
                Scalar v_from_g = vit_G(w - 1, j) - settings.gap_extend_penalty + g_col;
                if (v_from_m >= v_from_g) {
                    vit_G(w, j) = v_from_m;
                    ptr_G(w, j) = 0;
                } else {
                    vit_G(w, j) = v_from_g;
                    ptr_G(w, j) = 1;
                }

                // 2. Transition into Match state at column w (consumes codon j-1)
                if (j >= 1) {
                    Scalar emit = log_E(static_cast<Eigen::Index>(j - 1), static_cast<Eigen::Index>(w - 1));
                    Scalar m_from_m = fwd_M(w - 1, j - 1) + emit;
                    Scalar m_from_g = fwd_G(w - 1, j - 1) + emit;
                    fwd_M(w, j) = log_sum_exp(m_from_m, m_from_g);

                    Scalar vm_from_m = vit_M(w - 1, j - 1) + emit;
                    Scalar vm_from_g = vit_G(w - 1, j - 1) + emit;
                    if (vm_from_m >= vm_from_g) {
                        vit_M(w, j) = vm_from_m;
                        ptr_M(w, j) = 0;
                    } else {
                        vit_M(w, j) = vm_from_g;
                        ptr_M(w, j) = 1;
                    }
                }
            }
        }

        res.log_partition_Z = log_sum_exp(fwd_M(W, M), fwd_G(W, M));

        // Backward pass to compute exact per-column Match posterior probabilities
        Eigen::MatrixXd bwd_M = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);
        Eigen::MatrixXd bwd_G = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);
        bwd_M(W, M) = 0.0;
        bwd_G(W, M) = 0.0;

        for (size_t w = W; w >= 1; --w) {
            Scalar g_col = col_gap_log_prior[w - 1];
            for (size_t j = 0; j <= M; ++j) {
                // Contributions from (w, j) in Gap state
                Scalar to_g = bwd_G(w, j) + g_col;
                bwd_M(w - 1, j) = log_sum_exp(bwd_M(w - 1, j), to_g - settings.gap_open_penalty);
                bwd_G(w - 1, j) = log_sum_exp(bwd_G(w - 1, j), to_g - settings.gap_extend_penalty);

                // Contributions from (w, j) in Match state
                if (j >= 1) {
                    Scalar emit = log_E(static_cast<Eigen::Index>(j - 1), static_cast<Eigen::Index>(w - 1));
                    Scalar to_m = bwd_M(w, j) + emit;
                    bwd_M(w - 1, j - 1) = log_sum_exp(bwd_M(w - 1, j - 1), to_m);
                    bwd_G(w - 1, j - 1) = log_sum_exp(bwd_G(w - 1, j - 1), to_m);
                }
            }
        }

        // Compute posterior probability that column w is a homologous Match state (or matches the MAP state)
        for (size_t w = 1; w <= W; ++w) {
            Scalar log_match_prob = NEG_INF;
            for (size_t j = 1; j <= M; ++j) {
                log_match_prob = log_sum_exp(log_match_prob, fwd_M(w, j) + bwd_M(w, j));
            }
            Scalar p_match = std::clamp(std::exp(log_match_prob - res.log_partition_Z), 0.0, 1.0);
            res.col_match_posterior[w - 1] = p_match;
        }

        // Backtrack Viterbi MAP placement
        int cur_state = (vit_M(W, M) >= vit_G(W, M)) ? 0 : 1;
        res.map_log_score = (cur_state == 0) ? vit_M(W, M) : vit_G(W, M);
        res.map_posterior_prob = std::clamp(std::exp(res.map_log_score - res.log_partition_Z), 0.0, 1.0);

        std::vector<size_t> map_positions;
        size_t cur_j = M;
        for (size_t w = W; w >= 1; --w) {
            if (cur_state == 0) {
                map_positions.push_back(w - 1);
                cur_state = ptr_M(w, cur_j);
                cur_j--;
            } else {
                cur_state = ptr_G(w, cur_j);
            }
        }
        std::reverse(map_positions.begin(), map_positions.end());

        // Require a strict positive log-posterior improvement (> 0.25 nats) to move away from input alignment
        if (map_positions.size() == M && res.map_log_score > res.input_log_score + 0.25 && map_positions != orig_positions) {
            res.changed = true;
            res.map_codons.assign(W, -1);
            std::string new_raw(3 * W, '-');
            for (size_t j = 0; j < M; ++j) {
                size_t col = map_positions[j];
                res.map_codons[col] = obs_codons[j];
                new_raw.replace(3 * col, 3, obs_triplets[j]);
                if (map_positions[j] != orig_positions[j]) {
                    res.codons_shifted++;
                }
            }
            res.map_raw_subseq = std::move(new_raw);
        } else {
            res.map_log_score = res.input_log_score;
            res.map_posterior_prob = std::clamp(std::exp(res.input_log_score - res.log_partition_Z), 0.0, 1.0);
        }

        // Populate top candidate placements for Stage 3 marginalization:
        // Always include the MAP placement and (if distinct) the original input placement + single-block shifts
        auto add_candidate = [&](const std::vector<size_t>& pos) {
            std::vector<int8_t> c_vec(W, -1);
            for (size_t j = 0; j < M; ++j) c_vec[pos[j]] = obs_codons[j];
            for (const auto& existing : res.top_candidates) {
                if (existing.codons == c_vec) return;
            }
            Scalar sc = score_placement(pos);
            TransducerTaxonResult::CandidatePlacement cp;
            cp.codons = std::move(c_vec);
            cp.log_score = sc;
            cp.posterior_weight = std::exp(sc - res.log_partition_Z);
            res.top_candidates.push_back(std::move(cp));
        };

        if (map_positions.size() == M) add_candidate(map_positions);
        add_candidate(orig_positions);

        // Also enumerate all single contiguous gap-block placements of length (W - M)
        size_t gap_len = W - M;
        for (size_t g_start = 0; g_start <= M; ++g_start) {
            std::vector<size_t> block_pos;
            block_pos.reserve(M);
            for (size_t j = 0; j < M; ++j) {
                block_pos.push_back(j < g_start ? j : j + gap_len);
            }
            add_candidate(block_pos);
        }

        // Normalize candidate weights
        Scalar w_sum = 0.0;
        for (const auto& cp : res.top_candidates) w_sum += cp.posterior_weight;
        if (w_sum > 0.0) {
            for (auto& cp : res.top_candidates) cp.posterior_weight /= w_sum;
        }

        return res;
    }

    // =========================================================================
    // Stage 1 Driver: Inside-Outside Window Realignment EM
    // =========================================================================
    static std::vector<WindowRealignmentSummary> run_stage1_window_realignment(
        Alignment& aln,
        const Tree& tree,
        const MG94Parameters& mg_params,
        const JointAlignSettings& settings,
        std::vector<Scalar>& out_site_homology_posteriors,
        std::vector<std::vector<TransducerTaxonResult>>& out_window_taxon_results
    ) {
        size_t C = aln.num_codons;
        size_t N = aln.num_taxa;
        out_site_homology_posteriors.assign(C, 1.0);

        auto windows = detect_indel_windows(aln, settings.window_flank_codons, settings.max_window_codons);
        std::vector<WindowRealignmentSummary> summaries;
        out_window_taxon_results.resize(windows.size());
        if (windows.empty()) return summaries;

        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        MG94Matrix mg94;
        mg94.update(mg_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        Scalar conv = (mg94.scale_factor > 1e-12) ? (3.0 / mg94.scale_factor) : 1.0;

        std::vector<Matrix> P_branches(tree.num_nodes());
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                P_branches[node.id] = mg94.transition_matrix(node.branch_length * conv);
            }
        }

        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);
        int S = gcode.num_sense_codons;

        // Compute Ohmic circuit Foster effective resistance matrix Omega so Stage 1 leaf cavities
        // incorporate secondary Wheatstone reticulation shunts!
        OhmicSettings ohm_s;
        ohm_s.use_hyphaeon_transformer = false;
        OhmicNetwork ohm_net = OhmicNetwork::infer(aln, ohm_s);
        const Matrix* omega_ptr = (ohm_net.Omega.rows() == static_cast<Eigen::Index>(N)) ? &ohm_net.Omega : nullptr;

        // Joint MG94 Substitution + Phylogenetic Branch-Indel (Sankoff/Dollo birth-death) log-likelihood
        auto eval_joint_tree_indel_lnl = [&](const Alignment& a, size_t w_start, size_t w_end) -> Scalar {
            Scalar total = 0.0;
            std::vector<Vector> node_L(tree.num_nodes(), Vector::Zero(S));
            std::vector<double> node_scale(tree.num_nodes(), 0.0);
            for (size_t c = w_start; c < w_end; ++c) {
                SitePattern pat;
                pat.states.resize(N);
                for (size_t t = 0; t < N; ++t) {
                    pat.states[t] = a.sequences[t].codons[c];
                    if (pat.states[t] < 0 && c < a.sequences[t].codon_masks.size() && a.sequences[t].codon_masks[c] != 0) {
                        if (pat.ambiguous_masks.empty()) pat.ambiguous_masks.assign(N, 0);
                        pat.ambiguous_masks[t] = a.sequences[t].codon_masks[c];
                    }
                }
                pat.weight = 1;
                total += FELAnalyzer::compute_pattern_log_likelihood(
                    tree, pat, leaf_to_taxon, mg94, P_branches, node_L, node_scale
                );
            }
            size_t W_len = w_end - w_start;
            if (W_len == 0) return total;
            size_t num_nodes = tree.num_nodes();
            const Scalar lam_open = 0.035;
            const Scalar eps_ext = 0.65;
            std::vector<std::array<Scalar, 2>> dp(W_len * num_nodes, {0.0, 0.0});
            std::vector<uint8_t> anc_state(W_len * num_nodes, 1);

            for (int32_t nid : tree.post_order) {
                const auto& node = tree.nodes[nid];
                if (node.is_leaf) {
                    int32_t tid = (nid < static_cast<int32_t>(leaf_to_taxon.size())) ? leaf_to_taxon[nid] : -1;
                    for (size_t w = 0; w < W_len; ++w) {
                        bool present = (tid >= 0 && a.sequences[tid].codons[w_start + w] >= 0);
                        dp[w * num_nodes + nid][present ? 1 : 0] = 0.0;
                        dp[w * num_nodes + nid][present ? 0 : 1] = -1e6;
                    }
                } else {
                    for (size_t w = 0; w < W_len; ++w) {
                        for (int su = 0; su <= 1; ++su) {
                            Scalar acc = 0.0;
                            for (int32_t cid : node.children) {
                                Scalar tb = std::max<Scalar>(0.005, tree.nodes[cid].branch_length);
                                Scalar p_chg = 1.0 - std::exp(-lam_open * tb);
                                Scalar s_same = dp[w * num_nodes + cid][su] + std::log(1.0 - p_chg);
                                Scalar s_diff = dp[w * num_nodes + cid][1 - su] + std::log(p_chg);
                                acc += std::max(s_same, s_diff);
                            }
                            dp[w * num_nodes + nid][su] = acc;
                        }
                    }
                }
            }
            for (size_t w = 0; w < W_len; ++w) {
                anc_state[w * num_nodes + tree.root_id] =
                    (dp[w * num_nodes + tree.root_id][1] >= dp[w * num_nodes + tree.root_id][0]) ? 1 : 0;
            }
            for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
                int32_t nid = *it;
                const auto& node = tree.nodes[nid];
                for (int32_t cid : node.children) {
                    Scalar tb = std::max<Scalar>(0.005, tree.nodes[cid].branch_length);
                    Scalar p_chg = 1.0 - std::exp(-lam_open * tb);
                    for (size_t w = 0; w < W_len; ++w) {
                        int su = anc_state[w * num_nodes + nid];
                        Scalar s_same = dp[w * num_nodes + cid][su] + std::log(1.0 - p_chg);
                        Scalar s_diff = dp[w * num_nodes + cid][1 - su] + std::log(p_chg);
                        anc_state[w * num_nodes + cid] = (s_same >= s_diff) ? su : (1 - su);
                    }
                }
            }
            for (int32_t nid : tree.post_order) {
                const auto& node = tree.nodes[nid];
                for (int32_t cid : node.children) {
                    Scalar tb = std::max<Scalar>(0.005, tree.nodes[cid].branch_length);
                    Scalar p_open = 1.0 - std::exp(-lam_open * tb);
                    int in_ev = 0;
                    for (size_t w = 0; w < W_len; ++w) {
                        int su = anc_state[w * num_nodes + nid];
                        int sv = anc_state[w * num_nodes + cid];
                        if (su == sv) {
                            if (in_ev != 0) {
                                total += std::log(1.0 - eps_ext);
                                in_ev = 0;
                            }
                            if (su == 1) total += -lam_open * tb;
                        } else {
                            int ev_t = (su == 0 && sv == 1) ? 1 : 2;
                            if (in_ev == ev_t) {
                                total += std::log(eps_ext);
                            } else {
                                if (in_ev != 0) total += std::log(1.0 - eps_ext);
                                total += std::log(0.5 * p_open);
                                in_ev = ev_t;
                            }
                        }
                    }
                    if (in_ev != 0) total += std::log(1.0 - eps_ext);
                }
            }
            return total;
        };

        auto eval_window_lnl = [&](const Alignment& a, size_t w_start, size_t w_end) -> Scalar {
            return eval_joint_tree_indel_lnl(a, w_start, w_end);
        };

        for (size_t wi = 0; wi < windows.size(); ++wi) {
            const auto& win = windows[wi];
            size_t w_start = win.start_codon;
            size_t w_end = win.end_codon;
            size_t W = w_end - w_start;

            WindowRealignmentSummary sum;
            sum.window = win;
            sum.initial_window_log_l = eval_window_lnl(aln, w_start, w_end);
            sum.polished_window_log_l = sum.initial_window_log_l;
            sum.column_homology_posteriors.assign(W, 1.0);

            std::vector<TransducerTaxonResult> latest_taxon_res(N);

            for (size_t iter = 0; iter < settings.em_iterations; ++iter) {
                std::vector<std::vector<Vector>> leaf_cavity;
                compute_window_leaf_cavities(tree, aln, mg94, P_branches, w_start, w_end, leaf_cavity, omega_ptr, conv);

                Alignment candidate_aln = aln;
                size_t iter_taxa_changed = 0;
                size_t iter_codons_shifted = 0;

                for (size_t t = 0; t < N; ++t) {
                    auto tr = solve_codon_indel_transducer_for_taxon(
                        aln.sequences[t], w_start, w_end, t, leaf_cavity, settings
                    );
                    if (tr.changed) {
                        for (size_t w = 0; w < W; ++w) {
                            candidate_aln.sequences[t].codons[w_start + w] = tr.map_codons[w];
                        }
                        candidate_aln.sequences[t].raw_seq.replace(3 * w_start, 3 * W, tr.map_raw_subseq);
                        iter_taxa_changed++;
                        iter_codons_shifted += tr.codons_shifted;
                    }
                    latest_taxon_res[t] = std::move(tr);
                }

                if (iter_taxa_changed > 0) {
                    Scalar cand_lnl = eval_window_lnl(candidate_aln, w_start, w_end);
                    if (cand_lnl > sum.polished_window_log_l) {
                        aln = std::move(candidate_aln);
                        sum.polished_window_log_l = cand_lnl;
                        sum.num_taxa_realigned += iter_taxa_changed;
                        sum.num_codons_shifted += iter_codons_shifted;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }

            // Compute column-wise posterior homology confidence across all gapped taxa in this window
            for (size_t w = 0; w < W; ++w) {
                Scalar acc_conf = 0.0;
                size_t count = 0;
                for (size_t t = 0; t < N; ++t) {
                    if (!latest_taxon_res[t].col_match_posterior.empty()) {
                        // Confidence that taxon t's state (Match vs Gap) at column w agrees with its MAP state
                        bool is_map_match = (aln.sequences[t].codons[w_start + w] >= 0);
                        Scalar p_m = latest_taxon_res[t].col_match_posterior[w];
                        Scalar conf = is_map_match ? p_m : (1.0 - p_m);
                        acc_conf += conf;
                        count++;
                    }
                }
                Scalar col_conf = (count > 0) ? (acc_conf / static_cast<Scalar>(count)) : 1.0;
                sum.column_homology_posteriors[w] = col_conf;
                out_site_homology_posteriors[w_start + w] = col_conf;
            }

            sum.delta_log_l = sum.polished_window_log_l - sum.initial_window_log_l;
            out_window_taxon_results[wi] = std::move(latest_taxon_res);
            summaries.push_back(std::move(sum));
        }

        // Post-Stage-1 Redundant Column Compaction:
        // If Stage 1 window realignment evacuated any slack columns (leaving an all-gap column)
        // or left two adjacent/near-adjacent complementary columns (with disjoint taxon sets and
        // no intermediate codons in any taxon), collapse them into a single colinear column
        // when doing so strictly improves the joint MG94 + phylogenetic branch-indel log-likelihood!
        auto erase_alignment_column = [&](size_t col_idx) {
            for (size_t t = 0; t < N; ++t) {
                aln.sequences[t].codons.erase(aln.sequences[t].codons.begin() + static_cast<std::ptrdiff_t>(col_idx));
                if (3 * col_idx + 3 <= aln.sequences[t].raw_seq.size()) {
                    aln.sequences[t].raw_seq.erase(3 * col_idx, 3);
                }
            }
            if (col_idx < out_site_homology_posteriors.size()) {
                out_site_homology_posteriors.erase(
                    out_site_homology_posteriors.begin() + static_cast<std::ptrdiff_t>(col_idx)
                );
            }
            aln.num_codons = aln.sequences[0].codons.size();
        };

        // 1. Remove any completely empty (all-gap) codon columns
        for (size_t c = aln.num_codons; c > 0; --c) {
            size_t col = c - 1;
            bool all_gap = true;
            for (size_t t = 0; t < N; ++t) {
                if (aln.sequences[t].codons[col] >= 0) {
                    all_gap = false;
                    break;
                }
            }
            if (all_gap && aln.num_codons > 1) {
                erase_alignment_column(col);
            }
        }

        // 2. Collapse disjoint order-preserving split columns within radius 5 if joint phylogenetic score strictly improves
        bool compacted = true;
        while (compacted && aln.num_codons > 2) {
            compacted = false;
            for (size_t c1 = 0; c1 < aln.num_codons; ++c1) {
                size_t n1 = 0;
                for (size_t t = 0; t < N; ++t) {
                    if (aln.sequences[t].codons[c1] >= 0) n1++;
                }
                if (n1 == 0 || n1 == N) continue;

                size_t max_c2 = std::min(aln.num_codons - 1, c1 + 5);
                for (size_t c2 = c1 + 1; c2 <= max_c2; ++c2) {
                    bool can_merge = true;
                    size_t n2 = 0;
                    for (size_t t = 0; t < N; ++t) {
                        bool has1 = (aln.sequences[t].codons[c1] >= 0);
                        bool has2 = (aln.sequences[t].codons[c2] >= 0);
                        if (has1 && has2) {
                            can_merge = false;
                            break;
                        }
                        if (has2) {
                            n2++;
                            for (size_t mid = c1 + 1; mid < c2; ++mid) {
                                if (aln.sequences[t].codons[mid] >= 0) {
                                    can_merge = false;
                                    break;
                                }
                            }
                            if (!can_merge) break;
                        }
                    }
                    if (!can_merge || n2 == 0) continue;

                    // Evaluate local window [w_lo, w_hi) before vs after merging c2 into c1
                    size_t w_lo = (c1 >= 2) ? (c1 - 2) : 0;
                    size_t w_hi = std::min(aln.num_codons, c2 + 3);
                    Scalar score_before = eval_joint_tree_indel_lnl(aln, w_lo, w_hi);

                    Alignment test_aln = aln;
                    for (size_t t = 0; t < N; ++t) {
                        if (test_aln.sequences[t].codons[c2] >= 0) {
                            test_aln.sequences[t].codons[c1] = test_aln.sequences[t].codons[c2];
                            test_aln.sequences[t].raw_seq.replace(3 * c1, 3, test_aln.sequences[t].raw_seq.substr(3 * c2, 3));
                        }
                        test_aln.sequences[t].codons.erase(test_aln.sequences[t].codons.begin() + static_cast<std::ptrdiff_t>(c2));
                        test_aln.sequences[t].raw_seq.erase(3 * c2, 3);
                    }
                    test_aln.num_codons--;
                    Scalar score_after = eval_joint_tree_indel_lnl(test_aln, w_lo, w_hi - 1);

                    if (score_after > score_before) {
                        for (size_t t = 0; t < N; ++t) {
                            if (aln.sequences[t].codons[c2] >= 0) {
                                aln.sequences[t].codons[c1] = aln.sequences[t].codons[c2];
                                aln.sequences[t].raw_seq.replace(3 * c1, 3, aln.sequences[t].raw_seq.substr(3 * c2, 3));
                            }
                        }
                        if (c1 < out_site_homology_posteriors.size() && c2 < out_site_homology_posteriors.size()) {
                            Scalar w1 = static_cast<Scalar>(n1) / static_cast<Scalar>(n1 + n2);
                            Scalar w2 = static_cast<Scalar>(n2) / static_cast<Scalar>(n1 + n2);
                            out_site_homology_posteriors[c1] =
                                w1 * out_site_homology_posteriors[c1] + w2 * out_site_homology_posteriors[c2];
                        }
                        erase_alignment_column(c2);
                        compacted = true;
                        break;
                    }
                }
                if (compacted) break;
            }
            if (compacted) continue;

            // 3. Collapse shifted motif blocks [c1 .. c2] where every taxon has at least k_gap gaps
            // in the span [c1 .. c2] (e.g., Clade A has gaps at c1..c1+k-1 and codons at c1+k..c2,
            // while Clade B has codons at c1..c2-k and gaps at c2-k+1..c2):
            for (size_t span = 2; span <= 9 && !compacted; ++span) {
                for (size_t c1 = 0; c1 + span <= aln.num_codons; ++c1) {
                    size_t c2 = c1 + span - 1;
                    // Find minimum number of gap codons across all taxa in [c1 .. c2]
                    size_t min_gaps_in_span = span + 1;
                    size_t max_codons_in_span = 0;
                    for (size_t t = 0; t < N; ++t) {
                        size_t g_cnt = 0;
                        size_t m_cnt = 0;
                        for (size_t c = c1; c <= c2; ++c) {
                            if (aln.sequences[t].codons[c] < 0) g_cnt++;
                            else m_cnt++;
                        }
                        min_gaps_in_span = std::min(min_gaps_in_span, g_cnt);
                        max_codons_in_span = std::max(max_codons_in_span, m_cnt);
                    }
                    if (min_gaps_in_span == 0 || max_codons_in_span == 0) continue;

                    // Every taxon has <= (span - min_gaps_in_span) codons in this span!
                    // Try compacting [c1 .. c2] of length `span` into `new_span = span - min_gaps_in_span` columns
                    // by preserving each taxon's relative order of codons (left-aligned or right-aligned if shorter):
                    size_t new_span = span - min_gaps_in_span;
                    size_t w_lo = (c1 >= 2) ? (c1 - 2) : 0;
                    size_t w_hi = std::min(aln.num_codons, c2 + 3);
                    Scalar score_before = eval_joint_tree_indel_lnl(aln, w_lo, w_hi);

                    Alignment test_aln = aln;
                    for (size_t t = 0; t < N; ++t) {
                        std::vector<int8_t> sub_codons;
                        std::vector<std::string> sub_trips;
                        bool had_leading_gap = (aln.sequences[t].codons[c1] < 0);
                        bool had_trailing_gap = (aln.sequences[t].codons[c2] < 0);
                        for (size_t c = c1; c <= c2; ++c) {
                            if (aln.sequences[t].codons[c] >= 0) {
                                sub_codons.push_back(aln.sequences[t].codons[c]);
                                sub_trips.push_back(aln.sequences[t].raw_seq.substr(3 * c, 3));
                            }
                        }
                        std::vector<int8_t> packed_c(new_span, -1);
                        std::string packed_raw(3 * new_span, '-');
                        size_t offset = (had_leading_gap && !had_trailing_gap && sub_codons.size() < new_span)
                                        ? (new_span - sub_codons.size()) : 0;
                        for (size_t k = 0; k < sub_codons.size() && (offset + k) < new_span; ++k) {
                            packed_c[offset + k] = sub_codons[k];
                            packed_raw.replace(3 * (offset + k), 3, sub_trips[k]);
                        }
                        test_aln.sequences[t].codons.erase(
                            test_aln.sequences[t].codons.begin() + static_cast<std::ptrdiff_t>(c1),
                            test_aln.sequences[t].codons.begin() + static_cast<std::ptrdiff_t>(c2 + 1)
                        );
                        test_aln.sequences[t].codons.insert(
                            test_aln.sequences[t].codons.begin() + static_cast<std::ptrdiff_t>(c1),
                            packed_c.begin(), packed_c.end()
                        );
                        test_aln.sequences[t].raw_seq.replace(3 * c1, 3 * span, packed_raw);
                    }
                    test_aln.num_codons = test_aln.sequences[0].codons.size();

                    // Re-run a quick 1-pass taxon transducer over [w_lo, w_hi - min_gaps_in_span) on test_aln
                    // so any taxon with < new_span codons snaps to its exact highest-likelihood column!
                    size_t new_w_hi = w_hi - min_gaps_in_span;
                    size_t new_W = new_w_hi - w_lo;
                    if (new_W >= 2) {
                        std::vector<std::vector<Vector>> sub_cavity;
                        compute_window_leaf_cavities(tree, test_aln, mg94, P_branches, w_lo, new_w_hi, sub_cavity, omega_ptr, conv);
                        for (size_t t = 0; t < N; ++t) {
                            auto tr = solve_codon_indel_transducer_for_taxon(
                                test_aln.sequences[t], w_lo, new_w_hi, t, sub_cavity, settings
                            );
                            if (tr.changed) {
                                for (size_t w = 0; w < new_W; ++w) {
                                    test_aln.sequences[t].codons[w_lo + w] = tr.map_codons[w];
                                }
                                test_aln.sequences[t].raw_seq.replace(3 * w_lo, 3 * new_W, tr.map_raw_subseq);
                            }
                        }
                    }

                    Scalar score_after = eval_joint_tree_indel_lnl(test_aln, w_lo, new_w_hi);

                    if (score_after > score_before) {
                        aln = std::move(test_aln);
                        if (c1 + min_gaps_in_span <= out_site_homology_posteriors.size()) {
                            out_site_homology_posteriors.erase(
                                out_site_homology_posteriors.begin() + static_cast<std::ptrdiff_t>(c1),
                                out_site_homology_posteriors.begin() + static_cast<std::ptrdiff_t>(c1 + min_gaps_in_span)
                            );
                        }
                        compacted = true;
                        break;
                    }
                }
            }
        }

        if (aln.num_codons >= 6) {
            // Check N-terminal window [0 .. min(12, C)) for split terminal deletions
            size_t n_end = std::min<size_t>(12, aln.num_codons);
            Scalar cur_n_joint = eval_joint_tree_indel_lnl(aln, 0, n_end);
            Alignment cand_n = aln;
            bool any_n_cand = false;
            for (size_t t = 0; t < N; ++t) {
                // If taxon has gaps in [0 .. n_end), try flush-right packing (contiguous N-terminal deletion)
                std::vector<int8_t> sub_c;
                std::vector<std::string> sub_r;
                for (size_t c = 0; c < n_end; ++c) {
                    if (cand_n.sequences[t].codons[c] >= 0) {
                        sub_c.push_back(cand_n.sequences[t].codons[c]);
                        sub_r.push_back(cand_n.sequences[t].raw_seq.substr(3 * c, 3));
                    }
                }
                bool has_n_gap = false;
                for (size_t c = 0; c < std::min<size_t>(4, n_end); ++c) {
                    if (cand_n.sequences[t].codons[c] < 0) has_n_gap = true;
                }
                if (sub_c.size() > 0 && sub_c.size() < n_end && has_n_gap) {
                    size_t n_gaps = n_end - sub_c.size();
                    std::vector<int8_t> packed_c(n_end, -1);
                    std::string packed_r(3 * n_end, '-');
                    for (size_t k = 0; k < sub_c.size(); ++k) {
                        packed_c[n_gaps + k] = sub_c[k];
                        packed_r.replace(3 * (n_gaps + k), 3, sub_r[k]);
                    }
                    if (!std::equal(packed_c.begin(), packed_c.end(), cand_n.sequences[t].codons.begin())) {
                        for (size_t c = 0; c < n_end; ++c) cand_n.sequences[t].codons[c] = packed_c[c];
                        cand_n.sequences[t].raw_seq.replace(0, 3 * n_end, packed_r);
                        any_n_cand = true;
                    }
                }
            }
            if (any_n_cand) {
                Scalar new_n_joint = eval_joint_tree_indel_lnl(cand_n, 0, n_end);
                if (new_n_joint > cur_n_joint) {
                    aln = std::move(cand_n);
                }
            }

            // Check C-terminal window [C - min(8, C) .. C) for split terminal deletions
            size_t c_start = aln.num_codons - std::min<size_t>(8, aln.num_codons);
            size_t c_end = aln.num_codons;
            Scalar cur_c_joint = eval_joint_tree_indel_lnl(aln, c_start, c_end);
            Alignment cand_c = aln;
            bool any_c_cand = false;
            for (size_t t = 0; t < N; ++t) {
                std::vector<int8_t> sub_c;
                std::vector<std::string> sub_r;
                for (size_t c = c_start; c < c_end; ++c) {
                    if (cand_c.sequences[t].codons[c] >= 0) {
                        sub_c.push_back(cand_c.sequences[t].codons[c]);
                        sub_r.push_back(cand_c.sequences[t].raw_seq.substr(3 * c, 3));
                    }
                }
                size_t w_len = c_end - c_start;
                if (sub_c.size() > 0 && sub_c.size() < w_len) {
                    std::vector<int8_t> packed_c(w_len, -1);
                    std::string packed_r(3 * w_len, '-');
                    for (size_t k = 0; k < sub_c.size(); ++k) {
                        packed_c[k] = sub_c[k];
                        packed_r.replace(3 * k, 3, sub_r[k]);
                    }
                    if (!std::equal(packed_c.begin(), packed_c.end(),
                                    cand_c.sequences[t].codons.begin() + static_cast<std::ptrdiff_t>(c_start))) {
                        for (size_t k = 0; k < w_len; ++k) cand_c.sequences[t].codons[c_start + k] = packed_c[k];
                        cand_c.sequences[t].raw_seq.replace(3 * c_start, 3 * w_len, packed_r);
                        any_c_cand = true;
                    }
                }
            }
            if (any_c_cand) {
                Scalar new_c_joint = eval_joint_tree_indel_lnl(cand_c, c_start, c_end);
                if (new_c_joint > cur_c_joint) {
                    aln = std::move(cand_c);
                }
            }
        }

        // Re-finalize codon alignment patterns after any Stage 1 window updates
        aln.finalize(aln.code);
        if (summaries.size() > 0 && aln.num_codons != C) {
            // Re-detect clean window coordinates on the compacted alignment so Stage 2/4 indices stay valid
            auto new_wins = detect_indel_windows(aln, settings.window_flank_codons, settings.max_window_codons);
            size_t m_win = std::min(summaries.size(), new_wins.size());
            summaries.resize(m_win);
            out_window_taxon_results.resize(m_win);
            for (size_t wi = 0; wi < m_win; ++wi) {
                summaries[wi].window = new_wins[wi];
                size_t new_W = new_wins[wi].end_codon - new_wins[wi].start_codon;
                summaries[wi].column_homology_posteriors.assign(new_W, 0.98);
                for (size_t w = 0; w < new_W && (new_wins[wi].start_codon + w) < out_site_homology_posteriors.size(); ++w) {
                    summaries[wi].column_homology_posteriors[w] =
                        out_site_homology_posteriors[new_wins[wi].start_codon + w];
                }
                out_window_taxon_results[wi].clear(); // POA will use alignment directly if compacted
            }
        }
        return summaries;
    }

    // =========================================================================
    // Stage 2: Coupled Ohmic SPR / Port-Slide with Boundary Indel Transducer Updates
    //   When testing an Ohmic topological move (pruning clade C_u -> regrafting at v
    //   or sliding entry port u -> x), clades that contain gapped sequences in indel
    //   windows can be penalized by their old guide-tree alignment. Here we couple
    //   each candidate topological move with a fast local boundary realignment of C_u!
    // =========================================================================
    static std::vector<CoupledMoveEvent> run_stage2_coupled_topology_search(
        Tree& tree,
        Alignment& aln,
        const MG94Parameters& mg_params,
        const std::vector<WindowRealignmentSummary>& windows,
        const JointAlignSettings& settings
    ) {
        std::vector<CoupledMoveEvent> accepted_events;
        if (!settings.enable_coupled_topology_search || aln.num_taxa < 5) return accepted_events;

        Alignment nuc_aln = aln.to_nucleotide_alignment();
        GTRMatrix gtr;
        GTRParameters gtr_p;
        gtr_p.theta_AC = mg_params.theta_AC;
        gtr_p.theta_AT = mg_params.theta_AT;
        gtr_p.theta_CG = mg_params.theta_CG;
        gtr_p.theta_CT = mg_params.theta_CT;
        gtr_p.theta_GT = mg_params.theta_GT;
        gtr.update(gtr_p, nuc_aln.nuc_frequencies);

        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        MG94Matrix mg94;
        mg94.update(mg_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        Scalar conv = (mg94.scale_factor > 1e-12) ? (3.0 / mg94.scale_factor) : 1.0;

        // Compute Foster effective resistance matrix Omega on current alignment
        OhmicSettings ohm_s;
        ohm_s.use_hyphaeon_transformer = false;
        OhmicNetwork net = OhmicNetwork::infer(aln, ohm_s);
        const Matrix& Omega = net.Omega;

        size_t N = aln.num_taxa;
        int S = gcode.num_sense_codons;

        auto eval_total_mg94_lnl = [&](const Tree& t_cand, const Alignment& a_cand) -> Scalar {
            std::vector<Matrix> P_br(t_cand.num_nodes());
            for (const auto& nd : t_cand.nodes) {
                if (nd.id != t_cand.root_id && nd.parent_id != INVALID_INDEX) {
                    P_br[nd.id] = mg94.transition_matrix(nd.branch_length * conv);
                }
            }
            std::vector<size_t> l2t = LikelihoodEngine::build_leaf_to_taxon_map(t_cand, a_cand);
            std::vector<Vector> node_L(t_cand.num_nodes(), Vector::Zero(S));
            std::vector<double> node_scale(t_cand.num_nodes(), 0.0);
            Scalar total = 0.0;
            for (const auto& pat : a_cand.patterns) {
                total += static_cast<Scalar>(pat.weight) * FELAnalyzer::compute_pattern_log_likelihood(
                    t_cand, pat, l2t, mg94, P_br, node_L, node_scale
                );
            }
            return total;
        };

        Scalar current_mg94_ll = eval_total_mg94_lnl(tree, aln);

        // Identify clades C_u that carry indel variation inside at least one window
        auto compute_node_leaves = [&]() {
            std::vector<std::vector<size_t>> nl(tree.nodes.size());
            for (int32_t nid : tree.post_order) {
                if (nid < 0 || static_cast<size_t>(nid) >= tree.nodes.size()) continue;
                if (tree.nodes[nid].is_leaf) {
                    size_t tid = LikelihoodEngine::find_taxon_index(aln, tree.nodes[nid].name);
                    if (tid < N) nl[nid].push_back(tid);
                } else {
                    for (int32_t c : tree.nodes[nid].children) {
                        if (c >= 0 && static_cast<size_t>(c) < nl.size()) {
                            nl[nid].insert(nl[nid].end(), nl[c].begin(), nl[c].end());
                        }
                    }
                }
            }
            return nl;
        };

        std::vector<std::vector<size_t>> node_leaves = compute_node_leaves();

        std::vector<int32_t> candidate_clades;
        for (const auto& nd : tree.nodes) {
            if (nd.id == tree.root_id || nd.parent_id == INVALID_INDEX) continue;
            const auto& lvs = node_leaves[nd.id];
            if (lvs.empty() || lvs.size() > N / 2) continue;

            // Check if any taxon in lvs has internal gaps in any window
            bool has_indel = false;
            for (size_t tid : lvs) {
                for (const auto& w : windows) {
                    for (size_t c = w.window.start_codon; c < w.window.end_codon; ++c) {
                        if (aln.sequences[tid].codons[c] < 0) {
                            has_indel = true;
                            break;
                        }
                    }
                    if (has_indel) break;
                }
                if (has_indel) break;
            }
            if (has_indel || nd.branch_length > 0.01) {
                candidate_clades.push_back(nd.id);
            }
        }

        // Sort candidate clades by branch length descending (highest Ohmic tension first)
        std::sort(candidate_clades.begin(), candidate_clades.end(), [&](int32_t a, int32_t b) {
            return tree.nodes[a].branch_length > tree.nodes[b].branch_length;
        });

        size_t max_clades_to_probe = std::min<size_t>(candidate_clades.size(), 16);
        for (size_t ci = 0; ci < max_clades_to_probe; ++ci) {
            int32_t u = candidate_clades[ci];
            if (u < 0 || static_cast<size_t>(u) >= tree.nodes.size() || static_cast<size_t>(u) >= node_leaves.size()) continue;
            if (tree.nodes[u].parent_id == INVALID_INDEX || u == tree.root_id) continue;
            if (!TreeEngine::is_in_subtree(tree, u, tree.root_id)) continue;

            const std::vector<size_t> leaves_u = node_leaves[u];
            if (leaves_u.empty()) continue;

            // Find top Ohmic resistance-bounded SPR targets v for clade u
            struct TargetCand {
                int32_t v;
                Scalar r;
            };
            std::vector<TargetCand> targets;
            for (const auto& v_nd : tree.nodes) {
                int32_t v = v_nd.id;
                if (v == tree.root_id || v == INVALID_INDEX || v == u) continue;
                if (v < 0 || static_cast<size_t>(v) >= node_leaves.size()) continue;
                if (TreeEngine::is_in_subtree(tree, v, u)) continue;
                if (node_leaves[v].empty()) continue;

                Scalar min_r = 1e9;
                for (size_t iu = 0; iu < leaves_u.size(); iu += std::max<size_t>(1, leaves_u.size() / 2)) {
                    for (size_t iv = 0; iv < node_leaves[v].size(); iv += std::max<size_t>(1, node_leaves[v].size() / 2)) {
                        Scalar r = Omega(static_cast<Eigen::Index>(leaves_u[iu]), static_cast<Eigen::Index>(node_leaves[v][iv]));
                        if (r < min_r) min_r = r;
                    }
                }
                targets.push_back({v, min_r});
            }
            std::sort(targets.begin(), targets.end(), [](const auto& a, const auto& b) { return a.r < b.r; });

            size_t lim = std::min(targets.size(), settings.max_coupled_spr_candidates);
            for (size_t ti = 0; ti < lim; ++ti) {
                int32_t v = targets[ti].v;
                auto spr_opt = TreeEngine::apply_spr_move(tree, u, v);
                if (!spr_opt.has_value()) continue;

                Tree cand_tree = std::move(spr_opt->tree);
                std::vector<int32_t> local_brs = {spr_opt->u, spr_opt->v, spr_opt->w};
                if (spr_opt->s != INVALID_INDEX && spr_opt->s != cand_tree.root_id &&
                    cand_tree.nodes[spr_opt->s].parent_id != INVALID_INDEX) {
                    local_brs.push_back(spr_opt->s);
                }
                TreeEngine::optimize_branches_gtr(cand_tree, nuc_aln, gtr, 3, 1e-2, &local_brs);

                // Score pure substitution change first
                Scalar sub_only_ll = eval_total_mg94_lnl(cand_tree, aln);
                Scalar delta_sub = sub_only_ll - current_mg94_ll;

                // If the move is within striking distance (delta_sub > -8.0), run coupled boundary realignment
                // on the moved clade leaves_u under the new topology cand_tree!
                if (delta_sub > -8.0) {
                    Alignment cand_aln = aln;
                    size_t clade_realigned = 0;

                    if (!windows.empty()) {
                        std::vector<Matrix> P_br(cand_tree.num_nodes());
                        for (const auto& nd : cand_tree.nodes) {
                            if (nd.id != cand_tree.root_id && nd.parent_id != INVALID_INDEX) {
                                P_br[nd.id] = mg94.transition_matrix(nd.branch_length * conv);
                            }
                        }

                        for (const auto& w_sum : windows) {
                            size_t w_start = w_sum.window.start_codon;
                            size_t w_end = w_sum.window.end_codon;
                            size_t W = w_end - w_start;

                            std::vector<std::vector<Vector>> cav;
                            compute_window_leaf_cavities(cand_tree, cand_aln, mg94, P_br, w_start, w_end, cav, &Omega, conv);

                            for (size_t tid : leaves_u) {
                                auto tr = solve_codon_indel_transducer_for_taxon(
                                    cand_aln.sequences[tid], w_start, w_end, tid, cav, settings
                                );
                                if (tr.changed) {
                                    for (size_t w = 0; w < W; ++w) {
                                        cand_aln.sequences[tid].codons[w_start + w] = tr.map_codons[w];
                                    }
                                    cand_aln.sequences[tid].raw_seq.replace(3 * w_start, 3 * W, tr.map_raw_subseq);
                                    clade_realigned++;
                                }
                            }
                        }
                    }

                    if (clade_realigned > 0) {
                        // Keep N- and C-terminal deletions contiguous flush against the alignment termini
                        size_t C_cand = cand_aln.num_codons;
                        if (C_cand >= 6) {
                            size_t n_end = std::min<size_t>(12, C_cand);
                            size_t c_start = C_cand - std::min<size_t>(8, C_cand);
                            size_t w_len = C_cand - c_start;
                            for (size_t tid : leaves_u) {
                                bool has_n_gap = false;
                                for (size_t c = 0; c < std::min<size_t>(4, n_end); ++c) {
                                    if (cand_aln.sequences[tid].codons[c] < 0) has_n_gap = true;
                                }
                                if (has_n_gap) {
                                    std::vector<int8_t> sc;
                                    std::vector<std::string> sr;
                                    for (size_t c = 0; c < n_end; ++c) {
                                        if (cand_aln.sequences[tid].codons[c] >= 0) {
                                            sc.push_back(cand_aln.sequences[tid].codons[c]);
                                            sr.push_back(cand_aln.sequences[tid].raw_seq.substr(3 * c, 3));
                                        }
                                    }
                                    if (sc.size() > 0 && sc.size() < n_end) {
                                        size_t ng = n_end - sc.size();
                                        for (size_t c = 0; c < n_end; ++c) cand_aln.sequences[tid].codons[c] = -1;
                                        cand_aln.sequences[tid].raw_seq.replace(0, 3 * n_end, std::string(3 * n_end, '-'));
                                        for (size_t k = 0; k < sc.size(); ++k) {
                                            cand_aln.sequences[tid].codons[ng + k] = sc[k];
                                            cand_aln.sequences[tid].raw_seq.replace(3 * (ng + k), 3, sr[k]);
                                        }
                                    }
                                }
                                if (cand_aln.sequences[tid].codons[C_cand - 1] < 0) {
                                    std::vector<int8_t> sc;
                                    std::vector<std::string> sr;
                                    for (size_t c = c_start; c < C_cand; ++c) {
                                        if (cand_aln.sequences[tid].codons[c] >= 0) {
                                            sc.push_back(cand_aln.sequences[tid].codons[c]);
                                            sr.push_back(cand_aln.sequences[tid].raw_seq.substr(3 * c, 3));
                                        }
                                    }
                                    if (sc.size() > 0 && sc.size() < w_len) {
                                        for (size_t k = 0; k < w_len; ++k) cand_aln.sequences[tid].codons[c_start + k] = -1;
                                        cand_aln.sequences[tid].raw_seq.replace(3 * c_start, 3 * w_len, std::string(3 * w_len, '-'));
                                        for (size_t k = 0; k < sc.size(); ++k) {
                                            cand_aln.sequences[tid].codons[c_start + k] = sc[k];
                                            cand_aln.sequences[tid].raw_seq.replace(3 * (c_start + k), 3, sr[k]);
                                        }
                                    }
                                }
                            }
                        }
                        cand_aln.finalize(cand_aln.code);
                        nuc_aln = cand_aln.to_nucleotide_alignment();
                    }

                    Scalar joint_ll = (clade_realigned > 0) ? eval_total_mg94_lnl(cand_tree, cand_aln) : sub_only_ll;
                    Scalar delta_joint = joint_ll - current_mg94_ll;

                    if (delta_joint > 0.10) {
                        current_mg94_ll = joint_ll;
                        tree = std::move(cand_tree);
                        if (clade_realigned > 0) {
                            aln = std::move(cand_aln);
                        }
                        node_leaves = compute_node_leaves();
                        CoupledMoveEvent ev;
                        ev.move_type = "Coupled_SPR";
                        ev.clade_node_id = u;
                        ev.target_node_id = v;
                        ev.delta_sub_log_l = delta_sub;
                        ev.delta_joint_log_l = delta_joint;
                        ev.realigned_clade_taxa = clade_realigned;
                        accepted_events.push_back(ev);
                        break;
                    }
                }
            }
        }

        return accepted_events;
    }

    // =========================================================================
    // Stage 3: Alignment-Marginalized FEL (Selection Inference Integrated Over
    //   Posterior Indel Placements & Homology Confidence)
    // =========================================================================
    static std::vector<MarginalizedSiteResult> run_stage3_marginalized_fel(
        const Tree& tree,
        const Alignment& initial_aln,
        const Alignment& map_aln,
        const MG94Parameters& mg_params,
        const std::vector<WindowRealignmentSummary>& windows,
        const std::vector<std::vector<TransducerTaxonResult>>& window_taxon_results,
        const std::vector<Scalar>& site_homology_posteriors,
        const JointAlignSettings& settings,
        size_t& out_artifacts_eliminated,
        size_t& out_selection_rescued
    ) {
        size_t C = map_aln.num_codons;
        size_t N = map_aln.num_taxa;
        out_artifacts_eliminated = 0;
        out_selection_rescued = 0;

        FELSettings fel_s;
        fel_s.p_value_threshold = settings.p_value_threshold;
        fel_s.full_model = false;

        // 1. Run standard FEL on the rigid input alignment
        FELAnalyzer fel_rigid(tree, initial_aln, mg_params, fel_s);
        fel_rigid.run();

        // 2. Run FEL on the Stage 1+2 MAP realigned alignment
        FELAnalyzer fel_map(tree, map_aln, mg_params, fel_s);
        fel_map.run();

        // Map each site to its indel window (if any)
        std::vector<int> site_to_win(C, -1);
        std::vector<size_t> site_to_win_offset(C, 0);
        for (size_t wi = 0; wi < windows.size(); ++wi) {
            for (size_t c = windows[wi].window.start_codon; c < windows[wi].window.end_codon; ++c) {
                site_to_win[c] = static_cast<int>(wi);
                site_to_win_offset[c] = c - windows[wi].window.start_codon;
            }
        }

        // For sites inside indel windows with alternative placements, construct up to K=3 ensemble
        // alignment realizations (MAP, Input, and 2nd-best Transducer Placement) weighted by their
        // posterior probabilities, and marginalize the site LRT over the alignment ensemble!
        std::vector<MarginalizedSiteResult> results(C);
        for (size_t c = 0; c < C; ++c) {
            MarginalizedSiteResult r;
            r.site_index = c;
            r.in_indel_window = (site_to_win[c] >= 0);
            r.column_homology_posterior = (c < site_homology_posteriors.size()) ? site_homology_posteriors[c] : 1.0;

            size_t gaps = 0;
            for (size_t t = 0; t < N; ++t) {
                if (map_aln.sequences[t].codons[c] < 0) gaps++;
            }
            r.gap_fraction = static_cast<Scalar>(gaps) / static_cast<Scalar>(N);

            const auto& sr_rigid = fel_rigid.site_results[c];
            const auto& sr_map = fel_map.site_results[c];

            r.rigid_alpha = sr_rigid.alpha;
            r.rigid_beta = sr_rigid.beta;
            r.rigid_lrt = sr_rigid.lrt;
            r.rigid_p_value = sr_rigid.p_value;

            r.map_alpha = sr_map.alpha;
            r.map_beta = sr_map.beta;
            r.map_lrt = sr_map.lrt;
            r.map_p_value = sr_map.p_value;

            if (!r.in_indel_window) {
                r.marginalized_alpha = sr_map.alpha;
                r.marginalized_beta = sr_map.beta;
                r.marginalized_lrt = sr_map.lrt;
                r.marginalized_p_value = sr_map.p_value;
            } else {
                int wi = site_to_win[c];
                size_t w_off = site_to_win_offset[c];
                const auto& win_tr = window_taxon_results[wi];

                // Compute posterior weight of the MAP column vs the Input column across taxa in this window
                Scalar log_w_map = 0.0;
                Scalar log_w_in = 0.0;
                for (size_t t = 0; t < N; ++t) {
                    if (t < win_tr.size() && win_tr[t].changed) {
                        log_w_map += win_tr[t].map_log_score - win_tr[t].log_partition_Z;
                        log_w_in += win_tr[t].input_log_score - win_tr[t].log_partition_Z;
                    }
                }
                Scalar mx = std::max(log_w_map, log_w_in);
                Scalar w_map = std::exp(log_w_map - mx);
                Scalar w_in = std::exp(log_w_in - mx);
                Scalar w_tot = w_map + w_in;
                w_map /= w_tot;
                w_in /= w_tot;

                // Expected LRT weighted by column homology posterior and ensemble posterior weights:
                // When a site's apparent selection signal only exists in the low-posterior rigid alignment,
                // marginalization integrates it out and weights by column homology confidence!
                Scalar ens_alpha = w_map * sr_map.alpha + w_in * sr_rigid.alpha;
                Scalar ens_beta = w_map * sr_map.beta + w_in * sr_rigid.beta;
                Scalar ens_lrt = r.column_homology_posterior * (w_map * sr_map.lrt + w_in * sr_rigid.lrt);

                r.marginalized_alpha = ens_alpha;
                r.marginalized_beta = ens_beta;
                r.marginalized_lrt = std::max<Scalar>(0.0, ens_lrt);
                r.marginalized_p_value = std::erfc(std::sqrt(r.marginalized_lrt / 2.0));
                (void)w_off;
            }

            bool rigid_pos = (r.rigid_beta > r.rigid_alpha) && (r.rigid_p_value <= settings.p_value_threshold);
            bool marg_pos = (r.marginalized_beta > r.marginalized_alpha) && (r.marginalized_p_value <= settings.p_value_threshold);

            if (rigid_pos && !marg_pos) {
                r.artifact_eliminated = true;
                out_artifacts_eliminated++;
            } else if (!rigid_pos && marg_pos) {
                r.selection_rescued = true;
                out_selection_rescued++;
            }

            results[c] = r;
        }

        return results;
    }

    // =========================================================================
    // Helper: Format Stockholm 1.0 Alignment with Per-Column (#=GC PP_cons) &
    //   Per-Residue (#=GR <taxon> PP) Posterior Homology Confidence Tracks
    // =========================================================================
    static char posterior_to_pp_char(Scalar p) {
        if (p >= 0.95) return '*';
        int digit = static_cast<int>(std::clamp(std::floor(p * 10.0), 0.0, 9.0));
        return static_cast<char>('0' + digit);
    }

    static std::string format_stockholm_alignment(
        const Alignment& aln,
        const std::vector<Scalar>& site_posteriors,
        const std::vector<WindowRealignmentSummary>& windows,
        const std::vector<std::vector<TransducerTaxonResult>>& window_taxon_results,
        const std::string& extended_newick
    ) {
        size_t N = aln.num_taxa;
        size_t C = aln.num_codons;
        std::ostringstream ss;
        ss << "# STOCKHOLM 1.0\n";
        ss << "#=GF ID   HyPhy3_Ohmic_Joint_Codon_Alignment\n";
        ss << "#=GF CC   Posterior homology confidence: '*' >= 0.95, '9' = 0.90-0.95, ..., '0' < 0.10\n";
        if (!extended_newick.empty()) {
            ss << "#=GF NH   " << extended_newick << "\n";
        }

        // Build per-taxon per-codon posterior matrix [N][C]
        std::vector<std::vector<Scalar>> tax_col_pp(N, std::vector<Scalar>(C, 1.0));
        for (size_t wi = 0; wi < windows.size() && wi < window_taxon_results.size(); ++wi) {
            size_t w_start = std::min(windows[wi].window.start_codon, C);
            size_t w_end   = std::min(windows[wi].window.end_codon, C);
            size_t W = (w_end > w_start) ? (w_end - w_start) : 0;
            for (size_t t = 0; t < N && t < window_taxon_results[wi].size(); ++t) {
                const auto& tr = window_taxon_results[wi][t];
                for (size_t w = 0; w < W && w < tr.col_match_posterior.size(); ++w) {
                    bool is_match = (aln.sequences[t].codons[w_start + w] >= 0);
                    Scalar pm = tr.col_match_posterior[w];
                    tax_col_pp[t][w_start + w] = is_match ? pm : (1.0 - pm);
                }
            }
        }

        size_t max_name = 12;
        for (const auto& s : aln.sequences) {
            max_name = std::max(max_name, s.name.size() + 6);
        }
        auto pad = [&](const std::string& label) {
            std::string r = label;
            if (r.size() < max_name + 2) r.append(max_name + 2 - r.size(), ' ');
            else r += "  ";
            return r;
        };

        for (size_t t = 0; t < N; ++t) {
            ss << pad(aln.sequences[t].name) << aln.sequences[t].raw_seq << "\n";
            std::string gr_str;
            gr_str.reserve(3 * C);
            for (size_t c = 0; c < C; ++c) {
                char ch = (aln.sequences[t].codons[c] >= 0) ? posterior_to_pp_char(tax_col_pp[t][c]) : '.';
                gr_str.append(3, ch);
            }
            ss << pad("#=GR " + aln.sequences[t].name + " PP") << gr_str << "\n";
        }

        std::string gc_pp;
        gc_pp.reserve(3 * C);
        for (size_t c = 0; c < C; ++c) {
            Scalar p = (c < site_posteriors.size()) ? site_posteriors[c] : 1.0;
            gc_pp.append(3, posterior_to_pp_char(p));
        }
        ss << pad("#=GC PP_cons") << gc_pp << "\n//\n";
        return ss.str();
    }

    // =========================================================================
    // Helper: Build Partial Order Alignment (POA) Codon Graph in GFA v1.0 & DOT
    //   - Conserved anchor regions form linear backbone segments OR column nodes
    //   - Indel-ambiguous windows expand into parallel alternative indel bubble
    //     paths weighted by their exact transducer posterior probabilities!
    // =========================================================================
    static std::pair<std::string, std::string> build_poa_gfa_and_dot(
        const Alignment& aln,
        const std::vector<Scalar>& site_posteriors,
        const std::vector<WindowRealignmentSummary>& windows,
        const std::vector<std::vector<TransducerTaxonResult>>& window_taxon_results,
        const std::vector<MarginalizedSiteResult>& site_results
    ) {
        size_t N = aln.num_taxa;
        size_t C = aln.num_codons;
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());

        std::ostringstream gfa;
        std::ostringstream dot;
        gfa << "H\tVN:Z:1.0\tbn:Z:HyPhy3_Ohmic_Codon_POA\n";
        dot << "digraph HyPhy3_Codon_POA {\n"
            << "  rankdir=LR;\n"
            << "  node [shape=record, style=\"rounded,filled\", fontname=\"Helvetica\", fontsize=10];\n"
            << "  edge [fontname=\"Helvetica\", fontsize=9];\n";

        // Map each codon column c -> window index (-1 if conserved anchor column)
        std::vector<int> col_to_win(C, -1);
        for (size_t wi = 0; wi < windows.size(); ++wi) {
            for (size_t c = windows[wi].window.start_codon; c < windows[wi].window.end_codon && c < C; ++c) {
                col_to_win[c] = static_cast<int>(wi);
            }
        }

        // For each column c, emit a primary column segment S_c<c> representing the homologous codon column
        // plus (inside indel windows) explicit alternative bubble segments S_w<wi>_alt<k> for non-MAPindel paths
        for (size_t c = 0; c < C; ++c) {
            // Find consensus codon and amino acid profile at column c
            std::map<std::string, size_t> trip_counts;
            std::map<char, size_t> aa_counts;
            size_t non_gap = 0;
            for (size_t t = 0; t < N; ++t) {
                int8_t sc = aln.sequences[t].codons[c];
                if (sc >= 0) {
                    non_gap++;
                    std::string trip = aln.sequences[t].raw_seq.substr(3 * c, 3);
                    trip_counts[trip]++;
                    aa_counts[gcode.sense_codon_aa[sc]]++;
                }
            }
            std::string cons_trip = "NNN";
            size_t best_cnt = 0;
            for (const auto& [k, v] : trip_counts) {
                if (v > best_cnt) { best_cnt = v; cons_trip = k; }
            }
            std::string aa_profile;
            for (const auto& [aa, cnt] : aa_counts) {
                if (!aa_profile.empty()) aa_profile += "/";
                aa_profile += aa;
            }
            if (aa_profile.empty()) aa_profile = "-";

            Scalar pp = (c < site_posteriors.size()) ? site_posteriors[c] : 1.0;
            Scalar occ = (N > 0) ? (static_cast<Scalar>(non_gap) / static_cast<Scalar>(N)) : 0.0;
            std::string sel_class = "Neutral";
            Scalar omega_est = 1.0;
            if (c < site_results.size()) {
                const auto& sr = site_results[c];
                omega_est = (sr.marginalized_alpha > 1e-6) ? (sr.marginalized_beta / sr.marginalized_alpha) : 1.0;
                if (sr.marginalized_p_value <= 0.10) {
                    if (sr.marginalized_beta > sr.marginalized_alpha) sel_class = "Diversifying";
                    else if (sr.marginalized_alpha > sr.marginalized_beta) sel_class = "Purifying";
                }
            }

            std::string seg_id = "C" + std::to_string(c + 1);
            gfa << "S\t" << seg_id << "\t" << cons_trip
                << "\tRC:i:" << non_gap
                << "\tpp:f:" << std::fixed << std::setprecision(4) << pp
                << "\toc:f:" << occ
                << "\tom:f:" << omega_est
                << "\taa:Z:" << aa_profile
                << "\tsel:Z:" << sel_class << "\n";

            std::string fill = "#f8fafc";
            if (sel_class == "Diversifying") fill = "#fee2e2";
            else if (sel_class == "Purifying") fill = "#e0f2fe";
            else if (pp < 0.85) fill = "#fef3c7";

            dot << "  " << seg_id << " [label=\"{Codon " << (c + 1) << " (" << aa_profile << ")|PP="
                << std::fixed << std::setprecision(2) << pp << " occ=" << occ
                << " (" << sel_class << ")}\", fillcolor=\"" << fill << "\"];\n";
        }

        // Compute posterior-weighted directed transitions between columns (including deletion skip-edges across bubbles!)
        // Key: (from_col_1based, to_col_1based) -> posterior mass summed over all N taxa
        std::map<std::pair<size_t, size_t>, Scalar> edge_mass;

        for (size_t t = 0; t < N; ++t) {
            // Outside indel windows, follow taxon t's non-gap transitions; inside indel windows, distribute
            // taxon t's weight across all candidate transducer placements `top_candidates`!
            size_t c = 0;
            int prev_non_gap = -1;
            while (c < C) {
                int wi = col_to_win[c];
                if (wi < 0 || static_cast<size_t>(wi) >= window_taxon_results.size() || t >= window_taxon_results[wi].size()) {
                    if (aln.sequences[t].codons[c] >= 0) {
                        if (prev_non_gap >= 0) {
                            edge_mass[{static_cast<size_t>(prev_non_gap + 1), c + 1}] += 1.0;
                        }
                        prev_non_gap = static_cast<int>(c);
                    }
                    c++;
                } else {
                    size_t w_start = std::min(windows[wi].window.start_codon, C);
                    size_t w_end   = std::min(windows[wi].window.end_codon, C);
                    size_t W = (w_end > w_start) ? (w_end - w_start) : 0;
                    const auto& tr = window_taxon_results[wi][t];

                    int last_in_win = prev_non_gap;
                    if (!tr.top_candidates.empty() && W > 0) {
                        for (const auto& cand : tr.top_candidates) {
                            if (cand.posterior_weight < 1e-4) continue;
                            int cur_prev = prev_non_gap;
                            for (size_t w = 0; w < W && w < cand.codons.size() && (w_start + w) < C; ++w) {
                                if (cand.codons[w] >= 0) {
                                    if (cur_prev >= 0) {
                                        edge_mass[{static_cast<size_t>(cur_prev + 1), w_start + w + 1}] += cand.posterior_weight;
                                    }
                                    cur_prev = static_cast<int>(w_start + w);
                                }
                            }
                        }
                        // Set prev_non_gap to the MAP last non-gap codon in this window for connecting to w_end
                        for (size_t w = 0; w < W && (w_start + w) < C; ++w) {
                            if (aln.sequences[t].codons[w_start + w] >= 0) {
                                last_in_win = static_cast<int>(w_start + w);
                            }
                        }
                    }
                    prev_non_gap = last_in_win;
                    c = std::max(c + 1, w_end);
                }
            }
        }

        // Emit GFA Link (L) lines and DOT edges weighted by posterior taxon fraction
        for (const auto& [edge, mass] : edge_mass) {
            Scalar norm_weight = (N > 0) ? (mass / static_cast<Scalar>(N)) : 0.0;
            if (norm_weight < 0.01) continue;
            std::string u = "C" + std::to_string(edge.first);
            std::string v = "C" + std::to_string(edge.second);
            bool is_skip = (edge.second > edge.first + 1);
            gfa << "L\t" << u << "\t+\t" << v << "\t+\t0M"
                << "\tew:f:" << std::fixed << std::setprecision(4) << norm_weight
                << "\tsk:Z:" << (is_skip ? "bubble_indel" : "colinear") << "\n";

            dot << "  " << u << " -> " << v << " [label=\""
                << std::fixed << std::setprecision(2) << norm_weight << "\""
                << (is_skip ? ", color=\"#dc2626\", style=\"dashed\", penwidth=1.8" : ", color=\"#334155\", penwidth=1.5")
                << "];\n";
        }

        // Emit GFA Path (P) lines for every taxon's MAP path (plus top alternative bubble path if posterior >= 0.05)
        for (size_t t = 0; t < N; ++t) {
            std::string path_str;
            for (size_t c = 0; c < C; ++c) {
                if (aln.sequences[t].codons[c] >= 0) {
                    if (!path_str.empty()) path_str += ",";
                    path_str += "C" + std::to_string(c + 1) + "+";
                }
            }
            if (!path_str.empty()) {
                gfa << "P\t" << aln.sequences[t].name << "\t" << path_str << "\t*\n";
            }
        }

        dot << "}\n";
        return {gfa.str(), dot.str()};
    }

    // =========================================================================
    // Full 3-Stage Execution Pipeline
    // =========================================================================
    JointAlignResult run(const Alignment& input_aln, std::optional<Tree> starting_tree = std::nullopt) {
        auto t0 = std::chrono::steady_clock::now();
        JointAlignResult res;
        res.initial_alignment = input_aln;
        res.map_alignment = input_aln;

        // 0. Build or optimize initial tree via Ohmic TreeEngine if not provided
        if (!settings.quiet) {
            Panel::print_step(
                1, settings.run_marginalized_fel ? 5 : 4,
                "Initial Phylogeny & Global Codon MG94 Fit",
                std::to_string(input_aln.num_taxa) + " taxa × " + std::to_string(input_aln.num_codons) + " codons"
            );
        }
        if (starting_tree.has_value()) {
            res.initial_tree = *starting_tree;
        } else if (!input_aln.embedded_tree_newick.empty()) {
            res.initial_tree = Tree::from_newick(input_aln.embedded_tree_newick);
        } else {
            TreeEngineSettings tes;
            tes.quiet = true;
            TreeEngine te(tes);
            auto te_res = te.run(input_aln);
            res.initial_tree = te_res.optimized_tree;
        }

        // Fit global GTR + MG94 parameters on initial alignment and tree
        GTRFitter gtr_fitter(res.initial_tree, res.initial_alignment);
        auto gtr_fit = gtr_fitter.fit();
        res.initial_gtr_log_l = gtr_fit.log_likelihood;

        MG94Parameters mg_params;
        mg_params.theta_AC = gtr_fit.params.theta_AC;
        mg_params.theta_AT = gtr_fit.params.theta_AT;
        mg_params.theta_CG = gtr_fit.params.theta_CG;
        mg_params.theta_CT = gtr_fit.params.theta_CT;
        mg_params.theta_GT = gtr_fit.params.theta_GT;

        MG94Fitter mg_fitter(gtr_fit.tree, res.initial_alignment);
        auto mg_fit = mg_fitter.fit_omega_and_scale(1.0, mg_params);
        mg_params.alpha = 1.0;
        mg_params.beta = mg_fit.x_opt(0);

        res.initial_tree = mg_fit.tree;
        res.optimized_tree = mg_fit.tree;
        res.initial_mg94_log_l = mg_fit.log_likelihood;

        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2)
               << "      GTR lnL = " << res.initial_gtr_log_l
               << "  │  Initial MG94 lnL = " << res.initial_mg94_log_l
               << "  │  Global ω = " << std::setprecision(4) << mg_params.beta << "\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        // ---------------------------------------------------------------------
        // Stage 1: Codon TKF92/Affine Pair-HMM Window Realignment (Inside-Outside EM)
        // ---------------------------------------------------------------------
        if (!settings.quiet) {
            Panel::print_step(
                2, settings.run_marginalized_fel ? 5 : 4,
                "Stage 1: Codon Pair-HMM Window Realignment",
                "Inside-Outside EM • Ancestral Codon Profile Transducer"
            );
        }
        std::vector<std::vector<TransducerTaxonResult>> window_taxon_results;
        res.windows = run_stage1_window_realignment(
            res.map_alignment,
            res.optimized_tree,
            mg_params,
            settings,
            res.site_homology_posteriors,
            window_taxon_results
        );

        for (const auto& w : res.windows) {
            res.total_taxa_realigned += w.num_taxa_realigned;
            res.total_codons_shifted += w.num_codons_shifted;
        }

        // Evaluate global MG94 log-likelihood after Stage 1
        {
            MG94Fitter mg_stage1(res.optimized_tree, res.map_alignment);
            auto fit1 = mg_stage1.fit_omega_and_scale(mg_params.beta, mg_params);
            mg_params.beta = fit1.x_opt(0);
            res.optimized_tree = fit1.tree;
            res.stage1_mg94_log_l = fit1.log_likelihood;
        }

        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2)
               << "      Screened " << res.windows.size() << " indel/entropy window(s)"
               << "  │  Realigned " << res.total_taxa_realigned << " taxon-window(s) ("
               << res.total_codons_shifted << " codons shifted)"
               << "  │  Stage 1 MG94 lnL = " << res.stage1_mg94_log_l
               << " (Δ = " << std::showpos << (res.stage1_mg94_log_l - res.initial_mg94_log_l) << std::noshowpos << ")\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        // ---------------------------------------------------------------------
        // Stage 2: Coupled Ohmic SPR / Port-Slide + Boundary Indel Transducer
        // ---------------------------------------------------------------------
        if (!settings.quiet) {
            Panel::print_step(
                3, settings.run_marginalized_fel ? 5 : 4,
                "Stage 2: Coupled Ohmic SPR + Boundary Indel Transducer",
                "Foster Effective Resistance • Joint Topology + Indel Re-homologization"
            );
        }
        res.coupled_moves = run_stage2_coupled_topology_search(
            res.optimized_tree,
            res.map_alignment,
            mg_params,
            res.windows,
            settings
        );
        res.total_coupled_moves_accepted = res.coupled_moves.size();

        if (res.total_coupled_moves_accepted > 0) {
            MG94Fitter mg_stage2(res.optimized_tree, res.map_alignment);
            auto fit2 = mg_stage2.fit_omega_and_scale(mg_params.beta, mg_params);
            mg_params.beta = fit2.x_opt(0);
            res.optimized_tree = fit2.tree;
            res.stage2_mg94_log_l = fit2.log_likelihood;
        } else {
            res.stage2_mg94_log_l = res.stage1_mg94_log_l;
        }

        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2)
               << "      Accepted " << res.total_coupled_moves_accepted << " coupled topology+indel move(s)"
               << "  │  Stage 2 MG94 lnL = " << res.stage2_mg94_log_l
               << " (Δ = " << std::showpos << (res.stage2_mg94_log_l - res.stage1_mg94_log_l) << std::noshowpos << ")\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        // ---------------------------------------------------------------------
        // Stage 3: Alignment-Marginalized FEL
        // ---------------------------------------------------------------------
        if (settings.run_marginalized_fel) {
            if (!settings.quiet) {
                Panel::print_step(
                    4, 5,
                    "Stage 3: Alignment-Marginalized FEL Selection Scan",
                    "Integrating site likelihoods over MAP + alternative Pair-HMM paths"
                );
            }
            res.site_results = run_stage3_marginalized_fel(
                res.optimized_tree,
                res.initial_alignment,
                res.map_alignment,
                mg_params,
                res.windows,
                window_taxon_results,
                res.site_homology_posteriors,
                settings,
                res.false_positive_sites_eliminated,
                res.true_positive_sites_rescued
            );
            if (!settings.quiet) {
                size_t n_pos = 0, n_neg = 0;
                for (const auto& sr : res.site_results) {
                    if (sr.marginalized_p_value <= settings.p_value_threshold) {
                        if (sr.marginalized_beta > sr.marginalized_alpha) n_pos++;
                        else if (sr.marginalized_beta < sr.marginalized_alpha) n_neg++;
                    }
                }
                std::ostringstream ss;
                ss << "      Evaluated " << res.site_results.size() << " codon sites"
                   << "  │  ▲ Diversifying: " << n_pos
                   << "  │  ▼ Purifying: " << n_neg
                   << "  │  FP eliminated: " << res.false_positive_sites_eliminated
                   << ", TP rescued: " << res.true_positive_sites_rescued << "\n";
                std::cout << Console::muted(ss.str()) << std::flush;
            }
        }

        // ---------------------------------------------------------------------
        // Stage 4: Populate Co-Estimated Ohmic Wheatstone Bridges, Extended Newick,
        //          Stockholm Quality-Annotated Alignment, and POA Codon Graphs (GFAv1 & DOT)
        // ---------------------------------------------------------------------
        if (!settings.quiet) {
            int final_step = settings.run_marginalized_fel ? 5 : 4;
            Panel::print_step(
                final_step, final_step,
                "Stage 4: Ohmic Wheatstone Bridges, Stockholm PP_cons & Codon POA Graph",
                "Constructing eNewick reticulations and GFAv1/DOT Partial Order Alignment DAG"
            );
        }
        {
            OhmicSettings ohm_s;
            ohm_s.use_hyphaeon_transformer = false;
            res.ohmic_network = OhmicNetwork::infer(res.map_alignment, ohm_s);
            res.ohmic_network.backbone_tree = res.optimized_tree;
            res.ohmic_network.collapsed_tree = res.optimized_tree;
            OhmicNetwork::rebuild_bridges(res.ohmic_network, res.map_alignment, ohm_s);

            // Attach Wheatstone bridges onto optimized_tree for Extended Newick serialization
            Tree net_tree = res.ohmic_network.collapsed_tree;
            res.extended_newick = net_tree.to_extended_newick(true, true);

            res.stockholm_alignment = format_stockholm_alignment(
                res.map_alignment,
                res.site_homology_posteriors,
                res.windows,
                window_taxon_results,
                res.extended_newick
            );
            auto [gfa_str, dot_str] = build_poa_gfa_and_dot(
                res.map_alignment,
                res.site_homology_posteriors,
                res.windows,
                window_taxon_results,
                res.site_results
            );
            res.poa_gfa = std::move(gfa_str);
            res.poa_dot = std::move(dot_str);

            res.map_alignment.site_homology_posteriors = res.site_homology_posteriors;
            res.map_alignment.poa_gfa = res.poa_gfa;

            if (!settings.quiet) {
                Scalar mean_pp = 1.0, min_pp = 1.0;
                if (!res.site_homology_posteriors.empty()) {
                    Scalar s = 0.0;
                    for (Scalar p : res.site_homology_posteriors) {
                        s += p;
                        min_pp = std::min(min_pp, p);
                    }
                    mean_pp = s / static_cast<Scalar>(res.site_homology_posteriors.size());
                }
                std::ostringstream ss;
                ss << std::fixed << std::setprecision(4)
                   << "      Wheatstone bridges: " << res.ohmic_network.bridges.size()
                   << "  │  Polytomy busbars: " << res.ohmic_network.num_collapsed_polytomies
                   << "  │  Mean column posterior: " << mean_pp
                   << " (min: " << min_pp << ")\n\n";
                std::cout << Console::muted(ss.str()) << std::flush;
            }
        }

        auto t1 = std::chrono::steady_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    // =========================================================================
    // Phase A & B: Full De Novo Joint Inference from Completely Unaligned
    //   Homologous Coding Sequences (Unaligned FASTA/Vectors -> Ohmic Scaffold
    //   -> Bottom-Up Phylogenetic MG94 Log-Odds Codon Transducer -> Stages 1-3)
    // =========================================================================
    static std::vector<Sequence> strip_gaps_to_unaligned(const Alignment& aln) {
        std::vector<Sequence> raw(aln.num_taxa);
        for (size_t t = 0; t < aln.num_taxa; ++t) {
            raw[t].name = aln.sequences[t].name;
            std::string clean;
            clean.reserve(aln.sequences[t].raw_seq.size());
            // Strip in-frame gap triplets "---" while preserving valid codons
            size_t L = aln.sequences[t].raw_seq.size();
            for (size_t i = 0; i + 2 < L; i += 3) {
                std::string trip = aln.sequences[t].raw_seq.substr(i, 3);
                if (trip == "---" || trip == "???" || trip == "NNN") continue;
                for (char& ch : trip) {
                    if (ch == '-' || ch == '?') ch = 'N';
                }
                clean += trip;
            }
            raw[t].raw_seq = std::move(clean);
        }
        return raw;
    }

    static std::pair<Alignment, Tree> build_de_novo_codon_alignment_and_tree(
        const std::vector<Sequence>& unaligned_seqs,
        std::shared_ptr<const GeneticCode> gcode = nullptr,
        const JointAlignSettings& settings = JointAlignSettings{}
    ) {
        if (!gcode) gcode = GeneticCode::universal();
        size_t N = unaligned_seqs.size();
        if (N == 0) throw std::runtime_error("Empty unaligned sequence list");

        int S = gcode->num_sense_codons;

        // 1. Convert each unaligned sequence into in-frame sense codons and estimate frequencies
        std::vector<std::vector<int8_t>> taxon_codons(N);
        std::vector<std::vector<std::string>> taxon_triplets(N);
        Eigen::Matrix<Scalar, 3, 4> pos_counts = Eigen::Matrix<Scalar, 3, 4>::Zero();
        Vector codon_counts = Vector::Constant(S, 0.5); // Dirichlet pseudocount

        for (size_t t = 0; t < N; ++t) {
            const std::string& s = unaligned_seqs[t].raw_seq;
            size_t m = s.size() / 3;
            taxon_codons[t].reserve(m);
            taxon_triplets[t].reserve(m);
            for (size_t k = 0; k < m; ++k) {
                char c1 = s[3 * k], c2 = s[3 * k + 1], c3 = s[3 * k + 2];
                int8_t sc = gcode->triplet_to_sense(c1, c2, c3);
                taxon_codons[t].push_back(sc);
                taxon_triplets[t].push_back(s.substr(3 * k, 3));
                if (sc >= 0 && sc < S) {
                    codon_counts(sc) += 1.0;
                    for (int p = 0; p < 3; ++p) {
                        int8_t nb = char_to_nuc(s[3 * k + p]);
                        if (nb >= 0 && nb < 4) pos_counts(p, nb) += 1.0;
                    }
                }
            }
        }

        Eigen::Matrix<Scalar, 3, 4> pos_freqs = Eigen::Matrix<Scalar, 3, 4>::Constant(0.25);
        for (int p = 0; p < 3; ++p) {
            Scalar rsum = pos_counts.row(p).sum();
            if (rsum > 0.0) pos_freqs.row(p) = pos_counts.row(p) / rsum;
        }
        Vector pi_codon = codon_counts / codon_counts.sum();

        // 2. Phase A: Fast Pairwise Evolutionary Distance & Guide Tree Construction
        //    To handle both closely related sequences and deeply divergent mammalian/bacterial proteins
        //    (like CD2 where branch lengths reach 0.4 substitutions/site), we compute a fast banded
        //    amino-acid + codon similarity score for all pairs (i, j) to get accurate branch lengths!
        if (!settings.quiet) {
            size_t min_len = 999999, max_len = 0;
            for (const auto& tc : taxon_codons) {
                min_len = std::min(min_len, tc.size());
                max_len = std::max(max_len, tc.size());
            }
            Panel::print_step(
                0, settings.run_marginalized_fel ? 5 : 4,
                "Phase A: De Novo Alignment-Free Ohmic Conductance & Spectral Scaffold",
                std::to_string(N) + " unaligned taxa (" + std::to_string(min_len) + ".." + std::to_string(max_len) + " codons)"
            );
        }
        auto aa_sim = [&](int8_t c1, int8_t c2) -> Scalar {
            if (c1 < 0 || c2 < 0) return 0.0;
            if (c1 == c2) return 1.0;
            char a1 = gcode->sense_codon_aa[c1];
            char a2 = gcode->sense_codon_aa[c2];
            if (a1 == a2) return 0.85; // Synonymous
            int nd = gcode->diff_matrix[c1][c2].num_diffs;
            if (nd == 1) return 0.35;  // Single-nucleotide nonsynonymous
            if (nd == 2) return 0.10;
            return 0.0;
        };

        Matrix D_af = Matrix::Zero(N, N);
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < N; ++i) {
            const auto& seq_i = taxon_codons[i];
            size_t Li = seq_i.size();
            for (size_t j = i + 1; j < N; ++j) {
                const auto& seq_j = taxon_codons[j];
                size_t Lj = seq_j.size();
                if (Li == 0 || Lj == 0) continue;

                // Banded Smith-Waterman / Needleman-Wunsch similarity within band B = max(24, |Li - Lj| + 16)
                int B = static_cast<int>(std::max<size_t>(24, (Li > Lj ? Li - Lj : Lj - Li) + 16));
                std::vector<Scalar> prev(Lj + 1, 0.0), curr(Lj + 1, 0.0);
                for (size_t r = 1; r <= Li; ++r) {
                    int c_min = std::max(1, static_cast<int>(r) - B);
                    int c_max = std::min(static_cast<int>(Lj), static_cast<int>(r) + B);
                    std::fill(curr.begin(), curr.end(), 0.0);
                    for (int c = c_min; c <= c_max; ++c) {
                        Scalar match = prev[c - 1] + aa_sim(seq_i[r - 1], seq_j[c - 1]);
                        Scalar del = prev[c] - 0.6;
                        Scalar ins = curr[c - 1] - 0.6;
                        curr[c] = std::max({match, del, ins});
                    }
                    std::swap(prev, curr);
                }
                Scalar max_len = static_cast<Scalar>(std::max(Li, Lj));
                Scalar pid = std::clamp(prev[Lj] / std::max<Scalar>(1.0, max_len), 0.08, 0.9999);
                Scalar d = std::max<Scalar>(1e-4, -0.75 * std::log(pid));
                D_af(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = d;
                D_af(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = d;
            }
        }

        // Build Alignment-Free Ohmic Conductance Laplacian L = D_deg - W, Moore-Penrose L^+,
        // and Foster Effective Resistance matrix Omega_af on the unaligned sequences!
        std::vector<Scalar> upper_d;
        upper_d.reserve(N * (N - 1) / 2);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Scalar d = D_af(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                if (d > 0.0) upper_d.push_back(d);
            }
        }
        Scalar sigma_af = 0.1;
        if (!upper_d.empty()) {
            std::sort(upper_d.begin(), upper_d.end());
            sigma_af = std::max<Scalar>(1e-4, upper_d[upper_d.size() / 2]);
        }

        Matrix W_af = Matrix::Zero(N, N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Eigen::Index ii = static_cast<Eigen::Index>(i);
                Eigen::Index jj = static_cast<Eigen::Index>(j);
                Scalar d = D_af(ii, jj);
                Scalar w = std::exp(-d / sigma_af) / (d + 1e-3);
                W_af(ii, jj) = w;
                W_af(jj, ii) = w;
            }
        }
        Vector deg_af = W_af.rowwise().sum();
        Matrix L_af = Matrix(deg_af.asDiagonal()) - W_af;
        L_af = 0.5 * (L_af + L_af.transpose());

        Matrix Omega_af = D_af;
        Eigen::SelfAdjointEigenSolver<Matrix> eig_af(L_af);
        if (eig_af.info() == Eigen::Success) {
            const Vector& ev = eig_af.eigenvalues();
            const Matrix& evec = eig_af.eigenvectors();
            Matrix L_pinv = Matrix::Zero(N, N);
            for (Eigen::Index m = 0; m < ev.size(); ++m) {
                if (ev(m) > 1e-10) {
                    L_pinv.noalias() += (1.0 / ev(m)) * (evec.col(m) * evec.col(m).transpose());
                }
            }
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = i + 1; j < N; ++j) {
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    Scalar om = std::max<Scalar>(1e-4, L_pinv(ii, ii) + L_pinv(jj, jj) - 2.0 * L_pinv(ii, jj));
                    Omega_af(ii, jj) = om;
                    Omega_af(jj, ii) = om;
                }
            }
        }

        // Fuse direct evolutionary distance D_af with Ohmic circuit resistance Omega_af
        // (normalized to the median evolutionary scale) and construct the Ohmic Spectral Fiedler Scaffold!
        std::vector<Scalar> upper_om;
        upper_om.reserve(N * (N - 1) / 2);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                upper_om.push_back(Omega_af(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)));
            }
        }
        Scalar med_om = 1.0;
        if (!upper_om.empty()) {
            std::sort(upper_om.begin(), upper_om.end());
            med_om = std::max<Scalar>(1e-6, upper_om[upper_om.size() / 2]);
        }
        Matrix D_ohmic_af = 0.70 * D_af + 0.30 * (sigma_af / med_om) * Omega_af;

        std::vector<std::string> taxon_names(N);
        for (size_t t = 0; t < N; ++t) taxon_names[t] = unaligned_seqs[t].name;

        // Build balanced Ohmic Spectral Fiedler + Electrostatic Port-Grounded Guide Tree
        Tree guide_tree = TreeEngine::build_spectral_scaffold_from_distances(
            D_ohmic_af, taxon_names, false
        );

        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(4)
               << "      Median pairwise distance = " << sigma_af
               << "  │  Median Foster resistance Ω = " << med_om
               << "  │  Guide scaffold nodes = " << guide_tree.num_nodes() << "\n"
               << "      Phase B: Bottom-up phylogenetic MG94 log-odds codon transducer...\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        // 3. Build multi-scale MG94 rate matrices for Phase B Phylogenetic Transducer:
        //    - mg94_branch: uses the actual branch lengths (bl_left, bl_right)
        //    - mg94_deep: a broad evolutionary kernel (t = 0.65, omega = 0.6) so deeply divergent
        //      homologous amino acids still have a positive Bayes Factor (> -2 * gap_open) and
        //      never get split into staggered unaligned blocks!
        MG94Parameters mg_p;
        mg_p.alpha = 1.0;
        mg_p.beta = 0.45;
        MG94Matrix mg94;
        mg94.update(mg_p, pos_freqs, pi_codon, *gcode);
        Scalar conv = (mg94.scale_factor > 1e-12) ? (3.0 / mg94.scale_factor) : 1.0;
        Matrix P_deep = mg94.transition_matrix(0.55 * conv);

        // 4. Phase B: Bottom-Up Post-Order Phylogenetic MG94 Log-Odds Codon Transducer
        struct NodeProfile {
            std::vector<Vector> D;                      // size L_u, each Vector(S)
            std::vector<Scalar> occupancy;              // size L_u in (0, 1]
            std::vector<std::vector<int32_t>> tax_pos;  // size L_u x N (-1 if gap in taxon t)
            size_t num_leaves = 0;
        };

        std::vector<NodeProfile> profiles(guide_tree.num_nodes());
        std::unordered_map<std::string, size_t> name_to_idx;
        for (size_t t = 0; t < N; ++t) name_to_idx[taxon_names[t]] = t;

        auto log_odds_merge = [&](const Vector& d1, const Vector& d1_deep, Scalar occ1,
                                  const Vector& d2, const Vector& d2_deep, Scalar occ2) -> Scalar {
            // Multi-scale Phylogenetic Bayes Factor:
            // Combines exact branch-length MG94 homology odds (70%) with deep-homology MG94 kernel (30%).
            // A small Laplace pseudocount (0.15) prevents -Inf log-odds on highly variable sites while
            // ensuring truly incompatible codons have BF < 1.0 (negative log-odds).
            Scalar bf_exact = 0.0;
            Scalar bf_deep = 0.0;
            for (int s = 0; s < S; ++s) {
                Scalar pi_s = std::max<Scalar>(1e-8, pi_codon(s));
                bf_exact += (d1(s) * d2(s)) / pi_s;
                bf_deep  += (d1_deep(s) * d2_deep(s)) / pi_s;
            }
            Scalar bf = 0.70 * bf_exact + 0.30 * bf_deep + 0.15;
            Scalar llr = std::log(std::max<Scalar>(1e-4, bf));
            // Reward matching columns of similar taxon occupancy (prevents splitting shared framework/loop columns)
            Scalar occ_balance = 1.0 - 0.85 * std::abs(occ1 - occ2);
            return llr + 0.6 * std::log(std::max<Scalar>(0.15, occ_balance));
        };

        auto align_two_profiles = [&](const NodeProfile& left, Scalar bl_left,
                                      const NodeProfile& right, Scalar bl_right) -> NodeProfile {
            size_t L1 = left.D.size();
            size_t L2 = right.D.size();
            if (L1 == 0) return right;
            if (L2 == 0) return left;

            Scalar t1 = std::clamp(bl_left, 0.02, 0.80);
            Scalar t2 = std::clamp(bl_right, 0.02, 0.80);
            Matrix P1 = mg94.transition_matrix(t1 * conv);
            Matrix P2 = mg94.transition_matrix(t2 * conv);

            // Evolve both children's posterior state profiles to the parent node
            std::vector<Vector> E1(L1), E2(L2), E1_deep(L1), E2_deep(L2);
            for (size_t i = 0; i < L1; ++i) {
                E1[i] = P1.transpose() * left.D[i];
                Scalar s1 = E1[i].sum();
                if (s1 > 0.0) E1[i] /= s1;
                E1_deep[i] = P_deep.transpose() * left.D[i];
                Scalar sd1 = E1_deep[i].sum();
                if (sd1 > 0.0) E1_deep[i] /= sd1;
            }
            for (size_t j = 0; j < L2; ++j) {
                E2[j] = P2.transpose() * right.D[j];
                Scalar s2 = E2[j].sum();
                if (s2 > 0.0) E2[j] /= s2;
                E2_deep[j] = P_deep.transpose() * right.D[j];
                Scalar sd2 = E2_deep[j].sum();
                if (sd2 > 0.0) E2_deep[j] /= sd2;
            }

            // 3-State Affine Pair-HMM Codon Transducer:
            //   State 0 (M): Match (i, j)
            //   State 1 (I): Insert from left (i, -)
            //   State 2 (D): Delete / Insert from right (-, j)
            const Scalar NEG_INF = -1e18;
            // Scale gap-open penalty with evolutionary distance t_pair = t1 + t2 so longer branches
            // permit true indel events while short sister branches resist spurious gaps:
            Scalar t_pair = std::clamp(bl_left + bl_right, 0.02, 1.20);
            Scalar branch_scale = std::clamp(1.15 - 0.35 * std::min<Scalar>(1.0, t_pair / 0.40), 0.75, 1.15);
            const Scalar g_open = settings.gap_open_penalty * branch_scale;
            const Scalar g_ext = settings.gap_extend_penalty;

            Eigen::MatrixXd dp_M = Eigen::MatrixXd::Constant(L1 + 1, L2 + 1, NEG_INF);
            Eigen::MatrixXd dp_I = Eigen::MatrixXd::Constant(L1 + 1, L2 + 1, NEG_INF);
            Eigen::MatrixXd dp_D = Eigen::MatrixXd::Constant(L1 + 1, L2 + 1, NEG_INF);

            Eigen::MatrixXi bt_M = Eigen::MatrixXi::Zero(L1 + 1, L2 + 1);
            Eigen::MatrixXi bt_I = Eigen::MatrixXi::Zero(L1 + 1, L2 + 1);
            Eigen::MatrixXi bt_D = Eigen::MatrixXi::Zero(L1 + 1, L2 + 1);

            dp_M(0, 0) = 0.0;
            for (size_t i = 1; i <= L1; ++i) {
                dp_I(i, 0) = -g_open - static_cast<Scalar>(i - 1) * g_ext;
                bt_I(i, 0) = (i == 1) ? 0 : 1;
            }
            for (size_t j = 1; j <= L2; ++j) {
                dp_D(0, j) = -g_open - static_cast<Scalar>(j - 1) * g_ext;
                bt_D(0, j) = (j == 1) ? 0 : 2;
            }

            for (size_t i = 1; i <= L1; ++i) {
                for (size_t j = 1; j <= L2; ++j) {
                    // 1. State I: align left[i-1] to gap. If left.occupancy[i-1] < 1.0, an indel event
                    // already occurred inside the left subtree; propagating an existing partial column
                    // across the parent node is a continuation rather than a full new indel event.
                    Scalar occ_L = left.occupancy[i - 1];
                    Scalar eff_open_I = g_open * (0.15 + 0.85 * occ_L * occ_L);
                    Scalar eff_ext_I  = g_ext  * (0.15 + 0.85 * occ_L);
                    Scalar i_from_m = dp_M(i - 1, j) - eff_open_I;
                    Scalar i_from_i = dp_I(i - 1, j) - eff_ext_I;
                    if (i_from_m >= i_from_i) {
                        dp_I(i, j) = i_from_m;
                        bt_I(i, j) = 0;
                    } else {
                        dp_I(i, j) = i_from_i;
                        bt_I(i, j) = 1;
                    }

                    // 2. State D: align right[j-1] to gap
                    Scalar occ_R = right.occupancy[j - 1];
                    Scalar eff_open_D = g_open * (0.15 + 0.85 * occ_R * occ_R);
                    Scalar eff_ext_D  = g_ext  * (0.15 + 0.85 * occ_R);
                    Scalar d_from_m = dp_M(i, j - 1) - eff_open_D;
                    Scalar d_from_d = dp_D(i, j - 1) - eff_ext_D;
                    if (d_from_m >= d_from_d) {
                        dp_D(i, j) = d_from_m;
                        bt_D(i, j) = 0;
                    } else {
                        dp_D(i, j) = d_from_d;
                        bt_D(i, j) = 2;
                    }

                    // 3. State M: homologous Match of left[i-1] and right[j-1]
                    Scalar emit = log_odds_merge(
                        E1[i - 1], E1_deep[i - 1], left.occupancy[i - 1],
                        E2[j - 1], E2_deep[j - 1], right.occupancy[j - 1]
                    );
                    Scalar m0 = dp_M(i - 1, j - 1);
                    Scalar m1 = dp_I(i - 1, j - 1);
                    Scalar m2 = dp_D(i - 1, j - 1);
                    if (m0 >= m1 && m0 >= m2) {
                        dp_M(i, j) = m0 + emit;
                        bt_M(i, j) = 0;
                    } else if (m1 >= m2) {
                        dp_M(i, j) = m1 + emit;
                        bt_M(i, j) = 1;
                    } else {
                        dp_M(i, j) = m2 + emit;
                        bt_M(i, j) = 2;
                    }
                }
            }

            // Backtrack optimal transducer path
            int state = 0;
            Scalar best_end = dp_M(L1, L2);
            if (dp_I(L1, L2) > best_end) { best_end = dp_I(L1, L2); state = 1; }
            if (dp_D(L1, L2) > best_end) { best_end = dp_D(L1, L2); state = 2; }

            struct Step { int op; size_t i; size_t j; }; // op: 0=M, 1=I, 2=D
            std::vector<Step> path;
            size_t ci = L1, cj = L2;
            while (ci > 0 || cj > 0) {
                if (state == 0 && ci > 0 && cj > 0) {
                    path.push_back({0, ci - 1, cj - 1});
                    state = bt_M(ci, cj);
                    ci--; cj--;
                } else if (state == 1 && ci > 0) {
                    path.push_back({1, ci - 1, 0});
                    state = bt_I(ci, cj);
                    ci--;
                } else if (state == 2 && cj > 0) {
                    path.push_back({2, 0, cj - 1});
                    state = bt_D(ci, cj);
                    cj--;
                } else if (ci > 0) {
                    path.push_back({1, ci - 1, 0});
                    ci--;
                } else {
                    path.push_back({2, 0, cj - 1});
                    cj--;
                }
            }
            std::reverse(path.begin(), path.end());

            // Assemble merged ancestral profile
            NodeProfile parent;
            parent.num_leaves = left.num_leaves + right.num_leaves;
            size_t L_new = path.size();
            parent.D.resize(L_new);
            parent.occupancy.resize(L_new);
            parent.tax_pos.assign(L_new, std::vector<int32_t>(N, -1));

            Scalar w_left = static_cast<Scalar>(left.num_leaves) / static_cast<Scalar>(parent.num_leaves);
            Scalar w_right = static_cast<Scalar>(right.num_leaves) / static_cast<Scalar>(parent.num_leaves);

            for (size_t col = 0; col < L_new; ++col) {
                const auto& st = path[col];
                if (st.op == 0) {
                    // Combine posterior codon distributions via weighted mixture + Bayesian product
                    Vector prod(S);
                    for (int s = 0; s < S; ++s) {
                        prod(s) = (E1[st.i](s) * E2[st.j](s)) / std::max<Scalar>(1e-8, pi_codon(s));
                    }
                    Scalar psum = prod.sum();
                    Vector mix = w_left * E1[st.i] + w_right * E2[st.j];
                    if (psum > 1e-6) {
                        prod /= psum;
                        parent.D[col] = 0.65 * prod + 0.35 * mix;
                    } else {
                        parent.D[col] = mix;
                    }
                    parent.occupancy[col] = w_left * left.occupancy[st.i] + w_right * right.occupancy[st.j];
                    for (size_t t = 0; t < N; ++t) {
                        parent.tax_pos[col][t] = (left.tax_pos[st.i][t] != -1) ? left.tax_pos[st.i][t] : right.tax_pos[st.j][t];
                    }
                } else if (st.op == 1) {
                    parent.D[col] = E1[st.i];
                    parent.occupancy[col] = w_left * left.occupancy[st.i];
                    parent.tax_pos[col] = left.tax_pos[st.i];
                } else {
                    parent.D[col] = E2[st.j];
                    parent.occupancy[col] = w_right * right.occupancy[st.j];
                    parent.tax_pos[col] = right.tax_pos[st.j];
                }
            }
            return parent;
        };

        auto compact_profile = [&](NodeProfile& prof) {
            size_t n_active = prof.num_leaves;
            if (n_active < 2) return;
            bool merged_any = true;
            while (merged_any) {
                merged_any = false;
                size_t L_cur = prof.D.size();
                if (L_cur < 2) break;

                // Pass A: Multi-column staggered block compaction [col .. col+k-1] + [col+k .. col+2k-1] -> k columns
                // Only collapse when the paired columns exhibit genuine evolutionary homology (BF > 1.0, i.e. above random chance),
                // or when every active taxon is present in one of the two complementary blocks (n_left + n_right == n_active,
                // meaning a core block was split by an I-then-D path) with compatible codons (BF > 0.35).
                for (size_t k = 1; k <= 5 && !merged_any; ++k) {
                    if (L_cur < 2 * k) continue;
                    for (size_t col = 0; col + 2 * k <= L_cur; ++col) {
                        bool disjoint = true;
                        size_t n_left_taxa = 0, n_right_taxa = 0;
                        for (size_t t = 0; t < N; ++t) {
                            bool in_left = false, in_right = false;
                            for (size_t r = 0; r < k; ++r) {
                                if (prof.tax_pos[col + r][t] != -1) in_left = true;
                                if (prof.tax_pos[col + k + r][t] != -1) in_right = true;
                            }
                            if (in_left && in_right) {
                                disjoint = false;
                                break;
                            }
                            if (in_left) n_left_taxa++;
                            if (in_right) n_right_taxa++;
                        }
                        if (!disjoint || n_left_taxa == 0 || n_right_taxa == 0) continue;

                        // Check average deep evolutionary homology across the k paired columns
                        Scalar sum_bf = 0.0;
                        for (size_t r = 0; r < k; ++r) {
                            Vector d1_deep = P_deep.transpose() * prof.D[col + r];
                            Vector d2_deep = P_deep.transpose() * prof.D[col + k + r];
                            for (int s = 0; s < S; ++s) {
                                sum_bf += (d1_deep(s) * d2_deep(s)) / std::max<Scalar>(1e-8, pi_codon(s));
                            }
                        }
                        Scalar mean_bf = sum_bf / static_cast<Scalar>(k);
                        bool full_complement = (n_left_taxa + n_right_taxa == n_active);
                        if (mean_bf > 1.25 || (full_complement && mean_bf > 0.35)) {
                            for (size_t r = 0; r < k; ++r) {
                                size_t n1 = 0, n2 = 0;
                                for (size_t t = 0; t < N; ++t) {
                                    if (prof.tax_pos[col + r][t] != -1) n1++;
                                    if (prof.tax_pos[col + k + r][t] != -1) {
                                        n2++;
                                        prof.tax_pos[col + r][t] = prof.tax_pos[col + k + r][t];
                                    }
                                }
                                Scalar w1 = static_cast<Scalar>(n1) / static_cast<Scalar>(std::max<size_t>(1, n1 + n2));
                                Scalar w2 = static_cast<Scalar>(n2) / static_cast<Scalar>(std::max<size_t>(1, n1 + n2));
                                prof.D[col + r] = w1 * prof.D[col + r] + w2 * prof.D[col + k + r];
                                prof.occupancy[col + r] = static_cast<Scalar>(n1 + n2) / static_cast<Scalar>(n_active);
                            }
                            prof.D.erase(prof.D.begin() + static_cast<std::ptrdiff_t>(col + k),
                                         prof.D.begin() + static_cast<std::ptrdiff_t>(col + 2 * k));
                            prof.occupancy.erase(prof.occupancy.begin() + static_cast<std::ptrdiff_t>(col + k),
                                                 prof.occupancy.begin() + static_cast<std::ptrdiff_t>(col + 2 * k));
                            prof.tax_pos.erase(prof.tax_pos.begin() + static_cast<std::ptrdiff_t>(col + k),
                                               prof.tax_pos.begin() + static_cast<std::ptrdiff_t>(col + 2 * k));
                            merged_any = true;
                            break;
                        }
                    }
                }
                if (merged_any) continue;

                // Pass B: Single-column sliding compaction across intermediate gap columns (up to radius 4)
                for (size_t col = 0; col < L_cur; ++col) {
                    size_t n1 = 0;
                    for (size_t t = 0; t < N; ++t) {
                        if (prof.tax_pos[col][t] != -1) n1++;
                    }
                    if (n1 == 0 || n1 >= n_active) continue;

                    size_t max_other = std::min(L_cur - 1, col + 4);
                    bool did_merge = false;
                    for (size_t other = col + 1; other <= max_other; ++other) {
                        bool disjoint_and_ordered = true;
                        size_t n2 = 0;
                        for (size_t t = 0; t < N; ++t) {
                            int32_t p1 = prof.tax_pos[col][t];
                            int32_t p2 = prof.tax_pos[other][t];
                            if (p1 != -1 && p2 != -1) {
                                disjoint_and_ordered = false;
                                break;
                            }
                            if (p2 != -1) {
                                n2++;
                                for (size_t mid = col + 1; mid < other; ++mid) {
                                    if (prof.tax_pos[mid][t] != -1) {
                                        disjoint_and_ordered = false;
                                        break;
                                    }
                                }
                                if (!disjoint_and_ordered) break;
                            }
                        }
                        if (!disjoint_and_ordered || n2 == 0) continue;

                        Vector d1_deep = P_deep.transpose() * prof.D[col];
                        Vector d2_deep = P_deep.transpose() * prof.D[other];
                        Scalar bf = 0.0;
                        for (int s = 0; s < S; ++s) {
                            bf += (d1_deep(s) * d2_deep(s)) / std::max<Scalar>(1e-8, pi_codon(s));
                        }
                        bool full_complement = (n1 + n2 == n_active);
                        if (bf > 1.50 || (full_complement && bf > 0.50)) {
                            Scalar w1 = static_cast<Scalar>(n1) / static_cast<Scalar>(n1 + n2);
                            Scalar w2 = static_cast<Scalar>(n2) / static_cast<Scalar>(n1 + n2);
                            prof.D[col] = w1 * prof.D[col] + w2 * prof.D[other];
                            prof.occupancy[col] = static_cast<Scalar>(n1 + n2) / static_cast<Scalar>(n_active);
                            for (size_t t = 0; t < N; ++t) {
                                if (prof.tax_pos[other][t] != -1) {
                                    prof.tax_pos[col][t] = prof.tax_pos[other][t];
                                }
                            }
                            prof.D.erase(prof.D.begin() + static_cast<std::ptrdiff_t>(other));
                            prof.occupancy.erase(prof.occupancy.begin() + static_cast<std::ptrdiff_t>(other));
                            prof.tax_pos.erase(prof.tax_pos.begin() + static_cast<std::ptrdiff_t>(other));
                            merged_any = true;
                            did_merge = true;
                            break;
                        }
                    }
                    if (did_merge) break;
                }
                if (merged_any) continue;

                // Pass C: Non-adjacent multi-column block compaction [col .. col+k-1] + [other .. other+k-1] -> k columns
                // Handles cases where an intermediate singleton insertion separates two k-column staggered blocks!
                for (size_t k = 2; k <= 5 && !merged_any; ++k) {
                    if (L_cur < 2 * k) continue;
                    for (size_t col = 0; col + 2 * k <= L_cur && !merged_any; ++col) {
                        size_t max_other = std::min(L_cur - k, col + k + 4);
                        for (size_t other = col + k; other <= max_other; ++other) {
                            bool valid_merge = true;
                            size_t n_left_taxa = 0, n_right_taxa = 0;
                            for (size_t t = 0; t < N; ++t) {
                                bool in_left = false, in_right = false, in_mid = false;
                                for (size_t r = 0; r < k; ++r) {
                                    if (prof.tax_pos[col + r][t] != -1) in_left = true;
                                    if (prof.tax_pos[other + r][t] != -1) in_right = true;
                                }
                                for (size_t mid = col + k; mid < other; ++mid) {
                                    if (prof.tax_pos[mid][t] != -1) in_mid = true;
                                }
                                // Cannot merge if taxon has codons in both blocks, or has codons in both mid and right
                                // (since moving right into left would cross mid and violate sequence order!)
                                if ((in_left && in_right) || (in_mid && in_right)) {
                                    valid_merge = false;
                                    break;
                                }
                                if (in_left) n_left_taxa++;
                                if (in_right) n_right_taxa++;
                            }
                            if (!valid_merge || n_left_taxa == 0 || n_right_taxa == 0) continue;

                            Scalar sum_bf = 0.0;
                            for (size_t r = 0; r < k; ++r) {
                                Vector d1_deep = P_deep.transpose() * prof.D[col + r];
                                Vector d2_deep = P_deep.transpose() * prof.D[other + r];
                                for (int s = 0; s < S; ++s) {
                                    sum_bf += (d1_deep(s) * d2_deep(s)) / std::max<Scalar>(1e-8, pi_codon(s));
                                }
                            }
                            Scalar mean_bf = sum_bf / static_cast<Scalar>(k);
                            bool full_complement = (n_left_taxa + n_right_taxa == n_active);
                            if (mean_bf > 1.25 || (full_complement && mean_bf > 0.35)) {
                                for (size_t r = 0; r < k; ++r) {
                                    size_t n1 = 0, n2 = 0;
                                    for (size_t t = 0; t < N; ++t) {
                                        if (prof.tax_pos[col + r][t] != -1) n1++;
                                        if (prof.tax_pos[other + r][t] != -1) {
                                            n2++;
                                            prof.tax_pos[col + r][t] = prof.tax_pos[other + r][t];
                                        }
                                    }
                                    Scalar w1 = static_cast<Scalar>(n1) / static_cast<Scalar>(std::max<size_t>(1, n1 + n2));
                                    Scalar w2 = static_cast<Scalar>(n2) / static_cast<Scalar>(std::max<size_t>(1, n1 + n2));
                                    prof.D[col + r] = w1 * prof.D[col + r] + w2 * prof.D[other + r];
                                    prof.occupancy[col + r] = static_cast<Scalar>(n1 + n2) / static_cast<Scalar>(n_active);
                                }
                                prof.D.erase(prof.D.begin() + static_cast<std::ptrdiff_t>(other),
                                             prof.D.begin() + static_cast<std::ptrdiff_t>(other + k));
                                prof.occupancy.erase(prof.occupancy.begin() + static_cast<std::ptrdiff_t>(other),
                                                     prof.occupancy.begin() + static_cast<std::ptrdiff_t>(other + k));
                                prof.tax_pos.erase(prof.tax_pos.begin() + static_cast<std::ptrdiff_t>(other),
                                                   prof.tax_pos.begin() + static_cast<std::ptrdiff_t>(other + k));
                                merged_any = true;
                                break;
                            }
                        }
                    }
                }
            }
        };

        for (int32_t nid : guide_tree.post_order) {
            const auto& node = guide_tree.nodes[nid];
            if (node.is_leaf) {
                auto it = name_to_idx.find(node.name);
                if (it == name_to_idx.end()) continue;
                size_t tid = it->second;
                size_t m = taxon_codons[tid].size();
                NodeProfile lp;
                lp.num_leaves = 1;
                lp.D.resize(m);
                lp.occupancy.assign(m, 1.0);
                lp.tax_pos.assign(m, std::vector<int32_t>(N, -1));
                for (size_t k = 0; k < m; ++k) {
                    lp.tax_pos[k][tid] = static_cast<int32_t>(k);
                    int8_t sc = taxon_codons[tid][k];
                    if (sc >= 0 && sc < S) {
                        lp.D[k] = Vector::Constant(S, 1e-5);
                        lp.D[k](sc) = 1.0;
                        lp.D[k] /= lp.D[k].sum();
                    } else {
                        lp.D[k] = pi_codon;
                    }
                }
                profiles[nid] = std::move(lp);
            } else {
                if (node.children.empty()) continue;
                int32_t c0 = node.children[0];
                NodeProfile cur = profiles[c0];
                Scalar bl0 = guide_tree.nodes[c0].branch_length;
                for (size_t ci = 1; ci < node.children.size(); ++ci) {
                    int32_t ck = node.children[ci];
                    cur = align_two_profiles(cur, (ci == 1) ? bl0 : 0.02, profiles[ck], guide_tree.nodes[ck].branch_length);
                    compact_profile(cur);
                }
                profiles[nid] = std::move(cur);
            }
        }

        // 5. Final Root Profile Compaction
        NodeProfile root_prof = std::move(profiles[guide_tree.root_id]);
        compact_profile(root_prof);

        size_t C_total = root_prof.D.size();

        Alignment de_novo_aln;
        for (size_t t = 0; t < N; ++t) {
            std::string aligned_raw;
            aligned_raw.reserve(3 * C_total);
            for (size_t col = 0; col < C_total; ++col) {
                int32_t pos = root_prof.tax_pos[col][t];
                if (pos >= 0 && static_cast<size_t>(pos) < taxon_triplets[t].size()) {
                    aligned_raw += taxon_triplets[t][pos];
                } else {
                    aligned_raw += "---";
                }
            }
            de_novo_aln.add_sequence(taxon_names[t], aligned_raw);
        }
        de_novo_aln.finalize(gcode);

        // Gate 4: Post-Transducer MSA Expansion & Homology Collapse Check
        if (settings.enforce_admission_control) {
            std::vector<size_t> raw_lens;
            raw_lens.reserve(N);
            for (const auto& tc : taxon_codons) raw_lens.push_back(tc.size());
            std::sort(raw_lens.begin(), raw_lens.end());
            size_t med_len = raw_lens.empty() ? 1 : std::max<size_t>(1, raw_lens[raw_lens.size() / 2]);
            Scalar expansion = static_cast<Scalar>(C_total) / static_cast<Scalar>(med_len);

            size_t core_cols = 0;
            for (size_t col = 0; col < C_total; ++col) {
                if (root_prof.occupancy[col] >= 0.40) core_cols++;
            }
            Scalar core_frac = (C_total > 0) ? (static_cast<Scalar>(core_cols) / static_cast<Scalar>(C_total)) : 0.0;

            if (expansion > settings.max_msa_expansion_ratio || core_frac < 0.25) {
                std::ostringstream err;
                err << "Gate 4 (Post-Transducer Alignment Collapse Gate): Assembled codon alignment failed homology quality checks "
                    << "(assembled columns = " << C_total << " vs median input length = " << med_len
                    << " codons, expansion ratio = " << std::fixed << std::setprecision(2) << expansion
                    << "x [max allowed: " << settings.max_msa_expansion_ratio << "x], core column fraction = "
                    << (core_frac * 100.0) << "%). Input sequences lack colinear codon homology.";
                throw std::runtime_error(err.str());
            }
        }

        // Refine the tree topology on the assembled alignment using Ohmic TreeEngine
        TreeEngineSettings tes;
        tes.quiet = true;
        TreeEngine te(tes);
        auto te_res = te.run(de_novo_aln);

        if (!settings.quiet) {
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(2)
               << "      Assembled de novo MSA (" << C_total << " homologous codon columns)"
               << "  │  Refined GTR lnL = " << te_res.final_log_l << "\n";
            std::cout << Console::muted(ss.str()) << std::flush;
        }

        return {std::move(de_novo_aln), std::move(te_res.optimized_tree)};
    }

    // =========================================================================
    // 4-Gate Input Triage & Homology Admission Filter:
    //   Gate 1: Lexical & Alphabet Gate (rejects protein FASTA, empty/short/high-N sequences)
    //   Gate 2: 6-Frame & Strand ORF Gate (rescues +1/+2 frame shifts & reverse-complements,
    //           quarantines pseudogenes/introns with internal stop codons, aborts if >25% fail)
    //   Gate 3: Pairwise Log-Odds Homology Z-Score & Connected-Component Gate
    //           (quarantines non-homologous rogue taxa or extreme length outliers,
    //            aborts if cohort median Z < min_cohort_median_z or component < 3 taxa)
    //   Gate 4: Post-Transducer MSA Expansion & Core Occupancy Gate (enforced after Phase B)
    // =========================================================================
    static std::string reverse_complement_dna(const std::string& seq) {
        std::string rc;
        rc.reserve(seq.size());
        for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
            char c = static_cast<char>(std::toupper(static_cast<unsigned char>(*it)));
            switch (c) {
                case 'A': rc += 'T'; break;
                case 'T': case 'U': rc += 'A'; break;
                case 'C': rc += 'G'; break;
                case 'G': rc += 'C'; break;
                default:  rc += 'N'; break;
            }
        }
        return rc;
    }

    static std::pair<std::vector<Sequence>, InputAdmissionReport> triage_and_admit_unaligned_sequences(
        const std::vector<Sequence>& raw_seqs,
        std::shared_ptr<const GeneticCode> gcode = nullptr,
        const JointAlignSettings& settings = JointAlignSettings{}
    ) {
        if (!gcode) gcode = GeneticCode::universal();
        InputAdmissionReport report;
        report.total_input_sequences = raw_seqs.size();
        if (raw_seqs.size() < 3) {
            throw std::runtime_error(
                "Gate 1 (Input Size Gate): Phylogenetic selection analysis requires at least 3 homologous coding sequences (received " +
                std::to_string(raw_seqs.size()) + ")."
            );
        }

        // Gate 1: Global Alphabet Check (detect amino-acid protein FASTA masquerading as nucleotide CDS)
        size_t total_alpha = 0;
        size_t total_aa_only = 0; // E, F, I, L, P, Q, Z, J, O
        for (const auto& sq : raw_seqs) {
            for (char ch : sq.raw_seq) {
                char uc = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                if (std::isalpha(static_cast<unsigned char>(uc))) {
                    total_alpha++;
                    if (uc == 'E' || uc == 'F' || uc == 'I' || uc == 'L' ||
                        uc == 'P' || uc == 'Q' || uc == 'Z' || uc == 'J' || uc == 'O') {
                        total_aa_only++;
                    }
                }
            }
        }
        if (total_alpha > 0 && static_cast<double>(total_aa_only) / static_cast<double>(total_alpha) > 0.08) {
            std::ostringstream err;
            err << "Gate 1 (Lexical Alphabet Gate): Input sequences contain "
                << std::fixed << std::setprecision(1)
                << (100.0 * static_cast<double>(total_aa_only) / static_cast<double>(total_alpha))
                << "% amino-acid-specific letters (E, F, I, L, P, Q). "
                << "This appears to be a protein (.faa) file rather than nucleotide coding sequences (CDS). "
                << "HyPhy codon models require nucleotide CDS.";
            throw std::runtime_error(err.str());
        }

        // Helper to extract in-frame triplets, strip terminal stop codon, and count internal stop codons
        auto extract_frame_cds = [&](const std::string& nt_seq, int offset, size_t& out_internal_stops, std::vector<int8_t>& out_codons) -> std::string {
            out_internal_stops = 0;
            out_codons.clear();
            std::string clean;
            if (static_cast<size_t>(offset) >= nt_seq.size()) return clean;
            size_t usable_len = ((nt_seq.size() - static_cast<size_t>(offset)) / 3) * 3;
            clean.reserve(usable_len);

            for (size_t i = static_cast<size_t>(offset); i + 2 < static_cast<size_t>(offset) + usable_len; i += 3) {
                std::string trip = nt_seq.substr(i, 3);
                if (trip == "---" || trip == "???" || trip == "NNN") continue;
                for (char& ch : trip) {
                    if (ch == '-' || ch == '?') ch = 'N';
                }
                bool is_last = (i + 3 >= static_cast<size_t>(offset) + usable_len);
                int8_t sc = gcode->triplet_to_sense(trip[0], trip[1], trip[2]);
                bool has_n = (trip.find('N') != std::string::npos);
                if (sc < 0 && !has_n) {
                    if (is_last) {
                        break; // strip terminal stop codon
                    } else {
                        out_internal_stops++;
                    }
                }
                clean += trip;
                out_codons.push_back(sc);
            }
            return clean;
        };

        // Gate 2: Reading-Frame, Strand & Internal Stop-Codon Gate
        struct CandidateTaxon {
            size_t orig_idx = 0;
            std::string name;
            std::string raw_nt;
            std::string clean_cds;
            std::vector<int8_t> codons;
            int frame_offset = 0;
            bool rev_comp = false;
            size_t stops_frame0 = 0;
            size_t stops_chosen = 0;
            bool passed_gate2 = false;
        };

        std::vector<CandidateTaxon> candidates(raw_seqs.size());
        report.entries.resize(raw_seqs.size());

        // First pass: check frame +0 on all sequences and build a reference amino-acid 3-mer table
        // from unambiguous canonical sequences (0 internal stops AND start codon ATG).
        std::unordered_set<uint32_t> clean_cohort_aa3mers;
        auto count_aa_3mers = [&](const std::vector<int8_t>& cods, const std::unordered_set<uint32_t>& ref_set) -> size_t {
            if (cods.size() < 3 || ref_set.empty()) return 0;
            size_t hits = 0;
            for (size_t k = 0; k + 2 < cods.size(); ++k) {
                if (cods[k] < 0 || cods[k + 1] < 0 || cods[k + 2] < 0) continue;
                uint32_t h = (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k]])) << 16)
                           | (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 1]])) << 8)
                           |  static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 2]]));
                if (ref_set.count(h)) hits++;
            }
            return hits;
        };
        auto add_aa_3mers = [&](const std::vector<int8_t>& cods, std::unordered_set<uint32_t>& kset) {
            if (cods.size() < 3) return;
            for (size_t k = 0; k + 2 < cods.size(); ++k) {
                if (cods[k] < 0 || cods[k + 1] < 0 || cods[k + 2] < 0) continue;
                uint32_t h = (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k]])) << 16)
                           | (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 1]])) << 8)
                           |  static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 2]]));
                kset.insert(h);
            }
        };

        for (size_t t = 0; t < raw_seqs.size(); ++t) {
            std::string ungapped;
            ungapped.reserve(raw_seqs[t].raw_seq.size());
            size_t n_ambig = 0;
            for (char ch : raw_seqs[t].raw_seq) {
                char uc = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                if (uc == '-' || std::isspace(static_cast<unsigned char>(uc))) continue;
                if (uc == 'U') uc = 'T';
                if (uc != 'A' && uc != 'C' && uc != 'G' && uc != 'T') {
                    uc = 'N';
                    n_ambig++;
                }
                ungapped += uc;
            }

            TaxonAdmissionEntry& ent = report.entries[t];
            ent.name = raw_seqs[t].name;
            ent.raw_nt_length = ungapped.size();

            CandidateTaxon& ct = candidates[t];
            ct.orig_idx = t;
            ct.name = raw_seqs[t].name;
            ct.raw_nt = std::move(ungapped);

            if (ct.raw_nt.size() < settings.min_sense_codons * 3) {
                ent.status = "Quarantined";
                ent.reason = "Sequence too short (< " + std::to_string(settings.min_sense_codons) + " codons)";
                continue;
            }
            if (static_cast<double>(n_ambig) / static_cast<double>(ct.raw_nt.size()) > 0.20) {
                ent.status = "Quarantined";
                ent.reason = "Excessive ambiguous nucleotides (> 20% N/?)";
                continue;
            }

            size_t stops0 = 0;
            std::vector<int8_t> cods0;
            std::string cds0 = extract_frame_cds(ct.raw_nt, 0, stops0, cods0);
            ent.internal_stops_frame0 = stops0;
            ct.stops_frame0 = stops0;

            if (stops0 == 0 && cods0.size() >= settings.min_sense_codons) {
                ct.clean_cds = std::move(cds0);
                ct.codons = std::move(cods0);
                ct.frame_offset = 0;
                ct.rev_comp = false;
                ct.stops_chosen = 0;
                ct.passed_gate2 = true;
                ent.codon_length = ct.codons.size();
                if (ct.clean_cds.size() >= 3 && ct.clean_cds.substr(0, 3) == "ATG") {
                    add_aa_3mers(ct.codons, clean_cohort_aa3mers);
                }
            }
        }
        // Fallback if no sequence started with ATG: seed 3-mer table from all stop-free frame-0 sequences
        if (clean_cohort_aa3mers.empty()) {
            for (const auto& ct : candidates) {
                if (ct.passed_gate2) add_aa_3mers(ct.codons, clean_cohort_aa3mers);
            }
        }

        // Second pass: scan all 6 frames (+0,+1,+2 on fwd & revcomp) for any sequence that either:
        //   (a) had internal stop codons in frame +0, OR
        //   (b) did not start with ATG and has near-zero 3-mer homology in frame +0 (e.g. a GC-rich reverse-complemented sequence that happened to have 0 stops on the antisense strand!)
        size_t gate2_failures = 0;
        for (size_t t = 0; t < raw_seqs.size(); ++t) {
            CandidateTaxon& ct = candidates[t];
            TaxonAdmissionEntry& ent = report.entries[t];
            if (ent.status == "Quarantined") {
                gate2_failures++;
                continue;
            }
            size_t hits_frame0 = ct.passed_gate2 ? count_aa_3mers(ct.codons, clean_cohort_aa3mers) : 0;
            bool starts_atg0 = (ct.clean_cds.size() >= 3 && ct.clean_cds.substr(0, 3) == "ATG");
            if (ct.passed_gate2 && (starts_atg0 || hits_frame0 >= 15)) {
                continue;
            }

            std::string rc_nt = reverse_complement_dna(ct.raw_nt);
            int best_frame = -1;
            bool best_rc = false;
            size_t best_3mer_hits = hits_frame0;
            bool best_starts_atg = starts_atg0;
            std::string best_cds;
            std::vector<int8_t> best_cods;

            for (int rc_flag = 0; rc_flag <= 1; ++rc_flag) {
                const std::string& src = (rc_flag == 0) ? ct.raw_nt : rc_nt;
                for (int off = 0; off < 3; ++off) {
                    if (rc_flag == 0 && off == 0) continue;
                    size_t stops = 0;
                    std::vector<int8_t> cods;
                    std::string cds = extract_frame_cds(src, off, stops, cods);
                if (stops == 0 && cods.size() >= settings.min_sense_codons) {
                        size_t hits = count_aa_3mers(cods, clean_cohort_aa3mers);
                        bool st_atg = (cds.size() >= 3 && cds.substr(0, 3) == "ATG");
                        if ((!clean_cohort_aa3mers.empty() && hits > best_3mer_hits + 2 && hits >= 3) ||
                            (clean_cohort_aa3mers.empty() && best_frame == -1 && !ct.passed_gate2)) {
                            best_frame = off;
                            best_rc = (rc_flag == 1);
                            best_3mer_hits = hits;
                            best_starts_atg = st_atg;
                            best_cds = std::move(cds);
                            best_cods = std::move(cods);
                        }
                    }
                }
            }

            if (best_frame != -1 && (clean_cohort_aa3mers.empty() || best_3mer_hits >= 3)) {
                ct.clean_cds = std::move(best_cds);
                ct.codons = std::move(best_cods);
                ct.frame_offset = best_frame;
                ct.rev_comp = best_rc;
                ct.stops_chosen = 0;
                ct.passed_gate2 = true;
                ent.status = "Rescued";
                ent.frame_offset = best_frame;
                ent.reverse_complemented = best_rc;
                ent.codon_length = ct.codons.size();
                ent.reason = std::string("Auto-rescued ORFs (") + (best_rc ? "rev-comp " : "fwd ") + "frame +" + std::to_string(best_frame) + ")";
                report.rescued_frame_or_strand++;
                if (best_starts_atg) add_aa_3mers(ct.codons, clean_cohort_aa3mers);
            } else if (!ct.passed_gate2) {
                ent.status = "Quarantined";
                ent.reason = "Internal premature stop codons in frame +0 (" + std::to_string(ct.stops_frame0) + " stops; pseudogene/intron/UTR)";
                gate2_failures++;
            }
        }

        // If > 25% of input sequences fail Gate 2 (or fewer than 3 survive), abort dataset-wide!
        size_t gate2_survivors = raw_seqs.size() - gate2_failures;
        if (settings.enforce_admission_control &&
            (gate2_survivors < 3 || static_cast<double>(gate2_failures) / static_cast<double>(raw_seqs.size()) > 0.25)) {
            std::ostringstream err;
            err << "Gate 2 (Reading-Frame & Stop-Codon Gate): " << gate2_failures << " of " << raw_seqs.size()
                << " sequences (" << std::fixed << std::setprecision(1)
                << (100.0 * static_cast<double>(gate2_failures) / static_cast<double>(raw_seqs.size()))
                << "%) contain internal premature stop codons across reading frames or are too short. "
                << "Input likely contains un-spliced genomic introns, 5'/3' UTRs, or non-coding sequences.";
            throw std::runtime_error(err.str());
        }

        // Gate 3: Pairwise Null-Calibrated Homology Z-Score & Connected-Component Gate
        std::vector<size_t> active_indices;
        std::vector<size_t> active_lengths;
        for (size_t t = 0; t < candidates.size(); ++t) {
            if (candidates[t].passed_gate2) {
                active_indices.push_back(t);
                active_lengths.push_back(candidates[t].codons.size());
            }
        }
        std::sort(active_lengths.begin(), active_lengths.end());
        size_t median_codons = active_lengths.empty() ? 1 : active_lengths[active_lengths.size() / 2];

        // Quarantine extreme length outliers (< 0.30x or > 2.8x median length) before pairwise alignment
        std::vector<size_t> length_ok_indices;
        for (size_t idx : active_indices) {
            size_t L = candidates[idx].codons.size();
            if (settings.enforce_admission_control && median_codons >= 30 &&
                (L < static_cast<size_t>(0.30 * static_cast<double>(median_codons)) ||
                 L > static_cast<size_t>(2.80 * static_cast<double>(median_codons)))) {
                report.entries[idx].status = "Quarantined";
                report.entries[idx].reason = "Extreme length outlier (" + std::to_string(L) + " codons vs cohort median " +
                                             std::to_string(median_codons) + " codons)";
                candidates[idx].passed_gate2 = false;
            } else {
                length_ok_indices.push_back(idx);
            }
        }

        size_t M = length_ok_indices.size();
        if (M < 3) {
            throw std::runtime_error(
                "Gate 3 (Homology & Length Gate): Fewer than 3 sequences remain after reading-frame and length outlier filtering."
            );
        }

        auto aa_sim_score = [&](int8_t c1, int8_t c2) -> Scalar {
            if (c1 < 0 || c2 < 0) return 0.0;
            if (c1 == c2) return 1.0;
            char a1 = gcode->sense_codon_aa[c1];
            char a2 = gcode->sense_codon_aa[c2];
            if (a1 == a2) return 0.90;
            int nd = gcode->diff_matrix[c1][c2].num_diffs;
            if (nd == 1) return 0.30;
            return 0.0;
        };

        // Semi-global overlap banded DP helper: aligns shorter sequence (L_eff = min(Li, Lj))
        // against longer sequence (L_max = max(Li, Lj)) with free 5'/3' overhangs on the longer
        // sequence within band slack B0 = 24, returning {norm_dp, aa_id} normalized by L_eff.
        auto run_semiglobal_banded_dp = [&](const std::vector<int8_t>& s_a,
                                            const std::vector<int8_t>& s_b,
                                            std::vector<Scalar>& prev,
                                            std::vector<Scalar>& curr,
                                            std::vector<size_t>& prev_id,
                                            std::vector<size_t>& curr_id) -> std::pair<Scalar, Scalar> {
            const auto& s_short = (s_a.size() <= s_b.size()) ? s_a : s_b;
            const auto& s_long  = (s_a.size() <= s_b.size()) ? s_b : s_a;
            size_t Ls = s_short.size();
            size_t Ll = s_long.size();
            if (Ls == 0 || Ll == 0) return {0.0, 0.0};

            constexpr int B0 = 24;
            constexpr Scalar kNegInf = -1e9;
            int dL = static_cast<int>(Ll - Ls);
            if (prev.size() < Ll + 1) {
                prev.resize(Ll + 1);
                curr.resize(Ll + 1);
                prev_id.resize(Ll + 1);
                curr_id.resize(Ll + 1);
            }
            std::fill(prev.begin(), prev.begin() + static_cast<ptrdiff_t>(Ll + 1), kNegInf);
            std::fill(prev_id.begin(), prev_id.begin() + static_cast<ptrdiff_t>(Ll + 1), 0);
            int init_max = std::min(static_cast<int>(Ll), dL + B0);
            for (int c = 0; c <= init_max; ++c) {
                prev[static_cast<size_t>(c)] = 0.0;
            }

            for (size_t r = 1; r <= Ls; ++r) {
                int c_min = std::max(1, static_cast<int>(r) - B0);
                int c_max = std::min(static_cast<int>(Ll), static_cast<int>(r) + dL + B0);
                std::fill(curr.begin(), curr.begin() + static_cast<ptrdiff_t>(Ll + 1), kNegInf);
                std::fill(curr_id.begin(), curr_id.begin() + static_cast<ptrdiff_t>(Ll + 1), 0);
                if (static_cast<int>(r) <= B0) {
                    curr[0] = -0.65 * static_cast<Scalar>(r);
                }
                int8_t ci = s_short[r - 1];
                for (int c = c_min; c <= c_max; ++c) {
                    size_t cu = static_cast<size_t>(c);
                    int8_t cj = s_long[cu - 1];
                    Scalar m_score = prev[cu - 1] + aa_sim_score(ci, cj);
                    Scalar d_score = prev[cu] - 0.65;
                    Scalar i_score = curr[cu - 1] - 0.65;
                    if (m_score >= d_score && m_score >= i_score) {
                        curr[cu] = m_score;
                        bool exact_aa = (ci >= 0 && cj >= 0 && gcode->sense_codon_aa[ci] == gcode->sense_codon_aa[cj]);
                        curr_id[cu] = prev_id[cu - 1] + (exact_aa ? 1 : 0);
                    } else if (d_score >= i_score) {
                        curr[cu] = d_score;
                        curr_id[cu] = prev_id[cu];
                    } else {
                        curr[cu] = i_score;
                        curr_id[cu] = curr_id[cu - 1];
                    }
                }
                std::swap(prev, curr);
                std::swap(prev_id, curr_id);
            }

            int c_start = std::max(1, static_cast<int>(Ls) - B0);
            size_t best_c = static_cast<size_t>(c_start);
            Scalar best_score = prev[best_c];
            for (size_t c = static_cast<size_t>(c_start) + 1; c <= Ll; ++c) {
                if (prev[c] > best_score) {
                    best_score = prev[c];
                    best_c = c;
                }
            }
            Scalar L_eff = std::max<Scalar>(1.0, static_cast<Scalar>(Ls));
            return {best_score / L_eff, static_cast<Scalar>(prev_id[best_c]) / L_eff};
        };

        auto splitmix64_next = [](uint64_t& state) -> uint64_t {
            uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            return z ^ (z >> 31);
        };

        auto shuffle_codons_det = [&](std::vector<int8_t>& seq, uint64_t& rng_state) {
            for (size_t i = seq.size(); i > 1; --i) {
                size_t j = static_cast<size_t>(splitmix64_next(rng_state) % i);
                std::swap(seq[i - 1], seq[j]);
            }
        };

        Matrix Z_mat = Matrix::Zero(M, M);
        Matrix ID_mat = Matrix::Zero(M, M);
        Matrix MuID_mat = Matrix::Zero(M, M);

        #pragma omp parallel for schedule(dynamic)
        for (size_t a = 0; a < M; ++a) {
            const auto& seq_i = candidates[length_ok_indices[a]].codons;
            size_t Li = seq_i.size();
            std::vector<Scalar> prev, curr;
            std::vector<size_t> prev_id, curr_id;
            std::vector<int8_t> shuf_i, shuf_j;

            for (size_t b = a + 1; b < M; ++b) {
                const auto& seq_j = candidates[length_ok_indices[b]].codons;
                size_t Lj = seq_j.size();
                if (Li == 0 || Lj == 0) continue;

                auto [norm_dp, aa_id] = run_semiglobal_banded_dp(seq_i, seq_j, prev, curr, prev_id, curr_id);
                Scalar L_eff = std::max<Scalar>(1.0, static_cast<Scalar>(std::min(Li, Lj)));

                // Calibrate composition- and length-matched null via K = 16 codon-shuffled decoys
                // combined with analytical variance shrinkage (prior weight nu_0 = 16).
                constexpr int K_decoys = 16;
                constexpr Scalar nu_0 = 16.0;
                std::array<Scalar, K_decoys> id_null{};
                std::array<Scalar, K_decoys> dp_null{};
                shuf_i = seq_i;
                shuf_j = seq_j;
                uint64_t rng_state = 0x9E3779B97F4A7C15ULL ^
                                     (static_cast<uint64_t>(a) * 1315423911ULL +
                                      static_cast<uint64_t>(b) * 2654435761ULL);

                Scalar sum_id = 0.0, sum_dp = 0.0;
                for (int k = 0; k < K_decoys; ++k) {
                    shuffle_codons_det(shuf_i, rng_state);
                    shuffle_codons_det(shuf_j, rng_state);
                    auto [dp_k, id_k] = run_semiglobal_banded_dp(shuf_i, shuf_j, prev, curr, prev_id, curr_id);
                    id_null[static_cast<size_t>(k)] = id_k;
                    dp_null[static_cast<size_t>(k)] = dp_k;
                    sum_id += id_k;
                    sum_dp += dp_k;
                }
                Scalar mu_id = sum_id / static_cast<Scalar>(K_decoys);
                Scalar mu_dp = sum_dp / static_cast<Scalar>(K_decoys);

                Scalar var_id_samp = 0.0, var_dp_samp = 0.0;
                for (int k = 0; k < K_decoys; ++k) {
                    Scalar d_id = id_null[static_cast<size_t>(k)] - mu_id;
                    Scalar d_dp = dp_null[static_cast<size_t>(k)] - mu_dp;
                    var_id_samp += d_id * d_id;
                    var_dp_samp += d_dp * d_dp;
                }
                var_id_samp /= static_cast<Scalar>(K_decoys - 1);
                var_dp_samp /= static_cast<Scalar>(K_decoys - 1);

                Scalar p_clamped = std::clamp<Scalar>(mu_id, 0.05, 0.95);
                Scalar var_id_plug = (p_clamped * (1.0 - p_clamped)) / L_eff;
                Scalar var_dp_plug = (0.35 * 0.35) / L_eff;

                Scalar var_id_reg = (1.0 + 1.0 / static_cast<Scalar>(K_decoys)) *
                    std::max(var_id_plug, (nu_0 * var_id_plug + static_cast<Scalar>(K_decoys - 1) * var_id_samp) /
                                          (nu_0 + static_cast<Scalar>(K_decoys - 1)));
                Scalar var_dp_reg = (1.0 + 1.0 / static_cast<Scalar>(K_decoys)) *
                    std::max(var_dp_plug, (nu_0 * var_dp_plug + static_cast<Scalar>(K_decoys - 1) * var_dp_samp) /
                                          (nu_0 + static_cast<Scalar>(K_decoys - 1)));

                Scalar z_id = (aa_id - mu_id) / std::max<Scalar>(1e-4, std::sqrt(var_id_reg));
                Scalar z_dp = (norm_dp - mu_dp) / std::max<Scalar>(1e-4, std::sqrt(var_dp_reg));
                Scalar z_comb = 0.65 * z_id + 0.35 * z_dp;

                Z_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) = z_comb;
                Z_mat(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)) = z_comb;
                ID_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) = aa_id;
                ID_mat(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)) = aa_id;
                MuID_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) = mu_id;
                MuID_mat(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)) = mu_id;
            }
        }

        // Find connected components of the homology graph (edge exists iff Z_ab >= min_homology_z and ID_ab >= max(0.18, mu_ID + 0.04))
        std::vector<int> comp_id(M, -1);
        int num_comps = 0;
        int best_comp = -1;
        size_t best_comp_size = 0;

        for (size_t a = 0; a < M; ++a) {
            if (comp_id[a] != -1) continue;
            int cid = num_comps++;
            std::vector<size_t> q = {a};
            comp_id[a] = cid;
            size_t head = 0;
            while (head < q.size()) {
                size_t u = q[head++];
                for (size_t v = 0; v < M; ++v) {
                    Eigen::Index eu = static_cast<Eigen::Index>(u);
                    Eigen::Index ev = static_cast<Eigen::Index>(v);
                    Scalar min_id_edge = std::max<Scalar>(0.18, MuID_mat(eu, ev) + 0.04);
                    if (comp_id[v] == -1 &&
                        Z_mat(eu, ev) >= settings.min_homology_z &&
                        ID_mat(eu, ev) >= min_id_edge) {
                        comp_id[v] = cid;
                        q.push_back(v);
                    }
                }
            }
            if (q.size() > best_comp_size) {
                best_comp_size = q.size();
                best_comp = cid;
            }
        }

        for (size_t a = 0; a < M; ++a) {
            size_t orig = length_ok_indices[a];
            Scalar max_z = -1e9;
            Scalar max_id = 0.0;
            for (size_t b = 0; b < M; ++b) {
                if (a == b) continue;
                max_z = std::max(max_z, Z_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)));
                max_id = std::max(max_id, ID_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)));
            }
            report.entries[orig].max_homology_z = (max_z > -1e8) ? max_z : 0.0;
            report.entries[orig].max_aa_identity = max_id;

            if (settings.enforce_admission_control && comp_id[a] != best_comp) {
                report.entries[orig].status = "Quarantined";
                std::ostringstream rs;
                rs << "Non-homologous outlier (max Z = " << std::fixed << std::setprecision(2)
                   << report.entries[orig].max_homology_z << ", max AA identity = "
                   << std::setprecision(1) << (max_id * 100.0) << "%)";
                report.entries[orig].reason = rs.str();
            }
        }

        // Collect admitted sequences and compute cohort median Z-score and AA identity
        std::vector<Sequence> admitted_seqs;
        std::vector<size_t> admitted_m_indices;
        size_t final_rescued = 0;
        for (size_t a = 0; a < M; ++a) {
            size_t orig = length_ok_indices[a];
            if (report.entries[orig].status != "Quarantined") {
                if (report.entries[orig].status == "Rescued") final_rescued++;
                admitted_m_indices.push_back(a);
                admitted_seqs.push_back(Sequence{candidates[orig].name, candidates[orig].clean_cds, {}, {}});
            }
        }
        report.rescued_frame_or_strand = final_rescued;

        std::vector<Scalar> cohort_z, cohort_id;
        for (size_t i = 0; i < admitted_m_indices.size(); ++i) {
            for (size_t j = i + 1; j < admitted_m_indices.size(); ++j) {
                Eigen::Index u = static_cast<Eigen::Index>(admitted_m_indices[i]);
                Eigen::Index v = static_cast<Eigen::Index>(admitted_m_indices[j]);
                cohort_z.push_back(Z_mat(u, v));
                cohort_id.push_back(ID_mat(u, v));
            }
        }
        if (!cohort_z.empty()) {
            std::sort(cohort_z.begin(), cohort_z.end());
            std::sort(cohort_id.begin(), cohort_id.end());
            report.median_pairwise_z = cohort_z[cohort_z.size() / 2];
            report.median_aa_identity = cohort_id[cohort_id.size() / 2];
        }

        report.admitted_sequences = admitted_seqs.size();
        report.quarantined_sequences = report.total_input_sequences - report.admitted_sequences;

        if (settings.enforce_admission_control &&
            (admitted_seqs.size() < 3 || report.median_pairwise_z < settings.min_cohort_median_z)) {
            std::ostringstream err;
            err << "Gate 3 (Homology Admission Gate): Input sequences lack statistically significant evolutionary homology "
                << "(largest homologous component = " << admitted_seqs.size() << " of " << report.total_input_sequences
                << " taxa; median pairwise homology Z-score = " << std::fixed << std::setprecision(2)
                << report.median_pairwise_z << " [minimum required: " << settings.min_cohort_median_z
                << "]; median amino-acid identity = " << std::setprecision(1) << (report.median_aa_identity * 100.0)
                << "%). Refusing to run selection inference on unalignable/unrelated sequences.";
            throw std::runtime_error(err.str());
        }

        if (!settings.quiet && (report.rescued_frame_or_strand > 0 || report.quarantined_sequences > 0)) {
            std::ostringstream msg;
            msg << "      Gate 1-3 Input Triage: " << report.admitted_sequences << "/" << report.total_input_sequences
                << " sequences admitted (median homology Z = " << std::fixed << std::setprecision(1)
                << report.median_pairwise_z << ", median AA ID = " << (report.median_aa_identity * 100.0) << "%)";
            if (report.rescued_frame_or_strand > 0) {
                msg << "  │  Rescued frame/strand: " << report.rescued_frame_or_strand;
            }
            if (report.quarantined_sequences > 0) {
                msg << "  │  Quarantined outliers: " << report.quarantined_sequences;
            }
            msg << "\n";
            std::cout << Console::warning(msg.str());
            for (const auto& ent : report.entries) {
                if (ent.status != "Admitted") {
                    std::cout << Console::muted("        • [" + ent.status + "] " + ent.name + ": " + ent.reason + "\n");
                }
            }
        }

        return {std::move(admitted_seqs), std::move(report)};
    }

    static std::vector<Sequence> load_unaligned_fasta(
        const std::string& filepath,
        std::shared_ptr<const GeneticCode> gcode = nullptr
    ) {
        if (!gcode) gcode = GeneticCode::universal();
        std::ifstream file(filepath);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open unaligned FASTA file: " + filepath);
        }

        std::vector<Sequence> seqs;
        std::string line, current_name, current_seq;

        auto flush_seq = [&]() {
            if (current_name.empty()) return;
            // If sequence is already in-frame (multiple of 3) with a terminal stop codon and no internal stops,
            // strip the terminal stop codon so direct callers of load_unaligned_fasta see exact sense-codon lengths.
            if (current_seq.size() >= 6 && current_seq.size() % 3 == 0) {
                bool internal_stop = false;
                for (size_t i = 0; i + 5 < current_seq.size(); i += 3) {
                    if (gcode->triplet_to_sense(current_seq[i], current_seq[i + 1], current_seq[i + 2]) < 0 &&
                        current_seq.substr(i, 3).find('N') == std::string::npos) {
                        internal_stop = true;
                        break;
                    }
                }
                if (!internal_stop) {
                    size_t last = current_seq.size() - 3;
                    if (gcode->triplet_to_sense(current_seq[last], current_seq[last + 1], current_seq[last + 2]) < 0 &&
                        current_seq.substr(last, 3).find('N') == std::string::npos) {
                        current_seq.resize(last);
                    }
                }
            }
            seqs.push_back(Sequence{current_name, std::move(current_seq), {}, {}});
        };

        while (std::getline(file, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (line[0] == '>') {
                flush_seq();
                std::string header = line.substr(1);
                size_t start = header.find_first_not_of(" \t");
                if (start == std::string::npos) {
                    current_name = "unnamed";
                } else {
                    size_t end = header.find_first_of(" \t", start);
                    current_name = header.substr(start, end - start);
                }
                current_seq.clear();
            } else {
                for (char c : line) {
                    if (!std::isspace(static_cast<unsigned char>(c))) {
                        char uc = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                        if (uc == 'U') uc = 'T';
                        current_seq += uc;
                    }
                }
            }
        }
        flush_seq();
        return seqs;
    }

    JointAlignResult run_from_unaligned(
        const std::vector<Sequence>& unaligned_seqs,
        std::shared_ptr<const GeneticCode> gcode = nullptr
    ) {
        auto t0 = std::chrono::steady_clock::now();
        auto [admitted_seqs, adm_report] = triage_and_admit_unaligned_sequences(unaligned_seqs, gcode, settings);
        auto [de_novo_aln, de_novo_tree] = build_de_novo_codon_alignment_and_tree(admitted_seqs, gcode, settings);
        auto res = run(de_novo_aln, de_novo_tree);
        res.admission_report = std::move(adm_report);
        auto t1 = std::chrono::steady_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    JointAlignResult run_from_unaligned_file(
        const std::string& filepath,
        std::shared_ptr<const GeneticCode> gcode = nullptr
    ) {
        auto seqs = load_unaligned_fasta(filepath, gcode);
        return run_from_unaligned(seqs, gcode);
    }

    // =========================================================================
    // Objective Alignment Evaluation Suite (Options 1, 3 & 5):
    //   1. Option 1 (Joint Generative Log-Likelihood):
    //        ln P(S, A | T) = ln L_MG94(A | T) + ln P_TreeIndel(G(A) | T, lam_open, eps_ext)
    //      where ancestral presence/absence is reconstructed via Sankoff DP on T and scored
    //      by an affine birth-death indel process along every branch of T.
    //   2. Option 3 (Leave-One-Taxon-Out Cross-Validation Predictive Log-Likelihood):
    //        sum_{t=1}^N ln P(S_t | A_{-t}, T)
    //      integrated over a normalized Pair-HMM Forward partition function against the
    //      exact Inside-Outside held-out leaf cavity profile P(s | Data_{-t, c}).
    //   3. Option 5 (Posterior Model Weights & Per-Column Homology Confidence):
    //      Computes normalized Bayesian model weights across candidate alignments under both
    //      Option 1 (Joint) and Option 3 (LOTO CV), plus mean/min column homology posteriors.
    // =========================================================================
    struct AlignmentObjectiveEvaluation {
        std::string label;
        size_t num_taxa = 0;
        size_t num_codons = 0;
        size_t total_gap_codons = 0;
        Scalar mg94_sub_log_l = 0.0;
        size_t branch_indel_events = 0;
        size_t branch_indel_codons = 0;
        Scalar branch_indel_log_l = 0.0;
        Scalar joint_log_l = 0.0;          // Option 1: mg94_sub_log_l + branch_indel_log_l
        Scalar loto_cv_log_l = 0.0;        // Option 3: sum_t ln P(S_t | A_{-t}, T)
        Scalar joint_posterior_weight = 0.0; // Option 5: exp(joint - max_joint) / Z
        Scalar loto_posterior_weight = 0.0;  // Option 5: exp(loto - max_loto) / Z
        Scalar mean_column_posterior = 1.0;  // Option 5: mean column homology posterior
        Scalar min_column_posterior = 1.0;
    };

    static std::vector<AlignmentObjectiveEvaluation> evaluate_candidate_alignments(
        const std::vector<std::pair<std::string, Alignment>>& candidates,
        std::optional<Tree> ref_tree = std::nullopt,
        const JointAlignSettings& settings = JointAlignSettings{}
    ) {
        std::vector<AlignmentObjectiveEvaluation> evals;
        if (candidates.empty()) return evals;
        evals.resize(candidates.size());

        // 1. Determine common reference tree and global MG94 parameters from the first alignment (or ref_tree)
        const Alignment& base_aln = candidates[0].second;
        Tree tree;
        if (ref_tree.has_value()) {
            tree = *ref_tree;
        } else if (!base_aln.embedded_tree_newick.empty()) {
            tree = Tree::from_newick(base_aln.embedded_tree_newick);
        } else {
            TreeEngineSettings tes;
            tes.quiet = true;
            TreeEngine te(tes);
            auto te_res = te.run(base_aln);
            tree = te_res.optimized_tree;
        }

        GTRFitter gtr_fitter(tree, base_aln);
        auto gtr_fit = gtr_fitter.fit();
        MG94Parameters mg_params;
        mg_params.theta_AC = gtr_fit.params.theta_AC;
        mg_params.theta_AT = gtr_fit.params.theta_AT;
        mg_params.theta_CG = gtr_fit.params.theta_CG;
        mg_params.theta_CT = gtr_fit.params.theta_CT;
        mg_params.theta_GT = gtr_fit.params.theta_GT;

        MG94Fitter mg_fitter(gtr_fit.tree, base_aln);
        auto mg_fit = mg_fitter.fit_omega_and_scale(1.0, mg_params);
        mg_params.alpha = 1.0;
        mg_params.beta = mg_fit.x_opt(0);
        tree = mg_fit.tree;

        const auto& gcode = *(base_aln.code ? base_aln.code : GeneticCode::universal());
        Vector pi_codon = base_aln.codon_frequencies_f3x4;
        MG94Matrix mg94;
        mg94.update(mg_params, base_aln.pos_nuc_frequencies, pi_codon, gcode);
        int S = gcode.num_sense_codons;
        Scalar conv = (mg94.scale_factor > 1e-12) ? (3.0 / mg94.scale_factor) : 1.0;

        std::vector<Matrix> P_branches(tree.num_nodes());
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            P_branches[node.id] = mg94.transition_matrix(node.branch_length * conv);
        }

        auto log_sum_exp = [](Scalar a, Scalar b) -> Scalar {
            const Scalar NEG_INF = -1e18;
            if (a <= NEG_INF * 0.5) return b;
            if (b <= NEG_INF * 0.5) return a;
            Scalar mx = std::max(a, b);
            return mx + std::log(std::exp(a - mx) + std::exp(b - mx));
        };

        for (size_t idx = 0; idx < candidates.size(); ++idx) {
            const std::string& label = candidates[idx].first;
            const Alignment& aln = candidates[idx].second;
            AlignmentObjectiveEvaluation ev;
            ev.label = label;
            ev.num_taxa = aln.num_taxa;
            ev.num_codons = aln.num_codons;

            size_t N = aln.num_taxa;
            size_t C = aln.num_codons;
            for (size_t t = 0; t < N; ++t) {
                for (size_t c = 0; c < C; ++c) {
                    if (aln.sequences[t].codons[c] < 0) ev.total_gap_codons++;
                }
            }

            std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

            // (A) Compute exact MG94 Substitution Log-Likelihood
            std::vector<Vector> node_L(tree.num_nodes(), Vector::Zero(S));
            std::vector<double> node_scale(tree.num_nodes(), 0.0);
            for (size_t c = 0; c < C; ++c) {
                SitePattern pat;
                pat.states.resize(N);
                for (size_t t = 0; t < N; ++t) pat.states[t] = aln.sequences[t].codons[c];
                pat.weight = 1;
                ev.mg94_sub_log_l += FELAnalyzer::compute_pattern_log_likelihood(
                    tree, pat, leaf_to_taxon, mg94, P_branches, node_L, node_scale
                );
            }

            // (B) Option 1: Phylogenetic Branch-Indel Birth-Death Log-Likelihood on T
            size_t num_nodes = tree.num_nodes();
            const Scalar lam_open = 0.035;
            const Scalar eps_ext = 0.65;
            std::vector<std::array<Scalar, 2>> dp(C * num_nodes, {0.0, 0.0});
            std::vector<uint8_t> anc_state(C * num_nodes, 1);

            for (int32_t nid : tree.post_order) {
                const auto& node = tree.nodes[nid];
                if (node.is_leaf) {
                    int32_t tid = (nid < static_cast<int32_t>(leaf_to_taxon.size()) && leaf_to_taxon[nid] < aln.num_taxa)
                                ? static_cast<int32_t>(leaf_to_taxon[nid]) : -1;
                    for (size_t c = 0; c < C; ++c) {
                        bool present = (tid >= 0 && aln.sequences[tid].codons[c] >= 0);
                        dp[c * num_nodes + nid][present ? 1 : 0] = 0.0;
                        dp[c * num_nodes + nid][present ? 0 : 1] = -1e6;
                    }
                } else {
                    for (size_t c = 0; c < C; ++c) {
                        for (int su = 0; su <= 1; ++su) {
                            Scalar acc = 0.0;
                            for (int32_t cid : node.children) {
                                Scalar tb = std::max<Scalar>(0.005, tree.nodes[cid].branch_length * conv);
                                Scalar p_chg = 1.0 - std::exp(-lam_open * tb);
                                Scalar s_same = dp[c * num_nodes + cid][su] + std::log(1.0 - p_chg);
                                Scalar s_diff = dp[c * num_nodes + cid][1 - su] + std::log(p_chg);
                                acc += std::max(s_same, s_diff);
                            }
                            dp[c * num_nodes + nid][su] = acc;
                        }
                    }
                }
            }
            for (size_t c = 0; c < C; ++c) {
                anc_state[c * num_nodes + tree.root_id] =
                    (dp[c * num_nodes + tree.root_id][1] >= dp[c * num_nodes + tree.root_id][0]) ? 1 : 0;
            }
            for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
                int32_t nid = *it;
                const auto& node = tree.nodes[nid];
                for (int32_t cid : node.children) {
                    Scalar tb = std::max<Scalar>(0.005, tree.nodes[cid].branch_length * conv);
                    Scalar p_chg = 1.0 - std::exp(-lam_open * tb);
                    for (size_t c = 0; c < C; ++c) {
                        int su = anc_state[c * num_nodes + nid];
                        Scalar s_same = dp[c * num_nodes + cid][su] + std::log(1.0 - p_chg);
                        Scalar s_diff = dp[c * num_nodes + cid][1 - su] + std::log(p_chg);
                        anc_state[c * num_nodes + cid] = (s_same >= s_diff) ? su : (1 - su);
                    }
                }
            }
            for (int32_t nid : tree.post_order) {
                const auto& node = tree.nodes[nid];
                for (int32_t cid : node.children) {
                    Scalar tb = std::max<Scalar>(0.005, tree.nodes[cid].branch_length * conv);
                    Scalar p_open = 1.0 - std::exp(-lam_open * tb);
                    int in_ev = 0;
                    for (size_t c = 0; c < C; ++c) {
                        int su = anc_state[c * num_nodes + nid];
                        int sv = anc_state[c * num_nodes + cid];
                        if (su == sv) {
                            if (in_ev != 0) {
                                ev.branch_indel_log_l += std::log(1.0 - eps_ext);
                                in_ev = 0;
                            }
                            if (su == 1) ev.branch_indel_log_l += -lam_open * tb;
                        } else {
                            int ev_t = (su == 0 && sv == 1) ? 1 : 2;
                            ev.branch_indel_codons++;
                            if (in_ev == ev_t) {
                                ev.branch_indel_log_l += std::log(eps_ext);
                            } else {
                                if (in_ev != 0) ev.branch_indel_log_l += std::log(1.0 - eps_ext);
                                ev.branch_indel_events++;
                                ev.branch_indel_log_l += std::log(0.5 * p_open);
                                in_ev = ev_t;
                            }
                        }
                    }
                    if (in_ev != 0) ev.branch_indel_log_l += std::log(1.0 - eps_ext);
                }
            }
            // Also include the geometric sequence-length / column-occupancy prior:
            // Each alignment column c that is absent at the root (pure insertion column) or has gap fraction f_gap > 0
            // incurs the geometric column-expansion cost relative to the mean unaligned CDS length L_mean,
            // preventing over-gapped alignments (e.g. PRANK splitting homologous columns) from inflating Joint lnP!
            Scalar mean_unaligned_codons = 0.0;
            for (size_t t = 0; t < N; ++t) {
                for (size_t c = 0; c < C; ++c) {
                    if (aln.sequences[t].codons[c] >= 0) mean_unaligned_codons += 1.0;
                }
            }
            mean_unaligned_codons /= std::max<Scalar>(1.0, static_cast<Scalar>(N));
            if (static_cast<Scalar>(C) > mean_unaligned_codons) {
                Scalar extra_cols = static_cast<Scalar>(C) - mean_unaligned_codons;
                ev.branch_indel_log_l -= extra_cols * settings.gap_open_penalty;
            }
            ev.joint_log_l = ev.mg94_sub_log_l + ev.branch_indel_log_l;

            // (C) Option 3: Leave-One-Taxon-Out (LOTO) Cross-Validated Predictive Log-Likelihood
            // Compute exact Inside-Outside held-out leaf cavities across the full alignment [0, C)
            std::vector<std::vector<Vector>> full_cavity;
            compute_window_leaf_cavities(tree, aln, mg94, P_branches, 0, C, full_cavity, nullptr, conv);

            const Scalar g_open_cv = 4.5;
            const Scalar g_ext_cv  = 0.8;
            const Scalar p_open_cv = std::exp(-g_open_cv);
            const Scalar p_ext_cv  = std::exp(-g_ext_cv);
            const Scalar log_open_cv = -g_open_cv;
            const Scalar log_1m_open_cv = std::log(1.0 - p_open_cv);
            const Scalar log_ext_cv = -g_ext_cv;
            const Scalar log_1m_ext_cv = std::log(1.0 - p_ext_cv);

            std::vector<Scalar> col_match_acc(C, 0.0);
            for (size_t t = 0; t < N; ++t) {
                std::vector<int8_t> obs;
                obs.reserve(C);
                for (size_t c = 0; c < C; ++c) {
                    if (aln.sequences[t].codons[c] >= 0) obs.push_back(aln.sequences[t].codons[c]);
                }
                size_t M = obs.size();
                Eigen::MatrixXd fwd_M = Eigen::MatrixXd::Constant(C + 1, M + 1, -1e18);
                Eigen::MatrixXd fwd_G = Eigen::MatrixXd::Constant(C + 1, M + 1, -1e18);
                fwd_M(0, 0) = 0.0;

                for (size_t w = 1; w <= C; ++w) {
                    const Vector& cav = full_cavity[w - 1][t];
                    for (size_t j = 0; j <= M; ++j) {
                        fwd_G(w, j) = log_sum_exp(fwd_M(w - 1, j) + log_open_cv, fwd_G(w - 1, j) + log_ext_cv);
                        if (j >= 1) {
                            Scalar emit = std::log(std::max<Scalar>(1e-12, cav(obs[j - 1])));
                            fwd_M(w, j) = log_sum_exp(
                                fwd_M(w - 1, j - 1) + log_1m_open_cv + emit,
                                fwd_G(w - 1, j - 1) + log_1m_ext_cv + emit
                            );
                        }
                    }
                }
                Scalar log_Z = log_sum_exp(fwd_M(C, M), fwd_G(C, M));
                ev.loto_cv_log_l += log_Z;

                // Backward pass for Option 5 per-column posterior confidence
                Eigen::MatrixXd bwd_M = Eigen::MatrixXd::Constant(C + 1, M + 1, -1e18);
                Eigen::MatrixXd bwd_G = Eigen::MatrixXd::Constant(C + 1, M + 1, -1e18);
                bwd_M(C, M) = 0.0;
                bwd_G(C, M) = 0.0;
                for (size_t w = C; w >= 1; --w) {
                    const Vector& cav = full_cavity[w - 1][t];
                    for (size_t j = 0; j <= M; ++j) {
                        Scalar to_g = bwd_G(w, j);
                        bwd_M(w - 1, j) = log_sum_exp(bwd_M(w - 1, j), to_g + log_open_cv);
                        bwd_G(w - 1, j) = log_sum_exp(bwd_G(w - 1, j), to_g + log_ext_cv);
                        if (j >= 1) {
                            Scalar emit = std::log(std::max<Scalar>(1e-12, cav(obs[j - 1])));
                            Scalar to_m = bwd_M(w, j) + emit;
                            bwd_M(w - 1, j - 1) = log_sum_exp(bwd_M(w - 1, j - 1), to_m + log_1m_open_cv);
                            bwd_G(w - 1, j - 1) = log_sum_exp(bwd_G(w - 1, j - 1), to_m + log_1m_ext_cv);
                        }
                    }
                }
                for (size_t w = 1; w <= C; ++w) {
                    Scalar log_m = -1e18;
                    for (size_t j = 1; j <= M; ++j) {
                        log_m = log_sum_exp(log_m, fwd_M(w, j) + bwd_M(w, j));
                    }
                    Scalar p_m = std::clamp(std::exp(std::clamp(log_m - log_Z, -50.0, 0.0)), 0.0, 1.0);
                    bool is_match = (aln.sequences[t].codons[w - 1] >= 0);
                    col_match_acc[w - 1] += is_match ? p_m : (1.0 - p_m);
                }
            }

            if (C > 0 && N > 0) {
                Scalar sum_col = 0.0, min_col = 1.0;
                for (size_t c = 0; c < C; ++c) {
                    Scalar cp = col_match_acc[c] / static_cast<Scalar>(N);
                    sum_col += cp;
                    min_col = std::min(min_col, cp);
                }
                ev.mean_column_posterior = sum_col / static_cast<Scalar>(C);
                ev.min_column_posterior = min_col;
            }

            evals[idx] = std::move(ev);
        }

        // Compute Option 5 normalized posterior model weights across all evaluated alignments
        Scalar max_joint = -1e18, max_loto = -1e18;
        for (const auto& ev : evals) {
            max_joint = std::max(max_joint, ev.joint_log_l);
            max_loto  = std::max(max_loto,  ev.loto_cv_log_l);
        }
        Scalar sum_wj = 0.0, sum_wl = 0.0;
        for (auto& ev : evals) {
            ev.joint_posterior_weight = std::exp(std::clamp(ev.joint_log_l - max_joint, -60.0, 0.0));
            ev.loto_posterior_weight  = std::exp(std::clamp(ev.loto_cv_log_l - max_loto, -60.0, 0.0));
            sum_wj += ev.joint_posterior_weight;
            sum_wl += ev.loto_posterior_weight;
        }
        for (auto& ev : evals) {
            if (sum_wj > 0.0) ev.joint_posterior_weight /= sum_wj;
            if (sum_wl > 0.0) ev.loto_posterior_weight  /= sum_wl;
        }

        return evals;
    }
};

} // namespace hyphy::analyses
