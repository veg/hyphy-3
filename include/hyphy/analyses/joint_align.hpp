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
    Matrix admitted_dp_similarity;
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
    Scalar transducer_gap_open = 2.85;   // Phase 2 progressive transducer affine gap-open penalty (in nats)
    Scalar transducer_gap_extend = 0.65; // Phase 2 progressive transducer affine gap-extend penalty (in nats)
    Scalar transducer_laplace_floor = 0.22; // Phase 2 homology Bayes-factor Laplace floor \epsilon_0
    Scalar gap_open_penalty = 3.2;       // Phase 3 window realignment log-penalty for opening a codon deletion block
    Scalar gap_extend_penalty = 0.85;    // Phase 3 window realignment log-penalty per additional deleted codon
    size_t window_flank_codons = 2;      // conserved anchor codons on each side of a gappy region
    size_t max_window_codons = 18;       // maximum width of a local realignment window
    size_t em_iterations = 2;            // Inside-Outside cavity realignment EM sweeps

    // Subquadratic scaling settings for large alignments (N > 64 taxa or L > 200 codons)
    size_t subquadratic_taxon_threshold = 64;   // switch to O(N sqrt(N)) landmark + k-mer k-NN graph when N > threshold
    size_t subquadratic_length_threshold = 200; // switch to O(B * L) seed-chain anchor-banded profile DP when min(L1, L2) > threshold

    // Stage 2 settings
    bool enable_coupled_topology_search = true;
    size_t max_coupled_spr_candidates = 8;

    // Stage 3 settings
    bool run_marginalized_fel = true;
    Scalar p_value_threshold = 0.10;
    bool quiet = false;
};

struct TransducerTaxonResult {
    std::vector<int8_t> map_codons;          // length W
    std::vector<uint64_t> map_masks;         // length W: resolved IUPAC bitmasks
    std::string map_raw_subseq;              // length 3*W
    Scalar input_log_score = 0.0;
    Scalar map_log_score = 0.0;
    Scalar log_partition_Z = 0.0;
    Scalar map_posterior_prob = 1.0;
    std::vector<Scalar> col_match_posterior; // length W: posterior prob that column w is a homologous Match
    bool changed = false;
    size_t codons_shifted = 0;

    // Top K candidate placements with their exact Forward-Backward posterior weights
    struct CandidatePlacement {
        std::vector<int8_t> codons;
        Scalar log_score = -1e20;
        Scalar posterior_weight = 0.0;
    };
    std::vector<CandidatePlacement> top_candidates;
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
    std::vector<std::vector<TransducerTaxonResult>> window_taxon_results;
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
    // Empirical Two-Stratum (Core + Loop) Physicochemical & Multi-Hit Codon Model
    //   Shared across Stage 1 Phase 2 (progressive profile alignment) and
    //   Phase 3 (iterative Inside-Outside window refinement).
    // =========================================================================
    static Scalar aa_phys_dist(char a1, char a2) {
        if (a1 == a2) return 0.0;
        struct AAProps { char aa; Scalar pol, vol, chg, hyd; };
        static constexpr AAProps table[] = {
            {'A',  8.1,  31.0,  0.0,  1.8}, {'R', 10.5, 124.0,  1.0, -4.5},
            {'N', 11.6,  56.0,  0.0, -3.5}, {'D', 13.0,  54.0, -1.0, -3.5},
            {'C',  5.5,  55.0,  0.0,  2.5}, {'Q', 10.5,  85.0,  0.0, -3.5},
            {'E', 12.3,  83.0, -1.0, -3.5}, {'G',  9.0,   3.0,  0.0, -0.4},
            {'H', 10.4,  96.0,  0.5, -3.2}, {'I',  5.2, 111.0,  0.0,  4.5},
            {'L',  4.9, 111.0,  0.0,  3.8}, {'K', 11.3, 119.0,  1.0, -3.9},
            {'M',  5.7, 105.0,  0.0,  1.9}, {'F',  5.2, 132.0,  0.0,  2.8},
            {'P',  8.0,  32.5,  0.0, -1.6}, {'S',  9.2,  32.0,  0.0, -0.8},
            {'T',  8.6,  61.0,  0.0, -0.7}, {'W',  5.4, 170.0,  0.0, -0.9},
            {'Y',  6.2, 136.0,  0.0, -1.3}, {'V',  5.9,  84.0,  0.0,  4.2}
        };
        const AAProps *p1 = nullptr, *p2 = nullptr;
        for (const auto& e : table) { if (e.aa == a1) p1 = &e; if (e.aa == a2) p2 = &e; }
        if (!p1 || !p2) return 0.5;
        Scalar dp = (p1->pol - p2->pol) / 8.1;
        Scalar dv = (p1->vol - p2->vol) / 167.0;
        Scalar dc = (p1->chg - p2->chg) / 2.0;
        Scalar dh = (p1->hyd - p2->hyd) / 9.0;
        return std::clamp<Scalar>(std::sqrt((dp * dp + dv * dv + dc * dc + dh * dh) / 4.0), 0.0, 1.0);
    }

    struct EmpRegimeParams {
        Scalar omega_hat = 0.25;
        Scalar theta_tv  = 0.35;
        Scalar theta_ct  = 1.00;
        Scalar lambda_aa = 1.80;
        Scalar delta_2nt = 0.15;
        Scalar psi_3nt   = 0.04;
        Scalar obs_mean_daa = 0.25;
    };

    struct EmpRegimeCounts {
        Scalar sum_ag = 1.0, sum_ct = 1.0, sum_tv = 1.0;
        Scalar sum_syn1 = 1.0, sum_ns1 = 1.0, sum_ns1_daa = 0.25;
        Scalar sum_1nt = 1.0, sum_2nt = 0.5, sum_3nt = 0.1;

        void init_from_neutral(Scalar mu0_daa) {
            sum_ns1_daa = mu0_daa;
        }

        void add_codon_pair(int8_t c1, int8_t c2, const GeneticCode& code) {
            if (c1 < 0 || c2 < 0 || c1 == c2) return;
            const auto& diff = code.diff_matrix[c1][c2];
            Scalar d_aa = diff.synonymous ? 0.0 : aa_phys_dist(code.sense_codon_aa[c1], code.sense_codon_aa[c2]);
            if (diff.num_diffs == 1) {
                sum_1nt += 1.0;
                int b1 = diff.diff_from[0], b2 = diff.diff_to[0];
                bool p1 = (b1 == 0 || b1 == 2), p2 = (b2 == 0 || b2 == 2);
                if (p1 && p2) sum_ag += 1.0;
                else if (!p1 && !p2) sum_ct += 1.0;
                else sum_tv += 1.0;
                if (diff.synonymous) {
                    sum_syn1 += 1.0;
                } else {
                    sum_ns1 += 1.0;
                    sum_ns1_daa += d_aa;
                }
            } else if (diff.num_diffs == 2) {
                sum_2nt += 1.0;
            } else if (diff.num_diffs == 3) {
                sum_3nt += 1.0;
            }
        }

        EmpRegimeParams finalize(
            Scalar mu0_daa, Scalar var0_daa,
            Scalar syn_to_ns_opp_ratio,
            Scalar opp_1_to_2, Scalar opp_1_to_3
        ) const {
            EmpRegimeParams p;
            p.theta_ct  = std::clamp<Scalar>(sum_ct / sum_ag, 0.40, 2.50);
            p.theta_tv  = std::clamp<Scalar>(sum_tv / (sum_ag + sum_ct), 0.10, 0.85);
            Scalar raw_omega = (sum_ns1 / sum_syn1) * syn_to_ns_opp_ratio;
            p.omega_hat = std::clamp<Scalar>(raw_omega, 0.06, 1.50);
            p.obs_mean_daa = sum_ns1_daa / std::max<Scalar>(1.0, sum_ns1);
            Scalar raw_lambda = (mu0_daa - p.obs_mean_daa) / std::max<Scalar>(1e-4, var0_daa);
            p.lambda_aa = std::clamp<Scalar>(raw_lambda, 0.0, 3.50);
            p.delta_2nt = std::clamp<Scalar>((sum_2nt / sum_1nt) * opp_1_to_2, 0.06, 0.55);
            p.psi_3nt   = std::clamp<Scalar>((sum_3nt / sum_1nt) * opp_1_to_3, 0.01, 0.25);
            return p;
        }
    };

    struct EmpiricalSpectralGenerator {
        Matrix U_pi, V_pi;
        Vector eigenvalues;
        Scalar scale_factor = 1.0;

        void build(const EmpRegimeParams& st, const Eigen::Matrix<Scalar, 3, 4>& pf,
                   const Vector& pi, const GeneticCode& code) {
            int nS = code.num_sense_codons;
            Vector sqrt_pi(nS), inv_sqrt_pi(nS);
            for (int i = 0; i < nS; ++i) {
                sqrt_pi(i) = std::sqrt(std::max<Scalar>(pi(i), 1e-12));
                inv_sqrt_pi(i) = 1.0 / sqrt_pi(i);
            }
            Scalar r_nuc[4][4] = {
                {0.0, st.theta_tv, 1.0, st.theta_tv},
                {st.theta_tv, 0.0, st.theta_tv, st.theta_ct},
                {1.0, st.theta_tv, 0.0, st.theta_tv},
                {st.theta_tv, st.theta_ct, st.theta_tv, 0.0}
            };
            Scalar sum_aa_w = 0.0;
            int count_ns1 = 0;
            for (int i = 0; i < nS; ++i) {
                for (int j = i + 1; j < nS; ++j) {
                    const auto& diff = code.diff_matrix[i][j];
                    if (diff.num_diffs == 1 && !diff.synonymous) {
                        sum_aa_w += std::exp(-st.lambda_aa * aa_phys_dist(code.sense_codon_aa[i], code.sense_codon_aa[j]));
                        count_ns1++;
                    }
                }
            }
            Scalar mean_aa_w = (count_ns1 > 0) ? (sum_aa_w / static_cast<Scalar>(count_ns1)) : 1.0;

            Matrix Q = Matrix::Zero(nS, nS);
            for (int i = 0; i < nS; ++i) {
                for (int j = 0; j < nS; ++j) {
                    if (i == j) continue;
                    const auto& diff = code.diff_matrix[i][j];
                    Scalar d_aa = diff.synonymous ? 0.0 : aa_phys_dist(code.sense_codon_aa[i], code.sense_codon_aa[j]);
                    Scalar aa_mod = diff.synonymous ? 1.0 : (st.omega_hat * std::exp(-st.lambda_aa * d_aa) / mean_aa_w);
                    if (diff.num_diffs == 1) {
                        Q(i, j) = aa_mod * r_nuc[diff.diff_from[0]][diff.diff_to[0]] * pf(diff.diff_pos[0], diff.diff_to[0]);
                    } else if (diff.num_diffs == 2) {
                        Scalar nr = r_nuc[diff.diff_from[0]][diff.diff_to[0]] * r_nuc[diff.diff_from[1]][diff.diff_to[1]];
                        Scalar tf = pf(diff.diff_pos[0], diff.diff_to[0]) * pf(diff.diff_pos[1], diff.diff_to[1]);
                        Scalar mh = diff.synonymous ? 1.0 : st.delta_2nt;
                        Q(i, j) = aa_mod * mh * nr * tf;
                    } else if (diff.num_diffs == 3) {
                        Scalar nr = r_nuc[diff.diff_from[0]][diff.diff_to[0]] * r_nuc[diff.diff_from[1]][diff.diff_to[1]] * r_nuc[diff.diff_from[2]][diff.diff_to[2]];
                        Scalar tf = pf(diff.diff_pos[0], diff.diff_to[0]) * pf(diff.diff_pos[1], diff.diff_to[1]) * pf(diff.diff_pos[2], diff.diff_to[2]);
                        Scalar mh = diff.synonymous ? 1.0 : st.psi_3nt;
                        Q(i, j) = aa_mod * mh * nr * tf;
                    }
                }
                Q(i, i) = -Q.row(i).sum();
            }
            Scalar avg_rate = 0.0;
            for (int i = 0; i < nS; ++i) avg_rate -= pi(i) * Q(i, i);
            scale_factor = std::max<Scalar>(1e-12, avg_rate);

            Matrix S_mat = sqrt_pi.asDiagonal() * Q * inv_sqrt_pi.asDiagonal();
            S_mat = 0.5 * (S_mat + S_mat.transpose());
            Eigen::SelfAdjointEigenSolver<Matrix> solver(S_mat);
            eigenvalues = solver.eigenvalues().cwiseMin(0.0);
            Matrix eigvecs = solver.eigenvectors();
            U_pi = inv_sqrt_pi.asDiagonal() * eigvecs;
            V_pi = eigvecs.transpose() * sqrt_pi.asDiagonal();
        }

        Matrix transition_matrix(Scalar branch_length) const {
            Scalar conv = 3.0 / scale_factor;
            Vector exp_l = (eigenvalues * (branch_length * conv)).array().exp();
            return ((U_pi * exp_l.asDiagonal()) * V_pi).cwiseMax(0.0);
        }
    };

    struct EmpiricalTwoStratumModel {
        EmpRegimeParams emp_core, emp_loop;
        Scalar w_loop = 0.30;
        Scalar rho_loop_to_core = 4.0;
        EmpiricalSpectralGenerator q_core, q_loop;
        Vector pi;

        Matrix transition_matrix(Scalar bl) const {
            Scalar t = std::max<Scalar>(1e-4, bl);
            return (1.0 - w_loop) * q_core.transition_matrix(t)
                 + w_loop * q_loop.transition_matrix(rho_loop_to_core * t);
        }

        static EmpiricalTwoStratumModel calibrate(
            const std::vector<std::vector<int8_t>>& taxon_codons,
            const Eigen::Matrix<Scalar, 3, 4>& pos_freqs,
            const Vector& pi_codon,
            const GeneticCode& gcode,
            size_t subquadratic_taxon_threshold = 64
        ) {
            EmpiricalTwoStratumModel m;
            m.pi = pi_codon;
            int S = gcode.num_sense_codons;
            size_t N = taxon_codons.size();

            Scalar sum_d0 = 0.0, sum_d0_sq = 0.0;
            Scalar n_ns1_opp = 0.0, n_syn1_opp = 0.0;
            Scalar n_1nt_opp = 0.0, n_2nt_opp = 0.0, n_3nt_opp = 0.0;
            for (int i = 0; i < S; ++i) {
                for (int j = i + 1; j < S; ++j) {
                    const auto& diff = gcode.diff_matrix[i][j];
                    if (diff.num_diffs == 1) {
                        n_1nt_opp += 1.0;
                        if (diff.synonymous) {
                            n_syn1_opp += 1.0;
                        } else {
                            Scalar d = aa_phys_dist(gcode.sense_codon_aa[i], gcode.sense_codon_aa[j]);
                            sum_d0 += d;
                            sum_d0_sq += d * d;
                            n_ns1_opp += 1.0;
                        }
                    } else if (diff.num_diffs == 2) {
                        n_2nt_opp += 1.0;
                    } else if (diff.num_diffs == 3) {
                        n_3nt_opp += 1.0;
                    }
                }
            }
            Scalar mu0_daa = (n_ns1_opp > 0) ? (sum_d0 / n_ns1_opp) : 0.32;
            Scalar var0_daa = (n_ns1_opp > 0) ? std::max<Scalar>(1e-4, sum_d0_sq / n_ns1_opp - mu0_daa * mu0_daa) : 0.035;
            Scalar syn_to_ns_opp_ratio = (n_ns1_opp > 0) ? (n_syn1_opp / n_ns1_opp) : (0.25 / 0.75);
            Scalar opp_1_to_2 = (n_2nt_opp > 0) ? (n_1nt_opp / n_2nt_opp) : 0.35;
            Scalar opp_1_to_3 = (n_3nt_opp > 0) ? (n_1nt_opp / n_3nt_opp) : 0.32;

            EmpRegimeCounts cnt_all, cnt_core, cnt_loop;
            cnt_all.init_from_neutral(mu0_daa);
            cnt_core.init_from_neutral(mu0_daa);
            cnt_loop.init_from_neutral(mu0_daa);
            Scalar sum_r_core = 0.0, sum_r_loop = 0.0;
            size_t n_core = 0, n_loop = 0;

            const bool use_subquadratic_cal = (N > subquadratic_taxon_threshold);
            const size_t cal_stride = use_subquadratic_cal ? std::max<size_t>(1, N / 6) : 1;

            for (size_t u = 0; u < N; ++u) {
                const auto& su = taxon_codons[u];
                for (size_t v = u + 1; v < N; ++v) {
                    if (use_subquadratic_cal && (v - u) > 2 && ((v - u) % cal_stride != 0)) continue;
                    const auto& sv = taxon_codons[v];
                    if (su.empty() || sv.empty()) continue;
                    size_t Lmin = std::min(su.size(), sv.size());
                    size_t Lmax = std::max(su.size(), sv.size());
                    if (Lmin < 5 || Lmax > Lmin * 2) continue;

                    // Align (su, sv) via banded affine codon DP so internal indels never shift
                    // post-indel codons out of register during Stage-0 rate matrix calibration.
                    const auto& s_short = (su.size() <= sv.size()) ? su : sv;
                    const auto& s_long  = (su.size() <= sv.size()) ? sv : su;
                    int Ls = static_cast<int>(s_short.size());
                    int Ll = static_cast<int>(s_long.size());
                    int dL = Ll - Ls;
                    int B0 = std::max(24, dL + 8);
                    int band_w = 2 * B0 + dL + 3;
                    constexpr Scalar kNegInf = -1e9;
                    constexpr Scalar cal_g_open = 1.10;
                    constexpr Scalar cal_g_ext  = 0.30;

                    std::vector<Scalar> prev_M(Ll + 1, kNegInf), curr_M(Ll + 1, kNegInf);
                    std::vector<Scalar> prev_D(Ll + 1, kNegInf), curr_D(Ll + 1, kNegInf);
                    std::vector<uint8_t> bt_mat((Ls + 1) * band_w, 0); // 0=M, 1=D(up), 2=I(left)
                    auto b_idx = [&](int r, int c) -> size_t {
                        int off = c - (r - B0);
                        return static_cast<size_t>(r * band_w + std::clamp(off, 0, band_w - 1));
                    };

                    int init_max = std::min(Ll, dL + B0);
                    for (int c = 0; c <= init_max; ++c) {
                        prev_M[c] = -std::min<Scalar>(cal_g_open + cal_g_ext * c, 3.80);
                        bt_mat[b_idx(0, c)] = 2;
                    }
                    prev_M[0] = 0.0;

                    for (int r = 1; r <= Ls; ++r) {
                        int c_min = std::max(1, r - B0);
                        int c_max = std::min(Ll, r + dL + B0);
                        int prev_c_max = std::min(Ll, (r - 1) + dL + B0);
                        std::fill(curr_M.begin() + std::max(0, c_min - 1), curr_M.begin() + c_max + 1, kNegInf);
                        std::fill(curr_D.begin() + std::max(0, c_min - 1), curr_D.begin() + c_max + 1, kNegInf);
                        if (r <= B0) {
                            Scalar pen = -std::min<Scalar>(cal_g_open + cal_g_ext * r, 3.80);
                            curr_D[0] = pen;
                            curr_M[0] = pen;
                            bt_mat[b_idx(r, 0)] = 1;
                        }
                        int8_t ci = s_short[r - 1];
                        Scalar ins_sc = kNegInf;
                        for (int c = c_min; c <= c_max; ++c) {
                            int8_t cj = s_long[c - 1];
                            Scalar d_sc = (c <= prev_c_max)
                                        ? std::max(prev_M[c] - cal_g_open - cal_g_ext, prev_D[c] - cal_g_ext)
                                        : kNegInf;
                            curr_D[c] = d_sc;
                            bool left_in_band = (c > c_min) || (c == 1 && r <= B0);
                            Scalar i_op = left_in_band ? (curr_M[c - 1] - cal_g_open - cal_g_ext) : kNegInf;
                            Scalar i_ex = (c > c_min) ? (ins_sc - cal_g_ext) : kNegInf;
                            ins_sc = std::max(i_op, i_ex);
                            Scalar sim = 0.0;
                            if (ci >= 0 && cj >= 0) {
                                if (ci == cj) sim = 1.00;
                                else if (gcode.sense_codon_aa[ci] == gcode.sense_codon_aa[cj]) sim = 0.90;
                                else if (gcode.diff_matrix[ci][cj].num_diffs == 1) sim = 0.30;
                            }
                            Scalar m_sc = prev_M[c - 1] + sim;
                            Scalar best_cell = m_sc;
                            uint8_t best_bt = 0;
                            if (d_sc > best_cell) {
                                best_cell = d_sc;
                                best_bt = 1;
                            }
                            if (ins_sc > best_cell) {
                                best_cell = ins_sc;
                                best_bt = 2;
                            }
                            curr_M[c] = best_cell;
                            bt_mat[b_idx(r, c)] = best_bt;
                        }
                        std::swap(prev_M, curr_M);
                        std::swap(prev_D, curr_D);
                    }

                    int best_c = Ll;
                    Scalar best_sc = prev_M[Ll];
                    for (int c = std::max(1, Ls - B0); c <= Ll; ++c) {
                        Scalar sc = prev_M[c] - std::min<Scalar>(cal_g_open + cal_g_ext * (Ll - c), 3.80);
                        if (sc > best_sc) { best_sc = sc; best_c = c; }
                    }

                    std::vector<std::pair<int8_t, int8_t>> matched_pairs;
                    matched_pairs.reserve(Lmin);
                    int cr = Ls, cc = best_c;
                    int aa_matches = 0;
                    while (cr > 0 && cc > 0) {
                        uint8_t mv = bt_mat[b_idx(cr, cc)];
                        if (mv == 0) {
                            int8_t c1 = s_short[cr - 1], c2 = s_long[cc - 1];
                            if (c1 >= 0 && c2 >= 0) {
                                matched_pairs.push_back({c1, c2});
                                if (c1 == c2 || gcode.sense_codon_aa[c1] == gcode.sense_codon_aa[c2]) aa_matches++;
                            }
                            cr--; cc--;
                        } else if (mv == 1) {
                            cr--;
                        } else {
                            cc--;
                        }
                    }
                    if (matched_pairs.size() < Lmin / 2) continue;
                    Scalar pair_aa_id = static_cast<Scalar>(aa_matches) / static_cast<Scalar>(matched_pairs.size());
                    if (pair_aa_id < 0.35) continue;
                    std::reverse(matched_pairs.begin(), matched_pairs.end());

                    struct AlignedSite {
                        int8_t c1, c2;
                        Scalar d_nt;
                        Scalar d_aa;
                    };
                    std::vector<AlignedSite> sites;
                    sites.reserve(matched_pairs.size());
                    Scalar sum_dnt = 0.0, sum_daa = 0.0;
                    for (const auto& [c1, c2] : matched_pairs) {
                        Scalar d_nt = 0.0, d_aa = 0.0;
                        if (c1 != c2) {
                            const auto& diff = gcode.diff_matrix[c1][c2];
                            d_nt = static_cast<Scalar>(diff.num_diffs) / 3.0;
                            if (!diff.synonymous) {
                                d_aa = aa_phys_dist(gcode.sense_codon_aa[c1], gcode.sense_codon_aa[c2]);
                            }
                        }
                        sites.push_back({c1, c2, d_nt, d_aa});
                        sum_dnt += d_nt;
                        sum_daa += d_aa;
                        cnt_all.add_codon_pair(c1, c2, gcode);
                    }
                    if (sites.size() < 5 || sum_dnt <= 1e-9) continue;
                    Scalar mean_dnt = sum_dnt / static_cast<Scalar>(sites.size());
                    Scalar mean_daa = sum_daa / static_cast<Scalar>(sites.size());

                    int M = static_cast<int>(sites.size());
                    for (int k = 0; k < M; ++k) {
                        int w0 = std::max(0, k - 2);
                        int w1 = std::min(M - 1, k + 2);
                        Scalar win_dnt = 0.0, win_daa = 0.0;
                        for (int w = w0; w <= w1; ++w) {
                            win_dnt += sites[w].d_nt;
                            win_daa += sites[w].d_aa;
                        }
                        win_dnt /= static_cast<Scalar>(w1 - w0 + 1);
                        win_daa /= static_cast<Scalar>(w1 - w0 + 1);
                        bool is_loop = (mean_daa > 1e-9) ? (win_daa > mean_daa) : (win_dnt > mean_dnt);
                        Scalar r_site = (mean_daa > 1e-9)
                            ? 0.5 * (sites[k].d_nt / mean_dnt + sites[k].d_aa / mean_daa)
                            : (sites[k].d_nt / mean_dnt);
                        if (!is_loop) {
                            n_core++;
                            sum_r_core += r_site;
                            cnt_core.add_codon_pair(sites[k].c1, sites[k].c2, gcode);
                        } else {
                            n_loop++;
                            sum_r_loop += r_site;
                            cnt_loop.add_codon_pair(sites[k].c1, sites[k].c2, gcode);
                        }
                    }
                }
            }

            EmpRegimeParams emp_all = cnt_all.finalize(mu0_daa, var0_daa, syn_to_ns_opp_ratio, opp_1_to_2, opp_1_to_3);
            m.emp_core = cnt_core.finalize(mu0_daa, var0_daa, syn_to_ns_opp_ratio, opp_1_to_2, opp_1_to_3);
            m.emp_loop = cnt_loop.finalize(mu0_daa, var0_daa, syn_to_ns_opp_ratio, opp_1_to_2, opp_1_to_3);
            m.emp_core.theta_ct = emp_all.theta_ct;
            m.emp_core.theta_tv = emp_all.theta_tv;
            m.emp_loop.theta_ct = emp_all.theta_ct;
            m.emp_loop.theta_tv = emp_all.theta_tv;

            if (n_core > 0 && n_loop > 0) {
                Scalar w_l = static_cast<Scalar>(n_loop) / static_cast<Scalar>(n_core + n_loop);
                Scalar r_c = sum_r_core / static_cast<Scalar>(n_core);
                Scalar r_l = sum_r_loop / static_cast<Scalar>(n_loop);
                m.w_loop = std::clamp<Scalar>(w_l, 0.10, 0.60);
                m.rho_loop_to_core = std::clamp<Scalar>(r_l / std::max<Scalar>(0.28, r_c), 1.5, 11.0);
            }

            m.q_core.build(m.emp_core, pos_freqs, pi_codon, gcode);
            m.q_loop.build(m.emp_loop, pos_freqs, pi_codon, gcode);
            return m;
        }

        static EmpiricalTwoStratumModel from_alignment(const Alignment& aln) {
            const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
            std::vector<std::vector<int8_t>> tc(aln.num_taxa);
            for (size_t t = 0; t < aln.num_taxa; ++t) {
                tc[t].reserve(aln.num_codons);
                for (int8_t c : aln.sequences[t].codons) {
                    if (c >= 0) tc[t].push_back(c);
                }
            }
            return calibrate(tc, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        }
    };

    // =========================================================================
    // Helper: Detect Indel-Ambiguous Codon Windows
    // =========================================================================
    static inline bool is_gap_codon(const Sequence& seq, size_t c) {
        if (3 * c + 3 <= seq.raw_seq.size()) {
            return seq.raw_seq[3 * c] == '-' && seq.raw_seq[3 * c + 1] == '-' && seq.raw_seq[3 * c + 2] == '-';
        }
        return seq.codons[c] < 0;
    }

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
                if (!is_gap_codon(aln.sequences[t], c)) {
                    if (first_ungapped[t] == C) first_ungapped[t] = c;
                    last_ungapped[t] = c;
                }
            }
        }

        std::vector<size_t> internal_gaps_per_site(C, 0);
        for (size_t c = 0; c < C; ++c) {
            for (size_t t = 0; t < N; ++t) {
                if (is_gap_codon(aln.sequences[t], c)) {
                    bool is_long_leading_missing = (c < first_ungapped[t] && first_ungapped[t] > 8);
                    bool is_long_trailing_missing = (c > last_ungapped[t] && (C - 1 - last_ungapped[t]) > 8);
                    if (!is_long_leading_missing && !is_long_trailing_missing) {
                        internal_gaps_per_site[c]++;
                    }
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

        // Populate window statistics and tile wide regions with 50% overlapping half-stride windows
        std::vector<IndelWindow> final_windows;
        size_t stride = std::max<size_t>(1, max_window_len / 2);
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
                        if (!is_gap_codon(aln.sequences[t], k)) {
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
                if (cur_e >= w.end_codon) break;
                cur_s += stride;
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
        Scalar conv_factor = 1.0,
        const EmpiricalTwoStratumModel* two_stratum = nullptr
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
                    P_shunt[t] = (two_stratum != nullptr)
                        ? two_stratum->transition_matrix(d_sh)
                        : mg94.transition_matrix(d_sh * conv_factor);
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
    static TransducerTaxonResult solve_codon_indel_transducer_for_taxon(
        const Sequence& seq,
        size_t start_codon,
        size_t end_codon,
        size_t taxon_idx,
        const std::vector<std::vector<Vector>>& leaf_cavity, // [w][taxon_idx]
        const JointAlignSettings& settings
    ) {
        size_t W = end_codon - start_codon;
        int S_dim = static_cast<int>(leaf_cavity[0][taxon_idx].size());
        uint64_t all_sense_mask = (S_dim < 64) ? ((1ULL << S_dim) - 1ULL) : ~0ULL;

        TransducerTaxonResult res;
        res.map_codons.resize(W);
        res.map_masks.resize(W, all_sense_mask);
        res.col_match_posterior.assign(W, 1.0);

        // Extract unaligned non-gap codons (including partial IUPAC ambiguous codons) inside the window
        std::vector<int8_t> obs_codons;
        std::vector<uint64_t> obs_masks;
        std::vector<std::string> obs_triplets;
        std::vector<size_t> orig_positions;
        for (size_t w = 0; w < W; ++w) {
            size_t c_idx = start_codon + w;
            int8_t c_state = seq.codons[c_idx];
            uint64_t c_mask = (c_idx < seq.codon_masks.size() && seq.codon_masks[c_idx] != 0)
                ? seq.codon_masks[c_idx]
                : ((c_state >= 0) ? (1ULL << c_state) : all_sense_mask);
            res.map_codons[w] = c_state;
            res.map_masks[w] = c_mask;
            if (!is_gap_codon(seq, c_idx)) {
                obs_codons.push_back(c_state);
                obs_masks.push_back(c_mask);
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

        // Precompute exact Bayesian homology log-odds emission score E(j, w):
        //   E(j, w) = ln( P_cavity(obs_codons[j] at w) / (1 / S) )
        // At an all-gap / uninformative column (cav = 1/S), E(j, w) = 0.0 (neutral).
        // At an informative column, E(j, w) > 0 iff obs_codons[j] is more probable under the
        // phylogenetic leaf cavity than under random background emission (preventing 1-codon
        // nonsynonymous stepping-stone jumps into long deletions).
        Eigen::MatrixXd log_E(M, W);
        std::vector<Scalar> col_gap_log_prior(W, 0.0);
        Scalar uniform_p = 1.0 / static_cast<Scalar>(std::max(1, S_dim));
        Scalar log_uniform_p = std::log(uniform_p);
        for (size_t w = 0; w < W; ++w) {
            const Vector& cav = leaf_cavity[w][taxon_idx];
            Scalar max_p = cav.maxCoeff();
            Scalar col_informativeness = std::clamp((max_p - uniform_p) / std::max<Scalar>(1e-6, 1.0 - uniform_p), 0.0, 1.0);

            for (size_t j = 0; j < M; ++j) {
                Scalar p_emit = 1e-9;
                if (obs_codons[j] >= 0) {
                    p_emit = cav(obs_codons[j]);
                } else {
                    uint64_t msk = obs_masks[j];
                    Scalar sum_p = 0.0;
                    int cnt_m = 0;
                    for (int s = 0; s < S_dim; ++s) {
                        if ((msk >> s) & 1ULL) {
                            sum_p += cav(s);
                            cnt_m++;
                        }
                    }
                    p_emit = (cnt_m > 0) ? (sum_p / static_cast<Scalar>(cnt_m)) : uniform_p;
                }
                Scalar log_p = std::log(std::max<Scalar>(1e-9, p_emit));
                log_E(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(w)) = (log_p - log_uniform_p);
            }
            col_gap_log_prior[w] = -settings.gap_extend_penalty * col_informativeness;
        }

        bool has_left_codon = false;
        for (size_t c = 0; c < start_codon; ++c) {
            if (!is_gap_codon(seq, c)) { has_left_codon = true; break; }
        }
        bool has_right_codon = false;
        for (size_t c = end_codon; c < seq.codons.size(); ++c) {
            if (!is_gap_codon(seq, c)) { has_right_codon = true; break; }
        }
        const bool left_in_gap = (start_codon > 0 && has_left_codon && is_gap_codon(seq, start_codon - 1));
        const bool right_in_gap = (end_codon < seq.codons.size() && has_right_codon && is_gap_codon(seq, end_codon));
        const Scalar right_m_boundary_pen = right_in_gap
            ? std::max<Scalar>(0.0, settings.gap_open_penalty - settings.gap_extend_penalty)
            : 0.0;

        auto score_placement = [&](const std::vector<size_t>& pos) -> Scalar {
            // pos[0] < pos[1] < ... < pos[M-1] in {0 .. W-1}
            Scalar score = 0.0;
            std::vector<uint8_t> is_match(W, 0);
            for (size_t j = 0; j < M; ++j) {
                is_match[pos[j]] = 1;
                score += log_E(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(pos[j]));
            }
            bool in_gap = left_in_gap;
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
            if (!in_gap) {
                score -= right_m_boundary_pen;
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

        if (left_in_gap) {
            fwd_G(0, 0) = 0.0;
            vit_G(0, 0) = 0.0;
        } else {
            fwd_M(0, 0) = 0.0;
            vit_M(0, 0) = 0.0;
        }

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

        fwd_M(W, M) -= right_m_boundary_pen;
        vit_M(W, M) -= right_m_boundary_pen;
        res.log_partition_Z = log_sum_exp(fwd_M(W, M), fwd_G(W, M));

        // Backward pass to compute exact per-column Match posterior probabilities
        Eigen::MatrixXd bwd_M = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);
        Eigen::MatrixXd bwd_G = Eigen::MatrixXd::Constant(W + 1, M + 1, NEG_INF);
        bwd_M(W, M) = -right_m_boundary_pen;
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
            res.map_masks.assign(W, all_sense_mask);
            std::string new_raw(3 * W, '-');
            for (size_t j = 0; j < M; ++j) {
                size_t col = map_positions[j];
                res.map_codons[col] = obs_codons[j];
                res.map_masks[col] = obs_masks[j];
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

        // Sort alternative placements (index >= 1) by descending posterior_weight, keeping index 0 as MAP
        if (res.top_candidates.size() > 2) {
            std::sort(res.top_candidates.begin() + 1, res.top_candidates.end(),
                      [](const TransducerTaxonResult::CandidatePlacement& a,
                         const TransducerTaxonResult::CandidatePlacement& b) {
                          return a.posterior_weight > b.posterior_weight;
                      });
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

        // Unify Stage 1 Phase 2 and Phase 3 under the empirical two-stratum physicochemical + multi-hit model:
        // P_b = (1 - \hat{w}_loop) \exp(t_b \hat{Q}_core) + \hat{w}_loop \exp(\hat{\rho}_{loop/core} t_b \hat{Q}_loop)
        EmpiricalTwoStratumModel two_stratum = EmpiricalTwoStratumModel::from_alignment(aln);

        std::vector<Matrix> P_branches(tree.num_nodes());
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                P_branches[node.id] = two_stratum.transition_matrix(node.branch_length);
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

        // Calibrate empirical phylogenetic branch-indel parameters (lam_open, eps_ext, pi_indel) from (tree, aln)
        size_t cal_gap_codons = 0;
        size_t cal_gap_runs = 0;
        for (size_t t = 0; t < N; ++t) {
            bool in_g = false;
            for (size_t c = 0; c < C; ++c) {
                if (is_gap_codon(aln.sequences[t], c)) {
                    cal_gap_codons++;
                    if (!in_g) {
                        cal_gap_runs++;
                        in_g = true;
                    }
                } else {
                    in_g = false;
                }
            }
        }
        Scalar total_tree_bl = 0.0;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && node.parent_id != INVALID_INDEX) {
                total_tree_bl += std::max<Scalar>(0.005, node.branch_length);
            }
        }
        Scalar mean_gap_len = (cal_gap_runs > 0)
            ? (static_cast<Scalar>(cal_gap_codons) / static_cast<Scalar>(cal_gap_runs))
            : 2.5;
        const Scalar emp_eps_ext = std::clamp<Scalar>(1.0 - 1.0 / std::max<Scalar>(1.25, mean_gap_len), 0.30, 0.85);
        const Scalar emp_lam_open = (cal_gap_runs > 0 && C > 0)
            ? std::clamp<Scalar>(
                  static_cast<Scalar>(cal_gap_runs) / std::max<Scalar>(1.0, static_cast<Scalar>(C) * total_tree_bl),
                  0.015, 0.12)
            : 0.035;
        const Scalar emp_pi_1 = std::clamp<Scalar>(
            1.0 - static_cast<Scalar>(cal_gap_codons) / std::max<Scalar>(1.0, static_cast<Scalar>(N * C)),
            0.55, 0.95
        );
        const Scalar emp_pi_0 = 1.0 - emp_pi_1;

        // Joint Two-Stratum Substitution + 2-State Felsenstein Affine-Memory Branch-Indel log-likelihood
        auto eval_joint_tree_indel_lnl = [&](const Alignment& a, size_t w_start, size_t w_end) -> Scalar {
            Scalar total = 0.0;
            std::vector<Vector> node_L(tree.num_nodes(), Vector::Zero(S));
            std::vector<double> node_scale(tree.num_nodes(), 0.0);
            for (size_t c = w_start; c < w_end; ++c) {
                SitePattern pat;
                pat.states.resize(N);
                for (size_t t = 0; t < N; ++t) {
                    pat.states[t] = a.sequences[t].codons[c];
                    if (pat.states[t] < 0 && !is_gap_codon(a.sequences[t], c) &&
                        c < a.sequences[t].codon_masks.size() && a.sequences[t].codon_masks[c] != 0) {
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

            // 2-State Felsenstein Sum-Product Marginalization with Horizontal Affine Memory:
            // q_10[cid] = P(I_{pa, w-1} = 1, I_{cid, w-1} = 0 | I_{., 1:w-1}) (posterior deletion on branch cid)
            // q_01[cid] = P(I_{pa, w-1} = 0, I_{cid, w-1} = 1 | I_{., 1:w-1}) (posterior insertion on branch cid)
            std::vector<Scalar> q_10(num_nodes, 0.0);
            std::vector<Scalar> q_01(num_nodes, 0.0);
            std::vector<Scalar> p_00(num_nodes, 1.0), p_01(num_nodes, 0.0);
            std::vector<Scalar> p_10(num_nodes, 0.0), p_11(num_nodes, 1.0);
            std::vector<std::array<Scalar, 2>> d_in(num_nodes, {1.0, 1.0});
            std::vector<std::array<Scalar, 2>> m_up(num_nodes, {1.0, 1.0});
            std::vector<std::array<Scalar, 2>> v_out(num_nodes, {0.0, 0.0});

            for (size_t w = 0; w < W_len; ++w) {
                // Step A: Build 2x2 branch transition matrices conditioned on horizontal affine memory from w-1
                for (int32_t nid : tree.post_order) {
                    if (nid == tree.root_id || tree.nodes[nid].parent_id == INVALID_INDEX) continue;
                    Scalar tb = std::max<Scalar>(0.005, tree.nodes[nid].branch_length);
                    Scalar p_open = 1.0 - std::exp(-emp_lam_open * tb);
                    p_10[nid] = (1.0 - q_10[nid]) * (0.5 * p_open) + q_10[nid] * emp_eps_ext;
                    p_01[nid] = (1.0 - q_01[nid]) * (0.5 * p_open) + q_01[nid] * emp_eps_ext;
                    p_11[nid] = std::max<Scalar>(1e-12, 1.0 - p_10[nid]);
                    p_00[nid] = std::max<Scalar>(1e-12, 1.0 - p_01[nid]);
                }

                // Step B: 2-State Felsenstein Inside (post-order bottom-up sum-product) pass
                for (int32_t nid : tree.post_order) {
                    const auto& node = tree.nodes[nid];
                    if (node.is_leaf) {
                        int32_t tid = (nid < static_cast<int32_t>(leaf_to_taxon.size())) ? leaf_to_taxon[nid] : -1;
                        bool present = (tid >= 0 && !is_gap_codon(a.sequences[tid], w_start + w));
                        d_in[nid][0] = present ? 0.0 : 1.0;
                        d_in[nid][1] = present ? 1.0 : 0.0;
                    } else {
                        Scalar acc0 = 1.0, acc1 = 1.0;
                        for (int32_t cid : node.children) {
                            m_up[cid][0] = p_00[cid] * d_in[cid][0] + p_01[cid] * d_in[cid][1];
                            m_up[cid][1] = p_10[cid] * d_in[cid][0] + p_11[cid] * d_in[cid][1];
                            acc0 *= m_up[cid][0];
                            acc1 *= m_up[cid][1];
                        }
                        d_in[nid][0] = acc0;
                        d_in[nid][1] = acc1;
                    }
                }

                Scalar L_col = emp_pi_0 * d_in[tree.root_id][0] + emp_pi_1 * d_in[tree.root_id][1];
                total += std::log(std::max<Scalar>(1e-300, L_col));

                // Step C: 2-State Felsenstein Outside (pre-order top-down) pass to update q_10 and q_01
                v_out[tree.root_id][0] = emp_pi_0;
                v_out[tree.root_id][1] = emp_pi_1;
                for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
                    int32_t uid = *it;
                    const auto& node = tree.nodes[uid];
                    if (node.is_leaf) continue;
                    size_t k_ch = node.children.size();
                    for (size_t i = 0; i < k_ch; ++i) {
                        int32_t cid = node.children[i];
                        Scalar v_sib0 = v_out[uid][0];
                        Scalar v_sib1 = v_out[uid][1];
                        for (size_t j = 0; j < k_ch; ++j) {
                            if (i != j) {
                                int32_t sj = node.children[j];
                                v_sib0 *= m_up[sj][0];
                                v_sib1 *= m_up[sj][1];
                            }
                        }
                        v_out[cid][0] = v_sib0 * p_00[cid] + v_sib1 * p_10[cid];
                        v_out[cid][1] = v_sib0 * p_01[cid] + v_sib1 * p_11[cid];

                        Scalar w00 = v_sib0 * p_00[cid] * d_in[cid][0];
                        Scalar w01 = v_sib0 * p_01[cid] * d_in[cid][1];
                        Scalar w10 = v_sib1 * p_10[cid] * d_in[cid][0];
                        Scalar w11 = v_sib1 * p_11[cid] * d_in[cid][1];
                        Scalar w_sum = w00 + w01 + w10 + w11;
                        if (w_sum > 1e-300) {
                            q_01[cid] = w01 / w_sum;
                            q_10[cid] = w10 / w_sum;
                        } else {
                            q_01[cid] = 0.0;
                            q_10[cid] = 0.0;
                        }
                    }
                }
            }
            return total;
        };

        auto eval_window_lnl = [&](const Alignment& a, size_t w_start, size_t w_end) -> Scalar {
            size_t m_start = (w_start > 0) ? (w_start - 1) : 0;
            size_t m_end = std::min(a.num_codons, w_end + 1);
            return eval_joint_tree_indel_lnl(a, m_start, m_end);
        };

        summaries.resize(windows.size());
        for (size_t wi = 0; wi < windows.size(); ++wi) {
            const auto& win = windows[wi];
            size_t w_start = win.start_codon;
            size_t w_end = win.end_codon;
            summaries[wi].window = win;
            summaries[wi].initial_window_log_l = eval_window_lnl(aln, w_start, w_end);
            summaries[wi].polished_window_log_l = summaries[wi].initial_window_log_l;
            summaries[wi].column_homology_posteriors.assign(w_end - w_start, 1.0);
            out_window_taxon_results[wi].resize(N);
        }

        auto step_single_window = [&](size_t wi) -> bool {
            const auto& win = windows[wi];
            size_t w_start = win.start_codon;
            size_t w_end = win.end_codon;
            size_t W = w_end - w_start;

            Scalar cur_lnl = eval_window_lnl(aln, w_start, w_end);
            std::vector<std::vector<Vector>> leaf_cavity;
            compute_window_leaf_cavities(tree, aln, mg94, P_branches, w_start, w_end, leaf_cavity, omega_ptr, conv, &two_stratum);

            Alignment batch_aln = aln;
            size_t batch_taxa_changed = 0;
            size_t batch_codons_shifted = 0;
            std::vector<TransducerTaxonResult> cand_taxon_res(N);
            struct TaxonProposal {
                size_t taxon_idx;
                Scalar gain;
            };
            std::vector<TaxonProposal> proposals;

            for (size_t t = 0; t < N; ++t) {
                auto tr = solve_codon_indel_transducer_for_taxon(
                    aln.sequences[t], w_start, w_end, t, leaf_cavity, settings
                );
                if (tr.changed) {
                    for (size_t w = 0; w < W; ++w) {
                        batch_aln.sequences[t].codons[w_start + w] = tr.map_codons[w];
                        if (w_start + w < batch_aln.sequences[t].codon_masks.size()) {
                            batch_aln.sequences[t].codon_masks[w_start + w] = tr.map_masks[w];
                        }
                    }
                    batch_aln.sequences[t].raw_seq.replace(3 * w_start, 3 * W, tr.map_raw_subseq);
                    batch_taxa_changed++;
                    batch_codons_shifted += tr.codons_shifted;
                    proposals.push_back({t, tr.map_log_score - tr.input_log_score});
                }
                cand_taxon_res[t] = std::move(tr);
            }

            if (out_window_taxon_results[wi][0].col_match_posterior.empty()) {
                out_window_taxon_results[wi] = cand_taxon_res;
            }

            if (proposals.empty()) return false;

            // Sort proposing taxa in descending order of transducer log-odds gain
            std::sort(proposals.begin(), proposals.end(),
                      [](const TaxonProposal& a, const TaxonProposal& b) { return a.gain > b.gain; });

            // Coordinate-wise Gauss-Seidel updates in descending gain order (preventing sister swaps
            // and rejecting noisy moves), compared against the coordinated clade batch move:
            Alignment seq_aln = aln;
            Scalar seq_lnl = cur_lnl;
            size_t seq_taxa_changed = 0;
            size_t seq_codons_shifted = 0;
            std::vector<TransducerTaxonResult> seq_taxon_res = cand_taxon_res;

            for (size_t pi_idx = 0; pi_idx < proposals.size(); ++pi_idx) {
                size_t t = proposals[pi_idx].taxon_idx;
                if (seq_taxa_changed > 0) {
                    compute_window_leaf_cavities(tree, seq_aln, mg94, P_branches, w_start, w_end, leaf_cavity, omega_ptr, conv, &two_stratum);
                    seq_taxon_res[t] = solve_codon_indel_transducer_for_taxon(
                        seq_aln.sequences[t], w_start, w_end, t, leaf_cavity, settings
                    );
                }
                if (!seq_taxon_res[t].changed) continue;

                Alignment trial_aln = seq_aln;
                for (size_t w = 0; w < W; ++w) {
                    trial_aln.sequences[t].codons[w_start + w] = seq_taxon_res[t].map_codons[w];
                    if (w_start + w < trial_aln.sequences[t].codon_masks.size()) {
                        trial_aln.sequences[t].codon_masks[w_start + w] = seq_taxon_res[t].map_masks[w];
                    }
                }
                trial_aln.sequences[t].raw_seq.replace(3 * w_start, 3 * W, seq_taxon_res[t].map_raw_subseq);
                Scalar trial_lnl = eval_window_lnl(trial_aln, w_start, w_end);
                if (trial_lnl > seq_lnl) {
                    seq_aln = std::move(trial_aln);
                    seq_lnl = trial_lnl;
                    seq_taxa_changed++;
                    seq_codons_shifted += seq_taxon_res[t].codons_shifted;
                }
            }

            Scalar batch_lnl = (batch_taxa_changed > 1) ? eval_window_lnl(batch_aln, w_start, w_end) : -1e30;
            if (seq_lnl >= batch_lnl && seq_lnl > cur_lnl) {
                aln = std::move(seq_aln);
                summaries[wi].polished_window_log_l = seq_lnl;
                summaries[wi].num_taxa_realigned += seq_taxa_changed;
                summaries[wi].num_codons_shifted += seq_codons_shifted;
                out_window_taxon_results[wi] = std::move(seq_taxon_res);
                return true;
            } else if (batch_lnl > cur_lnl) {
                aln = std::move(batch_aln);
                summaries[wi].polished_window_log_l = batch_lnl;
                summaries[wi].num_taxa_realigned += batch_taxa_changed;
                summaries[wi].num_codons_shifted += batch_codons_shifted;
                out_window_taxon_results[wi] = std::move(cand_taxon_res);
                return true;
            }
            return false;
        };

        // Group contiguous overlapping half-stride tiles into chains and polish with
        // bidirectional (forward + backward) Gauss-Seidel sweeps
        size_t chain_start = 0;
        while (chain_start < windows.size()) {
            size_t chain_end = chain_start;
            while (chain_end + 1 < windows.size() &&
                   windows[chain_end + 1].start_codon < windows[chain_end].end_codon) {
                chain_end++;
            }

            for (size_t iter = 0; iter < settings.em_iterations; ++iter) {
                bool any_accepted = false;
                // Forward pass: chain_start -> ... -> chain_end
                for (size_t wi = chain_start; wi <= chain_end; ++wi) {
                    if (step_single_window(wi)) any_accepted = true;
                }
                // Backward pass across overlapping tiles: chain_end - 1 -> ... -> chain_start
                if (chain_end > chain_start) {
                    for (size_t wi = chain_end; wi > chain_start; --wi) {
                        if (step_single_window(wi - 1)) any_accepted = true;
                    }
                }
                if (!any_accepted) break;
            }

            for (size_t wi = chain_start; wi <= chain_end; ++wi) {
                const auto& win = windows[wi];
                size_t w_start = win.start_codon;
                size_t w_end = win.end_codon;
                size_t W = w_end - w_start;
                summaries[wi].polished_window_log_l = eval_window_lnl(aln, w_start, w_end);
                summaries[wi].delta_log_l = summaries[wi].polished_window_log_l - summaries[wi].initial_window_log_l;

                for (size_t w = 0; w < W; ++w) {
                    Scalar acc_conf = 0.0;
                    size_t count = 0;
                    for (size_t t = 0; t < N; ++t) {
                        if (!out_window_taxon_results[wi][t].col_match_posterior.empty()) {
                            bool is_map_match = (aln.sequences[t].codons[w_start + w] >= 0);
                            Scalar p_m = out_window_taxon_results[wi][t].col_match_posterior[w];
                            Scalar conf = is_map_match ? p_m : (1.0 - p_m);
                            acc_conf += conf;
                            count++;
                        }
                    }
                    Scalar col_conf = (count > 0) ? (acc_conf / static_cast<Scalar>(count)) : 1.0;
                    summaries[wi].column_homology_posteriors[w] = col_conf;
                    out_site_homology_posteriors[w_start + w] = col_conf;
                }
            }

            chain_start = chain_end + 1;
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

        // Remove any completely empty (all-gap) codon columns evacuated during window realignment
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

        // Compute Foster effective resistance network on current alignment
        OhmicSettings ohm_s;
        ohm_s.use_hyphaeon_transformer = false;
        OhmicNetwork net = OhmicNetwork::infer(aln, ohm_s);
        size_t N = aln.num_taxa;
        const Matrix* omega_ptr = (net.Omega.rows() == static_cast<Eigen::Index>(N) &&
                                   net.Omega.cols() == static_cast<Eigen::Index>(N))
                                ? &net.Omega : nullptr;

        int S = gcode.num_sense_codons;

        // Calibrate empirical affine indel parameters on current alignment (matching Stage 1 Phase 3)
        size_t C = aln.num_codons;
        size_t cal_gap_runs = 0, cal_gap_ext = 0, cal_gap_codons = 0;
        for (size_t t = 0; t < N; ++t) {
            bool in_g = false;
            for (size_t c = 0; c < C; ++c) {
                bool g = is_gap_codon(aln.sequences[t], c);
                if (g) {
                    cal_gap_codons++;
                    if (!in_g) { cal_gap_runs++; in_g = true; }
                    else { cal_gap_ext++; }
                } else {
                    in_g = false;
                }
            }
        }
        Scalar total_tree_len = 0.0;
        for (const auto& nd : tree.nodes) {
            if (nd.id != tree.root_id && nd.parent_id != INVALID_INDEX) {
                total_tree_len += std::max<Scalar>(0.005, nd.branch_length);
            }
        }
        const Scalar emp_eps_ext = (cal_gap_runs + cal_gap_ext > 0)
            ? std::clamp<Scalar>(static_cast<Scalar>(cal_gap_ext) / static_cast<Scalar>(cal_gap_runs + cal_gap_ext), 0.25, 0.92)
            : 0.65;
        const Scalar emp_lam_open = (C > 0 && total_tree_len > 1e-4)
            ? std::clamp<Scalar>(static_cast<Scalar>(cal_gap_runs) / ( static_cast<Scalar>(C) * total_tree_len), 0.005, 0.25)
            : 0.035;
        const Scalar emp_pi_1 = std::clamp<Scalar>(
            1.0 - static_cast<Scalar>(cal_gap_codons) / std::max<Scalar>(1.0, static_cast<Scalar>(N * C)),
            0.55, 0.95
        );
        const Scalar emp_pi_0 = 1.0 - emp_pi_1;

        auto eval_window_indel_lnl = [&](const Tree& t_cand, const Alignment& a_cand, const std::vector<size_t>& l2t, size_t w_start, size_t w_end) -> Scalar {
            size_t m_start = (w_start > 0) ? (w_start - 1) : 0;
            size_t m_end = std::min(a_cand.num_codons, w_end + 1);
            size_t W_len = (m_end > m_start) ? (m_end - m_start) : 0;
            if (W_len == 0) return 0.0;
            size_t num_nodes = t_cand.num_nodes();

            std::vector<Scalar> q_10(num_nodes, 0.0), q_01(num_nodes, 0.0);
            std::vector<Scalar> p_00(num_nodes, 1.0), p_01(num_nodes, 0.0);
            std::vector<Scalar> p_10(num_nodes, 0.0), p_11(num_nodes, 1.0);
            std::vector<std::array<Scalar, 2>> d_in(num_nodes, {1.0, 1.0});
            std::vector<std::array<Scalar, 2>> m_up(num_nodes, {1.0, 1.0});
            std::vector<std::array<Scalar, 2>> v_out(num_nodes, {0.0, 0.0});

            Scalar indel_ll = 0.0;
            for (size_t w = 0; w < W_len; ++w) {
                for (int32_t nid : t_cand.post_order) {
                    if (nid == t_cand.root_id || t_cand.nodes[nid].parent_id == INVALID_INDEX) continue;
                    Scalar tb = std::max<Scalar>(0.005, t_cand.nodes[nid].branch_length);
                    Scalar p_open = 1.0 - std::exp(-emp_lam_open * tb);
                    p_10[nid] = (1.0 - q_10[nid]) * (0.5 * p_open) + q_10[nid] * emp_eps_ext;
                    p_01[nid] = (1.0 - q_01[nid]) * (0.5 * p_open) + q_01[nid] * emp_eps_ext;
                    p_11[nid] = std::max<Scalar>(1e-12, 1.0 - p_10[nid]);
                    p_00[nid] = std::max<Scalar>(1e-12, 1.0 - p_01[nid]);
                }
                for (int32_t nid : t_cand.post_order) {
                    const auto& node = t_cand.nodes[nid];
                    if (node.is_leaf) {
                        int32_t tid = (nid < static_cast<int32_t>(l2t.size())) ? static_cast<int32_t>(l2t[nid]) : -1;
                        bool present = (tid >= 0 && !is_gap_codon(a_cand.sequences[tid], m_start + w));
                        d_in[nid][0] = present ? 0.0 : 1.0;
                        d_in[nid][1] = present ? 1.0 : 0.0;
                    } else {
                        Scalar acc0 = 1.0, acc1 = 1.0;
                        for (int32_t cid : node.children) {
                            m_up[cid][0] = p_00[cid] * d_in[cid][0] + p_01[cid] * d_in[cid][1];
                            m_up[cid][1] = p_10[cid] * d_in[cid][0] + p_11[cid] * d_in[cid][1];
                            acc0 *= m_up[cid][0];
                            acc1 *= m_up[cid][1];
                        }
                        d_in[nid][0] = acc0;
                        d_in[nid][1] = acc1;
                    }
                }
                Scalar L_col = emp_pi_0 * d_in[t_cand.root_id][0] + emp_pi_1 * d_in[t_cand.root_id][1];
                indel_ll += std::log(std::max<Scalar>(1e-300, L_col));

                v_out[t_cand.root_id][0] = emp_pi_0;
                v_out[t_cand.root_id][1] = emp_pi_1;
                for (auto it = t_cand.post_order.rbegin(); it != t_cand.post_order.rend(); ++it) {
                    int32_t uid = *it;
                    const auto& node = t_cand.nodes[uid];
                    if (node.is_leaf) continue;
                    size_t k_ch = node.children.size();
                    for (size_t i = 0; i < k_ch; ++i) {
                        int32_t cid = node.children[i];
                        Scalar v_sib0 = v_out[uid][0], v_sib1 = v_out[uid][1];
                        for (size_t j = 0; j < k_ch; ++j) {
                            if (i != j) {
                                int32_t sj = node.children[j];
                                v_sib0 *= m_up[sj][0];
                                v_sib1 *= m_up[sj][1];
                            }
                        }
                        v_out[cid][0] = v_sib0 * p_00[cid] + v_sib1 * p_10[cid];
                        v_out[cid][1] = v_sib0 * p_01[cid] + v_sib1 * p_11[cid];

                        Scalar w00 = v_sib0 * p_00[cid] * d_in[cid][0];
                        Scalar w01 = v_sib0 * p_01[cid] * d_in[cid][1];
                        Scalar w10 = v_sib1 * p_10[cid] * d_in[cid][0];
                        Scalar w11 = v_sib1 * p_11[cid] * d_in[cid][1];
                        Scalar w_sum = w00 + w01 + w10 + w11;
                        if (w_sum > 1e-300) {
                            q_01[cid] = w01 / w_sum;
                            q_10[cid] = w10 / w_sum;
                        } else {
                            q_01[cid] = 0.0;
                            q_10[cid] = 0.0;
                        }
                    }
                }
            }
            return indel_ll;
        };

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
            for (const auto& w_sum : windows) {
                total += eval_window_indel_lnl(t_cand, a_cand, l2t, w_sum.window.start_codon, w_sum.window.end_codon);
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

        auto compute_candidate_clades = [&]() {
            std::vector<int32_t> cands;
            for (const auto& nd : tree.nodes) {
                if (nd.id < 0 || static_cast<size_t>(nd.id) >= node_leaves.size()) continue;
                if (nd.id == tree.root_id || nd.parent_id == INVALID_INDEX) continue;
                const auto& lvs = node_leaves[nd.id];
                if (lvs.empty() || lvs.size() > N / 2) continue;

                // Check if any taxon in lvs has internal gaps in any window
                bool has_indel = false;
                for (size_t tid : lvs) {
                    if (tid >= aln.sequences.size()) continue;
                    for (const auto& w : windows) {
                        for (size_t c = w.window.start_codon; c < w.window.end_codon && c < aln.sequences[tid].codons.size(); ++c) {
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
                    cands.push_back(nd.id);
                }
            }
            // Sort candidate clades by branch length descending (highest Ohmic tension first)
            std::sort(cands.begin(), cands.end(), [&](int32_t a, int32_t b) {
                return tree.nodes[a].branch_length > tree.nodes[b].branch_length;
            });
            return cands;
        };

        std::vector<int32_t> candidate_clades = compute_candidate_clades();

        size_t max_clades_to_probe = std::min<size_t>(candidate_clades.size(), 16);
        for (size_t ci = 0; ci < max_clades_to_probe; ++ci) {
            if (ci >= candidate_clades.size()) break;
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
                        Scalar r = net.taxon_effective_resistance(leaves_u[iu], node_leaves[v][iv]);
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
                std::vector<int32_t> local_brs;
                for (int32_t bid : {spr_opt->u, spr_opt->v, spr_opt->w, spr_opt->s}) {
                    if (bid != INVALID_INDEX && bid >= 0 &&
                        static_cast<size_t>(bid) < cand_tree.nodes.size() &&
                        bid != cand_tree.root_id &&
                        cand_tree.nodes[bid].parent_id != INVALID_INDEX) {
                        local_brs.push_back(bid);
                    }
                }
                TreeEngine::optimize_branches_gtr(cand_tree, nuc_aln, gtr, 3, 1e-2, &local_brs);

                // Score uncoupled topological change (including indel likelihood on current alignment)
                Scalar sub_only_ll = eval_total_mg94_lnl(cand_tree, aln);
                Scalar delta_sub = sub_only_ll - current_mg94_ll;

                // If the move is within striking distance (delta_sub > -15.0), run coupled boundary realignment
                // on the moved clade leaves_u under the new topology cand_tree!
                if (delta_sub > -15.0) {
                    Alignment cand_aln = aln;
                    size_t clade_realigned = 0;

                    if (!windows.empty()) {
                        EmpiricalTwoStratumModel two_stratum = EmpiricalTwoStratumModel::from_alignment(cand_aln);
                        std::vector<Matrix> P_br(cand_tree.num_nodes());
                        for (const auto& nd : cand_tree.nodes) {
                            if (nd.id != cand_tree.root_id && nd.parent_id != INVALID_INDEX) {
                                P_br[nd.id] = two_stratum.transition_matrix(nd.branch_length);
                            }
                        }

                        for (const auto& w_sum : windows) {
                            size_t w_start = w_sum.window.start_codon;
                            size_t w_end = w_sum.window.end_codon;
                            size_t W = w_end - w_start;

                            std::vector<std::vector<Vector>> cav;
                            compute_window_leaf_cavities(cand_tree, cand_aln, mg94, P_br, w_start, w_end, cav, omega_ptr, conv, &two_stratum);

                            for (size_t tid : leaves_u) {
                                if (tid >= cand_aln.sequences.size()) continue;
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

                    Alignment cand_nuc_aln;
                    if (clade_realigned > 0) {
                        cand_aln.finalize(cand_aln.code);
                        cand_nuc_aln = cand_aln.to_nucleotide_alignment();
                        TreeEngine::optimize_branches_gtr(cand_tree, cand_nuc_aln, gtr, 3, 1e-2, &local_brs);
                    }

                    Scalar joint_ll = (clade_realigned > 0) ? eval_total_mg94_lnl(cand_tree, cand_aln) : sub_only_ll;
                    Scalar delta_joint = joint_ll - current_mg94_ll;

                    if (delta_joint > 0.10) {
                        current_mg94_ll = joint_ll;
                        tree = std::move(cand_tree);
                        if (clade_realigned > 0) {
                            aln = std::move(cand_aln);
                            nuc_aln = std::move(cand_nuc_aln);
                        }
                        node_leaves = compute_node_leaves();
                        candidate_clades = compute_candidate_clades();
                        max_clades_to_probe = std::min<size_t>(candidate_clades.size(), 16);
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

        // 1. Run FEL on the Stage 1+2 MAP realigned alignment
        FELAnalyzer fel_map(tree, map_aln, mg_params, fel_s);
        fel_map.run();

        // 2. Only run a separate FEL pass on initial_aln if Stage 1+2 modified at least one sequence
        bool any_seq_changed = (initial_aln.num_codons != map_aln.num_codons);
        if (!any_seq_changed) {
            for (size_t t = 0; t < N; ++t) {
                if (initial_aln.sequences[t].codons != map_aln.sequences[t].codons) {
                    any_seq_changed = true;
                    break;
                }
            }
        }
        std::optional<FELAnalyzer> fel_rigid_opt;
        if (any_seq_changed && initial_aln.num_codons == map_aln.num_codons) {
            fel_rigid_opt.emplace(tree, initial_aln, mg_params, fel_s);
            fel_rigid_opt->run();
        }
        const auto& rigid_site_results = fel_rigid_opt.has_value() ? fel_rigid_opt->site_results : fel_map.site_results;

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

            const auto& sr_rigid = rigid_site_results[c];
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
        if (input_aln.num_taxa < 3) {
            throw std::runtime_error(
                "Gate 1 (Input Size Gate): Phylogenetic selection analysis requires at least 3 homologous coding sequences (received " +
                std::to_string(input_aln.num_taxa) + ")."
            );
        }
        {
            std::unordered_set<std::string> seen_names;
            for (const auto& sq : input_aln.sequences) {
                if (!seen_names.insert(sq.name).second) {
                    throw std::runtime_error(
                        "Gate 1 (Duplicate Taxon Identifier): Input alignment contains duplicate sequence identifier '" +
                        sq.name + "'. Every sequence must have a unique taxon name."
                    );
                }
            }
        }

        JointAlignResult res;
        res.initial_alignment = input_aln;
        res.map_alignment = input_aln;

        const auto& gcode_adm = *(input_aln.code ? input_aln.code : GeneticCode::universal());
        auto is_raw_stop_triplet = [&](const std::string& raw, size_t c) -> bool {
            if (3 * c + 2 >= raw.size()) return false;
            char c1 = raw[3 * c], c2 = raw[3 * c + 1], c3 = raw[3 * c + 2];
            if (c1 == '-' || c2 == '-' || c3 == '-' || c1 == 'N' || c2 == 'N' || c3 == 'N' ||
                c1 == '?' || c2 == '?' || c3 == '?') {
                return false;
            }
            return gcode_adm.is_stop_triplet(c1, c2, c3);
        };

        // Check if aligned input has a native stop column in the trailing 15% followed by 3' UTR codons (E2)
        size_t utr_trim_col = res.map_alignment.num_codons;
        if (res.map_alignment.num_codons >= 50) {
            size_t tail_start = std::max<size_t>(40, (res.map_alignment.num_codons * 85) / 100);
            for (size_t c = tail_start; c < res.map_alignment.num_codons; ++c) {
                size_t stop_taxa = 0, non_gap_taxa = 0;
                for (size_t t = 0; t < res.map_alignment.num_taxa; ++t) {
                    if (3 * c + 2 < res.map_alignment.sequences[t].raw_seq.size()) {
                        std::string trip = res.map_alignment.sequences[t].raw_seq.substr(3 * c, 3);
                        if (trip != "---" && trip != "NNN" && trip != "???") {
                            non_gap_taxa++;
                            if (is_raw_stop_triplet(res.map_alignment.sequences[t].raw_seq, c)) stop_taxa++;
                        }
                    }
                }
                if (non_gap_taxa >= 3 && stop_taxa * 4 >= non_gap_taxa) {
                    utr_trim_col = c;
                    break;
                }
            }
            if (utr_trim_col < res.map_alignment.num_codons) {
                auto trim_aln_to = [&](Alignment& a, size_t new_c) {
                    for (auto& sq : a.sequences) {
                        if (sq.codons.size() > new_c) sq.codons.resize(new_c);
                        if (sq.codon_masks.size() > new_c) sq.codon_masks.resize(new_c);
                        if (sq.raw_seq.size() > 3 * new_c) sq.raw_seq.resize(3 * new_c);
                    }
                    a.num_codons = new_c;
                    a.finalize(a.code);
                };
                trim_aln_to(res.initial_alignment, utr_trim_col);
                trim_aln_to(res.map_alignment, utr_trim_col);
            }
        }

        // Populate admission audit summary and enforce Gate 2 reading-frame check for aligned inputs
        const Alignment& cur_aln = res.map_alignment;
        res.admission_report.total_input_sequences = cur_aln.num_taxa;
        res.admission_report.admitted_sequences = cur_aln.num_taxa;
        res.admission_report.entries.resize(cur_aln.num_taxa);
        size_t aligned_stop_failures = 0;
        for (size_t t = 0; t < cur_aln.num_taxa; ++t) {
            auto& ent = res.admission_report.entries[t];
            ent.name = cur_aln.sequences[t].name;
            ent.raw_nt_length = input_aln.sequences[t].raw_seq.size();
            size_t valid_cods = 0, stops = 0;
            for (size_t c = 0; c < cur_aln.num_codons; ++c) {
                int8_t sc = cur_aln.sequences[t].codons[c];
                if (sc >= 0) valid_cods++;
                else if (c + 1 < cur_aln.num_codons && is_raw_stop_triplet(cur_aln.sequences[t].raw_seq, c)) {
                    stops++;
                }
            }
            ent.codon_length = valid_cods;
            ent.internal_stops_frame0 = stops;
            if (stops > 0) {
                aligned_stop_failures++;
                ent.reason = "Aligned sequence (" + std::to_string(stops) + " internal stop(s) masked)";
            } else if (utr_trim_col < input_aln.num_codons) {
                ent.reason = "Aligned in-frame CDS (3' UTR trimmed at codon " + std::to_string(utr_trim_col + 1) + ")";
            } else if (cur_aln.num_codons > 0 && valid_cods * 2 < cur_aln.num_codons) {
                ent.reason = "Partial coverage (" + std::to_string(valid_cods) + "/" + std::to_string(cur_aln.num_codons) + " codons)";
            } else {
                ent.reason = "Clean in-frame CDS";
            }
        }

        if (settings.enforce_admission_control &&
            (cur_aln.num_taxa - aligned_stop_failures < 3 ||
             static_cast<double>(aligned_stop_failures) / static_cast<double>(cur_aln.num_taxa) > 0.25)) {
            std::string code_hint;
            const std::vector<std::string> alt_codes = {
                "Vertebrate-mtDNA", "Invertebrate-mtDNA", "Yeast-mtDNA", "Mold-Protozoan-mtDNA",
                "Ciliate-Nuclear", "Euplotid-Nuclear", "Echinoderm-mtDNA", "Ascidian-mtDNA",
                "Flatworm-mtDNA", "Blepharisma-Nuclear", "Chlorophycean-mtDNA", "Trematode-mtDNA"
            };
            size_t best_clean = cur_aln.num_taxa - aligned_stop_failures;
            std::string best_code_name;
            for (const auto& cname : alt_codes) {
                auto alt_gc = GeneticCode::from_name(cname);
                if (!alt_gc) continue;
                size_t clean_cnt = 0;
                for (size_t t = 0; t < cur_aln.num_taxa; ++t) {
                    bool has_st = false;
                    for (size_t c = 0; c + 1 < cur_aln.num_codons; ++c) {
                        if (3 * c + 2 >= cur_aln.sequences[t].raw_seq.size()) continue;
                        char c1 = cur_aln.sequences[t].raw_seq[3 * c];
                        char c2 = cur_aln.sequences[t].raw_seq[3 * c + 1];
                        char c3 = cur_aln.sequences[t].raw_seq[3 * c + 2];
                        if (c1 == '-' || c2 == '-' || c3 == '-' || c1 == 'N' || c2 == 'N' || c3 == 'N') continue;
                        if (alt_gc->is_stop_triplet(c1, c2, c3)) { has_st = true; break; }
                    }
                    if (!has_st) clean_cnt++;
                }
                if (clean_cnt > best_clean && clean_cnt >= 3 &&
                    static_cast<double>(cur_aln.num_taxa - clean_cnt) / static_cast<double>(cur_aln.num_taxa) <= 0.25) {
                    best_clean = clean_cnt;
                    best_code_name = cname;
                }
            }
            if (!best_code_name.empty()) {
                code_hint = " Hint: Re-running with '--code " + best_code_name + "' resolves internal stop codons in " +
                            std::to_string(best_clean) + "/" + std::to_string(cur_aln.num_taxa) + " sequences.";
            } else {
                code_hint = " Hint: If sequences are out-of-frame or unaligned, re-run with '-u' (--unaligned) for automatic reading-frame rescue.";
            }
            std::ostringstream err;
            err << "Gate 2 (Aligned Reading-Frame & Stop-Codon Gate): " << aligned_stop_failures << " of " << cur_aln.num_taxa
                << " aligned sequences (" << std::fixed << std::setprecision(1)
                << (100.0 * static_cast<double>(aligned_stop_failures) / static_cast<double>(cur_aln.num_taxa))
                << "%) contain internal premature stop codons in frame +0." << code_hint;
            throw std::runtime_error(err.str());
        }

        // Compute pairwise amino-acid identity and homology Z-scores across the aligned cohort for the audit report
        {
            std::vector<Scalar> pair_ids, pair_zs;
            for (size_t i = 0; i < cur_aln.num_taxa; ++i) {
                Scalar max_id = 0.0, max_z = 0.0;
                for (size_t j = 0; j < cur_aln.num_taxa; ++j) {
                    if (i == j) continue;
                    size_t shared = 0, same_aa = 0;
                    for (size_t c = 0; c < cur_aln.num_codons; ++c) {
                        int8_t ci = cur_aln.sequences[i].codons[c];
                        int8_t cj = cur_aln.sequences[j].codons[c];
                        if (ci >= 0 && cj >= 0) {
                            shared++;
                            if (gcode_adm.sense_codon_aa[ci] == gcode_adm.sense_codon_aa[cj]) same_aa++;
                        }
                    }
                    if (shared > 0) {
                        Scalar id_ij = static_cast<Scalar>(same_aa) / static_cast<Scalar>(shared);
                        Scalar p0 = 0.06;
                        Scalar sd0 = std::sqrt(p0 * (1.0 - p0) / static_cast<Scalar>(shared));
                        Scalar z_ij = (id_ij - p0) / std::max<Scalar>(1e-4, sd0);
                        max_id = std::max(max_id, id_ij);
                        max_z = std::max(max_z, z_ij);
                        if (j > i) {
                            pair_ids.push_back(id_ij);
                            pair_zs.push_back(z_ij);
                        }
                    }
                }
                res.admission_report.entries[i].max_aa_identity = max_id;
                res.admission_report.entries[i].max_homology_z = max_z;
            }
            if (!pair_ids.empty()) {
                std::sort(pair_ids.begin(), pair_ids.end());
                std::sort(pair_zs.begin(), pair_zs.end());
                res.admission_report.median_aa_identity = pair_ids[pair_ids.size() / 2];
                res.admission_report.median_pairwise_z = pair_zs[pair_zs.size() / 2];
            }
        }

        // 0. Build or optimize initial tree via Ohmic TreeEngine if not provided
        if (!settings.quiet) {
            Panel::print_step(
                1, settings.run_marginalized_fel ? 5 : 4,
                "Initial Phylogeny & Global Codon MG94 Fit",
                std::to_string(res.map_alignment.num_taxa) + " taxa × " + std::to_string(res.map_alignment.num_codons) + " codons"
            );
        }
        bool user_supplied_tree = false;
        if (starting_tree.has_value()) {
            res.initial_tree = *starting_tree;
            user_supplied_tree = true;
        } else if (!input_aln.embedded_tree_newick.empty()) {
            res.initial_tree = Tree::from_newick(input_aln.embedded_tree_newick);
            user_supplied_tree = true;
        } else {
            TreeEngineSettings tes;
            tes.quiet = true;
            if (res.map_alignment.num_taxa > settings.subquadratic_taxon_threshold) {
                tes.perform_spr = false;
                tes.perform_topological_gradients = false;
                tes.lbfgs_max_iterations = 15;
            }
            TreeEngine te(tes);
            auto te_res = te.run(res.map_alignment);
            res.initial_tree = te_res.optimized_tree;
        }

        // Sanitize numeric support values on internal nodes (e.g. FastTree "0.776") so branch names are unambiguous (E19)
        for (auto& nd : res.initial_tree.nodes) {
            if (!nd.is_leaf && !nd.name.empty()) {
                bool is_numeric = true;
                for (char ch : nd.name) {
                    if (!std::isdigit(static_cast<unsigned char>(ch)) && ch != '.' && ch != '-' && ch != '+') {
                        is_numeric = false;
                        break;
                    }
                }
                if (is_numeric) {
                    nd.name = "Node_" + std::to_string(nd.id);
                }
            }
        }

        // Validate leaf correspondence when a tree was explicitly supplied or embedded (T1-08 / T5-10 / #34)
        if (user_supplied_tree) {
            std::vector<std::string> unmatched_leaves;
            std::vector<std::string> matched_tree_leaf_names;
            std::vector<bool> taxon_matched(res.initial_alignment.num_taxa, false);
            // Pass 1: exact and case-insensitive matches
            for (const auto& nd : res.initial_tree.nodes) {
                if (!nd.is_leaf) continue;
                std::string clean_name = nd.name;
                if (clean_name.size() >= 2 &&
                    ((clean_name.front() == '\'' && clean_name.back() == '\'') ||
                     (clean_name.front() == '"' && clean_name.back() == '"'))) {
                    clean_name = clean_name.substr(1, clean_name.size() - 2);
                }
                size_t tid = static_cast<size_t>(INVALID_INDEX);
                auto it = res.initial_alignment.taxon_to_index.find(clean_name);
                if (it != res.initial_alignment.taxon_to_index.end()) {
                    tid = it->second;
                } else {
                    std::string lower = clean_name;
                    std::transform(lower.begin(), lower.end(), lower.begin(),
                                   [](unsigned char c) { return std::tolower(c); });
                    auto it_l = res.initial_alignment.lower_taxon_to_index.find(lower);
                    if (it_l != res.initial_alignment.lower_taxon_to_index.end()) {
                        tid = it_l->second;
                    }
                }
                if (tid != static_cast<size_t>(INVALID_INDEX) && tid < res.initial_alignment.num_taxa && !taxon_matched[tid]) {
                    taxon_matched[tid] = true;
                    matched_tree_leaf_names.push_back(nd.name);
                }
            }
            // Pass 2: unique prefix fallback only for remaining unmatched leaves and unmatched alignment taxa
            for (const auto& nd : res.initial_tree.nodes) {
                if (!nd.is_leaf) continue;
                if (std::find(matched_tree_leaf_names.begin(), matched_tree_leaf_names.end(), nd.name) !=
                    matched_tree_leaf_names.end()) {
                    continue;
                }
                size_t tid = LikelihoodEngine::find_taxon_index(res.initial_alignment, nd.name);
                if (tid == static_cast<size_t>(INVALID_INDEX) || tid >= res.initial_alignment.num_taxa || taxon_matched[tid]) {
                    unmatched_leaves.push_back(nd.name);
                } else {
                    taxon_matched[tid] = true;
                    matched_tree_leaf_names.push_back(nd.name);
                }
            }
            std::vector<std::string> missing_taxa;
            for (size_t t = 0; t < res.initial_alignment.num_taxa; ++t) {
                if (!taxon_matched[t]) missing_taxa.push_back(res.initial_alignment.sequences[t].name);
            }
            if (!missing_taxa.empty()) {
                std::ostringstream err;
                err << "Gate 1 (Tree-Alignment Taxon Mismatch): Supplied phylogenetic tree does not match alignment taxa."
                    << " Alignment taxa absent from tree (" << missing_taxa.size() << "/"
                    << res.initial_alignment.num_taxa << "): ";
                for (size_t i = 0; i < std::min<size_t>(missing_taxa.size(), 10); ++i) {
                    if (i > 0) err << ", ";
                    err << missing_taxa[i];
                }
                if (missing_taxa.size() > 10) {
                    err << ", ... (+" << (missing_taxa.size() - 10) << " more)";
                }
                err << ".";
                if (!unmatched_leaves.empty()) {
                    err << " Tree leaves absent from alignment (" << unmatched_leaves.size() << "): " << unmatched_leaves[0] << ".";
                }
                throw std::runtime_error(err.str());
            }
            if (!unmatched_leaves.empty()) {
                if (!settings.quiet) {
                    std::ostringstream warn;
                    warn << "      Note: Pruning " << unmatched_leaves.size()
                         << " extra tree tip(s) absent from the " << res.initial_alignment.num_taxa
                         << "-taxon alignment (e.g. " << unmatched_leaves[0] << ").\n";
                    std::cout << Console::warning(warn.str()) << std::flush;
                }
                res.initial_tree = res.initial_tree.prune_to_leaves(matched_tree_leaf_names);
            }
        }

        LikelihoodEngine::reconcile_tree_leaf_names(res.initial_tree, res.initial_alignment);

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
            res.window_taxon_results = std::move(window_taxon_results);

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

    // =========================================================================
    // Subquadratic O(N L) Dayhoff-6 & Amino-Acid k-mer Bit-Sketch Preprocessing
    // and Farthest-Point Landmark + k-NN Graph Construction for Large Alignments
    // =========================================================================
    struct LandmarkKmerGraph {
        size_t num_taxa = 0;
        size_t num_landmarks = 0;
        std::vector<size_t> landmarks;                  // K landmark indices in [0, N)
        std::vector<bool> is_landmark;                  // size N
        std::vector<std::array<Scalar, 20>> aa_freqs;   // size N: amino-acid composition vector f_u(A)
        std::vector<std::vector<size_t>> sparse_adj;    // size N: sorted neighbor indices v > u in E_sparse

        static int aa_to_index20(char aa) {
            switch (aa) {
                case 'A': return 0;  case 'R': return 1;  case 'N': return 2;  case 'D': return 3;
                case 'C': return 4;  case 'Q': return 5;  case 'E': return 6;  case 'G': return 7;
                case 'H': return 8;  case 'I': return 9;  case 'L': return 10; case 'K': return 11;
                case 'M': return 12; case 'F': return 13; case 'P': return 14; case 'S': return 15;
                case 'T': return 16; case 'W': return 17; case 'Y': return 18; case 'V': return 19;
                default:  return -1;
            }
        }

        static int aa_to_dayhoff6(char aa) {
            switch (aa) {
                case 'A': case 'G': case 'P': case 'S': case 'T': return 0;
                case 'D': case 'E': case 'N': case 'Q':           return 1;
                case 'H': case 'K': case 'R':                     return 2;
                case 'I': case 'L': case 'M': case 'V':           return 3;
                case 'F': case 'W': case 'Y':                     return 4;
                case 'C':                                         return 5;
                default:                                          return -1;
            }
        }

        static Scalar sketch_jaccard(
            const std::array<uint64_t, 128>& a,
            const std::array<uint64_t, 128>& b
        ) {
            size_t inter = 0, uni = 0;
            for (size_t w = 0; w < 128; ++w) {
                inter += static_cast<size_t>(__builtin_popcountll(a[w] & b[w]));
                uni   += static_cast<size_t>(__builtin_popcountll(a[w] | b[w]));
            }
            return (uni > 0) ? (static_cast<Scalar>(inter) / static_cast<Scalar>(uni)) : 0.0;
        }

        Scalar expected_aa_identity(size_t u, size_t v) const {
            Scalar dot = 0.0;
            for (int a = 0; a < 20; ++a) {
                dot += aa_freqs[u][a] * aa_freqs[v][a];
            }
            return std::clamp<Scalar>(dot, 0.04, 0.25);
        }

        static LandmarkKmerGraph build(
            const std::vector<std::vector<int8_t>>& taxon_codons,
            const GeneticCode& gcode,
            size_t knn_degree = 8
        ) {
            LandmarkKmerGraph g;
            size_t N = taxon_codons.size();
            g.num_taxa = N;
            g.is_landmark.assign(N, false);
            g.aa_freqs.resize(N);
            g.sparse_adj.resize(N);
            if (N == 0) return g;

            std::vector<std::array<uint64_t, 128>> sketches(N);
            std::vector<size_t> lengths(N, 0);
            for (size_t u = 0; u < N; ++u) {
                sketches[u].fill(0ULL);
                g.aa_freqs[u].fill(0.0);
                const auto& seq = taxon_codons[u];
                lengths[u] = seq.size();
                std::vector<int> aa_idx;
                std::vector<int> d6_idx;
                aa_idx.reserve(seq.size());
                d6_idx.reserve(seq.size());
                size_t valid_aa = 0;
                for (int8_t c : seq) {
                    if (c >= 0 && c < gcode.num_sense_codons) {
                        char aa = gcode.sense_codon_aa[c];
                        int ai = aa_to_index20(aa);
                        int di = aa_to_dayhoff6(aa);
                        aa_idx.push_back(ai);
                        d6_idx.push_back(di);
                        if (ai >= 0) {
                            g.aa_freqs[u][static_cast<size_t>(ai)] += 1.0;
                            valid_aa++;
                        }
                    } else {
                        aa_idx.push_back(-1);
                        d6_idx.push_back(-1);
                    }
                }
                if (valid_aa > 0) {
                    for (int a = 0; a < 20; ++a) {
                        g.aa_freqs[u][static_cast<size_t>(a)] /= static_cast<Scalar>(valid_aa);
                    }
                } else {
                    g.aa_freqs[u].fill(0.05);
                }

                // Words 0..63: exact AA 3-mer bits
                for (size_t i = 0; i + 2 < aa_idx.size(); ++i) {
                    if (aa_idx[i] >= 0 && aa_idx[i + 1] >= 0 && aa_idx[i + 2] >= 0) {
                        uint32_t code = static_cast<uint32_t>(aa_idx[i] * 400 + aa_idx[i + 1] * 20 + aa_idx[i + 2]);
                        uint32_t h = (code * 2654435761u) >> 20; // 0..4095
                        sketches[u][h >> 6] |= (1ULL << (h & 63));
                    }
                }
                // Words 64..127: Dayhoff-6 physicochemical 5-mer bits
                for (size_t i = 0; i + 4 < d6_idx.size(); ++i) {
                    if (d6_idx[i] >= 0 && d6_idx[i + 1] >= 0 && d6_idx[i + 2] >= 0 &&
                        d6_idx[i + 3] >= 0 && d6_idx[i + 4] >= 0) {
                        uint32_t code = static_cast<uint32_t>(
                            d6_idx[i] * 1296 + d6_idx[i + 1] * 216 + d6_idx[i + 2] * 36 + d6_idx[i + 3] * 6 + d6_idx[i + 4]
                        );
                        uint32_t h = (code * 2654435761u) >> 20; // 0..4095
                        sketches[u][64 + (h >> 6)] |= (1ULL << (h & 63));
                    }
                }
            }

            // Select K = min(N, max(32, 2 * ceil(sqrt(N)))) landmarks via farthest-point traversal
            size_t K = std::min<size_t>(N, std::max<size_t>(32, 2 * static_cast<size_t>(std::ceil(std::sqrt(static_cast<double>(N))))));
            std::vector<size_t> sorted_lens = lengths;
            std::sort(sorted_lens.begin(), sorted_lens.end());
            size_t med_len = sorted_lens[N / 2];

            size_t first_lm = 0;
            size_t best_len_diff = std::numeric_limits<size_t>::max();
            for (size_t u = 0; u < N; ++u) {
                size_t diff = (lengths[u] >= med_len) ? (lengths[u] - med_len) : (med_len - lengths[u]);
                if (diff < best_len_diff) {
                    best_len_diff = diff;
                    first_lm = u;
                }
            }

            g.landmarks.push_back(first_lm);
            g.is_landmark[first_lm] = true;
            std::vector<Scalar> max_sim_to_lm(N, 0.0);
            for (size_t u = 0; u < N; ++u) {
                max_sim_to_lm[u] = sketch_jaccard(sketches[u], sketches[first_lm]);
            }

            while (g.landmarks.size() < K) {
                size_t next_lm = N;
                Scalar min_max_sim = 2.0;
                for (size_t u = 0; u < N; ++u) {
                    if (g.is_landmark[u] || lengths[u] == 0) continue;
                    if (max_sim_to_lm[u] < min_max_sim) {
                        min_max_sim = max_sim_to_lm[u];
                        next_lm = u;
                    }
                }
                if (next_lm >= N) break;
                g.landmarks.push_back(next_lm);
                g.is_landmark[next_lm] = true;
                for (size_t u = 0; u < N; ++u) {
                    if (!g.is_landmark[u]) {
                        max_sim_to_lm[u] = std::max(max_sim_to_lm[u], sketch_jaccard(sketches[u], sketches[next_lm]));
                    }
                }
            }
            g.num_landmarks = g.landmarks.size();

            // Build symmetric sparse adjacency list (all taxon-landmark edges + top knn_degree k-mer neighbors per taxon)
            std::vector<std::vector<size_t>> knn_lists(N);
            #pragma omp parallel for schedule(static)
            for (size_t u = 0; u < N; ++u) {
                std::vector<std::pair<Scalar, size_t>> cand;
                cand.reserve(N);
                for (size_t v = 0; v < N; ++v) {
                    if (u == v) continue;
                    cand.push_back({sketch_jaccard(sketches[u], sketches[v]), v});
                }
                size_t k_take = std::min(knn_degree, cand.size());
                std::partial_sort(cand.begin(), cand.begin() + static_cast<ptrdiff_t>(k_take), cand.end(),
                    [](const auto& a, const auto& b) {
                        if (a.first != b.first) return a.first > b.first;
                        return a.second < b.second;
                    });
                knn_lists[u].reserve(k_take);
                for (size_t k = 0; k < k_take; ++k) {
                    knn_lists[u].push_back(cand[k].second);
                }
            }

            for (size_t u = 0; u < N; ++u) {
                if (g.is_landmark[u]) {
                    for (size_t v = u + 1; v < N; ++v) {
                        g.sparse_adj[u].push_back(v);
                    }
                } else {
                    std::vector<size_t> nbrs;
                    for (size_t lm : g.landmarks) {
                        if (lm > u) nbrs.push_back(lm);
                    }
                    for (size_t v : knn_lists[u]) {
                        if (v > u) nbrs.push_back(v);
                    }
                    for (size_t v = u + 1; v < N; ++v) {
                        for (size_t w : knn_lists[v]) {
                            if (w == u) { nbrs.push_back(v); break; }
                        }
                    }
                    std::sort(nbrs.begin(), nbrs.end());
                    nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
                    g.sparse_adj[u] = std::move(nbrs);
                }
            }
            return g;
        }
    };

    static std::pair<Scalar, Scalar> run_semiglobal_banded_dp_pair(
        const std::vector<int8_t>& s_a,
        const std::vector<int8_t>& s_b,
        const GeneticCode& gcode,
        std::vector<Scalar>& prev_M,
        std::vector<Scalar>& curr_M,
        std::vector<Scalar>& prev_D,
        std::vector<Scalar>& curr_D,
        std::vector<size_t>& prev_id_M,
        std::vector<size_t>& curr_id_M,
        std::vector<size_t>& prev_id_D,
        std::vector<size_t>& curr_id_D
    ) {
        const auto& s_short = (s_a.size() <= s_b.size()) ? s_a : s_b;
        const auto& s_long  = (s_a.size() <= s_b.size()) ? s_b : s_a;
        size_t Ls = s_short.size();
        size_t Ll = s_long.size();
        if (Ls == 0 || Ll == 0) return {0.0, 0.0};

        auto aa_sim_score = [&](int8_t c1, int8_t c2) -> Scalar {
            if (c1 < 0 || c2 < 0) return 0.0;
            if (c1 == c2) return 1.00;
            char a1 = gcode.sense_codon_aa[c1];
            char a2 = gcode.sense_codon_aa[c2];
            if (a1 == a2) return 0.90;
            int nd = gcode.diff_matrix[c1][c2].num_diffs;
            if (nd == 1) return 0.30;
            return 0.0;
        };

        constexpr int B0 = 24;
        constexpr Scalar kNegInf = -1e9;
        constexpr Scalar g_open = 1.10;
        constexpr Scalar g_ext  = 0.30;
        constexpr Scalar g_cap  = 3.80;
        int dL = static_cast<int>(Ll - Ls);
        if (prev_M.size() < Ll + 1) {
            prev_M.resize(Ll + 1);
            curr_M.resize(Ll + 1);
            prev_D.resize(Ll + 1);
            curr_D.resize(Ll + 1);
            prev_id_M.resize(Ll + 1);
            curr_id_M.resize(Ll + 1);
            prev_id_D.resize(Ll + 1);
            curr_id_D.resize(Ll + 1);
        }
        std::fill(prev_M.begin(), prev_M.begin() + static_cast<ptrdiff_t>(Ll + 1), kNegInf);
        std::fill(prev_D.begin(), prev_D.begin() + static_cast<ptrdiff_t>(Ll + 1), kNegInf);
        std::fill(prev_id_M.begin(), prev_id_M.begin() + static_cast<ptrdiff_t>(Ll + 1), 0);
        std::fill(prev_id_D.begin(), prev_id_D.begin() + static_cast<ptrdiff_t>(Ll + 1), 0);
        std::fill(curr_M.begin(), curr_M.begin() + static_cast<ptrdiff_t>(Ll + 1), kNegInf);
        std::fill(curr_D.begin(), curr_D.begin() + static_cast<ptrdiff_t>(Ll + 1), kNegInf);
        std::fill(curr_id_M.begin(), curr_id_M.begin() + static_cast<ptrdiff_t>(Ll + 1), 0);
        std::fill(curr_id_D.begin(), curr_id_D.begin() + static_cast<ptrdiff_t>(Ll + 1), 0);
        int init_max = std::min(static_cast<int>(Ll), dL + B0);
        for (int c = 0; c <= init_max; ++c) {
            prev_M[static_cast<size_t>(c)] = 0.0;
        }

        for (size_t r = 1; r <= Ls; ++r) {
            int c_min = std::max(1, static_cast<int>(r) - B0);
            int c_max = std::min(static_cast<int>(Ll), static_cast<int>(r) + dL + B0);
            int prev_c_max = std::min(static_cast<int>(Ll), static_cast<int>(r - 1) + dL + B0);
            int fill_lo = std::max(0, c_min - 2);
            int fill_hi = std::min(static_cast<int>(Ll), c_max + 2);
            curr_M[0] = kNegInf;
            curr_D[0] = kNegInf;
            curr_id_M[0] = 0;
            curr_id_D[0] = 0;
            std::fill(curr_M.begin() + fill_lo, curr_M.begin() + fill_hi + 1, kNegInf);
            std::fill(curr_D.begin() + fill_lo, curr_D.begin() + fill_hi + 1, kNegInf);
            std::fill(curr_id_M.begin() + fill_lo, curr_id_M.begin() + fill_hi + 1, 0);
            std::fill(curr_id_D.begin() + fill_lo, curr_id_D.begin() + fill_hi + 1, 0);
            if (static_cast<int>(r) <= B0) {
                Scalar pen = -std::min(g_open + g_ext * static_cast<Scalar>(r), g_cap);
                curr_D[0] = pen;
                curr_M[0] = pen;
            }
            int8_t ci = s_short[r - 1];
            Scalar ins_sc = kNegInf;
            size_t ins_id = 0;
            Scalar best_row_sc = (c_min == 1) ? curr_M[0] : kNegInf;
            size_t best_row_id = 0;

            for (int c = c_min; c <= c_max; ++c) {
                size_t cu = static_cast<size_t>(c);
                int8_t cj = s_long[cu - 1];

                if (c <= prev_c_max) {
                    Scalar d_op = prev_M[cu] - g_open - g_ext;
                    Scalar d_ex = prev_D[cu] - g_ext;
                    if (d_op >= d_ex) {
                        curr_D[cu] = d_op;
                        curr_id_D[cu] = prev_id_M[cu];
                    } else {
                        curr_D[cu] = d_ex;
                        curr_id_D[cu] = prev_id_D[cu];
                    }
                } else {
                    curr_D[cu] = kNegInf;
                    curr_id_D[cu] = 0;
                }

                bool left_in_band = (c > c_min) || (c == 1 && static_cast<int>(r) <= B0);
                Scalar i_op = left_in_band ? (curr_M[cu - 1] - g_open - g_ext) : kNegInf;
                Scalar i_ex = (c > c_min) ? (ins_sc - g_ext) : kNegInf;
                if (i_op > i_ex) {
                    ins_sc = i_op;
                    ins_id = curr_id_M[cu - 1];
                } else {
                    ins_sc = i_ex;
                }
                Scalar e_sc = best_row_sc - g_cap;

                Scalar m_sc = prev_M[cu - 1] + aa_sim_score(ci, cj);
                bool exact_aa = (ci >= 0 && cj >= 0 && gcode.sense_codon_aa[ci] == gcode.sense_codon_aa[cj]);
                size_t m_id = prev_id_M[cu - 1] + (exact_aa ? 1 : 0);

                Scalar best_sc = m_sc;
                size_t best_id = m_id;
                if (curr_D[cu] > best_sc) {
                    best_sc = curr_D[cu];
                    best_id = curr_id_D[cu];
                }
                if (ins_sc > best_sc) {
                    best_sc = ins_sc;
                    best_id = ins_id;
                }
                if (e_sc > best_sc) {
                    best_sc = e_sc;
                    best_id = best_row_id;
                }
                curr_M[cu] = best_sc;
                curr_id_M[cu] = best_id;
                if (best_sc > best_row_sc) {
                    best_row_sc = best_sc;
                    best_row_id = best_id;
                }
            }
            std::swap(prev_M, curr_M);
            std::swap(prev_D, curr_D);
            std::swap(prev_id_M, curr_id_M);
            std::swap(prev_id_D, curr_id_D);
        }

        int c_start = std::max(1, static_cast<int>(Ls) - B0);
        size_t best_c = static_cast<size_t>(c_start);
        Scalar best_score = prev_M[best_c];
        for (size_t c = static_cast<size_t>(c_start) + 1; c <= Ll; ++c) {
            if (prev_M[c] > best_score) {
                best_score = prev_M[c];
                best_c = c;
            }
        }
        Scalar L_eff = std::max<Scalar>(1.0, static_cast<Scalar>(Ls));
        return {best_score / L_eff, static_cast<Scalar>(prev_id_M[best_c]) / L_eff};
    }

    struct PairwiseCodonTracks {
        Scalar dp_similarity = 0.0;
        Scalar aa_identity = 0.0;
        Scalar d_tot = 0.0;
        Scalar d_syn = 0.0;
        Scalar d_ns = 0.0;
        Scalar d_tn93 = 0.0;
    };

    static PairwiseCodonTracks run_semiglobal_multitrack_dp_pair(
        const std::vector<int8_t>& s_a,
        const std::vector<int8_t>& s_b,
        const std::string& raw_a,
        const std::string& raw_b,
        const GeneticCode& gcode,
        const Vector4& pi
    ) {
        bool swapped = (s_a.size() > s_b.size());
        const auto& s_short = swapped ? s_b : s_a;
        const auto& s_long  = swapped ? s_a : s_b;
        const std::string& r_short = swapped ? raw_b : raw_a;
        const std::string& r_long  = swapped ? raw_a : raw_b;
        size_t Ls = s_short.size();
        size_t Ll = s_long.size();
        if (Ls == 0 || Ll == 0) return {};

        auto aa_sim_score = [&](int8_t c1, int8_t c2) -> Scalar {
            if (c1 < 0 || c2 < 0) return 0.0;
            if (c1 == c2) return 1.00;
            char a1 = gcode.sense_codon_aa[c1];
            char a2 = gcode.sense_codon_aa[c2];
            if (a1 == a2) return 0.90;
            int nd = gcode.diff_matrix[c1][c2].num_diffs;
            if (nd == 1) return 0.30;
            return 0.0;
        };

        constexpr int B0 = 24;
        constexpr Scalar kNegInf = -1e9;
        constexpr Scalar g_open = 1.10;
        constexpr Scalar g_ext  = 0.30;
        constexpr Scalar g_cap  = 3.80;
        int dL = static_cast<int>(Ll - Ls);

        struct TrackAcc {
            size_t exact_aa = 0;
            Scalar syn_sc = 0.0;
            Scalar ns_sc = 0.0;
            Scalar ag_ts = 0.0;
            Scalar ct_ts = 0.0;
            Scalar tv = 0.0;
            Scalar valid_nt = 0.0;
        };

        std::vector<Scalar> prev_M(Ll + 1, kNegInf), curr_M(Ll + 1, kNegInf);
        std::vector<Scalar> prev_D(Ll + 1, kNegInf), curr_D(Ll + 1, kNegInf);
        std::vector<TrackAcc> prev_tr_M(Ll + 1), curr_tr_M(Ll + 1);
        std::vector<TrackAcc> prev_tr_D(Ll + 1), curr_tr_D(Ll + 1);

        int init_max = std::min(static_cast<int>(Ll), dL + B0);
        for (int c = 0; c <= init_max; ++c) {
            prev_M[static_cast<size_t>(c)] = 0.0;
        }

        for (size_t r = 1; r <= Ls; ++r) {
            int c_min = std::max(1, static_cast<int>(r) - B0);
            int c_max = std::min(static_cast<int>(Ll), static_cast<int>(r) + dL + B0);
            int prev_c_max = std::min(static_cast<int>(Ll), static_cast<int>(r - 1) + dL + B0);
            int fill_lo = std::max(0, c_min - 2);
            int fill_hi = std::min(static_cast<int>(Ll), c_max + 2);
            curr_M[0] = kNegInf;
            curr_D[0] = kNegInf;
            curr_tr_M[0] = TrackAcc{};
            curr_tr_D[0] = TrackAcc{};
            std::fill(curr_M.begin() + fill_lo, curr_M.begin() + fill_hi + 1, kNegInf);
            std::fill(curr_D.begin() + fill_lo, curr_D.begin() + fill_hi + 1, kNegInf);
            std::fill(curr_tr_M.begin() + fill_lo, curr_tr_M.begin() + fill_hi + 1, TrackAcc{});
            std::fill(curr_tr_D.begin() + fill_lo, curr_tr_D.begin() + fill_hi + 1, TrackAcc{});
            if (static_cast<int>(r) <= B0) {
                Scalar pen = -std::min(g_open + g_ext * static_cast<Scalar>(r), g_cap);
                curr_D[0] = pen;
                curr_M[0] = pen;
            }
            int8_t ci = s_short[r - 1];
            Scalar ins_sc = kNegInf;
            TrackAcc ins_tr{};
            Scalar best_row_sc = (c_min == 1) ? curr_M[0] : kNegInf;
            TrackAcc best_row_tr{};

            for (int c = c_min; c <= c_max; ++c) {
                size_t cu = static_cast<size_t>(c);
                int8_t cj = s_long[cu - 1];

                if (c <= prev_c_max) {
                    Scalar d_op = prev_M[cu] - g_open - g_ext;
                    Scalar d_ex = prev_D[cu] - g_ext;
                    if (d_op >= d_ex) {
                        curr_D[cu] = d_op;
                        curr_tr_D[cu] = prev_tr_M[cu];
                    } else {
                        curr_D[cu] = d_ex;
                        curr_tr_D[cu] = prev_tr_D[cu];
                    }
                } else {
                    curr_D[cu] = kNegInf;
                    curr_tr_D[cu] = TrackAcc{};
                }

                bool left_in_band = (c > c_min) || (c == 1 && static_cast<int>(r) <= B0);
                Scalar i_op = left_in_band ? (curr_M[cu - 1] - g_open - g_ext) : kNegInf;
                Scalar i_ex = (c > c_min) ? (ins_sc - g_ext) : kNegInf;
                if (i_op > i_ex) {
                    ins_sc = i_op;
                    ins_tr = curr_tr_M[cu - 1];
                } else {
                    ins_sc = i_ex;
                }
                Scalar e_sc = best_row_sc - g_cap;

                Scalar m_sc = prev_M[cu - 1] + aa_sim_score(ci, cj);
                TrackAcc m_tr = prev_tr_M[cu - 1];
                if (ci >= 0 && cj >= 0) {
                    int nd = gcode.diff_matrix[ci][cj].num_diffs;
                    if (ci == cj) {
                        m_tr.exact_aa += 1;
                        m_tr.syn_sc += 1.0;
                        m_tr.ns_sc  += 1.0;
                    } else if (gcode.sense_codon_aa[ci] == gcode.sense_codon_aa[cj]) {
                        m_tr.exact_aa += 1;
                        m_tr.syn_sc += (nd == 1) ? 0.65 : 0.35;
                        m_tr.ns_sc  += 1.0;
                    } else {
                        m_tr.syn_sc += (nd == 1) ? 0.45 : 0.15;
                        m_tr.ns_sc  += (nd == 1) ? 0.45 : 0.10;
                    }
                }
                if (3 * (r - 1) + 2 < r_short.size() && 3 * (cu - 1) + 2 < r_long.size()) {
                    for (int p = 0; p < 3; ++p) {
                        int b1 = char_to_nuc(r_short[3 * (r - 1) + p]);
                        int b2 = char_to_nuc(r_long[3 * (cu - 1) + p]);
                        if (b1 >= 0 && b1 < 4 && b2 >= 0 && b2 < 4) {
                            m_tr.valid_nt += 1.0;
                            if (b1 != b2) {
                                bool p1 = (b1 == 0 || b1 == 2);
                                bool p2 = (b2 == 0 || b2 == 2);
                                if (p1 && p2) m_tr.ag_ts += 1.0;
                                else if (!p1 && !p2) m_tr.ct_ts += 1.0;
                                else m_tr.tv += 1.0;
                            }
                        }
                    }
                }

                Scalar best_sc = m_sc;
                TrackAcc best_tr = m_tr;
                if (curr_D[cu] > best_sc) {
                    best_sc = curr_D[cu];
                    best_tr = curr_tr_D[cu];
                }
                if (ins_sc > best_sc) {
                    best_sc = ins_sc;
                    best_tr = ins_tr;
                }
                if (e_sc > best_sc) {
                    best_sc = e_sc;
                    best_tr = best_row_tr;
                }
                curr_M[cu] = best_sc;
                curr_tr_M[cu] = best_tr;
                if (best_sc > best_row_sc) {
                    best_row_sc = best_sc;
                    best_row_tr = best_tr;
                }
            }
            std::swap(prev_M, curr_M);
            std::swap(prev_D, curr_D);
            std::swap(prev_tr_M, curr_tr_M);
            std::swap(prev_tr_D, curr_tr_D);
        }

        int c_start = std::max(1, static_cast<int>(Ls) - B0);
        size_t best_c = static_cast<size_t>(c_start);
        Scalar best_score = prev_M[best_c];
        for (size_t c = static_cast<size_t>(c_start) + 1; c <= Ll; ++c) {
            if (prev_M[c] > best_score) {
                best_score = prev_M[c];
                best_c = c;
            }
        }
        Scalar L_eff = std::max<Scalar>(1.0, static_cast<Scalar>(Ls));
        const auto& tr = prev_tr_M[best_c];
        PairwiseCodonTracks res;
        res.dp_similarity = best_score / L_eff;
        res.aa_identity   = static_cast<Scalar>(tr.exact_aa) / L_eff;
        Scalar pid_tot = std::clamp<Scalar>(res.dp_similarity, 0.08, 0.9999);
        Scalar pid_syn = std::clamp<Scalar>(tr.syn_sc / L_eff, 0.08, 0.9999);
        Scalar pid_ns  = std::clamp<Scalar>(tr.ns_sc / L_eff, 0.08, 0.9999);
        res.d_tot = std::max<Scalar>(1e-4, -0.75 * std::log(pid_tot));
        res.d_syn = std::max<Scalar>(1e-4, -0.75 * std::log(pid_syn));
        res.d_ns  = std::max<Scalar>(1e-4, -0.75 * std::log(pid_ns));

        if (tr.valid_nt >= 15.0) {
            Scalar gR = pi(0) + pi(2);
            Scalar gY = pi(1) + pi(3);
            Scalar k1_denom = 2.0 * pi(0) * pi(2) / gR;
            Scalar k2_denom = 2.0 * pi(1) * pi(3) / gY;
            Scalar k3_denom = 2.0 * gR * gY;
            Scalar P1 = tr.ag_ts / tr.valid_nt;
            Scalar P2 = tr.ct_ts / tr.valid_nt;
            Scalar Q  = tr.tv / tr.valid_nt;
            Scalar arg1 = std::max<Scalar>(0.02, 1.0 - P1 / k1_denom - Q / (2.0 * gR));
            Scalar arg2 = std::max<Scalar>(0.02, 1.0 - P2 / k2_denom - Q / (2.0 * gY));
            Scalar arg3 = std::max<Scalar>(0.02, 1.0 - Q / k3_denom);
            Scalar term1 = k1_denom * std::log(arg1);
            Scalar term2 = k2_denom * std::log(arg2);
            Scalar term3 = 2.0 * (gR * gY - (pi(0) * pi(2) * gY / gR) - (pi(1) * pi(3) * gR / gY)) * std::log(arg3);
            res.d_tn93 = std::max<Scalar>(1e-4, -(term1 + term2 + term3));
        } else {
            res.d_tn93 = res.d_tot;
        }
        return res;
    }

    static std::pair<Alignment, Tree> build_de_novo_codon_alignment_and_tree(
        const std::vector<Sequence>& unaligned_seqs,
        std::shared_ptr<const GeneticCode> gcode = nullptr,
        const JointAlignSettings& settings = JointAlignSettings{},
        const Matrix* precomputed_dp_similarity = nullptr
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
        Vector4 pi_nuc = pos_freqs.colwise().mean().transpose();
        for (int b = 0; b < 4; ++b) {
            if (pi_nuc(b) < 1e-4) pi_nuc(b) = 1e-4;
        }
        pi_nuc /= pi_nuc.sum();

        // 2. Phase A: Multi-Track Alignment-Free DP, Ohmic Conductance Graph & Grounded Green's Spectral Scaffold
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

        (void)precomputed_dp_similarity;
        const bool use_subquadratic_taxa = (N > settings.subquadratic_taxon_threshold);
        LandmarkKmerGraph lm_graph;
        if (use_subquadratic_taxa) {
            lm_graph = LandmarkKmerGraph::build(taxon_codons, *gcode, 8);
        }

        Matrix D_tot  = Matrix::Zero(N, N);
        Matrix D_syn  = Matrix::Zero(N, N);
        Matrix D_ns   = Matrix::Zero(N, N);
        Matrix D_tn93 = Matrix::Zero(N, N);
        Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> pair_eval =
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(N, N);

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < N; ++i) {
            const auto& seq_i = taxon_codons[i];
            size_t Li = seq_i.size();
            if (use_subquadratic_taxa) {
                for (size_t j : lm_graph.sparse_adj[i]) {
                    const auto& seq_j = taxon_codons[j];
                    size_t Lj = seq_j.size();
                    if (Li == 0 || Lj == 0) continue;
                    auto tr = run_semiglobal_multitrack_dp_pair(
                        seq_i, seq_j, unaligned_seqs[i].raw_seq, unaligned_seqs[j].raw_seq, *gcode, pi_nuc
                    );
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    D_tot(ii, jj)  = D_tot(jj, ii)  = tr.d_tot;
                    D_syn(ii, jj)  = D_syn(jj, ii)  = tr.d_syn;
                    D_ns(ii, jj)   = D_ns(jj, ii)   = tr.d_ns;
                    D_tn93(ii, jj) = D_tn93(jj, ii) = tr.d_tn93;
                    pair_eval(ii, jj) = pair_eval(jj, ii) = 1;
                }
            } else {
                for (size_t j = i + 1; j < N; ++j) {
                    const auto& seq_j = taxon_codons[j];
                    size_t Lj = seq_j.size();
                    if (Li == 0 || Lj == 0) continue;

                    auto tr = run_semiglobal_multitrack_dp_pair(
                        seq_i, seq_j, unaligned_seqs[i].raw_seq, unaligned_seqs[j].raw_seq, *gcode, pi_nuc
                    );
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    D_tot(ii, jj)  = D_tot(jj, ii)  = tr.d_tot;
                    D_syn(ii, jj)  = D_syn(jj, ii)  = tr.d_syn;
                    D_ns(ii, jj)   = D_ns(jj, ii)   = tr.d_ns;
                    D_tn93(ii, jj) = D_tn93(jj, ii) = tr.d_tn93;
                    pair_eval(ii, jj) = pair_eval(jj, ii) = 1;
                }
            }
        }

        // Depth-adaptive multi-track fusion:
        // Shallow splits are resolved by the synonymous-weighted track D_syn and transition/transversion D_tn93,
        // while deep splits (> 0.35 substitutions/site) blend in the non-synonymous track D_ns to prevent saturation rollover.
        Scalar sum_tn = 0.0, sum_syn = 0.0, sum_ns = 0.0;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Eigen::Index ii = static_cast<Eigen::Index>(i);
                Eigen::Index jj = static_cast<Eigen::Index>(j);
                if (!pair_eval(ii, jj)) continue;
                if (D_tn93(ii, jj) < 0.35) {
                    sum_tn  += D_tn93(ii, jj);
                    sum_syn += D_syn(ii, jj);
                    sum_ns  += D_ns(ii, jj);
                }
            }
        }
        Scalar scale_syn = (sum_syn > 1e-9 && sum_tn > 1e-9) ? (sum_tn / sum_syn) : 1.0;
        Scalar scale_ns  = (sum_ns > 1e-9 && sum_tn > 1e-9) ? (sum_tn / sum_ns) : 1.0;

        Matrix D_af = Matrix::Zero(N, N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i + 1; j < N; ++j) {
                Eigen::Index ii = static_cast<Eigen::Index>(i);
                Eigen::Index jj = static_cast<Eigen::Index>(j);
                if (!pair_eval(ii, jj)) continue;
                Scalar d0 = D_tn93(ii, jj);
                Scalar ds = D_syn(ii, jj) * scale_syn;
                Scalar dn = D_ns(ii, jj) * scale_ns;
                Scalar alpha_s = 0.20 * std::exp(-(d0 / 0.35) * (d0 / 0.35));
                Scalar d = alpha_s * ds + (1.0 - alpha_s) * (0.85 * d0 + 0.15 * dn);
                D_af(ii, jj) = D_af(jj, ii) = d;
            }
        }

        // Subquadratic Four-Point Landmark Quartet Triangulation for unevaluated pairs (i, j) not in E_sparse:
        // Given exact distances from all N taxa to the K landmarks P, for each taxon u we find its two closest
        // landmarks p1(u), p2(u) in P and reconstruct D_af(i, j) = D_af(i, q) + D_af(j, p) - D_af(p, q)
        // clamped to the exact metric triangle bounds across the anchor landmarks.
        std::vector<size_t> near_lm1(N, 0), near_lm2(N, 0);
        if (use_subquadratic_taxa && lm_graph.num_landmarks >= 2) {
            for (size_t u = 0; u < N; ++u) {
                size_t best1 = lm_graph.landmarks[0];
                size_t best2 = lm_graph.landmarks[1];
                Scalar d1 = D_af(static_cast<Eigen::Index>(u), static_cast<Eigen::Index>(best1));
                Scalar d2 = D_af(static_cast<Eigen::Index>(u), static_cast<Eigen::Index>(best2));
                if (d2 < d1) { std::swap(d1, d2); std::swap(best1, best2); }
                for (size_t k = 2; k < lm_graph.num_landmarks; ++k) {
                    size_t lm = lm_graph.landmarks[k];
                    Scalar d = D_af(static_cast<Eigen::Index>(u), static_cast<Eigen::Index>(lm));
                    if (d < d1) {
                        d2 = d1; best2 = best1;
                        d1 = d;  best1 = lm;
                    } else if (d < d2) {
                        d2 = d;  best2 = lm;
                    }
                }
                near_lm1[u] = best1;
                near_lm2[u] = best2;
            }

            for (size_t i = 0; i < N; ++i) {
                Eigen::Index ii = static_cast<Eigen::Index>(i);
                size_t p1 = near_lm1[i], p2 = near_lm2[i];
                Eigen::Index ep1 = static_cast<Eigen::Index>(p1);
                Eigen::Index ep2 = static_cast<Eigen::Index>(p2);
                for (size_t j = i + 1; j < N; ++j) {
                    Eigen::Index jj = static_cast<Eigen::Index>(j);
                    if (pair_eval(ii, jj)) continue;
                    size_t q1 = near_lm1[j], q2 = near_lm2[j];
                    size_t p = p1;
                    size_t q = (q1 != p) ? q1 : q2;
                    Eigen::Index ep = static_cast<Eigen::Index>(p);
                    Eigen::Index eq = static_cast<Eigen::Index>(q);
                    Eigen::Index eq1 = static_cast<Eigen::Index>(q1);
                    Eigen::Index eq2 = static_cast<Eigen::Index>(q2);

                    Scalar d_4pt = D_af(ii, eq) + D_af(jj, ep) - D_af(ep, eq);
                    std::array<Eigen::Index, 4> anchors = {ep1, ep2, eq1, eq2};
                    Scalar d_min = 1e-4;
                    Scalar d_max = 1e9;
                    for (Eigen::Index er : anchors) {
                        Scalar dip = D_af(ii, er);
                        Scalar djp = D_af(jj, er);
                        d_min = std::max(d_min, std::abs(dip - djp));
                        d_max = std::min(d_max, dip + djp);
                    }
                    if (d_max < d_min) d_max = d_min;
                    Scalar d_est = std::clamp(d_4pt, d_min, d_max);
                    D_af(ii, jj) = D_af(jj, ii) = d_est;
                }
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

        Matrix Omega_af = D_af;
        if (!use_subquadratic_taxa) {
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
        } else {
            // Subquadratic O(K^3 + N^2) Landmark Circuit Nystrom Resistance Extension:
            // Solve exact Moore-Penrose pseudoinverse on the K x K landmark sub-circuit L_lm,
            // then interpolate non-landmark port resistances via their two nearest landmark ports!
            size_t K = lm_graph.num_landmarks;
            Matrix W_lm = Matrix::Zero(K, K);
            for (size_t a = 0; a < K; ++a) {
                Eigen::Index ia = static_cast<Eigen::Index>(lm_graph.landmarks[a]);
                for (size_t b = a + 1; b < K; ++b) {
                    Eigen::Index ib = static_cast<Eigen::Index>(lm_graph.landmarks[b]);
                    Scalar d = D_af(ia, ib);
                    Scalar w = std::exp(-d / sigma_af) / (d + 1e-3);
                    W_lm(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) = w;
                    W_lm(static_cast<Eigen::Index>(b), static_cast<Eigen::Index>(a)) = w;
                }
            }
            Vector deg_lm = W_lm.rowwise().sum();
            Matrix L_lm = Matrix(deg_lm.asDiagonal()) - W_lm;
            L_lm = 0.5 * (L_lm + L_lm.transpose());
            Matrix Omega_lm = Matrix::Zero(K, K);
            Eigen::SelfAdjointEigenSolver<Matrix> eig_lm(L_lm);
            if (eig_lm.info() == Eigen::Success) {
                const Vector& ev = eig_lm.eigenvalues();
                const Matrix& evec = eig_lm.eigenvectors();
                Matrix L_pinv = Matrix::Zero(K, K);
                for (Eigen::Index m = 0; m < ev.size(); ++m) {
                    if (ev(m) > 1e-10) {
                        L_pinv.noalias() += (1.0 / ev(m)) * (evec.col(m) * evec.col(m).transpose());
                    }
                }
                for (size_t a = 0; a < K; ++a) {
                    for (size_t b = a + 1; b < K; ++b) {
                        Eigen::Index ea = static_cast<Eigen::Index>(a);
                        Eigen::Index eb = static_cast<Eigen::Index>(b);
                        Scalar om = std::max<Scalar>(1e-4, L_pinv(ea, ea) + L_pinv(eb, eb) - 2.0 * L_pinv(ea, eb));
                        Omega_lm(ea, eb) = Omega_lm(eb, ea) = om;
                    }
                }
                std::vector<size_t> taxon_to_lm_idx(N, 0);
                for (size_t a = 0; a < K; ++a) taxon_to_lm_idx[lm_graph.landmarks[a]] = a;
                std::vector<Scalar> attach_r(N, 0.0);
                for (size_t u = 0; u < N; ++u) {
                    if (lm_graph.is_landmark[u]) continue;
                    Eigen::Index eu = static_cast<Eigen::Index>(u);
                    Eigen::Index ep1 = static_cast<Eigen::Index>(near_lm1[u]);
                    Eigen::Index ep2 = static_cast<Eigen::Index>(near_lm2[u]);
                    Scalar d_stem = std::max<Scalar>(0.0, 0.5 * (D_af(eu, ep1) + D_af(eu, ep2) - D_af(ep1, ep2)));
                    attach_r[u] = d_stem;
                }
                Scalar med_d_lm = std::max<Scalar>(1e-4, sigma_af);
                Scalar med_om_lm = std::max<Scalar>(1e-4, Omega_lm.mean());
                Scalar r_scale = med_om_lm / med_d_lm;
                for (size_t i = 0; i < N; ++i) {
                    Eigen::Index ii = static_cast<Eigen::Index>(i);
                    Eigen::Index la = static_cast<Eigen::Index>(taxon_to_lm_idx[near_lm1[i]]);
                    for (size_t j = i + 1; j < N; ++j) {
                        Eigen::Index jj = static_cast<Eigen::Index>(j);
                        Eigen::Index lb = static_cast<Eigen::Index>(taxon_to_lm_idx[near_lm1[j]]);
                        Scalar om = (la == lb)
                            ? (r_scale * D_af(ii, jj))
                            : (Omega_lm(la, lb) + r_scale * (attach_r[i] + attach_r[j]));
                        Omega_af(ii, jj) = Omega_af(jj, ii) = std::max<Scalar>(1e-4, om);
                    }
                }
            }
        }

        // Fuse direct evolutionary distance D_af with Ohmic circuit resistance Omega_af
        // (normalized to the median evolutionary scale) and construct the Grounded Green's Spectral Scaffold!
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

        // Build balanced Grounded Dirichlet Green's Function Spectral Guide Tree
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

        // 3. Calibrate Stage-0 Empirical Codon Substitution Generator Q_emp (Core + Loop 2-rate mixture):
        EmpiricalTwoStratumModel two_stratum = EmpiricalTwoStratumModel::calibrate(
            taxon_codons, pos_freqs, pi_codon, *gcode, settings.subquadratic_taxon_threshold
        );

        // Empirical Two-Stratum Core + Loop transition probability matrix at true branch length t_b (t_0 = 0):
        //   P_{s,x}(t_b) = (1 - \hat{w}_loop) \exp(t_b \hat{Q}_core) + \hat{w}_loop \exp(\hat{\rho}_{loop/core} t_b \hat{Q}_loop)
        auto branch_transition_matrix = [&](Scalar bl) -> Matrix {
            return two_stratum.transition_matrix(bl);
        };

        // 4. Phase B: Bottom-Up Post-Order Phylogenetic Codon Transducer (Forward Felsenstein Pruning)
        //    Each profile column c at node v stores its standard Felsenstein partial conditional likelihood
        //    vector \ell_c^{(v)} \in R_+^S (\ell_{c,x}^{(v)} \propto P(Data_{v,c} | X_v = x)), normalized so
        //    that its equilibrium expectation is unity (\sum_{x=1}^S \pi_x \ell_{c,x}^{(v)} = 1).
        struct NodeProfile {
            std::vector<Vector> L;                      // Felsenstein partial likelihoods \ell_c^{(v)}, size L_v x S
            std::vector<Scalar> occupancy;              // clade column occupancy o_c^{(v)} in (0, 1]
            std::vector<std::vector<int32_t>> tax_pos;  // size L_v x N (-1 if gap in taxon t)
            size_t num_leaves = 0;
        };

        std::vector<NodeProfile> profiles(guide_tree.num_nodes());
        std::unordered_map<std::string, size_t> name_to_idx;
        for (size_t t = 0; t < N; ++t) name_to_idx[taxon_names[t]] = t;

        // Forward Felsenstein Homology Bayes Factor at parent u:
        // Given forward conditional likelihood vectors L_i^{(L)} = P(t_L) \ell_i^{(L)} and L_j^{(R)} = P(t_R) \ell_j^{(R)}
        // (each representing the probability that an ancestral codon s at parent u evolves forward down branch
        // t_L or t_R to produce the observed leaf codons in clade L or R, normalized so \pi^T L_i^{(L)} = 1),
        // the marginal likelihood ratio of homologous origin (single ancestor s ~ \pi at u) vs. independent
        // origins (s_L ~ \pi, s_R ~ \pi) is simply \sum_s \pi_s L_{i,s}^{(L)} L_{j,s}^{(R)}.
        const Scalar eps_floor = settings.transducer_laplace_floor;
        auto homology_bayes_factor = [&](const Vector& L_left, const Vector& L_right) -> Scalar {
            Scalar bf = 0.0;
            for (int s = 0; s < S; ++s) {
                bf += pi_codon(s) * L_left(s) * L_right(s);
            }
            return bf;
        };

        auto align_two_profiles = [&](const NodeProfile& left, Scalar bl_left,
                                      const NodeProfile& right, Scalar bl_right) -> NodeProfile {
            size_t L1 = left.L.size();
            size_t L2 = right.L.size();
            if (L1 == 0) return right;
            if (L2 == 0) return left;

            Matrix P1 = branch_transition_matrix(bl_left);
            Matrix P2 = branch_transition_matrix(bl_right);

            // Forward ancestor-to-child Felsenstein conditional likelihood vectors at parent u:
            //   L1_fwd[i](s) = \sum_{x=1}^S P_{s,x}(t_L) \ell_{i,x}^{(L)}
            //   L2_fwd[j](s) = \sum_{y=1}^S P_{s,y}(t_R) \ell_{j,y}^{(R)}
            std::vector<Vector> L1_fwd(L1), L2_fwd(L2);
            for (size_t i = 0; i < L1; ++i) {
                L1_fwd[i] = P1 * left.L[i];
                Scalar z1 = pi_codon.dot(L1_fwd[i]);
                if (z1 > 0.0) L1_fwd[i] /= z1;
            }
            for (size_t j = 0; j < L2; ++j) {
                L2_fwd[j] = P2 * right.L[j];
                Scalar z2 = pi_codon.dot(L2_fwd[j]);
                if (z2 > 0.0) L2_fwd[j] /= z2;
            }

            // 3-State Affine Pair-HMM Codon Transducer (Gotoh's algorithm):
            //   State 0 (M): Homologous match (i, j)
            //   State 1 (I): Left-clade insertion / right-clade deletion (i, -)
            //   State 2 (D): Right-clade insertion / left-clade deletion (-, j)
            const Scalar NEG_INF = -1e18;
            const Scalar g_open = settings.transducer_gap_open;
            const Scalar g_ext  = settings.transducer_gap_extend;

            // Subquadratic O(B * L) Collinear 3-AA Seed-Chain Anchor Corridor when min(L1, L2) > subquadratic_length_threshold
            const bool use_corridor = (std::min(L1, L2) > settings.subquadratic_length_threshold);
            std::vector<size_t> row_lo(L1 + 1, 0), row_hi(L1 + 1, L2);
            if (use_corridor) {
                auto modal_aa = [&](const std::vector<Vector>& L_prof) -> std::vector<int> {
                    std::vector<int> res(L_prof.size(), -1);
                    for (size_t k = 0; k < L_prof.size(); ++k) {
                        int best_s = -1;
                        Scalar best_w = -1.0;
                        for (int s = 0; s < S; ++s) {
                            Scalar w = pi_codon(s) * L_prof[k](s);
                            if (w > best_w) { best_w = w; best_s = s; }
                        }
                        if (best_s >= 0) {
                            res[k] = LandmarkKmerGraph::aa_to_index20(gcode->sense_codon_aa[best_s]);
                        }
                    }
                    return res;
                };
                auto aa1 = modal_aa(L1_fwd);
                auto aa2 = modal_aa(L2_fwd);

                // Extract low-multiplicity (<= 3 occurrences) 3-AA seed matches
                std::unordered_map<uint32_t, std::vector<size_t>> pos1, pos2;
                for (size_t i = 0; i + 2 < L1; ++i) {
                    if (aa1[i] >= 0 && aa1[i + 1] >= 0 && aa1[i + 2] >= 0) {
                        uint32_t kmer = static_cast<uint32_t>(aa1[i] * 400 + aa1[i + 1] * 20 + aa1[i + 2]);
                        pos1[kmer].push_back(i + 1);
                    }
                }
                for (size_t j = 0; j + 2 < L2; ++j) {
                    if (aa2[j] >= 0 && aa2[j + 1] >= 0 && aa2[j + 2] >= 0) {
                        uint32_t kmer = static_cast<uint32_t>(aa2[j] * 400 + aa2[j + 1] * 20 + aa2[j + 2]);
                        pos2[kmer].push_back(j + 1);
                    }
                }

                std::vector<std::pair<size_t, size_t>> seeds;
                for (const auto& [kmer, v1] : pos1) {
                    if (v1.size() > 3) continue;
                    auto it2 = pos2.find(kmer);
                    if (it2 == pos2.end() || it2->second.size() > 3) continue;
                    for (size_t r : v1) {
                        for (size_t c : it2->second) {
                            seeds.push_back({r, c});
                        }
                    }
                }
                std::sort(seeds.begin(), seeds.end(), [](const auto& a, const auto& b) {
                    if (a.first != b.first) return a.first < b.first;
                    return a.second > b.second; // reverse second for strict LIS
                });

                // Longest strictly increasing collinear chain in O(M log M)
                std::vector<std::pair<size_t, size_t>> anchors;
                anchors.push_back({0, 0});
                if (!seeds.empty()) {
                    std::vector<size_t> tail_val;
                    std::vector<size_t> tail_idx;
                    std::vector<ptrdiff_t> parent(seeds.size(), -1);
                    for (size_t m = 0; m < seeds.size(); ++m) {
                        size_t c = seeds[m].second;
                        auto it = std::lower_bound(tail_val.begin(), tail_val.end(), c);
                        size_t pos = static_cast<size_t>(std::distance(tail_val.begin(), it));
                        if (it == tail_val.end()) {
                            tail_val.push_back(c);
                            tail_idx.push_back(m);
                        } else {
                            *it = c;
                            tail_idx[pos] = m;
                        }
                        if (pos > 0) parent[m] = static_cast<ptrdiff_t>(tail_idx[pos - 1]);
                    }
                    std::vector<std::pair<size_t, size_t>> chain;
                    ptrdiff_t cur = static_cast<ptrdiff_t>(tail_idx.back());
                    while (cur >= 0) {
                        chain.push_back(seeds[static_cast<size_t>(cur)]);
                        cur = parent[static_cast<size_t>(cur)];
                    }
                    std::reverse(chain.begin(), chain.end());
                    for (const auto& pt : chain) anchors.push_back(pt);
                }
                anchors.push_back({L1, L2});

                int max_diag_excursion = 0;
                for (size_t k = 0; k + 1 < anchors.size(); ++k) {
                    int dr = static_cast<int>(anchors[k + 1].first) - static_cast<int>(anchors[k].first);
                    int dc = static_cast<int>(anchors[k + 1].second) - static_cast<int>(anchors[k].second);
                    max_diag_excursion = std::max(max_diag_excursion, std::abs(dc - dr));
                }
                int dL_abs = static_cast<int>(std::max(L1, L2) - std::min(L1, L2));
                int B_prof = std::max({64, dL_abs + 48, max_diag_excursion + 32});

                row_lo[0] = 0;
                row_hi[0] = std::min<size_t>(L2, static_cast<size_t>(B_prof));
                size_t seg = 0;
                for (size_t i = 1; i <= L1; ++i) {
                    while (seg + 2 < anchors.size() && anchors[seg + 1].first < i) {
                        seg++;
                    }
                    double r0 = static_cast<double>(anchors[seg].first);
                    double c0 = static_cast<double>(anchors[seg].second);
                    double r1 = static_cast<double>(anchors[seg + 1].first);
                    double c1 = static_cast<double>(anchors[seg + 1].second);
                    double frac = (r1 > r0) ? ((static_cast<double>(i) - r0) / (r1 - r0)) : 0.0;
                    int j_cen = static_cast<int>(std::llround(c0 + frac * (c1 - c0)));
                    row_lo[i] = (static_cast<int>(i) <= B_prof)
                              ? 0
                              : static_cast<size_t>(std::min<int>(static_cast<int>(L2), std::max(1, j_cen - B_prof)));
                    row_hi[i] = static_cast<size_t>(std::min<int>(static_cast<int>(L2), std::max(1, j_cen + B_prof)));
                }
                row_hi[L1] = L2;
                row_lo[L1] = std::min(row_lo[L1], L2);
                // Enforce monotonic corridor connectivity from (0, 0) to (L1, L2)
                for (size_t i = 1; i <= L1; ++i) {
                    row_lo[i] = std::max(row_lo[i], row_lo[i - 1]);
                    row_hi[i] = std::max(row_hi[i], row_hi[i - 1]);
                    if (row_lo[i] > row_hi[i - 1]) row_lo[i] = row_hi[i - 1];
                }
                for (size_t i = L1; i > 0; --i) {
                    if (row_hi[i - 1] < row_lo[i]) row_hi[i - 1] = row_lo[i];
                    if (row_lo[i - 1] > row_lo[i]) row_lo[i - 1] = row_lo[i];
                }
            }

            // Allocate compact O(B * L) corridor buffers instead of dense O(L1 * L2) matrices
            std::vector<size_t> row_offset(L1 + 2, 0);
            for (size_t i = 0; i <= L1; ++i) {
                row_offset[i + 1] = row_offset[i] + (row_hi[i] - row_lo[i] + 1);
            }
            const size_t total_cells = row_offset[L1 + 1];
            std::vector<Scalar> dp_M_buf(total_cells, NEG_INF);
            std::vector<Scalar> dp_I_buf(total_cells, NEG_INF);
            std::vector<Scalar> dp_D_buf(total_cells, NEG_INF);
            std::vector<uint8_t> bt_M_buf(total_cells, 0);
            std::vector<uint8_t> bt_I_buf(total_cells, 0);
            std::vector<uint8_t> bt_D_buf(total_cells, 0);

            auto in_band = [&](size_t r, size_t c) -> bool {
                return r <= L1 && c >= row_lo[r] && c <= row_hi[r];
            };
            auto cell_idx = [&](size_t r, size_t c) -> size_t {
                return row_offset[r] + (c - row_lo[r]);
            };
            auto get_M = [&](size_t r, size_t c) -> Scalar {
                return in_band(r, c) ? dp_M_buf[cell_idx(r, c)] : NEG_INF;
            };
            auto get_I = [&](size_t r, size_t c) -> Scalar {
                return in_band(r, c) ? dp_I_buf[cell_idx(r, c)] : NEG_INF;
            };
            auto get_D = [&](size_t r, size_t c) -> Scalar {
                return in_band(r, c) ? dp_D_buf[cell_idx(r, c)] : NEG_INF;
            };

            dp_M_buf[cell_idx(0, 0)] = 0.0;
            for (size_t i = 1; i <= L1; ++i) {
                if (!in_band(i, 0)) break;
                Scalar occ_L = left.occupancy[i - 1];
                size_t idx = cell_idx(i, 0);
                dp_I_buf[idx] = (i == 1) ? (-g_open * occ_L) : (get_I(i - 1, 0) - g_ext * occ_L);
                bt_I_buf[idx] = (i == 1) ? 0 : 1;
            }
            for (size_t j = 1; j <= row_hi[0]; ++j) {
                Scalar occ_R = right.occupancy[j - 1];
                size_t idx = cell_idx(0, j);
                dp_D_buf[idx] = (j == 1) ? (-g_open * occ_R) : (get_D(0, j - 1) - g_ext * occ_R);
                bt_D_buf[idx] = (j == 1) ? 0 : 2;
            }

            for (size_t i = 1; i <= L1; ++i) {
                size_t j_start = std::max<size_t>(1, row_lo[i]);
                size_t j_end   = row_hi[i];
                Scalar occ_L = std::max<Scalar>(0.50, left.occupancy[i - 1]);
                Scalar eff_open_I = g_open * occ_L;
                Scalar eff_ext_I  = g_ext  * occ_L;
                for (size_t j = j_start; j <= j_end; ++j) {
                    size_t idx = cell_idx(i, j);
                    // 1. State I: align left[i-1] to gap, weighted by clade occupancy o_i^(L) (floored at 0.50)
                    if (in_band(i - 1, j)) {
                        Scalar i_from_m = get_M(i - 1, j) - eff_open_I;
                        Scalar i_from_i = get_I(i - 1, j) - eff_ext_I;
                        if (i_from_m >= i_from_i) {
                            dp_I_buf[idx] = i_from_m;
                            bt_I_buf[idx] = 0;
                        } else {
                            dp_I_buf[idx] = i_from_i;
                            bt_I_buf[idx] = 1;
                        }
                    }

                    // 2. State D: align right[j-1] to gap, weighted by clade occupancy o_j^(R) (floored at 0.50)
                    if (in_band(i, j - 1)) {
                        Scalar occ_R = std::max<Scalar>(0.50, right.occupancy[j - 1]);
                        Scalar eff_open_D = g_open * occ_R;
                        Scalar eff_ext_D  = g_ext  * occ_R;
                        Scalar d_from_m = get_M(i, j - 1) - eff_open_D;
                        Scalar d_from_d = get_D(i, j - 1) - eff_ext_D;
                        if (d_from_m >= d_from_d) {
                            dp_D_buf[idx] = d_from_m;
                            bt_D_buf[idx] = 0;
                        } else {
                            dp_D_buf[idx] = d_from_d;
                            bt_D_buf[idx] = 2;
                        }
                    }

                    // 3. State M: homologous Match of left[i-1] and right[j-1]
                    if (in_band(i - 1, j - 1)) {
                        Scalar bf = homology_bayes_factor(L1_fwd[i - 1], L2_fwd[j - 1]);
                        Scalar emit = std::log(std::max<Scalar>(1e-4, bf + eps_floor));
                        Scalar m0 = get_M(i - 1, j - 1);
                        Scalar m1 = get_I(i - 1, j - 1);
                        Scalar m2 = get_D(i - 1, j - 1);
                        Scalar best_m = m0;
                        uint8_t best_bt = 0;
                        if (m1 > best_m) {
                            best_m = m1;
                            best_bt = 1;
                        }
                        if (m2 > best_m) {
                            best_m = m2;
                            best_bt = 2;
                        }
                        dp_M_buf[idx] = best_m + emit;
                        bt_M_buf[idx] = best_bt;
                    }
                }
            }

            // Backtrack optimal transducer path
            int state = 0;
            Scalar best_end = get_M(L1, L2);
            if (get_I(L1, L2) > best_end) { best_end = get_I(L1, L2); state = 1; }
            if (get_D(L1, L2) > best_end) { best_end = get_D(L1, L2); state = 2; }

            struct Step { int op; size_t i; size_t j; }; // op: 0=M, 1=I, 2=D
            std::vector<Step> path;
            size_t ci = L1, cj = L2;
            while (ci > 0 || cj > 0) {
                if (state == 0 && ci > 0 && cj > 0 && in_band(ci, cj)) {
                    path.push_back({0, ci - 1, cj - 1});
                    state = bt_M_buf[cell_idx(ci, cj)];
                    ci--; cj--;
                } else if (state == 1 && ci > 0 && in_band(ci, cj)) {
                    path.push_back({1, ci - 1, 0});
                    state = bt_I_buf[cell_idx(ci, cj)];
                    ci--;
                } else if (state == 2 && cj > 0 && in_band(ci, cj)) {
                    path.push_back({2, 0, cj - 1});
                    state = bt_D_buf[cell_idx(ci, cj)];
                    cj--;
                } else if (ci > 0 && (cj == 0 || !in_band(ci, cj - 1))) {
                    path.push_back({1, ci - 1, 0});
                    ci--;
                } else if (cj > 0) {
                    path.push_back({2, 0, cj - 1});
                    cj--;
                } else {
                    path.push_back({1, ci - 1, 0});
                    ci--;
                }
            }
            std::reverse(path.begin(), path.end());

            // Assemble merged parent profile using Felsenstein's pruning rule:
            //   - Match (M): \ell_{c,s}^{(u)} = L_{i,s}^{(L)} * L_{j,s}^{(R)} / BF(i, j)
            //   - Insertion (I, gap in R): \ell_{c,s}^{(u)} = L_{i,s}^{(L)}
            //   - Deletion (D, gap in L): \ell_{c,s}^{(u)} = L_{j,s}^{(R)}
            NodeProfile parent;
            parent.num_leaves = left.num_leaves + right.num_leaves;
            size_t L_new = path.size();
            parent.L.resize(L_new);
            parent.occupancy.resize(L_new);
            parent.tax_pos.assign(L_new, std::vector<int32_t>(N, -1));

            Scalar w_left = static_cast<Scalar>(left.num_leaves) / static_cast<Scalar>(parent.num_leaves);
            Scalar w_right = static_cast<Scalar>(right.num_leaves) / static_cast<Scalar>(parent.num_leaves);

            for (size_t col = 0; col < L_new; ++col) {
                const auto& st = path[col];
                if (st.op == 0) {
                    Vector prod = L1_fwd[st.i].cwiseProduct(L2_fwd[st.j]);
                    Scalar bf = pi_codon.dot(prod);
                    if (bf > 1e-12) {
                        parent.L[col] = prod / bf;
                    } else {
                        parent.L[col] = w_left * L1_fwd[st.i] + w_right * L2_fwd[st.j];
                    }
                    parent.occupancy[col] = w_left * left.occupancy[st.i] + w_right * right.occupancy[st.j];
                    for (size_t t = 0; t < N; ++t) {
                        parent.tax_pos[col][t] = (left.tax_pos[st.i][t] != -1) ? left.tax_pos[st.i][t] : right.tax_pos[st.j][t];
                    }
                } else if (st.op == 1) {
                    parent.L[col] = L1_fwd[st.i];
                    parent.occupancy[col] = w_left * left.occupancy[st.i];
                    parent.tax_pos[col] = left.tax_pos[st.i];
                } else {
                    parent.L[col] = L2_fwd[st.j];
                    parent.occupancy[col] = w_right * right.occupancy[st.j];
                    parent.tax_pos[col] = right.tax_pos[st.j];
                }
            }
            return parent;
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
                lp.L.resize(m);
                lp.occupancy.assign(m, 1.0);
                lp.tax_pos.assign(m, std::vector<int32_t>(N, -1));
                for (size_t k = 0; k < m; ++k) {
                    lp.tax_pos[k][tid] = static_cast<int32_t>(k);
                    int8_t sc = taxon_codons[tid][k];
                    if (sc >= 0 && sc < S) {
                        // Standard Felsenstein leaf conditional likelihood \ell_{k,x}^{(t)} = I[x == sc] / \pi_{sc}
                        // normalized by the marginal equilibrium probability \pi_{sc} so \sum_x \pi_x \ell_{k,x}^{(t)} = 1.
                        lp.L[k] = Vector::Zero(S);
                        lp.L[k](sc) = 1.0 / std::max<Scalar>(1e-8, pi_codon(sc));
                    } else {
                        lp.L[k] = Vector::Ones(S);
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
                    cur = align_two_profiles(cur, (ci == 1) ? bl0 : 1e-4, profiles[ck], guide_tree.nodes[ck].branch_length);
                }
                profiles[nid] = std::move(cur);
            }
        }

        // 5. Assemble Root Profile Alignment
        NodeProfile root_prof = std::move(profiles[guide_tree.root_id]);

        size_t C_total = root_prof.L.size();

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
        if (N > settings.subquadratic_taxon_threshold) {
            tes.perform_spr = false;
            tes.perform_topological_gradients = false;
            tes.lbfgs_max_iterations = 15;
        }
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
        {
            std::unordered_set<std::string> seen_names;
            for (const auto& sq : raw_seqs) {
                if (!seen_names.insert(sq.name).second) {
                    throw std::runtime_error(
                        "Gate 1 (Duplicate Taxon Identifier): Input contains duplicate sequence identifier '" +
                        sq.name + "'. Every sequence must have a unique taxon name."
                    );
                }
            }
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

        // Helper to extract in-frame triplets, strip terminal stop codon (and short 5'/3' UTR flanks), and count internal stop codons (E1)
        auto extract_frame_cds = [&](const std::string& nt_seq, int offset, size_t& out_internal_stops,
                                     std::vector<int8_t>& out_codons, bool* out_utr_trimmed = nullptr) -> std::string {
            out_internal_stops = 0;
            out_codons.clear();
            if (out_utr_trimmed) *out_utr_trimmed = false;
            std::string clean;
            if (static_cast<size_t>(offset) >= nt_seq.size()) return clean;
            size_t usable_len = ((nt_seq.size() - static_cast<size_t>(offset)) / 3) * 3;

            std::vector<std::string> trips;
            std::vector<int8_t> scs;
            std::vector<bool> is_stop;
            trips.reserve(usable_len / 3);
            scs.reserve(usable_len / 3);
            is_stop.reserve(usable_len / 3);

            for (size_t i = static_cast<size_t>(offset); i + 2 < static_cast<size_t>(offset) + usable_len; i += 3) {
                std::string trip = nt_seq.substr(i, 3);
                if (trip == "---" || trip == "???" || trip == "NNN") continue;
                for (char& ch : trip) {
                    if (ch == '-' || ch == '?') ch = 'N';
                }
                int8_t sc = gcode->triplet_to_sense(trip[0], trip[1], trip[2]);
                bool has_n = (trip.find('N') != std::string::npos);
                trips.push_back(std::move(trip));
                scs.push_back(sc);
                is_stop.push_back(sc < 0 && !has_n);
            }
            size_t L = trips.size();
            if (L == 0) return clean;

            // Check if a short 5' UTR leader (<= 20 codons) with an early stop precedes an ATG-initiated ORF
            size_t start_idx = 0;
            size_t first_stop = L;
            for (size_t k = 0; k < L; ++k) {
                if (is_stop[k]) { first_stop = k; break; }
            }
            if (first_stop < std::min<size_t>(20, L / 6)) {
                for (size_t k = first_stop + 1; k <= std::min<size_t>(25, L / 5); ++k) {
                    if (trips[k] == "ATG") {
                        size_t next_stop = L;
                        for (size_t j = k + 1; j < L; ++j) {
                            if (is_stop[j]) { next_stop = j; break; }
                        }
                        if (next_stop - k >= std::max<size_t>(settings.min_sense_codons, (L * 4) / 5)) {
                            start_idx = k;
                            first_stop = next_stop;
                            if (out_utr_trimmed) *out_utr_trimmed = true;
                            break;
                        }
                    }
                }
            }

            // If the first stop codon after start_idx is either the final codon OR lies within a trailing 3' UTR tail
            // (<= min(25 codons, 12.5% of L) from the 3' end of a >=50-codon ORF), treat it as the terminal stop codon
            // and trim any 3' UTR codons beyond it (E1: full-segment viral/GenBank records with 3' packaging tails).
            size_t end_idx = L;
            size_t max_3utr_codons = (L >= 50) ? std::min<size_t>(25, std::max<size_t>(1, L / 8)) : 0;
            if (first_stop >= start_idx && first_stop < L) {
                if (first_stop + 1 == L) {
                    end_idx = first_stop;
                } else if ((first_stop - start_idx) >= std::max<size_t>(settings.min_sense_codons, 50) &&
                           (L - 1 - first_stop) <= max_3utr_codons) {
                    end_idx = first_stop;
                    if (out_utr_trimmed) *out_utr_trimmed = true;
                }
            }

            clean.reserve((end_idx - start_idx) * 3);
            for (size_t k = start_idx; k < end_idx; ++k) {
                if (is_stop[k]) {
                    if (k + 1 == end_idx) break;
                    out_internal_stops++;
                }
                clean += trips[k];
                out_codons.push_back(scs[k]);
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
            bool utr_trimmed = false;
            size_t stops_frame0 = 0;
            size_t stops_chosen = 0;
            bool passed_gate2 = false;
        };

        std::vector<CandidateTaxon> candidates(raw_seqs.size());
        report.entries.resize(raw_seqs.size());

        // First pass: check frame +0 on all sequences and build a per-sequence amino-acid 3-mer frequency map
        // so that leave-one-out 3-mer matching can exclude a sequence's own 3-mers (fixing T1-03!).
        std::unordered_map<uint32_t, size_t> cohort_aa3mer_seq_counts;
        std::vector<bool> in_cohort_3mer_map(raw_seqs.size(), false);
        auto extract_unique_3mers = [&](const std::vector<int8_t>& cods) -> std::unordered_set<uint32_t> {
            std::unordered_set<uint32_t> u;
            if (cods.size() < 3) return u;
            u.reserve(cods.size());
            for (size_t k = 0; k + 2 < cods.size(); ++k) {
                if (cods[k] < 0 || cods[k + 1] < 0 || cods[k + 2] < 0) continue;
                uint32_t h = (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k]])) << 16)
                           | (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 1]])) << 8)
                           |  static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 2]]));
                u.insert(h);
            }
            return u;
        };
        auto add_seq_to_cohort_3mers = [&](size_t seq_idx, const std::vector<int8_t>& cods) {
            if (in_cohort_3mer_map[seq_idx]) return;
            auto u = extract_unique_3mers(cods);
            for (uint32_t h : u) cohort_aa3mer_seq_counts[h]++;
            in_cohort_3mer_map[seq_idx] = true;
        };
        auto count_other_cohort_3mers = [&](const std::vector<int8_t>& cods, const std::unordered_set<uint32_t>* self_set) -> size_t {
            if (cods.size() < 3 || cohort_aa3mer_seq_counts.empty()) return 0;
            size_t hits = 0;
            for (size_t k = 0; k + 2 < cods.size(); ++k) {
                if (cods[k] < 0 || cods[k + 1] < 0 || cods[k + 2] < 0) continue;
                uint32_t h = (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k]])) << 16)
                           | (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 1]])) << 8)
                           |  static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 2]]));
                auto it = cohort_aa3mer_seq_counts.find(h);
                if (it != cohort_aa3mer_seq_counts.end()) {
                    size_t self_contrib = (self_set && self_set->count(h)) ? 1 : 0;
                    if (it->second > self_contrib) hits++;
                }
            }
            return hits;
        };

        for (size_t t = 0; t < raw_seqs.size(); ++t) {
            std::string clean_raw;
            clean_raw.reserve(raw_seqs[t].raw_seq.size());
            bool has_gaps = false;
            for (char ch : raw_seqs[t].raw_seq) {
                if (std::isspace(static_cast<unsigned char>(ch))) continue;
                char uc = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
                if (uc == '-') has_gaps = true;
                if (uc == 'U') uc = 'T';
                clean_raw += uc;
            }
            std::string ungapped;
            ungapped.reserve(clean_raw.size());
            size_t n_ambig = 0;
            if (has_gaps && clean_raw.size() >= 3 && (clean_raw.size() % 3 == 0)) {
                for (size_t k = 0; k + 2 < clean_raw.size(); k += 3) {
                    char c1 = clean_raw[k], c2 = clean_raw[k + 1], c3 = clean_raw[k + 2];
                    if (c1 == '-' || c2 == '-' || c3 == '-') continue;
                    auto norm_nt = [&](char c) -> char {
                        if (c == 'A' || c == 'C' || c == 'G' || c == 'T') return c;
                        return 'N';
                    };
                    char n1 = norm_nt(c1), n2 = norm_nt(c2), n3 = norm_nt(c3);
                    if (n1 == 'N' && n2 == 'N' && n3 == 'N') continue;
                    if (n1 == 'N') n_ambig++;
                    if (n2 == 'N') n_ambig++;
                    if (n3 == 'N') n_ambig++;
                    ungapped += n1;
                    ungapped += n2;
                    ungapped += n3;
                }
            } else {
                for (char uc : clean_raw) {
                    if (uc == '-') continue;
                    if (uc != 'A' && uc != 'C' && uc != 'G' && uc != 'T') {
                        uc = 'N';
                        n_ambig++;
                    }
                    ungapped += uc;
                }
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
            bool utr_trim0 = false;
            std::string cds0 = extract_frame_cds(ct.raw_nt, 0, stops0, cods0, &utr_trim0);
            ent.internal_stops_frame0 = stops0;
            ct.stops_frame0 = stops0;
            ct.utr_trimmed = utr_trim0;

            if (stops0 == 0 && cods0.size() >= settings.min_sense_codons) {
                ct.clean_cds = std::move(cds0);
                ct.codons = std::move(cods0);
                ct.frame_offset = 0;
                ct.rev_comp = false;
                ct.stops_chosen = 0;
                ct.passed_gate2 = true;
                ent.codon_length = ct.codons.size();
                if (utr_trim0) {
                    ent.reason = "Clean in-frame CDS (UTR flank trimmed)";
                }
                if (ct.clean_cds.size() >= 3 && ct.clean_cds.substr(0, 3) == "ATG") {
                    add_seq_to_cohort_3mers(t, ct.codons);
                }
            }
        }
        // Fallback if fewer than 2 sequences started with ATG: seed 3-mer table from all stop-free frame-0 sequences
        size_t atg_seeded = 0;
        for (bool b : in_cohort_3mer_map) if (b) atg_seeded++;
        if (atg_seeded < 2) {
            for (size_t t = 0; t < candidates.size(); ++t) {
                if (candidates[t].passed_gate2) add_seq_to_cohort_3mers(t, candidates[t].codons);
            }
        }

        // Second pass: scan all 6 frames (+0,+1,+2 on fwd & revcomp) for any sequence that either:
        //   (a) had internal stop codons in frame +0, OR
        //   (b) did not start with ATG and has near-zero leave-one-out 3-mer homology in frame +0
        //       (e.g. a reverse-complemented sequence that happened to have 0 stops on the antisense strand, T1-03!)
        size_t gate2_failures = 0;
        for (size_t t = 0; t < raw_seqs.size(); ++t) {
            CandidateTaxon& ct = candidates[t];
            TaxonAdmissionEntry& ent = report.entries[t];
            if (ent.status == "Quarantined") {
                gate2_failures++;
                continue;
            }
            std::unordered_set<uint32_t> self_3mers;
            if (in_cohort_3mer_map[t]) self_3mers = extract_unique_3mers(ct.codons);
            const std::unordered_set<uint32_t>* self_ptr = in_cohort_3mer_map[t] ? &self_3mers : nullptr;

            size_t hits_frame0 = ct.passed_gate2 ? count_other_cohort_3mers(ct.codons, self_ptr) : 0;
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
                        size_t min_rescue_hits = std::max<size_t>(8, cods.size() / 8);
                        size_t hits = count_other_cohort_3mers(cods, self_ptr);
                        bool st_atg = (cds.size() >= 3 && cds.substr(0, 3) == "ATG");
                        if ((!cohort_aa3mer_seq_counts.empty() && hits > best_3mer_hits + 4 && hits >= min_rescue_hits) ||
                            (cohort_aa3mer_seq_counts.empty() && st_atg && best_frame == -1 && !ct.passed_gate2)) {
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

            size_t req_hits = std::max<size_t>(8, best_cods.size() / 8);
            if (best_frame != -1 && ((cohort_aa3mer_seq_counts.empty() && best_starts_atg) || best_3mer_hits >= req_hits)) {
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
                if (best_starts_atg) add_seq_to_cohort_3mers(t, ct.codons);
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
            std::string code_hint;
            const std::vector<std::string> alt_codes = {
                "Vertebrate-mtDNA", "Invertebrate-mtDNA", "Yeast-mtDNA", "Mold-Protozoan-mtDNA",
                "Ciliate-Nuclear", "Euplotid-Nuclear", "Echinoderm-mtDNA", "Ascidian-mtDNA",
                "Flatworm-mtDNA", "Blepharisma-Nuclear", "Chlorophycean-mtDNA", "Trematode-mtDNA"
            };
            size_t best_clean = gate2_survivors;
            std::string best_code_name;
            for (const auto& cname : alt_codes) {
                auto alt_gc = GeneticCode::from_name(cname);
                if (!alt_gc) continue;
                size_t clean_cnt = 0;
                for (const auto& ct : candidates) {
                    if (ct.raw_nt.size() < settings.min_sense_codons * 3) continue;
                    size_t usable = (ct.raw_nt.size() / 3) * 3;
                    bool has_st = false;
                    for (size_t i = 0; i + 5 < usable; i += 3) {
                        char c1 = ct.raw_nt[i], c2 = ct.raw_nt[i + 1], c3 = ct.raw_nt[i + 2];
                        if (c1 == 'N' || c2 == 'N' || c3 == 'N') continue;
                        if (alt_gc->is_stop_triplet(c1, c2, c3)) { has_st = true; break; }
                    }
                    if (!has_st) clean_cnt++;
                }
                if (clean_cnt > best_clean && clean_cnt >= 3 &&
                    static_cast<double>(raw_seqs.size() - clean_cnt) / static_cast<double>(raw_seqs.size()) <= 0.25) {
                    best_clean = clean_cnt;
                    best_code_name = cname;
                }
            }
            if (!best_code_name.empty()) {
                code_hint = " Hint: Re-running with '--code " + best_code_name + "' resolves internal stop codons in " +
                            std::to_string(best_clean) + "/" + std::to_string(raw_seqs.size()) +
                            " sequences (possible mitochondrial or non-standard genetic code).";
            }
            std::ostringstream err;
            err << "Gate 2 (Reading-Frame & Stop-Codon Gate): " << gate2_failures << " of " << raw_seqs.size()
                << " sequences (" << std::fixed << std::setprecision(1)
                << (100.0 * static_cast<double>(gate2_failures) / static_cast<double>(raw_seqs.size()))
                << "%) contain internal premature stop codons across reading frames or are too short. "
                << "Input likely contains un-spliced genomic introns, 5'/3' UTRs, or non-coding sequences."
                << code_hint;
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
        Matrix DP_mat = Matrix::Zero(M, M);
        Matrix MuID_mat = Matrix::Zero(M, M);
        Eigen::MatrixXi eval_mask = Eigen::MatrixXi::Zero(M, M);

        const bool use_subquadratic_gate3 = (M > settings.subquadratic_taxon_threshold);
        std::vector<std::vector<int8_t>> cand_codons(M);
        for (size_t a = 0; a < M; ++a) {
            cand_codons[a] = candidates[length_ok_indices[a]].codons;
        }
        LandmarkKmerGraph lm_graph;
        if (use_subquadratic_gate3) {
            lm_graph = LandmarkKmerGraph::build(cand_codons, *gcode, 8);
        }

        auto eval_pair_with_decoys = [&](size_t a, size_t b,
                                         std::vector<Scalar>& prev_M, std::vector<Scalar>& curr_M,
                                         std::vector<Scalar>& prev_D, std::vector<Scalar>& curr_D,
                                         std::vector<size_t>& prev_id_M, std::vector<size_t>& curr_id_M,
                                         std::vector<size_t>& prev_id_D, std::vector<size_t>& curr_id_D,
                                         std::vector<int8_t>& shuf_i, std::vector<int8_t>& shuf_j) {
            const auto& seq_i = cand_codons[a];
            const auto& seq_j = cand_codons[b];
            size_t Li = seq_i.size();
            size_t Lj = seq_j.size();
            if (Li == 0 || Lj == 0) return;

            auto [norm_dp, aa_id] = run_semiglobal_banded_dp_pair(
                seq_i, seq_j, *gcode, prev_M, curr_M, prev_D, curr_D, prev_id_M, curr_id_M, prev_id_D, curr_id_D);
            Scalar L_eff = std::max<Scalar>(1.0, static_cast<Scalar>(std::min(Li, Lj)));

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
                auto [dp_k, id_k] = run_semiglobal_banded_dp_pair(
                    shuf_i, shuf_j, *gcode, prev_M, curr_M, prev_D, curr_D, prev_id_M, curr_id_M, prev_id_D, curr_id_D);
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

            Eigen::Index ea = static_cast<Eigen::Index>(a);
            Eigen::Index eb = static_cast<Eigen::Index>(b);
            Z_mat(ea, eb) = Z_mat(eb, ea) = z_comb;
            ID_mat(ea, eb) = ID_mat(eb, ea) = aa_id;
            DP_mat(ea, eb) = DP_mat(eb, ea) = norm_dp;
            MuID_mat(ea, eb) = MuID_mat(eb, ea) = mu_id;
            eval_mask(ea, eb) = eval_mask(eb, ea) = 1;
        };

        if (!use_subquadratic_gate3) {
            #pragma omp parallel for schedule(dynamic)
            for (size_t a = 0; a < M; ++a) {
                std::vector<Scalar> prev_M, curr_M, prev_D, curr_D;
                std::vector<size_t> prev_id_M, curr_id_M, prev_id_D, curr_id_D;
                std::vector<int8_t> shuf_i, shuf_j;
                for (size_t b = a + 1; b < M; ++b) {
                    eval_pair_with_decoys(
                        a, b, prev_M, curr_M, prev_D, curr_D, prev_id_M, curr_id_M, prev_id_D, curr_id_D, shuf_i, shuf_j);
                }
            }
        } else {
            // Subquadratic Gate 3:
            // 1. Calibrate decoy excursion parameters on a compact landmark sample (<= 32 pairs)
            size_t K_lm = lm_graph.landmarks.size();
            std::vector<std::pair<size_t, size_t>> calib_pairs;
            for (size_t k = 0; k < K_lm && calib_pairs.size() < 32; ++k) {
                size_t u = lm_graph.landmarks[k];
                size_t v = lm_graph.landmarks[(k + 1) % K_lm];
                if (u != v) {
                    if (u > v) std::swap(u, v);
                    calib_pairs.push_back({u, v});
                }
            }
            Scalar sum_c_id = 0.0, sum_c_dp = 0.0;
            size_t n_calib = calib_pairs.size();
            {
                std::vector<Scalar> prev_M, curr_M, prev_D, curr_D;
                std::vector<size_t> prev_id_M, curr_id_M, prev_id_D, curr_id_D;
                std::vector<int8_t> shuf_i, shuf_j;
                for (const auto& [u, v] : calib_pairs) {
                    eval_pair_with_decoys(
                        u, v, prev_M, curr_M, prev_D, curr_D, prev_id_M, curr_id_M, prev_id_D, curr_id_D, shuf_i, shuf_j);
                    Scalar L_eff = std::max<Scalar>(
                        1.0, static_cast<Scalar>(std::min(cand_codons[u].size(), cand_codons[v].size())));
                    Scalar p0_uv = lm_graph.expected_aa_identity(u, v);
                    Scalar mu_id_uv = MuID_mat(static_cast<Eigen::Index>(u), static_cast<Eigen::Index>(v));
                    sum_c_id += (mu_id_uv - p0_uv) * std::sqrt(L_eff);

                    // Estimate decoy DP baseline from inverted Z_comb
                    Scalar aa_id_uv = ID_mat(static_cast<Eigen::Index>(u), static_cast<Eigen::Index>(v));
                    Scalar norm_dp_uv = DP_mat(static_cast<Eigen::Index>(u), static_cast<Eigen::Index>(v));
                    Scalar z_uv = Z_mat(static_cast<Eigen::Index>(u), static_cast<Eigen::Index>(v));
                    Scalar p_clamped = std::clamp<Scalar>(mu_id_uv, 0.05, 0.95);
                    Scalar sd_id = std::sqrt((17.0 / 16.0) * (p_clamped * (1.0 - p_clamped)) / L_eff);
                    Scalar sd_dp = std::sqrt((17.0 / 16.0) * (0.35 * 0.35) / L_eff);
                    Scalar z_id_est = (aa_id_uv - mu_id_uv) / std::max<Scalar>(1e-4, sd_id);
                    Scalar z_dp_est = (z_uv - 0.65 * z_id_est) / 0.35;
                    Scalar mu_dp_est = norm_dp_uv - z_dp_est * sd_dp;
                    sum_c_dp += mu_dp_est;
                }
            }
            Scalar c_id_hat = (n_calib > 0) ? (sum_c_id / static_cast<Scalar>(n_calib)) : 0.85;
            Scalar mu_dp_hat = (n_calib > 0) ? (sum_c_dp / static_cast<Scalar>(n_calib)) : -0.45;

            // 2. Evaluate only sparse landmark + k-NN edges in O(M sqrt(M)) using 1 DP per edge
            std::vector<std::pair<size_t, size_t>> sparse_edges;
            for (size_t a = 0; a < M; ++a) {
                for (size_t b : lm_graph.sparse_adj[a]) {
                    if (a < b && eval_mask(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) == 0) {
                        sparse_edges.push_back({a, b});
                    }
                }
            }

            #pragma omp parallel
            {
                std::vector<Scalar> prev_M, curr_M, prev_D, curr_D;
                std::vector<size_t> prev_id_M, curr_id_M, prev_id_D, curr_id_D;

                #pragma omp for schedule(dynamic, 16)
                for (size_t e = 0; e < sparse_edges.size(); ++e) {
                    size_t a = sparse_edges[e].first;
                    size_t b = sparse_edges[e].second;
                    const auto& seq_i = cand_codons[a];
                    const auto& seq_j = cand_codons[b];
                    size_t Li = seq_i.size();
                    size_t Lj = seq_j.size();
                    if (Li == 0 || Lj == 0) continue;

                    auto [norm_dp, aa_id] = run_semiglobal_banded_dp_pair(
                        seq_i, seq_j, *gcode, prev_M, curr_M, prev_D, curr_D, prev_id_M, curr_id_M, prev_id_D, curr_id_D);
                    Scalar L_eff = std::max<Scalar>(1.0, static_cast<Scalar>(std::min(Li, Lj)));

                    Scalar p0_ab = lm_graph.expected_aa_identity(a, b);
                    Scalar mu_id = std::clamp<Scalar>(p0_ab + c_id_hat / std::sqrt(L_eff), 0.04, 0.90);
                    Scalar mu_dp = mu_dp_hat;

                    Scalar p_clamped = std::clamp<Scalar>(mu_id, 0.05, 0.95);
                    Scalar var_id_reg = (17.0 / 16.0) * (p_clamped * (1.0 - p_clamped)) / L_eff;
                    Scalar var_dp_reg = (17.0 / 16.0) * (0.35 * 0.35) / L_eff;

                    Scalar z_id = (aa_id - mu_id) / std::max<Scalar>(1e-4, std::sqrt(var_id_reg));
                    Scalar z_dp = (norm_dp - mu_dp) / std::max<Scalar>(1e-4, std::sqrt(var_dp_reg));
                    Scalar z_comb = 0.65 * z_id + 0.35 * z_dp;

                    Eigen::Index ea = static_cast<Eigen::Index>(a);
                    Eigen::Index eb = static_cast<Eigen::Index>(b);
                    Z_mat(ea, eb) = Z_mat(eb, ea) = z_comb;
                    ID_mat(ea, eb) = ID_mat(eb, ea) = aa_id;
                    DP_mat(ea, eb) = DP_mat(eb, ea) = norm_dp;
                    MuID_mat(ea, eb) = MuID_mat(eb, ea) = mu_id;
                    eval_mask(ea, eb) = eval_mask(eb, ea) = 1;
                }
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
                Eigen::Index eu = static_cast<Eigen::Index>(u);
                if (use_subquadratic_gate3) {
                    for (size_t v : lm_graph.sparse_adj[u]) {
                        Eigen::Index ev = static_cast<Eigen::Index>(v);
                        if (eval_mask(eu, ev) == 0) continue;
                        Scalar min_id_edge = std::max<Scalar>(0.18, MuID_mat(eu, ev) + 0.04);
                        if (comp_id[v] == -1 &&
                            Z_mat(eu, ev) >= settings.min_homology_z &&
                            ID_mat(eu, ev) >= min_id_edge) {
                            comp_id[v] = cid;
                            q.push_back(v);
                        }
                    }
                } else {
                    for (size_t v = 0; v < M; ++v) {
                        Eigen::Index ev = static_cast<Eigen::Index>(v);
                        if (eval_mask(eu, ev) == 0) continue;
                        Scalar min_id_edge = std::max<Scalar>(0.18, MuID_mat(eu, ev) + 0.04);
                        if (comp_id[v] == -1 &&
                            Z_mat(eu, ev) >= settings.min_homology_z &&
                            ID_mat(eu, ev) >= min_id_edge) {
                            comp_id[v] = cid;
                            q.push_back(v);
                        }
                    }
                }
            }
            if (q.size() > best_comp_size) {
                best_comp_size = q.size();
                best_comp = cid;
            }
        }

        // Exact 16-decoy verification for any candidate disconnected from best_comp in subquadratic mode
        if (use_subquadratic_gate3 && best_comp >= 0) {
            std::vector<Scalar> prev_M, curr_M, prev_D, curr_D;
            std::vector<size_t> prev_id_M, curr_id_M, prev_id_D, curr_id_D;
            std::vector<int8_t> shuf_i, shuf_j;
            for (size_t a = 0; a < M; ++a) {
                if (comp_id[a] == best_comp) continue;
                size_t best_b = a;
                Scalar best_z = -1e9;
                for (size_t b : lm_graph.sparse_adj[a]) {
                    if (comp_id[b] == best_comp) {
                        Scalar z_ab = Z_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b));
                        if (z_ab > best_z) {
                            best_z = z_ab;
                            best_b = b;
                        }
                    }
                }
                if (best_b != a) {
                    eval_pair_with_decoys(
                        a, best_b, prev_M, curr_M, prev_D, curr_D, prev_id_M, curr_id_M, prev_id_D, curr_id_D, shuf_i, shuf_j);
                    Eigen::Index ea = static_cast<Eigen::Index>(a);
                    Eigen::Index eb = static_cast<Eigen::Index>(best_b);
                    Scalar min_id_edge = std::max<Scalar>(0.18, MuID_mat(ea, eb) + 0.04);
                    if (Z_mat(ea, eb) >= settings.min_homology_z && ID_mat(ea, eb) >= min_id_edge) {
                        comp_id[a] = best_comp;
                    }
                }
            }
        }

        for (size_t a = 0; a < M; ++a) {
            size_t orig = length_ok_indices[a];
            Scalar max_z = -1e9;
            Scalar max_id = 0.0;
            for (size_t b = 0; b < M; ++b) {
                if (a == b || eval_mask(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)) == 0) continue;
                max_z = std::max(max_z, Z_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)));
                max_id = std::max(max_id, ID_mat(static_cast<Eigen::Index>(a), static_cast<Eigen::Index>(b)));
            }
            report.entries[orig].max_homology_z = (max_z > -1e8) ? max_z : 0.0;
            report.entries[orig].max_aa_identity = max_id;

            if (settings.enforce_admission_control && comp_id[a] != best_comp) {
                report.entries[orig].status = "Quarantined";
                report.entries[orig].frame_offset = 0;
                report.entries[orig].reverse_complemented = false;
                std::ostringstream rs;
                if (candidates[orig].stops_frame0 > 0) {
                    rs << "Internal premature stop codons in frame +0 (" << candidates[orig].stops_frame0
                       << " stops; max Z = " << std::fixed << std::setprecision(2)
                       << report.entries[orig].max_homology_z << ")";
                } else {
                    rs << "Non-homologous outlier (max Z = " << std::fixed << std::setprecision(2)
                       << report.entries[orig].max_homology_z << ", max AA identity = "
                       << std::setprecision(1) << (max_id * 100.0) << "%)";
                }
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

        size_t N_adm = admitted_m_indices.size();
        report.admitted_dp_similarity = Matrix::Ones(N_adm, N_adm);
        std::vector<Scalar> cohort_z, cohort_id;
        for (size_t i = 0; i < N_adm; ++i) {
            for (size_t j = i + 1; j < N_adm; ++j) {
                Eigen::Index u = static_cast<Eigen::Index>(admitted_m_indices[i]);
                Eigen::Index v = static_cast<Eigen::Index>(admitted_m_indices[j]);
                if (eval_mask(u, v) != 0) {
                    cohort_z.push_back(Z_mat(u, v));
                    cohort_id.push_back(ID_mat(u, v));
                    report.admitted_dp_similarity(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = DP_mat(u, v);
                    report.admitted_dp_similarity(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = DP_mat(u, v);
                }
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
        std::string raw = Alignment::read_file_or_gzip(filepath);
        std::istringstream file(raw);

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
            } else if (line[0] == '(') {
                break;
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
        if (seqs.empty()) {
            // Fallback for PHYLIP / NEXUS / gzipped files passed to --unaligned
            Alignment loaded = Alignment::load(filepath, gcode);
            return strip_gaps_to_unaligned(loaded);
        }
        return seqs;
    }

    JointAlignResult run_from_unaligned(
        const std::vector<Sequence>& unaligned_seqs,
        std::shared_ptr<const GeneticCode> gcode = nullptr,
        std::optional<Tree> starting_tree = std::nullopt
    ) {
        auto t0 = std::chrono::steady_clock::now();
        auto [admitted_seqs, adm_report] = triage_and_admit_unaligned_sequences(unaligned_seqs, gcode, settings);
        auto [de_novo_aln, de_novo_tree] = build_de_novo_codon_alignment_and_tree(
            admitted_seqs, gcode, settings, &adm_report.admitted_dp_similarity);
        auto res = run(de_novo_aln, starting_tree.has_value() ? starting_tree : std::optional<Tree>(de_novo_tree));
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
        LikelihoodEngine::reconcile_tree_leaf_names(tree, base_aln);

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
            constexpr int B_cv = 24;
            constexpr int BAND_W = 2 * B_cv + 1;
            const Scalar NEG_INF_CV = -1e18;

            for (size_t t = 0; t < N; ++t) {
                std::vector<int8_t> obs;
                obs.reserve(C);
                std::vector<int> cum_obs(C + 1, 0);
                for (size_t c = 0; c < C; ++c) {
                    if (aln.sequences[t].codons[c] >= 0) obs.push_back(aln.sequences[t].codons[c]);
                    cum_obs[c + 1] = static_cast<int>(obs.size());
                }
                int M = static_cast<int>(obs.size());
                Eigen::MatrixXd fwd_M = Eigen::MatrixXd::Constant(C + 1, BAND_W, NEG_INF_CV);
                Eigen::MatrixXd fwd_G = Eigen::MatrixXd::Constant(C + 1, BAND_W, NEG_INF_CV);
                fwd_M(0, B_cv) = 0.0;

                for (size_t w = 1; w <= C; ++w) {
                    const Vector& cav = full_cavity[w - 1][t];
                    int center_prev = cum_obs[w - 1];
                    int center_curr = cum_obs[w];
                    int j_lo = std::max(0, center_curr - B_cv);
                    int j_hi = std::min(M, center_curr + B_cv);
                    for (int j = j_lo; j <= j_hi; ++j) {
                        int k_curr = j - center_curr + B_cv;
                        int k_prev_same = j - center_prev + B_cv;
                        if (k_prev_same >= 0 && k_prev_same < BAND_W) {
                            fwd_G(w, k_curr) = log_sum_exp(
                                fwd_M(w - 1, k_prev_same) + log_open_cv,
                                fwd_G(w - 1, k_prev_same) + log_ext_cv
                            );
                        }
                        int k_prev_diag = (j - 1) - center_prev + B_cv;
                        if (j >= 1 && k_prev_diag >= 0 && k_prev_diag < BAND_W) {
                            Scalar emit = std::log(std::max<Scalar>(1e-12, cav(obs[static_cast<size_t>(j - 1)])));
                            fwd_M(w, k_curr) = log_sum_exp(
                                fwd_M(w - 1, k_prev_diag) + log_1m_open_cv + emit,
                                fwd_G(w - 1, k_prev_diag) + log_1m_ext_cv + emit
                            );
                        }
                    }
                }
                Scalar log_Z = log_sum_exp(fwd_M(C, B_cv), fwd_G(C, B_cv));
                ev.loto_cv_log_l += log_Z;

                // Backward pass for Option 5 per-column posterior confidence
                Eigen::MatrixXd bwd_M = Eigen::MatrixXd::Constant(C + 1, BAND_W, NEG_INF_CV);
                Eigen::MatrixXd bwd_G = Eigen::MatrixXd::Constant(C + 1, BAND_W, NEG_INF_CV);
                bwd_M(C, B_cv) = 0.0;
                bwd_G(C, B_cv) = 0.0;
                for (size_t w = C; w >= 1; --w) {
                    const Vector& cav = full_cavity[w - 1][t];
                    int center_prev = cum_obs[w - 1];
                    int center_curr = cum_obs[w];
                    int j_lo = std::max(0, center_curr - B_cv);
                    int j_hi = std::min(M, center_curr + B_cv);
                    for (int j = j_lo; j <= j_hi; ++j) {
                        int k_curr = j - center_curr + B_cv;
                        int k_prev_same = j - center_prev + B_cv;
                        if (k_prev_same >= 0 && k_prev_same < BAND_W) {
                            Scalar to_g = bwd_G(w, k_curr);
                            bwd_M(w - 1, k_prev_same) = log_sum_exp(bwd_M(w - 1, k_prev_same), to_g + log_open_cv);
                            bwd_G(w - 1, k_prev_same) = log_sum_exp(bwd_G(w - 1, k_prev_same), to_g + log_ext_cv);
                        }
                        int k_prev_diag = (j - 1) - center_prev + B_cv;
                        if (j >= 1 && k_prev_diag >= 0 && k_prev_diag < BAND_W) {
                            Scalar emit = std::log(std::max<Scalar>(1e-12, cav(obs[static_cast<size_t>(j - 1)])));
                            Scalar to_m = bwd_M(w, k_curr) + emit;
                            bwd_M(w - 1, k_prev_diag) = log_sum_exp(bwd_M(w - 1, k_prev_diag), to_m + log_1m_open_cv);
                            bwd_G(w - 1, k_prev_diag) = log_sum_exp(bwd_G(w - 1, k_prev_diag), to_m + log_1m_ext_cv);
                        }
                    }
                }
                for (size_t w = 1; w <= C; ++w) {
                    Scalar log_m = NEG_INF_CV;
                    int center_curr = cum_obs[w];
                    int j_lo = std::max(1, center_curr - B_cv);
                    int j_hi = std::min(M, center_curr + B_cv);
                    for (int j = j_lo; j <= j_hi; ++j) {
                        int k_curr = j - center_curr + B_cv;
                        log_m = log_sum_exp(log_m, fwd_M(w, k_curr) + bwd_M(w, k_curr));
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
