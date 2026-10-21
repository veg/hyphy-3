#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "nlohmann/json.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/stats.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <numeric>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

struct ABSRELRateDistribution {
    std::vector<Scalar> rates;   // omega values
    std::vector<Scalar> weights; // mixture proportions summing to 1.0
};

struct ABSRELBranchResult {
    std::string branch_name;
    int32_t node_id = -1;
    Scalar baseline_omega = 1.0;
    Scalar baseline_branch_length = 0.0;
    int rate_classes = 1;
    ABSRELRateDistribution rate_distribution;
    Scalar full_branch_length = 0.0;
    Scalar full_es = 0.0; // synonymous subs/site
    Scalar full_en = 0.0; // nonsynonymous subs/site
    Scalar lrt = 0.0;
    Scalar uncorrected_p_value = 1.0;
    Scalar corrected_p_value = 1.0;
    bool is_tested = true;
    bool is_positive = false;
    int sites_ebf_100 = 0;
    std::vector<Scalar> site_ebf;
};

struct ABSRELFitSummary {
    Scalar log_likelihood = 0.0;
    size_t parameters = 0;
    Scalar aicc = 0.0;
};

struct ABSRELSettings {
    int max_rate_classes = 3;
    Scalar p_threshold = 0.05;
    std::string test_branches = "All"; // "All", "Internal", "Leaves"
    bool verbose = false;
};

struct ABSRELResult {
    ABSRELFitSummary gtr_fit;
    ABSRELFitSummary baseline_fit;
    ABSRELFitSummary full_adaptive_fit;
    GTRParameters gtr_params; // populated so the GTR fit can emit substitution rates
    std::unordered_map<std::string, ABSRELBranchResult> branches;
    std::vector<std::string> tested_branches;
    std::vector<std::string> positive_branches;
    Scalar p_threshold = 0.05;
    double runtime_seconds = 0.0;

nlohmann::json to_legacy_json(const Tree& tree, const Alignment& aln) const;

nlohmann::json to_modern_json(const Tree& tree, const Alignment& aln, const Provenance& prov = {}) const;
    nlohmann::json to_json(
        const Tree& tree,
        const Alignment& aln,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        if (format == JSONFormat::Legacy) {
            return to_legacy_json(tree, aln);
        }
        return to_modern_json(tree, aln, prov);
    }

    void save_json(
        const std::string& filepath,
        const Tree& tree,
        const Alignment& aln,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open file for writing: " + filepath);
        }
        nlohmann::json j = to_json(tree, aln, format, prov);
        out << j.dump(2) << "\n";
    }

};

class ABSRELAnalyzer {
public:
    Tree tree;
    Alignment aln;
    ABSRELSettings settings;

    const GeneticCode* gcode = nullptr;
    Eigen::Matrix<Scalar, 3, 4> pos_nuc_freqs;
    Vector codon_freqs;
    GTRParameters gtr_params;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;

    std::vector<size_t> leaf_to_taxon;

    // Per-branch parameters:
    // branch_syn_alpha[id]: synonymous rate multiplier alpha_b
    std::vector<Scalar> branch_alpha;
    // branch_omegas[id]: omega rate categories for branch id
    std::vector<std::vector<Scalar>> branch_omegas;
    // branch_weights[id]: mixture proportions for branch id
    std::vector<std::vector<Scalar>> branch_weights;

    // Scale factors and precomputed rate templates for Q(1, omega) = Q_syn + omega * Q_nonsyn
    Scalar syn_scale_factor = 1.0;
    Scalar nonsyn_scale_factor = 1.0;
    Matrix Q_syn_template;
    Matrix Q_nonsyn_template;
    Matrix S_syn_template;
    Matrix S_nonsyn_template;
    Vector sqrt_pi_vec;
    Vector inv_sqrt_pi_vec;
    Vector pattern_weights_vec;

    std::vector<bool> is_test_branch;

    ABSRELAnalyzer(Tree t, Alignment a, ABSRELSettings s = ABSRELSettings{})
        : tree(std::move(t)), aln(std::move(a)), settings(s) {
        gcode = aln.code ? aln.code.get() : GeneticCode::universal().get();
        pos_nuc_freqs = aln.pos_nuc_frequencies;
        codon_freqs = aln.codon_frequencies_f3x4;

        int S = gcode->num_sense_codons;
        sqrt_pi_vec.resize(S);
        inv_sqrt_pi_vec.resize(S);
        for (int i = 0; i < S; ++i) {
            sqrt_pi_vec(i) = std::sqrt(std::max(codon_freqs(i), 1e-12));
            inv_sqrt_pi_vec(i) = 1.0 / sqrt_pi_vec(i);
        }

        size_t num_patterns = aln.patterns.size();
        pattern_weights_vec.resize(num_patterns);
        for (size_t p = 0; p < num_patterns; ++p) {
            pattern_weights_vec(p) = aln.patterns[p].weight;
        }

        auto sel = tree.select_branches(settings.test_branches, "", {}, false);
        is_test_branch = std::move(sel.is_selected);

        size_t num_nodes = tree.num_nodes();
        leaf_to_taxon.assign(num_nodes, static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        branch_alpha.assign(num_nodes, 0.01);
        branch_omegas.resize(num_nodes);
        branch_weights.resize(num_nodes);
        for (size_t i = 0; i < num_nodes; ++i) {
            branch_omegas[i] = {0.5};
            branch_weights[i] = {1.0};
        }

        compute_syn_scale();
    }

    static ABSRELAnalyzer create(Tree t, Alignment a, ABSRELSettings s = ABSRELSettings{}) {
        return ABSRELAnalyzer(std::move(t), std::move(a), s);
    }

    void compute_syn_scale() {
        int S = gcode->num_sense_codons;
        Scalar r_nuc[4][4] = {
            {0.0, gtr_params.theta_AC, 1.0, gtr_params.theta_AT},
            {gtr_params.theta_AC, 0.0, gtr_params.theta_CG, gtr_params.theta_CT},
            {1.0, gtr_params.theta_CG, 0.0, gtr_params.theta_GT},
            {gtr_params.theta_AT, gtr_params.theta_CT, gtr_params.theta_GT, 0.0}
        };

        Q_syn_template = Matrix::Zero(S, S);
        Q_nonsyn_template = Matrix::Zero(S, S);

        for (int i = 0; i < S; ++i) {
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                const auto& diff = gcode->diff_matrix[i][j];
                if (diff.num_diffs == 1) {
                    Scalar rate = r_nuc[diff.diff_from[0]][diff.diff_to[0]] *
                                  pos_nuc_freqs(diff.diff_pos[0], diff.diff_to[0]);
                    if (diff.synonymous) {
                        Q_syn_template(i, j) = rate;
                    } else {
                        Q_nonsyn_template(i, j) = rate;
                    }
                }
            }
            Q_syn_template(i, i) = -Q_syn_template.row(i).sum();
            Q_nonsyn_template(i, i) = -Q_nonsyn_template.row(i).sum();
        }

        syn_scale_factor = -(codon_freqs.array() * Q_syn_template.diagonal().array()).sum();
        nonsyn_scale_factor = -(codon_freqs.array() * Q_nonsyn_template.diagonal().array()).sum();

        S_syn_template = sqrt_pi_vec.asDiagonal() * Q_syn_template * inv_sqrt_pi_vec.asDiagonal();
        S_syn_template = 0.5 * (S_syn_template + S_syn_template.transpose());
        S_nonsyn_template = sqrt_pi_vec.asDiagonal() * Q_nonsyn_template * inv_sqrt_pi_vec.asDiagonal();
        S_nonsyn_template = 0.5 * (S_nonsyn_template + S_nonsyn_template.transpose());
    }

    // Build rate matrix Q(1, omega) using pre-split synonymous and nonsynonymous templates
    MG94Matrix build_q_matrix(Scalar omega) const {
        MG94Matrix mat;
        mat.pi = codon_freqs;
        mat.sqrt_pi = sqrt_pi_vec;
        mat.inv_sqrt_pi = inv_sqrt_pi_vec;
        mat.Q = Q_syn_template + omega * Q_nonsyn_template;
        mat.scale_factor = syn_scale_factor + omega * nonsyn_scale_factor;

        Matrix S_mat = S_syn_template + omega * S_nonsyn_template;
        Eigen::SelfAdjointEigenSolver<Matrix> solver(S_mat);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("MG94 eigendecomposition failed");
        }

        mat.eigenvalues = solver.eigenvalues().cwiseMin(0.0);
        mat.eigenvectors = solver.eigenvectors();
        mat.inv_eigenvectors = mat.eigenvectors.transpose();
        mat.U_pi = inv_sqrt_pi_vec.asDiagonal() * mat.eigenvectors;
        mat.V_pi = mat.inv_eigenvectors * sqrt_pi_vec.asDiagonal();
        return mat;
    }

    // Evaluate transition matrix for branch node_id given its current parameters
    Matrix compute_branch_transition_matrix(int32_t node_id) const {
        int S = gcode->num_sense_codons;
        Scalar alpha = branch_alpha[node_id];
        const auto& omegas = branch_omegas[node_id];
        const auto& weights = branch_weights[node_id];
        size_t K = omegas.size();

        Matrix P_mix = Matrix::Zero(S, S);
        for (size_t k = 0; k < K; ++k) {
            MG94Matrix mat = build_q_matrix(omegas[k]);
            Matrix P_k = mat.transition_matrix(alpha);
            P_mix += weights[k] * P_k;
        }
        return P_mix;
    }

    // Precompute all branch transition matrices across the tree
    std::vector<Matrix> compute_all_branch_transition_matrices() const {
        return compute_all_branch_transition_matrices_with(branch_alpha, branch_omegas, branch_weights);
    }

    std::vector<Matrix> compute_all_branch_transition_matrices_with(
        const std::vector<Scalar>& alphas,
        const std::vector<std::vector<Scalar>>& omegas,
        const std::vector<std::vector<Scalar>>& weights
    ) const {
        size_t num_nodes = tree.num_nodes();
        int S = gcode->num_sense_codons;
        std::vector<Matrix> P_branches(num_nodes);
        std::vector<int32_t> active_branches;
        active_branches.reserve(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                active_branches.push_back(node.id);
            }
        }
        int n_threads = LikelihoodEngine::optimal_omp_threads(active_branches.size(), 4);
        #pragma omp parallel for num_threads(n_threads) schedule(static)
        for (size_t i = 0; i < active_branches.size(); ++i) {
            int32_t bid = active_branches[i];
            Scalar a = alphas[bid];
            P_branches[bid] = Matrix::Zero(S, S);
            for (size_t k = 0; k < omegas[bid].size(); ++k) {
                MG94Matrix mat = build_q_matrix(omegas[bid][k]);
                P_branches[bid] += weights[bid][k] * mat.transition_matrix(a);
            }
        }
        return P_branches;
    }

    // Compute full tree log-likelihood
    Scalar compute_tree_log_likelihood(const std::vector<Matrix>& P_branches) const {
        int S = gcode->num_sense_codons;
        return LikelihoodEngine::compute_tree_log_likelihood_from_P(
            tree, aln, P_branches, codon_freqs, S, &leaf_to_taxon
        );
    }

    // Compute expected branch length (substitutions per nucleotide site) for branch node_id
    Scalar compute_expected_branch_length(int32_t node_id) const {
        Scalar alpha = branch_alpha[node_id];
        const auto& omegas = branch_omegas[node_id];
        const auto& weights = branch_weights[node_id];
        Scalar total_sub_rate = 0.0;
        for (size_t k = 0; k < omegas.size(); ++k) {
            Scalar sf = syn_scale_factor + omegas[k] * nonsyn_scale_factor;
            total_sub_rate += weights[k] * (sf / 3.0);
        }
        return alpha * total_sub_rate;
    }

    // Compute synonymous and nonsynonymous component branch lengths
    std::pair<Scalar, Scalar> compute_es_en(int32_t node_id) const {
        Scalar alpha = branch_alpha[node_id];
        Scalar syn_rate = syn_scale_factor / 3.0;
        Scalar es = alpha * syn_rate;

        Scalar total_len = compute_expected_branch_length(node_id);
        Scalar en = std::max(0.0, total_len - es);
        return {es, en};
    }

    // AICc helper: N = num_sequences * num_codons
    Scalar compute_aicc(Scalar log_l, size_t num_params) const {
        return hyphy::core::stats::compute_aicc(log_l, num_params, aln.num_taxa * aln.num_codons);
    }

    // Asymptotic p-value for mixture distribution: 0.5 * [0.4 * chi2_1 + 0.6 * chi2_2]
    static Scalar compute_mixture_p_value(Scalar lrt) {
        return hyphy::core::stats::mixture_chi2_p_value(lrt, 0.20, 0.30);
    }

    // Vectorized evaluation of branch log-likelihood from packed [S x num_patterns] V_mat and D_mat
    Scalar eval_branch_ll_from_matrices(
        const Matrix& P_b,
        const Matrix& V_mat,
        const Matrix& D_mat,
        Scalar log_scale_offset,
        Matrix& PD_scratch
    ) const {
        PD_scratch.noalias() = P_b * D_mat;
        size_t num_patterns = aln.patterns.size();
        Scalar ll = log_scale_offset;
        Scalar c_ll = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            Scalar Lp = V_mat.col(p).dot(PD_scratch.col(p));
            Scalar term = (Lp > 0.0) ? (pattern_weights_vec(p) * std::log(Lp)) : (pattern_weights_vec(p) * (-1e20));
            LikelihoodEngine::neumaier_acc(ll, c_ll, term);
        }
        return ll + c_ll;
    }

    // Zero-allocation helper to compute V_mat [S x P], D_mat [S x P], and log_scale_offset for a single branch
    void compute_branch_context_matrices(
        int32_t bid,
        const std::vector<Matrix>& P_branches,
        Matrix& V_mat,
        Matrix& D_mat,
        Scalar& log_scale_offset
    ) const {
        int S = gcode->num_sense_codons;
        size_t num_patterns = aln.patterns.size();
        size_t num_nodes = tree.num_nodes();
        V_mat.resize(S, num_patterns);
        D_mat.resize(S, num_patterns);
        std::vector<Scalar> scale_offsets(num_patterns, 0.0);

        int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 16);
        constexpr size_t B_TILE = 256;
        size_t num_tiles = (num_patterns + B_TILE - 1) / B_TILE;

        #pragma omp parallel num_threads(n_threads)
        {
            LikelihoodEngine::BatchedInsideOutsideWorkspace ws;
            ws.init(num_nodes, S, std::min(num_patterns, B_TILE));

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * B_TILE;
                size_t p_end = std::min(p_start + B_TILE, num_patterns);
                size_t B = p_end - p_start;
                LikelihoodEngine::compute_inside_outside_batch_in_place(
                    tree, aln, p_start, p_end, leaf_to_taxon, P_branches, codon_freqs, S, ws
                );
                V_mat.middleCols(p_start, B) = ws.V[bid].leftCols(B);
                D_mat.middleCols(p_start, B) = ws.D[bid].leftCols(B);
                for (size_t b = 0; b < B; ++b) {
                    size_t p = p_start + b;
                    if (ws.likelihood(b) > 0.0) {
                        scale_offsets[p] = aln.patterns[p].weight * (ws.log_likelihood(b) - std::log(ws.likelihood(b)));
                    }
                }
            }
        }

        log_scale_offset = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            log_scale_offset += scale_offsets[p];
        }
    }

    // Phase 1: Fit GTR
    void fit_gtr() {
        GTRFitter gtr_fitter(tree, aln);
        auto gtr_res = gtr_fitter.fit();
        gtr_params = gtr_res.params;
        gtr_log_l = gtr_res.log_likelihood;
        gtr_aicc = gtr_res.aicc;

        compute_syn_scale();

        // Initialize branch alpha using nucleotide branch lengths
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != tree.root_id) {
                // Initial conversion from nucleotide substitutions to synonymous time
                Scalar bl = std::max(node.branch_length, 1e-4);
                branch_alpha[node.id] = bl * (3.0 / std::max(syn_scale_factor, 0.1));
                branch_omegas[node.id] = {0.3};
                branch_weights[node.id] = {1.0};
            }
        }
    }

    // Phase 2: Fit Baseline MG94xREV Model (each branch has 1 omega)
    ABSRELFitSummary fit_baseline(
        int max_cycles = 8,
        bool show_progress = false,
        bool force_progress = false,
        std::function<void(const std::string&, double)> progress_cb = nullptr
    ) {
        int S = gcode->num_sense_codons;
        size_t num_patterns = aln.patterns.size();
        size_t num_nodes = tree.num_nodes();
        size_t num_branches = num_nodes - 1;

        std::vector<int32_t> branch_ids;
        branch_ids.reserve(num_branches);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                branch_ids.push_back(node.id);
            }
        }

        std::unique_ptr<ProgressBar> pbar_p2;
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar_p2 = std::make_unique<ProgressBar>(
                max_cycles * num_branches,
                "[aBSREL] Baseline Model",
                "branches",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        std::vector<Matrix> P_branches = compute_all_branch_transition_matrices();
        Scalar prev_ll = compute_tree_log_likelihood(P_branches);

        std::vector<Matrix> V_all(num_nodes, Matrix::Zero(S, num_patterns));
        std::vector<Matrix> D_all(num_nodes, Matrix::Zero(S, num_patterns));

        for (int cycle = 0; cycle < max_cycles; ++cycle) {
            // 1. Batched GEMM Inside-Outside pass for all patterns across entire tree
            int io_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 16);
            constexpr size_t B_TILE = 256;
            size_t num_tiles = (num_patterns + B_TILE - 1) / B_TILE;

            #pragma omp parallel num_threads(io_threads)
            {
                LikelihoodEngine::BatchedInsideOutsideWorkspace ws;
                ws.init(num_nodes, S, std::min(num_patterns, B_TILE));

                #pragma omp for schedule(static)
                for (size_t t = 0; t < num_tiles; ++t) {
                    size_t p_start = t * B_TILE;
                    size_t p_end = std::min(p_start + B_TILE, num_patterns);
                    size_t B = p_end - p_start;
                    LikelihoodEngine::compute_inside_outside_batch_in_place(
                        tree, aln, p_start, p_end, leaf_to_taxon, P_branches, codon_freqs, S, ws
                    );
                    for (int32_t bid : branch_ids) {
                        V_all[bid].middleCols(p_start, B) = ws.V[bid].leftCols(B);
                        D_all[bid].middleCols(p_start, B) = ws.D[bid].leftCols(B);
                    }
                }
            }

            std::vector<Scalar> new_alpha(num_nodes);
            std::vector<Scalar> new_omega(num_nodes);

            // 2. Parallel branch optimization using cached Inside-Outside context matrices
            int br_threads = LikelihoodEngine::optimal_omp_threads(branch_ids.size(), 2);
            #pragma omp parallel num_threads(br_threads)
            {
                Matrix PD_scratch(S, num_patterns);
                #pragma omp for schedule(static)
                for (size_t bi = 0; bi < branch_ids.size(); ++bi) {
                    int32_t bid = branch_ids[bi];
                    const auto& node = tree.nodes[bid];

                    auto branch_obj = [&](Scalar log_alpha, Scalar log_omega) -> Scalar {
                        Scalar cur_alpha = std::exp(log_alpha);
                        Scalar cur_omega = std::exp(log_omega);
                        if (cur_alpha > 50.0 || cur_omega > 5000.0) return 1e20;

                        MG94Matrix mat = build_q_matrix(cur_omega);
                        Matrix P_b = mat.transition_matrix(cur_alpha);
                        return -eval_branch_ll_from_matrices(P_b, V_all[bid], D_all[bid], 0.0, PD_scratch);
                    };

                    Scalar init_log_alpha = std::log(std::clamp(branch_alpha[bid], 1e-4, 20.0));
                    Scalar init_log_omega = std::log(std::clamp(branch_omegas[bid][0], 1e-3, 50.0));

                    auto opt_res = NelderMead2D::minimize(branch_obj, init_log_alpha, init_log_omega, -10.0, 10.0, 1e-4, 80);
                    new_alpha[bid] = std::exp(opt_res.a);
                    new_omega[bid] = std::exp(opt_res.b);

                    if (pbar_p2) {
                        pbar_p2->tick();
                        std::ostringstream ss;
                        ss << node.name << " (c" << (cycle + 1) << "/" << max_cycles << " w=" << std::fixed << std::setprecision(2) << new_omega[bid] << ")";
                        pbar_p2->set_status(ss.str());
                    }
                }
            }

            // F-052: Monotonicity guard on Jacobi update (damped log-space backtracking, then Gauss-Seidel fallback)
            std::vector<Scalar> prev_alpha = branch_alpha;
            std::vector<Scalar> prev_omega_0(num_nodes, 0.5);
            for (int32_t bid : branch_ids) {
                prev_omega_0[bid] = branch_omegas[bid][0];
            }

            Scalar cur_ll = -std::numeric_limits<Scalar>::infinity();
            Scalar step = 1.0;
            bool accepted = false;
            for (int ls = 0; ls < 5; ++ls) {
                for (int32_t bid : branch_ids) {
                    Scalar la_prev = std::log(std::max(1e-8, prev_alpha[bid]));
                    Scalar lw_prev = std::log(std::max(1e-8, prev_omega_0[bid]));
                    Scalar la_new = std::log(std::max(1e-8, new_alpha[bid]));
                    Scalar lw_new = std::log(std::max(1e-8, new_omega[bid]));
                    branch_alpha[bid] = std::exp(la_prev + step * (la_new - la_prev));
                    branch_omegas[bid][0] = std::exp(lw_prev + step * (lw_new - lw_prev));
                    branch_weights[bid][0] = 1.0;
                }
                P_branches = compute_all_branch_transition_matrices();
                cur_ll = compute_tree_log_likelihood(P_branches);
                if (cur_ll >= prev_ll) {
                    accepted = true;
                    break;
                }
                step *= 0.5;
            }
            if (!accepted) {
                // Roll back Jacobi step and perform a serial Gauss-Seidel sweep so Phase 2 continues making monotone progress
                for (int32_t bid : branch_ids) {
                    branch_alpha[bid] = prev_alpha[bid];
                    branch_omegas[bid][0] = prev_omega_0[bid];
                    branch_weights[bid][0] = 1.0;
                }
                P_branches = compute_all_branch_transition_matrices();
                cur_ll = prev_ll;
                Matrix V_gs(S, num_patterns);
                Matrix D_gs(S, num_patterns);
                Matrix PD_gs(S, num_patterns);
                for (int32_t bid : branch_ids) {
                    Scalar gs_offset = 0.0;
                    compute_branch_context_matrices(bid, P_branches, V_gs, D_gs, gs_offset);
                    auto gs_obj = [&](Scalar log_alpha, Scalar log_omega) -> Scalar {
                        Scalar ca = std::exp(log_alpha);
                        Scalar cw = std::exp(log_omega);
                        if (ca > 50.0 || cw > 5000.0) return 1e20;
                        MG94Matrix mat = build_q_matrix(cw);
                        Matrix P_b = mat.transition_matrix(ca);
                        return -eval_branch_ll_from_matrices(P_b, V_gs, D_gs, gs_offset, PD_gs);
                    };
                    Scalar init_la = std::log(std::clamp(branch_alpha[bid], 1e-4, 20.0));
                    Scalar init_lw = std::log(std::clamp(branch_omegas[bid][0], 1e-3, 50.0));
                    auto opt_gs = NelderMead2D::minimize(gs_obj, init_la, init_lw, -10.0, 10.0, 1e-4, 60);
                    if (-opt_gs.f >= cur_ll) {
                        branch_alpha[bid] = std::exp(opt_gs.a);
                        branch_omegas[bid][0] = std::exp(opt_gs.b);
                        P_branches[bid] = build_q_matrix(branch_omegas[bid][0]).transition_matrix(branch_alpha[bid]);
                        cur_ll = -opt_gs.f;
                    }
                }
            }

            if (progress_cb) {
                double cycle_frac = static_cast<double>(cycle + 1) / static_cast<double>(max_cycles);
                progress_cb("Phase 2: Baseline Model (cycle " + std::to_string(cycle + 1) + "/" + std::to_string(max_cycles) + ")", 0.20 + 0.20 * cycle_frac);
            }

            if (std::abs(cur_ll - prev_ll) < 1e-3) {
                prev_ll = cur_ll;
                break;
            }
            prev_ll = cur_ll;
        }

        size_t baseline_params = 5 + 9 + 2 * num_branches;
        Scalar baseline_aicc = compute_aicc(prev_ll, baseline_params);

        if (pbar_p2) {
            pbar_p2->set_total(pbar_p2->completed());
            std::ostringstream summary;
            summary << "Baseline fit complete: lnL = " << std::fixed << std::setprecision(2) << prev_ll
                    << " (AICc = " << baseline_aicc << ")";
            pbar_p2->finish(summary.str());
        }

        return {prev_ll, baseline_params, baseline_aicc};
    }

    // Optimize a single branch for M rate classes using local Inside-Outside cache
    // Returns: {best_log_l, best_omegas, best_weights, best_alpha}
    struct LocalBranchFit {
        Scalar log_l = -1e20;
        std::vector<Scalar> omegas;
        std::vector<Scalar> weights;
        Scalar alpha = 0.01;
    };

    LocalBranchFit optimize_branch_mixture_mat(
        int32_t bid,
        int M,
        const Matrix& V_mat,
        const Matrix& D_mat,
        bool constrain_null = false,
        Scalar override_alpha = -1.0,
        const std::vector<Scalar>* override_omegas = nullptr,
        const std::vector<Scalar>* override_weights = nullptr,
        Scalar log_scale_offset = 0.0
    ) const {
        int S = gcode->num_sense_codons;
        size_t num_patterns = aln.patterns.size();
        Scalar cur_alpha = (override_alpha > 0.0) ? override_alpha : std::max(1e-5, branch_alpha[bid]);
        const auto& cur_omegas = (override_omegas != nullptr) ? *override_omegas : branch_omegas[bid];
        const auto& cur_weights = (override_weights != nullptr) ? *override_weights : branch_weights[bid];

        Matrix PD_scratch(S, num_patterns);

        if (M == 1) {
            MG94Matrix mat_w1 = build_q_matrix(1.0);
            auto obj_1d = [&](Scalar log_a, Scalar log_w) -> Scalar {
                Scalar a = std::exp(log_a);
                Matrix P_b = constrain_null ? mat_w1.transition_matrix(a) : build_q_matrix(std::exp(log_w)).transition_matrix(a);
                return -eval_branch_ll_from_matrices(P_b, V_mat, D_mat, log_scale_offset, PD_scratch);
            };

            LocalBranchFit res;
            if (constrain_null) {
                auto obj_a = [&](Scalar log_a) -> Scalar {
                    return obj_1d(log_a, 0.0);
                };
                Scalar scaled_alpha = cur_alpha;
                if (!cur_omegas.empty() && cur_omegas[0] > 1.0) {
                    Scalar cur_exp = syn_scale_factor + cur_omegas[0] * nonsyn_scale_factor;
                    Scalar null_exp = syn_scale_factor + nonsyn_scale_factor;
                    if (null_exp > 1e-12) {
                        scaled_alpha = cur_alpha * (cur_exp / null_exp);
                    }
                }
                Scalar init_log_a = std::log(std::clamp(scaled_alpha, 1e-6, 50.0));
                auto [best_log_a, best_nll] = Brent1D::minimize(obj_a, init_log_a - 6.0, init_log_a, init_log_a + 6.0, 1e-4, 80);
                res.log_l = -best_nll;
                res.omegas = {1.0};
                res.weights = {1.0};
                res.alpha = std::exp(best_log_a);
            } else {
                Scalar init_log_a = std::log(cur_alpha);
                Scalar init_log_w = std::log(std::max(1e-4, cur_omegas.empty() ? 0.5 : cur_omegas[0]));
                auto opt = NelderMead2D::minimize(obj_1d, init_log_a, init_log_w, -12.0, 10.0, 1e-4, 80);
                res.log_l = -opt.f;
                res.omegas = {std::exp(opt.b)};
                res.weights = {1.0};
                res.alpha = std::exp(opt.a);
            }
            return res;
        }

        // M == 2 classes:
        // Parameters:
        // x[0] = log(alpha)
        // x[1] = logit(p1): p1 = 1 / (1 + exp(-x[1]))
        // x[2] = logit(omega1): omega1 = 1 / (1 + exp(-x[2])) in [0, 1]
        // x[3] = log(omega2): omega2 = exp(x[3]) in [0, 10000]
        if (M == 2) {
            if (constrain_null) {
                MG94Matrix mat_w2_one = build_q_matrix(1.0);
                // 3D optimization: log(alpha), logit(p1), logit(omega1) with omega2 = 1.0 fixed
                auto obj_3d = [&](const std::array<Scalar, 3>& x) -> Scalar {
                    Scalar a = std::exp(x[0]);
                    Scalar p1 = 1.0 / (1.0 + std::exp(-x[1]));
                    Scalar p2 = 1.0 - p1;
                    Scalar w1 = 1.0 / (1.0 + std::exp(-x[2]));

                    MG94Matrix mat1 = build_q_matrix(w1);
                    Matrix P_mix = p1 * mat1.transition_matrix(a) + p2 * mat_w2_one.transition_matrix(a);
                    return -eval_branch_ll_from_matrices(P_mix, V_mat, D_mat, log_scale_offset, PD_scratch);
                };

                Scalar init_log_a = std::log(cur_alpha);
                std::array<Scalar, 3> lb = {init_log_a - 5.0, -6.0, -14.0};
                std::array<Scalar, 3> ub = {init_log_a + 6.0,  6.0,   6.0};

                // Warm-start from current branch state if available
                Scalar cur_p1 = cur_weights.empty() ? 0.90 : cur_weights[0];
                Scalar cur_w1 = cur_omegas.empty() ? 0.20 : cur_omegas[0];
                cur_p1 = std::clamp(cur_p1, 0.005, 0.995);
                cur_w1 = std::clamp(cur_w1, 1e-6, 0.995);
                std::array<Scalar, 3> best_pt = {
                    init_log_a,
                    std::log(cur_p1 / (1.0 - cur_p1)),
                    std::log(cur_w1 / (1.0 - cur_w1))
                };
                Scalar best_sc = obj_3d(best_pt);

                // Precompute eigendecompositions for the 6 w1 grid values, and PD = P(a)*D_mat for each (a_mult, w)
                constexpr std::array<Scalar, 6> a_grid = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
                constexpr std::array<Scalar, 5> p1_grid = {0.98, 0.90, 0.75, 0.50, 0.20};
                constexpr std::array<Scalar, 6> w1_grid = {0.0, 0.05, 0.20, 0.40, 0.60, 0.80};

                std::array<MG94Matrix, 6> mat_w1_grid;
                std::array<Scalar, 6> x2_grid;
                for (size_t iw1 = 0; iw1 < 6; ++iw1) {
                    Scalar w1_c = std::clamp(w1_grid[iw1], 1.0 / (1.0 + std::exp(14.0)), 1.0 / (1.0 + std::exp(-6.0)));
                    Scalar x2 = std::log(w1_c / (1.0 - w1_c));
                    x2_grid[iw1] = x2;
                    Scalar w1_exact = 1.0 / (1.0 + std::exp(-x2));
                    mat_w1_grid[iw1] = build_q_matrix(w1_exact);
                }

                // Precompute L_p(a, w) = V_mat.col(p).dot(P(a, w) * D_mat.col(p)) for each grid point!
                // Since Lp(p1) = p1 * L_p(a, w1) + (1 - p1) * L_p(a, 1.0), evaluating all 5 p1 values requires ZERO matrix multiplications!
                Vector L_w2(num_patterns);
                std::array<Vector, 6> L_w1;
                for (size_t iw1 = 0; iw1 < 6; ++iw1) L_w1[iw1].resize(num_patterns);

                for (Scalar a_mult : a_grid) {
                    Scalar cand_log_a = std::log(cur_alpha * a_mult);
                    Scalar a = std::exp(cand_log_a);

                    PD_scratch.noalias() = mat_w2_one.transition_matrix(a) * D_mat;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        L_w2(p) = V_mat.col(p).dot(PD_scratch.col(p));
                    }

                    for (size_t iw1 = 0; iw1 < 6; ++iw1) {
                        PD_scratch.noalias() = mat_w1_grid[iw1].transition_matrix(a) * D_mat;
                        for (size_t p = 0; p < num_patterns; ++p) {
                            L_w1[iw1](p) = V_mat.col(p).dot(PD_scratch.col(p));
                        }
                    }

                    for (Scalar p1_val : p1_grid) {
                        Scalar cand_x1 = std::log(p1_val / (1.0 - p1_val));
                        Scalar p1 = 1.0 / (1.0 + std::exp(-cand_x1));
                        Scalar p2 = 1.0 - p1;
                        for (size_t iw1 = 0; iw1 < 6; ++iw1) {
                            Scalar ll = log_scale_offset;
                            Scalar c_ll = 0.0;
                            for (size_t p = 0; p < num_patterns; ++p) {
                                Scalar Lp = p1 * L_w1[iw1](p) + p2 * L_w2(p);
                                Scalar term = (Lp > 0.0) ? (pattern_weights_vec(p) * std::log(Lp)) : (pattern_weights_vec(p) * (-1e20));
                                LikelihoodEngine::neumaier_acc(ll, c_ll, term);
                            }
                            Scalar sc = -(ll + c_ll);
                            if (sc < best_sc) {
                                best_sc = sc;
                                best_pt = {cand_log_a, cand_x1, x2_grid[iw1]};
                            }
                        }
                    }
                }

                auto opt = NelderMeadND<3>::minimize(obj_3d, best_pt, lb, ub, 1e-4, 120);
                Scalar opt_a = std::exp(opt.x[0]);
                Scalar opt_p1 = 1.0 / (1.0 + std::exp(-opt.x[1]));
                Scalar opt_w1 = 1.0 / (1.0 + std::exp(-opt.x[2]));

                LocalBranchFit res;
                res.log_l = -opt.f;
                res.alpha = opt_a;
                res.omegas = {opt_w1, 1.0};
                res.weights = {opt_p1, 1.0 - opt_p1};
                return res;
            } else {
                // 4D optimization: log(alpha), logit(p1), logit(omega1), log(omega2)
                auto obj_4d = [&](const std::array<Scalar, 4>& x) -> Scalar {
                    Scalar a = std::exp(x[0]);
                    Scalar p1 = 1.0 / (1.0 + std::exp(-x[1]));
                    Scalar p2 = 1.0 - p1;
                    Scalar w1 = 1.0 / (1.0 + std::exp(-x[2]));
                    Scalar w2 = std::exp(x[3]);

                    MG94Matrix mat1 = build_q_matrix(w1);
                    MG94Matrix mat2 = build_q_matrix(w2);
                    Matrix P_mix = p1 * mat1.transition_matrix(a) + p2 * mat2.transition_matrix(a);
                    return -eval_branch_ll_from_matrices(P_mix, V_mat, D_mat, log_scale_offset, PD_scratch);
                };

                Scalar init_log_a = std::log(cur_alpha);
                std::array<Scalar, 4> lb = {init_log_a - 5.0, -6.0, -14.0, -2.3};
                std::array<Scalar, 4> ub = {init_log_a + 8.0,  6.0,   6.0, 10.0};

                // Warm-start from current branch state if available
                Scalar cur_p1 = cur_weights.empty() ? 0.90 : cur_weights[0];
                Scalar cur_w1 = cur_omegas.empty() ? 0.20 : cur_omegas[0];
                Scalar cur_w2 = cur_omegas.size() > 1 ? cur_omegas[1] : 5.0;
                cur_p1 = std::clamp(cur_p1, 0.005, 0.995);
                cur_w1 = std::clamp(cur_w1, 1e-6, 0.995);
                cur_w2 = std::clamp(cur_w2, 0.1, 20000.0);

                std::array<Scalar, 4> best_pt = {
                    init_log_a,
                    std::log(cur_p1 / (1.0 - cur_p1)),
                    std::log(cur_w1 / (1.0 - cur_w1)),
                    std::log(cur_w2)
                };
                Scalar best_sc = obj_4d(best_pt);

                // F-053: Factorized 5 x 5 x 6 x 10 = 1500-point 4D grid search spanning HyPhy's aBSREL.bf:1068-1071 ranges:
                // 5 a_grid values, 5 p1 values ({0.98, 0.95, 0.90, 0.75, 0.50}),
                // 6 w1 values ({0.0, 0.05, 0.20, 0.40, 0.60, 0.80}), and 10 w2 values ({0.25, 0.5, 0.75, 1.0, 1.5, 4.0, 15.0, 80.0, 500.0, 2500.0}).
                constexpr std::array<Scalar, 5> a_grid = {0.25, 0.5, 1.0, 2.0, 6.0};
                constexpr std::array<Scalar, 5> p1_grid = {0.98, 0.95, 0.90, 0.75, 0.50};
                constexpr std::array<Scalar, 6> w1_grid = {0.0, 0.05, 0.20, 0.40, 0.60, 0.80};
                constexpr std::array<Scalar, 10> w2_grid = {0.25, 0.5, 0.75, 1.0, 1.5, 4.0, 15.0, 80.0, 500.0, 2500.0};

                std::array<MG94Matrix, 6> mat_w1_grid;
                std::array<Scalar, 6> x2_grid;
                for (size_t iw1 = 0; iw1 < 6; ++iw1) {
                    Scalar w1_c = std::clamp(w1_grid[iw1], 1.0 / (1.0 + std::exp(14.0)), 1.0 / (1.0 + std::exp(-6.0)));
                    Scalar x2 = std::log(w1_c / (1.0 - w1_c));
                    x2_grid[iw1] = x2;
                    Scalar w1_exact = 1.0 / (1.0 + std::exp(-x2));
                    mat_w1_grid[iw1] = build_q_matrix(w1_exact);
                }

                std::array<MG94Matrix, 10> mat_w2_grid;
                std::array<Scalar, 10> x3_grid;
                for (size_t iw2 = 0; iw2 < 10; ++iw2) {
                    Scalar x3 = std::log(w2_grid[iw2]);
                    x3_grid[iw2] = x3;
                    Scalar w2_exact = std::exp(x3);
                    mat_w2_grid[iw2] = build_q_matrix(w2_exact);
                }

                std::array<Vector, 6> L_w1;
                std::array<Vector, 10> L_w2;
                for (size_t iw1 = 0; iw1 < 6; ++iw1) L_w1[iw1].resize(num_patterns);
                for (size_t iw2 = 0; iw2 < 10; ++iw2) L_w2[iw2].resize(num_patterns);

                for (Scalar a_mult : a_grid) {
                    Scalar cand_log_a = std::log(cur_alpha * a_mult);
                    Scalar a = std::exp(cand_log_a);

                    for (size_t iw1 = 0; iw1 < 6; ++iw1) {
                        PD_scratch.noalias() = mat_w1_grid[iw1].transition_matrix(a) * D_mat;
                        for (size_t p = 0; p < num_patterns; ++p) {
                            L_w1[iw1](p) = V_mat.col(p).dot(PD_scratch.col(p));
                        }
                    }
                    for (size_t iw2 = 0; iw2 < 10; ++iw2) {
                        PD_scratch.noalias() = mat_w2_grid[iw2].transition_matrix(a) * D_mat;
                        for (size_t p = 0; p < num_patterns; ++p) {
                            L_w2[iw2](p) = V_mat.col(p).dot(PD_scratch.col(p));
                        }
                    }

                    for (Scalar p1_val : p1_grid) {
                        Scalar cand_x1 = std::log(p1_val / (1.0 - p1_val));
                        Scalar p1 = 1.0 / (1.0 + std::exp(-cand_x1));
                        Scalar p2 = 1.0 - p1;
                        for (size_t iw1 = 0; iw1 < 6; ++iw1) {
                            for (size_t iw2 = 0; iw2 < 10; ++iw2) {
                                Scalar ll = log_scale_offset;
                                Scalar c_ll = 0.0;
                                for (size_t p = 0; p < num_patterns; ++p) {
                                    Scalar Lp = p1 * L_w1[iw1](p) + p2 * L_w2[iw2](p);
                                    Scalar term = (Lp > 0.0) ? (pattern_weights_vec(p) * std::log(Lp)) : (pattern_weights_vec(p) * (-1e20));
                                    LikelihoodEngine::neumaier_acc(ll, c_ll, term);
                                }
                                Scalar sc = -(ll + c_ll);
                                if (sc < best_sc) {
                                    best_sc = sc;
                                    best_pt = {cand_log_a, cand_x1, x2_grid[iw1], x3_grid[iw2]};
                                }
                            }
                        }
                    }
                }

                auto opt = NelderMeadND<4>::minimize(obj_4d, best_pt, lb, ub, 1e-4, 150);
                Scalar opt_a = std::exp(opt.x[0]);
                Scalar opt_p1 = 1.0 / (1.0 + std::exp(-opt.x[1]));
                Scalar opt_w1 = 1.0 / (1.0 + std::exp(-opt.x[2]));
                Scalar opt_w2 = std::exp(opt.x[3]);

                LocalBranchFit res;
                res.log_l = -opt.f;
                res.alpha = opt_a;
                res.omegas = {opt_w1, opt_w2};
                res.weights = {opt_p1, 1.0 - opt_p1};
                return res;
            }
        }

        // Fallback for M >= 3
        LocalBranchFit res;
        res.log_l = -1e20;
        return res;
    }

    LocalBranchFit optimize_branch_mixture(
        int32_t bid,
        int M,
        const std::vector<Vector>& V_patterns,
        const std::vector<Vector>& D_patterns,
        bool constrain_null = false,
        Scalar override_alpha = -1.0,
        const std::vector<Scalar>* override_omegas = nullptr,
        const std::vector<Scalar>* override_weights = nullptr,
        Scalar log_scale_offset = 0.0
    ) const {
        int S = gcode->num_sense_codons;
        size_t num_patterns = aln.patterns.size();
        Matrix V_mat(S, num_patterns);
        Matrix D_mat(S, num_patterns);
        for (size_t p = 0; p < num_patterns; ++p) {
            V_mat.col(p) = V_patterns[p];
            D_mat.col(p) = D_patterns[p];
        }
        return optimize_branch_mixture_mat(
            bid, M, V_mat, D_mat, constrain_null, override_alpha, override_omegas, override_weights, log_scale_offset
        );
    }

    // Run complete aBSREL analysis
    ABSRELResult run(
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool show_progress = false,
        bool force_progress = false
    ) {
        auto t_start = std::chrono::high_resolution_clock::now();
        ABSRELResult result;
        result.p_threshold = settings.p_threshold;

        // Phase 1: GTR
        if (progress_cb) progress_cb("Phase 1: Fitting Nucleotide GTR Model", 0.05);
        fit_gtr();
        result.gtr_params = gtr_params;
        result.gtr_fit = {gtr_log_l, 39, gtr_aicc};

        // Phase 2: Baseline MG94xREV
        if (progress_cb) progress_cb("Phase 2: Fitting Baseline MG94xREV Model", 0.20);
        auto base_summary = fit_baseline(8, show_progress, force_progress, progress_cb);
        result.baseline_fit = base_summary;

        // Record baseline results per branch
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            ABSRELBranchResult bres;
            bres.branch_name = node.name;
            bres.node_id = node.id;
            bres.baseline_omega = branch_omegas[node.id][0];
            bres.baseline_branch_length = compute_expected_branch_length(node.id);
            bres.rate_classes = 1;
            bres.rate_distribution.rates = branch_omegas[node.id];
            bres.rate_distribution.weights = branch_weights[node.id];
            result.branches[node.name] = bres;
        }

        // Phase 3: Step-up complexity analysis
        if (progress_cb) progress_cb("Phase 3: Branch Complexity Model Selection", 0.40);
        size_t num_patterns = aln.patterns.size();
        int S = gcode->num_sense_codons;

        // Sort branches by baseline length descending (with deterministic tie-breaker by node id)
        std::vector<int32_t> sorted_branches;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                sorted_branches.push_back(node.id);
            }
        }
        std::sort(sorted_branches.begin(), sorted_branches.end(), [&](int32_t a, int32_t b) {
            Scalar la = result.branches[tree.nodes[a].name].baseline_branch_length;
            Scalar lb = result.branches[tree.nodes[b].name].baseline_branch_length;
            if (std::abs(la - lb) > 1e-12) return la > lb;
            return a < b;
        });

        Scalar current_best_aicc = base_summary.aicc;
        size_t current_params = base_summary.parameters;

        std::vector<Matrix> P_branches = compute_all_branch_transition_matrices();

        std::unique_ptr<ProgressBar> pbar_p3;
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar_p3 = std::make_unique<ProgressBar>(
                sorted_branches.size(),
                "[aBSREL] Complexity Selection",
                "branches",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        Matrix V_mat_p3(S, num_patterns);
        Matrix D_mat_p3(S, num_patterns);

        for (size_t b_idx = 0; b_idx < sorted_branches.size(); ++b_idx) {
            int32_t bid = sorted_branches[b_idx];
            const std::string& bname = tree.nodes[bid].name;

            if (pbar_p3) {
                pbar_p3->set_status(bname + ": testing complexity (M=1.." + std::to_string(settings.max_rate_classes) + ")");
            }

            // Compute local (V, D) context matrices for branch bid (zero-alloc, deterministic)
            Scalar log_scale_offset = 0.0;
            compute_branch_context_matrices(bid, P_branches, V_mat_p3, D_mat_p3, log_scale_offset);

            int current_classes = 1;
            LocalBranchFit best_fit;
            best_fit.omegas = branch_omegas[bid];
            best_fit.weights = branch_weights[bid];
            best_fit.alpha = branch_alpha[bid];

            for (int M = 2; M <= settings.max_rate_classes; ++M) {
                if (pbar_p3) {
                    pbar_p3->set_status(bname + ": evaluating M=" + std::to_string(M) + " rate classes...");
                }
                LocalBranchFit cand_fit = optimize_branch_mixture_mat(bid, M, V_mat_p3, D_mat_p3, false, -1.0, nullptr, nullptr, log_scale_offset);
                size_t cand_params = current_params + 2;
                Scalar cand_aicc = compute_aicc(cand_fit.log_l, cand_params);

                if (cand_aicc < current_best_aicc) {
                    current_best_aicc = cand_aicc;
                    current_params = cand_params;
                    current_classes = M;
                    best_fit = cand_fit;
                } else {
                    break;
                }
            }

            if (current_classes > 1) {
                branch_alpha[bid] = best_fit.alpha;
                branch_omegas[bid] = best_fit.omegas;
                branch_weights[bid] = best_fit.weights;
                P_branches[bid] = compute_branch_transition_matrix(bid);
                result.branches[bname].rate_classes = current_classes;
                result.branches[bname].rate_distribution.rates = best_fit.omegas;
                result.branches[bname].rate_distribution.weights = best_fit.weights;
            }

            if (pbar_p3) {
                std::string stat = bname + ": " + std::to_string(current_classes) + " rate class" + (current_classes > 1 ? "es" : "");
                pbar_p3->set_status(stat);
                pbar_p3->tick();
            }

            if (progress_cb) {
                double frac = 0.40 + 0.25 * (static_cast<double>(b_idx + 1) / static_cast<double>(sorted_branches.size()));
                std::string msg = "Phase 3: Branch Complexity (" + bname + " -> " + std::to_string(current_classes) + " rate class" + (current_classes > 1 ? "es" : "") + ")";
                progress_cb(msg, frac);
            }
        }

        if (pbar_p3) {
            size_t multi_rate_count = 0;
            for (const auto& pair : result.branches) {
                if (pair.second.rate_classes > 1) multi_rate_count++;
            }
            std::ostringstream summary;
            summary << "Model selection complete (" << multi_rate_count << " branch(es) with >1 rate class)";
            pbar_p3->finish(summary.str());
        }

        // Phase 4: Full Adaptive Model Fitting
        if (progress_cb) progress_cb("Phase 4: Full Adaptive Model Fitting", 0.65);

        // Joint full model optimization:
        // 1. Joint L-BFGS on all branch lengths using analytical Inside-Outside gradients
        // 2. GTR nucleotide exchangeability rate optimization (Brent 1D on 5 rates)
        // 3. Branch mixture parameter refinement (Nelder-Mead on local mixtures)
        std::vector<int32_t> all_branches;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                all_branches.push_back(node.id);
            }
        }
        size_t n_b = all_branches.size();

        auto run_joint_branch_length_lbfgs_with_params = [&](
            std::vector<Scalar>& cur_alpha,
            const std::vector<std::vector<Scalar>>& cur_omegas,
            const std::vector<std::vector<Scalar>>& cur_weights
        ) -> Scalar {
            // Precompute the MG94Matrix spectral decompositions for each branch's omega components once per L-BFGS call!
            // Inside BranchObjective::operator(), only alpha_ref[bid] changes, so we reuse cached eigendecompositions
            // instead of calling SelfAdjointEigenSolver<Matrix> on every line-search step!
            size_t num_nodes = tree.num_nodes();
            std::vector<std::vector<MG94Matrix>> cached_branch_mg(num_nodes);
            for (int32_t bid : all_branches) {
                size_t K_b = cur_omegas[bid].size();
                cached_branch_mg[bid].resize(K_b);
                for (size_t k = 0; k < K_b; ++k) {
                    cached_branch_mg[bid][k] = build_q_matrix(cur_omegas[bid][k]);
                }
            }

            struct BranchObjective {
                ABSRELAnalyzer& abs;
                const std::vector<int32_t>& nodes;
                int S_val;
                size_t n_pat;
                std::vector<Scalar>& alpha_ref;
                const std::vector<std::vector<Scalar>>& omegas_ref;
                const std::vector<std::vector<Scalar>>& weights_ref;
                const std::vector<std::vector<MG94Matrix>>& branch_mg;
                std::vector<Scalar> pattern_ll;
                Matrix pattern_dL;
                Scalar best_fx = std::numeric_limits<Scalar>::infinity();
                std::vector<Scalar> best_alpha;

                Scalar operator()(const Vector& y, Vector& grad) {
                    size_t n = nodes.size();
                    for (size_t i = 0; i < n; ++i) {
                        alpha_ref[nodes[i]] = std::exp(y(i));
                    }

                    size_t num_nodes = abs.tree.num_nodes();
                    std::vector<Matrix> P_b(num_nodes);
                    std::vector<Matrix> dP_branches(num_nodes);

                    int br_threads = LikelihoodEngine::optimal_omp_threads(n, 4);
                    #pragma omp parallel for num_threads(br_threads) schedule(static)
                    for (size_t i = 0; i < n; ++i) {
                        int32_t bid = nodes[i];
                        Scalar alpha = alpha_ref[bid];
                        P_b[bid] = Matrix::Zero(S_val, S_val);
                        dP_branches[bid] = Matrix::Zero(S_val, S_val);
                        for (size_t k = 0; k < omegas_ref[bid].size(); ++k) {
                            Scalar p = weights_ref[bid][k];
                            const MG94Matrix& mat = branch_mg[bid][k];
                            Matrix P_k = mat.transition_matrix(alpha);
                            P_b[bid].noalias() += p * P_k;
                            dP_branches[bid].noalias() += p * (mat.Q * P_k);
                        }
                    }

                    if (pattern_ll.size() != n_pat) pattern_ll.resize(n_pat);
                    if (pattern_dL.rows() != static_cast<Eigen::Index>(num_nodes) ||
                        pattern_dL.cols() != static_cast<Eigen::Index>(n_pat)) {
                        pattern_dL.resize(num_nodes, n_pat);
                    }
                    std::fill(pattern_ll.begin(), pattern_ll.end(), 0.0);
                    pattern_dL.setZero();

                    int pat_threads = LikelihoodEngine::optimal_omp_threads(n_pat, 16);
                    constexpr size_t B_TILE = 256;
                    size_t num_tiles = (n_pat + B_TILE - 1) / B_TILE;

                    #pragma omp parallel num_threads(pat_threads)
                    {
                        LikelihoodEngine::BatchedInsideOutsideWorkspace ws;
                        ws.init(num_nodes, S_val, std::min(n_pat, B_TILE));
                        Matrix tmp_mat(S_val, std::min(n_pat, B_TILE));
                        Eigen::RowVectorXd inv_L_row(std::min(n_pat, B_TILE));

                        #pragma omp for schedule(static)
                        for (size_t t = 0; t < num_tiles; ++t) {
                            size_t p_start = t * B_TILE;
                            size_t p_end = std::min(p_start + B_TILE, n_pat);
                            size_t B = p_end - p_start;
                            LikelihoodEngine::compute_inside_outside_batch_in_place(
                                abs.tree, abs.aln, p_start, p_end, abs.leaf_to_taxon, P_b, abs.codon_freqs, S_val, ws
                            );
                            for (size_t b = 0; b < B; ++b) {
                                size_t p = p_start + b;
                                Scalar L = ws.likelihood(b);
                                if (L > 0.0) {
                                    pattern_ll[p] = abs.aln.patterns[p].weight * ws.log_likelihood(b);
                                    inv_L_row(b) = abs.aln.patterns[p].weight / L;
                                } else {
                                    inv_L_row(b) = 0.0;
                                }
                            }
                            for (int32_t bid : nodes) {
                                auto dP_D = tmp_mat.leftCols(B);
                                dP_D.noalias() = dP_branches[bid] * ws.D[bid].leftCols(B);
                                pattern_dL.row(bid).segment(p_start, B).array() =
                                    (ws.V[bid].leftCols(B).cwiseProduct(dP_D).colwise().sum()).array() *
                                    inv_L_row.head(B).array();
                            }
                        }
                    }

                    Scalar total_ll = LikelihoodEngine::neumaier_sum(pattern_ll);

                    grad.resize(n);
                    for (size_t i = 0; i < n; ++i) {
                        int32_t bid = nodes[i];
                        Scalar dL_sum = 0.0;
                        Scalar c_dL = 0.0;
                        for (size_t p = 0; p < n_pat; ++p) {
                            LikelihoodEngine::neumaier_acc(dL_sum, c_dL, pattern_dL(bid, p));
                        }
                        Scalar alpha = alpha_ref[bid];
                        grad(i) = -alpha * (dL_sum + c_dL);
                    }

                    Scalar cur_fx = -total_ll;
                    if (std::isfinite(cur_fx) && cur_fx < best_fx) {
                        best_fx = cur_fx;
                        best_alpha = alpha_ref;
                    }

                    return cur_fx;
                }
            };

            BranchObjective b_obj{*this, all_branches, S, num_patterns, cur_alpha, cur_omegas, cur_weights, cached_branch_mg, {}, {}, std::numeric_limits<Scalar>::infinity(), cur_alpha};
            Vector y(n_b), lb(n_b), ub(n_b);
            const Scalar log_max_alpha = std::log(1000.0);
            for (size_t i = 0; i < n_b; ++i) {
                y(i) = std::log(std::clamp(cur_alpha[all_branches[i]], 1e-6, 1000.0));
                lb(i) = -14.0;
                ub(i) = log_max_alpha;
            }

            LBFGSpp::LBFGSBParam<Scalar> opt_param;
            opt_param.m = 6;
            opt_param.epsilon = 1e-4;
            opt_param.max_iterations = 35;

            LBFGSpp::LBFGSBSolver<Scalar> solver(opt_param);
            Scalar fx = 0.0;
            try {
                solver.minimize(b_obj, y, fx, lb, ub);
            } catch (...) {}

            // F-099: Always restore alpha_ref to best_alpha so a rejected line-search trial never corrupts state
            if (!b_obj.best_alpha.empty()) {
                cur_alpha = b_obj.best_alpha;
            }
            if (std::isfinite(b_obj.best_fx)) {
                return -b_obj.best_fx;
            }
            auto P_b = compute_all_branch_transition_matrices_with(cur_alpha, cur_omegas, cur_weights);
            return compute_tree_log_likelihood(P_b);
        };

        auto run_joint_branch_length_lbfgs = [&]() {
            return run_joint_branch_length_lbfgs_with_params(branch_alpha, branch_omegas, branch_weights);
        };

        auto run_gtr_optimization_with_params = [&](
            GTRParameters& gtr_p,
            const std::vector<Scalar>& cur_alpha,
            const std::vector<std::vector<Scalar>>& cur_omegas,
            const std::vector<std::vector<Scalar>>& cur_weights
        ) -> Scalar {
            const GTRParameters saved_member_gtr = gtr_params;
            gtr_params = gtr_p;
            compute_syn_scale();

            auto opt_rate = [&](Scalar& r_ref) {
                auto rate_obj = [&](Scalar val) -> Scalar {
                    Scalar old_val = r_ref;
                    r_ref = val;
                    compute_syn_scale();
                    auto P_b = compute_all_branch_transition_matrices_with(cur_alpha, cur_omegas, cur_weights);
                    Scalar lnl = compute_tree_log_likelihood(P_b);
                    r_ref = old_val;
                    return -lnl;
                };
                Scalar cur_val = r_ref;
                Scalar cur_nll = rate_obj(cur_val);
                auto [best_r, best_nll] = Brent1D::minimize(rate_obj, 1e-4, cur_val, 20.0, 1e-3, 15);
                if (best_nll < cur_nll) {
                    r_ref = best_r;
                } else {
                    r_ref = cur_val;
                }
                compute_syn_scale();
            };

            opt_rate(gtr_params.theta_AC);
            opt_rate(gtr_params.theta_AT);
            opt_rate(gtr_params.theta_CG);
            opt_rate(gtr_params.theta_CT);
            opt_rate(gtr_params.theta_GT);

            gtr_p = gtr_params;
            auto P_b = compute_all_branch_transition_matrices_with(cur_alpha, cur_omegas, cur_weights);
            Scalar lnl = compute_tree_log_likelihood(P_b);
            if (&gtr_p != &gtr_params) {
                gtr_params = saved_member_gtr;
                compute_syn_scale();
            }
            return lnl;
        };

        auto run_gtr_optimization = [&]() {
            return run_gtr_optimization_with_params(gtr_params, branch_alpha, branch_omegas, branch_weights);
        };

        auto run_mixture_refinement_with_params = [&](
            std::vector<Scalar>& cur_alpha,
            std::vector<std::vector<Scalar>>& cur_omegas,
            std::vector<std::vector<Scalar>>& cur_weights,
            int32_t constrained_bid = -1,
            bool update_result_branches = false
        ) -> Scalar {
            auto P_b = compute_all_branch_transition_matrices_with(cur_alpha, cur_omegas, cur_weights);
            Matrix V_mat_ref(S, num_patterns);
            Matrix D_mat_ref(S, num_patterns);
            for (int32_t bid_iter : all_branches) {
                int M = static_cast<int>(cur_omegas[bid_iter].size());
                bool is_constrained = (bid_iter == constrained_bid);
                if (is_constrained && M == 1) {
                    // For M=1 constrained null, omega=1 is fixed and alpha is handled by joint L-BFGS
                    continue;
                }
                Scalar log_scale_offset = 0.0;
                compute_branch_context_matrices(bid_iter, P_b, V_mat_ref, D_mat_ref, log_scale_offset);

                Matrix PD_scratch(S, num_patterns);
                Scalar prev_b_ll = eval_branch_ll_from_matrices(P_b[bid_iter], V_mat_ref, D_mat_ref, log_scale_offset, PD_scratch);

                auto refined_fit = optimize_branch_mixture_mat(
                    bid_iter, M, V_mat_ref, D_mat_ref, is_constrained,
                    cur_alpha[bid_iter], &cur_omegas[bid_iter], &cur_weights[bid_iter], log_scale_offset
                );
                if (refined_fit.log_l >= prev_b_ll) {
                    cur_alpha[bid_iter] = refined_fit.alpha;
                    cur_omegas[bid_iter] = refined_fit.omegas;
                    cur_weights[bid_iter] = refined_fit.weights;
                    P_b[bid_iter] = Matrix::Zero(S, S);
                    for (size_t k = 0; k < cur_omegas[bid_iter].size(); ++k) {
                        MG94Matrix mat = build_q_matrix(cur_omegas[bid_iter][k]);
                        P_b[bid_iter] += cur_weights[bid_iter][k] * mat.transition_matrix(cur_alpha[bid_iter]);
                    }
                    if (update_result_branches) {
                        result.branches[tree.nodes[bid_iter].name].rate_distribution.rates = refined_fit.omegas;
                        result.branches[tree.nodes[bid_iter].name].rate_distribution.weights = refined_fit.weights;
                    }
                }
            }
            return compute_tree_log_likelihood(P_b);
        };

        auto run_mixture_refinement = [&]() {
            return run_mixture_refinement_with_params(branch_alpha, branch_omegas, branch_weights, -1, true);
        };

        Scalar prev_full_ll = compute_tree_log_likelihood(P_branches);
        const int max_joint_passes = 15;
        std::unique_ptr<ProgressBar> pbar_p4;
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar_p4 = std::make_unique<ProgressBar>(
                max_joint_passes,
                "[aBSREL] Full Model Refinement",
                "passes",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        for (int iter = 0; iter < max_joint_passes; ++iter) {
            run_joint_branch_length_lbfgs();
            run_gtr_optimization();
            Scalar cur_full_ll = run_mixture_refinement();

            if (settings.verbose) {
                std::cout << "  [Phase 4 Joint Pass " << iter << "] Log(L) = " << cur_full_ll << " (delta = " << cur_full_ll - prev_full_ll << ")\n";
            }
            if (pbar_p4) {
                std::ostringstream stat;
                stat << "Pass " << (iter + 1) << "/" << max_joint_passes << ": lnL = "
                     << std::fixed << std::setprecision(2) << cur_full_ll;
                pbar_p4->set_status(stat.str());
                pbar_p4->tick();
            }
            if (std::abs(cur_full_ll - prev_full_ll) < 1e-3) {
                break;
            }
            prev_full_ll = cur_full_ll;
        }

        P_branches = compute_all_branch_transition_matrices();
        Scalar full_ll = compute_tree_log_likelihood(P_branches);
        Scalar full_aicc = compute_aicc(full_ll, current_params);
        result.full_adaptive_fit = {full_ll, current_params, full_aicc};

        if (pbar_p4) {
            std::ostringstream summary;
            summary << "Refinement complete: lnL = " << std::fixed << std::setprecision(2) << full_ll
                    << " (AICc = " << full_aicc << ")";
            pbar_p4->finish(summary.str());
        }

        auto update_full_branch_summaries = [&]() {
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int32_t bid_node = node.id;
                const std::string& bname_node = node.name;
                result.branches[bname_node].rate_distribution.rates = branch_omegas[bid_node];
                result.branches[bname_node].rate_distribution.weights = branch_weights[bid_node];
                result.branches[bname_node].full_branch_length = compute_expected_branch_length(bid_node);
                auto [es, en] = compute_es_en(bid_node);
                result.branches[bname_node].full_es = es;
                result.branches[bname_node].full_en = en;
            }
        };
        update_full_branch_summaries();

        // Phase 5: Testing Selected Branches for Positive Selection
        if (progress_cb) progress_cb("Phase 5: Hypothesis Testing for Selection", 0.80);
        std::vector<std::string> to_test;
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            bool test_this = (static_cast<size_t>(node.id) < is_test_branch.size()) ? is_test_branch[node.id] : true;
            if (test_this) {
                to_test.push_back(node.name);
            } else {
                result.branches[node.name].is_tested = false;
            }
        }
        result.tested_branches = to_test;

        std::unique_ptr<ProgressBar> pbar_p5;
        std::atomic<size_t> cand_pos_count{0};
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar_p5 = std::make_unique<ProgressBar>(
                to_test.size(),
                "[aBSREL] Selection Testing",
                "branches",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        // Precompute full-model Inside-Outside context ONCE across the tree for all branches!
        // (Previously this was recomputed from scratch for every tested branch with omega_max > 1).
        size_t num_nodes_total = tree.num_nodes();
        std::vector<Matrix> V_full_all(num_nodes_total, Matrix::Zero(S, num_patterns));
        std::vector<Matrix> D_full_all(num_nodes_total, Matrix::Zero(S, num_patterns));
        auto recompute_full_io_cache = [&]() {
            int io_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 16);
            constexpr size_t B_TILE = 256;
            size_t num_tiles = (num_patterns + B_TILE - 1) / B_TILE;

            #pragma omp parallel num_threads(io_threads)
            {
                LikelihoodEngine::BatchedInsideOutsideWorkspace ws;
                ws.init(num_nodes_total, S, std::min(num_patterns, B_TILE));

                #pragma omp for schedule(static)
                for (size_t t = 0; t < num_tiles; ++t) {
                    size_t p_start = t * B_TILE;
                    size_t p_end = std::min(p_start + B_TILE, num_patterns);
                    size_t B = p_end - p_start;
                    LikelihoodEngine::compute_inside_outside_batch_in_place(
                        tree, aln, p_start, p_end, leaf_to_taxon, P_branches, codon_freqs, S, ws
                    );
                    for (int32_t bid_node : all_branches) {
                        V_full_all[bid_node].middleCols(p_start, B) = ws.V[bid_node].leftCols(B);
                        D_full_all[bid_node].middleCols(p_start, B) = ws.D[bid_node].leftCols(B);
                    }
                }
            }
        };
        recompute_full_io_cache();

        Matrix V_null_mat(S, num_patterns);
        Matrix D_null_mat(S, num_patterns);

        const int max_p5_restarts = std::max<int>(10, static_cast<int>(to_test.size()) + 2);
        for (int p5_restart = 0; p5_restart <= max_p5_restarts; ++p5_restart) {
            bool full_model_improved = false;
            size_t t_idx = 0;
            for (const std::string& bname : to_test) {
                auto& bres = result.branches[bname];
                int32_t bid = bres.node_id;

                if (pbar_p5 && p5_restart == 0) {
                    size_t cpos = cand_pos_count.load(std::memory_order_relaxed);
                    std::string tag = (cpos > 0) ? (" (+" + std::to_string(cpos) + " pos)") : "";
                    pbar_p5->set_status(bname + tag);
                }

                Scalar max_omega = *std::max_element(branch_omegas[bid].begin(), branch_omegas[bid].end());
                if (max_omega <= 1.0) {
                    bres.lrt = 0.0;
                    bres.uncorrected_p_value = 1.0;
                    bres.sites_ebf_100 = 0;
                    if (pbar_p5 && p5_restart == 0) {
                        pbar_p5->tick();
                    }
                    t_idx++;
                    if (progress_cb && p5_restart == 0) {
                        double frac = 0.80 + 0.15 * (static_cast<double>(t_idx) / static_cast<double>(to_test.size()));
                        progress_cb("Phase 5: Hypothesis Testing (" + bname + ": w_max <= 1)", frac);
                    }
                    continue;
                }

                // Full tree optimization under Null Hypothesis (omega_max <= 1.0 on branch bid)
                // Hot start from global best-fit model
                std::vector<Scalar> null_alpha = branch_alpha;
                std::vector<std::vector<Scalar>> null_omegas = branch_omegas;
                std::vector<std::vector<Scalar>> null_weights = branch_weights;
                GTRParameters null_gtr = gtr_params;

                Scalar null_ll = -1e20;

                if (bres.rate_classes == 1) {
                    Scalar full_omega = null_omegas[bid][0];
                    null_omegas[bid] = {1.0};
                    null_weights[bid] = {1.0};

                    // Hot-start alpha to preserve expected substitutions per site
                    Scalar cur_exp = syn_scale_factor + full_omega * nonsyn_scale_factor;
                    Scalar null_exp = syn_scale_factor + nonsyn_scale_factor;
                    if (null_exp > 1e-12) {
                        null_alpha[bid] = std::clamp(branch_alpha[bid] * (cur_exp / null_exp), 1e-6, 1000.0);
                    }

                    // F-017: Multi-pass null optimization for 1-class branches (matching M >= 2)
                    const int max_null_passes = 6;
                    Scalar prev_nll = -1e20;
                    for (int n_pass = 0; n_pass < max_null_passes; ++n_pass) {
                        auto P_null = compute_all_branch_transition_matrices_with(null_alpha, null_omegas, null_weights);
                        Scalar null_scale_offset = 0.0;
                        compute_branch_context_matrices(bid, P_null, V_null_mat, D_null_mat, null_scale_offset);
                        LocalBranchFit fit = optimize_branch_mixture_mat(
                            bid, 1, V_null_mat, D_null_mat, true,
                            null_alpha[bid], &null_omegas[bid], &null_weights[bid], null_scale_offset
                        );
                        null_alpha[bid] = fit.alpha;
                        null_omegas[bid] = fit.omegas;
                        null_weights[bid] = fit.weights;

                        Scalar cur_nll = run_joint_branch_length_lbfgs_with_params(null_alpha, null_omegas, null_weights);
                        if (n_pass > 0 && std::abs(cur_nll - prev_nll) < 1e-3) {
                            null_ll = cur_nll;
                            break;
                        }
                        null_ll = cur_nll;
                        prev_nll = cur_nll;
                    }
                } else {
                    // Rate classes >= 2: constrain omega_max = 1.0 on branch bid
                    Scalar full_omega_max = null_omegas[bid].back();
                    null_omegas[bid].back() = 1.0;
                    if (full_omega_max > 2.0) {
                        Scalar cur_exp = syn_scale_factor + full_omega_max * nonsyn_scale_factor;
                        Scalar null_exp = syn_scale_factor + nonsyn_scale_factor;
                        Scalar p_max = null_weights[bid].back();
                        if (null_exp > 1e-12 && p_max > 0.05) {
                            Scalar scale_ratio = 1.0 + p_max * (cur_exp / null_exp - 1.0);
                            scale_ratio = std::clamp(scale_ratio, 1.0, 10.0);
                            null_alpha[bid] = std::clamp(branch_alpha[bid] * scale_ratio, 1e-6, 1000.0);
                        }
                    }

                    // F-016 & F-017: Up to 6 passes with 1e-3 threshold
                    const int max_null_passes = 6;
                    Scalar prev_nll = -1e20;

                    for (int n_pass = 0; n_pass < max_null_passes; ++n_pass) {
                        auto P_null = compute_all_branch_transition_matrices_with(null_alpha, null_omegas, null_weights);

                        Scalar null_scale_offset = 0.0;
                        compute_branch_context_matrices(bid, P_null, V_null_mat, D_null_mat, null_scale_offset);

                        LocalBranchFit fit = optimize_branch_mixture_mat(
                            bid, bres.rate_classes, V_null_mat, D_null_mat, true,
                            null_alpha[bid], &null_omegas[bid], &null_weights[bid], null_scale_offset
                        );
                        null_alpha[bid] = fit.alpha;
                        null_omegas[bid] = fit.omegas;
                        null_weights[bid] = fit.weights;

                        Scalar cur_nll = run_joint_branch_length_lbfgs_with_params(null_alpha, null_omegas, null_weights);
                        if (n_pass > 0 && std::abs(cur_nll - prev_nll) < 1e-3) {
                            null_ll = cur_nll;
                            break;
                        }
                        null_ll = cur_nll;
                        prev_nll = cur_nll;
                    }
                }

                // F-016 & F-017: Relax GTR rates and non-tested branches' mixtures inside the null loop
                // (holding only bid's top omega <= 1) whenever full_ll > null_ll
                if (std::isfinite(null_ll) && full_ll > null_ll) {
                    const GTRParameters saved_full_gtr = gtr_params;
                    Scalar prev_full_null_ll = null_ll;
                    const int max_full_null_passes = 6;
                    for (int n_pass = 0; n_pass < max_full_null_passes; ++n_pass) {
                        gtr_params = null_gtr;
                        compute_syn_scale();
                        run_gtr_optimization_with_params(null_gtr, null_alpha, null_omegas, null_weights);
                        gtr_params = null_gtr;
                        compute_syn_scale();
                        run_mixture_refinement_with_params(null_alpha, null_omegas, null_weights, bid, false);
                        Scalar cur_nll = run_joint_branch_length_lbfgs_with_params(null_alpha, null_omegas, null_weights);
                        if (cur_nll > null_ll) {
                            null_ll = cur_nll;
                        }
                        if (std::abs(cur_nll - prev_full_null_ll) < 1e-3) {
                            break;
                        }
                        prev_full_null_ll = cur_nll;
                    }
                    gtr_params = saved_full_gtr;
                    compute_syn_scale();
                }

                // F-080: If null_ll > full_ll + 1e-4, warm-start the full model from the null configuration
                // (and restart Phase 5 if a restart round remains) so all branches are tested against a single
                // consistent full model. Since the null configuration is itself a valid point in the unconstrained
                // Full Model space, if the unconstrained polish does not beat null_ll (or on the final round),
                // adopt the null configuration directly; if recomputed new_full_ll still fails to beat full_ll,
                // restore the saved full model state and transition matrices.
                if (p5_restart < max_p5_restarts && std::isfinite(null_ll) && null_ll > -1e19 && null_ll > full_ll + 1e-4) {
                    auto saved_alpha = branch_alpha;
                    auto saved_omegas = branch_omegas;
                    auto saved_weights = branch_weights;
                    GTRParameters saved_gtr = gtr_params;

                    branch_alpha = null_alpha;
                    branch_omegas = null_omegas;
                    branch_weights = null_weights;
                    gtr_params = null_gtr;
                    compute_syn_scale();

                    Scalar prev_ws_ll = null_ll;
                    for (int ws_pass = 0; ws_pass < 15; ++ws_pass) {
                        run_joint_branch_length_lbfgs();
                        run_gtr_optimization();
                        Scalar cur_ws_ll = run_mixture_refinement();
                        if (std::abs(cur_ws_ll - prev_ws_ll) < 1e-3) break;
                        prev_ws_ll = cur_ws_ll;
                    }
                    P_branches = compute_all_branch_transition_matrices();
                    Scalar new_full_ll = compute_tree_log_likelihood(P_branches);
                    if (new_full_ll < null_ll || !std::isfinite(new_full_ll)) {
                        branch_alpha = null_alpha;
                        branch_omegas = null_omegas;
                        branch_weights = null_weights;
                        gtr_params = null_gtr;
                        compute_syn_scale();
                        P_branches = compute_all_branch_transition_matrices();
                        new_full_ll = compute_tree_log_likelihood(P_branches);
                    }
                    if (std::isfinite(new_full_ll) && new_full_ll > full_ll) {
                        full_ll = new_full_ll;
                        full_aicc = compute_aicc(full_ll, current_params);
                        result.full_adaptive_fit = {full_ll, current_params, full_aicc};
                        update_full_branch_summaries();
                        recompute_full_io_cache();
                        full_model_improved = true;
                        break; // Restart Phase 5 so all branches are tested against the improved full model
                    } else {
                        branch_alpha = std::move(saved_alpha);
                        branch_omegas = std::move(saved_omegas);
                        branch_weights = std::move(saved_weights);
                        gtr_params = saved_gtr;
                        compute_syn_scale();
                        P_branches = compute_all_branch_transition_matrices();
                    }
                }

                const Matrix& V_branch_full = V_full_all[bid];
                const Matrix& D_branch_full = D_full_all[bid];

                const bool valid_lls = std::isfinite(null_ll) && std::isfinite(full_ll) &&
                                       null_ll > -1e19 && full_ll > -1e19;
                Scalar lrt = valid_lls ? 2.0 * (full_ll - null_ll) : 0.0;
                if (!valid_lls || lrt <= 0.0 || !std::isfinite(lrt)) {
                    lrt = 0.0;
                    bres.lrt = 0.0;
                    // F-080: When a tested branch had unconstrained omega_K > 1 and lrt == 0.0, set p_value = 0.5 (matching aBSREL.bf:938)
                    bres.uncorrected_p_value = 0.5;
                } else {
                    bres.lrt = lrt;
                    bres.uncorrected_p_value = compute_mixture_p_value(lrt);
                }

                // Compute Empirical Bayes Factors per site
                size_t K = bres.rate_classes;
                Scalar p_pos = bres.rate_distribution.weights.back();
                Scalar prior_odds = (p_pos > 0.0 && p_pos < 1.0) ? (p_pos / (1.0 - p_pos)) : 1.0;

                int sites_ebf_count = 0;
                bres.site_ebf.assign(aln.num_codons, 1.0);

                // Compute pattern EBFs
                std::vector<Scalar> pattern_ebf(num_patterns, 1.0);
                Scalar alpha = branch_alpha[bid];
                std::vector<Matrix> P_comps(K);
                for (size_t k = 0; k < K; ++k) {
                    MG94Matrix mat = build_q_matrix(bres.rate_distribution.rates[k]);
                    P_comps[k] = mat.transition_matrix(alpha);
                }

                for (size_t p = 0; p < num_patterns; ++p) {
                    std::vector<Scalar> comp_L(K);
                    Scalar denom = 0.0;
                    for (size_t k = 0; k < K; ++k) {
                        Scalar L_k = V_branch_full.col(p).dot(P_comps[k] * D_branch_full.col(p));
                        comp_L[k] = bres.rate_distribution.weights[k] * L_k;
                        denom += comp_L[k];
                    }
                    if (denom > 1e-300 && p_pos < 1.0) {
                        Scalar post_pos = comp_L.back() / denom;
                        Scalar one_minus_post = std::max<Scalar>(1e-300, 1.0 - post_pos);
                        Scalar post_odds = post_pos / one_minus_post;
                        pattern_ebf[p] = post_odds / prior_odds;
                    }
                }

                // Map pattern EBF to site EBF
                for (size_t s = 0; s < aln.num_codons; ++s) {
                    size_t p = aln.site_to_pattern[s];
                    Scalar ebf = pattern_ebf[p];
                    bres.site_ebf[s] = ebf;
                    if (ebf >= 100.0) {
                        sites_ebf_count++;
                    }
                }
                bres.sites_ebf_100 = sites_ebf_count;

                if (bres.uncorrected_p_value <= settings.p_threshold) {
                    cand_pos_count.fetch_add(1, std::memory_order_relaxed);
                }
                if (pbar_p5 && p5_restart == 0) {
                    size_t cpos = cand_pos_count.load(std::memory_order_relaxed);
                    if (cpos > 0) {
                        pbar_p5->set_status("\033[1;32m+" + std::to_string(cpos) + "\033[0m candidate positive (" + bname + ")");
                    } else {
                        pbar_p5->set_status(bname);
                    }
                    pbar_p5->tick();
                }

                t_idx++;
                if (progress_cb && p5_restart == 0) {
                    double frac = 0.80 + 0.15 * (static_cast<double>(t_idx) / static_cast<double>(to_test.size()));
                    progress_cb("Phase 5: Hypothesis Testing (" + bname + ": p=" + std::to_string(bres.uncorrected_p_value).substr(0, 6) + ")", frac);
                }
            }
            if (!full_model_improved) break;
        }

        // Phase 6: Holm-Bonferroni correction
        if (progress_cb) progress_cb("Phase 6: Multiple Testing Correction", 0.95);
        size_t m = to_test.size();
        std::vector<Scalar> raw_p_vals(m);
        for (size_t i = 0; i < m; ++i) {
            raw_p_vals[i] = result.branches[to_test[i]].uncorrected_p_value;
        }
        std::vector<Scalar> corr_p_vals = hyphy::core::stats::holm_bonferroni(raw_p_vals);
        for (size_t i = 0; i < m; ++i) {
            const std::string& bname = to_test[i];
            Scalar final_p = corr_p_vals[i];
            result.branches[bname].corrected_p_value = final_p;
            if (final_p <= settings.p_threshold) {
                result.branches[bname].is_positive = true;
                result.positive_branches.push_back(bname);
            }
        }

        if (pbar_p5) {
            std::ostringstream summary;
            summary << "\033[1;32m" << result.positive_branches.size()
                    << " branch(es) under episodic diversifying selection\033[0m (p ≤ "
                    << std::fixed << std::setprecision(2) << settings.p_threshold << ")";
            pbar_p5->finish(summary.str());
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        result.runtime_seconds = std::chrono::duration<double>(t_end - t_start).count();
        if (progress_cb) progress_cb("Analysis Complete", 1.0);

        return result;
    }
};

} // namespace hyphy::analyses
