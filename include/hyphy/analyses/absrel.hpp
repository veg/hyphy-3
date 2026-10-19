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
    bool do_srv = false;
    int syn_rate_classes = 3;
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

    // Scale factors for Q matrices
    Scalar syn_scale_factor = 1.0;

    std::vector<bool> is_test_branch;

    ABSRELAnalyzer(Tree t, Alignment a, ABSRELSettings s = ABSRELSettings{})
        : tree(std::move(t)), aln(std::move(a)), settings(s) {
        gcode = aln.code ? aln.code.get() : GeneticCode::universal().get();
        pos_nuc_freqs = aln.pos_nuc_frequencies;
        codon_freqs = aln.codon_frequencies_f3x4;

        auto sel = tree.select_branches(settings.test_branches, "", {}, false);
        is_test_branch = std::move(sel.is_selected);

        size_t num_nodes = tree.num_nodes();
        leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

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
        MG94Matrix mat;
        MG94Parameters p;
        p.alpha = 1.0;
        p.beta = 0.0; // synonymous only
        p.theta_AC = gtr_params.theta_AC;
        p.theta_AT = gtr_params.theta_AT;
        p.theta_CG = gtr_params.theta_CG;
        p.theta_CT = gtr_params.theta_CT;
        p.theta_GT = gtr_params.theta_GT;
        mat.update(p, pos_nuc_freqs, codon_freqs, *gcode);
        syn_scale_factor = mat.scale_factor;
    }

    // Build rate matrix Q(1, omega)
    MG94Matrix build_q_matrix(Scalar omega) const {
        MG94Matrix mat;
        MG94Parameters p;
        p.alpha = 1.0;
        p.beta = omega;
        p.theta_AC = gtr_params.theta_AC;
        p.theta_AT = gtr_params.theta_AT;
        p.theta_CG = gtr_params.theta_CG;
        p.theta_CT = gtr_params.theta_CT;
        p.theta_GT = gtr_params.theta_GT;
        mat.update(p, pos_nuc_freqs, codon_freqs, *gcode);
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
        size_t num_nodes = tree.num_nodes();
        std::vector<Matrix> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                P_branches[node.id] = compute_branch_transition_matrix(node.id);
            }
        }
        return P_branches;
    }

    std::vector<Matrix> compute_all_branch_transition_matrices_with(
        const std::vector<Scalar>& alphas,
        const std::vector<std::vector<Scalar>>& omegas,
        const std::vector<std::vector<Scalar>>& weights
    ) const {
        size_t num_nodes = tree.num_nodes();
        int S = gcode->num_sense_codons;
        std::vector<Matrix> P_branches(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                int32_t bid = node.id;
                Scalar a = alphas[bid];
                P_branches[bid] = Matrix::Zero(S, S);
                for (size_t k = 0; k < omegas[bid].size(); ++k) {
                    MG94Matrix mat = build_q_matrix(omegas[bid][k]);
                    P_branches[bid] += weights[bid][k] * mat.transition_matrix(a);
                }
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
            MG94Matrix mat = build_q_matrix(omegas[k]);
            total_sub_rate += weights[k] * (mat.scale_factor / 3.0);
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
        int max_cycles = 4,
        bool show_progress = false,
        bool force_progress = false,
        std::function<void(const std::string&, double)> progress_cb = nullptr
    ) {
        int S = gcode->num_sense_codons;
        size_t num_patterns = aln.patterns.size();
        size_t num_branches = tree.num_nodes() - 1;

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

        for (int cycle = 0; cycle < max_cycles; ++cycle) {
            // 1. Single Inside-Outside pass for all patterns across entire tree
            std::vector<std::vector<Vector>> V_all(tree.num_nodes(), std::vector<Vector>(num_patterns));
            std::vector<std::vector<Vector>> D_all(tree.num_nodes(), std::vector<Vector>(num_patterns));

            #pragma omp parallel for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                auto io = LikelihoodEngine::compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, codon_freqs, S
                );
                for (int32_t bid : branch_ids) {
                    V_all[bid][p] = std::move(io.V[bid]);
                    D_all[bid][p] = std::move(io.D[bid]);
                }
            }

            std::vector<Scalar> new_alpha(tree.num_nodes());
            std::vector<Scalar> new_omega(tree.num_nodes());

            // 2. Parallel branch optimization using cached Inside-Outside context
            #pragma omp parallel for schedule(dynamic)
            for (size_t bi = 0; bi < branch_ids.size(); ++bi) {
                int32_t bid = branch_ids[bi];
                const auto& node = tree.nodes[bid];

                auto branch_obj = [&](Scalar log_alpha, Scalar log_omega) -> Scalar {
                    Scalar cur_alpha = std::exp(log_alpha);
                    Scalar cur_omega = std::exp(log_omega);
                    if (cur_alpha > 50.0 || cur_omega > 5000.0) return 1e20;

                    MG94Matrix mat = build_q_matrix(cur_omega);
                    Matrix P_b = mat.transition_matrix(cur_alpha);

                    Scalar ll = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar Lp = V_all[bid][p].dot(P_b * D_all[bid][p]);
                        if (Lp > 0.0) {
                            ll += aln.patterns[p].weight * std::log(Lp);
                        } else {
                            ll += aln.patterns[p].weight * (-1e20);
                        }
                    }
                    return -ll;
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

            for (int32_t bid : branch_ids) {
                branch_alpha[bid] = new_alpha[bid];
                branch_omegas[bid][0] = new_omega[bid];
                branch_weights[bid][0] = 1.0;
                MG94Matrix updated_mat = build_q_matrix(new_omega[bid]);
                P_branches[bid] = updated_mat.transition_matrix(new_alpha[bid]);
            }

            if (progress_cb) {
                double cycle_frac = static_cast<double>(cycle + 1) / static_cast<double>(max_cycles);
                progress_cb("Phase 2: Baseline Model (cycle " + std::to_string(cycle + 1) + "/" + std::to_string(max_cycles) + ")", 0.20 + 0.20 * cycle_frac);
            }

            Scalar cur_ll = compute_tree_log_likelihood(P_branches);
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

    LocalBranchFit optimize_branch_mixture(
        int32_t bid,
        int M,
        const std::vector<Vector>& V_patterns,
        const std::vector<Vector>& D_patterns,
        bool constrain_null = false,
        Scalar override_alpha = -1.0,
        const std::vector<Scalar>* override_omegas = nullptr,
        const std::vector<Scalar>* override_weights = nullptr
    ) const {
        size_t num_patterns = aln.patterns.size();
        Scalar cur_alpha = (override_alpha > 0.0) ? override_alpha : std::max(1e-5, branch_alpha[bid]);
        const auto& cur_omegas = (override_omegas != nullptr) ? *override_omegas : branch_omegas[bid];
        const auto& cur_weights = (override_weights != nullptr) ? *override_weights : branch_weights[bid];

        if (M == 1) {
            auto obj_1d = [&](Scalar log_a, Scalar log_w) -> Scalar {
                Scalar a = std::exp(log_a);
                Scalar w = constrain_null ? 1.0 : std::exp(log_w);
                MG94Matrix mat = build_q_matrix(w);
                Matrix P_b = mat.transition_matrix(a);

                Scalar ll = 0.0;
                for (size_t p = 0; p < num_patterns; ++p) {
                    Scalar Lp = V_patterns[p].dot(P_b * D_patterns[p]);
                    if (Lp > 0.0) {
                        ll += aln.patterns[p].weight * std::log(Lp);
                    } else {
                        ll += aln.patterns[p].weight * (-1e20);
                    }
                }
                return -ll;
            };

            LocalBranchFit res;
            if (constrain_null) {
                auto obj_a = [&](Scalar log_a) -> Scalar {
                    return obj_1d(log_a, 0.0);
                };
                Scalar scaled_alpha = cur_alpha;
                if (!cur_omegas.empty() && cur_omegas[0] > 1.0) {
                    Scalar cur_exp = build_q_matrix(cur_omegas[0]).scale_factor;
                    Scalar null_exp = build_q_matrix(1.0).scale_factor;
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
                // 3D optimization: log(alpha), logit(p1), logit(omega1) with omega2 = 1.0 fixed
                auto obj_3d = [&](const std::array<Scalar, 3>& x) -> Scalar {
                    Scalar a = std::exp(x[0]);
                    Scalar p1 = 1.0 / (1.0 + std::exp(-x[1]));
                    Scalar p2 = 1.0 - p1;
                    Scalar w1 = 1.0 / (1.0 + std::exp(-x[2]));
                    Scalar w2 = 1.0;

                    MG94Matrix mat1 = build_q_matrix(w1);
                    MG94Matrix mat2 = build_q_matrix(w2);
                    Matrix P_mix = p1 * mat1.transition_matrix(a) + p2 * mat2.transition_matrix(a);

                    Scalar ll = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar Lp = V_patterns[p].dot(P_mix * D_patterns[p]);
                        if (Lp > 0.0) {
                            ll += aln.patterns[p].weight * std::log(Lp);
                        } else {
                            ll += aln.patterns[p].weight * (-1e20);
                        }
                    }
                    return -ll;
                };

                Scalar init_log_a = std::log(cur_alpha);
                std::array<Scalar, 3> lb = {init_log_a - 5.0, -6.0, -8.0};
                std::array<Scalar, 3> ub = {init_log_a + 6.0,  6.0,  5.0};

                // Warm-start from current branch state if available
                Scalar cur_p1 = cur_weights.empty() ? 0.90 : cur_weights[0];
                Scalar cur_w1 = cur_omegas.empty() ? 0.20 : cur_omegas[0];
                cur_p1 = std::clamp(cur_p1, 0.01, 0.99);
                cur_w1 = std::clamp(cur_w1, 1e-4, 0.99);
                std::array<Scalar, 3> best_pt = {
                    init_log_a,
                    std::log(cur_p1 / (1.0 - cur_p1)),
                    std::log(cur_w1 / (1.0 - cur_w1))
                };
                Scalar best_sc = obj_3d(best_pt);

                for (Scalar a_mult : {0.25, 0.5, 1.0, 2.0, 4.0, 8.0}) {
                    Scalar cand_log_a = std::log(cur_alpha * a_mult);
                    for (Scalar p1 : {0.98, 0.90, 0.75, 0.50, 0.20}) {
                        Scalar cand_x1 = std::log(p1 / (1.0 - p1));
                        for (Scalar w1 : {0.005, 0.05, 0.20, 0.50, 0.80}) {
                            Scalar cand_x2 = std::log(w1 / (1.0 - w1));
                            std::array<Scalar, 3> pt = {cand_log_a, cand_x1, cand_x2};
                            Scalar sc = obj_3d(pt);
                            if (sc < best_sc) {
                                best_sc = sc;
                                best_pt = pt;
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

                    Scalar ll = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Scalar Lp = V_patterns[p].dot(P_mix * D_patterns[p]);
                        if (Lp > 0.0) {
                            ll += aln.patterns[p].weight * std::log(Lp);
                        } else {
                            ll += aln.patterns[p].weight * (-1e20);
                        }
                    }
                    return -ll;
                };

                Scalar init_log_a = std::log(cur_alpha);
                std::array<Scalar, 4> lb = {init_log_a - 5.0, -6.0, -8.0, -2.3};
                std::array<Scalar, 4> ub = {init_log_a + 8.0,  6.0,  6.0, 10.0};

                // Warm-start from current branch state if available
                Scalar cur_p1 = cur_weights.empty() ? 0.90 : cur_weights[0];
                Scalar cur_w1 = cur_omegas.empty() ? 0.20 : cur_omegas[0];
                Scalar cur_w2 = cur_omegas.size() > 1 ? cur_omegas[1] : 5.0;
                cur_p1 = std::clamp(cur_p1, 0.005, 0.995);
                cur_w1 = std::clamp(cur_w1, 1e-4, 0.995);
                cur_w2 = std::clamp(cur_w2, 0.1, 20000.0);

                std::array<Scalar, 4> best_pt = {
                    init_log_a,
                    std::log(cur_p1 / (1.0 - cur_p1)),
                    std::log(cur_w1 / (1.0 - cur_w1)),
                    std::log(cur_w2)
                };
                Scalar best_sc = obj_4d(best_pt);

                for (Scalar a_mult : {0.5, 1.0, 2.0, 6.0}) {
                    Scalar cand_log_a = std::log(cur_alpha * a_mult);
                    for (Scalar p1 : {0.98, 0.92, 0.80, 0.50}) {
                        Scalar cand_x1 = std::log(p1 / (1.0 - p1));
                        for (Scalar w1 : {0.005, 0.05, 0.20}) {
                            Scalar cand_x2 = std::log(w1 / (1.0 - w1));
                            for (Scalar w2 : {2.0, 15.0, 80.0, 500.0, 2500.0}) {
                                Scalar cand_x3 = std::log(w2);
                                std::array<Scalar, 4> pt = {cand_log_a, cand_x1, cand_x2, cand_x3};
                                Scalar sc = obj_4d(pt);
                                if (sc < best_sc) {
                                    best_sc = sc;
                                    best_pt = pt;
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
        auto base_summary = fit_baseline(4, show_progress, force_progress, progress_cb);
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

        // Sort branches by baseline length descending
        std::vector<int32_t> sorted_branches;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                sorted_branches.push_back(node.id);
            }
        }
        std::sort(sorted_branches.begin(), sorted_branches.end(), [&](int32_t a, int32_t b) {
            return result.branches[tree.nodes[a].name].baseline_branch_length >
                   result.branches[tree.nodes[b].name].baseline_branch_length;
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

        for (size_t b_idx = 0; b_idx < sorted_branches.size(); ++b_idx) {
            int32_t bid = sorted_branches[b_idx];
            const std::string& bname = tree.nodes[bid].name;

            if (pbar_p3) {
                pbar_p3->set_status(bname + ": testing complexity (M=1.." + std::to_string(settings.max_rate_classes) + ")");
            }

            // Compute local (V, D) context for branch bid
            std::vector<Vector> V_patterns(num_patterns);
            std::vector<Vector> D_patterns(num_patterns);
            #pragma omp parallel for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                auto io = LikelihoodEngine::compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, codon_freqs, S
                );
                V_patterns[p] = io.V[bid];
                D_patterns[p] = io.D[bid];
            }

            int current_classes = 1;
            LocalBranchFit best_fit;
            best_fit.omegas = branch_omegas[bid];
            best_fit.weights = branch_weights[bid];
            best_fit.alpha = branch_alpha[bid];

            for (int M = 2; M <= settings.max_rate_classes; ++M) {
                if (pbar_p3) {
                    pbar_p3->set_status(bname + ": evaluating M=" + std::to_string(M) + " rate classes...");
                }
                LocalBranchFit cand_fit = optimize_branch_mixture(bid, M, V_patterns, D_patterns, false);
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
            struct BranchObjective {
                ABSRELAnalyzer& abs;
                const std::vector<int32_t>& nodes;
                int S_val;
                size_t n_pat;
                std::vector<Scalar>& alpha_ref;
                const std::vector<std::vector<Scalar>>& omegas_ref;
                const std::vector<std::vector<Scalar>>& weights_ref;

                Scalar operator()(const Vector& y, Vector& grad) {
                    size_t n = nodes.size();
                    for (size_t i = 0; i < n; ++i) {
                        alpha_ref[nodes[i]] = std::exp(y(i));
                    }

                    size_t num_nodes = abs.tree.num_nodes();
                    std::vector<Matrix> P_b(num_nodes);
                    std::vector<Matrix> dP_branches(num_nodes, Matrix::Zero(S_val, S_val));
                    for (int32_t bid : nodes) {
                        Scalar alpha = alpha_ref[bid];
                        P_b[bid] = Matrix::Zero(S_val, S_val);
                        for (size_t k = 0; k < omegas_ref[bid].size(); ++k) {
                            Scalar w = omegas_ref[bid][k];
                            Scalar p = weights_ref[bid][k];
                            MG94Matrix mat = abs.build_q_matrix(w);
                            Matrix P_k = mat.transition_matrix(alpha);
                            P_b[bid] += p * P_k;
                            dP_branches[bid] += p * (mat.Q * P_k);
                        }
                    }

                    Scalar total_ll = 0.0;
                    std::vector<Scalar> dL_da(num_nodes, 0.0);

                    #pragma omp parallel
                    {
                        Scalar local_ll = 0.0;
                        std::vector<Scalar> local_dL(num_nodes, 0.0);

                        #pragma omp for schedule(dynamic)
                        for (size_t p = 0; p < n_pat; ++p) {
                            const auto& pattern = abs.aln.patterns[p];
                            auto io = LikelihoodEngine::compute_inside_outside(
                                abs.tree, pattern, abs.leaf_to_taxon, P_b, abs.codon_freqs, S_val
                            );
                            Scalar L = io.likelihood;
                            if (L > 0.0) {
                                local_ll += pattern.weight * std::log(L);
                                Scalar inv_L = pattern.weight / L;
                                for (int32_t bid : nodes) {
                                    Scalar term = io.V[bid].dot(dP_branches[bid] * io.D[bid]);
                                    local_dL[bid] += inv_L * term;
                                }
                            }
                        }

                        #pragma omp critical
                        {
                            total_ll += local_ll;
                            for (int32_t bid : nodes) {
                                dL_da[bid] += local_dL[bid];
                            }
                        }
                    }

                    grad.resize(n);
                    for (size_t i = 0; i < n; ++i) {
                        int32_t bid = nodes[i];
                        Scalar alpha = alpha_ref[bid];
                        grad(i) = -alpha * dL_da[bid];
                    }

                    return -total_ll;
                }
            };

            BranchObjective b_obj{*this, all_branches, S, num_patterns, cur_alpha, cur_omegas, cur_weights};
            Vector y(n_b), lb(n_b), ub(n_b);
            for (size_t i = 0; i < n_b; ++i) {
                y(i) = std::log(std::clamp(cur_alpha[all_branches[i]], 1e-6, 50.0));
                lb(i) = -14.0;
                ub(i) =   5.0;
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
            return -fx;
        };

        auto run_joint_branch_length_lbfgs = [&]() {
            return run_joint_branch_length_lbfgs_with_params(branch_alpha, branch_omegas, branch_weights);
        };

        auto run_gtr_optimization = [&]() {
            auto opt_rate = [&](Scalar& r_ref) {
                auto rate_obj = [&](Scalar val) -> Scalar {
                    Scalar old_val = r_ref;
                    r_ref = val;
                    compute_syn_scale();
                    auto P_b = compute_all_branch_transition_matrices();
                    Scalar lnl = compute_tree_log_likelihood(P_b);
                    r_ref = old_val;
                    return -lnl;
                };
                auto [best_r, _] = Brent1D::minimize(rate_obj, 1e-4, r_ref, 20.0, 1e-3, 15);
                r_ref = best_r;
                compute_syn_scale();
            };

            opt_rate(gtr_params.theta_AC);
            opt_rate(gtr_params.theta_AT);
            opt_rate(gtr_params.theta_CG);
            opt_rate(gtr_params.theta_CT);
            opt_rate(gtr_params.theta_GT);

            auto P_b = compute_all_branch_transition_matrices();
            return compute_tree_log_likelihood(P_b);
        };

        auto run_mixture_refinement = [&]() {
            auto P_b = compute_all_branch_transition_matrices();
            for (int32_t bid : all_branches) {
                int M = branch_omegas[bid].size();
                std::vector<Vector> V_patterns(num_patterns);
                std::vector<Vector> D_patterns(num_patterns);
                #pragma omp parallel for schedule(dynamic)
                for (size_t p = 0; p < num_patterns; ++p) {
                    const auto& pattern = aln.patterns[p];
                    auto io = LikelihoodEngine::compute_inside_outside(
                        tree, pattern, leaf_to_taxon, P_b, codon_freqs, S
                    );
                    V_patterns[p] = io.V[bid];
                    D_patterns[p] = io.D[bid];
                }

                auto refined_fit = optimize_branch_mixture(bid, M, V_patterns, D_patterns, false);
                branch_alpha[bid] = refined_fit.alpha;
                branch_omegas[bid] = refined_fit.omegas;
                branch_weights[bid] = refined_fit.weights;
                P_b[bid] = compute_branch_transition_matrix(bid);

                result.branches[tree.nodes[bid].name].rate_distribution.rates = refined_fit.omegas;
                result.branches[tree.nodes[bid].name].rate_distribution.weights = refined_fit.weights;
            }
            return compute_tree_log_likelihood(P_b);
        };

        Scalar prev_full_ll = compute_tree_log_likelihood(P_branches);
        const int max_joint_passes = 6;
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
            if (std::abs(cur_full_ll - prev_full_ll) < 0.05) {
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

        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            int32_t bid = node.id;
            const std::string& bname = node.name;
            result.branches[bname].full_branch_length = compute_expected_branch_length(bid);
            auto [es, en] = compute_es_en(bid);
            result.branches[bname].full_es = es;
            result.branches[bname].full_en = en;
        }

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

        size_t t_idx = 0;
        for (const std::string& bname : to_test) {
            auto& bres = result.branches[bname];
            int32_t bid = bres.node_id;

            if (pbar_p5) {
                size_t cpos = cand_pos_count.load(std::memory_order_relaxed);
                std::string tag = (cpos > 0) ? (" (+" + std::to_string(cpos) + " pos)") : "";
                pbar_p5->set_status(bname + tag);
            }

            Scalar max_omega = *std::max_element(branch_omegas[bid].begin(), branch_omegas[bid].end());
            if (max_omega <= 1.0) {
                bres.lrt = 0.0;
                bres.uncorrected_p_value = 1.0;
                bres.sites_ebf_100 = 0;
                if (pbar_p5) {
                    pbar_p5->tick();
                }
                t_idx++;
                if (progress_cb) {
                    double frac = 0.80 + 0.15 * (static_cast<double>(t_idx) / static_cast<double>(to_test.size()));
                    progress_cb("Phase 5: Hypothesis Testing (" + bname + ": w_max <= 1)", frac);
                }
                continue;
            }

            // Compute local (V, D) context under Full model (used for Empirical Bayes Factors)
            std::vector<Vector> V_patterns(num_patterns);
            std::vector<Vector> D_patterns(num_patterns);
            #pragma omp parallel for schedule(dynamic)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                auto io = LikelihoodEngine::compute_inside_outside(
                    tree, pattern, leaf_to_taxon, P_branches, codon_freqs, S
                );
                V_patterns[p] = io.V[bid];
                D_patterns[p] = io.D[bid];
            }

            // Full tree optimization under Null Hypothesis (omega_max <= 1.0 on branch bid)
            // Hot start from global best-fit model
            std::vector<Scalar> null_alpha = branch_alpha;
            std::vector<std::vector<Scalar>> null_omegas = branch_omegas;
            std::vector<std::vector<Scalar>> null_weights = branch_weights;

            Scalar null_ll = -1e20;

            if (bres.rate_classes == 1) {
                Scalar full_omega = null_omegas[bid][0];
                null_omegas[bid] = {1.0};
                null_weights[bid] = {1.0};

                // Hot-start alpha to preserve expected substitutions per site
                Scalar cur_exp = build_q_matrix(full_omega).scale_factor;
                Scalar null_exp = build_q_matrix(1.0).scale_factor;
                if (null_exp > 1e-12) {
                    null_alpha[bid] = std::clamp(branch_alpha[bid] * (cur_exp / null_exp), 1e-6, 50.0);
                }

                // Full tree optimization using Fréchet derivative quasi-Newton (L-BFGS-B)
                null_ll = run_joint_branch_length_lbfgs_with_params(null_alpha, null_omegas, null_weights);
            } else {
                // Rate classes >= 2: constrain omega_max = 1.0 on branch bid
                Scalar full_omega_max = null_omegas[bid].back();
                null_omegas[bid].back() = 1.0;
                if (full_omega_max > 2.0) {
                    Scalar cur_exp = build_q_matrix(full_omega_max).scale_factor;
                    Scalar null_exp = build_q_matrix(1.0).scale_factor;
                    Scalar p_max = null_weights[bid].back();
                    if (null_exp > 1e-12 && p_max > 0.05) {
                        Scalar scale_ratio = 1.0 + p_max * (cur_exp / null_exp - 1.0);
                        scale_ratio = std::clamp(scale_ratio, 1.0, 10.0);
                        null_alpha[bid] = std::clamp(branch_alpha[bid] * scale_ratio, 1e-6, 50.0);
                    }
                }

                // Full tree optimization: alternating mixture refinement on branch bid
                // and Fréchet quasi-Newton branch length optimization on the entire tree
                const int max_null_passes = 3;
                Scalar prev_nll = -1e20;

                for (int n_pass = 0; n_pass < max_null_passes; ++n_pass) {
                    auto P_null = compute_all_branch_transition_matrices_with(null_alpha, null_omegas, null_weights);

                    std::vector<Vector> V_null(num_patterns);
                    std::vector<Vector> D_null(num_patterns);
                    #pragma omp parallel for schedule(dynamic)
                    for (size_t p = 0; p < num_patterns; ++p) {
                        const auto& pattern = aln.patterns[p];
                        auto io = LikelihoodEngine::compute_inside_outside(
                            tree, pattern, leaf_to_taxon, P_null, codon_freqs, S
                        );
                        V_null[p] = io.V[bid];
                        D_null[p] = io.D[bid];
                    }

                    LocalBranchFit fit = optimize_branch_mixture(
                        bid, bres.rate_classes, V_null, D_null, true,
                        null_alpha[bid], &null_omegas[bid], &null_weights[bid]
                    );
                    null_alpha[bid] = fit.alpha;
                    null_omegas[bid] = fit.omegas;
                    null_weights[bid] = fit.weights;

                    Scalar cur_nll = run_joint_branch_length_lbfgs_with_params(null_alpha, null_omegas, null_weights);
                    if (n_pass > 0 && std::abs(cur_nll - prev_nll) < 0.02) {
                        null_ll = cur_nll;
                        break;
                    }
                    null_ll = cur_nll;
                    prev_nll = cur_nll;
                }
            }

            Scalar lrt = 0.0;
            if (null_ll <= -1e19 || full_ll <= -1e19 || !std::isfinite(null_ll) || !std::isfinite(full_ll)) {
                lrt = 0.0;
                bres.lrt = 0.0;
                bres.uncorrected_p_value = 1.0;
            } else {
                lrt = 2.0 * (full_ll - null_ll);
                if (lrt < 0.0) lrt = 0.0;
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
                    Scalar L_k = V_patterns[p].dot(P_comps[k] * D_patterns[p]);
                    comp_L[k] = bres.rate_distribution.weights[k] * L_k;
                    denom += comp_L[k];
                }
                if (denom > 1e-300 && p_pos < 1.0) {
                    Scalar post_pos = comp_L.back() / denom;
                    if (post_pos < 1.0) {
                        Scalar post_odds = post_pos / (1.0 - post_pos);
                        pattern_ebf[p] = post_odds / prior_odds;
                    } else {
                        pattern_ebf[p] = 1e6;
                    }
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
            if (pbar_p5) {
                size_t cpos = cand_pos_count.load(std::memory_order_relaxed);
                if (cpos > 0) {
                    pbar_p5->set_status("\033[1;32m+" + std::to_string(cpos) + "\033[0m candidate positive (" + bname + ")");
                } else {
                    pbar_p5->set_status(bname);
                }
                pbar_p5->tick();
            }

            t_idx++;
            if (progress_cb) {
                double frac = 0.80 + 0.15 * (static_cast<double>(t_idx) / static_cast<double>(to_test.size()));
                progress_cb("Phase 5: Hypothesis Testing (" + bname + ": p=" + std::to_string(bres.uncorrected_p_value).substr(0, 6) + ")", frac);
            }
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
