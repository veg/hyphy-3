#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/analyses/busted_types.hpp"
#include "hyphy/analyses/mss.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/squarem.hpp"
#include "nlohmann/json.hpp"
#include "hyphy/core/provenance.hpp"

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

class BUSTEDAnalyzer {
public:
    Tree tree;
    Alignment aln;
    std::shared_ptr<const GeneticCode> code;
    mutable MG94Parameters base_params;
    mutable std::vector<Scalar> syn_pair_rates_storage; // Optional S*S flat array for BUSTED-MSS

    // Relative branch lengths tau_b (normalized from MG94 baseline)
    mutable std::vector<Scalar> branch_tau;
    mutable std::vector<std::array<Scalar, 3>> bridge_tau;
    mutable std::vector<bool> is_test_branch;
    // 0 = Test, 1 = Background/Comparison, 2 = Unclassified/Nuisance
    mutable std::vector<int> branch_partition;
    std::vector<size_t> leaf_to_taxon;
    mutable size_t num_test_branches = 0;
    mutable size_t num_background_branches = 0;
    mutable size_t num_unclassified_branches = 0;

    // Precomputed off-diagonal saturation replacement kernel P_sink(i, j) = pi_j / (1 - pi_i) for i != j, 0 for i == j
    Matrix P_sink;
    Vector inv_one_minus_pi;

    // GTR baseline results
    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    // Global MG94 baseline results
    Scalar mg94_log_l = 0.0;
    Scalar mg94_omega = 1.0;
    Scalar mg94_background_omega = -1.0;
    Scalar mg94_tree_length = 0.0;
    mutable size_t num_mss_free_params = 0;

    BUSTEDAnalyzer(
        Tree t,
        Alignment a,
        MG94Parameters base_p = MG94Parameters{}
    ) : tree(std::move(t)), aln(std::move(a)), base_params(base_p) {
        code = aln.code ? aln.code : GeneticCode::universal();
        size_t num_nodes = tree.num_nodes();
        branch_tau.assign(num_nodes, 0.0);
        is_test_branch.assign(num_nodes, true); // Default: all branches are Test
        branch_partition.assign(num_nodes, 0);
        if (tree.root_id >= 0 && static_cast<size_t>(tree.root_id) < num_nodes) {
            is_test_branch[tree.root_id] = false;
        }
        num_test_branches = (num_nodes > 0) ? (num_nodes - 1) : 0; // Non-root branches
        num_background_branches = 0;
        num_unclassified_branches = 0;
        leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);
        build_sink_kernel();
        prepare_baseline_lengths();
    }

    void select_branch_partitions(
        const std::string& test_spec,
        const std::string& comp_spec = "",
        bool auto_detect_tags = false
    ) const {
        size_t num_nodes = tree.num_nodes();
        auto sel_test = tree.select_branches(test_spec, "", {}, auto_detect_tags);
        is_test_branch = sel_test.is_selected;
        num_test_branches = sel_test.num_selected;
        branch_partition.assign(num_nodes, 1);
        num_background_branches = 0;
        num_unclassified_branches = 0;

        auto to_lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return s;
        };
        auto trim_braces = [](std::string s) {
            while (!s.empty() && (std::isspace(static_cast<unsigned char>(s.front())) || s.front() == '{')) s.erase(s.begin());
            while (!s.empty() && (std::isspace(static_cast<unsigned char>(s.back())) || s.back() == '}')) s.pop_back();
            return s;
        };
        auto is_bg_synonym = [](const std::string& tl) {
            return tl == "r" || tl == "ref" || tl == "reference" || tl == "bg" || tl == "background";
        };

        bool has_tagged_bg = false;
        if (comp_spec.empty() && auto_detect_tags) {
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && !is_test_branch[node.id] && !node.model_tag.empty()) {
                    if (is_bg_synonym(to_lower(trim_braces(node.model_tag)))) {
                        has_tagged_bg = true;
                        break;
                    }
                }
            }
        }

        if (!comp_spec.empty()) {
            auto sel_comp = tree.select_branches(comp_spec, "", {}, false);
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int nid = node.id;
                if (is_test_branch[nid]) {
                    branch_partition[nid] = 0;
                } else if (sel_comp.is_selected[nid]) {
                    branch_partition[nid] = 1;
                    num_background_branches++;
                } else {
                    branch_partition[nid] = 2;
                    num_unclassified_branches++;
                }
            }
            if (num_background_branches == 0) {
                throw std::runtime_error("No non-test branches matched '--comparison " + comp_spec + "'.");
            }
        } else if (has_tagged_bg) {
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int nid = node.id;
                if (is_test_branch[nid]) {
                    branch_partition[nid] = 0;
                } else if (is_bg_synonym(to_lower(trim_braces(node.model_tag)))) {
                    branch_partition[nid] = 1;
                    num_background_branches++;
                } else {
                    branch_partition[nid] = 2;
                    num_unclassified_branches++;
                }
            }
        } else {
            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int nid = node.id;
                if (is_test_branch[nid]) {
                    branch_partition[nid] = 0;
                } else {
                    branch_partition[nid] = 1;
                    num_background_branches++;
                }
            }
        }
    }

    void select_test_branches(const std::string& spec) const {
        select_branch_partitions(spec, "", false);
    }

    void build_sink_kernel() {
        int S = code->num_sense_codons;
        P_sink = Matrix::Zero(S, S);
        inv_one_minus_pi = Vector::Zero(S);
        for (int i = 0; i < S; ++i) {
            Scalar denom = std::max(1e-12, 1.0 - aln.codon_frequencies_f3x4(i));
            inv_one_minus_pi(i) = 1.0 / denom;
            for (int j = 0; j < S; ++j) {
                if (i != j) {
                    P_sink(i, j) = aln.codon_frequencies_f3x4(j) * inv_one_minus_pi(i);
                }
            }
        }
    }

    // O(S) evaluation of V^T * P_sink * D using the rank-1 off-diagonal identity
    inline Scalar eval_sink_bilinear(const Vector& V, const Vector& D) const {
        Scalar pi_dot_D = aln.codon_frequencies_f3x4.dot(D);
        int S = static_cast<int>(V.size());
        Scalar acc = 0.0;
        for (int i = 0; i < S; ++i) {
            Scalar p_sink_D_i = (pi_dot_D - aln.codon_frequencies_f3x4(i) * D(i)) * inv_one_minus_pi(i);
            acc += V(i) * std::max(0.0, p_sink_D_i);
        }
        return acc;
    }

    static BUSTEDAnalyzer create_and_fit(
        Tree raw_tree,
        Alignment aln,
        const std::vector<bool>* is_test_branch_ptr = nullptr
    ) {
        GTRFitter gtr_fitter(raw_tree, aln);
        auto gtr_res = gtr_fitter.fit();

        MG94Parameters base_p;
        base_p.theta_AC = gtr_res.params.theta_AC;
        base_p.theta_AT = gtr_res.params.theta_AT;
        base_p.theta_CG = gtr_res.params.theta_CG;
        base_p.theta_CT = gtr_res.params.theta_CT;
        base_p.theta_GT = gtr_res.params.theta_GT;

        MG94Fitter mg_fitter(gtr_res.tree, aln);
        Scalar bg_omega = -1.0;
        auto mg_res = mg_fitter.fit_full_model(
            1.0, base_p, nullptr, 100, 1e-4,
            is_test_branch_ptr, is_test_branch_ptr ? &bg_omega : nullptr
        );
        base_p = mg_res.params;

        BUSTEDAnalyzer analyzer(mg_res.tree, std::move(aln), base_p);
        analyzer.has_gtr_fit = true;
        analyzer.gtr_log_l = gtr_res.log_likelihood;
        analyzer.gtr_aicc = gtr_res.aicc;
        analyzer.gtr_rates = gtr_res.params;
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                analyzer.gtr_branch_lengths[node.name] = node.branch_length;
            }
        }

        analyzer.mg94_log_l = mg_res.log_likelihood;
        analyzer.mg94_omega = mg_res.params.beta;
        analyzer.mg94_background_omega = bg_omega;
        return analyzer;
    }

    void prepare_baseline_lengths() {
        MG94Matrix mat;
        mat.update(base_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        Scalar Q_scale = mat.scale_factor;
        Scalar conversion_factor = (Q_scale > 1e-12) ? (3.0 / Q_scale) : 1.0;

        mg94_tree_length = 0.0;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar scaled_len = node.branch_length * conversion_factor;
                branch_tau[node.id] = scaled_len;
                mg94_tree_length += scaled_len;
            }
        }
        bridge_tau.resize(tree.bridges.size());
        for (size_t k = 0; k < tree.bridges.size(); ++k) {
            for (int r = 0; r < 3; ++r) {
                bridge_tau[k][r] = tree.bridges[k].branch_lengths[r] * conversion_factor;
            }
        }
    }

    Scalar bridge_resolution_tau(
        size_t b_idx,
        int r,
        Scalar s,
        const std::vector<Scalar>* custom_branch_lengths = nullptr
    ) const {
        Scalar base_tau = (b_idx < bridge_tau.size()) ? bridge_tau[b_idx][r] : tree.bridges[b_idx].branch_lengths[r];
        if (custom_branch_lengths != nullptr && b_idx < tree.bridges.size()) {
            const auto& wb = tree.bridges[b_idx];
            int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
            if (bid >= 0 && static_cast<size_t>(bid) < custom_branch_lengths->size() &&
                static_cast<size_t>(bid) < branch_tau.size() && branch_tau[bid] > 1e-12) {
                base_tau *= ((*custom_branch_lengths)[bid] / branch_tau[bid]);
            }
        }
        return s * base_tau;
    }

    std::vector<std::array<Matrix, 3>> build_one_set_bridge_mix(
        Scalar s,
        Scalar alpha_m,
        const std::vector<MG94Matrix>& mg,
        const MG94Matrix& mg_bg,
        const std::vector<Scalar>& weights,
        Scalar eps_sink,
        Scalar bio_scale,
        const std::vector<Scalar>* custom_branch_lengths = nullptr
    ) const {
        int S = code->num_sense_codons;
        size_t K = mg.size();
        std::vector<std::array<Matrix, 3>> P_bridges(tree.bridges.size());
        for (size_t b_idx = 0; b_idx < tree.bridges.size(); ++b_idx) {
            const auto& wb = tree.bridges[b_idx];
            int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
            bool is_test = (bid >= 0 && static_cast<size_t>(bid) < is_test_branch.size()) ? is_test_branch[bid] : true;
            for (int r = 0; r < 3; ++r) {
                if (wb.weights[r] > 0.0 || r == 0) {
                    Scalar bl = bridge_resolution_tau(b_idx, r, s, custom_branch_lengths);
                    if (is_test) {
                        P_bridges[b_idx][r] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink > 0.0) ? bio_scale : 1.0;
                        for (size_t k = 0; k < K; ++k) {
                            P_bridges[b_idx][r] += (use_bio * weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink > 0.0) {
                            P_bridges[b_idx][r].noalias() += eps_sink * P_sink;
                        }
                    } else {
                        P_bridges[b_idx][r] = mg_bg.transition_matrix(alpha_m * bl);
                    }
                } else {
                    P_bridges[b_idx][r] = Matrix::Identity(S, S);
                }
            }
        }
        return P_bridges;
    }

    std::vector<std::array<Matrix, 3>> build_two_set_bridge_mix(
        Scalar s,
        Scalar alpha_m,
        const std::vector<MG94Matrix>& mg_t,
        const std::vector<MG94Matrix>& mg_b,
        const std::vector<MG94Matrix>& mg_u,
        const MG94Matrix& mg_unclass,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>* weights_unclass,
        bool has_unclass_dist,
        Scalar eps_sink_test,
        Scalar eps_sink_bg,
        Scalar bio_scale_t,
        Scalar bio_scale_b,
        const std::vector<Scalar>* custom_branch_lengths = nullptr
    ) const {
        int S = code->num_sense_codons;
        size_t K_t = mg_t.size();
        size_t K_b = mg_b.size();
        size_t K_u = mg_u.size();
        std::vector<std::array<Matrix, 3>> P_bridges(tree.bridges.size());
        for (size_t b_idx = 0; b_idx < tree.bridges.size(); ++b_idx) {
            const auto& wb = tree.bridges[b_idx];
            int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
            int part = (bid >= 0 && static_cast<size_t>(bid) < branch_partition.size())
                ? branch_partition[bid]
                : ((bid >= 0 && static_cast<size_t>(bid) < is_test_branch.size() && is_test_branch[bid]) ? 0 : 1);
            for (int r = 0; r < 3; ++r) {
                if (wb.weights[r] > 0.0 || r == 0) {
                    Scalar bl = bridge_resolution_tau(b_idx, r, s, custom_branch_lengths);
                    if (part == 0) {
                        P_bridges[b_idx][r] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_test > 0.0) ? bio_scale_t : 1.0;
                        for (size_t k = 0; k < K_t; ++k) {
                            P_bridges[b_idx][r] += (use_bio * weights_test[k]) * mg_t[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_test > 0.0) {
                            P_bridges[b_idx][r].noalias() += eps_sink_test * P_sink;
                        }
                    } else if (part == 1) {
                        P_bridges[b_idx][r] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_bg > 0.0) ? bio_scale_b : 1.0;
                        for (size_t k = 0; k < K_b; ++k) {
                            P_bridges[b_idx][r] += (use_bio * weights_bg[k]) * mg_b[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_bg > 0.0) {
                            P_bridges[b_idx][r].noalias() += eps_sink_bg * P_sink;
                        }
                    } else if (has_unclass_dist && weights_unclass != nullptr) {
                        P_bridges[b_idx][r] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_u; ++k) {
                            P_bridges[b_idx][r] += (*weights_unclass)[k] * mg_u[k].transition_matrix(alpha_m * bl);
                        }
                    } else {
                        P_bridges[b_idx][r] = mg_unclass.transition_matrix(alpha_m * bl);
                    }
                } else {
                    P_bridges[b_idx][r] = Matrix::Identity(S, S);
                }
            }
        }
        return P_bridges;
    }

    // Evaluates tree log-likelihood for arbitrary K rate categories, optional SRV (M classes), optional saturation error sink (eps), and custom branch lengths
    Scalar evaluate_log_l(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        std::vector<Scalar>* out_pattern_ll = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        size_t M = syn_rates.size();
        Scalar bio_scale = 1.0 - eps_sink;
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;

        // Build K rate matrices
        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        bool has_bg = (num_test_branches + 1 < num_nodes);
        MG94Matrix mg_bg;
        if (has_bg) {
            MG94Parameters p_bg = active_base;
            p_bg.alpha = 1.0;
            p_bg.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega;
            p_bg.delta = delta;
            p_bg.psi = psi;
            mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        // Precompute branch transition matrices P_mix[m][node_id]
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[node.id] : branch_tau[node.id]);
                    if (is_test_branch[node.id]) {
                        P_mix[m][node.id] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink > 0.0) ? bio_scale : 1.0;
                        for (size_t k = 0; k < K; ++k) {
                            P_mix[m][node.id] += (use_bio * weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink > 0.0) {
                            P_mix[m][node.id].noalias() += eps_sink * P_sink;
                        }
                    } else {
                        P_mix[m][node.id] = mg_bg.transition_matrix(alpha_m * bl);
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_one_set_bridge_mix(s, alpha_m, mg, mg_bg, weights, eps_sink, bio_scale, custom_branch_lengths);
            }
        }

        if (out_pattern_ll) {
            out_pattern_ll->assign(num_patterns, 0.0);
        }

        std::vector<Scalar> pat_ll_vec(num_patterns, -1e20);
        int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            LikelihoodEngine::BatchedInsideWorkspace ws;
            ws.init(num_nodes, S, tile_size);
            Eigen::MatrixXd terms_mat(M, tile_size);

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                if (M == 1) {
                    LikelihoodEngine::compute_inside_batch_in_place(
                        tree, aln, p_start, p_end, leaf_to_taxon, P_mix[0], aln.codon_frequencies_f3x4, S, ws,
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[0]
                    );
                    for (Eigen::Index b = 0; b < B; ++b) {
                        Scalar pll = ws.root_log_L(b);
                        pat_ll_vec[p_start + b] = pll;
                        if (out_pattern_ll) (*out_pattern_ll)[p_start + b] = pll;
                    }
                } else {
                    terms_mat.leftCols(B).setConstant(-1e20);
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] <= 1e-12) continue;
                        Scalar log_sw = std::log(syn_weights[m]);
                        LikelihoodEngine::compute_inside_batch_in_place(
                            tree, aln, p_start, p_end, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, ws,
                            tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                        );
                        for (Eigen::Index b = 0; b < B; ++b) {
                            if (ws.root_L(b) > 0.0) {
                                terms_mat(m, b) = log_sw + ws.root_log_L(b);
                            }
                        }
                    }
                    for (Eigen::Index b = 0; b < B; ++b) {
                        Scalar max_term = terms_mat.col(b).maxCoeff();
                        Scalar pll = -1e20;
                        if (max_term > -1e19) {
                            Scalar sum_exp = 0.0;
                            for (size_t m = 0; m < M; ++m) {
                                if (terms_mat(m, b) > -1e19) {
                                    sum_exp += std::exp(terms_mat(m, b) - max_term);
                                }
                            }
                            if (sum_exp > 0.0) {
                                pll = max_term + std::log(sum_exp);
                            }
                        }
                        pat_ll_vec[p_start + b] = pll;
                        if (out_pattern_ll) (*out_pattern_ll)[p_start + b] = pll;
                    }
                }
            }
        }
        Scalar total_log_l = 0.0, c_ll = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            LikelihoodEngine::neumaier_acc(total_log_l, c_ll, aln.patterns[p].weight * pat_ll_vec[p]);
        }
        return total_log_l + c_ll;
    }

    // Expectation Step: Computes new mixture weights (and optional saturation sink weight eps_sink) via batched inside-outside GEMM traversal
    std::vector<Scalar> compute_expected_weights(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0,
        Scalar* out_new_eps_sink = nullptr,
        Scalar max_eps_sink = 0.01,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        size_t K = omegas.size();
        if (K <= 1 && !out_new_eps_sink) return {1.0};

        size_t M = syn_rates.size();
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        Scalar bio_scale = (out_new_eps_sink != nullptr) ? (1.0 - eps_sink) : 1.0;
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;

        std::vector<std::vector<std::vector<Matrix>>> P_comp(M, std::vector<std::vector<Matrix>>(K, std::vector<Matrix>(num_nodes)));
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<MG94Matrix> mg_vec(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            mg_vec[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            for (size_t m = 0; m < M; ++m) {
                Scalar alpha_m = syn_rates[m];
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[node.id] : branch_tau[node.id]);
                        P_comp[m][k][node.id] = mg_vec[k].transition_matrix(alpha_m * bl);
                    }
                }
            }
        }
        MG94Matrix mg_bg_step;
        if (num_test_branches + 1 < num_nodes && !tree.bridges.empty()) {
            MG94Parameters p_bg = active_base;
            p_bg.alpha = 1.0;
            p_bg.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega;
            p_bg.delta = delta;
            p_bg.psi = psi;
            mg_bg_step.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P_mix[m][node.id] = Matrix::Zero(S, S);
                    Scalar use_bio = (out_new_eps_sink != nullptr && is_test_branch[node.id]) ? bio_scale : 1.0;
                    for (size_t k = 0; k < K; ++k) {
                        P_mix[m][node.id] += (use_bio * weights[k]) * P_comp[m][k][node.id];
                    }
                    if (out_new_eps_sink != nullptr && is_test_branch[node.id]) {
                        P_mix[m][node.id].noalias() += eps_sink * P_sink;
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_one_set_bridge_mix(s, syn_rates[m], mg_vec, mg_bg_step, weights, (out_new_eps_sink != nullptr) ? eps_sink : 0.0, bio_scale, custom_branch_lengths);
            }
        }

        Matrix pat_weight_accum = Matrix::Zero(K, num_patterns);
        std::vector<Scalar> pat_sink_accum(num_patterns, 0.0);
        std::vector<Scalar> pat_wt_accum(num_patterns, 0.0);
        int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<LikelihoodEngine::BatchedInsideOutsideWorkspace> io_m(M);
            for (size_t m = 0; m < M; ++m) io_m[m].init(num_nodes, S, tile_size);
            Eigen::MatrixXd post_m(M, tile_size);
            Matrix PD_tmp(S, tile_size);

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                post_m.leftCols(B).setZero();
                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    LikelihoodEngine::compute_inside_outside_batch_in_place(
                        tree, aln, p_start, p_end, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, io_m[m],
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                    );
                }

                for (Eigen::Index b = 0; b < B; ++b) {
                    Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] > 1e-12 && io_m[m].likelihood(b) > 0.0) {
                            Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood(b);
                            if (term > max_ll) max_ll = term;
                        }
                    }
                    if (max_ll > -1e19) {
                        Scalar sum_exp = 0.0;
                        for (size_t m = 0; m < M; ++m) {
                            if (syn_weights[m] > 1e-12 && io_m[m].likelihood(b) > 0.0) {
                                Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood(b) - max_ll);
                                post_m(m, b) = e;
                                sum_exp += e;
                            }
                        }
                        if (sum_exp > 0.0) {
                            post_m.col(b) /= sum_exp;
                        }
                    }
                }

                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && is_test_branch[node.id]) {
                        int nid = node.id;
                        auto PD_sub = PD_tmp.leftCols(B);
                        for (size_t k = 0; k < K; ++k) {
                            Scalar prior_k = bio_scale * weights[k];
                            for (size_t m = 0; m < M; ++m) {
                                if (syn_weights[m] <= 1e-12) continue;
                                // Batched GEMM: P_comp[m][k][nid] (S x S) * D[nid] (S x B) -> (S x B)
                                PD_sub.noalias() = P_comp[m][k][nid] * io_m[m].D[nid].leftCols(B);
                                auto V_sub = io_m[m].V[nid].leftCols(B);
                                for (Eigen::Index b = 0; b < B; ++b) {
                                    Scalar pm = post_m(m, b);
                                    Scalar L_b = io_m[m].likelihood(b);
                                    if (pm > 1e-12 && L_b > 0.0) {
                                        size_t p = p_start + b;
                                        Scalar eff_wt = aln.patterns[p].weight;
                                        Scalar ratio = V_sub.col(b).dot(PD_sub.col(b)) / L_b;
                                        pat_weight_accum(k, p) += eff_wt * prior_k * pm * ratio;
                                    }
                                }
                            }
                        }
                        for (Eigen::Index b = 0; b < B; ++b) {
                            if (post_m.col(b).sum() <= 0.0) continue;
                            size_t p = p_start + b;
                            Scalar eff_wt = aln.patterns[p].weight;
                            if (out_new_eps_sink != nullptr) {
                                Scalar sink_ratio = 0.0;
                                for (size_t m = 0; m < M; ++m) {
                                    Scalar pm = post_m(m, b);
                                    Scalar L_b = io_m[m].likelihood(b);
                                    if (pm > 1e-12 && L_b > 0.0) {
                                        sink_ratio += pm * (eval_sink_bilinear(io_m[m].V[nid].col(b), io_m[m].D[nid].col(b)) / L_b);
                                    }
                                }
                                pat_sink_accum[p] += eff_wt * eps_sink * sink_ratio;
                            }
                            pat_wt_accum[p] += eff_wt;
                        }
                    }
                }
            }
        }

        Vector weight_accum = Vector::Zero(K);
        Scalar sink_accum = 0.0;
        Scalar total_pattern_weight = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            weight_accum += pat_weight_accum.col(p);
            sink_accum += pat_sink_accum[p];
            total_pattern_weight += pat_wt_accum[p];
        }

        std::vector<Scalar> new_weights(K);
        if (total_pattern_weight > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t k = 0; k < K; ++k) {
                new_weights[k] = std::clamp(weight_accum(k) / total_pattern_weight, 1e-6, 1.0 - 1e-6);
                sum_w += new_weights[k];
            }
            for (size_t k = 0; k < K; ++k) {
                new_weights[k] /= sum_w;
            }
            if (out_new_eps_sink != nullptr) {
                Scalar raw_eps = sink_accum / std::max(1.0, total_pattern_weight);
                *out_new_eps_sink = std::clamp(raw_eps, 1e-7, max_eps_sink);
            }
        } else {
            new_weights = weights;
            if (out_new_eps_sink != nullptr) {
                *out_new_eps_sink = eps_sink;
            }
        }
        return new_weights;
    }

    // Expectation Step: Computes new synonymous rate weights via batched GEMM site posteriors (with log-underflow scaling)
    std::vector<Scalar> compute_expected_syn_weights(
        Scalar s,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates,
        const std::vector<Scalar>& syn_weights,
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        size_t M = syn_rates.size();
        if (M <= 1) return {1.0};

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        Scalar bio_scale = 1.0 - eps_sink;
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;

        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        bool has_bg = (num_test_branches + 1 < num_nodes);
        MG94Matrix mg_bg;
        if (has_bg) {
            MG94Parameters p_bg = active_base;
            p_bg.alpha = 1.0;
            p_bg.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega;
            p_bg.delta = delta;
            p_bg.psi = psi;
            mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[node.id] : branch_tau[node.id]);
                    if (is_test_branch[node.id]) {
                        P_mix[m][node.id] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink > 0.0) ? bio_scale : 1.0;
                        for (size_t k = 0; k < K; ++k) {
                            P_mix[m][node.id] += (use_bio * weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink > 0.0) {
                            P_mix[m][node.id].noalias() += eps_sink * P_sink;
                        }
                    } else {
                        P_mix[m][node.id] = mg_bg.transition_matrix(alpha_m * bl);
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_one_set_bridge_mix(s, alpha_m, mg, mg_bg, weights, eps_sink, bio_scale, custom_branch_lengths);
            }
        }

        Matrix pat_syn_accum = Matrix::Zero(M, num_patterns);
        std::vector<Scalar> pat_wt_accum(num_patterns, 0.0);
        int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            LikelihoodEngine::BatchedInsideWorkspace ws;
            ws.init(num_nodes, S, tile_size);
            Eigen::MatrixXd terms_mat(M, tile_size);

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                terms_mat.leftCols(B).setConstant(-std::numeric_limits<Scalar>::infinity());
                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    Scalar log_sw = std::log(syn_weights[m]);
                    LikelihoodEngine::compute_inside_batch_in_place(
                        tree, aln, p_start, p_end, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, ws,
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                    );
                    for (Eigen::Index b = 0; b < B; ++b) {
                        if (ws.root_L(b) > 0.0) {
                            terms_mat(m, b) = log_sw + ws.root_log_L(b);
                        }
                    }
                }

                for (Eigen::Index b = 0; b < B; ++b) {
                    Scalar max_term = terms_mat.col(b).maxCoeff();
                    if (max_term > -1e19) {
                        Scalar sum_exp = 0.0;
                        for (size_t m = 0; m < M; ++m) {
                            if (terms_mat(m, b) > -1e19) {
                                sum_exp += std::exp(terms_mat(m, b) - max_term);
                            }
                        }
                        if (sum_exp > 0.0) {
                            size_t p = p_start + b;
                            Scalar eff_wt = aln.patterns[p].weight;
                            for (size_t m = 0; m < M; ++m) {
                                if (terms_mat(m, b) > -1e19) {
                                    pat_syn_accum(m, p) = eff_wt * (std::exp(terms_mat(m, b) - max_term) / sum_exp);
                                }
                            }
                            pat_wt_accum[p] = eff_wt;
                        }
                    }
                }
            }
        }

        Vector syn_accum = Vector::Zero(M);
        Scalar total_pattern_weight = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            syn_accum += pat_syn_accum.col(p);
            total_pattern_weight += pat_wt_accum[p];
        }

        std::vector<Scalar> new_syn_weights(M);
        if (total_pattern_weight > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t m = 0; m < M; ++m) {
                new_syn_weights[m] = std::clamp(syn_accum(m) / total_pattern_weight, 1e-6, 1.0 - 1e-6);
                sum_w += new_syn_weights[m];
            }
            for (size_t m = 0; m < M; ++m) {
                new_syn_weights[m] /= sum_w;
            }
        } else {
            new_syn_weights = syn_weights;
        }
        return new_syn_weights;
    }

    // Evaluates log-likelihood and exact analytical branch gradients d ln L / d t_b
    // under the BUSTED mixture model using batched Inside-Outside GEMM adjoints (with optional SRV and error sink)
    std::pair<Scalar, std::vector<Scalar>> compute_busted_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink = 0.0,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = omegas.size();
        size_t M = syn_rates.size();
        Scalar bio_scale = 1.0 - eps_sink;
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;

        // 1. Build K rate matrices
        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        bool has_bg = (num_test_branches + 1 < num_nodes);
        MG94Matrix mg_bg;
        if (has_bg) {
            MG94Parameters p_bg = active_base;
            p_bg.alpha = 1.0;
            p_bg.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega;
            p_bg.delta = delta;
            p_bg.psi = psi;
            mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        // 2. Precompute P_mix[m] and dP_mix[m] for each synonymous rate category and branch
        // Note: d/dt_b (eps_sink * P_sink) == 0 because the saturation replacement kernel is branch-length invariant!
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> dP_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = branch_lengths[nid];
                    if (is_test_branch[nid]) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        dP_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink > 0.0) ? bio_scale : 1.0;
                        for (size_t k = 0; k < K; ++k) {
                            Matrix P_k = mg[k].transition_matrix(alpha_m * bl);
                            P_mix[m][nid] += (use_bio * weights[k]) * P_k;
                            dP_mix[m][nid] += (use_bio * weights[k]) * alpha_m * (mg[k].Q * P_k);
                        }
                        if (eps_sink > 0.0) {
                            P_mix[m][nid].noalias() += eps_sink * P_sink;
                        }
                    } else {
                        Matrix P_b = mg_bg.transition_matrix(alpha_m * bl);
                        P_mix[m][nid] = P_b;
                        dP_mix[m][nid] = alpha_m * (mg_bg.Q * P_b);
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_one_set_bridge_mix(1.0, alpha_m, mg, mg_bg, weights, eps_sink, bio_scale, &branch_lengths);
            }
        }

        // 3. Parallel batched inside-outside GEMM traversal (deterministic pattern-indexed accumulation)
        std::vector<Scalar> pat_ll_vec(num_patterns, 0.0);
        Matrix pat_grad_mat = Matrix::Zero(num_nodes, num_patterns);
        int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns, 32);
        const size_t tile_size = std::min<size_t>(64, std::max<size_t>(1, num_patterns));
        const size_t num_tiles = (num_patterns + tile_size - 1) / tile_size;

        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<LikelihoodEngine::BatchedInsideOutsideWorkspace> io_m(M);
            for (size_t m = 0; m < M; ++m) io_m[m].init(num_nodes, S, tile_size);
            Eigen::MatrixXd post_w(M, tile_size);
            Matrix dP_D(S, tile_size);

            #pragma omp for schedule(static)
            for (size_t t = 0; t < num_tiles; ++t) {
                size_t p_start = t * tile_size;
                size_t p_end = std::min(num_patterns, p_start + tile_size);
                Eigen::Index B = static_cast<Eigen::Index>(p_end - p_start);

                if (M == 1) {
                    LikelihoodEngine::compute_inside_outside_batch_in_place(
                        tree, aln, p_start, p_end, leaf_to_taxon, P_mix[0], aln.codon_frequencies_f3x4, S, io_m[0],
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[0]
                    );
                    for (Eigen::Index b = 0; b < B; ++b) {
                        size_t p = p_start + b;
                        Scalar eff_wt = aln.patterns[p].weight;
                        if (io_m[0].likelihood(b) > 0.0) {
                            pat_ll_vec[p] = eff_wt * io_m[0].log_likelihood(b);
                        }
                    }
                    auto dP_D_sub = dP_D.leftCols(B);
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            dP_D_sub.noalias() = dP_mix[0][vid] * io_m[0].D[vid].leftCols(B);
                            auto V_sub = io_m[0].V[vid].leftCols(B);
                            for (Eigen::Index b = 0; b < B; ++b) {
                                Scalar L_b = io_m[0].likelihood(b);
                                if (L_b > 0.0) {
                                    size_t p = p_start + b;
                                    Scalar eff_wt = aln.patterns[p].weight;
                                    pat_grad_mat(vid, p) = (eff_wt / L_b) * V_sub.col(b).dot(dP_D_sub.col(b));
                                }
                            }
                        }
                    }
                } else {
                    post_w.leftCols(B).setZero();
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] <= 1e-12) continue;
                        LikelihoodEngine::compute_inside_outside_batch_in_place(
                            tree, aln, p_start, p_end, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, io_m[m],
                            tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                        );
                    }
                    for (Eigen::Index b = 0; b < B; ++b) {
                        Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
                        for (size_t m = 0; m < M; ++m) {
                            if (syn_weights[m] > 1e-12 && io_m[m].likelihood(b) > 0.0) {
                                Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood(b);
                                if (term > max_ll) max_ll = term;
                            }
                        }
                        if (max_ll > -1e19) {
                            Scalar sum_exp = 0.0;
                            for (size_t m = 0; m < M; ++m) {
                                if (syn_weights[m] > 1e-12 && io_m[m].likelihood(b) > 0.0) {
                                    Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood(b) - max_ll);
                                    post_w(m, b) = e;
                                    sum_exp += e;
                                }
                            }
                            if (sum_exp > 0.0) {
                                size_t p = p_start + b;
                                Scalar eff_wt = aln.patterns[p].weight;
                                pat_ll_vec[p] = eff_wt * (max_ll + std::log(sum_exp));
                                post_w.col(b) /= sum_exp;
                            }
                        }
                    }
                    auto dP_D_sub = dP_D.leftCols(B);
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int32_t vid = node.id;
                            for (size_t m = 0; m < M; ++m) {
                                if (syn_weights[m] <= 1e-12) continue;
                                dP_D_sub.noalias() = dP_mix[m][vid] * io_m[m].D[vid].leftCols(B);
                                auto V_sub = io_m[m].V[vid].leftCols(B);
                                for (Eigen::Index b = 0; b < B; ++b) {
                                    Scalar pw = post_w(m, b);
                                    Scalar L_b = io_m[m].likelihood(b);
                                    if (pw > 1e-12 && L_b > 0.0) {
                                        size_t p = p_start + b;
                                        Scalar eff_wt = aln.patterns[p].weight;
                                        pat_grad_mat(vid, p) += eff_wt * pw * (V_sub.col(b).dot(dP_D_sub.col(b)) / L_b);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        Scalar total_log_l = LikelihoodEngine::neumaier_sum(pat_ll_vec);
        std::vector<Scalar> grad_b(num_nodes, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            Scalar s_acc = 0.0, c_acc = 0.0;
            for (size_t p = 0; p < num_patterns; ++p) {
                LikelihoodEngine::neumaier_acc(s_acc, c_acc, pat_grad_mat(i, p));
            }
            grad_b[i] = s_acc + c_acc;
        }

        return {total_log_l, grad_b};
    }

    // High-performance joint branch length refinement via exact Inside-Outside adjoint gradients and L-BFGS in log-branch-length space
    std::vector<Scalar> refine_branch_lengths(
        const std::vector<Scalar>& init_bl,
        const std::vector<Scalar>& omegas,
        const std::vector<Scalar>& weights,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        int max_iters = 25,
        Scalar eps_sink = 0.0,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        std::vector<Scalar> log_bl(num_nodes, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-8, 10000.0);
                log_bl[i] = std::log(bl[i]);
            }
        }

        // L-BFGS state in log-branch-length space y_i = log(t_i), where d ln L / d y_i = t_i * (d ln L / d t_i)
        const size_t m_history = 8;
        std::vector<std::vector<Scalar>> s_hist;
        std::vector<std::vector<Scalar>> y_hist;
        std::vector<Scalar> rho_hist;

        auto eval_log_space = [&](const std::vector<Scalar>& y_vec) -> std::tuple<Scalar, std::vector<Scalar>, std::vector<Scalar>> {
            std::vector<Scalar> cur_bl(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    cur_bl[i] = std::exp(std::clamp(y_vec[i], -18.420680744, 9.210340372));
                }
            }
            auto [ll, grad_t] = compute_busted_branch_gradients(cur_bl, omegas, weights, syn_rates, syn_weights, delta, psi, eps_sink, custom_base_params);
            std::vector<Scalar> grad_y(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    grad_y[i] = cur_bl[i] * grad_t[i];
                }
            }
            return {ll, grad_y, cur_bl};
        };

        auto [curr_ll, curr_grad, curr_bl] = eval_log_space(log_bl);
        bl = curr_bl;

        Scalar best_ll = curr_ll;
        std::vector<Scalar> best_bl = bl;

        for (int iter = 0; iter < max_iters; ++iter) {
            // Check convergence in log-space gradient
            Scalar max_g = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_g = std::max(max_g, std::abs(curr_grad[i]));
                }
            }
            if (max_g < 1e-4) break;

            // Two-loop L-BFGS recursion to compute search direction d = H * grad
            std::vector<Scalar> q = curr_grad;
            std::vector<Scalar> alpha(s_hist.size());

            for (int j = static_cast<int>(s_hist.size()) - 1; j >= 0; --j) {
                Scalar s_dot_q = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_q += s_hist[j][i] * q[i];
                }
                alpha[j] = rho_hist[j] * s_dot_q;
                for (size_t i = 0; i < num_nodes; ++i) {
                    q[i] -= alpha[j] * y_hist[j][i];
                }
            }

            // Initial Hessian scaling H0 = gamma * I
            Scalar gamma = 0.25;
            if (!s_hist.empty()) {
                Scalar s_dot_y = 0.0;
                Scalar y_dot_y = 0.0;
                const auto& s_last = s_hist.back();
                const auto& y_last = y_hist.back();
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_y += s_last[i] * y_last[i];
                    y_dot_y += y_last[i] * y_last[i];
                }
                if (y_dot_y > 1e-12) {
                    gamma = std::clamp(s_dot_y / y_dot_y, 1e-3, 5.0);
                }
            }

            std::vector<Scalar> r(num_nodes);
            for (size_t i = 0; i < num_nodes; ++i) {
                r[i] = gamma * q[i];
            }

            for (size_t j = 0; j < s_hist.size(); ++j) {
                Scalar y_dot_r = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    y_dot_r += y_hist[j][i] * r[i];
                }
                Scalar beta = rho_hist[j] * y_dot_r;
                for (size_t i = 0; i < num_nodes; ++i) {
                    r[i] += s_hist[j][i] * (alpha[j] - beta);
                }
            }

            // Direction d (for maximization, d points along r)
            std::vector<Scalar> d = r;
            d[tree.root_id] = 0.0;

            // Ensure ascent direction (d . grad > 0)
            Scalar dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    dir_dot_grad += d[i] * curr_grad[i];
                }
            }
            if (dir_dot_grad <= 0.0) {
                d = curr_grad;
                d[tree.root_id] = 0.0;
                dir_dot_grad = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        dir_dot_grad += d[i] * curr_grad[i];
                    }
                }
                s_hist.clear();
                y_hist.clear();
                rho_hist.clear();
            }

            // Clamp per-branch log-step to [-2.0, +2.0] (at most ~7.4x multiplier per step) without throttling other branches!
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    d[i] = std::clamp(d[i], -2.0, 2.0);
                }
            }
            dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    dir_dot_grad += d[i] * curr_grad[i];
                }
            }

            Scalar step_size = 1.0;
            std::vector<Scalar> new_log_bl(num_nodes, 0.0);
            std::vector<Scalar> new_bl(num_nodes, 0.0);
            Scalar new_ll = curr_ll;
            std::vector<Scalar> new_grad;

            bool line_search_ok = false;
            for (int ls = 0; ls < 12; ++ls) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) {
                        new_log_bl[i] = std::clamp(log_bl[i] + step_size * d[i], -18.420680744, 9.210340372);
                    } else {
                        new_log_bl[i] = 0.0;
                    }
                }
                auto [test_ll, test_grad, test_bl] = eval_log_space(new_log_bl);
                if (test_ll > curr_ll + 1e-4 * step_size * dir_dot_grad || (test_ll > curr_ll + 1e-5 && ls >= 3)) {
                    new_ll = test_ll;
                    new_grad = test_grad;
                    new_bl = test_bl;
                    line_search_ok = true;
                    break;
                }
                step_size *= 0.5;
            }

            if (!line_search_ok) {
                break;
            }

            // Update L-BFGS history in log-branch-length space
            std::vector<Scalar> s_k(num_nodes);
            std::vector<Scalar> y_k(num_nodes);
            Scalar s_dot_y = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                s_k[i] = new_log_bl[i] - log_bl[i];
                y_k[i] = curr_grad[i] - new_grad[i]; // for maximization
                s_dot_y += s_k[i] * y_k[i];
            }

            if (s_dot_y > 1e-8) {
                if (s_hist.size() >= m_history) {
                    s_hist.erase(s_hist.begin());
                    y_hist.erase(y_hist.begin());
                    rho_hist.erase(rho_hist.begin());
                }
                s_hist.push_back(s_k);
                y_hist.push_back(y_k);
                rho_hist.push_back(1.0 / s_dot_y);
            }

            log_bl = new_log_bl;
            bl = new_bl;
            curr_ll = new_ll;
            curr_grad = new_grad;

            if (curr_ll > best_ll) {
                best_ll = curr_ll;
                best_bl = bl;
            }
        }

        return best_bl;
    }

    // Fits BUSTED model for specified number of rate categories K >= 1, optional SRV (M >= 1), and optional saturation error sink
    BUSTEDFit fit_model(
        bool is_constrained,
        size_t K = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        const std::string& model_desc = "",
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01,
        const BUSTEDFit* warm_start_fit = nullptr
    ) const {
        K = std::max<size_t>(1, K);
        BUSTEDFit fit;
        fit.num_rate_classes = K;
        Scalar s = 1.0;

        bool allow_double = (multiple_hits == "Double" || multiple_hits == "Double+Triple");
        bool allow_triple = (multiple_hits == "Double+Triple");
        Scalar delta = warm_start_fit ? warm_start_fit->delta : 0.0;
        Scalar psi = warm_start_fit ? warm_start_fit->psi : 0.0;
        Scalar eps_sink = use_error_sink
            ? (warm_start_fit && warm_start_fit->error_sink_weight > 0.0
                ? std::min(warm_start_fit->error_sink_weight, max_eps_sink)
                : std::min(0.005, max_eps_sink * 0.5))
            : 0.0;

        MG94Parameters local_base = (warm_start_fit && warm_start_fit->has_gtr_params)
            ? warm_start_fit->gtr_params
            : base_params;

        // Initialize SRV synonymous rate distribution for arbitrary M >= 1
        size_t M = srv ? std::max<size_t>(1, num_syn_rates) : 1;
        std::vector<Scalar> syn_rates(M, 1.0);
        std::vector<Scalar> syn_weights(M, 1.0 / static_cast<Scalar>(M));

        if (srv && M == 2) {
            syn_rates = {0.35, 1.65};
            syn_weights = {0.50, 0.50};
        } else if (srv && M == 3) {
            syn_rates = {0.30, 1.00, 3.10};
            syn_weights = {0.30, 0.60, 0.10};
        } else if (srv && M >= 4) {
            Scalar sum_w = 0.0;
            for (size_t m = 0; m < M; ++m) {
                Scalar frac = static_cast<Scalar>(m) / static_cast<Scalar>(M - 1);
                syn_rates[m] = 0.15 * std::pow(25.0, frac);
                syn_weights[m] = 1.0 / static_cast<Scalar>(M);
                sum_w += syn_weights[m];
            }
            for (size_t m = 0; m < M; ++m) syn_weights[m] /= sum_w;
            Scalar mean_r = 0.0;
            for (size_t m = 0; m < M; ++m) mean_r += syn_weights[m] * syn_rates[m];
            if (mean_r > 1e-12) {
                for (size_t m = 0; m < M; ++m) syn_rates[m] /= mean_r;
            }
        }
        if (warm_start_fit && srv && warm_start_fit->test_distribution.syn_rates.size() == M &&
            warm_start_fit->test_distribution.syn_weights.size() == M) {
            syn_rates = warm_start_fit->test_distribution.syn_rates;
            syn_weights = warm_start_fit->test_distribution.syn_weights;
        }

        const std::vector<Scalar>* ecm_bl_ptr = nullptr;
        if (warm_start_fit && !warm_start_fit->branch_lengths.empty()) {
            ecm_bl_ptr = &warm_start_fit->branch_lengths;
            s = 1.0;
        }

        // Helper to sort (rates, weights) in ascending rate order without hard penalty barriers
        auto sort_rate_weight_pairs = [](std::vector<Scalar>& rates, std::vector<Scalar>& wts, size_t n_sort) {
            if (n_sort <= 1) return;
            std::vector<size_t> idx(n_sort);
            for (size_t i = 0; i < n_sort; ++i) idx[i] = i;
            std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return rates[a] < rates[b]; });
            std::vector<Scalar> r_tmp(n_sort), w_tmp(n_sort);
            for (size_t i = 0; i < n_sort; ++i) {
                r_tmp[i] = rates[idx[i]];
                w_tmp[i] = wts[idx[i]];
            }
            for (size_t i = 0; i < n_sort; ++i) {
                rates[i] = r_tmp[i];
                wts[i] = w_tmp[i];
            }
        };

        // Helper to renormalize SRV rates to exact unit mean and absorb the mean into s or bl
        auto renormalize_srv_unit_mean = [&](Scalar& scale_s, std::vector<Scalar>* bl_vec) {
            if (!srv || M < 2) return;
            sort_rate_weight_pairs(syn_rates, syn_weights, M);
            Scalar mean_r = 0.0;
            for (size_t m = 0; m < M; ++m) {
                mean_r += syn_weights[m] * syn_rates[m];
            }
            if (mean_r > 1e-8 && std::abs(mean_r - 1.0) > 1e-12) {
                for (size_t m = 0; m < M; ++m) {
                    syn_rates[m] /= mean_r;
                }
                if (bl_vec) {
                    for (size_t i = 0; i < bl_vec->size(); ++i) {
                        if (i != static_cast<size_t>(tree.root_id)) {
                            (*bl_vec)[i] = std::clamp((*bl_vec)[i] * mean_r, 1e-8, 10000.0);
                        }
                    }
                } else {
                    scale_s = std::clamp(scale_s * mean_r, 1e-4, 100.0);
                }
            }
        };

        // Helper: optimize SRV synonymous rates for any M >= 2 without hard ordering penalty (F-019, F-020, F-021)
        auto optimize_srv_rates = [&](Scalar& cur_s, const std::vector<Scalar>& cur_o, const std::vector<Scalar>& cur_w, std::vector<Scalar>* bl_vec, int max_it) {
            if (!srv || M < 2) return;
            const std::vector<Scalar>* b_ptr = bl_vec ? bl_vec : ecm_bl_ptr;
            Scalar eval_s = bl_vec ? 1.0 : cur_s;
            if (M == 2) {
                auto [r_opt, nll_opt] = Brent1D::minimize([&](Scalar ratio) {
                    Scalar q0 = syn_weights[0], q1 = syn_weights[1];
                    Scalar denom = std::max(1e-12, q0 * ratio + q1);
                    Scalar a0 = ratio / denom;
                    Scalar a1 = 1.0 / denom;
                    return -evaluate_log_l(eval_s, cur_o, cur_w, {a0, a1}, syn_weights, b_ptr, nullptr, delta, psi, eps_sink, &local_base);
                }, 1e-3, std::clamp(syn_rates[0] / std::max(1e-6, syn_rates[1]), 1e-3, 100.0), 100.0, 1e-4, max_it);
                Scalar denom = std::max(1e-12, syn_weights[0] * r_opt + syn_weights[1]);
                std::vector<Scalar> cand_r = {r_opt / denom, 1.0 / denom};
                if (evaluate_log_l(eval_s, cur_o, cur_w, cand_r, syn_weights, b_ptr, nullptr, delta, psi, eps_sink, &local_base) >=
                    evaluate_log_l(eval_s, cur_o, cur_w, syn_rates, syn_weights, b_ptr, nullptr, delta, psi, eps_sink, &local_base)) {
                    syn_rates = cand_r;
                }
            } else {
                size_t D_syn = M - 1;
                std::vector<Scalar> x0(D_syn), lb(D_syn, 1e-3), ub(D_syn, 100.0);
                Scalar ref_r = std::max(1e-6, syn_rates[M - 1]);
                for (size_t m = 0; m < D_syn; ++m) {
                    x0[m] = std::clamp(syn_rates[m] / ref_r, 1e-3, 100.0);
                }
                auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                    std::vector<Scalar> cand_r(M, 1.0);
                    Scalar mean_c = syn_weights[M - 1];
                    for (size_t m = 0; m < D_syn; ++m) {
                        cand_r[m] = std::max(1e-4, x[m]);
                        mean_c += syn_weights[m] * cand_r[m];
                    }
                    mean_c = std::max(1e-12, mean_c);
                    for (size_t m = 0; m < M; ++m) cand_r[m] /= mean_c;
                    return -evaluate_log_l(eval_s, cur_o, cur_w, cand_r, syn_weights, b_ptr, nullptr, delta, psi, eps_sink, &local_base);
                }, x0, lb, ub, 1e-4, max_it);
                std::vector<Scalar> cand_r(M, 1.0);
                Scalar mean_c = syn_weights[M - 1];
                for (size_t m = 0; m < D_syn; ++m) {
                    cand_r[m] = std::max(1e-4, opt.x[m]);
                    mean_c += syn_weights[m] * cand_r[m];
                }
                mean_c = std::max(1e-12, mean_c);
                for (size_t m = 0; m < M; ++m) cand_r[m] /= mean_c;
                if (evaluate_log_l(eval_s, cur_o, cur_w, cand_r, syn_weights, b_ptr, nullptr, delta, psi, eps_sink, &local_base) >=
                    evaluate_log_l(eval_s, cur_o, cur_w, syn_rates, syn_weights, b_ptr, nullptr, delta, psi, eps_sink, &local_base)) {
                    syn_rates = cand_r;
                }
            }
            renormalize_srv_unit_mean(cur_s, bl_vec);
        };

        // Helper: if a purifying rate class (k < K-1) has collapsed in weight (< 0.005) or merged with a neighbor,
        // propose a split of the purifying mass into (0.0, omega_mid, omega_3) and accept ONLY if it strictly improves lnL
        auto try_rescue_collapsed_classes = [&](
            Scalar& cur_s,
            std::vector<Scalar>& cur_o,
            std::vector<Scalar>& cur_w,
            const std::vector<Scalar>* bl_ptr
        ) {
            if (K != 3) return;
            bool collapsed = (cur_w[0] < 0.005) || (cur_w[1] < 0.005) || (std::abs(cur_o[1] - cur_o[0]) < 1e-3);
            if (!collapsed) return;

            Scalar base_ll = evaluate_log_l(cur_s, cur_o, cur_w, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
            size_t dom = (cur_w[0] >= cur_w[1]) ? 0 : 1;
            Scalar w_pur = std::clamp(cur_w[0] + cur_w[1], 0.04, 0.99);
            Scalar w_top = std::clamp(1.0 - w_pur, 0.01, 0.96);
            Scalar o_dom = std::clamp(cur_o[dom], 0.002, 0.80);
            Scalar o_mg = std::clamp(mg94_omega, 0.01, 0.80);
            Scalar o_top = is_constrained ? 1.0 : std::max(1.05, cur_o[2]);

            std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> rescue_cands = {
                {{0.0, std::min(0.95, std::max(0.02, o_dom * 2.2)), o_top}, {0.65 * w_pur, 0.35 * w_pur, w_top}},
                {{0.0, std::min(0.95, std::max(0.02, o_mg * 0.5)), o_top},  {0.65 * w_pur, 0.35 * w_pur, w_top}},
                {{0.0, std::min(0.95, o_mg), o_top},                        {0.70 * w_pur, 0.30 * w_pur, w_top}}
            };

            std::vector<Scalar> best_rc_o = cur_o;
            std::vector<Scalar> best_rc_w = cur_w;
            Scalar best_rc_ll = -1e20;
            for (const auto& [rc_o, rc_w] : rescue_cands) {
                Scalar ll = evaluate_log_l(cur_s, rc_o, rc_w, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                if (ll > best_rc_ll) {
                    best_rc_ll = ll;
                    best_rc_o = rc_o;
                    best_rc_w = rc_w;
                }
            }

            // Polish the best rescue candidate via 2 alternating (M-step on omegas, E-step on weights) rounds
            Scalar cand_s = cur_s;
            auto cand_o = best_rc_o;
            auto cand_w = best_rc_w;
            Scalar max_omega = use_error_sink ? 100.0 : 10000.0;
            for (int r_it = 0; r_it < 2; ++r_it) {
                if (is_constrained) {
                    auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                        Scalar w0 = std::min(x[1], x[2]);
                        Scalar w1 = std::max(x[1], x[2]);
                        return -evaluate_log_l(x[0], {w0, w1, 1.0}, cand_w, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                    }, {cand_s, cand_o[0], cand_o[1]}, {0.01, 0.0, 0.0}, {50.0, 1.0, 1.0}, 1e-4, 40);
                    cand_s = opt.x[0];
                    cand_o = {std::min(opt.x[1], opt.x[2]), std::max(opt.x[1], opt.x[2]), 1.0};
                } else {
                    auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                        Scalar w0 = std::min(x[1], x[2]);
                        Scalar w1 = std::max(x[1], x[2]);
                        return -evaluate_log_l(x[0], {w0, w1, x[3]}, cand_w, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                    }, {cand_s, cand_o[0], cand_o[1], cand_o[2]}, {0.01, 0.0, 0.0, 1.0}, {50.0, 1.0, 1.0, max_omega}, 1e-4, 45);
                    cand_s = opt.x[0];
                    cand_o = {std::min(opt.x[1], opt.x[2]), std::max(opt.x[1], opt.x[2]), opt.x[3]};
                }
                cand_w = compute_expected_weights(cand_s, cand_o, cand_w, syn_rates, syn_weights, bl_ptr, delta, psi, eps_sink, nullptr, max_eps_sink, &local_base);
            }
            Scalar final_rc_ll = evaluate_log_l(cand_s, cand_o, cand_w, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
            if (final_rc_ll > base_ll + 1e-4) {
                cur_s = cand_s;
                cur_o = cand_o;
                cur_w = cand_w;
            }
        };

        // Initialize rates and weights for arbitrary K >= 1
        std::vector<Scalar> omegas(K, 1.0);
        std::vector<Scalar> weights(K, 1.0 / static_cast<Scalar>(K));

        if (K == 1) {
            weights[0] = 1.0;
            omegas[0] = is_constrained ? std::min(1.0, mg94_omega) : mg94_omega;
            if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == 1) {
                omegas[0] = is_constrained
                    ? std::min(1.0, warm_start_fit->test_distribution.omegas[0])
                    : warm_start_fit->test_distribution.omegas[0];
            }
        } else if (K == 2) {
            std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> candidates;
            if (is_constrained) {
                candidates = {
                    {{0.10, 1.0}, {0.75, 0.25}},
                    {{0.02, 1.0}, {0.85, 0.15}},
                    {{0.25, 1.0}, {0.60, 0.40}},
                    {{0.0, 1.0},  {0.70, 0.30}}
                };
                if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == 2) {
                    Scalar w0 = std::min(1.0, warm_start_fit->test_distribution.omegas[0]);
                    Scalar p0 = std::clamp(warm_start_fit->test_distribution.weights[0], 0.02, 0.98);
                    candidates.push_back({{w0, 1.0}, {p0, 1.0 - p0}});
                }
            } else {
                candidates = {
                    {{0.10, 3.0}, {0.75, 0.25}},
                    {{0.02, 4.5}, {0.90, 0.10}},
                    {{0.05, 2.0}, {0.80, 0.20}},
                    {{0.02, 3.0}, {0.95, 0.05}}
                };
                if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == 2) {
                    Scalar w0 = std::min(1.0, warm_start_fit->test_distribution.omegas[0]);
                    Scalar w1 = std::max(1.0, warm_start_fit->test_distribution.omegas[1]);
                    Scalar p0 = std::clamp(warm_start_fit->test_distribution.weights[0], 0.02, 0.98);
                    candidates.push_back({{w0, w1}, {p0, 1.0 - p0}});
                    candidates.push_back({{w0, 1.35}, {p0, 1.0 - p0}});
                    candidates.push_back({{w0, 2.5}, {p0, 1.0 - p0}});
                    candidates.push_back({{w0, 3.5}, {0.95, 0.05}});
                    candidates.push_back({{w0, 5.0}, {std::max(0.5, p0), 1.0 - std::max(0.5, p0)}});
                }
            }
            Scalar best_start_ll = -1e20;
            for (const auto& [cand_o, cand_w] : candidates) {
                Scalar ll = evaluate_log_l(s, cand_o, cand_w, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                if (ll > best_start_ll) {
                    best_start_ll = ll;
                    omegas = cand_o;
                    weights = cand_w;
                }
            }
        } else if (K == 3) {
            Scalar w_mg = std::clamp(mg94_omega, 0.01, 0.90);
            if (is_constrained) {
                std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> candidates = {
                    {{0.05, 0.50, 1.0}, {0.45, 0.45, 0.10}},
                    {{0.01, 0.20, 1.0}, {0.30, 0.30, 0.40}},
                    {{0.10, 0.80, 1.0}, {0.20, 0.20, 0.60}},
                    {{0.001, 0.05, 1.0}, {0.15, 0.15, 0.70}},
                    {{0.0, 0.5 * w_mg, 1.0}, {0.65, 0.27, 0.08}},
                    {{0.0, w_mg, 1.0},       {0.65, 0.28, 0.07}},
                    {{0.0, 1.0, 1.0},        {0.40, 0.30, 0.30}},
                    {{0.05, 1.0, 1.0},       {0.35, 0.35, 0.30}}
                };
                if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == 3) {
                    const auto& co = warm_start_fit->test_distribution.omegas;
                    const auto& cw = warm_start_fit->test_distribution.weights;
                    Scalar o0 = std::clamp(co[0], 0.0, 1.0);
                    Scalar o1 = std::clamp(co[1], o0, 1.0);
                    Scalar p0 = std::clamp(cw[0], 0.001, 0.98);
                    Scalar p1 = std::clamp(cw[1], 0.001, 0.99 - p0);
                    Scalar p2 = std::max(0.001, 1.0 - p0 - p1);
                    Scalar sum_p = p0 + p1 + p2;
                    candidates.push_back({{o0, o1, 1.0}, {p0 / sum_p, p1 / sum_p, p2 / sum_p}});
                    Scalar p_pur = std::clamp(cw[0] + cw[1], 0.10, 0.97);
                    candidates.push_back({{0.0, std::min(0.95, std::max(0.02, o1)), 1.0}, {0.68 * p_pur, 0.32 * p_pur, 1.0 - p_pur}});
                    candidates.push_back({{o0, 1.0, 1.0}, {p_pur * 0.85, p_pur * 0.15, 1.0 - p_pur}});
                }
                Scalar best_start_ll = -1e20;
                bool best_use_warm_bl = (ecm_bl_ptr != nullptr);
                Scalar best_s = s;
                std::vector<Scalar> best_syn_r = syn_rates;
                std::vector<Scalar> best_syn_w = syn_weights;
                std::vector<Scalar> fresh_syn_r = (srv && M == 3) ? std::vector<Scalar>{0.30, 1.00, 3.10} : syn_rates;
                std::vector<Scalar> fresh_syn_w = (srv && M == 3) ? std::vector<Scalar>{0.30, 0.60, 0.10} : syn_weights;
                for (const auto& [cand_o, cand_w] : candidates) {
                    Scalar ll = evaluate_log_l(s, cand_o, cand_w, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                    if (ll > best_start_ll) {
                        best_start_ll = ll;
                        omegas = cand_o;
                        weights = cand_w;
                        best_use_warm_bl = (ecm_bl_ptr != nullptr);
                        best_s = s;
                        best_syn_r = syn_rates;
                        best_syn_w = syn_weights;
                    }
                    if (ecm_bl_ptr != nullptr) {
                        Scalar ll_fresh = evaluate_log_l(1.0, cand_o, cand_w, fresh_syn_r, fresh_syn_w, nullptr, nullptr, delta, psi, eps_sink, &local_base);
                        if (ll_fresh > best_start_ll) {
                            best_start_ll = ll_fresh;
                            omegas = cand_o;
                            weights = cand_w;
                            best_use_warm_bl = false;
                            best_s = 1.0;
                            best_syn_r = fresh_syn_r;
                            best_syn_w = fresh_syn_w;
                        }
                    }
                }
                if (!best_use_warm_bl) {
                    ecm_bl_ptr = nullptr;
                    s = best_s;
                    syn_rates = best_syn_r;
                    syn_weights = best_syn_w;
                }
            } else {
                std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> candidates = {
                    {{0.05, 0.50, 4.0}, {0.50, 0.40, 0.10}},
                    {{0.01, 0.50, 2.0}, {0.40, 0.10, 0.50}},
                    {{0.01, 1.00, 2.0}, {0.40, 0.05, 0.55}},
                    {{0.02, 0.20, 3.0}, {0.35, 0.35, 0.30}},
                    {{0.10, 0.90, 5.0}, {0.30, 0.40, 0.30}},
                    {{0.005, 0.10, 1.8}, {0.45, 0.15, 0.40}},
                    {{0.0, 0.5 * w_mg, 2.0}, {0.65, 0.28, 0.07}},
                    {{0.0, w_mg, 2.5},       {0.68, 0.27, 0.05}}
                };
                if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == 3) {
                    const auto& co = warm_start_fit->test_distribution.omegas;
                    const auto& cw = warm_start_fit->test_distribution.weights;
                    Scalar o0 = std::clamp(co[0], 0.0, 1.0);
                    Scalar o1 = std::clamp(co[1], o0, 1.0);
                    Scalar o2 = std::max(1.0, co[2]);
                    candidates.push_back({{o0, o1, o2}, cw});
                    Scalar p2 = std::max(0.02, cw[2] * 0.7);
                    Scalar p0 = std::clamp(cw[0], 0.02, 0.92);
                    Scalar p1 = std::max(0.02, 1.0 - p0 - p2);
                    Scalar sum_p = p0 + p1 + p2;
                    candidates.push_back({{o0, std::min(0.95, std::max(o0 + 0.01, o1)), 1.6}, {p0 / sum_p, p1 / sum_p, p2 / sum_p}});
                    candidates.push_back({{o0, std::min(0.95, std::max(o0 + 0.01, o1)), 2.5}, {p0 / sum_p, p1 / sum_p, p2 / sum_p}});
                    candidates.push_back({{o0, 1.0, 2.5}, {p0 / sum_p, p1 / sum_p, p2 / sum_p}});
                }
                Scalar best_start_ll = -1e20;
                Scalar best_start_eps = eps_sink;
                std::vector<Scalar> eps_cands = use_error_sink
                    ? std::vector<Scalar>{eps_sink, std::min(0.001, max_eps_sink), std::min(0.005, max_eps_sink), std::min(0.01, max_eps_sink)}
                    : std::vector<Scalar>{0.0};
                if (use_error_sink && warm_start_fit && warm_start_fit->error_sink_weight > 0.0) {
                    eps_cands.push_back(std::min(warm_start_fit->error_sink_weight, max_eps_sink));
                }
                for (Scalar cand_e : eps_cands) {
                    for (const auto& [cand_o, cand_w] : candidates) {
                        Scalar ll = evaluate_log_l(s, cand_o, cand_w, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, cand_e, &local_base);
                        if (ll > best_start_ll) {
                            best_start_ll = ll;
                            omegas = cand_o;
                            weights = cand_w;
                            best_start_eps = cand_e;
                        }
                    }
                }
                eps_sink = best_start_eps;
            }
        } else {
            // K >= 4: construct K-element candidate distributions (F-004)
            std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> candidates;
            for (Scalar top_o : (is_constrained ? std::vector<Scalar>{1.0} : std::vector<Scalar>{2.0, 4.0})) {
                for (Scalar p_top : {0.10, 0.25}) {
                    std::vector<Scalar> co(K), cw(K);
                    Scalar rem_w = 1.0 - p_top;
                    for (size_t k = 0; k + 1 < K; ++k) {
                        Scalar frac = (K > 2) ? (static_cast<Scalar>(k) / static_cast<Scalar>(K - 2)) : 0.5;
                        co[k] = 0.01 * std::pow(80.0, frac);
                        cw[k] = rem_w / static_cast<Scalar>(K - 1);
                    }
                    co[K - 1] = top_o;
                    cw[K - 1] = p_top;
                    candidates.push_back({co, cw});
                }
            }
            if (warm_start_fit && warm_start_fit->test_distribution.omegas.size() == K) {
                auto co = warm_start_fit->test_distribution.omegas;
                auto cw = warm_start_fit->test_distribution.weights;
                for (size_t k = 0; k + 1 < K; ++k) co[k] = std::clamp(co[k], 0.0, 1.0);
                co[K - 1] = is_constrained ? 1.0 : std::max(1.0, co[K - 1]);
                candidates.push_back({co, cw});
            }
            Scalar best_start_ll = -1e20;
            for (const auto& [cand_o, cand_w] : candidates) {
                Scalar ll = evaluate_log_l(s, cand_o, cand_w, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                if (ll > best_start_ll) {
                    best_start_ll = ll;
                    omegas = cand_o;
                    weights = cand_w;
                }
            }
        }

        if (allow_double) {
            Scalar best_grid_ll = -1e20;
            std::vector<Scalar> d_cands = {0.0, 0.02, 0.08, 0.20, 0.60};
            if (warm_start_fit && warm_start_fit->delta > 0.0) d_cands.push_back(warm_start_fit->delta);
            std::vector<Scalar> p_cands = allow_triple ? std::vector<Scalar>{0.0, 0.02, 0.08} : std::vector<Scalar>{0.0};
            if (allow_triple && warm_start_fit && warm_start_fit->psi > 0.0) p_cands.push_back(warm_start_fit->psi);
            for (Scalar d_c : d_cands) {
                for (Scalar p_c : p_cands) {
                    Scalar cand_ll = evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, d_c, p_c, eps_sink, &local_base);
                    if (cand_ll > best_grid_ll) {
                        best_grid_ll = cand_ll;
                        delta = d_c;
                        psi = p_c;
                    }
                }
            }
        }

        // General SQUAREM helper for arbitrary K >= 2simplex weights
        auto squarem_update_weights = [&]() {
            if (K == 1) {
                if (use_error_sink) {
                    weights = compute_expected_weights(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps_sink, &eps_sink, max_eps_sink, &local_base);
                }
                return;
            }
            for (int sq = 0; sq < 2; ++sq) {
                Scalar eps0 = eps_sink;
                auto w0 = weights;
                Scalar eps1 = eps0;
                auto w1 = compute_expected_weights(s, omegas, w0, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps0, use_error_sink ? &eps1 : nullptr, max_eps_sink, &local_base);
                Scalar eps2 = eps1;
                auto w2 = compute_expected_weights(s, omegas, w1, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps1, use_error_sink ? &eps2 : nullptr, max_eps_sink, &local_base);
                eps_sink = eps2;
                Scalar r2 = 0.0, v2 = 0.0;
                std::vector<Scalar> r_vec(K - 1), v_vec(K - 1);
                for (size_t k = 0; k + 1 < K; ++k) {
                    r_vec[k] = w1[k] - w0[k];
                    v_vec[k] = (w2[k] - w1[k]) - r_vec[k];
                    r2 += r_vec[k] * r_vec[k];
                    v2 += v_vec[k] * v_vec[k];
                }
                if (v2 > 1e-12) {
                    Scalar step = -std::sqrt(r2 / v2);
                    std::vector<Scalar> w_acc(K, 0.0);
                    Scalar sum_p = 0.0;
                    for (size_t k = 0; k + 1 < K; ++k) {
                        w_acc[k] = std::clamp(w0[k] - 2.0 * step * r_vec[k] + step * step * v_vec[k], 1e-5, 1.0 - 1e-5);
                        sum_p += w_acc[k];
                    }
                    if (sum_p > 1.0 - 1e-5) {
                        for (size_t k = 0; k + 1 < K; ++k) w_acc[k] = (w_acc[k] / sum_p) * 0.99;
                        sum_p = 0.99;
                    }
                    w_acc[K - 1] = 1.0 - sum_p;
                    Scalar eps_stab = eps2;
                    auto w_stab = compute_expected_weights(s, omegas, w_acc, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps2, use_error_sink ? &eps_stab : nullptr, max_eps_sink, &local_base);
                    if (evaluate_log_l(s, omegas, w_stab, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_stab, &local_base) >=
                        evaluate_log_l(s, omegas, w2, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps2, &local_base)) {
                        weights = w_stab;
                        eps_sink = eps_stab;
                    } else {
                        weights = w2;
                        eps_sink = eps2;
                    }
                } else {
                    weights = w2;
                    eps_sink = eps2;
                }
            }
            if (use_error_sink) {
                Scalar best_eps = eps_sink;
                Scalar best_eps_ll = evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                for (Scalar cand_e : {1e-7, 2e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2}) {
                    if (cand_e > max_eps_sink) continue;
                    Scalar ll_e = evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, cand_e, &local_base);
                    if (ll_e > best_eps_ll + 1e-4) {
                        best_eps_ll = ll_e;
                        best_eps = cand_e;
                    }
                }
                eps_sink = best_eps;
            }
        };

        auto squarem_update_syn_weights = [&]() {
            if (!srv || M < 2) return;
            for (int sq = 0; sq < 2; ++sq) {
                auto q0 = syn_weights;
                auto q1 = compute_expected_syn_weights(s, omegas, weights, syn_rates, q0, ecm_bl_ptr, delta, psi, eps_sink, &local_base);
                auto q2 = compute_expected_syn_weights(s, omegas, weights, syn_rates, q1, ecm_bl_ptr, delta, psi, eps_sink, &local_base);
                Scalar r2 = 0.0, v2 = 0.0;
                std::vector<Scalar> r_vec(M - 1), v_vec(M - 1);
                for (size_t m = 0; m + 1 < M; ++m) {
                    r_vec[m] = q1[m] - q0[m];
                    v_vec[m] = (q2[m] - q1[m]) - r_vec[m];
                    r2 += r_vec[m] * r_vec[m];
                    v2 += v_vec[m] * v_vec[m];
                }
                if (v2 > 1e-12) {
                    Scalar step = -std::sqrt(r2 / v2);
                    std::vector<Scalar> q_cand(M, 0.0);
                    Scalar sum_q = 0.0;
                    for (size_t m = 0; m + 1 < M; ++m) {
                        q_cand[m] = std::clamp(q0[m] - 2.0 * step * r_vec[m] + step * step * v_vec[m], 1e-5, 1.0 - 1e-5);
                        sum_q += q_cand[m];
                    }
                    if (sum_q > 1.0 - 1e-5) {
                        for (size_t m = 0; m + 1 < M; ++m) q_cand[m] = (q_cand[m] / sum_q) * 0.99;
                        sum_q = 0.99;
                    }
                    q_cand[M - 1] = 1.0 - sum_q;
                    auto q_stab = compute_expected_syn_weights(s, omegas, weights, syn_rates, q_cand, ecm_bl_ptr, delta, psi, eps_sink, &local_base);
                    if (evaluate_log_l(s, omegas, weights, syn_rates, q_stab, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base) >=
                        evaluate_log_l(s, omegas, weights, syn_rates, q2, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base)) {
                        syn_weights = q_stab;
                    } else {
                        syn_weights = q2;
                    }
                } else {
                    syn_weights = q2;
                }
            }
            renormalize_srv_unit_mean(s, nullptr);
        };

        // Helper: optimize omegas (and optional tree scale s) for arbitrary K >= 1
        // Note: all non-top classes k < K-1 are bounded in [0.0, 1.0] (range01, F-002),
        // and the top class K-1 is fixed to 1.0 (when constrained & K >= 2) or bounded in [1.0, max_omega] (unconstrained).
        auto optimize_omegas = [&](Scalar& cur_s, const std::vector<Scalar>* bl_ptr, bool opt_scale, int max_it) {
            Scalar max_omega = use_error_sink ? 100.0 : 10000.0;
            if (is_constrained) {
                if (K == 1) {
                    if (opt_scale) {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            return -evaluate_log_l(x[0], {x[1]}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {cur_s, std::clamp(omegas[0], 0.0, 1.0)}, {0.01, 0.0}, {50.0, 1.0}, 1e-4, max_it);
                        cur_s = opt.x[0];
                        omegas[0] = opt.x[1];
                    } else {
                        auto [w_opt, nll_opt] = Brent1D::minimize([&](Scalar w0) {
                            return -evaluate_log_l(1.0, {w0}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, 0.0, std::clamp(omegas[0], 0.0, 1.0), 1.0, 1e-4, max_it);
                        omegas[0] = w_opt;
                    }
                } else if (K == 2) {
                    if (opt_scale) {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            return -evaluate_log_l(x[0], {x[1], 1.0}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {cur_s, std::clamp(omegas[0], 0.0, 1.0)}, {0.01, 0.0}, {50.0, 1.0}, 1e-4, max_it);
                        cur_s = opt.x[0];
                        omegas[0] = opt.x[1];
                        omegas[1] = 1.0;
                    } else {
                        auto [w_opt, nll_opt] = Brent1D::minimize([&](Scalar w0) {
                            return -evaluate_log_l(1.0, {w0, 1.0}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, 0.0, std::clamp(omegas[0], 0.0, 1.0), 1.0, 1e-4, max_it);
                        omegas[0] = w_opt;
                        omegas[1] = 1.0;
                    }
                } else if (K == 3) {
                    if (opt_scale) {
                        auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                            Scalar w0 = std::min(x[1], x[2]);
                            Scalar w1 = std::max(x[1], x[2]);
                            return -evaluate_log_l(x[0], {w0, w1, 1.0}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {cur_s, std::clamp(omegas[0], 0.0, 1.0), std::clamp(omegas[1], 0.0, 1.0)}, {0.01, 0.0, 0.0}, {50.0, 1.0, 1.0}, 1e-4, max_it);
                        cur_s = opt.x[0];
                        omegas[0] = opt.x[1];
                        omegas[1] = opt.x[2];
                        omegas[2] = 1.0;
                    } else {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            Scalar w0 = std::min(x[0], x[1]);
                            Scalar w1 = std::max(x[0], x[1]);
                            return -evaluate_log_l(1.0, {w0, w1, 1.0}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {std::clamp(omegas[0], 0.0, 1.0), std::clamp(omegas[1], 0.0, 1.0)}, {0.0, 0.0}, {1.0, 1.0}, 1e-4, max_it);
                        omegas[0] = opt.x[0];
                        omegas[1] = opt.x[1];
                        omegas[2] = 1.0;
                    }
                    sort_rate_weight_pairs(omegas, weights, 2);
                } else {
                    size_t D = (K - 1) + (opt_scale ? 1 : 0);
                    std::vector<Scalar> x0(D), lb(D, 0.0), ub(D, 1.0);
                    size_t off = 0;
                    if (opt_scale) {
                        x0[0] = cur_s; lb[0] = 0.01; ub[0] = 50.0; off = 1;
                    }
                    for (size_t k = 0; k + 1 < K; ++k) {
                        x0[off + k] = std::clamp(omegas[k], 0.0, 1.0);
                    }
                    auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                        Scalar sc = opt_scale ? x[0] : 1.0;
                        std::vector<Scalar> cand_o(K, 1.0);
                        for (size_t k = 0; k + 1 < K; ++k) cand_o[k] = x[off + k];
                        return -evaluate_log_l(sc, cand_o, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                    }, x0, lb, ub, 1e-4, max_it);
                    if (opt_scale) cur_s = opt.x[0];
                    for (size_t k = 0; k + 1 < K; ++k) omegas[k] = opt.x[off + k];
                    omegas[K - 1] = 1.0;
                    sort_rate_weight_pairs(omegas, weights, K - 1);
                }
            } else {
                if (K == 1) {
                    if (opt_scale) {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            return -evaluate_log_l(x[0], {x[1]}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {cur_s, std::clamp(omegas[0], 0.0, max_omega)}, {0.01, 0.0}, {50.0, max_omega}, 1e-4, max_it);
                        cur_s = opt.x[0];
                        omegas[0] = opt.x[1];
                    } else {
                        auto [w_opt, nll_opt] = Brent1D::minimize([&](Scalar w0) {
                            return -evaluate_log_l(1.0, {w0}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, 0.0, std::clamp(omegas[0], 0.0, max_omega), max_omega, 1e-4, max_it);
                        omegas[0] = w_opt;
                    }
                } else if (K == 2) {
                    if (opt_scale) {
                        auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                            return -evaluate_log_l(x[0], {x[1], x[2]}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {cur_s, std::clamp(omegas[0], 0.0, 1.0), std::clamp(omegas[1], 1.0, max_omega)}, {0.01, 0.0, 1.0}, {50.0, 1.0, max_omega}, 1e-4, max_it);
                        cur_s = opt.x[0];
                        omegas[0] = opt.x[1];
                        omegas[1] = opt.x[2];
                    } else {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            return -evaluate_log_l(1.0, {x[0], x[1]}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {std::clamp(omegas[0], 0.0, 1.0), std::clamp(omegas[1], 1.0, max_omega)}, {0.0, 1.0}, {1.0, max_omega}, 1e-4, max_it);
                        omegas[0] = opt.x[0];
                        omegas[1] = opt.x[1];
                    }
                } else if (K == 3) {
                    if (opt_scale) {
                        auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                            Scalar w0 = std::min(x[1], x[2]);
                            Scalar w1 = std::max(x[1], x[2]);
                            return -evaluate_log_l(x[0], {w0, w1, x[3]}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {cur_s, std::clamp(omegas[0], 0.0, 1.0), std::clamp(omegas[1], 0.0, 1.0), std::clamp(omegas[2], 1.0, max_omega)},
                           {0.01, 0.0, 0.0, 1.0}, {50.0, 1.0, 1.0, max_omega}, 1e-4, max_it);
                        cur_s = opt.x[0];
                        omegas[0] = opt.x[1];
                        omegas[1] = opt.x[2];
                        omegas[2] = opt.x[3];
                    } else {
                        auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                            Scalar w0 = std::min(x[0], x[1]);
                            Scalar w1 = std::max(x[0], x[1]);
                            return -evaluate_log_l(1.0, {w0, w1, x[2]}, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                        }, {std::clamp(omegas[0], 0.0, 1.0), std::clamp(omegas[1], 0.0, 1.0), std::clamp(omegas[2], 1.0, max_omega)},
                           {0.0, 0.0, 1.0}, {1.0, 1.0, max_omega}, 1e-4, max_it);
                        omegas[0] = opt.x[0];
                        omegas[1] = opt.x[1];
                        omegas[2] = opt.x[2];
                    }
                    sort_rate_weight_pairs(omegas, weights, 2);
                } else {
                    size_t D = K + (opt_scale ? 1 : 0);
                    std::vector<Scalar> x0(D), lb(D, 0.0), ub(D, 1.0);
                    size_t off = 0;
                    if (opt_scale) {
                        x0[0] = cur_s; lb[0] = 0.01; ub[0] = 50.0; off = 1;
                    }
                    for (size_t k = 0; k + 1 < K; ++k) {
                        x0[off + k] = std::clamp(omegas[k], 0.0, 1.0);
                    }
                    x0[off + K - 1] = std::clamp(omegas[K - 1], 1.0, max_omega);
                    lb[off + K - 1] = 1.0;
                    ub[off + K - 1] = max_omega;
                    auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                        Scalar sc = opt_scale ? x[0] : 1.0;
                        std::vector<Scalar> cand_o(K);
                        for (size_t k = 0; k < K; ++k) cand_o[k] = x[off + k];
                        return -evaluate_log_l(sc, cand_o, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
                    }, x0, lb, ub, 1e-4, max_it);
                    if (opt_scale) cur_s = opt.x[0];
                    for (size_t k = 0; k < K; ++k) omegas[k] = opt.x[off + k];
                    sort_rate_weight_pairs(omegas, weights, K - 1);
                }
            }
        };

        auto optimize_mh = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr, int max_it) {
            if (allow_triple) {
                auto mh_opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                    return -evaluate_log_l(cur_s, omegas, weights, syn_rates, syn_weights, bl_ptr, nullptr, x[0], x[1], eps_sink, &local_base);
                }, { std::clamp(delta, 0.0, 1000.0), std::clamp(psi, 0.0, 1000.0) }, {0.0, 0.0}, {1000.0, 1000.0}, 1e-4, max_it);
                delta = std::max(0.0, mh_opt.x[0]);
                psi = std::max(0.0, mh_opt.x[1]);
            } else if (allow_double) {
                auto [d_opt, nll_opt] = Brent1D::minimize([&](Scalar d) {
                    return -evaluate_log_l(cur_s, omegas, weights, syn_rates, syn_weights, bl_ptr, nullptr, d, 0.0, eps_sink, &local_base);
                }, 0.0, std::clamp(delta, 0.0, 1000.0), 1000.0, 1e-4, max_it);
                delta = std::max(0.0, d_opt);
            }
        };

        // Helper: re-polish the 5 GTR nucleotide exchangeabilities (F-049)
        auto optimize_gtr_theta = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr) {
            Scalar cur_ll = evaluate_log_l(cur_s, omegas, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
            std::vector<Scalar> x0 = {
                std::log(std::clamp(local_base.theta_AC, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_AT, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_CG, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_CT, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_GT, 1e-4, 1e4))
            };
            std::vector<Scalar> lb(5, std::log(1e-4)), ub(5, std::log(1e4));
            auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                MG94Parameters cand = local_base;
                cand.theta_AC = std::exp(x[0]);
                cand.theta_AT = std::exp(x[1]);
                cand.theta_CG = std::exp(x[2]);
                cand.theta_CT = std::exp(x[3]);
                cand.theta_GT = std::exp(x[4]);
                return -evaluate_log_l(cur_s, omegas, weights, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink, &cand);
            }, x0, lb, ub, 1e-4, 40);
            if (-opt.f > cur_ll + 1e-5) {
                local_base.theta_AC = std::exp(opt.x[0]);
                local_base.theta_AT = std::exp(opt.x[1]);
                local_base.theta_CG = std::exp(opt.x[2]);
                local_base.theta_CT = std::exp(opt.x[3]);
                local_base.theta_GT = std::exp(opt.x[4]);
            }
        };

        // ECM Optimization Loop (convergence-terminated up to 20 iterations, F-001)
        int max_ecm = (K == 1 && !srv && !allow_double && !use_error_sink) ? 1 : 20;
        Scalar prev_ecm_ll = evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
        for (int ecm_iter = 0; ecm_iter < max_ecm; ++ecm_iter) {
            squarem_update_weights();
            squarem_update_syn_weights();
            optimize_omegas(s, ecm_bl_ptr, true, (K == 3 && !is_constrained) ? 80 : 50);

            if (ecm_iter == 1 || ecm_iter == 4) {
                try_rescue_collapsed_classes(s, omegas, weights, ecm_bl_ptr);
            }

            optimize_srv_rates(s, omegas, weights, nullptr, 40);
            optimize_mh(s, ecm_bl_ptr, 30);

            Scalar cur_ecm_ll = evaluate_log_l(s, omegas, weights, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink, &local_base);
            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << " (ECM " << (ecm_iter + 1) << "): lnL = "
                   << std::fixed << std::setprecision(2) << cur_ecm_ll;
                if (!is_constrained && !omegas.empty()) {
                    ss << ", ω" << K << " = " << std::setprecision(2) << omegas.back();
                }
                if (use_error_sink) {
                    ss << ", ε = " << std::setprecision(4) << eps_sink;
                }
                pbar->set_status(ss.str());
                pbar->tick();
            }
            if (ecm_iter >= 5 && std::abs(cur_ecm_ll - prev_ecm_ll) < 1e-3) {
                break;
            }
            prev_ecm_ll = cur_ecm_ll;
        }

        if (!refine_branches) {
            optimize_gtr_theta(s, ecm_bl_ptr);
        }

        // Branch length initialization and refinement (convergence-terminated up to 6 cycles, F-005)
        std::vector<Scalar> bl(tree.num_nodes(), 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                bl[node.id] = s * (ecm_bl_ptr ? (*ecm_bl_ptr)[node.id] : branch_tau[node.id]);
            }
        }

        if (refine_branches) {
            int max_cycles = (K == 1 && !srv && !allow_double && !use_error_sink) ? 2 : 6;
            Scalar prev_cycle_ll = evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink, &local_base);
            for (int cycle = 0; cycle < max_cycles; ++cycle) {
                bl = refine_branch_lengths(bl, omegas, weights, syn_rates, syn_weights, delta, psi, 25, eps_sink, &local_base);

                if (K >= 2 || use_error_sink) {
                    for (int em = 0; em < 2; ++em) {
                        weights = compute_expected_weights(1.0, omegas, weights, syn_rates, syn_weights, &bl, delta, psi, eps_sink, use_error_sink ? &eps_sink : nullptr, max_eps_sink, &local_base);
                        if (srv && M >= 2) {
                            syn_weights = compute_expected_syn_weights(1.0, omegas, weights, syn_rates, syn_weights, &bl, delta, psi, eps_sink, &local_base);
                        }
                    }
                    if (use_error_sink) {
                        Scalar best_eps = eps_sink;
                        Scalar best_eps_ll = evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink, &local_base);
                        for (Scalar cand_e : {1e-7, 2e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2}) {
                            if (cand_e > max_eps_sink) continue;
                            Scalar ll_e = evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, cand_e, &local_base);
                            if (ll_e > best_eps_ll + 1e-4) {
                                best_eps_ll = ll_e;
                                best_eps = cand_e;
                            }
                        }
                        eps_sink = best_eps;
                    }
                }

                // Polish omegas given refined branch lengths
                Scalar one_s = 1.0;
                optimize_omegas(one_s, &bl, false, 55);

                if (cycle == 0 || cycle == 2) {
                    try_rescue_collapsed_classes(one_s, omegas, weights, &bl);
                }

                // Polish SRV rates given bl and enforce unit mean for all M >= 2 (F-019, F-020, F-021)
                if (srv && M >= 2) {
                    optimize_srv_rates(one_s, omegas, weights, &bl, 40);
                }

                // Polish MH rates inside refinement loop so branch lengths and delta/psi co-adapt (F-006)
                if (allow_double) {
                    optimize_mh(1.0, &bl, 30);
                }

                // Re-polish the 5 GTR nucleotide exchangeabilities theta_AC..GT (F-049)
                optimize_gtr_theta(1.0, &bl);

                Scalar cur_cycle_ll = evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink, &local_base);
                if (cycle >= 2 && std::abs(cur_cycle_ll - prev_cycle_ll) < 1e-3) {
                    break;
                }
                prev_cycle_ll = cur_cycle_ll;
            }

            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Branch lengths refined";
                pbar->set_status(ss.str());
                pbar->tick();
            }
            if (allow_double && pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Multiple hits refined (δ=" << std::fixed << std::setprecision(3) << delta << ")";
                pbar->set_status(ss.str());
                pbar->tick();
            }
        }

        // Ensure ascending omega order in non-top classes and top class at K-1
        if (K >= 2) {
            sort_rate_weight_pairs(omegas, weights, K - 1);
        }
        if (srv && M >= 2) {
            Scalar dummy_s = s;
            renormalize_srv_unit_mean(dummy_s, refine_branches ? &bl : nullptr);
            if (!refine_branches) {
                s = dummy_s;
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        bl[node.id] = s * (ecm_bl_ptr ? (*ecm_bl_ptr)[node.id] : branch_tau[node.id]);
                    }
                }
            }
        }

        Scalar best_lnL = evaluate_log_l(1.0, omegas, weights, syn_rates, syn_weights, &bl, &fit.site_log_likelihoods, delta, psi, eps_sink, &local_base);
        fit.log_likelihood = best_lnL;
        fit.tree_scale = s;
        fit.branch_lengths = bl;
        fit.delta = delta;
        fit.psi = psi;
        fit.error_sink_weight = eps_sink;
        fit.has_gtr_params = true;
        fit.gtr_params = local_base;
        if (warm_start_fit && !warm_start_fit->mss_class_rates.empty()) {
            fit.mss_class_rates = warm_start_fit->mss_class_rates;
        }
        if (local_base.syn_pair_rates != nullptr) {
            fit.mss_syn_pair_rates = *local_base.syn_pair_rates;
            fit.gtr_params.syn_pair_rates = &fit.mss_syn_pair_rates;
        }

        // Compute substitution rates and fractions
        Scalar total_rate = 0.0;
        Scalar double_rate = 0.0;
        Scalar triple_rate = 0.0;
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = local_base;
            p.alpha = 1.0;
            p.beta = omegas[k];
            p.delta = delta;
            p.psi = psi;
            MG94Matrix mat;
            mat.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            auto sr = mat.compute_substitution_rates(*code);
            total_rate += weights[k] * sr.total_rate;
            double_rate += weights[k] * sr.rate_double;
            triple_rate += weights[k] * sr.rate_triple;
        }
        fit.frac_delta = (total_rate > 1e-12) ? (double_rate / total_rate) : 0.0;
        fit.frac_psi = (total_rate > 1e-12) ? (triple_rate / total_rate) : 0.0;

        fit.test_distribution.omegas = omegas;
        fit.test_distribution.weights = weights;
        fit.test_distribution.error_sink_weight = eps_sink;
        if (srv) {
            fit.test_distribution.syn_rates = syn_rates;
            fit.test_distribution.syn_weights = syn_weights;
        }

        // Detect collapsed rate classes. annotations[k] annotates class k
        // directly (empty string = unannotated); a pairwise collapse note is
        // attached to the first class of the pair.
        fit.test_distribution.annotations.assign(K, "");
        for (size_t i = 0; i < K; ++i) {
            if (weights[i] < 0.005) {
                fit.test_distribution.annotations[i] = "Rate class " + std::to_string(i) + " collapsed (weight < 0.5%)";
            }
        }
        for (size_t i = 0; i + 1 < K; ++i) {
            if (std::abs(omegas[i+1] - omegas[i]) < 0.015 && fit.test_distribution.annotations[i].empty()) {
                fit.test_distribution.annotations[i] = "Rate classes " + std::to_string(i) + " and " + std::to_string(i+1) + " collapsed into identical omega";
            }
        }

        size_t N_codons = aln.num_codons;
        size_t num_branches = tree.num_nodes() - 1;
        size_t num_mixture_params = (2 * K - 1) - (is_constrained && K >= 2 ? 1 : 0) + (use_error_sink ? 1 : 0);
        size_t num_syn_params = srv ? 2 * (M - 1) : 0;
        size_t num_mh_params = (allow_double ? 1 : 0) + (allow_triple ? 1 : 0);
        // 5 GTR exchangeabilities + 9 CF3x4 empirical frequency parameters (matching HyPhy's +9 AICc offset, F-054)
        size_t K_params = 14 + (refine_branches ? num_branches : 1) + num_mixture_params + num_syn_params + num_mh_params + num_mss_free_params;
        fit.aicc = -2.0 * fit.log_likelihood + 2.0 * K_params * N_codons / std::max(1.0, static_cast<double>(N_codons) - static_cast<double>(K_params) - 1.0);

        return fit;
    }

    // Computes branch-site Error Sink Bayes Factors in a single O(patterns * branches) Inside-Outside pass
    void compute_error_sink_diagnostics(
        const BUSTEDFit& fit,
        Scalar bf_threshold,
        std::vector<Scalar>& out_site_max_bf,
        std::vector<BUSTEDErrorEvent>& out_flagged_events
    ) const {
        size_t num_sites = aln.num_codons;
        out_site_max_bf.assign(num_sites, 0.0);
        out_flagged_events.clear();

        Scalar eps = fit.error_sink_weight;
        if (eps <= 1e-9 || eps >= 1.0 - 1e-9) return;
        Scalar eff_prior = std::max(eps, 0.0005);
        Scalar prior_odds = eff_prior / (1.0 - eff_prior);
        Scalar bio_scale = 1.0 - eps;

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K = fit.test_distribution.omegas.size();
        std::vector<Scalar> syn_rates = fit.test_distribution.syn_rates.empty() ? std::vector<Scalar>{1.0} : fit.test_distribution.syn_rates;
        std::vector<Scalar> syn_weights = fit.test_distribution.syn_weights.empty() ? std::vector<Scalar>{1.0} : fit.test_distribution.syn_weights;
        size_t M = syn_rates.size();

        const MG94Parameters& active_base = fit.has_gtr_params ? fit.gtr_params : base_params;
        std::vector<MG94Matrix> mg(K);
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = fit.test_distribution.omegas[k];
            p.delta = fit.delta;
            p.psi = fit.psi;
            mg[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M, std::vector<std::array<Matrix, 3>>(tree.bridges.size()));
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar bl = fit.branch_lengths[node.id];
                    P_mix[m][node.id] = Matrix::Zero(S, S);
                    Scalar use_bio = is_test_branch[node.id] ? bio_scale : 1.0;
                    for (size_t k = 0; k < K; ++k) {
                        P_mix[m][node.id] += (use_bio * fit.test_distribution.weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                    }
                    if (is_test_branch[node.id]) {
                        P_mix[m][node.id].noalias() += eps * P_sink;
                    }
                }
            }
            if (!tree.bridges.empty()) {
                for (size_t b = 0; b < tree.bridges.size(); ++b) {
                    int bid = tree.bridges[b].bridge_node_id;
                    bool is_test = (bid >= 0 && static_cast<size_t>(bid) < is_test_branch.size()) ? is_test_branch[bid] : true;
                    Scalar use_bio = is_test ? bio_scale : 1.0;
                    for (int r = 0; r < 3; ++r) {
                        Scalar bl = bridge_resolution_tau(b, r, 1.0, &fit.branch_lengths);
                        P_bridges_mix[m][b][r] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K; ++k) {
                            P_bridges_mix[m][b][r] += (use_bio * fit.test_distribution.weights[k]) * mg[k].transition_matrix(alpha_m * bl);
                        }
                        if (is_test) {
                            P_bridges_mix[m][b][r].noalias() += eps * P_sink;
                        }
                    }
                }
            }
        }

        // Per-pattern, per-node posterior probability of the saturation sink (underflow-safe via log-likelihood)
        std::vector<std::vector<Scalar>> pat_node_post(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        const int n_threads_sink = LikelihoodEngine::optimal_omp_threads(num_patterns);
        #pragma omp parallel num_threads(n_threads_sink)
        {
            std::vector<LikelihoodEngine::InsideOutsideWorkspace> io_m(M);
            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    LikelihoodEngine::compute_inside_outside_in_place(
                        tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, io_m[m],
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                    );
                    if (io_m[m].likelihood > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                        if (term > max_ll) max_ll = term;
                    }
                }
                if (max_ll > -1e19) {
                    Scalar sum_exp = 0.0;
                    std::vector<Scalar> post_m(M, 0.0);
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                            Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood - max_ll);
                            post_m[m] = e;
                            sum_exp += e;
                        }
                    }
                    if (sum_exp > 0.0) {
                        for (size_t m = 0; m < M; ++m) post_m[m] /= sum_exp;
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id && is_test_branch[node.id]) {
                                int nid = node.id;
                                Scalar sink_ratio = 0.0;
                                for (size_t m = 0; m < M; ++m) {
                                    if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                        sink_ratio += post_m[m] * (eval_sink_bilinear(io_m[m].V[nid], io_m[m].D[nid]) / io_m[m].likelihood);
                                    }
                                }
                                Scalar post = std::clamp(eps * sink_ratio, 0.0, 1.0 - 1e-9);
                                pat_node_post[p][nid] = post;
                            }
                        }
                    }
                }
            }
        }

        for (size_t s_idx = 0; s_idx < num_sites; ++s_idx) {
            size_t p_idx = aln.site_to_pattern[s_idx];
            Scalar max_bf = 0.0;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && is_test_branch[node.id]) {
                    Scalar post = pat_node_post[p_idx][node.id];
                    Scalar post_odds = post / std::max(1e-6, 1.0 - post);
                    Scalar bf = post_odds / prior_odds;
                    if (bf > max_bf) max_bf = bf;
                    if (bf >= bf_threshold && post >= 0.50) {
                        out_flagged_events.push_back({s_idx + 1, node.name, post, bf});
                    }
                }
            }
            out_site_max_bf[s_idx] = max_bf;
        }
    }

    // Fast 4-Gate Surprisal & Score Automatic Model Selector (BUSTED-Auto)
    // Diagnoses Error Sink (O(1) saturation + spatial adjacency), Multi-Hits (O(t_b) quarantined tangent score),
    // Synonymous Rate Variation (Neyman-Scott overdispersion + SQUAREM probe), and Omega classes K in a single baseline pass.
    BUSTEDAutoDiagnostics diagnose_and_select_model(BUSTEDSettings& settings) const {
        BUSTEDAutoDiagnostics diag;
        diag.enabled = true;

        size_t num_sites = aln.num_codons;
        size_t num_patterns = aln.patterns.size();
        size_t num_nodes = tree.num_nodes();
        int S = code->num_sense_codons;

        // Step 0: Fast K=2 baseline warm-up (M0: M=1, delta=0, psi=0, eps=0)
        std::vector<Scalar> omegas0 = {std::min(0.05, 0.25 * mg94_omega), 1.0};
        std::vector<Scalar> weights0 = {0.80, 0.20};
        Scalar s0 = 1.0;
        for (int it = 0; it < 2; ++it) {
            weights0 = compute_expected_weights(s0, omegas0, weights0);
            auto opt0 = NelderMead2D::minimize([&](Scalar w1, Scalar sc) {
                return -evaluate_log_l(sc, {w1, 1.0}, weights0);
            }, omegas0[0], s0, 1e-4, 1.0, 1e-3, 20);
            omegas0[0] = opt0.a;
            s0 = opt0.b;
        }
        weights0 = compute_expected_weights(s0, omegas0, weights0);
        Scalar base_ll = evaluate_log_l(s0, omegas0, weights0);
        Scalar mean_omega0 = weights0[0] * omegas0[0] + weights0[1] * omegas0[1];

        // Build baseline branch transition matrices P0[node_id]
        std::vector<MG94Matrix> mg0(2);
        for (size_t k = 0; k < 2; ++k) {
            MG94Parameters p = base_params;
            p.alpha = 1.0;
            p.beta = omegas0[k];
            mg0[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        std::vector<Matrix> P0(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar bl = s0 * branch_tau[node.id];
                P0[node.id] = weights0[0] * mg0[0].transition_matrix(bl) + weights0[1] * mg0[1].transition_matrix(bl);
            }
        }
        std::vector<std::array<Matrix, 3>> P0_bridges(tree.bridges.size());
        for (size_t b = 0; b < tree.bridges.size(); ++b) {
            for (int r = 0; r < 3; ++r) {
                Scalar bl = bridge_resolution_tau(b, r, s0);
                P0_bridges[b][r] = weights0[0] * mg0[0].transition_matrix(bl) + weights0[1] * mg0[1].transition_matrix(bl);
            }
        }

        // Build pure 2-hit (Q2_off), 3-hit (Q3_off), and synonymous (Q_syn) tangent operators at mean_omega0
        Scalar r_nuc[4][4] = {
            {0.0, base_params.theta_AC, 1.0, base_params.theta_AT},
            {base_params.theta_AC, 0.0, base_params.theta_CG, base_params.theta_CT},
            {1.0, base_params.theta_CG, 0.0, base_params.theta_GT},
            {base_params.theta_AT, base_params.theta_CT, base_params.theta_GT, 0.0}
        };
        Scalar q_scale = std::max(1e-6, mg0[1].scale_factor);
        Matrix Q2_off = Matrix::Zero(S, S);
        Matrix Q3_off = Matrix::Zero(S, S);
        Matrix Q_syn = Matrix::Zero(S, S);
        for (int i = 0; i < S; ++i) {
            Scalar syn_row_sum = 0.0;
            for (int j = 0; j < S; ++j) {
                if (i == j) continue;
                const auto& diff = code->diff_matrix[i][j];
                Scalar mod = diff.synonymous ? 1.0 : mean_omega0;
                if (diff.num_diffs == 1 && diff.synonymous) {
                    Scalar rate1 = r_nuc[diff.diff_from[0]][diff.diff_to[0]] * aln.pos_nuc_frequencies(diff.diff_pos[0], diff.diff_to[0]) / q_scale;
                    Q_syn(i, j) = rate1;
                    syn_row_sum += rate1;
                } else if (diff.num_diffs == 2) {
                    Scalar rate2 = mod * r_nuc[diff.diff_from[0]][diff.diff_to[0]] * r_nuc[diff.diff_from[1]][diff.diff_to[1]] *
                                   aln.pos_nuc_frequencies(diff.diff_pos[0], diff.diff_to[0]) * aln.pos_nuc_frequencies(diff.diff_pos[1], diff.diff_to[1]) / q_scale;
                    Q2_off(i, j) = rate2;
                } else if (diff.num_diffs == 3) {
                    Scalar rate3 = mod * r_nuc[diff.diff_from[0]][diff.diff_to[0]] * r_nuc[diff.diff_from[1]][diff.diff_to[1]] * r_nuc[diff.diff_from[2]][diff.diff_to[2]] *
                                   aln.pos_nuc_frequencies(diff.diff_pos[0], diff.diff_to[0]) * aln.pos_nuc_frequencies(diff.diff_pos[1], diff.diff_to[1]) * aln.pos_nuc_frequencies(diff.diff_pos[2], diff.diff_to[2]) / q_scale;
                    Q3_off(i, j) = rate3;
                }
            }
            Q_syn(i, i) = -syn_row_sum;
        }

        // Single parallel Inside-Outside pass on M0
        std::vector<std::vector<Scalar>> pat_R_sink(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        std::vector<std::vector<Scalar>> pat_T2(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        std::vector<std::vector<Scalar>> pat_T3(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        std::vector<Scalar> pat_syn_score(num_patterns, 0.0);
        std::vector<Scalar> pat_syn_info(num_patterns, 0.0);

        const int n_threads_auto = LikelihoodEngine::optimal_omp_threads(num_patterns);
        #pragma omp parallel num_threads(n_threads_auto)
        {
            LikelihoodEngine::InsideOutsideWorkspace io_ws;
            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                LikelihoodEngine::compute_inside_outside_in_place(
                    tree, pattern, leaf_to_taxon, P0, aln.codon_frequencies_f3x4, S, io_ws,
                    tree.bridges.empty() ? nullptr : &P0_bridges
                );
                Scalar L_p = io_ws.likelihood;
                if (L_p > 0.0) {
                    Scalar inv_L = 1.0 / L_p;
                    Scalar u_syn = 0.0;
                    Scalar v_syn = 0.0;
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id) {
                            int nid = node.id;
                            Scalar bl = s0 * branch_tau[nid];
                            if (is_test_branch[nid]) {
                                Scalar sink_val = eval_sink_bilinear(io_ws.V[nid], io_ws.D[nid]);
                                pat_R_sink[p][nid] = sink_val * inv_L;
                                pat_T2[p][nid] = bl * io_ws.V[nid].dot(Q2_off * io_ws.D[nid]) * inv_L;
                                pat_T3[p][nid] = bl * io_ws.V[nid].dot(Q3_off * io_ws.D[nid]) * inv_L;
                            }
                            Scalar branch_u = bl * io_ws.V[nid].dot(Q_syn * (P0[nid] * io_ws.D[nid])) * inv_L;
                            u_syn += branch_u;
                            v_syn += branch_u * branch_u;
                        }
                    }
                    pat_syn_score[p] = u_syn;
                    pat_syn_info[p] = std::max(1e-6, v_syn);
                }
            }
        }

        // Gate 1: Evaluate Error Sink Surprisal & Spatial Adjacency
        Scalar max_R_sink = 0.0;
        int max_sink_nid = -1;
        size_t candidate_count = 0;
        size_t adjacent_clusters = 0;
        Scalar eps_init = 0.002;
        Scalar sink_resp_sum = 0.0;
        Scalar total_branch_sites = static_cast<Scalar>(num_sites * std::max<size_t>(1, num_test_branches));

        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id && is_test_branch[node.id]) {
                int nid = node.id;
                int prev_high_site = -10;
                for (size_t s_idx = 0; s_idx < num_sites; ++s_idx) {
                    size_t p_idx = aln.site_to_pattern[s_idx];
                    Scalar r_val = pat_R_sink[p_idx][nid];
                    if (r_val > max_R_sink) {
                        max_R_sink = r_val;
                        max_sink_nid = nid;
                    }
                    Scalar bits = (r_val > 1.0) ? std::log2(r_val) : 0.0;
                    if (bits > 8.0) {
                        candidate_count++;
                    }
                    if (bits > 6.5) {
                        if (static_cast<int>(s_idx) - prev_high_site <= 2) {
                            adjacent_clusters++;
                        }
                        prev_high_site = static_cast<int>(s_idx);
                    }
                    Scalar resp = (eps_init * r_val) / ((1.0 - eps_init) + eps_init * r_val);
                    sink_resp_sum += resp;
                }
            }
        }

        diag.max_sink_bf = max_R_sink;
        diag.max_sink_surprisal_bits = (max_R_sink > 1.0) ? std::log2(max_R_sink) : 0.0;
        diag.sink_candidate_sites = candidate_count;
        diag.sink_adjacent_clusters = adjacent_clusters;

        Scalar eps_probe = std::clamp(sink_resp_sum / total_branch_sites, 1e-6, settings.error_sink_max_weight);
        Scalar sink_probe_ll = base_ll;
        if (diag.max_sink_surprisal_bits > 6.5) {
            sink_probe_ll = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, 0.0, 0.0, eps_probe);
            // If a single branch has extreme surprisal (which inflates its MG94 branch length), also check 1-step branch refinement under sink
            if (max_sink_nid >= 0 && diag.max_sink_surprisal_bits >= 10.0) {
                std::vector<Scalar> probe_bl(num_nodes, 0.0);
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) probe_bl[node.id] = s0 * branch_tau[node.id];
                }
                auto ref_bl = refine_branch_lengths(probe_bl, omegas0, weights0, {1.0}, {1.0}, 0.0, 0.0, 4, eps_probe);
                auto base_ref_bl = refine_branch_lengths(probe_bl, omegas0, weights0, {1.0}, {1.0}, 0.0, 0.0, 4, 0.0);
                Scalar ref_sink_ll = evaluate_log_l(1.0, omegas0, weights0, {1.0}, {1.0}, &ref_bl, nullptr, 0.0, 0.0, eps_probe);
                Scalar ref_base_ll = evaluate_log_l(1.0, omegas0, weights0, {1.0}, {1.0}, &base_ref_bl, nullptr, 0.0, 0.0, 0.0);
                if (ref_sink_ll - ref_base_ll > sink_probe_ll - base_ll) {
                    sink_probe_ll = base_ll + (ref_sink_ll - ref_base_ll);
                }
                // Recompute max_R_sink on the un-inflated branch length
                if (ref_bl[max_sink_nid] < probe_bl[max_sink_nid] * 0.5 && ref_bl[max_sink_nid] > 1e-8) {
                    Scalar ratio_gain = std::pow(probe_bl[max_sink_nid] / ref_bl[max_sink_nid], 2.0);
                    max_R_sink *= ratio_gain;
                    diag.max_sink_bf = max_R_sink;
                    diag.max_sink_surprisal_bits = (max_R_sink > 1.0) ? std::log2(max_R_sink) : 0.0;
                }
            }
        }
        diag.sink_delta_ll = std::max(0.0, sink_probe_ll - base_ll);

        // Activate Error Sink if Delta LL > 2.2 AND either extreme short-branch surprisal (>= 13.0 bits, BF >= 8192) or spatial adjacency cluster (>= 8.0 bits)
        if (diag.sink_delta_ll > 2.2 && (diag.max_sink_surprisal_bits >= 13.0 || (diag.sink_adjacent_clusters >= 1 && diag.max_sink_surprisal_bits >= 8.0))) {
            diag.activated_error_sink = true;
        }

        // Gate 2: Evaluate Multi-Hit Tangent Score & Probe (Quarantining high-R_sink branch-sites if Error Sink is active)
        Scalar eps_active = diag.activated_error_sink ? eps_probe : 0.0;
        Scalar s_delta = 0.0;
        Scalar s_psi = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            Scalar w_p = aln.patterns[p].weight;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && is_test_branch[node.id]) {
                    int nid = node.id;
                    Scalar r_val = pat_R_sink[p][nid];
                    Scalar r_sink = (eps_active > 0.0) ? ((eps_active * r_val) / ((1.0 - eps_active) + eps_active * r_val)) : 0.0;
                    Scalar clean_wt = w_p * (1.0 - r_sink);
                    s_delta += clean_wt * pat_T2[p][nid];
                    s_psi += clean_wt * pat_T3[p][nid];
                }
            }
        }
        diag.double_hit_score = s_delta;
        diag.triple_hit_score = s_psi;

        Scalar mh_base_ll = diag.activated_error_sink ? sink_probe_ll : base_ll;
        if (s_delta > 0.15) {
            Scalar ll_d1 = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, 0.06, 0.0, eps_active);
            Scalar ll_d2 = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, 0.20, 0.0, eps_active);
            Scalar best_d = (ll_d2 > ll_d1) ? 0.20 : 0.06;
            Scalar best_d_ll = std::max(ll_d1, ll_d2);
            diag.double_hit_delta_ll = std::max(0.0, best_d_ll - mh_base_ll);

            // Check if diffuse 2-hit evolution (with eps=0) strictly outperforms an isolated non-clustered error sink
            if (diag.activated_error_sink && diag.sink_adjacent_clusters == 0 && diag.max_sink_surprisal_bits < 16.0) {
                Scalar pure_d_ll = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, best_d, 0.0, 0.0);
                if (pure_d_ll > sink_probe_ll + 0.5 && (pure_d_ll - base_ll) > 2.2) {
                    diag.activated_error_sink = false;
                    eps_active = 0.0;
                    mh_base_ll = base_ll;
                    best_d_ll = pure_d_ll;
                    diag.double_hit_delta_ll = pure_d_ll - base_ll;
                }
            }

            if (diag.double_hit_delta_ll > 2.2) {
                diag.selected_multiple_hits = "Double";
                if (s_psi > 0.05) {
                    Scalar ll_p = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, best_d, 0.06, eps_active);
                    diag.triple_hit_delta_ll = std::max(0.0, ll_p - best_d_ll);
                    if (diag.triple_hit_delta_ll > 2.2) {
                        diag.selected_multiple_hits = "Double+Triple";
                    }
                }
            }
        }

        // Gate 3: Evaluate Synonymous Rate Variation (SRV) Overdispersion & Probe (F-022)
        Scalar num_disp = 0.0;
        Scalar den_disp = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            Scalar w_p = aln.patterns[p].weight;
            Scalar u2 = pat_syn_score[p] * pat_syn_score[p];
            Scalar v = pat_syn_info[p];
            num_disp += w_p * (u2 - v);
            den_disp += w_p * v * v;
        }
        diag.srv_overdispersion_z = (den_disp > 1e-12) ? (num_disp / std::sqrt(2.0 * den_disp)) : 0.0;

        Scalar active_d = (diag.selected_multiple_hits != "None") ? 0.10 : 0.0;
        Scalar active_p = (diag.selected_multiple_hits == "Double+Triple") ? 0.04 : 0.0;
        Scalar pre_srv_ll = evaluate_log_l(s0, omegas0, weights0, {1.0}, {1.0}, nullptr, nullptr, active_d, active_p, eps_active);

        std::vector<Scalar> syn_r3 = {0.30, 1.00, 3.10};
        std::vector<Scalar> syn_w3 = {0.30, 0.60, 0.10};
        Scalar s_srv3 = s0;
        for (int it = 0; it < 2; ++it) {
            syn_w3 = compute_expected_syn_weights(s_srv3, omegas0, weights0, syn_r3, syn_w3, nullptr, active_d, active_p, eps_active);
            Scalar mean_r3 = syn_w3[0] * syn_r3[0] + syn_w3[1] * syn_r3[1] + syn_w3[2] * syn_r3[2];
            if (mean_r3 > 1e-8) {
                for (Scalar& r : syn_r3) r /= mean_r3;
                s_srv3 *= mean_r3;
            }
        }
        auto syn3_opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
            Scalar r0 = std::max(1e-3, x[0]);
            Scalar r1 = std::max(1e-3, x[1]);
            Scalar denom = std::max(1e-12, syn_w3[0] * r0 + syn_w3[1] * r1 + syn_w3[2]);
            return -evaluate_log_l(s_srv3, omegas0, weights0, {r0 / denom, r1 / denom, 1.0 / denom}, syn_w3, nullptr, nullptr, active_d, active_p, eps_active);
        }, {syn_r3[0] / std::max(1e-6, syn_r3[2]), syn_r3[1] / std::max(1e-6, syn_r3[2])}, {1e-3, 1e-3}, {50.0, 50.0}, 1e-3, 20);
        {
            Scalar r0 = std::max(1e-3, syn3_opt.x[0]);
            Scalar r1 = std::max(1e-3, syn3_opt.x[1]);
            Scalar denom = std::max(1e-12, syn_w3[0] * r0 + syn_w3[1] * r1 + syn_w3[2]);
            std::vector<Scalar> cand_r3 = {r0 / denom, r1 / denom, 1.0 / denom};
            if (-syn3_opt.f > evaluate_log_l(s_srv3, omegas0, weights0, syn_r3, syn_w3, nullptr, nullptr, active_d, active_p, eps_active)) {
                syn_r3 = cand_r3;
            }
        }
        Scalar srv3_ll = evaluate_log_l(s_srv3, omegas0, weights0, syn_r3, syn_w3, nullptr, nullptr, active_d, active_p, eps_active);
        diag.srv_delta_ll = std::max(0.0, srv3_ll - pre_srv_ll);

        if ((diag.srv_delta_ll > 4.5 || (diag.srv_overdispersion_z >= 2.0 && diag.srv_delta_ll > 3.0)) &&
            syn_w3[0] >= 0.02 && syn_w3[1] >= 0.02 && syn_w3[2] >= 0.02) {
            diag.activated_srv = true;
            diag.selected_syn_rates = 3;
        } else if (diag.srv_delta_ll > 2.2 || (diag.srv_overdispersion_z >= 1.65 && diag.srv_delta_ll > 1.5)) {
            std::vector<Scalar> syn_r2 = {0.35, 1.65};
            std::vector<Scalar> syn_w2 = {0.50, 0.50};
            Scalar s_srv2 = s0;
            syn_w2 = compute_expected_syn_weights(s_srv2, omegas0, weights0, syn_r2, syn_w2, nullptr, active_d, active_p, eps_active);
            Scalar mean_r2 = syn_w2[0] * syn_r2[0] + syn_w2[1] * syn_r2[1];
            if (mean_r2 > 1e-8) {
                syn_r2[0] /= mean_r2; syn_r2[1] /= mean_r2; s_srv2 *= mean_r2;
            }
            auto [r2_opt, nll2] = Brent1D::minimize([&](Scalar ratio) {
                Scalar denom = std::max(1e-12, syn_w2[0] * ratio + syn_w2[1]);
                return -evaluate_log_l(s_srv2, omegas0, weights0, {ratio / denom, 1.0 / denom}, syn_w2, nullptr, nullptr, active_d, active_p, eps_active);
            }, 1e-3, syn_r2[0] / std::max(1e-6, syn_r2[1]), 50.0, 1e-3, 20);
            Scalar srv2_ll = std::max(-nll2, evaluate_log_l(s_srv2, omegas0, weights0, syn_r2, syn_w2, nullptr, nullptr, active_d, active_p, eps_active));
            if ((srv2_ll - pre_srv_ll > 2.2 || (diag.srv_overdispersion_z >= 1.65 && srv2_ll - pre_srv_ll > 1.5)) &&
                syn_w2[0] >= 0.02 && syn_w2[1] >= 0.02) {
                diag.activated_srv = true;
                diag.selected_syn_rates = 2;
                diag.srv_delta_ll = std::max(diag.srv_delta_ll, srv2_ll - pre_srv_ll);
            }
        }

        // Gate 4: Fast Unrefined Null Screen for Optimal Omega Classes K in {1, 2, 3} using AICc step-up (F-054)
        size_t max_screen_k = std::max<size_t>(3, settings.max_k);
        auto [opt_screen_k, screen_null_fits] = select_optimal_k(
            max_screen_k, false, diag.activated_srv, diag.selected_syn_rates,
            diag.selected_multiple_hits, nullptr, diag.activated_error_sink, settings.error_sink_max_weight
        );
        // For BUSTED mixture testing, require at least K=2 unless user explicitly set max_k=1
        size_t chosen_k = std::max<size_t>(2, opt_screen_k);
        if (chosen_k == 3 && screen_null_fits.size() >= 3) {
            const auto& nk3 = screen_null_fits[2];
            bool k3_collapsed = (std::abs(nk3.test_distribution.omegas[1] - nk3.test_distribution.omegas[0]) < 0.025) ||
                                (nk3.test_distribution.weights[0] < 0.015) ||
                                (nk3.test_distribution.weights[1] < 0.015);
            if (k3_collapsed) {
                chosen_k = 2;
                diag.k_rationale = "K=2 selected (pruned collapsed/low-weight 3rd background rate class)";
            } else {
                diag.k_rationale = "K=3 retained via AICc step-up (delta lnL = +" +
                                   std::to_string(nk3.log_likelihood - screen_null_fits[1].log_likelihood) + ")";
            }
        } else {
            diag.k_rationale = "K=" + std::to_string(chosen_k) + " selected via AICc step-up parsimony";
        }
        diag.selected_k = chosen_k;

        // Apply selected specification to settings
        settings.error_sink = diag.activated_error_sink;
        settings.multiple_hits = diag.selected_multiple_hits;
        settings.srv = diag.activated_srv;
        settings.num_syn_rate_classes = diag.selected_syn_rates;
        settings.num_rate_classes = diag.selected_k;
        settings.auto_select_k = false;

        return diag;
    }

    // Step-up model selection across K in {1, ..., max_k} using AICc on the Null model
    std::pair<size_t, std::vector<BUSTEDFit>> select_optimal_k(
        size_t max_k = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01
    ) const {
        std::vector<BUSTEDFit> null_fits;
        null_fits.reserve(max_k);
        for (size_t k = 1; k <= max_k; ++k) {
            std::string desc = "Null K=" + std::to_string(k);
            if (num_background_branches > 0) {
                null_fits.push_back(fit_two_set_model(true, false, k, refine_branches, srv, num_syn_rates, multiple_hits, pbar, desc, use_error_sink, max_eps_sink));
            } else {
                null_fits.push_back(fit_model(true, k, refine_branches, srv, num_syn_rates, multiple_hits, pbar, desc, use_error_sink, max_eps_sink));
            }
        }

        size_t optimal_k = 1;
        for (size_t k = 2; k <= max_k; ++k) {
            // Delta AICc threshold: improvement must be at least 2.0 to justify adding a rate class
            if (null_fits[k - 1].aicc < null_fits[optimal_k - 1].aicc - 2.0) {
                optimal_k = k;
            }
        }
        return {optimal_k, null_fits};
    }

    // Two-set likelihood evaluation for BUSTED with separate Test (partition 0), Background/Comparison (partition 1), and optional Unclassified/Nuisance (partition 2) K-class distributions
    Scalar evaluate_two_set_log_l(
        Scalar s,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        std::vector<Scalar>* out_pattern_ll = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0,
        const std::vector<Scalar>* omegas_unclass = nullptr,
        const std::vector<Scalar>* weights_unclass = nullptr,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        size_t M = syn_rates.size();
        Scalar bio_scale_t = 1.0 - eps_sink_test;
        Scalar bio_scale_b = 1.0 - eps_sink_bg;

        std::vector<MG94Matrix> mg_t(K_t);
        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas_test[k];
            p.delta = delta;
            p.psi = psi;
            mg_t[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<MG94Matrix> mg_b(K_b);
        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas_bg[k];
            p.delta = delta;
            p.psi = psi;
            mg_b[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        bool has_unclass_dist = (num_unclassified_branches > 0 && omegas_unclass != nullptr && weights_unclass != nullptr && !omegas_unclass->empty() && omegas_unclass->size() == weights_unclass->size());
        size_t K_u = has_unclass_dist ? omegas_unclass->size() : 0;
        std::vector<MG94Matrix> mg_u(K_u);
        MG94Matrix mg_unclass;
        if (has_unclass_dist) {
            for (size_t k = 0; k < K_u; ++k) {
                MG94Parameters p_u = active_base;
                p_u.alpha = 1.0;
                p_u.beta = (*omegas_unclass)[k];
                p_u.delta = delta;
                p_u.psi = psi;
                mg_u[k].update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
        } else if (num_unclassified_branches > 0) {
            MG94Parameters p_u = active_base;
            p_u.alpha = 1.0;
            p_u.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega;
            p_u.delta = delta;
            p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[nid] : branch_tau[nid]);
                    int part = (static_cast<size_t>(nid) < branch_partition.size()) ? branch_partition[nid] : (is_test_branch[nid] ? 0 : 1);
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_test > 0.0) ? bio_scale_t : 1.0;
                        for (size_t k = 0; k < K_t; ++k) {
                            P_mix[m][nid] += (use_bio * weights_test[k]) * mg_t[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_test > 0.0) {
                            P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                        }
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_bg > 0.0) ? bio_scale_b : 1.0;
                        for (size_t k = 0; k < K_b; ++k) {
                            P_mix[m][nid] += (use_bio * weights_bg[k]) * mg_b[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_bg > 0.0) {
                            P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                        }
                    } else if (has_unclass_dist) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_u; ++k) {
                            P_mix[m][nid] += (*weights_unclass)[k] * mg_u[k].transition_matrix(alpha_m * bl);
                        }
                    } else {
                        P_mix[m][nid] = mg_unclass.transition_matrix(alpha_m * bl);
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_two_set_bridge_mix(
                    s, alpha_m, mg_t, mg_b, mg_u, mg_unclass,
                    weights_test, weights_bg, weights_unclass, has_unclass_dist,
                    eps_sink_test, eps_sink_bg, bio_scale_t, bio_scale_b, custom_branch_lengths
                );
            }
        }

        if (out_pattern_ll) {
            out_pattern_ll->assign(num_patterns, 0.0);
        }

        std::vector<Scalar> pat_ll_buf(num_patterns, -1e20);
        const int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns);
        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));
            std::vector<double> node_log_scale(num_nodes, 0.0);
            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar pll = -1e20;
                if (M == 1) {
                    const auto* p_br = tree.bridges.empty() ? nullptr : &P_bridges_mix[0];
                    for (int32_t node_id : tree.post_order) {
                        const auto& node = tree.nodes[node_id];
                        if (node.is_leaf) {
                            size_t t_idx = leaf_to_taxon[node_id];
                            if (t_idx != static_cast<size_t>(-1)) {
                                int8_t state = pattern.states[t_idx];
                                if (state >= 0 && state < S) {
                                    node_L[node_id].setZero();
                                    node_L[node_id](state) = 1.0;
                                } else {
                                    node_L[node_id].setOnes();
                                }
                            } else {
                                node_L[node_id].setOnes();
                            }
                            node_log_scale[node_id] = 0.0;
                        } else {
                            double child_scale_sum = 0.0;
                            tree.eval_internal_node(node_id, P_mix[0], node_L, node_log_scale, node_L[node_id], child_scale_sum, p_br);
                            node_log_scale[node_id] = child_scale_sum;
                        }
                    }
                    Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    if (cond_L > 0.0) {
                        pll = std::log(cond_L) + node_log_scale[tree.root_id];
                    }
                } else {
                    Scalar max_term = -std::numeric_limits<Scalar>::infinity();
                    std::vector<Scalar> terms(M);
                    std::vector<bool> valid(M, false);

                    for (size_t m = 0; m < M; ++m) {
                        const auto* p_br = tree.bridges.empty() ? nullptr : &P_bridges_mix[m];
                        for (int32_t node_id : tree.post_order) {
                            const auto& node = tree.nodes[node_id];
                            if (node.is_leaf) {
                                size_t t_idx = leaf_to_taxon[node_id];
                                if (t_idx != static_cast<size_t>(-1)) {
                                    int8_t state = pattern.states[t_idx];
                                    if (state >= 0 && state < S) {
                                        node_L[node_id].setZero();
                                        node_L[node_id](state) = 1.0;
                                    } else {
                                        node_L[node_id].setOnes();
                                    }
                                } else {
                                    node_L[node_id].setOnes();
                                }
                                node_log_scale[node_id] = 0.0;
                            } else {
                                double child_scale_sum = 0.0;
                                tree.eval_internal_node(node_id, P_mix[m], node_L, node_log_scale, node_L[node_id], child_scale_sum, p_br);
                                node_log_scale[node_id] = child_scale_sum;
                            }
                        }
                        Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                        if (cond_L > 0.0) {
                            Scalar term = std::log(syn_weights[m]) + std::log(cond_L) + node_log_scale[tree.root_id];
                            terms[m] = term;
                            valid[m] = true;
                            if (term > max_term) max_term = term;
                        }
                    }

                    if (max_term > -1e19) {
                        Scalar sum_exp = 0.0;
                        for (size_t m = 0; m < M; ++m) {
                            if (valid[m]) sum_exp += std::exp(terms[m] - max_term);
                        }
                        if (sum_exp > 0.0) {
                            pll = max_term + std::log(sum_exp);
                        }
                    }
                }
                pat_ll_buf[p] = pll;
                if (out_pattern_ll) {
                    (*out_pattern_ll)[p] = pll;
                }
            }
        }
        Scalar total_log_l = 0.0, c_ll = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            LikelihoodEngine::neumaier_acc(total_log_l, c_ll, aln.patterns[p].weight * pat_ll_buf[p]);
        }
        return total_log_l + c_ll;
    }

    // Two-set Expectation Step: simultaneously updates Test (partition 0), Background (partition 1), and optional Unclassified (partition 2) mixture weights in one Inside-Outside pass
    std::pair<std::vector<Scalar>, std::vector<Scalar>> compute_expected_two_set_weights(
        Scalar s,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0,
        Scalar* out_new_eps_test = nullptr,
        Scalar* out_new_eps_bg = nullptr,
        Scalar max_eps_sink = 0.01,
        const std::vector<Scalar>* omegas_unclass = nullptr,
        const std::vector<Scalar>* weights_unclass = nullptr,
        std::vector<Scalar>* out_new_weights_unclass = nullptr,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        bool has_unclass_dist = (num_unclassified_branches > 0 && omegas_unclass != nullptr && weights_unclass != nullptr && !omegas_unclass->empty() && omegas_unclass->size() == weights_unclass->size());
        size_t K_u = has_unclass_dist ? omegas_unclass->size() : 0;
        size_t M = syn_rates.size();
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        Scalar bio_scale_t = (out_new_eps_test != nullptr) ? (1.0 - eps_sink_test) : 1.0;
        Scalar bio_scale_b = (out_new_eps_bg != nullptr) ? (1.0 - eps_sink_bg) : 1.0;

        std::vector<std::vector<std::vector<Matrix>>> P_comp_t(M, std::vector<std::vector<Matrix>>(K_t, std::vector<Matrix>(num_nodes)));
        std::vector<std::vector<std::vector<Matrix>>> P_comp_b(M, std::vector<std::vector<Matrix>>(K_b, std::vector<Matrix>(num_nodes)));
        std::vector<std::vector<std::vector<Matrix>>> P_comp_u(has_unclass_dist ? M : 0, std::vector<std::vector<Matrix>>(K_u, std::vector<Matrix>(has_unclass_dist ? num_nodes : 0)));
        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<MG94Matrix> mg_t_vec(K_t), mg_b_vec(K_b), mg_u_vec(K_u);

        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas_test[k];
            p.delta = delta;
            p.psi = psi;
            mg_t_vec[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            for (size_t m = 0; m < M; ++m) {
                Scalar alpha_m = syn_rates[m];
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && branch_partition[node.id] == 0) {
                        Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[node.id] : branch_tau[node.id]);
                        P_comp_t[m][k][node.id] = mg_t_vec[k].transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0;
            p.beta = omegas_bg[k];
            p.delta = delta;
            p.psi = psi;
            mg_b_vec[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            for (size_t m = 0; m < M; ++m) {
                Scalar alpha_m = syn_rates[m];
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id && branch_partition[node.id] == 1) {
                        Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[node.id] : branch_tau[node.id]);
                        P_comp_b[m][k][node.id] = mg_b_vec[k].transition_matrix(alpha_m * bl);
                    }
                }
            }
        }

        MG94Matrix mg_unclass;
        if (has_unclass_dist) {
            for (size_t k = 0; k < K_u; ++k) {
                MG94Parameters p_u = active_base;
                p_u.alpha = 1.0;
                p_u.beta = (*omegas_unclass)[k];
                p_u.delta = delta;
                p_u.psi = psi;
                mg_u_vec[k].update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
                for (size_t m = 0; m < M; ++m) {
                    Scalar alpha_m = syn_rates[m];
                    for (const auto& node : tree.nodes) {
                        if (node.id != tree.root_id && branch_partition[node.id] == 2) {
                            Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[node.id] : branch_tau[node.id]);
                            P_comp_u[m][k][node.id] = mg_u_vec[k].transition_matrix(alpha_m * bl);
                        }
                    }
                }
            }
        } else if (num_unclassified_branches > 0) {
            MG94Parameters p_u = active_base;
            p_u.alpha = 1.0;
            p_u.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega;
            p_u.delta = delta;
            p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    int part = branch_partition[nid];
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_t; ++k) {
                            P_mix[m][nid] += (bio_scale_t * weights_test[k]) * P_comp_t[m][k][nid];
                        }
                        if (out_new_eps_test != nullptr) {
                            P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                        }
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_b; ++k) {
                            P_mix[m][nid] += (bio_scale_b * weights_bg[k]) * P_comp_b[m][k][nid];
                        }
                        if (out_new_eps_bg != nullptr) {
                            P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                        }
                    } else if (has_unclass_dist) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_u; ++k) {
                            P_mix[m][nid] += (*weights_unclass)[k] * P_comp_u[m][k][nid];
                        }
                    } else {
                        Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[nid] : branch_tau[nid]);
                        P_mix[m][nid] = mg_unclass.transition_matrix(alpha_m * bl);
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_two_set_bridge_mix(
                    s, alpha_m, mg_t_vec, mg_b_vec, mg_u_vec, mg_unclass,
                    weights_test, weights_bg, weights_unclass, has_unclass_dist,
                    (out_new_eps_test != nullptr) ? eps_sink_test : 0.0,
                    (out_new_eps_bg != nullptr) ? eps_sink_bg : 0.0,
                    bio_scale_t, bio_scale_b, custom_branch_lengths
                );
            }
        }

        size_t K_u_safe = std::max<size_t>(1, K_u);
        std::vector<Vector> pat_accum_t(num_patterns, Vector::Zero(K_t));
        std::vector<Vector> pat_accum_b(num_patterns, Vector::Zero(K_b));
        std::vector<Vector> pat_accum_u(num_patterns, Vector::Zero(K_u_safe));
        std::vector<Scalar> pat_sink_t(num_patterns, 0.0), pat_sink_b(num_patterns, 0.0);
        std::vector<Scalar> pat_wt_t(num_patterns, 0.0), pat_wt_b(num_patterns, 0.0), pat_wt_u(num_patterns, 0.0);

        const int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns);
        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<LikelihoodEngine::InsideOutsideWorkspace> io_m(M);

            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar max_ll = -std::numeric_limits<Scalar>::infinity();

                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    LikelihoodEngine::compute_inside_outside_in_place(
                        tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, io_m[m],
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                    );
                    if (io_m[m].likelihood > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                        if (term > max_ll) max_ll = term;
                    }
                }

                if (max_ll > -1e19) {
                    Scalar sum_exp = 0.0;
                    std::vector<Scalar> post_m(M, 0.0);
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                            Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood - max_ll);
                            post_m[m] = e;
                            sum_exp += e;
                        }
                    }
                    if (sum_exp > 0.0) {
                        for (size_t m = 0; m < M; ++m) post_m[m] /= sum_exp;
                        Vector local_t = Vector::Zero(K_t);
                        Vector local_b = Vector::Zero(K_b);
                        Vector local_u = Vector::Zero(K_u_safe);
                        Scalar local_sink_t = 0.0, local_sink_b = 0.0;
                        Scalar local_wt_t = 0.0, local_wt_b = 0.0, local_wt_u = 0.0;

                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int nid = node.id;
                                int part = branch_partition[nid];
                                if (part == 0) {
                                    for (size_t k = 0; k < K_t; ++k) {
                                        Scalar ratio_k = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                ratio_k += post_m[m] * (io_m[m].V[nid].dot(P_comp_t[m][k][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_t(k) += pattern.weight * (bio_scale_t * weights_test[k]) * ratio_k;
                                    }
                                    if (out_new_eps_test != nullptr) {
                                        Scalar sink_ratio = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                sink_ratio += post_m[m] * (eval_sink_bilinear(io_m[m].V[nid], io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_sink_t += pattern.weight * eps_sink_test * sink_ratio;
                                    }
                                    local_wt_t += pattern.weight;
                                } else if (part == 1) {
                                    for (size_t k = 0; k < K_b; ++k) {
                                        Scalar ratio_k = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                ratio_k += post_m[m] * (io_m[m].V[nid].dot(P_comp_b[m][k][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_b(k) += pattern.weight * (bio_scale_b * weights_bg[k]) * ratio_k;
                                    }
                                    if (out_new_eps_bg != nullptr) {
                                        Scalar sink_ratio = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                sink_ratio += post_m[m] * (eval_sink_bilinear(io_m[m].V[nid], io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_sink_b += pattern.weight * eps_sink_bg * sink_ratio;
                                    }
                                    local_wt_b += pattern.weight;
                                } else if (has_unclass_dist && out_new_weights_unclass != nullptr) {
                                    for (size_t k = 0; k < K_u; ++k) {
                                        Scalar ratio_k = 0.0;
                                        for (size_t m = 0; m < M; ++m) {
                                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                                ratio_k += post_m[m] * (io_m[m].V[nid].dot(P_comp_u[m][k][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                            }
                                        }
                                        local_u(k) += pattern.weight * (*weights_unclass)[k] * ratio_k;
                                    }
                                    local_wt_u += pattern.weight;
                                }
                            }
                        }
                        pat_accum_t[p] = local_t;
                        pat_accum_b[p] = local_b;
                        pat_accum_u[p] = local_u;
                        pat_sink_t[p] = local_sink_t;
                        pat_sink_b[p] = local_sink_b;
                        pat_wt_t[p] = local_wt_t;
                        pat_wt_b[p] = local_wt_b;
                        pat_wt_u[p] = local_wt_u;
                    }
                }
            }
        }

        Vector accum_t = Vector::Zero(K_t);
        Vector accum_b = Vector::Zero(K_b);
        Vector accum_u = Vector::Zero(K_u_safe);
        Scalar sink_accum_t = 0.0, sink_accum_b = 0.0;
        Scalar total_wt_t = 0.0, total_wt_b = 0.0, total_wt_u = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            accum_t += pat_accum_t[p];
            accum_b += pat_accum_b[p];
            accum_u += pat_accum_u[p];
            sink_accum_t += pat_sink_t[p];
            sink_accum_b += pat_sink_b[p];
            total_wt_t += pat_wt_t[p];
            total_wt_b += pat_wt_b[p];
            total_wt_u += pat_wt_u[p];
        }

        std::vector<Scalar> new_wt_t = weights_test;
        if (total_wt_t > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t k = 0; k < K_t; ++k) {
                new_wt_t[k] = std::clamp(accum_t(k) / total_wt_t, 1e-6, 1.0 - 1e-6);
                sum_w += new_wt_t[k];
            }
            for (size_t k = 0; k < K_t; ++k) new_wt_t[k] /= sum_w;
            if (out_new_eps_test != nullptr) {
                *out_new_eps_test = std::clamp(sink_accum_t / total_wt_t, 1e-7, max_eps_sink);
            }
        }

        std::vector<Scalar> new_wt_b = weights_bg;
        if (total_wt_b > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t k = 0; k < K_b; ++k) {
                new_wt_b[k] = std::clamp(accum_b(k) / total_wt_b, 1e-6, 1.0 - 1e-6);
                sum_w += new_wt_b[k];
            }
            for (size_t k = 0; k < K_b; ++k) new_wt_b[k] /= sum_w;
            if (out_new_eps_bg != nullptr) {
                *out_new_eps_bg = std::clamp(sink_accum_b / total_wt_b, 1e-7, max_eps_sink);
            }
        }

        if (has_unclass_dist && out_new_weights_unclass != nullptr) {
            *out_new_weights_unclass = *weights_unclass;
            if (total_wt_u > 0.0) {
                Scalar sum_w = 0.0;
                for (size_t k = 0; k < K_u; ++k) {
                    (*out_new_weights_unclass)[k] = std::clamp(accum_u(k) / total_wt_u, 1e-6, 1.0 - 1e-6);
                    sum_w += (*out_new_weights_unclass)[k];
                }
                for (size_t k = 0; k < K_u; ++k) (*out_new_weights_unclass)[k] /= sum_w;
            }
        }

        return {new_wt_t, new_wt_b};
    }

    // Two-set synonymous rate E-step (with log-underflow scaling)
    std::vector<Scalar> compute_expected_two_set_syn_weights(
        Scalar s,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates,
        const std::vector<Scalar>& syn_weights,
        const std::vector<Scalar>* custom_branch_lengths = nullptr,
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0,
        const std::vector<Scalar>* omegas_unclass = nullptr,
        const std::vector<Scalar>* weights_unclass = nullptr,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;
        size_t M = syn_rates.size();
        if (M <= 1) return {1.0};

        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        bool has_unclass_dist = (num_unclassified_branches > 0 && omegas_unclass != nullptr && weights_unclass != nullptr && !omegas_unclass->empty() && omegas_unclass->size() == weights_unclass->size());
        size_t K_u = has_unclass_dist ? omegas_unclass->size() : 0;
        Scalar bio_scale_t = 1.0 - eps_sink_test;
        Scalar bio_scale_b = 1.0 - eps_sink_bg;

        std::vector<MG94Matrix> mg_t(K_t), mg_b(K_b), mg_u(K_u);
        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0; p.beta = omegas_test[k]; p.delta = delta; p.psi = psi;
            mg_t[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0; p.beta = omegas_bg[k]; p.delta = delta; p.psi = psi;
            mg_b[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        MG94Matrix mg_unclass;
        if (has_unclass_dist) {
            for (size_t k = 0; k < K_u; ++k) {
                MG94Parameters p_u = active_base;
                p_u.alpha = 1.0; p_u.beta = (*omegas_unclass)[k]; p_u.delta = delta; p_u.psi = psi;
                mg_u[k].update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
        } else if (num_unclassified_branches > 0) {
            MG94Parameters p_u = active_base;
            p_u.alpha = 1.0; p_u.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega; p_u.delta = delta; p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = s * (custom_branch_lengths ? (*custom_branch_lengths)[nid] : branch_tau[nid]);
                    int part = branch_partition[nid];
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_test > 0.0) ? bio_scale_t : 1.0;
                        for (size_t k = 0; k < K_t; ++k) {
                            P_mix[m][nid] += (use_bio * weights_test[k]) * mg_t[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_test > 0.0) P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_bg > 0.0) ? bio_scale_b : 1.0;
                        for (size_t k = 0; k < K_b; ++k) {
                            P_mix[m][nid] += (use_bio * weights_bg[k]) * mg_b[k].transition_matrix(alpha_m * bl);
                        }
                        if (eps_sink_bg > 0.0) P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                    } else if (has_unclass_dist) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_u; ++k) {
                            P_mix[m][nid] += (*weights_unclass)[k] * mg_u[k].transition_matrix(alpha_m * bl);
                        }
                    } else {
                        P_mix[m][nid] = mg_unclass.transition_matrix(alpha_m * bl);
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_two_set_bridge_mix(
                    s, alpha_m, mg_t, mg_b, mg_u, mg_unclass,
                    weights_test, weights_bg, weights_unclass, has_unclass_dist,
                    eps_sink_test, eps_sink_bg, bio_scale_t, bio_scale_b, custom_branch_lengths
                );
            }
        }

        std::vector<Vector> pat_syn_accum(num_patterns, Vector::Zero(M));
        std::vector<Scalar> pat_valid_wt(num_patterns, 0.0);
        const int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns);

        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<Vector> node_L(num_nodes, Vector::Zero(S));
            std::vector<double> node_log_scale(num_nodes, 0.0);

            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar max_term = -std::numeric_limits<Scalar>::infinity();
                std::vector<Scalar> terms(M, -std::numeric_limits<Scalar>::infinity());
                std::vector<bool> valid(M, false);

                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    const auto* p_br = tree.bridges.empty() ? nullptr : &P_bridges_mix[m];
                    for (int32_t node_id : tree.post_order) {
                        const auto& node = tree.nodes[node_id];
                        if (node.is_leaf) {
                            size_t t_idx = leaf_to_taxon[node_id];
                            if (t_idx != static_cast<size_t>(-1)) {
                                int8_t state = pattern.states[t_idx];
                                if (state >= 0 && state < S) {
                                    node_L[node_id].setZero();
                                    node_L[node_id](state) = 1.0;
                                } else {
                                    node_L[node_id].setOnes();
                                }
                            } else {
                                node_L[node_id].setOnes();
                            }
                            node_log_scale[node_id] = 0.0;
                        } else {
                            double child_scale_sum = 0.0;
                            tree.eval_internal_node(node_id, P_mix[m], node_L, node_log_scale, node_L[node_id], child_scale_sum, p_br);
                            node_log_scale[node_id] = child_scale_sum;
                        }
                    }
                    Scalar cond_L = aln.codon_frequencies_f3x4.dot(node_L[tree.root_id]);
                    if (cond_L > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + std::log(cond_L) + node_log_scale[tree.root_id];
                        terms[m] = term;
                        valid[m] = true;
                        if (term > max_term) max_term = term;
                    }
                }

                if (max_term > -1e19) {
                    Scalar sum_exp = 0.0;
                    std::vector<Scalar> rho(M, 0.0);
                    for (size_t m = 0; m < M; ++m) {
                        if (valid[m]) {
                            rho[m] = std::exp(terms[m] - max_term);
                            sum_exp += rho[m];
                        }
                    }
                    if (sum_exp > 0.0) {
                        for (size_t m = 0; m < M; ++m) {
                            pat_syn_accum[p](m) = pattern.weight * (rho[m] / sum_exp);
                        }
                        pat_valid_wt[p] = pattern.weight;
                    }
                }
            }
        }

        Vector syn_accum = Vector::Zero(M);
        Scalar total_pattern_weight = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            syn_accum += pat_syn_accum[p];
            total_pattern_weight += pat_valid_wt[p];
        }

        std::vector<Scalar> new_syn_weights(M);
        if (total_pattern_weight > 0.0) {
            Scalar sum_w = 0.0;
            for (size_t m = 0; m < M; ++m) {
                new_syn_weights[m] = std::clamp(syn_accum(m) / total_pattern_weight, 1e-6, 1.0 - 1e-6);
                sum_w += new_syn_weights[m];
            }
            for (size_t m = 0; m < M; ++m) new_syn_weights[m] /= sum_w;
        } else {
            new_syn_weights = syn_weights;
        }
        return new_syn_weights;
    }

    // Two-set analytical branch gradients
    std::pair<Scalar, std::vector<Scalar>> compute_two_set_branch_gradients(
        const std::vector<Scalar>& branch_lengths,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0,
        const std::vector<Scalar>* omegas_unclass = nullptr,
        const std::vector<Scalar>* weights_unclass = nullptr,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        const MG94Parameters& active_base = custom_base_params ? *custom_base_params : base_params;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K_t = omegas_test.size();
        size_t K_b = omegas_bg.size();
        bool has_unclass_dist = (num_unclassified_branches > 0 && omegas_unclass != nullptr && weights_unclass != nullptr && !omegas_unclass->empty() && omegas_unclass->size() == weights_unclass->size());
        size_t K_u = has_unclass_dist ? omegas_unclass->size() : 0;
        size_t M = syn_rates.size();
        Scalar bio_scale_t = 1.0 - eps_sink_test;
        Scalar bio_scale_b = 1.0 - eps_sink_bg;

        std::vector<MG94Matrix> mg_t(K_t), mg_b(K_b), mg_u(K_u);
        for (size_t k = 0; k < K_t; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0; p.beta = omegas_test[k]; p.delta = delta; p.psi = psi;
            mg_t[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        for (size_t k = 0; k < K_b; ++k) {
            MG94Parameters p = active_base;
            p.alpha = 1.0; p.beta = omegas_bg[k]; p.delta = delta; p.psi = psi;
            mg_b[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }
        MG94Matrix mg_unclass;
        if (has_unclass_dist) {
            for (size_t k = 0; k < K_u; ++k) {
                MG94Parameters p_u = active_base;
                p_u.alpha = 1.0; p_u.beta = (*omegas_unclass)[k]; p_u.delta = delta; p_u.psi = psi;
                mg_u[k].update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
        } else if (num_unclassified_branches > 0) {
            MG94Parameters p_u = active_base;
            p_u.alpha = 1.0; p_u.beta = (mg94_background_omega >= 0.0) ? mg94_background_omega : mg94_omega; p_u.delta = delta; p_u.psi = psi;
            mg_unclass.update(p_u, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> dP_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M);
        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    int nid = node.id;
                    Scalar bl = branch_lengths[nid];
                    int part = branch_partition[nid];
                    if (part == 0) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        dP_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_test > 0.0) ? bio_scale_t : 1.0;
                        for (size_t k = 0; k < K_t; ++k) {
                            Matrix P_k = mg_t[k].transition_matrix(alpha_m * bl);
                            P_mix[m][nid] += (use_bio * weights_test[k]) * P_k;
                            dP_mix[m][nid] += (use_bio * weights_test[k]) * alpha_m * (mg_t[k].Q * P_k);
                        }
                        if (eps_sink_test > 0.0) P_mix[m][nid].noalias() += eps_sink_test * P_sink;
                    } else if (part == 1) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        dP_mix[m][nid] = Matrix::Zero(S, S);
                        Scalar use_bio = (eps_sink_bg > 0.0) ? bio_scale_b : 1.0;
                        for (size_t k = 0; k < K_b; ++k) {
                            Matrix P_k = mg_b[k].transition_matrix(alpha_m * bl);
                            P_mix[m][nid] += (use_bio * weights_bg[k]) * P_k;
                            dP_mix[m][nid] += (use_bio * weights_bg[k]) * alpha_m * (mg_b[k].Q * P_k);
                        }
                        if (eps_sink_bg > 0.0) P_mix[m][nid].noalias() += eps_sink_bg * P_sink;
                    } else if (has_unclass_dist) {
                        P_mix[m][nid] = Matrix::Zero(S, S);
                        dP_mix[m][nid] = Matrix::Zero(S, S);
                        for (size_t k = 0; k < K_u; ++k) {
                            Matrix P_k = mg_u[k].transition_matrix(alpha_m * bl);
                            P_mix[m][nid] += (*weights_unclass)[k] * P_k;
                            dP_mix[m][nid] += (*weights_unclass)[k] * alpha_m * (mg_u[k].Q * P_k);
                        }
                    } else {
                        Matrix P_u = mg_unclass.transition_matrix(alpha_m * bl);
                        P_mix[m][nid] = P_u;
                        dP_mix[m][nid] = alpha_m * (mg_unclass.Q * P_u);
                    }
                }
            }
            if (!tree.bridges.empty()) {
                P_bridges_mix[m] = build_two_set_bridge_mix(
                    1.0, alpha_m, mg_t, mg_b, mg_u, mg_unclass,
                    weights_test, weights_bg, weights_unclass, has_unclass_dist,
                    eps_sink_test, eps_sink_bg, bio_scale_t, bio_scale_b, &branch_lengths
                );
            }
        }

        std::vector<Scalar> pat_ll_vec(num_patterns, 0.0);
        std::vector<std::vector<Scalar>> pat_grad(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        const int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns);

        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<LikelihoodEngine::InsideOutsideWorkspace> io_m(M);

            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                if (M == 1) {
                    LikelihoodEngine::compute_inside_outside_in_place(
                        tree, pattern, leaf_to_taxon, P_mix[0], aln.codon_frequencies_f3x4, S, io_m[0],
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[0]
                    );
                    if (io_m[0].likelihood > 0.0) {
                        pat_ll_vec[p] = pattern.weight * io_m[0].log_likelihood;
                        Scalar inv_L = pattern.weight / io_m[0].likelihood;
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id) {
                                int32_t vid = node.id;
                                Scalar dL_dt = io_m[0].V[vid].dot(dP_mix[0][vid] * io_m[0].D[vid]);
                                pat_grad[p][vid] = inv_L * dL_dt;
                            }
                        }
                    }
                } else {
                    Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
                    for (size_t m = 0; m < M; ++m) {
                        if (syn_weights[m] <= 1e-12) continue;
                        LikelihoodEngine::compute_inside_outside_in_place(
                            tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, io_m[m],
                            tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                        );
                        if (io_m[m].likelihood > 0.0) {
                            Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                            if (term > max_ll) max_ll = term;
                        }
                    }

                    if (max_ll > -1e19) {
                        Scalar sum_exp = 0.0;
                        std::vector<Scalar> post_w(M, 0.0);
                        for (size_t m = 0; m < M; ++m) {
                            if (syn_weights[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood - max_ll);
                                post_w[m] = e;
                                sum_exp += e;
                            }
                        }
                        if (sum_exp > 0.0) {
                            Scalar pat_ll = max_ll + std::log(sum_exp);
                            pat_ll_vec[p] = pattern.weight * pat_ll;
                            Scalar inv_sum = 1.0 / sum_exp;
                            for (size_t m = 0; m < M; ++m) post_w[m] *= inv_sum;

                            for (const auto& node : tree.nodes) {
                                if (node.id != tree.root_id) {
                                    int32_t vid = node.id;
                                    Scalar d_lnL_dt = 0.0;
                                    for (size_t m = 0; m < M; ++m) {
                                        if (post_w[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                            d_lnL_dt += post_w[m] * (io_m[m].V[vid].dot(dP_mix[m][vid] * io_m[m].D[vid]) / io_m[m].likelihood);
                                        }
                                    }
                                    pat_grad[p][vid] = pattern.weight * d_lnL_dt;
                                }
                            }
                        }
                    }
                }
            }
        }

        Scalar total_log_l = LikelihoodEngine::neumaier_sum(pat_ll_vec);
        std::vector<Scalar> grad_b(num_nodes, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            Scalar s_acc = 0.0, c_acc = 0.0;
            for (size_t p = 0; p < num_patterns; ++p) {
                LikelihoodEngine::neumaier_acc(s_acc, c_acc, pat_grad[p][i]);
            }
            grad_b[i] = s_acc + c_acc;
        }

        return {total_log_l, grad_b};
    }

    std::vector<Scalar> refine_two_set_branch_lengths(
        const std::vector<Scalar>& init_bl,
        const std::vector<Scalar>& omegas_test,
        const std::vector<Scalar>& weights_test,
        const std::vector<Scalar>& omegas_bg,
        const std::vector<Scalar>& weights_bg,
        const std::vector<Scalar>& syn_rates = {1.0},
        const std::vector<Scalar>& syn_weights = {1.0},
        Scalar delta = 0.0,
        Scalar psi = 0.0,
        int max_iters = 25,
        Scalar eps_sink_test = 0.0,
        Scalar eps_sink_bg = 0.0,
        const std::vector<Scalar>* omegas_unclass = nullptr,
        const std::vector<Scalar>* weights_unclass = nullptr,
        const MG94Parameters* custom_base_params = nullptr
    ) const {
        size_t num_nodes = tree.num_nodes();
        std::vector<Scalar> bl = init_bl;
        std::vector<Scalar> log_bl(num_nodes, 0.0);
        for (size_t i = 0; i < num_nodes; ++i) {
            if (i != static_cast<size_t>(tree.root_id)) {
                bl[i] = std::clamp(bl[i], 1e-8, 10000.0);
                log_bl[i] = std::log(bl[i]);
            }
        }

        const size_t m_history = 8;
        std::vector<std::vector<Scalar>> s_hist, y_hist;
        std::vector<Scalar> rho_hist;

        auto eval_log_space = [&](const std::vector<Scalar>& y_vec) -> std::tuple<Scalar, std::vector<Scalar>, std::vector<Scalar>> {
            std::vector<Scalar> cur_bl(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    cur_bl[i] = std::exp(std::clamp(y_vec[i], -18.420680744, 9.210340372));
                }
            }
            auto [ll, grad_t] = compute_two_set_branch_gradients(
                cur_bl, omegas_test, weights_test, omegas_bg, weights_bg, syn_rates, syn_weights, delta, psi, eps_sink_test, eps_sink_bg, omegas_unclass, weights_unclass, custom_base_params
            );
            std::vector<Scalar> grad_y(num_nodes, 0.0);
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    grad_y[i] = cur_bl[i] * grad_t[i];
                }
            }
            return {ll, grad_y, cur_bl};
        };

        auto [curr_ll, curr_grad, curr_bl] = eval_log_space(log_bl);
        bl = curr_bl;
        Scalar best_ll = curr_ll;
        std::vector<Scalar> best_bl = bl;

        for (int iter = 0; iter < max_iters; ++iter) {
            Scalar max_g = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    max_g = std::max(max_g, std::abs(curr_grad[i]));
                }
            }
            if (max_g < 1e-4) break;

            std::vector<Scalar> q = curr_grad;
            std::vector<Scalar> alpha(s_hist.size());
            for (int j = static_cast<int>(s_hist.size()) - 1; j >= 0; --j) {
                Scalar s_dot_q = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) s_dot_q += s_hist[j][i] * q[i];
                alpha[j] = rho_hist[j] * s_dot_q;
                for (size_t i = 0; i < num_nodes; ++i) q[i] -= alpha[j] * y_hist[j][i];
            }

            Scalar gamma = 0.25;
            if (!s_hist.empty()) {
                Scalar s_dot_y = 0.0, y_dot_y = 0.0;
                const auto& s_last = s_hist.back();
                const auto& y_last = y_hist.back();
                for (size_t i = 0; i < num_nodes; ++i) {
                    s_dot_y += s_last[i] * y_last[i];
                    y_dot_y += y_last[i] * y_last[i];
                }
                if (y_dot_y > 1e-12) gamma = std::clamp(s_dot_y / y_dot_y, 1e-3, 5.0);
            }

            std::vector<Scalar> r(num_nodes);
            for (size_t i = 0; i < num_nodes; ++i) r[i] = gamma * q[i];
            for (size_t j = 0; j < s_hist.size(); ++j) {
                Scalar y_dot_r = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) y_dot_r += y_hist[j][i] * r[i];
                Scalar beta = rho_hist[j] * y_dot_r;
                for (size_t i = 0; i < num_nodes; ++i) r[i] += s_hist[j][i] * (alpha[j] - beta);
            }

            std::vector<Scalar> d = r;
            d[tree.root_id] = 0.0;
            Scalar dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) dir_dot_grad += d[i] * curr_grad[i];
            }
            if (dir_dot_grad <= 0.0) {
                d = curr_grad;
                d[tree.root_id] = 0.0;
                dir_dot_grad = 0.0;
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (i != static_cast<size_t>(tree.root_id)) dir_dot_grad += d[i] * curr_grad[i];
                }
                s_hist.clear(); y_hist.clear(); rho_hist.clear();
            }

            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) {
                    d[i] = std::clamp(d[i], -2.0, 2.0);
                }
            }
            dir_dot_grad = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (i != static_cast<size_t>(tree.root_id)) dir_dot_grad += d[i] * curr_grad[i];
            }

            Scalar step_size = 1.0;
            std::vector<Scalar> new_log_bl(num_nodes, 0.0);
            std::vector<Scalar> new_bl(num_nodes, 0.0);
            Scalar new_ll = curr_ll;
            std::vector<Scalar> new_grad;
            bool line_search_ok = false;
            for (int ls = 0; ls < 12; ++ls) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    new_log_bl[i] = (i != static_cast<size_t>(tree.root_id))
                        ? std::clamp(log_bl[i] + step_size * d[i], -18.420680744, 9.210340372)
                        : 0.0;
                }
                auto [test_ll, test_grad, test_bl] = eval_log_space(new_log_bl);
                if (test_ll > curr_ll + 1e-4 * step_size * dir_dot_grad || (test_ll > curr_ll + 1e-5 && ls >= 3)) {
                    new_ll = test_ll;
                    new_grad = test_grad;
                    new_bl = test_bl;
                    line_search_ok = true;
                    break;
                }
                step_size *= 0.5;
            }
            if (!line_search_ok) break;

            std::vector<Scalar> s_k(num_nodes), y_k(num_nodes);
            Scalar s_dot_y = 0.0;
            for (size_t i = 0; i < num_nodes; ++i) {
                s_k[i] = new_log_bl[i] - log_bl[i];
                y_k[i] = curr_grad[i] - new_grad[i];
                s_dot_y += s_k[i] * y_k[i];
            }
            if (s_dot_y > 1e-8) {
                if (s_hist.size() >= m_history) {
                    s_hist.erase(s_hist.begin()); y_hist.erase(y_hist.begin()); rho_hist.erase(rho_hist.begin());
                }
                s_hist.push_back(s_k); y_hist.push_back(y_k); rho_hist.push_back(1.0 / s_dot_y);
            }

            log_bl = new_log_bl;
            bl = new_bl; curr_ll = new_ll; curr_grad = new_grad;
            if (curr_ll > best_ll) {
                best_ll = curr_ll; best_bl = bl;
            }
        }
        return best_bl;
    }

    // Fits two-set (or three-set when num_unclassified_branches > 0) BUSTED / BUSTED-PH model
    // with independent K-class omega distributions on Test, Background, and (if present) Unclassified branches.
    BUSTEDFit fit_two_set_model(
        bool constrain_test,
        bool constrain_bg,
        size_t K = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        const std::string& model_desc = "",
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01,
        const BUSTEDRateDistribution* init_test_dist = nullptr,
        const BUSTEDRateDistribution* init_bg_dist = nullptr,
        const BUSTEDRateDistribution* init_unclass_dist = nullptr,
        const std::vector<Scalar>* init_branch_lengths = nullptr,
        const MG94Parameters* init_base_params = nullptr,
        Scalar init_delta = 0.0,
        Scalar init_psi = 0.0
    ) const {
        K = std::max<size_t>(1, K);
        BUSTEDFit fit;
        fit.num_rate_classes = K;
        fit.has_background_distribution = true;
        bool has_u = (num_unclassified_branches > 0);
        fit.has_unclassified_distribution = has_u;
        Scalar s = 1.0;

        MG94Parameters local_base = init_base_params ? *init_base_params : base_params;

        bool allow_double = (multiple_hits == "Double" || multiple_hits == "Double+Triple");
        bool allow_triple = (multiple_hits == "Double+Triple");
        Scalar delta = allow_double ? std::max(0.0, init_delta) : 0.0;
        Scalar psi = allow_triple ? std::max(0.0, init_psi) : 0.0;
        Scalar eps_sink_t = use_error_sink
            ? (init_test_dist && init_test_dist->error_sink_weight > 0.0
                ? std::min(init_test_dist->error_sink_weight, max_eps_sink)
                : std::min(0.005, max_eps_sink * 0.5))
            : 0.0;
        Scalar eps_sink_b = use_error_sink
            ? (init_bg_dist && init_bg_dist->error_sink_weight > 0.0
                ? std::min(init_bg_dist->error_sink_weight, max_eps_sink)
                : std::min(0.005, max_eps_sink * 0.5))
            : 0.0;

        auto sort_pairs = [](std::vector<Scalar>& rates, std::vector<Scalar>& wts, size_t n_sort) {
            if (n_sort <= 1) return;
            std::vector<size_t> idx(n_sort);
            for (size_t i = 0; i < n_sort; ++i) idx[i] = i;
            std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return rates[a] < rates[b]; });
            std::vector<Scalar> r_tmp(n_sort), w_tmp(n_sort);
            for (size_t i = 0; i < n_sort; ++i) {
                r_tmp[i] = rates[idx[i]];
                w_tmp[i] = wts[idx[i]];
            }
            for (size_t i = 0; i < n_sort; ++i) {
                rates[i] = r_tmp[i];
                wts[i] = w_tmp[i];
            }
        };

        auto init_dist = [&](bool is_con, const BUSTEDRateDistribution* warm, std::vector<Scalar>& omegas, std::vector<Scalar>& weights) {
            if (warm && warm->omegas.size() == K && warm->weights.size() == K) {
                omegas = warm->omegas;
                weights = warm->weights;
                if (is_con) {
                    omegas.back() = (K == 1) ? std::min(1.0, omegas.back()) : 1.0;
                    for (size_t k = 0; k + 1 < K; ++k) omegas[k] = std::clamp(omegas[k], 0.0, 1.0);
                } else {
                    for (size_t k = 0; k + 1 < K; ++k) omegas[k] = std::clamp(omegas[k], 0.0, 1.0);
                    if (omegas.back() <= 1.001) {
                        omegas.back() = 2.0;
                    }
                }
                return;
            }
            omegas.resize(K);
            weights.resize(K);
            if (K == 1) {
                weights[0] = 1.0;
                omegas[0] = is_con ? std::min(1.0, mg94_omega) : mg94_omega;
            } else if (K == 2) {
                weights = {0.75, 0.25};
                omegas = is_con ? std::vector<Scalar>{0.10, 1.0} : std::vector<Scalar>{0.10, 3.0};
            } else if (K == 3) {
                if (is_con) {
                    omegas = {0.05, 0.50, 1.0};
                    weights = {0.45, 0.45, 0.10};
                } else {
                    omegas = {0.05, 0.50, 3.5};
                    weights = {0.50, 0.40, 0.10};
                }
            } else {
                Scalar p_top = 0.12;
                Scalar rem_w = 1.0 - p_top;
                for (size_t k = 0; k + 1 < K; ++k) {
                    Scalar frac = static_cast<Scalar>(k) / static_cast<Scalar>(K - 2);
                    omegas[k] = 0.01 * std::pow(80.0, frac);
                    weights[k] = rem_w / static_cast<Scalar>(K - 1);
                }
                omegas[K - 1] = is_con ? 1.0 : 3.5;
                weights[K - 1] = p_top;
            }
        };

        std::vector<Scalar> omegas_t, weights_t, omegas_b, weights_b, omegas_u, weights_u;
        init_dist(constrain_test, init_test_dist, omegas_t, weights_t);
        init_dist(constrain_bg, init_bg_dist, omegas_b, weights_b);
        if (has_u) {
            init_dist(false, init_unclass_dist, omegas_u, weights_u);
        }
        const std::vector<Scalar>* p_omegas_u = has_u ? &omegas_u : nullptr;
        const std::vector<Scalar>* p_weights_u = has_u ? &weights_u : nullptr;
        const std::vector<Scalar>* ecm_bl_ptr = (init_branch_lengths && !init_branch_lengths->empty()) ? init_branch_lengths : nullptr;

        size_t M = srv ? std::max<size_t>(1, num_syn_rates) : 1;
        std::vector<Scalar> syn_rates(M, 1.0);
        std::vector<Scalar> syn_weights(M, 1.0 / static_cast<Scalar>(M));
        if (srv && M == 2) {
            syn_rates = {0.35, 1.65};
            syn_weights = {0.50, 0.50};
        } else if (srv && M == 3) {
            syn_rates = {0.30, 1.00, 3.10};
            syn_weights = {0.30, 0.60, 0.10};
        } else if (srv && M >= 4) {
            Scalar mean_r = 0.0;
            for (size_t m = 0; m < M; ++m) {
                Scalar frac = static_cast<Scalar>(m) / static_cast<Scalar>(M - 1);
                syn_rates[m] = 0.15 * std::pow(25.0, frac);
                mean_r += syn_weights[m] * syn_rates[m];
            }
            for (size_t m = 0; m < M; ++m) syn_rates[m] /= std::max(1e-12, mean_r);
        }
        if (srv && init_test_dist && init_test_dist->syn_rates.size() == M && init_test_dist->syn_weights.size() == M) {
            syn_rates = init_test_dist->syn_rates;
            syn_weights = init_test_dist->syn_weights;
        }

        // Multi-candidate initial grid search for K >= 2 (F-008, F-059)
        if (K >= 2) {
            auto build_k_cands = [&](bool is_con, const BUSTEDRateDistribution* warm) {
                std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> cands;
                if (K == 2) {
                    if (is_con) {
                        cands = {
                            {{0.10, 1.0}, {0.75, 0.25}},
                            {{0.02, 1.0}, {0.85, 0.15}},
                            {{0.25, 1.0}, {0.60, 0.40}},
                            {{0.0, 1.0},  {0.70, 0.30}}
                        };
                    } else {
                        cands = {
                            {{0.10, 3.0}, {0.75, 0.25}},
                            {{0.02, 4.5}, {0.90, 0.10}},
                            {{0.05, 2.0}, {0.80, 0.20}},
                            {{0.02, 3.0}, {0.95, 0.05}}
                        };
                    }
                } else if (K == 3) {
                    if (is_con) {
                        cands = {
                            {{0.05, 0.50, 1.0}, {0.45, 0.45, 0.10}},
                            {{0.02, 0.30, 1.0}, {0.50, 0.35, 0.15}},
                            {{0.10, 0.70, 1.0}, {0.40, 0.40, 0.20}},
                            {{0.0, 1.0, 1.0},   {0.40, 0.30, 0.30}}
                        };
                    } else {
                        cands = {
                            {{0.05, 0.40, 2.0}, {0.55, 0.35, 0.10}},
                            {{0.02, 0.25, 8.0}, {0.60, 0.35, 0.05}},
                            {{0.10, 0.80, 4.0}, {0.45, 0.45, 0.10}},
                            {{0.005, 0.10, 1.8}, {0.45, 0.15, 0.40}}
                        };
                    }
                } else {
                    for (Scalar top_o : (is_con ? std::vector<Scalar>{1.0} : std::vector<Scalar>{2.0, 4.0})) {
                        for (Scalar p_top : {0.10, 0.25}) {
                            std::vector<Scalar> co(K), cw(K);
                            Scalar rem_w = 1.0 - p_top;
                            for (size_t k = 0; k + 1 < K; ++k) {
                                Scalar frac = (K > 2) ? (static_cast<Scalar>(k) / static_cast<Scalar>(K - 2)) : 0.5;
                                co[k] = 0.01 * std::pow(80.0, frac);
                                cw[k] = rem_w / static_cast<Scalar>(K - 1);
                            }
                            co[K - 1] = top_o;
                            cw[K - 1] = p_top;
                            cands.push_back({co, cw});
                        }
                    }
                }
                if (warm && warm->omegas.size() == K && warm->weights.size() == K) {
                    auto co = warm->omegas;
                    auto cw = warm->weights;
                    for (size_t k = 0; k + 1 < K; ++k) co[k] = std::clamp(co[k], 0.0, 1.0);
                    co[K - 1] = is_con ? 1.0 : std::max(1.0, co[K - 1]);
                    cands.push_back({co, cw});
                }
                return cands;
            };

            Scalar best_init_ll = evaluate_two_set_log_l(
                s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights,
                ecm_bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
            );
            auto cands_t = build_k_cands(constrain_test, init_test_dist);
            auto cands_b = build_k_cands(constrain_bg, init_bg_dist);
            std::vector<Scalar> eps_cands = use_error_sink
                ? std::vector<Scalar>{eps_sink_t, std::min(0.001, max_eps_sink), std::min(0.005, max_eps_sink), std::min(0.01, max_eps_sink)}
                : std::vector<Scalar>{0.0};

            for (Scalar cand_e : eps_cands) {
                for (const auto& [ot_c, wt_c] : cands_t) {
                    Scalar ll_c = evaluate_two_set_log_l(
                        s, ot_c, wt_c, omegas_b, weights_b, syn_rates, syn_weights,
                        ecm_bl_ptr, nullptr, delta, psi, cand_e, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                    );
                    if (ll_c > best_init_ll) {
                        best_init_ll = ll_c;
                        omegas_t = ot_c;
                        weights_t = wt_c;
                        eps_sink_t = cand_e;
                    }
                }
                for (const auto& [ob_c, wb_c] : cands_b) {
                    Scalar ll_c = evaluate_two_set_log_l(
                        s, omegas_t, weights_t, ob_c, wb_c, syn_rates, syn_weights,
                        ecm_bl_ptr, nullptr, delta, psi, eps_sink_t, cand_e, p_omegas_u, p_weights_u, &local_base
                    );
                    if (ll_c > best_init_ll) {
                        best_init_ll = ll_c;
                        omegas_b = ob_c;
                        weights_b = wb_c;
                        eps_sink_b = cand_e;
                    }
                }
            }
            if ((constrain_test || constrain_bg) && ecm_bl_ptr != nullptr) {
                std::vector<Scalar> fresh_syn_r = (srv && M == 3) ? std::vector<Scalar>{0.30, 1.00, 3.10} : syn_rates;
                std::vector<Scalar> fresh_syn_w = (srv && M == 3) ? std::vector<Scalar>{0.30, 0.60, 0.10} : syn_weights;
                for (const auto& [ot_c, wt_c] : cands_t) {
                    Scalar ll_fresh = evaluate_two_set_log_l(
                        1.0, ot_c, wt_c, omegas_b, weights_b, fresh_syn_r, fresh_syn_w,
                        nullptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                    );
                    if (ll_fresh > best_init_ll) {
                        best_init_ll = ll_fresh;
                        omegas_t = ot_c;
                        weights_t = wt_c;
                        ecm_bl_ptr = nullptr;
                        s = 1.0;
                        syn_rates = fresh_syn_r;
                        syn_weights = fresh_syn_w;
                    }
                }
            }
        }

        auto renormalize_two_set_srv = [&](Scalar& scale_s, std::vector<Scalar>* bl_vec) {
            if (!srv || M < 2) return;
            sort_pairs(syn_rates, syn_weights, M);
            Scalar mean_r = 0.0;
            for (size_t m = 0; m < M; ++m) mean_r += syn_weights[m] * syn_rates[m];
            if (mean_r > 1e-8 && std::abs(mean_r - 1.0) > 1e-12) {
                for (size_t m = 0; m < M; ++m) syn_rates[m] /= mean_r;
                if (bl_vec) {
                    for (size_t i = 0; i < bl_vec->size(); ++i) {
                        if (i != static_cast<size_t>(tree.root_id)) {
                            (*bl_vec)[i] = std::clamp((*bl_vec)[i] * mean_r, 1e-8, 10000.0);
                        }
                    }
                } else {
                    scale_s = std::clamp(scale_s * mean_r, 1e-4, 100.0);
                }
            }
        };

        auto optimize_two_set_srv_rates = [&](Scalar& cur_s, std::vector<Scalar>* bl_vec, int max_it) {
            if (!srv || M < 2) return;
            const std::vector<Scalar>* b_ptr = bl_vec ? bl_vec : ecm_bl_ptr;
            Scalar eval_s = bl_vec ? 1.0 : cur_s;
            if (M == 2) {
                auto [r_opt, nll_opt] = Brent1D::minimize([&](Scalar ratio) {
                    Scalar denom = std::max(1e-12, syn_weights[0] * ratio + syn_weights[1]);
                    return -evaluate_two_set_log_l(eval_s, omegas_t, weights_t, omegas_b, weights_b, {ratio / denom, 1.0 / denom}, syn_weights, b_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, 1e-3, std::clamp(syn_rates[0] / std::max(1e-6, syn_rates[1]), 1e-3, 100.0), 100.0, 1e-4, max_it);
                Scalar denom = std::max(1e-12, syn_weights[0] * r_opt + syn_weights[1]);
                std::vector<Scalar> cand_r = {r_opt / denom, 1.0 / denom};
                if (-nll_opt >= evaluate_two_set_log_l(eval_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, b_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base)) {
                    syn_rates = cand_r;
                }
            } else {
                size_t D_syn = M - 1;
                std::vector<Scalar> x0(D_syn), lb(D_syn, 1e-3), ub(D_syn, 100.0);
                Scalar ref_r = std::max(1e-6, syn_rates[M - 1]);
                for (size_t m = 0; m < D_syn; ++m) x0[m] = std::clamp(syn_rates[m] / ref_r, 1e-3, 100.0);
                auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                    std::vector<Scalar> cand_r(M, 1.0);
                    Scalar mean_c = syn_weights[M - 1];
                    for (size_t m = 0; m < D_syn; ++m) {
                        cand_r[m] = std::max(1e-4, x[m]);
                        mean_c += syn_weights[m] * cand_r[m];
                    }
                    mean_c = std::max(1e-12, mean_c);
                    for (size_t m = 0; m < M; ++m) cand_r[m] /= mean_c;
                    return -evaluate_two_set_log_l(eval_s, omegas_t, weights_t, omegas_b, weights_b, cand_r, syn_weights, b_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, x0, lb, ub, 1e-4, max_it);
                std::vector<Scalar> cand_r(M, 1.0);
                Scalar mean_c = syn_weights[M - 1];
                for (size_t m = 0; m < D_syn; ++m) {
                    cand_r[m] = std::max(1e-4, opt.x[m]);
                    mean_c += syn_weights[m] * cand_r[m];
                }
                mean_c = std::max(1e-12, mean_c);
                for (size_t m = 0; m < M; ++m) cand_r[m] /= mean_c;
                if (-opt.f >= evaluate_two_set_log_l(eval_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, b_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base)) {
                    syn_rates = cand_r;
                }
            }
            renormalize_two_set_srv(cur_s, bl_vec);
        };

        if (allow_double && init_delta <= 0.0 && init_psi <= 0.0) {
            Scalar best_grid_ll = -1e20;
            std::vector<Scalar> d_cands = {0.0, 0.02, 0.08, 0.20, 0.60};
            std::vector<Scalar> p_cands = allow_triple ? std::vector<Scalar>{0.0, 0.02, 0.08} : std::vector<Scalar>{0.0};
            for (Scalar d_c : d_cands) {
                for (Scalar p_c : p_cands) {
                    Scalar cand_ll = evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, ecm_bl_ptr, nullptr, d_c, p_c, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                    if (cand_ll > best_grid_ll) {
                        best_grid_ll = cand_ll;
                        delta = d_c;
                        psi = p_c;
                    }
                }
            }
        }

        Scalar max_omega = use_error_sink ? 100.0 : 10000.0;

        auto optimize_partition_omegas = [&](
            bool is_con,
            std::vector<Scalar>& target_o,
            std::vector<Scalar>& target_w,
            Scalar& cur_s,
            const std::vector<Scalar>* bl_ptr,
            bool opt_scale,
            int max_it
        ) {
            if (is_con) {
                if (K == 1) {
                    if (opt_scale) {
                        auto opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                            std::vector<Scalar> co = {x[1]};
                            const auto& ot = (&target_o == &omegas_t) ? co : omegas_t;
                            const auto& ob = (&target_o == &omegas_b) ? co : omegas_b;
                            const auto* ou = (&target_o == &omegas_u) ? &co : p_omegas_u;
                            return -evaluate_two_set_log_l(x[0], ot, weights_t, ob, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, ou, p_weights_u, &local_base);
                        }, {cur_s, std::clamp(target_o[0], 0.0, 1.0)}, {0.01, 0.0}, {50.0, 1.0}, 1e-4, max_it);
                        cur_s = opt.x[0]; target_o[0] = opt.x[1];
                    } else {
                        auto [w_opt, nll_opt] = Brent1D::minimize([&](Scalar w0) {
                            std::vector<Scalar> co = {w0};
                            const auto& ot = (&target_o == &omegas_t) ? co : omegas_t;
                            const auto& ob = (&target_o == &omegas_b) ? co : omegas_b;
                            const auto* ou = (&target_o == &omegas_u) ? &co : p_omegas_u;
                            return -evaluate_two_set_log_l(1.0, ot, weights_t, ob, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, ou, p_weights_u, &local_base);
                        }, 0.0, std::clamp(target_o[0], 0.0, 1.0), 1.0, 1e-4, max_it);
                        target_o[0] = w_opt;
                    }
                } else {
                    size_t D = (K - 1) + (opt_scale ? 1 : 0);
                    std::vector<Scalar> x0(D), lb(D, 0.0), ub(D, 1.0);
                    size_t off = 0;
                    if (opt_scale) {
                        x0[0] = cur_s; lb[0] = 0.01; ub[0] = 50.0; off = 1;
                    }
                    for (size_t k = 0; k + 1 < K; ++k) x0[off + k] = std::clamp(target_o[k], 0.0, 1.0);
                    auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                        Scalar sc = opt_scale ? x[0] : 1.0;
                        std::vector<Scalar> co(K, 1.0);
                        for (size_t k = 0; k + 1 < K; ++k) co[k] = x[off + k];
                        const auto& ot = (&target_o == &omegas_t) ? co : omegas_t;
                        const auto& ob = (&target_o == &omegas_b) ? co : omegas_b;
                        const auto* ou = (&target_o == &omegas_u) ? &co : p_omegas_u;
                        return -evaluate_two_set_log_l(sc, ot, weights_t, ob, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, ou, p_weights_u, &local_base);
                    }, x0, lb, ub, 1e-4, max_it);
                    if (opt_scale) cur_s = opt.x[0];
                    for (size_t k = 0; k + 1 < K; ++k) target_o[k] = opt.x[off + k];
                    target_o[K - 1] = 1.0;
                    sort_pairs(target_o, target_w, K - 1);
                }
            } else {
                size_t D = K + (opt_scale ? 1 : 0);
                std::vector<Scalar> x0(D), lb(D, 0.0), ub(D, 1.0);
                size_t off = 0;
                if (opt_scale) {
                    x0[0] = cur_s; lb[0] = 0.01; ub[0] = 50.0; off = 1;
                }
                for (size_t k = 0; k + 1 < K; ++k) x0[off + k] = std::clamp(target_o[k], 0.0, 1.0);
                x0[off + K - 1] = std::clamp(target_o[K - 1], (K == 1 ? 0.0 : 1.0), max_omega);
                lb[off + K - 1] = (K == 1) ? 0.0 : 1.0;
                ub[off + K - 1] = max_omega;
                auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                    Scalar sc = opt_scale ? x[0] : 1.0;
                    std::vector<Scalar> co(K);
                    for (size_t k = 0; k < K; ++k) co[k] = x[off + k];
                    const auto& ot = (&target_o == &omegas_t) ? co : omegas_t;
                    const auto& ob = (&target_o == &omegas_b) ? co : omegas_b;
                    const auto* ou = (&target_o == &omegas_u) ? &co : p_omegas_u;
                    return -evaluate_two_set_log_l(sc, ot, weights_t, ob, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, ou, p_weights_u, &local_base);
                }, x0, lb, ub, 1e-4, max_it);
                if (opt_scale) cur_s = opt.x[0];
                for (size_t k = 0; k < K; ++k) target_o[k] = opt.x[off + k];
                if (K >= 2) sort_pairs(target_o, target_w, K - 1);
            }
        };

        auto try_rescue_two_set_collapsed_classes = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr) {
            if (K != 3) return;
            auto rescue_partition = [&](bool is_con, std::vector<Scalar>& cur_o, std::vector<Scalar>& cur_w, bool is_test) {
                bool collapsed = (cur_w[0] < 0.005) || (cur_w[1] < 0.005) || (std::abs(cur_o[1] - cur_o[0]) < 1e-3);
                if (!collapsed) return;
                Scalar base_ll = evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                size_t dom = (cur_w[0] >= cur_w[1]) ? 0 : 1;
                Scalar w_pur = std::clamp(cur_w[0] + cur_w[1], 0.04, 0.99);
                Scalar w_top = std::clamp(1.0 - w_pur, 0.01, 0.96);
                Scalar o_dom = std::clamp(cur_o[dom], 0.002, 0.80);
                Scalar o_mg = std::clamp(mg94_omega, 0.01, 0.80);
                Scalar o_top = is_con ? 1.0 : std::max(1.05, cur_o[2]);
                std::vector<std::pair<std::vector<Scalar>, std::vector<Scalar>>> rescue_cands = {
                    {{0.0, std::min(0.95, std::max(0.02, o_dom * 2.2)), o_top}, {0.65 * w_pur, 0.35 * w_pur, w_top}},
                    {{0.0, std::min(0.95, std::max(0.02, o_mg * 0.5)), o_top},  {0.65 * w_pur, 0.35 * w_pur, w_top}},
                    {{0.0, std::min(0.95, o_mg), o_top},                        {0.70 * w_pur, 0.30 * w_pur, w_top}}
                };
                for (const auto& [rc_o, rc_w] : rescue_cands) {
                    const auto& ot = is_test ? rc_o : omegas_t;
                    const auto& wt = is_test ? rc_w : weights_t;
                    const auto& ob = is_test ? omegas_b : rc_o;
                    const auto& wb = is_test ? weights_b : rc_w;
                    Scalar ll = evaluate_two_set_log_l(cur_s, ot, wt, ob, wb, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                    if (ll > base_ll + 1e-4) {
                        base_ll = ll;
                        cur_o = rc_o;
                        cur_w = rc_w;
                    }
                }
            };
            rescue_partition(constrain_test, omegas_t, weights_t, true);
            rescue_partition(constrain_bg, omegas_b, weights_b, false);
        };

        auto optimize_two_set_mh = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr, int max_it) {
            if (allow_triple) {
                auto mh_opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                    return -evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, x[0], x[1], eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, { std::clamp(delta, 0.0, 1000.0), std::clamp(psi, 0.0, 1000.0) }, {0.0, 0.0}, {1000.0, 1000.0}, 1e-4, max_it);
                delta = std::max(0.0, mh_opt.x[0]); psi = std::max(0.0, mh_opt.x[1]);
            } else if (allow_double) {
                auto [d_opt, nll_opt] = Brent1D::minimize([&](Scalar d) {
                    return -evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, d, 0.0, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, 0.0, std::clamp(delta, 0.0, 1000.0), 1000.0, 1e-4, max_it);
                delta = std::max(0.0, d_opt);
            }
        };

        auto optimize_two_set_gtr_theta = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr) {
            Scalar cur_ll = evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
            std::vector<Scalar> x0 = {
                std::log(std::clamp(local_base.theta_AC, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_AT, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_CG, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_CT, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_GT, 1e-4, 1e4))
            };
            std::vector<Scalar> lb(5, std::log(1e-4)), ub(5, std::log(1e4));
            auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                MG94Parameters cand = local_base;
                cand.theta_AC = std::exp(x[0]);
                cand.theta_AT = std::exp(x[1]);
                cand.theta_CG = std::exp(x[2]);
                cand.theta_CT = std::exp(x[3]);
                cand.theta_GT = std::exp(x[4]);
                return -evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &cand);
            }, x0, lb, ub, 1e-4, 40);
            if (-opt.f > cur_ll + 1e-5) {
                local_base.theta_AC = std::exp(opt.x[0]);
                local_base.theta_AT = std::exp(opt.x[1]);
                local_base.theta_CG = std::exp(opt.x[2]);
                local_base.theta_CT = std::exp(opt.x[3]);
                local_base.theta_GT = std::exp(opt.x[4]);
            }
        };

        auto scan_two_set_eps_sink = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr) {
            if (!use_error_sink) return;
            Scalar best_et = eps_sink_t;
            Scalar best_ll = evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
            for (Scalar cand_e : {1e-7, 2e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2}) {
                if (cand_e > max_eps_sink) continue;
                Scalar ll_e = evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, cand_e, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                if (ll_e > best_ll + 1e-4) {
                    best_ll = ll_e;
                    best_et = cand_e;
                }
            }
            eps_sink_t = best_et;
            Scalar best_eb = eps_sink_b;
            for (Scalar cand_e : {1e-7, 2e-4, 5e-4, 1e-3, 2.5e-3, 5e-3, 1e-2}) {
                if (cand_e > max_eps_sink) continue;
                Scalar ll_e = evaluate_two_set_log_l(cur_s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, bl_ptr, nullptr, delta, psi, eps_sink_t, cand_e, p_omegas_u, p_weights_u, &local_base);
                if (ll_e > best_ll + 1e-4) {
                    best_ll = ll_e;
                    best_eb = cand_e;
                }
            }
            eps_sink_b = best_eb;
        };

        int max_ecm = (K == 1 && !srv && !allow_double && !use_error_sink) ? 1 : 15;
        Scalar prev_ecm_ll = evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
        for (int ecm_iter = 0; ecm_iter < max_ecm; ++ecm_iter) {
            for (int em = 0; em < 2; ++em) {
                std::vector<Scalar> wt_u;
                auto [wt_t, wt_b] = compute_expected_two_set_weights(
                    s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, ecm_bl_ptr, delta, psi,
                    eps_sink_t, eps_sink_b,
                    use_error_sink ? &eps_sink_t : nullptr,
                    use_error_sink ? &eps_sink_b : nullptr,
                    max_eps_sink,
                    p_omegas_u, p_weights_u, has_u ? &wt_u : nullptr, &local_base
                );
                weights_t = wt_t;
                weights_b = wt_b;
                if (has_u && wt_u.size() == K) {
                    weights_u = wt_u;
                }
            }
            scan_two_set_eps_sink(s, ecm_bl_ptr);

            if (srv && M >= 2) {
                for (int em = 0; em < 2; ++em) {
                    syn_weights = compute_expected_two_set_syn_weights(
                        s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, ecm_bl_ptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                    );
                }
                renormalize_two_set_srv(s, nullptr);
            }

            // Optimize Test omegas + scale s
            optimize_partition_omegas(constrain_test, omegas_t, weights_t, s, ecm_bl_ptr, true, 45);
            // Optimize Background omegas + scale s
            optimize_partition_omegas(constrain_bg, omegas_b, weights_b, s, ecm_bl_ptr, true, 45);
            // Optimize Unclassified omegas
            if (has_u) {
                optimize_partition_omegas(false, omegas_u, weights_u, s, ecm_bl_ptr, true, 35);
            }

            if (ecm_iter == 1 || ecm_iter == 4) {
                try_rescue_two_set_collapsed_classes(s, ecm_bl_ptr);
            }

            if (srv && M >= 2) {
                optimize_two_set_srv_rates(s, nullptr, 30);
            }
            if (allow_double) {
                optimize_two_set_mh(s, ecm_bl_ptr, 25);
            }

            Scalar cur_ecm_ll = evaluate_two_set_log_l(s, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, ecm_bl_ptr, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << " (ECM " << (ecm_iter + 1) << "): lnL = "
                   << std::fixed << std::setprecision(2) << cur_ecm_ll;
                pbar->set_status(ss.str());
                pbar->tick();
            }
            if (ecm_iter >= 5 && std::abs(cur_ecm_ll - prev_ecm_ll) < 1e-3) {
                break;
            }
            prev_ecm_ll = cur_ecm_ll;
        }

        if (!refine_branches) {
            optimize_two_set_gtr_theta(s, ecm_bl_ptr);
        }

        std::vector<Scalar> bl(tree.num_nodes(), 0.0);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                bl[node.id] = s * (ecm_bl_ptr ? (*ecm_bl_ptr)[node.id] : branch_tau[node.id]);
            }
        }

        if (refine_branches) {
            int max_cycles = (K == 1 && !srv && !allow_double && !use_error_sink) ? 2 : 5;
            Scalar prev_cycle_ll = evaluate_two_set_log_l(1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
            for (int cycle = 0; cycle < max_cycles; ++cycle) {
                bl = refine_two_set_branch_lengths(bl, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, delta, psi, 25, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);

                if (K >= 2 || use_error_sink) {
                    std::vector<Scalar> wt_u;
                    auto [wt_t, wt_b] = compute_expected_two_set_weights(
                        1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, delta, psi,
                        eps_sink_t, eps_sink_b,
                        use_error_sink ? &eps_sink_t : nullptr,
                        use_error_sink ? &eps_sink_b : nullptr,
                        max_eps_sink,
                        p_omegas_u, p_weights_u, has_u ? &wt_u : nullptr, &local_base
                    );
                    weights_t = wt_t;
                    weights_b = wt_b;
                    if (has_u && wt_u.size() == K) {
                        weights_u = wt_u;
                    }
                    scan_two_set_eps_sink(1.0, &bl);
                }
                if (srv && M >= 2) {
                    syn_weights = compute_expected_two_set_syn_weights(1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }

                Scalar one_s = 1.0;
                optimize_partition_omegas(constrain_test, omegas_t, weights_t, one_s, &bl, false, 45);
                optimize_partition_omegas(constrain_bg, omegas_b, weights_b, one_s, &bl, false, 45);
                if (has_u) {
                    optimize_partition_omegas(false, omegas_u, weights_u, one_s, &bl, false, 35);
                }
                if (cycle == 0 || cycle == 2) {
                    try_rescue_two_set_collapsed_classes(1.0, &bl);
                }
                if (srv && M >= 2) {
                    optimize_two_set_srv_rates(one_s, &bl, 30);
                }
                if (allow_double) {
                    optimize_two_set_mh(1.0, &bl, 25);
                }
                optimize_two_set_gtr_theta(1.0, &bl);

                Scalar cur_cycle_ll = evaluate_two_set_log_l(1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                if (cycle >= 2 && std::abs(cur_cycle_ll - prev_cycle_ll) < 1e-3) {
                    break;
                }
                prev_cycle_ll = cur_cycle_ll;
            }
            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Branch lengths refined";
                pbar->set_status(ss.str());
                pbar->tick();
            }
            if (allow_double && pbar) {
                pbar->tick();
            }
        }

        if (srv && M >= 2) {
            Scalar dummy_s = s;
            renormalize_two_set_srv(dummy_s, refine_branches ? &bl : nullptr);
            if (!refine_branches) {
                s = dummy_s;
                for (const auto& node : tree.nodes) {
                    if (node.id != tree.root_id) {
                        bl[node.id] = s * (ecm_bl_ptr ? (*ecm_bl_ptr)[node.id] : branch_tau[node.id]);
                    }
                }
            }
        }

        Scalar best_lnL = evaluate_two_set_log_l(
            1.0, omegas_t, weights_t, omegas_b, weights_b, syn_rates, syn_weights, &bl, &fit.site_log_likelihoods, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
        );
        fit.log_likelihood = best_lnL;
        fit.tree_scale = s;
        fit.branch_lengths = bl;
        fit.delta = delta;
        fit.psi = psi;
        fit.error_sink_weight = eps_sink_t;
        fit.has_gtr_params = true;
        fit.gtr_params = local_base;
        if (local_base.syn_pair_rates != nullptr) {
            fit.mss_syn_pair_rates = *local_base.syn_pair_rates;
            fit.gtr_params.syn_pair_rates = &fit.mss_syn_pair_rates;
        }

        // Compute substitution rates and fractions
        Scalar total_rate = 0.0;
        Scalar double_rate = 0.0;
        Scalar triple_rate = 0.0;
        for (size_t k = 0; k < K; ++k) {
            MG94Parameters p = local_base;
            p.alpha = 1.0;
            p.beta = omegas_t[k];
            p.delta = delta;
            p.psi = psi;
            MG94Matrix mat;
            mat.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            auto sr = mat.compute_substitution_rates(*code);
            total_rate += weights_t[k] * sr.total_rate;
            double_rate += weights_t[k] * sr.rate_double;
            triple_rate += weights_t[k] * sr.rate_triple;
        }
        fit.frac_delta = (total_rate > 1e-12) ? (double_rate / total_rate) : 0.0;
        fit.frac_psi = (total_rate > 1e-12) ? (triple_rate / total_rate) : 0.0;

        fit.test_distribution.omegas = omegas_t;
        fit.test_distribution.weights = weights_t;
        fit.test_distribution.error_sink_weight = eps_sink_t;
        fit.background_distribution.omegas = omegas_b;
        fit.background_distribution.weights = weights_b;
        fit.background_distribution.error_sink_weight = eps_sink_b;
        if (has_u) {
            fit.unclassified_distribution.omegas = omegas_u;
            fit.unclassified_distribution.weights = weights_u;
            fit.unclassified_distribution.error_sink_weight = 0.0;
        }
        if (srv) {
            fit.test_distribution.syn_rates = syn_rates;
            fit.test_distribution.syn_weights = syn_weights;
            fit.background_distribution.syn_rates = syn_rates;
            fit.background_distribution.syn_weights = syn_weights;
            if (has_u) {
                fit.unclassified_distribution.syn_rates = syn_rates;
                fit.unclassified_distribution.syn_weights = syn_weights;
            }
        }

        size_t N_codons = aln.num_codons;
        size_t num_branches = tree.num_nodes() - 1;
        size_t num_mixture_params = (2 * K - 1 - (constrain_test && K >= 2 ? 1 : 0))
                                  + (2 * K - 1 - (constrain_bg && K >= 2 ? 1 : 0))
                                  + (has_u ? (2 * K - 1) : 0)
                                  + (use_error_sink ? 2 : 0);
        size_t num_syn_params = srv ? 2 * (M - 1) : 0;
        size_t num_mh_params = (allow_double ? 1 : 0) + (allow_triple ? 1 : 0);
        size_t K_params = 14 + (refine_branches ? num_branches : 1) + num_mixture_params + num_syn_params + num_mh_params + num_mss_free_params;
        fit.aicc = -2.0 * fit.log_likelihood + 2.0 * K_params * N_codons / std::max(1.0, static_cast<double>(N_codons) - static_cast<double>(K_params) - 1.0);

        return fit;
    }

    // Fits coupled-k RELAX two-set model:
    //   Reference (Background, partition 1): K-class distribution (omegas_b[c], weights[c])
    //   Test (Foreground, partition 0):      K-class distribution (omegas_t[c] = (omegas_b[c])^k, weights[c])
    //   Unclassified (partition 2, if any):  Unconstrained K-class distribution (omegas_u[c], weights_u[c])
    // Warm-started from shared_null_fit (the exact k = 1 RELAX Null model), guaranteeing lnL_RELAX >= lnL_shared.
    BUSTEDFit fit_relax_coupled_model(
        const BUSTEDFit& shared_null_fit,
        size_t K = 3,
        bool refine_branches = true,
        bool srv = false,
        size_t num_syn_rates = 3,
        const std::string& multiple_hits = "None",
        ProgressBar* pbar = nullptr,
        const std::string& model_desc = "",
        bool use_error_sink = false,
        Scalar max_eps_sink = 0.01,
        Scalar* out_relax_k = nullptr
    ) const {
        BUSTEDFit fit = shared_null_fit;
        fit.num_rate_classes = K;
        fit.has_background_distribution = true;
        bool has_u = (num_unclassified_branches > 0);
        fit.has_unclassified_distribution = has_u;
        MG94Parameters local_base = shared_null_fit.has_gtr_params ? shared_null_fit.gtr_params : base_params;

        bool allow_double = (multiple_hits == "Double" || multiple_hits == "Double+Triple");
        bool allow_triple = (multiple_hits == "Double+Triple");
        Scalar max_omega = use_error_sink ? 100.0 : 10000.0;

        std::vector<Scalar> omegas_b = shared_null_fit.test_distribution.omegas;
        std::vector<Scalar> weights = shared_null_fit.test_distribution.weights;
        if (omegas_b.size() != K || weights.size() != K) {
            omegas_b.resize(K, 1.0);
            weights.assign(K, 1.0 / static_cast<Scalar>(K));
        }

        std::vector<Scalar> omegas_u, weights_u;
        if (has_u) {
            if (shared_null_fit.has_unclassified_distribution &&
                shared_null_fit.unclassified_distribution.omegas.size() == K &&
                shared_null_fit.unclassified_distribution.weights.size() == K) {
                omegas_u = shared_null_fit.unclassified_distribution.omegas;
                weights_u = shared_null_fit.unclassified_distribution.weights;
            } else {
                omegas_u = omegas_b;
                weights_u = weights;
            }
        }
        const std::vector<Scalar>* p_omegas_u = has_u ? &omegas_u : nullptr;
        const std::vector<Scalar>* p_weights_u = has_u ? &weights_u : nullptr;

        size_t M = srv ? num_syn_rates : 1;
        std::vector<Scalar> syn_rates = (srv && shared_null_fit.test_distribution.syn_rates.size() == M)
            ? shared_null_fit.test_distribution.syn_rates : std::vector<Scalar>(M, 1.0);
        std::vector<Scalar> syn_weights = (srv && shared_null_fit.test_distribution.syn_weights.size() == M)
            ? shared_null_fit.test_distribution.syn_weights : std::vector<Scalar>(M, 1.0 / static_cast<Scalar>(M));

        Scalar delta = shared_null_fit.delta;
        Scalar psi = shared_null_fit.psi;
        Scalar eps_sink_t = (shared_null_fit.test_distribution.error_sink_weight > 0.0)
            ? shared_null_fit.test_distribution.error_sink_weight
            : shared_null_fit.error_sink_weight;
        Scalar eps_sink_b = (shared_null_fit.has_background_distribution && shared_null_fit.background_distribution.error_sink_weight > 0.0)
            ? shared_null_fit.background_distribution.error_sink_weight
            : eps_sink_t;

        std::vector<Scalar> base_bl = shared_null_fit.branch_lengths;
        if (base_bl.size() != tree.num_nodes()) {
            base_bl.assign(tree.num_nodes(), 0.0);
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    base_bl[node.id] = shared_null_fit.tree_scale * branch_tau[node.id];
                }
            }
        }
        std::vector<Scalar> cur_bl = base_bl;
        auto get_scaled_bl = [&](Scalar scale_s) -> const std::vector<Scalar>* {
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    cur_bl[node.id] = scale_s * base_bl[node.id];
                }
            }
            return &cur_bl;
        };

        auto couple_omegas = [&](const std::vector<Scalar>& wb, Scalar kv) {
            std::vector<Scalar> wt(wb.size());
            for (size_t c = 0; c < wb.size(); ++c) {
                if (std::abs(wb[c] - 1.0) < 1e-12) {
                    wt[c] = 1.0;
                } else if (wb[c] <= 1e-12) {
                    wt[c] = (kv > 0.0) ? 1e-8 : 1.0;
                } else {
                    wt[c] = std::clamp(std::pow(wb[c], kv), 1e-8, max_omega);
                }
            }
            return wt;
        };

        Scalar s = 1.0;
        Scalar k_relax = 1.0;
        std::vector<Scalar> omegas_t = couple_omegas(omegas_b, k_relax);

        auto sort_rate_weight_pairs = [](std::vector<Scalar>& rates, std::vector<Scalar>& wts, size_t n_sort) {
            if (n_sort <= 1) return;
            std::vector<size_t> idx(n_sort);
            for (size_t i = 0; i < n_sort; ++i) idx[i] = i;
            std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return rates[a] < rates[b]; });
            std::vector<Scalar> r_tmp(n_sort), w_tmp(n_sort);
            for (size_t i = 0; i < n_sort; ++i) {
                r_tmp[i] = rates[idx[i]];
                w_tmp[i] = wts[idx[i]];
            }
            for (size_t i = 0; i < n_sort; ++i) {
                rates[i] = r_tmp[i];
                wts[i] = w_tmp[i];
            }
        };

        auto renormalize_relax_srv = [&](Scalar& cur_s, std::vector<Scalar>* bl_vec) {
            if (!srv || M < 2) return;
            sort_rate_weight_pairs(syn_rates, syn_weights, M);
            Scalar mean_syn = 0.0;
            for (size_t m = 0; m < M; ++m) mean_syn += syn_weights[m] * syn_rates[m];
            if (mean_syn > 1e-12) {
                for (size_t m = 0; m < M; ++m) syn_rates[m] /= mean_syn;
                if (bl_vec) {
                    for (auto& b_val : *bl_vec) b_val *= mean_syn;
                } else {
                    cur_s *= mean_syn;
                }
            }
        };

        auto optimize_relax_srv_rates = [&](Scalar& cur_s, const std::vector<Scalar>* bl_ptr, int max_it = 30) {
            if (!srv || M < 2) return;
            renormalize_relax_srv(cur_s, nullptr);
            if (M == 2) {
                Scalar w0 = std::max(1e-6, syn_weights[0]);
                Scalar w1 = std::max(1e-6, syn_weights[1]);
                Scalar r0_max = (1.0 - w1 * 1e-4) / w0;
                if (r0_max > 1e-3) {
                    Scalar r0_init = std::clamp(syn_rates[0], 1e-3, std::min(1.0, r0_max - 1e-4));
                    auto [r0_opt, nll_opt] = Brent1D::minimize([&](Scalar r0) {
                        Scalar r1 = std::max(1e-4, (1.0 - w0 * r0) / w1);
                        return -evaluate_two_set_log_l(1.0, omegas_t, weights, omegas_b, weights, {r0, r1}, syn_weights, bl_ptr ? bl_ptr : get_scaled_bl(cur_s), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                    }, 1e-3, r0_init, std::min(1.0, r0_max), 1e-3, max_it);
                    syn_rates[0] = r0_opt;
                    syn_rates[1] = std::max(1e-4, (1.0 - w0 * r0_opt) / w1);
                }
            } else {
                std::vector<Scalar> x0(M - 1), lb(M - 1, 1e-3), ub(M - 1, 10.0);
                for (size_t m = 0; m + 1 < M; ++m) x0[m] = std::clamp(syn_rates[m], 1e-3, 10.0);
                auto eval_srv = [&](const std::vector<Scalar>& r_free) -> Scalar {
                    Scalar rem = 1.0;
                    for (size_t m = 0; m + 1 < M; ++m) rem -= syn_weights[m] * r_free[m];
                    Scalar r_last = rem / std::max(1e-6, syn_weights[M - 1]);
                    if (r_last <= 1e-4 || r_last > 50.0) return 1e20;
                    std::vector<Scalar> cand = r_free;
                    cand.push_back(r_last);
                    return -evaluate_two_set_log_l(1.0, omegas_t, weights, omegas_b, weights, cand, syn_weights, bl_ptr ? bl_ptr : get_scaled_bl(cur_s), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                };
                auto opt_syn = NelderMeadDynamic::minimize(eval_srv, x0, lb, ub, 1e-3, max_it);
                if (opt_syn.f < 1e19) {
                    Scalar rem = 1.0;
                    for (size_t m = 0; m + 1 < M; ++m) {
                        syn_rates[m] = opt_syn.x[m];
                        rem -= syn_weights[m] * syn_rates[m];
                    }
                    syn_rates[M - 1] = std::max(1e-4, rem / std::max(1e-6, syn_weights[M - 1]));
                }
            }
            renormalize_relax_srv(cur_s, nullptr);
        };

        Scalar best_tracked_ll = evaluate_two_set_log_l(
            1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
        );

        auto optimize_relax_mh = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr, int max_it = 50) {
            if (!allow_double) return;
            const std::vector<Scalar>* active_bl = bl_ptr ? bl_ptr : get_scaled_bl(cur_s);
            if (allow_triple) {
                auto mh_opt = NelderMeadND<2>::minimize([&](const std::array<Scalar, 2>& x) {
                    return -evaluate_two_set_log_l(1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, active_bl, nullptr, x[0], x[1], eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, { std::clamp(delta, 0.0, 1000.0), std::clamp(psi, 0.0, 1000.0) }, {0.0, 0.0}, {1000.0, 1000.0}, 1e-4, max_it);
                if (-mh_opt.f >= best_tracked_ll) {
                    delta = std::max(0.0, mh_opt.x[0]);
                    psi = std::max(0.0, mh_opt.x[1]);
                    best_tracked_ll = -mh_opt.f;
                }
            } else {
                auto [d_opt, nll_opt] = Brent1D::minimize([&](Scalar d) {
                    return -evaluate_two_set_log_l(1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, active_bl, nullptr, d, 0.0, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, 0.0, std::clamp(delta, 0.0, 1000.0), 1000.0, 1e-4, max_it);
                if (-nll_opt >= best_tracked_ll) {
                    delta = std::max(0.0, d_opt);
                    best_tracked_ll = -nll_opt;
                }
            }
        };

        auto optimize_relax_gtr_theta = [&](Scalar cur_s, const std::vector<Scalar>* bl_ptr, int max_it = 35) {
            const std::vector<Scalar>* active_bl = bl_ptr ? bl_ptr : get_scaled_bl(cur_s);
            std::array<Scalar, 5> x0 = {
                std::log(std::clamp(local_base.theta_AC, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_AT, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_CG, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_CT, 1e-4, 1e4)),
                std::log(std::clamp(local_base.theta_GT, 1e-4, 1e4))
            };
            const Scalar lo = std::log(1e-4), hi = std::log(1e4);
            auto opt_gtr = NelderMeadND<5>::minimize([&](const std::array<Scalar, 5>& x) {
                MG94Parameters cand_p = local_base;
                cand_p.theta_AC = std::exp(x[0]);
                cand_p.theta_AT = std::exp(x[1]);
                cand_p.theta_CG = std::exp(x[2]);
                cand_p.theta_CT = std::exp(x[3]);
                cand_p.theta_GT = std::exp(x[4]);
                return -evaluate_two_set_log_l(1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, active_bl, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &cand_p);
            }, x0, {lo, lo, lo, lo, lo}, {hi, hi, hi, hi, hi}, 1e-4, max_it);
            if (-opt_gtr.f > best_tracked_ll) {
                local_base.theta_AC = std::exp(opt_gtr.x[0]);
                local_base.theta_AT = std::exp(opt_gtr.x[1]);
                local_base.theta_CG = std::exp(opt_gtr.x[2]);
                local_base.theta_CT = std::exp(opt_gtr.x[3]);
                local_base.theta_GT = std::exp(opt_gtr.x[4]);
                best_tracked_ll = -opt_gtr.f;
            }
        };

        // Initial 1D search on k_relax holding shared distribution fixed
        {
            std::vector<Scalar> k_grid = {0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0};
            Scalar best_k = k_relax;
            for (Scalar kc : k_grid) {
                auto wt_c = couple_omegas(omegas_b, kc);
                Scalar ll_c = evaluate_two_set_log_l(
                    1.0, wt_c, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                );
                if (ll_c > best_tracked_ll) {
                    best_tracked_ll = ll_c;
                    best_k = kc;
                }
            }
            auto [brent_k, brent_nll] = Brent1D::minimize([&](Scalar kc) {
                Scalar kc_clamped = std::clamp(kc, 0.0, 50.0);
                auto wt_c = couple_omegas(omegas_b, kc_clamped);
                return -evaluate_two_set_log_l(
                    1.0, wt_c, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                );
            }, 0.0, best_k, 50.0, 1e-4, 50);
            if (-brent_nll >= best_tracked_ll) {
                best_tracked_ll = -brent_nll;
                k_relax = std::clamp(brent_k, 0.0, 50.0);
            } else {
                k_relax = best_k;
            }
            omegas_t = couple_omegas(omegas_b, k_relax);
        }

        size_t B_t = std::max<size_t>(1, num_test_branches);
        size_t B_b = std::max<size_t>(1, num_background_branches);
        Scalar frac_t = static_cast<Scalar>(B_t) / static_cast<Scalar>(B_t + B_b);
        Scalar frac_b = static_cast<Scalar>(B_b) / static_cast<Scalar>(B_t + B_b);

        int max_ecm = (K == 1 && !srv && !allow_double && !use_error_sink) ? 1 : 15;
        Scalar prev_ecm_ll = best_tracked_ll;
        for (int ecm_iter = 0; ecm_iter < max_ecm; ++ecm_iter) {
            // 1. Shared weights E-step (F-098: update eps_sink_t and eps_sink_b independently)
            if (K >= 2 || use_error_sink) {
                Scalar cand_eps_t = eps_sink_t, cand_eps_b = eps_sink_b;
                std::vector<Scalar> wt_u_cand;
                auto [wt_t, wt_b] = compute_expected_two_set_weights(
                    1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), delta, psi,
                    eps_sink_t, eps_sink_b,
                    use_error_sink ? &cand_eps_t : nullptr,
                    use_error_sink ? &cand_eps_b : nullptr,
                    max_eps_sink,
                    p_omegas_u, p_weights_u, has_u ? &wt_u_cand : nullptr, &local_base
                );
                std::vector<Scalar> w_cand(K);
                Scalar sum_w = 0.0;
                for (size_t c = 0; c < K; ++c) {
                    w_cand[c] = std::clamp(frac_t * wt_t[c] + frac_b * wt_b[c], 1e-6, 1.0 - 1e-6);
                    sum_w += w_cand[c];
                }
                for (size_t c = 0; c < K; ++c) w_cand[c] /= sum_w;
                const std::vector<Scalar>* p_wu_cand = (has_u && wt_u_cand.size() == K) ? &wt_u_cand : p_weights_u;
                Scalar cand_ll = evaluate_two_set_log_l(
                    1.0, omegas_t, w_cand, omegas_b, w_cand, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, cand_eps_t, cand_eps_b, p_omegas_u, p_wu_cand, &local_base
                );
                if (cand_ll >= best_tracked_ll - 1e-9) {
                    weights = w_cand;
                    if (has_u && wt_u_cand.size() == K) {
                        weights_u = wt_u_cand;
                    }
                    eps_sink_t = cand_eps_t;
                    eps_sink_b = cand_eps_b;
                    best_tracked_ll = std::max(best_tracked_ll, cand_ll);
                }
            }

            // 2. SRV weights E-step
            if (srv && M >= 2) {
                auto q_cand = compute_expected_two_set_syn_weights(
                    1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                );
                Scalar cand_ll = evaluate_two_set_log_l(
                    1.0, omegas_t, weights, omegas_b, weights, syn_rates, q_cand, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                );
                if (cand_ll >= best_tracked_ll - 1e-9) {
                    syn_weights = q_cand;
                    best_tracked_ll = std::max(best_tracked_ll, cand_ll);
                }
            }

            // 3. Joint M-step on (s, omegas_b, k_relax)
            if (K == 1) {
                auto opt = NelderMeadND<3>::minimize([&](const std::array<Scalar, 3>& x) {
                    std::vector<Scalar> wb = {x[1]};
                    auto wt = couple_omegas(wb, x[2]);
                    return -evaluate_two_set_log_l(1.0, wt, weights, wb, weights, syn_rates, syn_weights, get_scaled_bl(x[0]), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, {s, omegas_b[0], k_relax}, {0.01, 0.0, 0.0}, {50.0, max_omega, 50.0}, 1e-4, 40);
                if (-opt.f >= best_tracked_ll) {
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; k_relax = opt.x[2];
                    best_tracked_ll = -opt.f;
                }
            } else if (K == 2) {
                std::vector<std::array<Scalar, 4>> extra_grid;
                if (ecm_iter == 0) {
                    extra_grid = {
                        {s, omegas_b[0], std::max(2.0, omegas_b[1]), 0.2},
                        {s, omegas_b[0], std::max(2.0, omegas_b[1]), 1.0},
                        {s, omegas_b[0], std::max(2.0, omegas_b[1]), 2.0}
                    };
                }
                auto opt = NelderMeadND<4>::minimize([&](const std::array<Scalar, 4>& x) {
                    std::vector<Scalar> wb = {x[1], x[2]};
                    auto wt = couple_omegas(wb, x[3]);
                    return -evaluate_two_set_log_l(1.0, wt, weights, wb, weights, syn_rates, syn_weights, get_scaled_bl(x[0]), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, {s, omegas_b[0], omegas_b[1], k_relax}, {0.01, 0.0, 1.0, 0.0}, {50.0, 1.0, max_omega, 50.0}, 1e-4, 45, extra_grid);
                if (-opt.f >= best_tracked_ll) {
                    s = opt.x[0]; omegas_b[0] = opt.x[1]; omegas_b[1] = opt.x[2]; k_relax = opt.x[3];
                    best_tracked_ll = -opt.f;
                }
            } else if (K == 3) {
                std::vector<std::array<Scalar, 5>> extra_grid;
                if (ecm_iter == 0) {
                    extra_grid = {
                        {s, omegas_b[0], omegas_b[1], std::max(2.0, omegas_b[2]), 0.2},
                        {s, omegas_b[0], omegas_b[1], std::max(2.0, omegas_b[2]), 1.0},
                        {s, omegas_b[0], omegas_b[1], std::max(2.0, omegas_b[2]), 2.0}
                    };
                }
                auto opt = NelderMeadND<5>::minimize([&](const std::array<Scalar, 5>& x) {
                    std::vector<Scalar> wb = {std::min(x[1], x[2]), std::max(x[1], x[2]), x[3]};
                    auto wt = couple_omegas(wb, x[4]);
                    return -evaluate_two_set_log_l(1.0, wt, weights, wb, weights, syn_rates, syn_weights, get_scaled_bl(x[0]), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, {s, omegas_b[0], omegas_b[1], omegas_b[2], k_relax}, {0.01, 0.0, 0.0, 1.0, 0.0}, {50.0, 1.0, 1.0, max_omega, 50.0}, 1e-4, 50, extra_grid);
                if (-opt.f >= best_tracked_ll) {
                    s = opt.x[0]; omegas_b[0] = std::min(opt.x[1], opt.x[2]); omegas_b[1] = std::max(opt.x[1], opt.x[2]); omegas_b[2] = opt.x[3]; k_relax = opt.x[4];
                    best_tracked_ll = -opt.f;
                }
            } else {
                std::vector<Scalar> x0(K + 2), lb(K + 2, 0.0), ub(K + 2, 1.0);
                x0[0] = s; lb[0] = 0.01; ub[0] = 50.0;
                for (size_t c = 0; c < K; ++c) x0[c + 1] = omegas_b[c];
                lb[K] = 1.0; ub[K] = max_omega;
                x0[K + 1] = k_relax; lb[K + 1] = 0.0; ub[K + 1] = 50.0;
                auto opt = NelderMeadDynamic::minimize([&](const std::vector<Scalar>& x) {
                    std::vector<Scalar> wb(x.begin() + 1, x.begin() + 1 + K);
                    std::sort(wb.begin(), wb.begin() + (K - 1));
                    auto wt = couple_omegas(wb, x[K + 1]);
                    return -evaluate_two_set_log_l(1.0, wt, weights, wb, weights, syn_rates, syn_weights, get_scaled_bl(x[0]), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base);
                }, x0, lb, ub, 1e-4, 55);
                if (-opt.f >= best_tracked_ll) {
                    s = opt.x[0];
                    for (size_t c = 0; c < K; ++c) omegas_b[c] = opt.x[c + 1];
                    std::sort(omegas_b.begin(), omegas_b.begin() + (K - 1));
                    k_relax = opt.x[K + 1];
                    best_tracked_ll = -opt.f;
                }
            }
            omegas_t = couple_omegas(omegas_b, k_relax);

            // 4. Optimize SRV rates and MH delta/psi (F-048)
            if (srv && M >= 2) {
                optimize_relax_srv_rates(s, nullptr, 30);
                best_tracked_ll = std::max(best_tracked_ll, evaluate_two_set_log_l(1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, get_scaled_bl(s), nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base));
            }
            if (allow_double) {
                optimize_relax_mh(s, nullptr, 50);
            }

            if (pbar && ecm_iter < 15) {
                std::ostringstream ss;
                ss << model_desc << " (ECM " << (ecm_iter + 1) << ", k="
                   << std::fixed << std::setprecision(2) << k_relax << "): lnL = " << best_tracked_ll;
                pbar->set_status(ss.str());
                pbar->tick();
            }
            if (ecm_iter >= 2 && std::abs(best_tracked_ll - prev_ecm_ll) < 1e-3) {
                break;
            }
            prev_ecm_ll = best_tracked_ll;
        }

        std::vector<Scalar> bl = *get_scaled_bl(s);

        if (refine_branches) {
            auto cand_bl = refine_two_set_branch_lengths(
                bl, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, delta, psi, 20, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
            );
            Scalar cand_bl_ll = evaluate_two_set_log_l(
                1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, &cand_bl, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
            );
            if (cand_bl_ll >= best_tracked_ll) {
                bl = cand_bl;
                best_tracked_ll = cand_bl_ll;
            }
            if (pbar) {
                std::ostringstream ss;
                ss << model_desc << ": Branch lengths refined";
                pbar->set_status(ss.str());
                pbar->tick();
            }

            // Re-optimize GTR exchangeability parameters (F-054)
            optimize_relax_gtr_theta(1.0, &bl, 35);

            // Final 1D polish on k_relax, SRV rates, and MH delta/psi with refined branch lengths (F-048)
            auto [brent_k, brent_nll] = Brent1D::minimize([&](Scalar kc) {
                Scalar kc_clamped = std::clamp(kc, 0.0, 50.0);
                auto wt_c = couple_omegas(omegas_b, kc_clamped);
                return -evaluate_two_set_log_l(
                    1.0, wt_c, weights, omegas_b, weights, syn_rates, syn_weights, &bl, nullptr, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
                );
            }, 0.0, k_relax, 50.0, 1e-4, 50);
            if (-brent_nll >= best_tracked_ll) {
                best_tracked_ll = -brent_nll;
                k_relax = std::clamp(brent_k, 0.0, 50.0);
                omegas_t = couple_omegas(omegas_b, k_relax);
            }
            if (srv && M >= 2) {
                Scalar one_s = 1.0;
                optimize_relax_srv_rates(one_s, &bl, 25);
                renormalize_relax_srv(one_s, &bl);
            }
            if (allow_double) {
                optimize_relax_mh(1.0, &bl, 50);
                if (pbar) pbar->tick();
            }
        } else {
            optimize_relax_gtr_theta(1.0, &bl, 35);
            if (srv && M >= 2) {
                Scalar dummy_s = 1.0;
                renormalize_relax_srv(dummy_s, &bl);
            }
        }

        Scalar best_lnL = evaluate_two_set_log_l(
            1.0, omegas_t, weights, omegas_b, weights, syn_rates, syn_weights, &bl, &fit.site_log_likelihoods, delta, psi, eps_sink_t, eps_sink_b, p_omegas_u, p_weights_u, &local_base
        );
        if (best_lnL < shared_null_fit.log_likelihood) {
            fit = shared_null_fit;
            k_relax = 1.0;
        } else {
            fit.log_likelihood = best_lnL;
            fit.tree_scale = s * shared_null_fit.tree_scale;
            fit.branch_lengths = bl;
            fit.delta = delta;
            fit.psi = psi;
            fit.error_sink_weight = eps_sink_t;
            fit.has_gtr_params = true;
            fit.gtr_params = local_base;
            fit.mss_class_rates = shared_null_fit.mss_class_rates;
            if (local_base.syn_pair_rates != nullptr) {
                fit.mss_syn_pair_rates = *local_base.syn_pair_rates;
                fit.gtr_params.syn_pair_rates = &fit.mss_syn_pair_rates;
            }

            if (delta > 0.0 || psi > 0.0) {
                MG94Parameters p_mh = local_base;
                p_mh.alpha = 1.0;
                Scalar sum_w = 0.0;
                Scalar mean_w = 0.0;
                for (size_t c = 0; c < K; ++c) {
                    mean_w += weights[c] * (frac_t * omegas_t[c] + frac_b * omegas_b[c]);
                    sum_w += weights[c];
                }
                p_mh.beta = (sum_w > 0.0) ? (mean_w / sum_w) : 1.0;
                p_mh.delta = delta;
                p_mh.psi = psi;
                MG94Matrix mg_mh;
                mg_mh.update(p_mh, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
                auto sr = mg_mh.compute_substitution_rates(*code);
                fit.frac_delta = (sr.total_rate > 0.0) ? (sr.rate_double / sr.total_rate) : 0.0;
                fit.frac_psi = (sr.total_rate > 0.0) ? (sr.rate_triple / sr.total_rate) : 0.0;
            } else {
                fit.frac_delta = 0.0;
                fit.frac_psi = 0.0;
            }

            fit.test_distribution.omegas = omegas_t;
            fit.test_distribution.weights = weights;
            fit.test_distribution.error_sink_weight = eps_sink_t;
            fit.background_distribution.omegas = omegas_b;
            fit.background_distribution.weights = weights;
            fit.background_distribution.error_sink_weight = eps_sink_b;
            if (has_u) {
                fit.unclassified_distribution.omegas = omegas_u;
                fit.unclassified_distribution.weights = weights_u;
                fit.unclassified_distribution.error_sink_weight = 0.0;
            }
            if (srv) {
                fit.test_distribution.syn_rates = syn_rates;
                fit.test_distribution.syn_weights = syn_weights;
                fit.background_distribution.syn_rates = syn_rates;
                fit.background_distribution.syn_weights = syn_weights;
                if (has_u) {
                    fit.unclassified_distribution.syn_rates = syn_rates;
                    fit.unclassified_distribution.syn_weights = syn_weights;
                }
            }
        }

        if (out_relax_k) *out_relax_k = k_relax;

        size_t N_codons = aln.num_codons;
        size_t num_branches = tree.num_nodes() - 1;
        size_t num_mixture_params = (2 * K - 1) + 1 + (has_u ? (2 * K - 1) : 0) + (use_error_sink ? 2 : 0); // shared weights + BG omegas + k_relax + (eps_sink_t, eps_sink_b)
        size_t num_syn_params = srv ? 2 * (M - 1) : 0;
        size_t num_mh_params = (allow_double ? 1 : 0) + (allow_triple ? 1 : 0);
        size_t K_params = 14 + (refine_branches ? num_branches : 1) + num_mixture_params + num_syn_params + num_mh_params + num_mss_free_params;
        fit.aicc = -2.0 * fit.log_likelihood + 2.0 * K_params * N_codons / std::max(1.0, static_cast<double>(N_codons) - static_cast<double>(K_params) - 1.0);

        return fit;
    }

    // Fits Shared Distribution model (D_T == D_B on Test U Background; independent K-class on Unclassified if present)
    BUSTEDFit fit_shared_distribution_model(
        size_t K,
        const BUSTEDSettings& settings,
        ProgressBar* pbar,
        const std::string& model_desc,
        const BUSTEDFit* warm_unc = nullptr
    ) const {
        if (num_unclassified_branches == 0) {
            auto saved_is_test = is_test_branch;
            auto saved_num_test = num_test_branches;
            num_test_branches = 0;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    bool in_shared = (branch_partition[node.id] <= 1);
                    is_test_branch[node.id] = in_shared;
                    if (in_shared) num_test_branches++;
                }
            }
            BUSTEDFit fit = fit_model(
                false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar, model_desc,
                settings.error_sink, settings.error_sink_max_weight,
                warm_unc
            );
            fit.has_background_distribution = true;
            fit.background_distribution = fit.test_distribution;
            is_test_branch = saved_is_test;
            num_test_branches = saved_num_test;
            return fit;
        }

        // When num_unclassified_branches > 0, (Test U Background) shares one K-class distribution (partition 0)
        // and Unclassified has its own unconstrained K-class distribution (partition 1).
        auto saved_part = branch_partition;
        auto saved_nt = num_test_branches;
        auto saved_nb = num_background_branches;
        auto saved_nu = num_unclassified_branches;

        num_test_branches = saved_nt + saved_nb;
        num_background_branches = saved_nu;
        num_unclassified_branches = 0;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                branch_partition[node.id] = (saved_part[node.id] <= 1) ? 0 : 1;
            }
        }

        const BUSTEDRateDistribution* w_t = warm_unc ? &warm_unc->test_distribution : nullptr;
        const BUSTEDRateDistribution* w_b = (warm_unc && warm_unc->has_unclassified_distribution)
            ? &warm_unc->unclassified_distribution : nullptr;
        const std::vector<Scalar>* w_bl = (warm_unc && !warm_unc->branch_lengths.empty())
            ? &warm_unc->branch_lengths : nullptr;
        const MG94Parameters* w_gtr = (warm_unc && warm_unc->has_gtr_params)
            ? &warm_unc->gtr_params : nullptr;
        Scalar w_d = warm_unc ? warm_unc->delta : -1.0;
        Scalar w_p = warm_unc ? warm_unc->psi : -1.0;

        BUSTEDFit fit = fit_two_set_model(
            false, false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
            settings.multiple_hits, pbar, model_desc,
            settings.error_sink, settings.error_sink_max_weight,
            w_t, w_b, nullptr, w_bl, w_gtr, w_d, w_p
        );
        fit.has_unclassified_distribution = true;
        fit.unclassified_distribution = fit.background_distribution;
        fit.has_background_distribution = true;
        fit.background_distribution = fit.test_distribution;

        branch_partition = saved_part;
        num_test_branches = saved_nt;
        num_background_branches = saved_nb;
        num_unclassified_branches = saved_nu;
        return fit;
    }

    // Computes per-branch, per-site posterior probabilities and EBFs for omega_K > 1 (F-023)
    // and per-branch 2H / 3H / 2H+3H evidence ratios when multiple_hits != "None" (F-060)
    // in a single O(patterns * branches) Inside-Outside pass.
    void compute_branch_site_diagnostics(
        const BUSTEDFit& fit,
        const BUSTEDSettings& settings,
        BUSTEDResult& res
    ) const {
        res.branch_posterior_probabilities.clear();
        res.branch_empirical_bayes_factors.clear();
        res.branch_er_2h.clear();
        res.branch_er_3h.clear();
        res.branch_er_23h.clear();

        size_t num_sites = aln.num_codons;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        int S = code->num_sense_codons;
        size_t K_t = fit.test_distribution.omegas.size();
        if (K_t == 0 || fit.branch_lengths.size() != num_nodes) return;

        bool has_bg = fit.has_background_distribution && !fit.background_distribution.omegas.empty();
        bool has_uncl = fit.has_unclassified_distribution && !fit.unclassified_distribution.omegas.empty();
        size_t K_b = has_bg ? fit.background_distribution.omegas.size() : K_t;
        size_t K_u = has_uncl ? fit.unclassified_distribution.omegas.size() : K_t;

        std::vector<Scalar> syn_rates = fit.test_distribution.syn_rates.empty() ? std::vector<Scalar>{1.0} : fit.test_distribution.syn_rates;
        std::vector<Scalar> syn_weights = fit.test_distribution.syn_weights.empty() ? std::vector<Scalar>{1.0} : fit.test_distribution.syn_weights;
        size_t M = syn_rates.size();

        const MG94Parameters& active_base = fit.has_gtr_params ? fit.gtr_params : base_params;
        Scalar eps_t = settings.error_sink ? fit.test_distribution.error_sink_weight : 0.0;
        Scalar eps_b = settings.error_sink ? (has_bg ? fit.background_distribution.error_sink_weight : eps_t) : 0.0;
        Scalar bio_t = 1.0 - eps_t;
        Scalar bio_b = 1.0 - eps_b;

        bool allow_double = (settings.multiple_hits == "Double" || settings.multiple_hits == "Double+Triple");
        bool allow_triple = (settings.multiple_hits == "Double+Triple");

        // Build per-class MG94 matrices for Test, Background, Unclassified (and 1H/2H/no2H variants if MH enabled)
        auto build_mg_vec = [&](const std::vector<Scalar>& omegas, Scalar d_val, Scalar p_val) {
            std::vector<MG94Matrix> mats(omegas.size());
            for (size_t k = 0; k < omegas.size(); ++k) {
                MG94Parameters p = active_base;
                p.alpha = 1.0;
                p.beta = omegas[k];
                p.delta = d_val;
                p.psi = p_val;
                mats[k].update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *code);
            }
            return mats;
        };

        auto mg_t = build_mg_vec(fit.test_distribution.omegas, fit.delta, fit.psi);
        auto mg_b = has_bg ? build_mg_vec(fit.background_distribution.omegas, fit.delta, fit.psi) : mg_t;
        auto mg_u = has_uncl ? build_mg_vec(fit.unclassified_distribution.omegas, fit.delta, fit.psi) : mg_t;

        std::vector<MG94Matrix> mg_t_1h, mg_b_1h, mg_u_1h;
        std::vector<MG94Matrix> mg_t_2h, mg_b_2h, mg_u_2h;
        std::vector<MG94Matrix> mg_t_no2h, mg_b_no2h, mg_u_no2h;
        if (allow_double) {
            mg_t_1h = build_mg_vec(fit.test_distribution.omegas, 0.0, 0.0);
            mg_b_1h = has_bg ? build_mg_vec(fit.background_distribution.omegas, 0.0, 0.0) : mg_t_1h;
            mg_u_1h = has_uncl ? build_mg_vec(fit.unclassified_distribution.omegas, 0.0, 0.0) : mg_t_1h;
            if (allow_triple) {
                mg_t_2h = build_mg_vec(fit.test_distribution.omegas, fit.delta, 0.0);
                mg_b_2h = has_bg ? build_mg_vec(fit.background_distribution.omegas, fit.delta, 0.0) : mg_t_2h;
                mg_u_2h = has_uncl ? build_mg_vec(fit.unclassified_distribution.omegas, fit.delta, 0.0) : mg_t_2h;
                mg_t_no2h = build_mg_vec(fit.test_distribution.omegas, 0.0, fit.psi);
                mg_b_no2h = has_bg ? build_mg_vec(fit.background_distribution.omegas, 0.0, fit.psi) : mg_t_no2h;
                mg_u_no2h = has_uncl ? build_mg_vec(fit.unclassified_distribution.omegas, 0.0, fit.psi) : mg_t_no2h;
            }
        }

        std::vector<std::vector<Matrix>> P_mix(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<std::array<Matrix, 3>>> P_bridges_mix(M, std::vector<std::array<Matrix, 3>>(tree.bridges.size()));
        std::vector<std::vector<Matrix>> P_top_t(M, std::vector<Matrix>(num_nodes));
        std::vector<std::vector<Matrix>> P_mix_1h, P_mix_2h, P_mix_no2h;
        if (allow_double) {
            P_mix_1h.assign(M, std::vector<Matrix>(num_nodes));
            if (allow_triple) {
                P_mix_2h.assign(M, std::vector<Matrix>(num_nodes));
                P_mix_no2h.assign(M, std::vector<Matrix>(num_nodes));
            }
        }

        for (size_t m = 0; m < M; ++m) {
            Scalar alpha_m = syn_rates[m];
            auto build_branch_mix_at = [&](int part, Scalar bl,
                                           const std::vector<MG94Matrix>& m_t,
                                           const std::vector<MG94Matrix>& m_b,
                                           const std::vector<MG94Matrix>& m_u) {
                Matrix P = Matrix::Zero(S, S);
                if (part == 0 || !has_bg) {
                    Scalar scale = (part == 0) ? bio_t : 1.0;
                    for (size_t k = 0; k < K_t; ++k) {
                        P += (scale * fit.test_distribution.weights[k]) * m_t[k].transition_matrix(alpha_m * bl);
                    }
                    if (part == 0 && eps_t > 0.0) P.noalias() += eps_t * P_sink;
                } else if (part == 1) {
                    for (size_t k = 0; k < K_b; ++k) {
                        P += (bio_b * fit.background_distribution.weights[k]) * m_b[k].transition_matrix(alpha_m * bl);
                    }
                    if (eps_b > 0.0) P.noalias() += eps_b * P_sink;
                } else if (has_uncl) {
                    for (size_t k = 0; k < K_u; ++k) {
                        P += fit.unclassified_distribution.weights[k] * m_u[k].transition_matrix(alpha_m * bl);
                    }
                }
                return P;
            };

            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                int nid = node.id;
                Scalar bl = fit.branch_lengths[nid];
                int part = (static_cast<size_t>(nid) < branch_partition.size()) ? branch_partition[nid] : (is_test_branch[nid] ? 0 : 1);

                P_mix[m][nid] = build_branch_mix_at(part, bl, mg_t, mg_b, mg_u);
                if (is_test_branch[nid]) {
                    P_top_t[m][nid] = mg_t[K_t - 1].transition_matrix(alpha_m * bl);
                }
                if (allow_double) {
                    P_mix_1h[m][nid] = build_branch_mix_at(part, bl, mg_t_1h, mg_b_1h, mg_u_1h);
                    if (allow_triple) {
                        P_mix_2h[m][nid] = build_branch_mix_at(part, bl, mg_t_2h, mg_b_2h, mg_u_2h);
                        P_mix_no2h[m][nid] = build_branch_mix_at(part, bl, mg_t_no2h, mg_b_no2h, mg_u_no2h);
                    }
                }
            }
            for (size_t b = 0; b < tree.bridges.size(); ++b) {
                int bid = tree.bridges[b].bridge_node_id;
                int part = (static_cast<size_t>(bid) < branch_partition.size()) ? branch_partition[bid] : (is_test_branch[bid] ? 0 : 1);
                for (int r = 0; r < 3; ++r) {
                    Scalar bl = bridge_resolution_tau(b, r, 1.0, &fit.branch_lengths);
                    P_bridges_mix[m][b][r] = build_branch_mix_at(part, bl, mg_t, mg_b, mg_u);
                }
            }
        }

        Scalar prior_pos = (fit.test_distribution.omegas.back() > 1.0)
            ? (bio_t * fit.test_distribution.weights.back())
            : 0.0;

        std::vector<std::vector<Scalar>> pat_post_pos(num_patterns, std::vector<Scalar>(num_nodes, 0.0));
        std::vector<std::vector<Scalar>> pat_er_2h, pat_er_3h, pat_er_23h;
        if (allow_double) {
            pat_er_2h.assign(num_patterns, std::vector<Scalar>(num_nodes, 1.0));
            if (allow_triple) {
                pat_er_3h.assign(num_patterns, std::vector<Scalar>(num_nodes, 1.0));
                pat_er_23h.assign(num_patterns, std::vector<Scalar>(num_nodes, 1.0));
            }
        }

        const int n_threads = LikelihoodEngine::optimal_omp_threads(num_patterns);
        #pragma omp parallel num_threads(n_threads)
        {
            std::vector<LikelihoodEngine::InsideOutsideWorkspace> io_m(M);
            #pragma omp for schedule(static)
            for (size_t p = 0; p < num_patterns; ++p) {
                const auto& pattern = aln.patterns[p];
                Scalar max_ll = -std::numeric_limits<Scalar>::infinity();
                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] <= 1e-12) continue;
                    LikelihoodEngine::compute_inside_outside_in_place(
                        tree, pattern, leaf_to_taxon, P_mix[m], aln.codon_frequencies_f3x4, S, io_m[m],
                        tree.bridges.empty() ? nullptr : &P_bridges_mix[m]
                    );
                    if (io_m[m].likelihood > 0.0) {
                        Scalar term = std::log(syn_weights[m]) + io_m[m].log_likelihood;
                        if (term > max_ll) max_ll = term;
                    }
                }
                if (max_ll <= -1e19) continue;

                Scalar sum_exp = 0.0;
                std::vector<Scalar> post_m(M, 0.0);
                for (size_t m = 0; m < M; ++m) {
                    if (syn_weights[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                        Scalar e = std::exp(std::log(syn_weights[m]) + io_m[m].log_likelihood - max_ll);
                        post_m[m] = e;
                        sum_exp += e;
                    }
                }
                if (sum_exp <= 0.0) continue;
                for (size_t m = 0; m < M; ++m) post_m[m] /= sum_exp;

                for (const auto& node : tree.nodes) {
                    if (node.id == tree.root_id) continue;
                    int nid = node.id;
                    if (is_test_branch[nid] && prior_pos > 0.0) {
                        Scalar ratio_top = 0.0;
                        for (size_t m = 0; m < M; ++m) {
                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                ratio_top += post_m[m] * (io_m[m].V[nid].dot(P_top_t[m][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                            }
                        }
                        pat_post_pos[p][nid] = std::clamp(prior_pos * ratio_top, 0.0, 1.0 - 1e-12);
                    }
                    if (allow_double) {
                        Scalar ratio_1h = 0.0;
                        for (size_t m = 0; m < M; ++m) {
                            if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                ratio_1h += post_m[m] * (io_m[m].V[nid].dot(P_mix_1h[m][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                            }
                        }
                        if (allow_triple) {
                            Scalar ratio_2h = 0.0;
                            Scalar ratio_no2h = 0.0;
                            for (size_t m = 0; m < M; ++m) {
                                if (post_m[m] > 1e-12 && io_m[m].likelihood > 0.0) {
                                    ratio_2h += post_m[m] * (io_m[m].V[nid].dot(P_mix_2h[m][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                    ratio_no2h += post_m[m] * (io_m[m].V[nid].dot(P_mix_no2h[m][nid] * io_m[m].D[nid]) / io_m[m].likelihood);
                                }
                            }
                            // Match HyPhy BUSTED.bf (F-060):
                            //   "2H":    L(delta, psi) / L(0, psi)
                            //   "3H":    L(delta, psi) / L(delta, 0)
                            //   "2H+3H": L(delta, psi) / L(0, 0)
                            pat_er_2h[p][nid] = 1.0 / std::max(1e-12, ratio_no2h);
                            pat_er_3h[p][nid] = 1.0 / std::max(1e-12, ratio_2h);
                            pat_er_23h[p][nid] = 1.0 / std::max(1e-12, ratio_1h);
                        } else {
                            pat_er_2h[p][nid] = 1.0 / std::max(1e-12, ratio_1h);
                        }
                    }
                }
            }
        }

        Scalar prior_odds = (prior_pos > 1e-9 && prior_pos < 1.0 - 1e-9) ? (prior_pos / (1.0 - prior_pos)) : 0.0;
        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id || node.name.empty()) continue;
            int nid = node.id;
            if (is_test_branch[nid]) {
                std::vector<Scalar> b_post(num_sites, 0.0);
                std::vector<Scalar> b_ebf(num_sites, 0.0);
                if (prior_odds > 0.0) {
                    for (size_t s = 0; s < num_sites; ++s) {
                        size_t p = aln.site_to_pattern[s];
                        Scalar post = pat_post_pos[p][nid];
                        b_post[s] = post;
                        Scalar post_odds = post / std::max(1e-12, 1.0 - post);
                        b_ebf[s] = post_odds / prior_odds;
                    }
                }
                res.branch_posterior_probabilities[node.name] = std::move(b_post);
                res.branch_empirical_bayes_factors[node.name] = std::move(b_ebf);
            }
            if (allow_double) {
                std::vector<Scalar> b_2h(num_sites, 1.0);
                for (size_t s = 0; s < num_sites; ++s) {
                    b_2h[s] = pat_er_2h[aln.site_to_pattern[s]][nid];
                }
                res.branch_er_2h[node.name] = std::move(b_2h);
                if (allow_triple) {
                    std::vector<Scalar> b_3h(num_sites, 1.0);
                    std::vector<Scalar> b_23h(num_sites, 1.0);
                    for (size_t s = 0; s < num_sites; ++s) {
                        b_3h[s] = pat_er_3h[aln.site_to_pattern[s]][nid];
                        b_23h[s] = pat_er_23h[aln.site_to_pattern[s]][nid];
                    }
                    res.branch_er_3h[node.name] = std::move(b_3h);
                    res.branch_er_23h[node.name] = std::move(b_23h);
                }
            }
        }
    }

    // Run complete BUSTED or BUSTED-PH analysis with options
    BUSTEDResult run(
        BUSTEDSettings settings = {},
        bool show_progress = false,
        bool force_progress = false
    ) const {
        auto t0 = std::chrono::high_resolution_clock::now();
        select_branch_partitions(settings.test_branches, settings.comparison_branches, settings.phenotype_mode || settings.relax_mode);
        if ((settings.phenotype_mode || settings.relax_mode) && num_background_branches == 0) {
            throw std::runtime_error(
                "BUSTED-PH / RELAX mode requires a proper subset of foreground/test branches and a non-empty background/reference set. "
                "Specify foreground branches via --branches <spec> or label branches with {TEST}/{Foreground} in the tree."
            );
        }

        BUSTEDResult res;
        bool has_mss = (!settings.mss_type.empty() && MSSPartition::to_upper(settings.mss_type) != "NONE" && MSSPartition::to_upper(settings.mss_type) != "NO");
        MSSPartition mss_part;
        MSSNormalization mss_norm = MSSNormalization::Flux;
        std::vector<Scalar> mss_class_rates;
        if (has_mss) {
            mss_part = MSSPartition::build(settings.mss_type, *code, settings.mss_neutral);
            MSSAnalyzer mss_helper(code);
            MSSAnalyzer::GeneContext gctx;
            gctx.filepath = "";
            gctx.aln = aln;
            gctx.tree = tree;
            gctx.leaf_to_taxon = leaf_to_taxon;
            gctx.gtr_params = gtr_rates;
            gctx.gtr_log_l = gtr_log_l;
            gctx.mg94_params = base_params;
            gctx.mg94_params.syn_pair_rates = nullptr;
            gctx.mg94_log_l = mg94_log_l;
            gctx.mg94_omega = mg94_omega;
            gctx.omega = mg94_omega;
            gctx.tree_scaler = 1.0;
            gctx.branch_tau = branch_tau;
            mss_helper.genes.push_back(std::move(gctx));

            std::string unorm = MSSPartition::to_upper(settings.mss_normalization);
            if (unorm == "UNWEIGHTED") mss_norm = MSSNormalization::Unweighted;
            else if (unorm == "NEUTRAL_REF" || unorm == "NEUTRAL" || mss_part.neutral_class_id >= 0) mss_norm = MSSNormalization::NeutralRef;

            mss_class_rates.assign(mss_part.num_classes(), 1.0);
            std::vector<Scalar> g_scalers = {1.0};
            std::vector<Scalar> g_omegas = {mg94_omega};
            mss_helper.fit_partition(mss_part, mss_norm, mss_class_rates, g_scalers, g_omegas, false, 5);

            syn_pair_rates_storage = mss_part.expand_to_matrix(code->num_sense_codons, mss_class_rates);
            base_params.syn_pair_rates = &syn_pair_rates_storage;
            num_mss_free_params = (mss_part.num_classes() > 1) ? (mss_part.num_classes() - 1) : 0;

            std::vector<size_t> pair_counts(mss_part.num_classes(), 0);
            for (const auto& sp : mss_part.pairs) {
                if (sp.class_id >= 0 && static_cast<size_t>(sp.class_id) < pair_counts.size()) {
                    pair_counts[sp.class_id]++;
                }
            }
            for (size_t k = 0; k < mss_part.num_classes(); ++k) {
                BUSTEDMSSClassRate cr;
                cr.class_name = mss_part.class_names[k];
                cr.rate = mss_class_rates[k];
                cr.num_pairs = pair_counts[k];
                cr.is_reference = (mss_norm == MSSNormalization::NeutralRef && static_cast<int>(k) == mss_part.neutral_class_id);
                res.mss_rates.push_back(std::move(cr));
            }
        } else {
            base_params.syn_pair_rates = nullptr;
            num_mss_free_params = 0;
        }
        if (settings.full_auto) {
            res.auto_diagnostics = diagnose_and_select_model(settings);
        }
        res.settings = settings;

        bool allow_double = (settings.multiple_hits == "Double" || settings.multiple_hits == "Double+Triple");
        auto count_steps = [&](size_t K) -> size_t {
            size_t max_ecm = (K == 1 && !settings.srv && !allow_double && !settings.error_sink) ? 1 : 20;
            return max_ecm + (settings.refine_branch_lengths ? 1 : 0) + ((settings.refine_branch_lengths && allow_double) ? 1 : 0);
        };
        auto count_relax_steps = [&](size_t K) -> size_t {
            size_t max_ecm = (K == 1 && !settings.srv && !allow_double && !settings.error_sink) ? 1 : 15;
            return max_ecm + (settings.refine_branch_lengths ? 1 : 0) + ((settings.refine_branch_lengths && allow_double) ? 1 : 0);
        };

        bool run_relax_in_ph = settings.phenotype_mode && (settings.include_relax_in_ph || settings.relax_mode);
        size_t total_steps = 0;
        if (settings.phenotype_mode) {
            total_steps = count_steps(settings.num_rate_classes) * 4
                        + (run_relax_in_ph ? count_relax_steps(settings.num_rate_classes) : 0);
        } else if (settings.auto_select_k) {
            for (size_t k = 1; k <= settings.max_k; ++k) {
                total_steps += count_steps(k);
            }
            total_steps += count_steps(settings.max_k);
            if (settings.relax_mode) {
                total_steps += count_steps(settings.max_k) + count_relax_steps(settings.max_k);
            }
        } else {
            total_steps = count_steps(settings.num_rate_classes) * 2
                        + (settings.relax_mode ? (count_steps(settings.num_rate_classes) + count_relax_steps(settings.num_rate_classes)) : 0);
        }

        std::unique_ptr<ProgressBar> pbar;
        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar = std::make_unique<ProgressBar>(
                total_steps,
                settings.phenotype_mode ? "[BUSTED-PH] Mixture Models" : "[BUSTED] Mixture Models",
                "steps",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        // Helper to re-polish MSS synonymous class rates under any BUSTED or BUSTED-PH mixture fit (F-026).
        // Each BUSTEDFit owns its own mss_class_rates and mss_syn_pair_rates buffer so polishing a constrained
        // or nested fit never aliases or overwrites another fit's gtr_params.syn_pair_rates.
        auto polish_mss_for_fit = [&](BUSTEDFit& target_fit, bool is_two_set, bool update_res_mss) {
            if (!has_mss || mss_part.num_classes() <= 1) return;
            const int S_codons = code->num_sense_codons;
            std::vector<Scalar> rates_vec;
            if (!target_fit.mss_class_rates.empty()) {
                rates_vec = target_fit.mss_class_rates;
            } else if (target_fit.has_gtr_params && target_fit.gtr_params.syn_pair_rates != nullptr &&
                       target_fit.gtr_params.syn_pair_rates->size() == static_cast<size_t>(S_codons * S_codons)) {
                rates_vec.assign(mss_part.num_classes(), 1.0);
                const auto& mat = *target_fit.gtr_params.syn_pair_rates;
                for (const auto& sp : mss_part.pairs) {
                    if (sp.class_id >= 0 && static_cast<size_t>(sp.class_id) < rates_vec.size()) {
                        rates_vec[sp.class_id] = mat[sp.u * S_codons + sp.v];
                    }
                }
            } else {
                rates_vec = mss_class_rates;
            }
            const MG94Parameters* active_base = target_fit.has_gtr_params ? &target_fit.gtr_params : &base_params;
            std::vector<Scalar> syn_r = target_fit.test_distribution.syn_rates.empty() ? std::vector<Scalar>{1.0} : target_fit.test_distribution.syn_rates;
            std::vector<Scalar> syn_w = target_fit.test_distribution.syn_weights.empty() ? std::vector<Scalar>{1.0} : target_fit.test_distribution.syn_weights;
            const std::vector<Scalar>* p_om_u = target_fit.has_unclassified_distribution ? &target_fit.unclassified_distribution.omegas : nullptr;
            const std::vector<Scalar>* p_wt_u = target_fit.has_unclassified_distribution ? &target_fit.unclassified_distribution.weights : nullptr;
            Scalar cur_mss_ll = target_fit.log_likelihood;
            bool mss_improved = false;
            for (size_t c = 0; c < mss_part.num_classes(); ++c) {
                if (mss_norm == MSSNormalization::NeutralRef && static_cast<int>(c) == mss_part.neutral_class_id) continue;
                if (mss_norm != MSSNormalization::NeutralRef && c == 0) continue; // fix class 0 as gauge reference during 1D polish
                Scalar init_r = rates_vec[c];
                auto [r_opt, nll_opt] = Brent1D::minimize([&](Scalar cand_r) {
                    std::vector<Scalar> cand_rates = rates_vec;
                    cand_rates[c] = cand_r;
                    std::vector<Scalar> cand_mat = mss_part.expand_to_matrix(S_codons, cand_rates);
                    MG94Parameters p_cand = *active_base;
                    p_cand.syn_pair_rates = &cand_mat;
                    if (is_two_set) {
                        return -evaluate_two_set_log_l(
                            1.0,
                            target_fit.test_distribution.omegas, target_fit.test_distribution.weights,
                            target_fit.background_distribution.omegas, target_fit.background_distribution.weights,
                            syn_r, syn_w, &target_fit.branch_lengths, nullptr,
                            target_fit.delta, target_fit.psi,
                            target_fit.test_distribution.error_sink_weight, target_fit.background_distribution.error_sink_weight,
                            p_om_u, p_wt_u, &p_cand
                        );
                    }
                    return -evaluate_log_l(
                        1.0, target_fit.test_distribution.omegas, target_fit.test_distribution.weights,
                        syn_r, syn_w, &target_fit.branch_lengths, nullptr,
                        target_fit.delta, target_fit.psi, target_fit.error_sink_weight, &p_cand
                    );
                }, std::max(0.05, init_r * 0.25), init_r, std::min(20.0, std::max(3.0, init_r * 4.0)), 1e-3, 15);
                if (-nll_opt > cur_mss_ll + 1e-3) {
                    rates_vec[c] = r_opt;
                    cur_mss_ll = -nll_opt;
                    mss_improved = true;
                }
            }
            target_fit.mss_class_rates = rates_vec;
            if (mss_improved || target_fit.mss_syn_pair_rates.empty()) {
                target_fit.mss_syn_pair_rates = mss_part.expand_to_matrix(S_codons, rates_vec);
            }
            target_fit.has_gtr_params = true;
            target_fit.gtr_params.syn_pair_rates = &target_fit.mss_syn_pair_rates;
            if (update_res_mss) {
                mss_class_rates = rates_vec;
                syn_pair_rates_storage = target_fit.mss_syn_pair_rates;
                base_params.syn_pair_rates = &syn_pair_rates_storage;
                for (size_t c = 0; c < mss_part.num_classes() && c < res.mss_rates.size(); ++c) {
                    res.mss_rates[c].rate = rates_vec[c];
                }
            }
            if (mss_improved) {
                Scalar old_ll = target_fit.log_likelihood;
                Scalar aicc_penalty = target_fit.aicc + 2.0 * old_ll;
                if (is_two_set) {
                    target_fit.log_likelihood = evaluate_two_set_log_l(
                        1.0,
                        target_fit.test_distribution.omegas, target_fit.test_distribution.weights,
                        target_fit.background_distribution.omegas, target_fit.background_distribution.weights,
                        syn_r, syn_w, &target_fit.branch_lengths, &target_fit.site_log_likelihoods,
                        target_fit.delta, target_fit.psi,
                        target_fit.test_distribution.error_sink_weight, target_fit.background_distribution.error_sink_weight,
                        p_om_u, p_wt_u, &target_fit.gtr_params
                    );
                } else {
                    target_fit.log_likelihood = evaluate_log_l(
                        1.0, target_fit.test_distribution.omegas, target_fit.test_distribution.weights,
                        syn_r, syn_w, &target_fit.branch_lengths, &target_fit.site_log_likelihoods,
                        target_fit.delta, target_fit.psi, target_fit.error_sink_weight, &target_fit.gtr_params
                    );
                }
                target_fit.aicc = -2.0 * target_fit.log_likelihood + aicc_penalty;
            }
        };

        if (settings.phenotype_mode) {
            size_t K = settings.num_rate_classes;
            res.optimal_k = K;
            res.has_background_test = true;
            res.shared_df = static_cast<int>(2 * K - 1);

            // 1. Fit Unconstrained Two-Set Model first (F-007: M_unconstrained: omega_K^(T) >= 1, omega_K^(B) >= 1)
            res.unconstrained = fit_two_set_model(
                false, false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Unconstrained Two-Set K=" + std::to_string(K),
                settings.error_sink, settings.error_sink_max_weight
            );
            polish_mss_for_fit(res.unconstrained, true, true);
            const BUSTEDRateDistribution* warm_u = res.unconstrained.has_unclassified_distribution
                ? &res.unconstrained.unclassified_distribution : nullptr;
            const std::vector<Scalar>* warm_bl = !res.unconstrained.branch_lengths.empty()
                ? &res.unconstrained.branch_lengths : nullptr;
            const MG94Parameters* warm_gtr = res.unconstrained.has_gtr_params
                ? &res.unconstrained.gtr_params : nullptr;

            // 2. Fit Constrained Test Model warm-started from Unconstrained MLE (F-007: M_null_T: omega_K^(T) = 1, omega_K^(B) >= 1)
            res.constrained = fit_two_set_model(
                true, false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Constrained Test K=" + std::to_string(K),
                settings.error_sink, settings.error_sink_max_weight,
                &res.unconstrained.test_distribution, &res.unconstrained.background_distribution, warm_u,
                warm_bl, warm_gtr, res.unconstrained.delta, res.unconstrained.psi
            );
            polish_mss_for_fit(res.constrained, true, false);

            // 3. Fit Constrained Background Model warm-started from Unconstrained MLE (F-007: M_null_B: omega_K^(T) >= 1, omega_K^(B) = 1)
            res.constrained_background = fit_two_set_model(
                false, true, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Constrained Background K=" + std::to_string(K),
                settings.error_sink, settings.error_sink_max_weight,
                &res.unconstrained.test_distribution, &res.unconstrained.background_distribution, warm_u,
                warm_bl, warm_gtr, res.unconstrained.delta, res.unconstrained.psi
            );
            polish_mss_for_fit(res.constrained_background, true, false);

            // 4. Fit Shared Distributions Model (M_shared: D_T == D_B on Test U Background; also exact RELAX Null k = 1)
            res.shared_distributions = fit_shared_distribution_model(
                K, settings, pbar.get(), "Shared Dist (RELAX Null k=1) K=" + std::to_string(K), &res.unconstrained
            );
            polish_mss_for_fit(res.shared_distributions, !res.shared_distributions.background_distribution.omegas.empty(), false);

            // 4b. Fit Coupled-k RELAX Alternative Model (M_RELAX: omega_FG = omega_BG^k, shared weights)
            if (run_relax_in_ph) {
                res.has_relax_test = true;
                res.relax_alternative = fit_relax_coupled_model(
                    res.shared_distributions, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, pbar.get(), "RELAX Coupled-k K=" + std::to_string(K),
                    settings.error_sink, settings.error_sink_max_weight, &res.relax_k
                );
                polish_mss_for_fit(res.relax_alternative, true, false);
                res.relax_lrt = std::max(0.0, 2.0 * (res.relax_alternative.log_likelihood - res.shared_distributions.log_likelihood));
                res.relax_p_value = busted_chi2_p_value(res.relax_lrt, 1.0);
                if (res.relax_p_value <= settings.p_value_threshold) {
                    res.relax_direction = (res.relax_k > 1.0) ? "Intensification" : "Relaxation";
                } else {
                    res.relax_direction = "None";
                }
            }

            // Ensure proper nesting of log-likelihoods (M_shared, M_null_T, M_null_B, M_RELAX c M_unconstrained):
            // If any nested model exceeded M_unconstrained, re-optimize M_unconstrained warm-started from the best nested fit!
            const BUSTEDFit* best_nested = &res.constrained;
            if (res.constrained_background.log_likelihood > best_nested->log_likelihood) best_nested = &res.constrained_background;
            if (res.shared_distributions.log_likelihood > best_nested->log_likelihood) best_nested = &res.shared_distributions;
            if (res.has_relax_test && res.relax_alternative.log_likelihood > best_nested->log_likelihood) best_nested = &res.relax_alternative;

            if (best_nested->log_likelihood > res.unconstrained.log_likelihood + 1e-5) {
                const BUSTEDRateDistribution* retry_u = best_nested->has_unclassified_distribution
                    ? &best_nested->unclassified_distribution : warm_u;
                const std::vector<Scalar>* retry_bl = !best_nested->branch_lengths.empty()
                    ? &best_nested->branch_lengths : warm_bl;
                const MG94Parameters* retry_gtr = best_nested->has_gtr_params
                    ? &best_nested->gtr_params : warm_gtr;
                BUSTEDFit reopt_unc = fit_two_set_model(
                    false, false, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, nullptr, "",
                    settings.error_sink, settings.error_sink_max_weight,
                    &best_nested->test_distribution, &best_nested->background_distribution, retry_u,
                    retry_bl, retry_gtr, best_nested->delta, best_nested->psi
                );
                reopt_unc.mss_class_rates = best_nested->mss_class_rates;
                polish_mss_for_fit(reopt_unc, true, false);
                if (reopt_unc.log_likelihood > res.unconstrained.log_likelihood) {
                    res.unconstrained = std::move(reopt_unc);
                    polish_mss_for_fit(res.unconstrained, true, true);
                }
            }
            if (res.unconstrained.log_likelihood < best_nested->log_likelihood) {
                Scalar saved_aicc_diff = res.unconstrained.aicc + 2.0 * res.unconstrained.log_likelihood;
                res.unconstrained = *best_nested;
                res.unconstrained.aicc = -2.0 * res.unconstrained.log_likelihood + saved_aicc_diff;
                polish_mss_for_fit(res.unconstrained, true, true);
            }

            bool unc_failed = (res.unconstrained.log_likelihood <= -1e19 || !std::isfinite(res.unconstrained.log_likelihood));
            bool con_failed = (res.constrained.log_likelihood <= -1e19 || !std::isfinite(res.constrained.log_likelihood));
            bool bg_failed = (res.constrained_background.log_likelihood <= -1e19 || !std::isfinite(res.constrained_background.log_likelihood));
            bool shared_failed = (res.shared_distributions.log_likelihood <= -1e19 || !std::isfinite(res.shared_distributions.log_likelihood));

            if (unc_failed || con_failed || bg_failed || shared_failed) {
                res.is_failed = true;
                res.error_message = "Optimization failed: model log-likelihood collapsed to penalty threshold (-1e20). The fit did not converge.";
                res.lrt = 0.0;
                res.p_value = 1.0;
                res.background_lrt = 0.0;
                res.background_p_value = 1.0;
                res.shared_lrt = 0.0;
                res.shared_p_value = 1.0;
                res.union_intersection_p_value = 1.0;
                res.phenotype_associated = false;
                res.phenotype_verdict = "Analysis failed: likelihood collapsed or did not converge.";
            } else {
                // Test 1: Foreground / Test selection (LRT_T)
                res.lrt = std::max(0.0, 2.0 * (res.unconstrained.log_likelihood - res.constrained.log_likelihood));
                res.p_value = (res.lrt > 0.0) ? (0.5 * std::exp(-res.lrt / 2.0)) : 0.5;

                // Test 2: Background / Comparison selection (LRT_B)
                res.background_lrt = std::max(0.0, 2.0 * (res.unconstrained.log_likelihood - res.constrained_background.log_likelihood));
                res.background_p_value = (res.background_lrt > 0.0) ? (0.5 * std::exp(-res.background_lrt / 2.0)) : 0.5;

                // Test 3: Shared vs Separate Distributions (LRT_shared ~ chi^2_{2K-1})
                res.shared_lrt = std::max(0.0, 2.0 * (res.unconstrained.log_likelihood - res.shared_distributions.log_likelihood));
                res.shared_p_value = busted_chi2_p_value(res.shared_lrt, static_cast<Scalar>(res.shared_df));
            }

            Scalar alpha = settings.p_value_threshold;
            res.union_intersection_p_value = std::max(res.p_value, res.shared_p_value);
            res.background_delta_aicc = res.constrained_background.aicc - res.unconstrained.aicc;
            res.background_neutral_supported = (res.background_delta_aicc <= 2.0) || (res.background_p_value > 0.068);
            res.ranking_score = -std::log10(std::max(1e-8, res.p_value))
                                - std::log10(std::max(1e-8, res.shared_p_value))
                                + std::log10(std::max(1e-8, res.background_p_value));

            if (res.union_intersection_p_value <= alpha) {
                if (res.background_neutral_supported) {
                    res.phenotype_associated = true;
                    res.phenotype_verdict =
                        "The composite null hypothesis (no FG selection or no FG/BG difference) is rejected (p_UIT <= alpha), "
                        "and the neutral model for background branches is sufficiently supported by AIC-c: "
                        "statistical evidence that selection is associated with the trait.";
                } else {
                    res.phenotype_associated = false;
                    res.phenotype_verdict =
                        "The composite null hypothesis (no FG selection or no FG/BG difference) is rejected (p_UIT <= alpha), "
                        "but the neutral model for background branches is not supported by AIC-c: "
                        "selection is acting broadly on the tree, not just on branches with the trait.";
                }
            } else {
                res.phenotype_associated = false;
                res.phenotype_verdict =
                    "The composite null hypothesis for no selection on foreground or no difference between foreground and background "
                    "could not be rejected (p_UIT = max(p_FG, p_diff) > alpha): no statistical evidence that selection is associated with the trait.";
            }
        } else if (settings.auto_select_k) {
            auto [opt_k, null_fits] = select_optimal_k(
                settings.max_k,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes,
                settings.multiple_hits,
                pbar.get(),
                settings.error_sink,
                settings.error_sink_max_weight
            );
            res.optimal_k = opt_k;
            res.k_null_fits = null_fits;
            res.constrained = null_fits[opt_k - 1];
            polish_mss_for_fit(res.constrained, num_background_branches > 0, false);
            if (num_background_branches > 0) {
                const BUSTEDRateDistribution* warm_u = res.constrained.has_unclassified_distribution
                    ? &res.constrained.unclassified_distribution : nullptr;
                const std::vector<Scalar>* warm_bl = !res.constrained.branch_lengths.empty()
                    ? &res.constrained.branch_lengths : nullptr;
                const MG94Parameters* warm_gtr = res.constrained.has_gtr_params
                    ? &res.constrained.gtr_params : nullptr;
                res.unconstrained = fit_two_set_model(
                    false, false, opt_k, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, pbar.get(), "Unconstrained K=" + std::to_string(opt_k),
                    settings.error_sink, settings.error_sink_max_weight,
                    &res.constrained.test_distribution, &res.constrained.background_distribution, warm_u,
                    warm_bl, warm_gtr, res.constrained.delta, res.constrained.psi
                );
                polish_mss_for_fit(res.unconstrained, true, true);
                if (opt_k >= 2) {
                    const BUSTEDRateDistribution* alt_u = res.unconstrained.has_unclassified_distribution
                        ? &res.unconstrained.unclassified_distribution : nullptr;
                    const std::vector<Scalar>* alt_bl = !res.unconstrained.branch_lengths.empty()
                        ? &res.unconstrained.branch_lengths : nullptr;
                    const MG94Parameters* alt_gtr = res.unconstrained.has_gtr_params
                        ? &res.unconstrained.gtr_params : nullptr;
                    BUSTEDFit null_from_alt = fit_two_set_model(
                        true, false, opt_k, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                        settings.multiple_hits, nullptr, "", settings.error_sink, settings.error_sink_max_weight,
                        &res.unconstrained.test_distribution, &res.unconstrained.background_distribution, alt_u,
                        alt_bl, alt_gtr, res.unconstrained.delta, res.unconstrained.psi
                    );
                    polish_mss_for_fit(null_from_alt, true, false);
                    if (null_from_alt.log_likelihood > res.constrained.log_likelihood) {
                        res.constrained = std::move(null_from_alt);
                    }
                }
            } else {
                res.unconstrained = fit_model(
                    false,
                    opt_k,
                    settings.refine_branch_lengths,
                    settings.srv,
                    settings.num_syn_rate_classes,
                    settings.multiple_hits,
                    pbar.get(),
                    "Unconstrained K=" + std::to_string(opt_k),
                    settings.error_sink,
                    settings.error_sink_max_weight,
                    &res.constrained
                );
                polish_mss_for_fit(res.unconstrained, false, true);
                if (opt_k >= 2) {
                    BUSTEDFit null_from_alt = fit_model(
                        true, opt_k, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                        settings.multiple_hits, nullptr, "", settings.error_sink, settings.error_sink_max_weight, &res.unconstrained
                    );
                    polish_mss_for_fit(null_from_alt, false, false);
                    if (null_from_alt.log_likelihood > res.constrained.log_likelihood) {
                        res.constrained = std::move(null_from_alt);
                    }
                }
            }
        } else if (num_background_branches > 0) {
            // Subset of branches tested in standard BUSTED (F-003):
            // Background branches receive their own unconstrained K-class BS-REL mixture in both Alt and Null.
            // Fit Unconstrained first, then warm-start Constrained Null from Unconstrained MLE (F-007).
            size_t K_fit = settings.num_rate_classes;
            res.optimal_k = K_fit;
            res.unconstrained = fit_two_set_model(
                false, false, K_fit, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Unconstrained K=" + std::to_string(K_fit),
                settings.error_sink, settings.error_sink_max_weight
            );
            polish_mss_for_fit(res.unconstrained, true, true);
            const BUSTEDRateDistribution* warm_u = res.unconstrained.has_unclassified_distribution
                ? &res.unconstrained.unclassified_distribution : nullptr;
            const std::vector<Scalar>* warm_bl = !res.unconstrained.branch_lengths.empty()
                ? &res.unconstrained.branch_lengths : nullptr;
            const MG94Parameters* warm_gtr = res.unconstrained.has_gtr_params
                ? &res.unconstrained.gtr_params : nullptr;
            res.constrained = fit_two_set_model(
                true, false, K_fit, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                settings.multiple_hits, pbar.get(), "Constrained K=" + std::to_string(K_fit),
                settings.error_sink, settings.error_sink_max_weight,
                &res.unconstrained.test_distribution, &res.unconstrained.background_distribution, warm_u,
                warm_bl, warm_gtr, res.unconstrained.delta, res.unconstrained.psi
            );
            polish_mss_for_fit(res.constrained, true, false);
            if (res.constrained.log_likelihood > res.unconstrained.log_likelihood + 1e-5) {
                const std::vector<Scalar>* con_bl = !res.constrained.branch_lengths.empty()
                    ? &res.constrained.branch_lengths : warm_bl;
                const MG94Parameters* con_gtr = res.constrained.has_gtr_params
                    ? &res.constrained.gtr_params : warm_gtr;
                BUSTEDFit alt_rescued = fit_two_set_model(
                    false, false, K_fit, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, nullptr, "", settings.error_sink, settings.error_sink_max_weight,
                    &res.constrained.test_distribution, &res.constrained.background_distribution, warm_u,
                    con_bl, con_gtr, res.constrained.delta, res.constrained.psi
                );
                alt_rescued.mss_class_rates = res.constrained.mss_class_rates;
                polish_mss_for_fit(alt_rescued, true, false);
                if (alt_rescued.log_likelihood > res.unconstrained.log_likelihood) {
                    res.unconstrained = std::move(alt_rescued);
                    polish_mss_for_fit(res.unconstrained, true, true);
                }
            }
        } else {
            // Standard all-branch BUSTED:
            // Fit Unconstrained (Alternative) model first, then warm-start Constrained (Null) from Unconstrained MLE (F-007).
            size_t K_fit = settings.num_rate_classes;
            res.optimal_k = K_fit;
            res.unconstrained = fit_model(
                false,
                K_fit,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes,
                settings.multiple_hits,
                pbar.get(),
                "Unconstrained K=" + std::to_string(K_fit),
                settings.error_sink,
                settings.error_sink_max_weight
            );
            polish_mss_for_fit(res.unconstrained, false, true);

            res.constrained = fit_model(
                true,
                K_fit,
                settings.refine_branch_lengths,
                settings.srv,
                settings.num_syn_rate_classes,
                settings.multiple_hits,
                pbar.get(),
                "Constrained K=" + std::to_string(K_fit),
                settings.error_sink,
                settings.error_sink_max_weight,
                &res.unconstrained
            );
            polish_mss_for_fit(res.constrained, false, false);

            if (settings.srv && (2.0 * (res.unconstrained.log_likelihood - res.constrained.log_likelihood)) > 11.5) {
                BUSTEDFit con_indep = fit_model(
                    true,
                    K_fit,
                    settings.refine_branch_lengths,
                    settings.srv,
                    settings.num_syn_rate_classes,
                    settings.multiple_hits,
                    nullptr,
                    "",
                    settings.error_sink,
                    settings.error_sink_max_weight,
                    nullptr
                );
                polish_mss_for_fit(con_indep, false, false);
                if (con_indep.log_likelihood > res.constrained.log_likelihood) {
                    res.constrained = std::move(con_indep);
                }
            }

            // Bidirectional cross-check (F-007): if Constrained Null exceeded Unconstrained, re-optimize Unconstrained warm-started from Constrained Null
            for (int retry = 0; retry < 2 && res.constrained.log_likelihood > res.unconstrained.log_likelihood + 1e-5; ++retry) {
                BUSTEDFit alt_rescued = fit_model(
                    false, K_fit, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, nullptr, "", settings.error_sink, settings.error_sink_max_weight, &res.constrained
                );
                alt_rescued.mss_class_rates = res.constrained.mss_class_rates;
                polish_mss_for_fit(alt_rescued, false, false);
                if (alt_rescued.log_likelihood > res.unconstrained.log_likelihood) {
                    res.unconstrained = std::move(alt_rescued);
                    polish_mss_for_fit(res.unconstrained, false, true);
                } else {
                    break;
                }
            }
        }

        if (!settings.phenotype_mode) {
            bool unc_failed = (res.unconstrained.log_likelihood <= -1e19 || !std::isfinite(res.unconstrained.log_likelihood));
            bool con_failed = (res.constrained.log_likelihood <= -1e19 || !std::isfinite(res.constrained.log_likelihood));
            if (unc_failed || con_failed) {
                res.is_failed = true;
                res.error_message = "Optimization failed: model log-likelihood collapsed to penalty threshold (-1e20). The fit did not converge.";
                res.lrt = 0.0;
                res.p_value = 1.0;
            } else if (res.optimal_k == 1) {
                // K=1 boundary hypothesis test: H0: omega_0 <= 1 vs H1: omega_0 > 1 (F-081)
                if (res.unconstrained.test_distribution.omegas[0] <= 1.0 ||
                    res.unconstrained.log_likelihood <= res.constrained.log_likelihood) {
                    res.constrained = res.unconstrained;
                    res.lrt = 0.0;
                    res.p_value = 0.5;
                } else {
                    res.lrt = 2.0 * (res.unconstrained.log_likelihood - res.constrained.log_likelihood);
                    res.p_value = 0.5 * busted_chi2_p_value(res.lrt, 1.0);
                }
            } else if (res.unconstrained.log_likelihood <= res.constrained.log_likelihood) {
                // Reconcile full unconstrained fit with constrained MLE so parameters, site_log_likelihoods, and log_likelihood are 100% consistent (F-007)
                Scalar saved_aicc_penalty = res.unconstrained.aicc + 2.0 * res.unconstrained.log_likelihood;
                res.unconstrained = res.constrained;
                res.unconstrained.aicc = -2.0 * res.unconstrained.log_likelihood + saved_aicc_penalty;
                res.lrt = 0.0;
                res.p_value = 0.5;
            } else {
                res.lrt = 2.0 * (res.unconstrained.log_likelihood - res.constrained.log_likelihood);
                // Asymptotic null: 0.5 * chi^2_0 + 0.5 * chi^2_2
                res.p_value = 0.5 * std::exp(-res.lrt / 2.0);
            }

            // If --relax was requested on standard BUSTED, fit Shared (k=1) and RELAX (k free)
            if (settings.relax_mode && !res.is_failed) {
                size_t K = res.optimal_k;
                res.has_relax_test = true;
                res.shared_distributions = fit_shared_distribution_model(
                    K, settings, pbar.get(), "RELAX Null (k=1) K=" + std::to_string(K), &res.unconstrained
                );
                res.relax_alternative = fit_relax_coupled_model(
                    res.shared_distributions, K, settings.refine_branch_lengths, settings.srv, settings.num_syn_rate_classes,
                    settings.multiple_hits, pbar.get(), "RELAX Coupled-k K=" + std::to_string(K),
                    settings.error_sink, settings.error_sink_max_weight, &res.relax_k
                );
                res.relax_lrt = std::max(0.0, 2.0 * (res.relax_alternative.log_likelihood - res.shared_distributions.log_likelihood));
                res.relax_p_value = busted_chi2_p_value(res.relax_lrt, 1.0);
                if (res.relax_p_value <= settings.p_value_threshold) {
                    res.relax_direction = (res.relax_k > 1.0) ? "Intensification" : "Relaxation";
                } else {
                    res.relax_direction = "None";
                }
            }
        }

        if (pbar) {
            std::ostringstream summary;
            if (settings.phenotype_mode) {
                summary << (res.phenotype_associated ? "\033[1;32m" : "")
                        << (res.phenotype_associated ? "Phenotype-associated selection detected" : "No phenotype-specific selection")
                        << (res.phenotype_associated ? "\033[0m" : "")
                        << " (p_UIT=" << std::fixed << std::setprecision(4) << res.union_intersection_p_value
                        << ", ΔAICc_BG=" << std::setprecision(2) << res.background_delta_aicc
                        << ", p_BG=" << std::setprecision(4) << res.background_p_value;
                if (res.has_relax_test) {
                    summary << ", RELAX k=" << std::setprecision(2) << res.relax_k
                            << " [p=" << std::setprecision(4) << res.relax_p_value << "]";
                }
                summary << ")";
            } else if (res.p_value <= settings.p_value_threshold) {
                summary << "\033[1;32mEvidence of episodic diversifying selection\033[0m: LRT = "
                        << std::fixed << std::setprecision(2) << res.lrt
                        << ", p = " << std::setprecision(4) << res.p_value;
            } else {
                summary << "No evidence of diversifying selection: LRT = "
                        << std::fixed << std::setprecision(2) << res.lrt
                        << ", p = " << std::setprecision(4) << res.p_value;
            }
            pbar->finish(summary.str());
        }

        // Compute per-site Evidence Ratios (F-055):
        // 1. "optimized null" (res.evidence_ratios): exp(ell_alt,s - ell_null,s) without flooring at 1.0
        // 2. "constrained" (res.constrained_evidence_ratios): exp(ell_alt,s - ell_clamped,s) where omega_K := 1.0
        //    is substituted into the unconstrained fit holding all other unconstrained parameters fixed.
        size_t num_sites = aln.num_codons;
        res.evidence_ratios.assign(num_sites, 1.0);
        res.constrained_evidence_ratios.assign(num_sites, 1.0);

        if (!res.is_failed && !res.unconstrained.site_log_likelihoods.empty() && !res.constrained.site_log_likelihoods.empty()) {
            std::vector<Scalar> clamped_pat_ll(aln.patterns.size(), 0.0);
            std::vector<Scalar> clamped_omegas_t = res.unconstrained.test_distribution.omegas;
            if (!clamped_omegas_t.empty()) {
                clamped_omegas_t.back() = std::min<Scalar>(1.0, clamped_omegas_t.back());
            }
            std::vector<Scalar> syn_r = res.unconstrained.test_distribution.syn_rates.empty() ? std::vector<Scalar>{1.0} : res.unconstrained.test_distribution.syn_rates;
            std::vector<Scalar> syn_w = res.unconstrained.test_distribution.syn_weights.empty() ? std::vector<Scalar>{1.0} : res.unconstrained.test_distribution.syn_weights;
            const MG94Parameters* alt_base = res.unconstrained.has_gtr_params ? &res.unconstrained.gtr_params : &base_params;

            if (res.unconstrained.has_background_distribution && !res.unconstrained.background_distribution.omegas.empty()) {
                const std::vector<Scalar>* p_ou = res.unconstrained.has_unclassified_distribution ? &res.unconstrained.unclassified_distribution.omegas : nullptr;
                const std::vector<Scalar>* p_wu = res.unconstrained.has_unclassified_distribution ? &res.unconstrained.unclassified_distribution.weights : nullptr;
                evaluate_two_set_log_l(
                    1.0, clamped_omegas_t, res.unconstrained.test_distribution.weights,
                    res.unconstrained.background_distribution.omegas, res.unconstrained.background_distribution.weights,
                    syn_r, syn_w, &res.unconstrained.branch_lengths, &clamped_pat_ll,
                    res.unconstrained.delta, res.unconstrained.psi,
                    res.unconstrained.test_distribution.error_sink_weight,
                    res.unconstrained.background_distribution.error_sink_weight,
                    p_ou, p_wu, alt_base
                );
            } else {
                evaluate_log_l(
                    1.0, clamped_omegas_t, res.unconstrained.test_distribution.weights,
                    syn_r, syn_w, &res.unconstrained.branch_lengths, &clamped_pat_ll,
                    res.unconstrained.delta, res.unconstrained.psi,
                    res.unconstrained.error_sink_weight, alt_base
                );
            }

            for (size_t s_idx = 0; s_idx < num_sites; ++s_idx) {
                size_t p_idx = aln.site_to_pattern[s_idx];
                Scalar null_ll = res.constrained.site_log_likelihoods[p_idx];
                Scalar alt_ll = res.unconstrained.site_log_likelihoods[p_idx];
                Scalar clamp_ll = clamped_pat_ll[p_idx];

                res.evidence_ratios[s_idx] = std::exp(std::clamp(alt_ll - null_ll, -200.0, 200.0));
                res.constrained_evidence_ratios[s_idx] = std::exp(std::clamp(alt_ll - clamp_ll, -200.0, 200.0));
            }

            compute_branch_site_diagnostics(res.unconstrained, settings, res);
        }

        if (settings.error_sink) {
            compute_error_sink_diagnostics(
                res.unconstrained,
                settings.error_sink_bf_threshold,
                res.site_error_bfs,
                res.flagged_error_events
            );
            if (!settings.output_filtered_alignment.empty()) {
                write_filtered_alignment(settings.output_filtered_alignment, res);
            }
        } else {
            res.site_error_bfs.assign(num_sites, 0.0);
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        res.runtime_seconds = std::chrono::duration<double>(t1 - t0).count();
        return res;
    }

    // Masks branch-site error events flagged by BUSTED-E (BF >= error_sink_bf_threshold)
    // as "---" in all descendant leaf sequences of the flagged branch at that codon site,
    // matching HyPhy 2.5's error-filter.bf workflow.
    Alignment filter_alignment_by_error_sink(const BUSTEDResult& res) const {
        Alignment filtered;
        filtered.code = aln.code;
        filtered.embedded_tree_newick = aln.embedded_tree_newick;
        for (const auto& seq : aln.sequences) {
            filtered.add_sequence(seq.name, seq.raw_seq);
        }

        if (res.flagged_error_events.empty()) {
            filtered.finalize();
            return filtered;
        }

        // Precompute descendant leaf taxon indices for each tree node
        size_t num_nodes = tree.num_nodes();
        std::vector<std::vector<size_t>> node_descendant_taxa(num_nodes);
        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                if (t_idx != static_cast<size_t>(-1) && t_idx < filtered.sequences.size()) {
                    node_descendant_taxa[nid].push_back(t_idx);
                }
            } else {
                for (int32_t cid : node.children) {
                    const auto& child_taxa = node_descendant_taxa[cid];
                    node_descendant_taxa[nid].insert(node_descendant_taxa[nid].end(), child_taxa.begin(), child_taxa.end());
                }
            }
        }

        std::unordered_map<std::string, int32_t> name_to_nid;
        for (const auto& node : tree.nodes) {
            if (!node.name.empty()) {
                name_to_nid[node.name] = node.id;
            }
        }

        for (const auto& ev : res.flagged_error_events) {
            if (ev.site == 0 || ev.site > aln.num_codons) continue;
            size_t codon_idx = ev.site - 1;
            size_t nt_offset = codon_idx * 3;
            auto it = name_to_nid.find(ev.branch_name);
            if (it != name_to_nid.end()) {
                int32_t nid = it->second;
                for (size_t t_idx : node_descendant_taxa[nid]) {
                    if (nt_offset + 3 <= filtered.sequences[t_idx].raw_seq.size()) {
                        filtered.sequences[t_idx].raw_seq.replace(nt_offset, 3, "---");
                    }
                }
            }
        }

        filtered.finalize();
        return filtered;
    }

    // Writes the BUSTED-E error-filtered codon alignment to a FASTA file
    void write_filtered_alignment(const std::string& filepath, const BUSTEDResult& res) const {
        Alignment filtered = filter_alignment_by_error_sink(res);
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open filtered alignment file for writing: " + filepath);
        }
        for (const auto& seq : filtered.sequences) {
            out << ">" << seq.name << "\n" << seq.raw_seq << "\n";
        }
    }

nlohmann::json to_legacy_json(const BUSTEDResult& res) const;

nlohmann::json to_modern_json(const BUSTEDResult& res, const Provenance& prov = {}) const;
    nlohmann::json to_json(
        const BUSTEDResult& res,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        if (format == JSONFormat::Legacy) {
            return to_legacy_json(res);
        }
        return to_modern_json(res, prov);
    }

    void save_json(
        const std::string& filepath,
        const BUSTEDResult& res,
        JSONFormat format = JSONFormat::ModernV3,
        const Provenance& prov = {}
    ) const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open file for writing: " + filepath);
        }
        nlohmann::json j = to_json(res, format, prov);
        out << j.dump(2) << "\n";
    }

};

} // namespace hyphy::analyses
