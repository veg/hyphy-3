#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/analyses/site_visualizer.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <atomic>

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

struct MEMESiteResult {
    size_t site_index = 0;
    Scalar alpha = 0.0;
    Scalar beta1 = 0.0;
    Scalar p1 = 1.0;
    Scalar beta_plus = 0.0;
    Scalar p_plus = 0.0;
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    int branches_under_selection = 0;
    Scalar total_branch_length = 0.0;
    Scalar log_l_meme = 0.0;
    Scalar log_l_fel = 0.0;
    Scalar p_value_meme_vs_fel = 1.0;
    Scalar fel_alpha = 0.0;
    Scalar fel_beta = 0.0;
    std::vector<Scalar> branch_ebf;
    SiteSubstitutions subs;
};

class MEMEAnalyzer {
public:
    Tree tree;
    Alignment aln;
    MG94Parameters base_params;
    FELAnalyzer fel_analyzer;

    Scalar p_value_threshold = 0.1;
    Scalar baseline_tree_length = 0.0;
    Scalar conversion_factor = 1.0;
    Scalar global_log_l = 0.0;
    Scalar global_aicc = 0.0;
    std::vector<MEMESiteResult> site_results;

    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    std::string input_filepath;
    std::string tree_string;
    std::string test_branches = "All";
    std::vector<bool> is_test_branch;
    size_t num_tested = 0;
    size_t num_background = 0;

    MEMEAnalyzer(Tree t, Alignment a, MG94Parameters params = MG94Parameters{}, std::string branches_spec = "All")
        : tree(t), aln(std::move(a)), base_params(params), fel_analyzer(t, aln, params), test_branches(std::move(branches_spec)) {
        prepare_branch_lengths();
        select_test_branches(test_branches);
    }

    void select_test_branches(const std::string& spec) {
        test_branches = spec;
        auto sel = tree.select_branches(spec, "", {}, false);
        is_test_branch = std::move(sel.is_selected);
        num_tested = sel.num_selected;
        num_background = sel.num_unselected;
        fel_analyzer.settings.test_branches = spec;
        fel_analyzer.partition_branches();
    }

    static MEMEAnalyzer create_and_fit(
        Tree input_tree,
        Alignment aln,
        Scalar pvalue_threshold = 0.1,
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool full_model = true,
        std::string branches_spec = "All"
    ) {
        // Validate branch selection upfront before fitting GTR/MG94
        input_tree.select_branches(branches_spec, "", {}, false);

        if (progress_cb) progress_cb("Phase 1: Fitting Nucleotide GTR Model", 0.1);
        GTRFitter gtr_fitter(input_tree, aln);
        auto gtr_res = gtr_fitter.fit();

        if (progress_cb) progress_cb("Phase 2: Refining under Global MG94 Model", 0.3);
        MG94Parameters base_p;
        base_p.theta_AC = gtr_res.params.theta_AC;
        base_p.theta_AT = gtr_res.params.theta_AT;
        base_p.theta_CG = gtr_res.params.theta_CG;
        base_p.theta_CT = gtr_res.params.theta_CT;
        base_p.theta_GT = gtr_res.params.theta_GT;

        MG94Fitter mg_fitter(gtr_res.tree, aln);
        FitResult mg_res;
        if (full_model) {
            mg_res = mg_fitter.fit_full_model(1.0, base_p, progress_cb);
            base_p = mg_res.params;
        } else {
            mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
            base_p.alpha = 1.0;
            base_p.beta = mg_res.x_opt(0);
        }

        MEMEAnalyzer analyzer(mg_res.tree, std::move(aln), base_p, std::move(branches_spec));
        analyzer.p_value_threshold = pvalue_threshold;
        analyzer.global_log_l = mg_res.log_likelihood;
        size_t K_mg = 5 + 1 + (mg_res.tree.num_nodes() - 1);
        size_t N_codons = analyzer.aln.num_codons;
        analyzer.global_aicc = -2.0 * mg_res.log_likelihood + 2.0 * K_mg * N_codons / (N_codons - K_mg - 1.0);

        analyzer.has_gtr_fit = true;
        analyzer.gtr_log_l = gtr_res.log_likelihood;
        analyzer.gtr_aicc = gtr_res.aicc;
        analyzer.gtr_rates = gtr_res.params;
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                analyzer.gtr_branch_lengths[node.name] = node.branch_length;
            }
        }
        return analyzer;
    }

    void prepare_branch_lengths() {
        MG94Matrix mat;
        mat.update(base_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, *(aln.code ? aln.code : GeneticCode::universal()));
        Scalar Q_scale = mat.scale_factor;
        conversion_factor = (Q_scale > 1e-12) ? (3.0 / Q_scale) : 1.0;

        baseline_tree_length = 0.0;
        for (auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                baseline_tree_length += node.branch_length;
                node.branch_length *= conversion_factor;
            }
        }
        tree.scale_bridge_lengths(conversion_factor);
    }

    MEMESiteResult analyze_pattern(
        size_t p_idx,
        const std::vector<size_t>& leaf_to_taxon
    ) const {
        const auto& pattern = aln.patterns[p_idx];
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();

        MEMESiteResult res;
        res.site_index = p_idx;
        res.branch_ebf.assign(num_nodes, 1.0);

        // Run FEL on this pattern as the starting baseline
        SiteResult fel_res = fel_analyzer.analyze_pattern(p_idx, leaf_to_taxon);
        res.fel_alpha = fel_res.alpha;
        res.fel_beta = fel_res.beta;
        res.log_l_fel = fel_res.log_l_alt;

        // Invariant pattern check
        int8_t first_state = -1;
        bool is_invariant = true;
        for (int8_t s : pattern.states) {
            if (s >= 0) {
                if (first_state == -1) {
                    first_state = s;
                } else if (s != first_state) {
                    is_invariant = false;
                    break;
                }
            }
        }

        if (is_invariant) {
            res.alpha = 0.0;
            res.beta1 = 0.0;
            res.p1 = 1.0;
            res.beta_plus = 0.0;
            res.p_plus = 0.0;
            res.lrt = 0.0;
            res.p_value = 1.0;
            res.branches_under_selection = 0;
            res.log_l_meme = 0.0;
            res.log_l_fel = 0.0;
            res.p_value_meme_vs_fel = 1.0;
            res.total_branch_length = 0.0;
            return res;
        }

        auto [thread_node_L, thread_node_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S);
        std::vector<Matrix> P_branches(num_nodes);

        bool has_bg = (num_background > 0);
        auto eval_pattern_meme_lnl = [&](Scalar a, Scalar w1, Scalar p1, Scalar b_plus) -> Scalar {
            Scalar b1 = a * w1;
            MG94Parameters p_1 = base_params;
            p_1.alpha = a;
            p_1.beta = b1;

            MG94Parameters p_p = base_params;
            p_p.alpha = a;
            p_p.beta = b_plus;

            MG94Matrix mg1, mg_plus, mg_bg;
            mg1.update(p_1, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            mg_plus.update(p_p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            if (has_bg) {
                MG94Parameters p_bg = base_params;
                p_bg.alpha = a;
                p_bg.beta = (fel_res.beta_nuisance > 0.0) ? fel_res.beta_nuisance : b1;
                mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            }

            Scalar p_plus_wt = 1.0 - p1;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    bool is_test = (static_cast<size_t>(node.id) < is_test_branch.size()) ? is_test_branch[node.id] : true;
                    if (is_test) {
                        P_branches[node.id] = p1 * mg1.transition_matrix(node.branch_length) +
                                              p_plus_wt * mg_plus.transition_matrix(node.branch_length);
                    } else {
                        P_branches[node.id] = mg_bg.transition_matrix(node.branch_length);
                    }
                }
            }

            return LikelihoodEngine::compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg1, P_branches, thread_node_L, thread_node_scale);
        };

        // Construct grid starting points matching HyPhy 2.5 MEME.bf
        Scalar alpha_fel = std::clamp(fel_res.alpha, 1e-4, 100.0);
        Scalar beta_fel = std::clamp(fel_res.beta, 1e-4, 1000.0);
        std::vector<std::array<Scalar, 4>> start_grid;

        Scalar w1_0 = std::clamp(beta_fel / std::max(alpha_fel, 1e-4), 0.0, 1.0);

        if (alpha_fel >= beta_fel) {
            // Conserved Site Regime: targeted grid points matching HyPhy
            start_grid.push_back({alpha_fel, w1_0, 0.75, beta_fel});
            start_grid.push_back({alpha_fel, 0.0, 0.75, beta_fel});
            start_grid.push_back({alpha_fel, 0.0, 0.95, alpha_fel * 0.5});
            start_grid.push_back({alpha_fel, 0.25, 0.85, alpha_fel * 0.9});
        } else {
            // Positive Selection Regime: structured grid points matching HyPhy
            start_grid.push_back({alpha_fel, w1_0, 0.75, beta_fel});
            start_grid.push_back({alpha_fel, 0.0, 0.75, beta_fel});
            start_grid.push_back({alpha_fel, 0.25, 0.50, beta_fel * 2.0});
            start_grid.push_back({alpha_fel, 0.25, 0.75, beta_fel * 5.0});
            start_grid.push_back({alpha_fel, 0.10, 0.90, beta_fel * 15.0});
            start_grid.push_back({alpha_fel, 0.00, 0.95, beta_fel * 30.0});
            start_grid.push_back({alpha_fel, 1.00, 0.01, std::max(beta_fel, alpha_fel * 1.5)});
        }

        // Optimize Alternative Model (4D Nelder-Mead: alpha, omega1, p1, beta_plus)
        auto obj_alt = [&](const std::array<Scalar, 4>& x) -> Scalar {
            return -eval_pattern_meme_lnl(x[0], x[1], x[2], x[3]);
        };

        std::array<Scalar, 4> lb_alt = {1e-4, 0.0, 1e-6, 1e-4};
        std::array<Scalar, 4> ub_alt = {100.0, 1.0, 1.0 - 1e-6, 2000.0};
        auto opt_alt = NelderMeadND<4>::minimize(obj_alt, start_grid[0], lb_alt, ub_alt, 1e-4, 400, start_grid);

        res.alpha = opt_alt.x[0];
        Scalar omega1 = opt_alt.x[1];
        res.p1 = opt_alt.x[2];
        res.beta_plus = opt_alt.x[3];
        res.beta1 = res.alpha * omega1;
        res.p_plus = 1.0 - res.p1;
        res.log_l_meme = -opt_alt.f;

        // Hypothesis Testing: Alternative vs. Null (beta_plus == alpha)
        if (res.beta_plus > res.alpha && res.p_plus > 1e-6) {
            auto obj_null = [&](const std::array<Scalar, 3>& x) -> Scalar {
                return -eval_pattern_meme_lnl(x[0], x[1], x[2], x[0]); // beta_plus = alpha
            };

            std::vector<std::array<Scalar, 3>> null_grid;
            null_grid.push_back({res.alpha, omega1, res.p1});
            null_grid.push_back({std::clamp(fel_res.alpha_null, 1e-4, 50.0), 1.0, 0.5});
            null_grid.push_back({std::clamp(fel_res.alpha_null, 1e-4, 50.0), omega1, res.p1});
            null_grid.push_back({std::clamp(res.alpha, 1e-4, 50.0), 1.0, 0.5});
            for (Scalar a_cand : {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0}) {
                null_grid.push_back({a_cand, 0.5, 0.5});
                null_grid.push_back({a_cand, 0.0, 0.9});
            }

            std::array<Scalar, 3> lb_null = {1e-4, 0.0, 1e-6};
            std::array<Scalar, 3> ub_null = {100.0, 1.0, 1.0 - 1e-6};
            auto opt_null = NelderMeadND<3>::minimize(obj_null, null_grid[0], lb_null, ub_null, 1e-4, 400, null_grid);

            Scalar log_l_null = -opt_null.f;
            if (log_l_null <= -1e10 || res.log_l_meme <= -1e10 || !std::isfinite(log_l_null) || !std::isfinite(res.log_l_meme)) {
                res.lrt = 0.0;
                res.p_value = 1.0;
                res.branches_under_selection = 0;
            } else if (log_l_null > res.log_l_meme) {
                res.log_l_meme = log_l_null;
                res.lrt = 0.0;
                res.p_value = 2.0 / 3.0;
                res.branches_under_selection = 0;
            } else {
                Scalar lrt = 2.0 * (res.log_l_meme - log_l_null);
                res.lrt = lrt;

                // Asymptotic mixture p-value: 2/3 * [0.45 * erfc(sqrt(LRT/2)) + 0.55 * exp(-LRT/2)]
                res.p_value = (2.0 / 3.0) * (0.45 * std::erfc(std::sqrt(lrt / 2.0)) + 0.55 * std::exp(-lrt / 2.0));
            }

            // Compute branch Empirical Bayes Factors (EBF)
            MG94Parameters p_1 = base_params;
            p_1.alpha = res.alpha;
            p_1.beta = res.beta1;

            MG94Parameters p_p = base_params;
            p_p.alpha = res.alpha;
            p_p.beta = res.beta_plus;

            MG94Matrix mg1, mg_plus, mg_bg;
            mg1.update(p_1, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            mg_plus.update(p_p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            if (has_bg) {
                MG94Parameters p_bg = base_params;
                p_bg.alpha = res.alpha;
                p_bg.beta = (fel_res.beta_nuisance > 0.0) ? fel_res.beta_nuisance : res.beta1;
                mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            }

            std::vector<Matrix> P1(num_nodes), P_plus(num_nodes), P_mix(num_nodes);
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    bool is_test = (static_cast<size_t>(node.id) < is_test_branch.size()) ? is_test_branch[node.id] : true;
                    P1[node.id] = mg1.transition_matrix(node.branch_length);
                    P_plus[node.id] = mg_plus.transition_matrix(node.branch_length);
                    if (is_test) {
                        P_mix[node.id] = res.p1 * P1[node.id] + res.p_plus * P_plus[node.id];
                    } else {
                        P_mix[node.id] = mg_bg.transition_matrix(node.branch_length);
                    }
                }
            }

            int branches_under_sel = 0;
            std::vector<Matrix> P_temp = P_mix;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    bool is_test = (static_cast<size_t>(node.id) < is_test_branch.size()) ? is_test_branch[node.id] : true;
                    if (!is_test) {
                        res.branch_ebf[node.id] = 1.0;
                        continue;
                    }
                    // Conditional on branch node.id being in class 1
                    P_temp[node.id] = P1[node.id];
                    Scalar log_L1 = LikelihoodEngine::compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg1, P_temp, thread_node_L, thread_node_scale);

                    // Conditional on branch node.id being in class +
                    P_temp[node.id] = P_plus[node.id];
                    Scalar log_L_plus = LikelihoodEngine::compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg1, P_temp, thread_node_L, thread_node_scale);

                    P_temp[node.id] = P_mix[node.id]; // restore

                    Scalar ebf = 1.0;
                    if (log_L1 > -1e10 && log_L_plus > -1e10) {
                        ebf = std::exp(std::clamp(log_L_plus - log_L1, -50.0, 50.0));
                    } else if (log_L1 <= -1e10 && log_L_plus > -1e10) {
                        ebf = 1e6;
                    } else {
                        ebf = 0.0;
                    }
                    res.branch_ebf[node.id] = ebf;
                    if (ebf >= 100.0) {
                        branches_under_sel++;
                    }
                }
            }
            res.branches_under_selection = branches_under_sel;
        } else {
            res.lrt = 0.0;
            res.p_value = 2.0 / 3.0;
            res.branches_under_selection = 0;
        }

        // Compute MEME vs FEL LRT p-value (df = 2)
        Scalar delta_fel = 2.0 * (res.log_l_meme - res.log_l_fel);
        if (delta_fel > 0.0) {
            res.p_value_meme_vs_fel = std::exp(-delta_fel / 2.0);
        } else {
            res.p_value_meme_vs_fel = 1.0;
        }

        res.total_branch_length = 0.0;

        // Ancestral state reconstruction and substitution history
        res.subs = fel_analyzer.reconstruct_ancestral_pattern(p_idx, leaf_to_taxon);
        FELAnalyzer::format_substitution_strings(res.subs, false);

        return res;
    }

    std::vector<MEMESiteResult> run(bool show_progress = false, bool force_progress = false) {
        std::vector<size_t> leaf_to_taxon(tree.num_nodes(), static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        size_t num_patterns = aln.patterns.size();
        std::vector<MEMESiteResult> pattern_results(num_patterns);

        SiteVisualizer viz(aln.num_codons, num_patterns, SiteVisualizerConfig{
            .method_name = "MEME",
            .p_value_threshold = p_value_threshold,
            .pos_label = "episodic",
            .neg_label = "purifying",
            .sig_label = "significant",
            .quiet = !show_progress,
            .force_tty = force_progress,
            .enable_chart = true
        });

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            pattern_results[p] = analyze_pattern(p, leaf_to_taxon);

            const auto& r = pattern_results[p];
            SiteSelectionType type = SiteSelectionType::Neutral;
            Scalar metric = 0.0;
            std::string detail;
            if (r.p_value <= p_value_threshold && r.beta_plus > r.alpha) {
                type = SiteSelectionType::Positive;
                metric = r.beta_plus / std::max(r.alpha, 1e-4);
                detail = "▲ β+ = " + SiteVisualizer::format_scalar(r.beta_plus, 2) + " (p+ = " + SiteVisualizer::format_scalar(r.p_plus, 2) + "), p = " + SiteVisualizer::format_pvalue(r.p_value);
            }
            viz.record_pattern(p, aln.patterns[p].sites, type, r.p_value, metric, detail);
        }

        viz.finish();

        site_results.resize(aln.num_codons);
        for (size_t s = 0; s < aln.num_codons; ++s) {
            size_t p_idx = aln.site_to_pattern[s];
            site_results[s] = pattern_results[p_idx];
            site_results[s].site_index = s;
        }

        return site_results;
    }

nlohmann::json to_legacy_json() const;

nlohmann::json to_modern_json(const Provenance& prov = {}) const;
    nlohmann::json to_json(JSONFormat format = JSONFormat::ModernV3, const Provenance& prov = {}) const {
        if (format == JSONFormat::Legacy) {
            return to_legacy_json();
        }
        return to_modern_json(prov);
    }

    void save_json(const std::string& filepath, JSONFormat format = JSONFormat::ModernV3, const Provenance& prov = {}) const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Failed to open output file: " + filepath);
        }
        out << to_json(format, prov).dump(2);
    }
};

} // namespace hyphy::analyses
