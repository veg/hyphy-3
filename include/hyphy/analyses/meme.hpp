#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
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


    MEMEAnalyzer(Tree t, Alignment a, MG94Parameters params = MG94Parameters{})
        : tree(t), aln(std::move(a)), base_params(params), fel_analyzer(t, aln, params) {
        prepare_branch_lengths();
    }

    static MEMEAnalyzer create_and_fit(
        Tree input_tree,
        Alignment aln,
        Scalar pvalue_threshold = 0.1,
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool full_model = true
    ) {
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

        MEMEAnalyzer analyzer(mg_res.tree, std::move(aln), base_p);
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

        std::vector<Vector> thread_node_L(num_nodes, Vector::Zero(S));
        std::vector<double> thread_node_scale(num_nodes, 0.0);
        std::vector<Matrix> P_branches(num_nodes);

        auto eval_pattern_meme_lnl = [&](Scalar a, Scalar w1, Scalar p1, Scalar b_plus) -> Scalar {
            Scalar b1 = a * w1;
            MG94Parameters p_1 = base_params;
            p_1.alpha = a;
            p_1.beta = b1;

            MG94Parameters p_p = base_params;
            p_p.alpha = a;
            p_p.beta = b_plus;

            MG94Matrix mg1, mg_plus;
            mg1.update(p_1, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            mg_plus.update(p_p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

            Scalar p_plus_wt = 1.0 - p1;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P_branches[node.id] = p1 * mg1.transition_matrix(node.branch_length) +
                                          p_plus_wt * mg_plus.transition_matrix(node.branch_length);
                }
            }

            return FELAnalyzer::compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg1, P_branches, thread_node_L, thread_node_scale);
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
            if (log_l_null > res.log_l_meme) {
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

            MG94Matrix mg1, mg_plus;
            mg1.update(p_1, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            mg_plus.update(p_p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

            std::vector<Matrix> P1(num_nodes), P_plus(num_nodes), P_mix(num_nodes);
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P1[node.id] = mg1.transition_matrix(node.branch_length);
                    P_plus[node.id] = mg_plus.transition_matrix(node.branch_length);
                    P_mix[node.id] = res.p1 * P1[node.id] + res.p_plus * P_plus[node.id];
                }
            }

            int branches_under_sel = 0;
            std::vector<Matrix> P_temp = P_mix;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    // Conditional on branch node.id being in class 1
                    P_temp[node.id] = P1[node.id];
                    Scalar log_L1 = FELAnalyzer::compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg1, P_temp, thread_node_L, thread_node_scale);

                    // Conditional on branch node.id being in class +
                    P_temp[node.id] = P_plus[node.id];
                    Scalar log_L_plus = FELAnalyzer::compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg1, P_temp, thread_node_L, thread_node_scale);

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

        std::unique_ptr<ProgressBar> pbar;
        std::atomic<size_t> num_episodic{0};

        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar = std::make_unique<ProgressBar>(
                num_patterns,
                "[MEME] Patterns",
                "patterns",
                force_progress,
                ProgressBar::Style::SmoothBlocks
            );
        }

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            pattern_results[p] = analyze_pattern(p, leaf_to_taxon);
            if (pbar) {
                const auto& r = pattern_results[p];
                size_t w = aln.patterns[p].weight;
                if (r.p_value <= p_value_threshold && r.beta_plus > r.alpha) {
                    num_episodic.fetch_add(w, std::memory_order_relaxed);
                }
                size_t epi = num_episodic.load(std::memory_order_relaxed);
                if (epi > 0) {
                    std::string stat = "\033[1;32m+" + std::to_string(epi) + "\033[0m episodic sites";
                    pbar->set_status(stat);
                }
                pbar->tick();
            }
        }

        if (pbar) {
            size_t epi = num_episodic.load();
            std::ostringstream summary;
            summary << "\033[1;32m" << epi << " site(s) under episodic diversifying selection\033[0m (p ≤ "
                    << std::fixed << std::setprecision(2) << p_value_threshold << ")";
            pbar->finish(summary.str());
        }

        site_results.resize(aln.num_codons);
        for (size_t s = 0; s < aln.num_codons; ++s) {
            size_t p_idx = aln.site_to_pattern[s];
            site_results[s] = pattern_results[p_idx];
            site_results[s].site_index = s;
        }

        return site_results;
    }

    nlohmann::json to_json() const {
        nlohmann::json j;

        // Analysis metadata
        j["analysis"] = {
            {"authors", "Sergei L. Kosakovsky Pond, Steven Weaver"},
            {"citation", "Detecting Individual Sites Subject to Episodic Diversifying Selection. PLoS Genet 8(7): e1002764."},
            {"contact", "spond@temple.edu"},
            {"info", "MEME (Mixed Effects Model of Evolution) estimates a site-wise synonymous (alpha) and a two-category mixture of non-synonymous (beta1 with p1, and beta+ with p+) rates, and uses a likelihood ratio test to determine if beta+ > alpha at a site."},
            {"requirements", "in-frame codon alignment and a phylogenetic tree"},
            {"settings", {{"p-value", p_value_threshold}, {"code", aln.code ? aln.code->name : "Universal"}}},
            {"version", "4.1 (Modern C++ hyphy-next)"}
        };

        // Input
        j["input"]["file name"] = input_filepath;
        j["input"]["genetic code"] = aln.code ? aln.code->name : "Universal";
        j["input"]["number of sequences"] = aln.num_taxa;
        j["input"]["number of sites"] = aln.num_codons;
        j["input"]["partition count"] = 1;
        if (!tree_string.empty()) {
            j["input"]["trees"]["0"] = tree_string;
        }


        // Fits
        nlohmann::json fits_json;
        if (has_gtr_fit) {
            fits_json["Nucleotide GTR"] = {
                {"AIC-c", gtr_aicc},
                {"Log Likelihood", gtr_log_l},
                {"Equilibrium frequencies", {
                    aln.nuc_frequencies(0),
                    aln.nuc_frequencies(1),
                    aln.nuc_frequencies(2),
                    aln.nuc_frequencies(3)
                }},
                {"Rate Distributions", {
                    {"Substitution rate from nucleotide A to nucleotide C", gtr_rates.theta_AC},
                    {"Substitution rate from nucleotide A to nucleotide G", 1.0},
                    {"Substitution rate from nucleotide A to nucleotide T", gtr_rates.theta_AT},
                    {"Substitution rate from nucleotide C to nucleotide G", gtr_rates.theta_CG},
                    {"Substitution rate from nucleotide C to nucleotide T", gtr_rates.theta_CT},
                    {"Substitution rate from nucleotide G to nucleotide T", gtr_rates.theta_GT}
                }}
            };
        }

        std::vector<Scalar> eq_freqs(aln.codon_frequencies_f3x4.data(), aln.codon_frequencies_f3x4.data() + aln.codon_frequencies_f3x4.size());
        fits_json["Global MG94xREV"] = {
            {"AIC-c", global_aicc},
            {"Log Likelihood", global_log_l},
            {"Equilibrium frequencies", eq_freqs},
            {"Rate Distributions", {
                {"non-synonymous/synonymous rate ratio", base_params.beta / std::max(base_params.alpha, 1e-12)}
            }}
        };
        j["fits"] = fits_json;

        // Branch Attributes
        nlohmann::json branch_attrs;
        nlohmann::json branch_dict;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                nlohmann::json b_info;
                b_info["Global MG94xREV"] = node.branch_length / conversion_factor;
                auto it = gtr_branch_lengths.find(node.name);
                if (it != gtr_branch_lengths.end()) {
                    b_info["Nucleotide GTR"] = it->second;
                }
                branch_dict[node.name] = b_info;
            }
        }
        branch_attrs["0"] = branch_dict;
        j["branch attributes"] = branch_attrs;

        // MLE headers and content
        j["MLE"]["headers"] = {
            {"&alpha;", "Synonymous substitution rate at a site"},
            {"&beta;<sup>1</sup>", "Non-synonymous substitution rate at a site for the negative/neutral evolution component 1"},
            {"p<sup>1</sup>", "Mixture distribution weight allocated to negative/neutral evolution component 1"},
            {"&beta;<sup>+</sup>", "Non-synonymous substitution rate at a site for the positive selection component"},
            {"p<sup>+</sup>", "Mixture distribution weight allocated to the positive selection component"},
            {"LRT", "Likelihood ratio test statistic for episodic diversification, i.e., p<sup>+</sup> > 0 and &beta;<sup>+</sup> > &alpha;"},
            {"p-value", "Asymptotic p-value for episodic diversification, i.e., p<sup>+</sup> > 0 and &beta;<sup>+</sup> > &alpha;"},
            {"# branches under selection", "The count of branches with empirical Bayes factor >= 100 for the &beta;<sup>+</sup> rate"},
            {"Total branch length", "The total length of branches contributing to inference at this site"},
            {"MEME LogL", "Site Log-likelihood under the MEME model"},
            {"FEL LogL", "Site Log-likelihood under the FEL model"},
            {"LRT MEME vs FEL", "Likelihood ratio test statistic p-value for MEME vs FEL"},
            {"FEL &alpha;", "Synonymous substitution rate at a site under the FEL model"},
            {"FEL &beta;", "Non-synonymous substitution rate at a site under the FEL model"}
        };

        std::vector<std::vector<Scalar>> content(site_results.size());
        for (size_t s = 0; s < site_results.size(); ++s) {
            const auto& r = site_results[s];
            content[s] = {
                r.alpha,
                r.beta1,
                r.p1,
                r.beta_plus,
                r.p_plus,
                r.lrt,
                r.p_value,
                static_cast<Scalar>(r.branches_under_selection),
                r.total_branch_length,
                r.log_l_meme,
                r.log_l_fel,
                r.p_value_meme_vs_fel,
                r.fel_alpha,
                r.fel_beta
            };
        }
        j["MLE"]["content"]["0"] = content;

        // Tested branches
        nlohmann::json tested_dict;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                tested_dict[node.name] = "test";
            }
        }
        j["tested"]["0"] = tested_dict;

        return j;
    }

    void save_json(const std::string& filepath) const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Failed to open output file: " + filepath);
        }
        out << to_json().dump(2);
    }
};

} // namespace hyphy::analyses
