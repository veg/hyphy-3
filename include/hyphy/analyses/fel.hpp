#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/core/provenance.hpp"
#include "LBFGSB.h"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "nlohmann/json.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <memory>

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

struct SiteResult {
    size_t site_index = 0;
    Scalar alpha = 0.0;
    Scalar beta = 0.0;
    Scalar alpha_null = 0.0;
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    Scalar total_branch_length = 0.0;
    Scalar log_l_alt = 0.0;
    Scalar log_l_null = 0.0;
};

class FELAnalyzer {
public:
    Tree tree;
    Alignment aln;
    MG94Parameters base_params;

    Scalar p_value_threshold = 0.1;
    Scalar baseline_tree_length = 0.0;
    Scalar conversion_factor = 1.0;
    Scalar global_log_l = 0.0;
    Scalar global_aicc = 0.0;
    std::vector<SiteResult> site_results;

    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;
    bool full_model = true;

    FELAnalyzer(Tree t, Alignment a, MG94Parameters params = MG94Parameters{})
        : tree(std::move(t)), aln(std::move(a)), base_params(params) {
        prepare_branch_lengths();
    }

    // Static factory that performs Phase 1 GTR and Phase 2 MG94 fitting automatically
    static FELAnalyzer create_and_fit(
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

        FELAnalyzer analyzer(mg_res.tree, std::move(aln), base_p);
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
        analyzer.full_model = full_model;
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

    // Evaluate single pattern log-likelihood with underflow protection
    static Scalar compute_pattern_log_likelihood(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const MG94Matrix& mg94_model,
        const std::vector<Matrix>& P_branches,
        std::vector<Vector>& node_L,
        std::vector<double>& node_log_scale
    ) {
        int S = mg94_model.pi.size();
        const double SCALER_THRESH = 1e-150;
        const double SCALER_UP = 1e150;
        const double LOG_SCALER = std::log(SCALER_UP);

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
                node_L[node_id].setOnes();
                double child_scale_sum = 0.0;
                for (int32_t child_id : node.children) {
                    node_L[node_id].array() *= (P_branches[child_id] * node_L[child_id]).array();
                    child_scale_sum += node_log_scale[child_id];
                }

                double max_val = node_L[node_id].maxCoeff();
                int boost_count = 0;
                while (max_val < SCALER_THRESH && max_val > 0.0) {
                    node_L[node_id] *= SCALER_UP;
                    max_val *= SCALER_UP;
                    boost_count++;
                }
                node_log_scale[node_id] = child_scale_sum - boost_count * LOG_SCALER;
            }
        }

        Scalar pat_l = mg94_model.pi.dot(node_L[tree.root_id]);
        if (pat_l > 0.0) {
            return std::log(pat_l) + node_log_scale[tree.root_id];
        } else {
            return -1e20;
        }
    }

    // Backward-compatible unscaled pattern likelihood
    static Scalar compute_pattern_likelihood(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const MG94Matrix& mg94_model,
        const std::vector<Matrix>& P_branches,
        std::vector<Vector>& node_L
    ) {
        std::vector<double> node_log_scale(tree.num_nodes(), 0.0);
        Scalar log_l = compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg94_model, P_branches, node_L, node_log_scale);
        return (log_l > -1e10) ? std::exp(log_l) : 0.0;
    }

    // Optimize single pattern
    SiteResult analyze_pattern(size_t p_idx, const std::vector<size_t>& leaf_to_taxon) const {
        const auto& pattern = aln.patterns[p_idx];
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;

        SiteResult res;
        res.site_index = p_idx;

        // Check if invariant (all non-gap states identical)
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
            res.beta = 0.0;
            res.alpha_null = 0.0;
            res.lrt = 0.0;
            res.p_value = 1.0;
            res.total_branch_length = 0.0;
            return res;
        }

        std::vector<Vector> thread_node_L(tree.num_nodes(), Vector::Zero(S));
        std::vector<double> thread_node_scale(tree.num_nodes(), 0.0);
        std::vector<Matrix> P_branches(tree.num_nodes());

        auto eval_pattern_lnl = [&](Scalar a, Scalar b) -> Scalar {
            MG94Parameters p = base_params;
            p.alpha = a;
            p.beta = b;

            MG94Matrix mg;
            mg.update(p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P_branches[node.id] = mg.transition_matrix(node.branch_length);
                }
            }

            return compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg, P_branches, thread_node_L, thread_node_scale);
        };

        // 1. Starting grid search (matching HyPhy FEL start grid)
        const std::vector<std::pair<Scalar, Scalar>> start_grid = {
            {0.01, 0.1}, {1.0, 0.1}, {1.0, 0.5}, {1.0, 1.0},
            {1.0, 5.0}, {10.0, 0.1}, {0.01, 0.5}, {0.01, 5.0},
            {10.0, 0.5}, {10.0, 1.0}, {10.0, 50.0}, {100.0, 1.0}
        };

        Scalar best_grid_lnl = -1e20;
        Scalar best_a = 1.0, best_b = 1.0;
        for (const auto& pt : start_grid) {
            Scalar lnl = eval_pattern_lnl(pt.first, pt.second);
            if (lnl > best_grid_lnl) {
                best_grid_lnl = lnl;
                best_a = pt.first;
                best_b = pt.second;
            }
        }

        // 2. Alternative Model Optimization: Nelder-Mead 2D
        auto obj_alt = [&](Scalar a, Scalar b) -> Scalar {
            return -eval_pattern_lnl(a, b);
        };
        auto pt_alt = NelderMead2D::minimize(obj_alt, best_a, best_b, 1e-4, 1000.0, 1e-4, 100);
        res.alpha = pt_alt.a;
        res.beta = pt_alt.b;
        res.log_l_alt = -pt_alt.f;

        // 3. Null Model Optimization: alpha = beta (Brent 1D)
        auto obj_null = [&](Scalar a) -> Scalar {
            return -eval_pattern_lnl(a, a);
        };

        Scalar init_null = (std::min(res.alpha, 100.0) + 3.0 * std::min(res.beta, 100.0)) / 4.0;
        Scalar best_null_lnl = eval_pattern_lnl(init_null, init_null);
        Scalar best_null_a = init_null;
        for (Scalar cand : {0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0}) {
            Scalar lnl = eval_pattern_lnl(cand, cand);
            if (lnl > best_null_lnl) {
                best_null_lnl = lnl;
                best_null_a = cand;
            }
        }

        Scalar ax = 1e-4;
        Scalar bx = std::clamp(best_null_a, 0.01, 45.0);
        Scalar cx = 50.0;
        auto [opt_null_a, opt_null_nll] = Brent1D::minimize(obj_null, ax, bx, cx, 1e-4, 60);

        res.alpha_null = opt_null_a;
        res.log_l_null = -opt_null_nll;

        // 4. Compute LRT and p-value
        Scalar lrt = 0.0;
        if (res.log_l_alt > -1e10 && res.log_l_null > -1e10) {
            lrt = 2.0 * (res.log_l_alt - res.log_l_null);
            if (lrt < 0.0) lrt = 0.0;
        }
        res.lrt = lrt;
        res.p_value = std::erfc(std::sqrt(lrt / 2.0));
        res.total_branch_length = res.alpha_null * baseline_tree_length;

        return res;
    }

    std::vector<SiteResult> run(bool show_progress = false, bool force_progress = false) {
        std::vector<size_t> leaf_to_taxon(tree.num_nodes(), static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                leaf_to_taxon[node.id] = LikelihoodEngine::find_taxon_index(aln, node.name);
            }
        }

        size_t num_patterns = aln.patterns.size();
        std::vector<SiteResult> pattern_results(num_patterns);

        std::unique_ptr<ProgressBar> pbar;
        std::atomic<size_t> num_positive{0};
        std::atomic<size_t> num_negative{0};

        if (show_progress && (force_progress || ProgressBar::is_terminal())) {
            pbar = std::make_unique<ProgressBar>(
                num_patterns,
                "[FEL] Patterns",
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
                if (r.p_value <= p_value_threshold) {
                    if (r.beta > r.alpha) {
                        num_positive.fetch_add(w, std::memory_order_relaxed);
                    } else if (r.alpha > r.beta) {
                        num_negative.fetch_add(w, std::memory_order_relaxed);
                    }
                }
                size_t pos = num_positive.load(std::memory_order_relaxed);
                size_t neg = num_negative.load(std::memory_order_relaxed);
                if (pos > 0 || neg > 0) {
                    std::string stat = "\033[1;32m+" + std::to_string(pos) + "\033[0m \033[1;31m-" + std::to_string(neg) + "\033[0m sites";
                    pbar->set_status(stat);
                }
                pbar->tick();
            }
        }

        if (pbar) {
            size_t pos = num_positive.load();
            size_t neg = num_negative.load();
            std::ostringstream summary;
            summary << "\033[1;32m" << pos << " positive\033[0m, \033[1;31m" << neg << " negative\033[0m sites (p ≤ "
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

    nlohmann::json to_legacy_json(const std::string& input_filepath = "", const std::string& tree_string = "") const {
        nlohmann::json j;

        // Analysis block
        j["analysis"]["authors"] = "Sergei L Kosakovsky Pond and Simon DW Frost";
        j["analysis"]["citation"] = "Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection (2005). _Mol Biol Evol_ 22 (5): 1208-1222";
        j["analysis"]["contact"] = "spond@temple.edu";
        j["analysis"]["info"] = "FEL (Fixed Effects Likelihood) estimates site-wise synonymous (alpha) and non-synonymous (beta) rates, and uses a likelihood ratio test to determine if beta != alpha at a site.";
        j["analysis"]["requirements"] = "in-frame codon alignment and a phylogenetic tree";
        j["analysis"]["settings"]["ci"] = 0;
        j["analysis"]["settings"]["multihit"] = "None";
        j["analysis"]["settings"]["pvalue"] = p_value_threshold;
        j["analysis"]["settings"]["resample"] = 0;
        j["analysis"]["settings"]["srv"] = 1;
        j["analysis"]["version"] = "3.0-next";

        // MLE block
        j["MLE"]["headers"] = nlohmann::json::array({
            nlohmann::json::array({"alpha", "Synonymous substitution rate at a site"}),
            nlohmann::json::array({"beta", "Non-synonymous substitution rate at a site"}),
            nlohmann::json::array({"alpha=beta", "The rate estimate under the neutral model"}),
            nlohmann::json::array({"LRT", "Likelihood ratio test statistic for beta = alpha, versus beta &neq; alpha"}),
            nlohmann::json::array({"p-value", "Asymptotic p-value for evidence of selection, i.e. beta &neq; alpha"}),
            nlohmann::json::array({"Total branch length", "The total length of branches contributing to inference at this site, and used to scale dN-dS"})
        });

        nlohmann::json content = nlohmann::json::array();
        for (const auto& r : site_results) {
            content.push_back({
                r.alpha,
                r.beta,
                r.alpha_null,
                r.lrt,
                r.p_value,
                r.total_branch_length
            });
        }
        j["MLE"]["content"]["0"] = content;

        // Input block
        j["input"]["file name"] = input_filepath;
        j["input"]["genetic code"] = aln.code ? aln.code->name : "Universal";
        j["input"]["number of sequences"] = aln.num_taxa;
        j["input"]["number of sites"] = aln.num_codons;
        j["input"]["partition count"] = 1;
        if (!tree_string.empty()) {
            j["input"]["trees"]["0"] = tree_string;
        }

        // Branch attributes
        nlohmann::json branch_attr;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar orig_bl = (conversion_factor > 0.0) ? (node.branch_length / conversion_factor) : node.branch_length;
                branch_attr[node.name]["Global MG94xREV"] = orig_bl;
                if (has_gtr_fit) {
                    branch_attr[node.name]["Nucleotide GTR"] = gtr_branch_lengths.count(node.name) ? gtr_branch_lengths.at(node.name) : orig_bl;
                }
                branch_attr[node.name]["original name"] = node.name;
            }
        }
        j["branch attributes"]["0"] = branch_attr;
        j["branch attributes"]["attributes"]["Global MG94xREV"]["attribute type"] = "branch length";
        j["branch attributes"]["attributes"]["Global MG94xREV"]["display order"] = 1;
        if (has_gtr_fit) {
            j["branch attributes"]["attributes"]["Nucleotide GTR"]["attribute type"] = "branch length";
            j["branch attributes"]["attributes"]["Nucleotide GTR"]["display order"] = 0;
        }

        // Fits
        if (has_gtr_fit) {
            j["fits"]["Nucleotide GTR"]["Log-likelihood"] = gtr_log_l;
            j["fits"]["Nucleotide GTR"]["AIC-c"] = gtr_aicc;
            j["fits"]["Nucleotide GTR"]["display order"] = 0;
            nlohmann::json gtr_ef = nlohmann::json::array();
            for (int i = 0; i < 4; ++i) {
                gtr_ef.push_back({aln.nuc_frequencies(i)});
            }
            j["fits"]["Nucleotide GTR"]["Equilibrium frequencies"] = gtr_ef;
            j["fits"]["Nucleotide GTR"]["Rate Distributions"] = {
                {"Substitution rate from nucleotide A to nucleotide C", gtr_rates.theta_AC},
                {"Substitution rate from nucleotide A to nucleotide G", 1.0},
                {"Substitution rate from nucleotide A to nucleotide T", gtr_rates.theta_AT},
                {"Substitution rate from nucleotide C to nucleotide G", gtr_rates.theta_CG},
                {"Substitution rate from nucleotide C to nucleotide T", gtr_rates.theta_CT},
                {"Substitution rate from nucleotide G to nucleotide T", gtr_rates.theta_GT}
            };
        }

        j["fits"]["Global MG94xREV"]["Log-likelihood"] = global_log_l;
        j["fits"]["Global MG94xREV"]["AIC-c"] = global_aicc;
        j["fits"]["Global MG94xREV"]["display order"] = 1;
        nlohmann::json mg_ef = nlohmann::json::array();
        for (int i = 0; i < aln.codon_frequencies_f3x4.size(); ++i) {
            mg_ef.push_back({aln.codon_frequencies_f3x4(i)});
        }
        j["fits"]["Global MG94xREV"]["Equilibrium frequencies"] = mg_ef;
        j["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"] = {
            {base_params.beta / base_params.alpha, 1.0}
        };

        return j;
    }

    nlohmann::json to_modern_json(const std::string& input_filepath = "",
                                  const std::string& tree_string = "",
                                  const Provenance& prov = {}) const {
        nlohmann::json j;

        j["$schema"] = "https://raw.githubusercontent.com/veg/hyphy-3/main/schemas/v3/fel.v3.schema.json";
        j["schema_version"] = "3.0.0";

        // Analysis block
        j["analysis"]["id"] = "fel";
        j["analysis"]["name"] = "Fixed Effects Likelihood";
        j["analysis"]["version"] = "3.0.0";
        j["analysis"]["category"] = "site_selection";
        j["analysis"]["description"] = "FEL (Fixed Effects Likelihood) estimates site-wise synonymous (alpha) and non-synonymous (beta) rates, and uses a likelihood ratio test to determine if beta != alpha at a site.";
        j["analysis"]["citations"] = nlohmann::json::array({
            {
                {"citation", "Kosakovsky Pond SL, Frost SDW (2005). Not So Different After All: A Comparison of Methods for Detecting Amino Acid Sites Under Selection. Mol Biol Evol 22(5): 1208-1222."},
                {"doi", "10.1093/molbev/msi105"},
                {"pmid", "15703277"}
            }
        });
        j["analysis"]["settings"]["pvalue_threshold"] = p_value_threshold;
        j["analysis"]["settings"]["srv"] = true;
        j["analysis"]["settings"]["ci"] = false;
        j["analysis"]["settings"]["full_model"] = full_model;

        // Software block
        j["software"]["name"] = "hyphy";
        j["software"]["version"] = "3.0.0";
        j["software"]["git_commit"] = "e1352b6";
        j["software"]["build_type"] = "Release";
        j["software"]["compiler"] = Provenance::detect_compiler();

        // Provenance block
        if (!prov.invocation.cli_command.empty() || !prov.inputs.empty()) {
            j["provenance"] = prov.to_json();
        } else {
            Provenance auto_prov;
            auto_prov.invocation.cli_command = "hyphy3 fel";
            auto_prov.invocation.working_directory = Provenance::get_cwd();
            if (!input_filepath.empty()) {
                size_t sz = 0;
                std::string hash = crypto::SHA256::hash_file(input_filepath, &sz);
                auto_prov.inputs["alignment"] = {input_filepath, "FASTA/NEXUS", hash, sz};
            }
            auto_prov.execution.start_time = Provenance::current_iso8601();
            auto_prov.execution.end_time = Provenance::current_iso8601();
            auto_prov.execution.wall_time_seconds = 0.0;
            auto_prov.execution.cpu_threads = 1;
            auto_prov.execution.hostname = Provenance::get_hostname();
            auto_prov.execution.os = Provenance::detect_os();
            auto_prov.execution.compiler = Provenance::detect_compiler();
            j["provenance"] = auto_prov.to_json();
        }

        // Dataset block
        j["dataset"]["taxa_count"] = aln.num_taxa;
        j["dataset"]["codon_sites"] = aln.num_codons;
        j["dataset"]["nucleotide_sites"] = aln.num_codons * 3;
        j["dataset"]["unique_patterns"] = aln.patterns.size();
        j["dataset"]["taxa"] = aln.taxon_names;

        nlohmann::json gcode;
        gcode["id"] = aln.code ? aln.code->name : "Universal";
        gcode["name"] = aln.code ? aln.code->name : "Universal";
        gcode["sense_codons"] = aln.code ? static_cast<int>(aln.code->sense_codons.size()) : 61;
        if (aln.code) {
            gcode["stop_codons"] = aln.code->stop_codons;
        } else {
            gcode["stop_codons"] = {"TAA", "TAG", "TGA"};
        }
        j["dataset"]["genetic_code"] = gcode;

        j["dataset"]["partitions"] = nlohmann::json::array({
            {
                {"id", "default"},
                {"name", "Full Alignment"},
                {"span", nlohmann::json::array({nlohmann::json::array({1, aln.num_codons})})},
                {"sites_count", aln.num_codons},
                {"patterns_count", aln.patterns.size()}
            }
        });

        // Phylogeny block
        std::string nwk = tree_string.empty() ? tree.to_newick() : tree_string;
        j["phylogeny"]["newick"] = nwk;

        nlohmann::json nodes_json;
        double sum_bl_mg94 = 0.0;
        double sum_bl_gtr = 0.0;

        for (const auto& node : tree.nodes) {
            nlohmann::json node_entry;
            node_entry["type"] = node.children.empty() ? "leaf" : "internal";

            Scalar mg_bl = (conversion_factor > 0.0) ? (node.branch_length / conversion_factor) : node.branch_length;
            node_entry["branch_lengths"]["codon_mg94"] = mg_bl;
            if (node.id != tree.root_id) sum_bl_mg94 += mg_bl;

            if (has_gtr_fit) {
                Scalar gtr_bl = gtr_branch_lengths.count(node.name) ? gtr_branch_lengths.at(node.name) : mg_bl;
                node_entry["branch_lengths"]["nucleotide_gtr"] = gtr_bl;
                if (node.id != tree.root_id) sum_bl_gtr += gtr_bl;
            }

            if (!node.model_tag.empty()) {
                node_entry["model_tags"] = nlohmann::json::array({node.model_tag});
            }

            nodes_json[node.name] = node_entry;
        }
        j["phylogeny"]["nodes"] = nodes_json;
        j["phylogeny"]["tree_lengths"]["codon_mg94"] = sum_bl_mg94;
        if (has_gtr_fit) {
            j["phylogeny"]["tree_lengths"]["nucleotide_gtr"] = sum_bl_gtr;
        }

        // Model fits block
        if (has_gtr_fit) {
            j["model_fits"]["nucleotide_gtr"]["description"] = "General Time Reversible nucleotide model with empirical frequencies";
            j["model_fits"]["nucleotide_gtr"]["log_likelihood"] = gtr_log_l;
            j["model_fits"]["nucleotide_gtr"]["aicc"] = gtr_aicc;
            j["model_fits"]["nucleotide_gtr"]["parameters_count"] = 5 + static_cast<int>(tree.nodes.size()) - 1;
            j["model_fits"]["nucleotide_gtr"]["equilibrium_frequencies"] = {
                {"A", aln.nuc_frequencies(0)},
                {"C", aln.nuc_frequencies(1)},
                {"G", aln.nuc_frequencies(2)},
                {"T", aln.nuc_frequencies(3)}
            };
            j["model_fits"]["nucleotide_gtr"]["substitution_rates"] = {
                {"AC", gtr_rates.theta_AC},
                {"AG", 1.0},
                {"AT", gtr_rates.theta_AT},
                {"CG", gtr_rates.theta_CG},
                {"CT", gtr_rates.theta_CT},
                {"GT", gtr_rates.theta_GT}
            };
        }

        j["model_fits"]["codon_mg94"]["description"] = "Muse-Gaut 1994 x GTR codon model with shared global omega";
        j["model_fits"]["codon_mg94"]["log_likelihood"] = global_log_l;
        j["model_fits"]["codon_mg94"]["aicc"] = global_aicc;
        j["model_fits"]["codon_mg94"]["parameters_count"] = 6 + static_cast<int>(tree.nodes.size()) - 1;
        j["model_fits"]["codon_mg94"]["global_parameters"] = {
            {"omega", (base_params.alpha > 0.0) ? (base_params.beta / base_params.alpha) : 1.0},
            {"tree_scale", 1.0}
        };

        // Statistical tests block
        auto summary = get_summary();
        j["statistical_tests"]["site_level_summary"] = {
            {"method", "Fixed Effects Likelihood (FEL)"},
            {"test_statistic", "Likelihood Ratio Test (LRT)"},
            {"distribution", "Asymptotic Chi-squared"},
            {"degrees_of_freedom", 1},
            {"null_hypothesis", "beta = alpha (Neutral evolution)"},
            {"alternative_hypothesis", "beta != alpha (Positive or negative selection)"},
            {"threshold", p_value_threshold},
            {"total_sites", static_cast<int>(summary.tested_sites)},
            {"diversifying_count", static_cast<int>(summary.positively_selected)},
            {"purifying_count", static_cast<int>(summary.negatively_selected)},
            {"neutral_count", static_cast<int>(summary.tested_sites - summary.positively_selected - summary.negatively_selected)}
        };

        // Site results (Columnar format)
        j["site_results"]["columns"] = {
            {"site", {{"type", "integer"}, {"unit", "1-based codon position"}, {"description", "Alignment codon site index (1-based)"}}},
            {"alpha", {{"type", "float"}, {"unit", "substitutions/site"}, {"description", "Synonymous substitution rate estimate (alpha)"}}},
            {"beta", {{"type", "float"}, {"unit", "substitutions/site"}, {"description", "Non-synonymous substitution rate estimate (beta)"}}},
            {"alpha_null", {{"type", "float"}, {"unit", "substitutions/site"}, {"description", "Constrained neutral rate (alpha = beta)"}}},
            {"lrt", {{"type", "float"}, {"description", "Likelihood ratio test statistic for beta != alpha"}}},
            {"p_value", {{"type", "float"}, {"description", "Asymptotic p-value from chi^2_1 distribution"}}},
            {"total_branch_length", {{"type", "float"}, {"unit", "substitutions/site"}, {"description", "Total length of branches contributing to inference at this site"}}},
            {"classification", {{"type", "string"}, {"description", "Selection classification at threshold (diversifying, purifying, neutral)"}}}
        };

        std::vector<int> col_site;
        std::vector<double> col_alpha, col_beta, col_alpha_null, col_lrt, col_pval, col_tbl;
        std::vector<std::string> col_class;

        col_site.reserve(site_results.size());
        col_alpha.reserve(site_results.size());
        col_beta.reserve(site_results.size());
        col_alpha_null.reserve(site_results.size());
        col_lrt.reserve(site_results.size());
        col_pval.reserve(site_results.size());
        col_tbl.reserve(site_results.size());
        col_class.reserve(site_results.size());

        for (size_t i = 0; i < site_results.size(); ++i) {
            const auto& r = site_results[i];
            col_site.push_back(static_cast<int>(i + 1));
            col_alpha.push_back(r.alpha);
            col_beta.push_back(r.beta);
            col_alpha_null.push_back(r.alpha_null);
            col_lrt.push_back(r.lrt);
            col_pval.push_back(r.p_value);
            col_tbl.push_back(r.total_branch_length);

            if (r.p_value <= p_value_threshold) {
                if (r.beta > r.alpha) {
                    col_class.push_back("diversifying");
                } else if (r.alpha > r.beta) {
                    col_class.push_back("purifying");
                } else {
                    col_class.push_back("neutral");
                }
            } else {
                col_class.push_back("neutral");
            }
        }

        j["site_results"]["data"]["site"] = col_site;
        j["site_results"]["data"]["alpha"] = col_alpha;
        j["site_results"]["data"]["beta"] = col_beta;
        j["site_results"]["data"]["alpha_null"] = col_alpha_null;
        j["site_results"]["data"]["lrt"] = col_lrt;
        j["site_results"]["data"]["p_value"] = col_pval;
        j["site_results"]["data"]["total_branch_length"] = col_tbl;
        j["site_results"]["data"]["classification"] = col_class;

        return j;
    }

    nlohmann::json to_json(const std::string& input_filepath = "",
                           const std::string& tree_string = "",
                           JSONFormat format = JSONFormat::ModernV3,
                           const Provenance& prov = {}) const {
        if (format == JSONFormat::Legacy) {
            return to_legacy_json(input_filepath, tree_string);
        }
        return to_modern_json(input_filepath, tree_string, prov);
    }

    void save_json(const std::string& output_filepath,
                   const std::string& input_filepath = "",
                   const std::string& tree_string = "",
                   JSONFormat format = JSONFormat::ModernV3,
                   const Provenance& prov = {}) const {
        auto j = to_json(input_filepath, tree_string, format, prov);
        std::ofstream out(output_filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open output JSON file: " + output_filepath);
        }
        out << j.dump(2) << std::endl;
    }

    struct Summary {
        size_t tested_sites = 0;
        size_t positively_selected = 0;
        size_t negatively_selected = 0;
    };

    Summary get_summary() const {
        Summary s;
        s.tested_sites = site_results.size();
        for (const auto& r : site_results) {
            if (r.p_value <= p_value_threshold) {
                if (r.beta > r.alpha) {
                    s.positively_selected++;
                } else if (r.alpha > r.beta) {
                    s.negatively_selected++;
                }
            }
        }
        return s;
    }
};

} // namespace hyphy::analyses
