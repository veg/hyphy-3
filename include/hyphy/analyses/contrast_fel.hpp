#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/core/ancestral.hpp"
#include "hyphy/analyses/site_visualizer.hpp"
#include "hyphy/opt/nelder_mead.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <regex>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

// NelderMeadDynamic is defined in hyphy::opt::NelderMeadDynamic

struct ContrastBranchGroup {
    int id = 0;                     // 0 .. K-1 (or -1 for background)
    std::string name;              // e.g. "Test", "Reference", "Group1"
    bool is_testable = true;
    size_t branch_count = 0;
    Scalar total_branch_length = 0.0;
};

struct PairwiseContrastResult {
    size_t group_a = 0;
    size_t group_b = 1;
    std::string name_a;
    std::string name_b;
    Scalar beta_a = 0.0;
    Scalar beta_b = 0.0;
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    Scalar p_value_adjusted = 1.0;  // Holm-Bonferroni step-down
    std::string direction;          // e.g. "A < B", "A > B", "A = B"
};

struct ContrastSubEvent {
    int32_t node_id = -1;
    std::string branch_name;
    std::string group_name;
    int group_id = -1;
    Scalar branch_length = 0.0;

    int8_t parent_codon_idx = -1;
    int8_t child_codon_idx = -1;
    std::string parent_codon;
    std::string child_codon;

    char parent_aa = '?';
    char child_aa = '?';

    int nt_distance = 0;
    bool is_synonymous = false;

    Scalar transition_prob = 0.0;
    Scalar surprise_bits = 0.0;
};

struct ContrastSiteSubstitutions {
    size_t site_index = 0;
    std::string root_codon;
    char root_aa = '?';
    std::vector<ContrastSubEvent> events;

    std::vector<size_t> syn_counts;      // per group (size K, +1 if background)
    std::vector<size_t> nonsyn_counts;   // per group (size K, +1 if background)

    size_t subs_1nt = 0;
    size_t subs_2nt = 0;
    size_t subs_3nt = 0;

    Scalar mean_surprise_bits = 0.0;
    Scalar total_surprise_bits = 0.0;
    Scalar mean_branch_length = 0.0;

    std::string summary_string = "Conserved";
    std::string summary_string_partitioned = "Conserved";
};

struct ContrastSiteResult {
    size_t site_index = 0;
    Scalar alpha = 0.0;
    std::vector<Scalar> beta;           // size K (testable groups)
    Scalar beta_background = 0.0;       // nuisance parameter if background present
    Scalar alpha_null = 0.0;
    Scalar beta_null = 0.0;             // shared non-synonymous rate under H0

    Scalar log_l_alt = 0.0;
    Scalar log_l_null = 0.0;
    Scalar penalized_log_l_alt = 0.0;
    Scalar raw_lrt = 0.0;
    Scalar penalized_lrt = 0.0;
    Scalar shrinkage_penalty = 0.0;

    Scalar p_value = 1.0;               // Asymptotic p-value from chi^2_{K-1}
    Scalar q_value = 1.0;               // Alignment-wide Benjamini-Hochberg FDR
    Scalar total_tested_branch_length = 0.0;

    std::vector<PairwiseContrastResult> pairwise;
    ContrastSiteSubstitutions subs;

    std::string contrast_tag;           // e.g. "G1 ≪ G2 (Adaptive in G2)"
    std::string classification;         // "different" (p <= threshold) or "neutral/equal"
};

struct ContrastFELSettings {
    Scalar p_value_threshold = 0.05;
    Scalar q_value_threshold = 0.20;
    Scalar shrinkage_lambda = 1.0;      // Empirical Bayes Ridge shrinkage parameter
    bool srv = true;                    // Synonymous rate variation (alpha varies per site)
    bool full_model = true;             // Full MG94 branch length re-optimization
    std::string branch_specs = "";      // Predefined branch sets or auto-detect
    std::string branch_regex = "";
    std::vector<std::string> explicit_groups = {};
};

class ContrastFELAnalyzer {
public:
    Tree tree;
    Alignment aln;
    MG94Parameters base_params;
    ContrastFELSettings settings;

    Scalar p_value_threshold = 0.05;
    Scalar q_value_threshold = 0.20;
    Scalar shrinkage_lambda = 1.0;
    Scalar baseline_tree_length = 0.0;
    Scalar conversion_factor = 1.0;
    Scalar global_log_l = 0.0;
    Scalar global_aicc = 0.0;
    std::vector<ContrastSiteResult> site_results;

    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    // Branch group definitions
    std::vector<ContrastBranchGroup> testable_groups;
    ContrastBranchGroup background_group;
    bool has_background = false;
    std::vector<int> node_group_id;     // node.id -> 0..K-1 (testable) or -1 (background)
    std::vector<Scalar> group_omegas;
    Scalar background_omega = -1.0;

    MG94Matrix base_mg94_matrix;
    std::vector<Matrix> branch_transition_matrices;

    ContrastFELAnalyzer(
        Tree t,
        Alignment a,
        MG94Parameters params = MG94Parameters{},
        ContrastFELSettings s = ContrastFELSettings{},
        std::vector<Scalar> grp_om = {},
        Scalar bg_om = -1.0
    )
        : tree(std::move(t)), aln(std::move(a)), base_params(params), settings(std::move(s)),
          group_omegas(std::move(grp_om)), background_omega(bg_om) {
        p_value_threshold = settings.p_value_threshold;
        q_value_threshold = settings.q_value_threshold;
        shrinkage_lambda = settings.shrinkage_lambda;
        prepare_branch_partitions();
        prepare_branch_lengths();
    }

    void prepare_branch_partitions() {
        size_t n = tree.num_nodes();
        node_group_id.assign(n, -1);
        testable_groups.clear();
        has_background = false;

        // 1. Gather all tags in tree
        std::map<std::string, std::vector<int32_t>> tag_to_nodes;
        std::vector<int32_t> untagged_nodes;

        for (const auto& node : tree.nodes) {
            if (node.id == tree.root_id) continue;
            if (!node.model_tag.empty()) {
                tag_to_nodes[node.model_tag].push_back(node.id);
            } else {
                untagged_nodes.push_back(node.id);
            }
        }

        auto to_lower = [](std::string s) {
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
            return s;
        };
        auto trim = [](std::string s) {
            while (!s.empty() && (std::isspace(static_cast<unsigned char>(s.front())) || s.front() == '{')) s.erase(s.begin());
            while (!s.empty() && (std::isspace(static_cast<unsigned char>(s.back())) || s.back() == '}')) s.pop_back();
            return s;
        };

        std::string mode_lower = to_lower(trim(settings.branch_specs));

        if (mode_lower == "internal,leaves" || mode_lower == "internal,terminal" || mode_lower == "internal_leaves") {
            // Partition Internal vs Leaves
            ContrastBranchGroup g_int;
            g_int.id = 0;
            g_int.name = "Internal";
            g_int.is_testable = true;

            ContrastBranchGroup g_leaf;
            g_leaf.id = 1;
            g_leaf.name = "Leaves";
            g_leaf.is_testable = true;

            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                if (!node.is_leaf) {
                    node_group_id[node.id] = 0;
                    g_int.branch_count++;
                } else {
                    node_group_id[node.id] = 1;
                    g_leaf.branch_count++;
                }
            }
            testable_groups.push_back(g_int);
            testable_groups.push_back(g_leaf);
        } else if (!mode_lower.empty() && mode_lower != "auto" && mode_lower != "all") {
            // User explicitly specified branch groups via --branches / --branch-sets
            std::vector<std::string> requested_specs;
            std::stringstream ss(settings.branch_specs);
            std::string item;
            while (std::getline(ss, item, ',')) {
                std::string cleaned = trim(item);
                if (!cleaned.empty()) {
                    requested_specs.push_back(cleaned);
                }
            }

            std::vector<bool> assigned(n, false);
            int grp_counter = 0;

            for (const auto& req : requested_specs) {
                std::string req_lower = to_lower(req);
                std::string canonical_name = req;
                std::vector<int32_t> matched_nodes;

                // Check tag_to_nodes (exact or case-insensitive)
                bool found_tag = false;
                for (const auto& [tag, nodes] : tag_to_nodes) {
                    if (to_lower(tag) == req_lower) {
                        canonical_name = tag;
                        matched_nodes = nodes;
                        found_tag = true;
                        break;
                    }
                }

                if (!found_tag) {
                    if (req_lower == "internal") {
                        canonical_name = "Internal";
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id && !node.is_leaf && !assigned[node.id]) {
                                matched_nodes.push_back(node.id);
                            }
                        }
                    } else if (req_lower == "leaves" || req_lower == "terminal" || req_lower == "tips") {
                        canonical_name = "Leaves";
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id && node.is_leaf && !assigned[node.id]) {
                                matched_nodes.push_back(node.id);
                            }
                        }
                    } else if (req_lower == "unlabeled" || req_lower == "untagged" || req_lower == "reference" || req_lower == "background") {
                        canonical_name = req;
                        for (int32_t nid : untagged_nodes) {
                            if (!assigned[nid]) matched_nodes.push_back(nid);
                        }
                    } else {
                        // Try matching branch name directly
                        for (const auto& node : tree.nodes) {
                            if (node.id != tree.root_id && to_lower(node.name) == req_lower && !assigned[node.id]) {
                                matched_nodes.push_back(node.id);
                            }
                        }
                    }
                }

                if (matched_nodes.empty()) {
                    throw std::runtime_error("Contrast-FEL: Branch group '" + req + "' matched 0 branches in the tree.");
                }

                ContrastBranchGroup grp;
                grp.id = grp_counter++;
                grp.name = canonical_name;
                grp.is_testable = true;
                grp.branch_count = matched_nodes.size();
                for (int32_t nid : matched_nodes) {
                    node_group_id[nid] = grp.id;
                    assigned[nid] = true;
                }
                testable_groups.push_back(grp);
            }

            // Collect all remaining unassigned non-root branches
            std::vector<int32_t> remaining_nodes;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && !assigned[node.id]) {
                    remaining_nodes.push_back(node.id);
                }
            }

            if (testable_groups.size() == 1) {
                if (remaining_nodes.empty()) {
                    throw std::runtime_error("Contrast-FEL requires at least 2 branch groups to compare (all branches were assigned to '" + testable_groups[0].name + "').");
                }
                ContrastBranchGroup ref_grp;
                ref_grp.id = 1;
                ref_grp.name = "Reference";
                ref_grp.is_testable = true;
                ref_grp.branch_count = remaining_nodes.size();
                for (int32_t nid : remaining_nodes) {
                    node_group_id[nid] = 1;
                }
                testable_groups.push_back(ref_grp);
            } else if (!remaining_nodes.empty()) {
                has_background = true;
                background_group.id = -1;
                background_group.name = "Background";
                background_group.is_testable = false;
                background_group.branch_count = remaining_nodes.size();
                for (int32_t nid : remaining_nodes) {
                    node_group_id[nid] = -1;
                }
            }
        } else if (!tag_to_nodes.empty()) {
            // If tree has tags, use them!
            int grp_counter = 0;
            for (const auto& [tag, nodes] : tag_to_nodes) {
                ContrastBranchGroup grp;
                grp.id = grp_counter++;
                grp.name = tag;
                grp.is_testable = true;
                grp.branch_count = nodes.size();
                for (int32_t nid : nodes) {
                    node_group_id[nid] = grp.id;
                }
                testable_groups.push_back(grp);
            }

            if (!untagged_nodes.empty()) {
                if (testable_groups.size() == 1) {
                    // Only 1 tag exists in tree -> untagged becomes Group 1 (Reference/Background)
                    ContrastBranchGroup ref_grp;
                    ref_grp.id = 1;
                    ref_grp.name = "Reference";
                    ref_grp.is_testable = true;
                    ref_grp.branch_count = untagged_nodes.size();
                    for (int32_t nid : untagged_nodes) {
                        node_group_id[nid] = 1;
                    }
                    testable_groups.push_back(ref_grp);
                } else {
                    // Multiple tags exist -> untagged is Background nuisance
                    has_background = true;
                    background_group.id = -1;
                    background_group.name = "Background";
                    background_group.is_testable = false;
                    background_group.branch_count = untagged_nodes.size();
                    for (int32_t nid : untagged_nodes) {
                        node_group_id[nid] = -1;
                    }
                }
            }
        } else {
            // Fallback: Partition by Internal vs Leaves if no tags found
            ContrastBranchGroup g_int;
            g_int.id = 0;
            g_int.name = "Internal";
            g_int.is_testable = true;

            ContrastBranchGroup g_leaf;
            g_leaf.id = 1;
            g_leaf.name = "Leaves";
            g_leaf.is_testable = true;

            for (const auto& node : tree.nodes) {
                if (node.id == tree.root_id) continue;
                if (!node.is_leaf) {
                    node_group_id[node.id] = 0;
                    g_int.branch_count++;
                } else {
                    node_group_id[node.id] = 1;
                    g_leaf.branch_count++;
                }
            }
            testable_groups.push_back(g_int);
            testable_groups.push_back(g_leaf);
        }
    }

    void prepare_branch_lengths() {
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        base_mg94_matrix.update(base_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
        Scalar Q_scale = base_mg94_matrix.scale_factor;
        conversion_factor = (Q_scale > 1e-12) ? (3.0 / Q_scale) : 1.0;

        std::vector<Scalar> group_conv(testable_groups.size(), conversion_factor);
        for (size_t g = 0; g < testable_groups.size(); ++g) {
            if (g < group_omegas.size() && group_omegas[g] > 0.0) {
                MG94Parameters pg = base_params;
                pg.alpha = 1.0;
                pg.beta = group_omegas[g];
                MG94Matrix mgg;
                mgg.update(pg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
                if (mgg.scale_factor > 1e-12) {
                    group_conv[g] = 3.0 / mgg.scale_factor;
                }
            }
        }
        Scalar bg_conv = conversion_factor;
        if (has_background && background_omega > 0.0) {
            MG94Parameters pbg = base_params;
            pbg.alpha = 1.0;
            pbg.beta = background_omega;
            MG94Matrix mgbg;
            mgbg.update(pbg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            if (mgbg.scale_factor > 1e-12) {
                bg_conv = 3.0 / mgbg.scale_factor;
            }
        }

        baseline_tree_length = 0.0;
        for (auto& grp : testable_groups) {
            grp.total_branch_length = 0.0;
        }
        background_group.total_branch_length = 0.0;

        for (auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                baseline_tree_length += node.branch_length;
                int gid = node_group_id[node.id];
                Scalar node_conv = conversion_factor;
                if (gid >= 0 && gid < static_cast<int>(testable_groups.size())) {
                    testable_groups[gid].total_branch_length += node.branch_length;
                    node_conv = group_conv[gid];
                } else if (has_background) {
                    background_group.total_branch_length += node.branch_length;
                    node_conv = bg_conv;
                }
                node.branch_length *= node_conv;
            }
        }
        tree.scale_bridge_lengths(conversion_factor);

        branch_transition_matrices.resize(tree.num_nodes());
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                branch_transition_matrices[node.id] = base_mg94_matrix.transition_matrix(node.branch_length);
            }
        }
    }

    // Static factory that performs Phase 1 GTR and Phase 2 MG94 fitting automatically
    static ContrastFELAnalyzer create_and_fit(
        Tree input_tree,
        Alignment aln,
        ContrastFELSettings settings = ContrastFELSettings{},
        std::function<void(const std::string&, double)> progress_cb = nullptr
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

        // Partition branches beforehand so Phase 2 MG94 fits partitioned group omegas
        ContrastFELAnalyzer part_probe(gtr_res.tree, aln, base_p, settings);
        std::vector<bool> is_group0(gtr_res.tree.num_nodes(), false);
        bool has_group0 = false, has_other = false;
        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id && node.parent_id != INVALID_INDEX) {
                if (part_probe.node_group_id[node.id] == 0) {
                    is_group0[node.id] = true;
                    has_group0 = true;
                } else {
                    has_other = true;
                }
            }
        }
        const std::vector<bool>* p_part = (has_group0 && has_other) ? &is_group0 : nullptr;

        MG94Fitter mg_fitter(gtr_res.tree, aln);
        FitResult mg_res;
        Scalar omega_other = -1.0;
        if (settings.full_model) {
            mg_res = mg_fitter.fit_full_model(1.0, base_p, progress_cb, 250, 1e-5, p_part, p_part ? &omega_other : nullptr);
            base_p = mg_res.params;
        } else {
            mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
            base_p.alpha = 1.0;
            base_p.beta = mg_res.x_opt(0);
            omega_other = base_p.beta;
        }

        std::vector<Scalar> grp_omegas(part_probe.testable_groups.size(), base_p.beta);
        if (!grp_omegas.empty()) {
            grp_omegas[0] = base_p.beta;
            for (size_t g = 1; g < grp_omegas.size(); ++g) {
                grp_omegas[g] = (omega_other > 0.0) ? omega_other : base_p.beta;
            }
        }
        Scalar bg_om = (omega_other > 0.0) ? omega_other : base_p.beta;

        ContrastFELAnalyzer analyzer(mg_res.tree, std::move(aln), base_p, settings, grp_omegas, bg_om);
        analyzer.global_log_l = mg_res.log_likelihood;
        size_t K_mg = 5 + (p_part ? 2 : 1) + (mg_res.tree.num_nodes() - 1);
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

    // Pattern log-likelihood evaluator with scaling for underflow protection
    static Scalar compute_pattern_log_likelihood(
        const Tree& tree,
        const SitePattern& pattern,
        const std::vector<size_t>& leaf_to_taxon,
        const MG94Matrix& mg94_model,
        const std::vector<Matrix>& P_branches,
        std::vector<Vector>& node_L,
        std::vector<double>& node_log_scale,
        const std::vector<std::array<Matrix, 3>>* P_bridges = nullptr
    ) {
        int S = mg94_model.pi.size();

        for (int32_t node_id : tree.post_order) {
            const auto& node = tree.nodes[node_id];

            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[node_id];
                LikelihoodEngine::init_leaf_likelihood(pattern, t_idx, S, node_L[node_id]);
                node_log_scale[node_id] = 0.0;
            } else {
                double child_scale_sum = 0.0;
                tree.eval_internal_node_with_model(node_id, P_branches, mg94_model, node_L, node_log_scale, node_L[node_id], child_scale_sum, 1.0, nullptr, P_bridges);
                node_log_scale[node_id] = child_scale_sum;
            }
        }

        Scalar pat_l = mg94_model.pi.dot(node_L[tree.root_id]);
        if (pat_l > 0.0) {
            return std::log(pat_l) + node_log_scale[tree.root_id];
        } else {
            return -1e20;
        }
    }

    // Single site optimization with Empirical Bayes Ridge Shrinkage
    ContrastSiteResult analyze_pattern(size_t p_idx, const std::vector<size_t>& leaf_to_taxon) const {
        const auto& pattern = aln.patterns[p_idx];
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t K = testable_groups.size();

        ContrastSiteResult res;
        res.site_index = p_idx;
        res.beta.assign(K, 0.0);

        // Check if invariant across observed taxa
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
            res.beta.assign(K, 0.0);
            res.beta_background = 0.0;
            res.alpha_null = 0.0;
            res.beta_null = 0.0;
            res.penalized_lrt = 0.0;
            res.raw_lrt = 0.0;
            res.p_value = 1.0;
            res.q_value = 1.0;
            res.contrast_tag = "Invariant";
            res.classification = "equal";
            return res;
        }

        std::vector<Vector> thread_node_L(tree.num_nodes(), Vector::Zero(S));
        std::vector<double> thread_node_scale(tree.num_nodes(), 0.0);
        std::vector<Matrix> P_branches(tree.num_nodes());
        std::vector<std::array<Matrix, 3>> P_bridges(tree.bridges.size());

        // Pre-allocate matrices per group
        std::vector<MG94Matrix> mg_groups(K);
        MG94Matrix mg_bg;

        auto eval_pattern_lnl = [&](Scalar a, const std::vector<Scalar>& b_rates, Scalar b_bg) -> Scalar {
            for (size_t k = 0; k < K; ++k) {
                MG94Parameters p_k = base_params;
                p_k.alpha = a;
                p_k.beta = b_rates[k];
                mg_groups[k].update(p_k, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            }

            if (has_background) {
                MG94Parameters p_bg = base_params;
                p_bg.alpha = a;
                p_bg.beta = b_bg;
                mg_bg.update(p_bg, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);
            }

            auto model_for_node = [&](int32_t nid) -> const MG94Matrix& {
                int gid = (nid >= 0 && static_cast<size_t>(nid) < node_group_id.size()) ? node_group_id[nid] : -1;
                if (gid >= 0 && gid < static_cast<int>(K)) {
                    return mg_groups[gid];
                } else if (has_background) {
                    return mg_bg;
                } else {
                    return mg_groups[0];
                }
            };

            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P_branches[node.id] = model_for_node(node.id).transition_matrix(node.branch_length);
                }
            }

            const std::vector<std::array<Matrix, 3>>* p_bridges_ptr = nullptr;
            if (!tree.bridges.empty()) {
                for (size_t k = 0; k < tree.bridges.size(); ++k) {
                    const auto& wb = tree.bridges[k];
                    int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
                    const MG94Matrix& mg_br = model_for_node(bid);
                    for (int r = 0; r < 3; ++r) {
                        if (wb.weights[r] > 0.0 || r == 0) {
                            P_bridges[k][r] = mg_br.transition_matrix(wb.branch_lengths[r]);
                        }
                    }
                }
                p_bridges_ptr = &P_bridges;
            }

            return compute_pattern_log_likelihood(tree, pattern, leaf_to_taxon, mg_groups[0], P_branches, thread_node_L, thread_node_scale, p_bridges_ptr);
        };

        // Ridge Shrinkage Penalty Calculation
        auto calc_shrinkage_penalty = [&](const std::vector<Scalar>& b_rates) -> Scalar {
            if (shrinkage_lambda <= 0.0 || K < 2) return 0.0;
            Scalar sum_log_b = 0.0;
            for (size_t k = 0; k < K; ++k) {
                sum_log_b += std::log(std::max(b_rates[k], 1e-5));
            }
            Scalar mean_log_b = sum_log_b / static_cast<Scalar>(K);
            Scalar sq_dev = 0.0;
            for (size_t k = 0; k < K; ++k) {
                Scalar diff = std::log(std::max(b_rates[k], 1e-5)) - mean_log_b;
                sq_dev += diff * diff;
            }
            return 0.5 * shrinkage_lambda * sq_dev;
        };

        // 1. Fit Alternative Model (Nelder-Mead ND)
        // Parameter vector x: [alpha, beta_0, beta_1, ..., beta_{K-1}, (beta_bg)]
        size_t alt_dim = 1 + K + (has_background ? 1 : 0);
        Scalar min_rate_bound = (shrinkage_lambda > 0.0) ? 1e-4 : 0.0;
        std::vector<Scalar> lb_alt(alt_dim, min_rate_bound);
        std::vector<Scalar> ub_alt(alt_dim, 10000.0);
        lb_alt[0] = 0.0; ub_alt[0] = 10000.0; // alpha bound

        auto obj_alt = [&](const std::vector<Scalar>& x) -> Scalar {
            Scalar a = x[0];
            std::vector<Scalar> b_rates(K);
            for (size_t k = 0; k < K; ++k) b_rates[k] = x[1 + k];
            Scalar b_bg = has_background ? x[1 + K] : 1.0;

            Scalar lnl = eval_pattern_lnl(a, b_rates, b_bg);
            Scalar penalty = calc_shrinkage_penalty(b_rates);
            return -(lnl - penalty);
        };

        // Grid search for initialization
        std::vector<std::vector<Scalar>> grid_pts;
        for (Scalar cand_a : {0.01, 0.1, 1.0, 5.0}) {
            for (Scalar cand_b : {0.1, 1.0, 3.0}) {
                std::vector<Scalar> pt(alt_dim, cand_b);
                pt[0] = cand_a;
                grid_pts.push_back(pt);
            }
        }
        // Asymmetric points for 2 groups
        if (K == 2) {
            grid_pts.push_back({1.0, 0.1, 3.0});
            grid_pts.push_back({1.0, 3.0, 0.1});
            grid_pts.push_back({1.0, 0.05, 8.0});
            grid_pts.push_back({1.0, 8.0, 0.05});
            grid_pts.push_back({0.1, 0.05, 5.0});
            grid_pts.push_back({0.1, 5.0, 0.05});
            grid_pts.push_back({0.001, 0.001, 5.0});
            grid_pts.push_back({0.001, 5.0, 0.001});
            grid_pts.push_back({0.001, 0.001, 12.0});
            grid_pts.push_back({0.001, 12.0, 0.001});
            if (has_background) {
                for (auto& g : grid_pts) {
                    if (g.size() < alt_dim) g.push_back(1.0);
                }
            }
        }

        std::vector<Scalar> init_alt(alt_dim, 1.0);
        Scalar best_grid_score = 1e20;
        for (const auto& g : grid_pts) {
            Scalar sc = obj_alt(g);
            if (sc < best_grid_score) {
                best_grid_score = sc;
                init_alt = g;
            }
        }

        auto pt_alt = NelderMeadDynamic::minimize(obj_alt, init_alt, lb_alt, ub_alt, 1e-5, 200, grid_pts);
        res.alpha = pt_alt.x[0];
        for (size_t k = 0; k < K; ++k) res.beta[k] = pt_alt.x[1 + k];
        if (has_background) res.beta_background = pt_alt.x[1 + K];

        res.penalized_log_l_alt = -pt_alt.f;
        res.shrinkage_penalty = calc_shrinkage_penalty(res.beta);
        res.log_l_alt = res.penalized_log_l_alt + res.shrinkage_penalty;

        // 2. Fit Overall Null Model (All tested groups share beta_null)
        // Parameter vector x: [alpha, beta_shared, (beta_bg)]
        size_t null_dim = 2 + (has_background ? 1 : 0);
        std::vector<Scalar> lb_null(null_dim, 0.0);
        std::vector<Scalar> ub_null(null_dim, 10000.0);
        lb_null[0] = 0.0; ub_null[0] = 10000.0;

        auto obj_null = [&](const std::vector<Scalar>& x) -> Scalar {
            Scalar a = x[0];
            std::vector<Scalar> b_rates(K, x[1]);
            Scalar b_bg = has_background ? x[2] : 1.0;
            return -eval_pattern_lnl(a, b_rates, b_bg);
        };

        Scalar avg_beta = 0.0;
        for (size_t k = 0; k < K; ++k) avg_beta += res.beta[k];
        avg_beta /= static_cast<Scalar>(K);

        std::vector<Scalar> init_null(null_dim, 1.0);
        init_null[0] = res.alpha;
        init_null[1] = avg_beta;
        if (has_background) init_null[2] = res.beta_background;

        Scalar best_null_score = obj_null(init_null);
        for (Scalar cand_a : {res.alpha, 0.001, 0.1, 1.0}) {
            for (Scalar cand_b : {avg_beta, 0.1, 1.0}) {
                std::vector<Scalar> cand = init_null;
                cand[0] = cand_a;
                cand[1] = cand_b;
                Scalar sc = obj_null(cand);
                if (sc < best_null_score) {
                    best_null_score = sc;
                    init_null = cand;
                }
            }
        }

        auto pt_null = NelderMeadDynamic::minimize(obj_null, init_null, lb_null, ub_null, 1e-5, 120);
        res.alpha_null = pt_null.x[0];
        res.beta_null = pt_null.x[1];
        res.log_l_null = -pt_null.f;

        // Likelihood Ratio Test Calculations
        bool valid_fits = (res.log_l_alt > -1e10 && res.log_l_null > -1e10 && res.penalized_log_l_alt > -1e10 &&
                           std::isfinite(res.log_l_alt) && std::isfinite(res.log_l_null) && std::isfinite(res.penalized_log_l_alt));
        Scalar raw_lrt = 0.0;
        Scalar penalized_lrt = 0.0;
        if (valid_fits) {
            raw_lrt = 2.0 * (res.log_l_alt - res.log_l_null);
            if (raw_lrt < 0.0) raw_lrt = 0.0;
            penalized_lrt = 2.0 * (res.penalized_log_l_alt - res.log_l_null);
            if (penalized_lrt < 0.0) penalized_lrt = 0.0;
        }
        res.raw_lrt = raw_lrt;
        res.penalized_lrt = penalized_lrt;

        // Asymptotic p-value from chi^2_{K-1}
        int df = static_cast<int>(K - 1);
        if (!valid_fits) {
            res.p_value = 1.0;
        } else if (df == 1) {
            res.p_value = std::erfc(std::sqrt(penalized_lrt / 2.0));
        } else {
            // Incomplete gamma chi^2 survival for df >= 2
            // Approx using Wilson-Hilferty transformation
            Scalar s_df = static_cast<Scalar>(df);
            Scalar z = std::pow(penalized_lrt / s_df, 1.0 / 3.0) - (1.0 - 2.0 / (9.0 * s_df));
            z /= std::sqrt(2.0 / (9.0 * s_df));
            res.p_value = 0.5 * std::erfc(z / std::sqrt(2.0));
            if (res.p_value < 0.0) res.p_value = 0.0;
            if (res.p_value > 1.0) res.p_value = 1.0;
        }

        // Total tested branch length
        Scalar tbl = 0.0;
        for (const auto& grp : testable_groups) {
            tbl += grp.total_branch_length;
        }
        res.total_tested_branch_length = (conversion_factor > 0.0) ? (tbl / conversion_factor) : tbl;

        // 3. Pairwise contrasts (for K == 2, reuse omnibus test; for K > 2, fit pairwise constrained models)
        if (K == 2) {
            PairwiseContrastResult pw;
            pw.group_a = 0;
            pw.group_b = 1;
            pw.name_a = testable_groups[0].name;
            pw.name_b = testable_groups[1].name;
            pw.beta_a = res.beta[0];
            pw.beta_b = res.beta[1];
            pw.lrt = res.penalized_lrt;
            pw.p_value = res.p_value;
            pw.p_value_adjusted = res.p_value;
            if (pw.beta_a < pw.beta_b) pw.direction = pw.name_a + " < " + pw.name_b;
            else if (pw.beta_a > pw.beta_b) pw.direction = pw.name_a + " > " + pw.name_b;
            else pw.direction = pw.name_a + " = " + pw.name_b;
            res.pairwise.push_back(pw);
        } else if (K > 2) {
            for (size_t i = 0; i < K; ++i) {
                for (size_t j = i + 1; j < K; ++j) {
                    PairwiseContrastResult pw;
                    pw.group_a = i;
                    pw.group_b = j;
                    pw.name_a = testable_groups[i].name;
                    pw.name_b = testable_groups[j].name;
                    pw.beta_a = res.beta[i];
                    pw.beta_b = res.beta[j];

                    // Null model where beta_i = beta_j
                    // Dimensions: (1 for alpha, K-1 non-synonymous rates, (1 for bg))
                    size_t pw_dim = 1 + (K - 1) + (has_background ? 1 : 0);
                    std::vector<Scalar> lb_pw(pw_dim, 1e-4);
                    std::vector<Scalar> ub_pw(pw_dim, 1000.0);
                    lb_pw[0] = 1e-4; ub_pw[0] = 100.0;

                    auto obj_pw = [&](const std::vector<Scalar>& x) -> Scalar {
                        Scalar a = x[0];
                        std::vector<Scalar> b_rates(K);
                        size_t x_idx = 1;
                        Scalar shared_ij = x[x_idx++];
                        for (size_t k = 0; k < K; ++k) {
                            if (k == i || k == j) {
                                b_rates[k] = shared_ij;
                            } else {
                                b_rates[k] = x[x_idx++];
                            }
                        }
                        Scalar b_bg = has_background ? x.back() : 1.0;
                        Scalar lnl = eval_pattern_lnl(a, b_rates, b_bg);
                        Scalar pen = calc_shrinkage_penalty(b_rates);
                        return -(lnl - pen);
                    };

                    std::vector<Scalar> init_pw(pw_dim, 1.0);
                    init_pw[0] = res.alpha;
                    init_pw[1] = 0.5 * (res.beta[i] + res.beta[j]);
                    size_t p_idx_init = 2;
                    for (size_t k = 0; k < K; ++k) {
                        if (k != i && k != j) {
                            init_pw[p_idx_init++] = res.beta[k];
                        }
                    }
                    if (has_background) init_pw.back() = res.beta_background;

                    auto pt_pw = NelderMeadDynamic::minimize(obj_pw, init_pw, lb_pw, ub_pw, 1e-4, 80);
                    Scalar pw_fit_ll = -pt_pw.f;
                    Scalar pw_lrt = 0.0;
                    if (valid_fits && pw_fit_ll > -1e10 && std::isfinite(pw_fit_ll)) {
                        pw_lrt = 2.0 * (res.penalized_log_l_alt - pw_fit_ll);
                        if (pw_lrt < 0.0) pw_lrt = 0.0;
                        pw.lrt = pw_lrt;
                        pw.p_value = std::erfc(std::sqrt(pw_lrt / 2.0));
                    } else {
                        pw.lrt = 0.0;
                        pw.p_value = 1.0;
                    }

                    if (pw.beta_a < pw.beta_b) {
                        pw.direction = pw.name_a + " < " + pw.name_b;
                    } else if (pw.beta_a > pw.beta_b) {
                        pw.direction = pw.name_a + " > " + pw.name_b;
                    } else {
                        pw.direction = pw.name_a + " = " + pw.name_b;
                    }

                    res.pairwise.push_back(pw);
                }
            }

            // Holm-Bonferroni correction across pairwise tests at this site
            std::vector<size_t> sort_order(res.pairwise.size());
            std::iota(sort_order.begin(), sort_order.end(), 0);
            std::sort(sort_order.begin(), sort_order.end(), [&](size_t a, size_t b) {
                return res.pairwise[a].p_value < res.pairwise[b].p_value;
            });

            size_t m = res.pairwise.size();
            Scalar running_max_adj = 0.0;
            for (size_t rank = 0; rank < m; ++rank) {
                size_t orig_idx = sort_order[rank];
                Scalar unadj_p = res.pairwise[orig_idx].p_value;
                Scalar adj_p = unadj_p * static_cast<Scalar>(m - rank);
                adj_p = std::min(1.0, std::max(running_max_adj, adj_p));
                running_max_adj = adj_p;
                res.pairwise[orig_idx].p_value_adjusted = adj_p;
            }
        }

        // Directional tag formulation
        if (K == 2) {
            Scalar b0 = res.beta[0];
            Scalar b1 = res.beta[1];
            Scalar a = std::max(res.alpha, 1e-4);
            Scalar w0 = b0 / a;
            Scalar w1 = b1 / a;
            std::string n0 = testable_groups[0].name;
            std::string n1 = testable_groups[1].name;

            if (res.p_value <= p_value_threshold) {
                res.classification = "different";
                if (b0 < b1) {
                    if (w0 < 1.0 && w1 > 1.0) {
                        res.contrast_tag = n0 + " ≪ " + n1 + " (Adaptive in " + n1 + ")";
                    } else if (w1 <= 1.0) {
                        res.contrast_tag = n0 + " ≪ " + n1 + " (Relaxed in " + n1 + ")";
                    } else {
                        res.contrast_tag = n0 + " ≪ " + n1 + " (Accelerated in " + n1 + ")";
                    }
                } else if (b0 > b1) {
                    if (w1 < 1.0 && w0 > 1.0) {
                        res.contrast_tag = n0 + " ≫ " + n1 + " (Adaptive in " + n0 + ")";
                    } else if (w0 <= 1.0) {
                        res.contrast_tag = n0 + " ≫ " + n1 + " (Relaxed in " + n0 + ")";
                    } else {
                        res.contrast_tag = n0 + " ≫ " + n1 + " (Accelerated in " + n0 + ")";
                    }
                } else {
                    res.contrast_tag = "Equal";
                }
            } else {
                res.classification = "equal";
                res.contrast_tag = "Equal";
            }
        } else {
            if (res.p_value <= p_value_threshold) {
                res.classification = "different";
                res.contrast_tag = "Heterogeneous (df=" + std::to_string(df) + ")";
            } else {
                res.classification = "equal";
                res.contrast_tag = "Homogeneous";
            }
        }

        return res;
    }

    // Ancestral state reconstruction with branch-set partitioned substitution counts
    ContrastSiteSubstitutions reconstruct_ancestral_pattern(size_t p, const std::vector<size_t>& leaf_to_taxon) const {
        ContrastSiteSubstitutions subs;
        const auto& pattern = aln.patterns[p];
        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        int S = code.num_sense_codons;
        size_t K = testable_groups.size();

        subs.syn_counts.assign(K + (has_background ? 1 : 0), 0);
        subs.nonsyn_counts.assign(K + (has_background ? 1 : 0), 0);

        auto vit = hyphy::core::ancestral::reconstruct_viterbi_states(
            tree, pattern, leaf_to_taxon, branch_transition_matrices, base_mg94_matrix.pi, S
        );
        subs.root_codon = code.sense_codons[vit.root_state];
        subs.root_aa = code.sense_codon_aa[vit.root_state];
        if (vit.is_constant) {
            subs.summary_string = "Conserved";
            subs.summary_string_partitioned = "Conserved";
            return subs;
        }
        const auto& reconstructed_states = vit.node_states;

        // 4. Collect substitution events
        Scalar total_surprise = 0.0;
        Scalar total_mut_branch_length = 0.0;
        std::map<std::string, std::map<std::string, size_t>> group_mut_counts;

        for (int32_t node_id : tree.post_order) {
            if (node_id == tree.root_id) continue;
            const auto& node = tree.nodes[node_id];
            int32_t p_id = node.parent_id;
            int8_t parent_c = reconstructed_states[p_id];
            int8_t child_c = reconstructed_states[node_id];

            if (parent_c != child_c && parent_c >= 0 && child_c >= 0) {
                ContrastSubEvent ev;
                ev.node_id = node_id;
                ev.branch_name = node.name;
                ev.branch_length = (conversion_factor > 0.0) ? (node.branch_length / conversion_factor) : node.branch_length;
                ev.parent_codon_idx = parent_c;
                ev.child_codon_idx = child_c;
                ev.parent_codon = code.sense_codons[parent_c];
                ev.child_codon = code.sense_codons[child_c];
                ev.parent_aa = code.sense_codon_aa[parent_c];
                ev.child_aa = code.sense_codon_aa[child_c];

                int dist = 0;
                for (int pos = 0; pos < 3; ++pos) {
                    if (ev.parent_codon[pos] != ev.child_codon[pos]) dist++;
                }
                ev.nt_distance = dist;
                if (dist == 1) subs.subs_1nt++;
                else if (dist == 2) subs.subs_2nt++;
                else if (dist >= 3) subs.subs_3nt++;

                ev.is_synonymous = (ev.parent_aa == ev.child_aa);

                int gid = node_group_id[node_id];
                if (gid >= 0 && gid < static_cast<int>(K)) {
                    ev.group_id = gid;
                    ev.group_name = testable_groups[gid].name;
                    if (ev.is_synonymous) subs.syn_counts[gid]++;
                    else subs.nonsyn_counts[gid]++;
                } else if (has_background) {
                    ev.group_id = -1;
                    ev.group_name = "Background";
                    if (ev.is_synonymous) subs.syn_counts[K]++;
                    else subs.nonsyn_counts[K]++;
                } else {
                    ev.group_id = 0;
                    ev.group_name = testable_groups[0].name;
                    if (ev.is_synonymous) subs.syn_counts[0]++;
                    else subs.nonsyn_counts[0]++;
                }

                Scalar t_prob = branch_transition_matrices[node_id](parent_c, child_c);
                ev.transition_prob = t_prob;
                ev.surprise_bits = (t_prob > 1e-15) ? -std::log2(t_prob) : 50.0;

                total_surprise += ev.surprise_bits;
                total_mut_branch_length += ev.branch_length;

                std::string mut_str;
                if (ev.is_synonymous) {
                    mut_str = "syn";
                } else {
                    mut_str = std::string(1, ev.parent_aa) + "→" + std::string(1, ev.child_aa);
                    if (ev.nt_distance > 1) mut_str += "[" + std::to_string(ev.nt_distance) + "nt]";
                    if (ev.surprise_bits >= 10.0) mut_str += "⚡";
                }
                group_mut_counts[ev.group_name][mut_str]++;

                subs.events.push_back(std::move(ev));
            }
        }

        if (!subs.events.empty()) {
            subs.total_surprise_bits = total_surprise;
            subs.mean_surprise_bits = total_surprise / subs.events.size();
            subs.mean_branch_length = total_mut_branch_length / subs.events.size();

            std::ostringstream ss_all, ss_part;
            bool first_grp = true;
            for (const auto& [grp_name, mut_map] : group_mut_counts) {
                if (!first_grp) ss_part << " | ";
                first_grp = false;
                ss_part << grp_name << ": ";
                bool first_mut = true;
                for (const auto& [mut_str, cnt] : mut_map) {
                    if (!first_mut) ss_part << ", ";
                    first_mut = false;
                    ss_part << mut_str;
                    if (cnt > 1) ss_part << " (" << cnt << "×)";
                }
            }
            subs.summary_string_partitioned = ss_part.str();
            subs.summary_string = subs.summary_string_partitioned;
        } else {
            subs.summary_string = "Conserved";
            subs.summary_string_partitioned = "Conserved";
        }

        return subs;
    }

    std::vector<ContrastSiteResult> run(bool show_progress = false, bool force_progress = false) {
        size_t num_patterns = aln.patterns.size();
        SiteVisualizer viz(aln.num_codons, num_patterns, SiteVisualizerConfig{
            .method_name = "Contrast-FEL",
            .p_value_threshold = p_value_threshold,
            .pos_label = "different",
            .neg_label = "equal",
            .sig_label = "different",
            .quiet = !show_progress,
            .force_tty = force_progress,
            .enable_chart = true
        });

        std::vector<size_t> leaf_to_taxon(tree.num_nodes(), static_cast<size_t>(-1));
        for (const auto& node : tree.nodes) {
            if (node.is_leaf) {
                auto it = aln.taxon_to_index.find(node.name);
                if (it != aln.taxon_to_index.end()) {
                    leaf_to_taxon[node.id] = it->second;
                }
            }
        }

        std::vector<ContrastSiteResult> pattern_results(num_patterns);

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            pattern_results[p] = analyze_pattern(p, leaf_to_taxon);
            pattern_results[p].subs = reconstruct_ancestral_pattern(p, leaf_to_taxon);

            const auto& r = pattern_results[p];
            SiteSelectionType type = SiteSelectionType::Neutral;
            Scalar metric = r.penalized_lrt;
            std::string detail;
            if (r.p_value <= p_value_threshold) {
                type = SiteSelectionType::Significant;
                detail = "★ LRT = " + SiteVisualizer::format_scalar(metric, 2) + ", p = " + SiteVisualizer::format_pvalue(r.p_value);
                if (!r.contrast_tag.empty() && r.contrast_tag != "Equal") {
                    detail += " (" + r.contrast_tag + ")";
                }
            }
            viz.record_pattern(p, aln.patterns[p].sites, type, r.p_value, metric, detail);
        }

        viz.finish();

        // Expand pattern results to site results
        site_results.resize(aln.num_codons);
        for (size_t c = 0; c < aln.num_codons; ++c) {
            size_t p_idx = aln.site_to_pattern[c];
            site_results[c] = pattern_results[p_idx];
            site_results[c].site_index = c;
            site_results[c].subs.site_index = c;
        }

        // Benjamini-Hochberg FDR correction on penalized p-values across all sites
        size_t M = site_results.size();
        std::vector<Scalar> raw_p(M);
        for (size_t i = 0; i < M; ++i) {
            raw_p[i] = site_results[i].p_value;
        }
        std::vector<Scalar> q_vals = hyphy::core::stats::benjamini_hochberg_fdr(raw_p);
        for (size_t i = 0; i < M; ++i) {
            site_results[i].q_value = q_vals[i];
        }

        return site_results;
    }

    struct Summary {
        size_t total_sites = 0;
        size_t significant_p_sites = 0;
        size_t significant_q_sites = 0;
        std::map<std::string, size_t> contrast_counts;
    };

    Summary get_summary() const {
        Summary s;
        s.total_sites = site_results.size();
        for (const auto& r : site_results) {
            if (r.p_value <= p_value_threshold) {
                s.significant_p_sites++;
                s.contrast_counts[r.contrast_tag]++;
            }
            if (r.q_value <= q_value_threshold) {
                s.significant_q_sites++;
            }
        }
        return s;
    }

    // Modern JSON Schema Serializer
nlohmann::json to_modern_json(
        const std::string& input_filepath = "",
        const std::string& tree_string = "",
        const Provenance& prov = {}
    ) const;

    // Legacy HyPhy 2.5 Contrast-FEL JSON Serializer
nlohmann::json to_legacy_json(
        const std::string& input_filepath = "",
        const std::string& tree_string = ""
    ) const;

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
                   const Provenance& prov = {},
                   JSONFormat format = JSONFormat::ModernV3) const {
        auto j = to_json(input_filepath, tree_string, format, prov);
        std::ofstream out(output_filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open output JSON file: " + output_filepath);
        }
        out << j.dump(2) << std::endl;
    }
};

} // namespace hyphy::analyses
