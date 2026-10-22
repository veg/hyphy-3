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
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <sstream>
#include <queue>
#include <map>
#include <unordered_map>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

// ============================================================================
// Mathematical Statistical Utilities (Regularized Incomplete Gamma & Chi2)
// ============================================================================

using hyphy::core::stats::regularized_gamma_p;
using hyphy::core::stats::regularized_gamma_q;
using hyphy::core::stats::chi2_p_value;

// ============================================================================
// Data Structures for Directional FEL
// ============================================================================

struct TaxonComposition {
    std::string name;
    size_t total_codons = 0;
    std::array<size_t, 4> counts_overall = {0, 0, 0, 0}; // A, C, G, T
    std::array<size_t, 4> counts_cp1 = {0, 0, 0, 0};
    std::array<size_t, 4> counts_cp2 = {0, 0, 0, 0};
    std::array<size_t, 4> counts_cp3 = {0, 0, 0, 0};
    Scalar gc_overall = 0.0;
    Scalar gc_cp3 = 0.0;
};

struct CompositionTestResult {
    Scalar chi2_overall = 0.0;
    int df_overall = 0;
    Scalar p_value_overall = 1.0;

    Scalar chi2_cp1 = 0.0;
    int df_cp1 = 0;
    Scalar p_value_cp1 = 1.0;

    Scalar chi2_cp2 = 0.0;
    int df_cp2 = 0;
    Scalar p_value_cp2 = 1.0;

    Scalar chi2_cp3 = 0.0;
    int df_cp3 = 0;
    Scalar p_value_cp3 = 1.0;

    bool is_heterogeneous = false;
    std::vector<TaxonComposition> taxa;
};

struct RootCandidate {
    size_t rank = 0;
    int32_t node_u = -1;
    int32_t node_v = -1;
    std::string branch_name;
    bool is_terminal = false;
    Scalar edge_length = 0.0;
    Scalar log_l = 0.0;
    Scalar delta_log_l = 0.0;
    Scalar weight = 0.0; // Posterior / Akaike root weight
};

struct AutoRootResult {
    Tree best_tree;
    size_t chosen_candidate_idx = 0;
    std::string chosen_branch_name;
    Scalar chosen_log_l = 0.0;
    Scalar chosen_weight = 0.0;
    std::string strategy = "auto";
    std::vector<RootCandidate> candidates;
};

struct DirectionalRegime {
    int id = 1; // 0 = AT-drift, 1 = Neutral, 2 = GC-drift
    std::string name = "Neutral";
    Scalar psi = 0.0; // Directional GC drift log-odds parameter
    Scalar mean_gc = 0.5;
    size_t branch_count = 0;
    Eigen::Matrix<Scalar, 3, 4> pos_nuc_freqs = Eigen::Matrix<Scalar, 3, 4>::Zero();
    Vector codon_freqs;
    MG94Matrix base_matrix;
};

struct DirectionalSubstitutionEvent {
    int32_t node_id = -1;
    std::string branch_name;
    int regime_id = 1;
    std::string regime_name = "Neutral";
    Scalar branch_length = 0.0;

    int8_t parent_codon_idx = -1;
    int8_t child_codon_idx = -1;
    std::string parent_codon;
    std::string child_codon;

    char parent_aa = '?';
    char child_aa = '?';

    int nt_distance = 0;
    bool is_synonymous = false;

    int gc_change = 0; // +1 (W -> S), -1 (S -> W), 0 (neutral)
    bool is_concordant_with_regime = false;

    Scalar transition_prob = 0.0;
    Scalar surprise_bits = 0.0;
};

struct DirectionalSiteSubstitutions {
    size_t site_index = 0;
    std::string root_codon;
    char root_aa = '?';
    std::vector<DirectionalSubstitutionEvent> events;

    size_t syn_count = 0;
    size_t nonsyn_count = 0;
    size_t nonsyn_drive_concordant = 0;
    size_t nonsyn_drive_discordant = 0;

    size_t subs_1nt = 0;
    size_t subs_2nt = 0;
    size_t subs_3nt = 0;

    Scalar mean_surprise_bits = 0.0;
    Scalar total_surprise_bits = 0.0;
    Scalar site_perplexity = 1.0;

    std::string summary_string = "Conserved";
};

struct DirectionalSiteResult {
    size_t site_index = 0;
    Scalar alpha = 0.0;
    Scalar beta = 0.0;
    Scalar alpha_null = 0.0;
    Scalar lrt = 0.0;
    Scalar p_value = 1.0;
    Scalar total_branch_length = 0.0;
    Scalar log_l_alt = 0.0;
    Scalar log_l_null = 0.0;
    DirectionalSiteSubstitutions subs;
};

struct DirectionalFELSettings {
    Scalar p_value_threshold = 0.1;
    Scalar heterogeneity_threshold = 0.05;
    std::string root_strategy = "auto"; // "auto", "keep", or taxon name
    int num_regimes = 0; // 0 = auto (1 or 3 based on test), 1 = force 1, 3 = force 3
    bool full_model = true;
};

struct DirectionalFELSummary {
    size_t tested_sites = 0;
    size_t positively_selected = 0;
    size_t negatively_selected = 0;
    size_t neutral_sites = 0;
    Scalar runtime_seconds = 0.0;
};

// ============================================================================
// Directional FEL Analyzer Engine
// ============================================================================

class DirectionalFELAnalyzer {
public:
    Tree tree;
    Alignment aln;
    MG94Parameters base_params;
    DirectionalFELSettings settings;

    CompositionTestResult comp_test;
    AutoRootResult root_result;

    std::vector<DirectionalRegime> regimes;
    std::vector<int> branch_regime_ids; // dominant regime id per tree node
    std::vector<std::vector<Scalar>> branch_regime_weights; // soft mixture weights [w_0, w_1, w_2] per tree node
    Vector root_frequencies;

    Scalar p_value_threshold = 0.1;
    Scalar baseline_tree_length = 0.0;
    Scalar conversion_factor = 1.0;
    Scalar global_log_l = 0.0;
    Scalar global_aicc = 0.0;
    std::vector<DirectionalSiteResult> site_results;

    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;

    // Cache of branch transition matrices for ancestral reconstruction
    std::vector<Matrix> branch_transition_matrices;

    DirectionalFELAnalyzer(
        Tree t,
        Alignment a,
        MG94Parameters params = MG94Parameters{},
        DirectionalFELSettings s = DirectionalFELSettings{},
        CompositionTestResult comp = CompositionTestResult{},
        AutoRootResult root_res = AutoRootResult{}
    ) : tree(std::move(t)), aln(std::move(a)), base_params(params), settings(std::move(s)),
        comp_test(std::move(comp)), root_result(std::move(root_res)) {
        p_value_threshold = settings.p_value_threshold;
        setup_regimes_and_matrices();
    }

    // Evaluate single pattern log-likelihood with underflow protection
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
                tree.eval_internal_node_with_model(
                    node_id, P_branches, mg94_model, node_L, node_log_scale,
                    node_L[node_id], child_scale_sum, 1.0, nullptr, P_bridges
                );
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

    // ========================================================================
    // 1. Upfront Compositional Heterogeneity Test
    // ========================================================================
    static CompositionTestResult run_composition_test(const Alignment& aln, Scalar threshold = 0.05) {
        CompositionTestResult res;
        size_t N = aln.num_taxa;
        if (N < 2) return res;

        res.taxa.resize(N);
        std::vector<std::string> valid_chars = {"A", "C", "G", "T"};

        auto char_to_nt = [](char c) -> int {
            char u = std::toupper(static_cast<unsigned char>(c));
            if (u == 'A') return 0;
            if (u == 'C') return 1;
            if (u == 'G') return 2;
            if (u == 'T' || u == 'U') return 3;
            return -1;
        };

        for (size_t t = 0; t < N; ++t) {
            auto& tc = res.taxa[t];
            tc.name = aln.taxon_names[t];
            tc.total_codons = aln.num_codons;
            const auto& raw_seq = aln.sequences[t].raw_seq;

            for (size_t c = 0; c < aln.num_codons; ++c) {
                if (c * 3 + 2 >= raw_seq.size()) break;
                int b0 = char_to_nt(raw_seq[c * 3]);
                int b1 = char_to_nt(raw_seq[c * 3 + 1]);
                int b2 = char_to_nt(raw_seq[c * 3 + 2]);

                if (b0 >= 0) { tc.counts_cp1[b0]++; tc.counts_overall[b0]++; }
                if (b1 >= 0) { tc.counts_cp2[b1]++; tc.counts_overall[b1]++; }
                if (b2 >= 0) { tc.counts_cp3[b2]++; tc.counts_overall[b2]++; }
            }

            size_t tot_all = tc.counts_overall[0] + tc.counts_overall[1] + tc.counts_overall[2] + tc.counts_overall[3];
            size_t gc_all = tc.counts_overall[1] + tc.counts_overall[2];
            tc.gc_overall = (tot_all > 0) ? (static_cast<Scalar>(gc_all) / tot_all) : 0.5;

            size_t tot_cp3 = tc.counts_cp3[0] + tc.counts_cp3[1] + tc.counts_cp3[2] + tc.counts_cp3[3];
            size_t gc_cp3 = tc.counts_cp3[1] + tc.counts_cp3[2];
            tc.gc_cp3 = (tot_cp3 > 0) ? (static_cast<Scalar>(gc_cp3) / tot_cp3) : 0.5;
        }

        // Helper to compute contingency table chi2
        auto calc_table_chi2 = [&](const auto& count_extractor) -> std::pair<Scalar, int> {
            std::vector<std::array<Scalar, 4>> obs(N);
            std::vector<Scalar> row_sums(N, 0.0);
            std::array<Scalar, 4> col_sums = {0.0, 0.0, 0.0, 0.0};
            Scalar grand_total = 0.0;

            for (size_t t = 0; t < N; ++t) {
                const auto& arr = count_extractor(res.taxa[t]);
                for (int b = 0; b < 4; ++b) {
                    obs[t][b] = static_cast<Scalar>(arr[b]);
                    row_sums[t] += obs[t][b];
                    col_sums[b] += obs[t][b];
                    grand_total += obs[t][b];
                }
            }

            if (grand_total <= 0.0) return {0.0, static_cast<int>(3 * (N - 1))};

            Scalar chi2 = 0.0;
            for (size_t t = 0; t < N; ++t) {
                if (row_sums[t] <= 0.0) continue;
                for (int b = 0; b < 4; ++b) {
                    if (col_sums[b] <= 0.0) continue;
                    Scalar expected = (row_sums[t] * col_sums[b]) / grand_total;
                    if (expected > 1e-9) {
                        Scalar diff = obs[t][b] - expected;
                        chi2 += (diff * diff) / expected;
                    }
                }
            }
            int df = static_cast<int>(3 * (N - 1));
            return {chi2, df};
        };

        auto [c_all, df_all] = calc_table_chi2([](const TaxonComposition& tc) { return tc.counts_overall; });
        res.chi2_overall = c_all;
        res.df_overall = df_all;
        res.p_value_overall = chi2_p_value(c_all, df_all);

        auto [c_1, df_1] = calc_table_chi2([](const TaxonComposition& tc) { return tc.counts_cp1; });
        res.chi2_cp1 = c_1;
        res.df_cp1 = df_1;
        res.p_value_cp1 = chi2_p_value(c_1, df_1);

        auto [c_2, df_2] = calc_table_chi2([](const TaxonComposition& tc) { return tc.counts_cp2; });
        res.chi2_cp2 = c_2;
        res.df_cp2 = df_2;
        res.p_value_cp2 = chi2_p_value(c_2, df_2);

        auto [c_3, df_3] = calc_table_chi2([](const TaxonComposition& tc) { return tc.counts_cp3; });
        res.chi2_cp3 = c_3;
        res.df_cp3 = df_3;
        res.p_value_cp3 = chi2_p_value(c_3, df_3);

        res.is_heterogeneous = (res.p_value_overall < threshold || res.p_value_cp3 < threshold);
        return res;
    }

    // ========================================================================
    // 2. Tree Re-rooting via Graph Contraction and BFS Orientation
    // ========================================================================

    struct InternalEdge {
        int32_t to;
        Scalar length;
    };

    static Tree reroot_on_edge(
        const Tree& orig_tree,
        int32_t node_u,
        int32_t node_v,
        Scalar split_ratio = 0.5
    ) {
        size_t M = orig_tree.num_nodes();
        std::vector<std::vector<InternalEdge>> adj(M);

        for (size_t i = 0; i < M; ++i) {
            if (static_cast<int32_t>(i) != orig_tree.root_id) {
                int32_t p = orig_tree.nodes[i].parent_id;
                if (p < 0) continue;
                Scalar len = orig_tree.nodes[i].branch_length;
                adj[p].push_back({static_cast<int32_t>(i), len});
                adj[i].push_back({p, len});
            }
        }

        // Contract root if degree 2
        int32_t old_root = orig_tree.root_id;
        if (old_root >= 0 && adj[old_root].size() == 2) {
            int32_t c0 = adj[old_root][0].to;
            int32_t c1 = adj[old_root][1].to;
            Scalar combined = adj[old_root][0].length + adj[old_root][1].length;

            for (auto& e : adj[c0]) {
                if (e.to == old_root) { e.to = c1; e.length = combined; }
            }
            for (auto& e : adj[c1]) {
                if (e.to == old_root) { e.to = c0; e.length = combined; }
            }
            adj[old_root].clear();
        }

        // Find the edge length between node_u and node_v
        Scalar edge_len = 0.0;
        bool found = false;
        for (const auto& e : adj[node_u]) {
            if (e.to == node_v) {
                edge_len = e.length;
                found = true;
                break;
            }
        }
        if (!found) {
            return orig_tree; // Fallback if edge not adjacent
        }

        // Remove edge between node_u and node_v
        auto remove_edge = [&](int32_t from, int32_t to) {
            auto& vec = adj[from];
            vec.erase(std::remove_if(vec.begin(), vec.end(), [to](const InternalEdge& e) { return e.to == to; }), vec.end());
        };
        remove_edge(node_u, node_v);
        remove_edge(node_v, node_u);

        // Build new rooted tree
        Tree new_tree;
        new_tree.root_id = 0;
        TreeNode root_node;
        root_node.id = 0;
        root_node.name = "Root";
        root_node.is_leaf = false;
        new_tree.nodes.push_back(root_node);

        Scalar len_u = (1.0 - split_ratio) * edge_len;
        Scalar len_v = split_ratio * edge_len;

        // BFS traversal
        struct BFSItem {
            int32_t orig_node;
            int32_t new_node_id;
            int32_t orig_parent;
        };

        std::queue<BFSItem> q;

        // Add node_u as child 1 of root
        {
            int32_t u_new_id = static_cast<int32_t>(new_tree.nodes.size());
            TreeNode u_node;
            u_node.id = u_new_id;
            u_node.name = orig_tree.nodes[node_u].name;
            u_node.is_leaf = orig_tree.nodes[node_u].is_leaf;
            u_node.model_tag = orig_tree.nodes[node_u].model_tag;
            u_node.parent_id = 0;
            u_node.branch_length = len_u;
            new_tree.nodes.push_back(u_node);
            new_tree.nodes[0].children.push_back(u_new_id);
            q.push({node_u, u_new_id, -1});
        }

        // Add node_v as child 2 of root
        {
            int32_t v_new_id = static_cast<int32_t>(new_tree.nodes.size());
            TreeNode v_node;
            v_node.id = v_new_id;
            v_node.name = orig_tree.nodes[node_v].name;
            v_node.is_leaf = orig_tree.nodes[node_v].is_leaf;
            v_node.model_tag = orig_tree.nodes[node_v].model_tag;
            v_node.parent_id = 0;
            v_node.branch_length = len_v;
            new_tree.nodes.push_back(v_node);
            new_tree.nodes[0].children.push_back(v_new_id);
            q.push({node_v, v_new_id, -1});
        }

        while (!q.empty()) {
            auto item = q.front();
            q.pop();

            for (const auto& edge : adj[item.orig_node]) {
                if (edge.to == item.orig_parent) continue;

                int32_t child_new_id = static_cast<int32_t>(new_tree.nodes.size());
                TreeNode ch_node;
                ch_node.id = child_new_id;
                ch_node.name = orig_tree.nodes[edge.to].name;
                ch_node.is_leaf = orig_tree.nodes[edge.to].is_leaf;
                ch_node.model_tag = orig_tree.nodes[edge.to].model_tag;
                ch_node.parent_id = item.new_node_id;
                ch_node.branch_length = edge.length;

                new_tree.nodes.push_back(ch_node);
                new_tree.nodes[item.new_node_id].children.push_back(child_new_id);

                q.push({edge.to, child_new_id, item.orig_node});
            }
        }

        new_tree.compute_post_order();
        new_tree.leaf_name_to_id.clear();
        for (const auto& node : new_tree.nodes) {
            if (node.is_leaf) {
                new_tree.leaf_name_to_id[node.name] = node.id;
            }
        }
        return new_tree;
    }

    // ========================================================================
    // 3. Candidate Edge Enumeration & Auto-Rooting Likelihood Scoring
    // ========================================================================

    static AutoRootResult auto_select_root(
        const Tree& input_tree,
        const Alignment& aln,
        const GTRParameters& gtr_rates,
        const DirectionalFELSettings& settings
    ) {
        AutoRootResult result;
        result.strategy = settings.root_strategy;
        if (settings.root_strategy == "keep" || input_tree.has_wheatstone_bridges()) {
            result.best_tree = input_tree;
            result.chosen_candidate_idx = 0;
            result.chosen_branch_name = input_tree.has_wheatstone_bridges()
                ? "Ohmic Network Root (Preserved)"
                : "User Input Root (Preserved)";
            result.chosen_log_l = 0.0;
            result.chosen_weight = 1.0;
            return result;
        }
        size_t M = input_tree.num_nodes();

        // Build adjacency to list unique undirected edges
        std::vector<std::vector<InternalEdge>> adj(M);
        for (size_t i = 0; i < M; ++i) {
            if (static_cast<int32_t>(i) != input_tree.root_id) {
                int32_t p = input_tree.nodes[i].parent_id;
                if (p < 0) continue;
                Scalar len = input_tree.nodes[i].branch_length;
                adj[p].push_back({static_cast<int32_t>(i), len});
                adj[i].push_back({p, len});
            }
        }

        // Contract old degree-2 root
        int32_t old_root = input_tree.root_id;
        if (old_root >= 0 && adj[old_root].size() == 2) {
            int32_t c0 = adj[old_root][0].to;
            int32_t c1 = adj[old_root][1].to;
            Scalar combined = adj[old_root][0].length + adj[old_root][1].length;

            for (auto& e : adj[c0]) {
                if (e.to == old_root) { e.to = c1; e.length = combined; }
            }
            for (auto& e : adj[c1]) {
                if (e.to == old_root) { e.to = c0; e.length = combined; }
            }
            adj[old_root].clear();
        }

        struct UniqueEdge {
            int32_t u;
            int32_t v;
            Scalar length;
            std::string label;
            bool is_terminal;
        };

        std::vector<UniqueEdge> unique_edges;
        for (int32_t u = 0; u < static_cast<int32_t>(M); ++u) {
            for (const auto& e : adj[u]) {
                if (u < e.to) {
                    int32_t v = e.to;
                    bool u_leaf = input_tree.nodes[u].is_leaf;
                    bool v_leaf = input_tree.nodes[v].is_leaf;
                    std::string label;
                    bool is_term = false;
                    if (v_leaf) {
                        label = input_tree.nodes[v].name;
                        is_term = true;
                    } else if (u_leaf) {
                        label = input_tree.nodes[u].name;
                        is_term = true;
                    } else {
                        label = input_tree.nodes[u].name + " - " + input_tree.nodes[v].name;
                        is_term = false;
                    }
                    unique_edges.push_back({u, v, e.length, label, is_term});
                }
            }
        }

        // Build base directional model for fast scoring
        MG94Parameters base_p;
        base_p.theta_AC = gtr_rates.theta_AC;
        base_p.theta_AT = gtr_rates.theta_AT;
        base_p.theta_CG = gtr_rates.theta_CG;
        base_p.theta_CT = gtr_rates.theta_CT;
        base_p.theta_GT = gtr_rates.theta_GT;
        base_p.alpha = 1.0;
        base_p.beta = 1.0;

        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        MG94Matrix base_matrix;
        base_matrix.update(base_p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

        // Precompute leaf GC contents
        std::vector<Scalar> leaf_gc(aln.num_taxa, 0.5);
        for (size_t t = 0; t < aln.num_taxa; ++t) {
            const auto& seq = aln.sequences[t].raw_seq;
            size_t gc = 0, tot = 0;
            for (char c : seq) {
                char u = std::toupper(static_cast<unsigned char>(c));
                if (u == 'G' || u == 'C') { gc++; tot++; }
                else if (u == 'A' || u == 'T' || u == 'U') { tot++; }
            }
            if (tot > 0) leaf_gc[t] = static_cast<Scalar>(gc) / tot;
        }

        // Evaluate each candidate edge
        size_t E = unique_edges.size();
        std::vector<RootCandidate> cand_list(E);

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < E; ++i) {
            const auto& ue = unique_edges[i];
            Tree cand_tree = reroot_on_edge(input_tree, ue.u, ue.v, 0.5);

            // Reconstruct node GC contents on cand_tree
            size_t num_nodes = cand_tree.num_nodes();
            std::vector<Scalar> node_gc(num_nodes, 0.5);
            auto leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(cand_tree, aln);

            for (int32_t nid : cand_tree.post_order) {
                const auto& node = cand_tree.nodes[nid];
                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[nid];
                    if (t_idx != static_cast<size_t>(-1)) {
                        node_gc[nid] = leaf_gc[t_idx];
                    }
                } else if (!node.children.empty()) {
                    Scalar sum = 0.0;
                    for (int32_t ch : node.children) sum += node_gc[ch];
                    node_gc[nid] = sum / node.children.size();
                }
            }

            // Directional transition matrices per branch
            std::vector<Matrix> P_branches(num_nodes);
            std::vector<MG94Matrix> branch_mgs(num_nodes);
            for (const auto& node : cand_tree.nodes) {
                if (node.id != cand_tree.root_id) {
                    int32_t p = node.parent_id;
                    Scalar delta_gc = node_gc[node.id] - node_gc[p];
                    Scalar psi = 2.0 * delta_gc; // Directional drift proxy

                    // Construct directional pos_nuc_freqs
                    Eigen::Matrix<Scalar, 3, 4> pos_freqs = aln.pos_nuc_frequencies;
                    for (int pos = 0; pos < 3; ++pos) {
                        pos_freqs(pos, 1) *= std::exp(0.5 * psi); // C
                        pos_freqs(pos, 2) *= std::exp(0.5 * psi); // G
                        pos_freqs(pos, 0) *= std::exp(-0.5 * psi); // A
                        pos_freqs(pos, 3) *= std::exp(-0.5 * psi); // T
                        Scalar row_sum = pos_freqs.row(pos).sum();
                        if (row_sum > 1e-9) pos_freqs.row(pos) /= row_sum;
                    }

                    Vector codon_f = aln.codon_frequencies_f3x4;
                    for (int c = 0; c < gcode.num_sense_codons; ++c) {
                        const auto& trip = gcode.sense_codons[c];
                        int b0 = (trip[0] == 'C' ? 1 : (trip[0] == 'G' ? 2 : (trip[0] == 'A' ? 0 : 3)));
                        int b1 = (trip[1] == 'C' ? 1 : (trip[1] == 'G' ? 2 : (trip[1] == 'A' ? 0 : 3)));
                        int b2 = (trip[2] == 'C' ? 1 : (trip[2] == 'G' ? 2 : (trip[2] == 'A' ? 0 : 3)));
                        codon_f(c) = pos_freqs(0, b0) * pos_freqs(1, b1) * pos_freqs(2, b2);
                    }
                    Scalar sum_cf = codon_f.sum();
                    if (sum_cf > 1e-9) codon_f /= sum_cf;

                    branch_mgs[node.id].update(base_p, pos_freqs, codon_f, gcode);
                    P_branches[node.id] = branch_mgs[node.id].transition_matrix(node.branch_length);
                }
            }
            std::vector<std::array<Matrix, 3>> P_bridges(cand_tree.bridges.size());
            for (size_t b = 0; b < cand_tree.bridges.size(); ++b) {
                const auto& wb = cand_tree.bridges[b];
                int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
                const MG94Matrix& m_br = (bid >= 0 && bid != cand_tree.root_id && static_cast<size_t>(bid) < branch_mgs.size())
                    ? branch_mgs[bid]
                    : base_matrix;
                for (int r = 0; r < 3; ++r) {
                    P_bridges[b][r] = m_br.transition_matrix(wb.branch_lengths[r]);
                }
            }

            // Pattern likelihood evaluation
            Scalar tree_lnl = 0.0;
            int S = gcode.num_sense_codons;
            std::vector<Vector> node_L(num_nodes, Vector(S));
            std::vector<double> node_log_scale(num_nodes, 0.0);

            for (size_t p = 0; p < aln.patterns.size(); ++p) {
                const auto& pat = aln.patterns[p];
                Scalar pat_l = compute_pattern_log_likelihood(
                    cand_tree, pat, leaf_to_taxon, base_matrix, P_branches, node_L, node_log_scale,
                    cand_tree.bridges.empty() ? nullptr : &P_bridges
                );
                tree_lnl += pat.weight * pat_l;
            }

            cand_list[i].node_u = ue.u;
            cand_list[i].node_v = ue.v;
            cand_list[i].branch_name = ue.label;
            cand_list[i].is_terminal = ue.is_terminal;
            cand_list[i].edge_length = ue.length;
            cand_list[i].log_l = tree_lnl;
        }

        // Rank by log-likelihood
        Scalar max_l = -1e20;
        for (const auto& c : cand_list) {
            if (c.log_l > max_l) max_l = c.log_l;
        }

        Scalar weight_sum = 0.0;
        for (auto& c : cand_list) {
            c.delta_log_l = max_l - c.log_l;
            c.weight = std::exp(-c.delta_log_l);
            weight_sum += c.weight;
        }
        if (weight_sum > 0.0) {
            for (auto& c : cand_list) c.weight /= weight_sum;
        }

        std::sort(cand_list.begin(), cand_list.end(), [](const RootCandidate& a, const RootCandidate& b) {
            return a.log_l > b.log_l;
        });

        for (size_t r = 0; r < cand_list.size(); ++r) {
            cand_list[r].rank = r + 1;
        }

        result.candidates = cand_list;

        // Select root based on strategy
        size_t chosen_idx = 0;
        if (settings.root_strategy == "keep") {
            // Keep input tree as is
            result.best_tree = input_tree;
            result.chosen_candidate_idx = 0;
            result.chosen_branch_name = "User Input Root (Preserved)";
            result.chosen_log_l = cand_list.empty() ? 0.0 : cand_list[0].log_l;
            result.chosen_weight = cand_list.empty() ? 1.0 : cand_list[0].weight;
            return result;
        } else if (settings.root_strategy != "auto" && !settings.root_strategy.empty()) {
            // Root at user-specified outgroup
            std::string target_name = settings.root_strategy;
            bool found_outgroup = false;
            for (size_t i = 0; i < cand_list.size(); ++i) {
                if (cand_list[i].branch_name == target_name) {
                    chosen_idx = i;
                    found_outgroup = true;
                    break;
                }
            }
            if (!found_outgroup) {
                chosen_idx = 0; // fallback to ML
            }
        }

        const auto& best_cand = cand_list[chosen_idx];
        result.chosen_candidate_idx = chosen_idx;
        result.chosen_branch_name = best_cand.branch_name;
        result.chosen_log_l = best_cand.log_l;
        result.chosen_weight = best_cand.weight;
        result.best_tree = reroot_on_edge(input_tree, best_cand.node_u, best_cand.node_v, 0.5);

        return result;
    }

    // ========================================================================
    // 4. Directional Regimes Setup & Matrix Preparation
    // ========================================================================

    void setup_regimes_and_matrices() {
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();

        // Determine number of regimes
        int K = 1;
        if (settings.num_regimes == 3) {
            K = 3;
        } else if (settings.num_regimes == 1) {
            K = 1;
        } else {
            // Auto: 3 if heterogeneous, 1 if homogeneous
            K = comp_test.is_heterogeneous ? 3 : 1;
        }

        // Calculate node GC contents to assign branches
        std::vector<Scalar> node_gc(num_nodes, 0.5);
        auto leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        for (int32_t nid : tree.post_order) {
            const auto& node = tree.nodes[nid];
            if (node.is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                if (t_idx != static_cast<size_t>(-1) && t_idx < comp_test.taxa.size()) {
                    node_gc[nid] = comp_test.taxa[t_idx].gc_cp3;
                }
            } else if (!node.children.empty()) {
                Scalar sum = 0.0;
                for (int32_t ch : node.children) sum += node_gc[ch];
                node_gc[nid] = sum / node.children.size();
            }
        }

        // Alignment background GC at CP3
        Scalar bg_gc = 0.5;
        if (!comp_test.taxa.empty()) {
            Scalar sum = 0.0;
            for (const auto& tc : comp_test.taxa) sum += tc.gc_cp3;
            bg_gc = sum / comp_test.taxa.size();
        }

        // Compute branch delta_gc relative to alignment background
        std::vector<Scalar> branch_deltas;
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                Scalar d = node_gc[node.id] - bg_gc;
                branch_deltas.push_back(d);
            }
        }

        Scalar mean_delta = 0.0, sd_delta = 0.0;
        if (!branch_deltas.empty()) {
            for (Scalar d : branch_deltas) mean_delta += d;
            mean_delta /= branch_deltas.size();
            Scalar sum_sq = 0.0;
            for (Scalar d : branch_deltas) sum_sq += (d - mean_delta) * (d - mean_delta);
            sd_delta = std::sqrt(sum_sq / branch_deltas.size());
        }

        Scalar sigma = std::max(0.015, sd_delta);
        Scalar tau = std::max(0.02, 0.4 * sigma);
        Scalar kappa = 3.5 / sigma;

        // Assign regimes via Softmax mixture
        branch_regime_ids.assign(num_nodes, 1); // Default Neutral
        branch_regime_weights.assign(num_nodes, {0.0, 1.0, 0.0});
        std::array<Scalar, 3> soft_regime_counts = {0.0, 0.0, 0.0};
        std::array<Scalar, 3> soft_regime_gc_sum = {0.0, 0.0, 0.0};

        if (K == 3) {
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    Scalar d = node_gc[node.id] - bg_gc;
                    // Softmax logits for 3 regimes: 0 (AT-drift), 1 (Neutral), 2 (GC-drift)
                    Scalar ell_0 = kappa * (-d - tau);
                    Scalar ell_1 = 0.0;
                    Scalar ell_2 = kappa * (d - tau);

                    Scalar max_l = std::max(ell_0, std::max(ell_1, ell_2));
                    Scalar e0 = std::exp(ell_0 - max_l);
                    Scalar e1 = std::exp(ell_1 - max_l);
                    Scalar e2 = std::exp(ell_2 - max_l);
                    Scalar sum_e = e0 + e1 + e2;
                    Scalar w0 = e0 / sum_e;
                    Scalar w1 = e1 / sum_e;
                    Scalar w2 = e2 / sum_e;

                    branch_regime_weights[node.id] = {w0, w1, w2};

                    int r_id = 1;
                    if (w2 > w1 && w2 > w0) r_id = 2;
                    else if (w0 > w1 && w0 > w2) r_id = 0;
                    branch_regime_ids[node.id] = r_id;

                    soft_regime_counts[0] += w0;
                    soft_regime_counts[1] += w1;
                    soft_regime_counts[2] += w2;

                    soft_regime_gc_sum[0] += w0 * node_gc[node.id];
                    soft_regime_gc_sum[1] += w1 * node_gc[node.id];
                    soft_regime_gc_sum[2] += w2 * node_gc[node.id];
                }
            }
        } else {
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    branch_regime_weights[node.id] = {0.0, 1.0, 0.0};
                    branch_regime_ids[node.id] = 1;
                    soft_regime_counts[1] += 1.0;
                    soft_regime_gc_sum[1] += node_gc[node.id];
                }
            }
        }

        // Build regimes
        regimes.clear();
        Scalar avg_gc_1 = (soft_regime_counts[1] > 1e-4) ? (soft_regime_gc_sum[1] / soft_regime_counts[1]) : bg_gc;
        Scalar avg_gc_2 = (soft_regime_counts[2] > 1e-4) ? (soft_regime_gc_sum[2] / soft_regime_counts[2]) : (bg_gc + 0.1);
        Scalar avg_gc_0 = (soft_regime_counts[0] > 1e-4) ? (soft_regime_gc_sum[0] / soft_regime_counts[0]) : (bg_gc - 0.1);

        avg_gc_1 = std::clamp(avg_gc_1, 0.05, 0.95);
        avg_gc_2 = std::clamp(avg_gc_2, 0.05, 0.95);
        avg_gc_0 = std::clamp(avg_gc_0, 0.05, 0.95);

        auto make_regime = [&](int id, const std::string& name, Scalar psi_val, size_t count, Scalar avg_gc) {
            DirectionalRegime reg;
            reg.id = id;
            reg.name = name;
            reg.psi = psi_val;
            reg.mean_gc = avg_gc;
            reg.branch_count = count;

            // Target pos_nuc_freqs
            reg.pos_nuc_freqs = aln.pos_nuc_frequencies;
            for (int p = 0; p < 3; ++p) {
                reg.pos_nuc_freqs(p, 1) *= std::exp(0.5 * psi_val); // C
                reg.pos_nuc_freqs(p, 2) *= std::exp(0.5 * psi_val); // G
                reg.pos_nuc_freqs(p, 0) *= std::exp(-0.5 * psi_val); // A
                reg.pos_nuc_freqs(p, 3) *= std::exp(-0.5 * psi_val); // T
                Scalar rsum = reg.pos_nuc_freqs.row(p).sum();
                if (rsum > 1e-9) reg.pos_nuc_freqs.row(p) /= rsum;
            }

            // Target codon freqs F3x4
            reg.codon_freqs.resize(S);
            for (int c = 0; c < S; ++c) {
                const auto& trip = gcode.sense_codons[c];
                int b0 = (trip[0] == 'C' ? 1 : (trip[0] == 'G' ? 2 : (trip[0] == 'A' ? 0 : 3)));
                int b1 = (trip[1] == 'C' ? 1 : (trip[1] == 'G' ? 2 : (trip[1] == 'A' ? 0 : 3)));
                int b2 = (trip[2] == 'C' ? 1 : (trip[2] == 'G' ? 2 : (trip[2] == 'A' ? 0 : 3)));
                reg.codon_freqs(c) = reg.pos_nuc_freqs(0, b0) * reg.pos_nuc_freqs(1, b1) * reg.pos_nuc_freqs(2, b2);
            }
            Scalar sum_cf = reg.codon_freqs.sum();
            if (sum_cf > 1e-9) reg.codon_freqs /= sum_cf;

            reg.base_matrix.update(base_params, reg.pos_nuc_freqs, reg.codon_freqs, gcode);
            return reg;
        };

        if (K == 3) {
            Scalar psi_gc = std::log((avg_gc_2 / (1.0 - avg_gc_2)) / (avg_gc_1 / (1.0 - avg_gc_1)));
            Scalar psi_at = std::log((avg_gc_0 / (1.0 - avg_gc_0)) / (avg_gc_1 / (1.0 - avg_gc_1)));
            psi_gc = std::clamp(psi_gc, 0.05, 3.0);
            psi_at = std::clamp(psi_at, -3.0, -0.05);

            regimes.push_back(make_regime(0, "AT-drift", psi_at, static_cast<size_t>(std::round(soft_regime_counts[0])), avg_gc_0));
            regimes.push_back(make_regime(1, "Neutral", 0.0, static_cast<size_t>(std::round(soft_regime_counts[1])), avg_gc_1));
            regimes.push_back(make_regime(2, "GC-drift", psi_gc, static_cast<size_t>(std::round(soft_regime_counts[2])), avg_gc_2));
        } else {
            regimes.push_back(make_regime(1, "Neutral", 0.0, static_cast<size_t>(std::round(soft_regime_counts[1])), avg_gc_1));
        }

        // Root frequencies: Neutral regime codon frequencies
        root_frequencies = (K == 3) ? regimes[1].codon_freqs : regimes[0].codon_freqs;

        // Tree length scale conversion
        Scalar Q_scale = regimes[K == 3 ? 1 : 0].base_matrix.scale_factor;
        conversion_factor = (Q_scale > 1e-12) ? (3.0 / Q_scale) : 1.0;

        baseline_tree_length = 0.0;
        for (auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                baseline_tree_length += node.branch_length;
                node.branch_length *= conversion_factor;
            }
        }
        tree.scale_bridge_lengths(conversion_factor);

        // Precompute branch transition matrices for reconstruction
        branch_transition_matrices.resize(num_nodes);
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                if (K == 3) {
                    const auto& w = branch_regime_weights[node.id];
                    if (w[1] > 0.99) {
                        branch_transition_matrices[node.id] = regimes[1].base_matrix.transition_matrix(node.branch_length);
                    } else if (w[2] > 0.99) {
                        branch_transition_matrices[node.id] = regimes[2].base_matrix.transition_matrix(node.branch_length);
                    } else if (w[0] > 0.99) {
                        branch_transition_matrices[node.id] = regimes[0].base_matrix.transition_matrix(node.branch_length);
                    } else {
                        Matrix P = w[1] * regimes[1].base_matrix.transition_matrix(node.branch_length);
                        if (w[2] > 1e-4) P += w[2] * regimes[2].base_matrix.transition_matrix(node.branch_length);
                        if (w[0] > 1e-4) P += w[0] * regimes[0].base_matrix.transition_matrix(node.branch_length);
                        branch_transition_matrices[node.id] = std::move(P);
                    }
                } else {
                    branch_transition_matrices[node.id] = regimes[0].base_matrix.transition_matrix(node.branch_length);
                }
            }
        }
    }

    // Factory: Phase 1 GTR + Composition Test + Auto-Root + Global MG94
    static DirectionalFELAnalyzer create_and_fit(
        Tree input_tree,
        Alignment aln,
        DirectionalFELSettings settings = DirectionalFELSettings{},
        std::function<void(const std::string&, double)> progress_cb = nullptr
    ) {
        if (progress_cb) progress_cb("Phase 1A: Compositional Heterogeneity Analysis", 0.1);
        auto comp_test = run_composition_test(aln, settings.heterogeneity_threshold);

        if (progress_cb) progress_cb("Phase 1B: Fitting Nucleotide GTR Model", 0.25);
        GTRFitter gtr_fitter(input_tree, aln);
        auto gtr_res = gtr_fitter.fit();

        if (progress_cb) progress_cb("Phase 1C: Directional Auto-Rooting Selection", 0.4);
        auto root_res = auto_select_root(gtr_res.tree, aln, gtr_res.params, settings);

        if (progress_cb) progress_cb("Phase 1D: Refining Global Directional MG94 Model", 0.55);
        MG94Parameters base_p;
        base_p.theta_AC = gtr_res.params.theta_AC;
        base_p.theta_AT = gtr_res.params.theta_AT;
        base_p.theta_CG = gtr_res.params.theta_CG;
        base_p.theta_CT = gtr_res.params.theta_CT;
        base_p.theta_GT = gtr_res.params.theta_GT;

        MG94Fitter mg_fitter(root_res.best_tree, aln);
        FitResult mg_res;
        if (settings.full_model) {
            mg_res = mg_fitter.fit_full_model(1.0, base_p, progress_cb);
            base_p = mg_res.params;
        } else {
            mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
            base_p.alpha = 1.0;
            base_p.beta = mg_res.x_opt(0);
        }

        DirectionalFELAnalyzer analyzer(mg_res.tree, std::move(aln), base_p, settings, comp_test, root_res);
        analyzer.global_log_l = mg_res.log_likelihood;
        size_t K_mg = 5 + 1 + (mg_res.tree.num_nodes() - 1) + (analyzer.regimes.size() > 1 ? 2 : 0);
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

    // ========================================================================
    // 5. Pattern-Level Site Selection Analysis (Phase 2)
    // ========================================================================

    DirectionalSiteResult analyze_pattern(size_t p, const std::vector<size_t>& leaf_to_taxon) {
        DirectionalSiteResult res;
        res.site_index = p;

        const auto& pattern = aln.patterns[p];
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        size_t K = regimes.size();

        // Check if invariant
        int8_t first_state = -1;
        bool is_invariant = true;
        for (size_t nid = 0; nid < num_nodes; ++nid) {
            if (tree.nodes[nid].is_leaf) {
                size_t t_idx = leaf_to_taxon[nid];
                if (t_idx != static_cast<size_t>(-1)) {
                    int8_t st = pattern.states[t_idx];
                    if (st >= 0 && st < S) {
                        if (first_state == -1) first_state = st;
                        else if (first_state != st) { is_invariant = false; break; }
                    }
                }
            }
        }

        if (is_invariant && first_state >= 0) {
            res.alpha = 0.0;
            res.beta = 0.0;
            res.alpha_null = 0.0;
            res.lrt = 0.0;
            res.p_value = 1.0;
            res.total_branch_length = 0.0;
            res.log_l_alt = 0.0;
            res.log_l_null = 0.0;
            return res;
        }

        // Thread-local likelihood scratchpad
        std::vector<Vector> thread_node_L(num_nodes, Vector(S));
        std::vector<double> thread_node_scale(num_nodes, 0.0);
        std::vector<Matrix> P_branches(num_nodes);
        std::vector<MG94Matrix> regime_matrices(K);

        auto eval_pattern_lnl = [&](Scalar a, Scalar b) -> Scalar {
            MG94Parameters p_curr = base_params;
            p_curr.alpha = a;
            p_curr.beta = b;

            for (size_t k = 0; k < K; ++k) {
                regime_matrices[k].update(p_curr, regimes[k].pos_nuc_freqs, regimes[k].codon_freqs, gcode);
            }

            auto eval_branch_P = [&](int32_t nid, Scalar bl) -> Matrix {
                if (K == 3 && nid >= 0 && static_cast<size_t>(nid) < branch_regime_weights.size()) {
                    const auto& w = branch_regime_weights[nid];
                    if (w[1] > 0.99) {
                        return regime_matrices[1].transition_matrix(bl);
                    } else if (w[2] > 0.99) {
                        return regime_matrices[2].transition_matrix(bl);
                    } else if (w[0] > 0.99) {
                        return regime_matrices[0].transition_matrix(bl);
                    } else {
                        Matrix P = w[1] * regime_matrices[1].transition_matrix(bl);
                        if (w[2] > 1e-4) P += w[2] * regime_matrices[2].transition_matrix(bl);
                        if (w[0] > 1e-4) P += w[0] * regime_matrices[0].transition_matrix(bl);
                        return P;
                    }
                }
                return regime_matrices[0].transition_matrix(bl);
            };

            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P_branches[node.id] = eval_branch_P(node.id, node.branch_length);
                }
            }
            std::vector<std::array<Matrix, 3>> P_bridges(tree.bridges.size());
            for (size_t b_idx = 0; b_idx < tree.bridges.size(); ++b_idx) {
                const auto& wb = tree.bridges[b_idx];
                int32_t bid = (wb.bridge_node_id >= 0 && wb.bridge_node_id != wb.hub_node_id) ? wb.bridge_node_id : wb.hub_node_id;
                for (int r = 0; r < 3; ++r) {
                    P_bridges[b_idx][r] = eval_branch_P(bid, wb.branch_lengths[r]);
                }
            }

            const auto& mg_root = (K == 3) ? regime_matrices[1] : regime_matrices[0];
            return compute_pattern_log_likelihood(
                tree, pattern, leaf_to_taxon, mg_root, P_branches, thread_node_L, thread_node_scale,
                tree.bridges.empty() ? nullptr : &P_bridges
            );
        };

        // Grid search initialization
        const std::vector<std::pair<Scalar, Scalar>> grid_points = {
            {0.01, 0.01}, {0.1, 0.1}, {0.5, 0.5}, {1.0, 1.0}, {2.0, 2.0},
            {0.1, 1.0}, {1.0, 0.1}, {0.5, 2.0}, {2.0, 0.5}, {5.0, 1.0}, {1.0, 5.0}
        };

        Scalar best_grid_lnl = -1e20;
        Scalar best_a = 1.0, best_b = 1.0;
        for (const auto& pt : grid_points) {
            Scalar lnl = eval_pattern_lnl(pt.first, pt.second);
            if (lnl > best_grid_lnl) {
                best_grid_lnl = lnl;
                best_a = pt.first;
                best_b = pt.second;
            }
        }

        // Alternative Model: Nelder-Mead 2D on (alpha, beta)
        auto obj_alt = [&](const std::array<Scalar, 2>& x) -> Scalar {
            return -eval_pattern_lnl(x[0], x[1]);
        };
        std::array<Scalar, 2> init_alt = {best_a, best_b};
        std::array<Scalar, 2> lb_alt = {1e-4, 1e-4};
        std::array<Scalar, 2> ub_alt = {1000.0, 1000.0};
        auto pt_alt = NelderMeadND<2>::minimize(obj_alt, init_alt, lb_alt, ub_alt, 1e-4, 80);

        res.alpha = pt_alt.x[0];
        res.beta = pt_alt.x[1];
        res.log_l_alt = -pt_alt.f;

        // Null Model: beta = alpha. Brent 1D search
        auto obj_null = [&](Scalar a) -> Scalar {
            return -eval_pattern_lnl(a, a);
        };

        Scalar best_null_lnl = -obj_null(best_a);
        Scalar best_null_a = best_a;
        for (Scalar cand_a : {0.01, 0.1, 0.5, 1.0, 2.0, 5.0}) {
            Scalar lnl = -obj_null(cand_a);
            if (lnl > best_null_lnl) {
                best_null_lnl = lnl;
                best_null_a = cand_a;
            }
        }

        Scalar ax = 1e-4;
        Scalar bx = std::clamp(best_null_a, 0.01, 45.0);
        Scalar cx = 50.0;
        auto [opt_null_a, opt_null_nll] = Brent1D::minimize(obj_null, ax, bx, cx, 1e-4, 60);

        res.alpha_null = opt_null_a;
        res.log_l_null = -opt_null_nll;

        // Likelihood Ratio Test
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

    // ========================================================================
    // 6. Ancestral State Reconstruction & Directed Substitutions
    // ========================================================================

    DirectionalSiteSubstitutions reconstruct_ancestral_pattern(size_t p, const std::vector<size_t>& leaf_to_taxon) const {
        DirectionalSiteSubstitutions subs;
        const auto& pattern = aln.patterns[p];
        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();

        auto vit = hyphy::core::ancestral::reconstruct_viterbi_states(
            tree, pattern, leaf_to_taxon, branch_transition_matrices, root_frequencies, S
        );
        subs.root_codon = gcode.sense_codons[vit.root_state];
        subs.root_aa = gcode.sense_codon_aa[vit.root_state];
        if (vit.is_constant) {
            subs.summary_string = "Conserved";
            return subs;
        }
        const auto& reconstructed = vit.node_states;

        // 5. Collect substitution events
        auto count_gc = [](const std::string& codon) -> int {
            int n = 0;
            for (char c : codon) {
                char u = std::toupper(static_cast<unsigned char>(c));
                if (u == 'G' || u == 'C') n++;
            }
            return n;
        };

        Scalar total_surprise = 0.0;
        for (size_t nid = 0; nid < num_nodes; ++nid) {
            if (static_cast<int32_t>(nid) == tree.root_id) continue;
            const auto& node = tree.nodes[nid];
            if (node.parent_id < 0) continue;
            int8_t p_st = reconstructed[node.parent_id];
            int8_t ch_st = reconstructed[nid];

            if (p_st != ch_st && p_st >= 0 && ch_st >= 0) {
                DirectionalSubstitutionEvent ev;
                ev.node_id = static_cast<int32_t>(nid);
                ev.branch_name = node.name;
                ev.regime_id = branch_regime_ids[nid];
                ev.regime_name = (ev.regime_id == 0 ? "AT-drift" : (ev.regime_id == 2 ? "GC-drift" : "Neutral"));
                ev.branch_length = (conversion_factor > 0.0) ? (node.branch_length / conversion_factor) : node.branch_length;
                ev.parent_codon_idx = p_st;
                ev.child_codon_idx = ch_st;
                ev.parent_codon = gcode.sense_codons[p_st];
                ev.child_codon = gcode.sense_codons[ch_st];
                ev.parent_aa = gcode.sense_codon_aa[p_st];
                ev.child_aa = gcode.sense_codon_aa[ch_st];

                const auto& diff = gcode.diff_matrix[p_st][ch_st];
                ev.nt_distance = diff.num_diffs;
                ev.is_synonymous = diff.synonymous;

                int gc_parent = count_gc(ev.parent_codon);
                int gc_child = count_gc(ev.child_codon);
                ev.gc_change = gc_child - gc_parent;

                if (ev.regime_id == 2 && ev.gc_change > 0) ev.is_concordant_with_regime = true;
                else if (ev.regime_id == 0 && ev.gc_change < 0) ev.is_concordant_with_regime = true;
                else ev.is_concordant_with_regime = false;

                const Matrix& P = branch_transition_matrices[nid];
                ev.transition_prob = std::max(1e-12, static_cast<double>(P(p_st, ch_st)));
                ev.surprise_bits = -std::log2(ev.transition_prob);
                total_surprise += ev.surprise_bits;

                if (ev.is_synonymous) {
                    subs.syn_count++;
                } else {
                    subs.nonsyn_count++;
                    if (ev.is_concordant_with_regime) subs.nonsyn_drive_concordant++;
                    else subs.nonsyn_drive_discordant++;
                }

                if (ev.nt_distance == 1) subs.subs_1nt++;
                else if (ev.nt_distance == 2) subs.subs_2nt++;
                else if (ev.nt_distance >= 3) subs.subs_3nt++;

                subs.events.push_back(std::move(ev));
            }
        }

        size_t total_events = subs.events.size();
        subs.total_surprise_bits = total_surprise;
        subs.mean_surprise_bits = (total_events > 0) ? (total_surprise / total_events) : 0.0;
        subs.site_perplexity = std::pow(2.0, subs.mean_surprise_bits);

        // Format summary string
        if (subs.events.empty()) {
            subs.summary_string = "Conserved";
        } else {
            std::ostringstream ss;
            bool first = true;
            for (const auto& ev : subs.events) {
                if (ev.is_synonymous) continue;
                if (!first) ss << ", ";
                ss << ev.parent_aa << " → " << ev.child_aa;
                if (ev.gc_change > 0) ss << " [GC↑]";
                else if (ev.gc_change < 0) ss << " [AT↑]";
                if (ev.surprise_bits >= 10.0) ss << "⚡";
                first = false;
            }
            if (subs.syn_count > 0) {
                ss << (first ? "" : " ") << "(+" << subs.syn_count << " syn)";
            }
            subs.summary_string = ss.str().empty() ? "Conserved (+syn)" : ss.str();
        }

        return subs;
    }

    // ========================================================================
    // 7. Full Execution Loop
    // ========================================================================

    std::vector<DirectionalSiteResult> run(bool show_progress = false, bool force_progress = false) {
        size_t num_patterns = aln.patterns.size();
        auto leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        std::vector<DirectionalSiteResult> pattern_results(num_patterns);

        SiteVisualizer viz(aln.num_codons, num_patterns, SiteVisualizerConfig{
            .method_name = "dFEL",
            .p_value_threshold = p_value_threshold,
            .pos_label = "positive",
            .neg_label = "negative",
            .sig_label = "significant",
            .quiet = !show_progress,
            .force_tty = force_progress,
            .enable_chart = true
        });

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            pattern_results[p] = analyze_pattern(p, leaf_to_taxon);
            pattern_results[p].subs = reconstruct_ancestral_pattern(p, leaf_to_taxon);

            const auto& r = pattern_results[p];
            SiteSelectionType type = SiteSelectionType::Neutral;
            Scalar metric = 0.0;
            std::string detail;
            if (r.p_value <= p_value_threshold) {
                if (r.beta > r.alpha) {
                    type = SiteSelectionType::Positive;
                    metric = r.beta / std::max(r.alpha, 1e-4);
                    detail = "▲ dN/dS = " + SiteVisualizer::format_scalar(metric, 2) + ", p = " + SiteVisualizer::format_pvalue(r.p_value);
                } else if (r.alpha > r.beta) {
                    type = SiteSelectionType::Negative;
                    metric = r.beta / std::max(r.alpha, 1e-4);
                    detail = "▼ dN/dS = " + SiteVisualizer::format_scalar(metric, 2) + ", p = " + SiteVisualizer::format_pvalue(r.p_value);
                }
            }
            viz.record_pattern(p, aln.patterns[p].sites, type, r.p_value, metric, detail);
        }

        viz.finish();

        site_results.resize(aln.num_codons);
        for (size_t s = 0; s < aln.num_codons; ++s) {
            size_t p_idx = aln.site_to_pattern[s];
            site_results[s] = pattern_results[p_idx];
            site_results[s].site_index = s;
            site_results[s].subs.site_index = s;
        }
        return site_results;
    }

    // ========================================================================
    // 8. JSON Export (Modern v3.0 & Legacy)
    // ========================================================================

nlohmann::json to_modern_json(const std::string& input_filepath = "", const std::string& tree_string = "") const;

nlohmann::json to_legacy_json(const std::string& input_filepath = "", const std::string& tree_string = "") const;

    void export_json(const std::string& filepath, JSONFormat format, const std::string& input_filepath = "", const std::string& tree_string = "") const {
        std::ofstream out(filepath);
        if (!out.is_open()) {
            throw std::runtime_error("Could not open output JSON file: " + filepath);
        }

        if (format == JSONFormat::Legacy) {
            out << to_legacy_json(input_filepath, tree_string).dump(2);
        } else {
            out << to_modern_json(input_filepath, tree_string).dump(2);
        }
    }

    DirectionalFELSummary get_summary() const {
        DirectionalFELSummary sum;
        sum.tested_sites = site_results.size();
        for (const auto& r : site_results) {
            if (r.p_value <= p_value_threshold) {
                if (r.beta > r.alpha) sum.positively_selected++;
                else if (r.alpha > r.beta) sum.negatively_selected++;
                else sum.neutral_sites++;
            } else {
                sum.neutral_sites++;
            }
        }
        return sum;
    }
};

} // namespace hyphy::analyses
