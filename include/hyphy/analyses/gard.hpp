#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/neighbor_joining.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/tree_engine.hpp"
#include "hyphy/analyses/ohmic.hpp"
#include "hyphy/core/wheatstone_likelihood.hpp"
#include "hyphy/analyses/gard_visualizer.hpp"
#include <nlohmann/json.hpp>

#include <vector>
#include <set>
#include <string>
#include <unordered_map>
#include <map>
#include <cmath>
#include <chrono>
#include <random>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <mutex>
#include <numeric>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;

enum class GARDTreeMethod {
    Ohmic,
    NeighborJoining
};

enum class GARDSRVMode {
    None,
    DiscreteGamma
};

enum class GARDTopologyTestMode {
    Both,
    SH,
    KH,
    None
};

struct BreakpointTopologyTest {
    size_t breakpoint_site = 0; // 1-based site position
    size_t left_partition_index = 0;
    size_t right_partition_index = 0;
    std::pair<size_t, size_t> left_interval;  // [start, end] 1-based
    std::pair<size_t, size_t> right_interval; // [start, end] 1-based

    Scalar left_tree_on_left_lnl = 0.0;
    Scalar right_tree_on_left_lnl = 0.0;
    Scalar delta_lnl_left = 0.0;

    Scalar right_tree_on_right_lnl = 0.0;
    Scalar left_tree_on_right_lnl = 0.0;
    Scalar delta_lnl_right = 0.0;

    Scalar kh_p_value_left = 1.0;
    Scalar kh_p_value_right = 1.0;

    Scalar sh_p_value_left = 1.0;
    Scalar sh_p_value_right = 1.0;

    bool statistically_significant = false;
};

struct PartitionTreeComparison {
    size_t partition_index = 0;
    std::pair<size_t, size_t> interval; // [start, end] 1-based
    size_t best_tree_index = 0;
    std::vector<Scalar> tree_log_likelihoods;
    std::vector<Scalar> kh_p_values; // vs partition best tree
    std::vector<Scalar> sh_p_values; // vs partition best tree
};

struct TopologyTestResult {
    bool run = false;
    std::string test_type = "both";
    size_t bootstrap_replicates = 1000;
    bool any_incongruence_significant = false;
    std::vector<BreakpointTopologyTest> breakpoint_tests;
    std::vector<PartitionTreeComparison> partition_matrix;
};

struct GARDMath {
    static inline double regularized_gamma_p(double a, double z) {
        return hyphy::core::stats::regularized_gamma_p(a, z);
    }
};

struct DiscreteGamma {
    Scalar alpha = 1.0;
    size_t num_classes = 4;
    std::vector<Scalar> rates;
    std::vector<Scalar> weights;

    static DiscreteGamma compute(Scalar alpha, size_t K = 4) {
        DiscreteGamma dg;
        dg.alpha = alpha;
        dg.num_classes = K;
        dg.weights.assign(K, 1.0 / static_cast<Scalar>(K));
        dg.rates.resize(K);

        if (K <= 1 || alpha <= 1e-4) {
            dg.rates.assign(K, 1.0);
            return dg;
        }

        std::vector<double> cutoffs(K + 1, 0.0);
        cutoffs[0] = 0.0;
        cutoffs[K] = 1e30;

        for (size_t k = 1; k < K; ++k) {
            double target_u = static_cast<double>(k) / static_cast<double>(K);
            double x_low = 1e-6;
            double x_high = std::max(5.0, 10.0 / alpha);
            while (GARDMath::regularized_gamma_p(alpha, alpha * x_high) < target_u && x_high < 1e6) {
                x_high *= 2.0;
            }

            double x_mid = 1.0;
            for (int iter = 0; iter < 40; ++iter) {
                x_mid = 0.5 * (x_low + x_high);
                double p = GARDMath::regularized_gamma_p(alpha, alpha * x_mid);
                if (std::abs(p - target_u) < 1e-10) break;
                if (p < target_u) {
                    x_low = x_mid;
                } else {
                    x_high = x_mid;
                }
            }
            cutoffs[k] = x_mid;
        }

        Scalar sum_rates = 0.0;
        for (size_t k = 0; k < K; ++k) {
            double z_lo = alpha * cutoffs[k];
            double z_hi = (k + 1 == K) ? 1e30 : (alpha * cutoffs[k + 1]);
            double p_hi = (k + 1 == K) ? 1.0 : GARDMath::regularized_gamma_p(alpha + 1.0, z_hi);
            double p_lo = (k == 0) ? 0.0 : GARDMath::regularized_gamma_p(alpha + 1.0, z_lo);
            double r = static_cast<double>(K) * (p_hi - p_lo);
            r = std::max(1e-4, r);
            dg.rates[k] = static_cast<Scalar>(r);
            sum_rates += dg.rates[k];
        }

        if (sum_rates > 0.0) {
            Scalar scale = static_cast<Scalar>(K) / sum_rates;
            for (size_t k = 0; k < K; ++k) {
                dg.rates[k] *= scale;
            }
        }

        return dg;
    }
};

struct GARDSettings {
    std::string input_file;
    std::string output_file;
    std::string run_mode = "Normal"; // "Normal" or "Faster"
    size_t max_breakpoints = 10000;
    size_t population_size = 32;
    Scalar mutation_rate = 0.35;               // Tuned: classic CHC 35% cataclysmic mutation rate (Eshelman 1991)
    Scalar rate_small_shifts = 0.80;            // 80% localized shifts, 20% global jumps
    Scalar caic_diversity_threshold = 0.25;    // Tuned: trigger cataclysm when fraction of unique models < 25%
    size_t max_stagnant_generations = 30;      // Tuned: 30 stagnant generations (Normal) vs 100 in 2005
    Scalar caic_improvement_threshold = 0.20;  // Tuned: 0.20 c-AIC meaningful improvement threshold
    bool quiet = false;
    bool dynamic_chart = true;
    bool force_tty = false;
    bool legacy_json = false;

    // Advanced tree search and reticulation settings
    GARDTreeMethod tree_method = GARDTreeMethod::Ohmic; // Fast O(N log N) Ohmic spectral trees
    GARDSRVMode srv_mode = GARDSRVMode::DiscreteGamma;  // Site-to-site rate variation: Discrete Gamma (4 classes)
    size_t gamma_rate_classes = 4;                      // Number of Gamma rate classes
    bool allow_bridges = true;                          // Wheatstone bridge screening & seeding

    // Post-hoc topological incongruence testing (KH & SH tests with RELL bootstrap)
    GARDTopologyTestMode topology_test_mode = GARDTopologyTestMode::Both;
    size_t topology_test_replicates = 1000;
};

struct GARDPartition {
    size_t start_site = 1; // 1-based start
    size_t end_site = 1;   // 1-based end
    core::Tree tree;
    Scalar log_likelihood = 0.0;
};

struct GARDModel {
    std::vector<size_t> breakpoints; // 0-based site coordinates
    Scalar log_likelihood = 0.0;
    Scalar aicc = 0.0;
    std::vector<GARDPartition> partitions;

    std::string to_key() const {
        if (breakpoints.empty()) return "{}";
        std::ostringstream ss;
        ss << "{";
        for (size_t i = 0; i < breakpoints.size(); ++i) {
            if (i > 0) ss << ", ";
            ss << breakpoints[i];
        }
        ss << "}";
        return ss.str();
    }

    static std::vector<size_t> parse_key(const std::string& key) {
        std::vector<size_t> res;
        if (key.size() <= 2) return res;
        std::string inner = key.substr(1, key.size() - 2);
        std::stringstream ss(inner);
        std::string token;
        while (std::getline(ss, token, ',')) {
            while (!token.empty() && std::isspace(static_cast<unsigned char>(token.front()))) token.erase(token.begin());
            while (!token.empty() && std::isspace(static_cast<unsigned char>(token.back()))) token.pop_back();
            if (!token.empty()) {
                res.push_back(std::stoull(token));
            }
        }
        std::sort(res.begin(), res.end());
        return res;
    }
};

struct GARDResult {
    Scalar baseline_aicc = 0.0;
    Scalar baseline_log_l = 0.0;
    core::Tree baseline_tree;
    core::GTRParameters baseline_rates;
    Vector4 baseline_frequencies = Vector4::Zero();

    // SRV
    GARDSRVMode srv_mode = GARDSRVMode::DiscreteGamma;
    Scalar gamma_alpha = 0.0;
    std::vector<Scalar> gamma_rates;
    std::vector<Scalar> gamma_weights;

    // Tree Engine
    GARDTreeMethod tree_method = GARDTreeMethod::Ohmic;

    // Ohmic & Wheatstone Bridges
    bool ohmic_bridges_detected = false;
    size_t num_bridges = 0;
    Scalar ohmic_log_l = 0.0;
    Scalar ohmic_aicc = 0.0;
    Scalar ohmic_delta_aicc = 0.0;
    std::vector<OhmicBridgeSummary> bridge_summaries;
    std::vector<size_t> bridge_candidate_breakpoints;
    std::vector<Scalar> site_bridge_delta_log_l;

    size_t num_sites = 0;
    size_t num_taxa = 0;
    size_t potential_breakpoints = 0;
    size_t min_partition_size = 0;
    std::vector<size_t> variable_sites;

    GARDModel best_model;
    Scalar single_tree_aicc = 0.0;
    bool topological_incongruence = false;
    TopologyTestResult topology_test_result;

    std::map<std::string, Scalar> site_breakpoint_support;
    std::vector<std::pair<size_t, Scalar>> improvements; // (num_bps, deltaAICc)
    std::vector<std::vector<size_t>> improvement_breakpoints; // per-step breakpoints
    size_t total_models_evaluated = 0;
    double time_elapsed = 0.0;
    std::string input_file;

    nlohmann::json to_json() const;
    nlohmann::json to_legacy_json() const;
};

/**
 * @brief Post-hoc topological incongruence testing for GARD.
 * 
 * Implements:
 * 1. Kishino-Hasegawa (KH 1989) test with asymptotic normal and RELL (Resampling of Estimated Log-Likelihoods) bootstrap.
 * 2. Shimodaira-Hasegawa (SH 1999) test with centered RELL bootstrap to properly correct for a posteriori tree selection.
 */
class GARDTopologyTester {
public:
    struct KHResult {
        Scalar delta_lnl = 0.0; // L1 - L2
        Scalar std_error = 0.0;
        Scalar p_asymptotic = 1.0;
        Scalar p_rell = 1.0;
    };

    static KHResult run_kh_test(
        const std::vector<Scalar>& lnl1,
        const std::vector<Scalar>& lnl2,
        size_t replicates = 1000,
        uint64_t seed = 42
    ) {
        size_t N = lnl1.size();
        KHResult res;
        if (N == 0 || lnl1.size() != lnl2.size()) return res;

        std::vector<Scalar> d(N);
        Scalar sum_d = 0.0;
        for (size_t i = 0; i < N; ++i) {
            d[i] = lnl1[i] - lnl2[i];
            sum_d += d[i];
        }
        res.delta_lnl = sum_d;
        Scalar mean_d = sum_d / static_cast<Scalar>(N);

        Scalar var_d = 0.0;
        for (size_t i = 0; i < N; ++i) {
            Scalar diff = d[i] - mean_d;
            var_d += diff * diff;
        }
        var_d = (N > 1) ? (var_d / static_cast<Scalar>(N - 1)) : 0.0;
        res.std_error = std::sqrt(static_cast<Scalar>(N) * var_d);

        // Asymptotic normal p-value (two-sided)
        if (res.std_error > 1e-12) {
            Scalar z = std::abs(res.delta_lnl) / res.std_error;
            res.p_asymptotic = std::erfc(z / std::sqrt(2.0));
        } else {
            res.p_asymptotic = (std::abs(res.delta_lnl) < 1e-9) ? 1.0 : 0.0;
        }

        // RELL bootstrap (two-sided)
        if (replicates > 0 && res.std_error > 1e-12) {
            std::vector<Scalar> c(N);
            for (size_t i = 0; i < N; ++i) {
                c[i] = d[i] - mean_d;
            }

            std::mt19937_64 rng(seed);
            std::uniform_int_distribution<size_t> dist(0, N - 1);

            size_t count_extreme = 0;
            Scalar abs_obs = std::abs(res.delta_lnl);

            for (size_t r = 0; r < replicates; ++r) {
                Scalar b_sum = 0.0;
                for (size_t j = 0; j < N; ++j) {
                    b_sum += c[dist(rng)];
                }
                if (std::abs(b_sum) >= abs_obs - 1e-9) {
                    count_extreme++;
                }
            }
            res.p_rell = static_cast<Scalar>(count_extreme) / static_cast<Scalar>(replicates);
        } else {
            res.p_rell = res.p_asymptotic;
        }

        return res;
    }

    struct SHResult {
        size_t best_index = 0;
        std::vector<Scalar> total_lnl;
        std::vector<Scalar> delta_lnl; // best - m >= 0
        std::vector<Scalar> p_values;  // SH p-values
    };

    static SHResult run_sh_test(
        const std::vector<std::vector<Scalar>>& site_lnl_matrix, // M models x N sites
        size_t replicates = 1000,
        uint64_t seed = 42
    ) {
        SHResult res;
        size_t M = site_lnl_matrix.size();
        if (M == 0) return res;
        size_t N = site_lnl_matrix[0].size();
        if (N == 0) return res;

        res.total_lnl.resize(M, 0.0);
        res.delta_lnl.resize(M, 0.0);
        res.p_values.resize(M, 1.0);

        size_t best_m = 0;
        Scalar max_l = -1e300;
        for (size_t m = 0; m < M; ++m) {
            Scalar sum = 0.0;
            for (size_t i = 0; i < N; ++i) {
                sum += site_lnl_matrix[m][i];
            }
            res.total_lnl[m] = sum;
            if (sum > max_l) {
                max_l = sum;
                best_m = m;
            }
        }
        res.best_index = best_m;

        for (size_t m = 0; m < M; ++m) {
            res.delta_lnl[m] = max_l - res.total_lnl[m];
        }

        if (M <= 1 || replicates == 0) {
            return res;
        }

        // Center site log-likelihoods for each model: c_{m, i} = v_{m, i} - mean_m
        std::vector<std::vector<Scalar>> centered(M, std::vector<Scalar>(N, 0.0));
        for (size_t m = 0; m < M; ++m) {
            Scalar mean_m = res.total_lnl[m] / static_cast<Scalar>(N);
            for (size_t i = 0; i < N; ++i) {
                centered[m][i] = site_lnl_matrix[m][i] - mean_m;
            }
        }

        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<size_t> dist(0, N - 1);

        std::vector<size_t> counts(M, 0);
        std::vector<Scalar> rep_sums(M, 0.0);

        for (size_t r = 0; r < replicates; ++r) {
            std::fill(rep_sums.begin(), rep_sums.end(), 0.0);
            for (size_t j = 0; j < N; ++j) {
                size_t idx = dist(rng);
                for (size_t m = 0; m < M; ++m) {
                    rep_sums[m] += centered[m][idx];
                }
            }

            Scalar rep_max = rep_sums[0];
            for (size_t m = 1; m < M; ++m) {
                if (rep_sums[m] > rep_max) rep_max = rep_sums[m];
            }

            for (size_t m = 0; m < M; ++m) {
                Scalar rep_delta = rep_max - rep_sums[m];
                if (rep_delta >= res.delta_lnl[m] - 1e-9) {
                    counts[m]++;
                }
            }
        }

        for (size_t m = 0; m < M; ++m) {
            if (m == best_m) {
                res.p_values[m] = 1.0;
            } else {
                res.p_values[m] = static_cast<Scalar>(counts[m]) / static_cast<Scalar>(replicates);
            }
        }

        return res;
    }
};

class GARDAnalyzer {
public:
    GARDSettings settings;
    const core::Alignment& alignment;

    struct PairHash {
        size_t operator()(const std::pair<size_t, size_t>& p) const {
            return std::hash<size_t>()(p.first) ^ (std::hash<size_t>()(p.second) << 16);
        }
    };

    struct PartitionFit {
        core::Tree tree;
        Scalar log_likelihood = 0.0;
    };

    static Scalar compute_gtr_srv_log_likelihood(
        const Tree& tree,
        const Alignment& aln,
        const GTRMatrix& gtr_model,
        const std::vector<Scalar>& srv_rates,
        const std::vector<Scalar>& srv_weights,
        const std::vector<size_t>* custom_leaf_to_taxon = nullptr,
        std::vector<Scalar>* out_site_log_l = nullptr
    ) {
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        size_t K_rates = srv_rates.empty() ? 1 : srv_rates.size();
        if (num_nodes == 0 || num_patterns == 0) return 0.0;

        std::vector<Scalar> effective_rates = srv_rates.empty() ? std::vector<Scalar>{1.0} : srv_rates;
        std::vector<Scalar> effective_weights = srv_weights.empty() ? std::vector<Scalar>{1.0} : srv_weights;

        std::vector<std::vector<Matrix4>> P_branches(K_rates, std::vector<Matrix4>(num_nodes, Matrix4::Identity()));
        for (size_t k = 0; k < K_rates; ++k) {
            Scalar r_k = effective_rates[k];
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    P_branches[k][node.id] = gtr_model.transition_matrix(node.branch_length * r_k);
                }
            }
        }

        std::vector<size_t> local_map;
        const std::vector<size_t>& leaf_to_taxon = custom_leaf_to_taxon
            ? *custom_leaf_to_taxon
            : (local_map = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln));

        std::vector<Scalar> pattern_log_l(num_patterns, 0.0);
        Scalar total_log_l = 0.0;

        auto eval_pattern = [&](size_t p, std::vector<Vector4>& node_L, std::vector<double>& node_log_scale) -> Scalar {
            const auto& pattern = aln.patterns[p];

            // Initialize leaf states once for this pattern
            for (int32_t node_id : tree.post_order) {
                const auto& node = tree.nodes[node_id];
                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[node_id];
                    if (t_idx != static_cast<size_t>(-1)) {
                        int8_t state = pattern.states[t_idx];
                        if (state >= 0 && state < 4) {
                            node_L[node_id].setZero();
                            node_L[node_id](state) = 1.0;
                        } else {
                            node_L[node_id].setOnes();
                        }
                    } else {
                        node_L[node_id].setOnes();
                    }
                    node_log_scale[node_id] = 0.0;
                }
            }

            Scalar max_log = -1e300;
            std::vector<Scalar> log_terms(K_rates, -1e300);

            for (size_t k = 0; k < K_rates; ++k) {
                Scalar w_k = effective_weights[k];
                const auto& P_k = P_branches[k];

                for (int32_t node_id : tree.post_order) {
                    const auto& node = tree.nodes[node_id];
                    if (!node.is_leaf) {
                        double child_scale_sum = 0.0;
                        tree.eval_internal_node_with_model(node_id, P_k, gtr_model, node_L, node_log_scale, node_L[node_id], child_scale_sum, effective_rates[k]);
                        node_log_scale[node_id] = child_scale_sum;
                    }
                }

                Scalar cond_l = gtr_model.pi.dot(node_L[tree.root_id]);
                if (cond_l > 0.0) {
                    log_terms[k] = std::log(w_k) + std::log(cond_l) + node_log_scale[tree.root_id];
                    if (log_terms[k] > max_log) max_log = log_terms[k];
                }
            }

            Scalar pat_lnl = -690.7755;
            if (max_log > -1e200) {
                Scalar sum_exp = 0.0;
                for (size_t k = 0; k < K_rates; ++k) {
                    if (log_terms[k] > -1e200) {
                        sum_exp += std::exp(log_terms[k] - max_log);
                    }
                }
                if (sum_exp > 0.0) {
                    pat_lnl = max_log + std::log(sum_exp);
                }
            }
            return pat_lnl;
        };

        if (omp_in_parallel()) {
            std::vector<Vector4> node_L(num_nodes);
            std::vector<double> node_log_scale(num_nodes, 0.0);
            for (size_t p = 0; p < num_patterns; ++p) {
                pattern_log_l[p] = eval_pattern(p, node_L, node_log_scale);
            }
        } else {
            #pragma omp parallel
            {
                std::vector<Vector4> node_L(num_nodes);
                std::vector<double> node_log_scale(num_nodes, 0.0);
                #pragma omp for schedule(static)
                for (size_t p = 0; p < num_patterns; ++p) {
                    pattern_log_l[p] = eval_pattern(p, node_L, node_log_scale);
                }
            }
        }

        Scalar total_comp = 0.0;
        for (size_t p = 0; p < num_patterns; ++p) {
            LikelihoodEngine::neumaier_acc(total_log_l, total_comp, aln.patterns[p].weight * pattern_log_l[p]);
        }
        total_log_l += total_comp;

        if (out_site_log_l != nullptr) {
            size_t S = aln.site_to_pattern.size();
            out_site_log_l->resize(S, 0.0);
            for (size_t s = 0; s < S; ++s) {
                (*out_site_log_l)[s] = pattern_log_l[aln.site_to_pattern[s]];
            }
        }

        return total_log_l;
    }

    static void sweep_gtr_srv_branches_cavity(
        Tree& curr_tree,
        const Alignment& nuc_aln,
        const GTRMatrix& gtr,
        const std::vector<Scalar>& srv_rates,
        const std::vector<Scalar>& srv_weights,
        const std::vector<size_t>& leaf_to_taxon,
        Scalar tol = 1e-3,
        int max_iter = 10
    ) {
        size_t num_nodes = curr_tree.num_nodes();
        size_t num_patterns = nuc_aln.patterns.size();
        size_t K_rates = srv_rates.empty() ? 1 : srv_rates.size();
        if (num_nodes == 0 || num_patterns == 0) return;

        std::vector<Scalar> effective_rates = srv_rates.empty() ? std::vector<Scalar>{1.0} : srv_rates;
        std::vector<Scalar> effective_weights = srv_weights.empty() ? std::vector<Scalar>{1.0} : srv_weights;

        std::vector<std::vector<std::vector<Vector4>>> D(
            K_rates, std::vector<std::vector<Vector4>>(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones())));
        std::vector<std::vector<std::vector<Vector4>>> P_D(
            K_rates, std::vector<std::vector<Vector4>>(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones())));
        std::vector<std::vector<std::vector<Vector4>>> U(
            K_rates, std::vector<std::vector<Vector4>>(num_nodes, std::vector<Vector4>(num_patterns, Vector4::Ones())));
        std::vector<std::vector<std::vector<double>>> D_log_scale(
            K_rates, std::vector<std::vector<double>>(num_nodes, std::vector<double>(num_patterns, 0.0)));
        std::vector<std::vector<std::vector<double>>> U_log_scale(
            K_rates, std::vector<std::vector<double>>(num_nodes, std::vector<double>(num_patterns, 0.0)));
        std::vector<std::vector<Matrix4>> P_branch(
            K_rates, std::vector<Matrix4>(num_nodes, Matrix4::Identity()));

        for (size_t k = 0; k < K_rates; ++k) {
            Scalar r_k = effective_rates[k];
            for (const auto& node : curr_tree.nodes) {
                if (node.id != curr_tree.root_id && (node.is_leaf || !node.children.empty())) {
                    P_branch[k][node.id] = gtr.transition_matrix(node.branch_length * r_k);
                }
            }
        }

        for (size_t k = 0; k < K_rates; ++k) {
            for (int32_t nid : curr_tree.post_order) {
                const auto& node = curr_tree.nodes[nid];
                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[nid];
                    for (size_t p = 0; p < num_patterns; ++p) {
                        int8_t state = (t_idx != static_cast<size_t>(-1)) ? nuc_aln.patterns[p].states[t_idx] : -1;
                        if (state >= 0 && state < 4) {
                            D[k][nid][p].setZero();
                            D[k][nid][p](state) = 1.0;
                        } else {
                            D[k][nid][p].setOnes();
                        }
                        D_log_scale[k][nid][p] = 0.0;
                        if (nid != curr_tree.root_id) {
                            P_D[k][nid][p] = P_branch[k][nid] * D[k][nid][p];
                        }
                    }
                } else {
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Vector4 d_vec = Vector4::Ones();
                        double sc = 0.0;
                        for (int32_t cid : node.children) {
                            d_vec.array() *= P_D[k][cid][p].array();
                            sc += D_log_scale[k][cid][p];
                            Scalar mx = d_vec.maxCoeff();
                            if (mx > 0.0 && mx < 1e-100) {
                                d_vec /= mx;
                                sc += std::log(mx);
                            }
                        }
                        D[k][nid][p] = d_vec;
                        D_log_scale[k][nid][p] = sc;
                        if (nid != curr_tree.root_id) {
                            P_D[k][nid][p] = P_branch[k][nid] * d_vec;
                        }
                    }
                }
            }
            for (size_t p = 0; p < num_patterns; ++p) {
                U[k][curr_tree.root_id][p] = gtr.pi;
                U_log_scale[k][curr_tree.root_id][p] = 0.0;
            }
        }

        std::vector<std::vector<Vector4>> c_spec(K_rates, std::vector<Vector4>(num_patterns));
        std::vector<std::vector<double>> c_log_scale(K_rates, std::vector<double>(num_patterns, 0.0));
        std::vector<double> max_c_scale(num_patterns, 0.0);

        std::function<void(int32_t)> dfs = [&](int32_t uid) {
            const auto& node = curr_tree.nodes[uid];
            if (uid != curr_tree.root_id && (node.is_leaf || !node.children.empty())) {
                for (size_t p = 0; p < num_patterns; ++p) {
                    double mx_sc = -1e300;
                    for (size_t k = 0; k < K_rates; ++k) {
                        Vector4 u_hat = gtr.U_pi.transpose() * U[k][uid][p];
                        Vector4 d_hat = gtr.V_pi * D[k][uid][p];
                        c_spec[k][p] = u_hat.cwiseProduct(d_hat);
                        c_log_scale[k][p] = U_log_scale[k][uid][p] + D_log_scale[k][uid][p];
                        if (c_log_scale[k][p] > mx_sc) mx_sc = c_log_scale[k][p];
                    }
                    max_c_scale[p] = mx_sc;
                }
                Scalar init_bl = curr_tree.nodes[uid].branch_length;
                auto bl_obj = [&](Scalar bl) -> Scalar {
                    std::vector<Vector4> exp_lam(K_rates);
                    for (size_t k = 0; k < K_rates; ++k) {
                        exp_lam[k] = (gtr.eigenvalues * (bl * effective_rates[k])).array().exp();
                    }
                    Scalar lnl = 0.0;
                    for (size_t p = 0; p < num_patterns; ++p) {
                        double ref_sc = max_c_scale[p];
                        Scalar pat_lk = 0.0;
                        for (size_t k = 0; k < K_rates; ++k) {
                            Scalar cond_lk = c_spec[k][p].dot(exp_lam[k]);
                            if (cond_lk > 0.0) {
                                pat_lk += effective_weights[k] * cond_lk * std::exp(c_log_scale[k][p] - ref_sc);
                            }
                        }
                        lnl += nuc_aln.patterns[p].weight * ((pat_lk > 1e-300) ? (ref_sc + std::log(pat_lk)) : -1e20);
                    }
                    return -lnl;
                };

                auto res_bl = opt::Brent1D::minimize(bl_obj, 1e-6, std::max<Scalar>(init_bl, 1e-4), 5.0, tol, max_iter);
                Scalar best_bl = res_bl.first;
                curr_tree.nodes[uid].branch_length = best_bl;
                for (size_t k = 0; k < K_rates; ++k) {
                    P_branch[k][uid] = gtr.transition_matrix(best_bl * effective_rates[k]);
                    for (size_t p = 0; p < num_patterns; ++p) {
                        P_D[k][uid][p] = P_branch[k][uid] * D[k][uid][p];
                    }
                }
            }

            if (!node.is_leaf && !node.children.empty()) {
                size_t n_ch = node.children.size();
                for (size_t i = 0; i < n_ch; ++i) {
                    int32_t cid = node.children[i];
                    for (size_t k = 0; k < K_rates; ++k) {
                        for (size_t p = 0; p < num_patterns; ++p) {
                            Vector4 u_c = (uid == curr_tree.root_id)
                                ? gtr.pi
                                : (P_branch[k][uid].transpose() * U[k][uid][p]);
                            double sc = (uid == curr_tree.root_id) ? 0.0 : U_log_scale[k][uid][p];
                            for (size_t j = 0; j < n_ch; ++j) {
                                if (j != i) {
                                    int32_t sib = node.children[j];
                                    u_c.array() *= P_D[k][sib][p].array();
                                    sc += D_log_scale[k][sib][p];
                                    Scalar mx_j = u_c.maxCoeff();
                                    if (mx_j > 0.0 && mx_j < 1e-100) {
                                        u_c /= mx_j;
                                        sc += std::log(mx_j);
                                    }
                                }
                            }
                            Scalar mx = u_c.maxCoeff();
                            if (mx > 0.0 && mx < 1e-100) {
                                u_c /= mx;
                                sc += std::log(mx);
                            }
                            U[k][cid][p] = u_c;
                            U_log_scale[k][cid][p] = sc;
                        }
                    }
                    dfs(cid);
                }
                for (size_t k = 0; k < K_rates; ++k) {
                    for (size_t p = 0; p < num_patterns; ++p) {
                        Vector4 d_vec = Vector4::Ones();
                        double sc = 0.0;
                        for (int32_t cid : node.children) {
                            d_vec.array() *= P_D[k][cid][p].array();
                            sc += D_log_scale[k][cid][p];
                            Scalar mx_c = d_vec.maxCoeff();
                            if (mx_c > 0.0 && mx_c < 1e-100) {
                                d_vec /= mx_c;
                                sc += std::log(mx_c);
                            }
                        }
                        D[k][uid][p] = d_vec;
                        D_log_scale[k][uid][p] = sc;
                        if (uid != curr_tree.root_id) {
                            P_D[k][uid][p] = P_branch[k][uid] * d_vec;
                        }
                    }
                }
            }
        };

        dfs(curr_tree.root_id);
    }

    static std::pair<Tree, Scalar> fit_gtr_segment_srv(
        Tree tree,
        const Alignment& slice,
        const GTRParameters& params,
        const Vector4& freqs,
        const std::vector<Scalar>& srv_rates,
        const std::vector<Scalar>& srv_weights,
        int passes = 1
    ) {
        if (srv_rates.size() <= 1 && tree.bridges.empty()) {
            return opt::GTRBranchLengthFitter::fit(std::move(tree), slice, params, freqs, passes);
        }

        std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, slice);
        GTRMatrix gtr;
        gtr.update(params, freqs);

        // Stage 1: Fast scaffold using unstratified cavity engine (O(N), ~0.1ms)
        auto fit_scaffold = opt::GTRBranchLengthFitter::fit(std::move(tree), slice, params, freqs, 1);
        Tree curr_tree = std::move(fit_scaffold.first);

        std::vector<int32_t> branch_node_ids;
        for (const auto& node : curr_tree.nodes) {
            if (node.id != curr_tree.root_id) {
                branch_node_ids.push_back(node.id);
            }
        }
        size_t num_branches = branch_node_ids.size();

        // Stage 2: Rescale tree under the discrete Gamma SRV mixture (1D Brent, ~10 evals)
        auto scale_obj = [&](Scalar s) -> Scalar {
            Tree t_scale = curr_tree;
            for (size_t b = 0; b < num_branches; ++b) {
                t_scale.nodes[branch_node_ids[b]].branch_length *= s;
            }
            t_scale.scale_bridge_lengths(s);
            return -compute_gtr_srv_log_likelihood(t_scale, slice, gtr, srv_rates, srv_weights, &leaf_to_taxon);
        };
        auto res_s = opt::Brent1D::minimize(scale_obj, 0.05, 1.0, 20.0, 1e-3, 10);
        Scalar opt_s = res_s.first;
        for (size_t b = 0; b < num_branches; ++b) {
            curr_tree.nodes[branch_node_ids[b]].branch_length *= opt_s;
        }
        curr_tree.scale_bridge_lengths(opt_s);

        // Stage 3: Closed-form spectral cavity sweep under SRV (if tree has no bridges and passes >= 1)
        if (curr_tree.bridges.empty()) {
            for (int p = 0; p < passes; ++p) {
                sweep_gtr_srv_branches_cavity(curr_tree, slice, gtr, srv_rates, srv_weights, leaf_to_taxon, 1e-3, 10);
            }
        } else {
            for (auto& br : curr_tree.bridges) {
                if (br.bridge_node_id >= 0 && static_cast<size_t>(br.bridge_node_id) < curr_tree.nodes.size() &&
                    br.bridge_node_id != br.hub_node_id) {
                    br.branch_lengths[0] = curr_tree.nodes[br.bridge_node_id].branch_length;
                }
            }
        }

        Scalar final_lnl = compute_gtr_srv_log_likelihood(curr_tree, slice, gtr, srv_rates, srv_weights, &leaf_to_taxon);
        return {curr_tree, final_lnl};
    }

    GARDAnalyzer(const core::Alignment& aln, GARDSettings s = GARDSettings())
        : settings(s), alignment(aln) {
        if (settings.run_mode == "Faster") {
            settings.max_stagnant_generations = std::min(settings.population_size / 2, static_cast<size_t>(15));
            settings.caic_improvement_threshold = 1.5;
        }
    }

    GARDResult run() {
        auto start_time = std::chrono::high_resolution_clock::now();
        GARDResult result;
        result.input_file = settings.input_file;
        result.num_sites = alignment.num_nucleotides;
        result.num_taxa = alignment.num_taxa;
        result.tree_method = settings.tree_method;
        result.srv_mode = settings.srv_mode;

        if (result.num_taxa < 3) {
            throw std::runtime_error("GARD requires at least 3 sequences for phylogenetic inference");
        }

        // 1. Identify variable sites
        result.variable_sites.clear();
        for (size_t site = 0; site < result.num_sites; ++site) {
            int8_t first_state = -1;
            bool is_variable = false;
            for (size_t t = 0; t < result.num_taxa; ++t) {
                int8_t s = core::char_to_nuc(alignment.sequences[t].raw_seq[site]);
                if (s >= 0 && s < 4) {
                    if (first_state == -1) {
                        first_state = s;
                    } else if (s != first_state) {
                        is_variable = true;
                        break;
                    }
                }
            }
            if (is_variable) {
                result.variable_sites.push_back(site);
            }
        }

        result.potential_breakpoints = result.variable_sites.empty() ? 0 : (result.variable_sites.size() - 1);
        result.min_partition_size = 2 * result.num_taxa - 3;

        size_t min_expected_sites = 5 + 2 * result.min_partition_size + 2;
        if (result.num_sites < min_expected_sites) {
            throw std::runtime_error("Alignment is too short to permit c-AIC based model comparison. Need at least " +
                                     std::to_string(min_expected_sites) + " sites for " +
                                     std::to_string(result.num_taxa) + " sequences.");
        }

        auto validate_breakpoints = [&](const std::vector<size_t>& bps) -> bool {
            size_t last_bp = 0;
            for (size_t bp : bps) {
                if (bp < last_bp || (bp - last_bp + 1) < result.min_partition_size) {
                    return false;
                }
                last_bp = bp + 1;
            }
            if (result.num_sites < last_bp || (result.num_sites - last_bp) < result.min_partition_size) {
                return false;
            }
            return true;
        };

        if (!settings.quiet) {
            std::cout << "[GARD] Loaded alignment with " << result.num_taxa << " sequences, "
                      << result.num_sites << " sites (" << result.variable_sites.size()
                      << " variable sites, " << result.potential_breakpoints << " potential breakpoints).\n";
            std::cout << "[GARD] Minimum partition size: " << result.min_partition_size << " sites.\n";
            std::cout << "[GARD] Tree Search Engine: "
                      << ((settings.tree_method == GARDTreeMethod::Ohmic) ? "Ohmic Spectral Scaffolding (O(N log N))" : "Neighbor-Joining (O(N^3))") << "\n";
            std::cout << "[GARD] Rate Variation: "
                      << ((settings.srv_mode == GARDSRVMode::DiscreteGamma) ? ("Discrete Gamma (" + std::to_string(settings.gamma_rate_classes) + " rate classes)") : "None (Constant)") << "\n";
        }

        // Wheatstone bridge screening (linear-time, no transformer)
        if (settings.allow_bridges) {
            if (!settings.quiet) {
                std::cout << "[GARD] Screening for Ohmic Wheatstone bridges (linear-time, no transformer)...\n";
            }
            try {
                OhmicAnalysisSettings ohm_settings;
                ohm_settings.network_settings.use_hyphaeon_transformer = false;
                ohm_settings.quiet = true;
                OhmicAnalyzer ohm_analyzer(alignment, ohm_settings);
                OhmicResult ohm_res = ohm_analyzer.run();

                result.num_bridges = ohm_res.num_bridges;
                result.ohmic_bridges_detected = (ohm_res.num_bridges > 0);
                result.ohmic_log_l = ohm_res.ohmic_log_l;
                result.ohmic_aicc = ohm_res.ohmic_aicc;
                result.bridge_summaries = ohm_res.bridge_summaries;
                result.site_bridge_delta_log_l = ohm_res.site_delta_log_l;

                if (result.ohmic_bridges_detected && !ohm_res.site_delta_log_l.empty()) {
                    size_t num_s = ohm_res.site_delta_log_l.size();
                    int win = std::max(5, static_cast<int>(result.min_partition_size / 2));
                    std::vector<Scalar> smoothed(num_s, 0.0);
                    for (int s = 0; s < static_cast<int>(num_s); ++s) {
                        int w_start = std::max(0, s - win);
                        int w_end = std::min(static_cast<int>(num_s), s + win + 1);
                        Scalar sum = 0.0;
                        for (int w = w_start; w < w_end; ++w) {
                            sum += ohm_res.site_delta_log_l[w];
                        }
                        smoothed[s] = sum / static_cast<Scalar>(w_end - w_start);
                    }

                    std::vector<std::pair<Scalar, size_t>> peak_cands;
                    for (size_t s = 1; s + 1 < num_s; ++s) {
                        if (smoothed[s] > 0.05 && smoothed[s] >= smoothed[s - 1] && smoothed[s] >= smoothed[s + 1]) {
                            peak_cands.push_back({smoothed[s], s});
                        }
                    }
                    std::sort(peak_cands.begin(), peak_cands.end(), [](const auto& a, const auto& b) {
                        return a.first > b.first;
                    });

                    for (const auto& [score, s] : peak_cands) {
                        auto it = std::lower_bound(result.variable_sites.begin(), result.variable_sites.end(), s);
                        size_t cand_bp = s;
                        if (it != result.variable_sites.end()) {
                            cand_bp = *it;
                        }
                        if (validate_breakpoints({cand_bp})) {
                            bool too_close = false;
                            for (size_t existing : result.bridge_candidate_breakpoints) {
                                if (cand_bp >= existing && (cand_bp - existing) < result.min_partition_size) {
                                    too_close = true; break;
                                }
                                if (existing >= cand_bp && (existing - cand_bp) < result.min_partition_size) {
                                    too_close = true; break;
                                }
                            }
                            if (!too_close) {
                                result.bridge_candidate_breakpoints.push_back(cand_bp);
                                if (result.bridge_candidate_breakpoints.size() >= 5) break;
                            }
                        }
                    }

                    if (!settings.quiet && !result.bridge_candidate_breakpoints.empty()) {
                        std::cout << "[GARD] Wheatstone bridge active: " << result.num_bridges
                                  << " bridge(s), Ohmic c-AIC = " << std::fixed << std::setprecision(2)
                                  << result.ohmic_aicc << ".\n";
                        std::cout << "[GARD] Bridge trajectory seeded "
                                  << result.bridge_candidate_breakpoints.size()
                                  << " candidate breakpoint(s): ";
                        for (size_t i = 0; i < result.bridge_candidate_breakpoints.size(); ++i) {
                            if (i > 0) std::cout << ", ";
                            std::cout << (result.bridge_candidate_breakpoints[i] + 1);
                        }
                        std::cout << "\n";
                    }
                }
            } catch (const std::exception& e) {
                if (!settings.quiet) {
                    std::cout << "[GARD] Note: Bridge screening skipped: " << e.what() << "\n";
                }
            }
        }

        // 2. Baseline fit on entire alignment
        if (!settings.quiet) {
            std::cout << "[GARD] Fitting baseline (single-partition; no breakpoints) model...\n";
        }

        core::Tree baseline_tree = (settings.tree_method == GARDTreeMethod::Ohmic && result.num_taxa > 50)
            ? TreeEngine::build_spectral_scaffold(alignment)
            : core::NeighborJoining::infer(alignment);
        opt::GTRFitter baseline_fitter(baseline_tree, alignment);
        auto base_fit = baseline_fitter.fit();

        result.baseline_tree = base_fit.tree;
        result.baseline_rates = base_fit.params;
        result.baseline_frequencies = alignment.nuc_frequencies;
        result.baseline_log_l = base_fit.log_likelihood;

        // Optimize discrete Gamma rate variation if enabled
        if (settings.srv_mode == GARDSRVMode::DiscreteGamma) {
            if (!settings.quiet) {
                std::cout << "[GARD] Optimizing discrete Gamma rate variation ("
                          << settings.gamma_rate_classes << " rate classes)...\n";
            }
            GTRMatrix base_gtr;
            base_gtr.update(result.baseline_rates, result.baseline_frequencies);
            auto alpha_obj = [&](Scalar a) -> Scalar {
                DiscreteGamma dg = DiscreteGamma::compute(a, settings.gamma_rate_classes);
                return -compute_gtr_srv_log_likelihood(
                    result.baseline_tree, alignment, base_gtr, dg.rates, dg.weights
                );
            };
            auto [opt_a, _] = opt::Brent1D::minimize(alpha_obj, 0.05, 1.0, 50.0, 1e-3, 20);
            DiscreteGamma best_dg = DiscreteGamma::compute(opt_a, settings.gamma_rate_classes);
            result.gamma_alpha = opt_a;
            result.gamma_rates = best_dg.rates;
            result.gamma_weights = best_dg.weights;
            result.baseline_log_l = compute_gtr_srv_log_likelihood(
                result.baseline_tree, alignment, base_gtr, best_dg.rates, best_dg.weights
            );
            if (!settings.quiet) {
                std::cout << "[GARD] Estimated Gamma shape alpha = " << std::fixed << std::setprecision(4)
                          << result.gamma_alpha << " (rates: ";
                for (size_t k = 0; k < result.gamma_rates.size(); ++k) {
                    if (k > 0) std::cout << ", ";
                    std::cout << std::fixed << std::setprecision(3) << result.gamma_rates[k];
                }
                std::cout << ")\n";
            }
        } else {
            result.gamma_alpha = 0.0;
            result.gamma_rates = {1.0};
            result.gamma_weights = {1.0};
        }

        size_t K_srv = (settings.srv_mode == GARDSRVMode::DiscreteGamma) ? 1 : 0;
        size_t K_base = 5 + result.min_partition_size + K_srv;
        result.baseline_aicc = -2.0 * result.baseline_log_l +
            (2.0 * result.num_sites * K_base) / (result.num_sites - K_base - 1.0);

        if (result.ohmic_bridges_detected) {
            result.ohmic_delta_aicc = result.baseline_aicc - result.ohmic_aicc;
        }

        if (!settings.quiet) {
            std::cout << "[GARD] Baseline fit: Log-L = " << std::fixed << std::setprecision(4)
                      << result.baseline_log_l << ", c-AIC = " << result.baseline_aicc << "\n";
        }

        // Set baseline as initial best model
        GARDModel best_overall_model;
        best_overall_model.breakpoints.clear();
        best_overall_model.log_likelihood = result.baseline_log_l;
        best_overall_model.aicc = result.baseline_aicc;
        GARDPartition base_part;
        base_part.start_site = 1;
        base_part.end_site = result.num_sites;
        base_part.tree = result.baseline_tree;
        base_part.log_likelihood = result.baseline_log_l;
        best_overall_model.partitions.push_back(std::move(base_part));

        result.best_model = best_overall_model;

        // Thread-safe partition cache
        std::unordered_map<std::pair<size_t, size_t>, PartitionFit, PairHash> partition_cache;
        std::mutex cache_mutex;

        // Pre-populate whole alignment in cache
        partition_cache[{0, result.num_sites}] = {result.baseline_tree, result.baseline_log_l};

        auto fit_segment = [&](size_t start, size_t end) -> PartitionFit {
            std::pair<size_t, size_t> key = {start, end};
            {
                std::lock_guard<std::mutex> lock(cache_mutex);
                auto it = partition_cache.find(key);
                if (it != partition_cache.end()) {
                    return it->second;
                }
            }

            core::Alignment slice = alignment.slice_nucleotide(start, end);
            core::Tree seg_tree = (settings.tree_method == GARDTreeMethod::Ohmic && result.num_taxa > 50)
                ? TreeEngine::build_spectral_scaffold(slice)
                : core::NeighborJoining::infer(slice);
            int passes = (settings.run_mode == "Faster") ? 1 : 2;
            auto [opt_tree, lnl] = fit_gtr_segment_srv(
                std::move(seg_tree), slice, result.baseline_rates, result.baseline_frequencies,
                result.gamma_rates, result.gamma_weights, passes
            );

            PartitionFit fit{std::move(opt_tree), lnl};
            {
                std::lock_guard<std::mutex> lock(cache_mutex);
                partition_cache[key] = fit;
            }
            return fit;
        };

        // Master list of evaluated models: canonical key -> c-AIC
        std::unordered_map<std::string, Scalar> master_list;
        std::mutex master_mutex;
        master_list[best_overall_model.to_key()] = result.baseline_aicc;

        auto evaluate_candidate_model = [&](const std::vector<size_t>& bps, bool build_partitions = false) -> GARDModel {
            GARDModel model;
            model.breakpoints = bps;
            std::sort(model.breakpoints.begin(), model.breakpoints.end());

            std::string key = model.to_key();
            if (!build_partitions) {
                std::lock_guard<std::mutex> lock(master_mutex);
                auto it = master_list.find(key);
                if (it != master_list.end()) {
                    model.aicc = it->second;
                    // If cached aicc already known, return stub model with aicc
                    return model;
                }
            }

            size_t num_bp = model.breakpoints.size();
            std::vector<std::pair<size_t, size_t>> segments;
            segments.reserve(num_bp + 1);

            size_t curr_start = 0;
            for (size_t b = 0; b < num_bp; ++b) {
                size_t curr_end = model.breakpoints[b] + 1;
                segments.push_back({curr_start, curr_end});
                curr_start = curr_end;
            }
            segments.push_back({curr_start, result.num_sites});

            Scalar total_lnl = 0.0;
            model.partitions.reserve(segments.size());
            for (const auto& [s, e] : segments) {
                PartitionFit pf = fit_segment(s, e);
                total_lnl += pf.log_likelihood;
                GARDPartition part;
                part.start_site = s + 1;
                part.end_site = e;
                part.tree = pf.tree;
                part.log_likelihood = pf.log_likelihood;
                model.partitions.push_back(std::move(part));
            }

            model.log_likelihood = total_lnl;
            size_t K = 5 + K_srv + (num_bp + 1) * result.min_partition_size;
            if (result.num_sites <= K + 1) {
                model.aicc = std::numeric_limits<Scalar>::infinity();
            } else {
                model.aicc = -2.0 * total_lnl + (2.0 * result.num_sites * K) / (result.num_sites - K - 1.0);
            }

            {
                std::lock_guard<std::mutex> lock(master_mutex);
                master_list[key] = model.aicc;
            }

            return model;
        };

        // 3. Exhaustive Single Breakpoint Analysis (B = 1)
        if (!settings.quiet) {
            std::cout << "[GARD] Performing exhaustive single breakpoint analysis ("
                      << result.potential_breakpoints << " candidate positions)...\n";
        }

        std::vector<size_t> candidate_single_bps;
        for (size_t idx = 0; idx < result.potential_breakpoints; ++idx) {
            size_t bp = result.variable_sites[idx];
            if (validate_breakpoints({bp})) {
                candidate_single_bps.push_back(bp);
            }
        }

        GARDVisualizer visualizer(result.num_sites, result.num_taxa, settings.quiet, settings.force_tty, settings.dynamic_chart);

        GARDModel best_single_bp_model;
        best_single_bp_model.aicc = std::numeric_limits<Scalar>::infinity();

        visualizer.start_single_bp(candidate_single_bps.size(), result.baseline_aicc);

        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < candidate_single_bps.size(); ++i) {
            size_t bp = candidate_single_bps[i];
            GARDModel m = evaluate_candidate_model({bp});

            visualizer.record_candidate({bp}, m.aicc);
            #pragma omp critical
            {
                if (m.aicc < best_single_bp_model.aicc) {
                    best_single_bp_model = m;
                    visualizer.update_best({bp}, m.aicc);
                }
            }
            visualizer.tick(1);
        }

        visualizer.finish_single_bp(best_single_bp_model.breakpoints, best_single_bp_model.aicc, result.baseline_aicc);

        bool single_bp_improves = (best_single_bp_model.aicc < result.baseline_aicc);
        if (single_bp_improves) {
            if (best_single_bp_model.partitions.empty()) {
                best_single_bp_model = evaluate_candidate_model(best_single_bp_model.breakpoints, true);
            }
            Scalar delta = result.baseline_aicc - best_single_bp_model.aicc;
            result.improvements.push_back({1, delta});
            result.improvement_breakpoints.push_back(best_single_bp_model.breakpoints);
            best_overall_model = best_single_bp_model;
            result.best_model = best_overall_model;
        }

        // 4. Multi-Breakpoint CHC Genetic Algorithm (B >= 2)
        bool adding_improves = single_bp_improves;
        size_t current_B = 1;

        std::mt19937 rng(42);

        while (adding_improves && current_B < settings.max_breakpoints &&
               (current_B + 2) * result.min_partition_size <= result.num_sites) {
            current_B += 1;

            visualizer.start_ga(current_B, settings.population_size, best_overall_model.aicc, best_overall_model.breakpoints);

            // Initialize population of models
            std::vector<std::vector<size_t>> population;
            population.reserve(settings.population_size);

            // 1. Seeds with bridge candidate breakpoints (if detected)
            if (!result.bridge_candidate_breakpoints.empty()) {
                for (size_t bp_seed : result.bridge_candidate_breakpoints) {
                    if (population.size() >= settings.population_size / 2) break;
                    std::vector<size_t> cand = best_overall_model.breakpoints;
                    if (std::find(cand.begin(), cand.end(), bp_seed) == cand.end()) {
                        cand.push_back(bp_seed);
                    }
                    std::sort(cand.begin(), cand.end());
                    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
                    int attempts = 0;
                    while (cand.size() < current_B && attempts < 50) {
                        size_t rand_idx = rng() % result.potential_breakpoints;
                        cand.push_back(result.variable_sites[rand_idx]);
                        std::sort(cand.begin(), cand.end());
                        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
                        attempts++;
                    }
                    if (cand.size() == current_B && validate_breakpoints(cand)) {
                        population.push_back(cand);
                    }
                }
            }

            // 2. Seeds with best overall model + 1 new random variable site
            size_t num_seeded = settings.population_size / 2;
            while (population.size() < num_seeded) {
                int attempts = 0;
                while (attempts < 100) {
                    std::vector<size_t> cand = best_overall_model.breakpoints;
                    size_t rand_idx = rng() % result.potential_breakpoints;
                    cand.push_back(result.variable_sites[rand_idx]);
                    std::sort(cand.begin(), cand.end());
                    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
                    if (cand.size() == current_B && validate_breakpoints(cand)) {
                        population.push_back(cand);
                        break;
                    }
                    attempts++;
                }
                if (attempts >= 100) break;
            }

            // 3. Other half purely random
            while (population.size() < settings.population_size) {
                int attempts = 0;
                while (attempts < 100) {
                    std::vector<size_t> cand;
                    cand.reserve(current_B);
                    while (cand.size() < current_B) {
                        size_t rand_idx = rng() % result.potential_breakpoints;
                        size_t bp = result.variable_sites[rand_idx];
                        if (std::find(cand.begin(), cand.end(), bp) == cand.end()) {
                            cand.push_back(bp);
                        }
                    }
                    std::sort(cand.begin(), cand.end());
                    if (validate_breakpoints(cand)) {
                        population.push_back(cand);
                        break;
                    }
                    attempts++;
                }
            }

            // Generational GA loop
            bool ga_terminated = false;
            size_t generation = 0;
            size_t generations_stagnant = 0;
            size_t generations_no_new_models = 0;
            size_t incest_distance = std::max<size_t>(1, current_B / 2);
            Scalar previous_best_caic = std::numeric_limits<Scalar>::infinity();
            GARDModel current_best_model_for_B;
            current_best_model_for_B.aicc = std::numeric_limits<Scalar>::infinity();

            auto count_diff = [&](const std::vector<size_t>& v1, const std::vector<size_t>& v2) -> size_t {
                size_t diff = 0, i1 = 0, i2 = 0;
                while (i1 < v1.size() && i2 < v2.size()) {
                    if (v1[i1] == v2[i2]) { i1++; i2++; }
                    else if (v1[i1] < v2[i2]) { diff++; i1++; }
                    else { diff++; i2++; }
                }
                return diff + (v1.size() - i1) + (v2.size() - i2);
            };

            while (!ga_terminated) {
                generation++;

                // Produce child models via recombination with incest prevention
                std::vector<std::vector<size_t>> child_models;
                child_models.reserve(settings.population_size);

                size_t mating_pairs_accepted = 0;
                for (size_t c = 0; c < settings.population_size; ++c) {
                    int attempts = 0;
                    while (attempts < 10) {
                        size_t p1 = rng() % population.size();
                        size_t p2 = rng() % population.size();
                        while (p2 == p1 && population.size() > 1) {
                            p2 = rng() % population.size();
                        }

                        if (count_diff(population[p1], population[p2]) < incest_distance) {
                            attempts++;
                            continue;
                        }

                        mating_pairs_accepted++;
                        std::vector<size_t> child_bp(current_B);
                        for (size_t b = 0; b < current_B; ++b) {
                            child_bp[b] = (rng() % 2 == 0) ? population[p1][b] : population[p2][b];
                        }
                        std::sort(child_bp.begin(), child_bp.end());
                        child_bp.erase(std::unique(child_bp.begin(), child_bp.end()), child_bp.end());

                        if (child_bp.size() == current_B && validate_breakpoints(child_bp)) {
                            child_models.push_back(std::move(child_bp));
                            break;
                        }
                        attempts++;
                    }
                }

                if (child_models.empty() || mating_pairs_accepted == 0) {
                    if (incest_distance > 0) incest_distance--;
                }

                // Intergenerational pool = parentModels + childModels
                std::vector<std::vector<size_t>> pool = population;
                for (auto& cm : child_models) {
                    pool.push_back(std::move(cm));
                }

                // Unique models
                std::sort(pool.begin(), pool.end());
                pool.erase(std::unique(pool.begin(), pool.end()), pool.end());

                // Evaluate unevaluated models
                std::vector<std::pair<std::vector<size_t>, Scalar>> evaluated_pool(pool.size());
                #pragma omp parallel for schedule(dynamic)
                for (size_t idx = 0; idx < pool.size(); ++idx) {
                    GARDModel m = evaluate_candidate_model(pool[idx]);
                    evaluated_pool[idx] = {pool[idx], m.aicc};
                    visualizer.record_candidate(pool[idx], m.aicc);
                    visualizer.tick(1);
                }

                // Elitist selection: top population_size models
                std::sort(evaluated_pool.begin(), evaluated_pool.end(),
                          [](const auto& a, const auto& b) { return a.second < b.second; });

                std::vector<std::vector<size_t>> selected;
                size_t num_keep = std::min(settings.population_size, evaluated_pool.size());
                selected.reserve(num_keep);
                for (size_t k = 0; k < num_keep; ++k) {
                    selected.push_back(evaluated_pool[k].first);
                }

                Scalar current_best_caic = evaluated_pool[0].second;
                if (current_best_caic < current_best_model_for_B.aicc) {
                    current_best_model_for_B = evaluate_candidate_model(evaluated_pool[0].first, true);
                }

                // Check convergence & diversity
                if (selected == population) {
                    generations_no_new_models++;
                    if (incest_distance > 0) incest_distance--;
                } else {
                    generations_no_new_models = 0;
                }

                std::set<std::vector<size_t>> unique_set(selected.begin(), selected.end());
                Scalar genotypic_diversity = static_cast<Scalar>(unique_set.size()) / static_cast<Scalar>(selected.size());

                size_t max_allowed_no_new = std::max(static_cast<size_t>(2), settings.max_stagnant_generations / 4);
                bool trigger_cataclysm = (genotypic_diversity < settings.caic_diversity_threshold) ||
                                         (incest_distance == 0) ||
                                         (generations_no_new_models > max_allowed_no_new);

                if (trigger_cataclysm) {
                    // Cataclysmic mutation: keep elite, mutate others
                    std::vector<std::vector<size_t>> mutated_pop;
                    mutated_pop.push_back(selected[0]);

                    int max_step = std::max(1, static_cast<int>(result.potential_breakpoints / 100));
                    std::poisson_distribution<int> poiss(max_step);

                    for (size_t i = 1; i < selected.size(); ++i) {
                        int attempts = 0;
                        while (attempts < 10) {
                            std::vector<size_t> mut_bp = selected[i];
                            for (size_t b = 0; b < current_B; ++b) {
                                Scalar r_val = static_cast<Scalar>(rng() % 10000) / 10000.0;
                                if (r_val < settings.mutation_rate) {
                                    Scalar shift_r = static_cast<Scalar>(rng() % 10000) / 10000.0;
                                    if (shift_r < settings.rate_small_shifts) {
                                        int step = std::max(1, poiss(rng));
                                        if (rng() % 2 == 0) step = -step;

                                        // Find index in variable_sites
                                        auto it = std::lower_bound(result.variable_sites.begin(), result.variable_sites.end(), mut_bp[b]);
                                        if (it != result.variable_sites.end()) {
                                            int orig_v_idx = static_cast<int>(std::distance(result.variable_sites.begin(), it));
                                            int new_v_idx = orig_v_idx + step;
                                            if (new_v_idx >= 0 && new_v_idx < static_cast<int>(result.variable_sites.size())) {
                                                mut_bp[b] = result.variable_sites[new_v_idx];
                                            }
                                        }
                                    } else {
                                        size_t rand_v = rng() % result.potential_breakpoints;
                                        mut_bp[b] = result.variable_sites[rand_v];
                                    }
                                }
                            }
                            std::sort(mut_bp.begin(), mut_bp.end());
                            mut_bp.erase(std::unique(mut_bp.begin(), mut_bp.end()), mut_bp.end());
                            if (mut_bp.size() == current_B && validate_breakpoints(mut_bp)) {
                                mutated_pop.push_back(std::move(mut_bp));
                                break;
                            }
                            attempts++;
                        }
                    }

                    // Reset incest distance threshold after cataclysmic restart
                    incest_distance = std::max<size_t>(1, current_B / 2);
                    generations_no_new_models = 0;

                    if (mutated_pop.size() <= 1) {
                        // Reinitialize with seeds
                        mutated_pop.clear();
                        mutated_pop.push_back(selected[0]);
                        while (mutated_pop.size() < settings.population_size) {
                            std::vector<size_t> cand = best_overall_model.breakpoints;
                            size_t rand_idx = rng() % result.potential_breakpoints;
                            cand.push_back(result.variable_sites[rand_idx]);
                            std::sort(cand.begin(), cand.end());
                            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
                            if (cand.size() == current_B && validate_breakpoints(cand)) {
                                mutated_pop.push_back(cand);
                            }
                        }
                    }
                    population = std::move(mutated_pop);
                } else {
                    population = std::move(selected);
                }

                // Check stagnation
                if (previous_best_caic - current_best_caic < settings.caic_improvement_threshold) {
                    generations_stagnant++;
                    if (generations_stagnant >= settings.max_stagnant_generations) {
                        ga_terminated = true;
                    }
                } else {
                    generations_stagnant = 0;
                }
                previous_best_caic = current_best_caic;

                visualizer.update_ga_generation(
                    generation,
                    settings.max_stagnant_generations,
                    generations_stagnant,
                    genotypic_diversity,
                    current_best_model_for_B.breakpoints,
                    current_best_model_for_B.aicc
                );
            } // end GA loop for current_B

            // Evaluate if this B improved over best overall
            bool improved = (current_best_model_for_B.aicc < best_overall_model.aicc);
            visualizer.finish_ga(current_B, current_best_model_for_B.breakpoints, current_best_model_for_B.aicc, best_overall_model.aicc, improved);

            if (improved) {
                if (current_best_model_for_B.partitions.empty()) {
                    current_best_model_for_B = evaluate_candidate_model(current_best_model_for_B.breakpoints, true);
                }
                Scalar delta = best_overall_model.aicc - current_best_model_for_B.aicc;
                result.improvements.push_back({current_B, delta});
                result.improvement_breakpoints.push_back(current_best_model_for_B.breakpoints);
                best_overall_model = current_best_model_for_B;
                result.best_model = best_overall_model;
            } else {
                adding_improves = false;
            }
        } // end multi-breakpoint loop

        // 5. Post-hoc Topological Incongruence Test
        // Fit best breakpoint placement constrained to baseline single tree topology
        if (!settings.quiet) {
            std::cout << "[GARD] Running post-hoc topological incongruence test...\n";
        }

        if (best_overall_model.breakpoints.empty()) {
            result.single_tree_aicc = result.baseline_aicc;
            result.topological_incongruence = false;
            result.topology_test_result.run = false;
        } else {
            size_t num_bp = best_overall_model.breakpoints.size();
            size_t num_parts = num_bp + 1;
            std::vector<std::pair<size_t, size_t>> segments;
            size_t curr_start = 0;
            for (size_t b = 0; b < num_bp; ++b) {
                size_t curr_end = best_overall_model.breakpoints[b] + 1;
                segments.push_back({curr_start, curr_end});
                curr_start = curr_end;
            }
            segments.push_back({curr_start, result.num_sites});

            Scalar single_tree_log_l = 0.0;
            int passes = (settings.run_mode == "Faster") ? 1 : 2;
            for (const auto& [s, e] : segments) {
                core::Alignment slice = alignment.slice_nucleotide(s, e);
                auto [opt_t, lnl] = fit_gtr_segment_srv(
                    result.baseline_tree, slice, result.baseline_rates, result.baseline_frequencies,
                    result.gamma_rates, result.gamma_weights, passes
                );
                single_tree_log_l += lnl;
            }

            size_t K_srv = (settings.srv_mode == GARDSRVMode::DiscreteGamma) ? 1 : 0;
            size_t K = 5 + K_srv + (num_bp + 1) * result.min_partition_size;
            result.single_tree_aicc = -2.0 * single_tree_log_l +
                (2.0 * result.num_sites * K) / (result.num_sites - K - 1.0);

            result.topological_incongruence = (best_overall_model.aicc < result.single_tree_aicc);

            if (settings.topology_test_mode != GARDTopologyTestMode::None) {
                result.topology_test_result.run = true;
                result.topology_test_result.test_type =
                    (settings.topology_test_mode == GARDTopologyTestMode::Both) ? "both" :
                    (settings.topology_test_mode == GARDTopologyTestMode::SH) ? "sh" : "kh";
                result.topology_test_result.bootstrap_replicates = settings.topology_test_replicates;

                GTRMatrix base_gtr;
                base_gtr.update(result.baseline_rates, result.baseline_frequencies);

                // Ensure partitions are fully built on best_overall_model
                if (best_overall_model.partitions.size() != num_parts) {
                    best_overall_model = evaluate_candidate_model(best_overall_model.breakpoints, true);
                    result.best_model = best_overall_model;
                }

                // site_log_l[p][t] = per-site log-likelihood vector of Tree t evaluated on Partition p
                std::vector<std::vector<std::vector<Scalar>>> site_log_l(
                    num_parts, std::vector<std::vector<Scalar>>(num_parts));

                for (size_t p = 0; p < num_parts; ++p) {
                    auto [s, e] = segments[p];
                    core::Alignment slice = alignment.slice_nucleotide(s, e);

                    for (size_t t = 0; t < num_parts; ++t) {
                        Tree tree_to_eval;
                        if (t == p) {
                            tree_to_eval = best_overall_model.partitions[p].tree;
                        } else {
                            auto [opt_t, _] = fit_gtr_segment_srv(
                                best_overall_model.partitions[t].tree, slice, result.baseline_rates, result.baseline_frequencies,
                                result.gamma_rates, result.gamma_weights, std::max(passes, 3)
                            );
                            tree_to_eval = std::move(opt_t);
                        }

                        std::vector<Scalar> s_lnl;
                        compute_gtr_srv_log_likelihood(
                            tree_to_eval, slice, base_gtr, result.gamma_rates, result.gamma_weights, nullptr, &s_lnl
                        );
                        site_log_l[p][t] = std::move(s_lnl);
                    }
                }

                // 1. Breakpoint-specific tests (Adjacent partitions across each breakpoint)
                bool any_sig = false;
                for (size_t b = 0; b < num_bp; ++b) {
                    size_t L = b;
                    size_t R = b + 1;
                    size_t bp_site = best_overall_model.breakpoints[b] + 1;

                    BreakpointTopologyTest bp_test;
                    bp_test.breakpoint_site = bp_site;
                    bp_test.left_partition_index = L;
                    bp_test.right_partition_index = R;
                    bp_test.left_interval = {segments[L].first + 1, segments[L].second};
                    bp_test.right_interval = {segments[R].first + 1, segments[R].second};

                    // On Left segment L: compare Tree L (MLE) vs Tree R
                    const auto& v_LL = site_log_l[L][L];
                    const auto& v_LR = site_log_l[L][R];
                    Scalar sum_LL = 0.0, sum_LR = 0.0;
                    for (Scalar x : v_LL) sum_LL += x;
                    for (Scalar x : v_LR) sum_LR += x;
                    bp_test.left_tree_on_left_lnl = sum_LL;
                    bp_test.right_tree_on_left_lnl = sum_LR;
                    bp_test.delta_lnl_left = sum_LL - sum_LR;

                    auto kh_L = GARDTopologyTester::run_kh_test(v_LL, v_LR, settings.topology_test_replicates, 42 + b * 2);
                    auto sh_L = GARDTopologyTester::run_sh_test({v_LL, v_LR}, settings.topology_test_replicates, 42 + b * 2);
                    bp_test.kh_p_value_left = kh_L.p_rell;
                    bp_test.sh_p_value_left = (sh_L.best_index == 0) ? sh_L.p_values[1] : sh_L.p_values[0];

                    // On Right segment R: compare Tree R (MLE) vs Tree L
                    const auto& v_RR = site_log_l[R][R];
                    const auto& v_RL = site_log_l[R][L];
                    Scalar sum_RR = 0.0, sum_RL = 0.0;
                    for (Scalar x : v_RR) sum_RR += x;
                    for (Scalar x : v_RL) sum_RL += x;
                    bp_test.right_tree_on_right_lnl = sum_RR;
                    bp_test.left_tree_on_right_lnl = sum_RL;
                    bp_test.delta_lnl_right = sum_RR - sum_RL;

                    auto kh_R = GARDTopologyTester::run_kh_test(v_RR, v_RL, settings.topology_test_replicates, 43 + b * 2);
                    auto sh_R = GARDTopologyTester::run_sh_test({v_RR, v_RL}, settings.topology_test_replicates, 43 + b * 2);
                    bp_test.kh_p_value_right = kh_R.p_rell;
                    bp_test.sh_p_value_right = (sh_R.best_index == 0) ? sh_R.p_values[1] : sh_R.p_values[0];

                    // Significance decision (Kosakovsky Pond et al. 2006: both flanking segments must
                    // significantly reject congruence after Bonferroni correction for 2 * num_bp tests)
                    Scalar bonf_alpha = 0.05 / static_cast<Scalar>(2 * std::max<size_t>(1, num_bp));
                    bool sig_kh = (bp_test.kh_p_value_left <= bonf_alpha && bp_test.kh_p_value_right <= bonf_alpha);
                    bool sig_sh = (bp_test.sh_p_value_left <= bonf_alpha && bp_test.sh_p_value_right <= bonf_alpha);
                    if (settings.topology_test_mode == GARDTopologyTestMode::KH) {
                        bp_test.statistically_significant = sig_kh;
                    } else if (settings.topology_test_mode == GARDTopologyTestMode::SH) {
                        bp_test.statistically_significant = sig_sh;
                    } else { // Both
                        bp_test.statistically_significant = sig_kh || sig_sh;
                    }

                    if (bp_test.statistically_significant) any_sig = true;
                    result.topology_test_result.breakpoint_tests.push_back(bp_test);
                }
                result.topology_test_result.any_incongruence_significant = any_sig;

                // 2. Full Partition Tree Matrix
                for (size_t p = 0; p < num_parts; ++p) {
                    PartitionTreeComparison comp;
                    comp.partition_index = p;
                    comp.interval = {segments[p].first + 1, segments[p].second};

                    auto sh_matrix = GARDTopologyTester::run_sh_test(
                        site_log_l[p], settings.topology_test_replicates, 100 + p
                    );
                    comp.best_tree_index = sh_matrix.best_index;
                    comp.tree_log_likelihoods = sh_matrix.total_lnl;
                    comp.sh_p_values = sh_matrix.p_values;

                    comp.kh_p_values.resize(num_parts, 1.0);
                    for (size_t t = 0; t < num_parts; ++t) {
                        if (t != comp.best_tree_index) {
                            auto kh_res = GARDTopologyTester::run_kh_test(
                                site_log_l[p][comp.best_tree_index], site_log_l[p][t],
                                settings.topology_test_replicates, 200 + p * num_parts + t
                            );
                            comp.kh_p_values[t] = kh_res.p_rell;
                        }
                    }
                    result.topology_test_result.partition_matrix.push_back(comp);
                }
            }
        }

        if (!settings.quiet) {
            std::cout << "[GARD] Single-tree c-AIC = " << std::fixed << std::setprecision(4)
                      << result.single_tree_aicc << " vs Partitioned c-AIC = " << best_overall_model.aicc << "\n";
            std::cout << "[GARD] Topological incongruence (recombination): "
                      << (result.topological_incongruence ? "SUPPORTED" : "NOT supported") << "\n";

            if (result.topology_test_result.run && !result.topology_test_result.breakpoint_tests.empty()) {
                std::cout << "[GARD] Breakpoint Topological Incongruence Tests (KH / SH with RELL bootstrap, "
                          << settings.topology_test_replicates << " reps):\n";
                for (const auto& bp_t : result.topology_test_result.breakpoint_tests) {
                    std::cout << "  • Breakpoint at site " << bp_t.breakpoint_site
                              << " (Partitions [" << bp_t.left_interval.first << ".." << bp_t.left_interval.second << "]"
                              << " vs [" << bp_t.right_interval.first << ".." << bp_t.right_interval.second << "]):\n"
                              << "      KH p-value: Left = " << std::fixed << std::setprecision(4) << bp_t.kh_p_value_left
                              << ", Right = " << bp_t.kh_p_value_right << "\n"
                              << "      SH p-value: Left = " << std::fixed << std::setprecision(4) << bp_t.sh_p_value_left
                              << ", Right = " << bp_t.sh_p_value_right << "\n"
                              << "      Verdict   : "
                              << (bp_t.statistically_significant ? "Significant Incongruence (Recombination)" : "Not Significant (Rate Shift)")
                              << "\n";
                }
            }
        }

        // 6. Compute Akaike weights for siteBreakPointSupport
        result.total_models_evaluated = master_list.size();
        Scalar min_master_caic = best_overall_model.aicc;

        std::unordered_map<size_t, Scalar> site_weights;
        Scalar total_aic_weight = 0.0;

        for (const auto& [key, score] : master_list) {
            if (std::isinf(score) || std::isnan(score)) continue;
            Scalar aicw = std::exp(-0.5 * (score - min_master_caic));
            if (aicw > 0.0) {
                total_aic_weight += aicw;
                std::vector<size_t> bps = GARDModel::parse_key(key);
                for (size_t bp : bps) {
                    site_weights[bp] += aicw;
                }
            }
        }

        result.site_breakpoint_support.clear();
        if (total_aic_weight > 0.0) {
            for (const auto& [bp, w] : site_weights) {
                result.site_breakpoint_support[std::to_string(bp)] = w / total_aic_weight;
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        result.time_elapsed = std::chrono::duration<double>(end_time - start_time).count();

        if (!settings.quiet) {
            std::cout << "[GARD] Total models evaluated: " << result.total_models_evaluated
                      << " in " << std::fixed << std::setprecision(2) << result.time_elapsed << " seconds.\n";
        }

        return result;
    }
};

} // namespace hyphy::analyses
