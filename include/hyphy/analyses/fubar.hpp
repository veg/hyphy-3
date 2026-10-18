#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/progress_bar.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/core/console.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <memory>
#include <random>
#include <chrono>
#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

enum class FUBARMethod {
    VariationalBayes,  // 0-th order Variational Bayes (VB0, recommended default)
    CollapsedGibbs,    // Collapsed Gibbs sampler (CG)
    MetropolisHastings // Full Metropolis-Hastings MCMC (MH)
};

inline std::string fubar_method_to_string(FUBARMethod m) {
    switch (m) {
        case FUBARMethod::VariationalBayes: return "Variational-Bayes";
        case FUBARMethod::CollapsedGibbs: return "Collapsed-Gibbs";
        case FUBARMethod::MetropolisHastings: return "Metropolis-Hastings";
    }
    return "Variational-Bayes";
}

inline FUBARMethod string_to_fubar_method(const std::string& str) {
    std::string s = str;
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    if (s == "variational-bayes" || s == "variational_bayes" || s == "vb" || s == "vb0") {
        return FUBARMethod::VariationalBayes;
    }
    if (s == "collapsed-gibbs" || s == "collapsed_gibbs" || s == "cg" || s == "gibbs") {
        return FUBARMethod::CollapsedGibbs;
    }
    if (s == "metropolis-hastings" || s == "metropolis_hastings" || s == "mh" || s == "mcmc") {
        return FUBARMethod::MetropolisHastings;
    }
    return FUBARMethod::VariationalBayes;
}

inline std::string format_float_val(Scalar val, int prec = 3) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(prec) << val;
    return ss.str();
}

struct FUBARSettings {
    int grid_points = 20;                     // Points per dimension (total G = D^2 = 400)
    FUBARMethod method = FUBARMethod::VariationalBayes;
    Scalar concentration = 0.5;              // Dirichlet prior hyperparameter alpha
    Scalar posterior_threshold = 0.9;         // Posterior cutoff for selection calls
    bool non_zero = false;                   // Enforce non-zero synonymous rate (min 1e-3)
    std::string tree_scaling = "grid";       // "grid" (default, HyPhy 2.5 standard) or "mg94"
    
    // MCMC settings
    int mcmc_chains = 1;
    int64_t mcmc_chain_length = 2000000;
    int64_t mcmc_burn_in = 1000000;
    int64_t mcmc_samples = 100;
    
    // Variational Bayes settings
    Scalar vb_tolerance = 1e-8;
    int vb_max_steps = 100000;
};

struct FUBARSiteResult {
    size_t site = 0;                          // 1-indexed codon site
    Scalar alpha = 0.0;                       // Mean posterior synonymous rate
    Scalar beta = 0.0;                        // Mean posterior non-synonymous rate
    Scalar beta_minus_alpha = 0.0;            // Mean posterior beta - alpha
    Scalar prob_neg = 0.0;                    // Posterior prob P(alpha > beta | D_s)
    Scalar prob_pos = 0.0;                    // Posterior prob P(beta > alpha | D_s)
    Scalar bayes_factor = 1.0;                // Empirical Bayes Factor for positive selection
    std::string classification = "neutral";   // "diversifying", "purifying", "neutral"
    std::vector<Scalar> posterior_weights;    // Size G: P(g | D_s)
};

struct FUBARResult {
    std::vector<FUBARSiteResult> site_results;
    std::vector<std::pair<Scalar, Scalar>> grid; // Size G: (alpha_g, beta_g)
    std::vector<Scalar> posterior_grid_weights;  // Size G: theta_g
    std::vector<size_t> positive_sites;         // 1-indexed sites with prob_pos >= threshold
    std::vector<size_t> purifying_sites;        // 1-indexed sites with prob_neg >= threshold
    Scalar expected_false_positives = 0.0;
    std::pair<int, int> false_positives_ci_95 = {0, 0};
    Scalar weight_non_positive = 0.0;           // Prior weight on alpha >= beta
    
    // Model fit info
    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;
    
    Scalar baseline_tree_length = 0.0;
    Scalar best_scaler_alpha = 1.0;
    Scalar best_scaler_beta = 1.0;
    Scalar scaled_tree_length = 0.0;
    
    // Timings (seconds)
    Scalar time_gtr_sec = 0.0;
    Scalar time_grid_sec = 0.0;
    Scalar time_posterior_sec = 0.0;
    Scalar time_total_sec = 0.0;
};

class FUBARAnalyzer {
public:
    Tree tree;
    Alignment aln;
    MG94Parameters base_params;
    FUBARSettings settings;
    
    bool has_gtr_fit = false;
    Scalar gtr_log_l = 0.0;
    Scalar gtr_aicc = 0.0;
    GTRParameters gtr_rates;
    std::unordered_map<std::string, Scalar> gtr_branch_lengths;
    Scalar conversion_factor = 1.0;
    Scalar baseline_tree_length = 0.0;
    Scalar time_gtr_sec = 0.0;

    FUBARAnalyzer(Tree t, Alignment a, MG94Parameters params = MG94Parameters{}, FUBARSettings s = FUBARSettings{})
        : tree(std::move(t)), aln(std::move(a)), base_params(params), settings(std::move(s)) {
        for (const auto& node : tree.nodes) {
            if (node.id != tree.root_id) {
                baseline_tree_length += node.branch_length;
            }
        }
    }

    // Static factory that performs Phase 1 GTR and branch-length conversion
    static FUBARAnalyzer create_and_fit(
        Tree input_tree,
        Alignment aln,
        FUBARSettings settings = FUBARSettings{},
        std::function<void(const std::string&, double)> progress_cb = nullptr
    ) {
        auto start_gtr = std::chrono::high_resolution_clock::now();
        if (progress_cb) progress_cb("Phase 1: Fitting Nucleotide GTR Model", 0.05);

        GTRFitter gtr_fitter(input_tree, aln);
        auto gtr_res = gtr_fitter.fit();
        auto end_gtr = std::chrono::high_resolution_clock::now();
        Scalar gtr_sec = std::chrono::duration<double>(end_gtr - start_gtr).count();

        MG94Parameters base_p;
        base_p.theta_AC = gtr_res.params.theta_AC;
        base_p.theta_AT = gtr_res.params.theta_AT;
        base_p.theta_CG = gtr_res.params.theta_CG;
        base_p.theta_CT = gtr_res.params.theta_CT;
        base_p.theta_GT = gtr_res.params.theta_GT;
        base_p.alpha = 1.0;
        base_p.beta = 1.0;

        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        MG94Matrix base_mg;
        base_mg.update(base_p, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

        // Convert nucleotide expected substitutions to codon expected substitutions
        // matching HyPhy 2.5 fubar.scalers.SetBranchLength
        Scalar codon_scale = base_mg.scale_factor;
        Scalar conv = (codon_scale > 1e-12) ? (3.0 / codon_scale) : 1.0;

        Tree codon_tree = gtr_res.tree;
        for (auto& node : codon_tree.nodes) {
            if (node.id != codon_tree.root_id) {
                node.branch_length *= conv;
            }
        }
        codon_tree.scale_bridge_lengths(conv);

        if (settings.tree_scaling == "mg94") {
            if (progress_cb) progress_cb("Phase 1b: Refining Tree Scale under MG94", 0.15);
            MG94Fitter mg_fitter(codon_tree, aln);
            auto mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
            codon_tree = mg_res.tree;
        }

        FUBARAnalyzer analyzer(codon_tree, std::move(aln), base_p, settings);
        analyzer.has_gtr_fit = true;
        analyzer.gtr_log_l = gtr_res.log_likelihood;
        analyzer.gtr_aicc = gtr_res.aicc;
        analyzer.gtr_rates = gtr_res.params;
        analyzer.conversion_factor = conv;
        analyzer.time_gtr_sec = gtr_sec;

        for (const auto& node : gtr_res.tree.nodes) {
            if (node.id != gtr_res.tree.root_id) {
                analyzer.gtr_branch_lengths[node.name] = node.branch_length;
            }
        }

        return analyzer;
    }

    // Generate bivariate (alpha, beta) grid matching HyPhy 2.5 fubar.DefineAlphaBetaGrid
    static std::vector<std::pair<Scalar, Scalar>> create_grid(int d, bool non_zero) {
        d = std::max(d, 5);
        int neg_points = static_cast<int>(std::floor(d * 0.7 + 0.5));
        int pos_points = static_cast<int>(std::floor((d - 1) * 0.3));
        if (neg_points + pos_points != d) {
            pos_points = d - neg_points;
        }

        std::vector<Scalar> u(d);
        Scalar neg_step = 1.0 / neg_points;
        for (int k = 0; k < neg_points; ++k) {
            u[k] = neg_step * k;
        }
        u[neg_points - 1] = 1.0;

        Scalar pos_step = std::pow(49.0, 1.0 / 3.0) / pos_points;
        for (int k = 1; k <= pos_points; ++k) {
            Scalar term = pos_step * k;
            u[neg_points + k - 1] = 1.0 + term * term * term;
        }

        Scalar min_val = non_zero ? std::max(1e-3, u[0]) : 0.0;
        std::vector<std::pair<Scalar, Scalar>> grid;
        grid.reserve(d * d);
        for (int r = 0; r < d; ++r) {
            for (int c = 0; c < d; ++c) {
                Scalar a = std::max(u[r], min_val);
                Scalar b = u[c];
                grid.push_back({a, b});
            }
        }
        return grid;
    }

    // Execute full FUBAR analysis
    FUBARResult run(
        std::function<void(const std::string&, double)> progress_cb = nullptr,
        bool verbose = true
    ) {
        auto start_total = std::chrono::high_resolution_clock::now();
        FUBARResult res;
        res.has_gtr_fit = has_gtr_fit;
        res.gtr_log_l = gtr_log_l;
        res.gtr_aicc = gtr_aicc;
        res.gtr_rates = gtr_rates;
        res.gtr_branch_lengths = gtr_branch_lengths;
        res.baseline_tree_length = baseline_tree_length;
        res.time_gtr_sec = time_gtr_sec;

        // 1. Grid definition
        int D = settings.grid_points;
        auto grid = create_grid(D, settings.non_zero);
        size_t G = grid.size();
        res.grid = grid;

        const auto& gcode = *(aln.code ? aln.code : GeneticCode::universal());
        int S_codons = gcode.num_sense_codons;
        size_t num_nodes = tree.num_nodes();
        size_t num_patterns = aln.patterns.size();
        size_t num_sites = aln.num_codons;

        auto leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        // 2. Phase 2: Compute conditional site likelihoods on the grid
        auto start_grid = std::chrono::high_resolution_clock::now();
        if (progress_cb) progress_cb("Phase 2: Computing Site Likelihoods on Rate Grid", 0.20);

        // Matrix of site log-likelihoods: [G x num_sites]
        std::vector<std::vector<Scalar>> site_log_l(G, std::vector<Scalar>(num_sites, -1e20));
        std::vector<Scalar> total_log_l_grid(G, 0.0);

        std::unique_ptr<ProgressBar> pbar_grid;
        if (verbose && ProgressBar::is_terminal()) {
            pbar_grid = std::make_unique<ProgressBar>(
                G, "[FUBAR] Rate Grid", "pts", false, ProgressBar::Style::SmoothBlocks
            );
        }

        #pragma omp parallel for schedule(dynamic)
        for (size_t g = 0; g < G; ++g) {
            Scalar a = grid[g].first;
            Scalar b = grid[g].second;

            MG94Parameters params = base_params;
            params.alpha = a;
            params.beta = b;

            MG94Matrix mg;
            mg.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

            auto P_branches = LikelihoodEngine::build_P_branches<Matrix>(tree, mg);

            auto [thread_node_L, thread_node_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S_codons);

            std::vector<Scalar> pattern_lnl(num_patterns, -1e20);
            Scalar total_ll = 0.0;

            for (size_t p = 0; p < num_patterns; ++p) {
                Scalar lnl = LikelihoodEngine::compute_pattern_log_likelihood(
                    tree, aln.patterns[p], leaf_to_taxon, mg, P_branches, thread_node_L, thread_node_scale
                );
                pattern_lnl[p] = lnl;
                total_ll += aln.patterns[p].weight * lnl;
            }

            total_log_l_grid[g] = total_ll;

            // Map pattern log-likelihoods to individual sites
            for (size_t s = 0; s < num_sites; ++s) {
                size_t p_idx = aln.site_to_pattern[s];
                site_log_l[g][s] = pattern_lnl[p_idx];
            }

            #pragma omp critical
            {
                if (pbar_grid) pbar_grid->tick();
            }
        }
        if (pbar_grid) pbar_grid->finish();

        // Check best scaling on the grid
        size_t best_grid_idx = 0;
        Scalar max_grid_ll = total_log_l_grid[0];
        for (size_t g = 1; g < G; ++g) {
            if (total_log_l_grid[g] > max_grid_ll) {
                max_grid_ll = total_log_l_grid[g];
                best_grid_idx = g;
            }
        }

        res.best_scaler_alpha = grid[best_grid_idx].first;
        res.best_scaler_beta = grid[best_grid_idx].second;

        // If best grid point differs substantially from (1.0, 1.0) and settings require scaling:
        bool needs_rescaling = (settings.tree_scaling == "grid") &&
                               (std::abs(res.best_scaler_alpha - 1.0) > 1e-4 || std::abs(res.best_scaler_beta - 1.0) > 1e-4);

        if (needs_rescaling) {
            if (progress_cb) progress_cb("Phase 2b: Re-evaluating Grid with Optimal Tree Scaling", 0.40);
            
            MG94Matrix base_mg;
            base_mg.update(base_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

            MG94Parameters opt_params = base_params;
            opt_params.alpha = res.best_scaler_alpha;
            opt_params.beta = res.best_scaler_beta;

            MG94Matrix opt_mg;
            opt_mg.update(opt_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

            Scalar tree_scale = (base_mg.scale_factor > 1e-12) ? (opt_mg.scale_factor / base_mg.scale_factor) : 1.0;

            for (auto& node : tree.nodes) {
                if (node.id != tree.root_id) {
                    node.branch_length *= tree_scale;
                }
            }
            tree.scale_bridge_lengths(tree_scale);
            baseline_tree_length *= tree_scale;
            res.baseline_tree_length = baseline_tree_length;

            std::unique_ptr<ProgressBar> pbar_rescale;
            if (verbose && ProgressBar::is_terminal()) {
                pbar_rescale = std::make_unique<ProgressBar>(
                    G, "[FUBAR] Rescaling Grid", "pts", false, ProgressBar::Style::SmoothBlocks
                );
            }

            #pragma omp parallel for schedule(dynamic)
            for (size_t g = 0; g < G; ++g) {
                Scalar a = grid[g].first;
                Scalar b = grid[g].second;

                MG94Parameters params = base_params;
                params.alpha = a;
                params.beta = b;

                MG94Matrix mg;
                mg.update(params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, gcode);

                auto P_branches = LikelihoodEngine::build_P_branches<Matrix>(tree, mg);

                auto [thread_node_L, thread_node_scale] = LikelihoodEngine::make_pattern_scratch<Vector>(num_nodes, S_codons);

                std::vector<Scalar> pattern_lnl(num_patterns, -1e20);

                for (size_t p = 0; p < num_patterns; ++p) {
                    Scalar lnl = LikelihoodEngine::compute_pattern_log_likelihood(
                        tree, aln.patterns[p], leaf_to_taxon, mg, P_branches, thread_node_L, thread_node_scale
                    );
                    pattern_lnl[p] = lnl;
                }

                for (size_t s = 0; s < num_sites; ++s) {
                    size_t p_idx = aln.site_to_pattern[s];
                    site_log_l[g][s] = pattern_lnl[p_idx];
                }

                #pragma omp critical
                {
                    if (pbar_rescale) pbar_rescale->tick();
                }
            }
            if (pbar_rescale) pbar_rescale->finish();
        }

        auto end_grid = std::chrono::high_resolution_clock::now();
        res.time_grid_sec = std::chrono::duration<double>(end_grid - start_grid).count();

        // Convert site log-likelihoods into normalized conditionals C (G x num_sites)
        // matching HyPhy 2.5 ConvertToConditionals
        Matrix C(G, num_sites);
        std::vector<Scalar> site_scalers(num_sites, 0.0);

        for (size_t s = 0; s < num_sites; ++s) {
            Scalar max_ll = site_log_l[0][s];
            for (size_t g = 1; g < G; ++g) {
                if (site_log_l[g][s] > max_ll) {
                    max_ll = site_log_l[g][s];
                }
            }

            Scalar norm_sum = 0.0;
            for (size_t g = 0; g < G; ++g) {
                Scalar diff = site_log_l[g][s] - max_ll;
                Scalar w = (diff > -700.0) ? std::exp(diff) : 0.0;
                C(g, s) = w;
                norm_sum += w;
            }

            Scalar inv_norm = (norm_sum > 0.0) ? (1.0 / norm_sum) : 1.0;
            for (size_t g = 0; g < G; ++g) {
                C(g, s) *= inv_norm;
            }

            site_scalers[s] = max_ll + std::log(std::max(norm_sum, 1e-300));
        }

        // 3. Phase 3: Posterior estimation over the rate grid
        auto start_post = std::chrono::high_resolution_clock::now();
        if (progress_cb) progress_cb("Phase 3: Estimating Posterior Grid Weights", 0.70);

        std::vector<Scalar> posterior_grid_weights;
        if (settings.method == FUBARMethod::VariationalBayes) {
            posterior_grid_weights = run_variational_bayes(
                C, settings.concentration, settings.vb_tolerance, settings.vb_max_steps
            );
        } else if (settings.method == FUBARMethod::CollapsedGibbs) {
            posterior_grid_weights = run_collapsed_gibbs(
                C, settings.concentration, settings.mcmc_chain_length, settings.mcmc_burn_in, settings.mcmc_samples
            );
        } else {
            posterior_grid_weights = run_metropolis_hastings(
                C, settings.concentration, settings.mcmc_chain_length, settings.mcmc_burn_in, settings.mcmc_samples
            );
        }

        auto end_post = std::chrono::high_resolution_clock::now();
        res.time_posterior_sec = std::chrono::duration<double>(end_post - start_post).count();
        res.posterior_grid_weights = posterior_grid_weights;

        // Prior weight assigned to negative or neutral regime (alpha >= beta)
        Scalar weight_non_pos = 0.0;
        for (size_t g = 0; g < G; ++g) {
            if (grid[g].first >= grid[g].second) {
                weight_non_pos += posterior_grid_weights[g];
            }
        }
        res.weight_non_positive = weight_non_pos;

        // 4. Phase 4: Site-level inference & tabulation
        if (progress_cb) progress_cb("Phase 4: Tabulating Site-Level Results", 0.90);
        res.site_results.resize(num_sites);

        std::vector<Scalar> p_err_positive;

        for (size_t s = 0; s < num_sites; ++s) {
            FUBARSiteResult& site_res = res.site_results[s];
            site_res.site = s + 1;
            site_res.posterior_weights.resize(G);

            Scalar denom = 0.0;
            for (size_t g = 0; g < G; ++g) {
                Scalar unnorm = posterior_grid_weights[g] * C(g, s);
                site_res.posterior_weights[g] = unnorm;
                denom += unnorm;
            }

            Scalar inv_denom = (denom > 0.0) ? (1.0 / denom) : (1.0 / G);
            Scalar mean_a = 0.0;
            Scalar mean_b = 0.0;
            Scalar p_neg = 0.0;
            Scalar p_pos = 0.0;

            for (size_t g = 0; g < G; ++g) {
                Scalar pp = site_res.posterior_weights[g] * inv_denom;
                site_res.posterior_weights[g] = pp;

                mean_a += grid[g].first * pp;
                mean_b += grid[g].second * pp;

                if (grid[g].first > grid[g].second) {
                    p_neg += pp;
                } else if (grid[g].first < grid[g].second) {
                    p_pos += pp;
                }
            }

            site_res.alpha = mean_a;
            site_res.beta = mean_b;
            site_res.beta_minus_alpha = mean_b - mean_a;
            site_res.prob_neg = p_neg;
            site_res.prob_pos = p_pos;

            // Empirical Bayes Factor for positive selection
            if (weight_non_pos > 0.0 && weight_non_pos < 1.0 && p_pos < 1.0) {
                site_res.bayes_factor = (p_pos / (1.0 - p_pos)) / ((1.0 - weight_non_pos) / weight_non_pos);
            } else if (p_pos >= 1.0) {
                site_res.bayes_factor = 1e6;
            } else {
                site_res.bayes_factor = 1.0;
            }

            // Classification
            if (p_pos >= settings.posterior_threshold) {
                site_res.classification = "diversifying";
                res.positive_sites.push_back(s + 1);
                p_err_positive.push_back(1.0 - p_pos);
            } else if (p_neg >= settings.posterior_threshold) {
                site_res.classification = "purifying";
                res.purifying_sites.push_back(s + 1);
            } else {
                site_res.classification = "neutral";
            }
        }

        // Expected number of false positives & 95% Confidence Interval
        res.expected_false_positives = 0.0;
        for (Scalar p_err : p_err_positive) {
            res.expected_false_positives += p_err;
        }
        res.false_positives_ci_95 = compute_enfp_ci(p_err_positive, 0.05);

        auto end_total = std::chrono::high_resolution_clock::now();
        res.time_total_sec = std::chrono::duration<double>(end_total - start_total).count();

        if (progress_cb) progress_cb("FUBAR Analysis Complete", 1.0);
        return res;
    }

    // 0-th order Variational Bayes procedure
    static std::vector<Scalar> run_variational_bayes(
        const Matrix& C,
        Scalar concentration,
        Scalar tolerance = 1e-8,
        int max_steps = 100000,
        std::function<void(int, Scalar)> progress_cb = nullptr
    ) {
        size_t G = C.rows();
        size_t S = C.cols();

        // Initial loading: row sums of conditionals
        Vector grid_weights = C.rowwise().sum();

        for (int step = 1; step <= max_steps; ++step) {
            Vector last_grid_weights = grid_weights;

            // Site normalizations: norm(s) = sum_g last_weights(g) * C(g, s)
            Vector norm = Vector::Zero(S);
            for (size_t s = 0; s < S; ++s) {
                Scalar sum_s = 0.0;
                for (size_t g = 0; g < G; ++g) {
                    sum_s += last_grid_weights(g) * C(g, s);
                }
                norm(s) = sum_s;
            }

            // Update weights: grid_weights(g) = sum_s (last_weights(g) * C(g, s) / norm(s)) + concentration
            grid_weights.setConstant(concentration);
            for (size_t s = 0; s < S; ++s) {
                Scalar inv_norm = (norm(s) > 1e-300) ? (1.0 / norm(s)) : 0.0;
                for (size_t g = 0; g < G; ++g) {
                    grid_weights(g) += (last_grid_weights(g) * C(g, s)) * inv_norm;
                }
            }

            Scalar total_weight = grid_weights.sum();
            grid_weights /= total_weight;

            Scalar err = (grid_weights - last_grid_weights).cwiseAbs().sum();
            if (progress_cb && (step % 10 == 0 || err < tolerance)) {
                progress_cb(step, err);
            }
            if (err < tolerance) {
                break;
            }
        }

        std::vector<Scalar> result(G);
        for (size_t g = 0; g < G; ++g) {
            result[g] = grid_weights(g);
        }
        return result;
    }

    // Collapsed Gibbs sampler
    static std::vector<Scalar> run_collapsed_gibbs(
        const Matrix& C,
        Scalar concentration,
        int64_t chain_length = 2000000,
        int64_t burn_in = 1000000,
        int64_t samples = 100,
        std::function<void(int64_t, int64_t)> progress_cb = nullptr
    ) {
        size_t G = C.rows();
        size_t S = C.cols();

        std::mt19937_64 rng(42);
        std::gamma_distribution<Scalar> init_gamma(concentration, 1.0);
        Vector current_weights(G);
        for (size_t g = 0; g < G; ++g) {
            current_weights(g) = init_gamma(rng);
        }
        current_weights /= current_weights.sum();

        int64_t sample_interval = std::max<int64_t>(1, (chain_length - burn_in) / samples);
        Vector sampled_weights_sum = Vector::Zero(G);
        int64_t samples_drawn = 0;

        for (int64_t step = 1; step <= chain_length; ++step) {
            Vector norm = Vector::Zero(S);
            for (size_t s = 0; s < S; ++s) {
                Scalar sum_s = 0.0;
                for (size_t g = 0; g < G; ++g) {
                    sum_s += current_weights(g) * C(g, s);
                }
                norm(s) = sum_s;
            }

            Vector dir_params = Vector::Constant(G, concentration);
            for (size_t s = 0; s < S; ++s) {
                Scalar inv_norm = (norm(s) > 1e-300) ? (1.0 / norm(s)) : 0.0;
                for (size_t g = 0; g < G; ++g) {
                    dir_params(g) += (current_weights(g) * C(g, s)) * inv_norm;
                }
            }

            for (size_t g = 0; g < G; ++g) {
                std::gamma_distribution<Scalar> gamma_dist(std::max<Scalar>(1e-4, dir_params(g)), 1.0);
                current_weights(g) = gamma_dist(rng);
            }
            Scalar sum_w = current_weights.sum();
            if (sum_w > 0.0) current_weights /= sum_w;

            if (step > burn_in && ((step - burn_in) % sample_interval == 0) && samples_drawn < samples) {
                sampled_weights_sum += current_weights;
                samples_drawn++;
                if (progress_cb) progress_cb(step, samples_drawn);
            }
        }

        if (samples_drawn > 0) {
            sampled_weights_sum /= static_cast<Scalar>(samples_drawn);
        } else {
            sampled_weights_sum = current_weights;
        }

        std::vector<Scalar> result(G);
        for (size_t g = 0; g < G; ++g) {
            result[g] = sampled_weights_sum(g);
        }
        return result;
    }

    // Full Metropolis-Hastings MCMC
    static std::vector<Scalar> run_metropolis_hastings(
        const Matrix& C,
        Scalar concentration,
        int64_t chain_length = 2000000,
        int64_t burn_in = 1000000,
        int64_t samples = 100,
        std::function<void(int64_t, int64_t)> progress_cb = nullptr
    ) {
        size_t G = C.rows();
        size_t S = C.cols();

        std::mt19937_64 rng(42);
        std::uniform_real_distribution<Scalar> unif(0.0, 1.0);
        std::uniform_int_distribution<size_t> rand_grid(0, G - 1);

        // Initial loading from row sums jiggled by +/- 20%
        Vector sum_by_site = C.rowwise().sum();
        Vector weights(G);
        for (size_t g = 0; g < G; ++g) {
            Scalar jitter = 0.8 + 0.4 * unif(rng);
            weights(g) = sum_by_site(g) * jitter;
        }
        weights /= weights.sum();

        Scalar baseline_step = std::max(std::min(0.001, 1.0 / static_cast<Scalar>(S)), 0.0001);

        Vector current_site_likelihoods = weights.transpose() * C;

        int64_t sample_interval = std::max<int64_t>(1, (chain_length - burn_in) / samples);
        Vector sampled_weights_sum = Vector::Zero(G);
        int64_t samples_drawn = 0;
        int64_t accepted_steps = 0;

        for (int64_t step = 1; step <= chain_length; ++step) {
            size_t idx1 = rand_grid(rng);
            size_t idx2 = rand_grid(rng);
            while (idx1 == idx2) {
                idx2 = rand_grid(rng);
            }

            if (step % 50000 == 0) {
                Scalar acc_rate = static_cast<Scalar>(accepted_steps) / static_cast<Scalar>(step);
                if (acc_rate < 0.25) baseline_step /= 1.6;
                else if (acc_rate > 0.5) baseline_step *= 1.6;
            }

            Scalar change = unif(rng) * baseline_step;
            if (weights(idx1) > change) {
                Scalar logLDiff = 0.0;
                for (size_t s = 0; s < S; ++s) {
                    Scalar diff_s = (C(idx2, s) - C(idx1, s)) * change;
                    Scalar new_lik = current_site_likelihoods(s) + diff_s;
                    logLDiff += std::log(std::max(new_lik, 1e-300)) - std::log(std::max(current_site_likelihoods(s), 1e-300));
                }

                Scalar diffPrior = (concentration - 1.0) * (
                    std::log((weights(idx1) - change) / weights(idx1)) +
                    std::log((weights(idx2) + change) / weights(idx2))
                );

                Scalar costOfMove = logLDiff + diffPrior;
                if (costOfMove >= 0.0 || unif(rng) <= std::exp(costOfMove)) {
                    for (size_t s = 0; s < S; ++s) {
                        current_site_likelihoods(s) += (C(idx2, s) - C(idx1, s)) * change;
                    }
                    weights(idx1) -= change;
                    weights(idx2) += change;
                    accepted_steps++;
                }
            }

            if (step > burn_in && ((step - burn_in) % sample_interval == 0) && samples_drawn < samples) {
                sampled_weights_sum += weights;
                samples_drawn++;
                if (progress_cb) progress_cb(step, samples_drawn);
            }
        }

        if (samples_drawn > 0) {
            sampled_weights_sum /= static_cast<Scalar>(samples_drawn);
        } else {
            sampled_weights_sum = weights;
        }

        std::vector<Scalar> result(G);
        for (size_t g = 0; g < G; ++g) {
            result[g] = sampled_weights_sum(g);
        }
        return result;
    }

    // Expected Number of False Positives 95% Confidence Interval
    // Convolves Bernoulli probabilities matching HyPhy 2.5 fubar.ComputeENFP_CI
    static std::pair<int, int> compute_enfp_ci(const std::vector<Scalar>& p_err, Scalar sig_level = 0.05) {
        size_t N = p_err.size();
        if (N == 0) return {0, 0};

        std::vector<Scalar> pdf = {1.0 - p_err[0], p_err[0]};
        for (size_t i = 1; i < N; ++i) {
            std::vector<Scalar> next_pdf(i + 2, 0.0);
            next_pdf[0] = pdf[0] * (1.0 - p_err[i]);
            for (size_t x = 1; x <= i; ++x) {
                next_pdf[x] = pdf[x] * (1.0 - p_err[i]) + pdf[x - 1] * p_err[i];
            }
            next_pdf[i + 1] = pdf[i] * p_err[i];
            pdf = std::move(next_pdf);
        }

        Scalar sum = 0.0;
        int lb = 0;
        for (size_t idx = 0; idx < pdf.size(); ++idx) {
            sum += pdf[idx];
            if (sum >= sig_level / 2.0) {
                lb = static_cast<int>(idx);
                break;
            }
        }

        sum = 0.0;
        int ub = static_cast<int>(pdf.size() - 1);
        for (size_t idx = 0; idx < pdf.size(); ++idx) {
            sum += pdf[idx];
            if (sum >= 1.0 - sig_level / 2.0) {
                ub = static_cast<int>(idx);
                break;
            }
        }
        return {lb, ub};
    }

    // Modern v3 JSON schema export
std::string to_modern_json(
        const FUBARResult& res,
        const std::string& cli_cmd = "hyphy3 fubar",
        const std::string& input_filepath = ""
    ) const;

std::string to_legacy_json(const FUBARResult& res, const std::string& input_filepath = "") const;

    // Unified JSON export: Modern v3.0 schema or legacy, with compatibility blocks
    std::string to_json(
        const FUBARResult& res,
        const std::string& fmt = "modern",
        const std::string& cli_cmd = "hyphy3 fubar",
        const std::string& input_filepath = ""
    ) const {
        if (fmt == "legacy") {
            return to_legacy_json(res, input_filepath);
        }
        nlohmann::json j = nlohmann::json::parse(to_modern_json(res, cli_cmd, input_filepath));
        nlohmann::json legacy = nlohmann::json::parse(to_legacy_json(res, input_filepath));

        // Inject compatibility blocks
        j["MLE"] = legacy["MLE"];
        j["grid"] = legacy["grid"];
        j["posterior"] = legacy["posterior"];
        j["settings"] = legacy["settings"];
        j["fits"] = legacy["fits"];
        j["input"] = legacy["input"];
        j["timers"] = legacy["timers"];
        j["runtime"] = legacy["runtime"];

        return j.dump(2);
    }

    // Print summary panel to console
    static void print_summary(const FUBARResult& res, const FUBARSettings& settings) {
        Panel::print_banner(
            "FUBAR RESULTS SUMMARY",
            "Fast, Unconstrained Bayesian AppRoximation",
            "Posterior distribution of substitution rates over bivariate (alpha, beta) grid",
            "HyPhy 3.0 • Site-level Selection Suite"
        );

        Table summary_tbl;
        summary_tbl.add_column("Metric");
        summary_tbl.add_column("Value");
        summary_tbl.add_column("Notes");

        summary_tbl.add_row({"Inference Method", fubar_method_to_string(settings.method), "Posterior weight optimization"});
        summary_tbl.add_row({"Grid Dimensions", std::to_string(settings.grid_points) + " x " + std::to_string(settings.grid_points), std::to_string(res.grid.size()) + " rate pairs"});
        summary_tbl.add_row({"Dirichlet Prior Alpha", format_float_val(settings.concentration, 2), "Hyperparameter"});
        summary_tbl.add_row({"Tree Scaling Factor", format_float_val(res.best_scaler_alpha, 3) + " / " + format_float_val(res.best_scaler_beta, 3), "Optimal grid scaling"});
        summary_tbl.add_row({"Total Tree Length", format_float_val(res.baseline_tree_length, 3), "Expected substitutions"});
        summary_tbl.add_row({"Positive Selection Threshold", "Prob[alpha < beta] >= " + format_float_val(settings.posterior_threshold, 2), "Credibility cutoff"});
        summary_tbl.add_row({"Diversifying Sites (Pos)", Console::bold(std::to_string(res.positive_sites.size())), "Posterior Prob >= " + format_float_val(settings.posterior_threshold, 2)});
        summary_tbl.add_row({"Purifying Sites (Neg)", Console::bold(std::to_string(res.purifying_sites.size())), "Posterior Prob >= " + format_float_val(settings.posterior_threshold, 2)});
        summary_tbl.add_row({"Expected False Positives", format_float_val(res.expected_false_positives, 2), "95% CI: [" + std::to_string(res.false_positives_ci_95.first) + " - " + std::to_string(res.false_positives_ci_95.second) + "]"});
        summary_tbl.add_row({"Total Runtime", format_float_val(res.time_total_sec, 2) + "s", "Grid: " + format_float_val(res.time_grid_sec, 2) + "s | Post: " + format_float_val(res.time_posterior_sec, 3) + "s"});

        summary_tbl.print();
        std::cout << "\n";

        if (res.positive_sites.empty()) {
            std::cout << Console::muted("  No sites detected under diversifying positive selection at posterior probability >= ") 
                      << settings.posterior_threshold << ".\n\n";
        } else {
            std::cout << Console::bold(Console::brand("  Detected " + std::to_string(res.positive_sites.size()) + 
                         " sites subject to diversifying positive selection (Prob[beta > alpha] >= " + 
                         format_float_val(settings.posterior_threshold, 2) + "):")) << "\n";
            std::cout << "  ";
            for (size_t i = 0; i < res.positive_sites.size(); ++i) {
                std::cout << Console::accent(std::to_string(res.positive_sites[i])) << (i + 1 < res.positive_sites.size() ? ", " : "\n");
            }
            std::cout << Console::muted("  Expected false positives: ") << format_float_val(res.expected_false_positives, 2)
                      << Console::muted(" (95% CI: ") << res.false_positives_ci_95.first << " - " << res.false_positives_ci_95.second << Console::muted(")\n\n");
        }
    }

    // Print table of significant sites
    static void print_results_table(const FUBARResult& res, const FUBARSettings& settings) {
        Table tbl;
        tbl.add_column("Codon");
        tbl.add_column("alpha");
        tbl.add_column("beta");
        tbl.add_column("beta - alpha");
        tbl.add_column("Prob[alpha>beta]");
        tbl.add_column("Prob[alpha<beta]");
        tbl.add_column("BayesFactor");
        tbl.add_column("Selection Mode");

        size_t shown_count = 0;
        for (const auto& sr : res.site_results) {
            bool is_pos = sr.prob_pos >= settings.posterior_threshold;
            bool is_neg = sr.prob_neg >= settings.posterior_threshold;

            if (is_pos || is_neg) {
                shown_count++;
                std::string mode_str;
                if (is_pos) {
                    mode_str = Console::danger("Diversifying (Prob=" + format_float_val(sr.prob_pos, 4) + ")");
                } else {
                    mode_str = Console::info("Purifying (Prob=" + format_float_val(sr.prob_neg, 4) + ")");
                }

                tbl.add_row({
                    std::to_string(sr.site),
                    format_float_val(sr.alpha, 3),
                    format_float_val(sr.beta, 3),
                    format_float_val(sr.beta_minus_alpha, 3),
                    format_float_val(sr.prob_neg, 4),
                    format_float_val(sr.prob_pos, 4),
                    format_float_val(sr.bayes_factor, 2),
                    mode_str
                });
            }
        }

        if (shown_count > 0) {
            std::cout << Console::bold("Significant Sites at Posterior Probability >= " + format_float_val(settings.posterior_threshold, 2) + ":\n");
            tbl.print();
            std::cout << "\n";
        }
    }

private:
};

} // namespace hyphy::analyses
