#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/core/stats.hpp"
#include "hyphy/core/provenance.hpp"
#include "hyphy/core/json_output.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "hyphy/analyses/prime.hpp"
#include "hyphy/analyses/mss.hpp"
#include "nlohmann/json.hpp"

#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hyphy::analyses {

using namespace hyphy::core;
using namespace hyphy::opt;

// ============================================================================
// Mutational Fitness & Evolutionary Trajectory Forecasting Engine
// ============================================================================
// Combines:
//   1. Internal-branch / lineage-partitioned FEL-SRV (alpha_j, beta_j)
//   2. Hierarchical G+S-PRIME Empirical Bayes biophysical exchangeabilities
//      (lambda_{j,k} shrunk toward global G-PRIME prior lambda_k^global)
//   3. Multi-class Synonymous Substitution (MSS / SynREVCodon) symmetric
//      codon-pair exchangeabilities Theta_syn(u, v)
//   4. Exact Inside-Outside marginal internal-node amino-acid posteriors
//      inducing a reversible site-specific stationary preference tilt phi_{j,a}
//
// Produces an exact reversible site-specific S x S codon generator Q_j with
// closed-form stationary distribution pi_j^*(u) proportional to
// pi_F3x4(u) * phi_{j, AA(u)}, enabling exact self-adjoint matrix exponentials
// P_j(Delta_t) = exp(Q_j * tau(Delta_t)), L x 20 DMS fitness matrices, and
// L x S codon rate & multi-step transition log-probability matrices.
// ============================================================================

struct ForecastSettings {
    std::string branches = "Internal";        // "Internal" (trunk default), "All", "Leaves", or tag
    std::string branch_regex = "";
    std::vector<std::string> branch_names = {};

    std::string property_set = "3PROP";       // "3PROP", "5PROP", "LCAP", "Atchley", "4PROP", "2PROP"
    std::string custom_property_file = "";
    Scalar shrinkage_lambda = 0.5;            // Ridge penalty weight gamma on site lambdas
    std::string shrinkage_prior = "global";   // "global" (G-PRIME-centered EB) or "zero"
    Scalar max_bound = 8.0;                   // Soft tanh saturation ceiling M

    std::string syn_model = "SynREVCodon";    // "SynREVCodon", "SynREV", "SynREV2g", "SynREV2", "None"/"MG94", or TSV
    Scalar syn_shrinkage = 0.5;               // Log-normal EB penalty on synonymous codon-pair rates

    Scalar ancestral_exponent = 0.5;          // Exponent eta >= 0 for internal-node AA preference tilt (0 disables)
    Scalar ancestral_pseudocount = 0.05;      // Dirichlet smoothing tau toward F3x4 baseline AA frequencies

    Scalar omega_prior_weight = 0.05;         // Empirical Bayes smoothing weight toward global omega for log(omega_eff)
    Scalar rate_prior_weight = 0.02;          // Small baseline rate floor so invariant sites retain finite-horizon dynamics
    Scalar horizon = 0.05;                    // Forecast evolutionary horizon Delta_t (expected codon substitutions/site)

    std::string reference_taxon = "";         // Optional taxon name defining focal codon/AA per site (empty = trunk consensus)
    Scalar p_value_threshold = 0.10;
};

struct ForecastSiteResult {
    size_t site_index = 0;                    // 0-based site index
    int focal_codon_index = 0;                // Sense codon index u_{0,j} in [0, S)
    std::string focal_codon;                  // 3-letter focal/wildtype codon
    char focal_aa = 'A';                      // 1-letter focal/wildtype amino acid
    int root_codon_index = 0;                 // ML marginal root sense codon index
    std::string root_codon;                   // 3-letter ML marginal root codon
    char root_aa = 'A';                       // 1-letter ML marginal root amino acid

    Scalar alpha = 0.0;                       // Site synonymous rate (SRV)
    Scalar beta = 0.0;                        // Site nonsynonymous rate on tested branches
    Scalar beta_background = 0.0;             // Site nonsynonymous rate on background branches
    Scalar effective_omega = 0.0;             // EB-smoothed omega_j^eff > 0
    Scalar omnibus_lrt = 0.0;                 // S-PRIME omnibus LRT
    Scalar omnibus_p_value = 1.0;
    Scalar omnibus_q_value = 1.0;

    Scalar expected_substitutions_horizon = 0.0; // 1 - P_j(Delta_t)[u_{0,j}, u_{0,j}]
    Scalar stationary_aa_entropy = 0.0;          // Shannon entropy (nats) of pi_{j,AA}^*
    Scalar effective_aa_count = 1.0;             // N_eff = exp(H(pi_{j,AA}^*)) in [1, 20]

    std::vector<Scalar> lambdas;              // Length K: hierarchical EB property lambdas
};

struct ForecastResult {
    ForecastSettings settings;
    std::string effective_branches = "Internal";
    std::string effective_syn_model = "SynREVCodon";
    std::string resolved_reference_taxon = "";

    std::vector<std::string> aa_alphabet;     // 20 canonical IUPAC amino acids ["A", ..., "Y"]
    std::vector<std::string> sense_codons;    // S sense codons in GeneticCode order
    std::vector<std::string> sense_codon_aa;  // S 1-letter translated AAs for each sense codon
    std::vector<std::string> property_names;  // K biophysical property names
    std::vector<AminoAcidProperty> properties;

    size_t total_sites = 0;
    size_t num_sense_codons = 61;
    size_t num_properties = 0;

    Scalar gtr_log_l = 0.0;
    Scalar prime_global_log_l = 0.0;
    Scalar prime_global_aicc = 0.0;
    Scalar global_omega = 1.0;
    std::vector<Scalar> global_lambdas;       // Length K: G-PRIME global prior center

    Scalar mss_baseline_log_l = 0.0;
    Scalar mss_log_l = 0.0;
    Scalar mss_omnibus_lrt = 0.0;
    Scalar mss_omnibus_p_value = 1.0;
    size_t mss_num_classes = 1;
    std::vector<std::string> mss_class_names;
    std::vector<Scalar> mss_class_rates;
    Matrix syn_pair_matrix;                   // (S x S) symmetric synonymous codon-pair multiplier matrix

    Scalar baseline_scale_factor = 1.0;       // Global codon rate matrix flux scale
    Scalar horizon_tau = 0.05;                // Scaled time parameter tau = horizon / baseline_scale_factor

    // Primary ML / Deep-Learning / DMS Target Tensors (Eigen::MatrixXd -> zero-copy NumPy in Python)
    Matrix aa_dms_matrix;                     // (L x 20) log mutational fitness shift relative to focal_aa (0 at focal_aa)
    Matrix aa_stationary_matrix;              // (L x 20) exact marginal stationary AA distribution pi_{j,AA}^*
    Matrix ancestral_aa_matrix;               // (L x 20) exact internal-node marginal AA posterior distribution
    Matrix site_lambda_matrix;                // (L x K)  hierarchical EB property constraints lambda_{j,k}
    Matrix codon_rate_matrix;                 // (L x S)  instantaneous off-diagonal substitution rates Q_j(u_{0,j}, v)
    Matrix codon_log_prob_matrix;             // (L x S)  finite-horizon log transition probabilities log P_j(Delta_t)_{u_{0,j}, v}

    std::vector<ForecastSiteResult> site_results;
    Scalar runtime_seconds = 0.0;

    nlohmann::json to_modern_json(
        const Tree& tree,
        const Alignment& aln,
        const Provenance& prov = {},
        const std::string& input_filepath = ""
    ) const;

    void write_dms_tsv(const std::string& filepath) const;
    void write_codon_tsv(const std::string& filepath) const;
};

class ForecastAnalyzer {
public:
    Tree tree;
    Alignment aln;
    ForecastSettings settings;

    ForecastAnalyzer(Tree t, Alignment a, ForecastSettings s = ForecastSettings{})
        : tree(std::move(t)), aln(std::move(a)), settings(std::move(s)) {}

    static ForecastAnalyzer create(Tree t, Alignment a, ForecastSettings s = ForecastSettings{}) {
        return ForecastAnalyzer(std::move(t), std::move(a), std::move(s));
    }

    // Build exact reversible site-specific S x S codon generator Q_j, stationary distribution pi_star,
    // and finite-horizon transition probability matrix P_horizon = exp(Q_j * tau).
    static void build_reversible_site_codon_model(
        int S,
        const GeneticCode& code,
        const GTRParameters& gtr,
        const Eigen::Matrix<Scalar, 3, 4>& pos_nuc_freqs,
        const Vector& f3x4_codon_freqs,
        const Matrix& syn_pair_matrix,
        const std::vector<Matrix>& codon_prop_diffs,
        Scalar alpha_gen,
        Scalar beta_gen,
        const std::vector<Scalar>& site_lambdas,
        const std::array<Scalar, 20>& phi_aa,
        Scalar max_bound,
        Scalar horizon_tau,
        Matrix& out_Q,
        Vector& out_pi_star,
        Matrix& out_P_horizon
    ) {
        const size_t K = site_lambdas.size();

        // 1. Exact closed-form stationary distribution pi_j^*(u) proportional to pi_F3x4(u) * phi_aa[AA(u)]
        out_pi_star.resize(S);
        Scalar z_sum = 0.0;
        for (int u = 0; u < S; ++u) {
            int aa_u = aa_char_to_index(code.sense_codon_aa[u]);
            Scalar phi_u = (aa_u >= 0 && aa_u < 20) ? std::max<Scalar>(phi_aa[aa_u], 1e-12) : 1.0;
            Scalar base_pi = std::max<Scalar>(f3x4_codon_freqs(u), 1e-12);
            out_pi_star(u) = base_pi * phi_u;
            z_sum += out_pi_star(u);
        }
        if (z_sum > 0.0) {
            out_pi_star /= z_sum;
        } else {
            out_pi_star.setConstant(1.0 / static_cast<Scalar>(S));
        }

        Vector sqrt_pi(S), inv_sqrt_pi(S);
        for (int u = 0; u < S; ++u) {
            sqrt_pi(u) = std::sqrt(std::max<Scalar>(out_pi_star(u), 1e-15));
            inv_sqrt_pi(u) = 1.0 / sqrt_pi(u);
        }

        const Scalar r_nuc[4][4] = {
            {0.0, gtr.theta_AC, 1.0, gtr.theta_AT},
            {gtr.theta_AC, 0.0, gtr.theta_CG, gtr.theta_CT},
            {1.0, gtr.theta_CG, 0.0, gtr.theta_GT},
            {gtr.theta_AT, gtr.theta_CT, gtr.theta_GT, 0.0}
        };

        out_Q.resize(S, S);
        out_Q.setZero();

        for (int u = 0; u < S; ++u) {
            int aa_u = aa_char_to_index(code.sense_codon_aa[u]);
            Scalar phi_u = (aa_u >= 0 && aa_u < 20) ? std::max<Scalar>(phi_aa[aa_u], 1e-12) : 1.0;

            for (int v = 0; v < S; ++v) {
                if (u == v) continue;
                const auto& diff = code.diff_matrix[u][v];
                if (diff.num_diffs != 1) continue;

                Scalar nuc_rate = r_nuc[diff.diff_from[0]][diff.diff_to[0]];
                Scalar target_nuc_freq = pos_nuc_freqs(diff.diff_pos[0], diff.diff_to[0]);

                Scalar rate_mod = 0.0;
                if (diff.synonymous) {
                    Scalar syn_mult = (syn_pair_matrix.rows() == S && syn_pair_matrix.cols() == S)
                        ? syn_pair_matrix(u, v)
                        : 1.0;
                    rate_mod = alpha_gen * syn_mult;
                } else {
                    int aa_v = aa_char_to_index(code.sense_codon_aa[v]);
                    Scalar phi_v = (aa_v >= 0 && aa_v < 20) ? std::max<Scalar>(phi_aa[aa_v], 1e-12) : 1.0;

                    Scalar eta = 0.0;
                    for (size_t k = 0; k < K; ++k) {
                        eta += site_lambdas[k] * codon_prop_diffs[k](u, v);
                    }
                    Scalar eta_eff = max_bound * std::tanh(eta / max_bound);
                    Scalar psi = std::exp(-eta_eff);
                    Scalar dir_tilt = std::sqrt(phi_v / phi_u);
                    rate_mod = beta_gen * psi * dir_tilt;
                }

                out_Q(u, v) = rate_mod * nuc_rate * target_nuc_freq;
            }
        }

        for (int u = 0; u < S; ++u) {
            out_Q(u, u) = -out_Q.row(u).sum();
        }

        // Exact self-adjoint eigendecomposition via detailed balance similarity transform
        Matrix S_mat = sqrt_pi.asDiagonal() * out_Q * inv_sqrt_pi.asDiagonal();
        S_mat = 0.5 * (S_mat + S_mat.transpose());

        Eigen::SelfAdjointEigenSolver<Matrix> solver(S_mat);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("ForecastAnalyzer: site codon generator eigendecomposition failed");
        }

        Vector evals = solver.eigenvalues().cwiseMin(0.0);
        const Matrix& evecs = solver.eigenvectors();
        Matrix U_pi = inv_sqrt_pi.asDiagonal() * evecs;
        Matrix V_pi = evecs.transpose() * sqrt_pi.asDiagonal();

        Vector exp_lam = (evals * std::max<Scalar>(horizon_tau, 0.0)).array().exp();
        out_P_horizon = ((U_pi * exp_lam.asDiagonal()) * V_pi).cwiseMax(0.0);

        // Enforce exact row stochasticity
        for (int u = 0; u < S; ++u) {
            Scalar rsum = out_P_horizon.row(u).sum();
            if (rsum > 0.0) {
                out_P_horizon.row(u) /= rsum;
            } else {
                out_P_horizon.row(u).setZero();
                out_P_horizon(u, u) = 1.0;
            }
        }
    }

    ForecastResult run(bool show_progress = false, bool force_progress = false) {
        auto t_start = std::chrono::high_resolution_clock::now();
        const auto& code = *(aln.code ? aln.code : GeneticCode::universal());
        const int S = code.num_sense_codons;
        const size_t L = aln.num_codons;
        const size_t num_patterns = aln.patterns.size();

        if (L == 0 || aln.num_taxa < 3) {
            throw std::runtime_error("ForecastAnalyzer requires a non-empty codon alignment with at least 3 taxa.");
        }

        LikelihoodEngine::reconcile_tree_leaf_names(tree, aln);

        // 1. Resolve branch partition (auto-fallback from Internal -> All on 3-taxon star trees)
        std::string eff_branches = settings.branches.empty() ? "Internal" : settings.branches;
        if (eff_branches == "Internal" && settings.branch_names.empty() && settings.branch_regex.empty()) {
            size_t internal_branch_count = 0;
            for (const auto& node : tree.nodes) {
                if (node.id != tree.root_id && !node.is_leaf) {
                    internal_branch_count++;
                }
            }
            if (internal_branch_count == 0) {
                eff_branches = "All";
            }
        }

        // Resolve optional focal reference taxon
        int ref_taxon_idx = -1;
        std::string resolved_ref_name;
        if (!settings.reference_taxon.empty()) {
            auto it = aln.taxon_to_index.find(settings.reference_taxon);
            if (it != aln.taxon_to_index.end()) {
                ref_taxon_idx = static_cast<int>(it->second);
                resolved_ref_name = aln.taxon_names[ref_taxon_idx];
            } else {
                std::string low = settings.reference_taxon;
                std::transform(low.begin(), low.end(), low.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                auto it_low = aln.lower_taxon_to_index.find(low);
                if (it_low != aln.lower_taxon_to_index.end()) {
                    ref_taxon_idx = static_cast<int>(it_low->second);
                    resolved_ref_name = aln.taxon_names[ref_taxon_idx];
                } else {
                    throw std::runtime_error("ForecastAnalyzer: reference taxon '" + settings.reference_taxon + "' not found in alignment.");
                }
            }
        }

        // 2. Fit Hierarchical G+S-PRIME (GTR + MG94 + G-PRIME + FEL-SRV + EB S-PRIME)
        PRIMESettings p_settings;
        p_settings.p_value_threshold = settings.p_value_threshold;
        p_settings.stage2_p_value = settings.p_value_threshold;
        p_settings.property_set = settings.property_set;
        p_settings.custom_property_file = settings.custom_property_file;
        p_settings.shrinkage_lambda = settings.shrinkage_lambda;
        p_settings.shrinkage_prior = settings.shrinkage_prior;
        p_settings.max_bound = settings.max_bound;
        p_settings.fit_global = true;
        p_settings.test_branches = eff_branches;
        p_settings.test_branch_regex = settings.branch_regex;
        p_settings.test_branch_names = settings.branch_names;

        PRIMEAnalyzer prime = PRIMEAnalyzer::create_and_fit(tree, aln, p_settings);
        PRIMEResult prime_res = prime.run(nullptr, show_progress, force_progress);
        tree = prime.tree; // Keep fitted codon-scaled branch lengths

        const size_t K = prime.K;
        const Scalar global_omega = (prime.global_fit.alpha > 1e-12)
            ? (prime.global_fit.beta / prime.global_fit.alpha)
            : std::max<Scalar>(prime.base_params.beta, 0.05);
        std::vector<Scalar> global_lambdas = (prime.global_fit.lambdas.size() == K)
            ? prime.global_fit.lambdas
            : std::vector<Scalar>(K, 0.0);

        // Compute global baseline scale factor for physical horizon calibration
        PRIMEMatrix pm_global;
        pm_global.update(
            1.0, global_omega, global_lambdas, prime.codon_prop_diffs,
            prime.gtr_rates, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4,
            code, settings.max_bound
        );
        const Scalar base_scale = std::max<Scalar>(pm_global.scale_factor, 1e-6);
        const Scalar horizon_tau = std::max<Scalar>(settings.horizon, 0.0) / base_scale;

        // 3. Fit Multi-class Synonymous Substitution (MSS / SynREVCodon) symmetric codon-pair matrix
        Matrix syn_pair_mat = Matrix::Ones(S, S);
        std::string syn_up = MSSPartition::to_upper(MSSPartition::trim(settings.syn_model));
        std::string eff_syn_model = "MG94";
        Scalar mss_base_ll = 0.0;
        Scalar mss_ll = 0.0;
        Scalar mss_lrt = 0.0;
        Scalar mss_pval = 1.0;
        size_t mss_nclasses = 1;
        std::vector<std::string> mss_cnames = {"NEUTRAL"};
        std::vector<Scalar> mss_crates = {1.0};

        if (!syn_up.empty() && syn_up != "NONE" && syn_up != "MG94" && syn_up != "UNIFORM") {
            MSSPartition part = MSSPartition::build(settings.syn_model, code);
            eff_syn_model = part.name;
            mss_nclasses = part.num_classes();
            mss_cnames = part.class_names;

            MSSAnalyzer mss(aln.code ? aln.code : GeneticCode::universal());
            MSSAnalyzer::GeneContext gctx;
            gctx.filepath = "";
            gctx.aln = aln;
            gctx.tree = prime.tree;
            gctx.leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(gctx.tree, gctx.aln);
            gctx.gtr_params = prime.gtr_rates;
            gctx.gtr_log_l = prime.gtr_log_l;
            gctx.mg94_params = prime.base_params;
            gctx.mg94_omega = prime.base_params.beta;
            gctx.omega = prime.base_params.beta;
            gctx.tree_scaler = 1.0;

            MG94Matrix mg_base;
            mg_base.update(gctx.mg94_params, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4, code);
            gctx.mg94_scale_factor = mg_base.scale_factor;
            gctx.branch_tau.assign(gctx.tree.num_nodes(), 0.0);
            for (const auto& node : gctx.tree.nodes) {
                if (node.id != gctx.tree.root_id) {
                    // prime.tree branch lengths are already in codon transition_matrix tau units
                    gctx.branch_tau[node.id] = node.branch_length;
                }
            }
            mss.genes.push_back(std::move(gctx));

            std::vector<Scalar> init_ones(mss_nclasses, 1.0);
            std::vector<Scalar> g_scalers = {1.0};
            std::vector<Scalar> g_omegas = {prime.base_params.beta};
            mss_base_ll = mss.evaluate_composite_log_l(part, init_ones, &g_scalers, &g_omegas);

            if (part.is_empirical_fixed) {
                mss_crates.assign(mss_nclasses, 1.0);
                std::vector<Scalar> flat = part.expand_to_matrix(S, mss_crates);
                for (int u = 0; u < S; ++u) {
                    for (int v = 0; v < S; ++v) syn_pair_mat(u, v) = flat[u * S + v];
                }
                mss_ll = mss.evaluate_composite_log_l(part, mss_crates, &g_scalers, &g_omegas);
            } else if (mss_nclasses > 1) {
                std::vector<Scalar> pair_flux = mss.compute_pair_flux_weights(part);
                std::vector<Scalar> class_flux(mss_nclasses, 0.0);
                for (size_t p = 0; p < part.pairs.size(); ++p) {
                    int cid = part.pairs[p].class_id;
                    if (cid >= 0 && static_cast<size_t>(cid) < mss_nclasses) {
                        class_flux[cid] += pair_flux[p];
                    }
                }
                int ref_k = static_cast<int>(std::max_element(class_flux.begin(), class_flux.end()) - class_flux.begin());
                mss_crates.assign(mss_nclasses, 1.0);

                const int n_cycles = (mss_nclasses > 30) ? 2 : 3;
                const int brent_iters = (mss_nclasses > 30) ? 16 : 24;
                const Scalar pen_w = std::max<Scalar>(settings.syn_shrinkage, 0.0);

                for (int cyc = 0; cyc < n_cycles; ++cyc) {
                    for (size_t k = 0; k < mss_nclasses; ++k) {
                        if (static_cast<int>(k) == ref_k || class_flux[k] <= 1e-12) continue;
                        Scalar init_r = mss_crates[k];
                        auto obj_k = [&](Scalar r_val) -> Scalar {
                            std::vector<Scalar> cand = mss_crates;
                            cand[k] = r_val;
                            MSSAnalyzer::normalize_class_rates(cand, MSSNormalization::Flux, -1, class_flux);
                            Scalar ll = mss.evaluate_composite_log_l(part, cand, &g_scalers, &g_omegas);
                            Scalar pen = 0.0;
                            if (pen_w > 0.0) {
                                for (size_t j = 0; j < mss_nclasses; ++j) {
                                    Scalar lr = std::log(std::max<Scalar>(cand[j], 1e-5));
                                    pen += 0.5 * pen_w * lr * lr;
                                }
                            }
                            return -(ll - pen);
                        };
                        Scalar lo = std::max<Scalar>(0.05, init_r * 0.2);
                        Scalar hi = std::min<Scalar>(20.0, std::max<Scalar>(5.0, init_r * 5.0));
                        auto [best_r, _] = Brent1D::minimize(obj_k, lo, init_r, hi, 2e-3, brent_iters);
                        mss_crates[k] = best_r;
                        MSSAnalyzer::normalize_class_rates(mss_crates, MSSNormalization::Flux, -1, class_flux);
                    }
                }

                mss_ll = mss.evaluate_composite_log_l(part, mss_crates, &g_scalers, &g_omegas);
                if (mss_ll < mss_base_ll && pen_w == 0.0) {
                    mss_crates.assign(mss_nclasses, 1.0);
                    mss_ll = mss_base_ll;
                }
                mss_lrt = std::max<Scalar>(0.0, 2.0 * (mss_ll - mss_base_ll));
                mss_pval = stats::chi2_p_value(mss_lrt, static_cast<Scalar>(mss_nclasses - 1));

                std::vector<Scalar> flat = part.expand_to_matrix(S, mss_crates);
                for (int u = 0; u < S; ++u) {
                    for (int v = 0; v < S; ++v) {
                        syn_pair_mat(u, v) = flat[u * S + v];
                    }
                }
            }
        }

        // 4. Baseline F3x4 Amino-Acid Equilibrium Frequencies pi_AA(b)
        std::array<Scalar, 20> f3x4_aa_freqs{};
        f3x4_aa_freqs.fill(0.0);
        for (int u = 0; u < S; ++u) {
            int aa_idx = aa_char_to_index(code.sense_codon_aa[u]);
            if (aa_idx >= 0 && aa_idx < 20) {
                f3x4_aa_freqs[aa_idx] += aln.codon_frequencies_f3x4(u);
            }
        }
        Scalar aa_sum = 0.0;
        for (int a = 0; a < 20; ++a) {
            f3x4_aa_freqs[a] = std::max<Scalar>(f3x4_aa_freqs[a], 1e-8);
            aa_sum += f3x4_aa_freqs[a];
        }
        for (int a = 0; a < 20; ++a) f3x4_aa_freqs[a] /= aa_sum;

        // 5. Exact Inside-Outside Marginal Internal-Node Posteriors per Unique Site Pattern
        const size_t num_nodes = tree.num_nodes();
        const bool has_bg = (prime.num_background > 0 && prime.num_tested > 0);
        const std::vector<size_t> leaf_to_taxon = LikelihoodEngine::build_leaf_to_taxon_map(tree, aln);

        // Representative site index for each pattern so we read pattern-fitted parameters from prime_res
        std::vector<size_t> pattern_rep_site(num_patterns, 0);
        for (size_t p = 0; p < num_patterns; ++p) {
            if (!aln.patterns[p].sites.empty()) {
                pattern_rep_site[p] = aln.patterns[p].sites[0];
            }
        }

        std::vector<Vector> pat_anc_codon_probs(num_patterns, Vector::Zero(S));
        std::vector<std::array<Scalar, 20>> pat_anc_aa_probs(num_patterns);
        std::vector<int> pat_root_codon(num_patterns, 0);

        #pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < num_patterns; ++p) {
            const size_t rep_site = pattern_rep_site[p];
            const auto& psr = prime_res.site_results[rep_site];

            Scalar a_pat = std::max<Scalar>(psr.alpha, 1e-4);
            Scalar b_pat = std::max<Scalar>(psr.beta, 1e-4 * global_omega);
            Scalar b_bg_pat = std::max<Scalar>(psr.beta_nuisance, 1e-4 * global_omega);
            std::vector<Scalar> lams_pat(K, 0.0);
            for (size_t k = 0; k < K && k < psr.properties.size(); ++k) {
                lams_pat[k] = psr.properties[k].lambda;
            }

            PRIMEMatrix pm_pat;
            pm_pat.update(
                a_pat, b_pat, lams_pat, prime.codon_prop_diffs,
                prime.gtr_rates, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4,
                code, settings.max_bound
            );

            PRIMEMatrix pm_bg_pat;
            if (has_bg) {
                std::vector<Scalar> zero_lams(K, 0.0);
                pm_bg_pat.update(
                    a_pat, b_bg_pat, zero_lams, prime.codon_prop_diffs,
                    prime.gtr_rates, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4,
                    code, settings.max_bound
                );
            }

            std::vector<Matrix> P_br(num_nodes, Matrix::Identity(S, S));
            for (size_t i = 0; i < num_nodes; ++i) {
                if (tree.nodes[i].id != tree.root_id) {
                    // When psr.alpha == psr.beta == 0 (invariant site), use tiny positive branch scaling
                    // so invariant codon identity propagates sharply while ties break by F3x4 prior
                    Scalar bl = (psr.alpha == 0.0 && psr.beta == 0.0)
                        ? (1e-4 * tree.nodes[i].branch_length)
                        : tree.nodes[i].branch_length;
                    if (prime.is_test_branch[i] || !has_bg) {
                        P_br[i] = pm_pat.transition_matrix(bl);
                    } else {
                        P_br[i] = pm_bg_pat.transition_matrix(bl);
                    }
                }
            }

            // Upward (Inside) pass with stored child messages M_child[cid] = P_br[cid] * L_in[cid]
            std::vector<Vector> L_in(num_nodes, Vector::Ones(S));
            std::vector<Vector> M_child(num_nodes, Vector::Ones(S));

            for (int32_t nid : tree.post_order) {
                const auto& node = tree.nodes[nid];
                if (node.is_leaf) {
                    size_t t_idx = leaf_to_taxon[nid];
                    LikelihoodEngine::init_leaf_likelihood(aln.patterns[p], t_idx, S, L_in[nid]);
                } else {
                    Vector prod = Vector::Ones(S);
                    for (int32_t cid : node.children) {
                        prod.array() *= M_child[cid].array();
                        Scalar mx = prod.maxCoeff();
                        if (mx > 0.0 && mx < 1e-100) prod /= mx;
                    }
                    Scalar mx = prod.maxCoeff();
                    if (mx > 0.0) prod /= mx;
                    L_in[nid] = prod;
                }

                if (nid != tree.root_id) {
                    M_child[nid] = P_br[nid] * L_in[nid];
                    Scalar mx = M_child[nid].maxCoeff();
                    if (mx > 0.0) M_child[nid] /= mx;
                }
            }

            // Downward (Outside) pass via division-free prefix/suffix sibling products
            std::vector<Vector> L_out(num_nodes, Vector::Ones(S));
            L_out[tree.root_id] = pm_pat.pi;

            for (auto it = tree.post_order.rbegin(); it != tree.post_order.rend(); ++it) {
                int32_t nid = *it;
                const auto& node = tree.nodes[nid];
                if (node.is_leaf || node.children.empty()) continue;

                const size_t deg = node.children.size();
                std::vector<Vector> prefix(deg + 1, Vector::Ones(S));
                std::vector<Vector> suffix(deg + 1, Vector::Ones(S));

                for (size_t r = 0; r < deg; ++r) {
                    prefix[r + 1] = prefix[r].cwiseProduct(M_child[node.children[r]]);
                    Scalar mx = prefix[r + 1].maxCoeff();
                    if (mx > 0.0) prefix[r + 1] /= mx;
                }
                for (size_t r = deg; r > 0; --r) {
                    suffix[r - 1] = suffix[r].cwiseProduct(M_child[node.children[r - 1]]);
                    Scalar mx = suffix[r - 1].maxCoeff();
                    if (mx > 0.0) suffix[r - 1] /= mx;
                }

                for (size_t r = 0; r < deg; ++r) {
                    int32_t cid = node.children[r];
                    Vector cavity = L_out[nid].cwiseProduct(prefix[r]).cwiseProduct(suffix[r + 1]);
                    Scalar mx_c = cavity.maxCoeff();
                    if (mx_c > 0.0) cavity /= mx_c;

                    L_out[cid] = P_br[cid].transpose() * cavity;
                    Scalar mx_o = L_out[cid].maxCoeff();
                    if (mx_o > 0.0) L_out[cid] /= mx_o;
                }
            }

            // Marginal posteriors at root and across internal nodes on the tested partition
            Vector root_post = L_out[tree.root_id].cwiseProduct(L_in[tree.root_id]);
            Scalar r_sum = root_post.sum();
            if (r_sum > 0.0) {
                root_post /= r_sum;
            } else {
                root_post = aln.codon_frequencies_f3x4;
            }
            Eigen::Index best_root_idx = 0;
            root_post.maxCoeff(&best_root_idx);
            pat_root_codon[p] = static_cast<int>(best_root_idx);

            Vector anc_codon_accum = Vector::Zero(S);
            size_t anc_node_count = 0;
            for (size_t nid = 0; nid < num_nodes; ++nid) {
                const auto& node = tree.nodes[nid];
                if (node.is_leaf) continue;
                if (nid != static_cast<size_t>(tree.root_id) && has_bg && !prime.is_test_branch[nid]) {
                    continue;
                }
                Vector post = L_out[nid].cwiseProduct(L_in[nid]);
                Scalar psum = post.sum();
                if (psum > 0.0) {
                    anc_codon_accum += (post / psum);
                    anc_node_count++;
                }
            }

            if (anc_node_count > 0) {
                anc_codon_accum /= static_cast<Scalar>(anc_node_count);
            } else {
                anc_codon_accum = root_post;
            }
            pat_anc_codon_probs[p] = anc_codon_accum;

            std::array<Scalar, 20> aa_p{};
            aa_p.fill(0.0);
            for (int u = 0; u < S; ++u) {
                int aa_idx = aa_char_to_index(code.sense_codon_aa[u]);
                if (aa_idx >= 0 && aa_idx < 20) {
                    aa_p[aa_idx] += anc_codon_accum(u);
                }
            }
            Scalar aa_tot = 0.0;
            for (int a = 0; a < 20; ++a) aa_tot += aa_p[a];
            if (aa_tot > 0.0) {
                for (int a = 0; a < 20; ++a) aa_p[a] /= aa_tot;
            } else {
                aa_p = f3x4_aa_freqs;
            }
            pat_anc_aa_probs[p] = aa_p;
        }

        // 6. Assemble Site-Specific Reversible Codon Generators, DMS Matrices & Codon Horizon Forecasts
        ForecastResult out;
        out.settings = settings;
        out.effective_branches = eff_branches;
        out.effective_syn_model = eff_syn_model;
        out.resolved_reference_taxon = resolved_ref_name;

        out.aa_alphabet.resize(20);
        for (int a = 0; a < 20; ++a) out.aa_alphabet[a] = std::string(1, aa_index_to_char(a));
        out.sense_codons = code.sense_codons;
        out.sense_codon_aa.resize(S);
        for (int u = 0; u < S; ++u) out.sense_codon_aa[u] = std::string(1, code.sense_codon_aa[u]);

        out.properties = prime.properties;
        out.property_names.resize(K);
        for (size_t k = 0; k < K; ++k) out.property_names[k] = prime.properties[k].name;

        out.total_sites = L;
        out.num_sense_codons = static_cast<size_t>(S);
        out.num_properties = K;

        out.gtr_log_l = prime.gtr_log_l;
        out.prime_global_log_l = prime.global_fit.log_likelihood;
        out.prime_global_aicc = prime.global_fit.aicc;
        out.global_omega = global_omega;
        out.global_lambdas = global_lambdas;

        out.mss_baseline_log_l = mss_base_ll;
        out.mss_log_l = mss_ll;
        out.mss_omnibus_lrt = mss_lrt;
        out.mss_omnibus_p_value = mss_pval;
        out.mss_num_classes = mss_nclasses;
        out.mss_class_names = std::move(mss_cnames);
        out.mss_class_rates = std::move(mss_crates);
        out.syn_pair_matrix = syn_pair_mat;

        out.baseline_scale_factor = base_scale;
        out.horizon_tau = horizon_tau;

        out.aa_dms_matrix = Matrix::Zero(L, 20);
        out.aa_stationary_matrix = Matrix::Zero(L, 20);
        out.ancestral_aa_matrix = Matrix::Zero(L, 20);
        out.site_lambda_matrix = Matrix::Zero(L, K);
        out.codon_rate_matrix = Matrix::Zero(L, S);
        out.codon_log_prob_matrix = Matrix::Zero(L, S);
        out.site_results.resize(L);

        const Scalar eps_omega = std::max<Scalar>(settings.omega_prior_weight, 1e-6);
        const Scalar eps_rate = std::max<Scalar>(settings.rate_prior_weight, 1e-6);
        const Scalar tau_pseudo = std::max<Scalar>(settings.ancestral_pseudocount, 1e-6);
        const Scalar eta_exp = std::max<Scalar>(settings.ancestral_exponent, 0.0);

        #pragma omp parallel for schedule(static)
        for (size_t site = 0; site < L; ++site) {
            const size_t p_idx = aln.site_to_pattern[site];
            const auto& psr = prime_res.site_results[site];

            ForecastSiteResult fsr;
            fsr.site_index = site;
            fsr.alpha = psr.alpha;
            fsr.beta = psr.beta;
            fsr.beta_background = psr.beta_nuisance;
            fsr.omnibus_lrt = psr.omnibus_lrt;
            fsr.omnibus_p_value = psr.omnibus_p_value;
            fsr.omnibus_q_value = psr.omnibus_q_value;

            fsr.effective_omega = (psr.beta + eps_omega * global_omega) / (psr.alpha + eps_omega);

            fsr.lambdas.resize(K, 0.0);
            for (size_t k = 0; k < K && k < psr.properties.size(); ++k) {
                fsr.lambdas[k] = psr.properties[k].lambda;
                out.site_lambda_matrix(site, k) = fsr.lambdas[k];
            }

            // Root & focal codon resolution
            int root_cod = pat_root_codon[p_idx];
            fsr.root_codon_index = root_cod;
            fsr.root_codon = code.sense_codons[root_cod];
            fsr.root_aa = code.sense_codon_aa[root_cod];

            int focal_cod = -1;
            if (ref_taxon_idx >= 0 && static_cast<size_t>(ref_taxon_idx) < aln.sequences.size()) {
                int8_t obs = aln.sequences[ref_taxon_idx].codons[site];
                if (obs >= 0 && obs < S) {
                    focal_cod = static_cast<int>(obs);
                }
            }
            if (focal_cod < 0) {
                Eigen::Index best_anc = 0;
                pat_anc_codon_probs[p_idx].maxCoeff(&best_anc);
                focal_cod = static_cast<int>(best_anc);
            }
            fsr.focal_codon_index = focal_cod;
            fsr.focal_codon = code.sense_codons[focal_cod];
            fsr.focal_aa = code.sense_codon_aa[focal_cod];
            const int focal_aa_idx = std::max(0, aa_char_to_index(fsr.focal_aa));

            // Ancestral AA probabilities & reversible directional preference tilt phi_{j,a}
            const auto& anc_aa = pat_anc_aa_probs[p_idx];
            std::array<Scalar, 20> phi_aa{};
            for (int a = 0; a < 20; ++a) {
                out.ancestral_aa_matrix(site, a) = anc_aa[a];
                if (eta_exp > 0.0) {
                    Scalar smoothed = (anc_aa[a] + tau_pseudo * f3x4_aa_freqs[a]) / (1.0 + tau_pseudo);
                    Scalar ratio = smoothed / std::max<Scalar>(f3x4_aa_freqs[a], 1e-8);
                    phi_aa[a] = std::pow(std::clamp<Scalar>(ratio, 1e-6, 1e6), eta_exp);
                } else {
                    phi_aa[a] = 1.0;
                }
            }

            // Build exact reversible site codon generator Q_j, stationary pi_j^*, and P_j(Delta_t)
            const Scalar alpha_gen = psr.alpha + eps_rate;
            const Scalar beta_gen = psr.beta + eps_rate * global_omega;

            Matrix Q_site;
            Vector pi_star;
            Matrix P_hor;
            build_reversible_site_codon_model(
                S, code, prime.gtr_rates, aln.pos_nuc_frequencies, aln.codon_frequencies_f3x4,
                syn_pair_mat, prime.codon_prop_diffs, alpha_gen, beta_gen, fsr.lambdas,
                phi_aa, settings.max_bound, horizon_tau, Q_site, pi_star, P_hor
            );

            // Marginal stationary AA distribution & entropy
            std::array<Scalar, 20> pi_aa_star{};
            pi_aa_star.fill(0.0);
            for (int u = 0; u < S; ++u) {
                int a_idx = aa_char_to_index(code.sense_codon_aa[u]);
                if (a_idx >= 0 && a_idx < 20) {
                    pi_aa_star[a_idx] += pi_star(u);
                }
            }
            Scalar H_aa = 0.0;
            for (int a = 0; a < 20; ++a) {
                out.aa_stationary_matrix(site, a) = pi_aa_star[a];
                if (pi_aa_star[a] > 1e-15) {
                    H_aa -= pi_aa_star[a] * std::log(pi_aa_star[a]);
                }
            }
            fsr.stationary_aa_entropy = H_aa;
            fsr.effective_aa_count = std::exp(H_aa);

            // 20-AA DMS mutational log-fitness matrix relative to focal_aa
            const Scalar log_omega_eff = std::log(std::max<Scalar>(fsr.effective_omega, 1e-12));
            const Scalar phi_focal = std::max<Scalar>(phi_aa[focal_aa_idx], 1e-12);

            for (int b = 0; b < 20; ++b) {
                if (b == focal_aa_idx) {
                    out.aa_dms_matrix(site, b) = 0.0;
                    continue;
                }
                Scalar eta = 0.0;
                for (size_t k = 0; k < K; ++k) {
                    Scalar dz = std::abs(prime.properties[k].z_scores[focal_aa_idx] - prime.properties[k].z_scores[b]);
                    eta += fsr.lambdas[k] * dz;
                }
                Scalar log_psi = -settings.max_bound * std::tanh(eta / settings.max_bound);
                Scalar log_dir = 0.5 * std::log(std::max<Scalar>(phi_aa[b], 1e-12) / phi_focal);
                out.aa_dms_matrix(site, b) = log_omega_eff + log_psi + log_dir;
            }

            // 61-codon instantaneous rates Q_site(focal_cod, v) & finite-horizon log P_hor(focal_cod, v)
            for (int v = 0; v < S; ++v) {
                out.codon_rate_matrix(site, v) = (v == focal_cod) ? 0.0 : std::max<Scalar>(0.0, Q_site(focal_cod, v));
                Scalar pv = std::max<Scalar>(P_hor(focal_cod, v), 1e-300);
                out.codon_log_prob_matrix(site, v) = std::log(pv);
            }
            fsr.expected_substitutions_horizon = std::clamp<Scalar>(1.0 - P_hor(focal_cod, focal_cod), 0.0, 1.0);

            out.site_results[site] = std::move(fsr);
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        out.runtime_seconds = std::chrono::duration<Scalar>(t_end - t_start).count();
        return out;
    }
};

} // namespace hyphy::analyses
