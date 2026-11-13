#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/forecast.hpp"

#include <cmath>
#include <fstream>
#include <cstdio>

using namespace hyphy::core;
using namespace hyphy::analyses;

inline std::string resolve_test_path(const std::string& path) {
    if (std::ifstream(path).good()) return path;
    if (std::ifstream("../" + path).good()) return "../" + path;
    if (std::ifstream("../../" + path).good()) return "../../" + path;
    return path;
}

TEST_CASE("FORECAST: Exact Detailed Balance, Stationary Convergence, & Row-Stochasticity") {
    auto code = GeneticCode::universal();
    const int S = code->num_sense_codons;

    Eigen::Matrix<Scalar, 3, 4> pos_nuc_freqs;
    pos_nuc_freqs << 0.30, 0.20, 0.28, 0.22,
                     0.26, 0.24, 0.19, 0.31,
                     0.18, 0.32, 0.30, 0.20;

    Vector codon_freqs(S);
    Scalar csum = 0.0;
    for (int u = 0; u < S; ++u) {
        const auto& s = code->sense_codons[u];
        auto nuc_i = [](char c) {
            return (c == 'A') ? 0 : (c == 'C') ? 1 : (c == 'G') ? 2 : 3;
        };
        codon_freqs(u) = pos_nuc_freqs(0, nuc_i(s[0])) *
                         pos_nuc_freqs(1, nuc_i(s[1])) *
                         pos_nuc_freqs(2, nuc_i(s[2]));
        csum += codon_freqs(u);
    }
    codon_freqs /= csum;

    GTRParameters gtr;
    gtr.theta_AC = 0.75;
    gtr.theta_AT = 0.45;
    gtr.theta_CG = 0.60;
    gtr.theta_CT = 1.35;
    gtr.theta_GT = 0.55;

    // Non-uniform symmetric SynREVCodon matrix
    MSSPartition part = MSSPartition::build("SynREVCodon", *code);
    std::vector<Scalar> syn_rates(part.num_classes(), 1.0);
    for (size_t k = 0; k < syn_rates.size(); ++k) {
        syn_rates[k] = 0.4 + 0.025 * static_cast<Scalar>(k % 40);
    }
    std::vector<Scalar> flat = part.expand_to_matrix(S, syn_rates);
    Matrix syn_mat(S, S);
    for (int u = 0; u < S; ++u) {
        for (int v = 0; v < S; ++v) syn_mat(u, v) = flat[u * S + v];
    }

    auto props = PropertyDatabase::get_3prop_properties();
    const size_t K = props.size();
    std::vector<Matrix> prop_diffs(K, Matrix::Zero(S, S));
    for (size_t k = 0; k < K; ++k) {
        for (int u = 0; u < S; ++u) {
            int au = aa_char_to_index(code->sense_codon_aa[u]);
            for (int v = 0; v < S; ++v) {
                int av = aa_char_to_index(code->sense_codon_aa[v]);
                prop_diffs[k](u, v) = std::abs(props[k].z_scores[au] - props[k].z_scores[av]);
            }
        }
    }

    // Non-uniform directional AA preference tilt phi_aa
    std::array<Scalar, 20> phi_aa{};
    for (int a = 0; a < 20; ++a) {
        phi_aa[a] = 0.25 + 0.15 * static_cast<Scalar>((a * 7 + 3) % 13);
    }

    std::vector<Scalar> site_lambdas = {1.8, -0.9, 1.2};
    Matrix Q, P_short, P_long;
    Vector pi_star;

    ForecastAnalyzer::build_reversible_site_codon_model(
        S, *code, gtr, pos_nuc_freqs, codon_freqs, syn_mat, prop_diffs,
        1.1, 0.35, site_lambdas, phi_aa, 8.0, 0.10, Q, pi_star, P_short
    );

    // 1. Exact Detailed Balance pi_star(u)*Q(u,v) == pi_star(v)*Q(v,u) across ALL 61x61 pairs
    CHECK(std::abs(pi_star.sum() - 1.0) < 1e-13);
    Scalar max_db_err = 0.0;
    for (int u = 0; u < S; ++u) {
        CHECK(std::abs(Q.row(u).sum()) < 1e-13);
        CHECK(std::abs(P_short.row(u).sum() - 1.0) < 1e-13);
        for (int v = 0; v < S; ++v) {
            Scalar f_uv = pi_star(u) * Q(u, v);
            Scalar f_vu = pi_star(v) * Q(v, u);
            max_db_err = std::max(max_db_err, std::abs(f_uv - f_vu));
            CHECK(P_short(u, v) >= 0.0);
        }
    }
    CHECK(max_db_err < 1e-13);

    // 2. Stationary invariance pi_star^T * P(Delta_t) == pi_star^T
    Vector pushed = P_short.transpose() * pi_star;
    CHECK((pushed - pi_star).cwiseAbs().maxCoeff() < 1e-12);

    // 3. Long-time ergodic convergence P(t -> inf) rows -> pi_star^T
    Matrix Q2;
    Vector pi_star2;
    ForecastAnalyzer::build_reversible_site_codon_model(
        S, *code, gtr, pos_nuc_freqs, codon_freqs, syn_mat, prop_diffs,
        1.1, 0.35, site_lambdas, phi_aa, 8.0, 2000.0, Q2, pi_star2, P_long
    );
    for (int u = 0; u < S; ++u) {
        CHECK((P_long.row(u).transpose() - pi_star).cwiseAbs().maxCoeff() < 1e-10);
    }
}

TEST_CASE("FORECAST: End-to-End DMS & Codon Trajectory Forecasting on cd2.fna") {
    std::string aln_path = resolve_test_path("benchmarks/data/cd2.fna");
    std::string tree_path = resolve_test_path("benchmarks/data/cd2.nwk");

    Alignment aln = Alignment::from_fasta(aln_path);
    Tree tree = Tree::from_newick_file(tree_path);

    ForecastSettings settings;
    settings.branches = "Internal";
    settings.property_set = "3PROP";
    settings.shrinkage_prior = "global";
    settings.shrinkage_lambda = 0.5;
    settings.syn_model = "SynREV";
    settings.ancestral_exponent = 0.5;
    settings.horizon = 0.05;
    settings.reference_taxon = aln.taxon_names.front();

    ForecastAnalyzer analyzer(tree, aln, settings);
    ForecastResult res = analyzer.run(false, false);

    const size_t L = aln.num_codons;
    const int S = static_cast<int>(res.num_sense_codons);
    REQUIRE(res.total_sites == L);
    REQUIRE(res.aa_dms_matrix.rows() == static_cast<Eigen::Index>(L));
    REQUIRE(res.aa_dms_matrix.cols() == 20);
    REQUIRE(res.aa_stationary_matrix.rows() == static_cast<Eigen::Index>(L));
    REQUIRE(res.aa_stationary_matrix.cols() == 20);
    REQUIRE(res.ancestral_aa_matrix.rows() == static_cast<Eigen::Index>(L));
    REQUIRE(res.ancestral_aa_matrix.cols() == 20);
    REQUIRE(res.codon_rate_matrix.rows() == static_cast<Eigen::Index>(L));
    REQUIRE(res.codon_rate_matrix.cols() == S);
    REQUIRE(res.codon_log_prob_matrix.rows() == static_cast<Eigen::Index>(L));
    REQUIRE(res.codon_log_prob_matrix.cols() == S);

    CHECK(res.effective_branches == "Internal");
    CHECK(res.effective_syn_model == "SynREV");
    CHECK(res.global_omega > 0.0);
    CHECK(res.global_lambdas.size() == 3);

    size_t invariant_sites_with_nonconstant_dms = 0;
    size_t sites_with_multistep_prob = 0;

    for (size_t s = 0; s < L; ++s) {
        const auto& sr = res.site_results[s];
        int focal_a = aa_char_to_index(sr.focal_aa);
        REQUIRE(focal_a >= 0);

        // Focal AA DMS shift is identically 0.0
        CHECK(std::abs(res.aa_dms_matrix(s, focal_a)) < 1e-12);

        // Stationary and ancestral AA distributions sum to 1.0
        CHECK(std::abs(res.aa_stationary_matrix.row(s).sum() - 1.0) < 1e-10);
        CHECK(std::abs(res.ancestral_aa_matrix.row(s).sum() - 1.0) < 1e-10);

        // Horizon codon transition probabilities sum to 1.0
        Scalar prob_sum = 0.0;
        for (int c = 0; c < S; ++c) {
            CHECK(std::isfinite(res.codon_log_prob_matrix(s, c)));
            prob_sum += std::exp(res.codon_log_prob_matrix(s, c));
        }
        CHECK(std::abs(prob_sum - 1.0) < 1e-9);

        // Check that invariant sites still receive non-constant within-site 20-AA DMS rankings
        // via hierarchical G-PRIME Empirical Bayes prior inheritance!
        if (sr.beta <= 1e-4) {
            Scalar min_dms = 1e30, max_dms = -1e30;
            for (int a = 0; a < 20; ++a) {
                if (a == focal_a) continue;
                min_dms = std::min(min_dms, res.aa_dms_matrix(s, a));
                max_dms = std::max(max_dms, res.aa_dms_matrix(s, a));
            }
            if (max_dms - min_dms > 1e-3) {
                invariant_sites_with_nonconstant_dms++;
            }
        }

        // Check 2-nt / 3-nt multi-step codons have strictly positive finite-horizon probability
        for (int c = 0; c < S; ++c) {
            if (aln.code->diff_matrix[sr.focal_codon_index][c].num_diffs >= 2) {
                if (res.codon_log_prob_matrix(s, c) > -50.0) {
                    sites_with_multistep_prob++;
                    break;
                }
            }
        }
    }

    CHECK(invariant_sites_with_nonconstant_dms > 0);
    CHECK(sites_with_multistep_prob == L);

    // JSON v3 & TSV serialization check
    nlohmann::json mj = res.to_modern_json(analyzer.tree, aln, Provenance{}, aln_path);
    CHECK(mj["schema_version"] == "3.0.0");
    CHECK(mj["analysis"]["id"] == "forecast");
    CHECK(mj["site_results"]["row_count"] == static_cast<int>(L));
    CHECK(mj["method_results"]["aa_dms_matrix"].size() == L);
    CHECK(mj["method_results"]["codon_log_prob_matrix"].size() == L);

    std::string tmp_dms = "/tmp/hyphy3_test_forecast_dms.tsv";
    std::string tmp_cod = "/tmp/hyphy3_test_forecast_codon.tsv";
    res.write_dms_tsv(tmp_dms);
    res.write_codon_tsv(tmp_cod);
    CHECK(std::ifstream(tmp_dms).good());
    CHECK(std::ifstream(tmp_cod).good());
    std::remove(tmp_dms.c_str());
    std::remove(tmp_cod.c_str());
}
