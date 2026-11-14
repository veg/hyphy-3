#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/analyses/prime.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace hyphy::core;
using namespace hyphy::analyses;

inline std::string resolve_test_path(const std::string& path) {
    if (std::ifstream(path).good()) return path;
    if (std::ifstream("../" + path).good()) return "../" + path;
    if (std::ifstream("../../" + path).good()) return "../../" + path;
    return path;
}

TEST_CASE("PRIME: Biophysical Property Standardization & Symmetry") {
    auto props5 = PropertyDatabase::get_5prop_properties();
    CHECK(props5.size() == 5);

    // Verify properties are standardized (mean=0, std=1)
    for (const auto& p : props5) {
        Scalar sum = 0.0;
        Scalar sum_sq = 0.0;
        for (int i = 0; i < 20; ++i) {
            sum += p.z_scores[i];
            sum_sq += p.z_scores[i] * p.z_scores[i];
        }
        Scalar mean_z = sum / 20.0;
        Scalar std_z = std::sqrt(sum_sq / 20.0 - mean_z * mean_z);

        CHECK(std::abs(mean_z) < 1e-12);
        CHECK(std::abs(std_z - 1.0) < 1e-12);
    }

    // Verify Atchley 5-factor representation
    auto atchley = PropertyDatabase::get_atchley_properties();
    CHECK(atchley.size() == 5);
    for (const auto& p : atchley) {
        Scalar sum = 0.0;
        Scalar sum_sq = 0.0;
        for (int i = 0; i < 20; ++i) {
            sum += p.z_scores[i];
            sum_sq += p.z_scores[i] * p.z_scores[i];
        }
        Scalar mean_z = sum / 20.0;
        Scalar std_z = std::sqrt(sum_sq / 20.0 - mean_z * mean_z);

        CHECK(std::abs(mean_z) < 1e-12);
        CHECK(std::abs(std_z - 1.0) < 1e-12);
    }
}

TEST_CASE("PRIME: PRIMEMatrix Eigendecomposition, Parity with MG94, & Detailed Balance") {
    auto gcode = GeneticCode::universal();
    int S = gcode->num_sense_codons;

    Eigen::Matrix<Scalar, 3, 4> pos_nuc_freqs;
    pos_nuc_freqs.fill(0.25);
    Vector codon_freqs = Vector::Constant(S, 1.0 / static_cast<Scalar>(S));

    GTRParameters gtr;
    gtr.theta_AC = 0.8;
    gtr.theta_AT = 0.4;
    gtr.theta_CG = 0.6;
    gtr.theta_CT = 1.2;
    gtr.theta_GT = 0.5;

    auto props = PropertyDatabase::get_5prop_properties();
    size_t K = props.size();

    std::vector<Matrix> codon_prop_diffs(K, Matrix::Zero(S, S));
    for (size_t k = 0; k < K; ++k) {
        for (int u = 0; u < S; ++u) {
            int idx_u = aa_char_to_index(gcode->sense_codon_aa[u]);
            for (int v = 0; v < S; ++v) {
                int idx_v = aa_char_to_index(gcode->sense_codon_aa[v]);
                codon_prop_diffs[k](u, v) = std::abs(props[k].z_scores[idx_u] - props[k].z_scores[idx_v]);
            }
        }
    }

    // 1. When all lambdas = 0, PRIMEMatrix must match MG94Matrix exactly
    PRIMEMatrix pmat_null;
    std::vector<Scalar> zero_lambdas(K, 0.0);
    pmat_null.update(1.2, 0.7, zero_lambdas, codon_prop_diffs, gtr, pos_nuc_freqs, codon_freqs, *gcode);

    MG94Parameters mg_params;
    mg_params.alpha = 1.2;
    mg_params.beta = 0.7;
    mg_params.theta_AC = gtr.theta_AC;
    mg_params.theta_AT = gtr.theta_AT;
    mg_params.theta_CG = gtr.theta_CG;
    mg_params.theta_CT = gtr.theta_CT;
    mg_params.theta_GT = gtr.theta_GT;

    MG94Matrix mg_mat;
    mg_mat.update(mg_params, pos_nuc_freqs, codon_freqs, *gcode);

    Scalar max_q_diff = (pmat_null.Q - mg_mat.Q).cwiseAbs().maxCoeff();
    CHECK(max_q_diff < 1e-12);

    Scalar max_eig_diff = (pmat_null.eigenvalues - mg_mat.eigenvalues).cwiseAbs().maxCoeff();
    CHECK(max_eig_diff < 1e-12);

    Matrix P_prime = pmat_null.transition_matrix(0.15);
    Matrix P_mg = mg_mat.transition_matrix(0.15);
    Scalar max_p_diff = (P_prime - P_mg).cwiseAbs().maxCoeff();
    CHECK(max_p_diff < 1e-12);

    // 2. Test Detailed Balance under non-zero property selection
    std::vector<Scalar> test_lambdas = {1.5, -0.8, 2.0, -1.2, 0.5};
    PRIMEMatrix pmat_alt;
    pmat_alt.update(1.0, 0.9, test_lambdas, codon_prop_diffs, gtr, pos_nuc_freqs, codon_freqs, *gcode);

    Scalar max_flux_asym = 0.0;
    for (int i = 0; i < S; ++i) {
        for (int j = 0; j < S; ++j) {
            Scalar flux_ij = codon_freqs(i) * pmat_alt.Q(i, j);
            Scalar flux_ji = codon_freqs(j) * pmat_alt.Q(j, i);
            max_flux_asym = std::max(max_flux_asym, std::abs(flux_ij - flux_ji));
        }
    }
    CHECK(max_flux_asym < 1e-12);

    // 3. Test Numerical Stability under Extreme Parameters
    std::vector<Scalar> extreme_lambdas = {100.0, -100.0, 50.0, -50.0, 0.0};
    PRIMEMatrix pmat_extreme;
    pmat_extreme.update(0.01, 10.0, extreme_lambdas, codon_prop_diffs, gtr, pos_nuc_freqs, codon_freqs, *gcode);

    for (int i = 0; i < S; ++i) {
        CHECK(!std::isnan(pmat_extreme.eigenvalues(i)));
        CHECK(!std::isinf(pmat_extreme.eigenvalues(i)));
        CHECK(pmat_extreme.eigenvalues(i) <= 1e-12);
    }
    Matrix P_extreme = pmat_extreme.transition_matrix(0.5);
    for (int i = 0; i < S; ++i) {
        Scalar row_sum = P_extreme.row(i).sum();
        CHECK(std::abs(row_sum - 1.0) < 1e-6);
        for (int j = 0; j < S; ++j) {
            CHECK(P_extreme(i, j) >= -1e-12);
            CHECK(P_extreme(i, j) <= 1.0 + 1e-6);
        }
    }
}

TEST_CASE("PRIME: Full S-PRIME Inference on cd2.fna with Guaranteed Nesting Monotonicity") {
    std::string aln_path = resolve_test_path("benchmarks/data/cd2.fna");
    std::string tree_path = resolve_test_path("benchmarks/data/cd2.nwk");

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick_file(tree_path);

    PRIMESettings settings;
    settings.property_set = "5PROP";
    settings.shrinkage_lambda = 0.5;
    settings.p_value_threshold = 0.10;
    settings.stage2_p_value = 0.10;
    settings.fit_global = true;

    auto analyzer = PRIMEAnalyzer::create_and_fit(tree, aln, settings);

    // Verify global G-PRIME results
    CHECK(analyzer.global_fit.lambdas.size() == 5);
    CHECK(analyzer.global_fit.log_likelihood < 0.0);
    std::cout << "G-PRIME baseline log-likelihood: " << analyzer.global_fit.log_likelihood << "\n";
    std::cout << "G-PRIME global beta: " << analyzer.global_fit.beta << "\n";
    std::cout << "G-PRIME lambdas: [";
    for (size_t k = 0; k < 5; ++k) {
        std::cout << (k > 0 ? ", " : "") << analyzer.global_fit.lambdas[k];
    }
    std::cout << "]\n";

    // Run S-PRIME
    auto results = analyzer.run(nullptr, false, false);
    CHECK(results.site_results.size() == aln.num_codons);
    CHECK(results.total_sites == aln.num_codons);

    // Verify Guaranteed Monotonicity & Nesting Axioms across all sites:
    // Axiom 1: LogL(Alt) >= LogL(FEL) by construction
    // Axiom 2: LogL(Alt) >= LogL(lambda_k = 0) by construction
    // Axiom 3: Omnibus LRT >= 0
    // Axiom 4: Property LRTs >= 0
    size_t invariant_count = 0;
    size_t variable_count = 0;

    for (const auto& sr : results.site_results) {
        if (sr.alpha == 0.0 && sr.beta == 0.0 && sr.omnibus_lrt == 0.0) {
            invariant_count++;
            CHECK(sr.omnibus_p_value == 1.0);
            CHECK(sr.omnibus_q_value == 1.0);
            continue;
        }
        variable_count++;

        // Axioms
        CHECK(sr.log_l_alt >= sr.log_l_fel - 1e-6);
        CHECK(sr.omnibus_lrt >= 0.0);
        CHECK(sr.omnibus_p_value >= 0.0);
        CHECK(sr.omnibus_p_value <= 1.0);
        CHECK(sr.omnibus_q_value >= 0.0);
        CHECK(sr.omnibus_q_value <= 1.0);

        for (const auto& p_inf : sr.properties) {
            CHECK(sr.log_l_alt >= p_inf.log_l_null - 1e-6);
            CHECK(p_inf.lrt >= 0.0);
            CHECK(p_inf.p_value >= 0.0);
            CHECK(p_inf.p_value <= 1.0);
            CHECK(p_inf.adjusted_p_value >= 0.0);
            CHECK(p_inf.adjusted_p_value <= 1.0);
        }
    }

    CHECK(invariant_count > 0);
    CHECK(variable_count > 0);
    std::cout << "Tested " << results.total_sites << " sites: "
              << invariant_count << " invariant, "
              << variable_count << " variable.\n";
    std::cout << "Significant sites under omnibus test: " << results.significant_sites << "\n";

    // Verify modern JSON serialization
    nlohmann::json modern_json = results.to_modern_json(tree, aln, "hyphy3 prime", aln_path);
    CHECK(modern_json["schema_version"] == "3.0.0");
    CHECK(modern_json["analysis"]["id"] == "prime");
    CHECK(modern_json["site_results"]["data"]["site"].size() == aln.num_codons);
    CHECK(modern_json["site_results"]["data"]["omnibus_q_value"].size() == aln.num_codons);

    // Verify legacy JSON serialization
    nlohmann::json legacy_json = results.to_legacy_json(tree, aln, aln_path);
    CHECK(legacy_json["analysis"]["version"] == "0.2.6");
    CHECK(legacy_json["MLE"]["content"]["0"].size() == aln.num_codons);
    CHECK(legacy_json["MLE"]["headers"].size() == 12 + 3 * 5); // 12 standard + 3*5 property columns = 27
}

TEST_CASE("PRIME site-level progress callbacks") {
    std::string aln_path = resolve_test_path("benchmarks/data/cd2.fna");
    std::string tree_path = resolve_test_path("benchmarks/data/cd2.nwk");
    Alignment aln = Alignment::from_fasta(aln_path);
    Tree tree = Tree::from_newick_file(tree_path);

    PRIMESettings settings;
    settings.property_set = "5PROP";
    settings.fit_global = false;

    auto analyzer = PRIMEAnalyzer::create_and_fit(tree, aln, settings);

    std::atomic<size_t> live_cb_count{0};
    size_t final_cb_count = 0;

    auto live_cb = [&](size_t site_idx, size_t cur_done, size_t total, const PRIMESiteResult& res) {
        live_cb_count.fetch_add(1, std::memory_order_relaxed);
        CHECK(site_idx < aln.num_codons);
        CHECK(cur_done <= total);
        CHECK(total == aln.num_codons);
        CHECK(res.omnibus_p_value >= 0.0);
    };

    auto final_cb = [&](size_t site_idx, const PRIMESiteResult& res) {
        final_cb_count++;
        CHECK(site_idx < aln.num_codons);
        CHECK(res.omnibus_q_value >= 0.0);
    };

    auto results = analyzer.run(nullptr, false, false, final_cb, live_cb);
    CHECK(live_cb_count.load() == aln.num_codons);
    CHECK(final_cb_count == aln.num_codons);
}

TEST_CASE("PRIME: LCAP Canonical A..Y Ordering & G-PRIME Global Empirical Bayes Prior") {
    auto lcap = PropertyDatabase::get_predefined_properties("LCAP");
    REQUIRE(lcap.size() == 5);
    for (const auto& p : lcap) {
        Scalar sum = 0.0, sum_sq = 0.0;
        for (int i = 0; i < 20; ++i) {
            sum += p.z_scores[i];
            sum_sq += p.z_scores[i] * p.z_scores[i];
        }
        CHECK(std::abs(sum / 20.0) < 1e-12);
        CHECK(std::abs(std::sqrt(sum_sq / 20.0) - 1.0) < 1e-12);
    }

    // Verify canonical IUPAC indexing (preventing the legacy FLIMVSPTAYHQNKDECWRG scramble):
    // 0: Chemical Composition -> C (idx 1) is maximum
    CHECK(lcap[0].z_scores[aa_char_to_index('C')] > 3.0);
    // 2: Molecular Volume -> W (idx 18) is maximum, G (idx 5) is minimum
    CHECK(lcap[2].z_scores[aa_char_to_index('W')] > 2.0);
    CHECK(lcap[2].z_scores[aa_char_to_index('G')] < -1.5);
    // 3: Iso-electric Point -> R (idx 14) is maximum, D (idx 2) is minimum
    CHECK(lcap[3].z_scores[aa_char_to_index('R')] > 2.4);
    CHECK(lcap[3].z_scores[aa_char_to_index('D')] < -1.7);

    // Test --prior global on Internal branches
    std::string aln_path = resolve_test_path("benchmarks/data/cd2.fna");
    std::string tree_path = resolve_test_path("benchmarks/data/cd2.nwk");
    Alignment aln = Alignment::from_fasta(aln_path);
    Tree tree = Tree::from_newick_file(tree_path);

    PRIMESettings settings;
    settings.property_set = "3PROP";
    settings.shrinkage_lambda = 0.5;
    settings.shrinkage_prior = "global";
    settings.test_branches = "Internal";

    auto analyzer = PRIMEAnalyzer::create_and_fit(tree, aln, settings);
    REQUIRE(analyzer.global_fit.lambdas.size() == 3);
    auto results = analyzer.run(nullptr, false, false);

    size_t inv_checked = 0;
    for (const auto& sr : results.site_results) {
        CHECK(sr.omnibus_lrt >= 0.0);
        if (sr.fel_beta <= 1e-4) {
            inv_checked++;
            CHECK(sr.omnibus_lrt == 0.0);
            CHECK(sr.omnibus_p_value == 1.0);
            // Invariant / zero-nonsynonymous sites inherit G-PRIME global lambdas under --prior global
            for (size_t k = 0; k < 3; ++k) {
                CHECK(std::abs(sr.properties[k].lambda - results.global_fit.lambdas[k]) < 1e-12);
            }
        }
    }
    CHECK(inv_checked > 0);

    nlohmann::json mj = results.to_modern_json(tree, aln, "hyphy3 prime --prior global", aln_path);
    CHECK(mj["statistical_tests"]["site_level_summary"]["shrinkage_prior"] == "global");
}

