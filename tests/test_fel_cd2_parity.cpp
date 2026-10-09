#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/analyses/fel.hpp"
#include "hyphy/opt/optimizer.hpp"
#include "nlohmann/json.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

using namespace hyphy::core;
using namespace hyphy::analyses;
using namespace hyphy::opt;

inline std::vector<double> compute_ranks(const std::vector<double>& v) {
    size_t n = v.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](size_t i1, size_t i2) {
        return v[i1] < v[i2];
    });

    std::vector<double> ranks(n);
    size_t i = 0;
    while (i < n) {
        size_t j = i + 1;
        while (j < n && std::abs(v[idx[j]] - v[idx[i]]) < 1e-12) {
            j++;
        }
        double rank = 0.5 * (i + j - 1) + 1.0;
        for (size_t k = i; k < j; ++k) {
            ranks[idx[k]] = rank;
        }
        i = j;
    }
    return ranks;
}

inline double pearson_correlation(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();
    if (n == 0 || n != y.size()) return 0.0;
    double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double num = 0.0, den_x = 0.0, den_y = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - mx;
        double dy = y[i] - my;
        num += dx * dy;
        den_x += dx * dx;
        den_y += dy * dy;
    }
    return (den_x > 0.0 && den_y > 0.0) ? (num / std::sqrt(den_x * den_y)) : 0.0;
}

inline double spearman_correlation(const std::vector<double>& x, const std::vector<double>& y) {
    return pearson_correlation(compute_ranks(x), compute_ranks(y));
}

TEST_CASE("FEL Parity on Full CD2 Dataset: Matched Branch Lengths") {
    std::string aln_path = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.fna";
    std::string tree_path = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2_fitted.nwk";
    std::string gt_json_path = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/ground_truth/cd2_fel.json";
    
    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick_file(tree_path);

    MG94Parameters base_p;
    base_p.theta_AC = 0.5521777318173139;
    base_p.theta_AT = 0.2653120245806623;
    base_p.theta_CG = 0.4932059218537349;
    base_p.theta_CT = 1.03133955537726;
    base_p.theta_GT = 0.3054313963126315;
    base_p.alpha = 1.0;
    base_p.beta = 1.012387352492807;

    FELAnalyzer fel(tree, aln, base_p);
    auto results = fel.run();
    CHECK(results.size() == 187);

    // Read HyPhy 2.5 ground truth
    std::ifstream gt_file(gt_json_path);
    REQUIRE(gt_file.is_open());
    nlohmann::json gt_j;
    gt_file >> gt_j;
    const auto& gt_mle = gt_j["MLE"]["content"]["0"];

    std::vector<double> pred_lrt(results.size()), gt_lrt(results.size());
    std::vector<double> pred_pval(results.size()), gt_pval(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        pred_lrt[i] = results[i].lrt;
        gt_lrt[i] = gt_mle[i][3].get<double>();
        pred_pval[i] = results[i].p_value;
        gt_pval[i] = gt_mle[i][4].get<double>();
    }

    double r_lrt = pearson_correlation(pred_lrt, gt_lrt);
    double rho_lrt = spearman_correlation(pred_lrt, gt_lrt);
    double r_pval = pearson_correlation(pred_pval, gt_pval);
    double rho_pval = spearman_correlation(pred_pval, gt_pval);

    std::cout << "[CD2 Matched Tree Parity]\n"
              << "  LRT Pearson r   : " << r_lrt << "\n"
              << "  LRT Spearman rho: " << rho_lrt << "\n"
              << "  p-val Pearson r : " << r_pval << "\n"
              << "  p-val Spearman  : " << rho_pval << "\n";

    CHECK(r_lrt > 0.998);
    CHECK(rho_lrt > 0.995);
    CHECK(r_pval > 0.993);
    CHECK(rho_pval > 0.995);
}

TEST_CASE("FEL Parity on Full CD2 Dataset: End-to-End from Unannotated Tree") {
    std::string aln_path = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.fna";
    std::string tree_path = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.nwk";
    std::string gt_json_path = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/ground_truth/cd2_fel.json";
    
    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick_file(tree_path);

    // Phase 1: Nucleotide GTR fitting
    GTRFitter gtr_fitter(tree, aln);
    auto gtr_res = gtr_fitter.fit();

    // Verify GTR parameter parity against HyPhy 2.5
    CHECK(doctest::Approx(gtr_res.params.theta_AC).epsilon(0.01) == 0.552);
    CHECK(doctest::Approx(gtr_res.params.theta_AT).epsilon(0.01) == 0.265);
    CHECK(doctest::Approx(gtr_res.params.theta_CG).epsilon(0.01) == 0.493);
    CHECK(doctest::Approx(gtr_res.params.theta_CT).epsilon(0.01) == 1.031);
    CHECK(doctest::Approx(gtr_res.params.theta_GT).epsilon(0.01) == 0.305);
    CHECK(gtr_res.log_likelihood > -3540.0);

    // Phase 2: Global MG94 refinement
    MG94Parameters base_p;
    base_p.theta_AC = gtr_res.params.theta_AC;
    base_p.theta_AT = gtr_res.params.theta_AT;
    base_p.theta_CG = gtr_res.params.theta_CG;
    base_p.theta_CT = gtr_res.params.theta_CT;
    base_p.theta_GT = gtr_res.params.theta_GT;

    MG94Fitter mg_fitter(gtr_res.tree, aln);
    auto mg_res = mg_fitter.fit_omega_and_scale(1.0, base_p);
    base_p.alpha = 1.0;
    base_p.beta = mg_res.x_opt(0);

    CHECK(doctest::Approx(mg_res.x_opt(0)).epsilon(0.05) == 1.012); // Global omega
    CHECK(mg_res.log_likelihood > -3475.0);

    // Phase 3: Site testing
    FELAnalyzer fel(mg_res.tree, aln, base_p);
    auto results = fel.run();
    CHECK(results.size() == 187);

    // Read HyPhy 2.5 ground truth
    std::ifstream gt_file(gt_json_path);
    REQUIRE(gt_file.is_open());
    nlohmann::json gt_j;
    gt_file >> gt_j;
    const auto& gt_mle = gt_j["MLE"]["content"]["0"];

    std::vector<double> pred_lrt(results.size()), gt_lrt(results.size());
    std::vector<double> pred_pval(results.size()), gt_pval(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        pred_lrt[i] = results[i].lrt;
        gt_lrt[i] = gt_mle[i][3].get<double>();
        pred_pval[i] = results[i].p_value;
        gt_pval[i] = gt_mle[i][4].get<double>();
    }

    double r_lrt = pearson_correlation(pred_lrt, gt_lrt);
    double rho_lrt = spearman_correlation(pred_lrt, gt_lrt);
    double r_pval = pearson_correlation(pred_pval, gt_pval);
    double rho_pval = spearman_correlation(pred_pval, gt_pval);

    std::cout << "[CD2 Unannotated Tree End-to-End Parity]\n"
              << "  LRT Pearson r   : " << r_lrt << "\n"
              << "  LRT Spearman rho: " << rho_lrt << "\n"
              << "  p-val Pearson r : " << r_pval << "\n"
              << "  p-val Spearman  : " << rho_pval << "\n";

    CHECK(r_lrt > 0.998);
    CHECK(rho_lrt > 0.995);
    CHECK(r_pval > 0.993);
    CHECK(rho_pval > 0.995);
}
