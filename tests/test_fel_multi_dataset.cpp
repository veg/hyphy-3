#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/core/likelihood.hpp"
#include "hyphy/analyses/fel.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>

using namespace hyphy::core;
using namespace hyphy::analyses;

inline std::vector<double> compute_ranks(const std::vector<double>& v) {
    size_t n = v.size();
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&v](size_t a, size_t b) {
        return v[a] < v[b];
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

TEST_CASE("Parity 1: CD2 (Universal Code, 10 taxa, 187 codons)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.fna");
    auto tree = Tree::from_newick_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.nwk");
    auto fel = FELAnalyzer::create_and_fit(tree, aln, 0.1);
    auto results = fel.run();

    std::ifstream gt_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/ground_truth/cd2_fel.json");
    REQUIRE(gt_file.is_open());
    nlohmann::json gt;
    gt_file >> gt;
    const auto& gt_mle = gt["MLE"]["content"]["0"];

    std::vector<double> p_lrt(results.size()), g_lrt(results.size());
    std::vector<double> p_pv(results.size()), g_pv(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        p_lrt[i] = results[i].lrt;
        g_lrt[i] = gt_mle[i][3].get<double>();
        p_pv[i] = results[i].p_value;
        g_pv[i] = gt_mle[i][4].get<double>();
    }

    double r_lrt = pearson_correlation(p_lrt, g_lrt);
    double rho_lrt = spearman_correlation(p_lrt, g_lrt);
    std::cout << "[Parity CD2] LRT r = " << r_lrt << ", rho = " << rho_lrt << std::endl;
    CHECK(r_lrt >= 0.998);
    CHECK(rho_lrt >= 0.995);
}

TEST_CASE("Parity 2: ADH (Universal Code, 23 taxa, 254 codons)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/tests/data/adh.nex");
    REQUIRE(!aln.embedded_tree_newick.empty());
    auto tree = Tree::from_newick(aln.embedded_tree_newick);
    auto fel = FELAnalyzer::create_and_fit(tree, aln, 0.1);
    auto results = fel.run();

    std::ifstream gt_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/ground_truth/adh_fel_rel.json");
    REQUIRE(gt_file.is_open());
    nlohmann::json gt;
    gt_file >> gt;
    const auto& gt_mle = gt["MLE"]["content"]["0"];

    std::vector<double> p_lrt(results.size()), g_lrt(results.size());
    std::vector<double> p_pv(results.size()), g_pv(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        p_lrt[i] = results[i].lrt;
        g_lrt[i] = gt_mle[i][3].get<double>();
        p_pv[i] = results[i].p_value;
        g_pv[i] = gt_mle[i][4].get<double>();
    }

    double r_lrt = pearson_correlation(p_lrt, g_lrt);
    double rho_lrt = spearman_correlation(p_lrt, g_lrt);
    std::cout << "[Parity ADH] LRT r = " << r_lrt << ", rho = " << rho_lrt << std::endl;
    CHECK(r_lrt >= 0.998);
    CHECK(rho_lrt >= 0.998);
}

TEST_CASE("Parity 3: BGLOBIN (Universal Code, 17 taxa, 144 codons)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/tests/data/bglobin.nex");
    REQUIRE(!aln.embedded_tree_newick.empty());
    auto tree = Tree::from_newick(aln.embedded_tree_newick);
    auto fel = FELAnalyzer::create_and_fit(tree, aln, 0.1);
    auto results = fel.run();

    std::ifstream gt_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/ground_truth/bglobin_fel_rel.json");
    REQUIRE(gt_file.is_open());
    nlohmann::json gt;
    gt_file >> gt;
    const auto& gt_mle = gt["MLE"]["content"]["0"];

    std::vector<double> p_lrt(results.size()), g_lrt(results.size());
    std::vector<double> p_pv(results.size()), g_pv(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        p_lrt[i] = results[i].lrt;
        g_lrt[i] = gt_mle[i][3].get<double>();
        p_pv[i] = results[i].p_value;
        g_pv[i] = gt_mle[i][4].get<double>();
    }

    double r_lrt = pearson_correlation(p_lrt, g_lrt);
    double rho_lrt = spearman_correlation(p_lrt, g_lrt);
    std::cout << "[Parity BGLOBIN] LRT r = " << r_lrt << ", rho = " << rho_lrt << std::endl;
    CHECK(r_lrt >= 0.998);
    CHECK(rho_lrt >= 0.998);
}

TEST_CASE("Parity 4: COXI (Vertebrate-mtDNA Non-Universal Code, 21 taxa, 510 codons)") {
    auto gcode = GeneticCode::vertebrate_mtdna();
    REQUIRE(gcode->num_sense_codons == 60);

    auto aln = Alignment::load("/Users/sergei/Development/hyphy/tests/hbltests/libv3/data/COXI.nex", gcode);
    REQUIRE(aln.num_codons == 510);
    REQUIRE(!aln.embedded_tree_newick.empty());

    auto tree = Tree::from_newick(aln.embedded_tree_newick);
    auto fel = FELAnalyzer::create_and_fit(tree, aln, 0.1);
    auto results = fel.run();

    std::ifstream gt_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/ground_truth/coxi_fel_rel.json");
    REQUIRE(gt_file.is_open());
    nlohmann::json gt;
    gt_file >> gt;
    const auto& gt_mle = gt["MLE"]["content"]["0"];

    std::vector<double> p_lrt(results.size()), g_lrt(results.size());
    std::vector<double> p_pv(results.size()), g_pv(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        p_lrt[i] = results[i].lrt;
        g_lrt[i] = gt_mle[i][3].get<double>();
        p_pv[i] = results[i].p_value;
        g_pv[i] = gt_mle[i][4].get<double>();
    }

    double r_lrt = pearson_correlation(p_lrt, g_lrt);
    double rho_lrt = spearman_correlation(p_lrt, g_lrt);
    double r_pv = pearson_correlation(p_pv, g_pv);
    double rho_pv = spearman_correlation(p_pv, g_pv);

    std::cout << "[Parity COXI (Vertebrate-mtDNA)]\n"
              << "  LRT Pearson r   : " << r_lrt << "\n"
              << "  LRT Spearman rho: " << rho_lrt << "\n"
              << "  p-val Pearson r : " << r_pv << "\n"
              << "  p-val Spearman  : " << rho_pv << "\n";

    CHECK(r_lrt >= 0.995);
    CHECK(rho_lrt >= 0.995);
    CHECK(r_pv >= 0.995);
    CHECK(rho_pv >= 0.995);

    double gt_gtr_logl = gt["fits"]["Nucleotide GTR"]["Log Likelihood"].get<double>();
    CHECK(doctest::Approx(fel.gtr_log_l).epsilon(1e-4) == gt_gtr_logl);

    double gt_w = gt["fits"]["Global MG94xREV"]["Rate Distributions"]["non-synonymous/synonymous rate ratio for *test*"][0][0].get<double>();
    double pred_w = fel.base_params.beta / fel.base_params.alpha;
    CHECK(std::abs(pred_w - gt_w) < 0.005);
}
