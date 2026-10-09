#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/core/types.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/rate_matrix.hpp"
#include "hyphy/analyses/fel.hpp"

using namespace hyphy::core;
using namespace hyphy::analyses;

TEST_CASE("FEL: Site-by-site dN/dS estimation on cd2_reduced with matched branch lengths") {
    std::string aln_path = "/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2_reduced.fna";
    auto aln = Alignment::from_fasta(aln_path);

    // Tree with exact fitted branch lengths from HyPhy 2.5 cd2_reduced_fel.json (Global MG94xREV)
    std::string nwk_fitted = "(Rat:0.0,Mouse:0.09748291588565174,((Human:0.0,Chimp:0.05545563815378509)Node12:0.1321590496511669,RhMonkey:0.0,Baboon:0.0)Node8:0.8462121433135308,Horse:0.95992364787302,Cat:1.197250478870998,Pig:0.4755389109237619,Cow:0.5939647347629947);";
    auto tree = Tree::from_newick(nwk_fitted);

    MG94Parameters base_p;
    // GTR substitution rates fitted by HyPhy 2.5 Nucleotide GTR phase:
    // AC: 0, AG: 1, AT: 0, CG: 1.939655, CT: 1.751827, GT: 0
    base_p.theta_AC = 0.0;
    base_p.theta_AT = 0.0;
    base_p.theta_CG = 1.939655;
    base_p.theta_CT = 1.751827;
    base_p.theta_GT = 0.0;
    base_p.alpha = 1.0;
    base_p.beta = 0.5817295780925678;

    FELAnalyzer fel(tree, aln, base_p);
    fel.global_log_l = -60.498956;
    fel.global_aicc = 222.2479;
    auto results = fel.run();

    std::string out_json = "/tmp/hyphy_next_fel_test.json";
    fel.save_json(out_json, aln_path, nwk_fitted);
    std::cout << "Saved FEL JSON to: " << out_json << "\n";

    std::cout << "\n=== FEL Site-by-Site Results (hyphy-next vs HyPhy 2.5) ===\n";
    // Ground truth from cd2_reduced_fel.json:
    // Site 0: alpha=0.0079, beta=8.3413, alpha=beta=2.9705, LRT=1.2225, pval=0.2689
    // Site 1: alpha=0.0,    beta=0.0,    alpha=beta=0.0,    LRT=0.0,    pval=1.0
    // Site 2: alpha=0.0,    beta=0.0,    alpha=beta=0.0,    LRT=0.0,    pval=1.0
    // Site 3: alpha=0.6517, beta=0.3637, alpha=beta=0.4357, LRT=0.1608, pval=0.6884
    // Site 4: alpha=1.4011, beta=0.6081, alpha=beta=0.7910, LRT=0.6365, pval=0.4250
    // Site 5: alpha=1.3208, beta=0.6669, alpha=beta=0.8145, LRT=0.3191, pval=0.5722

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        std::cout << "Site " << i 
                  << ": alpha=" << r.alpha 
                  << ", beta=" << r.beta 
                  << ", alpha=beta=" << r.alpha_null
                  << ", LRT=" << r.lrt 
                  << ", p-val=" << r.p_value 
                  << ", logL_alt=" << r.log_l_alt 
                  << ", logL_null=" << r.log_l_null << "\n";
    }

    // Invariant sites must be exactly 0
    CHECK(results[1].alpha == 0.0);
    CHECK(results[1].beta == 0.0);
    CHECK(results[1].p_value == 1.0);

    CHECK(results[2].alpha == 0.0);
    CHECK(results[2].beta == 0.0);
    CHECK(results[2].p_value == 1.0);

    // Variable sites should detect consistent rates and p-values
    CHECK(results[0].beta > results[0].alpha); // Pos 0 has non-syn change
    CHECK(results[3].p_value > 0.1);
    CHECK(results[4].p_value > 0.1);
    CHECK(results[5].p_value > 0.1);
}
