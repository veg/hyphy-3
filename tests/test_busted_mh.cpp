#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/busted.hpp"

#include <nlohmann/json.hpp>
#include <iostream>
#include <cmath>

using namespace hyphy::core;
using namespace hyphy::analyses;

TEST_CASE("BUSTED-MH: ADH Dataset (Double+Triple Multi-Hit Substitutions)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/tests/data/adh.nex");
    REQUIRE(!aln.embedded_tree_newick.empty());
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);

    BUSTEDSettings settings;
    settings.multiple_hits = "Double+Triple";
    settings.refine_branch_lengths = true;

    auto res = busted.run(settings);

    std::cout << "\n[BUSTED-MH ADH Double+Triple Results]\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n"
              << "  Unc Delta (Double) : " << res.unconstrained.delta << " (fraction: "
              << res.unconstrained.frac_delta * 100.0 << "%)\n"
              << "  Unc Psi (Triple)   : " << res.unconstrained.psi << " (fraction: "
              << res.unconstrained.frac_psi * 100.0 << "%)\n"
              << "  Con Delta (Double) : " << res.constrained.delta << " (fraction: "
              << res.constrained.frac_delta * 100.0 << "%)\n"
              << "  Con Psi (Triple)   : " << res.constrained.psi << " (fraction: "
              << res.constrained.frac_psi * 100.0 << "%)\n"
              << "  Runtime            : " << res.runtime_seconds << " seconds\n";

    CHECK(res.unconstrained.log_likelihood == doctest::Approx(-4703.2).epsilon(0.01));
    CHECK(res.constrained.log_likelihood == doctest::Approx(-4710.28).epsilon(0.01));
    CHECK(res.unconstrained.psi >= 0.0);
    CHECK(res.constrained.delta > 0.10);
    CHECK(res.constrained.frac_delta > 0.02);
    CHECK(res.constrained.psi > 0.05);
    CHECK(res.constrained.frac_psi > 0.002);
    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);

    // Validate Legacy JSON serialization
    auto json_out = busted.to_legacy_json(res);
    CHECK(json_out.contains("fits"));
    const auto& unc = json_out["fits"]["Unconstrained model"];
    CHECK(unc.contains("rate at which 2 nucleotides are changed instantly within a single codon"));
    CHECK(unc.contains("rate at which 3 nucleotides are changed instantly within a single codon"));
    CHECK(unc.contains("Fraction of subs rate at which 2 nucleotides are changed instantly within a single codon"));
    CHECK(unc.contains("Fraction of subs rate at which 3 nucleotides are changed instantly within a single codon"));

    const auto& con = json_out["fits"]["Constrained model"];
    CHECK(con.contains("rate at which 2 nucleotides are changed instantly within a single codon"));
    CHECK(con.contains("rate at which 3 nucleotides are changed instantly within a single codon"));

    // Validate Modern JSON serialization
    auto mod_out = busted.to_modern_json(res);
    CHECK(mod_out.contains("model_fits"));
    CHECK(mod_out["model_fits"]["unconstrained"].contains("multiple_hits"));
    CHECK(mod_out["model_fits"]["unconstrained"]["multiple_hits"].contains("double_hit_rate"));
    CHECK(mod_out["model_fits"]["unconstrained"]["multiple_hits"].contains("triple_hit_rate"));
}

TEST_CASE("BUSTED-MH: CD2 Dataset (Double vs Double+Triple)") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.fna");
    auto tree = Tree::from_newick_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.nwk");

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);

    // 1. Double only
    BUSTEDSettings settings_double;
    settings_double.multiple_hits = "Double";
    auto res_double = busted.run(settings_double);

    std::cout << "\n[BUSTED-MH CD2 Double Results]\n"
              << "  Unconstrained ln L : " << res_double.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res_double.constrained.log_likelihood << "\n"
              << "  Delta (Double)     : " << res_double.unconstrained.delta << " (fraction: "
              << res_double.unconstrained.frac_delta * 100.0 << "%)\n"
              << "  Psi (Triple)       : " << res_double.unconstrained.psi << "\n";

    CHECK(res_double.unconstrained.delta >= 0.0);
    CHECK(res_double.unconstrained.psi == 0.0);
    CHECK(res_double.unconstrained.frac_psi == 0.0);

    // 2. Double+Triple
    BUSTEDSettings settings_dt;
    settings_dt.multiple_hits = "Double+Triple";
    auto res_dt = busted.run(settings_dt);

    std::cout << "\n[BUSTED-MH CD2 Double+Triple Results]\n"
              << "  Unconstrained ln L : " << res_dt.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res_dt.constrained.log_likelihood << "\n"
              << "  Delta (Double)     : " << res_dt.unconstrained.delta << "\n"
              << "  Psi (Triple)       : " << res_dt.unconstrained.psi << " (fraction: "
              << res_dt.unconstrained.frac_psi * 100.0 << "%)\n";

    CHECK(res_dt.unconstrained.delta >= 0.0);
    CHECK(res_dt.unconstrained.psi >= 0.0);
    // Double+Triple likelihood should be at least as high as Double only (or within optimization tolerance)
    CHECK(res_dt.unconstrained.log_likelihood >= res_double.unconstrained.log_likelihood - 0.5);
}

TEST_CASE("BUSTED-MH: Invariance with Multiple Hits = None") {
    auto aln = Alignment::load("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.fna");
    auto tree = Tree::from_newick_file("/Users/sergei/Development/hyphy/hyphy-next/benchmarks/data/cd2.nwk");

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);

    BUSTEDSettings default_settings;
    auto res_default = busted.run(default_settings);

    BUSTEDSettings none_settings;
    none_settings.multiple_hits = "None";
    auto res_none = busted.run(none_settings);

    CHECK(res_none.unconstrained.log_likelihood == doctest::Approx(res_default.unconstrained.log_likelihood).epsilon(1e-6));
    CHECK(res_none.constrained.log_likelihood == doctest::Approx(res_default.constrained.log_likelihood).epsilon(1e-6));
    CHECK(res_none.unconstrained.delta == 0.0);
    CHECK(res_none.unconstrained.psi == 0.0);
}
