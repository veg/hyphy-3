#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/relax.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using namespace hyphy::core;
using namespace hyphy::analyses;

TEST_CASE("RELAX: Wertheim et al. 2015 Fig 4E Benchmark Parity (33 taxa, 286 codons)") {
    std::vector<std::string> candidates = {
        "tests/data/Fig4E.nex",
        "../tests/data/Fig4E.nex",
        "/Users/sergei/Development/hyphy-3/tests/data/Fig4E.nex",
        "/Users/sergei/Development/hyphy/tests/data/Fig4E.nex"
    };
    std::string path;
    for (const auto& c : candidates) {
        std::ifstream f(c);
        if (f.good()) { path = c; break; }
    }
    REQUIRE(!path.empty());
    auto aln = Alignment::load(path);
    REQUIRE(!aln.embedded_tree_newick.empty());
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    RELAXSettings settings;
    settings.p_value_threshold = 0.05;

    auto relax = RELAXAnalyzer::create(tree, aln, settings);

    // 1. Check branch classification from Newick {T} tags
    CHECK(relax.test_branch_names.size() == 17);
    CHECK(relax.ref_branch_names.size() == 47);

    // 2. Run RELAX analysis
    auto res = relax.run();

    std::cout << "\n=======================================================\n"
              << "   RELAX Fig4E Benchmark Validation                     \n"
              << "=======================================================\n"
              << "  GTR ln L              : " << res.gtr_log_likelihood << "\n"
              << "  MG94 Sep Rates ln L   : " << res.mg94_log_likelihood << "\n"
              << "    omega (Reference)   : " << res.mg94_omega_R << "\n"
              << "    omega (Test)        : " << res.mg94_omega_T << "\n"
              << "  RELAX Alternative ln L: " << res.alternative_fit.log_likelihood << "\n"
              << "  RELAX Null ln L       : " << res.null_fit.log_likelihood << "\n"
              << "  K parameter           : " << res.k << "\n"
              << "  LRT statistic         : " << res.lrt << "\n"
              << "  p-value               : " << res.p_value << "\n"
              << "  Runtime               : " << res.runtime_seconds << "s\n"
              << "=======================================================\n";

    // 3. Validation against HyPhy 2.5.101 (`hyphy relax --alignment Fig4E.nex --test T --reference R
    //    --models Minimal --rates 3 --srv No`; two independent runs, which differ slightly because
    //    HyPhy seeds its distribution grid randomly):
    //      GTR -5317.276 | MG94 -4985.166 | alt -4982.997 / -4982.983 | null -5044.587 / -5044.663
    //      LRT 123.180 / 123.361 | K 0.0000
    //    hyphy-3 (default: shared normalization, branch-site mixture, joint optimization) gives
    //    GTR -5317.18, MG94 -4959.81, alt -4957.65, null -5017.47, LRT 119.6, K 0.0000.
    //    The absolute codon log-likelihoods are ~25 units higher than HyPhy 2.5 on this alignment
    //    only (CD2 / adh agree within ~2.6): Fig4E has 49 partially-resolved codons (e.g. "AC-"),
    //    which HyPhy sums over compatible codons while hyphy-3's alignment encoding treats them as
    //    fully missing. Confirmed by masking those 49 codons to "---": HyPhy 2.5 then gives
    //    MG94 -4959.775 / null -5017.423 vs hyphy-3 -4959.720 / -5017.294. The offset cancels in
    //    likelihood differences, which are checked tightly.
    // GTR: HyPhy -5317.276, hyphy-3 -5317.18
    CHECK(std::abs(res.gtr_log_likelihood - (-5317.276)) < 0.5);

    // MG94 separate rates: HyPhy omega_R ~ 0.11, omega_T ~ 0.68
    CHECK(res.mg94_omega_R < 0.25);
    CHECK(res.mg94_omega_T > 0.40);
    CHECK(res.mg94_omega_T > res.mg94_omega_R);

    // Nesting: alternative >= null; alternative >= MG94 (K and the 3-class mixture nest it)
    CHECK(res.alternative_fit.log_likelihood >= res.null_fit.log_likelihood);
    CHECK(res.alternative_fit.log_likelihood >= res.mg94_log_likelihood - 0.5);

    // Gain of the RELAX alternative over MG94: HyPhy 2.5 = 2.18 (-4982.99 - -4985.17)
    CHECK(std::abs((res.alternative_fit.log_likelihood - res.mg94_log_likelihood) - 2.18) < 1.5);

    // Absolute logL: never worse than HyPhy 2.5 (see the partial-codon offset note above)
    CHECK(res.alternative_fit.log_likelihood > -4983.0 - 0.5);
    CHECK(res.null_fit.log_likelihood > -5044.66 - 0.5);
    CHECK(res.null_fit.log_likelihood < -5000.0);

    // Selection relaxation: HyPhy 2.5 infers K = 0.0000
    CHECK(res.k < 0.01);
    CHECK(res.is_relaxed);
    CHECK(res.is_significant);

    // Test statistic: HyPhy 2.5 LRT = 123.18-123.36 (hyphy-3: 119.6)
    CHECK(std::abs(res.lrt - 123.27) < 8.0);
    CHECK(res.p_value < 1e-15);

    // 4. Test Legacy JSON export
    auto legacy_out = res.to_legacy_json(tree, aln);
    CHECK(legacy_out.contains("analysis"));
    CHECK(legacy_out.contains("fits"));
    CHECK(legacy_out["fits"].contains("Nucleotide GTR"));
    CHECK(legacy_out["fits"].contains("MG94xREV with separate rates for branch sets"));
    CHECK(legacy_out["fits"].contains("RELAX alternative"));
    CHECK(legacy_out["fits"].contains("RELAX null"));
    CHECK(legacy_out.contains("test results"));
    CHECK(legacy_out["test results"]["LRT"] == res.lrt);
    CHECK(legacy_out["test results"]["p-value"] == res.p_value);
    CHECK(legacy_out["test results"]["relaxation or intensification parameter"] == res.k);
    CHECK(legacy_out.contains("tested"));
    CHECK(legacy_out.contains("branch attributes"));

    // Test Modern JSON export
    auto modern_out = res.to_modern_json(tree, aln);
    CHECK(modern_out.contains("statistical_tests"));
    CHECK(modern_out.contains("model_fits"));
    CHECK(modern_out["model_fits"].contains("alternative"));
    CHECK(modern_out["model_fits"].contains("null"));
    CHECK(modern_out["statistical_tests"]["hypothesis_test"]["statistic_value"] == res.lrt);
    CHECK(modern_out["statistical_tests"]["hypothesis_test"]["p_value"] == res.p_value);
    CHECK(modern_out["statistical_tests"]["hypothesis_test"]["k_parameter"] == res.k);
}

TEST_CASE("RELAX: User-specified test branch regex") {
    std::vector<std::string> candidates = {
        "tests/data/Fig4E.nex",
        "../tests/data/Fig4E.nex",
        "/Users/sergei/Development/hyphy-3/tests/data/Fig4E.nex",
        "/Users/sergei/Development/hyphy/tests/data/Fig4E.nex"
    };
    std::string path;
    for (const auto& c : candidates) {
        std::ifstream f(c);
        if (f.good()) { path = c; break; }
    }
    REQUIRE(!path.empty());
    auto aln = Alignment::load(path);
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    RELAXSettings settings;
    settings.test_branch_regex = ".*PSEUDOGENE.*";

    auto relax = RELAXAnalyzer::create(tree, aln, settings);

    // All pseudogenes should be marked as test
    CHECK(relax.test_branch_names.size() >= 12);
    for (const auto& name : relax.test_branch_names) {
        CHECK(name.find("PSEUDOGENE") != std::string::npos);
    }
}

TEST_CASE("Issue #4: RELAX case-insensitive {TEST}/{REFERENCE} tags and Unclassified partition") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    // Tag Human and Chimp as {TEST}, Baboon and RhMonkey as {REFERENCE}, leave others untagged (Unclassified)
    std::string nwk = "((((Pig:0.147969,Cow:0.213430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,((RhMonkey{REFERENCE}:0.002015,Baboon{REFERENCE}:0.003108):0.022733,(Human{TEST}:0.004349,Chimp{TEST}:0.000799):0.011873):0.101856):0.340802,Rat:0.050958,Mouse:0.097950);";

    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick(nwk);

    RELAXSettings settings;
    auto relax = RELAXAnalyzer::create(tree, aln, settings);

    CHECK(relax.test_branch_names.size() == 2);
    CHECK(relax.ref_branch_names.size() == 2);
    CHECK(!relax.unclassified_branch_names.empty());
}

TEST_CASE("RELAX: --test/--reference resolution (tag -> exact names -> regex)") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    // Group1 of hbltests CD2.nex as {T}: RhMonkey, Human, Chimp and the Human+Chimp clade (4 branches);
    // everything else {R}. Several names contain 't'/'T' (Cat, Rat, Horse, Mouse...).
    std::string nwk = "((((Pig{R}:0.147969,Cow{R}:0.213430){R}:0.085099,Horse{R}:0.165787,Cat{R}:0.264806){R}:0.058611,"
                      "((RhMonkey{T}:0.002015,Baboon{R}:0.003108){R}:0.022733,(Human{T}:0.004349,Chimp{T}:0.000799){T}:0.011873){R}:0.101856){R}:0.340802,"
                      "Rat{R}:0.050958,Mouse{R}:0.097950);";
    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick(nwk);

    auto count_test = [&](const std::string& spec, bool force_regex = false) {
        RELAXSettings s;
        s.test_branch_regex = spec;
        s.test_is_regex = force_regex;
        s.reference_branch_regex = "R";
        return RELAXAnalyzer::create(tree, aln, s);
    };

    // (a) tag present in the tree: exactly the 4 tagged branches, not every name containing 't'
    auto a = count_test("T");
    CHECK(a.test_branch_names.size() == 4);
    CHECK(a.ref_branch_names.size() == 12);
    CHECK(a.unclassified_branch_names.empty());
    for (const auto& nm : a.test_branch_names) {
        CHECK((nm == "RhMonkey" || nm == "Human" || nm == "Chimp" || a.branch_class.at(nm) == "Test"));
        CHECK(nm != "Cat");
        CHECK(nm != "Rat");
    }
    // braces and alias: {TEST} resolves to the {T} group
    CHECK(count_test("{T}").test_branch_names.size() == 4);
    CHECK(count_test("TEST").test_branch_names.size() == 4);

    // (b) comma-separated exact names
    auto b = count_test("Human,Chimp");
    CHECK(b.test_branch_names.size() == 2);

    // (c) explicit regex via prefix or flag
    CHECK(count_test("regex:^C").test_branch_names.size() == 3);   // Cow, Cat, Chimp
    CHECK(count_test("^C", true).test_branch_names.size() == 3);
    // fallback regex (with a warning) only when neither a tag nor an exact name matches
    CHECK(count_test("monkey").test_branch_names.size() == 1);      // RhMonkey (case-insensitive)
}

TEST_CASE("RELAX: Fig4E --test T selects the 17 {T} branches") {
    std::vector<std::string> candidates = {
        "tests/data/Fig4E.nex",
        "../tests/data/Fig4E.nex",
        "/Users/sergei/Development/hyphy-3/tests/data/Fig4E.nex"
    };
    std::string path;
    for (const auto& c : candidates) {
        std::ifstream f(c);
        if (f.good()) { path = c; break; }
    }
    REQUIRE(!path.empty());
    auto aln = Alignment::load(path);
    auto tree = Tree::from_newick(aln.embedded_tree_newick);
    RELAXSettings s;
    s.test_branch_regex = "T";
    s.reference_branch_regex = "R";
    auto relax = RELAXAnalyzer::create(tree, aln, s);
    CHECK(relax.test_branch_names.size() == 17);
    CHECK(relax.ref_branch_names.size() == 47);
}

TEST_CASE("RELAX: normalization / mixture-level invariants") {
    std::string aln_path = HYPHY3_SOURCE_DIR "/benchmarks/data/cd2_reduced.fna";
    std::string nwk = "((((Pig:0.147969,Cow:0.213430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,((RhMonkey:0.002015,Baboon:0.003108):0.022733,(Human{TEST}:0.004349,Chimp{TEST}:0.000799):0.011873):0.101856):0.340802,Rat:0.050958,Mouse:0.097950);";
    auto aln = Alignment::from_fasta(aln_path);
    auto tree = Tree::from_newick(nwk);

    // Defaults are the HyPhy 2.5 BS-REL semantics.
    RELAXSettings s_def;
    CHECK(s_def.normalization == RELAXNormalization::SharedMixture);
    CHECK(s_def.mixture_level == RELAXMixtureLevel::BranchSite);

    RELAXSettings s_pc;   // legacy: per-class normalization, site-level mixture
    s_pc.normalization = RELAXNormalization::PerClass;
    s_pc.mixture_level = RELAXMixtureLevel::Site;
    RELAXSettings s_sh;   // shared normalization, site-level mixture
    s_sh.normalization = RELAXNormalization::SharedMixture;
    s_sh.mixture_level = RELAXMixtureLevel::Site;
    auto r_pc = RELAXAnalyzer::create(tree, aln, s_pc);
    auto r_sh = RELAXAnalyzer::create(tree, aln, s_sh);

    // Single-class mixture: both normalizations coincide.
    Scalar l_pc1 = r_pc.evaluate_relax_mixture(r_pc.tree, {0.3, 0.3, 0.3}, {1.0, 0.0, 0.0}, 0.7, 1.0);
    Scalar l_sh1 = r_sh.evaluate_relax_mixture(r_sh.tree, {0.3, 0.3, 0.3}, {1.0, 0.0, 0.0}, 0.7, 1.0);
    CHECK(l_pc1 == doctest::Approx(l_sh1).epsilon(1e-12));

    // Genuine mixture: they differ.
    Scalar l_pc3 = r_pc.evaluate_relax_mixture(r_pc.tree, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7, 1.0);
    Scalar l_sh3 = r_sh.evaluate_relax_mixture(r_sh.tree, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7, 1.0);
    CHECK(std::isfinite(l_sh3));
    CHECK(std::abs(l_pc3 - l_sh3) > 1e-6);

    // Analytical branch gradient matches finite differences under SharedMixture.
    std::vector<Scalar> bl(r_sh.tree.num_nodes(), 0.0);
    for (const auto& n : r_sh.tree.nodes) if (n.id != r_sh.tree.root_id) bl[n.id] = n.branch_length;
    auto [ll, grad] = r_sh.compute_relax_branch_gradients(bl, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7);
    CHECK(ll == doctest::Approx(l_sh3).epsilon(1e-6) );
    int probe = r_sh.tree.nodes[0].id == r_sh.tree.root_id ? 1 : 0;
    for (const auto& n : r_sh.tree.nodes) {
        if (n.id == r_sh.tree.root_id) continue;
        if (n.name == "Human") { probe = n.id; break; }
    }
    const Scalar h = 1e-6;
    auto bp = bl; bp[probe] += h;
    auto bm = bl; bm[probe] -= h;
    Scalar lp = r_sh.compute_relax_branch_gradients(bp, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7).first;
    Scalar lm = r_sh.compute_relax_branch_gradients(bm, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7).first;
    CHECK(grad[probe] == doctest::Approx((lp - lm) / (2.0 * h)).epsilon(1e-3));

    // Branch-site mixture (HyPhy 2.5 BS-REL semantics) + shared normalization.
    RELAXSettings s_bs = s_sh;
    s_bs.mixture_level = RELAXMixtureLevel::BranchSite;
    auto r_bs = RELAXAnalyzer::create(tree, aln, s_bs);
    Scalar l_bs1 = r_bs.evaluate_relax_mixture(r_bs.tree, {0.3, 0.3, 0.3}, {1.0, 0.0, 0.0}, 0.7, 1.0);
    CHECK(l_bs1 == doctest::Approx(l_sh1).epsilon(1e-10));
    Scalar l_bs3 = r_bs.evaluate_relax_mixture(r_bs.tree, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7, 1.0);
    CHECK(std::isfinite(l_bs3));
    CHECK(std::abs(l_bs3 - l_sh3) > 1e-6);
    auto [ll_bs, grad_bs] = r_bs.compute_relax_branch_gradients(bl, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7);
    CHECK(ll_bs == doctest::Approx(l_bs3).epsilon(1e-6));
    Scalar lp_bs = r_bs.compute_relax_branch_gradients(bp, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7).first;
    Scalar lm_bs = r_bs.compute_relax_branch_gradients(bm, {0.05, 0.5, 4.0}, {0.5, 0.4, 0.1}, 0.7).first;
    CHECK(grad_bs[probe] == doctest::Approx((lp_bs - lm_bs) / (2.0 * h)).epsilon(1e-3));
}

