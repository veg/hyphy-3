#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/analyses/busted.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using namespace hyphy::core;
using namespace hyphy::analyses;

TEST_CASE("BUSTED: CD2 Dataset (10 taxa, 187 codons) & Ground-Truth Parity") {
    auto aln = Alignment::load(HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.fna");
    auto tree = Tree::from_newick_file(HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.nwk");

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    auto res = busted.run();

    std::cout << "\n[BUSTED CD2 Results (SRV=off)]\n"
              << "  Global MG94 ln L   : " << busted.mg94_log_l << " (omega = " << busted.mg94_omega << ")\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n"
              << "  Alt omegas         : [" << res.unconstrained.test_distribution.omegas[0] << ", "
                                          << res.unconstrained.test_distribution.omegas[1] << ", "
                                          << res.unconstrained.test_distribution.omegas[2] << "]\n"
              << "  Alt weights        : [" << res.unconstrained.test_distribution.weights[0] << ", "
                                          << res.unconstrained.test_distribution.weights[1] << ", "
                                          << res.unconstrained.test_distribution.weights[2] << "]\n"
              << "  Null omegas        : [" << res.constrained.test_distribution.omegas[0] << ", "
                                          << res.constrained.test_distribution.omegas[1] << ", "
                                          << res.constrained.test_distribution.omegas[2] << "]\n"
              << "  Runtime            : " << res.runtime_seconds << " seconds\n";

    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);
    CHECK(res.unconstrained.test_distribution.omegas[2] >= 1.0);
    CHECK(res.constrained.test_distribution.omegas[2] == 1.0);
    CHECK(res.evidence_ratios.size() == 187);

    // Load bundled HyPhy 2.5 ground-truth reference (cd2_busted.json has srv: "Yes")
    std::ifstream gt_file(HYPHY3_SOURCE_DIR "/benchmarks/ground_truth/cd2_busted.json");
    REQUIRE(gt_file.is_open());
    nlohmann::json gt;
    gt_file >> gt;

    Scalar gt_gtr_ll = gt["fits"]["Nucleotide GTR"]["Log Likelihood"].get<Scalar>();
    Scalar gt_mg94_ll = gt["fits"]["MG94xREV with separate rates for branch sets"]["Log Likelihood"].get<Scalar>();
    Scalar gt_unc_ll = gt["fits"]["Unconstrained model"]["Log Likelihood"].get<Scalar>();
    Scalar gt_con_ll = gt["fits"]["Constrained model"]["Log Likelihood"].get<Scalar>();
    Scalar gt_lrt = gt["test results"]["LRT"].get<Scalar>();
    Scalar gt_pval = gt["test results"]["p-value"].get<Scalar>();

    // Verify GTR and full Codon MG94xREV baseline parity against ground truth
    CHECK(std::abs(busted.gtr_log_l - gt_gtr_ll) < 1.0);
    CHECK(std::abs(busted.mg94_log_l - gt_mg94_ll) < 0.05);

    // Run with SRV=true to match cd2_busted.json settings ("srv": "Yes")
    BUSTEDSettings srv_settings;
    srv_settings.srv = true;
    srv_settings.num_rate_classes = 3;
    srv_settings.num_syn_rate_classes = 3;
    srv_settings.refine_branch_lengths = true;
    auto res_srv = busted.run(srv_settings);

    std::cout << "\n[BUSTED-S CD2 Ground-Truth Parity (SRV=on)]\n"
              << "  HyPhy 3 Unconstrained ln L : " << res_srv.unconstrained.log_likelihood
              << " (2.5 GT: " << gt_unc_ll << ")\n"
              << "  HyPhy 3 Constrained ln L   : " << res_srv.constrained.log_likelihood
              << " (2.5 GT: " << gt_con_ll << ")\n"
              << "  HyPhy 3 LRT                : " << res_srv.lrt
              << " (2.5 GT: " << gt_lrt << ")\n"
              << "  HyPhy 3 p-value            : " << res_srv.p_value
              << " (2.5 GT: " << gt_pval << ")\n";

    CHECK(std::abs(res_srv.unconstrained.log_likelihood - gt_unc_ll) < 2.5);
    CHECK(std::abs(res_srv.constrained.log_likelihood - gt_con_ll) < 2.5);
    CHECK(std::abs(res_srv.lrt - gt_lrt) < 1.5);
    CHECK(res_srv.p_value < 0.01);

    // Test Legacy JSON export
    auto legacy_out = busted.to_legacy_json(res);
    CHECK(legacy_out.contains("test results"));
    CHECK(legacy_out.contains("fits"));
    CHECK(legacy_out["fits"].contains("Unconstrained model"));
    CHECK(legacy_out["fits"].contains("Constrained model"));

    // Test Modern JSON export
    auto modern_out = busted.to_modern_json(res);
    CHECK(modern_out.contains("statistical_tests"));
    CHECK(modern_out.contains("model_fits"));
    CHECK(modern_out["model_fits"].contains("unconstrained"));
    CHECK(modern_out["model_fits"].contains("constrained"));
    CHECK(modern_out.contains("site_results"));
}

TEST_CASE("BUSTED: ADH Dataset (23 taxa, 254 codons)") {
    auto aln = Alignment::load(HYPHY3_SOURCE_DIR "/benchmarks/data/adh.nex");
    REQUIRE(!aln.embedded_tree_newick.empty());
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    auto res = busted.run();

    std::cout << "\n[BUSTED ADH Results]\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n"
              << "  Alt omegas         : [" << res.unconstrained.test_distribution.omegas[0] << ", "
                                          << res.unconstrained.test_distribution.omegas[1] << ", "
                                          << res.unconstrained.test_distribution.omegas[2] << "]\n"
              << "  Alt weights        : [" << res.unconstrained.test_distribution.weights[0] << ", "
                                          << res.unconstrained.test_distribution.weights[1] << ", "
                                          << res.unconstrained.test_distribution.weights[2] << "]\n"
              << "  Null omegas        : [" << res.constrained.test_distribution.omegas[0] << ", "
                                          << res.constrained.test_distribution.omegas[1] << ", "
                                          << res.constrained.test_distribution.omegas[2] << "]\n"
              << "  Runtime            : " << res.runtime_seconds << " seconds\n";

    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);
    CHECK(res.unconstrained.test_distribution.omegas[2] >= 1.0);
    CHECK(res.constrained.test_distribution.omegas[2] == 1.0);
    CHECK(res.evidence_ratios.size() == 254);
    CHECK(!res.unconstrained.branch_lengths.empty());
}

TEST_CASE("BUSTED: Automatic K Model Selection (AICc Step-Up)") {
    auto aln = Alignment::load(HYPHY3_SOURCE_DIR "/benchmarks/data/adh.nex");
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    BUSTEDSettings settings;
    settings.auto_select_k = true;
    settings.max_k = 3;
    settings.refine_branch_lengths = true;

    auto res = busted.run(settings);

    std::cout << "\n[BUSTED Auto-K Results]\n"
              << "  Optimal K          : " << res.optimal_k << "\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n";

    // For ADH, K=2 is strongly preferred over K=1 (delta AICc ~ -91), and K=3 is penalized (+15 AICc)
    CHECK(res.optimal_k == 2);
    CHECK(res.unconstrained.test_distribution.omegas.size() == 2);
    CHECK(res.constrained.test_distribution.omegas.size() == 2);
    CHECK(res.unconstrained.test_distribution.omegas[1] >= 1.0);
    CHECK(res.constrained.test_distribution.omegas[1] == 1.0);
    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value <= 0.05); // Positive selection still robustly detected under optimal K=2
}

TEST_CASE("BUSTED-S: Synonymous Rate Variation (SRV) on ADH") {
    auto aln = Alignment::load(HYPHY3_SOURCE_DIR "/benchmarks/data/adh.nex");
    auto tree = Tree::from_newick(aln.embedded_tree_newick);

    auto busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    BUSTEDSettings settings;
    settings.srv = true;
    settings.num_rate_classes = 3;
    settings.num_syn_rate_classes = 3;
    settings.refine_branch_lengths = true;

    auto res = busted.run(settings);

    std::cout << "\n[BUSTED-S (SRV) ADH Results]\n"
              << "  Unconstrained ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained ln L   : " << res.constrained.log_likelihood << "\n"
              << "  LRT                : " << res.lrt << "\n"
              << "  p-value            : " << res.p_value << "\n";
    if (!res.unconstrained.test_distribution.syn_rates.empty()) {
        std::cout << "  Syn rates          : ["
                  << res.unconstrained.test_distribution.syn_rates[0] << ", "
                  << res.unconstrained.test_distribution.syn_rates[1] << ", "
                  << res.unconstrained.test_distribution.syn_rates[2] << "]\n"
                  << "  Syn weights        : ["
                  << res.unconstrained.test_distribution.syn_weights[0] << ", "
                  << res.unconstrained.test_distribution.syn_weights[1] << ", "
                  << res.unconstrained.test_distribution.syn_weights[2] << "]\n";
    }

    CHECK(res.lrt >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);
    CHECK(res.unconstrained.test_distribution.syn_rates.size() == 3);
    CHECK(res.unconstrained.test_distribution.syn_weights.size() == 3);

    // Verify mean synonymous rate is ~ 1.0
    Scalar mean_syn = 0.0;
    for (size_t m = 0; m < 3; ++m) {
        mean_syn += res.unconstrained.test_distribution.syn_rates[m] * res.unconstrained.test_distribution.syn_weights[m];
    }
    CHECK(std::abs(mean_syn - 1.0) < 0.05);

    auto j = busted.to_legacy_json(res);
    CHECK(j["analysis"]["settings"]["srv"] == "Yes");
    CHECK(j["fits"]["Unconstrained model"]["Rate Distributions"].contains("Synonymous site-to-site rates"));

    auto j_mod = busted.to_modern_json(res);
    CHECK(j_mod["analysis"]["settings"]["srv"] == true);
    CHECK(j_mod["model_fits"]["unconstrained"]["rate_distributions"].contains("synonymous"));
}

TEST_CASE("BUSTED-E: Saturation Error Sink (Option A) Quarantines Misalignment Artifacts") {
    // 1. Construct a strictly purifying/neutral 6-taxon, 100-codon alignment on a tree with a very short sister pair (T1, T2)
    std::string nwk = "((T1:0.002,T2:0.002)N1:0.05,(T3:0.04,T4:0.04)N2:0.05,(T5:0.04,T6:0.04)N3:0.05);";
    auto tree = Tree::from_newick(nwk);

    // Base ancestral sequence of 300 4-fold degenerate codons (GCT=Ala, GTT=Val, CCT=Pro, GGT=Gly).
    // Note: 300 codons (298 invariant codons between sister tips T1 and T2) are required so that under
    // converged Constrained Null branch-length optimization (F-001, F-005, F-007, F-067), the Null model
    // cannot freely stretch terminal branch T1 to absorb the two 3-nt errors without paying a large
    // synonymous survival penalty across the invariant sites, ensuring the true MLE LRT of Standard BUSTED has p < 0.05.
    std::vector<std::string> base_codons(300);
    for (size_t s = 0; s < 300; ++s) {
        if (s % 4 == 0) base_codons[s] = "GCT";
        else if (s % 4 == 1) base_codons[s] = "GTT";
        else if (s % 4 == 2) base_codons[s] = "CCT";
        else base_codons[s] = "GGT";
    }

    auto make_seq = [&](int variant) {
        std::string seq;
        for (size_t s = 0; s < 300; ++s) {
            std::string c = base_codons[s];
            // Add realistic synonymous third-position variation across clades so omega < 1
            if (variant == 1 && (s % 5 == 0)) c[2] = 'C'; // T1, T2 share synonymous T->C
            if (variant == 2 && (s % 7 == 0)) c[2] = 'C'; // T3, T4 share synonymous T->C
            if (variant == 3 && (s % 3 == 0)) c[2] = 'A'; // T5, T6 share synonymous T->A
            seq += c;
        }
        return seq;
    };

    std::string s1 = make_seq(1);
    std::string s2 = make_seq(1);
    std::string s3 = make_seq(2);
    std::string s4 = make_seq(2);
    std::string s5 = make_seq(3);
    std::string s6 = make_seq(3);

    // Inject two 3-nucleotide non-synonymous frameshift/misalignment errors into T1 (branch length 0.002!)
    // At codon 10 (index 9): GCT (Ala) -> TGG (Trp, 3-nt change)
    // At codon 50 (index 49): GTT (Val) -> AAA (Lys, 3-nt change)
    s1.replace(9 * 3, 3, "TGG");
    s1.replace(49 * 3, 3, "AAA");

    Alignment corrupted_aln;
    corrupted_aln.add_sequence("T1", s1);
    corrupted_aln.add_sequence("T2", s2);
    corrupted_aln.add_sequence("T3", s3);
    corrupted_aln.add_sequence("T4", s4);
    corrupted_aln.add_sequence("T5", s5);
    corrupted_aln.add_sequence("T6", s6);
    corrupted_aln.finalize();

    // Run standard BUSTED (no error sink): the two 3-nt changes on branch T1 (t=0.002) force a massive omega_3 blow-up
    auto analyzer = BUSTEDAnalyzer::create_and_fit(tree, corrupted_aln);
    BUSTEDSettings std_settings;
    std_settings.error_sink = false;
    auto res_std = analyzer.run(std_settings);

    // Run BUSTED-E (--error-sink): the saturation replacement kernel absorbs the two 3-nt artifacts
    BUSTEDSettings sink_settings;
    sink_settings.error_sink = true;
    auto res_sink = analyzer.run(sink_settings);

    std::cout << "\n[BUSTED vs BUSTED-E on Corrupted Purifying Alignment]\n"
              << "  Standard BUSTED : LRT = " << res_std.lrt << ", p = " << res_std.p_value
              << ", omega_3 = " << res_std.unconstrained.test_distribution.omegas.back() << "\n"
              << "  BUSTED-E (Sink) : LRT = " << res_sink.lrt << ", p = " << res_sink.p_value
              << ", omega_3 = " << res_sink.unconstrained.test_distribution.omegas.back()
              << ", eps_sink = " << res_sink.unconstrained.error_sink_weight
              << ", flagged branch-sites = " << res_sink.flagged_error_events.size() << "\n";

    for (const auto& ev : res_sink.flagged_error_events) {
        std::cout << "    Flagged: codon " << ev.site << " on branch " << ev.branch_name
                  << " (BF = " << ev.bayes_factor << ", P = " << ev.posterior << ")\n";
    }

    // Standard BUSTED is fooled by the 3-nt misalignments on a short branch (p < 0.05)
    CHECK(res_std.p_value < 0.05);

    // BUSTED-E eliminates the false positive (p > 0.10) and flags both corrupted codons (1-based sites 10 and 50 on T1)
    CHECK(res_sink.p_value > 0.10);
    CHECK(res_sink.flagged_error_events.size() >= 2);
    bool found_c10 = false, found_c50 = false;
    for (const auto& ev : res_sink.flagged_error_events) {
        if (ev.site == 10 && ev.branch_name == "T1") found_c10 = true;
        if (ev.site == 50 && ev.branch_name == "T1") found_c50 = true;
    }
    CHECK(found_c10);
    CHECK(found_c50);

    // JSON export verification
    auto j_legacy = analyzer.to_legacy_json(res_sink);
    CHECK(j_legacy["Evidence Ratios"].contains("error_sink_max_bf"));
    auto j_modern = analyzer.to_modern_json(res_sink);
    CHECK(j_modern["analysis"]["settings"]["error_sink"] == true);
    CHECK(j_modern.contains("error_sink_diagnostics"));

    // Verify error-filter alignment masking (HyPhy 2.5 error-filter.bf equivalent)
    auto filtered_aln = analyzer.filter_alignment_by_error_sink(res_sink);
    CHECK(filtered_aln.sequences.size() == 6);
    CHECK(filtered_aln.sequences[0].raw_seq.substr(9 * 3, 3) == "---");
    CHECK(filtered_aln.sequences[0].raw_seq.substr(49 * 3, 3) == "---");
    CHECK(filtered_aln.sequences[1].raw_seq == s2);

    // Running standard BUSTED (no error sink) on the filtered alignment should no longer produce a false positive
    auto filtered_analyzer = BUSTEDAnalyzer::create_and_fit(tree, filtered_aln);
    auto res_filtered = filtered_analyzer.run(std_settings);
    CHECK(res_filtered.p_value > 0.10);
}

TEST_CASE("BUSTED-Auto: 4-Gate Surprisal & Score Model Selector") {
    // 1. Test on ADH: should skip Error Sink, skip Multi-Hits, activate SRV, select K=2, and detect positive selection
    {
        auto aln = Alignment::load(HYPHY3_SOURCE_DIR "/benchmarks/data/adh.nex");
        auto tree = Tree::from_newick(aln.embedded_tree_newick);
        auto analyzer = BUSTEDAnalyzer::create_and_fit(tree, aln);

        BUSTEDSettings auto_settings;
        auto_settings.full_auto = true;
        auto res = analyzer.run(auto_settings);

        const auto& d = res.auto_diagnostics;
        std::cout << "\n[BUSTED-Auto on ADH]\n"
                  << "  Gate 1 (Error Sink) : " << (d.activated_error_sink ? "ON" : "OFF")
                  << " (max surprisal = " << d.max_sink_surprisal_bits << " bits, delta_ll = " << d.sink_delta_ll << ")\n"
                  << "  Gate 2 (Multi-Hits) : " << d.selected_multiple_hits
                  << " (2-hit score = " << d.double_hit_score << ", delta_ll = " << d.double_hit_delta_ll << ")\n"
                  << "  Gate 3 (SRV)        : " << (d.activated_srv ? "ON" : "OFF")
                  << " (M = " << d.selected_syn_rates << ", Z = " << d.srv_overdispersion_z << ", delta_ll = " << d.srv_delta_ll << ")\n"
                  << "  Gate 4 (Omega K)    : K = " << d.selected_k << " (" << d.k_rationale << ")\n"
                  << "  Final LRT           : " << res.lrt << " (p = " << res.p_value << ")\n";

        CHECK(d.enabled == true);
        CHECK(d.activated_error_sink == false);
        CHECK(d.selected_multiple_hits == "None");
        CHECK(d.activated_srv == true);
        CHECK(d.selected_k == 2);
        CHECK(res.p_value < 0.05);

        auto j_mod = analyzer.to_modern_json(res);
        CHECK(j_mod.contains("auto_model_selection"));
    }

    // 2. Test on Corrupted Purifying Alignment: should automatically activate Error Sink, NOT activate Multi-Hits, and eliminate false positive
    {
        std::string nwk = "((T1:0.002,T2:0.002)N1:0.05,(T3:0.04,T4:0.04)N2:0.05,(T5:0.04,T6:0.04)N3:0.05);";
        auto tree = Tree::from_newick(nwk);

        std::vector<std::string> base_codons(100);
        for (size_t s = 0; s < 100; ++s) {
            if (s % 4 == 0) base_codons[s] = "GCT";
            else if (s % 4 == 1) base_codons[s] = "GTT";
            else if (s % 4 == 2) base_codons[s] = "CCT";
            else base_codons[s] = "GGT";
        }
        auto make_seq = [&](int variant) {
            std::string seq;
            for (size_t s = 0; s < 100; ++s) {
                std::string c = base_codons[s];
                if (variant == 1 && (s % 5 == 0)) c[2] = 'C';
                if (variant == 2 && (s % 7 == 0)) c[2] = 'C';
                if (variant == 3 && (s % 3 == 0)) c[2] = 'A';
                seq += c;
            }
            return seq;
        };
        std::string s1 = make_seq(1), s2 = make_seq(1), s3 = make_seq(2), s4 = make_seq(2), s5 = make_seq(3), s6 = make_seq(3);
        s1.replace(9 * 3, 3, "TGG");
        s1.replace(49 * 3, 3, "AAA");

        Alignment corrupted_aln;
        corrupted_aln.add_sequence("T1", s1);
        corrupted_aln.add_sequence("T2", s2);
        corrupted_aln.add_sequence("T3", s3);
        corrupted_aln.add_sequence("T4", s4);
        corrupted_aln.add_sequence("T5", s5);
        corrupted_aln.add_sequence("T6", s6);
        corrupted_aln.finalize();

        auto analyzer = BUSTEDAnalyzer::create_and_fit(tree, corrupted_aln);
        BUSTEDSettings auto_settings;
        auto_settings.full_auto = true;
        auto res = analyzer.run(auto_settings);

        const auto& d = res.auto_diagnostics;
        std::cout << "\n[BUSTED-Auto on Corrupted Purifying Alignment]\n"
                  << "  Gate 1 (Error Sink) : " << (d.activated_error_sink ? "ON" : "OFF")
                  << " (max surprisal = " << d.max_sink_surprisal_bits << " bits, delta_ll = " << d.sink_delta_ll << ")\n"
                  << "  Gate 2 (Multi-Hits) : " << d.selected_multiple_hits
                  << " (quarantined 2-hit score = " << d.double_hit_score << ", delta_ll = " << d.double_hit_delta_ll << ")\n"
                  << "  Gate 3 (SRV)        : " << (d.activated_srv ? "ON" : "OFF") << "\n"
                  << "  Gate 4 (Omega K)    : K = " << d.selected_k << "\n"
                  << "  Final LRT           : " << res.lrt << " (p = " << res.p_value << ")\n";

        CHECK(d.activated_error_sink == true);
        CHECK(d.selected_multiple_hits == "None");
        CHECK(res.p_value > 0.10);
        CHECK(res.flagged_error_events.size() == 2);
    }
}

TEST_CASE("BUSTED-PH: Phenotype-Associated Selection (4-Model Workflow)") {
    auto aln = Alignment::load(HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.fna");
    auto tree = Tree::from_newick_file(HYPHY3_SOURCE_DIR "/benchmarks/data/cd2.nwk");

    auto analyzer = BUSTEDAnalyzer::create_and_fit(tree, aln);
    BUSTEDSettings ph_settings;
    ph_settings.phenotype_mode = true;
    ph_settings.test_branches = "Human,Chimp,Baboon,RhMonkey,Node3,Node2";
    ph_settings.num_rate_classes = 3;

    auto res = analyzer.run(ph_settings);

    std::cout << "\n[BUSTED-PH on CD2 (Primates vs Non-Primates)]\n"
              << "  Unconstrained Two-Set ln L : " << res.unconstrained.log_likelihood << "\n"
              << "  Constrained Test ln L      : " << res.constrained.log_likelihood
              << " (LRT_T = " << res.lrt << ", p_T = " << res.p_value << ")\n"
              << "  Constrained BG ln L        : " << res.constrained_background.log_likelihood
              << " (LRT_B = " << res.background_lrt << ", p_B = " << res.background_p_value
              << ", delta_AICc = " << res.background_delta_aicc << ")\n"
              << "  Shared Dist ln L           : " << res.shared_distributions.log_likelihood
              << " (LRT_diff = " << res.shared_lrt << ", df = " << res.shared_df << ", p_diff = " << res.shared_p_value << ")\n"
              << "  Coupled-k RELAX ln L       : " << res.relax_alternative.log_likelihood
              << " (k = " << res.relax_k << ", LRT_RELAX = " << res.relax_lrt << ", p_RELAX = " << res.relax_p_value
              << ", dir = " << res.relax_direction << ")\n"
              << "  Union-Intersection p_UIT   : " << res.union_intersection_p_value << "\n"
              << "  BG Neutral Supported (AICc): " << (res.background_neutral_supported ? "YES" : "NO") << "\n"
              << "  Phenotype Associated       : " << (res.phenotype_associated ? "YES" : "NO")
              << " (" << res.phenotype_verdict << ")\n";

    CHECK(res.has_background_test == true);
    CHECK(res.has_relax_test == true);
    CHECK(res.shared_df == 5);
    CHECK(res.lrt >= 0.0);
    CHECK(res.background_lrt >= 0.0);
    CHECK(res.shared_lrt >= 0.0);
    CHECK(res.relax_lrt >= 0.0);
    CHECK(res.relax_k >= 0.0);
    CHECK(res.p_value >= 0.0);
    CHECK(res.p_value <= 1.0);
    CHECK(res.background_p_value >= 0.0);
    CHECK(res.background_p_value <= 1.0);
    CHECK(res.shared_p_value >= 0.0);
    CHECK(res.shared_p_value <= 1.0);
    CHECK(res.relax_p_value >= 0.0);
    CHECK(res.relax_p_value <= 1.0);
    CHECK(std::abs(res.union_intersection_p_value - std::max(res.p_value, res.shared_p_value)) < 1e-12);
    CHECK(std::abs(res.background_delta_aicc - (res.constrained_background.aicc - res.unconstrained.aicc)) < 1e-9);
    CHECK(res.unconstrained.has_background_distribution == true);
    CHECK(res.unconstrained.background_distribution.omegas.size() == 3);
    CHECK(res.constrained.test_distribution.omegas.back() == 1.0);
    CHECK(res.constrained_background.background_distribution.omegas.back() == 1.0);
    CHECK(!res.phenotype_verdict.empty());

    // Verify Legacy JSON
    auto j_leg = analyzer.to_legacy_json(res);
    CHECK(j_leg.contains("test results background"));
    CHECK(j_leg.contains("test results shared distributions"));
    CHECK(j_leg.contains("test results relaxation"));
    CHECK(j_leg.contains("BUSTED-PH"));
    CHECK(j_leg["BUSTED-PH"].contains("union_intersection_p_value"));
    CHECK(j_leg["BUSTED-PH"].contains("background_delta_aicc"));
    CHECK(j_leg.contains("phenotype association"));
    CHECK(j_leg["fits"].contains("Constrained model (background)"));
    CHECK(j_leg["fits"].contains("Unconstrained model with shared distributions"));
    CHECK(j_leg["fits"].contains("RELAX alternative (coupled-k)"));
    CHECK(j_leg["fits"]["Unconstrained model"]["Rate Distributions"].contains("Background"));

    // Verify Modern v3 JSON
    auto j_mod = analyzer.to_modern_json(res);
    CHECK(j_mod["analysis"]["id"] == "busted-ph");
    CHECK(j_mod["model_fits"].contains("constrained_background"));
    CHECK(j_mod["model_fits"].contains("shared_distributions"));
    CHECK(j_mod["model_fits"].contains("relax_alternative"));
    CHECK(j_mod["statistical_tests"].contains("background_selection"));
    CHECK(j_mod["statistical_tests"].contains("shared_distributions"));
    CHECK(j_mod["statistical_tests"].contains("relaxation"));
    CHECK(j_mod["statistical_tests"].contains("union_intersection"));
    CHECK(j_mod["statistical_tests"].contains("phenotype_association"));

    // Verify 3-way branch partitioning (Test, Background, Unclassified) via both --comparison and auto-detected {FG}/{BG} tags (Issue #2036)
    {
        // 1. Explicit --branches and --comparison leaving remaining branches in Unclassified
        BUSTEDSettings three_way_settings;
        three_way_settings.phenotype_mode = true;
        three_way_settings.test_branches = "Human,Chimp,Baboon,RhMonkey";
        three_way_settings.comparison_branches = "Mouse,Rat,Horse,Cat";
        three_way_settings.num_rate_classes = 3;

        auto res_3way = analyzer.run(three_way_settings);
        CHECK(analyzer.num_test_branches == 4);
        CHECK(analyzer.num_background_branches == 4);
        CHECK(analyzer.num_unclassified_branches > 0);
        CHECK(res_3way.unconstrained.has_unclassified_distribution == true);
        CHECK(res_3way.unconstrained.unclassified_distribution.omegas.size() == 3);
        CHECK(res_3way.constrained.has_unclassified_distribution == true);
        CHECK(res_3way.constrained_background.has_unclassified_distribution == true);
        CHECK(res_3way.shared_distributions.has_unclassified_distribution == true);

        auto j_leg_3 = analyzer.to_legacy_json(res_3way);
        CHECK(j_leg_3["fits"]["Unconstrained model"]["Rate Distributions"].contains("Unclassified"));
        CHECK(j_leg_3["tested"]["0"]["Cow"] == "unclassified");

        auto j_mod_3 = analyzer.to_modern_json(res_3way);
        CHECK(j_mod_3["model_fits"]["unconstrained"]["rate_distributions"].contains("unclassified"));

        // 2. Auto-detect {FG} and {BG} Newick tags when --branches and --comparison are omitted
        std::string tagged_nwk =
            "((((Human{FG}:0.03,Chimp{FG}:0.02)Node3{FG}:0.03,Baboon{FG}:0.05)Node2:0.08,RhMonkey:0.06)Node1:0.12,"
            "((Mouse{BG}:0.10,Rat{BG}:0.11)Node8{BG}:0.15,(Horse{BG}:0.09,Cat{BG}:0.08)Node11:0.06)Node7:0.05,"
            "(Cow:0.12,Pig:0.11)Node14:0.07);";
        auto tagged_tree = Tree::from_newick(tagged_nwk);
        BUSTEDAnalyzer tagged_analyzer(tagged_tree, aln);
        tagged_analyzer.select_branch_partitions("All", "", true);
        CHECK(tagged_analyzer.num_test_branches == 4);
        CHECK(tagged_analyzer.num_background_branches == 5);
        CHECK(tagged_analyzer.num_unclassified_branches > 0);
    }
}

TEST_CASE("BUSTED: High-Taxon (128 Taxa) Scaling and SRV Underflow Safety") {
    // Build a balanced 128-taxon binary tree (depth 7)
    auto build_balanced_newick = [](auto& self, int start_idx, int count, int depth) -> std::string {
        if (count == 1) {
            return "T" + std::to_string(start_idx) + ":0.03";
        }
        int half = count / 2;
        std::string left = self(self, start_idx, half, depth + 1);
        std::string right = self(self, start_idx + half, count - half, depth + 1);
        if (depth == 0) {
            return "(" + left + "," + right + ");";
        }
        return "(" + left + "," + right + "):0.02";
    };
    std::string nwk = build_balanced_newick(build_balanced_newick, 1, 128, 0);
    auto tree = Tree::from_newick(nwk);

    // Construct a 60-codon purifying/neutral alignment across 128 taxa with clade-correlated synonymous variation
    const std::vector<std::string> base_codons = {"GCT", "GTT", "CCT", "GGT", "ACT", "CTT"};
    const char syn_third[4] = {'T', 'C', 'A', 'G'};
    Alignment aln;
    for (int t = 1; t <= 128; ++t) {
        std::string seq;
        seq.reserve(60 * 3);
        int clade = (t - 1) / 8; // 16 clades of 8 taxa
        for (int s = 0; s < 60; ++s) {
            std::string c = base_codons[s % base_codons.size()];
            // Introduce synonymous variation across clades so sites are variable at 128 taxa
            if ((s + clade) % 3 == 0) {
                c[2] = syn_third[(clade + s) % 4];
            }
            // Occasional conservative nonsynonymous change at every 15th site in 2 clades
            if (s % 15 == 0 && (clade == 2 || clade == 11) && (t % 2 == 0)) {
                c = "GCC";
            }
            seq += c;
        }
        aln.add_sequence("T" + std::to_string(t), seq);
    }
    aln.finalize();

    auto analyzer = BUSTEDAnalyzer::create_and_fit(tree, aln);
    CHECK(std::isfinite(analyzer.gtr_log_l));
    CHECK(std::isfinite(analyzer.mg94_log_l));
    // Full MG94 optimization on 128 taxa should substantially improve over GTR baseline
    CHECK(analyzer.mg94_log_l > analyzer.gtr_log_l);

    // Run with SRV = true on 128 taxa to verify compute_expected_syn_weights and compute_expected_weights do not underflow
    BUSTEDSettings srv_settings;
    srv_settings.srv = true;
    srv_settings.num_syn_rate_classes = 3;
    auto res_srv = analyzer.run(srv_settings);

    std::cout << "\n[BUSTED 128-Taxon Scaling & SRV Underflow Test]\n"
              << "  GTR ln L           : " << analyzer.gtr_log_l << "\n"
              << "  MG94xREV ln L      : " << analyzer.mg94_log_l << "\n"
              << "  Constrained ln L   : " << res_srv.constrained.log_likelihood << "\n"
              << "  Unconstrained ln L : " << res_srv.unconstrained.log_likelihood << "\n"
              << "  LRT                : " << res_srv.lrt << " (p = " << res_srv.p_value << ")\n";

    CHECK(std::isfinite(res_srv.constrained.log_likelihood));
    CHECK(std::isfinite(res_srv.unconstrained.log_likelihood));
    CHECK(res_srv.constrained.log_likelihood >= analyzer.mg94_log_l - 1.0);
    CHECK(res_srv.unconstrained.log_likelihood >= res_srv.constrained.log_likelihood - 1e-6);
    CHECK(res_srv.unconstrained.test_distribution.syn_rates.size() == 3);
    double syn_weight_sum = 0.0;
    for (size_t m = 0; m < 3; ++m) {
        CHECK(std::isfinite(res_srv.unconstrained.test_distribution.syn_rates[m]));
        CHECK(std::isfinite(res_srv.unconstrained.test_distribution.syn_weights[m]));
        CHECK(res_srv.unconstrained.test_distribution.syn_weights[m] > 0.0);
        syn_weight_sum += res_srv.unconstrained.test_distribution.syn_weights[m];
    }
    CHECK(std::abs(syn_weight_sum - 1.0) < 1e-5);
    // Purifying/neutral alignment should not trigger false positive selection at 128 taxa
    CHECK(res_srv.p_value > 0.05);
}
