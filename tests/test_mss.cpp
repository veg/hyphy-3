#include "hyphy/analyses/mss.hpp"
#include "hyphy/analyses/busted.hpp"
#include "hyphy/core/alignment.hpp"
#include "hyphy/core/tree.hpp"
#include "hyphy/core/genetic_code.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#undef NDEBUG
#include <cassert>

using namespace hyphy::core;
using namespace hyphy::analyses;

#ifndef HYPHY3_SOURCE_DIR
#define HYPHY3_SOURCE_DIR "."
#endif

static std::string resolve_path(const std::string& rel) {
    std::ifstream f1(rel);
    if (f1.good()) return rel;
    std::string p2 = std::string(HYPHY3_SOURCE_DIR) + "/" + rel;
    std::ifstream f2(p2);
    if (f2.good()) return p2;
    return rel;
}

int main() {
    std::cout << "=== Running HyPhy 3 MSS (Multi-class Synonymous Substitution) Tests ===\n";

    auto code = GeneticCode::universal();

    // 1. Verify Universal Genetic Code Synonymous Pair Enumeration
    auto pairs = MSSPartition::enumerate_synonymous_pairs(*code);
    std::cout << "[Test 1] Universal genetic code single-step synonymous pairs: " << pairs.size() << "\n";
    assert(pairs.size() == 67 && "Universal code must have exactly 67 unordered single-nucleotide synonymous pairs");

    // 2. Verify Preset Partition Builders
    auto p_synrev = MSSPartition::build("SynREV", *code);
    std::cout << "[Test 2a] SynREV classes: " << p_synrev.class_names.size() << "\n";
    assert(p_synrev.class_names.size() == 18 && "SynREV must have 18 degenerate amino acid classes");

    auto p_codon = MSSPartition::build("SynREVCodon", *code);
    std::cout << "[Test 2b] SynREVCodon classes: " << p_codon.class_names.size() << "\n";
    assert(p_codon.class_names.size() == 67 && "SynREVCodon must have 67 codon-pair classes");

    auto p_rev2 = MSSPartition::build("SynREV2", *code);
    std::cout << "[Test 2c] SynREV2 classes: " << p_rev2.class_names.size() << "\n";
    assert(p_rev2.class_names.size() >= 25 && p_rev2.class_names.size() <= 45);

    auto p_rev2g = MSSPartition::build("SynREV2g", *code);
    std::cout << "[Test 2d] SynREV2g classes: " << p_rev2g.class_names.size() << "\n";
    assert(p_rev2g.class_names.size() >= 19 && p_rev2g.class_names.size() <= 35);

    // 3. Verify Deprecated/Broken HyPhy 2.5 Modes ("Full", "Random") are Rejected
    [[maybe_unused]] bool rejected_full = false, rejected_random = false;
    try {
        MSSPartition::build("Full", *code);
    } catch (const std::exception&) {
        rejected_full = true;
    }
    try {
        MSSPartition::build("Random", *code);
    } catch (const std::exception&) {
        rejected_random = true;
    }
    assert(rejected_full && "Deprecated 'Full' mode with Valine bug must be rejected");
    assert(rejected_random && "Deprecated 'Random' mode must be rejected");
    std::cout << "[Test 3] Deprecated/broken HyPhy 2.5 modes ('Full', 'Random') cleanly rejected.\n";

    // 4. Test Custom TSV Loading & Fitting on Drosophila adh.nex
    std::string adh_path = resolve_path("tests/data/adh.nex");
    auto aln = Alignment::from_nexus(adh_path, code);
    Tree tree = Tree::from_newick(aln.embedded_tree_newick);

    // Write a 3-column Codon TSV matching codons.tsv (SELECTED vs NEUTRAL)
    std::string tmp_codon_tsv = "/tmp/hyphy3_test_codons.tsv";
    {
        std::ofstream out(tmp_codon_tsv);
        for (int c = 0; c < code->num_sense_codons; ++c) {
            std::string cod = code->sense_codons[c];
            char aa = code->sense_codon_aa[c];
            // Mark C/G ending codons of F, L, V, A, G as SELECTED and others as NEUTRAL
            char cp3 = cod[2];
            bool pref = ((aa == 'F' || aa == 'L' || aa == 'V' || aa == 'A' || aa == 'G') && (cp3 == 'C' || cp3 == 'G'));
            out << aa << "\t" << cod << "\t" << (pref ? "SELECTED" : "NEUTRAL") << "\n";
        }
    }

    MSSAnalyzer analyzer(code);
    analyzer.add_dataset("adh.nex", tree, aln);

    MSSSettings tsv_settings;
    tsv_settings.model_type = tmp_codon_tsv;
    tsv_settings.neutral_ref = "NEUTRAL";
    tsv_settings.normalization = MSSNormalization::NeutralRef;
    tsv_settings.normalization_explicit = true;
    tsv_settings.compute_ci = true;
    tsv_settings.compute_lrt = true;

    MSSResult res_tsv = analyzer.run(tsv_settings);
    std::cout << "[Test 4] adh.nex 3-class Codon TSV fit: MG94 lnL = " << res_tsv.baseline_log_l
              << ", MSS lnL = " << res_tsv.mss_log_l
              << ", LRT = " << res_tsv.omnibus_lrt << " (df=" << res_tsv.omnibus_df << ")\n";
    assert(res_tsv.rates.size() == 3);
    assert(res_tsv.mss_log_l >= res_tsv.baseline_log_l - 1e-3 && "MSS fit must improve upon or match standard MG94");
    for (const auto& cr : res_tsv.rates) {
        std::cout << "  Class " << cr.class_name << ": rate = " << cr.rate
                  << " [" << cr.ci_low << ", " << cr.ci_high << "], ref=" << cr.is_reference << "\n";
        if (cr.is_reference) {
            assert(std::abs(cr.rate - 1.0) < 1e-9);
        } else {
            assert(cr.ci_low <= cr.rate + 1e-5 && cr.rate <= cr.ci_high + 1e-5);
        }
    }

    // 5. Test SynREV Preset with Flux Normalization
    MSSSettings synrev_settings;
    synrev_settings.model_type = "SynREV";
    synrev_settings.normalization = MSSNormalization::Flux;
    synrev_settings.normalization_explicit = true;
    synrev_settings.compute_ci = false;
    synrev_settings.compute_lrt = true;

    MSSResult res_synrev = analyzer.run(synrev_settings);
    std::cout << "[Test 5] adh.nex SynREV fit (Flux norm): MSS lnL = " << res_synrev.mss_log_l
              << ", Omnibus LRT = " << res_synrev.omnibus_lrt
              << " (df=" << res_synrev.omnibus_df << ", p=" << res_synrev.omnibus_p_value << ")\n";
    assert(res_synrev.rates.size() == 18);
    assert(res_synrev.mss_log_l > res_synrev.baseline_log_l + 5.0 && "SynREV on adh.nex must substantially improve log-likelihood");

    // Verify flux-weighted mean of class rates equals 1.0
    Scalar weighted_mean = 0.0;
    for (const auto& cr : res_synrev.rates) {
        weighted_mean += cr.equilibrium_flux_share * cr.rate;
    }
    std::cout << "  Flux-weighted mean of SynREV rates: " << weighted_mean << "\n";
    assert(std::abs(weighted_mean - 1.0) < 1e-4 && "Flux-weighted mean must equal 1.0");

    // 6. Test Fast Curvature-Guided DP Partition Selector (--select)
    MSSSettings sel_settings;
    sel_settings.model_type = "SynREV";
    sel_settings.select_partition = true;
    sel_settings.select_classes = 4;
    sel_settings.select_auto_k = true;
    sel_settings.ic_type = "BIC";
    sel_settings.compute_ci = false;
    sel_settings.compute_lrt = false;
    sel_settings.output_tsv = "/tmp/hyphy3_test_selected.tsv";

    MSSResult res_sel = analyzer.run(sel_settings);
    std::cout << "[Test 6] Fast DP Partition Selector on adh.nex: selected K = " << res_sel.selection.num_classes
              << ", lnL = " << res_sel.selection.best_partition_log_l
              << ", BIC = " << res_sel.selection.best_partition_ic << "\n";
    assert(res_sel.selection.enabled);
    assert(res_sel.selection.best_partition_log_l >= res_sel.baseline_log_l - 1e-3);

    // Reload the exported TSV and verify it parses cleanly as a 5-column processor-codon TSV
    auto reloaded_part = MSSPartition::load_from_tsv("/tmp/hyphy3_test_selected.tsv", *code);
    assert(reloaded_part.class_names.size() == res_sel.selection.num_classes);
    std::cout << "  Reloaded exported TSV with " << reloaded_part.class_names.size() << " classes successfully.\n";

    // 7. Test BUSTED-MSS Integration
    BUSTEDAnalyzer busted = BUSTEDAnalyzer::create_and_fit(tree, aln);
    BUSTEDSettings b_settings;
    b_settings.num_rate_classes = 2;
    b_settings.refine_branch_lengths = false;
    b_settings.mss_type = tmp_codon_tsv;
    b_settings.mss_neutral = "NEUTRAL";
    b_settings.mss_normalization = "neutral_ref";

    BUSTEDResult b_res = busted.run(b_settings, false, false);
    std::cout << "[Test 7] BUSTED-MSS on adh.nex: Unconstrained lnL = " << b_res.unconstrained.log_likelihood
              << ", LRT = " << b_res.lrt << ", MSS classes = " << b_res.mss_rates.size() << "\n";
    assert(!b_res.is_failed);
    assert(b_res.mss_rates.size() == 3);

    auto modern_j = busted.to_modern_json(b_res);
    assert(modern_j.contains("mss_rates"));
    assert(modern_j["mss_rates"].size() == 3);

    std::cout << "=== All HyPhy 3 MSS Tests Passed! ===\n";
    return 0;
}
