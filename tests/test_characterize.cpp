#include "hyphy/analyses/characterize.hpp"
#include "hyphy/analyses/joint_align.hpp"
#include "nlohmann/json.hpp"
#undef NDEBUG
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace hyphy::core;
using namespace hyphy::analyses;

int main() {
    // 1. Test Stage 0 Triage: Stop-free reverse-complement rescue (T1-03) & quarantined flag reset (T1-04)
    {
        std::string ifitm3_path = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa";
        auto base_seqs = JointAlignEngine::load_unaligned_fasta(ifitm3_path);
        assert(base_seqs.size() >= 6);

        // Strip leading ATG from all sequences (simulating partial CDS like bglobin where 0 sequences start with ATG)
        // and reverse-complement sequence 0 (choosing a stop-free slice if needed)
        std::vector<Sequence> no_atg_cohort(base_seqs.begin(), base_seqs.begin() + 6);
        for (auto& s : no_atg_cohort) {
            if (s.raw_seq.size() > 15) {
                s.raw_seq = s.raw_seq.substr(3); // remove initial ATG
            }
        }
        // Construct a stop-free 60-codon homologous sense sequence and reverse-complement it
        std::string sense_no_rc_stops;
        for (int k = 0; k < 60; ++k) {
            // GCT (Ala, revcomp AGC=Ser), GGT (Gly, revcomp ACC=Thr), CTG (Leu, revcomp CAG=Gln)
            const char* cods[] = {"GCT", "GGT", "CTG", "GTG", "TTC", "AAC"};
            sense_no_rc_stops += cods[(k * 3 + 1) % 6];
        }
        std::vector<Sequence> rc_cohort;
        for (int i = 0; i < 5; ++i) {
            std::string s = sense_no_rc_stops;
            // Introduce a few synonymous/conservative tweaks
            s.replace(3 * (i + 2), 3, "GCC");
            rc_cohort.push_back({"taxon_" + std::to_string(i), s, {}, {}});
        }
        // Taxon 0 is reverse-complemented (and stop-free in both orientations, with no leading ATG)
        rc_cohort[0].raw_seq = JointAlignEngine::reverse_complement_dna(rc_cohort[0].raw_seq);

        // Also add a stop-containing non-homologous outlier (T1-04: e06/e17) that has 0 stops on revcomp frame +1
        std::string stop_fwd_outlier = "ATGTAATAGTGACTA";
        for (int k = 0; k < 50; ++k) stop_fwd_outlier += "CAGCAGCAG"; // revcomp has CTG repeats
        rc_cohort.push_back({"outlier_with_stops", stop_fwd_outlier, {}, {}});

        JointAlignSettings js;
        js.quiet = true;
        auto [admitted, rep] = JointAlignEngine::triage_and_admit_unaligned_sequences(rc_cohort, nullptr, js);
        assert(rep.entries[0].status == "Rescued");
        assert(rep.entries[0].reverse_complemented == true);
        // Verify quarantined outlier has reverse_complemented == false and frame_offset == 0 (T1-04)
        assert(rep.entries[5].status == "Quarantined");
        assert(rep.entries[5].reverse_complemented == false);
        assert(rep.entries[5].frame_offset == 0);
        std::cout << "[PASS] Stage 0 Stop-Free Reverse-Complement Rescue (T1-03) & Quarantined Flag Reset (T1-04)\n";
    }

    // 2. Test Stage 0 Short 5'/3' UTR Trimming in -u Mode (E1) & Genetic Code Suggestion (T1-07)
    {
        std::string ifitm3_path = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa";
        auto seqs = JointAlignEngine::load_unaligned_fasta(ifitm3_path);
        // Add a 9-nt (3-codon) 5' UTR leader before ATG and a native TAA stop + 12-codon 3' UTR tail with stops to seqs[0]
        seqs[0].raw_seq = "ACACACACC" + seqs[0].raw_seq + "TAAGCTGACTGACTAATGACTGA";
        JointAlignSettings js;
        js.quiet = true;
        auto [admitted, rep] = JointAlignEngine::triage_and_admit_unaligned_sequences(seqs, nullptr, js);
        assert(admitted.size() == seqs.size());
        assert(rep.quarantined_sequences == 0);
        assert(rep.entries[0].status != "Quarantined");
        std::cout << "[PASS] Stage 0 Full-Segment 5'/3' UTR Automatic Leader/Tail Trimming (E1)\n";

        // Test T1-07: Vertebrate mtDNA sequences with TGA (Trp) codons rejected under Universal code with --code suggestion
        std::string mt_seq = "ATG";
        for (int k = 0; k < 12; ++k) mt_seq += "TGACTAAACGGTAACTTA";
        std::vector<Sequence> mt_cohort = {
            {"mt1", mt_seq, {}, {}},
            {"mt2", mt_seq, {}, {}},
            {"mt3", mt_seq, {}, {}}
        };
        bool caught_mt = false;
        try {
            JointAlignEngine::triage_and_admit_unaligned_sequences(mt_cohort, GeneticCode::universal(), js);
        } catch (const std::exception& e) {
            std::string msg = e.what();
            if (msg.find("--code Vertebrate-mtDNA") == std::string::npos) {
                std::cerr << "Unexpected T1-07 exception: " << msg << "\n";
            }
            assert(msg.find("--code Vertebrate-mtDNA") != std::string::npos);
            caught_mt = true;
        }
        assert(caught_mt);
        std::cout << "[PASS] Stage 0 Non-Universal Genetic Code Hint (--code Vertebrate-mtDNA) (T1-07)\n";
    }

    // 3. Test -a Mode Guards: Duplicate Taxon (T1-10), Tree-Alignment Mismatch (T1-08), Out-of-Frame Reject (G2-01), 3' Stop Trim (E2)
    {
        Alignment adh = Alignment::load(std::string(HYPHY3_SOURCE_DIR) + "/tests/data/adh.nex");
        Tree adh_tree = Tree::from_newick(adh.embedded_tree_newick);
        JointAlignSettings js;
        js.quiet = true;
        JointAlignEngine j_eng(js);

        // 3a. Duplicate taxon identifier in -a mode (T1-10)
        Alignment dup_aln = adh;
        dup_aln.sequences[1].name = dup_aln.sequences[0].name;
        bool caught_dup = false;
        try {
            j_eng.run(dup_aln);
        } catch (const std::exception& e) {
            assert(std::string(e.what()).find("Duplicate") != std::string::npos);
            caught_dup = true;
        }
        assert(caught_dup);

        // 3b. Tree-alignment taxon mismatch (T1-08 / T5-10)
        Tree bad_tree = adh_tree;
        for (auto& nd : bad_tree.nodes) {
            if (nd.is_leaf) {
                nd.name = "NOT_IN_ALIGNMENT";
                break;
            }
        }
        bool caught_mismatch = false;
        try {
            j_eng.run(adh, bad_tree);
        } catch (const std::exception& e) {
            assert(std::string(e.what()).find("Tree-Alignment Taxon Mismatch") != std::string::npos);
            caught_mismatch = true;
        }
        assert(caught_mismatch);

        // 3c. Out-of-frame alignment in -a mode (G2-01)
        Alignment shifted_aln = adh;
        for (auto& sq : shifted_aln.sequences) {
            sq.raw_seq = "A" + sq.raw_seq.substr(0, sq.raw_seq.size() - 1);
        }
        shifted_aln.finalize(adh.code);
        bool caught_shift = false;
        try {
            j_eng.run(shifted_aln);
        } catch (const std::exception& e) {
            assert(std::string(e.what()).find("Gate 2") != std::string::npos);
            caught_shift = true;
        }
        assert(caught_shift);

        // 3d. Trailing 3' UTR stop column auto-trimming in -a mode (E2)
        Alignment trailing_stop_aln = adh;
        for (auto& sq : trailing_stop_aln.sequences) {
            sq.raw_seq += "TAA";
        }
        trailing_stop_aln.finalize(adh.code);
        assert(trailing_stop_aln.num_codons == adh.num_codons + 1);
        auto trim_res = j_eng.run(trailing_stop_aln, adh_tree);
        assert(trim_res.map_alignment.num_codons == adh.num_codons);
        std::cout << "[PASS] -a Mode Guards: Duplicate Taxon (T1-10), Tree Mismatch (T1-08), Out-of-Frame (G2-01) & 3' Stop Trim (E2)\n";
    }

    // 4. Test Full CharacterizeEngine on adh.nex (-a mode) & JSON Schema Contract (T5-01, T5-03, E4, E6, E16)
    {
        Alignment adh = Alignment::load(std::string(HYPHY3_SOURCE_DIR) + "/tests/data/adh.nex");
        Tree adh_tree = Tree::from_newick(adh.embedded_tree_newick);

        CharacterizeSettings cs;
        cs.quiet = true;
        CharacterizeEngine c_eng(cs);
        auto res = c_eng.run_from_alignment(adh, adh_tree);

        assert(res.input_mode == "aligned");
        assert(res.tree_source == "user_supplied");
        assert(res.joint_align.admission_report.total_input_sequences == adh.num_taxa);
        assert(res.joint_align.admission_report.admitted_sequences == adh.num_taxa);
        assert(res.joint_align.admission_report.median_aa_identity > 0.50);
        assert(res.joint_align.admission_report.median_pairwise_z > 5.0);

        // Verify zero-AA-variation columns have beta_pervasive == 0 and meme_lrt == 0 (E6)
        for (size_t c = 0; c < res.site_atlas.size(); ++c) {
            uint32_t aa_mask = 0;
            size_t n_aa = 0;
            for (size_t t = 0; t < res.joint_align.map_alignment.num_taxa; ++t) {
                int8_t cod = res.joint_align.map_alignment.sequences[t].codons[c];
                if (cod >= 0 && cod < adh.code->num_sense_codons) {
                    char aa = adh.code->sense_codon_aa[cod];
                    int idx = (aa >= 'A' && aa <= 'Z') ? (aa - 'A') : -1;
                    if (idx >= 0 && (aa_mask & (1u << idx)) == 0) {
                        aa_mask |= (1u << idx);
                        n_aa++;
                    }
                }
            }
            if (n_aa <= 1) {
                assert(res.site_atlas[c].beta_pervasive == 0.0);
                assert(res.site_atlas[c].meme_lrt == 0.0);
                assert(res.site_atlas[c].meme_p_value == 1.0);
            }
        }

        // Verify JSON contract (T5-01)
        auto j = nlohmann::json::parse(CharacterizeEngine::to_json(res));
        assert(j.contains("$schema"));
        assert(j["analysis"]["schema_version"] == "3.0.0");
        assert(j["analysis"]["genetic_code"] == "Universal");
        assert(j["analysis"]["input_mode"] == "aligned");
        assert(j["analysis"]["tree_source"] == "user_supplied");
        assert(j["input_admission_audit"]["admitted_sequences"] == adh.num_taxa);
        assert(j["summary_counts"].contains("sites_simes_diversifying"));
        assert(j["summary_counts"].contains("sites_simes_purifying"));
        if (!j["clade_macro_regime"]["triggered"].get<bool>()) {
            assert(j["clade_macro_regime"]["clade_relax_p_value"].get<double>() == 1.0);
            assert(j["clade_macro_regime"]["clade_relax_lrt"].get<double>() == 0.0);
        }
        assert(res.gene_wide.sink_null_log_l_gain < 5.0);
        assert(res.gene_wide.busted_p_value < 0.05);
        std::cout << "[PASS] adh.nex BUSTED-E Error-Sink Keeps Selection Signal (gain="
                  << res.gene_wide.sink_null_log_l_gain << ", p=" << res.gene_wide.busted_p_value << ") (#31)\n";

        // 5. Test User-Specified Branch Tag (-b Foreground) with {Foreground} Newick annotations (#33)
        Tree tagged_tree = adh_tree;
        size_t expected_fg = 0;
        for (auto& nd : tagged_tree.nodes) {
            if (nd.id != tagged_tree.root_id && nd.is_leaf && expected_fg < 4) {
                nd.model_tag = "Foreground";
                expected_fg++;
            }
        }
        CharacterizeSettings cs_fg;
        cs_fg.quiet = true;
        cs_fg.foreground_tag = "Foreground";
        CharacterizeEngine c_eng_fg(cs_fg);
        auto res_fg = c_eng_fg.run_from_alignment(adh, tagged_tree);
        assert(res_fg.clade_macro_regime.triggered == true);
        assert(res_fg.clade_macro_regime.partition_source == "User {Foreground} Tags");
        assert(res_fg.clade_macro_regime.num_test_branches == expected_fg);
        size_t json_fg_count = 0;
        for (const auto& bp : res_fg.branch_profiles) {
            if (bp.macro_regime.rfind("Test", 0) == 0) {
                assert(bp.model_tag == "Foreground");
                json_fg_count++;
            }
        }
        assert(json_fg_count == expected_fg);
        std::cout << "[PASS] User-Specified Branch Tag (-b Foreground) Partition & model_tag Preservation (#33)\n";
    }

    // 6. Test that the BUSTED-E error sink shares the W_tensor scale (#31). A stale W_sink inflated
    // the Null lnL gain to 88.4 on bglobin and 27.6 on adh; HyPhy 2.5.101 BUSTED-E gains 8.0 and 0.5.
    {
        Alignment bglobin = Alignment::load(std::string(HYPHY3_SOURCE_DIR) + "/benchmarks/data/bglobin.nex");
        Tree bglobin_tree = Tree::from_newick(bglobin.embedded_tree_newick);
        CharacterizeSettings cs;
        cs.quiet = true;
        CharacterizeEngine c_eng(cs);
        auto res_bg = c_eng.run_from_alignment(bglobin, bglobin_tree);
        assert(res_bg.gene_wide.sink_null_log_l_gain > 0.0);
        assert(res_bg.gene_wide.sink_null_log_l_gain < 20.0);
        std::cout << "[PASS] bglobin.nex BUSTED-E Error-Sink Null Gain (gain="
                  << res_bg.gene_wide.sink_null_log_l_gain << ", p=" << res_bg.gene_wide.busted_p_value
                  << ", err_w=" << res_bg.gene_wide.error_sink_weight << ") (#31)\n";
    }

    // 7. Test Directional Tree-Alignment Taxon Validation on cd2.fna (#34):
    //    - Missing alignment taxon in supplied tree (9-tip tree missing Rat) -> rejected with Rat named
    //    - Extra tree tips in supplied superset tree (12-tip tree with ExtraSpeciesA/B) -> auto-pruned and accepted
    {
        Alignment cd2 = Alignment::load(std::string(HYPHY3_SOURCE_DIR) + "/benchmarks/data/cd2.fna");
        assert(cd2.num_taxa == 10);

        // 7a. 9-tip tree omitting 'Rat' must be rejected before model fitting
        Tree tree_missing_rat = Tree::from_newick(
            "((((Pig:0.147969,Cow:0.213430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,"
            "((RhMonkey:0.002015,Baboon:0.003108):0.022733,(Human:0.004349,Chimp:0.000799):0.011873):0.101856):0.340802,"
            "Mouse:0.097950);"
        );
        CharacterizeSettings cs;
        cs.quiet = true;
        CharacterizeEngine c_eng(cs);
        bool caught_missing_rat = false;
        try {
            c_eng.run_from_alignment(cd2, tree_missing_rat);
        } catch (const std::exception& e) {
            std::string msg = e.what();
            assert(msg.find("Tree-Alignment Taxon Mismatch") != std::string::npos);
            assert(msg.find("Rat") != std::string::npos);
            caught_missing_rat = true;
        }
        assert(caught_missing_rat);

        // 7b. 12-tip superset/species tree containing all 10 cd2 taxa + 2 extra tips (ExtraSpeciesA, ExtraSpeciesB)
        Tree tree_superset = Tree::from_newick(
            "((((Pig:0.147969,(Cow:0.200000,ExtraSpeciesA:0.050000):0.013430):0.085099,Horse:0.165787,Cat:0.264806):0.058611,"
            "((RhMonkey:0.002015,Baboon:0.003108):0.022733,((Human{Foreground}:0.004349,ExtraSpeciesB{Foreground}:0.003000):0.001000,Chimp:0.000799):0.011873):0.101856):0.340802,"
            "Rat:0.050958,Mouse:0.097950);"
        );
        CharacterizeSettings cs_sup;
        cs_sup.quiet = true;
        cs_sup.foreground_tag = "Foreground";
        CharacterizeEngine c_eng_sup(cs_sup);
        auto res_sup = c_eng_sup.run_from_alignment(cd2, tree_superset);
        size_t leaf_count = 0;
        bool found_human_fg = false;
        for (const auto& bp : res_sup.branch_profiles) {
            assert(bp.branch_name != "ExtraSpeciesA");
            assert(bp.branch_name != "ExtraSpeciesB");
            if (bp.is_leaf) leaf_count++;
            if (bp.branch_name == "Human" && bp.model_tag == "Foreground") {
                found_human_fg = true;
            }
        }
        assert(leaf_count == 10);
        assert(found_human_fg);
        std::cout << "[PASS] Directional Tree-Alignment Validation: Missing Taxon Rejected & Extra Tips Auto-Pruned (#34)\n";
    }

    // 8. Test Stage 2 Coupled SPR on >1,024-Taxon Alignment in ScalableWoodbury Mode (#36)
    {
        const size_t N_large = 1040;
        const size_t C_short = 15;
        Alignment large_aln;
        large_aln.code = GeneticCode::universal();
        large_aln.num_taxa = N_large;
        large_aln.num_codons = C_short;
        large_aln.num_nucleotides = 3 * C_short;
        large_aln.sequences.resize(N_large);
        large_aln.taxon_names.resize(N_large);

        const char* cod_pool[8] = {"GCT", "GCC", "GCA", "GCG", "GGT", "GGC", "GGA", "GGG"};
        std::ostringstream nwk;
        nwk << "(";
        for (size_t i = 0; i < N_large; ++i) {
            std::string name = "T" + std::to_string(i);
            large_aln.taxon_names[i] = name;
            std::string raw;
            for (size_t c = 0; c < C_short; ++c) {
                if (c == 7 && (i % 5 == 0)) {
                    raw += "---"; // internal gap window at codon 7
                } else {
                    size_t idx = ((i >> (c % 10)) ^ (c * 3 + i)) & 7;
                    raw += cod_pool[idx];
                }
            }
            large_aln.sequences[i].name = name;
            large_aln.sequences[i].raw_seq = raw;
            large_aln.sequences[i].codons.resize(C_short);
            large_aln.sequences[i].codon_masks.assign(C_short, 0);
            for (size_t c = 0; c < C_short; ++c) {
                std::string trip = raw.substr(3 * c, 3);
                large_aln.sequences[i].codons[c] = (trip == "---")
                    ? static_cast<int8_t>(-1)
                    : large_aln.code->triplet_to_sense(trip[0], trip[1], trip[2]);
            }
        }
        large_aln.finalize(large_aln.code);

        // Build balanced binary tree on N_large leaves
        auto build_sub_nwk = [&](auto& self, size_t lo, size_t hi) -> std::string {
            if (lo + 1 == hi) {
                return "T" + std::to_string(lo) + ":" + ((lo % 17 == 0) ? "0.025" : "0.005");
            }
            size_t mid = lo + (hi - lo) / 2;
            return "(" + self(self, lo, mid) + "," + self(self, mid, hi) + "):0.005";
        };
        size_t m1 = N_large / 3, m2 = (2 * N_large) / 3;
        std::string large_nwk = "(" + build_sub_nwk(build_sub_nwk, 0, m1) + ","
                              + build_sub_nwk(build_sub_nwk, m1, m2) + ","
                              + build_sub_nwk(build_sub_nwk, m2, N_large) + ");";
        Tree large_tree = Tree::from_newick(large_nwk);

        OhmicSettings ohm_s;
        ohm_s.use_hyphaeon_transformer = false;
        OhmicNetwork check_net = OhmicNetwork::infer(large_aln, ohm_s);
        assert(check_net.active_mode == OhmicScaleMode::ScalableWoodbury);
        assert(check_net.Omega.rows() == 0);
        assert(check_net.taxon_effective_resistance(0, 1039) >= 0.0);

        MG94Parameters mg_p;
        WindowRealignmentSummary w_sum;
        w_sum.window.start_codon = 5;
        w_sum.window.end_codon = 10;
        std::vector<WindowRealignmentSummary> wins = {w_sum};
        JointAlignSettings jas;
        jas.quiet = true;
        jas.max_coupled_spr_candidates = 2;

        auto moves = JointAlignEngine::run_stage2_coupled_topology_search(large_tree, large_aln, mg_p, wins, jas);
        (void)moves;
        std::cout << "[PASS] Stage 2 Coupled SPR on >1,024-Taxon Alignment in ScalableWoodbury Mode (#36)\n";
    }

    return 0;
}
