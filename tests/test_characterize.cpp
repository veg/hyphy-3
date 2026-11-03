#include "hyphy/analyses/characterize.hpp"
#include "hyphy/analyses/joint_align.hpp"
#include "nlohmann/json.hpp"
#undef NDEBUG
#include <cassert>
#include <iostream>
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
        std::cout << "[PASS] Full CharacterizeEngine on adh.nex & JSON v3 Schema Contract Verified\n";
    }

    return 0;
}
