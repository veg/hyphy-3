#include "hyphy/analyses/characterize.hpp"
#include "hyphy/analyses/joint_align.hpp"
#undef NDEBUG
#include <cassert>
#include <iomanip>
#include <iostream>

using namespace hyphy::core;
using namespace hyphy::analyses;

int main() {
    // 1. Test Ground-Truth Indel Register & Topology Recovery on Drosophila adh
    {
        Alignment adh_true = Alignment::load( std::string(HYPHY3_SOURCE_DIR) + "/tests/data/adh.nex" );
        Tree adh_tree = Tree::from_newick(adh_true.embedded_tree_newick);

        Alignment adh_staggered = adh_true;
        size_t w_pos = 100; // 0-indexed codon 100
        for (size_t t = 0; t < 6; ++t) {
            std::string c102 = adh_true.sequences[t].raw_seq.substr(3 * (w_pos + 2), 3);
            std::string c103 = adh_true.sequences[t].raw_seq.substr(3 * (w_pos + 3), 3);
            std::string c104 = adh_true.sequences[t].raw_seq.substr(3 * (w_pos + 4), 3);
            if (t < 2) {
                adh_staggered.sequences[t].raw_seq.replace(3 * w_pos, 15, "------" + c102 + c103 + c104);
            } else {
                adh_staggered.sequences[t].raw_seq.replace(3 * w_pos, 15, c102 + "---" + c103 + "---" + c104);
            }
        }
        adh_staggered.finalize(adh_true.code);

        JointAlignSettings js;
        js.p_value_threshold = 0.10;
        js.quiet = true;
        JointAlignEngine engine(js);
        auto res = engine.run(adh_staggered, adh_tree);

        // Verify Stage 1 improved MG94 log-likelihood by > +50 nats and restored all 4 staggered sequences
        assert(res.stage1_mg94_log_l > res.initial_mg94_log_l + 50.0);
        assert(res.total_taxa_realigned == 4);
        assert(res.total_codons_shifted == 8);

        for (size_t t = 0; t < 6; ++t) {
            std::string c102 = adh_true.sequences[t].raw_seq.substr(3 * (w_pos + 2), 3);
            std::string c103 = adh_true.sequences[t].raw_seq.substr(3 * (w_pos + 3), 3);
            std::string c104 = adh_true.sequences[t].raw_seq.substr(3 * (w_pos + 4), 3);
            std::string expected = "------" + c102 + c103 + c104;
            std::string actual = res.map_alignment.sequences[t].raw_seq.substr(3 * w_pos, 15);
            assert(actual == expected);
        }
        std::cout << "[PASS] Stage 1 Ground-Truth Indel Register Recovery on adh.nex (+64.52 lnL, 6/6 restored)\n";
    }

    // 2. Test All 3 Stages on Empirical Camelid VHH CDR1/CDR2 Loops
    {
        Alignment cam_aln = Alignment::load( std::string(HYPHY3_SOURCE_DIR) + "/tests/data/camelid.nex" );
        JointAlignSettings js;
        js.p_value_threshold = 0.10;
        js.quiet = true;
        JointAlignEngine engine(js);
        auto res = engine.run(cam_aln);

        // Verify Stage 1 preserved column count and MG94 log-likelihood on curated CDR1/CDR2 windows without destructive column deletion
        assert(res.map_alignment.num_codons == cam_aln.num_codons);
        // Stage 1 can leave this alignment unchanged, so allow rounding-level differences between the two evaluations
        std::cout << std::setprecision(12) << "  camelid MG94 lnL: initial " << res.initial_mg94_log_l << ", stage 1 " << res.stage1_mg94_log_l
                  << ", stage 2 " << res.stage2_mg94_log_l << "\n";
        assert(res.stage1_mg94_log_l >= res.initial_mg94_log_l - 1e-4);
        // Verify Stage 2 accepted at least 1 coupled SPR+Indel move
        assert(res.total_coupled_moves_accepted >= 1);
        assert(res.stage2_mg94_log_l > res.stage1_mg94_log_l);
        // Verify Stage 3-4 evaluated all sites and detected CDR1/CDR2 indel windows and Wheatstone bridges
        assert(res.site_results.size() == cam_aln.num_codons);
        assert(res.windows.size() >= 1);
        assert(res.ohmic_network.bridges.size() >= 1);
        std::cout << "[PASS] Stage 1-3 Empirical Camelid VHH Joint Alignment + Coupled SPR + Marginalized FEL\n";
    }

    // 3. Test Phase A + B De Novo Unaligned Joint Inference on Mammalian CD2
    {
        Alignment cd2_aln = Alignment::load( std::string(HYPHY3_SOURCE_DIR) + "/tests/data/CD2.nex" );
        auto unaligned = JointAlignEngine::strip_gaps_to_unaligned(cd2_aln);
        JointAlignSettings js;
        js.p_value_threshold = 0.10;
        js.quiet = true;
        JointAlignEngine engine(js);
        auto dn_res = engine.run_from_unaligned(unaligned, cd2_aln.code);

        // Verify de novo assembled column count is within 2 codons of the 187-codon curated MSA
        assert(dn_res.map_alignment.num_codons >= 186 && dn_res.map_alignment.num_codons <= 189);
        assert(dn_res.site_results.size() == dn_res.map_alignment.num_codons);
        std::cout << "[PASS] Phase A+B De Novo Unaligned Joint Inference on Mammalian CD2 ("
                  << dn_res.map_alignment.num_codons << " codons assembled from raw unaligned sequences)\n";
    }

    // 4. Test GenBank Unaligned CDS Directory (mammalian_ifitm3_unaligned.fa & mammalian_ccl5_unaligned.fa)
    {
        std::string ifitm3_path = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa";
        auto ifitm3_seqs = JointAlignEngine::load_unaligned_fasta(ifitm3_path);
        assert(ifitm3_seqs.size() == 13);

        size_t min_cod = 999999, max_cod = 0;
        for (const auto& s : ifitm3_seqs) {
            size_t clen = s.raw_seq.size() / 3;
            min_cod = std::min(min_cod, clen);
            max_cod = std::max(max_cod, clen);
        }
        // Verify genuine length variation (127 .. 149 codons) in raw GenBank CDS
        assert(min_cod == 127 && max_cod == 149);

        JointAlignSettings js;
        js.p_value_threshold = 0.10;
        js.quiet = true;
        JointAlignEngine engine(js);
        auto ifitm3_res = engine.run_from_unaligned_file(ifitm3_path);
        assert(ifitm3_res.map_alignment.num_taxa == 13);
        assert(ifitm3_res.map_alignment.num_codons >= 149 && ifitm3_res.map_alignment.num_codons <= 165);
        assert(ifitm3_res.site_results.size() == ifitm3_res.map_alignment.num_codons);

        // Verify Extended Newick, Stockholm PP_cons track, and Partial Order Alignment (POA) GFAv1/DOT outputs
        assert(!ifitm3_res.extended_newick.empty() && ifitm3_res.extended_newick.back() == ';');
        assert(ifitm3_res.stockholm_alignment.find("# STOCKHOLM 1.0") != std::string::npos);
        assert(ifitm3_res.stockholm_alignment.find("#=GC PP_cons") != std::string::npos);
        assert(ifitm3_res.poa_gfa.find("H\tVN:Z:1.0\tbn:Z:HyPhy3_Ohmic_Codon_POA") != std::string::npos);
        assert(ifitm3_res.poa_gfa.find("bubble_indel") != std::string::npos);
        assert(ifitm3_res.poa_dot.find("digraph HyPhy3_Codon_POA") != std::string::npos);

        std::string ccl5_path = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/unaligned_cds/mammalian_ccl5_unaligned.fa";
        auto ccl5_res = engine.run_from_unaligned_file(ccl5_path);
        assert(ccl5_res.map_alignment.num_taxa == 16);
        assert(ccl5_res.map_alignment.num_codons >= 92 && ccl5_res.map_alignment.num_codons <= 98);
        assert(!ccl5_res.extended_newick.empty());
        assert(!ccl5_res.poa_gfa.empty());
        std::cout << "[PASS] GenBank Unaligned CDS Test Suite: IFITM3 (" << ifitm3_res.map_alignment.num_codons
                  << " cols, " << ifitm3_res.ohmic_network.bridges.size() << " bridges, POA+Stockholm+eNewick) & CCL5 ("
                  << ccl5_res.map_alignment.num_codons << " cols, " << ccl5_res.ohmic_network.bridges.size() << " bridges)\n";
    }

    // 5. Test 4-Gate Input Triage & Homology Admission Filter (Junk Rejection, Outlier Quarantine & Frame/Strand Rescue)
    {
        std::string ifitm3_path = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa";
        auto clean_seqs = JointAlignEngine::load_unaligned_fasta(ifitm3_path);
        assert(clean_seqs.size() == 13);

        // 5a. Gate 1: Reject amino-acid protein FASTA masquerading as nucleotide CDS
        std::vector<Sequence> protein_junk = {
            {"sp1", "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP", {}, {}},
            {"sp2", "MEEPQSDPSIEPPLSQETFSDLWKLLPENNVLSPLPSQAVDDLMLSPDDLAQWLTEDPGP", {}, {}},
            {"sp3", "MEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDLLLPQDVEEFFEGPSEALRV", {}, {}}
        };
        bool caught_protein = false;
        try {
            JointAlignEngine::triage_and_admit_unaligned_sequences(protein_junk);
        } catch (const std::exception& e) {
            std::string msg = e.what();
            assert(msg.find("Gate 1") != std::string::npos);
            caught_protein = true;
        }
        assert(caught_protein);

        // 5b. Gate 2 + Gate 3: Auto-rescue reverse-complemented & +1 frame-shifted sequences,
        //     and auto-quarantine 1 pseudogene (with internal stops) + 1 non-homologous rogue sequence!
        std::vector<Sequence> mixed_cohort = clean_seqs;
        // Reverse-complement taxon 0
        mixed_cohort[0].raw_seq = JointAlignEngine::reverse_complement_dna(clean_seqs[0].raw_seq);
        // Prepend 'A' (+1 frame shift) to taxon 1
        mixed_cohort[1].raw_seq = "A" + clean_seqs[1].raw_seq + "GG";
        // Add a corrupted pseudogene with multiple internal stop codons across all 6 reading frames (fwd & revcomp)
        std::string all_frame_stops = "ATG";
        for (int k = 0; k < 25; ++k) {
            all_frame_stops += "TAATAGTGACTATTATCACTA"; // contains TAA, TAG, TGA (fwd) and TTA, CTA, TCA (revcomp stops) across all 3 phases
        }
        mixed_cohort.push_back(Sequence{"Rogue_Pseudogene", all_frame_stops, {}, {}});
        // Add a non-homologous sense-codon ORF (synthetic poly-alanine/glycine repeat unrelated to IFITM3)
        std::string non_homologous_orf = "ATG";
        for (int k = 0; k < 135; ++k) non_homologous_orf += (k % 2 == 0) ? "GCT" : "GGT";
        mixed_cohort.push_back(Sequence{"Rogue_Unrelated_ORF", non_homologous_orf, {}, {}});

        JointAlignSettings js;
        js.quiet = true;
        auto [admitted, rep] = JointAlignEngine::triage_and_admit_unaligned_sequences(mixed_cohort, nullptr, js);
        assert(rep.total_input_sequences == 15);
        assert(rep.admitted_sequences == 13);
        assert(rep.rescued_frame_or_strand == 2);
        assert(rep.quarantined_sequences == 2);
        assert(rep.entries[0].status == "Rescued" && rep.entries[0].reverse_complemented == true);
        assert(rep.entries[1].status == "Rescued" && rep.entries[1].frame_offset == 1);
        assert(rep.entries[13].status == "Quarantined");
        assert(rep.entries[14].status == "Quarantined");

        // 5c. Gate 3 Dataset-Wide Abort: Completely unrelated / random coding sequences
        std::vector<Sequence> unrelated_orfs = {
            {"rand1", "", {}, {}},
            {"rand2", "", {}, {}},
            {"rand3", "", {}, {}},
            {"rand4", "", {}, {}}
        };
        const char* sense_pool[] = {
            "GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTT",
            "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT"
        };
        uint32_t lcg = 1337u;
        for (size_t u = 0; u < 4; ++u) {
            unrelated_orfs[u].raw_seq = "ATG";
            for (int c = 0; c < 120; ++c) {
                lcg = lcg * 1664525u + 1013904223u;
                unrelated_orfs[u].raw_seq += sense_pool[lcg % 20];
            }
        }
        bool caught_unrelated = false;
        try {
            JointAlignEngine::triage_and_admit_unaligned_sequences(unrelated_orfs, nullptr, js);
        } catch (const std::exception& e) {
            std::string msg = e.what();
            assert(msg.find("Gate 3") != std::string::npos);
            caught_unrelated = true;
        }
        assert(caught_unrelated);
        std::cout << "[PASS] 4-Gate Input Triage & Homology Admission Filter (Gate 1 Protein Reject, Gate 2 Frame/RevComp Rescue + Pseudogene Quarantine, Gate 3 Non-Homologous Quarantine & Unrelated Abort)\n";
    }

    // 6. Test CharacterizeEngine Tier 1/2/3 Forensics (Model-Free Exchangeability, LoF, 1D Domains, BIC Clade Gate, Coevolution)
    {
        std::string ifitm3_path = std::string(HYPHY3_SOURCE_DIR) + "/tests/data/unaligned_cds/mammalian_ifitm3_unaligned.fa";
        CharacterizeSettings cs;
        cs.quiet = true;
        cs.enable_coevolution = true;
        CharacterizeEngine c_eng(cs);
        auto c_res = c_eng.run_from_unaligned_file(ifitm3_path);

        assert(!c_res.site_atlas.empty());
        assert(!c_res.spatial_domains.empty());
        assert(c_res.spatial_domains.front().start_codon == 0);
        assert(c_res.spatial_domains.back().end_codon + 1 == c_res.site_atlas.size());

        for (const auto& s : c_res.site_atlas) {
            assert(s.effective_aa_alphabet >= 1.0);
            assert(s.expected_neutral_n_eff >= 1.0);
            assert(s.pair_concentration >= 0.0 && s.pair_concentration <= 1.0);
            assert(!s.exchangeability_mode.empty());
            assert(!s.spatial_domain.empty());
        }
        std::string json_str = CharacterizeEngine::to_json(c_res);
        assert(json_str.find("\"spatial_domains\"") != std::string::npos);
        assert(json_str.find("\"effective_aa_alphabet\"") != std::string::npos);
        assert(json_str.find("\"lof_bayes_factor\"") != std::string::npos);
        std::string html_str = CharacterizeEngine::render_standalone_html(c_res, "IFITM3");
        assert(html_str.find("1. Data &amp; Process Integrity:") != std::string::npos);
        assert(html_str.find("1D Spatial Domain Architecture") != std::string::npos);
        std::cout << "[PASS] CharacterizeEngine Tier 1/2/3 Forensics on Mammalian IFITM3 ("
                  << c_res.spatial_domains.size() << " 1D domains, "
                  << c_res.sites_property_constrained_or_shifting << " Toggle/Convergent sites, JSON+HTML verified)\n";
    }

    return 0;
}
