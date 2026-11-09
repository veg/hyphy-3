#include "hyphy/analyses/characterize.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define HYPHY3_WASM_KEEPALIVE EMSCRIPTEN_KEEPALIVE
#else
#define HYPHY3_WASM_KEEPALIVE
#endif

namespace {
std::string g_last_html_report;
std::string g_last_json_report;
std::string g_last_stage0_json;
std::string g_last_stage1_json;
std::string g_last_stage2_json;
std::string g_last_stage3_json;
std::string g_last_error;

inline std::shared_ptr<const hyphy::core::GeneticCode> safe_genetic_code(const std::string& cname) {
    try {
        auto gc = hyphy::core::GeneticCode::from_name(cname);
        if (gc) return gc;
    } catch (...) {}
    return hyphy::core::GeneticCode::universal();
}

inline double safe_json_num(double v) {
    return std::isfinite(v) ? v : 0.0;
}

inline std::string escape_nexus_single_quotes(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == '\'') out += "''";
        else out += c;
    }
    return out;
}

struct ScopedTempFile {
    std::string path;
    explicit ScopedTempFile(std::string p) : path(std::move(p)) {}
    ~ScopedTempFile() {
        if (!path.empty()) std::remove(path.c_str());
    }
};

struct RawParsedSeq {
    std::string name;
    std::string raw_with_gaps;
    std::string ungapped_nt;
    size_t gap_count = 0;
    size_t ambig_count = 0;
    double gc_overall = 0.0;
    double gc_cp3 = 0.0;
};

std::vector<RawParsedSeq> parse_raw_sequences_preserving_gaps(const std::string& raw_text, bool& out_is_prealigned) {
    std::vector<RawParsedSeq> out;
    out_is_prealigned = false;

    std::istringstream iss(raw_text);
    std::string line;
    std::string cur_name;
    std::string cur_raw;

    auto finalize_record = [&](const std::string& nm, const std::string& rseq) {
        if (nm.empty() && rseq.empty()) return;
        RawParsedSeq rec;
        rec.name = nm.empty() ? ("Taxon_" + std::to_string(out.size() + 1)) : nm;
        rec.raw_with_gaps.reserve(rseq.size());
        rec.ungapped_nt.reserve(rseq.size());
        size_t gc_all = 0, acgt_all = 0;
        size_t gc_3 = 0, acgt_3 = 0;

        for (char c : rseq) {
            if (std::isspace(static_cast<unsigned char>(c))) continue;
            char uc = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            if (uc == 'U') uc = 'T';
            rec.raw_with_gaps += uc;
            if (uc == '-' || uc == '.') {
                rec.gap_count++;
                continue;
            }
            if (uc == 'A' || uc == 'C' || uc == 'G' || uc == 'T') {
                rec.ungapped_nt += uc;
                acgt_all++;
                if (uc == 'G' || uc == 'C') gc_all++;
                if (rec.ungapped_nt.size() % 3 == 0) {
                    acgt_3++;
                    if (uc == 'G' || uc == 'C') gc_3++;
                }
            } else {
                rec.ungapped_nt += 'N';
                rec.ambig_count++;
            }
        }
        rec.gc_overall = (acgt_all > 0) ? (static_cast<double>(gc_all) / static_cast<double>(acgt_all)) : 0.0;
        rec.gc_cp3 = (acgt_3 > 0) ? (static_cast<double>(gc_3) / static_cast<double>(acgt_3)) : 0.0;
        out.push_back(std::move(rec));
    };

    // Detect if NEXUS
    size_t p0 = 0;
    while (p0 < raw_text.size() && std::isspace(static_cast<unsigned char>(raw_text[p0]))) ++p0;
    if (p0 < raw_text.size() && raw_text[p0] == '#') {
        const std::string vpath = "/tmp/walrus_stage0_input.nex";
        ScopedTempFile tmp_guard(vpath);
        {
            std::ofstream vf(vpath, std::ios::binary | std::ios::trunc);
            vf << raw_text;
        }
        auto gcode = hyphy::core::GeneticCode::universal();
        try {
            auto aln = hyphy::core::Alignment::load(vpath, gcode);
            for (const auto& s : aln.sequences) {
                finalize_record(s.name, s.raw_seq);
            }
            out_is_prealigned = true;
        } catch (...) {
            auto seqs = hyphy::analyses::JointAlignEngine::load_unaligned_fasta(vpath, gcode);
            for (const auto& s : seqs) {
                finalize_record(s.name, s.raw_seq);
            }
            out_is_prealigned = false;
        }
        return out;
    }

    while (std::getline(iss, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur_name.empty() || !cur_raw.empty()) {
                finalize_record(cur_name, cur_raw);
                cur_raw.clear();
            }
            std::string hdr = line.substr(1);
            size_t sp = 0;
            while (sp < hdr.size() && std::isspace(static_cast<unsigned char>(hdr[sp]))) ++sp;
            hdr = hdr.substr(sp);
            size_t end_nm = 0;
            while (end_nm < hdr.size() && !std::isspace(static_cast<unsigned char>(hdr[end_nm]))) ++end_nm;
            cur_name = hdr.substr(0, end_nm);
        } else if (line[0] == '(' && out.size() >= 2) {
            // Trailing Newick tree line
            break;
        } else {
            cur_raw += line;
        }
    }
    if (!cur_name.empty() || !cur_raw.empty()) {
        finalize_record(cur_name, cur_raw);
    }

    if (out.size() >= 2) {
        size_t len0 = out[0].raw_with_gaps.size();
        bool all_equal = (len0 > 0 && len0 % 3 == 0);
        for (const auto& r : out) {
            if (r.raw_with_gaps.size() != len0) all_equal = false;
        }
        out_is_prealigned = all_equal;
    }
    return out;
}

bool looks_unaligned_or_needs_stage0(const std::string& raw_text) {
    bool is_prealigned = false;
    auto recs = parse_raw_sequences_preserving_gaps(raw_text, is_prealigned);
    return !is_prealigned;
}

void apply_nexus_translate_to_tree(hyphy::core::Tree& tree, const std::string& raw_text) {
    std::string upper = raw_text;
    for (char& c : upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    size_t tr_pos = upper.find("TRANSLATE");
    if (tr_pos == std::string::npos) return;
    if (tr_pos > 0 && !std::isspace(static_cast<unsigned char>(upper[tr_pos - 1])) && upper[tr_pos - 1] != ';') return;
    size_t semi_pos = raw_text.find(';', tr_pos + 9);
    if (semi_pos == std::string::npos) return;
    std::string block = raw_text.substr(tr_pos + 9, semi_pos - (tr_pos + 9));
    for (char& c : block) {
        if (c == ',') c = ' ';
    }
    std::unordered_map<std::string, std::string> tr_map;
    std::istringstream iss(block);
    std::string k, v;
    while (iss >> k >> v) {
        if (v.size() >= 2 && ((v.front() == '\'' && v.back() == '\'') || (v.front() == '"' && v.back() == '"'))) {
            v = v.substr(1, v.size() - 2);
        }
        tr_map[k] = v;
    }
    if (tr_map.empty()) return;
    for (auto& nd : tree.nodes) {
        if (nd.is_leaf) {
            auto it = tr_map.find(nd.name);
            if (it != tr_map.end()) {
                nd.name = it->second;
            }
        }
    }
    tree.leaf_name_to_id.clear();
    for (const auto& nd : tree.nodes) {
        if (nd.is_leaf) tree.leaf_name_to_id[nd.name] = nd.id;
    }
}

std::string extract_translated_embedded_newick(const std::string& raw_text) {
    std::string nwk;
    std::string upper = raw_text;
    for (char& c : upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    size_t tree_pos = upper.find("TREE ");
    if (tree_pos == std::string::npos) tree_pos = upper.find("TREE\t");
    if (tree_pos != std::string::npos) {
        size_t eq_pos = raw_text.find('=', tree_pos);
        if (eq_pos != std::string::npos) {
            size_t lparen = raw_text.find('(', eq_pos);
            if (lparen != std::string::npos) {
                int paren_depth = 0;
                size_t r_pos = lparen;
                for (; r_pos < raw_text.size(); ++r_pos) {
                    if (raw_text[r_pos] == '(') paren_depth++;
                    else if (raw_text[r_pos] == ')') {
                        paren_depth--;
                        if (paren_depth == 0) break;
                    }
                }
                if (r_pos < raw_text.size()) {
                    size_t semi = raw_text.find(';', r_pos);
                    if (semi != std::string::npos && semi - r_pos <= 64) {
                        nwk = raw_text.substr(lparen, semi - lparen + 1);
                    } else {
                        nwk = raw_text.substr(lparen, r_pos - lparen + 1) + ";";
                    }
                }
            }
        }
    }
    if (nwk.empty()) {
        std::istringstream iss(raw_text);
        std::string line;
        while (std::getline(iss, line)) {
            size_t p = 0;
            while (p < line.size() && std::isspace(static_cast<unsigned char>(line[p]))) ++p;
            if (p < line.size() && line[p] == '(' && line.find(':') != std::string::npos) {
                nwk = line.substr(p);
                break;
            }
        }
    }
    if (nwk.empty()) return "";
    try {
        hyphy::core::Tree t = hyphy::core::Tree::from_newick(nwk);
        apply_nexus_translate_to_tree(t, raw_text);
        return t.to_newick();
    } catch (...) {
        return nwk;
    }
}

struct FrameAnatomy {
    std::string label;
    int frame_offset = 0;
    bool rev_comp = false;
    size_t codon_count = 0;
    size_t internal_stops = 0;
    bool starts_with_atg = false;
    bool has_terminal_stop = false;
    size_t utr_5prime_trimmed = 0;
    size_t utr_3prime_trimmed = 0;
    size_t longest_orf_start = 0;
    size_t longest_orf_len = 0;
    size_t cohort_3mer_hits = 0;
    std::string aa_seq;
    std::string clean_cds;
    std::vector<int8_t> sense_codons;
    nlohmann::json stop_positions = nlohmann::json::array();
    nlohmann::json atg_positions = nlohmann::json::array();
};

FrameAnatomy analyze_single_reading_frame(
    const std::string& ungapped_nt,
    int offset,
    bool rev_comp,
    const hyphy::core::GeneticCode& gcode
) {
    FrameAnatomy fa;
    fa.frame_offset = offset;
    fa.rev_comp = rev_comp;
    fa.label = std::string(rev_comp ? "-" : "+") + std::to_string(offset);

    std::string src = rev_comp ? hyphy::analyses::JointAlignEngine::reverse_complement_dna(ungapped_nt) : ungapped_nt;
    if (static_cast<size_t>(offset) + 3 > src.size()) return fa;

    size_t usable_len = ((src.size() - static_cast<size_t>(offset)) / 3) * 3;
    size_t L = usable_len / 3;
    fa.codon_count = L;

    std::vector<std::string> trips;
    std::vector<int8_t> scs;
    std::vector<bool> is_stop;
    trips.reserve(L);
    scs.reserve(L);
    is_stop.reserve(L);

    size_t cur_orf_start = 0;
    size_t cur_orf_len = 0;

    for (size_t k = 0; k < L; ++k) {
        size_t nt_pos = static_cast<size_t>(offset) + k * 3;
        std::string trip = src.substr(nt_pos, 3);
        int8_t sc = gcode.triplet_to_sense(trip[0], trip[1], trip[2]);
        bool has_n = (trip.find('N') != std::string::npos);
        bool st = (sc < 0 && !has_n);

        trips.push_back(trip);
        scs.push_back(sc);
        is_stop.push_back(st);

        if (trip == "ATG" && fa.atg_positions.size() < 32) {
            fa.atg_positions.push_back(k + 1);
        }
        if (st) {
            fa.aa_seq += '*';
            if (cur_orf_len > fa.longest_orf_len) {
                fa.longest_orf_len = cur_orf_len;
                fa.longest_orf_start = cur_orf_start + 1;
            }
            cur_orf_start = k + 1;
            cur_orf_len = 0;
            if (k + 1 < L && fa.stop_positions.size() < 32) {
                fa.stop_positions.push_back({
                    {"codon", k + 1},
                    {"nt_pos", nt_pos + 1},
                    {"triplet", trip}
                });
            }
        } else {
            char aa = (sc >= 0 && static_cast<size_t>(sc) < gcode.sense_codon_aa.size())
                          ? gcode.sense_codon_aa[static_cast<size_t>(sc)]
                          : 'X';
            fa.aa_seq += aa;
            cur_orf_len++;
        }
    }
    if (cur_orf_len > fa.longest_orf_len) {
        fa.longest_orf_len = cur_orf_len;
        fa.longest_orf_start = cur_orf_start + 1;
    }

    if (L == 0) return fa;
    fa.starts_with_atg = (trips[0] == "ATG");
    fa.has_terminal_stop = is_stop[L - 1];

    // Check 5' and 3' UTR trimming rules identical to Gate 2
    size_t start_idx = 0;
    size_t first_stop = L;
    for (size_t k = 0; k < L; ++k) {
        if (is_stop[k]) { first_stop = k; break; }
    }
    if (first_stop < std::min<size_t>(20, L / 6)) {
        for (size_t k = first_stop + 1; k <= std::min<size_t>(25, L / 5) && k < L; ++k) {
            if (trips[k] == "ATG") {
                size_t next_stop = L;
                for (size_t j = k + 1; j < L; ++j) {
                    if (is_stop[j]) { next_stop = j; break; }
                }
                if (next_stop - k >= std::max<size_t>(15, (L * 4) / 5)) {
                    start_idx = k;
                    first_stop = next_stop;
                    fa.utr_5prime_trimmed = k;
                    break;
                }
            }
        }
    }

    size_t end_idx = L;
    size_t max_3utr_codons = (L >= 50) ? std::min<size_t>(25, std::max<size_t>(1, L / 8)) : 0;
    if (first_stop >= start_idx && first_stop < L) {
        if (first_stop + 1 == L) {
            end_idx = first_stop;
        } else if ((first_stop - start_idx) >= 50 && (L - 1 - first_stop) <= max_3utr_codons) {
            end_idx = first_stop;
            fa.utr_3prime_trimmed = L - 1 - first_stop;
        }
    }

    for (size_t k = start_idx; k < end_idx; ++k) {
        if (is_stop[k]) {
            if (k + 1 == end_idx) break;
            fa.internal_stops++;
        }
        fa.clean_cds += trips[k];
        if (scs[k] >= 0) fa.sense_codons.push_back(scs[k]);
    }
    return fa;
}

// Lightweight banded Needleman-Wunsch amino-acid alignment preview for the Click-Into Inspector
void build_pairwise_aa_preview(
    const std::string& q_aa,
    const std::string& t_aa,
    std::string& out_q_aln,
    std::string& out_mid,
    std::string& out_t_aln
) {
    const size_t n = std::min<size_t>(q_aa.size(), 600);
    const size_t m = std::min<size_t>(t_aa.size(), 600);
    if (n == 0 || m == 0) return;

    const int band = static_cast<int>(std::max<size_t>(32, std::abs(static_cast<int>(n) - static_cast<int>(m)) + 24));
    const int NEG_INF = -1000000;
    std::vector<int> dp((n + 1) * (m + 1), NEG_INF);
    std::vector<uint8_t> bt((n + 1) * (m + 1), 0); // 1: diag, 2: up (del in t), 3: left (ins in t)

    auto idx = [&](size_t i, size_t j) -> size_t { return i * (m + 1) + j; };
    dp[idx(0, 0)] = 0;
    for (size_t i = 1; i <= n && static_cast<int>(i) <= band; ++i) {
        dp[idx(i, 0)] = -2 * static_cast<int>(i);
        bt[idx(i, 0)] = 2;
    }
    for (size_t j = 1; j <= m && static_cast<int>(j) <= band; ++j) {
        dp[idx(0, j)] = -2 * static_cast<int>(j);
        bt[idx(0, j)] = 3;
    }

    for (size_t i = 1; i <= n; ++i) {
        int j_min = std::max(1, static_cast<int>(i) - band);
        int j_max = std::min(static_cast<int>(m), static_cast<int>(i) + band);
        for (int j_int = j_min; j_int <= j_max; ++j_int) {
            size_t j = static_cast<size_t>(j_int);
            char a = q_aa[i - 1];
            char b = t_aa[j - 1];
            int score_match = (a == b && a != 'X' && a != '*') ? 4 : ((a == '*' || b == '*') ? -4 : -1);
            int best = dp[idx(i - 1, j - 1)] + score_match;
            uint8_t dir = 1;
            if (dp[idx(i - 1, j)] - 3 > best) {
                best = dp[idx(i - 1, j)] - 3;
                dir = 2;
            }
            if (dp[idx(i, j - 1)] - 3 > best) {
                best = dp[idx(i, j - 1)] - 3;
                dir = 3;
            }
            dp[idx(i, j)] = best;
            bt[idx(i, j)] = dir;
        }
    }

    size_t i = n, j = m;
    std::string rq, rm, rt;
    while (i > 0 || j > 0) {
        uint8_t d = bt[idx(i, j)];
        if (i > 0 && j > 0 && (d == 1 || d == 0)) {
            char a = q_aa[i - 1];
            char b = t_aa[j - 1];
            rq += a;
            rt += b;
            rm += (a == b && a != 'X') ? '|' : ((a == '*' || b == '*') ? '!' : '.');
            --i;
            --j;
        } else if (i > 0 && (d == 2 || j == 0)) {
            rq += q_aa[i - 1];
            rt += '-';
            rm += ' ';
            --i;
        } else {
            rq += '-';
            rt += t_aa[j - 1];
            rm += ' ';
            --j;
        }
    }
    std::reverse(rq.begin(), rq.end());
    std::reverse(rm.begin(), rm.end());
    std::reverse(rt.begin(), rt.end());
    out_q_aln = std::move(rq);
    out_mid = std::move(rm);
    out_t_aln = std::move(rt);
}

} // namespace

extern "C" {

HYPHY3_WASM_KEEPALIVE
const char* walrus_get_last_json() {
    return g_last_json_report.c_str();
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_get_last_error() {
    return g_last_error.c_str();
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_stage0_run(
    const char* seq_content,
    const char* dataset_label,
    const char* code_name,
    int min_sense_codons,
    double min_homology_z,
    int strip_gaps_realign
) {
    try {
        auto t_start = std::chrono::steady_clock::now();
        g_last_error.clear();
        g_last_stage1_json.clear();
        g_last_stage2_json.clear();
        g_last_stage3_json.clear();
        g_last_json_report.clear();
        g_last_html_report.clear();
        std::string raw_seq = seq_content ? seq_content : "";
        if (raw_seq.empty()) {
            throw std::runtime_error("Empty sequence payload supplied to Stage 0 QA & Admission.");
        }
        std::string lbl = (dataset_label && dataset_label[0] != '\0') ? dataset_label : "DATASET";
        std::string cname = (code_name && code_name[0] != '\0') ? code_name : "Universal";
        auto gcode = safe_genetic_code(cname);

        bool input_is_prealigned = false;
        auto raw_recs = parse_raw_sequences_preserving_gaps(raw_seq, input_is_prealigned);
        if (raw_recs.empty()) {
            throw std::runtime_error("No valid FASTA or NEXUS sequences found in input payload.");
        }

        // Check duplicate names and lexical alphabet gate (amino-acid letters E,F,I,L,P,Q)
        std::unordered_set<std::string> seen_names;
        std::vector<std::string> duplicate_names;
        size_t total_alpha = 0, total_aa_only = 0;
        for (const auto& r : raw_recs) {
            if (!seen_names.insert(r.name).second) {
                duplicate_names.push_back(r.name);
            }
            for (char ch : r.raw_with_gaps) {
                if (std::isalpha(static_cast<unsigned char>(ch))) {
                    total_alpha++;
                    if (ch == 'E' || ch == 'F' || ch == 'I' || ch == 'L' || ch == 'P' || ch == 'Q' || ch == 'Z' || ch == 'J') {
                        total_aa_only++;
                    }
                }
            }
        }
        double aa_only_fraction = (total_alpha > 0) ? (static_cast<double>(total_aa_only) / static_cast<double>(total_alpha)) : 0.0;
        bool is_protein_error = (aa_only_fraction > 0.08);

        // Scan all 13 genetic codes to count how many sequences are stop-free in Frame +0
        const std::vector<std::string> all_codes = {
            "Universal", "Vertebrate-mtDNA", "Invertebrate-mtDNA", "Yeast-mtDNA",
            "Mold-Protozoan-mtDNA", "Ciliate-Nuclear", "Euplotid-Nuclear", "Echinoderm-mtDNA",
            "Ascidian-mtDNA", "Flatworm-mtDNA", "Blepharisma-Nuclear", "Chlorophycean-mtDNA", "Trematode-mtDNA"
        };
        nlohmann::json code_compat = nlohmann::json::array();
        size_t active_code_clean = 0;
        size_t best_alt_clean = 0;
        std::string recommended_code;

        for (const auto& alt_name : all_codes) {
            auto alt_gc = hyphy::core::GeneticCode::from_name(alt_name);
            if (!alt_gc) continue;
            size_t clean_f0 = 0;
            size_t total_stops_f0 = 0;
            for (const auto& r : raw_recs) {
                if (r.ungapped_nt.size() < 45) continue;
                auto fa0 = analyze_single_reading_frame(r.ungapped_nt, 0, false, *alt_gc);
                total_stops_f0 += fa0.internal_stops;
                if (fa0.internal_stops == 0) clean_f0++;
            }
            if (alt_name == cname) active_code_clean = clean_f0;
            if (clean_f0 > best_alt_clean) {
                best_alt_clean = clean_f0;
                recommended_code = alt_name;
            }
            code_compat.push_back({
                {"code", alt_name},
                {"clean_frame0_sequences", clean_f0},
                {"total_internal_stops_frame0", total_stops_f0}
            });
        }

        // Build 6-frame anatomy for every sequence under the selected genetic code
        std::vector<std::array<FrameAnatomy, 6>> all_six_frames(raw_recs.size());
        std::unordered_map<uint32_t, size_t> cohort_aa3mer_counts;

        auto extract_3mers = [&](const std::vector<int8_t>& cods) -> std::unordered_set<uint32_t> {
            std::unordered_set<uint32_t> u;
            if (cods.size() < 3) return u;
            for (size_t k = 0; k + 2 < cods.size(); ++k) {
                if (cods[k] < 0 || cods[k + 1] < 0 || cods[k + 2] < 0) continue;
                uint32_t h = (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k]])) << 16)
                           | (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 1]])) << 8)
                           |  static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cods[k + 2]]));
                u.insert(h);
            }
            return u;
        };

        std::vector<std::unordered_set<uint32_t>> seq_f0_3mers(raw_recs.size());
        for (size_t t = 0; t < raw_recs.size(); ++t) {
            for (int rc = 0; rc <= 1; ++rc) {
                for (int off = 0; off < 3; ++off) {
                    all_six_frames[t][rc * 3 + off] = analyze_single_reading_frame(raw_recs[t].ungapped_nt, off, rc == 1, *gcode);
                }
            }
            if (all_six_frames[t][0].internal_stops == 0 && all_six_frames[t][0].codon_count >= 15) {
                seq_f0_3mers[t] = extract_3mers(all_six_frames[t][0].sense_codons);
                for (uint32_t h : seq_f0_3mers[t]) cohort_aa3mer_counts[h]++;
            }
        }

        for (size_t t = 0; t < raw_recs.size(); ++t) {
            for (int f = 0; f < 6; ++f) {
                auto& fa = all_six_frames[t][f];
                size_t hits = 0;
                for (size_t k = 0; k + 2 < fa.sense_codons.size(); ++k) {
                    uint32_t h = (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[fa.sense_codons[k]])) << 16)
                               | (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[fa.sense_codons[k + 1]])) << 8)
                               |  static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[fa.sense_codons[k + 2]]));
                    auto it = cohort_aa3mer_counts.find(h);
                    if (it != cohort_aa3mer_counts.end()) {
                        size_t self_c = seq_f0_3mers[t].count(h) ? 1 : 0;
                        if (it->second > self_c) hits++;
                    }
                }
                fa.cohort_3mer_hits = hits;
            }
        }

        // Run C++ Stage 0 triage_and_admit_unaligned_sequences
        hyphy::analyses::JointAlignSettings jsettings;
        jsettings.min_sense_codons = static_cast<size_t>(std::max(5, min_sense_codons));
        jsettings.min_homology_z = std::max(0.5, min_homology_z);
        jsettings.quiet = true;

        std::vector<hyphy::core::Sequence> input_seqs;
        input_seqs.reserve(raw_recs.size());
        for (size_t t = 0; t < raw_recs.size(); ++t) {
            std::string nm = raw_recs[t].name;
            // Ensure unique names if duplicates exist so Stage 0 can still diagnose all records
            if ( duplicate_names.size() > 0 ) {
                nm = raw_recs[t].name + "#" + std::to_string(t + 1);
            }
            input_seqs.push_back(hyphy::core::Sequence{nm, raw_recs[t].ungapped_nt, {}, {}});
        }

        hyphy::analyses::InputAdmissionReport adm_report;
        std::string gate_abort_warning;
        if (!is_protein_error && input_seqs.size() >= 3) {
            try {
                jsettings.enforce_admission_control = true;
                auto [adm_seqs, rep] = hyphy::analyses::JointAlignEngine::triage_and_admit_unaligned_sequences(input_seqs, gcode, jsettings);
                adm_report = std::move(rep);
            } catch (const std::exception& e_strict) {
                gate_abort_warning = e_strict.what();
                try {
                    jsettings.enforce_admission_control = false;
                    auto [adm_seqs2, rep2] = hyphy::analyses::JointAlignEngine::triage_and_admit_unaligned_sequences(input_seqs, gcode, jsettings);
                    adm_report = std::move(rep2);
                } catch (...) {
                    adm_report.total_input_sequences = raw_recs.size();
                    adm_report.entries.resize(raw_recs.size());
                }
            }
        } else {
            adm_report.total_input_sequences = raw_recs.size();
            adm_report.entries.resize(raw_recs.size());
            if (is_protein_error) {
                gate_abort_warning = "Gate 1 (Lexical Alphabet Gate): Input contains " +
                    std::to_string(static_cast<int>( std::round(aa_only_fraction * 100.0) )) +
                    "% amino-acid-only letters (E, F, I, L, P, Q). HyPhy codon models require nucleotide coding sequences (CDS).";
            } else if (input_seqs.size() < 3) {
                gate_abort_warning = "Gate 1 (Input Size Gate): Phylogenetic selection analysis requires at least 3 homologous coding sequences.";
            }
        }

        // Compute median codon length of valid sequences
        std::vector<size_t> valid_lengths;
        for (size_t t = 0; t < raw_recs.size(); ++t) {
            size_t clen = (t < adm_report.entries.size() && adm_report.entries[t].codon_length > 0)
                              ? adm_report.entries[t].codon_length
                              : all_six_frames[t][0].codon_count;
            if (clen >= jsettings.min_sense_codons) valid_lengths.push_back(clen);
        }
        std::sort(valid_lengths.begin(), valid_lengths.end());
        size_t median_codons = valid_lengths.empty() ? 100 : valid_lengths[valid_lengths.size() / 2];
        size_t min_allowed_codons = std::max<size_t>(jsettings.min_sense_codons, static_cast<size_t>(std::floor(0.30 * static_cast<double>(median_codons))));
        size_t max_allowed_codons = static_cast<size_t>(std::ceil(2.80 * static_cast<double>(median_codons)));

        // Identify reference admitted sequences for nearest-neighbor homology preview
        std::vector<size_t> ref_indices;
        for (size_t t = 0; t < raw_recs.size() && t < adm_report.entries.size(); ++t) {
            if (adm_report.entries[t].status != "Quarantined" && all_six_frames[t][0].codon_count >= jsettings.min_sense_codons) {
                ref_indices.push_back(t);
            }
        }
        if (ref_indices.empty()) {
            for (size_t t = 0; t < std::min<size_t>(raw_recs.size(), 10); ++t) ref_indices.push_back(t);
        }

        nlohmann::json taxa_json = nlohmann::json::array();
        size_t count_admitted = 0, count_rescued = 0, count_quarantined = 0;

        for (size_t t = 0; t < raw_recs.size(); ++t) {
            const auto& r = raw_recs[t];
            hyphy::analyses::TaxonAdmissionEntry ent = (t < adm_report.entries.size()) ? adm_report.entries[t] : hyphy::analyses::TaxonAdmissionEntry{};
            ent.name = r.name;

            int chosen_f_idx = (ent.reverse_complemented ? 3 : 0) + std::clamp(ent.frame_offset, 0, 2);
            const auto& chosen_fa = all_six_frames[t][chosen_f_idx];
            const auto& f0_fa = all_six_frames[t][0];

            // Find best-frame candidate across all 6 frames even for quarantined sequences
            int best_alt_f = 0;
            for (int f = 1; f < 6; ++f) {
                if (all_six_frames[t][f].internal_stops < all_six_frames[t][best_alt_f].internal_stops ||
                    (all_six_frames[t][f].internal_stops == all_six_frames[t][best_alt_f].internal_stops &&
                     all_six_frames[t][f].cohort_3mer_hits > all_six_frames[t][best_alt_f].cohort_3mer_hits)) {
                    best_alt_f = f;
                }
            }

            // Find nearest neighbor among reference sequences by 3-mer overlap or amino acid identity
            size_t best_nn_idx = (ref_indices.empty()) ? 0 : ref_indices[0];
            if (best_nn_idx == t && ref_indices.size() > 1) best_nn_idx = ref_indices[1];
            const auto& q_aa_for_nn = (ent.status == "Quarantined") ? all_six_frames[t][best_alt_f].aa_seq : chosen_fa.aa_seq;
            std::unordered_set<uint32_t> q_3mers = extract_3mers(
                (ent.status == "Quarantined") ? all_six_frames[t][best_alt_f].sense_codons : chosen_fa.sense_codons
            );

            size_t best_shared_3mers = 0;
            for (size_t cand_idx : ref_indices) {
                if (cand_idx == t) continue;
                int cf = (adm_report.entries[cand_idx].reverse_complemented ? 3 : 0) + std::clamp(adm_report.entries[cand_idx].frame_offset, 0, 2);
                const auto& cand_cods = all_six_frames[cand_idx][cf].sense_codons;
                size_t shared = 0;
                for (size_t k = 0; k + 2 < cand_cods.size(); ++k) {
                    uint32_t h = (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cand_cods[k]])) << 16)
                               | (static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cand_cods[k + 1]])) << 8)
                               |  static_cast<uint32_t>(static_cast<uint8_t>(gcode->sense_codon_aa[cand_cods[k + 2]]));
                    if (q_3mers.count(h)) shared++;
                }
                if (shared > best_shared_3mers) {
                    best_shared_3mers = shared;
                    best_nn_idx = cand_idx;
                }
            }

            int nn_f = (best_nn_idx < adm_report.entries.size())
                           ? ((adm_report.entries[best_nn_idx].reverse_complemented ? 3 : 0) + std::clamp(adm_report.entries[best_nn_idx].frame_offset, 0, 2))
                           : 0;
            const std::string& nn_aa = all_six_frames[best_nn_idx][nn_f].aa_seq;
            std::string q_aln, mid_aln, t_aln;
            build_pairwise_aa_preview(q_aa_for_nn, nn_aa, q_aln, mid_aln, t_aln);

            // If quarantined before Gate 3, estimate AA identity and Z-score from the pairwise preview against best_nn_idx
            double disp_z = ent.max_homology_z;
            double disp_id = ent.max_aa_identity;
            if (disp_id <= 0.0 && !mid_aln.empty()) {
                size_t matches = 0, aligned_cols = 0;
                for (size_t k = 0; k < mid_aln.size(); ++k) {
                    if (q_aln[k] != '-' && t_aln[k] != '-') {
                        aligned_cols++;
                        if (mid_aln[k] == '|') matches++;
                    }
                }
                if (aligned_cols > 0) {
                    disp_id = static_cast<double>(matches) / static_cast<double>(aligned_cols);
                    double L_eff = static_cast<double>(aligned_cols);
                    disp_z = (disp_id - 0.08) / std::max(0.01, std::sqrt(0.08 * 0.92 / L_eff));
                }
            }

            // Build structured diagnostic verdict for the Click-Into Inspector
            std::string failed_gate = "None";
            std::string diag_title = "Admitted In-Frame Coding Sequence";
            std::string diag_explanation = "Sequence is in reading frame +0 on the forward strand with zero internal stop codons and strong evolutionary homology to the cohort.";
            std::string remediation_hint = "Ready for Stage 1 codon alignment.";

            double ambig_frac = (r.ungapped_nt.size() > 0) ? (static_cast<double>(r.ambig_count) / static_cast<double>(r.ungapped_nt.size())) : 0.0;
            size_t eff_codons = (ent.codon_length > 0) ? ent.codon_length : chosen_fa.codon_count;

            if (is_protein_error) {
                ent.status = "Quarantined";
                failed_gate = "Gate 1 (Lexical Alphabet)";
                diag_title = "Protein Amino-Acid Sequence Detected";
                diag_explanation = "This sequence contains amino-acid-specific characters (E, F, I, L, P, Q) rather than nucleotide triplets (A, C, G, T).";
                remediation_hint = "Upload the corresponding nucleotide coding sequence (CDS) FASTA file.";
            } else if (ent.status == "Rescued") {
                count_admitted++;
                count_rescued++;
                failed_gate = "Gate 2 (Auto-Rescued ORF)";
                diag_title = "Automatically Rescued Reading Frame / Strand";
                std::ostringstream exp;
                exp << "Raw sequence had " << f0_fa.internal_stops << " internal stop codon(s) or low 3-mer homology in forward frame +0, "
                    << "but scanning all 6 reading frames identified a clean, stop-free open reading frame on the "
                    << (ent.reverse_complemented ? "reverse-complement (-)" : "forward (+)")
                    << " strand at frame offset +" << ent.frame_offset
                    << " sharing " << chosen_fa.cohort_3mer_hits << " amino-acid 3-mers with the cohort.";
                if (chosen_fa.utr_5prime_trimmed > 0 || chosen_fa.utr_3prime_trimmed > 0) {
                    exp << " Trimmed " << chosen_fa.utr_5prime_trimmed << " 5'-UTR leader codon(s) and "
                        << chosen_fa.utr_3prime_trimmed << " 3'-UTR tail codon(s).";
                }
                diag_explanation = exp.str();
                remediation_hint = "Automatically corrected and admitted. Inspect the 6-frame barcode below to verify the rescued reading frame.";
            } else if (ent.status == "Quarantined") {
                count_quarantined++;
                if (r.ungapped_nt.size() < jsettings.min_sense_codons * 3) {
                    failed_gate = "Gate 1 (Minimum Length)";
                    diag_title = "Quarantined: Sequence Fragment Too Short";
                    diag_explanation = "Sequence contains only " + std::to_string(r.ungapped_nt.size() / 3) +
                        " ungapped codons (" + std::to_string(r.ungapped_nt.size()) +
                        " nt), which is below the minimum admission threshold of " +
                        std::to_string(jsettings.min_sense_codons) + " sense codons.";
                    remediation_hint = "Lower the 'Min Codons' threshold slider or click 'Force Admit' if this partial exon fragment is intentional.";
                } else if (ambig_frac > 0.20) {
                    failed_gate = "Gate 1 (Ambiguity Filter)";
                    diag_title = "Quarantined: Excessive Ambiguous Nucleotides (N/?)";
                    std::ostringstream exp;
                    exp << "Sequence contains " << r.ambig_count << " ambiguous nucleotides ("
                        << std::fixed << std::setprecision(1) << (ambig_frac * 100.0)
                        << "% of " << r.ungapped_nt.size() << " nt), exceeding the 20.0% maximum ambiguity limit.";
                    diag_explanation = exp.str();
                    remediation_hint = "High-N sequences degrade branch-length and site-level selection inference. Click 'Force Admit' only if missing regions are localized.";
                } else if (ent.reason.find("length outlier") != std::string::npos) {
                    failed_gate = "Gate 3 (Length Outlier)";
                    diag_title = "Quarantined: Extreme Sequence Length Outlier";
                    diag_explanation = "Sequence length (" + std::to_string(eff_codons) +
                        " codons) lies outside the permissible cohort length window [" +
                        std::to_string(min_allowed_codons) + ", " + std::to_string(max_allowed_codons) +
                        " codons] (0.30x to 2.80x cohort median of " + std::to_string(median_codons) + " codons).";
                    remediation_hint = "Often caused by un-spliced multi-exon concatamers or severe truncations. Click 'Force Admit' if this species has a genuine large indel.";
                } else if (ent.reason.find("Non-homologous") != std::string::npos) {
                    failed_gate = "Gate 3 (Homology Z-Score)";
                    diag_title = "Quarantined: Insufficient Evolutionary Homology";
                    std::ostringstream exp;
                    exp << "Sequence has zero internal stop codons in frame " << chosen_fa.label
                        << ", but its maximum null-calibrated pairwise homology Z-score against the cohort is Z = "
                        << std::fixed << std::setprecision(2) << disp_z << " (max AA identity = "
                        << std::setprecision(1) << (disp_id * 100.0) << "% vs '" << raw_recs[best_nn_idx].name
                        << "'), falling below the admission threshold Z >= " << std::setprecision(1) << jsettings.min_homology_z << ".";
                    diag_explanation = exp.str();
                    remediation_hint = "Inspect the Nearest-Neighbor Alignment below: if this is a fast-evolving divergent ortholog rather than a paralog/contaminant, lower Min Z or click 'Force Admit'.";
                } else {
                    failed_gate = "Gate 2 (Premature Stop Codons)";
                    // Check if a frameshift split the ORF across two frames!
                    size_t f0_longest = all_six_frames[t][0].longest_orf_len;
                    size_t f1_longest = all_six_frames[t][1].longest_orf_len;
                    size_t f2_longest = all_six_frames[t][2].longest_orf_len;
                    size_t second_best = std::max(f1_longest, f2_longest);
                    bool likely_frameshift = (f0_longest >= eff_codons / 3 && second_best >= eff_codons / 3 &&
                                              (f0_longest + second_best) >= (eff_codons * 4) / 5);
                    if (likely_frameshift) {
                        diag_title = "Quarantined: Mid-Sequence Reading-Frame Shift Detected";
                        std::ostringstream exp;
                        exp << "Frame +0 contains " << f0_fa.internal_stops << " premature stop codon(s), and no single reading frame is stop-free across the entire sequence. "
                            << "However, the 6-frame scan reveals a clean " << f0_longest << "-codon ORF in Frame +0 followed by a clean "
                            << second_best << "-codon ORF in Frame +" << (f1_longest >= f2_longest ? 1 : 2)
                            << " — the classic signature of a 1-bp or 2-bp sequencing/assembly indel shifting the reading frame mid-gene.";
                        diag_explanation = exp.str();
                        remediation_hint = "Stage 1's codon-aware frameshift transducer can repair 1-bp/2-bp frameshifts if you click 'Force Admit', or you can keep it quarantined as a pseudogene.";
                    } else {
                        diag_title = "Quarantined: Internal Premature Stop Codons in All 6 Frames";
                        std::ostringstream exp;
                        exp << "Forward Frame +0 contains " << f0_fa.internal_stops << " internal stop codon(s), and scanning all 6 reading frames (+0, +1, +2, -0, -1, -2) found no stop-free ORF "
                            << "(best candidate frame " << all_six_frames[t][best_alt_f].label << " still has "
                            << all_six_frames[t][best_alt_f].internal_stops << " internal stops).";
                        if (!recommended_code.empty() && recommended_code != cname && best_alt_clean > active_code_clean) {
                            exp << " Note: Switching genetic code to '" << recommended_code << "' resolves stop codons across the cohort!";
                        }
                        diag_explanation = exp.str();
                        remediation_hint = (!recommended_code.empty() && recommended_code != cname && best_alt_clean > active_code_clean)
                            ? ("Switch Genetic Code to '" + recommended_code + "' in the top bar.")
                            : "Likely a pseudogene, retained intron, or sequencing error. Inspect the red stop-codon ticks in the 6-frame barcode below.";
                    }
                }
            } else {
                count_admitted++;
                if (chosen_fa.utr_5prime_trimmed > 0 || chosen_fa.utr_3prime_trimmed > 0) {
                    diag_title = "Admitted In-Frame CDS (UTR Flanks Trimmed)";
                    diag_explanation = "Sequence is in Frame +0 with zero internal stops after automatically trimming " +
                        std::to_string(chosen_fa.utr_5prime_trimmed) + " 5'-UTR leader codon(s) and " +
                        std::to_string(chosen_fa.utr_3prime_trimmed) + " 3'-UTR trailing codon(s).";
                }
            }

            nlohmann::json six_frames_arr = nlohmann::json::array();
            for (int f = 0; f < 6; ++f) {
                const auto& fa = all_six_frames[t][f];
                six_frames_arr.push_back({
                    {"frame_label", fa.label},
                    {"frame_offset", fa.frame_offset},
                    {"reverse_complemented", fa.rev_comp},
                    {"is_selected", (f == chosen_f_idx && ent.status != "Quarantined")},
                    {"is_best_fallback", (f == best_alt_f)},
                    {"codon_count", fa.codon_count},
                    {"internal_stops", fa.internal_stops},
                    {"starts_with_atg", fa.starts_with_atg},
                    {"has_terminal_stop", fa.has_terminal_stop},
                    {"utr_5prime_trimmed", fa.utr_5prime_trimmed},
                    {"utr_3prime_trimmed", fa.utr_3prime_trimmed},
                    {"longest_orf_start", fa.longest_orf_start},
                    {"longest_orf_len", fa.longest_orf_len},
                    {"cohort_3mer_hits", fa.cohort_3mer_hits},
                    {"stop_positions", fa.stop_positions},
                    {"atg_positions", fa.atg_positions},
                    {"aa_seq", fa.aa_seq},
                    {"clean_cds", fa.clean_cds}
                });
            }

            taxa_json.push_back({
                {"index", t + 1},
                {"name", r.name},
                {"status", ent.status},
                {"failed_gate", failed_gate},
                {"reason", ent.reason},
                {"diagnostic_title", diag_title},
                {"diagnostic_explanation", diag_explanation},
                {"remediation_hint", remediation_hint},
                {"frame_offset", ent.frame_offset},
                {"reverse_complemented", ent.reverse_complemented},
                {"selected_frame_label", chosen_fa.label},
                {"raw_nt_length", r.raw_with_gaps.size()},
                {"ungapped_nt_length", r.ungapped_nt.size()},
                {"gap_count", r.gap_count},
                {"ambig_count", r.ambig_count},
                {"ambig_fraction", ambig_frac},
                {"gc_overall", r.gc_overall},
                {"gc_cp3", r.gc_cp3},
                {"codon_length", eff_codons},
                {"internal_stops_frame0", f0_fa.internal_stops},
                {"internal_stops_selected", chosen_fa.internal_stops},
                {"utr_5prime_trimmed", chosen_fa.utr_5prime_trimmed},
                {"utr_3prime_trimmed", chosen_fa.utr_3prime_trimmed},
                {"max_homology_z", safe_json_num(disp_z)},
                {"max_aa_identity", safe_json_num(disp_id)},
                {"clean_cds", (ent.status == "Quarantined") ? all_six_frames[t][best_alt_f].clean_cds : chosen_fa.clean_cds},
                {"raw_with_gaps", r.raw_with_gaps},
                {"translated_aa", (ent.status == "Quarantined") ? all_six_frames[t][best_alt_f].aa_seq : chosen_fa.aa_seq},
                {"six_frames", six_frames_arr},
                {"nearest_neighbor", {
                    {"target_name", (best_nn_idx < raw_recs.size()) ? raw_recs[best_nn_idx].name : ""},
                    {"query_aa_aligned", q_aln},
                    {"midline", mid_aln},
                    {"target_aa_aligned", t_aln}
                }}
            });
        }

        auto t_end = std::chrono::steady_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        nlohmann::json out_j = {
            {"stage", 0},
            {"stage_name", "Stage 0: Input Triage, Reading-Frame Rescue & Homology Admission"},
            {"dataset_label", lbl},
            {"genetic_code", cname},
            {"elapsed_ms", elapsed_ms},
            {"input_is_prealigned", input_is_prealigned},
            {"alignment_mode", (input_is_prealigned && strip_gaps_realign == 0) ? "preserve_alignment" : "de_novo_unaligned"},
            {"thresholds", {
                {"min_sense_codons", jsettings.min_sense_codons},
                {"min_homology_z", safe_json_num(jsettings.min_homology_z)},
                {"median_codons", median_codons},
                {"min_allowed_codons", min_allowed_codons},
                {"max_allowed_codons", max_allowed_codons}
            }},
            {"summary", {
                {"total_input_sequences", raw_recs.size()},
                {"admitted_sequences", count_admitted},
                {"rescued_frame_or_strand", count_rescued},
                {"quarantined_sequences", count_quarantined},
                {"median_pairwise_z", safe_json_num(adm_report.median_pairwise_z)},
                {"median_aa_identity", safe_json_num(adm_report.median_aa_identity)},
                {"aa_only_letter_fraction", safe_json_num(aa_only_fraction)},
                {"duplicate_names", duplicate_names},
                {"gate_abort_warning", gate_abort_warning},
                {"recommended_code", (best_alt_clean > active_code_clean) ? recommended_code : ""},
                {"recommended_code_clean_count", best_alt_clean},
                {"active_code_clean_count", active_code_clean},
                {"embedded_tree_newick", extract_translated_embedded_newick(raw_seq)}
            }},
            {"code_compatibility", code_compat},
            {"taxa", taxa_json}
        };

        g_last_stage0_json = out_j.dump();
        return g_last_stage0_json.c_str();
    } catch (const std::exception& e) {
        g_last_error = e.what();
        return "";
    } catch (...) {
        g_last_error = "Unhandled internal exception in WASM engine.";
        return "";
    }
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_stage1_run(
    const char* seq_content,
    const char* dataset_label,
    const char* code_name,
    int refine_indel_windows,
    int force_de_novo_realign
) {
    try {
        auto t_start = std::chrono::steady_clock::now();
        g_last_error.clear();
        std::string raw_seq = seq_content ? seq_content : "";
        if (raw_seq.empty()) {
            throw std::runtime_error("Empty sequence payload supplied to Stage 1 Codon Alignment.");
        }
        std::string lbl = (dataset_label && dataset_label[0] != '\0') ? dataset_label : "STAGE1_COHORT";
        std::string cname = (code_name && code_name[0] != '\0') ? code_name : "Universal";
        auto gcode = safe_genetic_code(cname);

        bool input_is_prealigned = false;
        auto raw_recs = parse_raw_sequences_preserving_gaps(raw_seq, input_is_prealigned);
        if (raw_recs.size() < 3) {
            throw std::runtime_error("Stage 1 requires at least 3 coding sequences.");
        }

        // 1. Detect and repair any 1-bp / 2-bp mid-gene reading-frame shifts in sequences that have internal stops in Frame +0
        //    (e.g., sequences force-admitted by the user from Stage 0).
        struct FrameshiftRepairRecord {
            size_t taxon_index = 0;
            std::string taxon_name;
            size_t split_codon = 0; // 1-based
            size_t split_nt = 0;    // 1-based
            int shift_bp = 1;       // 1 or 2 bp excised
            std::string excised_nt;
            size_t stops_before = 0;
            std::string before_aa;
            std::string after_aa;
            std::string headline;
            std::string explanation;
        };
        std::vector<FrameshiftRepairRecord> fs_repairs;
        std::unordered_map<size_t, size_t> taxon_to_fs_idx;

        std::vector<std::string> clean_ungapped(raw_recs.size());
        for (size_t t = 0; t < raw_recs.size(); ++t) {
            auto fa0 = analyze_single_reading_frame(raw_recs[t].ungapped_nt, 0, false, *gcode);
            if (fa0.internal_stops == 0) {
                clean_ungapped[t] = fa0.clean_cds.empty() ? raw_recs[t].ungapped_nt : fa0.clean_cds;
            } else {
                // Check if Frame +1 or Frame +2 has a clean suffix ORF that rescues the mid-gene frameshift
                auto fa1 = analyze_single_reading_frame(raw_recs[t].ungapped_nt, 1, false, *gcode);
                auto fa2 = analyze_single_reading_frame(raw_recs[t].ungapped_nt, 2, false, *gcode);
                int target_f = (fa1.longest_orf_len >= fa2.longest_orf_len) ? 1 : 2;
                const auto& fa_alt = (target_f == 1) ? fa1 : fa2;

                // Split point is just before the first stop codon in Frame +0 (or end of longest prefix ORF in +0)
                size_t first_stop_c = fa0.stop_positions.empty()
                    ? (fa0.codon_count / 2)
                    : (fa0.stop_positions[0].value("codon", fa0.codon_count / 2) - 1);
                if (first_stop_c >= 8 && fa_alt.longest_orf_len >= 8) {
                    // Search slightly upstream of first_stop_c for the split point where splicing (prefix of +0) + (suffix of +target_f)
                    // yields zero internal stops!
                    size_t best_split_c = first_stop_c;
                    std::string best_spliced_nt;
                    std::string best_excised;
                    for (size_t cand_c = (first_stop_c > 12 ? first_stop_c - 12 : 4); cand_c <= first_stop_c; ++cand_c) {
                        size_t nt_cut = 3 * cand_c;
                        if (nt_cut + static_cast<size_t>(target_f) >= raw_recs[t].ungapped_nt.size()) break;
                        std::string prefix = raw_recs[t].ungapped_nt.substr(0, nt_cut);
                        std::string excised = raw_recs[t].ungapped_nt.substr(nt_cut, static_cast<size_t>(target_f));
                        std::string suffix = raw_recs[t].ungapped_nt.substr(nt_cut + static_cast<size_t>(target_f));
                        size_t rem3 = suffix.size() % 3;
                        if (rem3 != 0) suffix.resize(suffix.size() - rem3);
                        std::string candidate = prefix + suffix;
                        auto check_fa = analyze_single_reading_frame(candidate, 0, false, *gcode);
                        if (check_fa.internal_stops == 0) {
                            best_split_c = cand_c;
                            best_spliced_nt = check_fa.clean_cds.empty() ? candidate : check_fa.clean_cds;
                            best_excised = excised;
                            break;
                        }
                    }
                    if (!best_spliced_nt.empty()) {
                        clean_ungapped[t] = best_spliced_nt;
                        auto repaired_fa = analyze_single_reading_frame(best_spliced_nt, 0, false, *gcode);
                        FrameshiftRepairRecord rec;
                        rec.taxon_index = t + 1;
                        rec.taxon_name = raw_recs[t].name;
                        rec.split_codon = best_split_c + 1;
                        rec.split_nt = 3 * best_split_c + 1;
                        rec.shift_bp = target_f;
                        rec.excised_nt = best_excised;
                        rec.stops_before = fa0.internal_stops;
                        rec.before_aa = fa0.aa_seq;
                        rec.after_aa = repaired_fa.aa_seq;
                        rec.headline = "Repaired " + std::to_string(target_f) + "-bp frameshift indel (+" + best_excised + ") near codon " + std::to_string(rec.split_codon);
                        std::ostringstream exp;
                        exp << "Sequence '" << raw_recs[t].name << "' contained " << fa0.internal_stops
                            << " premature stop codon(s) in Frame +0 due to a " << target_f << "-bp insertion (+"
                            << best_excised << ") near codon " << rec.split_codon << " (nt " << rec.split_nt
                            << ") that shifted translation into Frame +" << target_f
                            << ". Stage 1 excised the out-of-frame insertion and restored a clean "
                            << (best_spliced_nt.size() / 3) << "-codon open reading frame with 0 internal stops.";
                        rec.explanation = exp.str();
                        taxon_to_fs_idx[t] = fs_repairs.size();
                        fs_repairs.push_back(std::move(rec));
                    } else {
                        clean_ungapped[t] = raw_recs[t].ungapped_nt;
                    }
                } else {
                    clean_ungapped[t] = raw_recs[t].ungapped_nt;
                }
            }
            size_t r3 = clean_ungapped[t].size() % 3;
            if (r3 != 0) clean_ungapped[t].resize(clean_ungapped[t].size() - r3);
        }

        // 2. Build codon alignment & phylogenetic tree (either preserving pre-aligned input or running de novo Pair-HMM + POA)
        hyphy::analyses::JointAlignSettings jsettings;
        jsettings.quiet = true;
        jsettings.enforce_admission_control = false;

        hyphy::core::Alignment map_aln;
        hyphy::core::Tree tree;
        bool used_de_novo = (!input_is_prealigned || force_de_novo_realign != 0 || !fs_repairs.empty());

        if (!used_de_novo) {
            const std::string vpath = "/tmp/walrus_stage1_in.fa";
            ScopedTempFile tmp_guard(vpath);
            {
                std::ofstream vf(vpath, std::ios::binary | std::ios::trunc);
                for (const auto& r : raw_recs) {
                    vf << ">" << r.name << "\n" << r.raw_with_gaps << "\n";
                }
            }
            map_aln = hyphy::core::Alignment::load(vpath, gcode);
            hyphy::analyses::TreeEngineSettings tes;
            tes.quiet = true;
            hyphy::analyses::TreeEngine te(tes);
            tree = te.run(map_aln).optimized_tree;
        } else {
            std::vector<hyphy::core::Sequence> unaligned_seqs;
            unaligned_seqs.reserve(raw_recs.size());
            for (size_t t = 0; t < raw_recs.size(); ++t) {
                unaligned_seqs.push_back(hyphy::core::Sequence{raw_recs[t].name, clean_ungapped[t], {}, {}});
            }
            auto [dn_aln, dn_tree] = hyphy::analyses::JointAlignEngine::build_de_novo_codon_alignment_and_tree(
                unaligned_seqs, gcode, jsettings, nullptr);
            map_aln = std::move(dn_aln);
            tree = std::move(dn_tree);
        }

        hyphy::core::LikelihoodEngine::reconcile_tree_leaf_names(tree, map_aln);

        // 3. Fit initial MG94 model & optionally run Stage 1 Inside-Outside EM Window Realignment
        double initial_mg94_log_l = 0.0;
        double stage1_mg94_log_l = 0.0;
        double global_omega = 0.25;
        std::vector<hyphy::core::Scalar> site_homology_posteriors(map_aln.num_codons, 1.0);
        std::vector<hyphy::analyses::WindowRealignmentSummary> windows;
        std::vector<std::vector<hyphy::analyses::TransducerTaxonResult>> window_taxon_results;
        size_t total_taxa_realigned = 0;
        size_t total_codons_shifted = 0;

        {
            hyphy::analyses::GTRFitter gtr_fitter(tree, map_aln);
            auto gtr_fit = gtr_fitter.fit();
            hyphy::core::MG94Parameters mg_params;
            mg_params.theta_AC = gtr_fit.params.theta_AC;
            mg_params.theta_AT = gtr_fit.params.theta_AT;
            mg_params.theta_CG = gtr_fit.params.theta_CG;
            mg_params.theta_CT = gtr_fit.params.theta_CT;
            mg_params.theta_GT = gtr_fit.params.theta_GT;

            hyphy::analyses::MG94Fitter mg_fitter(gtr_fit.tree, map_aln);
            auto mg_fit = mg_fitter.fit_omega_and_scale(1.0, mg_params);
            mg_params.alpha = 1.0;
            mg_params.beta = mg_fit.x_opt(0);
            global_omega = mg_params.beta;
            tree = mg_fit.tree;
            initial_mg94_log_l = mg_fit.log_likelihood;
            stage1_mg94_log_l = initial_mg94_log_l;

            if (refine_indel_windows != 0) {
                auto backup_aln = map_aln;
                windows = hyphy::analyses::JointAlignEngine::run_stage1_window_realignment(
                    map_aln,
                    tree,
                    mg_params,
                    jsettings,
                    site_homology_posteriors,
                    window_taxon_results
                );
                for (const auto& w : windows) {
                    total_taxa_realigned += w.num_taxa_realigned;
                    total_codons_shifted += w.num_codons_shifted;
                }
                if (total_codons_shifted > 0) {
                    hyphy::analyses::MG94Fitter mg_stage1(tree, map_aln);
                    auto fit1 = mg_stage1.fit_omega_and_scale(mg_params.beta, mg_params);
                    if (!used_de_novo && fit1.log_likelihood < initial_mg94_log_l) {
                        map_aln = std::move(backup_aln);
                        site_homology_posteriors.assign(map_aln.num_codons, 1.0);
                        total_taxa_realigned = 0;
                        total_codons_shifted = 0;
                        stage1_mg94_log_l = initial_mg94_log_l;
                    } else {
                        global_omega = fit1.x_opt(0);
                        tree = fit1.tree;
                        stage1_mg94_log_l = fit1.log_likelihood;
                    }
                }
            } else {
                auto raw_wins = hyphy::analyses::JointAlignEngine::detect_indel_windows(
                    map_aln, jsettings.window_flank_codons, jsettings.max_window_codons);
                for (const auto& rw : raw_wins) {
                    hyphy::analyses::WindowRealignmentSummary ws;
                    ws.window = rw;
                    windows.push_back(ws);
                }
            }
        }

        // 4. Compute per-column consensus, gap fraction, AA conservation, and per-taxon alignment telemetry
        size_t N = map_aln.num_taxa;
        size_t C = map_aln.num_codons;

        std::string consensus_aa(C, '-');
        std::string consensus_nt(3 * C, '-');
        std::vector<double> col_gap_fraction(C, 0.0);
        std::vector<double> col_aa_conservation(C, 1.0);
        std::vector<double> col_posterior(C, 1.0);

        size_t gap_free_cols = 0;
        double sum_conservation = 0.0;
        double sum_posterior = 0.0;

        std::vector<std::string> taxon_aa(N, std::string(C, '-'));

        for (size_t c = 0; c < C; ++c) {
            size_t gap_cnt = 0;
            std::unordered_map<char, size_t> aa_counts;
            std::unordered_map<std::string, size_t> trip_counts;

            for (size_t t = 0; t < N; ++t) {
                int8_t sc = (c < map_aln.sequences[t].codons.size()) ? map_aln.sequences[t].codons[c] : -1;
                std::string trip = (3 * c + 2 < map_aln.sequences[t].raw_seq.size())
                                       ? map_aln.sequences[t].raw_seq.substr(3 * c, 3)
                                       : "---";
                if (sc >= 0 && sc < gcode->num_sense_codons) {
                    char aa = gcode->sense_codon_aa[sc];
                    taxon_aa[t][c] = aa;
                    aa_counts[aa]++;
                    trip_counts[trip]++;
                } else if (trip == "---") {
                    taxon_aa[t][c] = '-';
                    gap_cnt++;
                } else {
                    taxon_aa[t][c] = 'X';
                }
            }

            if (gap_cnt == 0) gap_free_cols++;
            col_gap_fraction[c] = (N > 0) ? (static_cast<double>(gap_cnt) / static_cast<double>(N)) : 0.0;

            char best_aa = '-';
            size_t best_aa_cnt = 0;
            for (const auto& kv : aa_counts) {
                if (kv.second > best_aa_cnt) {
                    best_aa_cnt = kv.second;
                    best_aa = kv.first;
                }
            }
            consensus_aa[c] = best_aa;

            std::string best_trip = "---";
            size_t best_trip_cnt = 0;
            for (const auto& kv : trip_counts) {
                if (kv.second > best_trip_cnt) {
                    best_trip_cnt = kv.second;
                    best_trip = kv.first;
                }
            }
            if (best_trip.size() == 3) {
                consensus_nt[3 * c] = best_trip[0];
                consensus_nt[3 * c + 1] = best_trip[1];
                consensus_nt[3 * c + 2] = best_trip[2];
            }

            size_t non_gap = N - gap_cnt;
            col_aa_conservation[c] = (non_gap > 0 && best_aa_cnt > 0)
                                         ? (static_cast<double>(best_aa_cnt) / static_cast<double>(non_gap))
                                         : 0.0;
            col_posterior[c] = (c < site_homology_posteriors.size()) ? site_homology_posteriors[c] : 1.0;
            sum_conservation += col_aa_conservation[c];
            sum_posterior += col_posterior[c];
        }

        // Build per-taxon telemetry and FASTA exports
        nlohmann::json taxa_json = nlohmann::json::array();
        nlohmann::json events_json = nlohmann::json::array();
        std::ostringstream aligned_cds_fasta, aligned_aa_fasta;

        // First add any Frameshift Repair events to events_json
        for (const auto& fs : fs_repairs) {
            events_json.push_back({
                {"event_type", "frameshift_repair"},
                {"badge", "Frameshift Repaired"},
                {"badge_style", "rescued"},
                {"title", fs.taxon_name},
                {"subtitle", "Codon " + std::to_string(fs.split_codon) + " • +" + std::to_string(fs.shift_bp) + " bp insertion excised"},
                {"summary", fs.headline},
                {"explanation", fs.explanation},
                {"start_codon", (fs.split_codon > 10 ? fs.split_codon - 10 : 1)},
                {"end_codon", std::min(C, fs.split_codon + 15)},
                {"focal_codon", fs.split_codon},
                {"taxon_index", fs.taxon_index},
                {"taxon_name", fs.taxon_name},
                {"before_aa", fs.before_aa},
                {"after_aa", fs.after_aa}
            });
        }

        // Next add Indel Window events to events_json
        for (size_t w_idx = 0; w_idx < windows.size(); ++w_idx) {
            const auto& w = windows[w_idx];
            size_t sc1 = w.window.start_codon + 1;
            size_t ec1 = w.window.end_codon;
            size_t width = (ec1 >= sc1) ? (ec1 - sc1 + 1) : 1;

            std::vector<std::string> gapped_names;
            double win_min_post = 1.0;
            double win_sum_post = 0.0;
            size_t win_cols = 0;
            for (size_t c = w.window.start_codon; c < w.window.end_codon && c < C; ++c) {
                win_min_post = std::min(win_min_post, col_posterior[c]);
                win_sum_post += col_posterior[c];
                win_cols++;
            }
            double win_mean_post = (win_cols > 0) ? (win_sum_post / static_cast<double>(win_cols)) : 1.0;

            for (size_t t = 0; t < N; ++t) {
                bool has_gap_in_win = false;
                for (size_t c = w.window.start_codon; c < w.window.end_codon && c < C; ++c) {
                    if (taxon_aa[t][c] == '-') {
                        has_gap_in_win = true;
                        break;
                    }
                }
                if (has_gap_in_win) gapped_names.push_back(map_aln.sequences[t].name);
            }

            // Extract competing POA / Pair-HMM candidate gap placements for taxa in this window
            nlohmann::json poa_paths_json = nlohmann::json::array();
            if (w_idx < window_taxon_results.size()) {
                const auto& tx_vec = window_taxon_results[w_idx];
                for (size_t t = 0; t < tx_vec.size() && t < N; ++t) {
                    const auto& tr = tx_vec[t];
                    nlohmann::json cands_j = nlohmann::json::array();
                    for (size_t k = 0; k < std::min<size_t>(4, tr.top_candidates.size()); ++k) {
                        const auto& cd = tr.top_candidates[k];
                        if (k > 0 && cd.posterior_weight < 0.005) continue;
                        std::string cand_aa;
                        for (int8_t sc : cd.codons) {
                            if (sc < 0) cand_aa += '-';
                            else if (sc < gcode->num_sense_codons) cand_aa += gcode->sense_codon_aa[sc];
                            else cand_aa += 'X';
                        }
                        cands_j.push_back({
                            {"rank", k + 1},
                            {"is_map", k == 0},
                            {"posterior_weight", cd.posterior_weight},
                            {"log_score", cd.log_score},
                            {"aa_str", cand_aa}
                        });
                    }
                    if (cands_j.size() >= 2 || (tr.changed && !cands_j.empty())) {
                        poa_paths_json.push_back({
                            {"taxon_index", t + 1},
                            {"taxon_name", map_aln.sequences[t].name},
                            {"map_posterior_prob", tr.map_posterior_prob},
                            {"changed", tr.changed},
                            {"codons_shifted", tr.codons_shifted},
                            {"candidates", cands_j}
                        });
                    }
                }
            }

            std::string badge_str = (w.num_codons_shifted > 0)
                ? "Polished Indel Region"
                : (win_min_post < 0.95 ? "POA Uncertainty Region" : "Indel Region");
            std::string badge_style = (w.num_codons_shifted > 0 || win_min_post < 0.95) ? "rescued" : "admitted";
            std::ostringstream sum_ss, exp_ss;
            sum_ss << "Codons " << sc1 << "–" << ec1 << " (" << width << " codons): "
                   << gapped_names.size() << " of " << N << " sequences gapped • Min homology P = "
                   << std::fixed << std::setprecision(0) << (win_min_post * 100.0) << "%";
            if (w.num_codons_shifted > 0) {
                sum_ss << " • " << w.num_taxa_realigned << " refined (ΔlnL = +"
                       << std::fixed << std::setprecision(1) << std::max(0.0, w.delta_log_l) << ")";
            } else if (!poa_paths_json.empty()) {
                sum_ss << " • " << poa_paths_json.size() << " sequence(s) with alternative POA gap paths";
            }

            exp_ss << "Local indel window spanning codons " << sc1 << " to " << ec1 << " (" << width
                   << " codons; mean column homology confidence = " << std::fixed << std::setprecision(1)
                   << (win_mean_post * 100.0) << "%, minimum = " << (win_min_post * 100.0) << "%). ";
            if (w.num_codons_shifted > 0) {
                exp_ss << "Inside-Outside Pair-HMM realignment shifted " << w.num_codons_shifted
                       << " codon(s) across " << w.num_taxa_realigned << " sequence(s), improving phylogenetic codon likelihood by ΔlnL = +"
                       << std::fixed << std::setprecision(2) << w.delta_log_l << ". ";
            }
            if (!poa_paths_json.empty()) {
                exp_ss << poa_paths_json.size() << " sequence(s) exhibit competing gap placements in the Forward-Backward transducer ensemble; Stage 3 marginalizes site likelihoods over these alternative paths.";
            } else {
                exp_ss << "Progressive POA gap placement has >99.5% posterior concentration under the phylogenetic codon transducer.";
            }

            events_json.push_back({
                {"event_type", "indel_window"},
                {"badge", badge_str},
                {"badge_style", badge_style},
                {"title", "Indel Region // Codons " + std::to_string(sc1) + "–" + std::to_string(ec1)},
                {"subtitle", std::to_string(width) + " codons • " + std::to_string(gapped_names.size()) + " gapped taxa • Min P = " + std::to_string(static_cast<int>( std::round(win_min_post * 100.0) )) + "%"},
                {"summary", sum_ss.str()},
                {"explanation", exp_ss.str()},
                {"start_codon", sc1},
                {"end_codon", ec1},
                {"focal_codon", sc1},
                {"affected_taxa", gapped_names},
                {"num_taxa_realigned", w.num_taxa_realigned},
                {"num_codons_shifted", w.num_codons_shifted},
                {"delta_log_l", w.delta_log_l},
                {"min_column_posterior", win_min_post},
                {"mean_column_posterior", win_mean_post},
                {"poa_uncertainty_paths", poa_paths_json}
            });
        }

        for (size_t t = 0; t < N; ++t) {
            const auto& sq = map_aln.sequences[t];
            aligned_cds_fasta << ">" << sq.name << "\n" << sq.raw_seq << "\n";
            aligned_aa_fasta << ">" << sq.name << "\n" << taxon_aa[t] << "\n";

            size_t gap_codons = 0;
            size_t indel_runs = 0;
            size_t match_cons = 0;
            size_t non_gap = 0;
            bool in_gap = false;

            for (size_t c = 0; c < C; ++c) {
                if (taxon_aa[t][c] == '-') {
                    gap_codons++;
                    if (!in_gap) {
                        indel_runs++;
                        in_gap = true;
                    }
                } else {
                    in_gap = false;
                    non_gap++;
                    if (taxon_aa[t][c] == consensus_aa[c]) match_cons++;
                }
            }

            double gap_frac = (C > 0) ? (static_cast<double>(gap_codons) / static_cast<double>(C)) : 0.0;
            double id_cons = (non_gap > 0) ? (static_cast<double>(match_cons) / static_cast<double>(non_gap)) : 0.0;
            bool has_fs = (taxon_to_fs_idx.count(t) > 0);

            std::string badge = "Aligned";
            std::string badge_style = "admitted";
            std::string summary_line = "Clean codon alignment (" + std::to_string(non_gap) + "/" + std::to_string(C) + " codons)";

            if (has_fs) {
                const auto& fs = fs_repairs[taxon_to_fs_idx[t]];
                badge = "Frameshift Repaired";
                badge_style = "rescued";
                summary_line = fs.headline;
            } else if (gap_codons > 0) {
                badge = std::to_string(gap_codons) + " Gap Codon" + (gap_codons > 1 ? "s" : "");
                badge_style = (gap_frac > 0.15) ? "rescued" : "admitted";
                summary_line = std::to_string(indel_runs) + " indel block(s) (" + std::to_string(gap_codons) + " gap codons)";
            }

            taxa_json.push_back({
                {"index", t + 1},
                {"name", sq.name},
                {"non_gap_codons", non_gap},
                {"gap_codons", gap_codons},
                {"gap_fraction", gap_frac},
                {"indel_runs", indel_runs},
                {"aa_identity_to_consensus", id_cons},
                {"frameshift_repaired", has_fs},
                {"badge", badge},
                {"badge_style", badge_style},
                {"summary_line", summary_line},
                {"aligned_nt", sq.raw_seq},
                {"aligned_aa", taxon_aa[t]}
            });
        }

        auto t_end = std::chrono::steady_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        nlohmann::json out_j = {
            {"stage", 1},
            {"stage_name", "Stage 1: Codon Alignment, Frameshift Repair & Indel Refinement"},
            {"dataset_label", lbl},
            {"genetic_code", cname},
            {"elapsed_ms", elapsed_ms},
            {"alignment_mode", used_de_novo ? "de_novo_poa" : "preserved_input"},
            {"summary", {
                {"num_taxa", N},
                {"num_codons", C},
                {"num_nt", 3 * C},
                {"gap_free_columns", gap_free_cols},
                {"gap_free_fraction", (C > 0 ? static_cast<double>(gap_free_cols) / static_cast<double>(C) : 1.0)},
                {"mean_aa_conservation", (C > 0 ? sum_conservation / static_cast<double>(C) : 1.0)},
                {"mean_column_posterior", (C > 0 ? sum_posterior / static_cast<double>(C) : 1.0)},
                {"indel_windows_count", windows.size()},
                {"frameshifts_repaired_count", fs_repairs.size()},
                {"taxa_realigned_in_polish", total_taxa_realigned},
                {"codons_shifted_in_polish", total_codons_shifted},
                {"initial_mg94_log_l", safe_json_num(initial_mg94_log_l)},
                {"stage1_mg94_log_l", safe_json_num(stage1_mg94_log_l)},
                {"delta_mg94_log_l", safe_json_num(stage1_mg94_log_l - initial_mg94_log_l)},
                {"global_omega", safe_json_num(global_omega)},
                {"tree_newick", tree.to_newick()}
            }},
            {"consensus", {
                {"consensus_aa", consensus_aa},
                {"consensus_nt", consensus_nt},
                {"column_gap_fraction", col_gap_fraction},
                {"column_aa_conservation", col_aa_conservation},
                {"column_posterior", col_posterior}
            }},
            {"exports", {
                {"aligned_cds_fasta", aligned_cds_fasta.str()},
                {"aligned_aa_fasta", aligned_aa_fasta.str()}
            }},
            {"events", events_json},
            {"taxa", taxa_json}
        };

        g_last_stage1_json = out_j.dump();
        return g_last_stage1_json.c_str();
    } catch (const std::exception& e) {
        g_last_error = e.what();
        return "";
    } catch (...) {
        g_last_error = "Unhandled internal exception in WASM engine.";
        return "";
    }
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_stage2_run(
    const char* seq_content,
    const char* tree_newick_in,
    const char* dataset_label,
    const char* code_name,
    int enable_coupled_spr,
    int enable_recomb_scan
) {
    try {
        auto t_start = std::chrono::steady_clock::now();
        g_last_error.clear();
        std::string raw_seq = seq_content ? seq_content : "";
        if (raw_seq.empty()) {
            throw std::runtime_error("Empty sequence payload supplied to Stage 2 Tree & Recombination.");
        }
        std::string nwk_in = tree_newick_in ? tree_newick_in : "";
        std::string lbl = (dataset_label && dataset_label[0] != '\0') ? dataset_label : "STAGE2_COHORT";
        std::string cname = (code_name && code_name[0] != '\0') ? code_name : "Universal";
        auto gcode = safe_genetic_code(cname);

        // 1. Load codon alignment (or build de novo if unaligned input was passed directly)
        hyphy::core::Alignment map_aln;
        hyphy::core::Tree work_tree;
        std::string tree_source = "inferred_nj_spr";

        size_t p0 = 0;
        while (p0 < raw_seq.size() && std::isspace(static_cast<unsigned char>(raw_seq[p0]))) ++p0;
        const bool is_nexus = (p0 < raw_seq.size() && raw_seq[p0] == '#');
        const std::string vpath = is_nexus ? "/tmp/walrus_stage2_in.nex" : "/tmp/walrus_stage2_in.fa";
        ScopedTempFile tmp_guard(vpath);
        {
            std::ofstream vf(vpath, std::ios::binary | std::ios::trunc);
            vf << raw_seq;
        }

        bool loaded_aligned = false;
        if (!looks_unaligned_or_needs_stage0(raw_seq)) {
            try {
                map_aln = hyphy::core::Alignment::load(vpath, gcode);
                loaded_aligned = true;
            } catch (...) {
                loaded_aligned = false;
            }
        }

        if (!loaded_aligned) {
            auto unaligned_seqs = hyphy::analyses::JointAlignEngine::load_unaligned_fasta(vpath, gcode);
            hyphy::analyses::JointAlignSettings jas;
            jas.quiet = true;
            auto [dn_aln, dn_tree] = hyphy::analyses::JointAlignEngine::build_de_novo_codon_alignment_and_tree(
                unaligned_seqs, gcode, jas, nullptr);
            map_aln = std::move(dn_aln);
            work_tree = std::move(dn_tree);
            tree_source = "de_novo_poa_tree";
        } else {
            if (!nwk_in.empty()) {
                try {
                    work_tree = hyphy::core::Tree::from_newick(nwk_in);
                    tree_source = "stage1_or_user_tree";
                } catch (...) {
                    nwk_in.clear();
                }
            }
            if (nwk_in.empty()) {
                std::string emb_nwk = extract_translated_embedded_newick(raw_seq);
                if (emb_nwk.empty()) emb_nwk = map_aln.embedded_tree_newick;
                if (!emb_nwk.empty()) {
                    try {
                        work_tree = hyphy::core::Tree::from_newick(emb_nwk);
                        tree_source = "embedded_in_alignment";
                    } catch (...) {
                        hyphy::analyses::TreeEngineSettings tes;
                        tes.quiet = true;
                        hyphy::analyses::TreeEngine te(tes);
                        work_tree = te.run(map_aln).optimized_tree;
                    }
                } else {
                    hyphy::analyses::TreeEngineSettings tes;
                    tes.quiet = true;
                    hyphy::analyses::TreeEngine te(tes);
                    work_tree = te.run(map_aln).optimized_tree;
                }
            }
        }

        apply_nexus_translate_to_tree(work_tree, raw_seq);
        hyphy::core::LikelihoodEngine::reconcile_tree_leaf_names(work_tree, map_aln);
        const size_t N = map_aln.num_taxa;
        const size_t C = map_aln.num_codons;
        const int S = static_cast<int>(gcode->num_sense_codons);

        // Validate leaf correspondence; prune extra tree leaves or rebuild if alignment taxa are missing from tree
        {
            auto check_l2t = hyphy::core::LikelihoodEngine::build_leaf_to_taxon_map(work_tree, map_aln);
            std::vector<bool> tax_seen(N, false);
            std::vector<std::string> matched_leaf_names;
            size_t total_tree_leaves = 0;
            for (const auto& nd : work_tree.nodes) {
                if (!nd.is_leaf) continue;
                total_tree_leaves++;
                if (nd.id >= 0 && static_cast<size_t>(nd.id) < check_l2t.size()) {
                    size_t tid = check_l2t[nd.id];
                    if (tid < N && !tax_seen[tid]) {
                        tax_seen[tid] = true;
                        matched_leaf_names.push_back(nd.name);
                    }
                }
            }
            if (matched_leaf_names.size() < N) {
                hyphy::analyses::TreeEngineSettings tes;
                tes.quiet = true;
                hyphy::analyses::TreeEngine te(tes);
                work_tree = te.run(map_aln).optimized_tree;
                tree_source = "inferred_nj_spr";
            } else if (total_tree_leaves > matched_leaf_names.size()) {
                work_tree = work_tree.prune_to_leaves(matched_leaf_names);
            }
            hyphy::core::LikelihoodEngine::reconcile_tree_leaf_names(work_tree, map_aln);
        }

        // 2. Fit GTR + MG94 branch lengths and optionally run Stage 2 Coupled SPR Topology Search
        double initial_mg94_log_l = 0.0;
        double final_mg94_log_l = 0.0;
        double global_omega = 0.25;
        std::vector<hyphy::analyses::CoupledMoveEvent> coupled_moves;

        hyphy::analyses::GTRFitter gtr_fitter(work_tree, map_aln);
        auto gtr_fit = gtr_fitter.fit();
        hyphy::core::MG94Parameters mg_params;
        mg_params.theta_AC = gtr_fit.params.theta_AC;
        mg_params.theta_AT = gtr_fit.params.theta_AT;
        mg_params.theta_CG = gtr_fit.params.theta_CG;
        mg_params.theta_CT = gtr_fit.params.theta_CT;
        mg_params.theta_GT = gtr_fit.params.theta_GT;

        hyphy::analyses::MG94Fitter mg_fitter(gtr_fit.tree, map_aln);
        auto mg_fit = mg_fitter.fit_omega_and_scale(1.0, mg_params);
        mg_params.alpha = 1.0;
        mg_params.beta = mg_fit.x_opt(0);
        global_omega = mg_params.beta;
        work_tree = mg_fit.tree;
        initial_mg94_log_l = mg_fit.log_likelihood;
        final_mg94_log_l = initial_mg94_log_l;

        hyphy::analyses::JointAlignSettings jsettings;
        jsettings.quiet = true;
        jsettings.enable_coupled_topology_search = (enable_coupled_spr != 0);

        if (enable_coupled_spr != 0 && N >= 5) {
            auto raw_wins = hyphy::analyses::JointAlignEngine::detect_indel_windows(
                map_aln, jsettings.window_flank_codons, jsettings.max_window_codons);
            std::vector<hyphy::analyses::WindowRealignmentSummary> win_summaries;
            for (const auto& rw : raw_wins) {
                hyphy::analyses::WindowRealignmentSummary ws;
                ws.window = rw;
                win_summaries.push_back(ws);
            }
            coupled_moves = hyphy::analyses::JointAlignEngine::run_stage2_coupled_topology_search(
                work_tree, map_aln, mg_params, win_summaries, jsettings);
            if (!coupled_moves.empty()) {
                hyphy::analyses::MG94Fitter mg_stage2(work_tree, map_aln);
                auto fit2 = mg_stage2.fit_omega_and_scale(mg_params.beta, mg_params);
                mg_params.beta = fit2.x_opt(0);
                global_omega = mg_params.beta;
                work_tree = fit2.tree;
                final_mg94_log_l = fit2.log_likelihood;
            }
        }

        // 3. Compositional Drift & Auto-Rooting Audit (dFEL)
        auto comp_test = hyphy::analyses::DirectionalFELAnalyzer::run_composition_test(map_aln, 0.05);
        bool directional_drift_active = false;
        std::string chosen_root_branch = "Stationary (Midpoint / Unrooted)";
        double chosen_root_weight = 1.0;

        if (comp_test.is_heterogeneous) {
            hyphy::analyses::DirectionalFELSettings dfs;
            dfs.p_value_threshold = 0.10;
            dfs.root_strategy = "auto";
            dfs.num_regimes = 0;
            auto root_res = hyphy::analyses::DirectionalFELAnalyzer::auto_select_root(
                work_tree, map_aln, gtr_fit.params, dfs);
            hyphy::analyses::DirectionalFELAnalyzer dfel_probe(
                root_res.best_tree, map_aln, mg_params, dfs, comp_test, root_res);
            comp_test = dfel_probe.comp_test;
            directional_drift_active = (dfel_probe.regimes.size() == 3);
            if (!root_res.chosen_branch_name.empty()) {
                chosen_root_branch = root_res.chosen_branch_name;
                chosen_root_weight = root_res.chosen_weight;
            }
        }

        // Ensure every internal node has a deterministic unique name (Node1, Node2, ...)
        {
            std::unordered_set<std::string> used_names;
            for (const auto& nd : work_tree.nodes) {
                if (!nd.name.empty()) used_names.insert(nd.name);
            }
            size_t next_idx = 1;
            for (auto& nd : work_tree.nodes) {
                if (nd.id == work_tree.root_id) {
                    if (nd.name.empty()) nd.name = "Root";
                    continue;
                }
                if (!nd.is_leaf && nd.name.empty()) {
                    std::string cand;
                    do {
                        cand = "Node" + std::to_string(next_idx++);
                    } while (used_names.count(cand));
                    nd.name = cand;
                    used_names.insert(cand);
                }
            }
        }

        auto leaf_to_taxon = hyphy::core::LikelihoodEngine::build_leaf_to_taxon_map(work_tree, map_aln);
        std::vector<std::vector<size_t>> wt_leaves(work_tree.num_nodes());
        for (int32_t nid : work_tree.post_order) {
            if (nid < 0 || static_cast<size_t>(nid) >= work_tree.num_nodes()) continue;
            if (work_tree.nodes[nid].is_leaf) {
                if (nid < static_cast<int32_t>(leaf_to_taxon.size()) && leaf_to_taxon[nid] < N) {
                    wt_leaves[nid].push_back(leaf_to_taxon[nid]);
                }
            } else {
                for (int32_t ch : work_tree.nodes[nid].children) {
                    if (ch >= 0 && static_cast<size_t>(ch) < wt_leaves.size()) {
                        wt_leaves[nid].insert(wt_leaves[nid].end(), wt_leaves[ch].begin(), wt_leaves[ch].end());
                    }
                }
            }
        }

        // Compute per-taxon GC3 and overall GC
        std::vector<double> taxon_gc3(N, 0.50);
        std::vector<double> taxon_gc_all(N, 0.50);
        double min_gc3 = 1.0, max_gc3 = 0.0, sum_gc3 = 0.0;
        for (size_t t = 0; t < N; ++t) {
            const std::string& rseq = map_aln.sequences[t].raw_seq;
            size_t gc3 = 0, v3 = 0, gca = 0, va = 0;
            for (size_t c = 0; c + 2 < rseq.size(); c += 3) {
                for (size_t k = 0; k < 3; ++k) {
                    char ch = static_cast<char>(std::toupper(static_cast<unsigned char>(rseq[c + k])));
                    if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
                        va++;
                        if (ch == 'G' || ch == 'C') gca++;
                        if (k == 2) {
                            v3++;
                            if (ch == 'G' || ch == 'C') gc3++;
                        }
                    }
                }
            }
            taxon_gc3[t] = v3 > 0 ? static_cast<double>(gc3) / static_cast<double>(v3) : 0.50;
            taxon_gc_all[t] = va > 0 ? static_cast<double>(gca) / static_cast<double>(va) : 0.50;
            if (N > 0) {
                min_gc3 = std::min(min_gc3, taxon_gc3[t]);
                max_gc3 = std::max(max_gc3, taxon_gc3[t]);
                sum_gc3 += taxon_gc3[t];
            }
        }
        double mean_gc3 = N > 0 ? sum_gc3 / static_cast<double>(N) : 0.50;

        // 4. Localized SPR Recombination & Reticulation Bridge Scan (with 2,000 site-order permutations)
        hyphy::core::MG94Matrix mg94_base;
        mg94_base.update(mg_params, map_aln.pos_nuc_frequencies, map_aln.codon_frequencies_f3x4, *gcode);
        const double conv = (mg94_base.scale_factor > 1e-12) ? (3.0 / mg94_base.scale_factor) : 1.0;

        struct Stage2RecombBridge {
            int32_t node_u = -1;
            int32_t node_v = -1;
            std::string clade_u_name;
            std::string target_v_name;
            std::vector<std::string> recipient_leaves;
            std::vector<std::string> donor_leaves;
            double conductance = 0.0;
            size_t start_codon = 1;
            size_t end_codon = 1;
            double delta_log_l = 0.0;
            double p_value = 1.0;
            bool routed = false;
            std::string spr_tree_newick;
            std::vector<double> site_delta_ll;
        };
        std::vector<Stage2RecombBridge> recomb_bridges;

        if (enable_recomb_scan != 0 && N >= 4 && C >= 10) {
            std::vector<hyphy::core::Matrix> P_base(work_tree.num_nodes());
            for (const auto& nd : work_tree.nodes) {
                if (nd.id != work_tree.root_id && nd.parent_id != hyphy::core::INVALID_INDEX) {
                    P_base[nd.id] = mg94_base.transition_matrix(nd.branch_length * conv);
                }
            }
            std::vector<hyphy::core::SitePattern> site_pats(C);
            for (size_t c = 0; c < C; ++c) {
                site_pats[c].states.resize(N);
                for (size_t t = 0; t < N; ++t) site_pats[c].states[t] = map_aln.sequences[t].codons[c];
                site_pats[c].weight = 1;
            }
            std::vector<hyphy::core::Scalar> base_site_ll(C, 0.0);
            std::vector<hyphy::core::Vector> node_L(work_tree.num_nodes(), hyphy::core::Vector::Zero(S));
            std::vector<double> node_sc(work_tree.num_nodes(), 0.0);
            for (size_t c = 0; c < C; ++c) {
                base_site_ll[c] = hyphy::analyses::FELAnalyzer::compute_pattern_log_likelihood(
                    work_tree, site_pats[c], leaf_to_taxon, mg94_base, P_base, node_L, node_sc);
            }

            struct BridgeProposal {
                int32_t u_id;
                int32_t v_id;
                double score;
            };
            std::vector<BridgeProposal> proposals;

            const size_t win_w = std::min<size_t>(30, std::max<size_t>(10, C / 6));
            const size_t win_step = std::max<size_t>(5, win_w / 2);
            const size_t max_u_clade = (N > 128) ? 8 : (N / 2);

            for (const auto& u_nd : work_tree.nodes) {
                int32_t u = u_nd.id;
                if (u == work_tree.root_id || u_nd.parent_id == hyphy::core::INVALID_INDEX) continue;
                if (wt_leaves[u].empty() || wt_leaves[u].size() > max_u_clade) continue;

                int32_t p_u = u_nd.parent_id;
                std::vector<size_t> sib_leaves;
                for (int32_t sib : work_tree.nodes[p_u].children) {
                    if (sib != u && sib >= 0 && static_cast<size_t>(sib) < wt_leaves.size()) {
                        sib_leaves.insert(sib_leaves.end(), wt_leaves[sib].begin(), wt_leaves[sib].end());
                    }
                }
                if (sib_leaves.empty()) continue;

                int32_t best_v = hyphy::core::INVALID_INDEX;
                double best_v_score = 0.0;

                for (const auto& v_nd : work_tree.nodes) {
                    int32_t v = v_nd.id;
                    if (v == work_tree.root_id || v == u || v == p_u || v_nd.parent_id == hyphy::core::INVALID_INDEX) continue;
                    if (v_nd.parent_id == p_u) continue;
                    if (hyphy::analyses::TreeEngine::is_in_subtree(work_tree, v, u) ||
                        hyphy::analyses::TreeEngine::is_in_subtree(work_tree, u, v)) continue;
                    if (wt_leaves[v].empty()) continue;

                    int32_t p_v = v_nd.parent_id;
                    std::vector<size_t> v_sib_leaves;
                    if (p_v >= 0 && static_cast<size_t>(p_v) < work_tree.num_nodes()) {
                        for (int32_t vsib : work_tree.nodes[p_v].children) {
                            if (vsib != v && vsib >= 0 && static_cast<size_t>(vsib) < wt_leaves.size()) {
                                v_sib_leaves.insert(v_sib_leaves.end(), wt_leaves[vsib].begin(), wt_leaves[vsib].end());
                            }
                        }
                    }

                    double max_win_gain = 0.0;
                    for (size_t ws = 0; ws + win_w <= C; ws += win_step) {
                        double min_leaf_gain = 1e9;
                        double sum_leaf_gain = 0.0;
                        bool all_leaves_valid = true;
                        for (size_t lu : wt_leaves[u]) {
                            double sum_uv = 0.0, cnt_uv = 0.0;
                            for (size_t lv : wt_leaves[v]) {
                                for (size_t c = ws; c < ws + win_w; ++c) {
                                    int cu = map_aln.sequences[lu].codons[c];
                                    int cv = map_aln.sequences[lv].codons[c];
                                    if (cu >= 0 && cv >= 0) {
                                        cnt_uv += 1.0;
                                        if (cu == cv) sum_uv += 1.0;
                                    }
                                }
                            }
                            double sum_us = 0.0, cnt_us = 0.0;
                            for (size_t ls : sib_leaves) {
                                for (size_t c = ws; c < ws + win_w; ++c) {
                                    int cu = map_aln.sequences[lu].codons[c];
                                    int cs = map_aln.sequences[ls].codons[c];
                                    if (cu >= 0 && cs >= 0) {
                                        cnt_us += 1.0;
                                        if (cu == cs) sum_us += 1.0;
                                    }
                                }
                            }
                            if (cnt_uv < 8.0 || cnt_us < 8.0) {
                                all_leaves_valid = false;
                                break;
                            }
                            double lg = (sum_uv / cnt_uv) - (sum_us / cnt_us);
                            min_leaf_gain = std::min(min_leaf_gain, lg);
                            sum_leaf_gain += lg;
                        }
                        if (!all_leaves_valid || min_leaf_gain <= 0.08) continue;

                        double mean_gain = sum_leaf_gain / static_cast<double>(wt_leaves[u].size());
                        double donor_bonus = 0.0;
                        if (!v_sib_leaves.empty()) {
                            size_t lv0 = wt_leaves[v][0];
                            size_t lvs0 = v_sib_leaves[0];
                            int m_vv = 0, v_cnt = 0;
                            for (size_t c = ws; c < ws + win_w; ++c) {
                                int cv = map_aln.sequences[lv0].codons[c];
                                int cvs = map_aln.sequences[lvs0].codons[c];
                                if (cv >= 0 && cvs >= 0) {
                                    v_cnt++;
                                    if (cv == cvs) m_vv++;
                                }
                            }
                            if (v_cnt >= 8) {
                                donor_bonus = 0.25 * (static_cast<double>(m_vv) / static_cast<double>(v_cnt));
                            }
                        }
                        max_win_gain = std::max(max_win_gain, mean_gain + donor_bonus);
                    }
                    if (max_win_gain > best_v_score) {
                        best_v_score = max_win_gain;
                        best_v = v;
                    }
                }
                if (best_v != hyphy::core::INVALID_INDEX && best_v_score > 0.15) {
                    proposals.push_back({u, best_v, best_v_score});
                }
            }

            std::sort(proposals.begin(), proposals.end(), [](const auto& a, const auto& b) {
                return a.score > b.score;
            });
            if (proposals.size() > 4) proposals.resize(4);
            const double bonf_factor = static_cast<double>(std::max<size_t>(1, proposals.size()));

            for (const auto& prop : proposals) {
                int32_t u_id = prop.u_id;
                int32_t v_id = prop.v_id;
                auto spr_opt = hyphy::analyses::TreeEngine::apply_spr_move(work_tree, u_id, v_id, false);
                if (!spr_opt.has_value()) continue;
                hyphy::core::Tree br_tree = std::move(spr_opt->tree);
                int32_t w_id = spr_opt->w;

                std::vector<hyphy::core::Matrix> P_base_short = P_base;
                if (u_id >= 0 && static_cast<size_t>(u_id) < work_tree.num_nodes() &&
                    u_id != work_tree.root_id && work_tree.nodes[u_id].parent_id != hyphy::core::INVALID_INDEX) {
                    P_base_short[u_id] = mg94_base.transition_matrix(0.20 * work_tree.nodes[u_id].branch_length * conv);
                    int32_t p_u = work_tree.nodes[u_id].parent_id;
                    if (p_u >= 0 && p_u != work_tree.root_id && work_tree.nodes[p_u].children.size() == 2) {
                        int32_t s_u = (work_tree.nodes[p_u].children[0] == u_id)
                                    ? work_tree.nodes[p_u].children[1] : work_tree.nodes[p_u].children[0];
                        double old_ts = work_tree.nodes[s_u].branch_length;
                        double old_tp = work_tree.nodes[p_u].branch_length;
                        P_base_short[s_u] = mg94_base.transition_matrix(0.15 * old_ts * conv);
                        P_base_short[p_u] = mg94_base.transition_matrix((old_tp + 0.85 * old_ts) * conv);
                    }
                }

                std::vector<hyphy::core::Matrix> P_br(br_tree.num_nodes());
                std::vector<hyphy::core::Matrix> P_br_short(br_tree.num_nodes());
                for (const auto& nd : br_tree.nodes) {
                    if (nd.id != br_tree.root_id && nd.parent_id != hyphy::core::INVALID_INDEX) {
                        P_br[nd.id] = mg94_base.transition_matrix(nd.branch_length * conv);
                        P_br_short[nd.id] = P_br[nd.id];
                    }
                }
                if (u_id >= 0 && static_cast<size_t>(u_id) < br_tree.num_nodes() &&
                    u_id != br_tree.root_id && br_tree.nodes[u_id].parent_id != hyphy::core::INVALID_INDEX) {
                    P_br_short[u_id] = mg94_base.transition_matrix(0.20 * br_tree.nodes[u_id].branch_length * conv);
                }
                if (v_id >= 0 && static_cast<size_t>(v_id) < work_tree.num_nodes() &&
                    w_id >= 0 && static_cast<size_t>(w_id) < br_tree.num_nodes()) {
                    double old_tv = std::max<double>(1e-4, work_tree.nodes[v_id].branch_length);
                    P_br_short[v_id] = mg94_base.transition_matrix(0.15 * old_tv * conv);
                    P_br_short[w_id] = mg94_base.transition_matrix(0.85 * old_tv * conv);
                }

                auto l2t_br = hyphy::core::LikelihoodEngine::build_leaf_to_taxon_map(br_tree, map_aln);
                std::vector<hyphy::core::Vector> nL_br(br_tree.num_nodes(), hyphy::core::Vector::Zero(S));
                std::vector<double> nS_br(br_tree.num_nodes(), 0.0);

                std::vector<double> site_diffs(C, 0.0);
                double best_sum = 0.0, cur_sum = 0.0;
                size_t best_s = 0, best_e = 0, cur_s = 0;
                for (size_t c = 0; c < C; ++c) {
                    double ll_b1 = base_site_ll[c];
                    double ll_b2 = hyphy::analyses::FELAnalyzer::compute_pattern_log_likelihood(
                        work_tree, site_pats[c], leaf_to_taxon, mg94_base, P_base_short, node_L, node_sc);
                    double ll_base_sym = std::max(ll_b1, ll_b2);

                    double ll_r1 = hyphy::analyses::FELAnalyzer::compute_pattern_log_likelihood(
                        br_tree, site_pats[c], l2t_br, mg94_base, P_br, nL_br, nS_br);
                    double ll_r2 = hyphy::analyses::FELAnalyzer::compute_pattern_log_likelihood(
                        br_tree, site_pats[c], l2t_br, mg94_base, P_br_short, nL_br, nS_br);
                    double ll_br_sym = std::max(ll_r1, ll_r2);

                    double diff = ll_br_sym - ll_base_sym;
                    site_diffs[c] = diff;
                    if (cur_sum + diff > 0.0) {
                        cur_sum += diff;
                        if (cur_sum > best_sum) {
                            best_sum = cur_sum;
                            best_s = cur_s;
                            best_e = c;
                        }
                    } else {
                        cur_sum = 0.0;
                        cur_s = c + 1;
                    }
                }

                size_t seg_len = (best_e >= best_s) ? (best_e - best_s + 1) : 0;
                size_t pos_sites = 0;
                for (size_t c = best_s; c <= best_e && c < C; ++c) {
                    if (site_diffs[c] > 0.05) pos_sites++;
                }

                if (best_sum > 3.0 && seg_len >= 8 && pos_sites >= 3) {
                    constexpr int kNumPerms = 2000;
                    int exceed_count = 0;
                    std::vector<double> perm_diffs = site_diffs;
                    uint64_t rng_state = 0x9E3779B97F4A7C15ULL ^
                        (static_cast<uint64_t>(std::max<int32_t>(0, u_id)) * 1315423911ULL +
                         static_cast<uint64_t>(std::max<int32_t>(0, v_id)) * 2654435761ULL);
                    auto next_u64 = [&]() -> uint64_t {
                        rng_state ^= rng_state >> 12;
                        rng_state ^= rng_state << 25;
                        rng_state ^= rng_state >> 27;
                        return rng_state * 2685821657736338717ULL;
                    };
                    for (int p_it = 0; p_it < kNumPerms; ++p_it) {
                        for (size_t i = C - 1; i > 0; --i) {
                            size_t j = static_cast<size_t>(next_u64() % (i + 1));
                            std::swap(perm_diffs[i], perm_diffs[j]);
                        }
                        double p_best = 0.0, p_cur = 0.0;
                        for (size_t c = 0; c < C; ++c) {
                            p_cur = std::max<double>(0.0, p_cur + perm_diffs[c]);
                            if (p_cur > p_best) p_best = p_cur;
                        }
                        if (p_best >= best_sum - 1e-12) exceed_count++;
                    }
                    double p_perm = static_cast<double>(exceed_count + 1) / static_cast<double>(kNumPerms + 1);
                    double p_corr = std::min<double>(1.0, bonf_factor * p_perm);
                    double bic_thresh = 3.0 * std::log(std::max<double>(10.0, static_cast<double>(C)));

                    Stage2RecombBridge rb;
                    rb.node_u = u_id;
                    rb.node_v = v_id;
                    rb.clade_u_name = work_tree.nodes[u_id].name;
                    rb.target_v_name = work_tree.nodes[v_id].name;
                    for (size_t lu : wt_leaves[u_id]) rb.recipient_leaves.push_back(map_aln.sequences[lu].name);
                    for (size_t lv : wt_leaves[v_id]) rb.donor_leaves.push_back(map_aln.sequences[lv].name);
                    rb.conductance = prop.score;
                    rb.start_codon = best_s + 1;
                    rb.end_codon = best_e + 1;
                    rb.delta_log_l = best_sum;
                    rb.p_value = p_corr;
                    rb.routed = (p_corr <= 0.01 && (2.0 * best_sum > bic_thresh));
                    rb.spr_tree_newick = br_tree.to_newick();
                    rb.site_delta_ll = std::move(site_diffs);
                    recomb_bridges.push_back(std::move(rb));
                }
            }
        }

        // 5. Compute root-to-tip depths, branch lengths, clade GC3, and identify pre-tagged or Auto-Discover focal clade
        std::vector<double> node_depth(work_tree.num_nodes(), 0.0);
        for (auto it = work_tree.post_order.rbegin(); it != work_tree.post_order.rend(); ++it) {
            int32_t nid = *it;
            if (nid < 0 || static_cast<size_t>(nid) >= work_tree.num_nodes()) continue;
            const auto& nd = work_tree.nodes[nid];
            if (nd.id == work_tree.root_id || nd.parent_id == hyphy::core::INVALID_INDEX) {
                node_depth[nid] = 0.0;
            } else {
                node_depth[nid] = node_depth[nd.parent_id] + std::max<double>(0.0, nd.branch_length);
            }
        }

        double total_tree_length = 0.0;
        double max_root_to_tip = 0.0;
        size_t num_branches = 0;
        std::vector<double> all_bls;
        bool has_embedded_foreground = false;
        std::string embedded_foreground_tag;

        for (const auto& nd : work_tree.nodes) {
            if (nd.id != work_tree.root_id && nd.parent_id != hyphy::core::INVALID_INDEX) {
                double bl = std::max<double>(0.0, nd.branch_length);
                total_tree_length += bl;
                num_branches++;
                all_bls.push_back(bl);
                if (!nd.model_tag.empty()) {
                    has_embedded_foreground = true;
                    if (embedded_foreground_tag.empty()) embedded_foreground_tag = nd.model_tag;
                }
            }
            if (nd.is_leaf) {
                max_root_to_tip = std::max(max_root_to_tip, node_depth[nd.id]);
            }
        }

        double median_bl = 0.05;
        if (!all_bls.empty()) {
            std::vector<double> s_bls = all_bls;
            std::sort(s_bls.begin(), s_bls.end());
            median_bl = s_bls[s_bls.size() / 2];
        }

        // Fast empirical nonsynonymous/synonymous branch ratio screen to identify the top Auto-Discover focal clade candidate
        std::vector<int> node_anc_codon(work_tree.num_nodes() * C, -1);
        for (size_t c = 0; c < C; ++c) {
            for (int32_t nid : work_tree.post_order) {
                if (nid < 0 || static_cast<size_t>(nid) >= work_tree.num_nodes()) continue;
                const auto& nd = work_tree.nodes[nid];
                if (nd.is_leaf) {
                    size_t tid = (nid < static_cast<int32_t>(leaf_to_taxon.size())) ? leaf_to_taxon[nid] : hyphy::core::INVALID_INDEX;
                    node_anc_codon[nid * C + c] = (tid < N) ? map_aln.sequences[tid].codons[c] : -1;
                } else {
                    std::array<int, 64> counts{};
                    int best_st = -1, best_cnt = 0;
                    for (int32_t ch : nd.children) {
                        int st = (ch >= 0 && static_cast<size_t>(ch) < work_tree.num_nodes()) ? node_anc_codon[ch * C + c] : -1;
                        if (st >= 0 && st < 64) {
                            counts[st]++;
                            if (counts[st] > best_cnt) {
                                best_cnt = counts[st];
                                best_st = st;
                            }
                        }
                    }
                    node_anc_codon[nid * C + c] = best_st;
                }
            }
        }

        std::vector<double> branch_emp_omega(work_tree.num_nodes(), global_omega);
        for (const auto& nd : work_tree.nodes) {
            if (nd.id == work_tree.root_id || nd.parent_id == hyphy::core::INVALID_INDEX) continue;
            int32_t pid = nd.parent_id;
            double n_syn = 0.5, n_non = 0.5 * std::max<double>(0.1, global_omega);
            for (size_t c = 0; c < C; ++c) {
                int c_p = node_anc_codon[pid * C + c];
                int c_u = node_anc_codon[nd.id * C + c];
                if (c_p >= 0 && c_u >= 0 && c_p != c_u && c_p < S && c_u < S) {
                    if (gcode->sense_codon_aa[c_p] == gcode->sense_codon_aa[c_u]) n_syn += 1.0;
                    else n_non += 1.0;
                }
            }
            branch_emp_omega[nd.id] = std::clamp<double>((n_non / (2.2 * n_syn)), 0.02, 5.0);
        }

        int32_t auto_clade_node_id = -1;
        std::string auto_clade_name = "None (Homogeneous)";
        double auto_clade_z = 0.0;
        double auto_clade_k_est = 1.0;
        for (const auto& nd : work_tree.nodes) {
            if (nd.id == work_tree.root_id || nd.is_leaf) continue;
            if (wt_leaves[nd.id].size() < 2 || wt_leaves[nd.id].size() + 2 > N) continue;
            std::unordered_set<int32_t> sub_ids;
            std::vector<int32_t> stk = {nd.id};
            while (!stk.empty()) {
                int32_t cur = stk.back();
                stk.pop_back();
                sub_ids.insert(cur);
                for (int32_t ch : work_tree.nodes[cur].children) stk.push_back(ch);
            }
            double s1 = 0.0, s2 = 0.0;
            size_t n1 = 0, n2 = 0;
            for (const auto& b_nd : work_tree.nodes) {
                if (b_nd.id == work_tree.root_id || b_nd.parent_id == hyphy::core::INVALID_INDEX) continue;
                double lw = std::log(branch_emp_omega[b_nd.id]);
                if (sub_ids.count(b_nd.id)) { s1 += lw; n1++; }
                else { s2 += lw; n2++; }
            }
            if (n1 >= 2 && n2 >= 2) {
                double z = std::abs(s1 / n1 - s2 / n2) * std::sqrt(static_cast<double>(n1 * n2) / static_cast<double>(n1 + n2));
                if (z > auto_clade_z) {
                    auto_clade_z = z;
                    auto_clade_node_id = nd.id;
                    auto_clade_name = nd.name;
                    auto_clade_k_est = std::clamp<double>(std::exp((s1 / n1) - (s2 / n2)), 0.20, 4.0);
                }
            }
        }

        std::unordered_set<int32_t> auto_clade_subtree;
        if (auto_clade_node_id >= 0) {
            std::vector<int32_t> stk = {auto_clade_node_id};
            while (!stk.empty()) {
                int32_t cur = stk.back();
                stk.pop_back();
                auto_clade_subtree.insert(cur);
                for (int32_t ch : work_tree.nodes[cur].children) stk.push_back(ch);
            }
        }

        // Build nodes JSON array (in pre-order / layout-ready order)
        nlohmann::json nodes_json = nlohmann::json::array();
        for (const auto& nd : work_tree.nodes) {
            double clade_gc3 = mean_gc3;
            std::vector<std::string> clade_leaf_names;
            if (!wt_leaves[nd.id].empty()) {
                double s_gc = 0.0;
                for (size_t tid : wt_leaves[nd.id]) {
                    s_gc += taxon_gc3[tid];
                    clade_leaf_names.push_back(map_aln.sequences[tid].name);
                }
                clade_gc3 = s_gc / static_cast<double>(wt_leaves[nd.id].size());
            }
            std::string comp_regime = "Balanced";
            if (clade_gc3 >= mean_gc3 + 0.08) comp_regime = "GC-Rich";
            else if (clade_gc3 <= mean_gc3 - 0.08) comp_regime = "AT-Rich";

            bool in_embedded_fg = !nd.model_tag.empty();
            bool in_auto_fg = (auto_clade_subtree.count(nd.id) > 0);

            nodes_json.push_back({
                {"node_id", nd.id},
                {"parent_id", nd.parent_id},
                {"name", nd.name},
                {"model_tag", nd.model_tag},
                {"is_leaf", nd.is_leaf},
                {"is_root", (nd.id == work_tree.root_id)},
                {"branch_length", (nd.id == work_tree.root_id ? 0.0 : std::max<double>(0.0, nd.branch_length))},
                {"depth_subs", node_depth[nd.id]},
                {"children", nd.children},
                {"leaf_count", wt_leaves[nd.id].size()},
                {"leaf_names", clade_leaf_names},
                {"gc3", clade_gc3},
                {"compositional_regime", comp_regime},
                {"empirical_omega", branch_emp_omega[nd.id]},
                {"in_embedded_foreground", in_embedded_fg},
                {"in_auto_foreground", in_auto_fg},
                {"is_auto_clade_root", (nd.id == auto_clade_node_id)}
            });
        }

        // Build structured Events list (Exceptions-First progressive disclosure)
        nlohmann::json events_json = nlohmann::json::array();
        size_t routed_recomb_count = 0;
        for (size_t ri = 0; ri < recomb_bridges.size(); ++ri) {
            const auto& rb = recomb_bridges[ri];
            if (rb.routed) routed_recomb_count++;
            std::ostringstream exp;
            exp << "Codons " << rb.start_codon << "–" << rb.end_codon << " ("
                << (rb.end_codon - rb.start_codon + 1) << " codons) support an alternative local SPR topology moving clade '"
                << rb.clade_u_name << "' (" << rb.recipient_leaves.size() << " sequence(s)) to branch '"
                << rb.target_v_name << "' with ΔlnL = +" << std::fixed << std::setprecision(1) << rb.delta_log_l
                << " (permutation p = " << std::setprecision(4) << rb.p_value << "). ";
            if (rb.routed) {
                exp << "Stage 2 automatically routes this segment to its local SPR topology so downstream selection tests are not confounded by recombination.";
            } else {
                exp << "Signal is below the strict routing threshold (p <= 0.01 & Schwarz BIC) and is retained on the primary backbone tree.";
            }
            events_json.push_back({
                {"event_id", "recomb_" + std::to_string(ri + 1)},
                {"event_type", "recombination_bridge"},
                {"badge", rb.routed ? "Recombinant Segment (Routed)" : "Candidate Reticulation"},
                {"badge_style", rb.routed ? "quarantined" : "rescued"},
                {"title", "Codons " + std::to_string(rb.start_codon) + "–" + std::to_string(rb.end_codon) + ": " + rb.clade_u_name + " ↔ " + rb.target_v_name},
                {"subtitle", "ΔlnL = +" + (std::ostringstream() << std::fixed << std::setprecision(1) << rb.delta_log_l).str() + " • p = " + (std::ostringstream() << std::fixed << std::setprecision(4) << rb.p_value).str()},
                {"summary", (rb.routed ? "Routed to local SPR tree: " : "Localized discordant segment: ") + rb.clade_u_name + " groups with " + rb.target_v_name + " across codons " + std::to_string(rb.start_codon) + "–" + std::to_string(rb.end_codon)},
                {"explanation", exp.str()},
                {"node_u", rb.node_u},
                {"node_v", rb.node_v},
                {"clade_u_name", rb.clade_u_name},
                {"target_v_name", rb.target_v_name},
                {"recipient_leaves", rb.recipient_leaves},
                {"donor_leaves", rb.donor_leaves},
                {"start_codon", rb.start_codon},
                {"end_codon", rb.end_codon},
                {"delta_log_l", rb.delta_log_l},
                {"p_value", rb.p_value},
                {"routed", rb.routed},
                {"spr_tree_newick", rb.spr_tree_newick},
                {"site_delta_ll", rb.site_delta_ll}
            });
        }

        for (size_t mi = 0; mi < coupled_moves.size(); ++mi) {
            const auto& mv = coupled_moves[mi];
            std::string u_nm = (mv.clade_node_id >= 0 && static_cast<size_t>(mv.clade_node_id) < work_tree.nodes.size())
                             ? work_tree.nodes[mv.clade_node_id].name : ("Node" + std::to_string(mv.clade_node_id));
            std::string v_nm = (mv.target_node_id >= 0 && static_cast<size_t>(mv.target_node_id) < work_tree.nodes.size())
                             ? work_tree.nodes[mv.target_node_id].name : ("Node" + std::to_string(mv.target_node_id));
            std::ostringstream s_dl;
            s_dl << std::fixed << std::setprecision(1) << mv.delta_joint_log_l;
            events_json.push_back({
                {"event_id", "spr_" + std::to_string(mi + 1)},
                {"event_type", "coupled_spr"},
                {"badge", "Topology Upgraded"},
                {"badge_style", "rescued"},
                {"title", mv.move_type + ": " + u_nm + " → " + v_nm},
                {"subtitle", "ΔlnL = +" + s_dl.str()},
                {"summary", "Accepted coupled SPR rearrangement improving joint codon likelihood by +" + s_dl.str()},
                {"explanation", "Stage 2 coupled topology search relocated clade '" + u_nm + "' to '" + v_nm + "', improving the MG94 codon log-likelihood by +" + s_dl.str() + "."},
                {"node_u", mv.clade_node_id},
                {"node_v", mv.target_node_id}
            });
        }

        if (comp_test.is_heterogeneous) {
            std::ostringstream exp;
            exp << "Third-codon-position GC content (GC3%) ranges from "
                << std::fixed << std::setprecision(1) << (min_gc3 * 100.0) << "% to " << (max_gc3 * 100.0)
                << "% across lineages (χ² = " << comp_test.chi2_cp3 << ", p = " << std::scientific << std::setprecision(2) << comp_test.p_value_cp3
                << "). Stage 2 selected root branch '" << chosen_root_branch << "' (confidence weight "
                << std::fixed << std::setprecision(2) << chosen_root_weight << ") to account for directional compositional drift.";
            events_json.push_back({
                {"event_id", "comp_drift"},
                {"event_type", "compositional_drift"},
                {"badge", directional_drift_active ? "3-Regime GC3 Drift" : "GC3 Heterogeneity"},
                {"badge_style", "rescued"},
                {"title", "GC3 Compositional Drift (" + std::to_string(static_cast<int>(std::round(min_gc3 * 100))) + "% – " + std::to_string(static_cast<int>(std::round(max_gc3 * 100))) + "% GC3)"},
                {"subtitle", "Root: " + chosen_root_branch},
                {"summary", "Significant GC3 compositional variation across lineages; auto-rooted on " + chosen_root_branch},
                {"explanation", exp.str()}
            });
        }

        // Check for extreme long branches (> 0.60 subs/codon and > 4x median)
        for (const auto& nd : work_tree.nodes) {
            if (nd.id == work_tree.root_id || nd.parent_id == hyphy::core::INVALID_INDEX) continue;
            if (nd.branch_length > 0.60 && nd.branch_length > 4.0 * std::max<double>(0.05, median_bl)) {
                std::ostringstream bl_s;
                bl_s << std::fixed << std::setprecision(2) << nd.branch_length;
                events_json.push_back({
                    {"event_id", "long_branch_" + std::to_string(nd.id)},
                    {"event_type", "long_branch"},
                    {"badge", "Long Branch"},
                    {"badge_style", "rescued"},
                    {"title", "Divergent Branch: " + nd.name},
                    {"subtitle", bl_s.str() + " subs/codon"},
                    {"summary", "Branch '" + nd.name + "' spans " + bl_s.str() + " substitutions/codon (>4× median branch length)"},
                    {"explanation", "High evolutionary divergence on branch '" + nd.name + "' (" + bl_s.str() + " subs/codon). Stage 3 automatically activates its Markov error-sink HMM and Empirical Bayes synonymous rate shrinkage to prevent long-branch false positives."},
                    {"node_u", nd.id}
                });
            }
        }

        // Build partitioned NEXUS export (with trees block for backbone + any routed recombinant segments)
        std::ostringstream part_nexus;
        part_nexus << "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=" << N << " NCHAR=" << (3 * C) << ";\n"
                   << "  FORMAT DATATYPE=DNA GAP=- MISSING=?;\n  MATRIX\n";
        for (size_t t = 0; t < N; ++t) {
            part_nexus << "    '" << escape_nexus_single_quotes(map_aln.sequences[t].name) << "'  " << map_aln.sequences[t].raw_seq << "\n";
        }
        part_nexus << "  ;\nEND;\n\nBEGIN TREES;\n"
                   << "  TREE backbone_tree = " << work_tree.to_newick() << "\n";
        for (size_t ri = 0; ri < recomb_bridges.size(); ++ri) {
            part_nexus << "  TREE recomb_segment_" << (ri + 1) << "_codons_"
                       << recomb_bridges[ri].start_codon << "_" << recomb_bridges[ri].end_codon
                       << " = " << recomb_bridges[ri].spr_tree_newick << "\n";
        }
        part_nexus << "END;\n";

        auto t_end = std::chrono::steady_clock::now();
        double elapsed_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

        nlohmann::json out_j = {
            {"stage", 2},
            {"stage_name", "Stage 2: Phylogeny, Rooting & Recombination"},
            {"dataset_label", lbl},
            {"genetic_code", cname},
            {"elapsed_ms", elapsed_ms},
            {"tree_source", tree_source},
            {"summary", {
                {"num_taxa", N},
                {"num_codons", C},
                {"num_branches", num_branches},
                {"total_tree_length", safe_json_num(total_tree_length)},
                {"max_root_to_tip", safe_json_num(max_root_to_tip)},
                {"median_branch_length", safe_json_num(median_bl)},
                {"global_omega", safe_json_num(global_omega)},
                {"initial_mg94_log_l", safe_json_num(initial_mg94_log_l)},
                {"final_mg94_log_l", safe_json_num(final_mg94_log_l)},
                {"delta_spr_log_l", safe_json_num(final_mg94_log_l - initial_mg94_log_l)},
                {"coupled_spr_moves_count", coupled_moves.size()},
                {"recombination_bridges_count", recomb_bridges.size()},
                {"routed_recombination_count", routed_recomb_count},
                {"tree_newick", work_tree.to_newick()},
                {"has_embedded_foreground", has_embedded_foreground},
                {"embedded_foreground_tag", embedded_foreground_tag},
                {"auto_clade_node_id", auto_clade_node_id},
                {"auto_clade_name", auto_clade_name},
                {"auto_clade_z", safe_json_num(auto_clade_z)},
                {"auto_clade_k_est", safe_json_num(auto_clade_k_est)}
            }},
            {"compositional_audit", {
                {"is_heterogeneous", comp_test.is_heterogeneous},
                {"directional_drift_active", directional_drift_active},
                {"chi2_overall", safe_json_num(comp_test.chi2_overall)},
                {"p_value_overall", safe_json_num(comp_test.p_value_overall)},
                {"chi2_cp3", safe_json_num(comp_test.chi2_cp3)},
                {"p_value_cp3", safe_json_num(comp_test.p_value_cp3)},
                {"min_gc3", safe_json_num(min_gc3)},
                {"max_gc3", safe_json_num(max_gc3)},
                {"mean_gc3", safe_json_num(mean_gc3)},
                {"chosen_root_branch", chosen_root_branch},
                {"chosen_root_weight", safe_json_num(chosen_root_weight)}
            }},
            {"nodes", nodes_json},
            {"events", events_json},
            {"exports", {
                {"tree_newick", work_tree.to_newick()},
                {"partitioned_nexus", part_nexus.str()}
            }}
        };

        g_last_stage2_json = out_j.dump();
        return g_last_stage2_json.c_str();
    } catch (const std::exception& e) {
        g_last_error = e.what();
        return "";
    } catch (...) {
        g_last_error = "Unhandled internal exception in WASM engine.";
        return "";
    }
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_get_last_html() {
    return g_last_html_report.c_str();
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_render_html_from_json(const char* raw_json, const char* dataset_label) {
    try {
        g_last_error.clear();
        std::string j_str = raw_json ? raw_json : "{}";
        std::string lbl = (dataset_label && dataset_label[0] != '\0') ? dataset_label : "WALRUS_SONAR_TARGET";
        g_last_json_report = j_str;
        g_last_html_report = hyphy::analyses::CharacterizeEngine::render_standalone_html_from_json_string(j_str, lbl);
        return g_last_html_report.c_str();
    } catch (const std::exception& e) {
        g_last_error = e.what();
        return "";
    } catch (...) {
        g_last_error = "Unhandled internal exception in WASM engine.";
        return "";
    }
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_stage3_run(
    const char* seq_content,
    const char* tree_newick_in,
    const char* dataset_label,
    const char* code_name
) {
    try {
        auto t_start = std::chrono::steady_clock::now();
        g_last_error.clear();
        std::string raw_seq = seq_content ? seq_content : "";
        if (raw_seq.empty()) {
            throw std::runtime_error("Empty sequence payload supplied to Stage 3 Selection & Forensics.");
        }
        std::string nwk_in = tree_newick_in ? tree_newick_in : "";
        std::string lbl = (dataset_label && dataset_label[0] != '\0') ? dataset_label : "STAGE3_COHORT";
        std::string cname = (code_name && code_name[0] != '\0') ? code_name : "Universal";

        size_t p0 = 0;
        while (p0 < raw_seq.size() && std::isspace(static_cast<unsigned char>(raw_seq[p0]))) ++p0;
        if (p0 < raw_seq.size() && raw_seq[p0] == '{') {
            g_last_json_report = raw_seq;
            g_last_stage3_json = raw_seq;
            g_last_html_report = hyphy::analyses::CharacterizeEngine::render_standalone_html_from_json_string(raw_seq, lbl);
            return g_last_stage3_json.c_str();
        }

        const bool is_nexus = (p0 < raw_seq.size() && raw_seq[p0] == '#');
        const std::string vpath = is_nexus ? "/tmp/walrus_stage3_in.nex" : "/tmp/walrus_stage3_in.fa";
        ScopedTempFile tmp_guard(vpath);
        {
            std::ofstream vf(vpath, std::ios::binary | std::ios::trunc);
            vf << raw_seq;
        }

        auto gcode = safe_genetic_code(cname);

        hyphy::analyses::CharacterizeSettings settings;
        settings.quiet = true;
        settings.genetic_code_name = cname;
        settings.enable_coevolution = true;

        std::optional<hyphy::core::Tree> opt_tree = std::nullopt;
        std::string effective_nwk = nwk_in;
        if (effective_nwk.empty()) {
            effective_nwk = extract_translated_embedded_newick(raw_seq);
        }
        if (!effective_nwk.empty()) {
            try {
                hyphy::core::Tree ut = hyphy::core::Tree::from_newick(effective_nwk);
                apply_nexus_translate_to_tree(ut, raw_seq);
                for (const auto& nd : ut.nodes) {
                    if (!nd.model_tag.empty()) {
                        settings.foreground_tag = nd.model_tag;
                        break;
                    }
                }
                opt_tree = std::move(ut);
                settings.enable_coupled_spr = false; // Stage 2 or embedded tree already finalized topology & foreground tags
            } catch (...) {
                opt_tree = std::nullopt;
            }
        }

        hyphy::analyses::CharacterizeEngine engine(settings);
        hyphy::analyses::CharacterizeResult res;

        if (looks_unaligned_or_needs_stage0(raw_seq)) {
            auto seqs = hyphy::analyses::JointAlignEngine::load_unaligned_fasta(vpath, gcode);
            res = engine.run_from_unaligned(seqs, gcode, opt_tree);
        } else {
            try {
                hyphy::core::Alignment aln = hyphy::core::Alignment::load(vpath, gcode);
                if (!effective_nwk.empty()) aln.embedded_tree_newick = effective_nwk;
                res = engine.run_from_alignment(aln, opt_tree);
            } catch (...) {
                auto seqs = hyphy::analyses::JointAlignEngine::load_unaligned_fasta(vpath, gcode);
                res = engine.run_from_unaligned(seqs, gcode, opt_tree);
            }
        }

        g_last_json_report = hyphy::analyses::CharacterizeEngine::to_json(res);
        g_last_html_report = hyphy::analyses::CharacterizeEngine::render_standalone_html_from_json_string(g_last_json_report, lbl);

        // Enrich Stage 3 JSON with tree_nodes containing ancestral & terminal codon_seq and aa_seq
        // so the browser can render visual codon-substitution trees for ANY codon site!
        nlohmann::json out_j = nlohmann::json::parse(g_last_json_report);
        const auto& map_aln = res.joint_align.map_alignment;
        hyphy::core::Tree work_tree = res.joint_align.optimized_tree;
        hyphy::core::LikelihoodEngine::reconcile_tree_leaf_names(work_tree, map_aln);

        const size_t N = map_aln.num_taxa;
        const size_t C = map_aln.num_codons;
        const int S = static_cast<int>(gcode->num_sense_codons);

        // Ensure deterministic node names
        {
            std::unordered_set<std::string> used_names;
            for (const auto& nd : work_tree.nodes) {
                if (!nd.name.empty()) used_names.insert(nd.name);
            }
            size_t next_idx = 1;
            for (auto& nd : work_tree.nodes) {
                if (nd.id == work_tree.root_id) {
                    if (nd.name.empty()) nd.name = "Root";
                    continue;
                }
                if (!nd.is_leaf && nd.name.empty()) {
                    std::string cand;
                    do {
                        cand = "Node" + std::to_string(next_idx++);
                    } while (used_names.count(cand));
                    nd.name = cand;
                    used_names.insert(cand);
                }
            }
        }

        auto leaf_to_taxon = hyphy::core::LikelihoodEngine::build_leaf_to_taxon_map(work_tree, map_aln);
        std::vector<double> node_depth(work_tree.num_nodes(), 0.0);
        for (int32_t i = static_cast<int32_t>(work_tree.post_order.size()) - 1; i >= 0; --i) {
            int32_t nid = work_tree.post_order[i];
            const auto& nd = work_tree.nodes[nid];
            if (nd.id != work_tree.root_id && nd.parent_id != hyphy::core::INVALID_INDEX &&
                static_cast<size_t>(nd.parent_id) < node_depth.size()) {
                node_depth[nid] = node_depth[nd.parent_id] + std::max<double>(0.0, nd.branch_length);
            }
        }

        // Empirical Bayes MAP ancestral codon state reconstruction across all C codons (with Fitch fallback)
        std::vector<int> node_anc_codon =
            (res.node_ancestral_codons.size() == work_tree.num_nodes() * C)
                ? res.node_ancestral_codons
                : std::vector<int>(work_tree.num_nodes() * C, -1);
        for (size_t c = 0; c < C; ++c) {
            std::vector<uint64_t> mask(work_tree.num_nodes(), 0ULL);
            for (int32_t nid : work_tree.post_order) {
                if (work_tree.nodes[nid].is_leaf) {
                    size_t tid = (nid < static_cast<int32_t>(leaf_to_taxon.size())) ? leaf_to_taxon[nid] : N;
                    int st = (tid < N && c < map_aln.sequences[tid].codons.size()) ? map_aln.sequences[tid].codons[c] : -1;
                    if (st >= 0 && st < S) {
                        mask[nid] = (1ULL << st);
                        node_anc_codon[nid * C + c] = st;
                    } else {
                        mask[nid] = (S >= 64) ? ~0ULL : ((1ULL << S) - 1ULL);
                    }
                } else {
                    uint64_t inter = ~0ULL;
                    uint64_t uni = 0ULL;
                    for (int32_t ch : work_tree.nodes[nid].children) {
                        inter &= mask[ch];
                        uni |= mask[ch];
                    }
                    mask[nid] = (inter != 0ULL) ? inter : uni;
                }
            }
            for (int32_t i = static_cast<int32_t>(work_tree.post_order.size()) - 1; i >= 0; --i) {
                int32_t nid = work_tree.post_order[i];
                if (node_anc_codon[nid * C + c] >= 0) continue;
                int32_t pid = work_tree.nodes[nid].parent_id;
                int p_st = (pid != hyphy::core::INVALID_INDEX && pid >= 0) ? node_anc_codon[pid * C + c] : -1;
                if (p_st >= 0 && (mask[nid] & (1ULL << p_st))) {
                    node_anc_codon[nid * C + c] = p_st;
                } else {
                    int best_st = -1;
                    for (int s = 0; s < S; ++s) {
                        if (mask[nid] & (1ULL << s)) { best_st = s; break; }
                    }
                    node_anc_codon[nid * C + c] = best_st;
                }
            }
        }

        nlohmann::json tree_nodes_json = nlohmann::json::array();
        for (const auto& nd : work_tree.nodes) {
            std::string codon_seq;
            std::string aa_seq;
            codon_seq.reserve(3 * C);
            aa_seq.reserve(C);

            size_t tid = (nd.is_leaf && nd.id >= 0 && static_cast<size_t>(nd.id) < leaf_to_taxon.size())
                             ? leaf_to_taxon[nd.id]
                             : N;
            for (size_t c = 0; c < C; ++c) {
                if (nd.is_leaf && tid < N && 3 * c + 2 < map_aln.sequences[tid].raw_seq.size()) {
                    std::string trip = map_aln.sequences[tid].raw_seq.substr(3 * c, 3);
                    codon_seq += trip;
                    int8_t sc = (c < map_aln.sequences[tid].codons.size()) ? map_aln.sequences[tid].codons[c] : -1;
                    if (sc >= 0 && sc < S) {
                        aa_seq += gcode->sense_codon_aa[sc];
                    } else if (trip == "---") {
                        aa_seq += '-';
                    } else {
                        aa_seq += 'X';
                    }
                } else {
                    int st = node_anc_codon[nd.id * C + c];
                    if (st >= 0 && st < S && static_cast<size_t>(st) < gcode->sense_codons.size()) {
                        codon_seq += gcode->sense_codons[st];
                        aa_seq += gcode->sense_codon_aa[st];
                    } else {
                        codon_seq += "---";
                        aa_seq += '-';
                    }
                }
            }

            tree_nodes_json.push_back({
                {"node_id", nd.id},
                {"parent_id", nd.parent_id},
                {"name", nd.name},
                {"model_tag", nd.model_tag},
                {"is_leaf", nd.is_leaf},
                {"is_root", (nd.id == work_tree.root_id)},
                {"branch_length", (nd.id == work_tree.root_id ? 0.0 : std::max<double>(0.0, nd.branch_length))},
                {"depth_subs", node_depth[nd.id]},
                {"children", nd.children},
                {"codon_seq", codon_seq},
                {"aa_seq", aa_seq}
            });
        }

        auto t_end = std::chrono::steady_clock::now();
        out_j["stage"] = 3;
        out_j["stage_name"] = "Stage 3: Selection & Evolutionary Forensics";
        out_j["dataset_label"] = lbl;
        out_j["elapsed_ms"] = std::chrono::duration<double, std::milli>(t_end - t_start).count();
        out_j["tree_nodes"] = std::move(tree_nodes_json);

        g_last_stage3_json = out_j.dump();
        return g_last_stage3_json.c_str();
    } catch (const std::exception& e) {
        g_last_error = e.what();
        return "";
    } catch (...) {
        g_last_error = "Unhandled internal exception in WASM engine.";
        return "";
    }
}

HYPHY3_WASM_KEEPALIVE
const char* walrus_ping_run(const char* seq_content, const char* dataset_label, const char* code_name) {
    try {
        g_last_error.clear();
        std::string raw_seq = seq_content ? seq_content : "";
        if (raw_seq.empty()) {
            throw std::runtime_error("Empty sequence payload supplied to ONE PING ONLY.");
        }
        std::string lbl = (dataset_label && dataset_label[0] != '\0') ? dataset_label : "RED_OCTOBER_TARGET";
        std::string cname = (code_name && code_name[0] != '\0') ? code_name : "Universal";

        // Trim leading whitespace to see if the payload is already a .CHARACTERIZE.json document
        size_t p0 = 0;
        while (p0 < raw_seq.size() && std::isspace(static_cast<unsigned char>(raw_seq[p0]))) ++p0;
        if (p0 < raw_seq.size() && raw_seq[p0] == '{') {
            g_last_json_report = raw_seq;
            g_last_html_report = hyphy::analyses::CharacterizeEngine::render_standalone_html_from_json_string(raw_seq, lbl);
            return g_last_html_report.c_str();
        }

        // Write input payload to Emscripten MEMFS virtual file
        const bool is_nexus = (p0 < raw_seq.size() && raw_seq[p0] == '#');
        const std::string vpath = is_nexus ? "/tmp/walrus_input.nex" : "/tmp/walrus_input.fa";
        ScopedTempFile tmp_guard(vpath);
        {
            std::ofstream vf(vpath, std::ios::binary | std::ios::trunc);
            vf << raw_seq;
        }

        auto gcode = safe_genetic_code(cname);
        hyphy::analyses::CharacterizeSettings settings;
        settings.quiet = false;

        hyphy::analyses::CharacterizeEngine engine(settings);
        hyphy::analyses::CharacterizeResult res;

        if (looks_unaligned_or_needs_stage0(raw_seq)) {
            res = engine.run_from_unaligned_file(vpath, gcode);
        } else {
            try {
                hyphy::core::Alignment aln = hyphy::core::Alignment::load(vpath, gcode);
                std::string tr_nwk = extract_translated_embedded_newick(raw_seq);
                if (!tr_nwk.empty()) aln.embedded_tree_newick = tr_nwk;
                res = engine.run_from_alignment(aln);
            } catch (...) {
                // Fallback to Stage 0/1 unaligned transducer if strict alignment load fails
                res = engine.run_from_unaligned_file(vpath, gcode);
            }
        }

        g_last_json_report = hyphy::analyses::CharacterizeEngine::to_json(res);
        g_last_html_report = hyphy::analyses::CharacterizeEngine::render_standalone_html_from_json_string(g_last_json_report, lbl);
        return g_last_html_report.c_str();
    } catch (const std::exception& e) {
        g_last_error = e.what();
        std::cerr << "[SONAR FAULT] " << g_last_error << "\n";
        return "";
    } catch (...) {
        g_last_error = "Unhandled internal exception in WASM engine.";
        return "";
    }
}

} // extern "C"
