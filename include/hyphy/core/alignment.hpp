#pragma once

#include "hyphy/core/types.hpp"
#include "hyphy/core/genetic_code.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <memory>

namespace hyphy::core {

struct Sequence {
    std::string name;
    std::string raw_seq;
    std::vector<int8_t> codons; // sense codon indices (or -1 for gap/stop/ambiguous)
    std::vector<uint64_t> codon_masks; // resolved sense-codon bitmask for each codon site
};

struct SitePattern {
    std::vector<int8_t> states; // State per taxon (size = num_taxa)
    std::vector<uint64_t> ambiguous_masks; // Bitmask of compatible states when states[t] < 0
    size_t weight = 0;          // Number of alignment sites matching this pattern
    std::vector<size_t> sites;  // Original site indices
};

class Alignment {
public:
    std::vector<std::string> taxon_names;
    std::unordered_map<std::string, size_t> taxon_to_index;
    std::unordered_map<std::string, size_t> lower_taxon_to_index;
    std::vector<Sequence> sequences;
    size_t num_taxa = 0;
    size_t num_codons = 0;
    size_t num_nucleotides = 0;

    // Genetic code
    std::shared_ptr<const GeneticCode> code;

    // Embedded Newick tree if parsed from NEXUS
    std::string embedded_tree_newick;

    // Pattern compression
    std::vector<SitePattern> patterns;
    std::vector<size_t> site_to_pattern;

    // Base frequencies
    Vector4 nuc_frequencies = Vector4::Zero();

    // Positional base frequencies (for F3x4: 3 positions x 4 bases)
    Eigen::Matrix<Scalar, 3, 4> pos_nuc_frequencies = Eigen::Matrix<Scalar, 3, 4>::Zero();

    // Codon frequencies (sense codons)
    Vector codon_frequencies_f1x4;
    Vector codon_frequencies_f3x4;

    // Optional alignment quality & Partial Order Alignment (POA) metadata from JointAlignEngine
    std::vector<Scalar> site_homology_posteriors; // Posterior column confidence in [0, 1] per codon site
    std::string poa_gfa;                          // GFAv1 Partial Order Alignment codon DAG representation

    static std::string read_file_or_gzip(const std::string& filepath) {
        std::ifstream test_f(filepath, std::ios::binary);
        if (!test_f.is_open()) {
            throw std::runtime_error("Could not open alignment file: " + filepath);
        }
        unsigned char magic[2] = {0, 0};
        test_f.read(reinterpret_cast<char*>(magic), 2);
        bool is_gz = (magic[0] == 0x1f && magic[1] == 0x8b) ||
                     (filepath.size() >= 3 && filepath.substr(filepath.size() - 3) == ".gz");
        test_f.close();

        if (is_gz) {
            std::string cmd = "gzip -dc \"" + filepath + "\"";
            FILE* pipe = popen(cmd.c_str(), "r");
            if (!pipe) {
                throw std::runtime_error("Failed to decompress gzip alignment file: " + filepath);
            }
            std::string out;
            char buf[65536];
            size_t nread = 0;
            while ((nread = std::fread(buf, 1, sizeof(buf), pipe)) > 0) {
                out.append(buf, nread);
            }
            pclose(pipe);
            if (out.empty()) {
                throw std::runtime_error("Decompressed gzip alignment is empty: " + filepath);
            }
            return out;
        }

        std::ifstream file(filepath);
        std::stringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }

    static Alignment from_fasta(const std::string& filepath, std::shared_ptr<const GeneticCode> gcode = nullptr, bool as_nucleotide = false) {
        Alignment aln;
        std::string raw = read_file_or_gzip(filepath);
        std::istringstream file(raw);

        std::string line;
        std::string current_name;
        std::string current_seq;

        while (std::getline(file, line)) {
            // Trim trailing \r or whitespace
            while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
                line.pop_back();
            }
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!current_name.empty()) {
                    aln.add_sequence(current_name, current_seq);
                    current_seq.clear();
                }
                current_name = line.substr(1);
                // Strip spaces
                auto first_space = current_name.find_first_of(" \t");
                if (first_space != std::string::npos) {
                    current_name = current_name.substr(0, first_space);
                }
            } else {
                current_seq += line;
            }
        }
        if (!current_name.empty()) {
            aln.add_sequence(current_name, current_seq);
        }

        if (as_nucleotide) {
            aln.finalize_nucleotide();
        } else {
            aln.finalize(gcode);
        }
        return aln;
    }

    static Alignment from_nexus(const std::string& filepath, std::shared_ptr<const GeneticCode> gcode = nullptr, bool as_nucleotide = false) {
        Alignment aln;
        std::string raw = read_file_or_gzip(filepath);

        // 1. Remove comments [ ... ]
        std::string text;
        text.reserve(raw.size());
        int comment_depth = 0;
        for (char c : raw) {
            if (c == '[') {
                comment_depth++;
            } else if (c == ']') {
                if (comment_depth > 0) comment_depth--;
            } else if (comment_depth == 0) {
                text += c;
            }
        }

        // Case-insensitive search helper
        auto find_case_insensitive = [](const std::string& haystack, const std::string& needle, size_t start_pos = 0) -> size_t {
            auto it = std::search(
                haystack.begin() + start_pos, haystack.end(),
                needle.begin(), needle.end(),
                [](char ch1, char ch2) { return std::toupper(static_cast<unsigned char>(ch1)) == std::toupper(static_cast<unsigned char>(ch2)); }
            );
            return (it != haystack.end()) ? std::distance(haystack.begin(), it) : std::string::npos;
        };

        // 2. Find MATRIX block
        size_t mat_pos = find_case_insensitive(text, "MATRIX");
        if (mat_pos == std::string::npos) {
            throw std::runtime_error("NEXUS file does not contain a MATRIX block: " + filepath);
        }

        // Parse optional TAXLABELS block and FORMAT NOLABELS / DIMENSIONS NCHAR before MATRIX
        std::vector<std::string> taxlabels;
        size_t tl_pos = find_case_insensitive(text, "TAXLABELS");
        if (tl_pos != std::string::npos && tl_pos < mat_pos) {
            size_t tl_end = text.find(';', tl_pos + 9);
            if (tl_end != std::string::npos && tl_end < mat_pos) {
                std::string tl_body = text.substr(tl_pos + 9, tl_end - (tl_pos + 9));
                size_t ti = 0, tn = tl_body.size();
                while (ti < tn) {
                    while (ti < tn && std::isspace(static_cast<unsigned char>(tl_body[ti]))) ti++;
                    if (ti >= tn) break;
                    std::string name;
                    if (tl_body[ti] == '\'' || tl_body[ti] == '"') {
                        char q = tl_body[ti++];
                        while (ti < tn && tl_body[ti] != q) name += tl_body[ti++];
                        if (ti < tn) ti++;
                    } else {
                        while (ti < tn && !std::isspace(static_cast<unsigned char>(tl_body[ti]))) {
                            name += tl_body[ti++];
                        }
                    }
                    if (!name.empty()) taxlabels.push_back(name);
                }
            }
        }

        bool nolabels = false;
        size_t fmt_pos = find_case_insensitive(text, "FORMAT");
        while (fmt_pos != std::string::npos && fmt_pos < mat_pos) {
            size_t fmt_end = text.find(';', fmt_pos);
            if (fmt_end != std::string::npos && fmt_end <= mat_pos) {
                std::string fmt_block = text.substr(fmt_pos, fmt_end - fmt_pos);
                if (find_case_insensitive(fmt_block, "NOLABELS") != std::string::npos) {
                    nolabels = true;
                    break;
                }
            }
            fmt_pos = find_case_insensitive(text, "FORMAT", fmt_pos + 6);
        }

        mat_pos += 6; // skip "MATRIX"

        size_t mat_end = text.find(';', mat_pos);
        if (mat_end == std::string::npos) {
            throw std::runtime_error("NEXUS MATRIX block is not terminated with a semicolon: " + filepath);
        }

        std::string mat_content = text.substr(mat_pos, mat_end - mat_pos);

        // Parse sequences in MATRIX
        size_t idx = 0;
        size_t n = mat_content.size();
        size_t nolabel_row_idx = 0;
        while (idx < n) {
            // skip whitespace
            while (idx < n && std::isspace(static_cast<unsigned char>(mat_content[idx]))) {
                idx++;
            }
            if (idx >= n) break;

            std::string taxon;
            if (nolabels) {
                if (taxlabels.empty()) {
                    throw std::runtime_error("NEXUS FORMAT specifies NOLABELS, but no TAXLABELS were found: " + filepath);
                }
                taxon = taxlabels[nolabel_row_idx % taxlabels.size()];
                nolabel_row_idx++;
            } else {
                if (mat_content[idx] == '\'' || mat_content[idx] == '"') {
                    char quote = mat_content[idx++];
                    while (idx < n && mat_content[idx] != quote) {
                        taxon += mat_content[idx++];
                    }
                    if (idx < n) idx++; // skip closing quote
                } else {
                    while (idx < n && !std::isspace(static_cast<unsigned char>(mat_content[idx]))) {
                        taxon += mat_content[idx++];
                    }
                }

                // skip whitespace to sequence
                while (idx < n && std::isspace(static_cast<unsigned char>(mat_content[idx]))) {
                    idx++;
                }
                if (idx >= n) break;
            }

            // read sequence characters
            std::string seq;
            while (idx < n) {
                char c = mat_content[idx];
                if (std::isalpha(static_cast<unsigned char>(c)) || c == '-' || c == '?' || c == '.') {
                    seq += (c == '.' ? '-' : c);
                    idx++;
                } else if (c == '\n' || c == '\r') {
                    idx++;
                    break;
                } else if (std::isspace(static_cast<unsigned char>(c))) {
                    if (nolabels) {
                        // In NOLABELS mode, horizontal spaces within a line are ignored; newline ends row
                        idx++;
                        continue;
                    }
                    // Check if next non-space character is the start of a new taxon
                    size_t lookahead = idx;
                    while (lookahead < n && std::isspace(static_cast<unsigned char>(mat_content[lookahead]))) lookahead++;
                    if (lookahead < n && (mat_content[lookahead] == '\'' || mat_content[lookahead] == '"')) {
                        break;
                    }
                    // In single-line contiguous NEXUS, sequence is contiguous until next taxon
                    break;
                } else {
                    idx++;
                }
            }

            if (!taxon.empty() && !seq.empty()) {
                auto it = aln.taxon_to_index.find(taxon);
                if (it == aln.taxon_to_index.end()) {
                    aln.add_sequence(taxon, seq);
                } else {
                    aln.sequences[it->second].raw_seq += seq;
                }
            }
        }

        // 3. Find embedded tree in BEGIN TREES; TREE <name> = <newick>;
        size_t tree_pos = find_case_insensitive(text, "TREE ");
        if (tree_pos == std::string::npos) {
            tree_pos = find_case_insensitive(text, "TREE\t");
        }
        if (tree_pos != std::string::npos) {
            size_t eq_pos = text.find('=', tree_pos);
            if (eq_pos != std::string::npos) {
                size_t lparen = text.find('(', eq_pos);
                if (lparen != std::string::npos) {
                    int paren_depth = 0;
                    size_t r_pos = lparen;
                    for (; r_pos < text.size(); ++r_pos) {
                        if (text[r_pos] == '(') paren_depth++;
                        else if (text[r_pos] == ')') {
                            paren_depth--;
                            if (paren_depth == 0) break;
                        }
                    }
                    if (r_pos < text.size()) {
                        aln.embedded_tree_newick = text.substr(lparen, r_pos - lparen + 1) + ";";
                    }
                }
            }
        }

        if (as_nucleotide) {
            aln.finalize_nucleotide();
        } else {
            aln.finalize(gcode);
        }
        return aln;
    }

    static Alignment from_fasta_nucleotide(const std::string& filepath) {
        return from_fasta(filepath, nullptr, true);
    }

    static Alignment from_nexus_nucleotide(const std::string& filepath) {
        return from_nexus(filepath, nullptr, true);
    }

    static Alignment load(const std::string& filepath, std::shared_ptr<const GeneticCode> gcode = nullptr) {
        std::string raw = read_file_or_gzip(filepath);
        std::istringstream file(raw);
        std::string line;
        while (std::getline(file, line)) {
            while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front()))) line.erase(line.begin());
            if (line.empty()) continue;
            if (line[0] == '#') {
                return from_nexus(filepath, gcode);
            }
            if (line[0] == '>') {
                return from_fasta(filepath, gcode);
            }
            std::string upper = line;
            for (auto& c : upper) c = std::toupper(static_cast<unsigned char>(c));
            if (upper.find("BEGIN") != std::string::npos) {
                return from_nexus(filepath, gcode);
            }
        }
        return from_fasta(filepath, gcode);
    }

    static Alignment load_nucleotide(const std::string& filepath) {
        std::string raw = read_file_or_gzip(filepath);
        std::istringstream file(raw);
        std::string line;
        while (std::getline(file, line)) {
            while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front()))) line.erase(line.begin());
            if (line.empty()) continue;
            if (line[0] == '#') {
                return from_nexus(filepath, nullptr, true);
            }
            if (line[0] == '>') {
                return from_fasta(filepath, nullptr, true);
            }
            std::string upper = line;
            for (auto& c : upper) c = std::toupper(static_cast<unsigned char>(c));
            if (upper.find("BEGIN") != std::string::npos) {
                return from_nexus(filepath, nullptr, true);
            }
        }
        return from_fasta(filepath, nullptr, true);
    }

    void add_sequence(const std::string& name, const std::string& seq) {
        Sequence s;
        s.name = name;
        s.raw_seq = seq;
        taxon_to_index[name] = taxon_names.size();
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
        lower_taxon_to_index[lower] = taxon_names.size();
        taxon_names.push_back(name);
        sequences.push_back(std::move(s));
    }

    void finalize(std::shared_ptr<const GeneticCode> gcode = nullptr) {
        if (gcode) {
            code = gcode;
        } else if (!code) {
            code = GeneticCode::universal();
        }

        if (sequences.empty()) {
            throw std::runtime_error("Empty alignment");
        }
        num_taxa = sequences.size();
        num_nucleotides = sequences[0].raw_seq.size();

        for (const auto& s : sequences) {
            if (s.raw_seq.size() != num_nucleotides) {
                throw std::runtime_error("Sequences have unequal lengths in alignment");
            }
        }

        if (num_nucleotides % 3 != 0) {
            throw std::runtime_error("Nucleotide alignment length is not a multiple of 3 (in-frame codons required)");
        }
        num_codons = num_nucleotides / 3;

        // Convert raw sequence to codons using genetic code
        for (auto& s : sequences) {
            s.codons.resize(num_codons, -1);
            s.codon_masks.resize(num_codons, 0);
            for (size_t c = 0; c < num_codons; ++c) {
                char c1 = s.raw_seq[c * 3];
                char c2 = s.raw_seq[c * 3 + 1];
                char c3 = s.raw_seq[c * 3 + 2];
                int8_t st = code->triplet_to_sense(c1, c2, c3);
                s.codons[c] = st;
                if (st >= 0) {
                    s.codon_masks[c] = (uint64_t(1) << st);
                } else {
                    s.codon_masks[c] = code->resolve_triplet_mask(c1, c2, c3);
                }
            }
        }

        compute_frequencies(code);
        compress_patterns();
    }

    void finalize_nucleotide() {
        if (sequences.empty()) {
            throw std::runtime_error("Empty alignment");
        }
        num_taxa = sequences.size();
        num_nucleotides = sequences[0].raw_seq.size();

        for (const auto& s : sequences) {
            if (s.raw_seq.size() != num_nucleotides) {
                throw std::runtime_error("Sequences have unequal lengths in alignment");
            }
        }
        num_codons = 0;

        Vector4 nuc_counts = Vector4::Zero();
        for (const auto& s : sequences) {
            for (size_t i = 0; i < num_nucleotides; ++i) {
                int8_t n = char_to_nuc(s.raw_seq[i]);
                if (n >= 0 && n < 4) {
                    nuc_counts(n) += 1.0;
                }
            }
        }
        Scalar total_nucs = nuc_counts.sum();
        if (total_nucs > 0.0) {
            nuc_frequencies = nuc_counts / total_nucs;
        } else {
            nuc_frequencies.setConstant(0.25);
        }

        patterns.clear();
        site_to_pattern.resize(num_nucleotides);
        std::unordered_map<std::string, size_t> pattern_map;

        for (size_t site = 0; site < num_nucleotides; ++site) {
            std::string col(num_taxa, 'N');
            std::vector<int8_t> states(num_taxa, -1);
            std::vector<uint64_t> masks(num_taxa, 0);

            for (size_t t = 0; t < num_taxa; ++t) {
                char nuc = sequences[t].raw_seq[site];
                uint8_t m = char_to_iupac_mask(nuc);
                col[t] = static_cast<char>('A' + m);
                states[t] = char_to_nuc(nuc);
                masks[t] = static_cast<uint64_t>(m);
            }

            auto it = pattern_map.find(col);
            if (it != pattern_map.end()) {
                size_t p_idx = it->second;
                patterns[p_idx].weight++;
                patterns[p_idx].sites.push_back(site);
                site_to_pattern[site] = p_idx;
            } else {
                size_t p_idx = patterns.size();
                pattern_map[col] = p_idx;

                SitePattern p;
                p.states = std::move(states);
                p.ambiguous_masks = std::move(masks);
                p.weight = 1;
                p.sites.push_back(site);
                patterns.push_back(std::move(p));
                site_to_pattern[site] = p_idx;
            }
        }
    }

    Alignment slice_nucleotide(size_t start, size_t end) const {
        if (start >= end || end > num_nucleotides) {
            throw std::runtime_error("Invalid nucleotide slice range: [" +
                                     std::to_string(start) + ", " + std::to_string(end) +
                                     ") for alignment of length " + std::to_string(num_nucleotides));
        }

        Alignment slice;
        slice.taxon_names = taxon_names;
        slice.taxon_to_index = taxon_to_index;
        slice.lower_taxon_to_index = lower_taxon_to_index;
        slice.embedded_tree_newick = embedded_tree_newick;
        slice.code = code;

        size_t len = end - start;
        slice.sequences.resize(num_taxa);
        for (size_t t = 0; t < num_taxa; ++t) {
            slice.sequences[t].name = sequences[t].name;
            slice.sequences[t].raw_seq = sequences[t].raw_seq.substr(start, len);
        }

        slice.finalize_nucleotide();
        return slice;
    }

    Alignment to_nucleotide_alignment() const {
        if (num_codons == 0 && !patterns.empty()) {
            return *this;
        }
        Alignment nuc_aln;
        nuc_aln.taxon_names = taxon_names;
        nuc_aln.taxon_to_index = taxon_to_index;
        nuc_aln.lower_taxon_to_index = lower_taxon_to_index;
        nuc_aln.num_taxa = num_taxa;
        nuc_aln.num_nucleotides = num_nucleotides;
        nuc_aln.num_codons = 0;
        nuc_aln.code = code;
        nuc_aln.embedded_tree_newick = embedded_tree_newick;
        nuc_aln.nuc_frequencies = nuc_frequencies;

        nuc_aln.sequences.resize(num_taxa);
        for (size_t t = 0; t < num_taxa; ++t) {
            nuc_aln.sequences[t].name = sequences[t].name;
            nuc_aln.sequences[t].raw_seq = sequences[t].raw_seq;
        }

        std::unordered_map<std::string, size_t> pattern_map;
        nuc_aln.site_to_pattern.resize(num_nucleotides);

        for (size_t site = 0; site < num_nucleotides; ++site) {
            std::string col(num_taxa, 'N');
            std::vector<int8_t> states(num_taxa, -1);
            std::vector<uint64_t> masks(num_taxa, 0);

            for (size_t t = 0; t < num_taxa; ++t) {
                char nuc = sequences[t].raw_seq[site];
                uint8_t m = char_to_iupac_mask(nuc);
                col[t] = static_cast<char>('A' + m);
                states[t] = char_to_nuc(nuc);
                masks[t] = static_cast<uint64_t>(m);
            }

            auto it = pattern_map.find(col);
            if (it != pattern_map.end()) {
                size_t p_idx = it->second;
                nuc_aln.patterns[p_idx].weight++;
                nuc_aln.patterns[p_idx].sites.push_back(site);
                nuc_aln.site_to_pattern[site] = p_idx;
            } else {
                size_t p_idx = nuc_aln.patterns.size();
                pattern_map[col] = p_idx;

                SitePattern p;
                p.states = std::move(states);
                p.ambiguous_masks = std::move(masks);
                p.weight = 1;
                p.sites.push_back(site);
                nuc_aln.patterns.push_back(std::move(p));
                nuc_aln.site_to_pattern[site] = p_idx;
            }
        }

        return nuc_aln;
    }

private:
    void compute_frequencies(std::shared_ptr<const GeneticCode> gcode = nullptr) {
        if (!gcode) gcode = code ? code : GeneticCode::universal();
        code = gcode;

        Vector4 nuc_counts = Vector4::Zero();
        Eigen::Matrix<Scalar, 3, 4> pos_counts = Eigen::Matrix<Scalar, 3, 4>::Zero();

        for (const auto& s : sequences) {
            for (size_t c = 0; c < num_codons; ++c) {
                for (int pos = 0; pos < 3; ++pos) {
                    int8_t n = char_to_nuc(s.raw_seq[c * 3 + pos]);
                    if (n >= 0) {
                        nuc_counts(n) += 1.0;
                        pos_counts(pos, n) += 1.0;
                    }
                }
            }
        }

        Scalar total_nucs = nuc_counts.sum();
        if (total_nucs > 0.0) {
            nuc_frequencies = nuc_counts / total_nucs;
        }
        bool padded_nuc = false;
        for (int b = 0; b < 4; ++b) {
            if (nuc_frequencies(b) < 1e-8) {
                nuc_frequencies(b) = 1e-8;
                padded_nuc = true;
            }
        }
        if (padded_nuc) {
            nuc_frequencies /= nuc_frequencies.sum();
        }

        for (int pos = 0; pos < 3; ++pos) {
            Scalar pos_total = pos_counts.row(pos).sum();
            if (pos_total > 0.0) {
                pos_nuc_frequencies.row(pos) = pos_counts.row(pos) / pos_total;
            }
            bool padded_pos = false;
            for (int b = 0; b < 4; ++b) {
                if (pos_nuc_frequencies(pos, b) < 1e-8) {
                    pos_nuc_frequencies(pos, b) = 1e-8;
                    padded_pos = true;
                }
            }
            if (padded_pos) {
                pos_nuc_frequencies.row(pos) /= pos_nuc_frequencies.row(pos).sum();
            }
        }

        int S = code->num_sense_codons;
        codon_frequencies_f1x4.resize(S);
        codon_frequencies_f3x4.resize(S);

        // Compute F1x4
        Scalar sum_f1x4 = 0.0;
        for (int idx = 0; idx < S; ++idx) {
            const auto& codon = code->sense_codons[idx];
            int8_t n1 = char_to_nuc(codon[0]);
            int8_t n2 = char_to_nuc(codon[1]);
            int8_t n3 = char_to_nuc(codon[2]);
            codon_frequencies_f1x4(idx) = nuc_frequencies(n1) * nuc_frequencies(n2) * nuc_frequencies(n3);
            sum_f1x4 += codon_frequencies_f1x4(idx);
        }
        if (sum_f1x4 > 0.0) {
            codon_frequencies_f1x4 /= sum_f1x4;
        }

        // Compute Corrected F3x4 (CF3x4) using fixed-point iteration matching HyPhy
        Eigen::Matrix<Scalar, 3, 4> p_cur = pos_nuc_frequencies;
        for (int it = 0; it < 200; ++it) {
            Eigen::Matrix<Scalar, 3, 4> p_prev = p_cur;
            Scalar sum_cf3x4 = 0.0;
            for (int idx = 0; idx < S; ++idx) {
                const auto& codon = code->sense_codons[idx];
                int8_t n1 = char_to_nuc(codon[0]);
                int8_t n2 = char_to_nuc(codon[1]);
                int8_t n3 = char_to_nuc(codon[2]);
                codon_frequencies_f3x4(idx) = p_cur(0, n1) * p_cur(1, n2) * p_cur(2, n3);
                sum_cf3x4 += codon_frequencies_f3x4(idx);
            }
            if (sum_cf3x4 > 0.0) {
                codon_frequencies_f3x4 /= sum_cf3x4;
            }

            Eigen::Matrix<Scalar, 3, 4> q = Eigen::Matrix<Scalar, 3, 4>::Zero();
            for (int idx = 0; idx < S; ++idx) {
                const auto& codon = code->sense_codons[idx];
                int8_t n1 = char_to_nuc(codon[0]);
                int8_t n2 = char_to_nuc(codon[1]);
                int8_t n3 = char_to_nuc(codon[2]);
                q(0, n1) += codon_frequencies_f3x4(idx);
                q(1, n2) += codon_frequencies_f3x4(idx);
                q(2, n3) += codon_frequencies_f3x4(idx);
            }

            for (int pos = 0; pos < 3; ++pos) {
                for (int b = 0; b < 4; ++b) {
                    p_cur(pos, b) *= (pos_nuc_frequencies(pos, b) / std::max(q(pos, b), 1e-12));
                }
                Scalar s = p_cur.row(pos).sum();
                if (s > 0.0) {
                    p_cur.row(pos) /= s;
                }
            }

            if ((p_cur - p_prev).cwiseAbs().maxCoeff() < 1e-14) {
                break;
            }
        }

        // Final CF3x4 normalize
        Scalar sum_cf3x4 = 0.0;
        for (int idx = 0; idx < S; ++idx) {
            const auto& codon = code->sense_codons[idx];
            int8_t n1 = char_to_nuc(codon[0]);
            int8_t n2 = char_to_nuc(codon[1]);
            int8_t n3 = char_to_nuc(codon[2]);
            codon_frequencies_f3x4(idx) = p_cur(0, n1) * p_cur(1, n2) * p_cur(2, n3);
            sum_cf3x4 += codon_frequencies_f3x4(idx);
        }
        if (sum_cf3x4 > 0.0) {
            codon_frequencies_f3x4 /= sum_cf3x4;
        }
        pos_nuc_frequencies = p_cur;
    }

    void compress_patterns() {
        patterns.clear();
        site_to_pattern.resize(num_codons);

        struct PatternKey {
            std::vector<int8_t> states;
            std::vector<uint64_t> masks;
            bool operator==(const PatternKey& o) const {
                return states == o.states && masks == o.masks;
            }
        };

        struct PatternHash {
            size_t operator()(const PatternKey& k) const {
                size_t h = k.states.size();
                for (size_t i = 0; i < k.states.size(); ++i) {
                    size_t v = (k.states[i] >= 0)
                        ? static_cast<size_t>(k.states[i] + 1)
                        : static_cast<size_t>(k.masks[i] ^ (k.masks[i] >> 32));
                    h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
                }
                return h;
            }
        };

        std::unordered_map<PatternKey, size_t, PatternHash> pattern_map;

        PatternKey key;
        key.states.resize(num_taxa);
        key.masks.resize(num_taxa);

        for (size_t c = 0; c < num_codons; ++c) {
            for (size_t t = 0; t < num_taxa; ++t) {
                key.states[t] = sequences[t].codons[c];
                key.masks[t] = (sequences[t].codon_masks.size() > c)
                    ? sequences[t].codon_masks[c]
                    : ((key.states[t] >= 0) ? (uint64_t(1) << key.states[t]) : ~uint64_t(0));
            }

            auto it = pattern_map.find(key);
            if (it == pattern_map.end()) {
                size_t new_idx = patterns.size();
                pattern_map[key] = new_idx;
                SitePattern p;
                p.states = key.states;
                p.ambiguous_masks = key.masks;
                p.weight = 1;
                p.sites.push_back(c);
                patterns.push_back(std::move(p));
                site_to_pattern[c] = new_idx;
            } else {
                size_t p_idx = it->second;
                patterns[p_idx].weight += 1;
                patterns[p_idx].sites.push_back(c);
                site_to_pattern[c] = p_idx;
            }
        }
    }
};

} // namespace hyphy::core
