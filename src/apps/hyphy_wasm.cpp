#include "hyphy/analyses/characterize.hpp"
#include <fstream>
#include <iostream>
#include <string>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define HYPHY3_WASM_KEEPALIVE EMSCRIPTEN_KEEPALIVE
#else
#define HYPHY3_WASM_KEEPALIVE
#endif

namespace {
std::string g_last_html_report;
std::string g_last_json_report;
std::string g_last_error;

bool looks_unaligned_or_needs_stage0(const std::string& raw_text) {
    // Trim leading whitespace
    size_t pos = 0;
    while (pos < raw_text.size() && std::isspace(static_cast<unsigned char>(raw_text[pos]))) {
        ++pos;
    }
    if (pos < raw_text.size() && raw_text[pos] == '#') {
        // NEXUS format is pre-aligned
        return false;
    }
    // Parse FASTA sequence lengths and check if all lengths are equal and divisible by 3
    std::istringstream iss(raw_text);
    std::string line;
    std::string cur_seq;
    size_t first_len = 0;
    bool have_first = false;
    bool has_gaps = false;
    size_t seq_count = 0;

    auto check_seq = [&](const std::string& s) -> bool {
        if (s.empty()) return false;
        seq_count++;
        for (char c : s) {
            if (c == '-') has_gaps = true;
        }
        if (s.size() % 3 != 0) return true; // Needs Stage 0/1 unaligned transducer
        if (!have_first) {
            first_len = s.size();
            have_first = true;
        } else if (s.size() != first_len) {
            return true; // Unequal lengths -> unaligned FASTA
        }
        return false;
    };

    while (std::getline(iss, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur_seq.empty()) {
                if (check_seq(cur_seq)) return true;
                cur_seq.clear();
            }
        } else if (line[0] == '(') {
            break;
        } else {
            cur_seq += line;
        }
    }
    if (!cur_seq.empty()) {
        if (check_seq(cur_seq)) return true;
    }
    // If there are no gaps at all in a multi-sequence FASTA, run full Stage 0/1 joint alignment
    return !has_gaps;
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
        {
            std::ofstream vf(vpath, std::ios::binary | std::ios::trunc);
            vf << raw_seq;
        }

        auto gcode = hyphy::core::GeneticCode::from_name(cname);
        hyphy::analyses::CharacterizeSettings settings;
        settings.quiet = false;

        hyphy::analyses::CharacterizeEngine engine(settings);
        hyphy::analyses::CharacterizeResult res;

        if (looks_unaligned_or_needs_stage0(raw_seq)) {
            res = engine.run_from_unaligned_file(vpath, gcode);
        } else {
            try {
                hyphy::core::Alignment aln = hyphy::core::Alignment::load(vpath, gcode);
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
    }
}

} // extern "C"
