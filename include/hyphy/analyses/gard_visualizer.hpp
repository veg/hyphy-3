#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <mutex>
#include <atomic>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "hyphy/core/types.hpp"
#include "hyphy/core/console.hpp"

namespace hyphy::analyses {

using namespace hyphy::core;

/**
 * @brief Real-time responsive visualizer and dynamic terminal chart for GARD.
 * 
 * Provides:
 * - 1:1 proportional alignment coordinate ruler with site labels.
 * - Dynamic colored breakpoint track with arrow markers (▲) for partition boundaries.
 * - Live exploration density sparkline showing where candidate models cluster across sites.
 * - Partition segment breakdown with segment lengths in matching colors.
 * - Smooth in-place ANSI terminal updates (~25 FPS) with lock-free atomic counters.
 * - Clean, non-destructive plain-text fallback for non-TTY / pipe / cluster log environments.
 */
class GARDVisualizer {
public:
    enum class Stage {
        None,
        SingleBP,
        MultiBP_GA,
        TopologicalTest
    };

    size_t num_sites = 0;
    size_t num_taxa = 0;
    bool quiet = false;
    bool is_tty = false;
    bool dynamic_chart = true;

    int chart_width = 65;
    std::vector<std::atomic<size_t>> density_bins;
    std::atomic<size_t> max_density{0};

    Stage current_stage = Stage::None;
    size_t current_B = 1;

    std::atomic<size_t> models_evaluated{0};
    size_t total_candidates = 0;
    std::chrono::steady_clock::time_point stage_start_time;
    std::chrono::steady_clock::time_point last_render_time;
    std::chrono::steady_clock::time_point last_log_time;

    mutable std::mutex best_mutex;
    std::vector<size_t> current_best_bps;
    Scalar current_best_caic = std::numeric_limits<Scalar>::infinity();
    Scalar baseline_caic = 0.0;
    std::string ga_status_details;
    std::string ga_status_compact;

    std::mutex render_mutex;
    int lines_rendered = 0;

    static constexpr const char* kSparklineGlyphs[9] = {
        " ", " ", "▂", "▃", "▄", "▅", "▆", "▇", "█"
    };

    GARDVisualizer(size_t sites, size_t taxa, bool q = false, bool force_tty = false, bool enable_chart = true)
        : num_sites(sites), num_taxa(taxa), quiet(q), dynamic_chart(enable_chart) {
        is_tty = force_tty || core::Console::is_terminal();
        int term_w = core::Console::terminal_width(80);
        chart_width = std::clamp(term_w - 14, 30, 75);
        density_bins = std::vector<std::atomic<size_t>>(chart_width);
        for (auto& b : density_bins) b.store(0);
    }

    void reset_density() {
        for (auto& b : density_bins) b.store(0);
        max_density.store(0);
    }

    void record_candidate(const std::vector<size_t>& bps, Scalar /* aicc */) {
        if (num_sites <= 1 || chart_width <= 1) return;
        for (size_t b : bps) {
            int bin = std::clamp(
                static_cast<int>(std::round(static_cast<double>(b) / (num_sites - 1) * (chart_width - 1))),
                0, chart_width - 1
            );
            size_t val = density_bins[bin].fetch_add(1, std::memory_order_relaxed) + 1;
            size_t cur = max_density.load(std::memory_order_relaxed);
            while (val > cur && !max_density.compare_exchange_weak(cur, val, std::memory_order_relaxed));
        }
    }

    void update_best(const std::vector<size_t>& bps, Scalar aicc) {
        bool improved = false;
        {
            std::lock_guard<std::mutex> lock(best_mutex);
            if (aicc < current_best_caic) {
                current_best_caic = aicc;
                current_best_bps = bps;
                std::sort(current_best_bps.begin(), current_best_bps.end());
                improved = true;
            }
        }
        if (improved) {
            if (!is_tty && !quiet) {
                print_non_tty_improvement();
            } else {
                try_render(true);
            }
        }
    }

    void tick(size_t count = 1) {
        models_evaluated.fetch_add(count, std::memory_order_relaxed);
        try_render(false);
    }

    void start_single_bp(size_t candidates, Scalar base_caic) {
        current_stage = Stage::SingleBP;
        current_B = 1;
        total_candidates = candidates;
        baseline_caic = base_caic;
        current_best_caic = std::numeric_limits<Scalar>::infinity();
        current_best_bps.clear();
        models_evaluated.store(0);
        reset_density();
        stage_start_time = std::chrono::steady_clock::now();
        last_render_time = stage_start_time;
        last_log_time = stage_start_time;
        lines_rendered = 0;

        if (!quiet) {
            std::cout << core::Console::accent("◆ [GARD B=1] ")
                      << core::Console::bold("Exhaustive Single Breakpoint Scan ")
                      << core::Console::muted("(" + std::to_string(candidates) + " candidate positions)\n")
                      << std::flush;
            if (is_tty && dynamic_chart) {
                try_render(true);
            }
        }
    }

    void finish_single_bp(const std::vector<size_t>& best_bps, Scalar best_caic, Scalar base_caic) {
        {
            std::lock_guard<std::mutex> lock(best_mutex);
            current_best_bps = best_bps;
            current_best_caic = best_caic;
            baseline_caic = base_caic;
        }

        if (!quiet) {
            if (dynamic_chart) {
                force_render_final();
            }
            lines_rendered = 0;

            if (best_caic < base_caic && !best_bps.empty()) {
                Scalar delta = base_caic - best_caic;
                std::cout << core::Console::success("  ✔ Best single breakpoint: site ")
                          << core::Console::bold(core::Console::brand(std::to_string(best_bps[0] + 1)))
                          << core::Console::success(" • c-AIC: ")
                          << core::Console::bold(core::Console::number(format_scalar(best_caic, 4)))
                          << core::Console::success(" (ΔAICc: ")
                          << core::Console::bold(core::Console::success("-" + format_scalar(delta, 4)))
                          << core::Console::success(")\n\n")
                          << std::flush;
            } else {
                std::cout << core::Console::warning("  ✖ No single breakpoint improved over baseline c-AIC (")
                          << format_scalar(base_caic, 4) << ").\n\n"
                          << std::flush;
            }
        }
    }

    void start_ga(size_t B, size_t pop_size, Scalar prev_best_caic, const std::vector<size_t>& initial_best_bps) {
        current_stage = Stage::MultiBP_GA;
        current_B = B;
        total_candidates = pop_size;
        baseline_caic = prev_best_caic;
        current_best_caic = prev_best_caic;
        current_best_bps = initial_best_bps;
        models_evaluated.store(0);
        reset_density();
        ga_status_details = "Initializing population...";
        ga_status_compact = "Init pop...";
        stage_start_time = std::chrono::steady_clock::now();
        last_render_time = stage_start_time;
        last_log_time = stage_start_time;
        lines_rendered = 0;

        if (!quiet) {
            std::cout << core::Console::accent("◆ [GARD B=" + std::to_string(B) + "] ")
                      << core::Console::bold("CHC Genetic Algorithm ")
                      << core::Console::muted("(searching " + std::to_string(B) + " breakpoints; pop = " + std::to_string(pop_size) + ")\n")
                      << std::flush;
            if (is_tty && dynamic_chart) {
                try_render(true);
            }
        }
    }

    void update_ga_generation(
        size_t generation,
        size_t max_stagnant,
        size_t stagnant,
        Scalar diversity,
        const std::vector<size_t>& best_bps_for_B,
        Scalar best_caic_for_B
    ) {
        std::ostringstream ss;
        ss << "Gen " << generation
           << " • Diversity: " << std::fixed << std::setprecision(3) << diversity
           << " • Stagnant: " << stagnant << "/" << max_stagnant;
        ga_status_details = ss.str();

        std::ostringstream css;
        css << "G" << generation
            << " (stg " << stagnant << "/" << max_stagnant
            << ", div " << std::fixed << std::setprecision(2) << diversity << ")";
        ga_status_compact = css.str();

        update_best(best_bps_for_B, best_caic_for_B);

        if (!is_tty && !quiet) {
            auto now = std::chrono::steady_clock::now();
            auto sec = std::chrono::duration_cast<std::chrono::seconds>(now - last_log_time).count();
            if (sec >= 5 || generation == 1) {
                last_log_time = now;
                std::cout << "[GARD B=" << current_B << "] " << ga_status_details
                          << " | Models: " << models_evaluated.load()
                          << " | Best c-AIC: " << std::fixed << std::setprecision(4) << current_best_caic
                          << "\n" << std::flush;
            }
        } else {
            try_render(true);
        }
    }

    void finish_ga(size_t B, const std::vector<size_t>& best_bps, Scalar best_caic, Scalar prev_caic, bool improved) {
        {
            std::lock_guard<std::mutex> lock(best_mutex);
            current_best_bps = best_bps;
            current_best_caic = best_caic;
            baseline_caic = prev_caic;
        }

        if (!quiet) {
            if (dynamic_chart) {
                force_render_final();
            }
            lines_rendered = 0;

            if (improved) {
                Scalar delta = prev_caic - best_caic;
                std::cout << core::Console::success("  ✔ " + std::to_string(B) + " breakpoints improved c-AIC to ")
                          << core::Console::bold(core::Console::number(format_scalar(best_caic, 4)))
                          << core::Console::success(" (ΔAICc: ")
                          << core::Console::bold(core::Console::success("-" + format_scalar(delta, 4)))
                          << core::Console::success(")\n\n")
                          << std::flush;
            } else {
                std::cout << core::Console::muted("  ℹ " + std::to_string(B) + " breakpoints did not improve c-AIC (")
                          << format_scalar(best_caic, 4) << " vs best " << format_scalar(prev_caic, 4)
                          << "). Concluding breakpoint search.\n\n"
                          << std::flush;
            }
        }
    }

    void try_render(bool force = false) {
        if (quiet || !dynamic_chart || !is_tty) return;

        auto now = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_render_time).count();
        if (!force && elapsed_ms < 40) {
            return;
        }

        if (!render_mutex.try_lock()) {
            return;
        }

        last_render_time = now;
        render_in_place();
        render_mutex.unlock();
    }

    void force_render_final() {
        std::lock_guard<std::mutex> lock(render_mutex);
        render_in_place();
        std::cout << "\n" << std::flush;
    }

    // Helper to safely truncate a styled ANSI line to at most max_w printable columns
    static std::string truncate_ansi_line(std::string_view s, size_t max_w) {
        if (core::Console::printable_len(s) <= max_w) {
            return std::string(s);
        }

        std::string result;
        result.reserve(s.size());
        size_t cur_w = 0;
        size_t i = 0;
        bool has_active_ansi = false;

        while (i < s.size()) {
            if (s[i] == '\033' && i + 1 < s.size() && s[i + 1] == '[') {
                size_t start_esc = i;
                i += 2;
                while (i < s.size() && s[i] != 'm') i++;
                if (i < s.size()) i++;
                std::string_view esc = s.substr(start_esc, i - start_esc);
                result.append(esc);
                if (esc == "\033[0m") {
                    has_active_ansi = false;
                } else {
                    has_active_ansi = true;
                }
            } else if ((static_cast<unsigned char>(s[i]) & 0xC0) == 0x80) {
                // UTF-8 continuation byte
                result.push_back(s[i]);
                i++;
            } else {
                // Start of a character / code point
                if (cur_w >= max_w) {
                    break;
                }
                result.push_back(s[i]);
                cur_w++;
                i++;
            }
        }

        if (has_active_ansi) {
            result.append("\033[0m");
        }
        return result;
    }

    // Formats the partition site intervals with adaptive multi-tier resolution detail
    std::string format_partitions_line(const std::vector<size_t>& bps, int avail_w) const {
        static const auto get_seg_style = [](int seg_idx, const std::string& ch) -> std::string {
            switch (seg_idx % 6) {
                case 0: return core::Console::brand(ch);    // Cyan
                case 1: return core::Console::warning(ch);  // Amber
                case 2: return core::Console::accent(ch);   // Purple
                case 3: return core::Console::success(ch);  // Green
                case 4: return core::Console::info(ch);     // Sky blue
                default: return core::Console::number(ch);  // Gold
            }
        };

        if (bps.empty()) {
            return get_seg_style(0, "[1.." + std::to_string(num_sites) + "] (" + std::to_string(num_sites) + " bp)");
        }

        struct PartInfo {
            size_t start;
            size_t end;
            size_t len;
        };

        std::vector<PartInfo> parts;
        parts.reserve(bps.size() + 1);
        size_t prev_site = 1;
        for (size_t b : bps) {
            size_t bp_site = b + 1;
            if (bp_site >= prev_site) {
                parts.push_back({prev_site, bp_site, bp_site - prev_site + 1});
                prev_site = bp_site + 1;
            }
        }
        if (prev_site <= num_sites) {
            parts.push_back({prev_site, num_sites, num_sites - prev_site + 1});
        }

        if (parts.empty()) {
            return get_seg_style(0, "[1.." + std::to_string(num_sites) + "]");
        }

        size_t min_len = parts[0].len;
        size_t max_len = parts[0].len;
        for (const auto& p : parts) {
            if (p.len < min_len) min_len = p.len;
            if (p.len > max_len) max_len = p.len;
        }

        // --- Tier 1: Full partition intervals with lengths [start..end] (len bp) ---
        {
            std::ostringstream ss;
            for (size_t p = 0; p < parts.size(); ++p) {
                if (p > 0) ss << core::Console::muted(" ── ");
                ss << get_seg_style(static_cast<int>(p),
                    "[" + std::to_string(parts[p].start) + ".." + std::to_string(parts[p].end) + "] (" + std::to_string(parts[p].len) + "bp)");
            }
            std::string t1 = ss.str();
            if (core::Console::printable_len(t1) <= static_cast<size_t>(avail_w)) {
                return t1;
            }
        }

        // --- Tier 2: Compact partition intervals without lengths [start..end] ---
        {
            std::ostringstream ss;
            for (size_t p = 0; p < parts.size(); ++p) {
                if (p > 0) ss << core::Console::muted(" ── ");
                ss << get_seg_style(static_cast<int>(p),
                    "[" + std::to_string(parts[p].start) + ".." + std::to_string(parts[p].end) + "]");
            }
            std::string t2 = ss.str();
            if (core::Console::printable_len(t2) <= static_cast<size_t>(avail_w)) {
                return t2;
            }
        }

        // --- Tier 3A: Breakpoint chain with min/max partition stats ---
        // e.g.: 9 parts (min 43, max 165bp) • 1 ── 86 ── 159 ── ... ── 762
        {
            std::string prefix = std::to_string(parts.size()) + " parts (min " +
                                 std::to_string(min_len) + ", max " + std::to_string(max_len) + "bp) • ";
            std::ostringstream ss;
            ss << core::Console::muted(prefix);
            ss << core::Console::bold("1");
            for (size_t p = 0; p < bps.size(); ++p) {
                ss << get_seg_style(static_cast<int>(p), " ── ");
                ss << core::Console::bold(core::Console::brand(std::to_string(bps[p] + 1)));
            }
            ss << get_seg_style(static_cast<int>(bps.size()), " ── ");
            ss << core::Console::bold(std::to_string(num_sites));
            std::string t3a = ss.str();
            if (core::Console::printable_len(t3a) <= static_cast<size_t>(avail_w)) {
                return t3a;
            }
        }

        // --- Tier 3B: Breakpoint chain with compact prefix ---
        // e.g.: 9 parts • 1 ── 86 ── 159 ── ... ── 762
        {
            std::string prefix = std::to_string(parts.size()) + " parts • ";
            std::ostringstream ss;
            ss << core::Console::muted(prefix);
            ss << core::Console::bold("1");
            for (size_t p = 0; p < bps.size(); ++p) {
                ss << get_seg_style(static_cast<int>(p), " ── ");
                ss << core::Console::bold(core::Console::brand(std::to_string(bps[p] + 1)));
            }
            ss << get_seg_style(static_cast<int>(bps.size()), " ── ");
            ss << core::Console::bold(std::to_string(num_sites));
            std::string t3b = ss.str();
            if (core::Console::printable_len(t3b) <= static_cast<size_t>(avail_w)) {
                return t3b;
            }
        }

        // --- Tier 4: Head & Tail omission summary with partition stats ---
        // e.g.: 10 parts (min 35, max 180bp) • [1..86] ── … (8 more) … ── [641..762]
        {
            std::string prefix = std::to_string(parts.size()) + " parts (min " +
                                 std::to_string(min_len) + ", max " + std::to_string(max_len) + "bp) • ";
            std::string short_prefix = std::to_string(parts.size()) + " parts • ";

            for (const auto& pref : {prefix, short_prefix}) {
                std::ostringstream ss;
                ss << core::Console::muted(pref);
                ss << get_seg_style(0, "[" + std::to_string(parts.front().start) + ".." + std::to_string(parts.front().end) + "]");
                size_t omitted = parts.size() > 2 ? parts.size() - 2 : 0;
                ss << core::Console::muted(" ── … (" + std::to_string(omitted) + " more) … ── ");
                ss << get_seg_style(static_cast<int>(parts.size() - 1),
                    "[" + std::to_string(parts.back().start) + ".." + std::to_string(parts.back().end) + "]");
                std::string t4 = ss.str();
                if (core::Console::printable_len(t4) <= static_cast<size_t>(avail_w)) {
                    return t4;
                }
            }

            // Minimalist fallback: 10 parts • [1..86] … [641..762]
            std::ostringstream ss_min;
            ss_min << core::Console::muted(std::to_string(parts.size()) + " parts • ")
                   << get_seg_style(0, "[" + std::to_string(parts.front().start) + ".." + std::to_string(parts.front().end) + "]")
                   << core::Console::muted(" … ")
                   << get_seg_style(static_cast<int>(parts.size() - 1),
                        "[" + std::to_string(parts.back().start) + ".." + std::to_string(parts.back().end) + "]");
            return ss_min.str();
        }
    }

    // Formats the 5-line dynamic dashboard:
    // Line 0: Status & Speed
    // Line 1: Ruler coordinates
    // Line 2: Breakpoint track
    // Line 3: Density sparkline
    // Line 4: Partition site intervals
    std::vector<std::string> format_dashboard_lines() const {
        std::vector<size_t> bps;
        Scalar caic = 0.0;
        Scalar base_caic = 0.0;
        {
            std::lock_guard<std::mutex> lock(best_mutex);
            bps = current_best_bps;
            caic = current_best_caic;
            base_caic = baseline_caic;
        }

        size_t eval_count = models_evaluated.load(std::memory_order_relaxed);
        auto now = std::chrono::steady_clock::now();
        double elapsed_sec = std::chrono::duration<double>(now - stage_start_time).count();
        double speed = (elapsed_sec > 1e-3) ? (eval_count / elapsed_sec) : 0.0;

        int term_w = core::Console::terminal_width(80);
        int max_line_w = std::max(30, term_w - 1);
        int chart_w = std::clamp(term_w - 14, 25, chart_width);

        // --- Line 0: Status line ---
        std::ostringstream status_ss;
        status_ss << core::Console::bold("Status:  ");
        if (current_stage == Stage::SingleBP) {
            double pct = (total_candidates > 0) ? (100.0 * eval_count / total_candidates) : 0.0;
            status_ss << core::Console::accent("[B=1] ")
                      << eval_count << "/" << total_candidates << " (" << std::fixed << std::setprecision(1) << pct << "%)"
                      << " [" << format_duration(elapsed_sec) << ", " << std::fixed << std::setprecision(0) << speed << " mod/s]";
            if (!std::isinf(caic) && !bps.empty()) {
                Scalar delta = base_caic - caic;
                status_ss << " • Best: " << core::Console::number(format_scalar(caic, 2))
                          << " (site " << core::Console::brand(std::to_string(bps[0] + 1))
                          << ", Δ: " << (delta >= 0 ? core::Console::success("-" + format_scalar(delta, 2)) : format_scalar(delta, 2)) << ")";
            }
        } else if (current_stage == Stage::MultiBP_GA) {
            status_ss << core::Console::accent("[B=" + std::to_string(current_B) + "] ")
                      << eval_count << " models [" << format_duration(elapsed_sec) << ", " << std::fixed << std::setprecision(0) << speed << " mod/s]";
            if (!ga_status_details.empty()) {
                status_ss << " • " << ga_status_details;
            }
            if (!std::isinf(caic)) {
                Scalar delta = base_caic - caic;
                status_ss << " • Best: " << core::Console::number(format_scalar(caic, 2));
                if (delta > 0) {
                    status_ss << " (Δ: " << core::Console::success("-" + format_scalar(delta, 2)) << ")";
                }
            }
        } else {
            status_ss << core::Console::muted("GARD Analysis in progress...");
        }

        std::string status_str = status_ss.str();
        if (core::Console::printable_len(status_str) > static_cast<size_t>(max_line_w)) {
            // Build compact status line
            std::ostringstream css;
            css << core::Console::bold("Status:  ");
            if (current_stage == Stage::SingleBP) {
                double pct = (total_candidates > 0) ? (100.0 * eval_count / total_candidates) : 0.0;
                css << core::Console::accent("[B=1] ")
                    << eval_count << "/" << total_candidates << " (" << std::fixed << std::setprecision(0) << pct << "%)"
                    << " [" << format_duration(elapsed_sec) << "]";
                if (!std::isinf(caic) && !bps.empty()) {
                    Scalar delta = base_caic - caic;
                    css << " • Best: " << core::Console::number(format_scalar(caic, 2))
                        << " (site " << core::Console::brand(std::to_string(bps[0] + 1))
                        << ", Δ: " << (delta >= 0 ? core::Console::success("-" + format_scalar(delta, 1)) : format_scalar(delta, 1)) << ")";
                }
            } else if (current_stage == Stage::MultiBP_GA) {
                css << core::Console::accent("[B=" + std::to_string(current_B) + "] ")
                    << eval_count << " mod [" << format_duration(elapsed_sec) << ", " << std::fixed << std::setprecision(0) << speed << "/s]";
                if (!ga_status_compact.empty()) {
                    css << " • " << ga_status_compact;
                }
                if (!std::isinf(caic)) {
                    Scalar delta = base_caic - caic;
                    css << " • Best: " << core::Console::number(format_scalar(caic, 2));
                    if (delta > 0) {
                        css << " (Δ: " << core::Console::success("-" + format_scalar(delta, 1)) << ")";
                    }
                }
            }
            std::string c_str = css.str();
            if (core::Console::printable_len(c_str) <= static_cast<size_t>(max_line_w)) {
                status_str = c_str;
            }
        }

        // --- Line 1: Ruler line ---
        std::vector<std::string> ruler_chars(chart_w, " ");
        std::vector<bool> ruler_occupied(chart_w, false);

        ruler_chars[0] = core::Console::muted("1");
        ruler_occupied[0] = true;
        if (chart_w > 1) ruler_occupied[1] = true;

        std::string end_str = std::to_string(num_sites) + " bp";
        int end_start = chart_w - static_cast<int>(end_str.length());
        if (end_start > 2) {
            for (int i = 0; i < static_cast<int>(end_str.length()); ++i) {
                if (end_start + i < chart_w && end_start + i >= 0) {
                    ruler_chars[end_start + i] = core::Console::muted(std::string(1, end_str[i]));
                    ruler_occupied[end_start + i] = true;
                }
            }
            if (end_start - 1 >= 0) ruler_occupied[end_start - 1] = true;
        }

        std::vector<int> bp_cols;
        bp_cols.reserve(bps.size());
        for (size_t b : bps) {
            int col = (num_sites > 1)
                ? std::clamp(static_cast<int>(std::round(static_cast<double>(b) / (num_sites - 1) * (chart_w - 1))), 0, chart_w - 1)
                : 0;
            bp_cols.push_back(col);
        }

        // Place breakpoint numbers on ruler with collision avoidance
        for (size_t i = 0; i < bps.size(); ++i) {
            std::string b_str = std::to_string(bps[i] + 1);
            int col = bp_cols[i];
            int b_start = col - static_cast<int>(b_str.length()) / 2;
            int b_end = b_start + static_cast<int>(b_str.length()) - 1;

            if (b_start >= 1 && b_end < chart_w) {
                int check_start = std::max(0, b_start - 1);
                int check_end = std::min(chart_w - 1, b_end + 1);
                bool collision = false;
                for (int c = check_start; c <= check_end; ++c) {
                    if (ruler_occupied[c]) {
                        collision = true;
                        break;
                    }
                }
                if (!collision) {
                    for (size_t j = 0; j < b_str.length(); ++j) {
                        ruler_chars[b_start + j] = core::Console::bold(core::Console::brand(std::string(1, b_str[j])));
                        ruler_occupied[b_start + j] = true;
                    }
                }
            }
        }

        std::ostringstream ruler_ss;
        ruler_ss << core::Console::muted("Ruler:   ");
        for (const auto& ch : ruler_chars) ruler_ss << ch;

        // --- Line 2: Breakpoint track line ---
        static const auto get_seg_style = [](int seg_idx, const std::string& ch) -> std::string {
            switch (seg_idx % 6) {
                case 0: return core::Console::brand(ch);    // Cyan
                case 1: return core::Console::warning(ch);  // Amber
                case 2: return core::Console::accent(ch);   // Purple
                case 3: return core::Console::success(ch);  // Green
                case 4: return core::Console::info(ch);     // Sky blue
                default: return core::Console::number(ch);  // Gold
            }
        };

        std::ostringstream track_ss;
        track_ss << core::Console::border("Track:   ") << core::Console::border("├");

        int current_seg = 0;
        size_t next_bp_idx = 0;

        for (int col = 1; col < chart_w - 1; ++col) {
            if (next_bp_idx < bp_cols.size() && col == bp_cols[next_bp_idx]) {
                track_ss << core::Console::bold(core::Console::danger("▲"));
                current_seg++;
                next_bp_idx++;
                while (next_bp_idx < bp_cols.size() && bp_cols[next_bp_idx] <= col) {
                    next_bp_idx++;
                }
            } else {
                track_ss << get_seg_style(current_seg, "─");
            }
        }
        track_ss << core::Console::border("┤");

        // --- Line 3: Density sparkline ---
        std::ostringstream density_ss;
        density_ss << core::Console::info("Density: ");
        size_t max_d = max_density.load(std::memory_order_relaxed);

        for (int col = 0; col < chart_w; ++col) {
            int bin = (chart_w > 1 && !density_bins.empty())
                ? std::clamp(static_cast<int>(std::round(static_cast<double>(col) / (chart_w - 1) * (density_bins.size() - 1))), 0, static_cast<int>(density_bins.size() - 1))
                : 0;
            size_t val = (bin < static_cast<int>(density_bins.size()))
                ? density_bins[bin].load(std::memory_order_relaxed) : 0;
            if (val == 0 || max_d == 0) {
                density_ss << core::Console::muted("·");
            } else {
                int lvl = std::clamp(static_cast<int>(std::round(static_cast<double>(val) / max_d * 8.0)), 1, 8);
                std::string g = kSparklineGlyphs[lvl];
                if (lvl <= 3) {
                    density_ss << core::Console::info(g);
                } else if (lvl <= 6) {
                    density_ss << core::Console::brand(g);
                } else {
                    density_ss << core::Console::warning(g);
                }
            }
        }

        // --- Line 4: Partition breakdown ---
        std::ostringstream parts_ss;
        parts_ss << core::Console::accent("Parts:   ");
        int avail_parts_w = std::max(20, max_line_w - 9);
        parts_ss << format_partitions_line(bps, avail_parts_w);

        return {
            truncate_ansi_line(status_str, max_line_w),
            truncate_ansi_line(ruler_ss.str(), max_line_w),
            truncate_ansi_line(track_ss.str(), max_line_w),
            truncate_ansi_line(density_ss.str(), max_line_w),
            truncate_ansi_line(parts_ss.str(), max_line_w)
        };
    }

private:
    void render_in_place() {
        auto lines = format_dashboard_lines();
        if (lines_rendered > 0) {
            std::cout << "\033[" << lines_rendered << "A";
        }
        for (const auto& line : lines) {
            std::cout << "\r\033[K" << line << "\n";
        }
        std::cout << std::flush;
        lines_rendered = static_cast<int>(lines.size());
    }

    void print_non_tty_improvement() {
        std::ostringstream ss;
        ss << "  [GARD B=" << current_B << "] Improved c-AIC to "
           << std::fixed << std::setprecision(4) << current_best_caic
           << " at site(s): ";
        for (size_t i = 0; i < current_best_bps.size(); ++i) {
            if (i > 0) ss << ", ";
            ss << (current_best_bps[i] + 1);
        }
        std::cout << ss.str() << "\n" << std::flush;
    }

    static std::string format_scalar(Scalar val, int precision = 4) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(precision) << val;
        return ss.str();
    }

    static std::string format_duration(double sec) {
        int s = static_cast<int>(sec);
        int m = s / 60;
        s %= 60;
        int h = m / 60;
        m %= 60;
        std::ostringstream ss;
        if (h > 0) {
            ss << std::setfill('0') << std::setw(2) << h << ":"
               << std::setfill('0') << std::setw(2) << m << ":"
               << std::setfill('0') << std::setw(2) << s;
        } else {
            ss << std::setfill('0') << std::setw(2) << m << ":"
               << std::setfill('0') << std::setw(2) << s;
        }
        return ss.str();
    }
};

} // namespace hyphy::analyses
