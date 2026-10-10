#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <atomic>
#include <mutex>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#define isatty _isatty
#define fileno _fileno
#define STDOUT_FILENO 1
#else
#include <unistd.h>
#include <sys/ioctl.h>
#endif

namespace hyphy::core {

/**
 * @brief High-performance, thread-safe modern terminal progress bar.
 * 
 * Features:
 * - Sub-character Unicode block smoothing (▏▎▍▌▋▊▉█).
 * - Real-time speed (items/sec), elapsed time, and ETA calculations.
 * - Dynamic status message tagging (e.g. tracking positive/negative sites).
 * - Dynamic terminal width adaptation (compact on <= 80 columns, expanded on wide terminals).
 * - Lock-free rate limiting (~30 FPS) for zero thread contention under OpenMP.
 * - Completely silent when stdout is not a TTY (preserving clean logs/pipes).
 */
class ProgressBar {
public:
    enum class Style {
        Blocks,
        SmoothBlocks,
        Lines
    };

    ProgressBar(
        size_t total,
        std::string task_name = "Processing",
        std::string unit = "it",
        bool force_tty = false,
        Style style = Style::SmoothBlocks
    ) : total_(total),
        task_name_(std::move(task_name)),
        unit_(std::move(unit)),
        force_tty_(force_tty),
        style_(style),
        start_time_(std::chrono::steady_clock::now()),
        last_render_time_(std::chrono::steady_clock::now())
    {
        is_tty_ = force_tty_ || is_terminal();
        if (is_tty_) {
            render(0);
        }
    }

    ~ProgressBar() {
        if (!is_finished_.load()) {
            finish();
        }
    }

    // Advance progress by count (thread-safe, non-blocking for OpenMP)
    void tick(size_t count = 1) {
        if (!is_tty_) return;
        size_t current = completed_.fetch_add(count, std::memory_order_relaxed) + count;
        try_render(current);
    }

    // Set absolute progress
    void update(size_t current) {
        if (!is_tty_) return;
        completed_.store(current, std::memory_order_relaxed);
        try_render(current);
    }

    // Update dynamic status label displayed next to stats
    void set_status(std::string status) {
        {
            std::lock_guard<std::mutex> lock(status_mutex_);
            status_ = std::move(status);
        }
        if (is_tty_) {
            try_render(completed_.load(std::memory_order_relaxed));
        }
    }

    // Conclude progress bar with an optional completion summary line
    void finish(const std::string& summary_message = "") {
        if (is_finished_.exchange(true)) return;
        if (!is_tty_) return;

        std::lock_guard<std::mutex> lock(render_mutex_);
        auto now = std::chrono::steady_clock::now();
        double total_sec = std::chrono::duration<double>(now - start_time_).count();
        double speed = (total_sec > 1e-4) ? (completed_.load() / total_sec) : 0.0;

        // Clear active progress bar line
        std::cout << "\r\033[K";
        
        // Print clean checkmark completion line
        std::cout << "\033[1;32m✓\033[0m "
                  << "\033[1;36m" << task_name_ << "\033[0m: "
                  << completed_.load() << "/" << total_ << " " << unit_ << " in "
                  << std::fixed << std::setprecision(2) << total_sec << "s "
                  << "(" << std::setprecision(1) << speed << " " << unit_ << "/s)";
        
        if (!summary_message.empty()) {
            std::cout << " • " << summary_message;
        } else {
            std::string stat;
            {
                std::lock_guard<std::mutex> slock(status_mutex_);
                stat = status_;
            }
            if (!stat.empty()) {
                std::cout << " • " << stat;
            }
        }
        std::cout << "\n" << std::flush;
    }

    static bool is_terminal() {
#if defined(_WIN32)
        return isatty(STDOUT_FILENO) != 0;
#else
        const char* term = std::getenv("TERM");
        if (term && std::string(term) == "dumb") return false;
        return isatty(STDOUT_FILENO) != 0;
#endif
    }

    static int terminal_width() {
#if defined(_WIN32)
        return 80;
#else
        if (!is_terminal()) return 80;
        struct winsize w;
        if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 20) {
            return static_cast<int>(w.ws_col);
        }
        return 80;
#endif
    }

private:
    void try_render(size_t current) {
        if (is_finished_.load()) return;

        auto now = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_render_time_).count();
        // Limit rendering to ~30 frames/sec unless complete
        if (current < total_ && elapsed_ms < 33) {
            return;
        }

        // Try lock ensures worker threads never block scientific computations
        if (!render_mutex_.try_lock()) {
            return;
        }

        last_render_time_ = now;
        render(current);
        render_mutex_.unlock();
    }

    void render(size_t current) {
        if (total_ == 0) return;
        current = std::min(current, total_);

        double fraction = static_cast<double>(current) / static_cast<double>(total_);
        fraction = std::clamp(fraction, 0.0, 1.0);
        double percent = fraction * 100.0;

        auto now = std::chrono::steady_clock::now();
        double elapsed_sec = std::chrono::duration<double>(now - start_time_).count();
        double speed = (elapsed_sec > 1e-3) ? (current / elapsed_sec) : 0.0;
        double remaining_sec = (speed > 1e-4) ? ((total_ - current) / speed) : 0.0;

        int term_w = terminal_width();
        bool compact = (term_w < 105);

        // Header / Task name
        std::string display_name = task_name_;
        if (compact && display_name.length() > 16) {
            display_name = display_name.substr(0, 16);
        }

        // Status text
        std::string stat;
        {
            std::lock_guard<std::mutex> slock(status_mutex_);
            stat = status_;
        }

        std::ostringstream stats_ss;
        stats_ss << " " << std::setw(5) << std::fixed << std::setprecision(1) << percent << "%"
                 << " [" << current << "/" << total_ << "]";

        if (!compact) {
            std::string time_str = format_time(elapsed_sec) + "<" + format_time(remaining_sec);
            stats_ss << " [" << time_str << ", " << std::setprecision(1) << speed << " " << unit_ << "/s]";
        } else {
            stats_ss << " [" << format_time(elapsed_sec) << ", " << std::setprecision(0) << speed << "/s]";
        }

        if (!stat.empty()) {
            stats_ss << " (" << stat << ")";
        }

        // Measure plain printable lengths (excluding ANSI escape codes)
        std::string stats_raw = stats_ss.str();
        size_t ansi_chars = 0;
        for (size_t i = 0; i + 1 < stats_raw.size(); ++i) {
            if (stats_raw[i] == '\033' && stats_raw[i+1] == '[') {
                size_t j = i + 2;
                while (j < stats_raw.size() && stats_raw[j] != 'm') j++;
                ansi_chars += (j - i + 1);
                i = j;
            }
        }
        int stats_printable_len = static_cast<int>(stats_raw.length() - ansi_chars);
        int prefix_printable_len = static_cast<int>(display_name.length()) + 1; // name + space

        int bar_width = term_w - prefix_printable_len - stats_printable_len - 2;
        if (bar_width < 8) bar_width = 8;
        if (bar_width > 35) bar_width = 35;

        std::ostringstream out;
        out << "\r\033[K";

        // Task name
        out << "\033[1;36m" << display_name << "\033[0m ";

        // Render bar
        if (style_ == Style::Lines) {
            int filled = static_cast<int>(fraction * bar_width);
            out << "\033[38;5;45m";
            for (int i = 0; i < filled; ++i) out << "━";
            if (filled < bar_width) {
                out << "\033[1;37m╸\033[38;5;238m";
                for (int i = filled + 1; i < bar_width; ++i) out << "━";
            }
            out << "\033[0m";
        } else if (style_ == Style::Blocks) {
            int filled = static_cast<int>(fraction * bar_width);
            out << "\033[38;5;48m";
            for (int i = 0; i < filled; ++i) out << "█";
            out << "\033[38;5;238m";
            for (int i = filled; i < bar_width; ++i) out << "░";
            out << "\033[0m";
        } else {
            // Style::SmoothBlocks (default)
            double total_blocks = fraction * bar_width;
            int full_blocks = static_cast<int>(total_blocks);
            double partial = total_blocks - full_blocks;
            int partial_idx = static_cast<int>(partial * 8.0);

            const char* sub_blocks[] = {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"};

            out << "\033[38;5;45m";
            for (int i = 0; i < full_blocks; ++i) {
                out << "█";
            }
            if (full_blocks < bar_width) {
                if (partial_idx > 0) {
                    out << sub_blocks[partial_idx];
                } else {
                    out << "\033[38;5;238m░\033[38;5;45m";
                }
                out << "\033[38;5;238m";
                for (int i = full_blocks + 1; i < bar_width; ++i) {
                    out << "░";
                }
            }
            out << "\033[0m";
        }

        // Stats text
        out << "\033[1;37m" << stats_raw << "\033[0m";
        std::cout << out.str() << std::flush;
    }

    std::string format_time(double sec) {
        if (std::isnan(sec) || std::isinf(sec) || sec < 0) sec = 0;
        int s = static_cast<int>(sec);
        int m = s / 60;
        s = s % 60;
        std::ostringstream ss;
        ss << std::setfill('0') << std::setw(2) << m << ":" << std::setw(2) << s;
        return ss.str();
    }

    size_t total_ = 0;
    std::string task_name_;
    std::string unit_;
    bool force_tty_ = false;
    Style style_ = Style::SmoothBlocks;
    bool is_tty_ = false;

    std::atomic<size_t> completed_{0};
    std::atomic<bool> is_finished_{false};

    std::chrono::steady_clock::time_point start_time_;
    std::chrono::steady_clock::time_point last_render_time_;

    std::mutex render_mutex_;
    std::string status_;
    std::mutex status_mutex_;
};

} // namespace hyphy::core
