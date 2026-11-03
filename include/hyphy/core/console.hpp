#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <utility>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#if defined(_WIN32)
#include <io.h>
#define isatty _isatty
#define STDOUT_FILENO 1
#else
#include <unistd.h>
#include <sys/ioctl.h>
#endif

namespace hyphy::core {

class Console {
public:
    static bool is_terminal() {
#if defined(_WIN32)
        return isatty(STDOUT_FILENO) != 0;
#else
        const char* term = std::getenv("TERM");
        if (term && std::string(term) == "dumb") return false;
        return isatty(STDOUT_FILENO) != 0;
#endif
    }

    static bool color_enabled() {
        if (!is_terminal()) return false;
        if (std::getenv("NO_COLOR") != nullptr) return false;
        return true;
    }

    static int terminal_width(int default_w = 80) {
#if defined(_WIN32)
        return default_w;
#else
        if (!is_terminal()) return default_w;
        struct winsize w;
        if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col >= 40) {
            return static_cast<int>(w.ws_col);
        }
        return default_w;
#endif
    }

    // Strip ANSI codes to measure exact printable character length
    static size_t printable_len(std::string_view s) {
        size_t len = 0;
        size_t i = 0;
        while (i < s.size()) {
            if (s[i] == '\033' && i + 1 < s.size() && s[i + 1] == '[') {
                i += 2;
                while (i < s.size() && s[i] != 'm') i++;
                if (i < s.size()) i++;
            } else if ((static_cast<unsigned char>(s[i]) & 0xC0) == 0x80) {
                // UTF-8 continuation byte: skip without incrementing printable length
                i++;
            } else {
                len++;
                i++;
            }
        }
        return len;
    }

    // Styling primitives
    static std::string wrap(std::string_view code, std::string_view text) {
        if (!color_enabled() || text.empty()) return std::string(text);
        return "\033[" + std::string(code) + "m" + std::string(text) + "\033[0m";
    }

    // Vibrant 256-color palette
    static std::string bold(std::string_view text)    { return wrap("1", text); }
    static std::string dim(std::string_view text)     { return wrap("2", text); }
    static std::string italic(std::string_view text)  { return wrap("3", text); }
    static std::string underline(std::string_view text){ return wrap("4", text); }

    static std::string brand(std::string_view text)   { return wrap("1;38;5;45", text); }   // Bold Electric Cyan
    static std::string accent(std::string_view text)  { return wrap("1;38;5;141", text); }  // Bold Soft Purple
    static std::string success(std::string_view text) { return wrap("1;38;5;48", text); }   // Bold Emerald Green
    static std::string warning(std::string_view text) { return wrap("1;38;5;214", text); }  // Bold Warm Amber
    static std::string danger(std::string_view text)  { return wrap("1;38;5;203", text); }  // Bold Coral Red
    static std::string info(std::string_view text)    { return wrap("1;38;5;75", text); }   // Bold Sky Blue
    static std::string muted(std::string_view text)   { return wrap("38;5;244", text); }    // Slate Grey
    static std::string border(std::string_view text)  { return wrap("38;5;239", text); }    // Border Dark Slate
    static std::string number(std::string_view text)  { return wrap("38;5;222", text); }    // Warm Gold Number
    static std::string badge(std::string_view text, std::string_view color_code = "48;5;24;38;5;231;1") {
        if (!color_enabled()) return "[" + std::string(text) + "]";
        return "\033[" + std::string(color_code) + "m " + std::string(text) + " \033[0m";
    }
};

// Elegant Unicode Box & Panel Renderers
class Panel {
public:
    static void print_banner(
        const std::string& tool_badge,
        const std::string& tool_title,
        const std::string& subtitle,
        const std::string& citation = "",
        int width = 76
    ) {
        std::string title_line = " " + Console::brand("⬢ " + tool_badge) + " " +
                                 Console::bold(tool_title);
        std::string sub_line = subtitle.empty() ? "" : (" " + Console::info(subtitle));
        std::string cit_line = citation.empty() ? "" : (" " + Console::muted(citation));

        size_t max_content = Console::printable_len(title_line);
        if (!sub_line.empty()) max_content = std::max(max_content, Console::printable_len(sub_line));
        if (!cit_line.empty()) max_content = std::max(max_content, Console::printable_len(cit_line));

        int w = std::min(width, Console::terminal_width(width));
        if (w < 40) w = 40;
        w = std::max(w, static_cast<int>(max_content) + 4);
        int inner_w = w - 4;

        std::string top = Console::border("╭");
        for (int i = 0; i < w - 2; ++i) top += Console::border("─");
        top += Console::border("╮");

        std::string bot = Console::border("╰");
        for (int i = 0; i < w - 2; ++i) bot += Console::border("─");
        bot += Console::border("╯");

        std::cout << "\n" << top << "\n";

        print_panel_line(title_line, inner_w);
        if (!sub_line.empty()) {
            print_panel_line(sub_line, inner_w);
        }
        if (!cit_line.empty()) {
            print_panel_line(cit_line, inner_w);
        }

        std::cout << bot << "\n\n" << std::flush;
    }

    static void print_step(
        int current,
        int total,
        const std::string& title,
        const std::string& subtitle = ""
    ) {
        std::cout << Console::accent("◆ [") 
                  << Console::brand(std::to_string(current) + "/" + std::to_string(total)) 
                  << Console::accent("] ")
                  << Console::bold(title);
        if (!subtitle.empty()) {
            std::cout << Console::muted(" (" + subtitle + ")");
        }
        std::cout << "\n" << std::flush;
    }

    static void print_card(
        const std::string& title,
        const std::vector<std::pair<std::string, std::string>>& items,
        int width = 76
    ) {
        int term_w = Console::terminal_width(100);
        int max_w = std::max(width, std::min(term_w, 100));

        // Find max printable key width
        size_t max_key = 0;
        for (const auto& item : items) {
            max_key = std::max(max_key, Console::printable_len(item.first));
        }

        // First pass: check natural max width before wrapping
        size_t natural_max = 0;
        for (const auto& item : items) {
            size_t line_len = 1 + max_key + 3 + Console::printable_len(item.second);
            natural_max = std::max(natural_max, line_len);
        }

        std::string header_title = " " + title + " ";
        int header_len = static_cast<int>(Console::printable_len(header_title));
        int w = std::max({width, std::min(static_cast<int>(natural_max) + 4, max_w), header_len + 6});
        if (w < 40) w = 40;
        int inner_w = w - 4;
        int val_col_start = 1 + static_cast<int>(max_key) + 3; // " " + key + " : "
        int val_max_w = std::max(16, inner_w - val_col_start);

        std::vector<std::string> formatted_lines;
        size_t actual_max = 0;
        for (const auto& item : items) {
            size_t key_len = Console::printable_len(item.first);
            size_t pad = (max_key > key_len) ? (max_key - key_len) : 0;
            std::vector<std::string> wrapped_vals = wrap_text(item.second, val_max_w);
            for (size_t i = 0; i < wrapped_vals.size(); ++i) {
                std::string line;
                if (i == 0) {
                    line = " " + item.first + std::string(pad, ' ') +
                           Console::border(" : ") + wrapped_vals[i];
                } else {
                    line = std::string(val_col_start, ' ') + wrapped_vals[i];
                }
                actual_max = std::max(actual_max, Console::printable_len(line));
                formatted_lines.push_back(std::move(line));
            }
        }

        w = std::max({w, static_cast<int>(actual_max) + 4, header_len + 6});
        inner_w = w - 4;

        // Top line with title embedded
        int line_left = 2;
        int line_right = w - 2 - line_left - header_len;
        if (line_right < 2) line_right = 2;

        std::string top = Console::border("┌");
        for (int i = 0; i < line_left; ++i) top += Console::border("─");
        top += Console::bold(header_title);
        for (int i = 0; i < line_right; ++i) top += Console::border("─");
        top += Console::border("┐");

        std::string bot = Console::border("└");
        for (int i = 0; i < w - 2; ++i) bot += Console::border("─");
        bot += Console::border("┘");

        std::cout << top << "\n";

        for (const auto& line : formatted_lines) {
            print_box_line(line, inner_w);
        }

        std::cout << bot << "\n\n" << std::flush;
    }

    static void print_summary_card(
        const std::string& title,
        const std::vector<std::pair<std::string, std::string>>& items,
        const std::string& conclusion = "",
        bool is_success = true,
        int width = 76
    ) {
        int term_w = Console::terminal_width(100);
        int max_w = std::max(width, std::min(term_w, 100));

        size_t max_key = 0;
        for (const auto& item : items) {
            max_key = std::max(max_key, Console::printable_len(item.first));
        }

        size_t natural_max = 0;
        for (const auto& item : items) {
            size_t line_len = 1 + max_key + 3 + Console::printable_len(item.second);
            natural_max = std::max(natural_max, line_len);
        }

        std::string header_title = " " + title + " ";
        int header_len = static_cast<int>(Console::printable_len(header_title));
        int w = std::max({width, std::min(static_cast<int>(natural_max) + 4, max_w), header_len + 6});
        if (w < 40) w = 40;
        int inner_w = w - 4;
        int val_col_start = 1 + static_cast<int>(max_key) + 3;
        int val_max_w = std::max(16, inner_w - val_col_start);

        std::vector<std::string> formatted_lines;
        size_t actual_max = 0;
        for (const auto& item : items) {
            size_t key_len = Console::printable_len(item.first);
            size_t pad = (max_key > key_len) ? (max_key - key_len) : 0;
            std::vector<std::string> wrapped_vals = wrap_text(item.second, val_max_w);
            for (size_t i = 0; i < wrapped_vals.size(); ++i) {
                std::string line;
                if (i == 0) {
                    line = " " + item.first + std::string(pad, ' ') +
                           Console::border(" : ") + wrapped_vals[i];
                } else {
                    line = std::string(val_col_start, ' ') + wrapped_vals[i];
                }
                actual_max = std::max(actual_max, Console::printable_len(line));
                formatted_lines.push_back(std::move(line));
            }
        }

        w = std::max({w, static_cast<int>(actual_max) + 4, header_len + 6});
        inner_w = w - 4;

        std::vector<std::string> conclusion_lines;
        if (!conclusion.empty()) {
            int conc_max_w = std::max(20, inner_w - 3); // " ✦ " is 3 printable chars
            std::vector<std::string> wrapped_conc = wrap_text(conclusion, conc_max_w);
            for (size_t i = 0; i < wrapped_conc.size(); ++i) {
                std::string prefix = (i == 0) ? "✦ " : "  ";
                std::string c_line = " " + (is_success ? Console::success(prefix + wrapped_conc[i])
                                                       : Console::warning(prefix + wrapped_conc[i]));
                actual_max = std::max(actual_max, Console::printable_len(c_line));
                conclusion_lines.push_back(std::move(c_line));
            }
            w = std::max(w, static_cast<int>(actual_max) + 4);
            inner_w = w - 4;
        }

        int line_left = 2;
        int line_right = w - 2 - line_left - header_len;
        if (line_right < 2) line_right = 2;

        std::string top = Console::border("╭");
        for (int i = 0; i < line_left; ++i) top += Console::border("─");
        top += Console::brand(header_title);
        for (int i = 0; i < line_right; ++i) top += Console::border("─");
        top += Console::border("╮");

        std::string bot = Console::border("╰");
        for (int i = 0; i < w - 2; ++i) bot += Console::border("─");
        bot += Console::border("╯");

        std::cout << "\n" << top << "\n";

        for (const auto& line : formatted_lines) {
            print_panel_line(line, inner_w);
        }

        if (!conclusion_lines.empty()) {
            std::cout << Console::border("│") << std::string(w - 2, ' ') << Console::border("│") << "\n";
            for (const auto& c_line : conclusion_lines) {
                print_panel_line(c_line, inner_w);
            }
        }

        std::cout << bot << "\n\n" << std::flush;
    }

private:
    static std::vector<std::string> wrap_text(const std::string& text, int max_width) {
        if (max_width <= 0 || static_cast<int>(Console::printable_len(text)) <= max_width) {
            return {text};
        }
        bool has_ansi = (text.find('\033') != std::string::npos);
        std::vector<std::string> lines;
        std::istringstream words(text);
        std::string word;
        std::string current;

        auto push_long_token = [&](const std::string& token) {
            if (has_ansi || static_cast<int>(Console::printable_len(token)) <= max_width) {
                current = token;
                return;
            }
            std::string rem = token;
            while (static_cast<int>(Console::printable_len(rem)) > max_width) {
                // Prefer splitting after '/' within [max_width/2, max_width]
                size_t split_pos = std::string::npos;
                size_t search_limit = std::min(rem.size(), static_cast<size_t>(max_width));
                size_t slash = rem.rfind('/', search_limit - 1);
                if (slash != std::string::npos && slash >= static_cast<size_t>(max_width / 3)) {
                    split_pos = slash + 1;
                } else {
                    split_pos = search_limit;
                }
                lines.push_back(rem.substr(0, split_pos));
                rem = rem.substr(split_pos);
            }
            current = rem;
        };

        while (words >> word) {
            if (current.empty()) {
                push_long_token(word);
            } else if (static_cast<int>(Console::printable_len(current) + 1 + Console::printable_len(word)) <= max_width) {
                current += " " + word;
            } else {
                lines.push_back(current);
                current.clear();
                push_long_token(word);
            }
        }
        if (!current.empty()) {
            lines.push_back(current);
        }
        if (lines.empty()) {
            lines.push_back(text);
        }
        return lines;
    }

    static void print_panel_line(const std::string& content, int inner_w) {
        size_t p_len = Console::printable_len(content);
        std::cout << Console::border("│ ") << content;
        if (static_cast<int>(p_len) < inner_w) {
            std::cout << std::string(inner_w - p_len, ' ');
        }
        std::cout << Console::border(" │\n");
    }

    static void print_box_line(const std::string& content, int inner_w) {
        size_t p_len = Console::printable_len(content);
        std::cout << Console::border("│ ") << content;
        if (static_cast<int>(p_len) < inner_w) {
            std::cout << std::string(inner_w - p_len, ' ');
        }
        std::cout << Console::border(" │\n");
    }
};

// Modern Unicode Table Generator
class Table {
public:
    enum class Align { Left, Right, Center };

    struct Column {
        std::string header;
        Align align = Align::Left;
        int min_width = 0;
        int computed_width = 0;
    };

    struct Row {
        std::vector<std::string> cells;
        std::string highlight_color; // Optional ANSI escape prefix
    };

    void add_column(std::string header, Align align = Align::Left, int min_width = 0) {
        columns_.push_back({std::move(header), align, min_width, 0});
    }

    void add_row(std::vector<std::string> cells, std::string highlight_color = "") {
        rows_.push_back({std::move(cells), std::move(highlight_color)});
    }

    size_t num_rows() const {
        return rows_.size();
    }

    void print() const {
        if (columns_.empty()) return;

        // 1. Compute widths
        std::vector<int> col_w(columns_.size(), 0);
        for (size_t c = 0; c < columns_.size(); ++c) {
            col_w[c] = std::max(static_cast<int>(Console::printable_len(columns_[c].header)), columns_[c].min_width);
        }

        for (const auto& row : rows_) {
            for (size_t c = 0; c < std::min(row.cells.size(), columns_.size()); ++c) {
                col_w[c] = std::max(col_w[c], static_cast<int>(Console::printable_len(row.cells[c])));
            }
        }

        // 2. Render borders
        auto render_horizontal = [&](const char* left, const char* mid, const char* cross, const char* right) {
            std::ostringstream ss;
            ss << left;
            for (size_t c = 0; c < col_w.size(); ++c) {
                for (int i = 0; i < col_w[c] + 2; ++i) ss << mid;
                if (c + 1 < col_w.size()) ss << cross;
            }
            ss << right << "\n";
            return Console::border(ss.str());
        };

        std::cout << render_horizontal("┌", "─", "┬", "┐");

        // 3. Render Header
        std::cout << Console::border("│");
        for (size_t c = 0; c < columns_.size(); ++c) {
            std::cout << " " << format_cell(columns_[c].header, col_w[c], Align::Center, true) << " " << Console::border("│");
        }
        std::cout << "\n";

        std::cout << render_horizontal("├", "─", "┼", "┤");

        // 4. Render Rows
        for (const auto& row : rows_) {
            std::cout << Console::border("│");
            for (size_t c = 0; c < columns_.size(); ++c) {
                std::string cell_val = (c < row.cells.size()) ? row.cells[c] : "";
                std::string formatted = format_cell(cell_val, col_w[c], columns_[c].align, false);
                if (!row.highlight_color.empty() && Console::color_enabled()) {
                    std::cout << " " << row.highlight_color << formatted << "\033[0m " << Console::border("│");
                } else {
                    std::cout << " " << formatted << " " << Console::border("│");
                }
            }
            std::cout << "\n";
        }

        std::cout << render_horizontal("└", "─", "┴", "┘") << std::flush;
    }

private:
    static std::string format_cell(const std::string& text, int width, Align align, bool is_header) {
        int text_len = static_cast<int>(Console::printable_len(text));
        int pad = std::max(0, width - text_len);

        std::ostringstream ss;
        if (align == Align::Right) {
            ss << std::string(pad, ' ') << (is_header ? Console::bold(text) : text);
        } else if (align == Align::Center) {
            int left_pad = pad / 2;
            int right_pad = pad - left_pad;
            ss << std::string(left_pad, ' ') << (is_header ? Console::bold(text) : text) << std::string(right_pad, ' ');
        } else {
            ss << (is_header ? Console::bold(text) : text) << std::string(pad, ' ');
        }
        return ss.str();
    }

    std::vector<Column> columns_;
    std::vector<Row> rows_;
};

} // namespace hyphy::core
