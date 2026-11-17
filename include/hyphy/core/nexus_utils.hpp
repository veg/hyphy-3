#pragma once

#include <algorithm>
#include <cctype>
#include <string>
#include <string_view>
#include <unordered_map>

namespace hyphy::core {

/// Strips NEXUS `[ ... ]` comments (supporting nested brackets) outside single- and
/// double-quoted strings while preserving line breaks and doubled `''` apostrophes.
inline std::string strip_nexus_comments(std::string_view raw) {
    std::string text;
    text.reserve(raw.size());
    int comment_depth = 0;
    char active_quote = '\0';
    for (size_t ri = 0; ri < raw.size(); ++ri) {
        char c = raw[ri];
        if (comment_depth > 0) {
            if (c == '[') {
                comment_depth++;
            } else if (c == ']') {
                comment_depth--;
            } else if (c == '\n' || c == '\r') {
                text += c;
            }
            continue;
        }
        if (active_quote != '\0') {
            text += c;
            if (c == active_quote) {
                if (c == '\'' && ri + 1 < raw.size() && raw[ri + 1] == '\'') {
                    text += raw[++ri];
                } else {
                    active_quote = '\0';
                }
            } else if (c == '\n' || c == '\r') {
                active_quote = '\0';
            }
            continue;
        }
        if (c == '[') {
            comment_depth++;
        } else {
            if (c == '\'' || c == '"') {
                active_quote = c;
            }
            text += c;
        }
    }
    return text;
}

namespace detail {

inline bool is_nexus_delim(char c) {
    return std::isspace(static_cast<unsigned char>(c)) || c == ';';
}

inline size_t find_nexus_keyword_ci(
    std::string_view text,
    std::string_view kw,
    size_t start_pos,
    size_t end_pos
) {
    if (end_pos > text.size()) end_pos = text.size();
    if (start_pos >= end_pos || kw.empty() || end_pos - start_pos < kw.size()) {
        return std::string::npos;
    }

    char active_quote = '\0';
    for (size_t i = start_pos; i + kw.size() <= end_pos; ++i) {
        char c = text[i];
        if (active_quote != '\0') {
            if (c == active_quote) {
                if (c == '\'' && i + 1 < end_pos && text[i + 1] == '\'') {
                    ++i;
                } else {
                    active_quote = '\0';
                }
            }
            continue;
        }
        if (c == '\'' || c == '"') {
            active_quote = c;
            continue;
        }

        if (i > start_pos && !is_nexus_delim(text[i - 1])) {
            continue;
        }
        bool match = true;
        for (size_t j = 0; j < kw.size(); ++j) {
            if (std::toupper(static_cast<unsigned char>(text[i + j])) !=
                std::toupper(static_cast<unsigned char>(kw[j]))) {
                match = false;
                break;
            }
        }
        if (!match) continue;

        size_t after = i + kw.size();
        if (after == end_pos || is_nexus_delim(text[after]) ||
            (kw == "TREE" && text[after] == '*') ||
            (kw == "UTREE" && text[after] == '*')) {
            return i;
        }
    }
    return std::string::npos;
}

inline size_t find_nexus_begin_trees(std::string_view text, size_t start_pos, size_t end_pos, size_t& out_stmt_end) {
    size_t pos = start_pos;
    while (pos < end_pos) {
        size_t b = find_nexus_keyword_ci(text, "BEGIN", pos, end_pos);
        if (b == std::string::npos) break;
        size_t p = b + 5;
        while (p < end_pos && std::isspace(static_cast<unsigned char>(text[p]))) ++p;
        if (p + 5 <= end_pos) {
            const std::string_view trees_kw = "TREES";
            bool ok = true;
            for (size_t j = 0; j < 5; ++j) {
                if (std::toupper(static_cast<unsigned char>(text[p + j])) != trees_kw[j]) {
                    ok = false;
                    break;
                }
            }
            size_t after = p + 5;
            if (ok && (after == end_pos || is_nexus_delim(text[after]))) {
                size_t semi = text.find(';', after);
                out_stmt_end = (semi != std::string::npos && semi < end_pos) ? (semi + 1) : after;
                return b;
            }
        }
        pos = b + 5;
    }
    return std::string::npos;
}

inline size_t find_quoted_end(std::string_view text, size_t quote_pos, size_t limit) {
    char q = text[quote_pos];
    size_t i = quote_pos + 1;
    while (i < limit) {
        if (text[i] == q) {
            if (q == '\'' && i + 1 < limit && text[i + 1] == '\'') {
                i += 2;
                continue;
            }
            return i;
        }
        ++i;
    }
    return limit;
}

inline std::string decode_quoted_token(std::string_view text, size_t quote_pos, size_t& out_next) {
    char q = text[quote_pos];
    std::string val;
    size_t i = quote_pos + 1;
    while (i < text.size()) {
        if (text[i] == q) {
            if (q == '\'' && i + 1 < text.size() && text[i + 1] == '\'') {
                val += '\'';
                i += 2;
                continue;
            }
            ++i;
            break;
        }
        val += text[i++];
    }
    out_next = i;
    return val;
}

inline std::string format_newick_leaf_label(const std::string& label) {
    if (label.empty()) return "''";
    bool needs_quotes = false;
    for (unsigned char uc : label) {
        if (std::isspace(uc) || uc == '(' || uc == ')' || uc == '[' || uc == ']' ||
            uc == '{' || uc == '}' || uc == ',' || uc == ':' || uc == ';' ||
            uc == '\'' || uc == '"') {
            needs_quotes = true;
            break;
        }
    }
    if (!needs_quotes) return label;
    std::string out = "'";
    out.reserve(label.size() + 4);
    for (char c : label) {
        if (c == '\'') out += "''";
        else out += c;
    }
    out += "'";
    return out;
}

inline std::unordered_map<std::string, std::string> parse_nexus_translate_map(
    std::string_view text,
    size_t scope_start,
    size_t scope_end
) {
    std::unordered_map<std::string, std::string> tr_map;
    size_t tr_pos = find_nexus_keyword_ci(text, "TRANSLATE", scope_start, scope_end);
    if (tr_pos == std::string::npos) return tr_map;

    size_t i = tr_pos + 9;
    // Scan TRANSLATE ... ; outside quotes
    size_t semi_pos = scope_end;
    {
        char q = '\0';
        for (size_t k = i; k < scope_end; ++k) {
            char c = text[k];
            if (q != '\0') {
                if (c == q) {
                    if (q == '\'' && k + 1 < scope_end && text[k + 1] == '\'') ++k;
                    else q = '\0';
                }
                continue;
            }
            if (c == '\'' || c == '"') q = c;
            else if (c == ';') {
                semi_pos = k;
                break;
            }
        }
    }

    auto read_tok = [&](std::string& out_tok) -> bool {
        while (i < semi_pos && (std::isspace(static_cast<unsigned char>(text[i])) || text[i] == ',')) {
            ++i;
        }
        if (i >= semi_pos) return false;
        if (text[i] == '\'' || text[i] == '"') {
            size_t next_i = i;
            out_tok = decode_quoted_token(text.substr(0, semi_pos), i, next_i);
            i = next_i;
            return true;
        }
        size_t start = i;
        while (i < semi_pos && !std::isspace(static_cast<unsigned char>(text[i])) && text[i] != ',') {
            ++i;
        }
        out_tok.assign(text.substr(start, i - start));
        return !out_tok.empty();
    };

    std::string key, label;
    while (read_tok(key) && read_tok(label)) {
        tr_map[key] = label;
    }
    return tr_map;
}

inline std::string apply_translate_to_newick(
    std::string_view nwk,
    const std::unordered_map<std::string, std::string>& tr_map
) {
    if (tr_map.empty() || nwk.empty()) return std::string(nwk);

    std::string out;
    out.reserve(nwk.size() + 64);

    bool expect_leaf = false;
    size_t i = 0;
    const size_t n = nwk.size();
    while (i < n) {
        char c = nwk[i];
        if (expect_leaf) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                out += c;
                ++i;
                continue;
            }
            if (c == '(' || c == ',') {
                out += c;
                ++i;
                expect_leaf = true;
                continue;
            }
            if (c == ')' || c == ':' || c == '{' || c == ';') {
                out += c;
                ++i;
                expect_leaf = false;
                continue;
            }

            // Read leaf token (either quoted or unquoted up to Newick delimiter)
            if (c == '\'' || c == '"') {
                size_t tok_start = i;
                size_t next_i = i;
                std::string raw_lbl = decode_quoted_token(nwk, i, next_i);
                auto it = tr_map.find(raw_lbl);
                if (it != tr_map.end()) {
                    out += format_newick_leaf_label(it->second);
                } else {
                    out.append(nwk.substr(tok_start, next_i - tok_start));
                }
                i = next_i;
                expect_leaf = false;
                continue;
            }

            size_t tok_start = i;
            while (i < n) {
                char tc = nwk[i];
                if (std::isspace(static_cast<unsigned char>(tc)) ||
                    tc == '(' || tc == ')' || tc == ':' || tc == ',' ||
                    tc == ';' || tc == '{' || tc == '[') {
                    break;
                }
                ++i;
            }
            std::string raw_tok(nwk.substr(tok_start, i - tok_start));
            auto it = tr_map.find(raw_tok);
            if (it != tr_map.end()) {
                out += format_newick_leaf_label(it->second);
            } else {
                out += raw_tok;
            }
            expect_leaf = false;
            continue;
        }

        // Not expecting a leaf token (may still be inside quotes on an internal label)
        if (c == '\'' || c == '"') {
            size_t q_end = find_quoted_end(nwk, i, n);
            size_t take = (q_end < n) ? (q_end + 1 - i) : (n - i);
            out.append(nwk.substr(i, take));
            i += take;
            continue;
        }
        if (c == '{') {
            size_t close_b = nwk.find('}', i + 1);
            if (close_b != std::string::npos) {
                out.append(nwk.substr(i, close_b - i + 1));
                i = close_b + 1;
                continue;
            }
        }
        out += c;
        if (c == '(' || c == ',') {
            expect_leaf = true;
        }
        ++i;
    }
    return out;
}

} // namespace detail

/// Extracts the first embedded tree Newick string from comment-stripped NEXUS text,
/// preferring the explicit `BEGIN TREES ... END;` block (if present) and resolving
/// any `TRANSLATE` leaf mappings directly onto leaf tokens.
inline std::string extract_nexus_tree_newick(std::string_view clean_text, size_t search_from = 0) {
    if (search_from > clean_text.size()) search_from = clean_text.size();

    size_t scope_start = search_from;
    size_t scope_end = clean_text.size();

    size_t trees_stmt_end = 0;
    size_t trees_begin = detail::find_nexus_begin_trees(clean_text, search_from, clean_text.size(), trees_stmt_end);
    if (trees_begin == std::string::npos && search_from > 0) {
        trees_begin = detail::find_nexus_begin_trees(clean_text, 0, search_from, trees_stmt_end);
    }
    if (trees_begin != std::string::npos) {
        scope_start = trees_stmt_end;
        size_t end_block = detail::find_nexus_keyword_ci(clean_text, "ENDBLOCK", scope_start, clean_text.size());
        size_t end_kw = detail::find_nexus_keyword_ci(clean_text, "END", scope_start, clean_text.size());
        if (end_block != std::string::npos && (end_kw == std::string::npos || end_block < end_kw)) {
            scope_end = end_block;
        } else if (end_kw != std::string::npos) {
            scope_end = end_kw;
        }
    }

    auto tr_map = detail::parse_nexus_translate_map(clean_text, scope_start, scope_end);

    size_t tree_pos = detail::find_nexus_keyword_ci(clean_text, "TREE", scope_start, scope_end);
    size_t utree_pos = detail::find_nexus_keyword_ci(clean_text, "UTREE", scope_start, scope_end);
    if (utree_pos != std::string::npos && (tree_pos == std::string::npos || utree_pos < tree_pos)) {
        tree_pos = utree_pos;
    }
    if (tree_pos == std::string::npos) {
        return "";
    }

    size_t eq_pos = clean_text.find('=', tree_pos);
    if (eq_pos == std::string::npos || eq_pos >= scope_end) {
        return "";
    }
    size_t lparen = clean_text.find('(', eq_pos);
    if (lparen == std::string::npos || lparen >= scope_end) {
        return "";
    }

    int paren_depth = 0;
    size_t r_pos = lparen;
    char quote_c = '\0';
    for (; r_pos < scope_end; ++r_pos) {
        char c = clean_text[r_pos];
        if (quote_c != '\0') {
            if (c == quote_c) {
                if (quote_c == '\'' && r_pos + 1 < scope_end && clean_text[r_pos + 1] == '\'') {
                    ++r_pos;
                } else {
                    quote_c = '\0';
                }
            }
            continue;
        }
        if (c == '\'' || c == '"') {
            quote_c = c;
            continue;
        }
        if (c == '(') {
            paren_depth++;
        } else if (c == ')') {
            paren_depth--;
            if (paren_depth == 0) break;
        }
    }
    if (r_pos >= scope_end || paren_depth != 0) {
        return "";
    }

    std::string raw_nwk;
    size_t semi = std::string::npos;
    for (size_t k = r_pos + 1; k < scope_end && (k - r_pos) <= 64; ++k) {
        char c = clean_text[k];
        if (quote_c != '\0') {
            if (c == quote_c) {
                if (quote_c == '\'' && k + 1 < scope_end && clean_text[k + 1] == '\'') ++k;
                else quote_c = '\0';
            }
            continue;
        }
        if (c == '\'' || c == '"') quote_c = c;
        else if (c == ';') {
            semi = k;
            break;
        }
    }
    if (semi != std::string::npos) {
        raw_nwk.assign(clean_text.substr(lparen, semi - lparen + 1));
    } else {
        raw_nwk.assign(clean_text.substr(lparen, r_pos - lparen + 1));
        raw_nwk += ';';
    }

    return detail::apply_translate_to_newick(raw_nwk, tr_map);
}

inline bool looks_like_nexus_tree_input(std::string_view sv) {
    size_t p = 0;
    while (p < sv.size() && std::isspace(static_cast<unsigned char>(sv[p]))) ++p;
    if (p + 6 <= sv.size()) {
        const std::string_view hdr = "#NEXUS";
        bool hdr_match = true;
        for (size_t j = 0; j < 6; ++j) {
            if (std::toupper(static_cast<unsigned char>(sv[p + j])) != hdr[j]) {
                hdr_match = false;
                break;
            }
        }
        if (hdr_match) return true;
    }
    size_t dummy_end = 0;
    return detail::find_nexus_begin_trees(sv, 0, sv.size(), dummy_end) != std::string::npos;
}

} // namespace hyphy::core
