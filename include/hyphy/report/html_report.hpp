#pragma once
// Self-contained HTML reports: the results viewer (viz/fel-viewer, embedded at build time by
// cmake/EmbedViewer.cmake), its stylesheet and the analysis JSON are inlined into one file that
// opens from file:// with no network access, server, Node or Python.

#include <hyphy/report/viewer_assets.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace hyphy::report {

// How the report presents one analysis: the nav badge (HTML) and the viewer element it mounts.
struct ReportViewer {
    std::string method;
    std::string badge;
    std::string tag;
};

// Analyses whose results the embedded viewer renders correctly. Extend as viewers are added.
inline const std::vector<ReportViewer>& report_viewers() {
    static const std::vector<ReportViewer> viewers = {
        {"fel", "FEL &bull; Pervasive Selection", "hyphy-fel-viewer"},
        {"meme", "MEME &bull; Episodic Selection", "hyphy-meme-viewer"},
    };
    return viewers;
}

inline const ReportViewer* find_viewer(const std::string& method) {
    const auto& v = report_viewers();
    auto it = std::find_if(v.begin(), v.end(), [&](const ReportViewer& r) { return r.method == method; });
    return it == v.end() ? nullptr : &*it;
}

inline const std::vector<std::string>& supported_analyses() {
    static const std::vector<std::string> methods = [] {
        std::vector<std::string> m;
        for (const auto& v : report_viewers()) m.push_back(v.method);
        return m;
    }();
    return methods;
}

inline bool has_viewer(const std::string& method) { return find_viewer(method) != nullptr; }

namespace detail {

inline std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

inline std::string replace_all(std::string s, const std::string& from, const std::string& to) {
    for (size_t pos = s.find(from); pos != std::string::npos; pos = s.find(from, pos + to.size())) {
        s.replace(pos, from.size(), to);
    }
    return s;
}

inline std::string html_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            case '\'': out += "&#39;"; break;
            default: out += c;
        }
    }
    return out;
}

// '<' only occurs inside JSON strings, where < is an equivalent escape. This keeps
// sequences such as "</script>" in taxon names from terminating the inline <script> block.
inline std::string json_for_script(const std::string& json_text) {
    return replace_all(json_text, "<", "\\u003c");
}

// "</script" can only appear inside JS string or regex literals, where "<\/script" is equivalent
inline std::string js_for_script(const std::string& js) {
    return replace_all(replace_all(js, "</script", "<\\/script"), "</SCRIPT", "<\\/SCRIPT");
}

inline std::string css_for_style(const std::string& css) {
    return replace_all(css, "</style", "<\\/style");
}

inline std::string basename(const std::string& path) {
    auto slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

}  // namespace detail

// Identify the analysis that produced a result JSON (hyphy3 modern v3 or HyPhy 2.5 legacy).
// Returns a lowercase method id ("fel", "meme", "contrast-fel", ...) or "" when unknown.
inline std::string detect_analysis(const nlohmann::json& j) {
    const auto& analysis = j.contains("analysis") && j["analysis"].is_object() ? j["analysis"] : nlohmann::json::object();

    // Modern v3 output carries an explicit id
    if (analysis.contains("id") && analysis["id"].is_string()) {
        return detail::lower(analysis["id"].get<std::string>());
    }

    // Legacy output: the info text begins with the method name ("FEL (Fixed Effects ...", "GARD : ...")
    if (analysis.contains("info") && analysis["info"].is_string()) {
        std::string info = detail::lower(analysis["info"].get<std::string>());
        size_t start = info.find_first_not_of(" \t\r\n");
        if (start != std::string::npos) {
            size_t end = info.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789-", start);
            std::string lead = info.substr(start, end == std::string::npos ? std::string::npos : end - start);
            static const std::vector<std::string> known = {"fel", "meme", "slac", "fubar", "busted", "absrel",
                                                           "relax", "prime", "gard", "contrast-fel", "bgm", "fade"};
            if (std::find(known.begin(), known.end(), lead) != known.end()) return lead;
        }
        if (info.find("(fubar") != std::string::npos) return "fubar";
        if (start != std::string::npos && info.compare(start, 15, "directional fel") == 0) return "dfel";
    }

    // Older legacy files without usable info: recognize FEL/MEME from their site-table headers
    if (j.contains("MLE") && j["MLE"].contains("headers") && j["MLE"]["headers"].is_array()) {
        bool beta_plus = false, alpha_eq_beta = false;
        for (const auto& h : j["MLE"]["headers"]) {
            std::string name = detail::lower(h.is_array() && !h.empty() && h[0].is_string() ? h[0].get<std::string>()
                                             : h.is_string() ? h.get<std::string>() : "");
            if (name.find("beta") != std::string::npos &&
                (name.find('+') != std::string::npos || name.find("plus") != std::string::npos)) {
                beta_plus = true;
            }
            if (name == "alpha=beta") alpha_eq_beta = true;
        }
        if (beta_plus) return "meme";
        if (alpha_eq_beta) return "fel";
    }
    return "";
}

inline std::string supported_analyses_text() {
    std::string s;
    for (std::string m : supported_analyses()) {
        std::transform(m.begin(), m.end(), m.begin(), [](unsigned char c) { return std::toupper(c); });
        if (!s.empty()) s += ", ";
        s += m;
    }
    return s;
}

// Render a complete standalone HTML report for the given result JSON text.
// Throws std::runtime_error for invalid JSON or analyses without an HTML viewer.
inline std::string render_html(const std::string& json_text, const std::string& source_name) {
    nlohmann::json j;
    try {
        j = nlohmann::json::parse(json_text);
    } catch (const std::exception& e) {
        throw std::runtime_error("Cannot create HTML report: input is not valid JSON (" + std::string(e.what()) + ")");
    }
    const std::string method = detect_analysis(j);
    const ReportViewer* viewer = find_viewer(method);
    if (!viewer) {
        throw std::runtime_error("Cannot create HTML report: no HTML viewer for " +
                                 (method.empty() ? std::string("this result file") : "'" + method + "' results") +
                                 " yet (available for: " + supported_analyses_text() + ")");
    }

    const std::string& badge = viewer->badge;
    const std::string& tag = viewer->tag;
    std::string viewer_tags;
    for (const auto& v : report_viewers()) viewer_tags += (viewer_tags.empty() ? "" : ", ") + v.tag;
    const std::string name = detail::html_escape(detail::basename(source_name));

    std::ostringstream html;
    html << R"(<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>HyPhy )" << badge << " &mdash; " << name << R"(</title>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    html, body {
      width: 100%; height: 100%; overflow: hidden;
      background-color: #F2F2EF; color: #1A1A18;
      font-family: -apple-system, BlinkMacSystemFont, "Helvetica Neue", Helvetica, Arial, sans-serif;
      -webkit-font-smoothing: antialiased;
    }
    #app-container { width: 100vw; height: 100vh; display: flex; flex-direction: column; }
    .app-nav {
      height: 42px; background-color: #FFFFFF; border-bottom: 1px solid #D1D1CB;
      display: flex; align-items: center; justify-content: space-between;
      padding: 0 20px; font-size: 11px; z-index: 100;
    }
    .app-brand {
      display: flex; align-items: center; gap: 10px; font-weight: 700; color: #1A1A18;
      letter-spacing: 1.2px; text-transform: uppercase;
    }
    .app-brand-badge {
      background-color: #1A1A18; color: #FFFFFF; padding: 2px 6px; border-radius: 2px;
      font-size: 10px; font-weight: 700; letter-spacing: 0.8px;
    }
    .app-brand-dot { color: #EB5E28; font-size: 12px; }
    .app-nav-title {
      color: #767672; font-family: "SF Mono", Monaco, monospace; font-size: 11px;
      background-color: #EDEDE9; border: 1px solid #D1D1CB; border-radius: 2px; padding: 2px 8px;
    }
    #viewer-mount {
      flex: 1 1 0%; min-height: 0; display: flex; position: relative;
      overflow-y: auto; overflow-x: hidden;
    }
    )" << viewer_tags << R"( { width: 100%; height: 100%; display: block; }
    /* Viewer stylesheet */
)" << detail::css_for_style(assets::viewer_css()) << R"(
  </style>
</head>
<body>
  <div id="app-container">
    <nav class="app-nav">
      <div class="app-brand">
        <span class="app-brand-badge">HYPHY 3</span>
        <span>)" << badge << R"(</span>
        <span class="app-brand-dot">&#9679;</span>
      </div>
      <div class="app-nav-title"><span>)" << name << R"(</span></div>
    </nav>
    <div id="viewer-mount">
      <)" << tag << R"( id="hyphy-viewer" theme="light"></)" << tag << R"(>
    </div>
  </div>
  <script>
)" << detail::js_for_script(assets::viewer_js()) << R"(
  </script>
  <script>
    const payload = )" << detail::json_for_script(json_text) << R"(;
    function startViewer() {
      const viewer = document.getElementById('hyphy-viewer');
      customElements.whenDefined(')" << tag << R"(').then(function () { viewer.load(payload); });
    }
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', startViewer);
    } else {
      startViewer();
    }
  </script>
</body>
</html>
)";
    return html.str();
}

// Read a result JSON file and write the standalone HTML report next to it (or wherever asked).
inline void write_html_report(const std::string& json_path, const std::string& html_path) {
    std::ifstream in(json_path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot create HTML report: unable to read " + json_path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string html = render_html(buffer.str(), json_path);
    std::ofstream out(html_path, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot create HTML report: unable to write " + html_path);
    out << html;
    if (!out) throw std::runtime_error("Cannot create HTML report: write failed for " + html_path);
}

}  // namespace hyphy::report
