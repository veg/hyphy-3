#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/report/html_report.hpp"

#include <stdexcept>
#include <string>

using hyphy::report::detect_analysis;
using hyphy::report::render_html;
using nlohmann::json;

static size_t count(const std::string& haystack, const std::string& needle) {
    size_t n = 0;
    for (size_t pos = haystack.find(needle); pos != std::string::npos; pos = haystack.find(needle, pos + 1)) ++n;
    return n;
}

TEST_CASE("Report: analysis detection for modern and legacy result JSON") {
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"id": "fel"}})J")) == "fel");
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"id": "MEME"}})J")) == "meme");
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"id": "contrast-fel"}})J")) == "contrast-fel");

    // Legacy HyPhy 2.5 files name the method at the start of analysis.info
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"info": "FEL (Fixed Effects Likelihood)\n estimates"}})J")) == "fel");
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"info": "MEME (Mixed Effects Model of Evolution) estimates"}})J")) == "meme");
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"info": "SLAC (Single Likelihood Ancestor Counting)"}})J")) == "slac");
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"info": "GARD : Genetic Algorithms for Recombination"}})J")) == "gard");
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"info": "Perform a Fast Unbiased AppRoximate Bayesian (FUBAR) analysis"}})J")) == "fubar");
    CHECK(detect_analysis(json::parse(R"J({"analysis": {"info": "Directional FEL detects site-level selection"}})J")) == "dfel");

    // Without usable info, FEL and MEME are recognized from their site-table headers
    CHECK(detect_analysis(json::parse(R"J({"MLE": {"headers": [["alpha", ""], ["beta", ""], ["alpha=beta", ""]]}})J")) == "fel");
    CHECK(detect_analysis(json::parse(R"J({"MLE": {"headers": [["&alpha;", ""], ["&beta;<sup>+</sup>", ""]]}})J")) == "meme");

    CHECK(detect_analysis(json::parse(R"J({"something": "else"})J")) == "");
}

TEST_CASE("Report: analyses without a viewer are refused instead of rendered as FEL") {
    CHECK_THROWS_AS(render_html(R"J({"analysis": {"id": "busted"}})J", "x.json"), std::runtime_error);
    CHECK_THROWS_AS(render_html(R"J({"analysis": {}})J", "x.json"), std::runtime_error);
    CHECK_THROWS_AS(render_html("not json", "x.json"), std::runtime_error);
    try {
        render_html(R"J({"analysis": {"id": "busted"}})J", "x.json");
    } catch (const std::runtime_error& e) {
        CHECK(std::string(e.what()).find("busted") != std::string::npos);
        CHECK(std::string(e.what()).find("FEL") != std::string::npos);
    }
}

TEST_CASE("Report: FEL and MEME reports embed the viewer, stylesheet and data") {
    std::string fel = render_html(R"J({"analysis": {"id": "fel"}, "dataset": {}})J", "/runs/cd2.FEL.json");
    CHECK(fel.find("<hyphy-fel-viewer") != std::string::npos);
    CHECK(fel.find("window.HyPhyFELViewer") != std::string::npos);
    CHECK(fel.find("window.HyPhyViewer") != std::string::npos);     // neutral alias of the same API
    CHECK(fel.find(".hyphy-fel-container") != std::string::npos);  // viewer stylesheet inlined
    CHECK(fel.find("FEL &bull; Pervasive Selection") != std::string::npos);
    CHECK(fel.find(">cd2.FEL.json<") != std::string::npos);        // file name only, not the path
    // No external resources: everything the page needs is inline
    CHECK(fel.find("<script src") == std::string::npos);
    CHECK(fel.find("<link") == std::string::npos);
    CHECK(fel.find("src=\"http") == std::string::npos);
    CHECK(fel.find("href=\"http") == std::string::npos);
    CHECK(fel.find("@import") == std::string::npos);

    std::string meme = render_html(R"J({"analysis": {"id": "meme"}})J", "m.json");
    CHECK(meme.find("<hyphy-meme-viewer") != std::string::npos);
    CHECK(meme.find("MEME &bull; Episodic Selection") != std::string::npos);
    CHECK(meme.find("whenDefined('hyphy-meme-viewer')") != std::string::npos);
    CHECK(meme.find("hyphy-fel-viewer, hyphy-meme-viewer { width: 100%;") != std::string::npos);
}

TEST_CASE("Report: the viewer table drives the supported analyses") {
    CHECK(supported_analyses() == std::vector<std::string>{"fel", "meme"});
    CHECK(supported_analyses_text() == "FEL, MEME");
    REQUIRE(find_viewer("meme") != nullptr);
    CHECK(find_viewer("meme")->tag == "hyphy-meme-viewer");
    CHECK(find_viewer("busted") == nullptr);
    CHECK(find_viewer("") == nullptr);
}

TEST_CASE("Report: data cannot break out of the inline script") {
    // A taxon name containing a closing script tag must not terminate the data block
    std::string html = render_html(R"J({"analysis": {"id": "fel"}, "taxa": ["</script><script>alert(1)</script>"]})J",
                                   "<evil>&.json");
    CHECK(count(html, "</script>") == 2);  // exactly the two real script blocks
    CHECK(html.find("\\u003c/script>\\u003cscript>alert(1)\\u003c/script>") != std::string::npos);
    CHECK(html.find("<evil>") == std::string::npos);  // file name is HTML-escaped
    CHECK(html.find("&lt;evil&gt;&amp;.json") != std::string::npos);
}
