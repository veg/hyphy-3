#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "hyphy/report/html_report.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

using hyphy::report::detect_analysis;
using hyphy::report::find_viewer;
using hyphy::report::render_html;
using hyphy::report::supported_analyses;
using hyphy::report::supported_analyses_text;
using nlohmann::json;

static size_t count(const std::string& haystack, const std::string& needle) {
    size_t n = 0;
    for (size_t pos = haystack.find(needle); pos != std::string::npos; pos = haystack.find(needle, pos + 1)) ++n;
    return n;
}

// tests/data/detect_cases.json is shared with the viewer's detectMethod test, so
// both detectors are held to the same answers
TEST_CASE("Report: analysis detection agrees with the shared test vectors") {
    std::ifstream in(HYPHY3_SOURCE_DIR "/tests/data/detect_cases.json");
    REQUIRE(in.good());
    const json cases = json::parse(in)["cases"];
    REQUIRE(cases.size() > 0);
    for (const auto& c : cases) {
        const std::string expected = c["method"].is_null() ? "" : c["method"].get<std::string>();
        INFO(c["name"].get<std::string>());
        CHECK(detect_analysis(c["json"]) == expected);
    }
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
    CHECK(meme.find("hyphy-fel-viewer, hyphy-meme-viewer, hyphy-viewer, hyphy-ohmic-viewer { width: 100%;") != std::string::npos);

    // FUBAR mounts the neutral element, which picks its viewer from the result itself
    std::string fubar = render_html(R"J({"analysis": {"id": "fubar"}})J", "f.json");
    CHECK(fubar.find("<hyphy-viewer id=\"hyphy-viewer\"") != std::string::npos);
    CHECK(fubar.find("FUBAR &bull; Bayesian Pervasive Selection") != std::string::npos);
    CHECK(fubar.find("whenDefined('hyphy-viewer')") != std::string::npos);

    std::string slac = render_html(R"J({"analysis": {"id": "slac"}})J", "s.json");
    CHECK(slac.find("<hyphy-viewer id=\"hyphy-viewer\"") != std::string::npos);
    CHECK(slac.find("SLAC &bull; Ancestral Counting") != std::string::npos);
    // Viewers that share an element name it once in the sizing rule
    CHECK(count(slac, "hyphy-viewer, hyphy-viewer") == 0);

    std::string dfel = render_html(R"J({"analysis": {"id": "dfel"}})J", "d.json");
    CHECK(dfel.find("<hyphy-viewer id=\"hyphy-viewer\"") != std::string::npos);
    CHECK(dfel.find("dFEL &bull; Directional Selection") != std::string::npos);
    // Other spellings of the dFEL id get the same report
    std::string dfel_alias = render_html(R"J({"analysis": {"id": "directional_fel"}})J", "d.json");
    CHECK(dfel_alias.find("<hyphy-viewer id=\"hyphy-viewer\"") != std::string::npos);
    CHECK(dfel_alias.find("dFEL &bull; Directional Selection") != std::string::npos);

    // Ohmic mounts its own viewer from its own bundle
    std::string ohmic = render_html(R"J({"analysis": {"id": "ohmic"}})J", "o.json");
    CHECK(ohmic.find("<hyphy-ohmic-viewer") != std::string::npos);
    CHECK(ohmic.find("window.HyPhyOhmicViewer") != std::string::npos);
    CHECK(ohmic.find("window.HyPhyFELViewer") == std::string::npos);
    CHECK(ohmic.find("OHMIC &bull; Phylogenetic Network") != std::string::npos);

    std::string gard = render_html(R"J({"analysis": {"id": "gard"}})J", "g.json");
    CHECK(gard.find("<hyphy-viewer id=\"hyphy-viewer\"") != std::string::npos);
    CHECK(gard.find("GARD &bull; Recombination Breakpoints") != std::string::npos);
}

TEST_CASE("Report: the viewer table drives the supported analyses") {
    CHECK(supported_analyses() == std::vector<std::string>{"fel", "meme", "fubar", "slac", "dfel", "gard", "ohmic"});
    CHECK(supported_analyses_text() == "FEL, MEME, FUBAR, SLAC, DFEL, GARD, OHMIC");
    REQUIRE(find_viewer("meme") != nullptr);
    CHECK(find_viewer("meme")->tag == "hyphy-meme-viewer");
    REQUIRE(find_viewer("fubar") != nullptr);
    CHECK(find_viewer("fubar")->tag == "hyphy-viewer");
    REQUIRE(find_viewer("slac") != nullptr);
    CHECK(find_viewer("slac")->tag == "hyphy-viewer");
    REQUIRE(find_viewer("dfel") != nullptr);
    CHECK(find_viewer("dfel")->tag == "hyphy-viewer");
    REQUIRE(find_viewer("gard") != nullptr);
    CHECK(find_viewer("gard")->tag == "hyphy-viewer");
    REQUIRE(find_viewer("ohmic") != nullptr);
    CHECK(find_viewer("ohmic")->tag == "hyphy-ohmic-viewer");
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

