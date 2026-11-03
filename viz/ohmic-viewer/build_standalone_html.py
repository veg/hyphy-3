#!/usr/bin/env python3
"""
Generates a 100% self-contained, single-file HTML report for HyPhy OHMIC results.
Inlines the JS engine, CSS stylesheet, and JSON dataset directly into the HTML.
Zero network access required; works directly via file://, double-click in Finder, or inside Cocoa WKWebView.
"""

import sys
import os
import json
import re


def build_standalone_html(json_path, output_html_path, title=None):
    dir_path = os.path.dirname(os.path.abspath(__file__))
    dist_js_path = os.path.join(dir_path, "dist", "hyphy-ohmic-viewer.umd.cjs")
    dist_css_path = os.path.join(dir_path, "dist", "ohmic-viewer.css")

    if not os.path.exists(dist_js_path) or not os.path.exists(dist_css_path):
        raise FileNotFoundError("Please build the library first via 'npm run build -- --mode lib'")

    with open(dist_js_path, "r", encoding="utf-8") as f:
        js_code = f.read()

    with open(dist_css_path, "r", encoding="utf-8") as f:
        css_code = f.read()

    # Strip sourcemaps so browser doesn't try to fetch local file:// sourcemap URLs
    js_code = re.sub(r'//# sourceMappingURL=.*', '', js_code)
    css_code = re.sub(r'/\*# sourceMappingURL=.* \*/', '', css_code)

    with open(json_path, "r", encoding="utf-8") as f:
        data_json_str = f.read()

    method_badge = "OHMIC • Electrical Network & Wheatstone Pruning"
    viewer_tag = "hyphy-ohmic-viewer"
    doc_title = title or f"HyPhy {method_badge} — {os.path.basename(json_path)}"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{doc_title}</title>
  <style>
    * {{
      box-sizing: border-box;
      margin: 0;
      padding: 0;
    }}
    html, body {{
      width: 100%;
      height: 100%;
      overflow: hidden;
      background-color: #F2F2EF;
      color: #1A1A18;
      font-family: -apple-system, BlinkMacSystemFont, "Helvetica Neue", Helvetica, Arial, sans-serif;
      -webkit-font-smoothing: antialiased;
    }}
    #app-container {{
      width: 100vw;
      height: 100vh;
      display: flex;
      flex-direction: column;
    }}
    /* Dieter Rams Braun Industrial Header */
    .app-nav {{
      height: 42px;
      background-color: #FFFFFF;
      border-bottom: 1px solid #D1D1CB;
      display: flex;
      align-items: center;
      justify-content: space-between;
      padding: 0 20px;
      font-size: 11px;
      z-index: 100;
    }}
    .app-brand {{
      display: flex;
      align-items: center;
      gap: 10px;
      font-weight: 700;
      color: #1A1A18;
      letter-spacing: 1.2px;
      text-transform: uppercase;
    }}
    .app-brand-badge {{
      background-color: #1A1A18;
      color: #FFFFFF;
      padding: 2px 6px;
      border-radius: 2px;
      font-size: 10px;
      font-weight: 700;
      letter-spacing: 0.8px;
    }}
    .app-brand-dot {{
      color: #EB5E28;
      font-size: 12px;
    }}
    .app-nav-title {{
      color: #767672;
      font-family: "SF Mono", Monaco, monospace;
      font-size: 11px;
      background-color: #EDEDE9;
      border: 1px solid #D1D1CB;
      border-radius: 2px;
      padding: 2px 8px;
    }}
    #viewer-mount {{
      flex: 1 1 0%;
      min-height: 0;
      display: flex;
      position: relative;
      overflow-y: auto;
      overflow-x: hidden;
    }}
    hyphy-ohmic-viewer, hyphy-ohmic-widget {{
      width: 100%;
      height: 100%;
      display: block;
    }}
    /* Inlined Component Styles */
    {css_code}
  </style>
</head>
<body>
  <div id="app-container">
    <nav class="app-nav">
      <div class="app-brand">
        <span class="app-brand-badge">HYPHY 3</span>
        <span>{method_badge}</span>
        <span class="app-brand-dot">●</span>
      </div>
      <div class="app-nav-title">
        <span>{os.path.basename(json_path)}</span>
      </div>
    </nav>
    <div id="viewer-mount">
      <{viewer_tag} id="ohmic-viewer" theme="light"></{viewer_tag}>
    </div>
  </div>

  <script>
    // Inlined Standalone Engine
    {js_code}

    const payload = {data_json_str};

    function startViewer() {{
      const viewer = document.getElementById('ohmic-viewer');
      if (window.customElements && customElements.whenDefined) {{
        customElements.whenDefined('{viewer_tag}').then(function() {{
          if (viewer && typeof viewer.load === 'function') {{
            viewer.load(payload);
          }}
        }});
      }} else if (viewer && typeof viewer.load === 'function') {{
        viewer.load(payload);
      }} else if (typeof HyPhyOhmicViewer !== 'undefined' && HyPhyOhmicViewer.HyPhyOhmicViewer) {{
        const mount = document.getElementById('viewer-mount');
        mount.innerHTML = '';
        const v = new HyPhyOhmicViewer.HyPhyOhmicViewer(mount, {{ theme: 'light' }});
        v.load(payload);
      }}
    }}

    if (document.readyState === 'loading') {{
      document.addEventListener('DOMContentLoaded', startViewer);
    }} else {{
      startViewer();
    }}
  </script>
</body>
</html>
"""

    with open(output_html_path, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"Generated standalone HTML report: {output_html_path} ({os.path.getsize(output_html_path):,} bytes)")


if __name__ == "__main__":
    dir_path = os.path.dirname(os.path.abspath(__file__))
    if len(sys.argv) < 3:
        for name in ["cd2_ohmic", "bglobin_ohmic", "adh_ohmic"]:
            json_p = os.path.join(dir_path, "public", "samples", f"{name}.json")
            if os.path.exists(json_p):
                build_standalone_html(json_p, os.path.join(dir_path, "dist", f"{name}_report.html"))
    else:
        build_standalone_html(sys.argv[1], sys.argv[2])
