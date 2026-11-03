# `@hyphy/ohmic-viewer`

Interactive HTML5 / Web Component visualizer for **HyPhy 3 OHMIC** (*Ohmic Phylogenetic Networks & Linear-Time Wheatstone Pruning*).

Built with the **Dieter Rams / Braun Industrial Minimalist Design System** matching `@hyphy/fel-viewer`, zero runtime dependencies, and full support for standalone offline HTML reports and native macOS `WKWebView` embedding.

---

## Key Features

1. **5-Panel Braun Industrial KPI Bar (`OhmicKpiBar`)**:
   - **Alignment & Electrical Solver**: Taxa, codons/sites, electrical super-nodes, solver mode (`ExactDense L⁺` vs. `Woodbury k-NN`), and `HYPHAEON 384D` transformer embedding badge.
   - **Foster's Theorem Invariant**: Displays $\sum_{e} c_e \Omega_e$ vs. theoretical $N - 1$ with exact trace residual check.
   - **Wheatstone 4-Way Bridges ($K$)**: Total series-parallel bridges, ML-retained shunts ($\hat{w}_{k,2} > 0.01$), and collapsed zero-conductance polytomies.
   - **Network $\Delta \ln L$ & $\Delta \text{AICc}$**: Log-likelihood gain (in nats) and $\Delta \text{AICc}$ improvement over the single binary tree.
   - **Interactive Shunt Stress Slider ($\tau_{\min}^{\text{Ohm}}$)**: Real-time filtering of Wheatstone shunts across the network canvas and bridge registry table.

2. **Phylogenetic Electrical Network & Wheatstone Diamond Canvas (`NetworkTreeCanvas`)**:
   - Switch between **Backbone + Shunts** and **Collapsed Polytomy Hubs**, and between **Phylogram** and **Cladogram** branch scaling.
   - Draws interactive **Wheatstone 4-Way Diamonds (`◆`)** at every bridge hub and dashed Braun-orange bypass shunt arcs connecting Clade $A \leftrightarrow$ Clade $C$.
   - Clicking any bridge highlights its 4 quartet clades ($A$: Steel Blue, $B$: Sister Teal, $C$: Shunt Partner Braun Orange, $D$: Outgroup) on the tree leaves.
   - High-DPI Retina **PNG** and vector **SVG** export.

3. **Pairwise Electrical Network Matrix (`FosterMatrixCanvas`)**:
   - Interactive $N \times N$ heatmap with 4-way selector:
     - **Foster Edge Probability ($P_{ij} = c_{ij}\Omega_{ij}$)**
     - **Effective Resistance ($\Omega_{ij}$)**
     - **Electrical Conductance ($c_{ij}$)**
     - **Metric Distance ($D_{ij}$)**
   - Synchronized Quartet Block highlighting ($(A \times B)$ vs. $(A \times C)$) when a Wheatstone bridge is selected.

4. **Site-Level Network Support Profile (`SiteReticulationPlot`)**:
   - Per-site log-likelihood difference $\Delta \ln L_s = \ln L_s^{\text{Ohm}} - \ln L_s^{\text{Tree}}$ across the alignment, pinpointing specific codons/sites that drive reticulate shunt support ($\Delta \ln L_s > 0$, Braun Orange) vs. the single tree ($\Delta \ln L_s < 0$, Steel Blue).

5. **Wheatstone Bridge Registry Table & Slide-Over Circuit/Simplex Inspector (`BridgesTable` & `BridgeDrawer`)**:
   - Sortable table of all $K$ bridges with CSV export.
   - Slide-over **Wheatstone Circuit & Quartet Simplex Inspector** featuring:
     - **4-Terminal Wheatstone Bridge Electrical Schematic SVG** showing the series-parallel pruning operator $D_{\text{hub}}^{\text{Ohm}} = \sum_{r=1}^3 w_{k,r} [P(t_{k,r})(M_{a_r} \odot M_{b_r})] \odot M_{c_r}$.
     - **Equilateral Barycentric Quartet Simplex SVG** plotting the trajectory from Initial Ohmic Weights $\mathbf{w}_k^{(0)}$ to ML-Fitted Weights $\hat{\mathbf{w}}_k$.

---

## Usage

### 1. Web Component (`<hyphy-ohmic-viewer>`)

```html
<link rel="stylesheet" href="dist/ohmic-viewer.css" />
<script type="module" src="dist/hyphy-ohmic-viewer.js"></script>

<hyphy-ohmic-viewer
  id="ohmic-viewer"
  theme="light"
  src="adh_ohmic.json">
</hyphy-ohmic-viewer>
```

### 2. Programmatic JavaScript API

```js
import { HyPhyOhmicViewer } from '@hyphy/ohmic-viewer';
import '@hyphy/ohmic-viewer/dist/ohmic-viewer.css';

const viewer = new HyPhyOhmicViewer('#mount', { theme: 'light' });
viewer.load(ohmicJsonPayload);
```

### 3. Generate Self-Contained Single-File HTML Report

```bash
npm run build:all
# Or build a standalone HTML report for any OHMIC JSON output:
python3 build_standalone_html.py path/to/output.ohmic.json path/to/report.html
```
