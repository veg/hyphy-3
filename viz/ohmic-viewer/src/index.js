import './styles/ohmic-viewer.css';
import { OhmicDataModel } from './data/ohmic-data.js';
import { OhmicKpiBar } from './components/ohmic-kpi-bar.js';
import { NetworkTreeCanvas } from './components/network-tree-canvas.js';
import { FosterMatrixCanvas } from './components/foster-matrix-canvas.js';
import { SiteReticulationPlot } from './components/site-reticulation-plot.js';
import { BridgesTable } from './components/bridges-table.js';
import { BridgeDrawer } from './components/bridge-drawer.js';

export class HyPhyOhmicViewer {
  constructor(container, options = {}) {
    this.container = typeof container === 'string' ? document.querySelector(container) : container;
    this.options = options;
    this.model = null;
    this.tauThreshold = options.initialTau ?? 0.0;
    this.selectedBridgeIndex = null;

    this._initLayout();
  }

  _initLayout() {
    this.container.classList.add('ohmic-viewer-root');
    if (this.options.theme) {
      this.container.setAttribute('theme', this.options.theme);
    }

    this.container.innerHTML = `
      <div class="ohmic-shell">
        <div id="ohmic-kpi-mount"></div>
        <div class="ohmic-main-split">
          <div class="ohmic-panel" id="ohmic-tree-mount"></div>
          <div class="ohmic-panel" id="ohmic-matrix-mount"></div>
        </div>
        <div class="ohmic-panel" id="ohmic-site-mount"></div>
        <div class="ohmic-panel" id="ohmic-bridges-mount"></div>
        <div class="ohmic-drawer" id="ohmic-drawer-mount"></div>
      </div>
    `;

    this.kpiBar = new OhmicKpiBar(this.container.querySelector('#ohmic-kpi-mount'), {
      initialTau: this.tauThreshold,
      onTauThresholdChange: (tau) => {
        this.tauThreshold = tau;
        this.treeCanvas.setTauThreshold(tau);
        this.bridgesTable.setTauThreshold(tau);
      }
    });

    this.treeCanvas = new NetworkTreeCanvas(this.container.querySelector('#ohmic-tree-mount'), {
      onBridgeSelect: (bridge) => this.selectBridge(bridge)
    });

    this.matrixCanvas = new FosterMatrixCanvas(this.container.querySelector('#ohmic-matrix-mount'));
    this.sitePlot = new SiteReticulationPlot(this.container.querySelector('#ohmic-site-mount'));

    this.bridgesTable = new BridgesTable(this.container.querySelector('#ohmic-bridges-mount'), {
      onBridgeSelect: (bridge) => this.selectBridge(bridge)
    });

    this.bridgeDrawer = new BridgeDrawer(this.container.querySelector('#ohmic-drawer-mount'), {
      onClose: () => {
        this.selectedBridgeIndex = null;
        this.treeCanvas.setSelectedBridge(null);
        this.matrixCanvas.setSelectedBridge(null);
        this.bridgesTable.setSelectedBridge(null);
      }
    });
  }

  selectBridge(bridge) {
    if (!bridge) return;
    this.selectedBridgeIndex = bridge.index;
    this.treeCanvas.setSelectedBridge(bridge.index);
    this.matrixCanvas.setSelectedBridge(bridge.index);
    this.bridgesTable.setSelectedBridge(bridge.index);
    this.bridgeDrawer.open(bridge);
  }

  load(jsonOrObject) {
    const raw = typeof jsonOrObject === 'string' ? JSON.parse(jsonOrObject) : jsonOrObject;
    this.model = new OhmicDataModel(raw);
    this.selectedBridgeIndex = this.model.bridges.length > 0 ? this.model.bridges[0].index : null;

    this.kpiBar.render(this.model, this.tauThreshold);
    this.treeCanvas.setModel(this.model, this.tauThreshold, this.selectedBridgeIndex);
    this.matrixCanvas.setModel(this.model, this.selectedBridgeIndex);
    this.sitePlot.setModel(this.model);
    this.bridgesTable.setModel(this.model, this.tauThreshold, this.selectedBridgeIndex);
    this.bridgeDrawer.close();
  }

  setTheme(theme) {
    this.container.setAttribute('theme', theme);
    this.treeCanvas.draw();
    this.matrixCanvas.draw();
    this.sitePlot.draw();
  }
}

class HyPhyOhmicElement extends HTMLElement {
  static get observedAttributes() {
    return ['src', 'theme'];
  }

  connectedCallback() {
    if (!this._viewer) {
      this._viewer = new HyPhyOhmicViewer(this, {
        theme: this.getAttribute('theme') || 'light'
      });
    }
    const src = this.getAttribute('src');
    if (src) {
      this.loadFromUrl(src);
    }
  }

  attributeChangedCallback(name, oldVal, newVal) {
    if (oldVal === newVal || !this._viewer) return;
    if (name === 'theme') {
      this._viewer.setTheme(newVal);
    } else if (name === 'src' && newVal) {
      this.loadFromUrl(newVal);
    }
  }

  async loadFromUrl(url) {
    const resp = await fetch(url);
    const json = await resp.json();
    this.load(json);
  }

  load(data) {
    if (!this._viewer) {
      this._viewer = new HyPhyOhmicViewer(this, {
        theme: this.getAttribute('theme') || 'light'
      });
    }
    this._viewer.load(data);
  }
}

if (typeof window !== 'undefined' && window.customElements) {
  if (!window.customElements.get('hyphy-ohmic-viewer')) {
    window.customElements.define('hyphy-ohmic-viewer', HyPhyOhmicElement);
  }
  if (!window.customElements.get('hyphy-ohmic-widget')) {
    window.customElements.define('hyphy-ohmic-widget', class extends HyPhyOhmicElement {});
  }
}

export { OhmicDataModel };
