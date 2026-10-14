// Main Entry Point: HyPhyFELViewer Class & <hyphy-fel-viewer> Web Component
// First-principles redesign: Hero Evidence Plot, Diversifying-Only Default, Embeddable Everywhere

import { FELData } from './data/fel-data.js';
import { KpiBar } from './components/kpi-bar.js';
import { EvidencePlot } from './components/evidence-plot.js';
import { CandidatesTable } from './components/candidates-table.js';
import { SiteDrawer } from './components/site-drawer.js';
import { formatOmegaInf } from './core/format.js';

export class HyPhyFELViewer {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      theme: 'light',
      threshold: 0.1,
      filterType: 'diversifying', // Diversifying only by default!
      onSiteSelect: null,
      onThresholdChange: null
    }, options);

    this.data = null;
    this.threshold = this.options.threshold;
    this.selectedSite = null;

    this.initDOM();
  }

  initDOM() {
    this.container.classList.add('hyphy-fel-container');
    if (this.options.theme === 'dark') {
      this.container.classList.add('hyphy-theme-dark');
    }

    // Main layout wrapper
    this.layout = document.createElement('div');
    this.layout.className = 'hyphy-fel-layout';
    this.container.appendChild(this.layout);

    // 1. KPI Summary Bar & Reactive Threshold Controller
    this.kpiBar = new KpiBar(this.layout, {
      threshold: this.threshold,
      onThresholdChange: (p) => this.setThreshold(p)
    });

    // 2. Hero Visual: Selection Evidence Landscape Plot (-log10 p)
    this.plot = new EvidencePlot(this.layout, {
      threshold: this.threshold,
      onSiteSelect: (site, siteResult) => this.selectSite(site)
    });

    // 3. Focused Candidate Sites Table (Diversifying Only by Default)
    this.table = new CandidatesTable(this.layout, {
      threshold: this.threshold,
      filterType: this.options.filterType,
      onSiteSelect: (site, siteResult, action) => this.selectSite(site, action)
    });

    // 4. Collapsible Diagnostics Section (Model Fits & Tree)
    this.initDiagnosticsSection();

    // 5. Slide-over Site Inspector Drawer (Progressive Disclosure)
    this.drawer = new SiteDrawer(this.container, {
      threshold: this.threshold,
      onClose: () => {
        this.selectedSite = null;
        this.plot.selectSite(null);
        this.table.selectSite(null);
      }
    });

    // Enable direct drag & drop onto the container
    this.initDragAndDrop();
  }

  initDiagnosticsSection() {
    this.diagWrapper = document.createElement('div');
    this.diagWrapper.className = 'hyphy-diagnostics-card';
    this.diagWrapper.innerHTML = `
      <button type="button" class="diagnostics-toggle-btn" id="diag-toggle-btn">
        <span>Model Fits &amp; Phylogeny Diagnostics</span>
        <span id="diag-toggle-icon">▼</span>
      </button>
      <div class="diagnostics-content" id="diag-content">
        <div class="diag-grid" id="diag-grid">
          <!-- Dynamic fit metrics -->
        </div>
        <div class="diag-item" style="margin-top: 8px;">
          <span class="diag-lbl">Phylogenetic Tree (Newick)</span>
          <div class="newick-box" id="diag-newick">No tree loaded</div>
        </div>
      </div>
    `;

    this.layout.appendChild(this.diagWrapper);

    const toggleBtn = this.diagWrapper.querySelector('#diag-toggle-btn');
    const content = this.diagWrapper.querySelector('#diag-content');
    const icon = this.diagWrapper.querySelector('#diag-toggle-icon');

    toggleBtn.addEventListener('click', () => {
      const isOpen = content.classList.toggle('active');
      icon.textContent = isOpen ? '▲' : '▼';
    });
  }

  updateDiagnostics() {
    if (!this.data) return;
    const grid = this.diagWrapper.querySelector('#diag-grid');
    const newickBox = this.diagWrapper.querySelector('#diag-newick');

    const isMeme = this.data.isMEME;
    const m = this.data.analysisType ? this.data.analysisType.toUpperCase() : 'FEL';

    let html = `
      <div class="diag-item">
        <span class="diag-lbl">Analysis Method</span>
        <span class="diag-val">${isMeme ? 'MEME (Mixed Effects Model of Episodic Selection)' : `${m} (Fixed Effects Likelihood)`}</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">${isMeme ? 'Hypothesis Test' : 'Codon Sites'}</span>
        <span class="diag-val font-mono">${isMeme ? 'H₀: β⁺ ≤ α vs H₁: β⁺ > α' : this.data.numSites}</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">${isMeme ? 'Codon Sites' : 'Taxa / Sequences'}</span>
        <span class="diag-val font-mono">${isMeme ? this.data.numSites : this.data.numTaxa}</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">${isMeme ? 'Taxa / Sequences' : 'Global Mean dN/dS (ω)'}</span>
        <span class="diag-val font-mono">${isMeme ? this.data.numTaxa : formatOmegaInf(this.data.globalOmega, 3)}</span>
      </div>
    `;

    // Legacy or modern fit details
    const raw = this.data.raw || {};
    const codonFit = raw.model_fits?.codon_mg94 || raw.fits?.['Global MG94xREV'];
    if (codonFit) {
      const logL = codonFit.log_likelihood ?? codonFit['Log Likelihood'];
      const aicc = codonFit.aicc ?? codonFit['AIC-c'];
      if (logL != null) {
        html += `
          <div class="diag-item">
            <span class="diag-lbl">MG94 Log-Likelihood</span>
            <span class="diag-val font-mono">${Number(logL).toFixed(2)}</span>
          </div>
        `;
      }
      if (aicc != null) {
        html += `
          <div class="diag-item">
            <span class="diag-lbl">MG94 AICc</span>
            <span class="diag-val font-mono">${Number(aicc).toFixed(2)}</span>
          </div>
        `;
      }
    }

    grid.innerHTML = html;
    newickBox.textContent = this.data.treeNewick || 'No Newick string available in JSON';
  }

  initDragAndDrop() {
    this.container.addEventListener('dragover', (e) => {
      e.preventDefault();
      e.stopPropagation();
    });

    this.container.addEventListener('drop', (e) => {
      e.preventDefault();
      e.stopPropagation();
      const file = e.dataTransfer?.files?.[0];
      if (file) {
        this.loadFile(file);
      }
    });
  }

  // Load a File object directly
  loadFile(file) {
    const reader = new FileReader();
    reader.onload = (e) => {
      try {
        const json = JSON.parse(e.target.result);
        this.load(json);
      } catch (err) {
        alert('Invalid JSON file format: ' + err.message);
      }
    };
    reader.readAsText(file);
  }

  // Load from a remote or local URL
  async loadFromUrl(url) {
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(`HTTP ${res.status}: ${res.statusText}`);
      const json = await res.json();
      this.load(json);
    } catch (err) {
      console.error(`Failed to load FEL JSON from ${url}:`, err);
    }
  }

  // Core API: Load a parsed HyPhy FEL JSON object
  load(json) {
    this.data = new FELData(json);

    // Propagate data to all subcomponents
    this.kpiBar.setData(this.data);
    this.plot.setData(this.data);
    this.table.setData(this.data);
    this.drawer.setData(this.data);
    this.updateDiagnostics();

    // Notify listeners
    this.container.dispatchEvent(new CustomEvent('data-loaded', {
      bubbles: true,
      composed: true,
      detail: {
        numSites: this.data.numSites,
        numTaxa: this.data.numTaxa,
        stats: this.data.getStats(this.threshold)
      }
    }));
  }

  // Reactive Threshold Controller API
  setThreshold(p) {
    this.threshold = p;
    this.kpiBar.setThreshold(p);
    this.plot.setThreshold(p);
    this.table.setThreshold(p);
    this.drawer.setThreshold(p);
    this.updateDiagnostics();

    this.container.dispatchEvent(new CustomEvent('threshold-change', {
      bubbles: true,
      composed: true,
      detail: { threshold: p }
    }));

    if (typeof this.options.onThresholdChange === 'function') {
      this.options.onThresholdChange(p);
    }
  }

  // Site Selection & Drawer Inspector API
  selectSite(siteNumber, activeTab = null) {
    this.selectedSite = siteNumber;
    this.plot.selectSite(siteNumber);
    this.table.selectSite(siteNumber);

    if (siteNumber != null) {
      this.drawer.inspectSite(siteNumber, activeTab);

      const r = this.data ? this.data.getSite(siteNumber) : null;
      this.container.dispatchEvent(new CustomEvent('site-select', {
        bubbles: true,
        composed: true,
        detail: { site: siteNumber, raw: r, tab: activeTab }
      }));

      if (typeof this.options.onSiteSelect === 'function') {
        this.options.onSiteSelect(siteNumber, r, activeTab);
      }
    } else {
      this.drawer.close();
    }
  }

  // Table Filter API
  setFilter(filterType) {
    this.table.setFilterType(filterType);
  }

  // Theme API
  setTheme(theme) {
    this.options.theme = theme;
    if (theme === 'dark') {
      this.container.classList.add('hyphy-theme-dark');
    } else {
      this.container.classList.remove('hyphy-theme-dark');
    }
    this.plot.render();
  }
}

// Custom Element / Web Component Definition: <hyphy-fel-viewer> & <hyphy-fel-widget>
export class HyPhyFELViewerElement extends HTMLElement {
  constructor() {
    super();
    this.viewer = null;
  }

  connectedCallback() {
    if (!this.viewer) {
      const theme = this.getAttribute('theme') || 'light';
      const threshold = parseFloat(this.getAttribute('threshold')) || 0.1;
      const filter = this.getAttribute('filter') || 'diversifying';

      this.viewer = new HyPhyFELViewer(this, {
        theme,
        threshold,
        filterType: filter
      });

      const src = this.getAttribute('src');
      if (src) {
        this.viewer.loadFromUrl(src);
      }
    }
  }

  static get observedAttributes() {
    return ['src', 'threshold', 'theme', 'filter'];
  }

  attributeChangedCallback(name, oldVal, newVal) {
    if (!this.viewer || oldVal === newVal) return;
    if (name === 'src' && newVal) {
      this.viewer.loadFromUrl(newVal);
    } else if (name === 'threshold') {
      this.viewer.setThreshold(parseFloat(newVal) || 0.1);
    } else if (name === 'theme') {
      this.viewer.setTheme(newVal);
    } else if (name === 'filter') {
      this.viewer.setFilter(newVal);
    }
  }

  // Passthrough methods
  load(data) {
    if (this.viewer) this.viewer.load(data);
  }

  loadFromUrl(url) {
    if (this.viewer) this.viewer.loadFromUrl(url);
  }

  setThreshold(p) {
    if (this.viewer) this.viewer.setThreshold(p);
  }

  selectSite(s) {
    if (this.viewer) this.viewer.selectSite(s);
  }

  setFilter(f) {
    if (this.viewer) this.viewer.setFilter(f);
  }
}

export class HyPhyFELWidgetElement extends HyPhyFELViewerElement {}

export class HyPhyMEMEViewer extends HyPhyFELViewer {}
export class HyPhyMEMEViewerElement extends HyPhyFELViewerElement {}
export class HyPhyMEMEWidgetElement extends HyPhyFELViewerElement {}

if (typeof customElements !== 'undefined') {
  if (!customElements.get('hyphy-fel-viewer')) {
    customElements.define('hyphy-fel-viewer', HyPhyFELViewerElement);
  }

  if (!customElements.get('hyphy-fel-widget')) {
    customElements.define('hyphy-fel-widget', HyPhyFELWidgetElement);
  }

  if (!customElements.get('hyphy-meme-viewer')) {
    customElements.define('hyphy-meme-viewer', HyPhyMEMEViewerElement);
  }

  if (!customElements.get('hyphy-meme-widget')) {
    customElements.define('hyphy-meme-widget', HyPhyMEMEWidgetElement);
  }
}
