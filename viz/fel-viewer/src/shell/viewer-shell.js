// ViewerShell: the layout every method shares, filled in by a method config
// KPI bar, hero, results table and collapsed diagnostics stack in the layout; the
// site drawer slides over it. The shell keeps the threshold and the focused site,
// fans them out to the mounted views, and handles loading and drag & drop.

import { methodConfig, configForResult } from '../methods/registry.js';

export class ViewerShell {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      method: null, // config to mount before any data arrives; otherwise the first load picks one
      theme: 'light',
      threshold: 0.1,
      filterType: null, // the method's default filter when unset
      onSiteSelect: null,
      onThresholdChange: null
    }, options);

    this.data = null;
    this.config = null;
    this.views = [];
    this.state = { threshold: this.options.threshold, focus: null };

    this.initDOM();
    const config = methodConfig(this.options.method);
    if (config) this.mount(config);
  }

  get threshold() {
    return this.state.threshold;
  }

  get selectedSite() {
    return this.state.focus;
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

    // Enable direct drag & drop onto the container
    this.initDragAndDrop();
  }

  // Build a method's views into the layout, replacing those of the previous method
  mount(config) {
    this.unmount();
    this.config = config;

    const hooks = {
      threshold: this.state.threshold,
      filterType: this.options.filterType || config.defaultFilter,
      onThresholdChange: (p) => this.setThreshold(p),
      onSiteSelect: (site, tab) => this.selectSite(site, tab),
      onClose: () => this.setState({ focus: null })
    };

    this.kpiBar = config.kpis(this.layout, hooks);
    this.hero = config.hero(this.layout, hooks);
    this.table = config.table(this.layout, hooks);
    this.initDiagnosticsSection();
    this.drawer = config.drawer(this.container, hooks);
    this.views = [this.kpiBar, this.hero, this.table, this.drawer];
  }

  unmount() {
    this.views.forEach(view => view.destroy?.());
    this.views = [];
    this.layout.replaceChildren();
    this.state.focus = null;
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

    grid.innerHTML = this.config.diagnostics(this.data);
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
      console.error(`Failed to load HyPhy JSON from ${url}:`, err);
    }
  }

  // Core API: load a parsed HyPhy result JSON. A method without a viewer of its
  // own falls back to the mounted config, or to FEL when nothing is mounted yet.
  load(json) {
    const config = configForResult(json) || this.config || methodConfig('fel');
    if (config !== this.config) this.mount(config);
    this.data = new config.Parser(json);

    // Propagate data to all subcomponents
    this.views.forEach(view => view.setData(this.data));
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

  // Shared view state fan-out: { threshold, focus } reach every mounted view
  setState(patch) {
    Object.assign(this.state, patch);
    if ('threshold' in patch) {
      this.views.forEach(view => view.setThreshold?.(this.state.threshold));
      this.updateDiagnostics();
    }
    if ('focus' in patch) {
      this.views.forEach(view => view.selectSite?.(this.state.focus));
    }
  }

  // Reactive Threshold Controller API
  setThreshold(p) {
    this.setState({ threshold: p });

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
    this.setState({ focus: siteNumber });

    if (siteNumber != null) {
      this.drawer?.inspectSite(siteNumber, activeTab);

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
      this.drawer?.close();
    }
  }

  // Table Filter API
  setFilter(filterType) {
    this.options.filterType = filterType;
    this.table?.setFilterType(filterType);
  }

  // Theme API
  setTheme(theme) {
    this.options.theme = theme;
    if (theme === 'dark') {
      this.container.classList.add('hyphy-theme-dark');
    } else {
      this.container.classList.remove('hyphy-theme-dark');
    }
    this.hero?.render();
  }
}
