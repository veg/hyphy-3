// ViewerShell: the layout every method shares, filled in by a method config
// KPI bar, hero, results table and collapsed diagnostics stack in the layout; the
// drawer slides over it. The shell keeps the threshold, the statistic it applies to
// ('p' | 'q') and the focused row, fans them out to the mounted views, and handles
// loading and drag & drop.

import { METHOD_CONFIGS, methodConfig, resolveConfig } from '../methods/registry.js';
import { unsupportedNotice } from '../components/unsupported-notice.js';

// The cutoff a report starts at: the caller's own, else the cutoff the run used,
// else the method controller's seed
export function effectiveThreshold(explicit, data, config) {
  return explicit ?? data?.runThreshold ?? config?.controller?.seed ?? null;
}

// A threshold attribute as an explicit cutoff, or null when it is absent or not a
// positive number, so the report falls back to the run's own cutoff
export function thresholdAttribute(value) {
  const t = parseFloat(value);
  return t > 0 ? t : null;
}

export class ViewerShell {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      method: null, // config to mount before any data arrives; otherwise the first load picks one
      theme: 'light',
      threshold: null, // the run's own cutoff, else the method's seed
      statistic: null, // the method controller's statistic when unset
      filterType: null, // the method's default filter when unset
      onSiteSelect: null,
      onThresholdChange: null
    }, options);

    this.data = null;
    this.config = null;
    this.views = [];
    this.state = { threshold: null, statistic: null, focus: null };

    this.initDOM();
    const config = methodConfig(this.options.method);
    if (config) this.mount(config);
  }

  // The threshold and statistic a newly mounted config's views start from
  seedState(config, data) {
    this.state.threshold = effectiveThreshold(this.options.threshold, data, config);
    this.state.statistic = this.options.statistic ?? config.controller?.statistic ?? null;
  }

  get threshold() {
    return this.state.threshold;
  }

  get statistic() {
    return this.state.statistic;
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
  mount(config, data = null) {
    this.unmount();
    this.config = config;
    this.seedState(config, data);

    const hooks = {
      threshold: this.state.threshold,
      statistic: this.state.statistic,
      filterType: this.options.filterType || config.defaultFilter,
      onThresholdChange: (p) => this.setThreshold(p),
      onStatisticChange: (statistic) => this.setStatistic(statistic),
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
    this.kpiBar = this.hero = this.table = this.drawer = null;
    this.layout.replaceChildren();
    this.state.focus = null;
  }

  // Replace the report with a notice for a recognized method that has no viewer
  showUnsupported(method) {
    this.unmount();
    this.config = null;
    this.data = null;
    this.layout.innerHTML = unsupportedNotice(method, Object.keys(METHOD_CONFIGS));
    this.container.dispatchEvent(new CustomEvent('data-unsupported', {
      bubbles: true,
      composed: true,
      detail: { method }
    }));
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

  // Core API: load a parsed HyPhy result JSON. A result whose method is not
  // recognized falls back to the mounted config, or to FEL when nothing is mounted
  // yet; a recognized method without a viewer gets a notice instead of a report.
  // Every load reseeds the threshold from the caller, the run or the method.
  load(json) {
    const { method, config } = resolveConfig(json, this.config || methodConfig('fel'));
    if (!config) {
      this.showUnsupported(method);
      return;
    }
    const data = new config.Parser(json);
    if (config !== this.config) {
      this.mount(config, data);
    } else {
      const threshold = effectiveThreshold(this.options.threshold, data, config);
      if (threshold !== this.state.threshold) this.setState({ threshold });
    }
    this.data = data;

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

  // Shared view state fan-out: { threshold, statistic, focus } reach every mounted view
  setState(patch) {
    Object.assign(this.state, patch);
    if ('threshold' in patch) {
      this.views.forEach(view => view.setThreshold?.(this.state.threshold));
      this.updateDiagnostics();
    }
    if ('statistic' in patch) {
      this.views.forEach(view => view.setStatistic?.(this.state.statistic));
    }
    if ('focus' in patch) {
      this.views.forEach(view => view.selectSite?.(this.state.focus));
    }
  }

  // The cutoff a load would start from now: the caller's, the run's or the method's
  seedThreshold() {
    return effectiveThreshold(this.options.threshold, this.data, this.config);
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

  // Which statistic the cutoff applies to, for methods that report both p and q
  setStatistic(statistic) {
    this.setState({ statistic });

    this.container.dispatchEvent(new CustomEvent('statistic-change', {
      bubbles: true,
      composed: true,
      detail: { statistic }
    }));
  }

  // Row Selection & Drawer Inspector API: the key is a site, branch or segment,
  // looked up through the method config's row(data, key)
  selectSite(siteNumber, activeTab = null) {
    this.setState({ focus: siteNumber });

    if (siteNumber != null) {
      this.drawer?.inspectSite(siteNumber, activeTab);

      const r = this.data ? this.config.row(this.data, siteNumber) : null;
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
