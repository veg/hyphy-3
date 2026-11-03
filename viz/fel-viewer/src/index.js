// Main Entry Point: viewer classes & web components (<hyphy-viewer>, <hyphy-fel-viewer>, ...)
// ViewerShell owns the layout; the method registry picks what fills it for each result

import { ViewerShell, thresholdAttribute } from './shell/viewer-shell.js';

// Neutral viewer: mounts whichever method the loaded result (or `method` option) names
export class HyPhyViewer extends ViewerShell {}

export class HyPhyFELViewer extends ViewerShell {
  constructor(container, options = {}) {
    super(container, Object.assign({ method: 'fel' }, options));
  }
}

export class HyPhyMEMEViewer extends HyPhyFELViewer {
  constructor(container, options = {}) {
    super(container, Object.assign({ method: 'meme' }, options));
  }
}

// Custom Element / Web Component Definition: <hyphy-viewer method="...">
export class HyPhyViewerElement extends HTMLElement {
  constructor() {
    super();
    this.viewer = null;
  }

  createViewer(options) {
    return new HyPhyViewer(this, Object.assign({ method: this.getAttribute('method') }, options));
  }

  connectedCallback() {
    if (!this.viewer) {
      const theme = this.getAttribute('theme') || 'light';
      const threshold = thresholdAttribute(this.getAttribute('threshold'));
      const filter = this.getAttribute('filter') || 'diversifying';

      this.viewer = this.createViewer({
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
      this.viewer.options.threshold = thresholdAttribute(newVal);
      const threshold = this.viewer.seedThreshold();
      if (threshold != null) this.viewer.setThreshold(threshold);
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

// <hyphy-fel-viewer> & <hyphy-fel-widget>
export class HyPhyFELViewerElement extends HyPhyViewerElement {
  createViewer(options) {
    return new HyPhyFELViewer(this, options);
  }
}

export class HyPhyFELWidgetElement extends HyPhyFELViewerElement {}

// <hyphy-meme-viewer> & <hyphy-meme-widget>
export class HyPhyMEMEViewerElement extends HyPhyFELViewerElement {
  createViewer(options) {
    return new HyPhyMEMEViewer(this, options);
  }
}

export class HyPhyMEMEWidgetElement extends HyPhyMEMEViewerElement {}

if (typeof customElements !== 'undefined') {
  if (!customElements.get('hyphy-viewer')) {
    customElements.define('hyphy-viewer', HyPhyViewerElement);
  }

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
