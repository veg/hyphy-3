// Overview track: a site track with an overview brush above its canvas for long
// alignments. Past brush.minSites the brush shows the whole alignment and the track
// zooms to the brushed window; shorter alignments draw the whole track and hide the
// brush. PREV/NEXT SIG and clicks on the brush select a site, scrolling the window to it.
//
// options.track(container) builds the SiteTrackCanvas; options.brush holds the
// OverviewBrush options (value, tone, rule, toneLabels, maxValue, minSites); window is
// the number of sites the track shows when zoomed, or window(numSites) for one that scales
// with the alignment; onSiteSelect(site) reports a site
// picked on either part.

import { OverviewBrush } from './overview-brush.js';

// The window of span sites (0-based start, end exclusive) that keeps site in view: the
// current one when it already holds the site, else one centred on it within the alignment
export function revealWindow(viewport, site, span, numSites) {
  if (site - 1 >= viewport.start && site - 1 < viewport.end) return viewport;
  const start = Math.max(0, Math.min(numSites - span, site - 1 - Math.floor(span / 2)));
  return { start, end: start + span };
}

export class OverviewTrack {
  constructor(container, options = {}) {
    this.options = Object.assign({
      track: null,
      brush: {},
      window: 150,
      onSiteSelect: null
    }, options);

    this.track = this.options.track(container);
    this.brushHost = document.createElement('div');
    this.brushHost.className = 'overview-track-brush';
    this.track.canvasContainer.before(this.brushHost);
    this.brush = new OverviewBrush(this.brushHost, Object.assign({}, this.options.brush, {
      threshold: this.track.threshold,
      onSiteSelect: (site) => this.pick(site),
      onViewportChange: (start, end) => this.track.setViewport(start + 1, end)
    }));
    this.span = 0;

    // The brush fills the host's content box, inside its padding
    this.resizeObserver = new ResizeObserver(([entry]) => this.brush.resize(Math.floor(entry.contentRect.width)));
    this.resizeObserver.observe(this.brushHost);
  }

  get zoomed() {
    return this.span > 0;
  }

  setData(data) {
    this.track.setData(data);
    const numSites = this.track.numSites;
    this.brush.setData(this.track.rows, numSites);
    const { window: size } = this.options;
    this.span = numSites >= this.brush.options.minSites ? Math.min(numSites, typeof size === 'function' ? size(numSites) : size) : 0;
    this.brushHost.style.display = this.zoomed ? '' : 'none';
    if (this.zoomed) {
      this.brush.setViewport(0, this.span);
      this.track.setViewport(1, this.span);
    } else {
      this.track.setViewport(null);
    }
  }

  setThreshold(t) {
    this.track.setThreshold(t);
    this.brush.setThreshold(t);
  }

  // Focus from the shell: mark the site on both parts and bring it into the window
  selectSite(site) {
    this.track.selectSite(site);
    this.brush.selectSite(site);
    if (site != null && this.zoomed) this.reveal(site);
  }

  reveal(site) {
    const next = revealWindow(this.brush.viewport, site, this.span, this.track.numSites);
    if (next === this.brush.viewport) return;
    this.brush.setViewport(next.start, next.end);
    this.track.setViewport(next.start + 1, next.end);
  }

  // A site picked on the brush
  pick(site) {
    this.reveal(site);
    this.track.selectSite(site);
    if (typeof this.options.onSiteSelect === 'function') this.options.onSiteSelect(site);
  }

  render() {
    this.track.render();
    this.brush.dirtyCache = true;
    this.brush.render();
  }

  destroy() {
    this.resizeObserver.disconnect();
    this.track.destroy();
    this.brush.destroy();
  }
}
