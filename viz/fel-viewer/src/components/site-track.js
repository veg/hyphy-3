// Site track: a HiDPI canvas with one codon x-scale shared by a main band and any
// bands stacked under it. The track owns the plumbing (resize, codon axis, hit testing,
// tooltip, selection ring, hover hairline, mode control, PNG/SVG export); bands own their
// marks.
//
// A band implements:
//   drawBackground(ctx, frame, track)   y grid under the codon axis (main band only)
//   drawMarks(ctx, frame, track)        rules and marks; returns the hit points
//                                       [{ site, x, y, raw, series }] (main band only)
//   draw(ctx, frame, track)             the whole band (stacked bands); may return hit
//                                       regions [{ x0, x1, y0, y1, key, kind }]
//   svg(frame, track)                   SVG elements for the export
//   tooltip(row, track, point)          optional tooltip HTML (main band only)
//   regionTooltip(region, track)        optional tooltip HTML for a hit region
//   select(region, track)               optional: a region of this band was clicked
//   measure(track)                      optional: sizes the band before each layout
//   height                              pixel height (stacked bands)
// A frame is { left, top, width, height, scaleX(site) }.
//
// A region of kind 'site' (key = site) hovers and selects that site like a main-band
// point; any other kind is the track's hoveredRegion/selectedRegion and goes to
// onRegionHover/onRegionSelect.

import { PALETTE } from '../data/palette.js';
import { downloadPNG, downloadSVG } from '../core/export.js';

export const TRACK_PADDING = { left: 52, right: 24, top: 24, bottom: 38, axis: 24, gap: 6 };

// Frames for the main band and each stacked band. The main band fills what the stack
// leaves; stacked bands sit under the codon axis labels, one gap apart.
export function stackBandFrames(width, height, stack, left = TRACK_PADDING.left) {
  const { right, top, bottom, axis, gap } = TRACK_PADDING;
  const plotW = width - left - right;
  const stackH = stack.reduce((sum, band) => sum + band.height + gap, 0);
  const main = { left, top, width: plotW, height: height - top - bottom - stackH };
  let y = top + main.height + axis;
  const frames = stack.map(band => {
    const frame = { left, top: y, width: plotW, height: band.height };
    y += band.height + gap;
    return frame;
  });
  return [main, ...frames];
}

// What the pointer at (mx, my) is over. Inside the frame of a stacked band that returned
// hit regions, only that band's regions count; elsewhere it is the main band's nearest
// point by x within the tolerance, and points sharing that x (a mirrored band's two
// series) are told apart by y. Returns { point } or { region }, or null.
export function resolveTrackHit({ points, regions, hitFrames }, mx, my, tolerance) {
  const frame = hitFrames.find(f => my >= f.top && my < f.top + f.height);
  if (frame) {
    const region = regions.find(r => mx >= r.x0 && mx < r.x1 && my >= r.y0 && my < r.y1);
    return region ? { region } : null;
  }

  let closest = null;
  let minDx = Infinity;
  let minDy = Infinity;
  for (const pt of points) {
    const dx = Math.abs(pt.x - mx);
    const dy = Math.abs(pt.y - my);
    if (dx < minDx || (dx === minDx && dy < minDy)) {
      minDx = dx;
      minDy = dy;
      closest = pt;
    }
  }
  return closest && minDx <= tolerance ? { point: closest } : null;
}

// Codon tick spacing for the shared x-axis
export function codonTickStep(numSites) {
  return numSites > 300 ? 50 : (numSites > 100 ? 20 : 10);
}

// Dashed horizontal threshold rule across a frame, with an optional label at its left end
export function drawThresholdRule(ctx, frame, y, color, label) {
  ctx.save();
  ctx.setLineDash([4, 4]);
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(frame.left, y);
  ctx.lineTo(frame.left + frame.width, y);
  ctx.stroke();
  if (label) {
    ctx.fillStyle = color;
    ctx.font = 'bold 9px "SF Mono", Monaco, monospace';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.fillText(label, frame.left + 6, y - 6);
  }
  ctx.restore();
}

// Amber reticle around the selected site's mark
export function drawSelectionRing(ctx, x, y) {
  ctx.strokeStyle = PALETTE.amberVivid;
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.arc(x, y, 11, 0, Math.PI * 2);
  ctx.stroke();
}

function trackModeMarkup(modes) {
  if (modes.length < 2) return '';
  const buttons = modes.map((mode, i) =>
    `            <button type="button" class="seg-btn${i === 0 ? ' active' : ''}" data-mode="${mode.id}">${mode.label}</button>`);
  return `<div class="segmented-control plot-mode-control">\n${buttons.join('\n')}\n          </div>\n          `;
}

export class SiteTrackCanvas {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      title: '',
      subtitle: '',
      modes: [],         // [{ id, label, subtitle(data), band }]; the active mode picks the main band
      band: null,        // the main band when there are no modes
      stack: [],         // bands drawn under the main one, on the same x-scale
      threshold: 0.1,
      sites: (data) => data.siteResults,
      siteCount: (data) => data.numSites,
      gutter: TRACK_PADDING.left, // left margin for y labels and stacked-band row labels
      exportName: 'hyphy_site_track',
      onSiteSelect: null,
      onSiteHover: null,
      onRegionSelect: null, // (region) for a clicked region that is not a site
      onRegionHover: null   // (region | null)
    }, options);

    this.threshold = this.options.threshold;
    this.modes = this.options.modes;
    this.mode = this.modes[0]?.id ?? null;
    this.stack = this.options.stack;
    this.data = null;
    this.selectedSite = null;
    this.selectedSeries = null; // which series of a mirrored band the selection came from
    this.hoveredSite = null;
    this.selectedRegion = null;
    this.hoveredRegion = null;
    this.points = [];    // Cached screen coordinates of the main band's marks
    this.regions = [];   // Hit regions the stacked bands returned, each with its band
    this.hitFrames = []; // Frames of the stacked bands that returned regions

    this.initDOM();
    this.initEvents();
    this.initResizeObserver();
  }

  get band() {
    return this.modes.find(m => m.id === this.mode)?.band ?? this.options.band;
  }

  get rows() {
    return this.options.sites(this.data);
  }

  get numSites() {
    return this.options.siteCount(this.data);
  }

  initDOM() {
    this.wrapper = document.createElement('div');
    this.wrapper.className = 'hyphy-evidence-plot-card';
    this.wrapper.innerHTML = `
      <div class="evidence-plot-header">
        <div class="plot-title-group">
          <span class="plot-title">${this.options.title}</span>
          <span class="plot-subtitle" id="plot-sub-label">${this.options.subtitle}</span>
        </div>
        <div class="plot-actions">
          ${trackModeMarkup(this.modes)}<div class="plot-export-group">
            <button type="button" class="plot-tool-btn" id="btn-export-png" title="Export High-Res PNG">PNG</button>
            <button type="button" class="plot-tool-btn" id="btn-export-svg" title="Export Vector SVG">SVG</button>
          </div>
        </div>
      </div>
      <div class="evidence-canvas-container" id="canvas-container">
        <canvas class="evidence-canvas" id="plot-canvas"></canvas>
        <div class="plot-tooltip" id="plot-tooltip"></div>
      </div>
    `;

    this.container.appendChild(this.wrapper);

    this.canvasContainer = this.wrapper.querySelector('#canvas-container');
    this.canvas = this.wrapper.querySelector('#plot-canvas');
    this.ctx = this.canvas.getContext('2d');
    this.tooltip = this.wrapper.querySelector('#plot-tooltip');
    this.subLabel = this.wrapper.querySelector('#plot-sub-label');

    // Bind mode buttons
    const modeBtns = this.wrapper.querySelectorAll('.plot-mode-control .seg-btn');
    modeBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        modeBtns.forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        this.setMode(btn.dataset.mode);
      });
    });

    // Bind export buttons
    this.wrapper.querySelector('#btn-export-png').addEventListener('click', () => this.exportPNG());
    this.wrapper.querySelector('#btn-export-svg').addEventListener('click', () => this.exportSVG());
  }

  initResizeObserver() {
    this.resizeObserver = new ResizeObserver(() => {
      this.resize();
    });
    this.resizeObserver.observe(this.canvasContainer);
  }

  destroy() {
    this.resizeObserver.disconnect();
  }

  resize() {
    if (!this.canvasContainer) return;
    const rect = this.canvasContainer.getBoundingClientRect();
    const width = Math.max(300, Math.floor(rect.width));
    const height = Math.max(220, Math.floor(rect.height || 280));

    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = width * dpr;
    this.canvas.height = height * dpr;
    this.canvas.style.width = `${width}px`;
    this.canvas.style.height = `${height}px`;

    this.width = width;
    this.height = height;
    this.dpr = dpr;

    this.render();
  }

  setData(data) {
    this.data = data;
    this.resize();
  }

  setThreshold(p) {
    this.threshold = p;
    this.render();
  }

  setMode(id) {
    this.mode = id;
    const mode = this.modes.find(m => m.id === id);
    if (mode?.subtitle) this.subLabel.textContent = mode.subtitle(this.data);
    this.render();
  }

  selectSite(siteNumber) {
    this.selectedSite = siteNumber;
    this.selectedSeries = null;
    this.render();
  }

  initEvents() {
    this.canvas.addEventListener('mousemove', (e) => this.handleMouseMove(e));
    this.canvas.addEventListener('mouseleave', () => this.handleMouseLeave());
    this.canvas.addEventListener('click', (e) => this.handleClick(e));
  }

  // The point or region under the pointer; points match by x within the given pixel
  // tolerance (wider when sites are sparse)
  hitAt(e, minTolerance) {
    if (!this.data) return null;
    const rect = this.canvas.getBoundingClientRect();
    const tolerance = Math.max(minTolerance, this.plotWidth / Math.max(1, this.numSites) * 2);
    return resolveTrackHit(this, e.clientX - rect.left, e.clientY - rect.top, tolerance);
  }

  // The site row behind a site region
  regionRow(region) {
    return region.raw ?? this.rows.find(r => r.site === region.key) ?? null;
  }

  setHoveredRegion(region) {
    const changed = region?.kind !== this.hoveredRegion?.kind || region?.key !== this.hoveredRegion?.key;
    this.hoveredRegion = region;
    if (changed && typeof this.options.onRegionHover === 'function') {
      this.options.onRegionHover(region);
    }
  }

  setHoveredSite(site) {
    const changed = site !== this.hoveredSite;
    this.hoveredSite = site;
    this.render();
    if (changed && typeof this.options.onSiteHover === 'function') {
      this.options.onSiteHover(site);
    }
  }

  handleMouseMove(e) {
    if (!this.data) return;
    const hit = this.hitAt(e, 12);
    const { point, region } = hit || {};
    this.setHoveredRegion(region ?? null);
    if (point) {
      this.setHoveredSite(point.site);
      this.showTooltip(this.band.tooltip?.(point.raw, this, point), e.clientX, e.clientY);
    } else if (region) {
      this.setHoveredSite(region.kind === 'site' ? region.key : null);
      this.showTooltip(region.band.regionTooltip?.(region, this), e.clientX, e.clientY);
    } else {
      this.setHoveredSite(null);
      this.hideTooltip();
    }
  }

  handleMouseLeave() {
    this.setHoveredRegion(null);
    this.setHoveredSite(null);
    this.hideTooltip();
  }

  handleClick(e) {
    const hit = this.hitAt(e, 16);
    if (!hit) return;
    const { point, region } = hit;
    if (point || region.kind === 'site') {
      const site = point ? point.site : region.key;
      this.selectedSite = site;
      this.selectedSeries = point?.series ?? null;
      this.render();
      if (typeof this.options.onSiteSelect === 'function') {
        this.options.onSiteSelect(site, point ? point.raw : this.regionRow(region));
      }
      return;
    }
    this.selectedRegion = region;
    region.band.select?.(region, this);
    this.render();
    if (typeof this.options.onRegionSelect === 'function') {
      this.options.onRegionSelect(region);
    }
  }

  showTooltip(html, clientX, clientY) {
    if (!html) {
      this.hideTooltip();
      return;
    }
    this.tooltip.innerHTML = html;
    this.tooltip.style.display = 'block';

    // Position tooltip safely within container
    const containerRect = this.canvasContainer.getBoundingClientRect();
    let left = clientX - containerRect.left + 16;
    let top = clientY - containerRect.top + 16;

    const tipRect = this.tooltip.getBoundingClientRect();
    if (left + tipRect.width > containerRect.width - 12) {
      left = clientX - containerRect.left - tipRect.width - 16;
    }
    if (top + tipRect.height > containerRect.height - 12) {
      top = clientY - containerRect.top - tipRect.height - 16;
    }

    this.tooltip.style.left = `${Math.max(8, left)}px`;
    this.tooltip.style.top = `${Math.max(8, top)}px`;
  }

  hideTooltip() {
    this.tooltip.style.display = 'none';
  }

  // Frames for the current size, each with the shared codon x-scale
  layout() {
    const numSites = Math.max(1, this.numSites);
    this.stack.forEach(b => b.measure?.(this));
    return stackBandFrames(this.width, this.height, this.stack, this.options.gutter).map(frame => Object.assign(frame, {
      scaleX: (site) => frame.left + ((site - 1) / Math.max(1, numSites - 1)) * frame.width
    }));
  }

  // Vertical codon gridlines through the main frame, labelled underneath. The gridlines
  // keep the line width the band's background left (its zero rule widens them).
  drawCodonAxis(ctx, frame) {
    const numSites = Math.max(1, this.numSites);
    const bottom = frame.top + frame.height;
    ctx.font = '10px "SF Mono", Monaco, monospace';
    ctx.fillStyle = PALETTE.stone;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';

    for (let s = 1; s <= numSites; s += codonTickStep(numSites)) {
      const x = frame.scaleX(s);
      ctx.strokeStyle = PALETTE.gridLine;
      ctx.beginPath();
      ctx.moveTo(x, frame.top);
      ctx.lineTo(x, bottom);
      ctx.stroke();

      ctx.fillText(String(s), x, bottom + 8);
    }
    // Last codon tick if not overlapping
    ctx.fillText(String(numSites), frame.scaleX(numSites), bottom + 8);
  }

  render() {
    if (!this.data || !this.width || !this.height) return;
    const ctx = this.ctx;

    ctx.save();
    ctx.scale(this.dpr, this.dpr);
    ctx.clearRect(0, 0, this.width, this.height);

    this.frames = this.layout();
    const [main, ...stacked] = this.frames;
    this.plotWidth = main.width;
    this.plotHeight = main.height;

    const band = this.band;
    band.drawBackground(ctx, main, this);
    this.drawCodonAxis(ctx, main);
    this.points = band.drawMarks(ctx, main, this);
    this.regions = [];
    this.hitFrames = [];
    this.stack.forEach((b, i) => {
      const regions = b.draw(ctx, stacked[i], this);
      if (!Array.isArray(regions)) return;
      this.hitFrames.push(stacked[i]);
      for (const region of regions) this.regions.push(Object.assign(region, { band: b }));
    });

    // Hover hairline through every band
    const hPt = this.hoveredSite ? this.points.find(p => p.site === this.hoveredSite) : null;
    if (hPt) {
      const last = this.frames[this.frames.length - 1];
      ctx.save();
      ctx.setLineDash([2, 2]);
      ctx.strokeStyle = PALETTE.ink;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(hPt.x, main.top);
      ctx.lineTo(hPt.x, last.top + last.height);
      ctx.stroke();
      ctx.restore();
    }

    ctx.restore();
  }

  exportPNG() {
    downloadPNG(this.canvas, `${this.options.exportName}.png`);
  }

  // Every band goes into the export, on the frames of the last render
  exportSVG() {
    if (!this.data || !this.frames) return;
    const w = this.width;
    const h = this.height;
    const [main, ...stacked] = this.frames;

    let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">\n`;
    svg += `  <style>\n    text { font-family: -apple-system, sans-serif; font-size: 10px; fill: ${PALETTE.stone}; }\n  </style>\n`;
    svg += `  <rect width="100%" height="100%" fill="${PALETTE.panel}"/>\n`;
    svg += this.band.svg(main, this);
    this.stack.forEach((b, i) => { svg += b.svg(stacked[i], this); });
    svg += `</svg>`;
    downloadSVG(svg, `${this.options.exportName}.svg`);
  }
}
