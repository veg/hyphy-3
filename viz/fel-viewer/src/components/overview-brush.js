// Overview brush: a whole-alignment strip above long site tracks. It draws a 1D profile
// with one bar per site, a draggable brush over the window the main track shows, a needle
// on the selected site, and PREV/NEXT buttons that jump between significant sites.
//
// Rows are { site, ... }; the options say how to read them:
//   value(row)              bar height on a 0..maxValue scale (default -log10 p_value)
//   tone(row, threshold)    'accent' or 'purifying' when the site passes, otherwise null
//   rule(threshold)         value of the dashed threshold rule, or null for none
// A site is significant when tone() returns a tone.

import { PALETTE, withAlpha } from '../data/palette.js';

// The next significant site after current (direction > 0) or before it, wrapping around.
// sites is sorted ascending; returns null when it is empty.
export function nextSignificantSite(sites, current, direction) {
  if (!sites.length) return null;
  if (direction > 0) return sites.find(site => site > current) ?? sites[0];
  return [...sites].reverse().find(site => site < current) ?? sites[sites.length - 1];
}

export class OverviewBrush {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      height: 56,
      maxValue: 5,                // -log10 p of 1e-5 fills the profile
      value: row => -Math.log10(Math.max(row.p_value, 1e-15)),
      tone: (row, threshold) => (row.p_value <= threshold ? 'accent' : null),
      rule: threshold => -Math.log10(Math.max(threshold, 1e-15)),
      toneLabels: { accent: 'Diversifying ▲', purifying: 'Purifying ▼' },
      unit: 'Sites',
      minSites: 300,              // shorter alignments fit the main track, so the brush hides
      threshold: 0.1,
      onSiteSelect: null,         // (site, row)
      onViewportChange: null      // (start, end), 0-based and end-exclusive
    }, options);

    this.rows = [];
    this.numSites = 0;
    this.significant = [];
    this.selectedSite = null;
    this.hoveredSite = null;
    this.viewport = { start: 0, end: 0 };

    // The profile only changes with the data, threshold or width, so it is cached offscreen
    this.offscreenCanvas = document.createElement('canvas');
    this.offscreenCtx = this.offscreenCanvas.getContext('2d');
    this.dirtyCache = true;

    this.drag = null;

    this.initDOM();
    this.setupEvents();
  }

  initDOM() {
    this.wrapper = document.createElement('div');
    this.wrapper.className = 'hyphy-overview-panel';
    this.wrapper.innerHTML = `
      <div class="hyphy-overview-header">
        <div class="overview-left-group">
          <span class="overview-title">Overview</span>
          <span class="overview-status"></span>
        </div>
        <div class="overview-nav-btns">
          <button type="button" class="plot-tool-btn ov-prev" title="Jump to the previous significant site">◀ PREV SIG</button>
          <button type="button" class="plot-tool-btn ov-next" title="Jump to the next significant site">NEXT SIG ▶</button>
        </div>
      </div>
      <div class="hyphy-overview-canvas-container">
        <canvas class="hyphy-overview-canvas"></canvas>
      </div>
    `;
    this.status = this.wrapper.querySelector('.overview-status');
    this.canvas = this.wrapper.querySelector('canvas');
    this.ctx = this.canvas.getContext('2d');
    this.container.appendChild(this.wrapper);

    this.wrapper.querySelector('.ov-prev').addEventListener('click', () => this.jumpSignificant(-1));
    this.wrapper.querySelector('.ov-next').addEventListener('click', () => this.jumpSignificant(1));
  }

  setData(rows, numSites) {
    this.rows = rows || [];
    this.numSites = numSites || this.rows.length;
    this.wrapper.style.display = this.numSites >= this.options.minSites ? '' : 'none';
    this.refresh();
  }

  setThreshold(threshold) {
    this.options.threshold = threshold;
    this.refresh();
  }

  selectSite(site) {
    this.selectedSite = site;
    this.render();
  }

  // The window the main track shows, in 0-based sites with end exclusive
  setViewport(start, end) {
    this.viewport = { start, end };
    this.updateStatus();
    this.render();
  }

  refresh() {
    const t = this.options.threshold;
    this.significant = this.rows.filter(row => this.options.tone(row, t)).map(row => row.site).sort((a, b) => a - b);
    this.dirtyCache = true;
    this.updateStatus();
    this.render();
  }

  updateStatus() {
    const total = this.numSites;
    if (!total) {
      this.status.textContent = '';
      return;
    }
    const start = Math.min(total, Math.max(1, this.viewport.start + 1));
    const end = Math.min(total, Math.max(start, this.viewport.end));
    const pct = Math.round(((end - start + 1) / total) * 100);
    const counts = {};
    for (const row of this.rows) {
      const tone = this.options.tone(row, this.options.threshold);
      if (tone) counts[tone] = (counts[tone] || 0) + 1;
    }
    const parts = [`${this.options.unit} ${start}–${end} of ${total} (${pct}%)`];
    for (const [tone, label] of Object.entries(this.options.toneLabels)) {
      if (counts[tone]) parts.push(`${counts[tone]} ${label}`);
    }
    this.status.textContent = parts.join(' • ');
  }

  jumpSignificant(direction) {
    const site = nextSignificantSite(this.significant, this.selectedSite ?? 0, direction);
    if (site === null) return;
    this.selectSite(site);
    if (this.options.onSiteSelect) {
      this.options.onSiteSelect(site, this.rows.find(row => row.site === site));
    }
  }

  resize(width) {
    if (!width || width <= 0) return;
    const dpr = window.devicePixelRatio || 1;
    this.width = width;
    this.height = this.options.height;
    for (const canvas of [this.canvas, this.offscreenCanvas]) {
      canvas.width = width * dpr;
      canvas.height = this.height * dpr;
    }
    this.canvas.style.width = width + 'px';
    this.canvas.style.height = this.height + 'px';
    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    this.ctx.scale(dpr, dpr);
    this.offscreenCtx.setTransform(1, 0, 0, 1, 0, 0);
    this.offscreenCtx.scale(dpr, dpr);
    this.dirtyCache = true;
    this.render();
  }

  renderOffscreenCache() {
    const ctx = this.offscreenCtx;
    const w = this.width;
    const h = this.height;
    const total = this.numSites;
    const { maxValue, threshold } = this.options;
    const profileH = h - 14;

    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = PALETTE.panel;
    ctx.fillRect(0, 0, w, h);
    ctx.fillStyle = PALETTE.overviewBg;
    ctx.fillRect(0, 0, w, profileH);

    // 1. Significance threshold dashed hairline
    const ruleValue = this.options.rule ? this.options.rule(threshold) : null;
    if (ruleValue !== null) {
      const ruleY = profileH - (Math.min(ruleValue, maxValue) / maxValue) * profileH;
      ctx.strokeStyle = withAlpha(PALETTE.accent, 0.45);
      ctx.lineWidth = 1;
      ctx.setLineDash([2, 3]);
      ctx.beginPath();
      ctx.moveTo(0, ruleY);
      ctx.lineTo(w, ruleY);
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // 2. One bar per site: significant sites in their tone with a landmark pip, the rest as low grey sticks
    const siteW = Math.max(1, w / total);
    for (const row of this.rows) {
      const value = this.options.value(row);
      if (!Number.isFinite(value)) continue;
      const x = ((row.site - 1) / total) * w;
      const barH = Math.max(1.5, (Math.min(Math.max(value, 0), maxValue) / maxValue) * profileH);
      const barY = profileH - barH;
      const tone = this.options.tone(row, threshold);

      if (tone) {
        ctx.fillStyle = tone === 'purifying' ? PALETTE.purifying : PALETTE.accent;
        ctx.fillRect(x, barY, Math.max(1.5, siteW), barH);
        ctx.beginPath();
        ctx.arc(x + siteW / 2, barY - 2, tone === 'purifying' ? 2 : 2.5, 0, Math.PI * 2);
        ctx.fill();
      } else {
        ctx.fillStyle = PALETTE.overviewEmpty;
        ctx.fillRect(x, profileH - Math.max(1, barH * 0.4), Math.max(1, siteW * 0.7), Math.max(1, barH * 0.4));
      }
    }

    // 3. Profile baseline rule
    ctx.strokeStyle = PALETTE.hairline;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, profileH);
    ctx.lineTo(w, profileH);
    ctx.stroke();

    // 4. Site ruler along the bottom
    const rulerY = h - 3;
    ctx.font = '8px "SF Mono", Monaco, sans-serif';
    ctx.fillStyle = PALETTE.stone;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    ctx.strokeStyle = PALETTE.dotMuted;
    ctx.lineWidth = 0.75;
    const step = total > 500 ? 100 : (total > 200 ? 50 : 25);
    for (let s = 1; s <= total; s += step) {
      const x = ((s - 0.5) / total) * w;
      ctx.fillText(s.toString(), x, rulerY);
      ctx.beginPath();
      ctx.moveTo(x, rulerY - 10);
      ctx.lineTo(x, rulerY - 7);
      ctx.stroke();
    }

    this.dirtyCache = false;
  }

  // Brush extent in pixels; never narrower than its two grips
  brushExtent() {
    const total = this.numSites;
    const startX = (this.viewport.start / total) * this.width;
    const endX = (Math.min(total, this.viewport.end) / total) * this.width;
    return { startX, lensW: Math.max(6, endX - startX) };
  }

  render() {
    if (!this.ctx || !this.width || !this.numSites) return;
    const ctx = this.ctx;
    const w = this.width;
    const h = this.height;
    const total = this.numSites;

    if (this.dirtyCache) this.renderOffscreenCache();
    ctx.clearRect(0, 0, w, h);
    ctx.drawImage(this.offscreenCanvas, 0, 0, w, h);

    // 1. Brush body, border and grips
    const { startX, lensW } = this.brushExtent();
    ctx.fillStyle = withAlpha(PALETTE.accent, 0.10);
    ctx.fillRect(startX, 0, lensW, h);
    ctx.strokeStyle = PALETTE.accent;
    ctx.lineWidth = 1.5;
    ctx.strokeRect(startX + 0.5, 0.5, lensW - 1, h - 1);
    ctx.fillStyle = PALETTE.accent;
    ctx.fillRect(startX - 1, 0, 2.5, h);
    ctx.fillRect(startX + lensW - 1.5, 0, 2.5, h);

    // 2. Selected site needle with a pointer cap
    if (this.selectedSite >= 1 && this.selectedSite <= total) {
      const needleX = ((this.selectedSite - 0.5) / total) * w;
      ctx.strokeStyle = PALETTE.accent;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(needleX, 0);
      ctx.lineTo(needleX, h);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(needleX - 3.5, 0);
      ctx.lineTo(needleX + 3.5, 0);
      ctx.lineTo(needleX, 5);
      ctx.closePath();
      ctx.fill();
    }

    // 3. Hovered site guide
    if (this.hoveredSite >= 1 && this.hoveredSite <= total && !this.drag) {
      const hx = ((this.hoveredSite - 0.5) / total) * w;
      ctx.strokeStyle = PALETTE.ink;
      ctx.lineWidth = 1;
      ctx.setLineDash([1, 2]);
      ctx.beginPath();
      ctx.moveTo(hx, 0);
      ctx.lineTo(hx, h);
      ctx.stroke();
      ctx.setLineDash([]);
    }
  }

  siteAt(x) {
    const frac = Math.max(0, Math.min(1, x / this.width));
    return Math.min(this.numSites, Math.max(1, Math.floor(frac * this.numSites) + 1));
  }

  moveViewport(start) {
    const span = this.viewport.end - this.viewport.start;
    const clamped = Math.max(0, Math.min(this.numSites - span, start));
    this.setViewport(clamped, clamped + span);
    if (this.options.onViewportChange) this.options.onViewportChange(clamped, clamped + span);
  }

  setupEvents() {
    const pointerX = (e) => e.clientX - this.canvas.getBoundingClientRect().left;

    this.canvas.addEventListener('mousemove', (e) => {
      if (!this.numSites) return;
      const mx = pointerX(e);
      if (this.drag) {
        this.moveViewport(this.drag.start + Math.round(((mx - this.drag.x) / this.width) * this.numSites));
        return;
      }
      const site = this.siteAt(mx);
      if (this.hoveredSite !== site) {
        this.hoveredSite = site;
        this.render();
      }
    });

    this.canvas.addEventListener('mouseleave', () => {
      if (!this.drag) {
        this.hoveredSite = null;
        this.render();
      }
    });

    // Inside the brush starts a drag; anywhere else centres the brush there and selects the site
    this.canvas.addEventListener('mousedown', (e) => {
      if (!this.numSites) return;
      const mx = pointerX(e);
      const { startX, lensW } = this.brushExtent();
      if (mx >= startX && mx <= startX + lensW) {
        this.drag = { x: mx, start: this.viewport.start };
        return;
      }
      const site = this.siteAt(mx);
      this.moveViewport(site - Math.floor((this.viewport.end - this.viewport.start) / 2));
      this.selectSite(site);
      if (this.options.onSiteSelect) {
        this.options.onSiteSelect(site, this.rows.find(row => row.site === site));
      }
    });

    this.endDrag = () => {
      if (this.drag) {
        this.drag = null;
        this.render();
      }
    };
    window.addEventListener('mouseup', this.endDrag);
  }

  destroy() {
    window.removeEventListener('mouseup', this.endDrag);
    this.wrapper.remove();
  }
}
