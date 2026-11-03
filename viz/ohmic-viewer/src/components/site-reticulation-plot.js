/**
 * Site-Level Network Support Profile (Delta ln L_s = ln L_s^Ohm - ln L_s^Tree).
 * Shows which alignment sites favor the Wheatstone 4-way network shunts (Braun Orange)
 * versus the single binary tree backbone (Steel Blue).
 */

export class SiteReticulationPlot {
  constructor(container) {
    this.container = container;
    this.model = null;
    this.hoveredSite = null;

    this._initDOM();
  }

  _initDOM() {
    this.container.innerHTML = `
      <div class="ohmic-panel-header">
        <div class="ohmic-panel-title-group">
          <span class="ohmic-panel-title">SITE-LEVEL WHEATSTONE NETWORK SUPPORT PROFILE (Δ LN L_s)</span>
          <span class="ohmic-panel-subtitle" id="ohmic-site-subtitle">Per-site log-likelihood gain: ln L_s(Ohmic Network) − ln L_s(Single Tree)</span>
        </div>
        <div class="ohmic-panel-controls">
          <span class="ohmic-legend-item"><span class="ohmic-swatch ohmic-swatch-orange"></span>Favors Wheatstone Network (ΔlnL_s &gt; 0)</span>
          <span class="ohmic-legend-item"><span class="ohmic-swatch ohmic-swatch-blue"></span>Favors Single Tree (ΔlnL_s &lt; 0)</span>
          <button class="ohmic-icon-btn" id="ohmic-site-export-png" title="Export Site Profile as PNG">PNG</button>
        </div>
      </div>
      <div class="ohmic-canvas-wrap ohmic-site-canvas-wrap" id="ohmic-site-canvas-wrap">
        <canvas id="ohmic-site-canvas"></canvas>
        <div class="ohmic-tooltip" id="ohmic-site-tooltip" style="display:none;"></div>
      </div>
    `;

    this.canvasWrap = this.container.querySelector('#ohmic-site-canvas-wrap');
    this.canvas = this.container.querySelector('#ohmic-site-canvas');
    this.tooltip = this.container.querySelector('#ohmic-site-tooltip');

    this.container.querySelector('#ohmic-site-export-png')?.addEventListener('click', () => {
      if (!this.canvas) return;
      const a = document.createElement('a');
      a.href = this.canvas.toDataURL('image/png');
      a.download = 'ohmic_site_delta_logl.png';
      a.click();
    });

    this.canvas.addEventListener('mousemove', (e) => this._onMouseMove(e));
    this.canvas.addEventListener('mouseleave', () => {
      this.hoveredSite = null;
      this.tooltip.style.display = 'none';
      this.draw();
    });

    this._lastObservedWidth = 0;
    this._resizeObserver = new ResizeObserver((entries) => {
      const w = Math.round(entries[0]?.contentRect?.width || 0);
      if (w > 0 && w !== this._lastObservedWidth) {
        this._lastObservedWidth = w;
        this.draw();
      }
    });
    this._resizeObserver.observe(this.canvasWrap);
  }

  setModel(model) {
    this.model = model;
    this.draw();
  }

  draw() {
    if (!this.canvas || !this.model) return;
    const sites = this.model.siteResults || [];
    const rect = this.canvasWrap.getBoundingClientRect();
    const width = Math.max(320, Math.round(rect.width || 800));
    const height = 165;

    const dpr = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(width * dpr);
    this.canvas.height = Math.round(height * dpr);
    this.canvas.style.width = `${width}px`;
    this.canvas.style.height = `${height}px`;

    const ctx = this.canvas.getContext('2d');
    ctx.save();
    ctx.scale(dpr, dpr);

    const isDark = this.container.closest('[theme="dark"]') !== null;
    const bg = isDark ? '#181816' : '#FFFFFF';
    const ink = isDark ? '#F2F2EF' : '#1A1A18';
    const muted = isDark ? '#8E8E88' : '#767672';
    const grid = isDark ? '#2A2A26' : '#EAEAE4';

    ctx.fillStyle = bg;
    ctx.fillRect(0, 0, width, height);

    if (sites.length === 0) {
      ctx.fillStyle = muted;
      ctx.font = '500 11px "SF Mono", Menlo, monospace';
      ctx.textAlign = 'center';
      ctx.fillText('No site-level log-likelihood profile available', width / 2, height / 2);
      ctx.restore();
      return;
    }

    const padLeft = 52;
    const padRight = 20;
    const padTop = 18;
    const padBottom = 26;
    const plotW = Math.max(100, width - padLeft - padRight);
    const plotH = Math.max(60, height - padTop - padBottom);
    const midY = padTop + plotH * 0.5;

    const maxAbs = Math.max(0.05, this.model.maxSiteDeltaLogL * 1.08);
    this.plotGeom = { padLeft, padRight, padTop, padBottom, plotW, plotH, midY, maxAbs, count: sites.length };

    // Horizontal reference lines (+maxAbs/2, 0, -maxAbs/2)
    ctx.strokeStyle = grid;
    ctx.lineWidth = 1;
    [-1, -0.5, 0, 0.5, 1].forEach((frac) => {
      const y = midY - frac * (plotH * 0.46);
      ctx.beginPath();
      ctx.moveTo(padLeft, y);
      ctx.lineTo(padLeft + plotW, y);
      ctx.stroke();

      if (frac === -1 || frac === 0 || frac === 1) {
        ctx.fillStyle = muted;
        ctx.font = '500 9.5px "SF Mono", Menlo, monospace';
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        const val = frac * maxAbs;
        ctx.fillText((val > 0 ? '+' : '') + val.toFixed(2), padLeft - 6, y);
      }
    });

    // Zero baseline
    ctx.strokeStyle = isDark ? '#585852' : '#9E9E98';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    ctx.moveTo(padLeft, midY);
    ctx.lineTo(padLeft + plotW, midY);
    ctx.stroke();

    // Bars for each site
    const S = sites.length;
    const stepX = plotW / S;
    const barW = Math.max(1.2, Math.min(8, stepX * 0.78));

    for (let i = 0; i < S; i++) {
      const d = sites[i].deltaLogL;
      const x = padLeft + (i + 0.5) * stepX;
      const barH = (Math.abs(d) / maxAbs) * (plotH * 0.46);

      ctx.fillStyle = d >= 0 ? '#EB5E28' : '#1E40AF';
      if (this.hoveredSite && this.hoveredSite.site === i + 1) {
        ctx.fillStyle = ink;
      }

      if (d >= 0) {
        ctx.fillRect(x - barW * 0.5, midY - barH, barW, Math.max(1, barH));
      } else {
        ctx.fillRect(x - barW * 0.5, midY, barW, Math.max(1, barH));
      }
    }

    // X-axis ticks
    ctx.fillStyle = muted;
    ctx.font = '500 9.5px "SF Mono", Menlo, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    const tickCount = Math.min(10, S);
    for (let t = 0; t <= tickCount; t++) {
      const siteIdx = Math.min(S, Math.max(1, Math.round(1 + (t / tickCount) * (S - 1))));
      const tx = padLeft + ((siteIdx - 0.5) / S) * plotW;
      ctx.fillText(`Site ${siteIdx}`, tx, padTop + plotH + 6);
    }

    ctx.restore();
  }

  _onMouseMove(e) {
    if (!this.plotGeom || !this.model) return;
    const sites = this.model.siteResults || [];
    if (sites.length === 0) return;

    const rect = this.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const { padLeft, plotW, count } = this.plotGeom;

    const idx = Math.floor(((mx - padLeft) / plotW) * count);
    if (idx >= 0 && idx < count) {
      const s = sites[idx];
      this.hoveredSite = s;
      this.draw();

      this.tooltip.style.display = 'block';
      this.tooltip.style.left = `${Math.min(mx + 14, this.canvas.clientWidth - 230)}px`;
      this.tooltip.style.top = `${Math.max(8, my - 20)}px`;
      this.tooltip.innerHTML = `
        <div class="ohmic-tt-title">SITE ${s.site}</div>
        <div class="ohmic-tt-row"><span>Δ ln L_s (Ohmic − Tree):</span><strong style="color:${s.deltaLogL >= 0 ? '#EB5E28' : '#1E40AF'}">${s.deltaLogL >= 0 ? '+' : ''}${s.deltaLogL.toFixed(4)} nats</strong></div>
        <div class="ohmic-tt-row"><span>Ohmic ln L_s:</span><strong>${s.ohmicLogL.toFixed(4)}</strong></div>
        <div class="ohmic-tt-row"><span>Single-Tree ln L_s:</span><strong>${s.singleTreeLogL.toFixed(4)}</strong></div>
      `;
    } else {
      this.hoveredSite = null;
      this.tooltip.style.display = 'none';
      this.draw();
    }
  }
}
