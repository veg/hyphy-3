// High-DPI Interactive Evidence Plot (Manhattan -log10 p)
// Pure HTML5 Canvas with Retina scaling, responsive resize, threshold guide, hover tooltips, and click inspection

import { PALETTE } from '../data/palette.js';
import { formatPValue, formatThreshold, formatOmegaInf } from '../core/format.js';
import { downloadPNG, downloadSVG } from '../core/export.js';

export class EvidencePlot {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      threshold: 0.1,
      plotMode: 'manhattan', // 'manhattan' | 'signed' | 'diff'
      onSiteSelect: null,
      onSiteHover: null
    }, options);

    this.threshold = this.options.threshold;
    this.plotMode = this.options.plotMode;
    this.data = null;
    this.selectedSite = null;
    this.hoveredSite = null;
    this.points = []; // Cached screen coordinates

    this.initDOM();
    this.initEvents();
    this.initResizeObserver();
  }

  initDOM() {
    this.wrapper = document.createElement('div');
    this.wrapper.className = 'hyphy-evidence-plot-card';
    this.wrapper.innerHTML = `
      <div class="evidence-plot-header">
        <div class="plot-title-group">
          <span class="plot-title">Selection Evidence Landscape</span>
          <span class="plot-subtitle" id="plot-sub-label">−log₁₀(p) vs Codon Position</span>
        </div>
        <div class="plot-actions">
          <div class="segmented-control plot-mode-control">
            <button type="button" class="seg-btn active" data-mode="manhattan">Evidence (−log₁₀ p)</button>
            <button type="button" class="seg-btn" data-mode="signed">Signed Evidence</button>
            <button type="button" class="seg-btn" data-mode="diff">Rate Diff (β − α)</button>
          </div>
          <div class="plot-export-group">
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
        this.setPlotMode(btn.dataset.mode);
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

  setPlotMode(mode) {
    this.plotMode = mode;
    if (mode === 'manhattan') {
      this.subLabel.textContent = this.data?.isMEME
        ? '−log₁₀(p) Evidence for Episodic Diversifying Selection'
        : '−log₁₀(p) vs Codon Position';
    } else if (mode === 'signed') {
      this.subLabel.textContent = this.data?.isMEME
        ? 'Episodic Rate Ratio: β⁺ / α'
        : 'Signed Evidence: +Diversifying / −Purifying';
    } else {
      this.subLabel.textContent = this.data?.isMEME
        ? 'Episodic Rate Difference: β⁺ − α'
        : 'Absolute Rate Difference: β (dN) − α (dS)';
    }
    this.render();
  }

  selectSite(siteNumber) {
    this.selectedSite = siteNumber;
    this.render();
  }

  initEvents() {
    this.canvas.addEventListener('mousemove', (e) => this.handleMouseMove(e));
    this.canvas.addEventListener('mouseleave', () => this.handleMouseLeave());
    this.canvas.addEventListener('click', (e) => this.handleClick(e));
  }

  handleMouseMove(e) {
    if (!this.data || this.points.length === 0) return;
    const rect = this.canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    // Find closest site by X coordinate
    let closest = null;
    let minDx = Infinity;

    for (const pt of this.points) {
      const dx = Math.abs(pt.x - mouseX);
      if (dx < minDx) {
        minDx = dx;
        closest = pt;
      }
    }

    if (closest && minDx <= Math.max(12, this.plotWidth / Math.max(1, this.data.numSites) * 2)) {
      this.hoveredSite = closest.site;
      this.render();
      this.showTooltip(closest, e.clientX, e.clientY);
    } else {
      this.hoveredSite = null;
      this.render();
      this.hideTooltip();
    }
  }

  handleMouseLeave() {
    this.hoveredSite = null;
    this.render();
    this.hideTooltip();
  }

  handleClick(e) {
    if (!this.data || this.points.length === 0) return;
    const rect = this.canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;

    let closest = null;
    let minDx = Infinity;

    for (const pt of this.points) {
      const dx = Math.abs(pt.x - mouseX);
      if (dx < minDx) {
        minDx = dx;
        closest = pt;
      }
    }

    if (closest && minDx <= Math.max(16, this.plotWidth / Math.max(1, this.data.numSites) * 2)) {
      this.selectedSite = closest.site;
      this.render();
      if (typeof this.options.onSiteSelect === 'function') {
        this.options.onSiteSelect(closest.site, closest.raw);
      }
    }
  }

  showTooltip(pt, clientX, clientY) {
    const r = pt.raw;
    const isSig = r.pValue <= this.threshold;
    const isDiv = this.data.isMEME ? isSig : (isSig && r.beta > r.alpha);
    const isPur = !this.data.isMEME && isSig && r.alpha > r.beta;

    let selBadge = '<span class="tip-badge tip-neutral">Neutral</span>';
    if (isDiv) {
      selBadge = this.data.isMEME
        ? '<span class="tip-badge tip-div">▲ Episodic Diversifying</span>'
        : '<span class="tip-badge tip-div">▲ Diversifying Selection</span>';
    } else if (isPur) {
      selBadge = '<span class="tip-badge tip-pur">▼ Purifying Selection</span>';
    }

    const pStr = formatPValue(r.pValue);
    const subStr = r.history ? `<div class="tip-row tip-sub"><span class="tip-lbl">Substitutions:</span> <span class="tip-val-sub">${r.history}</span></div>` : '';

    let statsRows = '';
    if (this.data.isMEME) {
      const betaPlus = r.beta_plus ?? r.beta;
      const pPlusPct = Math.round((r.p_plus ?? 0) * 100);
      const branches = r.branchesUnderSelection || (r.branchEBF ? Object.values(r.branchEBF).filter(v => v >= 100).length : 0);
      statsRows = `
        <div class="tip-row"><span class="tip-lbl">p-value:</span> <span class="tip-val font-mono">${pStr}</span></div>
        <div class="tip-row"><span class="tip-lbl">LRT:</span> <span class="tip-val font-mono">${r.lrt.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">α (dS):</span> <span class="tip-val font-mono">${r.alpha.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">β⁺ (dN):</span> <span class="tip-val font-mono">${betaPlus.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">p⁺ (% branches):</span> <span class="tip-val font-mono">${pPlusPct}%</span></div>
        <div class="tip-row"><span class="tip-lbl">Branches:</span> <span class="tip-val font-mono">${branches} (EBF≥100)</span></div>
      `;
    } else {
      const omegaStr = formatOmegaInf(r.omega);
      statsRows = `
        <div class="tip-row"><span class="tip-lbl">p-value:</span> <span class="tip-val font-mono">${pStr}</span></div>
        <div class="tip-row"><span class="tip-lbl">LRT:</span> <span class="tip-val font-mono">${r.lrt.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">α (dS):</span> <span class="tip-val font-mono">${r.alpha.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">β (dN):</span> <span class="tip-val font-mono">${r.beta.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">dN/dS (ω):</span> <span class="tip-val font-mono">${omegaStr}</span></div>
      `;
    }

    this.tooltip.innerHTML = `
      <div class="tip-header">
        <span class="tip-title">Codon #${r.site}</span>
        ${selBadge}
      </div>
      <div class="tip-grid">
        ${statsRows}
      </div>
      ${subStr}
      <div class="tip-hint">Click point to inspect site tree &amp; branches</div>
    `;

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

  getYValue(r) {
    if (this.plotMode === 'manhattan') {
      return r.negLogP;
    } else if (this.plotMode === 'signed') {
      if (this.data?.isMEME) {
        const bp = r.beta_plus ?? r.beta;
        return r.alpha > 1e-4 ? (bp / r.alpha) : (bp > 0 ? 10 : 1);
      }
      const sign = r.beta >= r.alpha ? 1 : -1;
      return sign * r.negLogP;
    } else {
      // diff: beta - alpha
      const b = this.data?.isMEME ? (r.beta_plus ?? r.beta) : r.beta;
      return b - r.alpha;
    }
  }

  render() {
    if (!this.data || !this.width || !this.height) return;
    const ctx = this.ctx;
    const dpr = this.dpr;
    const w = this.width;
    const h = this.height;

    ctx.save();
    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, w, h);

    // Margins
    const padLeft = 52;
    const padRight = 24;
    const padTop = 24;
    const padBottom = 38;

    const plotW = w - padLeft - padRight;
    const plotH = h - padTop - padBottom;
    this.plotWidth = plotW;
    this.plotHeight = plotH;

    const numSites = Math.max(1, this.data.numSites);
    const results = this.data.siteResults;

    // Determine Y range
    let yMin = 0;
    let yMax = 3.0; // minimum ceiling for -log10 p

    if (this.plotMode === 'manhattan') {
      for (const r of results) {
        if (r.negLogP > yMax) yMax = r.negLogP;
      }
      yMax = Math.ceil(yMax * 1.15);
      yMin = 0;
    } else if (this.plotMode === 'signed') {
      let maxAbs = 2.0;
      for (const r of results) {
        const v = Math.abs(r.negLogP);
        if (v > maxAbs) maxAbs = v;
      }
      maxAbs = Math.ceil(maxAbs * 1.15);
      yMax = maxAbs;
      yMin = -maxAbs;
    } else {
      // Rate diff
      let minD = 0;
      let maxD = 1.0;
      for (const r of results) {
        const diff = r.beta - r.alpha;
        if (diff < minD) minD = diff;
        if (diff > maxD) maxD = diff;
      }
      const margin = Math.max(1.0, (maxD - minD) * 0.1);
      yMin = Math.floor(minD - margin);
      yMax = Math.ceil(maxD + margin);
    }

    const scaleX = (siteIdx) => padLeft + ((siteIdx - 1) / Math.max(1, numSites - 1)) * plotW;
    const scaleY = (val) => padTop + (1 - (val - yMin) / (yMax - yMin)) * plotH;

    // 1. Grid Lines & Ticks (Dieter Rams subtle gray)
    ctx.strokeStyle = PALETTE.gridLine;
    ctx.lineWidth = 1;
    ctx.font = '10px "SF Mono", Monaco, monospace';
    ctx.fillStyle = PALETTE.stone;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';

    // Horizontal grid ticks
    const numYTicks = 5;
    for (let i = 0; i <= numYTicks; ++i) {
      const val = yMin + (i / numYTicks) * (yMax - yMin);
      const y = scaleY(val);
      ctx.beginPath();
      ctx.moveTo(padLeft, y);
      ctx.lineTo(padLeft + plotW, y);
      ctx.stroke();

      ctx.fillText(val.toFixed(this.plotMode === 'diff' ? 1 : 1), padLeft - 8, y);
    }

    // Zero line for signed or diff mode
    if (yMin < 0 && yMax > 0) {
      const zeroY = scaleY(0);
      ctx.strokeStyle = PALETTE.hairline;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(padLeft, zeroY);
      ctx.lineTo(padLeft + plotW, zeroY);
      ctx.stroke();
    }

    // Vertical X-axis ticks (Codon positions)
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    const xStep = Math.pow(10, Math.floor(Math.log10(numSites))) / 2 || 20;
    const actualStep = numSites > 300 ? 50 : (numSites > 100 ? 20 : 10);

    for (let s = 1; s <= numSites; s += actualStep) {
      const x = scaleX(s);
      ctx.strokeStyle = PALETTE.gridLine;
      ctx.beginPath();
      ctx.moveTo(x, padTop);
      ctx.lineTo(x, padTop + plotH);
      ctx.stroke();

      ctx.fillText(String(s), x, padTop + plotH + 8);
    }
    // Last codon tick if not overlapping
    const lastX = scaleX(numSites);
    ctx.fillText(String(numSites), lastX, padTop + plotH + 8);

    // 2. Significance Threshold Line (Dashed Alert Line)
    if (this.plotMode === 'manhattan' || this.plotMode === 'signed') {
      const threshNegLog = -Math.log10(Math.max(1e-15, this.threshold));
      const threshY = scaleY(threshNegLog);

      ctx.save();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = PALETTE.accent; // Braun Signal Orange
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(padLeft, threshY);
      ctx.lineTo(padLeft + plotW, threshY);
      ctx.stroke();

      // Threshold label badge
      ctx.fillStyle = PALETTE.accent;
      ctx.font = 'bold 9px "SF Mono", Monaco, monospace';
      ctx.textAlign = 'left';
      ctx.fillText(`p = ${formatThreshold(this.threshold)}`, padLeft + 6, threshY - 6);

      // If signed, draw negative threshold line for purifying
      if (this.plotMode === 'signed') {
        const negThreshY = scaleY(-threshNegLog);
        ctx.strokeStyle = PALETTE.purifying;
        ctx.beginPath();
        ctx.moveTo(padLeft, negThreshY);
        ctx.lineTo(padLeft + plotW, negThreshY);
        ctx.stroke();
      }
      ctx.restore();
    }

    // 3. Compute screen points & render glyphs
    this.points = [];
    const pointsToRender = [];

    for (let i = 0; i < results.length; ++i) {
      const r = results[i];
      const yVal = this.getYValue(r);
      const x = scaleX(r.site);
      const y = scaleY(yVal);

      const isSig = r.pValue <= this.threshold;
      const isDiv = this.data.isMEME ? isSig : (isSig && r.beta > r.alpha);
      const isPur = !this.data.isMEME && isSig && r.alpha > r.beta;

      const ptObj = {
        site: r.site,
        x,
        y,
        yVal,
        raw: r,
        isSig,
        isDiv,
        isPur
      };

      this.points.push(ptObj);
      pointsToRender.push(ptObj);
    }

    // Sort render order: neutral first, purifying second, diversifying on top!
    pointsToRender.sort((a, b) => {
      const scoreA = a.isDiv ? 2 : (a.isPur ? 1 : 0);
      const scoreB = b.isDiv ? 2 : (b.isPur ? 1 : 0);
      return scoreA - scoreB;
    });

    // Draw Points / Glyphs
    for (const pt of pointsToRender) {
      const isSelected = (this.selectedSite === pt.site);
      const isHovered = (this.hoveredSite === pt.site);

      if (pt.isDiv) {
        // Diversifying: ▲ Upward Triangle (Braun Signal Orange)
        const size = isSelected ? 8 : (isHovered ? 7 : 5.5);
        ctx.fillStyle = PALETTE.accent;
        ctx.strokeStyle = PALETTE.accentDeep;
        ctx.lineWidth = isSelected ? 2 : 1;

        ctx.beginPath();
        ctx.moveTo(pt.x, pt.y - size * 1.2);
        ctx.lineTo(pt.x - size, pt.y + size * 0.8);
        ctx.lineTo(pt.x + size, pt.y + size * 0.8);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
      } else if (pt.isPur) {
        // Purifying: ▼ Downward Triangle (Technical Cobalt)
        const size = isSelected ? 7 : (isHovered ? 6 : 4.5);
        ctx.fillStyle = PALETTE.purifying;
        ctx.strokeStyle = PALETTE.purifyingDeep;
        ctx.lineWidth = isSelected ? 2 : 1;

        ctx.beginPath();
        ctx.moveTo(pt.x, pt.y + size * 1.2);
        ctx.lineTo(pt.x - size, pt.y - size * 0.8);
        ctx.lineTo(pt.x + size, pt.y - size * 0.8);
        ctx.closePath();
        ctx.fill();
        ctx.stroke();
      } else {
        // Neutral: subtle gray dot
        const rSize = isSelected ? 5 : (isHovered ? 4 : 2.5);
        ctx.fillStyle = isHovered ? PALETTE.slate : PALETTE.slateLight;
        ctx.beginPath();
        ctx.arc(pt.x, pt.y, rSize, 0, Math.PI * 2);
        ctx.fill();
      }

      // Draw Selected Reticle Ring
      if (isSelected) {
        ctx.strokeStyle = PALETTE.amberVivid; // Vivid Amber
        ctx.lineWidth = 2;
        ctx.beginPath();
        ctx.arc(pt.x, pt.y, 11, 0, Math.PI * 2);
        ctx.stroke();
      }
    }

    // 4. Hover Hairline Crosshair
    if (this.hoveredSite) {
      const hPt = this.points.find(p => p.site === this.hoveredSite);
      if (hPt) {
        ctx.save();
        ctx.setLineDash([2, 2]);
        ctx.strokeStyle = PALETTE.ink;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(hPt.x, padTop);
        ctx.lineTo(hPt.x, padTop + plotH);
        ctx.stroke();
        ctx.restore();
      }
    }

    ctx.restore();
  }

  exportPNG() {
    downloadPNG(this.canvas, 'hyphy_fel_evidence_plot.png');
  }

  exportSVG() {
    if (!this.data) return;
    const w = this.width;
    const h = this.height;
    const numSites = this.data.numSites;
    const padLeft = 52;
    const padRight = 24;
    const padTop = 24;
    const padBottom = 38;
    const plotW = w - padLeft - padRight;
    const plotH = h - padTop - padBottom;

    let yMax = 3.0;
    for (const r of this.data.siteResults) {
      if (r.negLogP > yMax) yMax = r.negLogP;
    }
    yMax = Math.ceil(yMax * 1.15);

    const scaleX = (s) => padLeft + ((s - 1) / Math.max(1, numSites - 1)) * plotW;
    const scaleY = (v) => padTop + (1 - v / yMax) * plotH;

    let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">\n`;
    svg += `  <style>\n    text { font-family: -apple-system, sans-serif; font-size: 10px; fill: ${PALETTE.stone}; }\n  </style>\n`;
    svg += `  <rect width="100%" height="100%" fill="${PALETTE.panel}"/>\n`;

    // Threshold line
    const threshY = scaleY(-Math.log10(this.threshold));
    svg += `  <line x1="${padLeft}" y1="${threshY}" x2="${padLeft + plotW}" y2="${threshY}" stroke="${PALETTE.accent}" stroke-dasharray="4,4" stroke-width="1.5"/>\n`;
    svg += `  <text x="${padLeft + 6}" y="${threshY - 4}" fill="${PALETTE.accent}" font-weight="bold">p = ${this.threshold}</text>\n`;

    // Points
    for (const pt of this.points) {
      if (pt.isDiv) {
        svg += `  <polygon points="${pt.x},${pt.y - 6} ${pt.x - 5},${pt.y + 4} ${pt.x + 5},${pt.y + 4}" fill="${PALETTE.accent}" stroke="${PALETTE.accentDeep}"/>\n`;
      } else if (pt.isPur) {
        svg += `  <polygon points="${pt.x},${pt.y + 6} ${pt.x - 5},${pt.y - 4} ${pt.x + 5},${pt.y - 4}" fill="${PALETTE.purifying}" stroke="${PALETTE.purifyingDeep}"/>\n`;
      } else {
        svg += `  <circle cx="${pt.x}" cy="${pt.y}" r="2.5" fill="${PALETTE.slateLight}"/>\n`;
      }
    }

    svg += `</svg>`;
    downloadSVG(svg, 'hyphy_fel_evidence_plot.svg');
  }
}
