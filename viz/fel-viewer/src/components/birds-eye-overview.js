// Bird's Eye Overview: Whole-Gene Selection Landscape & 2D Alignment Minimap
// Dieter Rams Aesthetic: Calibrated frequency tuner dial with dual 1D profile and 2D matrix modes

import { ZAPPO_COLORS } from '../data/genetic-code.js';
import { PALETTE, withAlpha } from '../data/palette.js';

export class BirdsEyeOverview {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      mode: '1d',                   // '1d' (56px tuner profile) | '2d' (136px macro matrix)
      height1D: 56,
      height2D: 136,
      pThreshold: 0.1,
      displayMode: 'amino_acid',
      colorScheme: 'difference',
      onSiteSelect: null,
      onViewportChange: null        // (scrollLeft, scrollTop)
    }, options);

    this.data = null;
    this.taxaOrder = [];
    this.selectedSite = 1;
    this.hoveredSite = -1;
    this.visibleColRange = { start: 0, end: 10 };
    this.visibleRowRange = { start: 0, end: 10 };
    this.totalCols = 1;
    this.totalRows = 1;

    // Offscreen cache for high performance
    this.offscreenCanvas = document.createElement('canvas');
    this.offscreenCtx = this.offscreenCanvas.getContext('2d');
    this.dirtyCache = true;

    // Interaction state
    this.isDragging = false;
    this.dragStartX = 0;
    this.dragStartY = 0;
    this.dragStartCol = 0;
    this.dragStartRow = 0;

    this.initDOM();
    this.setupEvents();
  }

  initDOM() {
    this.wrapper = document.createElement('div');
    this.wrapper.className = 'hyphy-overview-panel';

    // 1. Control & Status Bar
    this.header = document.createElement('div');
    this.header.className = 'hyphy-overview-header';
    this.header.innerHTML = `
      <div class="overview-left-group">
        <span class="overview-title">Bird's Eye Overview</span>
        <span class="overview-status" id="ov-status">Loading...</span>
      </div>
      <div class="overview-right-group">
        <div class="overview-nav-btns">
          <button class="ov-btn" id="ov-prev" title="Jump to previous significant site (p <= threshold)">◀ PREV SIG</button>
          <button class="ov-btn" id="ov-next" title="Jump to next significant site (p <= threshold)">NEXT SIG ▶</button>
        </div>
        <div class="segmented-control ov-mode-toggle">
          <button class="seg-btn active" data-mode="1d" title="1D Gene Selection Profile & Tuner">⚊ 1D PROFILE</button>
          <button class="seg-btn" data-mode="2d" title="2D Alignment Macro Matrix">☷ 2D MATRIX</button>
        </div>
      </div>
    `;
    this.wrapper.appendChild(this.header);

    // 2. Interactive Canvas
    this.canvasContainer = document.createElement('div');
    this.canvasContainer.className = 'hyphy-overview-canvas-container';
    this.canvas = document.createElement('canvas');
    this.canvas.className = 'hyphy-overview-canvas';
    this.ctx = this.canvas.getContext('2d');
    this.canvasContainer.appendChild(this.canvas);
    this.wrapper.appendChild(this.canvasContainer);

    this.container.appendChild(this.wrapper);

    // Mode Toggle Buttons
    this.header.querySelectorAll('.ov-mode-toggle .seg-btn').forEach(btn => {
      btn.addEventListener('click', (e) => {
        this.header.querySelectorAll('.ov-mode-toggle .seg-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        this.setMode(btn.dataset.mode);
      });
    });

    // Prev / Next Sig Site Buttons
    this.header.querySelector('#ov-prev').addEventListener('click', () => this.jumpSignificant(-1));
    this.header.querySelector('#ov-next').addEventListener('click', () => this.jumpSignificant(1));
  }

  setMode(mode) {
    if (this.options.mode === mode) return;
    this.options.mode = mode;
    this.dirtyCache = true;
    this.resize(this.width);
    if (this.options.onModeChange) {
      this.options.onModeChange(mode);
    }
  }

  setData(felData, taxaOrder = []) {
    this.data = felData;
    this.taxaOrder = taxaOrder.length > 0 ? taxaOrder : (felData?.taxa || []);
    this.totalCols = felData?.numSites || 1;
    this.totalRows = this.taxaOrder.length || 1;
    this.dirtyCache = true;
    this.updateStatus();
    this.render();
  }

  setTaxaOrder(taxaOrder) {
    this.taxaOrder = taxaOrder;
    this.totalRows = this.taxaOrder.length || 1;
    this.dirtyCache = true;
    this.render();
  }

  setColorScheme(scheme) {
    this.options.colorScheme = scheme;
    this.dirtyCache = true;
    this.render();
  }

  setDisplayMode(mode) {
    this.options.displayMode = mode;
    this.dirtyCache = true;
    this.render();
  }

  setPThreshold(th) {
    this.options.pThreshold = th;
    this.dirtyCache = true;
    this.updateStatus();
    this.render();
  }

  setSelectedSite(site) {
    this.selectedSite = site;
    this.render();
  }

  // Update visible range coordinates from main alignment canvas
  setViewport(visibleColStart, visibleColEnd, visibleRowStart, visibleRowEnd, totalCols, totalRows) {
    this.visibleColRange = { start: visibleColStart, end: visibleColEnd };
    this.visibleRowRange = { start: visibleRowStart, end: visibleRowEnd };
    if (totalCols) this.totalCols = totalCols;
    if (totalRows) this.totalRows = totalRows;
    this.updateStatus();
    this.render();
  }

  updateStatus() {
    const statusEl = this.header.querySelector('#ov-status');
    if (!statusEl || !this.data) return;

    const totalSites = this.data.numSites;
    const th = this.options.pThreshold;
    const sigSites = this.data.siteResults.filter(r => r.pValue <= th);
    const startSite = Math.min(totalSites, Math.max(1, this.visibleColRange.start + 1));
    const endSite = Math.min(totalSites, Math.max(startSite, this.visibleColRange.end));
    const pct = totalSites > 0 ? Math.round(((endSite - startSite + 1) / totalSites) * 100) : 100;

    if (this.data.isMEME) {
      const epiCount = sigSites.length;
      statusEl.innerHTML = `
        <span class="status-range">Sites ${startSite}–${endSite} of ${totalSites} (${pct}%)</span>
        <span class="status-divider">•</span>
        <span class="status-div">${epiCount} Episodic Diversifying ▲</span>
      `;
    } else {
      const divCount = sigSites.filter(r => r.beta > r.alpha).length;
      const purCount = sigSites.filter(r => r.alpha > r.beta).length;
      statusEl.innerHTML = `
        <span class="status-range">Sites ${startSite}–${endSite} of ${totalSites} (${pct}%)</span>
        <span class="status-divider">•</span>
        <span class="status-div">${divCount} Diversifying ▲</span>
        <span class="status-divider">•</span>
        <span class="status-pur">${purCount} Purifying ▼</span>
      `;
    }
  }

  jumpSignificant(direction) {
    if (!this.data) return;
    const th = this.options.pThreshold;
    const curSite = this.selectedSite;
    const sites = this.data.siteResults;

    let target = null;
    if (direction > 0) {
      // Find first sig site > curSite
      target = sites.find(r => r.site > curSite && r.pValue <= th);
      if (!target) {
        // Wrap around to first sig site
        target = sites.find(r => r.pValue <= th);
      }
    } else {
      // Find last sig site < curSite
      for (let i = curSite - 2; i >= 0; --i) {
        if (sites[i] && sites[i].pValue <= th) {
          target = sites[i];
          break;
        }
      }
      if (!target) {
        // Wrap around to last sig site
        for (let i = sites.length - 1; i >= 0; --i) {
          if (sites[i] && sites[i].pValue <= th) {
            target = sites[i];
            break;
          }
        }
      }
    }

    if (target && this.options.onSiteSelect) {
      this.options.onSiteSelect(target.site, target);
    }
  }

  resize(width) {
    if (!width || width <= 0) return;
    const dpr = window.devicePixelRatio || 1;
    this.width = width;
    this.height = this.options.mode === '2d' ? this.options.height2D : this.options.height1D;

    this.canvas.width = width * dpr;
    this.canvas.height = this.height * dpr;
    this.canvas.style.width = width + 'px';
    this.canvas.style.height = this.height + 'px';

    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    this.ctx.scale(dpr, dpr);

    this.offscreenCanvas.width = width * dpr;
    this.offscreenCanvas.height = this.height * dpr;
    this.offscreenCtx.setTransform(1, 0, 0, 1, 0, 0);
    this.offscreenCtx.scale(dpr, dpr);

    this.dirtyCache = true;
    this.render();
  }

  renderOffscreenCache() {
    if (!this.data || !this.width || !this.height) return;
    const ctx = this.offscreenCtx;
    const w = this.width;
    const h = this.height;
    const totalSites = this.data.numSites;
    if (totalSites <= 0) return;

    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = PALETTE.panel;
    ctx.fillRect(0, 0, w, h);

    const is2D = this.options.mode === '2d';
    const profileH = is2D ? 40 : h - 14;
    const matrixH = is2D ? h - profileH - 12 : 0;
    const maxNegLogP = 5.0; // scale -log10(p) up to 1e-5

    // 1. Draw Whole-Gene Mutation Divergence Silhouette (Underlay)
    ctx.beginPath();
    ctx.moveTo(0, profileH);
    for (let s = 0; s < totalSites; ++s) {
      const x = (s / totalSites) * w;
      const refAA = this.data.consensusAA[s] || '-';
      let mutatedCount = 0;
      let validCount = 0;
      for (const t of this.taxaOrder) {
        const aa = this.data.translated[t]?.[s] || this.data.ancestralTranslated?.[t]?.[s];
        if (aa && aa !== '-' && aa !== '?') {
          validCount++;
          if (aa !== refAA) mutatedCount++;
        }
      }
      const divRate = validCount > 0 ? (mutatedCount / validCount) : 0;
      const y = profileH - divRate * (profileH * 0.75);
      ctx.lineTo(x, y);
    }
    ctx.lineTo(w, profileH);
    ctx.closePath();
    ctx.fillStyle = PALETTE.overviewBg;
    ctx.fill();
    ctx.strokeStyle = PALETTE.overviewRule;
    ctx.lineWidth = 0.75;
    ctx.stroke();

    // 2. Significance Threshold Dashed Hairline
    const th = this.options.pThreshold;
    const threshNegLog = -Math.log10(Math.max(th, 1e-15));
    const threshY = profileH - (Math.min(threshNegLog, maxNegLogP) / maxNegLogP) * profileH;
    ctx.strokeStyle = withAlpha(PALETTE.accent, 0.45);
    ctx.lineWidth = 1;
    ctx.setLineDash([2, 3]);
    ctx.beginPath();
    ctx.moveTo(0, threshY);
    ctx.lineTo(w, threshY);
    ctx.stroke();
    ctx.setLineDash([]);

    // 3. Selection Landscape Sticks & Landmark Pins
    const siteW = Math.max(1, w / totalSites);
    for (let s = 0; s < totalSites; ++s) {
      const r = this.data.siteResults[s];
      if (!r) continue;
      const x = (s / totalSites) * w;
      const isSig = r.pValue <= th;
      const val = Math.min(r.negLogP, maxNegLogP);
      const barH = Math.max(1.5, (val / maxNegLogP) * profileH);
      const barY = profileH - barH;

      if (isSig) {
        if (this.data.isMEME || r.beta > r.alpha) {
          // Diversifying / Episodic: Braun Signal Orange
          ctx.fillStyle = PALETTE.accent;
          ctx.fillRect(x, barY, Math.max(1.5, siteW), barH);
          // Landmark Pip
          ctx.beginPath();
          ctx.arc(x + siteW / 2, barY - 2, 2.5, 0, Math.PI * 2);
          ctx.fill();
        } else {
          // Purifying: Technical Cobalt Blue
          ctx.fillStyle = PALETTE.purifying;
          ctx.fillRect(x, barY, Math.max(1.5, siteW), barH);
          ctx.beginPath();
          ctx.arc(x + siteW / 2, barY - 2, 2, 0, Math.PI * 2);
          ctx.fill();
        }
      } else {
        // Neutral baseline stick
        ctx.fillStyle = PALETTE.overviewEmpty;
        ctx.fillRect(x, profileH - Math.max(1, barH * 0.4), Math.max(1, siteW * 0.7), Math.max(1, barH * 0.4));
      }
    }

    // 4. Subtle Profile Baseline Rule
    ctx.strokeStyle = PALETTE.hairline;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, profileH);
    ctx.lineTo(w, profileH);
    ctx.stroke();

    // 5. If 2D Mode: Draw Compressed Multi-Sequence Alignment Thumbnail
    if (is2D && matrixH > 10) {
      const matrixTop = profileH + 2;
      const numTaxa = this.taxaOrder.length;
      const rowStep = matrixH / Math.max(1, numTaxa);

      // Light backdrop
      ctx.fillStyle = PALETTE.overviewRow;
      ctx.fillRect(0, matrixTop, w, matrixH);

      const scheme = this.options.colorScheme;
      for (let r = 0; r < numTaxa; ++r) {
        const taxon = this.taxaOrder[r];
        const isAncestor = this.data.ancestralSequences?.[taxon] !== undefined;
        const aaSeq = isAncestor 
          ? (this.data.ancestralTranslated[taxon] || [])
          : (this.data.translated[taxon] || []);
        const y = matrixTop + r * rowStep;

        for (let s = 0; s < totalSites; ++s) {
          const aa = aaSeq[s] || '-';
          const refAA = scheme === 'difference_root'
            ? (this.data.rootAA?.[s] || this.data.consensusAA[s] || '-')
            : (this.data.consensusAA[s] || '-');

          if (aa !== refAA && aa !== '-') {
            const x = (s / totalSites) * w;
            ctx.fillStyle = ZAPPO_COLORS[aa] || PALETTE.accent;
            ctx.fillRect(x, y, Math.max(1, siteW), Math.max(1, rowStep));
          }
        }
      }

      // Matrix bottom hairline
      ctx.strokeStyle = PALETTE.hairline;
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(0, matrixTop + matrixH);
      ctx.lineTo(w, matrixTop + matrixH);
      ctx.stroke();
    }

    // 6. Calibrated Site Ruler (at the very bottom)
    const rulerY = h - 3;
    ctx.font = '8px "SF Mono", Monaco, sans-serif';
    ctx.fillStyle = PALETTE.stone;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';

    const step = totalSites > 500 ? 100 : (totalSites > 200 ? 50 : 25);
    for (let s = 1; s <= totalSites; s += step) {
      const x = ((s - 0.5) / totalSites) * w;
      ctx.fillText(s.toString(), x, rulerY);
      ctx.strokeStyle = PALETTE.dotMuted;
      ctx.lineWidth = 0.75;
      ctx.beginPath();
      ctx.moveTo(x, rulerY - 10);
      ctx.lineTo(x, rulerY - 7);
      ctx.stroke();
    }

    this.dirtyCache = false;
  }

  render() {
    if (!this.ctx || !this.width || !this.height) return;
    const ctx = this.ctx;
    const w = this.width;
    const h = this.height;
    const totalSites = this.data?.numSites || 1;
    const totalRows = this.totalRows || 1;
    const is2D = this.options.mode === '2d';

    if (this.dirtyCache) {
      this.renderOffscreenCache();
    }

    // Draw cached static landscape/matrix
    ctx.clearRect(0, 0, w, h);
    ctx.drawImage(this.offscreenCanvas, 0, 0);

    // 1. Draw Viewport Lens / Caliper (Interactive Layer)
    const siteW = w / totalSites;
    const startX = (this.visibleColRange.start / totalSites) * w;
    const endX = (Math.min(totalSites, this.visibleColRange.end) / totalSites) * w;
    const lensW = Math.max(6, endX - startX);

    let lensY = 0;
    let lensH = h;
    if (is2D) {
      const profileH = 40;
      const matrixH = h - profileH - 12;
      const rowStep = matrixH / Math.max(1, totalRows);
      const startY = profileH + 2 + (this.visibleRowRange.start / totalRows) * matrixH;
      const endY = profileH + 2 + (Math.min(totalRows, this.visibleRowRange.end) / totalRows) * matrixH;
      lensY = startY;
      lensH = Math.max(10, endY - startY);

      // In 2D mode, also draw horizontal lens stripe across the 1D profile
      ctx.fillStyle = withAlpha(PALETTE.accent, 0.08);
      ctx.fillRect(startX, 0, lensW, profileH);
      ctx.strokeStyle = withAlpha(PALETTE.accent, 0.4);
      ctx.lineWidth = 1;
      ctx.strokeRect(startX + 0.5, 0.5, lensW - 1, profileH - 1);
    }

    // Caliper Lens Body
    ctx.fillStyle = withAlpha(PALETTE.accent, 0.10);
    ctx.fillRect(startX, lensY, lensW, lensH);

    // Caliper Machined Border (Braun Signal Orange)
    ctx.strokeStyle = PALETTE.accent;
    ctx.lineWidth = 1.5;
    ctx.strokeRect(startX + 0.5, lensY + 0.5, lensW - 1, lensH - 1);

    // Caliper Left & Right Gripper Handles
    ctx.fillStyle = PALETTE.accent;
    ctx.fillRect(startX - 1, lensY, 2.5, lensH);
    ctx.fillRect(startX + lensW - 1.5, lensY, 2.5, lensH);

    // 2. Active Site Needle (Precision hair needle with pointer cap)
    if (this.selectedSite >= 1 && this.selectedSite <= totalSites) {
      const needleX = ((this.selectedSite - 0.5) / totalSites) * w;
      ctx.strokeStyle = PALETTE.accent;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(needleX, 0);
      ctx.lineTo(needleX, h);
      ctx.stroke();

      // Top triangle pointer
      ctx.fillStyle = PALETTE.accent;
      ctx.beginPath();
      ctx.moveTo(needleX - 3.5, 0);
      ctx.lineTo(needleX + 3.5, 0);
      ctx.lineTo(needleX, 5);
      ctx.closePath();
      ctx.fill();
    }

    // 3. Hovered Site Guide & Readout
    if (this.hoveredSite >= 1 && this.hoveredSite <= totalSites && !this.isDragging) {
      const hx = ((this.hoveredSite - 0.5) / totalSites) * w;
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

  setupEvents() {
    const getSiteFromX = (x) => {
      const totalSites = this.data?.numSites || 1;
      const frac = Math.max(0, Math.min(1, x / this.width));
      return Math.min(totalSites, Math.max(1, Math.floor(frac * totalSites) + 1));
    };

    const getRowFromY = (y) => {
      if (this.options.mode !== '2d') return 0;
      const profileH = 40;
      const matrixH = this.height - profileH - 12;
      const relY = y - (profileH + 2);
      const frac = Math.max(0, Math.min(1, relY / matrixH));
      return Math.min(this.totalRows - 1, Math.max(0, Math.floor(frac * this.totalRows)));
    };

    this.canvas.addEventListener('mousemove', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      if (this.isDragging) {
        // Dragging the viewport lens
        const totalSites = this.data?.numSites || 1;
        const dx = mx - this.dragStartX;
        const dCol = Math.round((dx / this.width) * totalSites);
        const newCol = Math.max(0, Math.min(totalSites - 1, this.dragStartCol + dCol));

        let newRow = this.dragStartRow;
        if (this.options.mode === '2d') {
          const profileH = 40;
          const matrixH = this.height - profileH - 12;
          const dy = my - this.dragStartY;
          const dRow = Math.round((dy / matrixH) * this.totalRows);
          newRow = Math.max(0, Math.min(this.totalRows - 1, this.dragStartRow + dRow));
        }

        if (this.options.onViewportChange) {
          this.options.onViewportChange(newCol, newRow);
        }
        return;
      }

      const site = getSiteFromX(mx);
      if (this.hoveredSite !== site) {
        this.hoveredSite = site;
        this.render();
      }
    });

    this.canvas.addEventListener('mouseleave', () => {
      if (!this.isDragging) {
        this.hoveredSite = -1;
        this.render();
      }
    });

    this.canvas.addEventListener('mousedown', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const totalSites = this.data?.numSites || 1;

      const startX = (this.visibleColRange.start / totalSites) * this.width;
      const endX = (Math.min(totalSites, this.visibleColRange.end) / totalSites) * this.width;
      const lensW = Math.max(6, endX - startX);

      // Check if clicked inside current lens
      if (mx >= startX && mx <= startX + lensW) {
        this.isDragging = true;
        this.dragStartX = mx;
        this.dragStartY = my;
        this.dragStartCol = this.visibleColRange.start;
        this.dragStartRow = this.visibleRowRange.start;
      } else {
        // Clicked outside lens: seek directly to clicked site
        const site = getSiteFromX(mx);
        const row = getRowFromY(my);
        const visibleCols = this.visibleColRange.end - this.visibleColRange.start;
        const newStartCol = Math.max(0, Math.min(totalSites - visibleCols, site - Math.floor(visibleCols / 2)));

        if (this.options.onViewportChange) {
          this.options.onViewportChange(newStartCol, row);
        }
        if (this.options.onSiteSelect && this.data) {
          this.options.onSiteSelect(site, this.data.siteResults[site - 1]);
        }
      }
    });

    window.addEventListener('mouseup', () => {
      if (this.isDragging) {
        this.isDragging = false;
        this.render();
      }
    });
  }
}
