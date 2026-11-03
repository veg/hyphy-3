/**
 * Interactive Pairwise Electrical Network Heatmap Canvas (N x N).
 * Supports Foster Edge Probability (P_ij = c_ij * Omega_ij), Effective Resistance (Omega_ij),
 * Electrical Conductance (c_ij), and Pairwise Metric Distance (D_ij).
 * Highlights Quartet Clade blocks (A x B vs A x C) when a Wheatstone bridge is selected.
 */

export class FosterMatrixCanvas {
  constructor(container) {
    this.container = container;
    this.model = null;
    this.metricMode = 'foster'; // 'foster' | 'omega' | 'conductance' | 'distance'
    this.selectedBridgeIndex = null;
    this.hoveredCell = null;

    this._initDOM();
  }

  _initDOM() {
    this.container.innerHTML = `
      <div class="ohmic-panel-header">
        <div class="ohmic-panel-title-group">
          <span class="ohmic-panel-title">ELECTRICAL NETWORK MATRIX (L⁺ SPECTRUM)</span>
          <span class="ohmic-panel-subtitle" id="ohmic-matrix-subtitle">Foster Edge Probability P_ij = c_ij Ω_ij</span>
        </div>
        <div class="ohmic-panel-controls">
          <div class="ohmic-seg-group">
            <button class="ohmic-seg-btn active" data-metric="foster" title="Foster Edge Probability P_ij = c_ij * Omega_ij">FOSTER P_ij</button>
            <button class="ohmic-seg-btn" data-metric="omega" title="Effective Resistance Omega_ij = (e_i - e_j)^T L^+ (e_i - e_j)">RESISTANCE Ω_ij</button>
            <button class="ohmic-seg-btn" data-metric="conductance" title="Pairwise Conductance c_ij">CONDUCTANCE c_ij</button>
            <button class="ohmic-seg-btn" data-metric="distance" title="Pairwise Evolutionary Distance D_ij">DISTANCE D_ij</button>
          </div>
          <button class="ohmic-icon-btn" id="ohmic-matrix-export-png" title="Export Heatmap as PNG">PNG</button>
        </div>
      </div>
      <div class="ohmic-matrix-colorbar" id="ohmic-matrix-colorbar">
        <span class="ohmic-cb-label" id="ohmic-cb-min">0.000</span>
        <div class="ohmic-cb-gradient" id="ohmic-cb-grad"></div>
        <span class="ohmic-cb-label" id="ohmic-cb-max">1.000</span>
        <span class="ohmic-cb-note" id="ohmic-cb-note">∑ P_ij (i&lt;j) = N − 1</span>
      </div>
      <div class="ohmic-canvas-wrap" id="ohmic-matrix-canvas-wrap">
        <canvas id="ohmic-matrix-canvas"></canvas>
        <div class="ohmic-tooltip" id="ohmic-matrix-tooltip" style="display:none;"></div>
      </div>
    `;

    this.canvasWrap = this.container.querySelector('#ohmic-matrix-canvas-wrap');
    this.canvas = this.container.querySelector('#ohmic-matrix-canvas');
    this.tooltip = this.container.querySelector('#ohmic-matrix-tooltip');
    this.subtitleEl = this.container.querySelector('#ohmic-matrix-subtitle');
    this.cbMinEl = this.container.querySelector('#ohmic-cb-min');
    this.cbMaxEl = this.container.querySelector('#ohmic-cb-max');
    this.cbNoteEl = this.container.querySelector('#ohmic-cb-note');

    this.container.querySelectorAll('[data-metric]').forEach((btn) => {
      btn.addEventListener('click', () => {
        this.container.querySelectorAll('[data-metric]').forEach((b) => b.classList.remove('active'));
        btn.classList.add('active');
        this.metricMode = btn.getAttribute('data-metric');
        this.draw();
      });
    });

    this.container.querySelector('#ohmic-matrix-export-png')?.addEventListener('click', () => {
      if (!this.canvas) return;
      const a = document.createElement('a');
      a.href = this.canvas.toDataURL('image/png');
      a.download = `ohmic_matrix_${this.metricMode}.png`;
      a.click();
    });

    this.canvas.addEventListener('mousemove', (e) => this._onMouseMove(e));
    this.canvas.addEventListener('mouseleave', () => {
      this.hoveredCell = null;
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

  setModel(model, selectedBridgeIndex = null) {
    this.model = model;
    this.selectedBridgeIndex = selectedBridgeIndex;
    this.draw();
  }

  setSelectedBridge(bridgeIndex) {
    this.selectedBridgeIndex = bridgeIndex;
    this.draw();
  }

  _getActiveMatrix() {
    if (!this.model) return null;
    if (this.metricMode === 'foster') return this.model.pairwiseFoster;
    if (this.metricMode === 'omega') return this.model.pairwiseOmega;
    if (this.metricMode === 'conductance') return this.model.pairwiseConductance;
    if (this.metricMode === 'distance') return this.model.pairwiseDistances;
    return this.model.pairwiseFoster;
  }

  _colorForValue(t, mode, isDark) {
    const clamped = Math.max(0, Math.min(1, t));
    if (mode === 'foster' || mode === 'conductance') {
      // Braun Ivory (#F4F4F0) -> Warm Amber -> Braun Orange (#EB5E28) -> Deep Charcoal Red (#9A2B0E)
      if (isDark) {
        const r = Math.round(28 + clamped * 207);
        const g = Math.round(28 + Math.pow(clamped, 1.2) * 66);
        const b = Math.round(26 + (1 - clamped) * 14);
        return `rgb(${r},${g},${b})`;
      }
      const r = Math.round(245 - clamped * 10);
      const g = Math.round(244 - clamped * 150);
      const b = Math.round(239 - clamped * 199);
      return `rgb(${r},${g},${b})`;
    } else {
      // Resistance / Distance: Ivory -> Steel Blue (#1E40AF)
      if (isDark) {
        const r = Math.round(26 + clamped * 30);
        const g = Math.round(28 + clamped * 84);
        const b = Math.round(34 + clamped * 175);
        return `rgb(${r},${g},${b})`;
      }
      const r = Math.round(245 - clamped * 215);
      const g = Math.round(245 - clamped * 181);
      const b = Math.round(242 - clamped * 67);
      return `rgb(${r},${g},${b})`;
    }
  }

  draw() {
    if (!this.canvas || !this.model) return;
    const mat = this._getActiveMatrix();
    const rect = this.canvasWrap.getBoundingClientRect();
    const width = Math.max(300, Math.round(rect.width || 460));
    const numLeaves = Math.max(8, this.model.taxaCount || 10);
    const height = Math.max(390, Math.min(640, numLeaves * 21 + 64));

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

    ctx.fillStyle = bg;
    ctx.fillRect(0, 0, width, height);

    if (!mat || mat.length === 0) {
      ctx.fillStyle = muted;
      ctx.font = '500 12px "SF Mono", Menlo, monospace';
      ctx.textAlign = 'center';
      ctx.fillText('Pairwise matrix omitted in Woodbury k-NN mode (N > 128)', width / 2, height / 2);
      ctx.restore();
      return;
    }

    const N = mat.length;
    const taxa = this.model.taxa.length === N
      ? this.model.taxa
      : Array.from({ length: N }, (_, i) => `T${i + 1}`);

    // Compute min/max off-diagonal
    let minVal = Infinity;
    let maxVal = -Infinity;
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        if (i === j) continue;
        const v = Number(mat[i][j] ?? 0);
        if (v < minVal) minVal = v;
        if (v > maxVal) maxVal = v;
      }
    }
    if (!Number.isFinite(minVal)) minVal = 0;
    if (!Number.isFinite(maxVal) || maxVal <= minVal) maxVal = minVal + 1;

    // Update header & colorbar labels
    const titles = {
      foster: 'Foster Edge Probability P_ij = c_ij Ω_ij',
      omega: 'Effective Resistance Distance Ω_ij (Laplacian L⁺)',
      conductance: 'Pairwise Electrical Conductance c_ij = 1 / R_ij',
      distance: 'Pairwise Evolutionary Distance D_ij'
    };
    if (this.subtitleEl) this.subtitleEl.textContent = titles[this.metricMode] || '';
    if (this.cbMinEl) this.cbMinEl.textContent = minVal.toFixed(4);
    if (this.cbMaxEl) this.cbMaxEl.textContent = maxVal.toFixed(4);
    if (this.cbNoteEl) {
      this.cbNoteEl.textContent = this.metricMode === 'foster'
        ? `∑ P_ij (i<j) = ${this.model.fosterSum.toFixed(2)}`
        : `${N} × ${N} super-nodes`;
    }

    const labelPadLeft = N <= 30 ? 82 : 36;
    const labelPadTop = N <= 30 ? 68 : 28;
    const padRight = 18;
    const padBottom = 18;

    const availW = Math.max(120, width - labelPadLeft - padRight);
    const availH = Math.max(120, height - labelPadTop - padBottom);
    const gridSize = Math.min(availW, availH);
    const cellSize = gridSize / N;
    const offsetX = labelPadLeft + (availW - gridSize) * 0.5;
    const offsetY = labelPadTop + (availH - gridSize) * 0.5;

    this.gridGeom = { N, taxa, cellSize, offsetX, offsetY, gridSize, mat };

    // Determine Quartet role per taxon if a bridge is selected
    const roleMap = new Map();
    if (this.selectedBridgeIndex !== null) {
      const br = this.model.bridges.find((b) => b.index === this.selectedBridgeIndex);
      if (br && Array.isArray(br.quartetClades)) {
        ['A', 'B', 'C', 'D'].forEach((role, rIdx) => {
          (br.quartetClades[rIdx] || []).forEach((t) => roleMap.set(t, role));
        });
      }
    }

    // Draw cells
    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        const x = offsetX + j * cellSize;
        const y = offsetY + i * cellSize;
        if (i === j) {
          ctx.fillStyle = isDark ? '#252522' : '#E6E6E1';
          ctx.fillRect(x, y, cellSize, cellSize);
          continue;
        }
        const v = Number(mat[i][j] ?? 0);
        const norm = (v - minVal) / (maxVal - minVal);
        ctx.fillStyle = this._colorForValue(norm, this.metricMode, isDark);
        ctx.fillRect(x, y, cellSize, cellSize);

        // Quartet block highlight: (A, B) in Steel Blue, (A, C) in Braun Orange
        const rI = roleMap.get(taxa[i]);
        const rJ = roleMap.get(taxa[j]);
        if ((rI === 'A' && rJ === 'C') || (rI === 'C' && rJ === 'A')) {
          ctx.strokeStyle = '#EB5E28';
          ctx.lineWidth = Math.max(1.2, cellSize * 0.12);
          ctx.strokeRect(x + 0.5, y + 0.5, cellSize - 1, cellSize - 1);
        } else if ((rI === 'A' && rJ === 'B') || (rI === 'B' && rJ === 'A')) {
          ctx.strokeStyle = '#1E40AF';
          ctx.lineWidth = Math.max(1.0, cellSize * 0.1);
          ctx.strokeRect(x + 0.5, y + 0.5, cellSize - 1, cellSize - 1);
        }
      }
    }

    // Subtle cell grid if N <= 28
    if (N <= 28) {
      ctx.strokeStyle = isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.08)';
      ctx.lineWidth = 0.5;
      for (let k = 0; k <= N; k++) {
        ctx.beginPath();
        ctx.moveTo(offsetX + k * cellSize, offsetY);
        ctx.lineTo(offsetX + k * cellSize, offsetY + gridSize);
        ctx.stroke();

        ctx.beginPath();
        ctx.moveTo(offsetX, offsetY + k * cellSize);
        ctx.lineTo(offsetX + gridSize, offsetY + k * cellSize);
        ctx.stroke();
      }
    }

    // Hover crosshair
    if (this.hoveredCell) {
      const { i, j } = this.hoveredCell;
      ctx.strokeStyle = '#1A1A18';
      ctx.lineWidth = 1.8;
      ctx.strokeRect(offsetX + j * cellSize, offsetY + i * cellSize, cellSize, cellSize);
    }

    // Row & Column Taxon Labels if N <= 30
    if (N <= 30) {
      const fontSize = Math.min(10, Math.max(7.5, cellSize * 0.55));
      ctx.font = `500 ${fontSize.toFixed(1)}px "SF Mono", Menlo, monospace`;
      ctx.textBaseline = 'middle';

      for (let i = 0; i < N; i++) {
        const name = taxa[i];
        const role = roleMap.get(name);
        ctx.fillStyle = role === 'A' ? '#1E40AF' : role === 'C' ? '#EB5E28' : role === 'B' ? '#0D9488' : ink;

        // Left row label
        ctx.textAlign = 'right';
        const shortName = name.length > 10 ? name.slice(0, 9) + '…' : name;
        ctx.fillText(shortName, offsetX - 5, offsetY + (i + 0.5) * cellSize);

        // Top rotated column label
        ctx.save();
        ctx.translate(offsetX + (i + 0.5) * cellSize, offsetY - 5);
        ctx.rotate(-Math.PI / 3);
        ctx.textAlign = 'left';
        ctx.fillText(shortName, 0, 0);
        ctx.restore();
      }
    }

    ctx.restore();
  }

  _onMouseMove(e) {
    if (!this.gridGeom || !this.model) return;
    const rect = this.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    const { N, taxa, cellSize, offsetX, offsetY } = this.gridGeom;

    const j = Math.floor((mx - offsetX) / cellSize);
    const i = Math.floor((my - offsetY) / cellSize);

    if (i >= 0 && i < N && j >= 0 && j < N) {
      this.hoveredCell = { i, j };
      this.draw();

      const fVal = this.model.pairwiseFoster?.[i]?.[j] ?? 0;
      const oVal = this.model.pairwiseOmega?.[i]?.[j] ?? 0;
      const cVal = this.model.pairwiseConductance?.[i]?.[j] ?? 0;
      const dVal = this.model.pairwiseDistances?.[i]?.[j] ?? 0;

      this.tooltip.style.display = 'block';
      this.tooltip.style.left = `${Math.min(mx + 14, this.canvas.clientWidth - 230)}px`;
      this.tooltip.style.top = `${Math.max(10, my - 16)}px`;
      this.tooltip.innerHTML = `
        <div class="ohmic-tt-title">${taxa[i]} ↔ ${taxa[j]}</div>
        <div class="ohmic-tt-row"><span>Foster P_ij (c_ij Ω_ij):</span><strong>${Number(fVal).toFixed(4)}</strong></div>
        <div class="ohmic-tt-row"><span>Effective Resistance Ω_ij:</span><strong>${Number(oVal).toFixed(4)}</strong></div>
        <div class="ohmic-tt-row"><span>Conductance c_ij:</span><strong>${Number(cVal).toFixed(3)}</strong></div>
        <div class="ohmic-tt-row"><span>Evolutionary Dist D_ij:</span><strong>${Number(dVal).toFixed(4)}</strong></div>
      `;
    } else {
      if (this.hoveredCell) {
        this.hoveredCell = null;
        this.draw();
      }
      this.tooltip.style.display = 'none';
    }
  }
}
