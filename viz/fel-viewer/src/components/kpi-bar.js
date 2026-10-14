// KPI Summary Bar & Reactive Threshold Controller
// Dieter Rams / Braun Industrial aesthetic: clean cards, tactile presets, reactive slider

import { formatThreshold, formatOmegaInf } from '../core/format.js';

export class KpiBar {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      threshold: 0.1,
      onThresholdChange: null
    }, options);

    this.threshold = this.options.threshold;
    this.data = null;
    this.initDOM();
  }

  initDOM() {
    this.el = document.createElement('div');
    this.el.className = 'hyphy-kpi-bar';
    this.el.innerHTML = `
      <div class="kpi-cards-grid">
        <!-- 1. Alignment Card -->
        <div class="kpi-card" id="kpi-alignment">
          <div class="kpi-label">Alignment</div>
          <div class="kpi-value" id="kpi-val-sites">—</div>
          <div class="kpi-sub" id="kpi-sub-taxa">— Taxa</div>
        </div>

        <!-- 2. Global Mean dN/dS Card -->
        <div class="kpi-card" id="kpi-omega">
          <div class="kpi-label">Global Mean dN/dS (ω)</div>
          <div class="kpi-value" id="kpi-val-omega">—</div>
          <div class="kpi-sub" id="kpi-sub-omega">Gene-wide Pressure</div>
        </div>

        <!-- 3. Diversifying Sites (Adaptive) Card -->
        <div class="kpi-card kpi-card-div" id="kpi-div">
          <div class="kpi-label">Diversifying (Positive)</div>
          <div class="kpi-value kpi-text-orange" id="kpi-val-div">—</div>
          <div class="kpi-sub" id="kpi-sub-div">Sites (β > α)</div>
        </div>

        <!-- 4. Purifying Sites (Negative) Card -->
        <div class="kpi-card kpi-card-pur" id="kpi-pur">
          <div class="kpi-label">Purifying (Constraint)</div>
          <div class="kpi-value kpi-text-blue" id="kpi-val-pur">—</div>
          <div class="kpi-sub" id="kpi-sub-pur">Sites (α > β)</div>
        </div>

        <!-- 5. Interactive Threshold Controller -->
        <div class="kpi-card kpi-card-control" id="kpi-controller">
          <div class="kpi-control-header">
            <span class="kpi-label">Significance Cutoff</span>
            <span class="kpi-thresh-display" id="kpi-thresh-display">p ≤ 0.10</span>
          </div>
          <div class="kpi-slider-row">
            <input type="range" class="kpi-slider" id="kpi-thresh-slider" min="0.001" max="0.20" step="0.005" value="0.10" />
          </div>
          <div class="kpi-presets-row">
            <button type="button" class="preset-pill active" data-p="0.10">p ≤ 0.10</button>
            <button type="button" class="preset-pill" data-p="0.05">p ≤ 0.05</button>
            <button type="button" class="preset-pill" data-p="0.01">p ≤ 0.01</button>
          </div>
        </div>
      </div>
    `;

    this.container.appendChild(this.el);

    // Bind slider & preset listeners
    this.slider = this.el.querySelector('#kpi-thresh-slider');
    this.display = this.el.querySelector('#kpi-thresh-display');
    this.presetBtns = this.el.querySelectorAll('.preset-pill');

    this.slider.addEventListener('input', (e) => {
      const val = parseFloat(e.target.value);
      this.setThreshold(val, true);
    });

    this.presetBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        const val = parseFloat(btn.dataset.p);
        this.slider.value = val;
        this.setThreshold(val, true);
      });
    });
  }

  setThreshold(p, triggerCallback = false) {
    this.threshold = p;
    this.display.textContent = `p ≤ ${formatThreshold(p)}`;

    // Update active preset button highlight
    this.presetBtns.forEach(btn => {
      const btnP = parseFloat(btn.dataset.p);
      if (Math.abs(btnP - p) < 0.002) {
        btn.classList.add('active');
      } else {
        btn.classList.remove('active');
      }
    });

    this.updateStats();

    if (triggerCallback && typeof this.options.onThresholdChange === 'function') {
      this.options.onThresholdChange(p);
    }
  }

  setData(data) {
    this.data = data;
    this.renderInitial();
    this.updateStats();
  }

  renderInitial() {
    if (!this.data) return;
    const sitesVal = this.el.querySelector('#kpi-val-sites');
    const taxaSub = this.el.querySelector('#kpi-sub-taxa');
    const omegaVal = this.el.querySelector('#kpi-val-omega');
    const omegaSub = this.el.querySelector('#kpi-sub-omega');

    sitesVal.textContent = `${this.data.numSites.toLocaleString()} Codons`;
    taxaSub.textContent = `${this.data.numTaxa} Taxa`;

    const omega = this.data.globalOmega;
    if (this.data.isMEME) {
      const omegaLabel = this.el.querySelector('#kpi-omega .kpi-label');
      omegaLabel.textContent = 'Analysis Model';
      omegaVal.textContent = 'MEME';
      omegaSub.textContent = 'Episodic Diversifying';

      const divLabel = this.el.querySelector('#kpi-div .kpi-label');
      divLabel.textContent = 'Episodically Diversifying';

      const purLabel = this.el.querySelector('#kpi-pur .kpi-label');
      purLabel.textContent = 'Neutral / Pervasive';
    } else {
      omegaVal.textContent = formatOmegaInf(omega);
      if (omega < 0.2) {
        omegaSub.textContent = 'Strong Purifying Constraint';
      } else if (omega < 0.7) {
        omegaSub.textContent = 'Pervasive Purifying Selection';
      } else if (omega <= 1.2) {
        omegaSub.textContent = 'Nearly Neutral Evolution';
      } else {
        omegaSub.textContent = 'Diversifying Selection';
      }
    }
  }

  updateStats() {
    if (!this.data) return;
    const stats = this.data.getStats(this.threshold);

    const divVal = this.el.querySelector('#kpi-val-div');
    const divSub = this.el.querySelector('#kpi-sub-div');
    const purVal = this.el.querySelector('#kpi-val-pur');
    const purSub = this.el.querySelector('#kpi-sub-pur');

    const pctDiv = ((stats.diversifying / Math.max(1, stats.totalSites)) * 100).toFixed(1);
    const pctPur = ((stats.purifying / Math.max(1, stats.totalSites)) * 100).toFixed(1);

    if (this.data.isMEME) {
      divVal.textContent = `▲ ${stats.diversifying}`;
      divSub.textContent = `${pctDiv}% of sites (p ≤ ${formatThreshold(this.threshold)})`;

      purVal.textContent = `${stats.neutral}`;
      purVal.className = 'kpi-value';
      purSub.textContent = `Conserved or Neutral sites`;
    } else {
      divVal.textContent = `▲ ${stats.diversifying}`;
      divSub.textContent = `${pctDiv}% of sites (β > α)`;

      purVal.textContent = `▼ ${stats.purifying}`;
      purVal.className = 'kpi-value kpi-text-blue';
      purSub.textContent = `${pctPur}% of sites (α > β)`;
    }
  }
}
