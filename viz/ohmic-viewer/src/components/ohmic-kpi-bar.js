/**
 * Dieter Rams / Braun Industrial KPI Header Bar for HyPhy OHMIC Viewer.
 */

export class OhmicKpiBar {
  constructor(container, options = {}) {
    this.container = container;
    this.onTauThresholdChange = options.onTauThresholdChange || (() => {});
    this.tauThreshold = options.initialTau ?? 0.0;
    this.model = null;
  }

  render(model, tauThreshold = this.tauThreshold) {
    this.model = model;
    this.tauThreshold = tauThreshold;
    if (!model) {
      this.container.innerHTML = '';
      return;
    }

    const activeCount = model.bridges.filter((b) => b.shuntStress >= this.tauThreshold).length;
    const mlRetainedCount = model.bridges.filter((b) => b.activeShuntWeight > 0.01).length;
    const fosterOk = model.fosterResidual < 1e-6;
    const deltaPositive = model.deltaLogL >= 0;
    const aiccPositive = model.deltaAICc >= 0;

    const solverShort = model.solver.includes('Woodbury')
      ? 'WOODBURY k-NN'
      : 'EXACT DENSE L⁺';

    const hyphaeonBadge = model.usedHyphaeon
      ? `<span class="ohmic-badge ohmic-badge-accent" title="384-D HyphAeon phylogenetic transformer embeddings fused with TN93 conductances">HYPHAEON ${model.embeddingDim}D</span>`
      : `<span class="ohmic-badge" title="${model.solver}">${solverShort}</span>`;

    const omegaText = model.mg94Params && Number.isFinite(model.mg94Params.omega_ratio)
      ? `ω = ${Number(model.mg94Params.omega_ratio).toFixed(3)}`
      : model.datatype.toUpperCase();

    this.container.innerHTML = `
      <div class="ohmic-kpi-grid">
        <!-- Card 1: Alignment & Solver -->
        <div class="ohmic-kpi-card">
          <div class="ohmic-kpi-header">
            <span class="ohmic-kpi-label">ALIGNMENT &amp; SOLVER</span>
            ${hyphaeonBadge}
          </div>
          <div class="ohmic-kpi-value-row">
            <span class="ohmic-kpi-main">${model.taxaCount}</span>
            <span class="ohmic-kpi-unit">TAXA</span>
            <span class="ohmic-kpi-sep">/</span>
            <span class="ohmic-kpi-main">${model.sitesCount}</span>
            <span class="ohmic-kpi-unit">${model.datatype === 'codon' ? 'CODONS' : 'SITES'}</span>
          </div>
          <div class="ohmic-kpi-sub">
            <span>${model.superNodesCount} electrical super-nodes • ${model.patternsCount} patterns • ${omegaText}</span>
          </div>
        </div>

        <!-- Card 2: Foster's Theorem Invariant -->
        <div class="ohmic-kpi-card">
          <div class="ohmic-kpi-header">
            <span class="ohmic-kpi-label">FOSTER SUM INVARIANT (∑ cₑΩₑ)</span>
            <span class="ohmic-badge ${fosterOk ? 'ohmic-badge-ok' : 'ohmic-badge-warn'}">
              ${fosterOk ? '● EXACT N−1' : 'Δ ' + model.fosterResidual.toExponential(1)}
            </span>
          </div>
          <div class="ohmic-kpi-value-row">
            <span class="ohmic-kpi-main">${model.fosterSum.toFixed(5)}</span>
            <span class="ohmic-kpi-unit">EXPECTED ${model.expectedFosterSum}</span>
          </div>
          <div class="ohmic-kpi-sub">
            <span>Effective resistance trace tr(L⁺L) = ${model.superNodesCount} − 1</span>
          </div>
        </div>

        <!-- Card 3: Wheatstone Bridges -->
        <div class="ohmic-kpi-card">
          <div class="ohmic-kpi-header">
            <span class="ohmic-kpi-label">WHEATSTONE 4-WAY BRIDGES</span>
            <span class="ohmic-badge ohmic-badge-blue">O(K) SERIES-PARALLEL</span>
          </div>
          <div class="ohmic-kpi-value-row">
            <span class="ohmic-kpi-main ohmic-text-accent">${activeCount}</span>
            <span class="ohmic-kpi-unit">/ ${model.bridges.length} SHUNTS</span>
            <span class="ohmic-kpi-sep">•</span>
            <span class="ohmic-kpi-main">${mlRetainedCount}</span>
            <span class="ohmic-kpi-unit">ML-ACTIVE</span>
          </div>
          <div class="ohmic-kpi-sub">
            <span>${model.numCollapsedPolytomies} collapsed zero-conductance polytomies</span>
          </div>
        </div>

        <!-- Card 4: Likelihood & AICc Improvement -->
        <div class="ohmic-kpi-card">
          <div class="ohmic-kpi-header">
            <span class="ohmic-kpi-label">NETWORK Δ LN(L) &amp; Δ AICc</span>
            <span class="ohmic-badge ${aiccPositive ? 'ohmic-badge-accent' : ''}">
              ${aiccPositive ? '▲ NETWORK PREFERRED' : 'SINGLE TREE'}
            </span>
          </div>
          <div class="ohmic-kpi-value-row">
            <span class="ohmic-kpi-main ${deltaPositive ? 'ohmic-text-accent' : ''}">
              ${deltaPositive ? '+' : ''}${model.deltaLogL.toFixed(2)}
            </span>
            <span class="ohmic-kpi-unit">NATS</span>
            <span class="ohmic-kpi-sep">/</span>
            <span class="ohmic-kpi-main ${aiccPositive ? 'ohmic-text-blue' : ''}">
              ${aiccPositive ? '+' : ''}${model.deltaAICc.toFixed(2)}
            </span>
            <span class="ohmic-kpi-unit">ΔAICc</span>
          </div>
          <div class="ohmic-kpi-sub">
            <span>Tree: ${model.singleTreeLogL.toFixed(1)} → Ohmic: ${model.ohmicLogL.toFixed(1)}</span>
          </div>
        </div>

        <!-- Card 5: Interactive Shunt Stress Filter -->
        <div class="ohmic-kpi-card ohmic-kpi-control-card">
          <div class="ohmic-kpi-header">
            <span class="ohmic-kpi-label">MIN SHUNT STRESS (τ_Ohm)</span>
            <span class="ohmic-kpi-mono-val" id="ohmic-tau-val">≥ ${this.tauThreshold.toFixed(2)}</span>
          </div>
          <div class="ohmic-slider-wrap">
            <input
              type="range"
              id="ohmic-tau-slider"
              class="ohmic-slider"
              min="0"
              max="1"
              step="0.02"
              value="${this.tauThreshold}"
            />
          </div>
          <div class="ohmic-kpi-sub">
            <span>Showing ${activeCount} of ${model.bridges.length} bridges • ${model.networkFavoringSitesCount} sites ΔlnL &gt; 0.05</span>
          </div>
        </div>
      </div>
    `;

    const slider = this.container.querySelector('#ohmic-tau-slider');
    const valLabel = this.container.querySelector('#ohmic-tau-val');
    if (slider) {
      slider.addEventListener('input', (e) => {
        const v = parseFloat(e.target.value);
        this.tauThreshold = v;
        if (valLabel) valLabel.textContent = `≥ ${v.toFixed(2)}`;
        this.onTauThresholdChange(v);
      });
    }
  }
}
