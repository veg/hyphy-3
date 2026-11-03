/**
 * Sortable Industrial Table of 4-Way Wheatstone Bridges.
 */

export class BridgesTable {
  constructor(container, options = {}) {
    this.container = container;
    this.onBridgeSelect = options.onBridgeSelect || (() => {});
    this.model = null;
    this.tauThreshold = 0.0;
    this.selectedBridgeIndex = null;
    this.sortKey = 'shuntStress';
    this.sortDesc = true;

    this._initDOM();
  }

  _initDOM() {
    this.container.innerHTML = `
      <div class="ohmic-panel-header">
        <div class="ohmic-panel-title-group">
          <span class="ohmic-panel-title">WHEATSTONE 4-WAY BRIDGE REGISTRY (SERIES-PARALLEL SHUNTS)</span>
          <span class="ohmic-panel-subtitle" id="ohmic-bridges-subtitle">Click any bridge row to highlight its quartet clades (A, B, C, D) and open the Wheatstone Circuit Inspector</span>
        </div>
        <div class="ohmic-panel-controls">
          <button class="ohmic-icon-btn" id="ohmic-bridges-export-csv" title="Export Bridges Table as CSV">EXPORT CSV</button>
        </div>
      </div>
      <div class="ohmic-table-wrap">
        <table class="ohmic-table" id="ohmic-bridges-table">
          <thead>
            <tr>
              <th data-sort="index">BRIDGE</th>
              <th data-sort="hubNode">HUB / SPLIT NODE</th>
              <th>QUARTET TERMINALS (A | B ∥ C)</th>
              <th data-sort="shuntStress">SHUNT STRESS (τ_Ohm)</th>
              <th data-sort="ohmicSupport">OHMIC SUPPORT (π_Ohm)</th>
              <th>INITIAL WEIGHTS (w₁⁰, w₂⁰)</th>
              <th data-sort="activeShuntWeight">ML FITTED WEIGHTS (ŵ₁, ŵ₂)</th>
              <th>BRIDGE LENGTHS (t̂₁, t̂₂)</th>
              <th>CIRCUIT</th>
            </tr>
          </thead>
          <tbody id="ohmic-bridges-tbody"></tbody>
        </table>
      </div>
    `;

    this.tbody = this.container.querySelector('#ohmic-bridges-tbody');

    this.container.querySelectorAll('th[data-sort]').forEach((th) => {
      th.addEventListener('click', () => {
        const key = th.getAttribute('data-sort');
        if (this.sortKey === key) {
          this.sortDesc = !this.sortDesc;
        } else {
          this.sortKey = key;
          this.sortDesc = true;
        }
        this.renderRows();
      });
    });

    this.container.querySelector('#ohmic-bridges-export-csv')?.addEventListener('click', () => this.exportCSV());
  }

  setModel(model, tauThreshold = 0.0, selectedBridgeIndex = null) {
    this.model = model;
    this.tauThreshold = tauThreshold;
    this.selectedBridgeIndex = selectedBridgeIndex;
    this.renderRows();
  }

  setTauThreshold(tau) {
    this.tauThreshold = tau;
    this.renderRows();
  }

  setSelectedBridge(bridgeIndex) {
    this.selectedBridgeIndex = bridgeIndex;
    this.renderRows();
  }

  renderRows() {
    if (!this.tbody || !this.model) return;

    const filtered = this.model.bridges
      .filter((b) => b.shuntStress >= this.tauThreshold)
      .slice()
      .sort((a, b) => {
        const va = a[this.sortKey];
        const vb = b[this.sortKey];
        if (typeof va === 'string') {
          return this.sortDesc ? vb.localeCompare(va) : va.localeCompare(vb);
        }
        return this.sortDesc ? vb - va : va - vb;
      });

    if (filtered.length === 0) {
      this.tbody.innerHTML = `
        <tr>
          <td colspan="9" class="ohmic-table-empty">No Wheatstone bridges match current Shunt Stress threshold (τ_Ohm ≥ ${this.tauThreshold.toFixed(2)}).</td>
        </tr>
      `;
      return;
    }

    this.tbody.innerHTML = filtered
      .map((b) => {
        const isSelected = this.selectedBridgeIndex === b.index;
        const cA = b.quartetClades?.[0] || [b.children[0]];
        const cB = b.quartetClades?.[1] || [b.children[1]];
        const cC = b.quartetClades?.[2] || [b.children[2]];
        const tauPct = Math.min(100, Math.max(2, b.shuntStress * 100));
        const mlShuntActive = b.fittedWeights[1] > 0.01;

        return `
          <tr class="${isSelected ? 'selected' : ''}" data-bridge-idx="${b.index}">
            <td class="ohmic-mono"><strong>B${b.index + 1}</strong></td>
            <td class="ohmic-mono">${b.hubNode} <span class="ohmic-muted">/ ${b.bridgeNode}</span></td>
            <td class="ohmic-mono">
              <span class="ohmic-pill ohmic-pill-a" title="Clade A: ${cA.join(', ')}">A: ${b.children[0]} (${cA.length})</span>
              <span class="ohmic-pill ohmic-pill-b" title="Clade B: ${cB.join(', ')}">B: ${b.children[1]} (${cB.length})</span>
              <span class="ohmic-pill ohmic-pill-c" title="Clade C: ${cC.join(', ')}">C: ${b.children[2]} (${cC.length})</span>
            </td>
            <td>
              <div class="ohmic-bar-cell">
                <span class="ohmic-mono">${b.shuntStress.toFixed(4)}</span>
                <div class="ohmic-mini-bar"><div class="ohmic-mini-bar-fill" style="width:${tauPct.toFixed(1)}%"></div></div>
              </div>
            </td>
            <td class="ohmic-mono">${b.ohmicSupport.toFixed(4)}</td>
            <td class="ohmic-mono">(${b.initialWeights[0].toFixed(3)}, <span class="ohmic-text-accent">${b.initialWeights[1].toFixed(3)}</span>)</td>
            <td class="ohmic-mono">
              (${b.fittedWeights[0].toFixed(3)}, <strong class="${mlShuntActive ? 'ohmic-text-accent' : 'ohmic-muted'}">${b.fittedWeights[1].toFixed(4)}</strong>)
            </td>
            <td class="ohmic-mono">(${b.fittedBranchLengths[0].toFixed(4)}, ${b.fittedBranchLengths[1].toFixed(4)})</td>
            <td>
              <button class="ohmic-inspect-btn" data-inspect-idx="${b.index}">INSPECT ◆</button>
            </td>
          </tr>
        `;
      })
      .join('');

    this.tbody.querySelectorAll('tr[data-bridge-idx]').forEach((tr) => {
      tr.addEventListener('click', () => {
        const idx = parseInt(tr.getAttribute('data-bridge-idx'), 10);
        const br = this.model.bridges.find((x) => x.index === idx);
        if (br) {
          this.selectedBridgeIndex = idx;
          this.renderRows();
          this.onBridgeSelect(br);
        }
      });
    });
  }

  exportCSV() {
    if (!this.model) return;
    const headers = [
      'bridge_index',
      'hub_node',
      'bridge_node',
      'child_A',
      'child_B',
      'child_C',
      'shunt_stress_tau',
      'ohmic_support_pi',
      'init_w1',
      'init_w2',
      'fitted_w1',
      'fitted_w2',
      'fitted_t1',
      'fitted_t2'
    ];
    const rows = this.model.bridges.map((b) => [
      b.index + 1,
      b.hubNode,
      b.bridgeNode,
      b.children[0],
      b.children[1],
      b.children[2],
      b.shuntStress.toFixed(6),
      b.ohmicSupport.toFixed(6),
      b.initialWeights[0].toFixed(6),
      b.initialWeights[1].toFixed(6),
      b.fittedWeights[0].toFixed(6),
      b.fittedWeights[1].toFixed(6),
      b.fittedBranchLengths[0].toFixed(6),
      b.fittedBranchLengths[1].toFixed(6)
    ]);
    const csv = [headers.join(','), ...rows.map((r) => r.join(','))].join('\n');
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'ohmic_wheatstone_bridges.csv';
    a.click();
    URL.revokeObjectURL(url);
  }
}
