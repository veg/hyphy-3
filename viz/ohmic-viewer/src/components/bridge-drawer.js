/**
 * Slide-over Wheatstone Circuit & Quartet Simplex Inspector Drawer.
 * Renders an interactive 4-terminal Wheatstone diamond circuit diagram SVG,
 * a barycentric quartet simplex (w_1, w_2, w_3) comparing Initial Ohmic vs ML-Fitted weights,
 * and full Quartet Clade taxon breakdowns (A, B, C, D).
 */

export class BridgeDrawer {
  constructor(container, options = {}) {
    this.container = container;
    this.onClose = options.onClose || (() => {});
    this.bridge = null;
  }

  open(bridge) {
    this.bridge = bridge;
    if (!bridge) {
      this.close();
      return;
    }
    this.container.classList.add('open');
    this._render();
  }

  close() {
    this.bridge = null;
    this.container.classList.remove('open');
    this.container.innerHTML = '';
    this.onClose();
  }

  _barycentricToXY(w1, w2, w3, cx, cy, radius) {
    // Top vertex: T1 (AB|CD, w1=1)
    // Bottom-right vertex: T2 (AC|BD, w2=1)
    // Bottom-left vertex: T3 (BC|AD, w3=1)
    const sum = Math.max(1e-9, w1 + w2 + w3);
    const p1 = w1 / sum;
    const p2 = w2 / sum;
    const p3 = w3 / sum;

    const v1 = { x: cx, y: cy - radius };
    const v2 = { x: cx + radius * Math.cos(Math.PI / 6), y: cy + radius * Math.sin(Math.PI / 6) };
    const v3 = { x: cx - radius * Math.cos(Math.PI / 6), y: cy + radius * Math.sin(Math.PI / 6) };

    return {
      x: p1 * v1.x + p2 * v2.x + p3 * v3.x,
      y: p1 * v1.y + p2 * v2.y + p3 * v3.y,
      v1,
      v2,
      v3
    };
  }

  _renderCircuitSVG(b) {
    const w1 = b.fittedWeights[0].toFixed(3);
    const w2 = b.fittedWeights[1].toFixed(3);
    const t1 = b.fittedBranchLengths[0].toFixed(4);
    const t2 = b.fittedBranchLengths[1].toFixed(4);

    return `
      <svg class="ohmic-circuit-svg" viewBox="0 0 340 190" xmlns="http://www.w3.org/2000/svg">
        <!-- Diamond frame -->
        <line x1="170" y1="24" x2="64" y2="95" stroke="#4A4A46" stroke-width="1.8" />
        <line x1="170" y1="24" x2="276" y2="95" stroke="#4A4A46" stroke-width="1.8" />
        <line x1="64" y1="95" x2="170" y2="166" stroke="#1E40AF" stroke-width="3.0" />
        <line x1="276" y1="95" x2="170" y2="166" stroke="#0D9488" stroke-width="2.0" />

        <!-- Central Wheatstone Shunt Bridge (A <-> C) -->
        <line x1="64" y1="95" x2="276" y2="95" stroke="#EB5E28" stroke-width="2.8" stroke-dasharray="5,3" />
        <rect x="136" y="83" width="68" height="24" rx="2" fill="#FFFFFF" stroke="#EB5E28" stroke-width="1.5"/>
        <text x="170" y="94" text-anchor="middle" font-family="SF Mono, monospace" font-size="8.5" font-weight="700" fill="#EB5E28">SHUNT w₂=${w2}</text>
        <text x="170" y="103" text-anchor="middle" font-family="SF Mono, monospace" font-size="7.5" fill="#767672">t₂=${t2}</text>

        <!-- Primary Split Label on A-B arm -->
        <rect x="74" y="124" width="68" height="22" rx="2" fill="#FFFFFF" stroke="#1E40AF" stroke-width="1.2"/>
        <text x="108" y="134" text-anchor="middle" font-family="SF Mono, monospace" font-size="8.5" font-weight="700" fill="#1E40AF">SPLIT w₁=${w1}</text>
        <text x="108" y="142" text-anchor="middle" font-family="SF Mono, monospace" font-size="7.5" fill="#767672">t₁=${t1}</text>

        <!-- 4 Terminals: D (Top Hub), A (Left), C (Right), B (Bottom) -->
        <circle cx="170" cy="24" r="11" fill="#1A1A18" />
        <text x="170" y="28" text-anchor="middle" font-family="SF Mono, monospace" font-size="10" font-weight="700" fill="#FFFFFF">D</text>
        <text x="190" y="27" font-family="SF Mono, monospace" font-size="9.5" fill="#1A1A18">${b.hubNode} (Outgroup)</text>

        <circle cx="64" cy="95" r="12" fill="#1E40AF" />
        <text x="64" y="99" text-anchor="middle" font-family="SF Mono, monospace" font-size="10" font-weight="700" fill="#FFFFFF">A</text>
        <text x="46" y="98" text-anchor="end" font-family="SF Mono, monospace" font-size="9.5" font-weight="700" fill="#1E40AF">${b.children[0]}</text>

        <circle cx="276" cy="95" r="12" fill="#EB5E28" />
        <text x="276" y="99" text-anchor="middle" font-family="SF Mono, monospace" font-size="10" font-weight="700" fill="#FFFFFF">C</text>
        <text x="294" y="98" text-anchor="start" font-family="SF Mono, monospace" font-size="9.5" font-weight="700" fill="#EB5E28">${b.children[2]}</text>

        <circle cx="170" cy="166" r="12" fill="#0D9488" />
        <text x="170" y="170" text-anchor="middle" font-family="SF Mono, monospace" font-size="10" font-weight="700" fill="#FFFFFF">B</text>
        <text x="190" y="170" font-family="SF Mono, monospace" font-size="9.5" font-weight="700" fill="#0D9488">${b.children[1]}</text>
      </svg>
    `;
  }

  _renderSimplexSVG(b) {
    const cx = 170;
    const cy = 112;
    const radius = 76;

    const initPt = this._barycentricToXY(
      b.initialWeights[0],
      b.initialWeights[1],
      b.initialWeights[2],
      cx,
      cy,
      radius
    );
    const fitPt = this._barycentricToXY(
      b.fittedWeights[0],
      b.fittedWeights[1],
      b.fittedWeights[2],
      cx,
      cy,
      radius
    );
    const { v1, v2, v3 } = initPt;

    return `
      <svg class="ohmic-simplex-svg" viewBox="0 0 340 190" xmlns="http://www.w3.org/2000/svg">
        <!-- Centroid medians -->
        <line x1="${v1.x}" y1="${v1.y}" x2="${(v2.x + v3.x) / 2}" y2="${(v2.y + v3.y) / 2}" stroke="#E2E2DC" stroke-width="1" stroke-dasharray="2,2" />
        <line x1="${v2.x}" y1="${v2.y}" x2="${(v1.x + v3.x) / 2}" y2="${(v1.y + v3.y) / 2}" stroke="#E2E2DC" stroke-width="1" stroke-dasharray="2,2" />
        <line x1="${v3.x}" y1="${v3.y}" x2="${(v1.x + v2.x) / 2}" y2="${(v1.y + v2.y) / 2}" stroke="#E2E2DC" stroke-width="1" stroke-dasharray="2,2" />

        <!-- Equilateral Triangle -->
        <polygon points="${v1.x},${v1.y} ${v2.x},${v2.y} ${v3.x},${v3.y}" fill="rgba(242,242,239,0.5)" stroke="#1A1A18" stroke-width="1.6" />

        <!-- Vertex labels -->
        <text x="${v1.x}" y="${v1.y - 10}" text-anchor="middle" font-family="SF Mono, monospace" font-size="9.5" font-weight="700" fill="#1E40AF">T₁: (A,B)|C [w₁=1]</text>
        <text x="${v2.x + 8}" y="${v2.y + 14}" text-anchor="middle" font-family="SF Mono, monospace" font-size="9.5" font-weight="700" fill="#EB5E28">T₂: (A,C)|B [w₂=1]</text>
        <text x="${v3.x - 8}" y="${v3.y + 14}" text-anchor="middle" font-family="SF Mono, monospace" font-size="9.5" font-weight="600" fill="#767672">T₃: (B,C)|A [w₃=1]</text>

        <!-- Trajectory line from Initial (Ohmic) to Fitted (ML) -->
        <line x1="${initPt.x}" y1="${initPt.y}" x2="${fitPt.x}" y2="${fitPt.y}" stroke="#1A1A18" stroke-width="1.5" stroke-dasharray="3,2" />

        <!-- Initial Ohmic Weight Point -->
        <circle cx="${initPt.x}" cy="${initPt.y}" r="5.5" fill="#FFFFFF" stroke="#1E40AF" stroke-width="2.2" />
        <text x="${initPt.x + 9}" y="${initPt.y - 5}" font-family="SF Mono, monospace" font-size="8.5" fill="#1E40AF">Ohmic w⁰</text>

        <!-- Fitted ML Weight Point -->
        <circle cx="${fitPt.x}" cy="${fitPt.y}" r="6" fill="#EB5E28" stroke="#FFFFFF" stroke-width="1.5" />
        <text x="${fitPt.x + 9}" y="${fitPt.y + 10}" font-family="SF Mono, monospace" font-size="8.5" font-weight="700" fill="#EB5E28">ML ŵ</text>
      </svg>
    `;
  }

  _render() {
    const b = this.bridge;
    if (!b) return;

    const clades = b.quartetClades || [[], [], [], []];
    const roleMeta = [
      { key: 'A', label: `Terminal A (${b.children[0]})`, color: '#1E40AF', taxa: clades[0] || [] },
      { key: 'B', label: `Terminal B (${b.children[1]})`, color: '#0D9488', taxa: clades[1] || [] },
      { key: 'C', label: `Shunt Partner C (${b.children[2]})`, color: '#EB5E28', taxa: clades[2] || [] },
      { key: 'D', label: `Hub Outgroup D (${b.hubNode})`, color: '#767672', taxa: clades[3] || [] }
    ];

    this.container.innerHTML = `
      <div class="ohmic-drawer-inner">
        <div class="ohmic-drawer-header">
          <div>
            <div class="ohmic-drawer-badge">WHEATSTONE 4-WAY BRIDGE #${b.index + 1}</div>
            <div class="ohmic-drawer-title">${b.hubNode} → (${b.children.join(' • ')})</div>
          </div>
          <button class="ohmic-drawer-close" id="ohmic-drawer-close-btn" title="Close Inspector">✕</button>
        </div>

        <div class="ohmic-drawer-body">
          <!-- Section 1: 4-Terminal Circuit Schematic -->
          <div class="ohmic-drawer-section">
            <div class="ohmic-drawer-sec-title">4-TERMINAL SERIES-PARALLEL WHEATSTONE CIRCUIT</div>
            ${this._renderCircuitSVG(b)}
            <div class="ohmic-formula-box">
              D_hub^Ohm = w₁ [P(t₁)(M_A ⊙ M_B)] ⊙ M_C + w₂ [P(t₂)(M_A ⊙ M_C)] ⊙ M_B
            </div>
          </div>

          <!-- Section 2: Quartet Barycentric Simplex -->
          <div class="ohmic-drawer-section">
            <div class="ohmic-drawer-sec-title">BARYCENTRIC QUARTET SIMPLEX (INITIAL OHMIC → ML FITTED)</div>
            ${this._renderSimplexSVG(b)}
          </div>

          <!-- Section 3: Electrical & ML Parameters -->
          <div class="ohmic-drawer-section">
            <div class="ohmic-drawer-sec-title">ELECTRICAL &amp; MAXIMUM LIKELIHOOD DIAGNOSTICS</div>
            <div class="ohmic-kv-grid">
              <div class="ohmic-kv-item">
                <span>Shunt Stress (τ_k^Ohm)</span>
                <strong>${b.shuntStress.toFixed(4)}</strong>
              </div>
              <div class="ohmic-kv-item">
                <span>Ohmic Support (π_k^Ohm)</span>
                <strong>${b.ohmicSupport.toFixed(4)}</strong>
              </div>
              <div class="ohmic-kv-item">
                <span>Initial Ohmic (w₁⁰, w₂⁰, w₃⁰)</span>
                <strong>(${b.initialWeights.map((x) => x.toFixed(3)).join(', ')})</strong>
              </div>
              <div class="ohmic-kv-item">
                <span>ML Fitted (ŵ₁, ŵ₂, ŵ₃)</span>
                <strong class="ohmic-text-accent">(${b.fittedWeights.map((x) => x.toFixed(4)).join(', ')})</strong>
              </div>
              <div class="ohmic-kv-item">
                <span>Primary Branch Length (t̂₁)</span>
                <strong>${b.fittedBranchLengths[0].toFixed(5)}</strong>
              </div>
              <div class="ohmic-kv-item">
                <span>Shunt Branch Length (t̂₂)</span>
                <strong>${b.fittedBranchLengths[1].toFixed(5)}</strong>
              </div>
            </div>
          </div>

          <!-- Section 4: Quartet Clade Taxon Partition -->
          <div class="ohmic-drawer-section">
            <div class="ohmic-drawer-sec-title">QUARTET CLADE TAXON PARTITION (A | B ∥ C | D)</div>
            <div class="ohmic-clade-list">
              ${roleMeta
                .map(
                  (rm) => `
                <div class="ohmic-clade-card">
                  <div class="ohmic-clade-card-hdr" style="border-left: 3px solid ${rm.color};">
                    <strong>${rm.label}</strong>
                    <span>${rm.taxa.length} taxa</span>
                  </div>
                  <div class="ohmic-clade-taxa">${rm.taxa.join(', ') || '—'}</div>
                </div>
              `
                )
                .join('')}
            </div>
          </div>
        </div>
      </div>
    `;

    this.container.querySelector('#ohmic-drawer-close-btn')?.addEventListener('click', () => this.close());
  }
}
