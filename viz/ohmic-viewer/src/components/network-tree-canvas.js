/**
 * Interactive Phylogenetic Electrical Network & Wheatstone Shunt Canvas.
 * Renders the backbone or collapsed tree with interactive 4-way Wheatstone
 * diamond shunt overlays, cross-branch shunt arcs (A <-> C), and quartet clade leaf tags.
 */

export class NetworkTreeCanvas {
  constructor(container, options = {}) {
    this.container = container;
    this.onBridgeSelect = options.onBridgeSelect || (() => {});
    this.model = null;
    this.treeMode = 'backbone'; // 'backbone' | 'collapsed'
    this.layoutMode = 'phylo'; // 'phylo' | 'cladogram'
    this.weightDisplayMode = 'initial'; // 'initial' | 'fitted'
    this.tauThreshold = 0.0;
    this.selectedBridgeIndex = null;
    this.hoveredBridgeIndex = null;
    this.hitRegions = [];

    this._initDOM();
  }

  _initDOM() {
    this.container.innerHTML = `
      <div class="ohmic-panel-header">
        <div class="ohmic-panel-title-group">
          <span class="ohmic-panel-title">WHEATSTONE PHYLOGENETIC NETWORK TOPOLOGY</span>
          <span class="ohmic-panel-subtitle" id="ohmic-tree-subtitle">Backbone tree + 4-way series-parallel shunts</span>
        </div>
        <div class="ohmic-panel-controls">
          <div class="ohmic-seg-group" title="Switch between binary backbone tree with shunt overlays vs collapsed 3-child hub topology">
            <button class="ohmic-seg-btn active" data-tree-mode="backbone">BACKBONE + SHUNTS</button>
            <button class="ohmic-seg-btn" data-tree-mode="collapsed">COLLAPSED HUBS</button>
          </div>
          <div class="ohmic-seg-group" title="Branch length scaling">
            <button class="ohmic-seg-btn active" data-layout-mode="phylo">PHYLO</button>
            <button class="ohmic-seg-btn" data-layout-mode="cladogram">CLADO</button>
          </div>
          <button class="ohmic-icon-btn" id="ohmic-tree-export-svg" title="Export Network as SVG">SVG</button>
          <button class="ohmic-icon-btn" id="ohmic-tree-export-png" title="Export Network as PNG">PNG</button>
        </div>
      </div>
      <div class="ohmic-tree-legend">
        <span class="ohmic-legend-item"><span class="ohmic-swatch ohmic-swatch-blue"></span>Primary Split (w₁: A|B)</span>
        <span class="ohmic-legend-item"><span class="ohmic-swatch ohmic-swatch-orange"></span>Wheatstone Shunt (w₂: A|C)</span>
        <span class="ohmic-legend-item"><span class="ohmic-swatch ohmic-swatch-teal"></span>Clade B Sister</span>
        <span class="ohmic-legend-item"><span class="ohmic-diamond-icon">◆</span>Click Bridge Hub to Inspect Circuit &amp; Simplex</span>
      </div>
      <div class="ohmic-canvas-wrap" id="ohmic-tree-canvas-wrap">
        <canvas id="ohmic-tree-canvas"></canvas>
        <div class="ohmic-tooltip" id="ohmic-tree-tooltip" style="display:none;"></div>
      </div>
    `;

    this.canvasWrap = this.container.querySelector('#ohmic-tree-canvas-wrap');
    this.canvas = this.container.querySelector('#ohmic-tree-canvas');
    this.tooltip = this.container.querySelector('#ohmic-tree-tooltip');

    this.container.querySelectorAll('[data-tree-mode]').forEach((btn) => {
      btn.addEventListener('click', () => {
        this.container.querySelectorAll('[data-tree-mode]').forEach((b) => b.classList.remove('active'));
        btn.classList.add('active');
        this.treeMode = btn.getAttribute('data-tree-mode');
        this.draw();
      });
    });

    this.container.querySelectorAll('[data-layout-mode]').forEach((btn) => {
      btn.addEventListener('click', () => {
        this.container.querySelectorAll('[data-layout-mode]').forEach((b) => b.classList.remove('active'));
        btn.classList.add('active');
        this.layoutMode = btn.getAttribute('data-layout-mode');
        this.draw();
      });
    });

    this.container.querySelector('#ohmic-tree-export-png')?.addEventListener('click', () => this.exportPNG());
    this.container.querySelector('#ohmic-tree-export-svg')?.addEventListener('click', () => this.exportSVG());

    this.canvas.addEventListener('mousemove', (e) => this._onMouseMove(e));
    this.canvas.addEventListener('mouseleave', () => {
      this.hoveredBridgeIndex = null;
      this.tooltip.style.display = 'none';
      this.draw();
    });
    this.canvas.addEventListener('click', (e) => this._onClick(e));

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

  setModel(model, tauThreshold = 0.0, selectedBridgeIndex = null) {
    this.model = model;
    this.tauThreshold = tauThreshold;
    this.selectedBridgeIndex = selectedBridgeIndex;
    this.draw();
  }

  setTauThreshold(tau) {
    this.tauThreshold = tau;
    this.draw();
  }

  setSelectedBridge(bridgeIndex) {
    this.selectedBridgeIndex = bridgeIndex;
    this.draw();
  }

  _computeLayout(treeObj, width, height) {
    const root = treeObj?.root;
    if (!root) return { nodes: [], leaves: [], maxDist: 1 };

    const padLeft = 36;
    const padRight = 165;
    const padTop = 28;
    const padBottom = 28;
    const plotW = Math.max(120, width - padLeft - padRight);
    const plotH = Math.max(120, height - padTop - padBottom);

    const leaves = treeObj.leaves;
    const N = Math.max(1, leaves.length);

    // Assign Y to leaves evenly
    leaves.forEach((leaf, idx) => {
      leaf.y = padTop + (N === 1 ? plotH / 2 : (idx / (N - 1)) * plotH);
    });

    // Post-order Y for internal nodes
    function assignY(node) {
      if (node.isLeaf || node.children.length === 0) return node.y;
      let sum = 0;
      for (const c of node.children) {
        sum += assignY(c);
      }
      node.y = sum / node.children.length;
      return node.y;
    }
    assignY(root);

    // Pre-order X distance from root
    let maxDist = 0;
    function assignDist(node, accDist, depth) {
      const step = node === root ? 0 : (node.branchLength || 0.01);
      node.distFromRoot = accDist + step;
      node.cladoDepth = depth;
      if (node.distFromRoot > maxDist) maxDist = node.distFromRoot;
      for (const c of node.children) {
        assignDist(c, node.distFromRoot, depth + 1);
      }
    }
    assignDist(root, 0, 0);

    // Find max cladogram depth
    let maxDepth = 1;
    for (const n of treeObj.nodes) {
      if (n.cladoDepth > maxDepth) maxDepth = n.cladoDepth;
    }

    for (const n of treeObj.nodes) {
      if (this.layoutMode === 'cladogram') {
        const frac = n.isLeaf ? 1.0 : n.cladoDepth / Math.max(1, maxDepth);
        n.x = padLeft + frac * plotW;
      } else {
        const frac = maxDist > 1e-9 ? n.distFromRoot / maxDist : n.cladoDepth / Math.max(1, maxDepth);
        n.x = padLeft + frac * plotW;
      }
    }

    return { nodes: treeObj.nodes, leaves, maxDist, padLeft, plotW };
  }

  _getQuartetRoleMap() {
    const roleMap = new Map();
    if (this.selectedBridgeIndex === null || !this.model) return roleMap;
    const br = this.model.bridges.find((b) => b.index === this.selectedBridgeIndex);
    if (!br || !Array.isArray(br.quartetClades)) return roleMap;

    const roles = ['A', 'B', 'C', 'D'];
    br.quartetClades.forEach((clade, idx) => {
      const role = roles[idx] || 'D';
      if (Array.isArray(clade)) {
        clade.forEach((taxon) => roleMap.set(taxon, role));
      }
    });
    return roleMap;
  }

  draw() {
    if (!this.canvas || !this.model) return;
    const rect = this.canvasWrap.getBoundingClientRect();
    const width = Math.max(320, Math.round(rect.width || 560));
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
    const gridLine = isDark ? '#2C2C28' : '#ECECE7';
    const accentOrange = '#EB5E28';
    const steelBlue = '#1E40AF';
    const tealSister = '#0D9488';

    ctx.fillStyle = bg;
    ctx.fillRect(0, 0, width, height);

    const treeObj = this.treeMode === 'collapsed' ? this.model.collapsedTree : this.model.backboneTree;
    const layout = this._computeLayout(treeObj, width, height);
    this.lastLayout = { layout, treeObj, width, height };
    this.hitRegions = [];

    const roleMap = this._getQuartetRoleMap();

    // Draw subtle vertical distance grid
    ctx.strokeStyle = gridLine;
    ctx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const gx = layout.padLeft + (i / 4) * layout.plotW;
      ctx.beginPath();
      ctx.moveTo(gx, 16);
      ctx.lineTo(gx, height - 16);
      ctx.stroke();
    }

    // 1. Draw standard tree branches (rectangular phylogram)
    for (const node of layout.nodes) {
      if (!node.parent) continue;
      const p = node.parent;

      // Check if this branch is a bridge internal node
      const brByBridge = this.model.bridgeByBridgeNodeName.get(node.name);
      const isActiveBridgeBranch = brByBridge && brByBridge.shuntStress >= this.tauThreshold;

      ctx.beginPath();
      ctx.moveTo(p.x, p.y);
      ctx.lineTo(p.x, node.y);
      ctx.lineTo(node.x, node.y);

      if (isActiveBridgeBranch) {
        ctx.strokeStyle = steelBlue;
        ctx.lineWidth = 2.4;
      } else {
        ctx.strokeStyle = isDark ? '#666660' : '#4A4A46';
        ctx.lineWidth = 1.5;
      }
      ctx.stroke();
    }

    // 2. Draw Wheatstone Bridge Shunt Arcs & Hub Diamonds
    for (const br of this.model.bridges) {
      if (br.shuntStress < this.tauThreshold) continue;

      const hubNode = treeObj.nodeByName.get(br.hubNode);
      const bridgeNode = treeObj.nodeByName.get(br.bridgeNode) || hubNode;
      const childA = treeObj.nodeByName.get(br.children[0]);
      const childB = treeObj.nodeByName.get(br.children[1]);
      const childC = treeObj.nodeByName.get(br.children[2]);

      const anchor = bridgeNode || hubNode;
      if (!anchor) continue;

      const isSelected = this.selectedBridgeIndex === br.index;
      const isHovered = this.hoveredBridgeIndex === br.index;

      // Draw curved Wheatstone shunt bypass chord connecting Clade A and Clade C
      if (childA && childC) {
        const ax = (anchor.x + childA.x) * 0.5;
        const ay = childA.y;
        const cx = ( (hubNode ? hubNode.x : anchor.x) + childC.x ) * 0.5;
        const cy = childC.y;
        const midX = Math.max(anchor.x, (ax + cx) * 0.5) + 18;
        const midY = (ay + cy) * 0.5;

        ctx.save();
        ctx.beginPath();
        ctx.moveTo(ax, ay);
        ctx.quadraticCurveTo(midX, midY, cx, cy);
        ctx.strokeStyle = accentOrange;
        ctx.lineWidth = isSelected ? 3.0 : isHovered ? 2.5 : Math.max(1.5, 3.2 * br.shuntStress);
        ctx.setLineDash([4, 3]);
        ctx.globalAlpha = isSelected || isHovered ? 1.0 : 0.82;
        ctx.stroke();
        ctx.restore();

        // Small resistor/shunt terminals at (ax, ay) and (cx, cy)
        [ [ax, ay], [cx, cy] ].forEach(([tx, ty]) => {
          ctx.beginPath();
          ctx.arc(tx, ty, 3, 0, Math.PI * 2);
          ctx.fillStyle = accentOrange;
          ctx.fill();
        });
      }

      // Draw Wheatstone Diamond at the hub/bridge anchor
      const dx = anchor.x;
      const dy = anchor.y;
      const r = isSelected ? 8.5 : isHovered ? 7.5 : 6.5;

      if (isSelected || isHovered) {
        ctx.beginPath();
        ctx.arc(dx, dy, r + 5, 0, Math.PI * 2);
        ctx.fillStyle = isSelected ? 'rgba(235, 94, 40, 0.22)' : 'rgba(30, 64, 175, 0.16)';
        ctx.fill();
      }

      ctx.beginPath();
      ctx.moveTo(dx, dy - r);
      ctx.lineTo(dx + r, dy);
      ctx.lineTo(dx, dy + r);
      ctx.lineTo(dx - r, dy);
      ctx.closePath();

      ctx.fillStyle = isSelected ? accentOrange : (br.activeShuntWeight > 0.01 ? accentOrange : steelBlue);
      ctx.fill();
      ctx.lineWidth = 1.5;
      ctx.strokeStyle = bg;
      ctx.stroke();

      // Badge label next to diamond
      ctx.font = '600 10px "SF Mono", Menlo, monospace';
      ctx.fillStyle = isSelected ? accentOrange : ink;
      const badgeText = `B${br.index + 1} τ=${br.shuntStress.toFixed(2)}`;
      ctx.fillText(badgeText, dx + r + 4, dy - 3);

      this.hitRegions.push({
        bridgeIndex: br.index,
        x: dx,
        y: dy,
        radius: 14,
        bridge: br
      });
    }

    // 3. Draw Leaf Nodes & Taxon Labels with Quartet Role Badges
    ctx.font = '500 11px "SF Mono", Menlo, monospace';
    ctx.textBaseline = 'middle';

    for (const leaf of layout.leaves) {
      const role = roleMap.get(leaf.name);
      let dotColor = isDark ? '#A0A09A' : '#3A3A36';
      let textColor = ink;
      let roleBadgeColor = null;

      if (role === 'A') {
        dotColor = steelBlue;
        textColor = steelBlue;
        roleBadgeColor = steelBlue;
      } else if (role === 'B') {
        dotColor = tealSister;
        textColor = tealSister;
        roleBadgeColor = tealSister;
      } else if (role === 'C') {
        dotColor = accentOrange;
        textColor = accentOrange;
        roleBadgeColor = accentOrange;
      } else if (role === 'D') {
        dotColor = muted;
        textColor = muted;
      }

      ctx.beginPath();
      ctx.arc(leaf.x, leaf.y, role ? 4.0 : 2.8, 0, Math.PI * 2);
      ctx.fillStyle = dotColor;
      ctx.fill();

      let labelX = leaf.x + 8;
      if (role && role !== 'D') {
        // Draw mini quartet terminal pill [A], [B], or [C]
        ctx.fillStyle = roleBadgeColor;
        ctx.fillRect(labelX, leaf.y - 6, 13, 12);
        ctx.fillStyle = '#FFFFFF';
        ctx.font = '700 9px "SF Mono", Menlo, monospace';
        ctx.fillText(role, labelX + 3, leaf.y + 0.5);
        labelX += 17;
      }

      ctx.font = role && role !== 'D' ? '700 11px "SF Mono", Menlo, monospace' : '500 11px "SF Mono", Menlo, monospace';
      ctx.fillStyle = textColor;
      ctx.fillText(leaf.name, labelX, leaf.y);
    }

    ctx.restore();
  }

  _findHit(e) {
    const rect = this.canvas.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;
    for (const h of this.hitRegions) {
      const dist = Math.hypot(mx - h.x, my - h.y);
      if (dist <= h.radius) return { hit: h, mx, my };
    }
    return { hit: null, mx, my };
  }

  _onMouseMove(e) {
    const { hit, mx, my } = this._findHit(e);
    if (hit) {
      this.canvas.style.cursor = 'pointer';
      if (this.hoveredBridgeIndex !== hit.bridgeIndex) {
        this.hoveredBridgeIndex = hit.bridgeIndex;
        this.draw();
      }
      const b = hit.bridge;
      this.tooltip.style.display = 'block';
      this.tooltip.style.left = `${Math.min(mx + 14, this.canvas.clientWidth - 250)}px`;
      this.tooltip.style.top = `${Math.max(10, my - 18)}px`;
      this.tooltip.innerHTML = `
        <div class="ohmic-tt-title">WHEATSTONE BRIDGE #${b.index + 1} (${b.hubNode})</div>
        <div class="ohmic-tt-row"><span>Triplet (A, B, C):</span><strong>${b.children.join(' • ')}</strong></div>
        <div class="ohmic-tt-row"><span>Shunt Stress (τ_Ohm):</span><strong>${b.shuntStress.toFixed(4)}</strong></div>
        <div class="ohmic-tt-row"><span>Ohmic Support (π_Ohm):</span><strong>${b.ohmicSupport.toFixed(4)}</strong></div>
        <div class="ohmic-tt-row"><span>Initial Weights (w₁, w₂):</span><strong>(${b.initialWeights[0].toFixed(3)}, ${b.initialWeights[1].toFixed(3)})</strong></div>
        <div class="ohmic-tt-row"><span>Fitted ML Weights:</span><strong>(${b.fittedWeights[0].toFixed(3)}, ${b.fittedWeights[1].toFixed(3)})</strong></div>
        <div class="ohmic-tt-hint">Click to inspect 4-terminal Wheatstone circuit &amp; quartet simplex</div>
      `;
    } else {
      this.canvas.style.cursor = 'default';
      if (this.hoveredBridgeIndex !== null) {
        this.hoveredBridgeIndex = null;
        this.draw();
      }
      this.tooltip.style.display = 'none';
    }
  }

  _onClick(e) {
    const { hit } = this._findHit(e);
    if (hit) {
      this.selectedBridgeIndex = hit.bridgeIndex;
      this.draw();
      this.onBridgeSelect(hit.bridge);
    }
  }

  exportPNG() {
    if (!this.canvas) return;
    const url = this.canvas.toDataURL('image/png');
    const a = document.createElement('a');
    a.href = url;
    a.download = `ohmic_network_${this.treeMode}.png`;
    a.click();
  }

  exportSVG() {
    if (!this.lastLayout || !this.model) return;
    const { layout, treeObj, width, height } = this.lastLayout;
    const lines = [];
    lines.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`);
    lines.push(`<rect width="100%" height="100%" fill="#FFFFFF"/>`);

    // Branches
    for (const node of layout.nodes) {
      if (!node.parent) continue;
      const p = node.parent;
      lines.push(`<path d="M ${p.x.toFixed(1)} ${p.y.toFixed(1)} L ${p.x.toFixed(1)} ${node.y.toFixed(1)} L ${node.x.toFixed(1)} ${node.y.toFixed(1)}" fill="none" stroke="#4A4A46" stroke-width="1.5"/>`);
    }

    // Bridges
    for (const br of this.model.bridges) {
      if (br.shuntStress < this.tauThreshold) continue;
      const hubNode = treeObj.nodeByName.get(br.hubNode);
      const bridgeNode = treeObj.nodeByName.get(br.bridgeNode) || hubNode;
      const childA = treeObj.nodeByName.get(br.children[0]);
      const childC = treeObj.nodeByName.get(br.children[2]);
      const anchor = bridgeNode || hubNode;
      if (!anchor) continue;

      if (childA && childC) {
        const ax = (anchor.x + childA.x) * 0.5;
        const ay = childA.y;
        const cx = ((hubNode ? hubNode.x : anchor.x) + childC.x) * 0.5;
        const cy = childC.y;
        const midX = Math.max(anchor.x, (ax + cx) * 0.5) + 18;
        const midY = (ay + cy) * 0.5;
        lines.push(`<path d="M ${ax.toFixed(1)} ${ay.toFixed(1)} Q ${midX.toFixed(1)} ${midY.toFixed(1)} ${cx.toFixed(1)} ${cy.toFixed(1)}" fill="none" stroke="#EB5E28" stroke-width="2.2" stroke-dasharray="4,3"/>`);
      }

      const r = 6.5;
      lines.push(`<polygon points="${anchor.x},${anchor.y - r} ${anchor.x + r},${anchor.y} ${anchor.x},${anchor.y + r} ${anchor.x - r},${anchor.y}" fill="#EB5E28" stroke="#FFFFFF" stroke-width="1.2"/>`);
    }

    // Leaves
    for (const leaf of layout.leaves) {
      lines.push(`<circle cx="${leaf.x.toFixed(1)}" cy="${leaf.y.toFixed(1)}" r="3" fill="#1A1A18"/>`);
      lines.push(`<text x="${(leaf.x + 8).toFixed(1)}" y="${(leaf.y + 3.5).toFixed(1)}" font-family="SF Mono, monospace" font-size="11" fill="#1A1A18">${leaf.name}</text>`);
    }

    lines.push(`</svg>`);
    const blob = new Blob([lines.join('\n')], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `ohmic_network_${this.treeMode}.svg`;
    a.click();
    URL.revokeObjectURL(url);
  }
}
