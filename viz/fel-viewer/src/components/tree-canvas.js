// High-performance phylogenetic tree canvas renderer with 1-to-1 leaf row registration

import { PALETTE, withAlpha } from '../data/palette.js';

export class TreeCanvas {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      rowHeight: 22,
      font: '11px -apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif',
      branchColor: PALETTE.ink,
      testedBranchColor: PALETTE.ink,
      backgroundBranchColor: PALETTE.branchGrey,
      leafTextColor: PALETTE.ink,
      highlightColor: PALETTE.accent,
      activeBranchColor: PALETTE.accent,
      cladogram: false,
      onBranchSelect: null,
      onLeafHover: null
    }, options);

    this.canvas = document.createElement('canvas');
    this.canvas.className = 'hyphy-tree-canvas';
    this.ctx = this.canvas.getContext('2d');
    this.container.appendChild(this.canvas);

    this.root = null;
    this.leaves = [];             // In top-to-bottom visual order
    this.leafOrder = [];          // Array of leaf names in display order
    this.nodesByName = {};
    this.maxDistance = 1.0;
    this.activeSiteEvents = [];    // Substitution events for selected site
    this.activeSiteEBF = {};       // Branch name -> Empirical Bayes Factor for selected site
    this.hoveredNode = null;
    this.selectedNode = null;
    this.scrollTop = 0;
    this.showAncestors = false;
    this.ancestralTaxa = [];

    this.setupEvents();
  }

  setNewick(newickStr, nodeMetadata = {}) {
    if (!newickStr) return;
    this.root = this.parseNewick(newickStr);
    this.applyMetadata(this.root, nodeMetadata);
    this.ladderize(this.root);
    this.leaves = [];
    this.nodesByName = {};
    this.collectLeaves(this.root);
    this.leafOrder = this.leaves.map(l => l.name);
    this.layout();
    this.render();
  }

  parseNewick(str) {
    let s = str.trim();
    if (s.endsWith(';')) s = s.slice(0, -1);

    const root = { children: [] };
    const stack = [root];
    let i = 0;

    while (i < s.length) {
      const ch = s[i];
      if (ch === '(') {
        const newNode = { children: [] };
        stack[stack.length - 1].children.push(newNode);
        stack.push(newNode);
        i++;
      } else if (ch === ',') {
        i++;
      } else if (ch === ')') {
        const finishedNode = stack.pop();
        i++;
        // Optional name and branch length on internal node: )Node5:0.2824
        let token = '';
        while (i < s.length && s[i] !== ',' && s[i] !== ')' && s[i] !== '(') {
          token += s[i];
          i++;
        }
        token = token.trim();
        if (token.length > 0) {
          let name = '';
          let bl = 0.05;
          if (token.includes(':')) {
            const parts = token.split(':');
            name = parts[0].trim().replace(/^['"]|['"]$/g, '');
            bl = parseFloat(parts[1]) || 0.05;
          } else {
            name = token.replace(/^['"]|['"]$/g, '');
          }
          if (name) finishedNode.name = name;
          finishedNode.branch_length = bl;
        }
      } else {
        // Leaf node: read name and branch length
        let token = '';
        while (i < s.length && s[i] !== ',' && s[i] !== ')' && s[i] !== '(') {
          token += s[i];
          i++;
        }
        token = token.trim();
        if (token.length > 0) {
          let name = '';
          let bl = 0.05;
          if (token.includes(':')) {
            const parts = token.split(':');
            name = parts[0].trim().replace(/^['"]|['"]$/g, '');
            bl = parseFloat(parts[1]) || 0.05;
          } else {
            name = token.replace(/^['"]|['"]$/g, '');
          }
          const leafNode = { name, branch_length: bl, children: [] };
          stack[stack.length - 1].children.push(leafNode);
        }
      }
    }

    return root.children[0] || root;
  }

  applyMetadata(node, meta) {
    if (!node) return;
    if (node.name && meta[node.name]) {
      const m = meta[node.name];
      node.tested = m.tested !== undefined ? m.tested : true;
      if (m.branch_lengths?.codon_mg94 !== undefined) {
        node.branch_length = m.branch_lengths.codon_mg94;
      }
    } else {
      node.tested = true;
    }
    for (const child of node.children) {
      this.applyMetadata(child, meta);
    }
  }

  ladderize(node) {
    if (!node.children || node.children.length === 0) {
      node.leafCount = 1;
      return;
    }
    let total = 0;
    for (const c of node.children) {
      this.ladderize(c);
      total += c.leafCount;
    }
    node.leafCount = total;
    node.children.sort((a, b) => a.leafCount - b.leafCount);
  }

  collectLeaves(node) {
    if (node.name) this.nodesByName[node.name] = node;
    if (!node.children || node.children.length === 0) {
      node.isLeaf = true;
      this.leaves.push(node);
      return;
    }
    node.isLeaf = false;
    for (const c of node.children) {
      c.parent = node;
      this.collectLeaves(c);
    }
  }

  layout() {
    if (!this.root || this.leaves.length === 0) return;

    // 1. Assign Y coordinates (1-to-1 registration with sequence rows)
    const rowH = this.options.rowHeight;
    this.leaves.forEach((leaf, idx) => {
      leaf.y = idx * rowH + rowH / 2;
    });

    const assignInternalY = (node) => {
      if (node.isLeaf) return;
      for (const c of node.children) assignInternalY(c);
      const firstY = node.children[0].y;
      const lastY = node.children[node.children.length - 1].y;
      node.y = (firstY + lastY) / 2;
    };
    assignInternalY(this.root);

    // 2. Assign X coordinates based on cumulative branch lengths
    const assignDistance = (node, currentDist) => {
      node.dist = currentDist + (node.branch_length || 0);
      if (node.children) {
        for (const c of node.children) {
          assignDistance(c, node.dist);
        }
      }
    };
    this.root.dist = 0;
    if (this.root.children) {
      for (const c of this.root.children) assignDistance(c, 0);
    }

    let maxDist = 0;
    const findMaxDist = (node) => {
      if (node.dist > maxDist) maxDist = node.dist;
      if (node.children) {
        for (const c of node.children) findMaxDist(c);
      }
    };
    findMaxDist(this.root);
    this.maxDistance = Math.max(maxDist, 1e-4);
  }

  setActiveSiteEvents(events) {
    this.activeSiteEvents = events || [];
    this.render();
  }

  setActiveSiteEBF(ebfMap) {
    this.activeSiteEBF = ebfMap || {};
    this.render();
  }

  setShowAncestors(visible, ancestralTaxa = []) {
    this.showAncestors = visible;
    this.ancestralTaxa = ancestralTaxa;
    this.render();
  }

  setScrollTop(st) {
    const totalCount = this.leaves.length + (this.showAncestors ? (this.ancestralTaxa?.length || 0) : 0);
    const totalH = totalCount * this.options.rowHeight;
    const maxSt = Math.max(0, totalH - (this.height || totalH));
    this.scrollTop = Math.max(0, Math.min(st, maxSt));
    this.render();
  }

  resize(width, height) {
    const dpr = window.devicePixelRatio || 1;
    this.width = width;
    this.height = height;
    this.canvas.width = width * dpr;
    this.canvas.height = height * dpr;
    this.canvas.style.width = width + 'px';
    this.canvas.style.height = height + 'px';
    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    this.ctx.scale(dpr, dpr);
    this.render();
  }

  render() {
    if (!this.root || !this.ctx || !this.width) return;
    const ctx = this.ctx;
    const w = this.width;
    const h = this.height;
    const padX = 16;
    const labelPad = 130; // Max width for leaf labels on the right
    const treeW = Math.max(w - padX - labelPad, 60);

    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = PALETTE.panel;
    ctx.fillRect(0, 0, w, h);
    ctx.save();
    ctx.translate(0, -this.scrollTop);

    const rowH = this.options.rowHeight;
    const totalH = this.leaves.length * rowH;

    // Map cumulative distance to canvas X
    const getX = (dist) => {
      return padX + (dist / this.maxDistance) * treeW;
    };

    // Index active site events by branch
    const eventsByBranch = {};
    for (const ev of this.activeSiteEvents) {
      eventsByBranch[ev.branch] = ev;
    }

    // Recursive drawing
    const drawBranch = (node) => {
      const px = node.parent ? getX(node.parent.dist) : padX;
      const nx = getX(node.dist);
      const py = node.parent ? node.parent.y : node.y;
      const ny = node.y;

      // Determine branch styling
      const isTested = node.tested !== false;
      const ev = eventsByBranch[node.name];
      const isHovered = this.hoveredNode === node;
      const ebf = (node.name && this.activeSiteEBF) ? (this.activeSiteEBF[node.name] || 0) : 0;
      const isEpisodic = ebf >= 100.0;
      const hasModerateEBF = ebf >= 5.0 && ebf < 100.0;

      // If branch is under episodic selection (EBF >= 100), draw an incandescent glowing aura first
      if (isEpisodic) {
        ctx.save();
        ctx.strokeStyle = withAlpha(PALETTE.accent, 0.28);
        ctx.lineWidth = 7.0;
        ctx.lineCap = 'round';
        ctx.beginPath();
        ctx.moveTo(px, ny);
        ctx.lineTo(nx, ny);
        ctx.stroke();
        ctx.restore();
      }

      ctx.lineWidth = isHovered ? 2.5 : (isEpisodic ? 2.25 : (isTested ? 1.25 : 1.0));
      if (isHovered) {
        ctx.strokeStyle = PALETTE.accent; // Iconic Braun Orange on hover focus
      } else if (isEpisodic) {
        ctx.strokeStyle = PALETTE.accent; // Braun Signal Orange for episodic branch
      } else if (hasModerateEBF) {
        ctx.strokeStyle = PALETTE.amber; // Amber for moderate EBF (5..100)
      } else if (ev) {
        ctx.strokeStyle = ev.is_synonymous ? PALETTE.purifying : PALETTE.accent; // Cobalt for syn, Braun Orange for nonsyn
      } else if (isTested) {
        ctx.strokeStyle = this.options.testedBranchColor;
      } else {
        ctx.strokeStyle = this.options.backgroundBranchColor;
      }

      // Horizontal edge from parent's X to node's X
      ctx.beginPath();
      ctx.moveTo(px, ny);
      ctx.lineTo(nx, ny);
      ctx.stroke();

      // If episodic branch (EBF >= 100), draw a distinctive episodic diamond badge
      if (isEpisodic) {
        ctx.save();
        const midX = (px + nx) / 2;
        ctx.fillStyle = PALETTE.accent;
        ctx.beginPath();
        ctx.moveTo(midX, ny - 4.5);
        ctx.lineTo(midX + 4.5, ny);
        ctx.lineTo(midX, ny + 4.5);
        ctx.lineTo(midX - 4.5, ny);
        ctx.closePath();
        ctx.fill();
        ctx.strokeStyle = PALETTE.panel;
        ctx.lineWidth = 1.0;
        ctx.stroke();

        // EBF label tag above diamond
        ctx.font = 'bold 8px "SF Mono", Monaco, monospace';
        ctx.fillStyle = PALETTE.accent;
        ctx.textAlign = 'center';
        ctx.textBaseline = 'bottom';
        const ebfLabel = ebf >= 1000 ? `EBF>${Math.round(ebf)}` : `EBF:${Math.round(ebf)}`;
        ctx.fillText(ebfLabel, midX, ny - 5);
        ctx.restore();
      }

      // If active substitution event on this branch, draw a precision rivet dot
      if (ev) {
        ctx.save();
        const midX = isEpisodic ? ((px + nx) / 2 + 8) : ((px + nx) / 2);
        ctx.fillStyle = ev.is_synonymous ? PALETTE.purifying : PALETTE.accent;
        ctx.beginPath();
        ctx.arc(midX, ny, 4.0, 0, Math.PI * 2);
        ctx.fill();
        ctx.strokeStyle = PALETTE.panel;
        ctx.lineWidth = 1.0;
        ctx.stroke();

        // High perplexity lightning badge
        if (ev.surprise_bits >= 10.0) {
          ctx.fillStyle = PALETTE.amberOchre;
          ctx.font = 'bold 10px sans-serif';
          ctx.fillText('⚡', midX - 4, ny - 6);
        }
        ctx.restore();
      }

      // If internal node, draw vertical bar connecting children
      if (node.children && node.children.length > 0) {
        const topY = node.children[0].y;
        const botY = node.children[node.children.length - 1].y;

        ctx.strokeStyle = isHovered ? PALETTE.accent : (isTested ? this.options.testedBranchColor : this.options.backgroundBranchColor);
        ctx.lineWidth = isHovered ? 2.5 : (isTested ? 1.75 : 1.25);
        ctx.beginPath();
        ctx.moveTo(nx, topY);
        ctx.lineTo(nx, botY);
        ctx.stroke();

        if (this.showAncestors && node.name) {
          ctx.font = isHovered ? 'bold 9px "SF Mono", Monaco, monospace' : '9px "SF Mono", Monaco, monospace';
          ctx.fillStyle = isHovered ? PALETTE.accent : (isEpisodic ? PALETTE.accent : PALETTE.branchGrey);
          ctx.textAlign = 'left';
          ctx.textBaseline = 'bottom';
          ctx.fillText(node.name, nx + 3, ny - 2);
        }

        for (const c of node.children) {
          drawBranch(c);
        }
      } else {
        // Leaf label
        ctx.font = (isHovered || isEpisodic) ? `bold ${this.options.font}` : this.options.font;
        ctx.fillStyle = (isHovered || isEpisodic) ? PALETTE.accent : this.options.leafTextColor;
        ctx.textBaseline = 'middle';
        ctx.textAlign = 'left';
        ctx.fillText(node.name || '', nx + 8, ny);
      }
    };

    drawBranch(this.root);

    // Draw Ancestral Node rows if enabled
    if (this.showAncestors && this.ancestralTaxa && this.ancestralTaxa.length > 0) {
      const startRow = this.leaves.length;
      for (let i = 0; i < this.ancestralTaxa.length; ++i) {
        const rowIdx = startRow + i;
        const y = rowIdx * rowH;
        const ancName = this.ancestralTaxa[i];
        const isAncHovered = this.hoveredNode?.name === ancName;

        // Subtle row background wash
        ctx.fillStyle = PALETTE.cellBg;
        ctx.fillRect(0, y, w, rowH);

        // Subtle hairline divider
        ctx.strokeStyle = PALETTE.gridLine;
        ctx.lineWidth = 0.5;
        ctx.beginPath();
        ctx.moveTo(0, y + rowH);
        ctx.lineTo(w, y + rowH);
        ctx.stroke();

        // Ancestral node row label
        ctx.font = isAncHovered ? 'bold 11px -apple-system, sans-serif' : 'italic 11px -apple-system, sans-serif';
        ctx.fillStyle = isAncHovered ? PALETTE.accent : PALETTE.stone;
        ctx.textBaseline = 'middle';
        ctx.textAlign = 'left';
        const isRoot = !this.nodesByName[ancName]?.parent;
        const badge = isRoot ? ' [Root]' : '';
        ctx.fillText(`⮑ ${ancName}${badge}`, padX + 8, y + rowH / 2);
      }
    }

    ctx.restore();
  }

  setupEvents() {
    this.canvas.addEventListener('mousemove', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const my = e.clientY - rect.top + this.scrollTop;
      const rowH = this.options.rowHeight;
      const leafIdx = Math.floor(my / rowH);

      if (leafIdx >= 0 && leafIdx < this.leaves.length) {
        const leaf = this.leaves[leafIdx];
        if (this.hoveredNode !== leaf) {
          this.hoveredNode = leaf;
          this.render();
          if (this.options.onLeafHover) {
            this.options.onLeafHover(leaf.name, leafIdx);
          }
        }
      } else if (this.showAncestors && leafIdx >= this.leaves.length && leafIdx < this.leaves.length + (this.ancestralTaxa?.length || 0)) {
        const ancName = this.ancestralTaxa[leafIdx - this.leaves.length];
        const ancNode = this.nodesByName[ancName];
        if (this.hoveredNode !== ancNode) {
          this.hoveredNode = ancNode;
          this.render();
          if (this.options.onLeafHover) {
            this.options.onLeafHover(ancName, leafIdx);
          }
        }
      } else if (this.hoveredNode) {
        this.hoveredNode = null;
        this.render();
        if (this.options.onLeafHover) {
          this.options.onLeafHover(null, -1);
        }
      }
    });

    this.canvas.addEventListener('mouseleave', () => {
      if (this.hoveredNode) {
        this.hoveredNode = null;
        this.render();
        if (this.options.onLeafHover) {
          this.options.onLeafHover(null, -1);
        }
      }
    });

    this.canvas.addEventListener('click', (e) => {
      if (this.hoveredNode && this.options.onBranchSelect) {
        this.options.onBranchSelect(this.hoveredNode.name);
      }
    });
  }

  highlightBranch(branchName) {
    if (!branchName) {
      this.hoveredNode = null;
    } else {
      this.hoveredNode = this.nodesByName[branchName] || null;
    }
    this.render();
  }
}

