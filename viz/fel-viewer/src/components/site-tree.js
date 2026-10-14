// Interactive Vector Site Tree Viewer (Phylogeny Annotated for Focal Codon Site)
// Maps ancestral root state, branch mutations (synonymous & nonsynonymous), and extant leaf states

import { translateCodon, SELECTION_COLORS } from '../data/genetic-code.js';
import { PALETTE, withAlpha } from '../data/palette.js';
import { formatEBF } from '../core/format.js';
import { downloadSVG } from '../core/export.js';
import { parseNewick } from '../core/newick.js';

export class SiteTree {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      mode: 'cladogram', // 'cladogram' | 'phylogram'
      labelMode: 'both', // 'both' | 'aa' | 'codon'
      rowHeight: 22,
      onBranchClick: null
    }, options);

    this.mode = this.options.mode;
    this.labelMode = this.options.labelMode;
    this.rowHeight = this.options.rowHeight;

    this.treeNewick = '';
    this.treeRoot = null;
    this.site = 1;
    this.siteResult = null;
    this.sequences = {};

    this.initDOM();
  }

  initDOM() {
    this.wrapper = document.createElement('div');
    this.wrapper.className = 'hyphy-site-tree-card';
    this.wrapper.innerHTML = `
      <div class="site-tree-header">
        <div class="site-tree-title-group">
          <span class="site-tree-title">Site Phylogeny</span>
          <span class="site-tree-subtitle" id="site-tree-sub">Codon #—</span>
        </div>
        <div class="site-tree-controls">
          <div class="segmented-control">
            <button type="button" class="seg-btn active" data-mode="cladogram">Cladogram</button>
            <button type="button" class="seg-btn" data-mode="phylogram">Phylogram</button>
          </div>
          <div class="segmented-control">
            <button type="button" class="seg-btn active" data-label="both">Codon + AA</button>
            <button type="button" class="seg-btn" data-label="aa">AA Only</button>
          </div>
          <button type="button" class="plot-tool-btn" id="btn-export-tree-svg" title="Export Vector SVG">SVG</button>
        </div>
      </div>

      <div class="site-tree-legend">
        <div class="tree-legend-item">
          <span class="tree-legend-line" style="background-color: ${PALETTE.treeLine};"></span>
          <span>Ancestral State Maintained</span>
        </div>
        <div class="tree-legend-item">
          <span class="tree-legend-line" style="background-color: ${PALETTE.accent}; height: 3px;"></span>
          <span style="color: ${PALETTE.accent}; font-weight: 700;">▲ Non-synonymous Mutation</span>
        </div>
        <div class="tree-legend-item">
          <span class="tree-legend-line" style="background-color: ${PALETTE.synonymous}; height: 3px;"></span>
          <span style="color: ${PALETTE.synonymous}; font-weight: 700;">Synonymous Mutation</span>
        </div>
        <div class="tree-legend-item">
          <span class="tree-legend-line" style="background-color: ${PALETTE.accent}; height: 3.5px; border-top: 1px solid ${PALETTE.accentPale};"></span>
          <span style="color: ${PALETTE.accent}; font-weight: 700;">⚡ Episodic Selection (EBF ≥ 100)</span>
        </div>
      </div>

      <div class="site-tree-scroll" id="tree-scroll">
        <div class="site-tree-svg-container" id="svg-container">
          <!-- SVG injected here -->
        </div>
      </div>
    `;

    this.container.appendChild(this.wrapper);

    this.svgContainer = this.wrapper.querySelector('#svg-container');
    this.subLabel = this.wrapper.querySelector('#site-tree-sub');

    // Controls
    const modeBtns = this.wrapper.querySelectorAll('[data-mode]');
    modeBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        modeBtns.forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        this.mode = btn.dataset.mode;
        this.render();
      });
    });

    const labelBtns = this.wrapper.querySelectorAll('[data-label]');
    labelBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        labelBtns.forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        this.labelMode = btn.dataset.label;
        this.render();
      });
    });

    this.wrapper.querySelector('#btn-export-tree-svg').addEventListener('click', () => this.exportSVG());
  }

  setTree(newickStr) {
    this.treeNewick = newickStr;
    if (newickStr) {
      // Substitution events and EBFs are matched against the labels as written,
      // {Group} tags included, and zero lengths draw as 0.05, as they always have.
      this.treeRoot = parseNewick(newickStr, { zeroLength: 0.05, tags: 'keep' });
    }
  }

  setSequences(seqs) {
    this.sequences = seqs || {};
  }

  setData(data) {
    this.setTree(data.treeNewick);
    this.setSequences(data.sequences);
  }

  setSite(siteNumber, siteResult) {
    this.site = siteNumber;
    this.siteResult = siteResult;
    this.render();
  }

  // Clone tree topology so layout coordinates and states don't mutate base tree
  cloneTree(node) {
    if (!node) return null;
    const n = {
      name: node.name,
      branch_length: node.branch_length || 0.05,
      children: (node.children || []).map(c => this.cloneTree(c))
    };
    return n;
  }

  assignStates(node, eventMap, rootCodon, rootAA) {
    const traverse = (n, parentCodon, parentAA) => {
      n.parentCodon = parentCodon;
      n.parentAA = parentAA;
      n.codon = parentCodon;
      n.aa = parentAA;
      n.hasMutation = false;
      n.event = null;

      if (n.name && eventMap.has(n.name)) {
        const ev = eventMap.get(n.name);
        n.codon = ev.child_codon;
        n.aa = ev.child_aa;
        n.hasMutation = true;
        n.event = ev;
      }

      // Check leaf alignment sequence if available
      if ((!n.children || n.children.length === 0) && n.name && this.sequences[n.name]) {
        const seqCodon = this.sequences[n.name].substr((this.site - 1) * 3, 3).toUpperCase();
        if (seqCodon && seqCodon.length === 3 && seqCodon !== '---') {
          n.codon = seqCodon;
          n.aa = translateCodon(seqCodon);
        }
      }

      for (const c of n.children || []) {
        traverse(c, n.codon, n.aa);
      }
    };

    traverse(node, rootCodon, rootAA);
  }

  layoutTree(root, mode, availableWidth) {
    // 1. Calculate Y positions (cladogram leaves)
    const yCursor = { val: 24 };
    const rowH = this.rowHeight;

    const computeY = (n) => {
      if (!n.children || n.children.length === 0) {
        n.y = yCursor.val;
        yCursor.val += rowH;
        n.leafCount = 1;
        return;
      }
      let sumY = 0;
      let leafSum = 0;
      for (const c of n.children) {
        computeY(c);
        sumY += c.y;
        leafSum += c.leafCount;
      }
      n.y = sumY / n.children.length;
      n.leafCount = leafSum;
    };
    computeY(root);

    const totalHeight = yCursor.val + 24;

    // 2. Calculate X positions
    let maxDist = 0;
    const computeDist = (n, currDist = 0) => {
      const d = currDist + (mode === 'phylogram' ? (n.branch_length || 0.05) : 1);
      n.dist = d;
      if (d > maxDist) maxDist = d;
      for (const c of n.children || []) {
        computeDist(c, d);
      }
    };
    root.dist = 0;
    for (const c of root.children || []) computeDist(c, 0);

    const leftPad = 32;
    const labelArea = 180; // space for leaf text on right
    const treeW = Math.max(160, availableWidth - leftPad - labelArea);

    const computeX = (n) => {
      n.x = leftPad + (maxDist > 0 ? (n.dist / maxDist) * treeW : 0);
      for (const c of n.children || []) computeX(c);
    };
    root.x = leftPad;
    for (const c of root.children || []) computeX(c);

    return { totalWidth: availableWidth, totalHeight };
  }

  render() {
    if (!this.treeRoot || !this.siteResult) {
      this.svgContainer.innerHTML = '<div class="empty-events-notice">No tree available for this dataset.</div>';
      return;
    }

    const r = this.siteResult;
    const rootCodon = r.rootCodon || '---';
    const rootAA = r.rootAA || translateCodon(rootCodon);

    this.subLabel.innerHTML = `Codon #${this.site} • Ancestral Root: <strong>${rootCodon}</strong> (${rootAA})`;

    const events = r.events || [];
    const eventMap = new Map();
    for (const ev of events) {
      if (ev.branch) eventMap.set(ev.branch, ev);
    }

    // Clone tree
    const tree = this.cloneTree(this.treeRoot);
    this.assignStates(tree, eventMap, rootCodon, rootAA);

    // Compute layout
    const containerW = Math.max(500, this.svgContainer.getBoundingClientRect().width || 600);
    const { totalWidth, totalHeight } = this.layoutTree(tree, this.mode, containerW);

    // Render SVG
    let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${totalWidth}" height="${totalHeight}" viewBox="0 0 ${totalWidth} ${totalHeight}">\n`;
    svg += `  <style>\n`;
    svg += `    .tree-line { fill: none; stroke: ${PALETTE.treeLine}; stroke-width: 1.5; stroke-linecap: round; stroke-linejoin: round; }\n`;
    svg += `    .tree-mut-nonsyn { stroke: ${PALETTE.accent} !important; stroke-width: 2.5 !important; }\n`;
    svg += `    .tree-mut-syn { stroke: ${PALETTE.synonymous} !important; stroke-width: 2.5 !important; }\n`;
    svg += `    .tree-branch-ebf { stroke: ${PALETTE.accent} !important; stroke-width: 3.5px !important; }\n`;
    svg += `    .tree-leaf-text { font-family: -apple-system, sans-serif; font-size: 11px; fill: ${PALETTE.ink}; font-weight: 500; dominant-baseline: middle; }\n`;
    svg += `    .tree-leaf-mono { font-family: "SF Mono", Monaco, monospace; font-size: 10px; font-weight: 700; }\n`;
    svg += `    .tree-badge-bg { rx: 3px; ry: 3px; }\n`;
    svg += `  </style>\n`;

    // 1. Draw Root Node
    svg += `  <circle cx="${tree.x}" cy="${tree.y}" r="4" fill="${PALETTE.ink}"/>\n`;
    svg += `  <text x="${tree.x - 8}" y="${tree.y}" text-anchor="end" class="tree-leaf-mono" fill="${PALETTE.stone}" dominant-baseline="middle">Root: ${rootCodon} (${rootAA})</text>\n`;

    // 2. Draw Branches & Nodes (Recursive)
    const drawBranches = (n) => {
      if (!n.children || n.children.length === 0) {
        // Leaf Node
        const isDiff = (n.aa !== rootAA);
        const codonLabel = n.codon || '---';
        const aaLabel = n.aa || '-';

        let badgeFill = PALETTE.inset;
        let badgeText = PALETTE.slate;
        if (isDiff) {
          badgeFill = withAlpha(PALETTE.accent, 0.15);
          badgeText = PALETTE.accent;
        } else if (n.hasMutation) {
          badgeFill = withAlpha(PALETTE.synonymous, 0.12);
          badgeText = PALETTE.synonymous;
        }

        let displayTag = `${codonLabel} [${aaLabel}]`;
        let tagW = 64;
        if (this.labelMode === 'aa') {
          displayTag = `[${aaLabel}]`;
          tagW = 32;
        }

        // Leaf dot
        svg += `  <circle cx="${n.x}" cy="${n.y}" r="3" fill="${isDiff ? PALETTE.accent : PALETTE.stone}"/>\n`;

        // Taxon name
        svg += `  <text x="${n.x + 8}" y="${n.y}" class="tree-leaf-text">${n.name || 'Leaf'}</text>\n`;

        // Codon / AA Pill Badge
        const tagX = n.x + 130;
        svg += `  <rect x="${tagX}" y="${n.y - 8}" width="${tagW}" height="16" fill="${badgeFill}" class="tree-badge-bg"/>\n`;
        svg += `  <text x="${tagX + tagW / 2}" y="${n.y}" class="tree-leaf-mono" fill="${badgeText}" text-anchor="middle" dominant-baseline="middle">${displayTag}</text>\n`;
        return;
      }

      // Vertical connector between children
      const minY = Math.min(...n.children.map(c => c.y));
      const maxY = Math.max(...n.children.map(c => c.y));
      svg += `  <line x1="${n.x}" y1="${minY}" x2="${n.x}" y2="${maxY}" class="tree-line"/>\n`;

      // Horizontal branch to each child
      for (const c of n.children) {
        let branchClass = 'tree-line';
        let strokeColor = PALETTE.treeLine;
        let badgeHtml = '';

        const ebf = (r.branchEBF && c.name && r.branchEBF[c.name]) || null;
        const isEBF = (ebf != null && ebf >= 100);

        if (c.hasMutation && c.event) {
          const isSyn = c.event.is_synonymous;
          branchClass = isSyn ? 'tree-line tree-mut-syn' : 'tree-line tree-mut-nonsyn';
          if (isEBF) branchClass += ' tree-branch-ebf';

          // Mutation Badge on branch
          const midX = (n.x + c.x) / 2;
          const mutLabel = isSyn
            ? `${c.event.parent_aa}→${c.event.child_aa} (syn)`
            : `${c.event.parent_aa}→${c.event.child_aa}${c.event.nt_distance > 1 ? '⚡' : ''}`;
          const label = isEBF ? `${mutLabel} • ⚡ EBF ${formatEBF(ebf, '≥10k')}` : mutLabel;
          const bW = Math.max(48, label.length * 6.5 + 10);

          badgeHtml = `
            <g class="branch-mutation-badge">
              <rect x="${midX - bW / 2}" y="${c.y - 15}" width="${bW}" height="13" rx="2" fill="${isSyn ? PALETTE.synonymous : PALETTE.accent}"/>
              <text x="${midX}" y="${c.y - 8}" fill="${PALETTE.knockout}" font-family="'SF Mono', Monaco, monospace" font-size="8.5px" font-weight="700" text-anchor="middle" dominant-baseline="middle">${label}</text>
            </g>
          `;
        } else if (isEBF) {
          branchClass += ' tree-branch-ebf';
          const midX = (n.x + c.x) / 2;
          const label = `⚡ EBF ${formatEBF(ebf, '≥10k')}`;
          const bW = Math.max(48, label.length * 6.5 + 10);

          badgeHtml = `
            <g class="branch-ebf-badge">
              <rect x="${midX - bW / 2}" y="${c.y - 15}" width="${bW}" height="13" rx="2" fill="${PALETTE.accent}"/>
              <text x="${midX}" y="${c.y - 8}" fill="${PALETTE.knockout}" font-family="'SF Mono', Monaco, monospace" font-size="8.5px" font-weight="700" text-anchor="middle" dominant-baseline="middle">${label}</text>
            </g>
          `;
        }

        svg += `  <line x1="${n.x}" y1="${c.y}" x2="${c.x}" y2="${c.y}" class="${branchClass}"/>\n`;
        if (badgeHtml) svg += badgeHtml;

        drawBranches(c);
      }
    };

    drawBranches(tree);
    svg += `</svg>\n`;

    this.svgContainer.innerHTML = svg;
    this.currentSVG = svg;
  }

  exportSVG() {
    if (!this.currentSVG) return;
    downloadSVG(this.currentSVG, `hyphy_site_${this.site}_phylogeny.svg`);
  }
}
