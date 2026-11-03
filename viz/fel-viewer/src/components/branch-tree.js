// Branch tree: the whole-tree renderer for branch-level methods. It lays a Newick tree out
// as a rectangular phylogram with one row per leaf, takes branch lengths from v3 node
// metadata under a chosen key, and leaves every colour decision to a per-branch style
// callback. The canvas and the SVG export draw the same list of primitives.
//
// style(node, meta) may return any of:
//   stroke, width, dash      the horizontal branch (dash is a setLineDash array)
//   segments                 [{ fraction, stroke }] painted in order from the parent end,
//                            e.g. aBSREL omega classes weighted by their proportions
//   label                    text centred above the branch midpoint
//   labelColor               defaults to the stroke
// node is the parsed Newick node ({ name, tag, children, parent }), meta is
// nodeMeta[node.name] and colors the palette in use (the page's themed tokens on the
// canvas, PALETTE in the export). Unset fields fall back to a plain ink branch.
//
// The hovered branch and the selected branch both draw in the hover colour. Hover follows
// the pointer (or highlight(name) from a linked view); the selection is set by a click or
// select(name) and stays until the next one. The tree redraws whenever its container's
// width changes, so a tree built inside a hidden tab draws once the tab is shown.

import { PALETTE } from '../data/palette.js';
import { parseNewick, lengthFrom } from '../core/newick.js';
import { downloadPNG, downloadSVG, escapeSvgText } from '../core/export.js';
import { trackColors } from './site-track.js';

const BRANCH_TREE_FONT = '-apple-system, BlinkMacSystemFont, "Helvetica Neue", Arial, sans-serif';

const BRANCH_TREE_SIZES = {
  full: { rowHeight: 22, font: 11, labelFont: 9, width: 1.25, scale: 1, padX: 16, padY: 12, rootRadius: 3 },
  mini: { rowHeight: 12, font: 9, labelFont: 8, width: 1, scale: 0.7, padX: 8, padY: 6, rootRadius: 2 }
};

// Parents, ladderized order, leaf rows and depths. Depth sums the lengths from nodeMeta
// under lengthKey, or the Newick lengths when lengthKey is null; a tree without any length
// (or cladogram: true) is laid out by edge count instead. Rows are leaf indices, with
// internal nodes centred on their first and last child.
export function branchTreeLayout(root, options = {}) {
  const { nodeMeta = {}, lengthKey = null, ladderize = true, cladogram = false } = options;

  const count = (node, parent) => {
    node.parent = parent;
    node.leafCount = node.children.length ? 0 : 1;
    for (const child of node.children) node.leafCount += count(child, node);
    if (ladderize) node.children.sort((a, b) => a.leafCount - b.leafCount);
    return node.leafCount;
  };
  count(root, null);

  const nodes = [];
  const leaves = [];
  const length = node => (lengthKey ? lengthFrom(nodeMeta[node.name], lengthKey) : node.branch_length) ?? 0;
  const place = (node, depth, edges) => {
    nodes.push(node);
    node.depth = depth;
    node.edges = edges;
    for (const child of node.children) place(child, depth + length(child), edges + 1);
    if (node.children.length) {
      node.row = (node.children[0].row + node.children[node.children.length - 1].row) / 2;
    } else {
      node.row = leaves.length;
      leaves.push(node);
    }
  };
  place(root, 0, 0);

  const deepest = () => nodes.reduce((max, n) => Math.max(max, n.depth), 0);
  let maxDepth = deepest();
  const byEdges = cladogram || !(maxDepth > 0);
  if (byEdges) {
    for (const node of nodes) node.depth = node.edges;
    maxDepth = Math.max(1, deepest());
  }
  return { root, nodes, leaves, maxDepth, cladogram: byEdges };
}

// The corners of a diamond marker of half-diagonal r, clockwise from the top
export function diamondPoints({ x, y, r }) {
  return [[x, y - r], [x + r, y], [x, y + r], [x - r, y]];
}

export class BranchTree {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      lengthKey: null,          // nodeMeta[name].branch_lengths key; null uses the Newick lengths
      style: null,              // (node, meta) => branch style, see the top of this file
      hoverColor: 'accent',     // a colour, or a palette key resolved in the theme at paint time
      rootShape: 'dot',         // 'dot' | 'diamond' (a root the method placed itself)
      mini: false,              // compact rows and type, for drawers
      ladderize: true,
      cladogram: false,
      internalLabels: false,    // draw internal node names next to their connector
      leafLabel: null,          // (node, meta) => leaf text; defaults to the node name
      exportName: 'branch-tree',
      onSelect: null,           // (name, node) on click
      onHover: null             // (name, node), or (null, null) when the pointer leaves
    }, options);
    this.size = this.options.mini ? BRANCH_TREE_SIZES.mini : BRANCH_TREE_SIZES.full;

    this.wrapper = document.createElement('div');
    this.wrapper.className = this.options.mini ? 'hyphy-branch-tree mini' : 'hyphy-branch-tree';
    this.legend = document.createElement('div');
    this.legend.className = 'branch-tree-legend';
    this.legend.style.display = 'none';
    this.canvas = document.createElement('canvas');
    this.canvas.className = 'hyphy-tree-canvas';
    this.ctx = this.canvas.getContext('2d');
    this.wrapper.append(this.legend, this.canvas);
    this.container.appendChild(this.wrapper);

    this.layout = null;
    this.nodeMeta = {};
    this.nodesByName = {};
    this.edges = [];
    this.hovered = null;
    this.selected = null;
    this.width = 0;

    this.setupEvents();
    this.resizeObserver = new ResizeObserver(() => {
      const width = this.container.clientWidth;
      if (width !== this.width) this.resize(width);
    });
    this.resizeObserver.observe(this.container);
  }

  destroy() {
    this.resizeObserver.disconnect();
  }

  setTree(newick, nodeMeta = {}) {
    const root = parseNewick(newick, { zeroLength: 'keep' });
    this.nodeMeta = nodeMeta;
    this.layout = root ? branchTreeLayout(root, {
      nodeMeta,
      lengthKey: this.options.lengthKey,
      ladderize: this.options.ladderize,
      cladogram: this.options.cladogram
    }) : null;
    this.nodesByName = {};
    for (const node of this.layout?.nodes || []) {
      if (node.name) this.nodesByName[node.name] = node;
    }
    this.hovered = null;
    this.selected = this.nodesByName[this.selected?.name] || null;
    this.resize(this.width || this.container.clientWidth);
  }

  // Swatches above the tree: [{ label, stroke, dash }]
  setLegend(items = []) {
    this.legend.replaceChildren(...items.map(item => {
      const entry = document.createElement('span');
      entry.className = 'branch-tree-legend-item';
      const swatch = document.createElement('span');
      swatch.className = 'branch-tree-swatch';
      swatch.style.borderTop = `2px ${item.dash ? 'dashed' : 'solid'} ${item.stroke}`;
      entry.append(swatch, item.label);
      return entry;
    }));
    this.legend.style.display = items.length ? '' : 'none';
  }

  setStyle(style) {
    this.options.style = style;
    this.render();
  }

  // Transient emphasis from a linked view; the next pointer move replaces it
  highlight(name) {
    this.hovered = (name && this.nodesByName[name]) || null;
    this.render();
  }

  // Persistent emphasis: survives hover until the next select or click
  select(name) {
    this.selected = (name && this.nodesByName[name]) || null;
    this.render();
  }

  leafLabel(node) {
    return this.options.leafLabel ? this.options.leafLabel(node, this.nodeMeta[node.name]) : node.name;
  }

  // The canvas is as tall as the tree needs; the container scrolls if it must
  resize(width) {
    if (!width || width <= 0) return;
    const { rowHeight, padY } = this.size;
    const dpr = window.devicePixelRatio || 1;
    this.width = width;
    this.height = padY * 2 + (this.layout?.leaves.length || 0) * rowHeight;
    this.canvas.width = width * dpr;
    this.canvas.height = this.height * dpr;
    this.canvas.style.width = width + 'px';
    this.canvas.style.height = this.height + 'px';
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.render();
  }

  font(px, bold) {
    return `${bold ? 'bold ' : ''}${px}px ${BRANCH_TREE_FONT}`;
  }

  // Lines, texts and dots for the current width, hover and selection. Also refreshes the
  // hit-test edges: one per branch, reaching past the leaf label for leaves.
  scene(colors = PALETTE) {
    const { nodes, leaves, maxDepth, root } = this.layout;
    const { rowHeight, font, labelFont, width: baseWidth, scale, padX, padY, rootRadius } = this.size;
    const hover = colors[this.options.hoverColor] ?? this.options.hoverColor;

    this.ctx.font = this.font(font, true);
    const labelW = leaves.reduce((max, leaf) => Math.max(max, this.ctx.measureText(this.leafLabel(leaf)).width), 0) + 8;
    const treeW = Math.max(this.width - padX * 2 - labelW, 40);
    const x = depth => padX + (depth / maxDepth) * treeW;
    const y = row => padY + row * rowHeight + rowHeight / 2;

    const lines = [];
    const texts = [];
    const rootR = this.options.rootShape === 'diamond' ? rootRadius + 2 : rootRadius;
    const dots = [{ x: x(0), y: y(root.row), r: rootR, fill: colors.ink, shape: this.options.rootShape }];
    this.edges = [];

    for (const node of nodes) {
      const hovered = node === this.hovered || node === this.selected;
      const ny = y(node.row);
      const nx = x(node.depth);

      if (node.parent) {
        const st = (this.options.style && this.options.style(node, this.nodeMeta[node.name], colors)) || {};
        const stroke = hovered ? hover : (st.stroke || colors.ink);
        const width = (st.width ? st.width * scale : baseWidth) + (hovered ? 1.25 : 0);
        const px = x(node.parent.depth);
        if (st.segments?.length && !hovered) {
          let from = px;
          for (const seg of st.segments) {
            const to = from + seg.fraction * (nx - px);
            lines.push({ x1: from, y1: ny, x2: to, y2: ny, stroke: seg.stroke, width, dash: st.dash });
            from = to;
          }
        } else {
          lines.push({ x1: px, y1: ny, x2: nx, y2: ny, stroke, width, dash: st.dash });
        }
        if (st.label) {
          texts.push({ x: (px + nx) / 2, y: ny - width / 2 - 2, text: st.label, size: labelFont, bold: hovered,
                       fill: hovered ? hover : (st.labelColor || stroke), align: 'center', baseline: 'bottom' });
        }
        this.edges.push({ node, x0: px, x1: node.children.length ? nx : nx + labelW, y: ny });
      }

      if (node.children.length) {
        lines.push({ x1: nx, y1: y(node.children[0].row), x2: nx, y2: y(node.children[node.children.length - 1].row),
                     stroke: hovered ? hover : colors.ink, width: baseWidth + (hovered ? 1.25 : 0) });
        if (this.options.internalLabels && node.name && node.parent) {
          texts.push({ x: nx + 3, y: ny - 2, text: node.name, size: labelFont, bold: hovered,
                       fill: hovered ? hover : colors.stone, align: 'left', baseline: 'bottom' });
        }
      } else {
        texts.push({ x: nx + 6, y: ny, text: this.leafLabel(node), size: font, bold: hovered,
                     fill: hovered ? hover : colors.ink, align: 'left', baseline: 'middle' });
      }
    }
    return { lines, texts, dots };
  }

  render() {
    if (!this.ctx || !this.width) return;
    const ctx = this.ctx;
    ctx.clearRect(0, 0, this.width, this.height);
    const colors = trackColors(this.canvas);
    ctx.fillStyle = colors.panel;
    ctx.fillRect(0, 0, this.width, this.height);
    if (!this.layout) return;

    const { lines, texts, dots } = this.scene(colors);
    ctx.lineCap = 'butt';
    for (const line of lines) {
      ctx.strokeStyle = line.stroke;
      ctx.lineWidth = line.width;
      ctx.setLineDash(line.dash || []);
      ctx.beginPath();
      ctx.moveTo(line.x1, line.y1);
      ctx.lineTo(line.x2, line.y2);
      ctx.stroke();
    }
    ctx.setLineDash([]);
    for (const dot of dots) {
      ctx.fillStyle = dot.fill;
      ctx.beginPath();
      if (dot.shape === 'diamond') {
        diamondPoints(dot).forEach(([px, py], i) => (i ? ctx.lineTo(px, py) : ctx.moveTo(px, py)));
        ctx.closePath();
      } else {
        ctx.arc(dot.x, dot.y, dot.r, 0, Math.PI * 2);
      }
      ctx.fill();
    }
    for (const text of texts) {
      ctx.font = this.font(text.size, text.bold);
      ctx.fillStyle = text.fill;
      ctx.textAlign = text.align;
      ctx.textBaseline = text.baseline;
      ctx.fillText(text.text, text.x, text.y);
    }
  }

  svg() {
    const w = this.width;
    const h = this.height;
    let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}">\n`;
    svg += `  <rect width="100%" height="100%" fill="${PALETTE.panel}"/>\n`;
    if (this.layout) {
      const { lines, texts, dots } = this.scene();
      for (const line of lines) {
        const dash = line.dash?.length ? ` stroke-dasharray="${line.dash.join(' ')}"` : '';
        svg += `  <line x1="${line.x1}" y1="${line.y1}" x2="${line.x2}" y2="${line.y2}" stroke="${line.stroke}" stroke-width="${line.width}"${dash}/>\n`;
      }
      for (const dot of dots) {
        svg += dot.shape === 'diamond'
          ? `  <polygon points="${diamondPoints(dot).map(p => p.join(',')).join(' ')}" fill="${dot.fill}"/>\n`
          : `  <circle cx="${dot.x}" cy="${dot.y}" r="${dot.r}" fill="${dot.fill}"/>\n`;
      }
      for (const text of texts) {
        const anchor = text.align === 'center' ? 'middle' : 'start';
        const baseline = text.baseline === 'middle' ? 'middle' : 'text-after-edge';
        const weight = text.bold ? ' font-weight="bold"' : '';
        svg += `  <text x="${text.x}" y="${text.y}" font-family="${escapeSvgText(BRANCH_TREE_FONT)}" font-size="${text.size}"${weight} fill="${text.fill}" text-anchor="${anchor}" dominant-baseline="${baseline}">${escapeSvgText(text.text)}</text>\n`;
      }
    }
    return svg + `</svg>`;
  }

  exportSVG() {
    downloadSVG(this.svg(), `${this.options.exportName}.svg`);
  }

  exportPNG() {
    downloadPNG(this.canvas, `${this.options.exportName}.png`);
  }

  // The branch whose horizontal edge (or leaf label) is nearest the point, within half a row
  nodeAt(mx, my) {
    let best = null;
    let bestDist = this.size.rowHeight / 2;
    for (const edge of this.edges) {
      if (mx < edge.x0 - 4 || mx > edge.x1 + 4) continue;
      const dist = Math.abs(my - edge.y);
      if (dist <= bestDist) {
        best = edge.node;
        bestDist = dist;
      }
    }
    return best;
  }

  setHovered(node) {
    if (this.hovered === node) return;
    this.hovered = node;
    this.render();
    if (this.options.onHover) this.options.onHover(node ? node.name : null, node);
  }

  setupEvents() {
    this.canvas.addEventListener('mousemove', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      this.setHovered(this.nodeAt(e.clientX - rect.left, e.clientY - rect.top));
    });

    this.canvas.addEventListener('mouseleave', () => this.setHovered(null));

    this.canvas.addEventListener('click', () => {
      if (!this.hovered) return;
      this.selected = this.hovered;
      this.render();
      if (this.options.onSelect) this.options.onSelect(this.hovered.name, this.hovered);
    });
  }
}
