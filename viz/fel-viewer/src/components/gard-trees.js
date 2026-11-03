// GARD's tree views: two trees side by side with the taxa that moved between them marked
// (GardTreePair, used by the breakpoint drawer and the tree card), and GardTreeCard, the
// card under the table with the step-up model selection, the partition-by-tree
// likelihood matrix and a tree comparison over any two segments.

import { BranchTree } from './branch-tree.js';
import { COMPOSITION_CSS_COLORS } from './composition-panel.js';
import { PALETTE } from '../data/palette.js';
import { GARD_RUN_ALPHA, gardCutoff, gardFormat, gardP } from '../data/gard-data.js';
import { downloadSVG, downloadSVGAsPNG, escapeSvgText } from '../core/export.js';

// Past this many leaves a tree is drawn as a cladogram: rows are too tight for lengths to read
export const GARD_CLADOGRAM_TAXA = 60;

const GARD_COLORS = Object.assign({}, COMPOSITION_CSS_COLORS, { accent: 'var(--accent)' });
const GARD_FILE_COLORS = Object.assign({}, PALETTE);

const GARD_MOVED_LEGEND = [{ label: 'Moved taxon (sister clade differs)', stroke: 'var(--accent)' }];
const GARD_MOVED_SHOWN = 24;

// 'S1 · 1–279', the segment a select offers
export function gardSegmentLabel(segment, index) {
  return `S${index + 1} · ${segment.start}–${segment.end}`;
}

// The segments a comparison starts on: the two around the selected breakpoint, else the
// first two; a single segment stands alone
export function gardDefaultPair(data, site) {
  const brk = data.getBreak(site);
  if (brk) return [brk.index - 1, brk.index];
  return data.segments.length > 1 ? [0, 1] : [0, null];
}

// A leaf's label, marked when its taxon moved
export function gardLeafLabel(name, moved) {
  return moved.has(name) ? `◆ ${name}` : name;
}

// The line under a tree pair: the RF distance and the moved taxa
export function gardTreeSummary(diff) {
  if (!diff) return '';
  const shown = diff.moved.slice(0, GARD_MOVED_SHOWN).map(escapeSvgText).join(', ');
  const more = diff.moved.length > GARD_MOVED_SHOWN ? ` and ${diff.moved.length - GARD_MOVED_SHOWN} more` : '';
  const moved = diff.moved.length ? `${diff.moved.length} moved: ${shown}${more}` : 'no taxon moved';
  return `Robinson–Foulds distance <strong>${diff.rf}</strong> of ${diff.max} · ${moved}`;
}

export class GardTreePair {
  constructor(container, options = {}) {
    this.options = Object.assign({ exportName: 'hyphy_gard_tree' }, options);
    this.data = null;
    this.moved = new Set();
    this.indexes = [null, null];

    this.el = document.createElement('div');
    this.el.className = 'gard-tree-pair';
    this.el.innerHTML = `
      <div class="gard-tree-cols">
        <div class="gard-tree-col" data-side="left">
          <div class="gard-tree-head"><span class="gard-tree-title" data-title="left"></span><button type="button" class="plot-tool-btn" data-svg="left" title="Export Vector SVG">SVG</button></div>
          <div class="gard-tree-host" data-tree="left"></div>
        </div>
        <div class="gard-tree-col" data-side="right">
          <div class="gard-tree-head"><span class="gard-tree-title" data-title="right"></span><button type="button" class="plot-tool-btn" data-svg="right" title="Export Vector SVG">SVG</button></div>
          <div class="gard-tree-host" data-tree="right"></div>
        </div>
      </div>
      <div class="gard-tree-summary" data-summary></div>
    `;
    container.appendChild(this.el);

    this.sides = ['left', 'right'].map((side, i) => ({
      col: this.el.querySelector(`[data-side="${side}"]`),
      title: this.el.querySelector(`[data-title="${side}"]`),
      tree: new BranchTree(this.el.querySelector(`[data-tree="${side}"]`), {
        ladderize: true,
        hoverColor: 'stone',
        style: (node, meta, colors) => (!node.children.length && this.moved.has(node.name) ? { stroke: colors.accent, width: 2 } : null),
        leafLabel: (node) => gardLeafLabel(node.name, this.moved),
        exportName: `${this.options.exportName}_${side}`,
        // The same taxon lights up in the other tree
        onHover: (name, node) => this.sides[1 - i]?.tree.highlight(node && !node.children.length ? name : null)
      })
    }));
    this.sides.forEach(({ tree }, i) => {
      this.el.querySelector(`[data-svg="${i ? 'right' : 'left'}"]`).addEventListener('click', () => tree.exportSVG());
    });
    this.summary = this.el.querySelector('[data-summary]');
  }

  get trees() {
    return this.sides.map(side => side.tree);
  }

  destroy() {
    this.trees.forEach(tree => tree.destroy());
  }

  render() {
    this.trees.forEach(tree => tree.render());
  }

  // Draw segments i and j (indexes into data.segments); j = null draws one tree
  setSegments(data, i, j = null) {
    this.data = data;
    this.indexes = [i, j];
    const diff = j == null ? null : data.compare(i, j);
    this.moved = new Set(diff?.moved ?? []);
    const cladogram = data.numTaxa > GARD_CLADOGRAM_TAXA;
    this.indexes.forEach((index, side) => {
      const { col, title, tree } = this.sides[side];
      col.style.display = index == null ? 'none' : '';
      if (index == null) return;
      const segment = data.segments[index];
      title.textContent = gardSegmentLabel(segment, index);
      tree.options.cladogram = cladogram;
      tree.setTree(segment.tree);
      tree.setLegend(side === 0 && this.moved.size ? GARD_MOVED_LEGEND : []);
      tree.options.exportName = `${this.options.exportName}_s${index + 1}`;
    });
    const scale = cladogram ? ` · cladogram above ${GARD_CLADOGRAM_TAXA} taxa (branch lengths not to scale)` : '';
    this.summary.innerHTML = j == null ? 'A single tree fits best: no breakpoints.'
      : (diff ? `${gardTreeSummary(diff)}${scale}` : `Trees cannot be compared: a tree is missing or its taxon names are not unique${scale}`);
  }
}

// Model selection: the step-up procedure's gain in AICc per breakpoint added, annotated
// with the sites each step added

// The display sites step i added to the previous step's set; none when a step lists no sites
export function gardAddedSites(steps, i) {
  const bps = steps[i]?.bps;
  if (!bps) return [];
  const before = new Set(steps[i - 1]?.bps ?? []);
  return bps.filter(b => !before.has(b));
}

// '+279', '+120, 410', '+12, 34, 56 …'
export function gardAddedLabel(sites) {
  if (!sites.length) return '';
  return `+${sites.slice(0, 3).join(', ')}${sites.length > 3 ? ' …' : ''}`;
}

// Bars of ΔAICc per step, the best model's step in ink. Returns { svg, width, height };
// colors is GARD_COLORS (CSS variables) for the page and GARD_FILE_COLORS for a download.
export function gardStepsSVG(steps, { best = steps.length, width = 480, colors = GARD_COLORS } = {}) {
  const pad = { left: 44, right: 16, top: 28, bottom: 44 };
  const height = 190;
  const plotW = width - pad.left - pad.right;
  const plotH = height - pad.top - pad.bottom;
  const max = Math.max(1e-9, ...steps.map(s => s.deltaAicc ?? 0));
  const slot = plotW / Math.max(1, steps.length);
  const barW = Math.min(56, slot * 0.6);
  const fill = (c) => `style="fill: ${c}"`;
  const text = (x, y, value, anchor, color, extra = '') =>
    `  <text x="${x}" y="${y}" text-anchor="${anchor}" font-size="10"${extra} ${fill(color)}>${escapeSvgText(value)}</text>\n`;

  let svg = `<svg xmlns="http://www.w3.org/2000/svg" class="gard-steps" role="img" aria-label="ΔAICc gained as each breakpoint is added to the best model" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="-apple-system, sans-serif">\n`;
  svg += `  <rect width="100%" height="100%" ${fill(colors.panel)}/>\n`;
  svg += `  <line x1="${pad.left}" y1="${pad.top + plotH}" x2="${width - pad.right}" y2="${pad.top + plotH}" style="stroke: ${colors.hairline}"/>\n`;
  svg += text(pad.left - 8, pad.top + 4, gardFormat(max, 0), 'end', colors.stone);
  svg += text(pad.left - 8, pad.top + plotH, '0', 'end', colors.stone);
  svg += text(4, pad.top - 14, 'ΔAICc gained', 'start', colors.stone, ' font-size="9"');
  steps.forEach((step, i) => {
    const x = pad.left + slot * i + (slot - barW) / 2;
    const h = Math.max(1, ((step.deltaAicc ?? 0) / max) * plotH);
    const isBest = step.k === best;
    const sites = gardAddedSites(steps, i);
    const tip = sites.length ? `<title>${step.k} BP: +${escapeSvgText(sites.join(', '))}</title>` : '';
    svg += `  <rect x="${x}" y="${pad.top + plotH - h}" width="${barW}" height="${h}" ${fill(isBest ? colors.ink : colors.slateLight)}>${tip}</rect>\n`;
    svg += text(x + barW / 2, pad.top + plotH - h - 5, gardFormat(step.deltaAicc, 1), 'middle', colors.ink, isBest ? ' font-weight="bold"' : '');
    svg += text(x + barW / 2, pad.top + plotH + 14, `${step.k} BP`, 'middle', isBest ? colors.ink : colors.stone, isBest ? ' font-weight="bold"' : '');
    // Fit the label to its slot (about 5px a character); the full list is in the bar's tooltip
    const room = Math.floor((slot - 4) / 5);
    let added = gardAddedLabel(sites);
    if (added.length > room) added = room >= 4 ? `${added.slice(0, room - 1).replace(/[ ,]+$/, '')}…` : '';
    if (added) svg += text(x + barW / 2, pad.top + plotH + 27, added, 'middle', colors.stone, ' font-size="9"');

  });
  svg += '</svg>';
  return { svg, width, height };
}

// Partition by tree: each partition's likelihood under every partition's tree. mode is
// 'logL', 'kh' or 'sh'; the best tree of a row is bold and a tree the row's KH or SH test
// rejects at cutoff (alpha / 2B) is shaded.
export const GARD_MATRIX_MODES = [
  { id: 'logL', label: 'log L' },
  { id: 'kh', label: 'KH p' },
  { id: 'sh', label: 'SH p' }
];

export function gardMatrixMarkup(matrix, { mode = 'logL', cutoff = 0 } = {}) {
  const trees = Math.max(...matrix.map(m => m.logL.length));
  const head = Array.from({ length: trees }, (_, j) => `<th title="The tree of partition ${j + 1}">Tree ${j + 1}</th>`).join('');
  const body = matrix.map(m => {
    const cells = Array.from({ length: trees }, (_, j) => {
      const rejected = (m.kh[j] != null && m.kh[j] <= cutoff) || (m.sh[j] != null && m.sh[j] <= cutoff);
      const classes = ['font-mono', j === m.bestTree && 'gard-matrix-best', rejected && 'gard-matrix-rejected'].filter(Boolean).join(' ');
      const value = mode === 'logL' ? gardFormat(m.logL[j], 1) : gardP((mode === 'kh' ? m.kh : m.sh)[j]);
      return `<td class="${classes}">${value}</td>`;
    }).join('');
    return `
              <tr><th>P${m.partition + 1}</th><td class="font-mono">${m.interval[0]}–${m.interval[1]}</td>${cells}</tr>`;
  }).join('');
  return `<table class="gard-matrix" aria-label="Log-likelihood or topology-test p of each partition under every partition's tree">
            <thead><tr><th>Partition</th><th>Sites</th>${head}</tr></thead>
            <tbody>${body}
            </tbody>
          </table>`;
}

function gardPanelHead(title, exportId, extra = '') {
  const tools = exportId ? `<span class="plot-export-group">
              <button type="button" class="plot-tool-btn" data-export="${exportId}" data-format="png" title="Export High-Res PNG">PNG</button>
              <button type="button" class="plot-tool-btn" data-export="${exportId}" data-format="svg" title="Export Vector SVG">SVG</button>
            </span>` : '';
  return `<div class="gard-panel-head"><span class="gard-panel-title">${title}</span>${extra}${tools}</div>`;
}

export class GardTreeCard {
  constructor(container, options = {}) {
    this.container = container;
    this.data = null;
    this.alpha = options.threshold ?? GARD_RUN_ALPHA;
    this.matrixMode = 'logL';
    this.pair = null;
    this.indexes = [0, null];

    this.el = document.createElement('div');
    this.el.className = 'hyphy-gard-card';
    this.el.innerHTML = `
      <div class="gard-card-row">
        <section class="gard-panel" id="gard-steps-panel">
          ${gardPanelHead('Model selection', 'steps')}
          <div class="gard-figure" id="gard-steps"></div>
          <div class="gard-caption" id="gard-steps-caption"></div>
        </section>
        <section class="gard-panel" id="gard-matrix-panel">
          ${gardPanelHead('Partition × tree', null, `<span class="segmented-control gard-matrix-modes">${GARD_MATRIX_MODES.map((m, i) =>
            `<button type="button" class="seg-btn${i ? '' : ' active'}" data-mode="${m.id}">${m.label}</button>`).join('')}</span>`)}
          <div class="gard-matrix-host" id="gard-matrix"></div>
          <div class="gard-caption" id="gard-matrix-caption"></div>
        </section>
      </div>
      <section class="gard-panel gard-compare">
        ${gardPanelHead('Tree comparison', null, `<span class="gard-selects" id="gard-selects">
              <select class="gard-select" id="gard-select-left" aria-label="First segment"></select>
              <span class="gard-vs">vs</span>
              <select class="gard-select" id="gard-select-right" aria-label="Second segment"></select>
            </span>`)}
        <div id="gard-pair"></div>
      </section>
    `;
    container.appendChild(this.el);

    this.selects = [this.el.querySelector('#gard-select-left'), this.el.querySelector('#gard-select-right')];
    this.selects.forEach((select, side) => select.addEventListener('change', () => {
      this.indexes[side] = Number(select.value);
      this.renderTrees();
    }));
    this.el.querySelectorAll('[data-export]').forEach(btn => {
      btn.addEventListener('click', () => this.export(btn.dataset.format));
    });
    this.el.querySelectorAll('[data-mode]').forEach(btn => {
      btn.addEventListener('click', () => {
        this.matrixMode = btn.dataset.mode;
        this.el.querySelectorAll('[data-mode]').forEach(b => b.classList.toggle('active', b === btn));
        this.renderMatrix();
      });
    });
    this.pair = new GardTreePair(this.el.querySelector('#gard-pair'), { exportName: 'hyphy_gard_compare' });
    this.resizeObserver = typeof ResizeObserver === 'function' ? new ResizeObserver(() => this.renderSteps()) : null;
    this.resizeObserver?.observe(this.el.querySelector('#gard-steps'));
  }

  destroy() {
    this.resizeObserver?.disconnect();
    this.pair.destroy();
  }

  // Repaint the canvas trees in the current theme; the SVG and the table follow the CSS tokens
  render() {
    this.pair.render();
  }

  setData(data) {
    this.data = data;
    this.indexes = gardDefaultPair(data, null);
    this.el.querySelector('#gard-matrix-panel').style.display = data.capabilities.matrix ? '' : 'none';
    // A matrix of many trees takes the full row
    this.el.querySelector('.gard-card-row').classList.toggle('gard-card-row-wide', data.matrix.length > 4);
    this.el.querySelector('#gard-selects').style.display = data.segments.length > 1 ? '' : 'none';
    this.selects.forEach(select => {
      select.innerHTML = data.segments.map((s, i) => `<option value="${i}">${gardSegmentLabel(s, i)}</option>`).join('');
    });
    this.renderSteps();
    this.renderMatrix();
    this.renderTrees();
  }

  setThreshold(alpha) {
    this.alpha = alpha;
    this.renderMatrix();
  }

  // The tree comparison follows the focused breakpoint to the two segments around it
  selectSite(site) {
    const brk = this.data?.getBreak(site);
    if (!brk) return;
    this.indexes = [brk.index - 1, brk.index];
    this.renderTrees();
  }

  figureWidth() {
    const width = this.el.querySelector('#gard-steps').clientWidth;
    return width > 0 ? width : 480;
  }

  stepsFigure(colors) {
    return gardStepsSVG(this.data.steps, { best: this.data.breaks.length, width: this.figureWidth(), colors });
  }

  renderSteps() {
    if (!this.data) return;
    const { steps } = this.data;
    this.el.querySelector('#gard-steps').innerHTML = steps.length ? this.stepsFigure(GARD_COLORS).svg : '';
    this.el.querySelector('#gard-steps-caption').textContent = steps.length
      ? 'Step-up procedure: the improvement in AICc as each breakpoint is added to the best model; sites under a bar are the breakpoints that step added. The best model is in ink.'
      : 'No breakpoint improved on the single tree (ΔAICc 0), so the step-up procedure stopped at the baseline.';
  }

  // The cutoff a partition's KH and SH tests must reach to reject a tree: alpha / 2B
  matrixCutoff() {
    return (this.alpha ?? 0) / (2 * Math.max(1, this.data.breaks.length));
  }

  renderMatrix() {
    if (!this.data?.capabilities.matrix) return;
    const cutoff = this.matrixCutoff();
    this.el.querySelector('#gard-matrix').innerHTML = gardMatrixMarkup(this.data.matrix, { mode: this.matrixMode, cutoff });
    this.el.querySelector('#gard-matrix-caption').textContent = `Each partition's likelihood under every partition's tree. Bold: the partition's best tree. Shaded: a tree its KH or SH test rejects at α / 2B = ${gardCutoff(cutoff)}.`;
  }

  renderTrees() {
    if (!this.data) return;
    const [i, j] = this.indexes;
    this.selects[0].value = String(i);
    this.selects[1].value = String(j ?? 0);
    this.pair.setSegments(this.data, i, j);
  }

  export(format) {
    if (!this.data?.steps.length) return;
    const figure = this.stepsFigure(GARD_FILE_COLORS);
    if (format === 'png') downloadSVGAsPNG(figure.svg, figure.width, figure.height, 'hyphy_gard_model_selection.png');
    else downloadSVG(figure.svg, 'hyphy_gard_model_selection.svg');
  }
}
