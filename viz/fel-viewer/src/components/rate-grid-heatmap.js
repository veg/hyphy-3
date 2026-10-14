// Rate grid heatmap: weights over an N x N alpha/beta grid (FUBAR), drawn as an ORDINAL
// grid of equal cells because the grid values are not evenly spaced. Tick labels give the
// real alpha (x) and beta (y) values at a subset of indices, fill is a sqrt greyscale of
// the weight with its own ramp legend, a dashed diagonal marks beta = alpha and the two
// half planes carry faint 'β > α' / 'α > β' labels. Hovering a cell reads it out under
// the canvas; PNG and SVG export sit in the header.
//
// options: { title, caption, size (plot side in px), grid: { alpha[], beta[], dims },
//   weights[] (alpha-major: g = iAlpha * dims + iBeta), mode: outline the heaviest cell,
//   mean: { alpha, beta } marked by interpolation on the ordinal axes, exportName }

import { PALETTE, withAlpha } from '../data/palette.js';
import { downloadPNG, downloadSVG, escapeSvgText } from '../core/export.js';

export const RATE_GRID_MARGIN = { left: 44, right: 58, top: 8, bottom: 30 };

// The distinct axis values: alpha changes along the major index, beta along the minor
export function rateGridAxes(grid) {
  const n = grid.dims;
  return {
    alpha: Array.from({ length: n }, (_, i) => grid.alpha[i * n]),
    beta: Array.from({ length: n }, (_, j) => grid.beta[j])
  };
}

// Fractional ordinal index of value v on ascending axis values, by linear interpolation
// between neighbours (clamped to the ends)
export function ordinalIndex(values, v) {
  if (v <= values[0]) return 0;
  for (let i = 1; i < values.length; i++) {
    if (v <= values[i]) return i - 1 + (v - values[i - 1]) / (values[i] - values[i - 1]);
  }
  return values.length - 1;
}

// Grey for a weight on a sqrt ramp from the panel (0) to near-ink (the maximum)
export function rateGridShade(w, max) {
  const level = Math.round(248 - Math.sqrt(Math.max(0, w) / (max || 1)) * 214);
  return `rgb(${level}, ${level}, ${level})`;
}

// Indices to label on ascending axis values. Always the ends, the neutral rate 1 (or
// the value nearest it) and the middle of the range above it, where the posterior mass
// of most sites sits; then, left to right, any index at least minPx from every label.
export function rateGridTicks(values, cell, minPx = 28) {
  const n = values.length;
  const every = Math.max(1, Math.ceil(minPx / cell));
  let one = 0;
  values.forEach((v, i) => {
    if (Math.abs(v - 1) < Math.abs(values[one] - 1)) one = i;
  });
  const ticks = new Set([0, n - 1, one, Math.round((one + n - 1) / 2)]);
  for (let i = 1; i < n - 1; i++) {
    if ([...ticks].every(t => Math.abs(t - i) >= every)) ticks.add(i);
  }
  return [...ticks].sort((a, b) => a - b);
}

function rateGridValue(v) {
  if (v === 0) return '0';
  return v >= 10 ? v.toFixed(0) : v.toPrecision(2);
}

function rateGridPercent(w) {
  return `${(w * 100).toFixed(w >= 0.01 ? 1 : 2)}%`;
}

export class RateGridHeatmap {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      title: '',
      caption: '',
      size: 200,
      grid: null,
      weights: [],
      mode: false,
      mean: null,
      exportName: 'hyphy_rate_grid'
    }, options);

    const { size } = this.options;
    const m = RATE_GRID_MARGIN;
    this.width = m.left + size + m.right;
    this.height = m.top + size + m.bottom;
    this.n = this.options.grid.dims;
    this.cell = size / this.n;
    this.axes = rateGridAxes(this.options.grid);
    this.max = Math.max(...this.options.weights);
    this.modeIndex = this.options.weights.indexOf(this.max);

    this.initDOM();
    this.draw();
  }

  initDOM() {
    this.el = document.createElement('div');
    this.el.className = 'rate-grid-heatmap';
    this.el.style.width = `${this.width}px`;
    this.el.innerHTML = `
      <div class="rate-grid-head">
        <span class="rate-grid-title">${this.options.title}</span>
        <span class="plot-export-group">
          <button type="button" class="plot-tool-btn" data-export="png" title="Export High-Res PNG">PNG</button>
          <button type="button" class="plot-tool-btn" data-export="svg" title="Export Vector SVG">SVG</button>
        </span>
      </div>
      <canvas class="rate-grid-canvas"></canvas>
      <div class="rate-grid-readout"></div>
    `;
    this.container.appendChild(this.el);

    this.canvas = this.el.querySelector('canvas');
    this.readout = this.el.querySelector('.rate-grid-readout');
    this.readout.textContent = this.options.caption;
    const dpr = window.devicePixelRatio || 1;
    this.dpr = dpr;
    this.canvas.width = this.width * dpr;
    this.canvas.height = this.height * dpr;
    this.canvas.style.width = `${this.width}px`;
    this.canvas.style.height = `${this.height}px`;

    this.el.querySelector('[data-export="png"]').addEventListener('click', () => downloadPNG(this.canvas, `${this.options.exportName}.png`));
    this.el.querySelector('[data-export="svg"]').addEventListener('click', () => downloadSVG(this.svg(), `${this.options.exportName}.svg`));
    this.canvas.addEventListener('mousemove', (e) => this.hover(e));
    this.canvas.addEventListener('mouseleave', () => {
      this.readout.textContent = this.options.caption;
    });
  }

  // Top-left corner of cell (i alpha, j beta); beta grows upward
  cellAt(i, j) {
    const m = RATE_GRID_MARGIN;
    return { x: m.left + i * this.cell, y: m.top + (this.n - 1 - j) * this.cell };
  }

  // The marks every rendering shares, in plot coordinates
  marks() {
    const m = RATE_GRID_MARGIN;
    const { size } = this.options;
    const n = this.n;
    const cells = this.options.weights.map((w, g) => {
      const i = Math.floor(g / n);
      const j = g % n;
      return Object.assign(this.cellAt(i, j), { fill: rateGridShade(w, this.max) });
    });

    const mean = this.options.mean && {
      x: m.left + (ordinalIndex(this.axes.alpha, this.options.mean.alpha) + 0.5) * this.cell,
      y: m.top + (n - 0.5 - ordinalIndex(this.axes.beta, this.options.mean.beta)) * this.cell
    };
    const mode = this.options.mode && this.cellAt(Math.floor(this.modeIndex / n), this.modeIndex % n);
    return {
      cells,
      mode,
      mean,
      xTicks: rateGridTicks(this.axes.alpha, this.cell).map(i => ({ x: m.left + (i + 0.5) * this.cell, label: rateGridValue(this.axes.alpha[i]) })),
      yTicks: rateGridTicks(this.axes.beta, this.cell).map(j => ({ y: m.top + (n - 0.5 - j) * this.cell, label: rateGridValue(this.axes.beta[j]) })),
      plot: { left: m.left, top: m.top, right: m.left + size, bottom: m.top + size },
      ramp: { x: m.left + size + 10, top: m.top, bottom: m.top + size, width: 8 }
    };
  }

  draw() {
    const ctx = this.canvas.getContext('2d');
    const k = this.marks();
    const { plot, ramp } = k;
    ctx.save();
    ctx.scale(this.dpr, this.dpr);
    ctx.clearRect(0, 0, this.width, this.height);

    for (const c of k.cells) {
      ctx.fillStyle = c.fill;
      ctx.fillRect(c.x, c.y, this.cell + 0.5, this.cell + 0.5);
    }
    ctx.strokeStyle = PALETTE.hairline;
    ctx.lineWidth = 1;
    ctx.strokeRect(plot.left, plot.top, plot.right - plot.left, plot.bottom - plot.top);

    // beta = alpha through the diagonal cells' centres
    ctx.save();
    ctx.setLineDash([4, 3]);
    ctx.strokeStyle = PALETTE.stone;
    ctx.beginPath();
    ctx.moveTo(plot.left, plot.bottom);
    ctx.lineTo(plot.right, plot.top);
    ctx.stroke();
    ctx.restore();

    ctx.font = 'bold 10px "SF Mono", Monaco, monospace';
    ctx.textBaseline = 'top';
    ctx.textAlign = 'left';
    ctx.fillStyle = withAlpha(PALETTE.accent, 0.7);
    ctx.fillText('β > α', plot.left + 6, plot.top + 6);
    ctx.textAlign = 'right';
    ctx.textBaseline = 'bottom';
    ctx.fillStyle = withAlpha(PALETTE.purifying, 0.7);
    ctx.fillText('α > β', plot.right - 6, plot.bottom - 6);

    // Mode and mean in ink on a knockout halo: neutral marks that read on either half
    if (k.mode) {
      ctx.strokeStyle = PALETTE.knockout;
      ctx.lineWidth = 3.5;
      ctx.strokeRect(k.mode.x, k.mode.y, this.cell, this.cell);
      ctx.strokeStyle = PALETTE.ink;
      ctx.lineWidth = 1.5;
      ctx.strokeRect(k.mode.x, k.mode.y, this.cell, this.cell);
    }
    if (k.mean) {
      ctx.strokeStyle = PALETTE.knockout;
      ctx.lineWidth = 3.5;
      this.cross(ctx, k.mean);
      ctx.strokeStyle = PALETTE.ink;
      ctx.lineWidth = 1.5;
      this.cross(ctx, k.mean);
    }

    ctx.font = '10px "SF Mono", Monaco, monospace';
    ctx.fillStyle = PALETTE.stone;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const t of k.xTicks) ctx.fillText(t.label, t.x, plot.bottom + 4);
    ctx.fillText('α (synonymous)', (plot.left + plot.right) / 2, plot.bottom + 17);
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    for (const t of k.yTicks) ctx.fillText(t.label, plot.left - 4, t.y);
    ctx.save();
    ctx.translate(6, (plot.top + plot.bottom) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center';
    ctx.fillText('β (non-syn.)', 0, 0);
    ctx.restore();

    // Ramp legend: 0 at the bottom, the heaviest cell at the top
    const steps = 32;
    const h = (ramp.bottom - ramp.top) / steps;
    for (let s = 0; s < steps; s++) {
      ctx.fillStyle = rateGridShade(((s + 0.5) / steps) ** 2 * this.max, this.max);
      ctx.fillRect(ramp.x, ramp.bottom - (s + 1) * h, ramp.width, h + 0.5);
    }
    ctx.strokeStyle = PALETTE.hairline;
    ctx.strokeRect(ramp.x, ramp.top, ramp.width, ramp.bottom - ramp.top);
    ctx.fillStyle = PALETTE.stone;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    ctx.fillText(rateGridPercent(this.max), ramp.x + ramp.width + 3, ramp.top);
    ctx.textBaseline = 'bottom';
    ctx.fillText('0', ramp.x + ramp.width + 3, ramp.bottom);
    ctx.restore();
  }

  // A diagonal ✕, as the caption names it
  cross(ctx, { x, y }) {
    ctx.beginPath();
    ctx.moveTo(x - 4, y - 4);
    ctx.lineTo(x + 4, y + 4);
    ctx.moveTo(x - 4, y + 4);
    ctx.lineTo(x + 4, y - 4);
    ctx.stroke();
  }

  // The cell under the pointer, read out as its alpha, beta and weight
  hover(e) {
    const rect = this.canvas.getBoundingClientRect();
    const m = RATE_GRID_MARGIN;
    const i = Math.floor((e.clientX - rect.left - m.left) / this.cell);
    const j = this.n - 1 - Math.floor((e.clientY - rect.top - m.top) / this.cell);
    if (i < 0 || j < 0 || i >= this.n || j >= this.n) {
      this.readout.textContent = this.options.caption;
      return;
    }
    const w = this.options.weights[i * this.n + j];
    this.readout.textContent = `α = ${rateGridValue(this.axes.alpha[i])}, β = ${rateGridValue(this.axes.beta[j])}: ${rateGridPercent(w)} of the mass`;
  }

  svg() {
    const k = this.marks();
    const { plot, ramp } = k;
    const text = (x, y, s, attrs = '') => `  <text x="${x}" y="${y}"${attrs}>${escapeSvgText(s)}</text>\n`;
    let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${this.width}" height="${this.height}" viewBox="0 0 ${this.width} ${this.height}">\n`;
    svg += `  <style>\n    text { font-family: -apple-system, sans-serif; font-size: 10px; fill: ${PALETTE.stone}; }\n  </style>\n`;
    svg += `  <rect width="100%" height="100%" fill="${PALETTE.panel}"/>\n`;
    for (const c of k.cells) {
      svg += `  <rect x="${c.x}" y="${c.y}" width="${this.cell}" height="${this.cell}" fill="${c.fill}"/>\n`;
    }
    svg += `  <rect x="${plot.left}" y="${plot.top}" width="${plot.right - plot.left}" height="${plot.bottom - plot.top}" fill="none" stroke="${PALETTE.hairline}"/>\n`;
    svg += `  <line x1="${plot.left}" y1="${plot.bottom}" x2="${plot.right}" y2="${plot.top}" stroke="${PALETTE.stone}" stroke-dasharray="4,3"/>\n`;
    svg += text(plot.left + 6, plot.top + 14, 'β > α', ` style="fill: ${PALETTE.accent}; font-weight: bold;" opacity="0.7"`);
    svg += text(plot.right - 6, plot.bottom - 6, 'α > β', ` style="fill: ${PALETTE.purifying}; font-weight: bold;" opacity="0.7" text-anchor="end"`);
    if (k.mode) {
      const rect = `x="${k.mode.x}" y="${k.mode.y}" width="${this.cell}" height="${this.cell}" fill="none"`;
      svg += `  <rect ${rect} stroke="${PALETTE.knockout}" stroke-width="3.5"/>\n`;
      svg += `  <rect ${rect} stroke="${PALETTE.ink}" stroke-width="1.5"/>\n`;
    }
    if (k.mean) {
      const { x, y } = k.mean;
      const d = `M${x - 4},${y - 4}L${x + 4},${y + 4}M${x - 4},${y + 4}L${x + 4},${y - 4}`;
      svg += `  <path d="${d}" stroke="${PALETTE.knockout}" stroke-width="3.5"/>\n`;
      svg += `  <path d="${d}" stroke="${PALETTE.ink}" stroke-width="1.5"/>\n`;
    }
    for (const t of k.xTicks) svg += text(t.x, plot.bottom + 12, t.label, ' text-anchor="middle"');
    svg += text((plot.left + plot.right) / 2, plot.bottom + 25, 'α (synonymous)', ' text-anchor="middle"');
    for (const t of k.yTicks) svg += text(plot.left - 4, t.y + 3, t.label, ' text-anchor="end"');
    svg += text(0, 0, 'β (non-syn.)', ` text-anchor="middle" transform="translate(12, ${(plot.top + plot.bottom) / 2}) rotate(-90)"`);
    svg += `  <defs><linearGradient id="rate-grid-ramp" x1="0" y1="1" x2="0" y2="0">`
      + [0, 0.25, 0.5, 0.75, 1].map(f => `<stop offset="${f}" stop-color="${rateGridShade(f * f * this.max, this.max)}"/>`).join('')
      + `</linearGradient></defs>\n`;
    svg += `  <rect x="${ramp.x}" y="${ramp.top}" width="${ramp.width}" height="${ramp.bottom - ramp.top}" fill="url(#rate-grid-ramp)" stroke="${PALETTE.hairline}"/>\n`;
    svg += text(ramp.x + ramp.width + 3, ramp.top + 8, rateGridPercent(this.max));
    svg += text(ramp.x + ramp.width + 3, ramp.bottom, '0');
    svg += `</svg>`;
    return svg;
  }
}
