// Composition panel: the Pearson chi-square screen for base-composition homogeneity across
// taxa (overall and by codon position) and a per-taxon GC dot strip. Both are plain markup,
// so the strip is drawn once inline (themed through CSS custom properties) and once with
// literal colours for export.
//
// colors: { ink, stone, hairline, gridLine, slateLight, panel }, either CSS var() strings
// (COMPOSITION_CSS_COLORS) for the page or PALETTE values for a downloaded file.

import { escapeSvgText } from '../core/export.js';

export const COMPOSITION_CSS_COLORS = {
  ink: 'var(--ink)',
  stone: 'var(--stone)',
  hairline: 'var(--hairline)',
  gridLine: 'var(--grid-line)',
  slateLight: 'var(--slate-light)',
  panel: 'var(--panel)'
};

// p-values as the screen prints them: three decimals, exponential below 0.001
export function compositionP(p) {
  if (p == null) return '—';
  return p < 0.001 ? p.toExponential(1) : p.toFixed(3);
}

// The four screen rows. Overall and CP3 decide whether regimes are activated (•); a p
// below the heterogeneity threshold is bold.
export function compositionTableMarkup(composition, threshold) {
  const rows = [
    ['Overall', composition.overall, true],
    ['CP1', composition.cp1, false],
    ['CP2', composition.cp2, false],
    ['CP3', composition.cp3, true]
  ].filter(([, test]) => test).map(([label, test, decides]) => {
    const strong = threshold != null && test.p != null && test.p < threshold;
    return `
              <tr><td>${label}${decides ? ' •' : ''}</td><td class="font-mono">${test.chi2 == null ? '—' : test.chi2.toFixed(2)}</td><td class="font-mono">${test.df ?? '—'}</td><td class="font-mono${strong ? ' font-bold' : ''}">${compositionP(test.p)}</td></tr>`;
  }).join('');
  return `
          <table class="diag-fits-table composition-table">
            <thead><tr><th>Screen</th><th><span class="case-literal">χ²</span></th><th>df</th><th>p</th></tr></thead>
            <tbody>${rows}
            </tbody>
          </table>`;
}

function compositionMean(values) {
  return values.reduce((sum, v) => sum + v, 0) / Math.max(1, values.length);
}

// Taxa whose GC3 stands out: beyond 2 SD of the mean, or among the two highest and two
// lowest when there are more than four taxa (every taxon tied with the second value too)
export function gcOutliers(taxa) {
  const values = taxa.map(t => t.gcCP3);
  const mean = compositionMean(values);
  const sd = Math.sqrt(compositionMean(values.map(v => (v - mean) ** 2)));
  const sorted = [...values].sort((a, b) => a - b);
  const end = (v) => taxa.length > 4 && (v <= sorted[1] || v >= sorted[sorted.length - 2]);
  return new Set(taxa.filter(t => Math.abs(t.gcCP3 - mean) > 2 * sd || end(t.gcCP3)).map(t => t.taxon));
}

// The axis: every plotted value with 0.02 either side
export function gcDomain(taxa, marks = []) {
  const values = [...taxa.flatMap(t => [t.gcCP3, t.gcOverall]), ...marks].filter(v => v != null);
  return [Math.min(...values) - 0.02, Math.max(...values) + 0.02];
}

// One row per taxon, sorted by GC3: a filled dot at GC3 and a hollow dot at overall GC,
// joined by a thin line; a dashed rule at the mean GC3; short labelled hairline ticks above
// the rows for regime mean GC (marks: [{ value, label }]). Returns { svg, width, height }.
export function gcStripSVG(taxa, { width = 360, marks = [], colors = COMPOSITION_CSS_COLORS } = {}) {
  const rows = taxa.filter(t => t.gcCP3 != null).sort((a, b) => a.gcCP3 - b.gcCP3);
  const compact = rows.length > 40;
  const rowH = compact ? 6 : 16;
  const pad = { left: compact ? 12 : 84, right: 44, top: marks.length ? 30 : 14, bottom: 26 };
  const height = pad.top + rows.length * rowH + pad.bottom;
  const [lo, hi] = gcDomain(rows, marks.map(m => m.value));
  const plotW = width - pad.left - pad.right;
  const x = (v) => pad.left + ((v - lo) / (hi - lo)) * plotW;
  const y = (i) => pad.top + i * rowH + rowH / 2;
  const fill = (c) => `style="fill: ${c}"`;
  const stroke = (c) => `style="stroke: ${c}"`;
  const text = (tx, ty, value, anchor, color, extra = '') =>
    `  <text x="${tx}" y="${ty}" text-anchor="${anchor}" dominant-baseline="middle" font-size="10"${extra} ${fill(color)}>${escapeSvgText(value)}</text>\n`;
  const outliers = gcOutliers(rows);
  const mean = compositionMean(rows.map(t => t.gcCP3));
  const bottom = pad.top + rows.length * rowH;

  let svg = `<svg xmlns="http://www.w3.org/2000/svg" class="gc-strip" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="-apple-system, sans-serif">\n`;
  svg += `  <rect width="100%" height="100%" ${fill(colors.panel)}/>\n`;
  const step = hi - lo > 0.3 ? 0.1 : 0.05;
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) {
    svg += `  <line x1="${x(v)}" y1="${pad.top}" x2="${x(v)}" y2="${bottom}" ${stroke(colors.gridLine)}/>\n`;
    svg += text(x(v), bottom + 10, v.toFixed(2), 'middle', colors.stone);
  }
  // A label that would overlap the previous one moves up a line
  let last = { x: -Infinity, half: 0, raised: true };
  for (const m of [...marks].sort((a, b) => a.value - b.value)) {
    const mx = x(m.value);
    const half = m.label.length * 2.8;
    const raised = mx - half < last.x + last.half + 4 && !last.raised;
    last = { x: mx, half, raised };
    svg += `  <line x1="${mx}" y1="${pad.top - 8}" x2="${mx}" y2="${pad.top - 2}" ${stroke(colors.stone)}/>\n`;
    svg += text(mx, pad.top - (raised ? 24 : 14), m.label, 'middle', colors.stone);
  }
  svg += `  <line x1="${x(mean)}" y1="${pad.top}" x2="${x(mean)}" y2="${bottom}" stroke-dasharray="3,3" ${stroke(colors.stone)}/>\n`;
  rows.forEach((t, i) => {
    const r = compact ? 2 : 3.5;
    if (t.gcOverall != null) {
      svg += `  <line x1="${x(t.gcOverall)}" y1="${y(i)}" x2="${x(t.gcCP3)}" y2="${y(i)}" ${stroke(colors.hairline)}/>\n`;
      svg += `  <circle cx="${x(t.gcOverall)}" cy="${y(i)}" r="${r}" stroke-width="1.25" style="fill: ${colors.panel}; stroke: ${colors.slateLight}"/>\n`;
    }
    svg += `  <circle cx="${x(t.gcCP3)}" cy="${y(i)}" r="${r}" ${fill(colors.ink)}/>\n`;
    if (!compact) svg += text(pad.left - 8, y(i), t.taxon, 'end', colors.ink);
    if (outliers.has(t.taxon)) {
      const label = compact ? `${t.taxon} ${t.gcCP3.toFixed(2)}` : t.gcCP3.toFixed(2);
      svg += text(x(t.gcCP3) + 7, y(i), label, 'start', colors.ink, ' font-weight="bold"');
    }
  });
  svg += text(pad.left + plotW / 2, bottom + 22, 'GC content', 'middle', colors.stone);
  svg += '</svg>';
  return { svg, width, height };
}
