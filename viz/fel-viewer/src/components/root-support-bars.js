// Root support bars: candidate root placements ranked by Akaike weight, one horizontal bar
// per candidate edge. The chosen placement is matched by name (it need not be rank 1 when
// an outgroup was requested) and drawn in ink, the rest in grey; a list cut short still ends
// with the chosen candidate, under its rank. A dashed rule follows the row where the
// cumulative weight first reaches 95%. Plain SVG markup like the composition strip, themed
// through the same colour sets.

import { escapeSvgText } from '../core/export.js';
import { COMPOSITION_CSS_COLORS } from './composition-panel.js';

// 'Node1 - Node2' -> 'Node1 – Node2'
export function rootEdgeLabel(branch) {
  return String(branch ?? '').replace(/ - /g, ' – ');
}

// How many top-ranked candidates it takes for the cumulative weight to reach 95%
export function rootCredibleCount(candidates, level = 0.95) {
  let sum = 0;
  for (const [i, c] of candidates.entries()) {
    sum += c.weight ?? 0;
    if (sum >= level - 1e-9) return i + 1;
  }
  return candidates.length;
}

export function rootWeightPercent(w) {
  return `${(w * 100).toFixed(1)}%`;
}

// The rows drawn for the first `limit` candidates: [{ c, rank }] in rank order, with the
// chosen candidate appended (after a gap) when it ranks below the limit
export function rootSupportRows(candidates, chosen, limit = candidates.length) {
  const rows = candidates.slice(0, limit).map((c, rank) => ({ c, rank }));
  const rank = candidates.findIndex(c => c.branch === chosen);
  if (rank >= limit) rows.push({ c: candidates[rank], rank });
  return rows;
}

// The bars for the first `limit` candidates (rank order) and the chosen one; chosen is the
// root_branch name. Returns { svg, width, height }.
export function rootSupportSVG(candidates, { chosen, limit = candidates.length, width = 360, colors = COMPOSITION_CSS_COLORS } = {}) {
  const shown = rootSupportRows(candidates, chosen, limit);
  const rowH = 20;
  const gap = shown.length > limit ? rowH / 2 : 0;
  const pad = { left: 118, right: 128, top: 6, bottom: 6 };
  const height = pad.top + shown.length * rowH + gap + pad.bottom;
  const plotW = width - pad.left - pad.right;
  const maxW = Math.max(...candidates.map(c => c.weight ?? 0), 1e-9);
  const credible = rootCredibleCount(candidates);
  const fill = (c) => `style="fill: ${c}"`;

  let svg = `<svg xmlns="http://www.w3.org/2000/svg" class="root-support" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" font-family="-apple-system, sans-serif">\n`;
  svg += `  <rect width="100%" height="100%" ${fill(colors.panel)}/>\n`;
  shown.forEach(({ c, rank }, i) => {
    const appended = i >= limit;
    const y = pad.top + i * rowH + (appended ? gap : 0);
    const mid = y + rowH / 2;
    if (appended) {
      svg += `  <line x1="4" y1="${y - gap / 2}" x2="${width - 4}" y2="${y - gap / 2}" stroke-dasharray="1,3" style="stroke: ${colors.stone}"/>\n`;
      svg += `  <text x="4" y="${mid}" dominant-baseline="middle" font-size="9" ${fill(colors.stone)}>#${rank + 1}</text>\n`;
    }
    const isChosen = c.branch === chosen;
    const ink = isChosen ? colors.ink : colors.stone;
    const weight = isChosen ? ' font-weight="bold"' : '';
    let labelX = pad.left - 8;
    if (c.terminal) {
      svg += `  <circle cx="${pad.left - 10}" cy="${mid}" r="2.5" stroke-width="1" style="fill: none; stroke: ${ink}"/>\n`;
      labelX -= 8;
    }
    svg += `  <text x="${labelX}" y="${mid}" text-anchor="end" dominant-baseline="middle" font-size="10"${weight} ${fill(ink)}>${escapeSvgText(rootEdgeLabel(c.branch))}</text>\n`;
    const barW = Math.max(1, ((c.weight ?? 0) / maxW) * plotW);
    svg += `  <rect x="${pad.left}" y="${y + 5}" width="${barW}" height="${rowH - 10}" ${fill(isChosen ? colors.ink : colors.slateLight)}/>\n`;
    const delta = c.deltaLogL == null ? '' : ` · ΔlnL ${c.deltaLogL.toFixed(2)}`;
    svg += `  <text x="${pad.left + barW + 6}" y="${mid}" dominant-baseline="middle" font-size="10"${weight} ${fill(ink)}>w = ${rootWeightPercent(c.weight ?? 0)}${delta}</text>\n`;
    if (!appended && rank + 1 === credible && credible < candidates.length) {
      const ry = y + rowH;
      svg += `  <line x1="4" y1="${ry}" x2="${width - 4}" y2="${ry}" stroke-dasharray="3,3" style="stroke: ${colors.stone}"/>\n`;
      svg += `  <text x="4" y="${ry - 2}" font-size="9" ${fill(colors.stone)}>95%</text>\n`;
    }
  });
  svg += '</svg>';
  return { svg, width, height };
}
