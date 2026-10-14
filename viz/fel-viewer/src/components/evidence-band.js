// Evidence band: one mark per site on a y-axis, coloured by the site's call at the live
// threshold. Serves p-style evidence (rule at -log10 t), posterior-style evidence (rule
// at t, direction 'ge') and a mirrored two-sided layout (SLAC: one series above zero,
// one below).

import { PALETTE } from '../data/palette.js';
import { formatThreshold } from '../core/format.js';
import { drawThresholdRule, drawSelectionRing, drawRuleLabel } from './site-track.js';

// Render order: neutral first, purifying second, diversifying on top
const EVIDENCE_TONE_ORDER = { div: 2, pur: 1 };

// Height of the threshold rule for a cutoff t
export function evidenceRuleValue(direction, t) {
  return direction === 'ge' ? t : -Math.log10(Math.max(1e-15, t));
}

// [min, max] padded around the plotted values and always including zero;
// symmetric about zero for a two-sided band
export function autoEvidenceDomain(values, symmetric) {
  let lo = 0;
  let hi = 0;
  for (const v of values) {
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  const margin = (hi - lo) * 0.1 || 1;
  if (symmetric) {
    const extent = Math.ceil(Math.max(-lo, hi) + margin);
    return [-extent, extent];
  }
  return [lo < 0 ? Math.floor(lo - margin) : 0, Math.ceil(hi + margin)];
}

export class EvidenceBand {
  constructor(options = {}) {
    this.options = Object.assign({
      value: (r) => r.negLogP,   // (row, track) => y; null leaves the site out
      tone: () => null,          // (row, track) => 'div' | 'pur' | null at the live threshold
      mirror: null,              // { value, tone }: a second series drawn downward from zero
      domain: 'auto',            // 'auto' | [min, max] | (rows, track) => [min, max]
      ticks: 5,                  // number of tick intervals, the tick values, or (min, max) => values
      tickFormat: (v) => v.toFixed(1),
      rule: 'le',                // 'le': p-style rule at -log10 t; 'ge': posterior-style rule at t; null: none
      mirrorRule: false,         // also draw the rule below zero, in the purifying hue
      ruleLabel: null,           // t => text; defaults to 'p = t' or 'P ≥ t'
      ruleLabelAt: 'start',      // 'start': on the rule at the left; 'end': at the right, above
                                 // the +t rule and below a mirrored one, drawn over the marks
      label: null,               // y caption, or [above zero, below zero] for a two-sided band
      shape: 'triangle',         // 'triangle' | 'dot'
      tooltip: null              // (row, track, point) => HTML; point.series is 1 below zero
    }, options);
    this.points = [];
  }

  get mirrored() {
    return Boolean(this.options.mirror);
  }

  domain(rows, track) {
    const { domain } = this.options;
    if (Array.isArray(domain)) return domain;
    if (typeof domain === 'function') return domain(rows, track);
    const values = [];
    for (const r of rows) {
      const v = this.options.value(r, track);
      if (Number.isFinite(v)) values.push(v);
      if (this.mirrored) {
        const m = this.options.mirror.value(r, track);
        if (Number.isFinite(m)) values.push(-m);
      }
    }
    return autoEvidenceDomain(values, this.mirrored);
  }

  tickValues() {
    const { ticks } = this.options;
    if (Array.isArray(ticks)) return ticks;
    if (typeof ticks === 'function') return ticks(this.yMin, this.yMax);
    return Array.from({ length: ticks + 1 }, (_, i) => this.yMin + (i / ticks) * (this.yMax - this.yMin));
  }

  ruleLabel(t) {
    if (this.options.ruleLabel) return this.options.ruleLabel(t);
    return this.options.rule === 'ge' ? `P ≥ ${formatThreshold(t)}` : `p = ${formatThreshold(t)}`;
  }

  // Sets the y-scale on the frame and draws the y grid, zero line and caption
  drawBackground(ctx, frame, track) {
    [this.yMin, this.yMax] = this.domain(track.rows, track);
    const { yMin, yMax } = this;
    frame.scaleY = (val) => frame.top + (1 - (val - yMin) / (yMax - yMin)) * frame.height;
    const right = frame.left + frame.width;

    ctx.strokeStyle = PALETTE.gridLine;
    ctx.lineWidth = 1;
    ctx.font = '10px "SF Mono", Monaco, monospace';
    ctx.fillStyle = PALETTE.stone;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';

    for (const val of this.tickValues()) {
      const y = frame.scaleY(val);
      ctx.beginPath();
      ctx.moveTo(frame.left, y);
      ctx.lineTo(right, y);
      ctx.stroke();

      ctx.fillText(this.options.tickFormat(val), frame.left - 8, y);
    }

    if (yMin < 0 && yMax > 0) {
      const zeroY = frame.scaleY(0);
      ctx.strokeStyle = PALETTE.hairline;
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(frame.left, zeroY);
      ctx.lineTo(right, zeroY);
      ctx.stroke();
    }

    this.drawCaption(ctx, frame);
  }

  // Rotated y caption in the left gutter; a pair splits at zero
  drawCaption(ctx, frame) {
    const { label } = this.options;
    if (!label) return;
    const zeroY = Math.min(Math.max(frame.scaleY(0), frame.top), frame.top + frame.height);
    const spans = Array.isArray(label)
      ? [[label[0], frame.top, zeroY], [label[1], zeroY, frame.top + frame.height]]
      : [[label, frame.top, frame.top + frame.height]];

    ctx.save();
    ctx.fillStyle = PALETTE.stone;
    ctx.font = '9px "SF Mono", Monaco, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'top';
    for (const [text, from, to] of spans) {
      ctx.save();
      ctx.translate(4, (from + to) / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.fillText(text, 0, 0);
      ctx.restore();
    }
    ctx.restore();
  }

  // Screen points for every plotted value, in row order (a mirrored row gives two:
  // series 0 above zero, series 1 below)
  computePoints(frame, track) {
    const points = [];
    const series = [[this.options, 1]];
    if (this.mirrored) series.push([this.options.mirror, -1]);
    for (const r of track.rows) {
      for (const [index, [s, sign]] of series.entries()) {
        const v = s.value(r, track);
        if (v == null) continue;
        const yVal = sign * v;
        points.push({
          site: r.site,
          x: frame.scaleX(r.site),
          y: frame.scaleY(yVal),
          yVal,
          raw: r,
          series: index,
          tone: s.tone(r, track)
        });
      }
    }
    return points;
  }

  // The one point the selection ring goes on: the selected site's point in the series the
  // selection came from, or its first point
  selectedPoint(track) {
    const candidates = this.points.filter(p => p.site === track.selectedSite);
    return candidates.find(p => p.series === track.selectedSeries) ?? candidates[0] ?? null;
  }

  get hasMirrorRule() {
    return this.options.mirrorRule || this.mirrored;
  }

  drawRules(ctx, frame, track) {
    if (!this.options.rule) return;
    const at = evidenceRuleValue(this.options.rule, track.threshold);
    const startLabel = this.options.ruleLabelAt === 'start' ? this.ruleLabel(track.threshold) : null;
    drawThresholdRule(ctx, frame, frame.scaleY(at), PALETTE.accent, startLabel);
    if (this.hasMirrorRule) {
      drawThresholdRule(ctx, frame, frame.scaleY(-at), PALETTE.purifying);
    }
  }

  // The 'end' labels: [y, colour, side of the rule] for the rule and its mirror. The
  // mirror's label goes below its rule unless that would run past the plot's floor.
  endRuleLabels(frame, track) {
    if (!this.options.rule || this.options.ruleLabelAt !== 'end') return [];
    const at = evidenceRuleValue(this.options.rule, track.threshold);
    const labels = [[frame.scaleY(at), PALETTE.accent, 'above']];
    if (this.hasMirrorRule) {
      const y = frame.scaleY(-at);
      labels.push([y, PALETTE.purifying, y + 17 <= frame.top + frame.height ? 'below' : 'above']);
    }
    return labels;
  }

  drawMarks(ctx, frame, track) {
    this.drawRules(ctx, frame, track);
    this.points = this.computePoints(frame, track);

    const order = this.points.slice().sort((a, b) =>
      (EVIDENCE_TONE_ORDER[a.tone] ?? 0) - (EVIDENCE_TONE_ORDER[b.tone] ?? 0));
    const zeroY = frame.scaleY(0);
    const selected = this.selectedPoint(track);

    for (const pt of order) {
      const isSelected = (pt === selected);
      const isHovered = (track.hoveredSite === pt.site);
      if (this.mirrored) this.drawStem(ctx, pt, zeroY);
      if (this.options.shape === 'dot') {
        this.drawDot(ctx, pt, isSelected, isHovered);
      } else {
        this.drawTriangle(ctx, pt, isSelected, isHovered);
      }
      if (isSelected) drawSelectionRing(ctx, pt.x, pt.y);
    }
    const text = this.ruleLabel(track.threshold);
    for (const [y, color, side] of this.endRuleLabels(frame, track)) {
      drawRuleLabel(ctx, frame, y, color, text, side);
    }
    return this.points;
  }

  drawStem(ctx, pt, zeroY) {
    ctx.strokeStyle = pt.tone === 'div' ? PALETTE.accent : (pt.tone === 'pur' ? PALETTE.purifying : PALETTE.gridLine);
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(pt.x, zeroY);
    ctx.lineTo(pt.x, pt.y);
    ctx.stroke();
  }

  // ▲ diversifying, ▼ purifying, small grey dot otherwise
  drawTriangle(ctx, pt, isSelected, isHovered) {
    if (pt.tone === 'div') {
      const size = isSelected ? 8 : (isHovered ? 7 : 5.5);
      ctx.fillStyle = PALETTE.accent;
      ctx.strokeStyle = PALETTE.accentDeep;
      ctx.lineWidth = isSelected ? 2 : 1;

      ctx.beginPath();
      ctx.moveTo(pt.x, pt.y - size * 1.2);
      ctx.lineTo(pt.x - size, pt.y + size * 0.8);
      ctx.lineTo(pt.x + size, pt.y + size * 0.8);
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    } else if (pt.tone === 'pur') {
      const size = isSelected ? 7 : (isHovered ? 6 : 4.5);
      ctx.fillStyle = PALETTE.purifying;
      ctx.strokeStyle = PALETTE.purifyingDeep;
      ctx.lineWidth = isSelected ? 2 : 1;

      ctx.beginPath();
      ctx.moveTo(pt.x, pt.y + size * 1.2);
      ctx.lineTo(pt.x - size, pt.y - size * 0.8);
      ctx.lineTo(pt.x + size, pt.y - size * 0.8);
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    } else {
      this.drawNeutralDot(ctx, pt, isSelected, isHovered);
    }
  }

  // Filled dot in the call's hue (posterior-style plots)
  drawDot(ctx, pt, isSelected, isHovered) {
    if (!pt.tone) {
      this.drawNeutralDot(ctx, pt, isSelected, isHovered);
      return;
    }
    ctx.fillStyle = pt.tone === 'div' ? PALETTE.accent : PALETTE.purifying;
    ctx.beginPath();
    ctx.arc(pt.x, pt.y, isSelected ? 6 : (isHovered ? 5 : 4), 0, Math.PI * 2);
    ctx.fill();
  }

  drawNeutralDot(ctx, pt, isSelected, isHovered) {
    const rSize = isSelected ? 5 : (isHovered ? 4 : 2.5);
    ctx.fillStyle = isHovered ? PALETTE.slate : PALETTE.slateLight;
    ctx.beginPath();
    ctx.arc(pt.x, pt.y, rSize, 0, Math.PI * 2);
    ctx.fill();
  }

  tooltip(row, track, point) {
    return this.options.tooltip?.(row, track, point) ?? '';
  }

  svg(frame, track) {
    let svg = '';
    const line = (y, color) =>
      `  <line x1="${frame.left}" y1="${y}" x2="${frame.left + frame.width}" y2="${y}" stroke="${color}" stroke-dasharray="4,4" stroke-width="1.5"/>\n`;

    if (this.options.rule) {
      const at = evidenceRuleValue(this.options.rule, track.threshold);
      const threshY = frame.scaleY(at);
      svg += line(threshY, PALETTE.accent);
      if (this.options.ruleLabelAt === 'start') {
        svg += `  <text x="${frame.left + 6}" y="${threshY - 4}" fill="${PALETTE.accent}" font-weight="bold">${this.ruleLabel(track.threshold)}</text>\n`;
      }
      if (this.hasMirrorRule) svg += line(frame.scaleY(-at), PALETTE.purifying);
    }

    for (const pt of this.points) {
      if (pt.tone === 'div' && this.options.shape === 'triangle') {
        svg += `  <polygon points="${pt.x},${pt.y - 6} ${pt.x - 5},${pt.y + 4} ${pt.x + 5},${pt.y + 4}" fill="${PALETTE.accent}" stroke="${PALETTE.accentDeep}"/>\n`;
      } else if (pt.tone === 'pur' && this.options.shape === 'triangle') {
        svg += `  <polygon points="${pt.x},${pt.y + 6} ${pt.x - 5},${pt.y - 4} ${pt.x + 5},${pt.y - 4}" fill="${PALETTE.purifying}" stroke="${PALETTE.purifyingDeep}"/>\n`;
      } else if (pt.tone) {
        svg += `  <circle cx="${pt.x}" cy="${pt.y}" r="4" fill="${pt.tone === 'div' ? PALETTE.accent : PALETTE.purifying}"/>\n`;
      } else {
        svg += `  <circle cx="${pt.x}" cy="${pt.y}" r="2.5" fill="${PALETTE.slateLight}"/>\n`;
      }
    }
    const text = this.ruleLabel(track.threshold);
    for (const [y, color, side] of this.endRuleLabels(frame, track)) {
      const x = frame.left + frame.width - 6;
      svg += `  <text x="${x}" y="${side === 'above' ? y - 5 : y + 13}" fill="${color}" font-weight="bold" text-anchor="end" paint-order="stroke" stroke="${PALETTE.panel}" stroke-width="3">${text}</text>\n`;
    }
    return svg;
  }
}
