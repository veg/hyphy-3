// Tick strip: a short greyscale band stacked under a site track (dFEL's GC/AT direction
// strip). Each site gets an upward dark tick and a downward light tick sized by two
// counts, and an optional black dot. Orange and blue stay reserved for selection.

import { PALETTE } from '../data/palette.js';

// Ceiling for tick heights: the given quantile of the nonzero counts, at least 1
export function tickStripCeiling(counts, quantile) {
  const nonzero = counts.filter(c => c > 0).sort((a, b) => a - b);
  if (nonzero.length === 0) return 1;
  return Math.max(1, nonzero[Math.ceil(quantile * nonzero.length) - 1]);
}

export class TickStripBand {
  constructor(options = {}) {
    this.options = Object.assign({
      height: 20,
      up: () => 0,          // row => count drawn as an upward dark tick
      down: () => 0,        // row => count drawn as a downward light tick
      flag: () => false,    // row => true draws a 3px black dot on the baseline
      clamp: 0.95,          // ticks saturate at this quantile of the nonzero counts
      label: ''             // short label in the left gutter
    }, options);
    this.height = this.options.height;
  }

  // One mark per site with something to draw; tick lengths in pixels
  marks(frame, track) {
    const { up, down, flag, clamp } = this.options;
    const rows = track.rows;
    const counts = rows.flatMap(r => [up(r), down(r)]);
    const ceiling = tickStripCeiling(counts, clamp);
    const half = frame.height / 2;
    const length = (count) => Math.min(1, count / ceiling) * half;

    return rows.map(r => ({
      site: r.site,
      x: frame.scaleX(r.site),
      up: length(up(r)),
      down: length(down(r)),
      flag: Boolean(flag(r))
    })).filter(m => m.up > 0 || m.down > 0 || m.flag);
  }

  // Drawn in the track's themed colours, so the strip follows the dark theme
  draw(ctx, frame, track) {
    const colors = track.colors ?? PALETTE;
    const mid = frame.top + frame.height / 2;
    ctx.save();
    ctx.strokeStyle = colors.gridLine;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(frame.left, mid);
    ctx.lineTo(frame.left + frame.width, mid);
    ctx.stroke();

    if (this.options.label) {
      ctx.fillStyle = colors.stone;
      ctx.font = '9px "SF Mono", Monaco, monospace';
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      ctx.fillText(this.options.label, frame.left - 8, mid);
    }

    ctx.lineWidth = 1.5;
    for (const m of this.marks(frame, track)) {
      if (m.up > 0) {
        ctx.strokeStyle = PALETTE.slate;
        ctx.beginPath();
        ctx.moveTo(m.x, mid);
        ctx.lineTo(m.x, mid - m.up);
        ctx.stroke();
      }
      if (m.down > 0) {
        ctx.strokeStyle = PALETTE.slateLight;
        ctx.beginPath();
        ctx.moveTo(m.x, mid);
        ctx.lineTo(m.x, mid + m.down);
        ctx.stroke();
      }
      if (m.flag) {
        ctx.fillStyle = colors.ink;
        ctx.beginPath();
        ctx.arc(m.x, mid, 1.5, 0, Math.PI * 2);
        ctx.fill();
      }
    }
    ctx.restore();
  }

  svg(frame, track) {
    const mid = frame.top + frame.height / 2;
    let svg = `  <line x1="${frame.left}" y1="${mid}" x2="${frame.left + frame.width}" y2="${mid}" stroke="${PALETTE.gridLine}"/>\n`;
    if (this.options.label) {
      svg += `  <text x="${frame.left - 8}" y="${mid}" text-anchor="end" dominant-baseline="middle">${this.options.label}</text>\n`;
    }
    for (const m of this.marks(frame, track)) {
      if (m.up > 0) svg += `  <line x1="${m.x}" y1="${mid}" x2="${m.x}" y2="${mid - m.up}" stroke="${PALETTE.slate}" stroke-width="1.5"/>\n`;
      if (m.down > 0) svg += `  <line x1="${m.x}" y1="${mid}" x2="${m.x}" y2="${mid + m.down}" stroke="${PALETTE.slateLight}" stroke-width="1.5"/>\n`;
      if (m.flag) svg += `  <circle cx="${m.x}" cy="${mid}" r="1.5" fill="${PALETTE.ink}"/>\n`;
    }
    return svg;
  }
}
