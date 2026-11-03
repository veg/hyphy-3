// Segment band: alignment segments as adjacent blocks stacked under a site track (GARD's
// best-model partitions), in alternating light and mid grey, each labelled with its
// range when it fits. A breakpoint sits between consecutive segments, half a site past
// the left segment's end, and gets a solid accent rule through the band and, with
// rules: true, through the main band too. breakTone(bp, track) recolours a break's rule and
// label ('accent' | 'ink' | 'muted'; GARD's incongruent, branch-length-only and untested
// breaks), and breakLabel(bp, track) rewrites its label; bp is { index, at, left, right }. Each block is a hit region of kind 'segment' (key = segment.key).
//
// A segment is { key, start, end, label? } in the track's site coordinates (1-based,
// inclusive), in alignment order.

import { PALETTE } from '../data/palette.js';
import { escapeSvgText } from '../core/export.js';

// Breakpoints between consecutive segments: [{ index (1-based), at (site coordinate) }]
export function segmentBreaks(segments) {
  return segments.slice(0, -1).map((s, i) => ({ index: i + 1, at: s.end + 0.5 }));
}

function segmentTitle(segment, index) {
  return segment.label ?? `Segment ${index + 1} · ${segment.start}–${segment.end}`;
}

// The colours a break tone names, from the track's themed palette when it has one
const SEGMENT_TONES = {
  accent: (colors) => PALETTE.accent,
  ink: (colors) => colors.ink,
  muted: (colors) => colors.stone
};

export class SegmentBand {
  constructor(options = {}) {
    this.options = Object.assign({
      segments: () => [],    // (data, track) => segments
      height: 28,
      rules: true,           // carry the breakpoint rules up through the main band
      breakTone: null,       // (bp, track) => 'accent' | 'ink' | 'muted'; accent when unset
      breakLabel: null,      // (bp, track) => label text; 'BP n' when unset
      emptyLabel: 'No recombination detected — single tree fits best',
      tooltip: null          // (segment, track) => HTML
    }, options);
    this.height = this.options.height;
  }

  // Blocks clipped to the frame: [{ x0, x1, segment, index, title }]
  blocks(frame, track) {
    const segments = this.options.segments(track.data, track);
    const right = frame.left + frame.width;
    const single = segments.length === 1;
    return segments.map((segment, index) => ({
      x0: Math.max(frame.left, frame.scaleX(segment.start - 0.5)),
      x1: Math.min(right, frame.scaleX(segment.end + 0.5)),
      segment,
      index,
      title: single ? this.options.emptyLabel : segmentTitle(segment, index)
    }));
  }

  // Breaks between the blocks' segments, with their tone and label
  breaks(blocks, track, colors) {
    const segments = blocks.map(b => b.segment);
    return segmentBreaks(segments).map((bp, i) => {
      const full = Object.assign(bp, { left: segments[i], right: segments[i + 1] });
      const tone = this.options.breakTone?.(full, track) ?? 'accent';
      return Object.assign(full, {
        color: SEGMENT_TONES[tone](colors),
        label: this.options.breakLabel?.(full, track) ?? `BP ${full.index}`
      });
    });
  }

  isHit(region, block) {
    return region?.kind === 'segment' && region.key === block.segment.key;
  }

  draw(ctx, frame, track) {
    const blocks = this.blocks(frame, track);
    const bottom = frame.top + frame.height;
    const main = track.frames[0];
    ctx.save();
    ctx.font = '10px "SF Mono", Monaco, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    for (const b of blocks) {
      const selected = this.isHit(track.selectedRegion, b);
      const hovered = this.isHit(track.hoveredRegion, b);
      ctx.fillStyle = b.index % 2 ? PALETTE.overviewEmpty : PALETTE.inset;
      ctx.fillRect(b.x0, frame.top, b.x1 - b.x0, frame.height);
      if (selected || hovered) {
        ctx.strokeStyle = PALETTE.ink;
        ctx.lineWidth = selected ? 2 : 1;
        ctx.strokeRect(b.x0 + 1, frame.top + 1, b.x1 - b.x0 - 2, frame.height - 2);
      }
      const label = ctx.measureText(b.title).width <= b.x1 - b.x0 - 8 ? b.title : String(b.index + 1);
      ctx.fillStyle = selected || hovered ? PALETTE.ink : PALETTE.ancestorInk;
      ctx.fillText(label, (b.x0 + b.x1) / 2, frame.top + frame.height / 2);
    }

    ctx.lineWidth = 1.5;
    ctx.font = 'bold 9px "SF Mono", Monaco, monospace';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    const top = this.options.rules ? main.top : frame.top;
    for (const bp of this.breaks(blocks, track, track.colors ?? PALETTE)) {
      const x = frame.scaleX(bp.at);
      ctx.strokeStyle = bp.color;
      ctx.fillStyle = bp.color;
      ctx.beginPath();
      ctx.moveTo(x, top);
      ctx.lineTo(x, bottom);
      ctx.stroke();
      if (this.options.rules) ctx.fillText(bp.label, x + 3, main.top + 2);
    }
    ctx.restore();

    return blocks.map(b => ({
      x0: b.x0, x1: b.x1, y0: frame.top, y1: bottom, key: b.segment.key, kind: 'segment', segment: b.segment, index: b.index
    }));
  }

  regionTooltip(region, track) {
    if (this.options.tooltip) return this.options.tooltip(region.segment, track);
    const { start, end } = region.segment;
    return `
      <div class="tip-header"><span class="tip-title">${escapeSvgText(segmentTitle(region.segment, region.index))}</span></div>
      <div class="tip-grid">
        <div class="tip-row"><span class="tip-lbl">Sites:</span> <span class="tip-val font-mono">${end - start + 1}</span></div>
      </div>
      <div class="tip-hint">Click to inspect this segment</div>
    `;
  }

  svg(frame, track) {
    const blocks = this.blocks(frame, track);
    const bottom = frame.top + frame.height;
    const top = this.options.rules ? track.frames[0].top : frame.top;
    let svg = '';
    for (const b of blocks) {
      svg += `  <rect x="${b.x0}" y="${frame.top}" width="${b.x1 - b.x0}" height="${frame.height}" fill="${b.index % 2 ? PALETTE.overviewEmpty : PALETTE.inset}"/>\n`;
      svg += `  <text x="${(b.x0 + b.x1) / 2}" y="${frame.top + frame.height / 2}" text-anchor="middle" dominant-baseline="middle">${escapeSvgText(b.title)}</text>\n`;
    }
    for (const bp of this.breaks(blocks, track, PALETTE)) {
      const x = frame.scaleX(bp.at);
      svg += `  <line x1="${x}" y1="${top}" x2="${x}" y2="${bottom}" stroke="${bp.color}" stroke-width="1.5"/>\n`;
      if (this.options.rules) svg += `  <text x="${x + 3}" y="${top + 10}" fill="${bp.color}" font-weight="bold">${escapeSvgText(bp.label)}</text>\n`;
    }
    return svg;
  }
}
