// Matrix band: a rows x sites diverging matrix stacked under a site track (PRIME's
// property x site lambda, Contrast-FEL's group x site relative beta, aBSREL's branch x site
// EBF). Each cell is a signed magnitude: positive in one hue, negative in its mirror,
// opacity scaled by |value| over the band's extent and by the cell's emphasis tier.
// When sites are narrower than a pixel, each pixel column shows its strongest member.
// Row labels sit in the track's left gutter; clicking one focuses that row (clicking it
// again clears the focus).
//
// A matrix row is { key, label, title?, muted? }; a muted row (Contrast-FEL's background
// set) is drawn at reduced saturation. Cell hit regions are kind 'site' with the row
// attached; label hit regions are kind 'row'.

import { PALETTE, withAlpha } from '../data/palette.js';
import { escapeSvgText } from '../core/export.js';

// [max |value|] over every non-blank cell; 1 when there is none
export function matrixExtent(sites, rows, value, track) {
  let extent = 0;
  for (const site of sites) {
    for (const row of rows) {
      const v = value(site, row, track);
      if (v != null && Math.abs(v) > extent) extent = Math.abs(v);
    }
  }
  return extent || 1;
}

// Cell columns across the frame: one per site, spanning half a site either side, or one
// per pixel column once sites are narrower than a pixel. Each holds the site rows it covers.
export function matrixColumns(sites, frame, numSites) {
  const right = frame.left + frame.width;
  const step = frame.width / Math.max(1, numSites - 1);
  if (step >= 1) {
    return sites.map(site => {
      const x = frame.scaleX(site.site);
      return { x0: Math.max(frame.left, x - step / 2), x1: Math.min(right, x + step / 2), members: [site] };
    });
  }
  const byPixel = new Map();
  for (const site of sites) {
    const px = Math.min(Math.floor(frame.scaleX(site.site)), Math.ceil(right) - 1);
    if (!byPixel.has(px)) byPixel.set(px, { x0: px, x1: px + 1, members: [] });
    byPixel.get(px).members.push(site);
  }
  return [...byPixel.values()];
}

export class MatrixBand {
  constructor(options = {}) {
    this.options = Object.assign({
      rows: [],                    // matrix rows, or (data, track) => matrix rows
      value: () => null,           // (site, row, track) => signed magnitude; null leaves the cell blank
      emphasis: () => 0,           // (site, row, track) => tier 0 | 1 | 2 at the live threshold
      extent: 'auto',              // 'auto' (max |value| in the file) | number | track => number
      positive: PALETTE.accent,    // hue for value > 0
      negative: PALETTE.purifying, // hue for value < 0
      tierOpacity: [0.15, 0.5, 1], // opacity by tier; tier 2 also gets an outline
      mutedOpacity: 0.4,
      rowHeight: 14,
      tooltip: null,               // (site, row, track) => HTML for a cell
      rowTooltip: null,            // (row, track) => HTML for a row label
      onFocus: null                // (row key | null) when a label click changes the focus
    }, options);
    this.focused = null;
    this.rowList = Array.isArray(this.options.rows) ? this.options.rows : [];
  }

  get height() {
    return this.rowList.length * this.options.rowHeight;
  }

  measure(track) {
    const { rows } = this.options;
    this.rowList = typeof rows === 'function' ? rows(track.data, track) : rows;
  }

  focus(key) {
    this.focused = key;
  }

  select(region) {
    if (region.kind !== 'row') return;
    this.focus(region.key === this.focused ? null : region.key);
    this.options.onFocus?.(this.focused);
  }

  extent(track) {
    const { extent, value } = this.options;
    if (typeof extent === 'number') return extent;
    if (typeof extent === 'function') return extent(track);
    return matrixExtent(track.rows, this.rowList, value, track);
  }

  // Every non-blank cell, after binning: { x0, x1, y0, y1, site, row, value, tier, fill, alpha }
  cells(frame, track) {
    const { value, emphasis, positive, negative, tierOpacity, mutedOpacity, rowHeight } = this.options;
    const extent = this.extent(track);
    const columns = matrixColumns(track.rows, frame, track.numSites);
    const cells = [];
    this.rowList.forEach((row, k) => {
      const y0 = frame.top + k * rowHeight;
      for (const column of columns) {
        let best = null;
        for (const site of column.members) {
          const v = value(site, row, track);
          if (v == null) continue;
          const tier = emphasis(site, row, track);
          if (!best || tier > best.tier || (tier === best.tier && Math.abs(v) > Math.abs(best.value))) {
            best = { site, value: v, tier };
          }
        }
        if (!best) continue;
        const magnitude = Math.min(1, Math.abs(best.value) / extent);
        cells.push({
          x0: column.x0, x1: column.x1, y0, y1: y0 + rowHeight,
          site: best.site, row, value: best.value, tier: best.tier,
          fill: best.value >= 0 ? positive : negative,
          alpha: magnitude * tierOpacity[best.tier] * (row.muted ? mutedOpacity : 1)
        });
      }
    });
    return cells;
  }

  labelColor(row, track) {
    const hovered = track.hoveredRegion?.kind === 'row' && track.hoveredRegion.key === row.key;
    return row.key === this.focused || hovered ? PALETTE.ink : PALETTE.stone;
  }

  draw(ctx, frame, track) {
    const { rowHeight } = this.options;
    const right = frame.left + frame.width;
    const cells = this.cells(frame, track);
    ctx.save();
    ctx.fillStyle = PALETTE.cellBg;
    ctx.fillRect(frame.left, frame.top, frame.width, frame.height);

    for (const c of cells) {
      ctx.fillStyle = withAlpha(c.fill, c.alpha);
      ctx.fillRect(c.x0, c.y0, c.x1 - c.x0, rowHeight);
    }
    ctx.strokeStyle = PALETTE.ink;
    ctx.lineWidth = 1;
    for (const c of cells) {
      if (c.tier === 2) ctx.strokeRect(c.x0 + 0.5, c.y0 + 0.5, Math.max(0, c.x1 - c.x0 - 1), rowHeight - 1);
    }

    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    const regions = [];
    this.rowList.forEach((row, k) => {
      const y0 = frame.top + k * rowHeight;
      const focused = row.key === this.focused;
      ctx.font = `${focused ? 'bold ' : ''}9px "SF Mono", Monaco, monospace`;
      ctx.fillStyle = this.labelColor(row, track);
      ctx.fillText(row.label, frame.left - 8, y0 + rowHeight / 2);
      if (focused) {
        ctx.strokeStyle = PALETTE.ink;
        ctx.strokeRect(frame.left - 0.5, y0 + 0.5, frame.width + 1, rowHeight - 1);
      }
      regions.push({ x0: 0, x1: frame.left - 2, y0, y1: y0 + rowHeight, key: row.key, kind: 'row', row });
    });

    ctx.strokeStyle = PALETTE.hairline;
    ctx.strokeRect(frame.left - 0.5, frame.top - 0.5, right - frame.left + 1, frame.height + 1);
    ctx.restore();

    for (const c of cells) {
      regions.push({ x0: c.x0, x1: c.x1, y0: c.y0, y1: c.y1, key: c.site.site, kind: 'site', raw: c.site, row: c.row });
    }
    return regions;
  }

  regionTooltip(region, track) {
    if (region.kind === 'row') return this.options.rowTooltip?.(region.row, track) ?? '';
    return this.options.tooltip?.(region.raw, region.row, track) ?? '';
  }

  svg(frame, track) {
    const { rowHeight } = this.options;
    let svg = `  <rect x="${frame.left}" y="${frame.top}" width="${frame.width}" height="${frame.height}" fill="${PALETTE.cellBg}" stroke="${PALETTE.hairline}"/>\n`;
    for (const c of this.cells(frame, track)) {
      const outline = c.tier === 2 ? ` stroke="${PALETTE.ink}"` : '';
      svg += `  <rect x="${c.x0}" y="${c.y0}" width="${c.x1 - c.x0}" height="${rowHeight}" fill="${c.fill}" fill-opacity="${c.alpha}"${outline}/>\n`;
    }
    this.rowList.forEach((row, k) => {
      const weight = row.key === this.focused ? ' font-weight="bold"' : '';
      svg += `  <text x="${frame.left - 8}" y="${frame.top + (k + 0.5) * rowHeight}" text-anchor="end" dominant-baseline="middle"${weight}>${escapeSvgText(row.label)}</text>\n`;
    });
    return svg;
  }
}
