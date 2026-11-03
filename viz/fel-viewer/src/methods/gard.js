// GARD method config: recombination breakpoints found by a genetic-algorithm search over
// alignment partitions, each partition fitted its own tree. The hero is the support each
// site has for a breakpoint after it (the Akaike weight over the models the search tried)
// with the best model's segments under it; the table has one row per breakpoint, called at
// a family-wise alpha. A breakpoint is incongruent (recombination) when the KH or SH test
// of each flanking segment's tree against its neighbour's rejects at alpha / 2B, with B the
// number of breakpoints; a tree that differs only in its branch lengths is rate
// heterogeneity, not a topology change. The Ohmic row marks candidate breakpoints from
// the spectral tree engine's seeded search.

import { GardResults, GARD_RUN_ALPHA, gardCutoff, gardFormat, gardP } from '../data/gard-data.js';
import { KpiBar } from '../components/kpi-bar.js';
import { SiteTrackCanvas } from '../components/site-track.js';
import { EvidenceBand } from '../components/evidence-band.js';
import { OverviewTrack } from '../components/overview-track.js';
import { SegmentBand } from '../components/segment-band.js';
import { evidenceTipRow } from '../components/evidence-plot.js';
import { ResultsTable } from '../components/results-table.js';
import { Drawer } from '../components/drawer.js';
import { GardTreePair, GardTreeCard, gardSegmentLabel } from '../components/gard-trees.js';
import { PALETTE } from '../data/palette.js';
import { formatThreshold } from '../core/format.js';
import { escapeSvgText } from '../core/export.js';

// The family-wise alpha for the topology tests, on a log slider
export const GARD_ALPHA = {
  title: 'Family-wise α (KH/SH)',
  statistic: 'α',
  direction: 'le',
  min: 0.001,
  max: 0.2,
  step: 0.001,
  scale: 'log',
  presets: [0.01, 0.05, 0.1],
  seed: GARD_RUN_ALPHA,
  format: formatThreshold,
  caption: ({ threshold, data }) => (data && data.breaks.length
    ? `Bonferroni per test α/2B = ${gardCutoff(threshold / (2 * data.breaks.length))}, B = ${data.breaks.length}`
    : 'Bonferroni per test α/2B'),
  enabled: (data) => data.capabilities.topology,
  disabledCaption: 'No topology test'
};

const GARD_CALLS = {
  incongruent: { tone: 'div', glyph: '✱', label: 'Incongruent', title: 'Topological incongruence (recombination)', breakTone: 'accent' },
  branch: { tone: 'neutral', glyph: '≈', label: 'Branch lengths only', title: 'Rate heterogeneity: the topology is not rejected', breakTone: 'ink' },
  untested: { tone: 'neutral', glyph: '–', label: 'Untested', title: 'No topology test for this breakpoint', breakTone: 'muted' }
};

export function gardCall(data, brk, alpha) {
  return GARD_CALLS[data.callOf(brk, alpha)];
}

function gardPlural(n, noun) {
  return `${n} ${noun}${n === 1 ? '' : 's'}`;
}

// KPI cards

function gardAlignmentSub(data) {
  const { variableSites, potentialBreakpoints } = data.dataset;
  const parts = [variableSites != null && `${variableSites} variable`, potentialBreakpoints != null && `${potentialBreakpoints} candidate breakpoints`];
  return parts.filter(Boolean).join(' · ') || '—';
}

export function gardTopologyValue(stats, { data }) {
  if (!data.breaks.length) return 'Not applicable';
  if (!data.capabilities.topology) {
    const supported = data.topology?.supported;
    return supported == null ? '–' : (supported ? 'Supported' : 'Not supported');
  }
  return `${stats.incongruent} of ${stats.breakpoints} incongruent`;
}

export function gardTopologySub(stats, { data }) {
  if (!data.breaks.length) return 'No breakpoints to test';
  const delta = data.topology?.deltaAicc;
  const text = delta == null ? '' : ` · ΔAICc_topology ${gardFormat(delta, 1)}`;
  if (!data.capabilities.topology) return `Topology tests not in this file${text}`;
  if (data.topology?.supported === false) return `Branch-length heterogeneity only${text}`;
  return `ΔAICc_topology ${gardFormat(delta, 1)} (topology, not just branch lengths)`;
}

export const GARD_KPI_CARDS = [
  {
    id: 'alignment',
    valueId: 'sites',
    subId: 'taxa',
    cardClass: 'kpi-card-text',
    label: 'Alignment',
    placeholder: '— Taxa',
    value: (stats, { data }) => `${data.numTaxa} taxa · ${data.numSites.toLocaleString()} sites`,
    sub: (stats, { data }) => gardAlignmentSub(data),
    tone: 'ink'
  },
  {
    id: 'breaks',
    cardClass: 'kpi-card-text',
    label: 'Breakpoints',
    value: (stats) => (stats.breakpoints ? gardPlural(stats.breakpoints, 'breakpoint') : 'No recombination detected'),
    sub: (stats, { data }) => (stats.breakpoints
      ? `ΔAICc ${gardFormat(data.deltaAicc, 1)} vs single tree`
      : 'A single tree fits best (ΔAICc 0)'),
    tone: (stats) => (stats.breakpoints ? 'ink' : 'muted')
  },
  {
    id: 'topology',
    cardClass: 'kpi-card-text',
    label: 'Topology · incongruent at α',
    value: gardTopologyValue,
    sub: gardTopologySub,
    tone: (stats, { data }) => {
      if (!data.capabilities.topology) return 'muted';
      return stats.incongruent > 0 ? 'accent' : 'ink';
    }
  },
  {
    id: 'ohmic',
    cardClass: 'kpi-card-text',
    label: 'Ohmic network',
    available: (data) => data.capabilities.ohmic,
    value: (stats, { data }) => gardPlural(data.ohmic.count, 'bridge'),
    sub: (stats, { data }) => `${gardPlural(data.ohmic.candidates.length, 'candidate breakpoint')} · ΔAICc vs baseline ${gardFormat(data.ohmic.deltaAicc, 1)}`,
    tone: (stats, { data }) => (data.ohmic.count ? 'ink' : 'muted')
  }
];

// Hero: breakpoint support per site, the best model's segments and the Ohmic candidates

// Support on its log axis, held at 1e-6 (the writer's values underflow to 1e-320)
export const GARD_SUPPORT_FLOOR = -6;

export function gardLogSupport(p) {
  return Math.max(GARD_SUPPORT_FLOOR, Math.log10(p));
}

// Ticks at whole decades: '≤1e-6' at the floor, 0.1, 1e-2 ... and 1
export function gardSupportTick(v) {
  if (v <= GARD_SUPPORT_FLOOR) return `≤1e-${-GARD_SUPPORT_FLOOR}`;
  if (v === 0) return '1';
  return v === -1 ? '0.1' : `1e${v}`;
}

// The support as a number to read: three decimals, exponential below 0.001
export function gardSupportText(p) {
  return p < 0.001 ? p.toExponential(2) : p.toFixed(3);
}

const gardSupportTone = (r) => (r.p >= 0.5 ? 'div' : null);

function gardSupportTooltip(r, track) {
  const brk = track.data.getBreak(r.site);
  const call = brk ? gardCall(track.data, brk, track.threshold) : null;
  const badge = brk ? `<span class="tip-badge tip-${call.tone}">BP ${brk.index} · ${call.label}</span>` : '';
  return `
      <div class="tip-header">
        <span class="tip-title">Breakpoint after site ${r.site}</span>${badge}
      </div>
      <div class="tip-grid">${evidenceTipRow('Support (Akaike weight)', gardSupportText(r.p))}</div>
      ${brk ? '<div class="tip-hint">Click to inspect the breakpoint</div>' : ''}
    `;
}

function gardSupportBand(options) {
  return new EvidenceBand(Object.assign({
    tone: gardSupportTone,
    shape: 'dot',
    rule: null,
    neutralAlpha: 0.7,
    tooltip: gardSupportTooltip
  }, options));
}

// Fresh bands per plot, since a band caches the points it last drew
export function gardTrackModes() {
  return [
    {
      id: 'log',
      label: 'Log',
      subtitle: () => 'Support for a breakpoint after each site (Akaike weight, log scale)',
      band: gardSupportBand({
        value: (r) => gardLogSupport(r.p),
        domain: [GARD_SUPPORT_FLOOR, 0],
        ticks: [-6, -5, -4, -3, -2, -1, 0],
        tickFormat: gardSupportTick,
        label: 'support'
      })
    },
    {
      id: 'linear',
      label: 'Linear',
      subtitle: () => 'Support for a breakpoint after each site (Akaike weight)',
      band: gardSupportBand({
        value: (r) => r.p,
        domain: [0, 1],
        ticks: [0, 0.25, 0.5, 0.75, 1],
        tickFormat: (v) => String(v),
        label: 'support'
      })
    }
  ];
}

// The overview brush: support on the same log axis, a site called supported from 0.5
export const GARD_BRUSH = {
  minSites: 1000,
  maxValue: -GARD_SUPPORT_FLOOR,
  unit: 'Sites',
  toneLabels: { accent: 'with support ≥ 0.5' },
  value: (r) => gardLogSupport(r.p) - GARD_SUPPORT_FLOOR,
  tone: (r) => (r.p >= 0.5 ? 'accent' : null),
  rule: () => null
};

// A zoomed window of a tenth of the alignment, at least 300 sites
export function gardWindow(numSites) {
  return Math.max(300, Math.ceil(numSites / 10));
}

// A break's rule and label in the segment band, by its call at the live alpha
export function gardBreakTone(bp, track) {
  return gardCall(track.data, track.data.breaks[bp.index - 1], track.threshold).breakTone;
}

export function gardBreakLabel(bp, track) {
  const call = gardCall(track.data, track.data.breaks[bp.index - 1], track.threshold);
  return call === GARD_CALLS.incongruent ? `BP ${bp.index} ✱` : `BP ${bp.index}`;
}

function gardSegmentTooltip(segment, track) {
  const { data } = track;
  const index = data.segments.indexOf(segment);
  const rows = [
    evidenceTipRow('Sites', `${segment.start}–${segment.end} (${segment.end - segment.start + 1})`),
    ...(segment.logL == null ? [] : [evidenceTipRow('log L', gardFormat(segment.logL, 1))])
  ];
  const hint = data.breaks.length ? '<div class="tip-hint">Click to inspect its breakpoint</div>' : '';
  return `
      <div class="tip-header"><span class="tip-title">${escapeSvgText(gardSegmentLabel(segment, index))}</span></div>
      <div class="tip-grid">${rows.join('')}</div>${hint}
    `;
}

// The breakpoint a clicked segment opens: the one on its right, else (the last segment) on its left
export function gardSegmentBreak(data, index) {
  return data.breaks[index] ?? data.breaks[index - 1] ?? null;
}

// A strip of hollow triangles at the Ohmic candidate breakpoints, under the segment band;
// it takes no room for a result without the Ohmic block
export class GardTickBand {
  constructor(options = {}) {
    this.options = Object.assign({ size: 5, hit: 5, label: 'Ohmic' }, options);
    this.height = 0;
  }

  measure(track) {
    this.height = track.data?.capabilities.ohmic ? 20 : 0;
  }

  // The candidates in view: [{ site, x }]
  marks(frame, track) {
    const [from, to] = track.siteRange();
    return track.data.ohmic.candidates.filter(s => s >= from && s <= to).map(site => ({ site, x: frame.scaleX(site) }));
  }

  triangle(x, mid) {
    const { size } = this.options;
    return [[x, mid - size], [x - size, mid + size * 0.8], [x + size, mid + size * 0.8]];
  }

  draw(ctx, frame, track) {
    if (!this.height) return null;
    const colors = track.colors ?? PALETTE;
    const mid = frame.top + frame.height / 2;
    const marks = this.marks(frame, track);
    ctx.save();
    ctx.strokeStyle = colors.gridLine;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(frame.left, mid);
    ctx.lineTo(frame.left + frame.width, mid);
    ctx.stroke();
    ctx.fillStyle = colors.stone;
    ctx.font = '9px "SF Mono", Monaco, monospace';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(this.options.label, frame.left - 8, mid);
    ctx.fillStyle = colors.panel;
    ctx.strokeStyle = PALETTE.slate;
    ctx.lineWidth = 1.25;
    for (const m of marks) {
      const [first, ...rest] = this.triangle(m.x, mid);
      ctx.beginPath();
      ctx.moveTo(...first);
      rest.forEach(([px, py]) => ctx.lineTo(px, py));
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    }
    ctx.restore();
    const { hit } = this.options;
    return marks.map(m => ({ x0: m.x - hit, x1: m.x + hit, y0: frame.top, y1: frame.top + frame.height, key: m.site, kind: 'candidate' }));
  }

  regionTooltip(region) {
    return `
      <div class="tip-header"><span class="tip-title">After site ${region.key}</span></div>
      <div class="tip-grid">${evidenceTipRow('Ohmic bridge candidate', 'seeded GA')}</div>
    `;
  }

  svg(frame, track) {
    if (!this.height) return '';
    const mid = frame.top + frame.height / 2;
    let svg = `  <line x1="${frame.left}" y1="${mid}" x2="${frame.left + frame.width}" y2="${mid}" stroke="${PALETTE.gridLine}"/>\n`;
    svg += `  <text x="${frame.left - 8}" y="${mid}" text-anchor="end" dominant-baseline="middle">${this.options.label}</text>\n`;
    for (const m of this.marks(frame, track)) {
      svg += `  <polygon points="${this.triangle(m.x, mid).map(p => p.join(',')).join(' ')}" fill="${PALETTE.panel}" stroke="${PALETTE.slate}" stroke-width="1.25"/>\n`;
    }
    return svg;
  }
}

export function gardHero(parent, hooks) {
  return new OverviewTrack(parent, {
    track: (container) => {
      const track = new SiteTrackCanvas(container, {
        title: 'Breakpoint Support',
        subtitle: 'Support for a breakpoint after each site (Akaike weight, log scale)',
        modes: gardTrackModes(),
        stack: [
          new SegmentBand({
            segments: (data) => data.segments,
            breakTone: gardBreakTone,
            breakLabel: gardBreakLabel,
            tooltip: gardSegmentTooltip
          }),
          new GardTickBand()
        ],
        sites: (data) => data.support,
        threshold: hooks.threshold,
        exportName: 'hyphy_gard_support_plot',
        onSiteSelect: (site) => hooks.onSiteSelect(site),
        onRegionSelect: (region) => {
          const brk = region.kind === 'segment' ? gardSegmentBreak(track.data, region.index) : null;
          if (brk) hooks.onSiteSelect(brk.site);
        }
      });
      return track;
    },
    brush: GARD_BRUSH,
    window: gardWindow,
    onSiteSelect: (site) => hooks.onSiteSelect(site)
  });
}

// Results table: one row per breakpoint

function gardCallBadge(data, r, alpha) {
  const call = gardCall(data, r, alpha);
  const key = data.callOf(r, alpha);
  return `<span class="omega-badge gard-call gard-call-${key}" title="${call.title}">${call.glyph} ${call.label}</span>`;
}

const gardInterval = (segment) => `${segment.start}–${segment.end}`;

// 'a / b' for a flank pair, dashes where the break has no test
function gardPair(test, left, right, format) {
  return test ? `${format(test[left])} / ${format(test[right])}` : '–';
}

// A KH or SH p pair: bold in the accent when both flanks reject at alpha / 2B
function gardPCell(r, { data, threshold }, left, right) {
  const cutoff = threshold / (2 * data.breaks.length);
  const rejects = r.test && r.test[left] != null && r.test[right] != null && r.test[left] <= cutoff && r.test[right] <= cutoff;
  return `<td class="font-mono${rejects ? ' font-bold text-orange' : ''}">${gardPair(r.test, left, right, gardP)}</td>`;
}

const gardCsv = (v) => v ?? '';

export const GARD_COLUMNS = [
  {
    id: 'bp',
    header: 'BP',
    width: 50,
    alwaysShown: true,
    cell: (r) => `<td class="font-mono">${r.index}</td>`,
    csv: [{ header: 'Breakpoint', value: (r) => r.index }],
    sort: (r) => r.index
  },
  {
    id: 'site',
    header: 'After site',
    width: 80,
    alwaysShown: true,
    title: 'The breakpoint falls between this site and the next',
    cell: (r) => `<td class="font-mono">${r.site}</td>`,
    csv: [{ header: 'After site', value: (r) => r.site }],
    sort: (r) => r.site
  },
  {
    id: 'left',
    header: 'Left segment',
    width: 100,
    cell: (r) => `<td class="font-mono">${gardInterval(r.left)}</td>`,
    csv: [{ header: 'Left start', value: (r) => r.left.start }, { header: 'Left end', value: (r) => r.left.end }],
    sort: (r) => r.left.start
  },
  {
    id: 'right',
    header: 'Right segment',
    width: 100,
    cell: (r) => `<td class="font-mono">${gardInterval(r.right)}</td>`,
    csv: [{ header: 'Right start', value: (r) => r.right.start }, { header: 'Right end', value: (r) => r.right.end }],
    sort: (r) => r.right.start
  },
  {
    id: 'call',
    header: 'Call',
    width: 150,
    alwaysShown: true,
    title: 'The call at the current α: both flanks must reject at α / 2B',
    cell: (r, { data, threshold }) => `<td>${gardCallBadge(data, r, threshold)}</td>`,
    csv: [{ header: 'Call', value: (r, { data, threshold }) => gardCall(data, r, threshold).label }],
    sort: (r) => Math.max(r.test?.khLeft ?? 1, r.test?.khRight ?? 1)
  },
  {
    id: 'dlogl',
    header: 'ΔlogL L / R',
    width: 110,
    availableWhen: (caps) => caps.topology,
    title: 'Log-likelihood gain of each flank\'s own tree over its neighbour\'s tree',
    cell: (r) => `<td class="font-mono">${gardPair(r.test, 'dLogLLeft', 'dLogLRight', (v) => gardFormat(v, 1))}</td>`,
    csv: [{ header: 'dlogL left', value: (r) => gardCsv(r.test?.dLogLLeft) }, { header: 'dlogL right', value: (r) => gardCsv(r.test?.dLogLRight) }],
    sort: (r) => Math.min(r.test?.dLogLLeft ?? -Infinity, r.test?.dLogLRight ?? -Infinity)
  },
  {
    id: 'kh',
    header: 'KH p L / R',
    width: 130,
    availableWhen: (caps) => caps.topology,
    cell: (r, view) => gardPCell(r, view, 'khLeft', 'khRight'),
    csv: [{ header: 'KH p left', value: (r) => gardCsv(r.test?.khLeft) }, { header: 'KH p right', value: (r) => gardCsv(r.test?.khRight) }],
    sort: (r) => Math.max(r.test?.khLeft ?? 1, r.test?.khRight ?? 1)
  },
  {
    id: 'sh',
    header: 'SH p L / R',
    width: 130,
    availableWhen: (caps) => caps.topology,
    cell: (r, view) => gardPCell(r, view, 'shLeft', 'shRight'),
    csv: [{ header: 'SH p left', value: (r) => gardCsv(r.test?.shLeft) }, { header: 'SH p right', value: (r) => gardCsv(r.test?.shRight) }],
    sort: (r) => Math.max(r.test?.shLeft ?? 1, r.test?.shRight ?? 1)
  },
  {
    id: 'rf',
    header: 'RF',
    width: 80,
    title: 'Robinson–Foulds distance between the trees either side, of its maximum',
    cell: (r) => `<td class="font-mono">${r.diff ? `${r.diff.rf} / ${r.diff.max}` : '–'}</td>`,
    csv: [{ header: 'RF', value: (r) => gardCsv(r.diff?.rf) }, { header: 'RF maximum', value: (r) => gardCsv(r.diff?.max) }],
    sort: (r) => r.diff?.rf ?? -1
  },
  {
    id: 'moved',
    header: 'Moved taxa',
    width: 90,
    title: 'Taxa whose sister clade differs between the trees either side',
    cell: (r) => `<td class="font-mono">${r.diff?.moved.length ?? '–'}</td>`,
    csv: [{ header: 'Moved taxa', value: (r) => gardCsv(r.diff?.moved.length) }],
    sort: (r) => r.diff?.moved.length ?? -1
  },
  {
    id: 'actions',
    header: 'Inspect',
    width: 80,
    align: 'right',
    cell: (r) => `<td style="text-align: right;"><button type="button" class="btn-inspect" data-action="trees" data-site="${r.site}">Inspect</button></td>`
  }
];

export const GARD_FILTERS = [
  {
    id: 'incongruent',
    key: 'inc',
    glyph: '✱',
    label: 'Incongruent',
    unit: 'incongruent breakpoint',
    count: (stats) => stats.incongruent,
    test: (r, t, { data }) => data.callOf(r, t) === 'incongruent',
    subtitle: (n, t) => `Showing ${n} breakpoints with significant topological incongruence (α ≤ ${formatThreshold(t)})`
  },
  {
    id: 'all',
    key: 'all',
    label: 'All Breakpoints',
    unit: 'breakpoint',
    count: (stats) => stats.breakpoints,
    test: () => true,
    subtitle: (n) => `Showing all ${n} breakpoints`
  }
];

// No breakpoints at all says so; an empty incongruent filter offers the full list
export function gardEmptyState(t, { filter, stats }) {
  if (!stats.breakpoints) {
    return { icon: '○', title: 'No breakpoints — a single tree fits best (ΔAICc 0)' };
  }
  const state = {
    icon: filter.glyph,
    title: `No ${filter.unit} at α ≤ ${formatThreshold(t)}`,
    desc: 'No breakpoint has both flanks rejecting the neighbouring tree at α / 2B.'
  };
  state.jump = { filter: 'all', label: `View All Breakpoints (${stats.breakpoints})` };
  return state;
}

export const GARD_TABLE = {
  title: 'Breakpoints',
  placeholder: 'Called at the family-wise α',
  rowKey: 'site',
  rows: (data) => data.breaks,
  columns: GARD_COLUMNS,
  filters: GARD_FILTERS,
  defaultFilter: 'all',
  defaultSort: { column: 'bp', dir: 'asc' },
  emptyState: gardEmptyState,
  exportName: 'gard-breakpoints'
};

// Drawer: the trees either side of a breakpoint and the topology test behind its call

function gardDrawerHeader(r, { data, threshold }) {
  const call = gardCall(data, r, threshold);
  return { tag: `${call.label} · after site ${r.site}`, tone: call.tone, badge: `BP ${r.index}` };
}

const GARD_TREES_TAB = {
  id: 'trees',
  label: 'Trees',
  comment: 'The trees either side of the breakpoint, moved taxa marked',
  mount: (panel) => {
    panel.innerHTML = `
        <div class="drawer-section">
          <div class="drawer-section-title">Trees either side of the breakpoint</div>
          <div class="gard-drawer-pair"></div>
        </div>`;
    return { pair: new GardTreePair(panel.querySelector('.gard-drawer-pair'), { exportName: 'hyphy_gard_breakpoint' }) };
  },
  render: (panel, r, { data, view }) => {
    view.pair.setSegments(data, r.index - 1, r.index);
  }
};

function gardStatCard(label, value, sub, hue = '') {
  return `
            <div class="drawer-stat-card">
              <div class="dstat-label">${label}</div>
              <div class="dstat-value font-mono font-bold ${hue}">${value}</div>
              <div class="dstat-sub">${sub}</div>
            </div>`;
}

const GARD_MODE_NAMES = { kh: 'KH', sh: 'SH', both: 'KH or SH' };

// The plain-language reading of a break's call
export function gardTestSentence(data, r, alpha) {
  const cutoff = gardCutoff(alpha / (2 * data.breaks.length));
  const mode = GARD_MODE_NAMES[data.testMode];
  return data.callOf(r, alpha) === 'incongruent'
    ? `Each flank's own tree fits significantly better than its neighbour's (${mode} at α / 2B = ${cutoff} on both flanks): the topologies differ, evidence of recombination.`
    : `At least one flank does not reject its neighbour's tree at α / 2B = ${cutoff} (${mode}): the trees may differ in branch lengths, which is rate heterogeneity rather than recombination.`;
}

const GARD_TEST_TAB = {
  id: 'test',
  label: 'Test',
  comment: 'KH and SH tests of each flank against its neighbour\'s tree',
  availableWhen: (data) => data.capabilities.topology,
  unavailable: 'Topology tests are not in this file.',
  render: (panel, r, { data, threshold }) => {
    const t = r.test;
    const cutoff = threshold / (2 * data.breaks.length);
    const passes = (p) => p != null && p <= cutoff;
    const hue = (p) => (passes(p) ? 'text-orange' : '');
    const cards = t ? [
      gardStatCard('KH p, left', gardP(t.khLeft), `left flank ${t.leftInterval?.join('–') ?? gardInterval(r.left)}`, hue(t.khLeft)),
      gardStatCard('KH p, right', gardP(t.khRight), `right flank ${t.rightInterval?.join('–') ?? gardInterval(r.right)}`, hue(t.khRight)),
      gardStatCard('SH p, left', gardP(t.shLeft), 'Shimodaira–Hasegawa', hue(t.shLeft)),
      gardStatCard('SH p, right', gardP(t.shRight), 'Shimodaira–Hasegawa', hue(t.shRight)),
      gardStatCard('ΔlogL, left', gardFormat(t.dLogLLeft, 2), 'own tree over the neighbour\'s'),
      gardStatCard('ΔlogL, right', gardFormat(t.dLogLRight, 2), 'own tree over the neighbour\'s')
    ].join('') : '';
    panel.innerHTML = t ? `<div class="drawer-stack">
        <div class="drawer-section">
          <div class="drawer-section-title">Topology test at α ≤ ${formatThreshold(threshold)}</div>
          <div class="drawer-stats-grid">${cards}
          </div>
          <div class="count-sentence">${gardTestSentence(data, r, threshold)}</div>
          <div class="empty-events-notice">Each flank's tree is tested against its neighbour's on the flank's sites (${data.topology?.replicates ?? '–'} resampling replicates); both flanks must pass the Bonferroni cutoff α / 2B = ${gardCutoff(cutoff)} for B = ${data.breaks.length} breakpoints. Orange: passes.</div>
        </div>
      </div>` : '<div class="empty-events-notice">No topology test for this breakpoint.</div>';
  }
};

export const GARD_DRAWER = {
  emptyTag: 'Breakpoint —',
  row: (data, site) => data.getBreak(site),
  header: gardDrawerHeader,
  tabs: [GARD_TREES_TAB, GARD_TEST_TAB]
};

// Diagnostics: the run, its model fits and the engines behind it

function gardDiagItem(label, value, valueClass = '', itemClass = '') {
  return value == null || value === '' ? '' : `
      <div class="diag-item${itemClass && ` ${itemClass}`}">
        <span class="diag-lbl">${label}</span>
        <span class="diag-val ${valueClass}">${value}</span>
      </div>`;
}

// The file's name, never its path
function gardBasename(path) {
  return path ? escapeSvgText(String(path).split(/[\\/]/).pop()) : null;
}

function gardSeconds(s) {
  if (s == null) return null;
  return s < 10 ? `${s.toFixed(2)} s` : `${Math.round(s)} s`;
}

function gardRateVariation(rv) {
  if (!rv) return null;
  return rv.mode === 'discrete_gamma'
    ? `Discrete gamma, ${rv.classes} classes, α = ${gardFormat(rv.shape_alpha, 3)}`
    : 'None (one rate class)';
}

function gardTreeEngine(engine) {
  if (!engine) return null;
  return `${engine.method === 'ohmic_spectral' ? 'Ohmic spectral' : 'Neighbor joining'} (${engine.complexity})`;
}

function gardTestRun(data) {
  const t = data.topology;
  if (!t?.testType) return null;
  const replicates = t.replicates == null ? '' : `, ${t.replicates} replicates`;
  return `${t.testType.toUpperCase()}${replicates}`;
}

function gardVersion(data) {
  const s = data.software;
  return [s.version && `hyphy ${s.version}`, s.git_commit && `commit ${s.git_commit}`].filter(Boolean).join(' · ') || null;
}

export function gardDiagnostics(data) {
  return [
    gardDiagItem('Analysis Method', 'GARD (Genetic Algorithm for Recombination Detection)'),
    gardDiagItem('Version', gardVersion(data)),
    gardDiagItem('Output format', data.formatLabel),
    gardDiagItem('Sites', data.numSites),
    gardDiagItem('Taxa / Sequences', data.numTaxa),
    gardDiagItem('Alignment', gardBasename(data.dataset.file ?? data.provenance.alignmentFile), 'diag-val-wrap'),
    gardDiagItem('Minimum partition size', data.dataset.minPartition),
    gardDiagItem('Rate variation', gardRateVariation(data.rateVariation)),
    gardDiagItem('Tree engine', gardTreeEngine(data.treeEngine)),
    gardDiagItem('Topology test', gardTestRun(data)),
    gardDiagItem('Models evaluated', data.modelsEvaluated?.toLocaleString()),
    gardDiagItem('Run time', gardSeconds(data.timeElapsed)),
    gardDiagItem('Alignment SHA-256', data.provenance.alignmentSha256, 'diag-val-wrap font-mono'),
    gardDiagItem('Method', data.description, 'diag-val-prose', 'diag-item-wide'),
    gardDiagItem('Citation', data.citation, 'diag-val-prose', 'diag-item-wide')
  ].filter(Boolean).join('');
}

// The three AICc values the verdicts rest on, and the Ohmic network's when it ran
export function gardFitsTable(data) {
  const rows = [
    ['Baseline (no breakpoints)', data.baseline.logL, data.baselineAicc],
    ['Single tree, all partitions', null, data.singleTreeAicc],
    [`Best model (${gardPlural(data.breaks.length, 'breakpoint')})`, data.segments.every(s => s.logL != null) ? data.segments.reduce((sum, s) => sum + s.logL, 0) : null, data.bestAicc]
  ];
  if (data.ohmic) rows.push(['Ohmic network', data.ohmic.logL, data.ohmic.aicc]);
  const body = rows.map(([name, logL, aicc]) => `
              <tr><td>${name}</td><td class="font-mono">${gardFormat(logL, 2)}</td><td class="font-mono">${gardFormat(aicc, 2)}</td></tr>`).join('');
  return `
        <div class="diag-item">
          <span class="diag-lbl">Model fits</span>
          <table class="diag-fits-table">
            <thead><tr><th>Model</th><th>log L</th><th><span class="case-literal">AICc</span></th></tr></thead>
            <tbody>${body}
            </tbody>
          </table>
        </div>`;
}

function gardFormatNotes(data) {
  const notes = [
    'Breakpoint support is the Akaike weight of the models in which a breakpoint falls after the site; the best model is the one with the lowest AICc.',
    'A breakpoint is incongruent when each flank\'s tree rejects its neighbour\'s by KH or SH at α / 2B, B being the number of breakpoints; calls are recomputed at the live α. A change in branch lengths alone is rate heterogeneity.',
    'Robinson–Foulds distances compare leaf-name splits, so they do not depend on rooting; a moved taxon is one whose sister clade differs between the two trees.'
  ];
  if (data.legacy) notes.push(`Topology and Ohmic blocks are unavailable or partial in ${data.formatLabel} output; the tree-by-partition matrix is not reported.`);
  else if (!data.capabilities.topology) notes.push('This file carries no per-breakpoint topology tests, so every breakpoint is shown untested.');
  return `<ul class="diag-notes">${notes.map(n => `<li>${n}</li>`).join('')}</ul>`;
}

export function gardDiagnosticsPanel(el, data) {
  el.insertAdjacentHTML('beforeend', `${gardFitsTable(data)}
        <div class="diag-item">
          <span class="diag-lbl">Format notes</span>
          ${gardFormatNotes(data)}
        </div>`);
}

export const GARD_CONFIG = {
  id: 'gard',
  Parser: GardResults,
  controller: GARD_ALPHA,
  row: (data, site) => data.getBreak(site),
  defaultFilter: GARD_TABLE.defaultFilter,
  formatLabel: (data) => data.formatLabel,
  kpis: (parent, hooks) => new KpiBar(parent, {
    cards: GARD_KPI_CARDS,
    controller: GARD_ALPHA,
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    onThresholdChange: hooks.onThresholdChange,
    onStatisticChange: hooks.onStatisticChange
  }),
  hero: gardHero,
  table: (parent, hooks) => new ResultsTable(parent, Object.assign({}, GARD_TABLE, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    filterType: hooks.filterType,
    onSelect: (site, row, action) => hooks.onSiteSelect(site, action)
  })),
  secondary: (parent) => new GardTreeCard(parent),
  drawer: (container, hooks) => new Drawer(container, Object.assign({}, GARD_DRAWER, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    siblings: hooks.siblings,
    onStep: hooks.onStep,
    onClose: hooks.onClose
  })),
  diagnosticsTitle: 'Model fits &amp; run details',
  diagnostics: gardDiagnostics,
  newick: (data) => data.treeNewick,
  diagnosticsPanel: gardDiagnosticsPanel
};
