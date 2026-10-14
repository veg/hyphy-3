// SLAC method config: site-level selection from substitutions counted over the maximum-
// likelihood ancestral reconstruction. Each codon gets two one-sided binomial p-values on
// its fractional counts: p(dN>dS) for too few synonymous substitutions given the
// synonymous opportunity, p(dN<dS) for too many. A site is diversifying when
// p(dN>dS) ≤ t, otherwise purifying when p(dN<dS) ≤ t; the two never both pass below
// t = 0.5 because they sum to at least 1. The p-values are uncorrected.

import { SLACData, slacSiteClass, SLAC_LENGTH_KEY } from '../data/slac-data.js';
import { ALIGNMENT_KPI, CODON_COLUMN, SITE_P_THRESHOLD, sitePercent, siteRow } from './fel.js';
import { KpiBar } from '../components/kpi-bar.js';
import { SiteTrackCanvas } from '../components/site-track.js';
import { EvidenceBand } from '../components/evidence-band.js';
import { OverviewTrack } from '../components/overview-track.js';
import { evidenceTipRow } from '../components/evidence-plot.js';
import { roundedDomain, roundedTicks } from '../components/posterior-plot.js';
import { ResultsTable, closestCalls } from '../components/results-table.js';
import { Drawer } from '../components/drawer.js';
import { BranchTree } from '../components/branch-tree.js';
import { formatThreshold, formatSigned } from '../core/format.js';

// The slider runs on round stops (0.005, 0.010, ... 0.200), so every preset and the
// maximum are reachable and the cutoff applied is the one the label prints
export const SLAC_P_THRESHOLD = Object.assign({}, SITE_P_THRESHOLD, {
  min: 0.005,
  caption: 'Uncorrected per-site p-values'
});

const SLAC_CALLS = {
  diversifying: { tone: 'div', glyph: '▲', label: 'Diversifying', badge: '▲ Diversifying selection', badgeClass: 'omega-div' },
  purifying: { tone: 'pur', glyph: '▼', label: 'Purifying', badge: '▼ Purifying selection', badgeClass: 'omega-pur' },
  neutral: { tone: 'neutral', glyph: '—', label: 'Neutral', badge: 'Neutral', badgeClass: '' },
  invariant: { tone: 'neutral', glyph: '○', label: 'Invariant', badge: 'Invariant (no substitutions)', badgeClass: 'omega-invariant' }
};

export function slacCall(r, t) {
  return SLAC_CALLS[slacSiteClass(r, t)];
}

// Notation inside uppercased labels keeps its case (dS, dN, p(dN>dS))
function slacLiteral(text) {
  return `<span class="case-literal">${text}</span>`;
}

// Inferred counts are fractional (ties split between paths): always one decimal
export function slacCount(n) {
  return n.toFixed(1);
}

// A p-value as the method reports it: four decimals, exponential below 1e-3. An exact 0
// is the degenerate P_S = 0 case, not underflow.
export function slacP(p) {
  if (p === 0) return '0';
  return p < 0.001 ? p.toExponential(2) : p.toFixed(4);
}

// Rates that are undefined where ES = 0 read as a dash, never 0
function slacRate(v, digits = 3) {
  return v == null ? '—' : v.toFixed(digits);
}

function slacSignedOrDash(v) {
  return v == null ? '—' : formatSigned(v);
}

// The degenerate call: no synonymous opportunity, yet non-synonymous changes
export function slacDegenerate(r) {
  return r.dSUndefined && r.N > 0;
}

// The two tests' names, escaped for markup
const SLAC_P_POS = 'p(dN&gt;dS)';
const SLAC_P_NEG = 'p(dN&lt;dS)';

const SLAC_DS_UNDEFINED = 'dS undefined: no synonymous opportunity at this site (ES = 0)';

export const SLAC_KPI_CARDS = [
  Object.assign({}, ALIGNMENT_KPI, {
    sub: (stats, { data }) => `${data.numTaxa} Taxa · ΣS ${slacCount(data.totalS)} · ΣN ${slacCount(data.totalN)}`
  }),
  {
    id: 'div',
    label: 'Diversifying (Positive)',
    value: (stats) => `▲ ${stats.diversifying}`,
    sub: (stats) => `${sitePercent(stats.diversifying, stats)}% of sites (dN > dS)`,
    tone: 'accent'
  },
  {
    id: 'pur',
    label: 'Purifying (Constraint)',
    value: (stats) => `▼ ${stats.purifying}`,
    sub: (stats) => `${sitePercent(stats.purifying, stats)}% of sites (dN < dS)`,
    tone: 'purifying'
  },
  {
    id: 'omega',
    label: 'Global dN/dS (ω)',
    value: (stats, { data }) => (data.globalOmega == null ? '—' : data.globalOmega.toFixed(3)),
    sub: (stats, { data }) => (data.globalOmega == null
      ? `Not reported in ${data.formatLabel} output`
      : 'MG94xREV fit · context, not a test'),
    tone: (stats, { data }) => (data.globalOmega == null ? 'muted' : 'ink')
  }
];

// Hero: the mirrored evidence landscape (p(dN>dS) up, p(dN<dS) down, ticks read as p) and
// the scaled dN − dS effect bars

export function slacNegLog10(p) {
  return -Math.log10(Math.max(p, 1e-15));
}

// Symmetric about zero, at least down to p = 0.001 either way. A degenerate p(dN>dS) = 0
// carries no magnitude, so it does not set the scale; its mark is pinned at the top edge.
export function slacEvidenceDomain(rows) {
  let max = 0;
  for (const r of rows) {
    max = Math.max(max, slacNegLog10(r.pNeg));
    if (!slacDegenerate(r)) max = Math.max(max, slacNegLog10(r.pPos));
  }
  const extent = Math.max(3, Math.ceil(max * 1.1));
  return [-extent, extent];
}

// Whole decades, thinned to about four either side of zero
export function slacEvidenceTicks(min, max) {
  const step = max <= 4 ? 1 : Math.ceil(max / 4);
  const ticks = [];
  for (let v = -Math.floor(max / step) * step; v <= max; v += step) ticks.push(v);
  return ticks;
}

// A tick at ±k reads as the p-value 10^-k it stands for: decimals while the axis stops
// at 0.001, powers of ten throughout once it runs further, so one axis never mixes them
export function slacPTick(v, max) {
  const k = Math.round(Math.abs(v));
  if (k === 0) return '1';
  return max > 3 ? `1e-${k}` : String(Number((10 ** -k).toPrecision(1)));
}

const slacDivTone = (r, track) => (slacSiteClass(r, track.threshold) === 'diversifying' ? 'div' : null);
const slacPurTone = (r, track) => (slacSiteClass(r, track.threshold) === 'purifying' ? 'pur' : null);

function slacEffectTone(r, track) {
  const tone = slacCall(r, track.threshold).tone;
  return tone === 'neutral' ? null : tone;
}

// A site with no inferred substitutions is plotted as an empty ring on the baseline; the
// degenerate ES = 0 call carries a warning beside its upward mark, which is held at the
// top edge since slacEvidenceDomain leaves it out
export function slacMarker(r, track, series) {
  if (r.noSubs) return 'hollow';
  return slacDegenerate(r) && series === 0 ? 'flag' : null;
}

// Both p-values on one line, the one that passes at t in its direction's hue
export function slacTipPValues(r, t) {
  const cls = slacSiteClass(r, t);
  const shown = (p, passes, hue) => (passes ? `<span class="${hue}">${slacP(p)}</span>` : slacP(p));
  return `${shown(r.pPos, cls === 'diversifying', 'tip-pass-div')} / ${shown(r.pNeg, cls === 'purifying', 'tip-pass-pur')}`;
}

function slacTooltip(r, track) {
  const call = slacCall(r, track.threshold);
  const rows = [
    evidenceTipRow('S / N observed', `${slacCount(r.S)} / ${slacCount(r.N)}`),
    evidenceTipRow('ES / EN expected', `${r.ES.toFixed(2)} / ${r.EN.toFixed(2)}`),
    evidenceTipRow('dS', slacRate(r.dS)),
    evidenceTipRow('dN', slacRate(r.dN)),
    evidenceTipRow('dN − dS (scaled)', slacSignedOrDash(r.dNmdS)),
    evidenceTipRow(`${SLAC_P_POS} / ${SLAC_P_NEG}`, slacTipPValues(r, track.threshold))
  ];
  const notes = [];
  if (r.noSubs) notes.push('No substitutions inferred, so not informative');
  if (r.dSUndefined) notes.push(`⚠ ${SLAC_DS_UNDEFINED}`);
  if (slacDegenerate(r)) notes.push('p(dN&gt;dS) = 0 is degenerate (P_S = 0); shown as the method reports it');
  const aa = r.consensusAA ? ` · ${r.consensusAA}` : '';
  return `
      <div class="tip-header">
        <span class="tip-title">Codon #${r.site}${aa}</span>
        <span class="tip-badge tip-${call.tone}">${call.badge}</span>
      </div>
      <div class="tip-grid">${rows.join('')}</div>
      ${notes.map(n => `<div class="tip-row tip-sub"><span class="tip-val-sub">${n}</span></div>`).join('')}
      <div class="tip-hint">Click to inspect the counting test</div>
    `;
}

const SLAC_EVIDENCE_SUBTITLE = 'Diversifying p-values up, purifying down, per codon (−log₁₀ scale)';

// Fresh bands per plot, since a band caches the points it last drew
export function slacTrackModes() {
  return [
    {
      id: 'evidence',
      label: 'Evidence',
      subtitle: () => SLAC_EVIDENCE_SUBTITLE,
      band: new EvidenceBand({
        value: (r) => slacNegLog10(r.pPos),
        tone: slacDivTone,
        mirror: { value: (r) => (r.noSubs ? null : slacNegLog10(r.pNeg)), tone: slacPurTone },
        marker: slacMarker,
        domain: slacEvidenceDomain,
        ticks: slacEvidenceTicks,
        tickFormat: slacPTick,
        ruleLabel: (t) => `p ≤ ${formatThreshold(t)}`,
        ruleLabelAt: 'end',
        label: ['more non-synonymous ▲', 'more synonymous ▼'],
        shape: 'dot',
        neutralAlpha: 0.45,
        toneRing: true,
        tooltip: slacTooltip
      })
    },
    {
      id: 'effect',
      label: 'Effect: scaled dN − dS',
      subtitle: () => 'dN − dS ÷ branch length contributing at the site; coloured where significant at the cutoff',
      band: new EvidenceBand({
        value: (r) => r.dNmdS,
        tone: slacEffectTone,
        marker: (r) => (r.noSubs ? 'hollow' : null),
        domain: (rows) => roundedDomain(rows.map(r => r.dNmdS).filter(v => v != null)),
        ticks: roundedTicks,
        tickFormat: (v) => String(v),
        rule: null,
        label: 'scaled dN − dS',
        shape: 'bar',
        toneRing: true,
        tooltip: slacTooltip
      })
    }
  ];
}

// The overview brush's profile: the stronger of the two tails, in the call's hue
export const SLAC_BRUSH = {
  minSites: 400,
  maxValue: 3,
  unit: 'Codons',
  value: (r) => Math.max(slacNegLog10(r.pPos), slacNegLog10(r.pNeg)),
  tone: (r, t) => ({ div: 'accent', pur: 'purifying' })[slacCall(r, t).tone] ?? null,
  rule: (t) => slacNegLog10(t)
};

export function slacHero(parent, hooks) {
  return new OverviewTrack(parent, {
    track: (container) => new SiteTrackCanvas(container, {
      title: 'Selection Evidence Landscape',
      subtitle: SLAC_EVIDENCE_SUBTITLE,
      modes: slacTrackModes(),
      threshold: hooks.threshold,
      exportName: 'hyphy_slac_evidence_plot',
      onSiteSelect: (site) => hooks.onSiteSelect(site)
    }),
    brush: SLAC_BRUSH,
    onSiteSelect: (site) => hooks.onSiteSelect(site)
  });
}

// Results table

function slacCallBadge(r, t) {
  const call = slacCall(r, t);
  const flag = r.dSUndefined ? ` <span class="slac-flag" title="${SLAC_DS_UNDEFINED}">⚠</span>` : '';
  return `<span class="omega-badge ${call.badgeClass}">${call.glyph} ${call.label}</span>${flag}`;
}

// A p cell: bold in the call's hue when that direction passes
function slacPCell(p, passes, hue) {
  return `<td class="font-mono ${passes ? `font-bold ${hue}` : ''}">${slacP(p)}</td>`;
}

// The inline diverging bar for scaled dN − dS, drawn against the largest |value| in the
// result; it takes the call's hue only when the site is called at t
export function slacEffectBar(r, t, extent) {
  if (r.dNmdS == null || !(extent > 0)) return '<span class="effect-bar"></span>';
  const half = Math.min(1, Math.abs(r.dNmdS) / extent) * 50;
  const side = r.dNmdS >= 0 ? `left: 50%; width: ${half}%;` : `left: ${50 - half}%; width: ${half}%;`;
  const tone = slacCall(r, t).tone;
  const hue = tone === 'neutral' ? '' : ` effect-${tone}`;
  return `<span class="effect-bar"><span class="effect-fill${hue}" style="${side}"></span></span>`;
}

const SLAC_SAMPLED_CSV = ['S', 'N', 'dS', 'dN'].flatMap(k => [
  { header: `${k} median`, value: (r) => r.sampled?.[k]?.median ?? '' },
  { header: `${k} 2.5%`, value: (r) => r.sampled?.[k]?.lo ?? '' },
  { header: `${k} 97.5%`, value: (r) => r.sampled?.[k]?.hi ?? '' }
]);

const slacCsvValue = (v) => v ?? '';

export const SLAC_COLUMNS = [
  Object.assign({}, CODON_COLUMN, { alwaysShown: true }),
  {
    id: 'partition',
    header: 'Partition',
    width: 70,
    availableWhen: (caps) => caps.multiPartition,
    cell: (r) => `<td class="font-mono">${r.partition}</td>`,
    csv: [{ header: 'Partition', value: (r) => r.partition }],
    sort: (r) => r.partition
  },
  {
    id: 'class',
    header: 'Selection',
    width: 130,
    alwaysShown: true,
    title: 'The call at the current cutoff; sorts by the smaller p-value',
    cell: (r, { threshold }) => `<td>${slacCallBadge(r, threshold)}</td>`,
    csv: [
      { header: 'Selection', value: (r, { threshold }) => slacCall(r, threshold).label },
      { header: 'dS undefined', value: (r) => (r.dSUndefined ? 'yes' : '') }
    ],
    sort: (r) => Math.min(r.pPos, r.pNeg)
  },
  {
    id: 'aa',
    header: 'AA',
    width: 50,
    title: 'Consensus amino acid over the leaves',
    cell: (r) => `<td class="font-mono">${r.consensusAA ?? '—'}</td>`,
    csv: [{ header: 'Consensus AA', value: (r) => slacCsvValue(r.consensusAA) }],
    sort: (r) => r.consensusAA ?? ''
  },
  {
    id: 'sn',
    header: 'S / N',
    width: 90,
    title: 'Inferred synonymous / non-synonymous substitutions',
    cell: (r) => `<td class="font-mono">${slacCount(r.S)} / ${slacCount(r.N)}</td>`,
    csv: [{ header: 'S', value: (r) => r.S }, { header: 'N', value: (r) => r.N }],
    sort: (r) => r.S + r.N
  },
  {
    id: 'esen',
    header: 'ES / EN',
    width: 100,
    title: 'Expected synonymous / non-synonymous sites',
    hiddenByDefault: true,
    cell: (r) => `<td class="font-mono">${r.ES.toFixed(2)} / ${r.EN.toFixed(2)}</td>`,
    csv: [{ header: 'ES', value: (r) => r.ES }, { header: 'EN', value: (r) => r.EN }],
    sort: (r) => r.ES / Math.max(1e-12, r.ES + r.EN)
  },
  {
    id: 'ds',
    header: slacLiteral('dS'),
    width: 70,
    title: 'Inferred synonymous substitution rate',
    cell: (r) => `<td class="font-mono">${slacRate(r.dS)}</td>`,
    csv: [{ header: 'dS', value: (r) => slacCsvValue(r.dS) }],
    sort: (r) => r.dS ?? -Infinity
  },
  {
    id: 'dn',
    header: slacLiteral('dN'),
    width: 70,
    title: 'Inferred non-synonymous substitution rate',
    cell: (r) => `<td class="font-mono">${slacRate(r.dN)}</td>`,
    csv: [{ header: 'dN', value: (r) => slacCsvValue(r.dN) }],
    sort: (r) => r.dN ?? -Infinity
  },
  {
    id: 'effect',
    header: `${slacLiteral('dN − dS')} (scaled)`,
    width: 150,
    title: 'dN − dS divided by the length of branches contributing at the site',
    cell: (r, { threshold, data }) => `<td class="font-mono"><div class="effect-cell"><span>${slacSignedOrDash(r.dNmdS)}</span>${slacEffectBar(r, threshold, data.effectExtent)}</div></td>`,
    csv: [{ header: 'dN-dS (scaled)', value: (r) => slacCsvValue(r.dNmdS) }],
    sort: (r) => r.dNmdS ?? -Infinity
  },
  {
    id: 'ppos',
    header: slacLiteral(SLAC_P_POS),
    width: 90,
    cell: (r, { threshold }) => slacPCell(r.pPos, slacSiteClass(r, threshold) === 'diversifying', 'text-orange'),
    csv: [{ header: 'p(dN>dS)', value: (r) => r.pPos }],
    sort: (r) => r.pPos
  },
  {
    id: 'pneg',
    header: slacLiteral(SLAC_P_NEG),
    width: 90,
    cell: (r, { threshold }) => slacPCell(r.pNeg, slacSiteClass(r, threshold) === 'purifying', 'text-blue'),
    csv: [{ header: 'p(dN<dS)', value: (r) => r.pNeg }],
    sort: (r) => r.pNeg
  },
  {
    id: 'detail',
    csv: [
      { header: 'P_S', value: (r) => slacCsvValue(r.P_S) },
      { header: 'Total branch length', value: (r) => slacCsvValue(r.tbl) }
    ]
  },
  {
    id: 'sampled',
    availableWhen: (caps) => caps.sampling,
    csv: SLAC_SAMPLED_CSV
  },
  {
    id: 'actions',
    header: 'Inspect',
    width: 80,
    align: 'right',
    cell: (r) => `<td style="text-align: right;"><button type="button" class="btn-inspect" data-action="summary" data-site="${r.site}">Inspect</button></td>`
  }
];

const slacCalled = (cls) => (r, t) => slacSiteClass(r, t) === cls;

export const SLAC_FILTERS = [
  {
    id: 'diversifying',
    key: 'div',
    glyph: '▲',
    label: 'Diversifying',
    unit: 'diversifying site',
    evidence: 'more non-synonymous substitutions',
    statistic: SLAC_P_POS,
    p: (r) => r.pPos,
    opposite: 'purifying',
    sort: { column: 'ppos', dir: 'asc' },
    count: (stats) => stats.diversifying,
    test: slacCalled('diversifying'),
    subtitle: (n, t) => `Showing ${n} diversifying sites (p(dN>dS) ≤ ${formatThreshold(t)})`
  },
  {
    id: 'purifying',
    key: 'pur',
    glyph: '▼',
    label: 'Purifying',
    unit: 'purifying site',
    evidence: 'fewer non-synonymous substitutions',
    statistic: SLAC_P_NEG,
    p: (r) => r.pNeg,
    opposite: 'diversifying',
    sort: { column: 'pneg', dir: 'asc' },
    count: (stats) => stats.purifying,
    test: slacCalled('purifying'),
    subtitle: (n, t) => `Showing ${n} purifying sites (p(dN<dS) ≤ ${formatThreshold(t)})`
  },
  {
    id: 'significant',
    key: 'sig',
    label: 'All Significant',
    unit: 'site',
    statistic: 'p',
    p: (r) => Math.min(r.pPos, r.pNeg),
    sort: { column: 'class', dir: 'asc' },
    count: (stats) => stats.significant,
    test: (r, t) => ['diversifying', 'purifying'].includes(slacSiteClass(r, t)),
    subtitle: (n, t) => `Showing ${n} significant sites (either p ≤ ${formatThreshold(t)})`
  },
  {
    id: 'all',
    key: 'all',
    label: 'All Sites',
    unit: 'site',
    count: (stats) => stats.totalSites,
    test: () => true,
    subtitle: (n) => `Showing all ${n} codon sites`
  }
];

// What the empty filter means, a jump to the mirror direction when it has sites, and the
// sites closest to passing
export function slacEmptyState(t, { filter, stats, rows }) {
  const title = `No ${filter.unit} reaches ${slacLiteral(filter.statistic)} ≤ ${formatThreshold(t)}`;
  const state = { icon: filter.glyph, title, fallback: closestCalls(rows, filter.p) };
  if (filter.evidence) {
    state.desc = `No codon has significantly ${filter.evidence} than its synonymous opportunity predicts at this cutoff.`;
  }
  const other = SLAC_FILTERS.find(f => f.id === filter.opposite);
  const otherCount = other ? other.count(stats) : 0;
  if (otherCount > 0) {
    state.jump = { filter: other.id, label: `View ${other.label} Sites (${other.glyph} ${otherCount})` };
  }
  return state;
}

export const SLAC_TABLE = {
  title: 'Candidate Sites',
  placeholder: 'Filtered by selection call',
  rowKey: 'site',
  rows: (data) => data.siteResults,
  columns: SLAC_COLUMNS,
  filters: SLAC_FILTERS,
  defaultFilter: 'diversifying',
  defaultSort: { column: 'codon', dir: 'asc' },
  search: {
    placeholder: 'Search codon or AA...',
    match: (r, q) => String(r.site).includes(q) || (r.consensusAA || '').toLowerCase() === q,
    noMatch: { title: 'No Matching Sites', desc: 'No codon sites match your current filter and search query.' }
  },
  emptyState: slacEmptyState,
  exportName: 'hyphy_slac_sites'
};

// Drawer: the call at the live cutoff, the counting test behind it, the site's counts and
// rates, their spread over sampled ancestors, and the leaves' codons

// The header badge is the table's short call ('○ Invariant'); the counting test spells
// out what an invariant site lacks
function slacDrawerHeader(r, { threshold }) {
  const call = slacCall(r, threshold);
  const aa = r.consensusAA ? ` · ${r.consensusAA}` : '';
  const flag = r.dSUndefined ? ' · ⚠ dS undefined' : '';
  const badge = call === SLAC_CALLS.invariant ? `${call.glyph} ${call.label}` : call.badge;
  return { tag: `Codon #${r.site}${aa}`, tone: call.tone, badge: `${badge}${flag}` };
}

function slacPercent(x) {
  return `${Math.round(x * 100)}%`;
}

// The expected synonymous share: P_S, or ES / (ES + EN) when the file lacks it
function slacExpectedSyn(r) {
  return r.P_S ?? r.ES / Math.max(1e-12, r.ES + r.EN);
}

// The plain-language reading of the test, quoting the p of the direction the counts lean
export function slacCountingSentence(r, t) {
  if (r.noSubs) return 'No substitutions inferred at this site; SLAC cannot test it.';
  if (slacDegenerate(r)) return 'No synonymous opportunity (ES = 0), so p(dN&gt;dS) = 0 is degenerate; the call is shown as the method reports it.';
  const total = r.S + r.N;
  const expectedN = 1 - slacExpectedSyn(r);
  const leansN = r.N / total > expectedN;
  const [p, label, hue] = leansN ? [r.pPos, SLAC_P_POS, 'text-orange'] : [r.pNeg, SLAC_P_NEG, 'text-blue'];
  const shown = p <= t ? `<strong class="${hue}">${slacP(p)}</strong>` : `<strong>${slacP(p)}</strong>`;
  return leansN
    ? `Of ${slacCount(total)} inferred substitutions, ${slacCount(r.N)} are non-synonymous (${slacPercent(r.N / total)}); ${slacPercent(expectedN)} would be expected if neutral → ${label} = ${shown}.`
    : `Of ${slacCount(total)} inferred substitutions, ${slacCount(r.S)} are synonymous (${slacPercent(r.S / total)}); ${slacPercent(1 - expectedN)} would be expected if neutral → ${label} = ${shown}.`;
}

// A synonymous share and its complement as whole percentages that add up to 100
export function slacShares(syn) {
  const s = Math.round(syn * 100);
  return `S ${s}% · N ${100 - s}%`;
}

// Two 100% bars, expected and observed synonymous | non-synonymous shares. The side where
// the observed share exceeds the expected one is outlined in its direction's hue.
function slacCountingSection(r, t) {
  const expS = slacExpectedSyn(r);
  const total = r.S + r.N;
  const obsS = total > 0 ? r.S / total : null;
  const bar = (syn, synClass = '', nonClass = '') => `<span class="count-bar">`
    + `<span class="count-syn${synClass}" style="width: ${syn * 100}%;"></span>`
    + `<span class="count-non${nonClass}" style="width: ${(1 - syn) * 100}%;"></span></span>`;
  const observed = obsS == null
    ? '<span class="count-bar count-bar-empty"></span>'
    : bar(obsS, obsS > expS ? ' count-excess-pur' : '', obsS < expS ? ' count-excess-div' : '');
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Counting test</div>
          <div class="count-test">
            <span class="count-label">Expected</span>${bar(expS)}<span class="count-val font-mono">${slacShares(expS)}</span>
            <span class="count-label">Observed</span>${observed}<span class="count-val font-mono">${obsS == null ? 'none' : slacShares(obsS)}</span>
          </div>
          <div class="count-sentence">${slacCountingSentence(r, t)}</div>
          <div class="empty-events-notice">Binomial tails on fractional counts (HyPhy's extended binomial): n = S + N, success = synonymous with probability P_S = ES / (ES + EN). Uncorrected for multiple sites.</div>
        </div>`;
}

function slacStatCard(label, value, sub, hue = '') {
  return `
            <div class="drawer-stat-card">
              <div class="dstat-label">${label}</div>
              <div class="dstat-value font-mono font-bold ${hue}">${value}</div>
              <div class="dstat-sub">${sub}</div>
            </div>`;
}

// A p card's caption: the call when this direction passes, else that it does not
function slacPCaption(passes, call, t) {
  return passes ? `${call} at this cutoff` : `not significant at p ≤ ${formatThreshold(t)}`;
}

function slacStatsSection(r, t) {
  const cls = slacSiteClass(r, t);
  const dsSub = r.dSUndefined ? '⚠ undefined: ES = 0' : 'synonymous rate';
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Counts and rates</div>
          <div class="drawer-stats-grid">${[
            slacStatCard('S (observed)', slacCount(r.S), 'inferred synonymous substitutions'),
            slacStatCard('N (observed)', slacCount(r.N), 'inferred non-synonymous substitutions'),
            slacStatCard('ES (expected)', r.ES.toFixed(3), 'synonymous sites'),
            slacStatCard('EN (expected)', r.EN.toFixed(3), 'non-synonymous sites'),
            slacStatCard(slacLiteral('dS'), slacRate(r.dS), dsSub),
            slacStatCard(slacLiteral('dN'), slacRate(r.dN), 'non-synonymous rate'),
            slacStatCard(`${slacLiteral('dN − dS')} (scaled)`, slacSignedOrDash(r.dNmdS), 'divided by the total branch length'),
            slacStatCard('Total branch length', slacRate(r.tbl, 4), 'length of branches contributing at this site, used to scale dN − dS'),
            slacStatCard(slacLiteral(SLAC_P_POS), slacP(r.pPos), slacPCaption(cls === 'diversifying', 'diversifying', t), cls === 'diversifying' ? 'text-orange' : ''),
            slacStatCard(slacLiteral(SLAC_P_NEG), slacP(r.pNeg), slacPCaption(cls === 'purifying', 'purifying', t), cls === 'purifying' ? 'text-blue' : '')
          ].join('')}
          </div>
        </div>`;
}

// One dot-and-whisker row: the 2.5–97.5% range as a line, the median as a tick, the point
// estimate as a dot, on a 0..max scale shared within a pair (S with N, dS with dN)
export function slacWhisker(point, q, max, width = 180) {
  const x = (v) => 4 + (Math.max(0, v) / max) * (width - 8);
  if (!q || q.lo == null || q.hi == null || point == null) return '';
  return `<svg class="whisker" width="${width}" height="14" viewBox="0 0 ${width} 14">`
    + `<line x1="4" y1="7" x2="${width - 4}" y2="7" class="whisker-axis"/>`
    + `<line x1="${x(q.lo)}" y1="7" x2="${x(q.hi)}" y2="7" class="whisker-range"/>`
    + (q.median == null ? '' : `<line x1="${x(q.median)}" y1="2" x2="${x(q.median)}" y2="12" class="whisker-median"/>`)
    + `<circle cx="${x(point)}" cy="7" r="3" class="whisker-point"/></svg>`;
}

function slacUncertaintySection(r, data) {
  if (!r.sampled) {
    return data.settings.samples === 0
      ? `
        <div class="drawer-section">
          <div class="empty-events-notice">No ancestral reconstructions were sampled in this run (samples = 0), so there is no uncertainty range.</div>
        </div>`
      : '';
  }
  const pairMax = (a, b) => Math.max(1e-9, ...[a, b].flatMap(k => [r[k], r.sampled[k]?.hi]).filter(v => v != null));
  const max = { S: pairMax('S', 'N'), N: pairMax('S', 'N'), dS: pairMax('dS', 'dN'), dN: pairMax('dS', 'dN') };
  const digits = { S: 1, N: 1, dS: 3, dN: 3 };
  const rows = ['S', 'N', 'dS', 'dN'].map(k => {
    const q = r.sampled[k];
    const fmt = (v) => (v == null ? '—' : v.toFixed(digits[k]));
    const text = q ? `${fmt(r[k])} · median ${fmt(q.median)} · 95% ${fmt(q.lo)}–${fmt(q.hi)}` : '—';
    return `
            <span class="whisker-label">${slacLiteral(k)}</span>${slacWhisker(r[k], q, max[k])}<span class="whisker-text font-mono">${text}</span>`;
  }).join('');
  const samples = data.settings.samples ? ` (${data.settings.samples} samples)` : '';
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Ancestral uncertainty${samples}</div>
          <div class="whisker-grid">${rows}
          </div>
          <div class="empty-events-notice">Dot: the estimate from the maximum-likelihood ancestors. Tick and line: median and 2.5–97.5% range over sampled ancestral reconstructions.</div>
        </div>`;
}

// The leaves' residues at the site, most common first
function slacResiduesSection(r, data) {
  const counts = new Map();
  for (const { codon, aa } of data.siteCodons(r.site)) {
    if (codon === '---') continue;
    const entry = counts.get(codon) || { codon, aa, n: 0 };
    entry.n++;
    counts.set(codon, entry);
  }
  if (!counts.size) return '';
  const chips = [...counts.values()].sort((a, b) => b.n - a.n || a.codon.localeCompare(b.codon)).map(c => `
            <div class="extant-chip">
              <span class="extant-codon font-mono">${c.codon} ${c.aa}</span>
              <span class="extant-count font-mono">×${c.n}</span>
            </div>`).join('');
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Observed extant residues</div>
          <div class="extant-chips-grid">${chips}
          </div>
        </div>`;
}

const SLAC_SUMMARY_TAB = {
  id: 'summary',
  label: 'Summary',
  comment: 'Counting test, counts, sampled ancestors and residues',
  render: (panel, r, { data, threshold }) => {
    panel.innerHTML = `<div class="drawer-stack">${slacCountingSection(r, threshold)}${slacStatsSection(r, threshold)}${slacUncertaintySection(r, data)}${slacResiduesSection(r, data)}
      </div>`;
  }
};

// 'Human  GAG · E', the leaf's codon at the site after its name
export function slacLeafLabel(name, entry) {
  return entry ? `${name}  ${entry.codon} · ${entry.aa}` : name;
}

const SLAC_TREE_TAB = {
  id: 'tree',
  label: 'Tree',
  comment: 'The tree with each leaf\'s codon at the site',
  availableWhen: (data) => data.capabilities.tree && data.capabilities.leafCodons,
  mount: (panel) => {
    panel.innerHTML = `
        <div class="drawer-section">
          <div class="drawer-section-title"><span>Leaf codons</span><button type="button" class="plot-tool-btn slac-tree-svg" title="Export Vector SVG">SVG</button></div>
          <div class="slac-tree-host"></div>
          <div class="empty-events-notice">Branch lengths from the MG94xREV fit the ancestors were reconstructed on.</div>
        </div>`;
    const view = { codons: new Map() };
    view.tree = new BranchTree(panel.querySelector('.slac-tree-host'), {
      lengthKey: SLAC_LENGTH_KEY,
      leafLabel: (node) => slacLeafLabel(node.name, view.codons.get(node.name)),
      exportName: 'hyphy_slac_site_tree'
    });
    view.setData = (data) => view.tree.setTree(data.treeNewick, data.nodeMeta);
    panel.querySelector('.slac-tree-svg').addEventListener('click', () => view.tree.exportSVG());
    return view;
  },
  render: (panel, r, { data, view }) => {
    view.codons = new Map(data.siteCodons(r.site).map(c => [c.taxon, c]));
    view.tree.options.exportName = `hyphy_slac_site${r.site}_tree`;
    view.tree.resize(view.tree.width || view.tree.container.clientWidth);
  }
};

export const SLAC_DRAWER = {
  emptyTag: 'Codon #—',
  row: siteRow,
  header: slacDrawerHeader,
  tabs: [SLAC_SUMMARY_TAB, SLAC_TREE_TAB]
};

// Diagnostics: the run, the test setup and the fits the reconstruction rests on

function slacDiagItem(label, value, valueClass = '', itemClass = '') {
  return value == null || value === '' ? '' : `
      <div class="diag-item${itemClass && ` ${itemClass}`}">
        <span class="diag-lbl">${label}</span>
        <span class="diag-val ${valueClass}">${value}</span>
      </div>`;
}

function slacSeconds(s) {
  if (s == null) return null;
  return s < 10 ? `${s.toFixed(2)} s` : `${Math.round(s)} s`;
}

function slacTestedItem(data) {
  const { tested, background, testedLength } = data.testing;
  if (tested == null) return null;
  const parts = [`${tested} tested · ${background} background`];
  if (data.settings.testBranches) parts.unshift(data.settings.testBranches);
  if (testedLength != null) parts.push(`tested tree length ${testedLength.toFixed(4)}`);
  return parts.join(' · ');
}

// The run's own cutoff, with the counts hyphy3 wrote at it when the file records them
function slacRunItem(data) {
  const run = data.runSummary;
  if (run) {
    return slacDiagItem(`Run cutoff, as written by ${slacLiteral('hyphy3')}`,
      `p ≤ ${formatThreshold(run.threshold)}: ${run.diversifying} diversifying · ${run.purifying} purifying · ${run.neutral} neutral`);
  }
  return slacDiagItem('Run cutoff', data.runThreshold == null ? null : `p ≤ ${formatThreshold(data.runThreshold)}`);
}

export function slacDiagnostics(data) {
  const product = data.format === 'hyphy25' ? 'HyPhy' : 'hyphy3';
  const versions = data.format === 'v3'
    ? `analysis ${data.version.analysis} · hyphy ${data.version.software}`
    : `analysis ${data.version.analysis} · ${product} ${data.version.software}`;
  const citation = data.citation && (data.citation.doi ? `${data.citation.text} doi:${data.citation.doi}` : data.citation.text);
  return [
    slacDiagItem('Analysis Method', 'SLAC (Single-Likelihood Ancestor Counting)'),
    slacDiagItem('Version', versions),
    slacDiagItem('Output format', data.formatLabel),
    slacDiagItem('Codon Sites', data.numSites),
    slacDiagItem('Taxa / Sequences', data.numTaxa),
    slacDiagItem(`Global dN/dS (${slacLiteral('ω')})`, data.globalOmega == null ? null : data.globalOmega.toFixed(3)),
    slacDiagItem('Branches', slacTestedItem(data)),
    slacDiagItem('Ancestor samples', data.settings.samples),
    slacDiagItem('Genetic code', data.settings.geneticCode),
    slacDiagItem('Ambiguities', 'Per-site values use resolved ambiguities (as hyphy3 v3 does)'),
    slacRunItem(data),
    slacDiagItem('Run time', slacSeconds(data.provenance.wallSeconds)),
    slacDiagItem('Command', data.provenance.command, 'diag-val-wrap'),
    slacDiagItem('Alignment SHA-256', data.provenance.alignmentSha256, 'diag-val-wrap font-mono'),
    slacDiagItem('Method', data.description, 'diag-val-prose', 'diag-item-wide'),
    slacDiagItem('Citation', citation, 'diag-val-prose', 'diag-item-wide')
  ].filter(Boolean).join('');
}

// The tree as the diagnostics print it: without {tags}, without the zero lengths v3
// writes (lengths live in the fits) and without v3's 'END' root label
export function slacNewickText(newick) {
  return newick.replace(/:0(?:\.0+)?(?=[,);])/g, '').replace(/\)END;?$/, ');');
}

// Every fit the file reports, with an em dash for a field it lacks
export function slacFitsTable(fits, formatLabel) {
  const cell = (v, digits) => (v == null ? '—' : (digits == null ? String(v) : v.toFixed(digits)));
  const rows = fits.map(f => `
              <tr><td>${f.name}</td><td class="font-mono">${cell(f.logL, 3)}</td><td class="font-mono">${cell(f.params)}</td><td class="font-mono">${cell(f.aicc, 2)}</td></tr>`).join('');
  return `
        <div class="diag-item">
          <span class="diag-lbl">Model fits (${slacLiteral(formatLabel)})</span>
          <table class="diag-fits-table">
            <thead><tr><th>Model</th><th>log L</th><th>Parameters</th><th>${slacLiteral('AICc')}</th></tr></thead>
            <tbody>${rows}
            </tbody>
          </table>
        </div>`;
}

// Notes on what each number means, so none is read across formats
function slacFormatNotes() {
  const notes = [
    'p(dN&gt;dS) and p(dN&lt;dS) are separate one-sided binomial tails on fractional counts; they sum to at least 1, so no site is called both ways. Calls and counts are recomputed at the live cutoff, with no multiple-testing correction.',
    'dS and scaled dN − dS are undefined where ES = 0 and are left blank, never 0. If such a site has non-synonymous substitutions its p(dN&gt;dS) is 0 by construction and is flagged ⚠.',
    'Counting on a reconstruction is conservative: a site with few inferred substitutions (small S + N) cannot reach a small p in either direction.',
    'Scaled dN − dS divides by the length of the branches contributing at the site; that length, the AICc values and the parameter counts are this file\'s own and are not comparable across formats.'
  ];
  return `<ul class="diag-notes">${notes.map(n => `<li>${n}</li>`).join('')}</ul>`;
}

export function slacDiagnosticsPanel(el, data) {
  if (data.fits.length) el.insertAdjacentHTML('beforeend', slacFitsTable(data.fits, data.formatLabel));
  el.insertAdjacentHTML('beforeend', `
        <div class="diag-item">
          <span class="diag-lbl">Format notes</span>
          ${slacFormatNotes()}
        </div>`);
}

export const SLAC_CONFIG = {
  id: 'slac',
  Parser: SLACData,
  controller: SLAC_P_THRESHOLD,
  row: siteRow,
  defaultFilter: SLAC_TABLE.defaultFilter,
  formatLabel: (data) => data.formatLabel,
  kpis: (parent, hooks) => new KpiBar(parent, {
    cards: SLAC_KPI_CARDS,
    controller: SLAC_P_THRESHOLD,
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    onThresholdChange: hooks.onThresholdChange,
    onStatisticChange: hooks.onStatisticChange
  }),
  hero: slacHero,
  table: (parent, hooks) => new ResultsTable(parent, Object.assign({}, SLAC_TABLE, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    filterType: hooks.filterType,
    onSelect: (site, row, action) => hooks.onSiteSelect(site, action)
  })),
  drawer: (container, hooks) => new Drawer(container, Object.assign({}, SLAC_DRAWER, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    siblings: hooks.siblings,
    onStep: hooks.onStep,
    onClose: hooks.onClose
  })),
  diagnosticsTitle: 'Model fits, tree &amp; run details',
  diagnostics: slacDiagnostics,
  newick: (data) => slacNewickText(data.treeNewick),
  diagnosticsPanel: slacDiagnosticsPanel
};
