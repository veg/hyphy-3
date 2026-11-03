// dFEL method config: FEL's per-site alpha vs beta test, fitted under a directional model
// that lets base composition drift along lineages (K regimes) on a likelihood-chosen root.
// A site is diversifying when p ≤ t and β > α, purifying when p ≤ t and α > β; the p-values
// are uncorrected. Under the site answer sits the model the test rests on: the composition
// screen, the K regimes actually fitted (the screen picks K only in auto mode) and the root.

import { DFELData, dfelSiteClass, dfelAtFloor, DFEL_RATE_CEILING, DFEL_REGIMES } from '../data/dfel-data.js';
import { ALIGNMENT_KPI, CODON_COLUMN, LRT_COLUMN, SITE_P_THRESHOLD, sitePercent, siteRow, felEvidenceTone } from './fel.js';
import { KpiBar } from '../components/kpi-bar.js';
import { EvidencePlot, evidenceTipRow } from '../components/evidence-plot.js';
import { TickStripBand } from '../components/tick-strip-band.js';
import { ResultsTable, closestCalls } from '../components/results-table.js';
import { Drawer } from '../components/drawer.js';
import { BranchTree } from '../components/branch-tree.js';
import { compositionTableMarkup, compositionP, gcStripSVG } from '../components/composition-panel.js';
import { rootSupportSVG, rootEdgeLabel, rootWeightPercent, rootCredibleCount } from '../components/root-support-bars.js';
import { PALETTE } from '../data/palette.js';
import { formatThreshold, formatPValue } from '../core/format.js';
import { csvQuote, downloadSVG, downloadSVGAsPNG } from '../core/export.js';

// The slider spans 0.001-0.20 in steps of 0.001, so every preset and any run cutoff written
// to three decimals sits on a stop
export const DFEL_P_THRESHOLD = Object.assign({}, SITE_P_THRESHOLD, {
  step: 0.001,
  caption: 'Uncorrected per-site p-values'
});

const DFEL_CALLS = {
  diversifying: { tone: 'div', glyph: '▲', label: 'Diversifying', badge: '▲ Diversifying selection', badgeClass: 'omega-div' },
  purifying: { tone: 'pur', glyph: '▼', label: 'Purifying', badge: '▼ Purifying selection', badgeClass: 'omega-pur' },
  neutral: { tone: 'neutral', glyph: '—', label: 'Neutral', badge: 'Neutral', badgeClass: '' },
  invariant: { tone: 'neutral', glyph: '○', label: 'Invariant', badge: 'Invariant (no substitutions)', badgeClass: 'omega-invariant' }
};

export function dfelCall(r, t) {
  return DFEL_CALLS[dfelSiteClass(r, t)];
}

// Greek and bracketed notation inside uppercased labels keeps its case
function dfelLiteral(text) {
  return `<span class="case-literal">${text}</span>`;
}

// A rate with its optimiser bounds spelled out: '0' at the floor, '≥1000 (bound)' at the cap
export function dfelRate(v) {
  if (dfelAtFloor(v)) return '0';
  return v >= DFEL_RATE_CEILING ? '≥1000 (bound)' : v.toFixed(2);
}

// omega = β/α: '∞' when α sits at the floor, '0' when β does, '—' when both do
export function dfelOmegaText(omega) {
  if (omega == null) return '—';
  if (omega === Infinity) return '∞';
  if (omega === 0) return '0';
  return omega >= 1000 ? '>1000' : omega.toFixed(2);
}

const DFEL_ROOT_FIXED = new Set(['auto', 'keep']);

// An outgroup was requested by name; whether it was found is root.branch === root.strategy
function dfelOutgroup(root) {
  return Boolean(root.strategy) && !DFEL_ROOT_FIXED.has(root.strategy);
}

// The requested outgroup is missing from the tree and the ML root stands in for it
export function dfelOutgroupMissing(root) {
  return Boolean(root?.branch) && dfelOutgroup(root) && root.branch !== root.strategy;
}

// 'w = 8.1% (weak)': a chosen root under half the Akaike weight is weakly supported
function dfelRootWeight(root) {
  if (root.weight == null) return null;
  return `w = ${rootWeightPercent(root.weight)}${root.weight < 0.5 ? ' (weak)' : ''}`;
}

// The root as one phrase: the chosen edge and its weight in auto mode, the outgroup when one
// was requested (the ML root when it was not found, which the model card warns about), or
// the input root under 'keep', whose written weight and log L are placeholders and never shown
export function dfelRootPhrase(root) {
  if (!root?.branch) return null;
  if (root.strategy === 'keep') return 'root: input (kept)';
  const edge = dfelOutgroup(root) && !dfelOutgroupMissing(root)
    ? `root: outgroup ${root.branch}`
    : `${dfelOutgroupMissing(root) ? 'ML root' : 'root'} ${rootEdgeLabel(root.branch)}`;
  return [edge, dfelRootWeight(root)].filter(Boolean).join(', ');
}

// The root in a few words for the KPI sub line; the edge name stays in the model card
export function dfelRootToken(root) {
  if (!root?.branch) return null;
  if (root.strategy === 'keep') return 'root kept';
  if (dfelOutgroupMissing(root)) return '⚠ outgroup not found';
  return [dfelOutgroup(root) ? 'outgroup root' : 'root', dfelRootWeight(root)].filter(Boolean).join(' ');
}

// The screen and the fit, as two facts: K follows the screen only in auto mode, so a K the
// screen would not pick is marked forced
export function dfelModelVerdict(data) {
  const comp = data.composition;
  if (!comp) return null;
  const K = data.K;
  const screen = comp.heterogeneous ? 'Heterogeneous' : 'Homogeneous';
  const forced = K != null && K !== (comp.heterogeneous ? 3 : 1);
  const fit = `K = ${K ?? '—'}${forced ? ' (forced)' : ''}`;
  return { forced, screen: `${screen} screen`, fit, summary: `${screen} screen · ${fit}` };
}

// The deciding screen test: Overall or CP3, whichever has the smaller p (regimes are
// activated when either falls under the heterogeneity threshold)
export function dfelDecidingTest(comp) {
  const tests = [['Overall', comp.overall], ['CP3', comp.cp3]].filter(([, t]) => t?.p != null);
  if (!tests.length) return null;
  const [label, test] = tests.reduce((a, b) => (b[1].p < a[1].p ? b : a));
  return `screen p = ${compositionP(test.p)} (${label})`;
}

export const DFEL_KPI_CARDS = [
  Object.assign({}, ALIGNMENT_KPI, {
    sub: (stats, { data }) => [`${data.numTaxa} Taxa`, data.patterns != null && `${data.patterns} patterns`].filter(Boolean).join(' · ')
  }),
  {
    id: 'div',
    label: 'Diversifying (Positive)',
    value: (stats) => `▲ ${stats.diversifying}`,
    sub: (stats) => `${sitePercent(stats.diversifying, stats)}% of sites (β > α)`,
    tone: 'accent'
  },
  {
    id: 'pur',
    label: 'Purifying (Constraint)',
    value: (stats) => `▼ ${stats.purifying}`,
    sub: (stats) => `${sitePercent(stats.purifying, stats)}% of sites (α > β)`,
    tone: 'purifying'
  },
  {
    // The screen verdict over the fitted K; a --legacy-json export has neither, so the card
    // is left out rather than shown empty
    id: 'model',
    label: 'Directional model',
    cardClass: 'kpi-card-text',
    available: (data) => data.capabilities.directionalModel,
    value: (stats, { data }) => {
      const verdict = dfelModelVerdict(data);
      return `${verdict.screen}\n${verdict.fit}`;
    },
    sub: (stats, { data }) => [dfelDecidingTest(data.composition), dfelRootToken(data.root)].filter(Boolean).join(' · '),
    tone: 'ink'
  }
];

// Hero: FEL's evidence landscape with the GC/AT direction strip of the reconstructed
// non-synonymous changes under it

// 'nonsyn: 3 GC↑ / 1 AT↑ / 0 ⚡, +1 syn'
export function dfelChangeSummary(r) {
  if (!r.tally) return null;
  if (r.subs.conserved) return 'none (invariant)';
  const t = r.tally;
  return `nonsyn: ${t.gc} GC↑ / ${t.at} AT↑ / ${t.surprising} ⚡, +${t.syn} syn`;
}

export const DFEL_EVIDENCE = {
  tone: felEvidenceTone,
  alpha: (r) => dfelRate(r.alpha),
  badges: { div: '▲ Diversifying Selection', pur: '▼ Purifying Selection' },
  rows: (r) => [
    evidenceTipRow('β (dN)', dfelRate(r.beta)),
    evidenceTipRow('ω (β/α)', dfelOmegaText(r.omega)),
    ...(r.tally ? [evidenceTipRow('Changes', dfelChangeSummary(r))] : [])
  ],
  signed: (r) => (r.beta >= r.alpha ? 1 : -1) * r.negLogP,
  diff: (r) => r.beta - r.alpha,
  subtitles: {
    manhattan: '−log₁₀(p) vs Codon Position (directional model)',
    signed: 'Signed Evidence: +Diversifying / −Purifying',
    diff: 'Absolute Rate Difference: β (dN) − α (dS)'
  },
  hint: 'Click point to inspect the site'
};

export const DFEL_STRIP = {
  label: 'GC/AT direction',
  title: 'Show or hide the direction of each site\'s reconstructed non-synonymous changes',
  available: (data) => data.capabilities.substitutions,
  // On by default only when the model fitted directional regimes
  shown: (data) => (data.K ?? 1) > 1,
  note: (data) => {
    const legend = 'Strip: dark ticks up count GC-ward non-synonymous changes, light ticks down AT-ward; a dot marks a ⚡ surprising change.';
    return (data.K ?? 1) > 1 ? legend : `Direction of reconstructed changes (single Neutral regime fitted). ${legend}`;
  }
};

export function dfelHero(parent, hooks) {
  return new EvidencePlot(parent, {
    preset: DFEL_EVIDENCE,
    threshold: hooks.threshold,
    exportName: 'hyphy_dfel_evidence_plot',
    stack: [new TickStripBand({
      up: (r) => r.tally?.gc ?? 0,
      down: (r) => r.tally?.at ?? 0,
      flag: (r) => (r.tally?.surprising ?? 0) > 0,
      label: 'GC/AT'
    })],
    stackToggle: DFEL_STRIP,
    onSiteSelect: (site) => hooks.onSiteSelect(site)
  });
}

// Results table

// The rate-balance bar (α grey, β orange) with ω beside it, the ω badge in the call's hue at t
function dfelOmegaCell(r, t) {
  const call = dfelCall(r, t);
  const total = Math.max(0.1, r.alpha + r.beta);
  const alphaPct = Math.round((r.alpha / total) * 100);
  const betaPct = Math.round((r.beta / total) * 100);
  return `<td class="omega-cell"><div class="dfel-omega-cell">
              <div class="micro-balance-bar" title="α ${dfelRate(r.alpha)} vs β ${dfelRate(r.beta)}">
                <div class="bar-alpha dfel-bar-alpha" style="width: ${alphaPct}%;"></div>
                <div class="bar-beta" style="width: ${betaPct}%;"></div>
              </div>
              <span class="omega-badge ${call.badgeClass}">${dfelOmegaText(r.omega)}</span>
            </div></td>`;
}

// '3 / 1' and a two-segment greyscale bar on the result's largest GC + AT count; a muted
// dash for a site with no GC- or AT-ward change
export function dfelGcAtCell(r, extent) {
  const t = r.tally;
  if (!t) return '<td></td>';
  if (t.gc + t.at === 0) return '<td class="font-mono"><span class="dfel-muted">—</span></td>';
  const width = (n) => (extent > 0 ? (n / extent) * 100 : 0);
  return `<td class="font-mono"><div class="gcat-cell"><span>${t.gc} / ${t.at}</span>`
    + `<span class="gcat-bar" title="${t.gc} GC-ward, ${t.at} AT-ward non-synonymous changes">`
    + `<span class="gcat-gc" style="width: ${width(t.gc)}%;"></span><span class="gcat-at" style="width: ${width(t.at)}%;"></span></span></div></td>`;
}

const DFEL_CHANGES_MAX = 60;

// The reconstructed changes on one line: whole tokens (⚡ inline) up to about 60 characters,
// the synonymous count moved out as a muted suffix, an invariant site said in words
export function dfelChangesCell(r) {
  const subs = r.subs;
  if (!subs) return '<td></td>';
  if (subs.conserved) return '<td class="sub-cell"><span class="dfel-muted">Invariant (no substitutions)</span></td>';
  const tokens = subs.events.map(e => `${e.from} → ${e.to}${e.dir ? ` [${e.dir}↑]` : ''}${e.surprising ? '⚡' : ''}`);
  let shown = tokens.join(', ');
  if (shown.length > DFEL_CHANGES_MAX) {
    // Whole tokens only, as many as fit
    let kept = tokens[0];
    for (const token of tokens.slice(1)) {
      if (kept.length + 2 + token.length > DFEL_CHANGES_MAX - 2) break;
      kept += `, ${token}`;
    }
    shown = `${kept}, …`;
  }
  const syn = subs.synCount ? `<span class="dfel-muted">+${subs.synCount} syn</span>` : '';
  return `<td class="sub-cell" title="${r.substitutions}"><div class="sub-text-wrapper">${shown}${shown && syn ? ' ' : ''}${syn}</div></td>`;
}

const dfelCsv = (v) => v ?? '';

export const DFEL_COLUMNS = [
  CODON_COLUMN,
  {
    id: 'alpha',
    header: dfelLiteral('α (dS)'),
    width: 80,
    title: 'Synonymous rate; 0 at the lower bound',
    cell: (r) => `<td class="font-mono">${dfelRate(r.alpha)}</td>`,
    csv: [{ header: 'alpha', value: (r) => r.alpha }],
    sort: (r) => r.alpha
  },
  {
    id: 'beta',
    header: dfelLiteral('β (dN)'),
    width: 80,
    title: 'Non-synonymous rate; 0 at the lower bound',
    cell: (r) => `<td class="font-mono">${dfelRate(r.beta)}</td>`,
    csv: [{ header: 'beta', value: (r) => r.beta }],
    sort: (r) => r.beta
  },
  {
    id: 'omega',
    header: dfelLiteral('ω = β/α'),
    width: 130,
    title: 'Rate balance (α grey, β orange) and ω; ∞ where α is at its lower bound',
    cell: (r, { threshold }) => dfelOmegaCell(r, threshold),
    csv: [{ header: 'omega', value: (r) => (r.omega === Infinity ? 'Infinity' : dfelCsv(r.omega)) }],
    sort: (r) => r.omega ?? -1
  },
  LRT_COLUMN,
  {
    id: 'p',
    header: 'p-value',
    width: 90,
    cell: (r, { threshold }) => `<td class="font-mono pval-cell ${r.pValue <= threshold ? 'pval-sig font-bold' : ''}">${formatPValue(r.pValue)}</td>`,
    csv: [{ header: 'p-value', value: (r) => r.pValue }],
    sort: (r) => r.pValue
  },
  {
    id: 'gcat',
    header: 'GC↑ / AT↑',
    width: 110,
    title: 'Non-synonymous changes that raise / lower the codon\'s GC count; sorts by GC − AT',
    availableWhen: (caps) => caps.substitutions,
    cell: (r, { data }) => dfelGcAtCell(r, data.gcAtExtent),
    csv: [
      { header: 'GC-ward nonsyn', value: (r) => r.tally.gc },
      { header: 'AT-ward nonsyn', value: (r) => r.tally.at },
      { header: 'No net GC change nonsyn', value: (r) => r.tally.flat }
    ],
    sort: (r) => (r.tally ? r.tally.gc - r.tally.at : 0)
  },
  {
    id: 'changes',
    header: 'Reconstructed changes',
    title: 'Non-synonymous changes over the Viterbi reconstruction; ⚡ marks a surprising change',
    availableWhen: (caps) => caps.substitutions,
    cell: (r) => dfelChangesCell(r),
    csv: [
      { header: 'Surprising (⚡)', value: (r) => r.tally.surprising },
      { header: 'Synonymous', value: (r) => r.tally.syn },
      { header: 'Substitutions', value: (r) => csvQuote(r.substitutions) }
    ]
  },
  {
    id: 'detail',
    csv: [
      { header: 'alpha=beta', value: (r) => dfelCsv(r.alphaNull) },
      { header: 'Call', value: (r, { threshold }) => dfelCall(r, threshold).label }
    ]
  },
  {
    id: 'actions',
    header: 'Inspect',
    width: 80,
    align: 'right',
    cell: (r) => `<td style="text-align: right;"><button type="button" class="btn-inspect" data-action="summary" data-site="${r.site}">Inspect</button></td>`
  }
];

const dfelCalled = (cls) => (r, t) => dfelSiteClass(r, t) === cls;

export const DFEL_FILTERS = [
  {
    id: 'diversifying',
    key: 'div',
    glyph: '▲',
    label: 'Diversifying',
    unit: 'diversifying site',
    evidence: 'β > α',
    opposite: 'purifying',
    leans: (r) => r.beta > r.alpha,
    sort: { column: 'p', dir: 'asc' },
    count: (stats) => stats.diversifying,
    test: dfelCalled('diversifying'),
    subtitle: (n, t) => `Showing ${n} diversifying sites (β > α, p ≤ ${formatThreshold(t)})`
  },
  {
    id: 'purifying',
    key: 'pur',
    glyph: '▼',
    label: 'Purifying',
    unit: 'purifying site',
    evidence: 'α > β',
    opposite: 'diversifying',
    leans: (r) => r.alpha > r.beta,
    sort: { column: 'p', dir: 'asc' },
    count: (stats) => stats.purifying,
    test: dfelCalled('purifying'),
    subtitle: (n, t) => `Showing ${n} purifying sites (α > β, p ≤ ${formatThreshold(t)})`
  },
  {
    id: 'significant',
    key: 'sig',
    label: 'All Significant',
    unit: 'site',
    sort: { column: 'p', dir: 'asc' },
    count: (stats) => stats.significant,
    test: (r, t) => ['diversifying', 'purifying'].includes(dfelSiteClass(r, t)),
    subtitle: (n, t) => `Showing ${n} significant sites (p ≤ ${formatThreshold(t)})`
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

// Jump to the opposite direction when it has sites, else this direction's closest calls
export function dfelEmptyState(t, { filter, stats, rows }) {
  const title = `No ${filter.unit} reaches p ≤ ${formatThreshold(t)}`;
  const desc = filter.evidence && `No codon site has ${filter.evidence} with p ≤ ${formatThreshold(t)} under the directional model.`;
  const other = DFEL_FILTERS.find(f => f.id === filter.opposite);
  const otherCount = other ? other.count(stats) : 0;
  if (otherCount > 0) {
    return { icon: filter.glyph, title, desc, jump: { filter: other.id, label: `View ${other.label} Sites (${other.glyph} ${otherCount})` } };
  }
  const leaning = filter.leans ? rows.filter(filter.leans) : rows;
  return { icon: filter.glyph, title, desc, fallback: closestCalls(leaning, r => r.pValue) };
}

export const DFEL_TABLE = {
  title: 'Candidate Sites',
  placeholder: 'Filtered by selective regime',
  rowKey: 'site',
  rows: (data) => data.siteResults,
  columns: DFEL_COLUMNS,
  filters: DFEL_FILTERS,
  defaultFilter: 'diversifying',
  defaultSort: { column: 'codon', dir: 'asc' },
  search: {
    placeholder: 'Search codon or amino acid...',
    match: (r, q) => String(r.site).includes(q)
      || (r.subs?.events || []).some(e => e.from.toLowerCase() === q || e.to.toLowerCase() === q),
    noMatch: { title: 'No Matching Sites', desc: 'No codon sites match your current filter and search query.' }
  },
  emptyState: dfelEmptyState,
  exportName: 'hyphy_dfel_sites'
};

// Drawer: the call at the live cutoff, the rates and test, then the reconstructed changes

function dfelDrawerHeader(r, { threshold }) {
  const call = dfelCall(r, threshold);
  return { tag: `Codon #${r.site}`, tone: call.tone, badge: call.badge };
}

function dfelStatCard(label, value, sub, hue = '') {
  return `
            <div class="drawer-stat-card">
              <div class="dstat-label">${label}</div>
              <div class="dstat-value font-mono font-bold ${hue}">${value}</div>
              <div class="dstat-sub">${sub}</div>
            </div>`;
}

function dfelOmegaSub(r) {
  if (r.omega == null) return 'both rates at the lower bound';
  if (r.omega === Infinity) return 'α at its lower bound';
  return r.omega === 0 ? 'β at its lower bound' : 'β / α';
}

// The drawer's p carries one digit more than the table's, as in FEL's and MEME's drawers
function dfelRatesSection(r, t) {
  const cls = dfelSiteClass(r, t);
  const significant = cls === 'diversifying' || cls === 'purifying';
  const hue = { diversifying: 'text-orange', purifying: 'text-blue' }[cls] ?? '';
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Maximum likelihood rates</div>
          <div class="drawer-stats-grid">${[
            dfelStatCard(`Synonymous rate (${dfelLiteral('α')})`, dfelRate(r.alpha), 'dS estimate'),
            dfelStatCard(`Non-synonymous rate (${dfelLiteral('β')})`, dfelRate(r.beta), 'dN estimate'),
            dfelStatCard(`Ratio (${dfelLiteral('ω')})`, dfelOmegaText(r.omega), dfelOmegaSub(r)),
            dfelStatCard(`Null rate (${dfelLiteral('α=β')})`, r.alphaNull == null ? '—' : dfelRate(r.alphaNull), 'shared rate under the null'),
            dfelStatCard('LRT', r.lrt.toFixed(2), '2ΔlnL vs α = β'),
            dfelStatCard('p-value', formatPValue(r.pValue, 3), significant ? `${DFEL_CALLS[cls].label.toLowerCase()} at p ≤ ${formatThreshold(t)}` : `not significant at p ≤ ${formatThreshold(t)}`, hue)
          ].join('')}
          </div>
        </div>`;
}

// 'n non-synonymous: g GC-ward, a AT-ward, u no net GC change; k high-surprise (...)'
export function dfelTallySentence(t) {
  return `${t.nonsyn} non-synonymous: ${t.gc} GC-ward, ${t.at} AT-ward, ${t.flat} no net GC change; `
    + `${t.surprising} high-surprise (⚡ = −log₂ P ≥ 10 bits under the branch's transition matrix)`;
}

// Shown only for a diversifying site whose amino-acid changes mostly share one GC direction,
// when the model fitted directional regimes
export function dfelDriftNote(r, t, K) {
  const tally = r.tally;
  if (!tally || !tally.nonsyn || !(K > 1) || dfelSiteClass(r, t) !== 'diversifying') return null;
  if (Math.max(tally.gc, tally.at) / tally.nonsyn < 0.75) return null;
  return 'Most amino-acid changes here move in the same base-composition direction; the directional model accounts for regime-level drift, but inspect before attributing to adaptation.';
}

function dfelChangesSection(r, t, data) {
  if (!r.subs) return '';
  const subs = r.subs;
  const title = '<div class="drawer-section-title">Reconstructed substitutions (Viterbi, directional model)</div>';
  if (subs.conserved) {
    return `
        <div class="drawer-section">${title}
          <div class="empty-events-notice">No substitutions reconstructed (invariant site).</div>
        </div>`;
  }
  if (subs.synOnly) {
    return `
        <div class="drawer-section">${title}
          <div class="empty-events-notice">Synonymous changes only (${subs.synCount}).</div>
        </div>`;
  }
  const chips = subs.events.map(e => {
    const dir = e.dir ? `<span class="dfel-dir dfel-dir-${e.dir.toLowerCase()}">${e.dir}↑</span>` : '';
    return `<span class="dfel-chip"><span class="font-mono">${e.from} → ${e.to}</span>${dir}${e.surprising ? '<span class="dfel-flash" title="Surprising change: ≥ 10 bits">⚡</span>' : ''}</span>`;
  }).join('');
  const note = dfelDriftNote(r, t, data.K);
  return `
        <div class="drawer-section">${title}
          <div class="dfel-chips">${chips}</div>
          ${subs.synCount ? `<div class="dfel-muted dfel-syn-line">+${subs.synCount} synonymous change${subs.synCount === 1 ? '' : 's'}</div>` : ''}
          <div class="count-sentence">${dfelTallySentence(r.tally)}</div>
          <div class="empty-events-notice">Repeated changes are separate branch events. Only non-synonymous changes are listed; GC↑ / AT↑ is the sign of the codon's GC-count change.</div>${note ? `
          <div class="dfel-drift-note">${note}</div>` : ''}
        </div>`;
}

const DFEL_SUMMARY_TAB = {
  id: 'summary',
  label: 'Summary',
  comment: 'Rates, test and reconstructed substitutions',
  render: (panel, r, { data, threshold }) => {
    const missing = data.capabilities.substitutions ? '' : `
        <div class="empty-events-notice">Reconstructed substitutions are not reported in ${data.formatLabel} output.</div>`;
    panel.innerHTML = `<div class="drawer-stack">${dfelRatesSection(r, threshold)}${dfelChangesSection(r, threshold, data)}${missing}
      </div>`;
  }
};

export const DFEL_DRAWER = {
  emptyTag: 'Codon #—',
  row: siteRow,
  header: dfelDrawerHeader,
  tabs: [DFEL_SUMMARY_TAB]
};

// The directional model card under the table: the composition screen, GC by taxon, the
// rooted tree styled by regime, root support and (K > 1) the regime table

// Branch styles by regime id; with a single regime every branch is a plain ink branch
export function dfelBranchStyle(regimeId, K, colors) {
  if (!(K > 1) || regimeId == null) return {};
  if (regimeId === 2) return { stroke: colors.ink, width: 2.5 };
  if (regimeId === 0) return { stroke: colors.ink, width: 2.5, dash: [5, 3] };
  return { stroke: colors.slateLight, width: 1 };
}

// 'AT 0% · Neutral 100% · GC 0%', the soft regime weights of one branch
export function dfelWeightsText(weights) {
  return DFEL_REGIMES.map(g => `${g.short} ${Math.round((weights[g.id] ?? 0) * 100)}%`).join(' · ');
}

function dfelWeightsBar(weights) {
  return `<span class="regime-weights-bar">${DFEL_REGIMES.map(g =>
    `<span class="regime-w regime-w-${g.id}" style="width: ${(weights[g.id] ?? 0) * 100}%;"></span>`).join('')}</span>`;
}

// What the diamond marks, under the tree
export function dfelRootCaption(root) {
  if (!root?.branch) return '';
  if (root.strategy === 'keep') return '◆ input root, kept from the input tree (root_strategy = keep)';
  const w = root.weight == null ? '' : ` (w = ${rootWeightPercent(root.weight)})`;
  return `◆ ${dfelOutgroup(root) && !dfelOutgroupMissing(root) ? 'outgroup root' : 'ML root'} on ${rootEdgeLabel(root.branch)}${w}`;
}

// Notes under the root bars: the weighting, the outgroup request and how flat the profile is
// (a missing outgroup is warned about in the card header)
export function dfelRootNotes(root) {
  const notes = ['Akaike weights ∝ exp(−ΔlnL); a flat profile means the root is not identifiable from these data.'];
  const n = root.candidates.length;
  if (n) notes.push(`95% of the weight needs ${rootCredibleCount(root.candidates)} of ${n} placements.`);
  if (dfelOutgroup(root)) {
    notes.unshift(`Outgroup requested: ${root.strategy}${dfelOutgroupMissing(root) ? ' (not found)' : ''}.`);
  }
  return notes;
}

// The header warning when the requested outgroup is not in the tree
export function dfelRootWarning(root) {
  return dfelOutgroupMissing(root) ? `⚠ Outgroup ${root.strategy} not found in the tree; the ML root is used.` : null;
}

const DFEL_ROOT_TOP = 10;

// Whether there are candidate placements to rank: none under 'keep', whose root is fixed
export function dfelRootRanked(root) {
  return root.strategy !== 'keep' && root.candidates.length > 0;
}

function dfelPanelHead(title, exportId) {
  const tools = exportId ? `<span class="plot-export-group">
              <button type="button" class="plot-tool-btn" data-export="${exportId}" data-format="png" title="Export High-Res PNG">PNG</button>
              <button type="button" class="plot-tool-btn" data-export="${exportId}" data-format="svg" title="Export Vector SVG">SVG</button>
            </span>` : '';
  return `<div class="dfel-panel-head"><span class="dfel-panel-title">${title}</span>${tools}</div>`;
}

export class DirectionalModelCard {
  constructor(container) {
    this.container = container;
    this.data = null;
    this.expanded = false;
    this.showAllRoots = false;

    this.el = document.createElement('div');
    this.el.className = 'hyphy-dfel-model-card';
    this.el.innerHTML = `
      <button type="button" class="dfel-model-toggle" id="dfel-model-toggle">
        <span class="dfel-model-title">Directional model</span>
        <span class="dfel-model-summary" id="dfel-model-summary">—</span>
        <span class="dfel-model-icon" id="dfel-model-icon">▼</span>
      </button>
      <div class="dfel-model-warning" id="dfel-model-warning" style="display: none;"></div>
      <div class="dfel-model-body" id="dfel-model-body" style="display: none;">
        <div class="dfel-model-row">
          <section class="dfel-panel dfel-panel-narrow">
            ${dfelPanelHead('Composition screen')}
            <div id="dfel-composition"></div>
          </section>
          <section class="dfel-panel dfel-panel-wide">
            ${dfelPanelHead('GC by taxon', 'gc')}
            <div class="dfel-figure" id="dfel-gc"></div>
            <div class="dfel-caption" id="dfel-gc-caption"></div>
          </section>
        </div>
        <div class="dfel-model-row">
          <section class="dfel-panel dfel-panel-wide">
            ${dfelPanelHead('Rooted tree by regime', 'tree')}
            <div class="dfel-tree-host" id="dfel-tree"></div>
            <div class="dfel-caption" id="dfel-tree-hover">&nbsp;</div>
            <div class="dfel-caption" id="dfel-tree-caption"></div>
            <div id="dfel-regimes"></div>
          </section>
          <section class="dfel-panel" id="dfel-root-panel">
            ${dfelPanelHead('Root support', 'root')}
            <div class="dfel-figure" id="dfel-root"></div>
            <button type="button" class="plot-tool-btn dfel-root-more" id="dfel-root-more" style="display: none;"></button>
            <div class="dfel-caption" id="dfel-root-caption"></div>
          </section>
        </div>
      </div>
    `;
    this.container.appendChild(this.el);

    this.toggleBtn = this.el.querySelector('#dfel-model-toggle');
    this.body = this.el.querySelector('#dfel-model-body');
    this.toggleBtn.addEventListener('click', () => this.setExpanded(!this.expanded));
    this.el.querySelector('#dfel-root-more').addEventListener('click', () => {
      this.showAllRoots = !this.showAllRoots;
      this.renderRoot();
    });
    this.el.querySelectorAll('[data-export]').forEach(btn => {
      btn.addEventListener('click', () => this.export(btn.dataset.export, btn.dataset.format));
    });

    this.tree = new BranchTree(this.el.querySelector('#dfel-tree'), {
      hoverColor: 'stone',
      rootShape: 'diamond',
      exportName: 'hyphy_dfel_regime_tree',
      onHover: (name) => this.showBranch(name)
    });

    this.resizeObserver = typeof ResizeObserver === 'function'
      ? new ResizeObserver(() => this.renderFigures())
      : null;
    this.resizeObserver?.observe(this.body);
  }

  destroy() {
    this.resizeObserver?.disconnect();
    this.tree.destroy();
  }

  // Repaint the canvas tree in the current theme; the SVG figures follow the CSS tokens
  render() {
    this.tree.render();
  }

  setExpanded(expanded) {
    this.expanded = expanded;
    this.body.style.display = expanded ? '' : 'none';
    this.el.querySelector('#dfel-model-icon').textContent = expanded ? '▲' : '▼';
    if (expanded) this.renderFigures();
  }

  setData(data) {
    this.data = data;
    this.showAllRoots = false;
    const verdict = dfelModelVerdict(data);
    // The legacy export carries none of the model, which the diagnostics' format notes say
    this.el.style.display = verdict ? '' : 'none';
    if (!verdict) {
      this.setExpanded(false);
      return;
    }
    this.el.querySelector('#dfel-model-summary').textContent = [verdict.summary, dfelRootPhrase(data.root)].filter(Boolean).join(' · ');
    const warning = this.el.querySelector('#dfel-model-warning');
    warning.textContent = dfelRootWarning(data.root) ?? '';
    warning.style.display = warning.textContent ? '' : 'none';
    this.el.querySelector('#dfel-root-panel').style.display = dfelRootRanked(data.root) ? '' : 'none';

    this.renderComposition();
    this.renderTree();
    this.renderRegimes();
    this.setExpanded(data.K > 1 || data.composition.heterogeneous);
  }

  renderComposition() {
    const { composition, settings, numTaxa } = this.data;
    const thr = settings.heterogeneityThreshold;
    const df = composition.overall?.df;
    this.el.querySelector('#dfel-composition').innerHTML = `${compositionTableMarkup(composition, thr)}
          <div class="dfel-caption">• Regimes are activated (auto mode) if Overall or CP3 p &lt; ${thr ?? '—'}; CP1 and CP2 are context only.</div>
          <div class="dfel-caption">Pearson ${dfelLiteral('χ²')} homogeneity across taxa, df = 3(N−1)${df != null ? ` = ${df} for N = ${numTaxa}` : ''}.</div>`;
  }

  // The regime means the GC strip ticks, when more than one regime was fitted
  regimeMarks() {
    return this.data.K > 1
      ? this.data.regimes.filter(g => g.meanGC != null)
        .map(g => ({ value: g.meanGC, label: DFEL_REGIMES.find(r => r.id === g.id)?.short ?? g.name }))
      : [];
  }

  figureWidth(id, fallback) {
    const width = this.el.querySelector(id).clientWidth;
    return width > 0 ? width : fallback;
  }

  // The SVG figures follow their panels' widths; drawn only while the card is open
  renderFigures() {
    if (!this.data?.composition || !this.expanded) return;
    const gc = gcStripSVG(this.data.composition.taxa, { width: this.figureWidth('#dfel-gc', 360), marks: this.regimeMarks() });
    this.el.querySelector('#dfel-gc').innerHTML = gc.svg;
    const marks = this.regimeMarks();
    const regimes = marks.length ? ` Ticks above: regime mean GC (${marks.map(m => `${m.label} ${m.value.toFixed(3)}`).join(', ')}).` : '';
    this.el.querySelector('#dfel-gc-caption').textContent =
      `Filled: GC at third codon positions; hollow: GC over all positions. Dashed: mean GC3.${regimes}`;
    this.renderRoot();
  }

  // Without ranked candidates the panel is hidden and the tree caption names the root
  renderRoot() {
    const root = this.data.root;
    if (!dfelRootRanked(root)) return;
    const host = this.el.querySelector('#dfel-root');
    const more = this.el.querySelector('#dfel-root-more');
    const caption = this.el.querySelector('#dfel-root-caption');
    const n = root.candidates.length;
    const limit = this.showAllRoots ? n : DFEL_ROOT_TOP;
    host.innerHTML = rootSupportSVG(root.candidates, { chosen: root.branch, limit, width: this.figureWidth('#dfel-root', 360) }).svg;
    more.style.display = n > DFEL_ROOT_TOP ? '' : 'none';
    more.textContent = this.showAllRoots ? `Show top ${DFEL_ROOT_TOP}` : `Show all ${n}`;
    caption.textContent = dfelRootNotes(root).join(' ');
  }

  renderTree() {
    const { data } = this;
    const K = data.K;
    this.tree.setStyle((node, meta, colors) => dfelBranchStyle(data.branchRegimes[node.name]?.regimeId, K, colors));
    this.tree.setTree(data.treeNewick);
    this.tree.setLegend(K > 1 ? [
      { label: 'GC-drift', stroke: 'var(--ink)' },
      { label: 'AT-drift', stroke: 'var(--ink)', dash: true },
      { label: 'Neutral', stroke: 'var(--slate-light)' }
    ] : []);
    const branches = Object.keys(data.branchRegimes).length || Math.max(0, (this.tree.layout?.nodes.length ?? 1) - 1);
    const rootCaption = dfelRootCaption(data.root);
    const parts = [rootCaption && `${rootCaption}.`];
    if (!(K > 1)) parts.push(`Single Neutral regime: all ${branches} branches.`);
    parts.push('Branch lengths as the result\'s tree writes them (no unit).');
    this.el.querySelector('#dfel-tree-caption').textContent = parts.filter(Boolean).join(' ');
  }

  // The hovered branch's regime and its soft weights
  showBranch(name) {
    const line = this.el.querySelector('#dfel-tree-hover');
    const b = name ? this.data?.branchRegimes[name] : null;
    if (!b) {
      line.innerHTML = '&nbsp;';
      return;
    }
    const regime = this.data.regime(b.regimeId)?.name ?? DFEL_REGIMES.find(g => g.id === b.regimeId)?.name ?? '—';
    line.innerHTML = `<strong>${name}</strong> · ${regime} ${dfelWeightsBar(b.weights)} ${dfelWeightsText(b.weights)}`;
  }

  renderRegimes() {
    const { data } = this;
    const host = this.el.querySelector('#dfel-regimes');
    if (!(data.K > 1)) {
      host.innerHTML = '';
      return;
    }
    const cell = (v, digits) => (v == null ? '—' : v.toFixed(digits));
    const rows = [...data.regimes].sort((a, b) => a.id - b.id).map(g => `
              <tr><td>${g.name}</td><td class="font-mono">${cell(g.psi, 3)}</td><td class="font-mono">${cell(g.meanGC, 3)}</td><td class="font-mono">${g.branchCount == null ? '—' : Math.round(g.branchCount)}</td></tr>`).join('');
    host.innerHTML = `
          <table class="diag-fits-table dfel-regime-table">
            <thead><tr><th>Regime</th><th><span class="case-literal">ψ</span></th><th>Mean GC</th><th>Branches</th></tr></thead>
            <tbody>${rows}
            </tbody>
          </table>
          <div class="dfel-caption">ψ &lt; 0 drifts toward AT, 0 is neutral, &gt; 0 toward GC; branches is the rounded soft count.</div>`;
  }

  export(id, format) {
    if (id === 'tree') {
      if (format === 'png') this.tree.exportPNG();
      else this.tree.exportSVG();
      return;
    }
    const figure = id === 'gc'
      ? gcStripSVG(this.data.composition.taxa, { width: this.figureWidth('#dfel-gc', 360), marks: this.regimeMarks(), colors: PALETTE })
      : rootSupportSVG(this.data.root.candidates, {
        chosen: this.data.root.branch,
        limit: this.showAllRoots ? this.data.root.candidates.length : DFEL_ROOT_TOP,
        width: this.figureWidth('#dfel-root', 360),
        colors: PALETTE
      });
    const name = id === 'gc' ? 'hyphy_dfel_gc_by_taxon' : 'hyphy_dfel_root_support';
    if (format === 'png') downloadSVGAsPNG(figure.svg, figure.width, figure.height, `${name}.png`);
    else downloadSVG(figure.svg, `${name}.svg`);
  }
}

// Diagnostics: fits (never compared across the nucleotide and codon models), the run's
// settings and cross-check, the dataset, provenance and format notes

function dfelDiagItem(label, value, valueClass = '', itemClass = '') {
  return value == null || value === '' ? '' : `
      <div class="diag-item${itemClass && ` ${itemClass}`}">
        <span class="diag-lbl">${label}</span>
        <span class="diag-val ${valueClass}">${value}</span>
      </div>`;
}

function dfelBasename(path) {
  return String(path).split(/[\\/]/).pop();
}

// The stored classification at the run's cutoff, as a parse check against the live call
function dfelRunItem(data) {
  const run = data.runCheck();
  if (!run) return dfelDiagItem('Run cutoff', data.runThreshold == null ? null : `p ≤ ${formatThreshold(data.runThreshold)}`);
  const agree = run.mismatches ? `; ${run.mismatches} sites differ from the live call` : '; matches the live call';
  return dfelDiagItem(`Run cutoff, as written by ${dfelLiteral('hyphy3')}`,
    `p ≤ ${formatThreshold(run.threshold)}: ${run.written.diversifying} diversifying · ${run.written.purifying} purifying · ${run.written.neutral} neutral${agree}`);
}

export function dfelDiagnostics(data) {
  const s = data.settings;
  const versions = [data.version.analysis && `analysis ${data.version.analysis}`, data.version.software && `hyphy ${data.version.software}`]
    .filter(Boolean).join(' · ');
  const file = data.dataset.file && `<span title="${data.dataset.file}">${dfelBasename(data.dataset.file)}</span>`;
  return [
    dfelDiagItem('Analysis Method', 'dFEL (Directional Fixed Effects Likelihood)'),
    dfelDiagItem('Version', versions),
    dfelDiagItem('Output format', data.formatLabel),
    dfelDiagItem('Codon Sites', data.numSites),
    dfelDiagItem('Taxa / Sequences', data.numTaxa),
    dfelDiagItem('Site patterns', data.patterns),
    dfelDiagItem('Alignment', file, 'diag-val-wrap'),
    dfelDiagItem('Genetic code', data.dataset.geneticCode),
    dfelRunItem(data),
    dfelDiagItem('Heterogeneity threshold', s.heterogeneityThreshold),
    dfelDiagItem('Root strategy', s.rootStrategy),
    dfelDiagItem('Regimes fitted (K)', s.numRegimes),
    dfelDiagItem('Tree length', data.treeLength == null ? null : `${data.treeLength.toFixed(4)} (baseline units)`),
    dfelDiagItem('Created', data.provenance?.createdAt),
    dfelDiagItem('Method', data.description, 'diag-val-prose', 'diag-item-wide'),
    dfelDiagItem('Citation', data.citation, 'diag-val-prose', 'diag-item-wide')
  ].filter(Boolean).join('');
}

// Each fit on its own row; the root log L belongs to a third model and is hidden under
// 'keep', where the writer emits a placeholder 0
export function dfelFitsTable(data) {
  const cell = (v, digits) => (v == null ? '—' : v.toFixed(digits));
  const rows = [
    ['Nucleotide GTR', data.fits.gtr],
    ['Directional MG94 (codon)', data.fits.directional]
  ].filter(([, fit]) => fit).map(([name, fit]) => [name, fit.logL, fit.aicc]);
  if (data.root?.logL != null && data.root.strategy !== 'keep') {
    rows.push(['Directional nucleotide model at chosen root', data.root.logL, null]);
  }
  if (!rows.length) return '';
  const body = rows.map(([name, logL, aicc]) => `
              <tr><td>${name}</td><td class="font-mono">${cell(logL, 2)}</td><td class="font-mono">${cell(aicc, 2)}</td></tr>`).join('');
  return `
        <div class="diag-item">
          <span class="diag-lbl">Model fits</span>
          <table class="diag-fits-table">
            <thead><tr><th>Model</th><th>log L</th><th>${dfelLiteral('AICc')}</th></tr></thead>
            <tbody>${body}
            </tbody>
          </table>
        </div>`;
}

function dfelFormatNotes(data) {
  const notes = [
    'p is the per-site 1-df LRT of α = β under the directional model, uncorrected for multiple sites; calls and counts are recomputed at the live cutoff from p, α and β.'
  ];
  if (data.capabilities.directionalModel) {
    notes.push(
      'The models are not compared: GTR is a nucleotide model and directional MG94 a codon model, and the root log L comes from a third (directional nucleotide) model.',
      'K is the number of regimes fitted. The composition screen picks it only in auto mode; --regimes 1 or 3 forces it, which the model card marks as forced.',
      'The tree is drawn with the lengths its Newick text carries, rescaled from the baseline fit, so its total differs from the tree length above.'
    );
  } else {
    notes.push(`The composition screen, rooting, regimes and reconstructed substitutions are not reported in ${data.formatLabel} output.`);
  }
  return `<ul class="diag-notes">${notes.map(n => `<li>${n}</li>`).join('')}</ul>`;
}

export function dfelDiagnosticsPanel(el, data) {
  el.insertAdjacentHTML('beforeend', `${dfelFitsTable(data)}
        <div class="diag-item">
          <span class="diag-lbl">Format notes</span>
          ${dfelFormatNotes(data)}
        </div>`);
}

export const DFEL_CONFIG = {
  id: 'dfel',
  Parser: DFELData,
  controller: DFEL_P_THRESHOLD,
  row: siteRow,
  defaultFilter: DFEL_TABLE.defaultFilter,
  formatLabel: (data) => data.formatLabel,
  kpis: (parent, hooks) => new KpiBar(parent, {
    cards: DFEL_KPI_CARDS,
    controller: DFEL_P_THRESHOLD,
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    onThresholdChange: hooks.onThresholdChange,
    onStatisticChange: hooks.onStatisticChange
  }),
  hero: dfelHero,
  table: (parent, hooks) => new ResultsTable(parent, Object.assign({}, DFEL_TABLE, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    filterType: hooks.filterType,
    onSelect: (site, row, action) => hooks.onSiteSelect(site, action)
  })),
  secondary: (parent) => new DirectionalModelCard(parent),
  drawer: (container, hooks) => new Drawer(container, Object.assign({}, DFEL_DRAWER, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    siblings: hooks.siblings,
    onStep: hooks.onStep,
    onClose: hooks.onClose
  })),
  diagnosticsTitle: 'Model fits &amp; run details',
  diagnostics: dfelDiagnostics,
  newick: (data) => data.treeNewick,
  diagnosticsPanel: dfelDiagnosticsPanel
};
