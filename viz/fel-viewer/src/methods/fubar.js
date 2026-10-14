// FUBAR method config: pervasive selection per codon from posterior probabilities over
// an alpha x beta rate grid. A site is called when its posterior reaches the cutoff
// (P ≥ t, so a higher t is stricter): diversifying when P[β>α] ≥ t, otherwise purifying
// when P[α>β] ≥ t. The purifying class is a hyphy3 convention, applied to every format.

import { FUBARData, fubarSiteClass } from '../data/fubar-data.js';
import { ALIGNMENT_KPI, CODON_COLUMN, siteRow } from './fel.js';
import { KpiBar } from '../components/kpi-bar.js';
import { PosteriorPlot } from '../components/posterior-plot.js';
import { evidenceTipRow } from '../components/evidence-plot.js';
import { ResultsTable, closestCalls } from '../components/results-table.js';
import { Drawer } from '../components/drawer.js';
import { RateGridHeatmap } from '../components/rate-grid-heatmap.js';
import { formatPosteriorCutoff, formatPosterior, formatEBF } from '../core/format.js';

export const FUBAR_POSTERIOR_CUTOFF = {
  title: 'Posterior Cutoff',
  statistic: 'P',
  direction: 'ge',
  min: 0.5,
  max: 0.999,
  step: 0.001,
  presets: [0.9, 0.95, 0.99],
  seed: 0.9,
  format: formatPosteriorCutoff
};

const FUBAR_CALLS = {
  diversifying: { tone: 'div', glyph: '▲', label: 'Diversifying', badge: '▲ Diversifying selection', badgeClass: 'omega-div' },
  purifying: { tone: 'pur', glyph: '▼', label: 'Purifying', badge: '▼ Purifying selection', badgeClass: 'omega-pur' },
  neutral: { tone: 'neutral', glyph: '—', label: 'Neutral', badge: 'Neutral', badgeClass: '' }
};

function fubarCall(r, t) {
  return FUBAR_CALLS[fubarSiteClass(r, t)];
}

// Greek and bracketed notation inside uppercased labels: 'α' must not read as 'Α'
function fubarLiteral(text) {
  return `<span class="case-literal">${text}</span>`;
}

function fubarPercent(p) {
  return `${(p * 100).toFixed(1)}%`;
}

export const FUBAR_KPI_CARDS = [
  ALIGNMENT_KPI,
  {
    id: 'div',
    label: 'Diversifying sites',
    value: (stats) => `▲ ${stats.diversifying}`,
    sub: (stats, { threshold }) => `P[β>α] ≥ ${formatPosteriorCutoff(threshold)}`,
    tone: 'accent'
  },
  {
    id: 'pur',
    label: 'Purifying sites',
    value: (stats) => `▼ ${stats.purifying}`,
    sub: (stats, { threshold }) => `P[α>β] ≥ ${formatPosteriorCutoff(threshold)}`,
    tone: 'purifying'
  },
  {
    id: 'efp',
    label: 'Expected false positives',
    value: (stats) => stats.expectedFalse.toFixed(2),
    sub: (stats) => `95% CI ${stats.falseCI[0]}–${stats.falseCI[1]} of ${stats.diversifying} diversifying calls`,
    tone: 'ink'
  }
];

// Posterior mass split three ways: P[α>β] | P[α=β] | P[β>α]
export function posteriorMassBar(r) {
  return `<span class="mass-bar" title="P[α>β] ${fubarPercent(r.probNeg)} · P[α=β] ${fubarPercent(r.probEq)} · P[β>α] ${fubarPercent(r.probPos)}">`
    + `<span class="mass-neg" style="width: ${r.probNeg * 100}%;"></span>`
    + `<span class="mass-eq" style="width: ${r.probEq * 100}%;"></span>`
    + `<span class="mass-pos" style="width: ${r.probPos * 100}%;"></span></span>`;
}

// A posterior cell: bold in the call's hue when it passes the cutoff
function posteriorCell(p, t, hue) {
  return p >= t ? `font-bold ${hue}` : '';
}

export const FUBAR_COLUMNS = [
  CODON_COLUMN,
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
    header: 'Class',
    width: 110,
    cell: (r, { threshold }) => {
      const call = fubarCall(r, threshold);
      return `<td><span class="omega-badge ${call.badgeClass}">${call.glyph} ${call.label}</span></td>`;
    },
    csv: [{ header: 'Class', value: (r, { threshold }) => fubarCall(r, threshold).label }],
    sort: (r) => r.probPos - r.probNeg
  },
  {
    id: 'alpha',
    header: fubarLiteral('α (dS)'),
    width: 70,
    title: 'Posterior mean synonymous rate',
    cell: (r) => `<td class="font-mono">${r.alpha.toFixed(2)}</td>`,
    csv: [{ header: 'alpha', value: (r) => r.alpha }],
    sort: (r) => r.alpha
  },
  {
    id: 'beta',
    header: fubarLiteral('β (dN)'),
    width: 70,
    title: 'Posterior mean non-synonymous rate',
    cell: (r) => `<td class="font-mono">${r.beta.toFixed(2)}</td>`,
    csv: [{ header: 'beta', value: (r) => r.beta }],
    sort: (r) => r.beta
  },
  {
    id: 'diff',
    header: fubarLiteral('β − α'),
    width: 70,
    title: 'Posterior mean β minus posterior mean α',
    cell: (r) => `<td class="font-mono">${r.diff >= 0 ? '+' : '−'}${Math.abs(r.diff).toFixed(2)}</td>`,
    csv: [{ header: 'beta-alpha', value: (r) => r.diff }],
    sort: (r) => r.diff
  },
  {
    id: 'ppos',
    header: fubarLiteral('P[β>α]'),
    width: 130,
    cell: (r, { threshold }) => `<td class="font-mono"><div class="posterior-cell"><span class="${posteriorCell(r.probPos, threshold, 'text-orange')}">${formatPosterior(r.probPos)}</span>${posteriorMassBar(r)}</div></td>`,
    csv: [{ header: 'P[beta>alpha]', value: (r) => r.probPos }],
    sort: (r) => r.probPos
  },
  {
    id: 'pneg',
    header: fubarLiteral('P[α>β]'),
    width: 80,
    cell: (r, { threshold }) => `<td class="font-mono ${posteriorCell(r.probNeg, threshold, 'text-blue')}">${formatPosterior(r.probNeg)}</td>`,
    csv: [{ header: 'P[alpha>beta]', value: (r) => r.probNeg }],
    sort: (r) => r.probNeg
  },
  {
    id: 'ebf',
    header: fubarLiteral('EBF[β>α]'),
    width: 80,
    title: 'Empirical Bayes factor: posterior odds ÷ prior odds of β > α (informational; the cutoff applies to the posterior)',
    cell: (r) => `<td class="font-mono">${formatEBF(r.ebf)}</td>`,
    csv: [{ header: 'EBF[beta>alpha]', value: (r) => r.ebf }],
    sort: (r) => r.ebf
  },
  {
    id: 'actions',
    header: '',
    width: 80,
    align: 'right',
    cell: (r) => `<td style="text-align: right;"><button type="button" class="btn-inspect" data-action="posterior" data-site="${r.site}">Inspect</button></td>`
  }
];

export const FUBAR_FILTERS = [
  {
    id: 'diversifying',
    key: 'div',
    glyph: '▲',
    label: 'Diversifying',
    unit: 'diversifying site',
    statistic: 'P[β>α]',
    posterior: (r) => r.probPos,
    opposite: 'purifying',
    count: (stats) => stats.diversifying,
    test: (r, t) => fubarSiteClass(r, t) === 'diversifying',
    subtitle: (n, t) => `Showing ${n} diversifying sites (P[β>α] ≥ ${formatPosteriorCutoff(t)})`
  },
  {
    id: 'purifying',
    key: 'pur',
    glyph: '▼',
    label: 'Purifying',
    unit: 'purifying site',
    statistic: 'P[α>β]',
    posterior: (r) => r.probNeg,
    opposite: 'diversifying',
    count: (stats) => stats.purifying,
    test: (r, t) => fubarSiteClass(r, t) === 'purifying',
    subtitle: (n, t) => `Showing ${n} purifying sites (P[α>β] ≥ ${formatPosteriorCutoff(t)})`
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

// Jump to the mirror class when it has sites, else show the sites closest to the cutoff
export function fubarEmptyState(t, { filter, stats, rows }) {
  const title = `No ${filter.unit} reaches ${fubarLiteral(filter.statistic)} ≥ ${formatPosteriorCutoff(t)}`;
  const other = FUBAR_FILTERS.find(f => f.id === filter.opposite);
  const otherCount = other.count(stats);
  if (otherCount > 0) {
    return { icon: filter.glyph, title, jump: { filter: other.id, label: `View ${other.label} Sites (${other.glyph} ${otherCount})` } };
  }
  return { icon: filter.glyph, title, fallback: closestCalls(rows, r => -filter.posterior(r)) };
}

export const FUBAR_TABLE = {
  title: 'Candidate Sites',
  placeholder: 'Filtered by posterior call',
  rowKey: 'site',
  rows: (data) => data.siteResults,
  columns: FUBAR_COLUMNS,
  filters: FUBAR_FILTERS,
  defaultFilter: 'diversifying',
  defaultSort: { column: 'ppos', dir: 'desc' },
  search: {
    placeholder: 'Search codon...',
    match: (r, q) => String(r.site).includes(q),
    noMatch: { title: 'No Matching Sites', desc: 'No codon sites match your current filter and search query.' }
  },
  emptyState: fubarEmptyState,
  exportName: 'hyphy_fubar_sites'
};

function fubarTooltip(r, track) {
  const call = fubarCall(r, track.threshold);
  const rows = [
    evidenceTipRow('α (dS)', r.alpha.toFixed(2)),
    evidenceTipRow('β (dN)', r.beta.toFixed(2)),
    evidenceTipRow('P[β>α]', formatPosterior(r.probPos)),
    evidenceTipRow('P[α>β]', formatPosterior(r.probNeg)),
    evidenceTipRow('EBF[β>α]', formatEBF(r.ebf))
  ];
  return `
      <div class="tip-header">
        <span class="tip-title">Codon #${r.site}</span>
        <span class="tip-badge tip-${call.tone}">${call.tone === 'neutral' ? 'Neutral' : call.badge}</span>
      </div>
      <div class="tip-grid">${rows.join('')}</div>
      <div class="tip-hint">Click point to inspect the site posterior</div>
    `;
}

// '20×20 grid, Variational-Bayes', from what the result says it ran with
function fubarRunNote(data) {
  const { gridPoints, method } = data.settings;
  return [gridPoints && `${gridPoints}×${gridPoints} grid`, method].filter(Boolean).join(', ');
}

export const FUBAR_EVIDENCE = {
  tone: (r, track) => {
    const tone = fubarCall(r, track.threshold).tone;
    return tone === 'neutral' ? null : tone;
  },
  pos: (r) => r.probPos,
  neg: (r) => r.probNeg,
  diff: (r) => r.diff,
  tooltip: fubarTooltip,
  ruleLabel: (t) => `P ≥ ${formatPosteriorCutoff(t)}`,
  labels: { pos: 'P[β>α]', signed: ['P[β>α]', 'P[α>β]'], diff: 'β − α' },
  subtitles: {
    pos: (data) => `Posterior probability of β > α per codon (FUBAR, ${fubarRunNote(data)})`,
    signed: (data) => `+P[β>α] or −P[α>β], whichever is larger (FUBAR, ${fubarRunNote(data)})`,
    diff: () => 'Posterior mean β − α per codon'
  }
};

// Drawer: the call at the live cutoff, the posterior mass split, the posterior
// summaries and where the site's posterior sits on the rate grid

function fubarDrawerHeader(r, { threshold }) {
  const call = fubarCall(r, threshold);
  return { tag: `Codon #${r.site}`, tone: call.tone, badge: call.badge };
}

function fubarMassSection(r, t) {
  const label = (cls, text, p) => `<span class="mass-label ${cls}">${text} <strong>${fubarPercent(p)}</strong></span>`;
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Posterior mass</div>
          <div class="mass-bar mass-bar-full">
            <span class="mass-neg" style="width: ${r.probNeg * 100}%;"></span>
            <span class="mass-eq" style="width: ${r.probEq * 100}%;"></span>
            <span class="mass-pos" style="width: ${r.probPos * 100}%;"></span>
            <span class="mass-tick mass-tick-neg" style="left: ${t * 100}%;" title="P[α>β] ≥ ${formatPosteriorCutoff(t)} past here"></span>
            <span class="mass-tick mass-tick-pos" style="left: ${(1 - t) * 100}%;" title="P[β>α] ≥ ${formatPosteriorCutoff(t)} past here"></span>
          </div>
          <div class="mass-labels">
            ${label('text-blue', 'P[α>β]', r.probNeg)}
            ${label('', 'P[α=β]', r.probEq)}
            ${label('text-orange', 'P[β>α]', r.probPos)}
          </div>
          <div class="empty-events-notice">Ticks mark the cutoff ${formatPosteriorCutoff(t)} from each end; a site is called when its colour reaches the tick.</div>
        </div>`;
}

function fubarStatCard(label, value, sub, hue = '') {
  return `
            <div class="drawer-stat-card">
              <div class="dstat-label">${fubarLiteral(label)}</div>
              <div class="dstat-value font-mono font-bold ${hue}">${value}</div>
              <div class="dstat-sub">${sub}</div>
            </div>`;
}

function fubarStatsSection(r, t, data) {
  const prior = data.priorPositive;
  const priorNote = prior == null ? '' : ` (prior P[β>α] = ${fubarPercent(prior)})`;
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Posterior summaries</div>
          <div class="drawer-stats-grid">${[
            fubarStatCard('Synonymous rate (α)', r.alpha.toFixed(2), 'posterior mean'),
            fubarStatCard('Non-synonymous rate (β)', r.beta.toFixed(2), 'posterior mean'),
            fubarStatCard('β − α', `${r.diff >= 0 ? '+' : '−'}${Math.abs(r.diff).toFixed(2)}`, 'posterior mean β minus posterior mean α'),
            fubarStatCard('P[β>α]', formatPosterior(r.probPos), r.probPos >= t ? 'diversifying at this cutoff' : 'below the cutoff', r.probPos >= t ? 'text-orange' : ''),
            fubarStatCard('P[α>β]', formatPosterior(r.probNeg), r.probNeg >= t ? 'purifying at this cutoff' : 'below the cutoff', r.probNeg >= t ? 'text-blue' : ''),
            fubarStatCard('EBF[β>α]', formatEBF(r.ebf), `posterior odds ÷ prior odds${priorNote}`)
          ].join('')}
          </div>
        </div>`;
}

const FUBAR_POSTERIOR_TAB = {
  id: 'posterior',
  label: 'Site posterior',
  comment: 'Posterior mass, summaries and rate grid',
  render: (panel, r, { data, threshold }) => {
    const weights = data.sitePosterior(r.site);
    panel.innerHTML = `<div class="drawer-stack">${fubarMassSection(r, threshold)}${fubarStatsSection(r, threshold, data)}
        <div class="drawer-section">
          <div class="drawer-section-title">Where this site sits on the rate grid</div>
          <div class="rate-grid-pair"></div>
          ${weights ? '' : '<div class="empty-events-notice">Per-site posterior grid not included in this file.</div>'}
        </div>
      </div>`;
    if (!data.grid) return;
    const pair = panel.querySelector('.rate-grid-pair');
    if (weights) {
      new RateGridHeatmap(pair, {
        title: `Codon #${r.site} posterior`,
        caption: 'Mode outlined; ✕ marks the posterior mean (ᾱ, β̄)',
        size: 160,
        grid: data.grid,
        weights,
        mode: true,
        mean: { alpha: r.alpha, beta: r.beta },
        exportName: `hyphy_fubar_site${r.site}_posterior`
      });
    }
    new RateGridHeatmap(pair, {
      title: 'Estimated grid weights',
      caption: 'The empirical prior every site shares',
      size: 170,
      grid: data.grid,
      weights: data.grid.weight,
      exportName: 'hyphy_fubar_grid_weights'
    });
  }
};

export const FUBAR_DRAWER = {
  emptyTag: 'Codon #—',
  row: siteRow,
  header: fubarDrawerHeader,
  tabs: [FUBAR_POSTERIOR_TAB]
};

// Diagnostics: the run, the grid and the nucleotide GTR fit the grid is built on

// valueClass 'diag-val-wrap' lets long unbroken text (paths) wrap inside its column
function fubarDiagItem(label, value, valueClass = '') {
  return value == null || value === '' ? '' : `
      <div class="diag-item">
        <span class="diag-lbl">${label}</span>
        <span class="diag-val ${valueClass}">${value}</span>
      </div>`;
}

function fubarSeconds(s) {
  if (s == null) return null;
  return s < 10 ? `${s.toFixed(2)} s` : `${Math.round(s)} s`;
}

// Chain settings matter only to the sampling methods
function fubarChainItem(settings) {
  if (!settings.method || settings.method === 'Variational-Bayes' || settings.chains == null) return '';
  const n = (v) => Number(v).toLocaleString();
  return fubarDiagItem('Chains', `${settings.chains} × ${n(settings.chainLength)} steps, burn-in ${n(settings.burnIn)}, ${n(settings.samples)} samples`);
}

export function fubarDiagnostics(data) {
  const s = data.settings;
  const gtr = data.gtr || {};
  const run = data.runSummary;
  const versions = data.format === 'v3'
    ? `analysis ${data.version.analysis} · hyphy ${data.version.software}`
    : `analysis ${data.version.analysis} · ${data.format === 'hyphy25' ? 'HyPhy' : 'hyphy3'} ${data.version.software}`;
  const citation = data.citation && (data.citation.doi
    ? `${data.citation.text} doi:${data.citation.doi}`
    : data.citation.text);
  return [
    fubarDiagItem('Analysis Method', 'FUBAR (Fast, Unconstrained Bayesian AppRoximation)'),
    fubarDiagItem('Version', versions),
    fubarDiagItem('Output format', data.formatLabel),
    fubarDiagItem('Codon Sites', data.numSites),
    fubarDiagItem('Taxa / Sequences', data.numTaxa),
    fubarDiagItem('Posterior estimation', s.method),
    fubarChainItem(s),
    fubarDiagItem('Rate grid', s.gridPoints && `${s.gridPoints} × ${s.gridPoints} (${s.totalGridPoints} points)`),
    fubarDiagItem('Dirichlet concentration', s.concentration),
    fubarDiagItem(fubarLiteral('Non-zero α constraint'), s.nonZero == null ? null : (s.nonZero ? 'yes' : 'no')),
    fubarDiagItem('Grid scalers', s.bestScalerAlpha == null ? null : `α ${s.bestScalerAlpha.toFixed(3)} · β ${s.bestScalerBeta.toFixed(3)}`),
    fubarDiagItem('Tree scaling', s.treeScaling),
    fubarDiagItem('Run cutoff', data.runThreshold == null ? null : `P ≥ ${formatPosteriorCutoff(data.runThreshold)}`),
    run && fubarDiagItem(`As written by hyphy3 at ${fubarLiteral('P')} ≥ ${formatPosteriorCutoff(run.threshold)}`,
      `${run.diversifying} diversifying · ${run.purifying} purifying · EFP ${run.expected.toFixed(2)} (95% CI ${run.ci[0]}–${run.ci[1]})`),
    fubarDiagItem(fubarLiteral('Prior P[β>α]'), data.priorPositive == null ? null : fubarPercent(data.priorPositive)),
    fubarDiagItem('GTR log L', gtr.logL == null ? null : gtr.logL.toFixed(3)),
    fubarDiagItem(`GTR AICc (${data.formatLabel})`, gtr.aicc == null ? null : gtr.aicc.toFixed(2)),
    fubarDiagItem('GTR parameters', gtr.params),
    fubarDiagItem('GTR relative rates', gtr.rates && gtr.rates.map(([k, v]) => `${k} ${Number(v.toFixed(3))}`).join(' · ')),
    fubarDiagItem('GTR tree length', gtr.treeLength == null ? null : `${gtr.treeLength.toFixed(3)} subs/site`),
    fubarDiagItem('Run time', fubarSeconds(data.timers.overall)),
    fubarDiagItem('Command', data.timers.command, 'diag-val-wrap'),
    fubarDiagItem('Citation', citation)
  ].filter(Boolean).join('');
}

// Notes on what each format means, so no number is read across formats
function fubarFormatNotes(data) {
  const notes = [
    'The purifying class (P[α>β] ≥ t, when not diversifying) is a hyphy3 convention, applied here to every format; HyPhy 2.5 reports positive sites only.',
    'Calls, counts and expected false positives are recomputed at the live cutoff; P[α=β] is the grid mass on the diagonal, so P[β>α] + P[α>β] < 1.',
    'AICc is shown as this file\'s producer wrote it and is not comparable across formats.',
    'GTR tree length sums the fitted GTR branch lengths, not the input tree.'
  ];
  if (data.format !== 'v3') {
    notes.push('HyPhy 2.5 prints the expected false positives to the console only; the viewer computes them from the site posteriors.');
  }
  return `<ul class="diag-notes">${notes.map(n => `<li>${n}</li>`).join('')}</ul>`;
}

export function fubarDiagnosticsPanel(el, data) {
  if (data.grid) {
    const caption = data.priorPositive == null
      ? 'Estimated grid weights (empirical prior)'
      : `Estimated grid weights (empirical prior). Mass with β > α: ${fubarPercent(data.priorPositive)}`;
    new RateGridHeatmap(el, {
      title: 'Global rate grid',
      caption,
      size: 320,
      grid: data.grid,
      weights: data.grid.weight,
      exportName: 'hyphy_fubar_rate_grid'
    });
  }
  el.insertAdjacentHTML('beforeend', `
        <div class="diag-item">
          <span class="diag-lbl">Format notes</span>
          ${fubarFormatNotes(data)}
        </div>`);
}

export const FUBAR_CONFIG = {
  id: 'fubar',
  Parser: FUBARData,
  controller: FUBAR_POSTERIOR_CUTOFF,
  row: siteRow,
  defaultFilter: FUBAR_TABLE.defaultFilter,
  formatLabel: (data) => data.formatLabel,
  kpis: (parent, hooks) => new KpiBar(parent, {
    cards: FUBAR_KPI_CARDS,
    controller: FUBAR_POSTERIOR_CUTOFF,
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    onThresholdChange: hooks.onThresholdChange,
    onStatisticChange: hooks.onStatisticChange
  }),
  hero: (parent, hooks) => new PosteriorPlot(parent, {
    preset: FUBAR_EVIDENCE,
    threshold: hooks.threshold,
    exportName: 'hyphy_fubar_posterior_plot',
    onSiteSelect: (site) => hooks.onSiteSelect(site)
  }),
  table: (parent, hooks) => new ResultsTable(parent, Object.assign({}, FUBAR_TABLE, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    filterType: hooks.filterType,
    onSelect: (site, row, action) => hooks.onSiteSelect(site, action)
  })),
  drawer: (container, hooks) => new Drawer(container, Object.assign({}, FUBAR_DRAWER, {
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    siblings: hooks.siblings,
    onStep: hooks.onStep,
    onClose: hooks.onClose
  })),
  diagnosticsTitle: 'Model fits, grid &amp; run details',
  diagnostics: fubarDiagnostics,
  diagnosticsPanel: fubarDiagnosticsPanel
};
