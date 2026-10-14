// MEME method config: episodic diversifying selection (beta+ vs alpha) per codon
// Shares FEL's site components and parser; the KPI cards, table, drawer and diagnostics differ here

import { FEL_CONFIG, mg94FitItems, SITE_P_THRESHOLD, FEL_KPI_CARDS, ALIGNMENT_KPI, sitePercent } from './fel.js';
import { CODON_COLUMN, LRT_COLUMN, SUBSTITUTIONS_COLUMN, SITE_ACTIONS_COLUMN, FEL_COLUMNS, FEL_FILTERS, FEL_TABLE, siteResultsTable } from './fel.js';
import { FEL_DRAWER, SITE_TREE_TAB, siteOverviewTab, siteDrawer } from './fel.js';
import { KpiBar } from '../components/kpi-bar.js';
import { closestCalls } from '../components/results-table.js';
import { PALETTE, withAlpha } from '../data/palette.js';
import { formatThreshold, formatPValue, formatOmegaInf, formatEBF } from '../core/format.js';

// The bar draws FEL's blank cards until a result loads, then relabels them
export const MEME_KPI_CARDS = [
  ALIGNMENT_KPI,
  {
    id: 'omega',
    label: 'Analysis Model',
    value: () => 'MEME',
    sub: () => 'Episodic Diversifying',
    tone: 'ink'
  },
  {
    id: 'div',
    label: 'Episodically Diversifying',
    value: (stats) => `▲ ${stats.diversifying}`,
    sub: (stats, { threshold }) => `${sitePercent(stats.diversifying, stats)}% of sites (p ≤ ${formatThreshold(threshold)})`,
    tone: 'accent'
  },
  {
    id: 'pur',
    label: 'Neutral / Pervasive',
    value: (stats) => `${stats.neutral}`,
    sub: () => 'Conserved or Neutral sites',
    tone: 'ink'
  }
];

// Branches with EBF >= 100 at a site, as the parser counted them or from the EBF map
function memeBranchCount(r) {
  return r.branchesUnderSelection || (r.branchEBF ? Object.values(r.branchEBF).filter(v => v >= 100).length : 0);
}

const memeBetaPlus = (r) => r.beta_plus ?? r.beta;

export const MEME_COLUMNS = [
  CODON_COLUMN,
  {
    id: 'rates',
    header: 'Rates (α : β⁺)',
    width: 140,
    cell: (r) => `<td class="font-mono rates-cell">
              <span class="rate-alpha font-mono" title="Synonymous rate (α)">α=${r.alpha.toFixed(2)}</span>
              <span class="rate-sep">|</span>
              <span class="rate-beta font-mono" title="Episodic non-synonymous rate (β⁺)">β⁺=${memeBetaPlus(r).toFixed(2)}</span>
            </td>`,
    csv: [
      { header: 'Alpha (dS)', value: (r) => r.alpha },
      { header: 'Beta+ (dN+)', value: memeBetaPlus }
    ],
    sort: (r) => memeBetaPlus(r) - r.alpha
  },
  {
    id: 'weight',
    header: 'Weight (p⁺)',
    width: 140,
    cell: (r) => {
      const pPlusPct = Math.round((r.p_plus ?? 0) * 100);
      return `<td class="balance-cell">
              <div class="meme-weight-bar-wrapper" title="Weight allocated to episodic selection: ${pPlusPct}%">
                <div class="micro-balance-bar">
                  <div class="bar-meme-weight" style="width: ${pPlusPct}%;"></div>
                </div>
                <span class="meme-weight-label">${pPlusPct}%</span>
              </div>
            </td>`;
    },
    csv: [{ header: 'Weight p+', value: (r) => r.p_plus ?? 0 }],
    sort: (r) => r.p_plus ?? 0
  },
  {
    id: 'branches',
    header: 'Branches (EBF≥100)',
    width: 110,
    cell: (r) => {
      const count = memeBranchCount(r);
      const badge = count > 0
        ? `<span class="ebf-branch-badge ebf-active" title="${count} branch(es) with Empirical Bayes Factor ≥ 100">⚡ ${count}</span>`
        : `<span class="font-mono" style="color: ${PALETTE.neutral}; font-size: 11px;">0</span>`;
      return `<td class="branches-cell">
              ${badge}
            </td>`;
    },
    csv: [{ header: 'Branches (EBF>=100)', value: memeBranchCount }],
    sort: memeBranchCount
  },
  LRT_COLUMN,
  {
    id: 'p',
    header: 'p-value',
    width: 90,
    cell: (r, { threshold }) => {
      const isSig = r.pValue <= threshold;
      return `<td class="font-mono pval-cell ${isSig ? 'pval-sig font-bold' : ''}" style="${isSig ? `color: ${PALETTE.accent};` : ''}">${formatPValue(r.pValue)}</td>`;
    },
    csv: [{ header: 'p-value', value: (r) => r.pValue }],
    sort: (r) => r.pValue
  },
  {
    id: 'regime',
    csv: [{ header: 'Regime', value: (r, { threshold }) => r.pValue <= threshold ? 'Episodic Diversifying' : 'Neutral' }]
  },
  SUBSTITUTIONS_COLUMN,
  SITE_ACTIONS_COLUMN
];

// MEME tests one direction, so FEL's purifying and significant pills are hidden
export const MEME_FILTERS = [
  Object.assign({}, FEL_FILTERS[0], {
    label: 'Episodic Diversifying',
    test: (r, t) => r.pValue <= t,
    subtitle: (n, t) => `Showing ${n} episodically diversifying sites (p ≤ ${formatThreshold(t)})`
  }),
  FEL_FILTERS[3]
];

function memeEmptyState(t, { rows }) {
  return {
    icon: '▲',
    title: `No episodically diversifying site reaches p ≤ ${formatThreshold(t)}`,
    desc: 'No codon site showed evidence of episodic diversifying selection at the current threshold.',
    fallback: closestCalls(rows, r => r.pValue)
  };
}

// The table draws FEL's blank columns and pills until a result loads
export const MEME_TABLE = Object.assign({}, FEL_TABLE, {
  columns: MEME_COLUMNS,
  filters: MEME_FILTERS,
  blank: { columns: FEL_COLUMNS, filters: FEL_FILTERS },
  emptyState: memeEmptyState,
  exportName: 'hyphy_meme_candidates'
});

function memeDrawerHeader(r, { threshold }) {
  return r.pValue <= threshold
    ? { tag: `Codon #${r.site}`, tone: 'div', badge: '▲ Episodic Diversifying Selection' }
    : { tag: `Codon #${r.site}`, tone: 'neutral', badge: 'Neutral / Constrained' };
}

function memeRatesGrid(r, threshold) {
  const isSig = r.pValue <= threshold;
  const omegaStr = formatOmegaInf(r.omega);
  const pStr = formatPValue(r.pValue, 3);
  const pPlusPct = Math.round((r.p_plus ?? 0) * 100);
  const p1Pct = Math.round((r.p1 ?? (1 - (r.p_plus ?? 0))) * 100);
  return `
        <div class="drawer-section">
          <div class="drawer-section-title">Maximum Likelihood Rates (MEME Mixture Model)</div>
          <div class="drawer-stats-grid">
            <div class="drawer-stat-card">
              <div class="dstat-label">Synonymous Rate (α)</div>
              <div class="dstat-value font-mono font-bold">${r.alpha.toFixed(2)}</div>
              <div class="dstat-sub">dS estimate</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Constrained Rate (β₁)</div>
              <div class="dstat-value font-mono font-bold">${(r.beta1 ?? 0).toFixed(2)}</div>
              <div class="dstat-sub">p₁ = ${p1Pct}% of branches</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Episodic Rate (β⁺)</div>
              <div class="dstat-value font-mono font-bold text-orange">${memeBetaPlus(r).toFixed(2)}</div>
              <div class="dstat-sub">p⁺ = ${pPlusPct}% under selection</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Ratio (β⁺ / α)</div>
              <div class="dstat-value font-mono font-bold">${omegaStr}</div>
              <div class="dstat-sub">${isSig ? 'Episodic Diversification' : 'Unconstrained'}</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">LRT Statistic</div>
              <div class="dstat-value font-mono font-bold">${r.lrt.toFixed(2)}</div>
              <div class="dstat-sub">2ΔlnL (vs β⁺ ≤ α)</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">p-value</div>
              <div class="dstat-value font-mono font-bold ${isSig ? 'text-orange' : ''}">${pStr}</div>
              <div class="dstat-sub">${isSig ? 'Episodically Diversifying' : 'Not Significant'}</div>
            </div>
          </div>
        </div>
      `;
}

// Branches with EBF >= 100 at the site, matched to their reconstructed substitution
function memeBranchSection(r, { events, rootCodon, rootAA }) {
  const supportedBranches = Object.entries(r.branchEBF || {})
    .filter(([_, ebf]) => ebf >= 100)
    .sort((a, b) => b[1] - a[1]);

  if (supportedBranches.length === 0) {
    return `
          <div class="drawer-section">
            <div class="drawer-section-title">Branches Under Episodic Selection</div>
            <div class="empty-events-notice">
              No individual branches reached the strict EBF ≥ 100 evidentiary cutoff at this codon site.
              ${r.p_plus > 0 ? `The episodic signal is distributed across branches (overall p⁺ = ${(r.p_plus * 100).toFixed(1)}%).` : ''}
            </div>
          </div>
        `;
  }

  return `
          <div class="drawer-section">
            <div class="drawer-section-title">
              <span>Branches Under Episodic Selection (EBF ≥ 100)</span>
              <span class="root-badge" style="background: ${withAlpha(PALETTE.accent, 0.15)}; color: ${PALETTE.accent}; border-color: ${withAlpha(PALETTE.accent, 0.3)};">
                <strong>${supportedBranches.length}</strong> branch(es) supported
              </span>
            </div>
            <div class="events-table-wrapper">
              <table class="events-table">
                <thead>
                  <tr>
                    <th>Branch</th>
                    <th>Empirical Bayes Factor (EBF)</th>
                    <th>Reconstructed Substitution</th>
                    <th>Type</th>
                  </tr>
                </thead>
                <tbody>
                  ${supportedBranches.map(([bName, ebf]) => {
                    const ev = events.find(e => e.branch === bName);
                    const subStr = ev
                      ? `<strong>${ev.parent_codon || rootCodon}</strong> (${ev.parent_aa || rootAA}) → <strong>${ev.child_codon}</strong> (${ev.child_aa})`
                      : `<span style="color: ${PALETTE.neutral}; font-style: italic;">No discrete mutation reconstructed</span>`;
                    const typeHtml = ev
                      ? `<span class="ev-tag ${ev.is_synonymous ? 'ev-syn' : 'ev-nonsyn'}">${ev.is_synonymous ? 'Synonymous' : 'Non-synonymous'}</span>`
                      : '—';
                    return `
                      <tr>
                        <td class="font-mono font-bold">${bName}</td>
                        <td class="font-mono"><span class="ebf-branch-badge ebf-active">⚡ EBF ${formatEBF(ebf)}</span></td>
                        <td class="font-mono">${subStr}</td>
                        <td>${typeHtml}</td>
                      </tr>
                    `;
                  }).join('')}
                </tbody>
              </table>
            </div>
          </div>
        `;
}

export const MEME_DRAWER = Object.assign({}, FEL_DRAWER, {
  header: memeDrawerHeader,
  tabs: [siteOverviewTab(memeRatesGrid, memeBranchSection), SITE_TREE_TAB]
});

function memeDiagnostics(data) {
  return `
      <div class="diag-item">
        <span class="diag-lbl">Analysis Method</span>
        <span class="diag-val">MEME (Mixed Effects Model of Episodic Selection)</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">Hypothesis Test</span>
        <span class="diag-val font-mono">H₀: β⁺ ≤ α vs H₁: β⁺ > α</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">Codon Sites</span>
        <span class="diag-val font-mono">${data.numSites}</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">Taxa / Sequences</span>
        <span class="diag-val font-mono">${data.numTaxa}</span>
      </div>
    ` + mg94FitItems(data.raw || {});
}

export const MEME_CONFIG = Object.assign({}, FEL_CONFIG, {
  id: 'meme',
  kpis: (parent, hooks) => new KpiBar(parent, {
    cards: MEME_KPI_CARDS,
    blank: FEL_KPI_CARDS,
    controller: SITE_P_THRESHOLD,
    threshold: hooks.threshold,
    statistic: hooks.statistic,
    onThresholdChange: hooks.onThresholdChange,
    onStatisticChange: hooks.onStatisticChange
  }),
  table: (parent, hooks) => siteResultsTable(parent, hooks, MEME_TABLE),
  drawer: (container, hooks) => siteDrawer(container, hooks, MEME_DRAWER),
  diagnostics: memeDiagnostics
});
