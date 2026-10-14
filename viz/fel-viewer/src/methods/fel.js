// FEL method config: site-level pervasive selection, alpha vs beta per codon
// The site components below also serve MEME; the hero and drawer read data.isMEME for its wording

import { FELData } from '../data/fel-data.js';
import { KpiBar } from '../components/kpi-bar.js';
import { EvidencePlot } from '../components/evidence-plot.js';
import { ResultsTable, closestCalls } from '../components/results-table.js';
import { SiteDrawer } from '../components/site-drawer.js';
import { formatOmegaInf, formatThreshold, formatPValue } from '../core/format.js';
import { csvQuote } from '../core/export.js';

// Diagnostics rows for the global MG94 fit, from v3 model_fits or legacy fits
export function mg94FitItems(raw) {
  const codonFit = raw.model_fits?.codon_mg94 || raw.fits?.['Global MG94xREV'];
  if (!codonFit) return '';
  const logL = codonFit.log_likelihood ?? codonFit['Log Likelihood'];
  const aicc = codonFit.aicc ?? codonFit['AIC-c'];
  let html = '';
  if (logL != null) {
    html += `
          <div class="diag-item">
            <span class="diag-lbl">MG94 Log-Likelihood</span>
            <span class="diag-val font-mono">${Number(logL).toFixed(2)}</span>
          </div>
        `;
  }
  if (aicc != null) {
    html += `
          <div class="diag-item">
            <span class="diag-lbl">MG94 AICc</span>
            <span class="diag-val font-mono">${Number(aicc).toFixed(2)}</span>
          </div>
        `;
  }
  return html;
}

function felDiagnostics(data) {
  const m = data.analysisType ? data.analysisType.toUpperCase() : 'FEL';
  return `
      <div class="diag-item">
        <span class="diag-lbl">Analysis Method</span>
        <span class="diag-val">${m} (Fixed Effects Likelihood)</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">Codon Sites</span>
        <span class="diag-val font-mono">${data.numSites}</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">Taxa / Sequences</span>
        <span class="diag-val font-mono">${data.numTaxa}</span>
      </div>
      <div class="diag-item">
        <span class="diag-lbl">Global Mean dN/dS (ω)</span>
        <span class="diag-val font-mono">${formatOmegaInf(data.globalOmega, 3)}</span>
      </div>
    ` + mg94FitItems(data.raw || {});
}

// Site p-value cutoff shared by FEL and MEME
export const SITE_P_THRESHOLD = {
  title: 'Significance Cutoff',
  statistic: 'p',
  direction: 'le',
  min: 0.001,
  max: 0.2,
  step: 0.005,
  presets: [0.1, 0.05, 0.01],
  seed: 0.1,
  format: formatThreshold
};

// Share of all sites, as the KPI sub lines print it
export function sitePercent(count, stats) {
  return ((count / Math.max(1, stats.totalSites)) * 100).toFixed(1);
}

export const ALIGNMENT_KPI = {
  id: 'alignment',
  valueId: 'sites',
  subId: 'taxa',
  label: 'Alignment',
  placeholder: '— Taxa',
  value: (stats, { data }) => `${data.numSites.toLocaleString()} Codons`,
  sub: (stats, { data }) => `${data.numTaxa} Taxa`,
  tone: 'ink'
};

function omegaPressure(omega) {
  if (omega < 0.2) return 'Strong Purifying Constraint';
  if (omega < 0.7) return 'Pervasive Purifying Selection';
  if (omega <= 1.2) return 'Nearly Neutral Evolution';
  return 'Diversifying Selection';
}

export const FEL_KPI_CARDS = [
  ALIGNMENT_KPI,
  {
    id: 'omega',
    comment: 'Global Mean dN/dS Card',
    label: 'Global Mean dN/dS (ω)',
    placeholder: 'Gene-wide Pressure',
    value: (stats, { data }) => formatOmegaInf(data.globalOmega),
    sub: (stats, { data }) => omegaPressure(data.globalOmega),
    tone: 'ink'
  },
  {
    id: 'div',
    comment: 'Diversifying Sites (Adaptive) Card',
    cardClass: 'kpi-card-div',
    label: 'Diversifying (Positive)',
    placeholder: 'Sites (β > α)',
    value: (stats) => `▲ ${stats.diversifying}`,
    sub: (stats) => `${sitePercent(stats.diversifying, stats)}% of sites (β > α)`,
    tone: 'accent'
  },
  {
    id: 'pur',
    comment: 'Purifying Sites (Negative) Card',
    cardClass: 'kpi-card-pur',
    label: 'Purifying (Constraint)',
    placeholder: 'Sites (α > β)',
    value: (stats) => `▼ ${stats.purifying}`,
    sub: (stats) => `${sitePercent(stats.purifying, stats)}% of sites (α > β)`,
    tone: 'purifying'
  }
];

// Results table columns shared with MEME
export const CODON_COLUMN = {
  id: 'codon',
  header: 'Codon',
  width: 70,
  cell: (r) => `<td class="font-mono font-bold site-cell">
              <span class="site-number">#${r.site}</span>
            </td>`,
  csv: [{ header: 'Codon', value: (r) => r.site }],
  sort: (r) => r.site
};

export const LRT_COLUMN = {
  id: 'lrt',
  header: 'LRT',
  width: 80,
  cell: (r) => `<td class="font-mono lrt-cell">${r.lrt.toFixed(2)}</td>`,
  csv: [{ header: 'LRT', value: (r) => r.lrt }],
  sort: (r) => r.lrt
};

function siteHistory(r) {
  return r.history || (r.events?.length ? `${r.events.length} events` : '—');
}

export const SUBSTITUTIONS_COLUMN = {
  id: 'substitutions',
  header: 'Substitutions Reconstructed',
  cell: (r) => `<td class="sub-cell" title="${siteHistory(r)}">
              <div class="sub-text-wrapper">${siteHistory(r)}</div>
            </td>`,
  csv: [{ header: 'Substitutions', value: (r) => csvQuote(r.history) }]
};

export const SITE_ACTIONS_COLUMN = {
  id: 'actions',
  header: 'Inspect / Tree',
  width: 140,
  align: 'right',
  cell: (r) => `<td style="text-align: right; white-space: nowrap;">
              <button type="button" class="btn-inspect" data-action="summary" data-site="${r.site}">Inspect</button>
              <button type="button" class="btn-tree" data-action="tree" data-site="${r.site}" title="View Site Phylogeny Tree">🌳 Tree</button>
            </td>`
};

const felDiversifying = (r, t) => r.pValue <= t && r.beta > r.alpha;
const felPurifying = (r, t) => r.pValue <= t && r.alpha > r.beta;

function felRegime(r, t) {
  if (felDiversifying(r, t)) return 'Diversifying';
  return felPurifying(r, t) ? 'Purifying' : 'Neutral';
}

export const FEL_COLUMNS = [
  CODON_COLUMN,
  {
    id: 'rates',
    header: 'Rates (α : β)',
    width: 130,
    cell: (r) => `<td class="font-mono rates-cell">
              <span class="rate-alpha font-mono" title="Synonymous rate (α)">α=${r.alpha.toFixed(2)}</span>
              <span class="rate-sep">|</span>
              <span class="rate-beta font-mono" title="Non-synonymous rate (β)">β=${r.beta.toFixed(2)}</span>
            </td>`,
    csv: [
      { header: 'Alpha (dS)', value: (r) => r.alpha },
      { header: 'Beta (dN)', value: (r) => r.beta }
    ],
    sort: (r) => r.beta - r.alpha
  },
  {
    id: 'balance',
    header: 'Balance (dS vs dN)',
    width: 140,
    cell: (r) => {
      const maxRate = Math.max(0.1, r.alpha + r.beta);
      const alphaPct = Math.round((r.alpha / maxRate) * 100);
      const betaPct = Math.round((r.beta / maxRate) * 100);
      return `<td class="balance-cell">
              <div class="micro-balance-bar" title="dS: ${r.alpha.toFixed(2)} vs dN: ${r.beta.toFixed(2)}">
                <div class="bar-alpha" style="width: ${alphaPct}%;"></div>
                <div class="bar-beta" style="width: ${betaPct}%;"></div>
              </div>
            </td>`;
    }
  },
  {
    id: 'omega',
    header: 'ω (dN/dS)',
    width: 90,
    cell: (r, { threshold }) => {
      const badge = { Diversifying: 'omega-div', Purifying: 'omega-pur', Neutral: '' }[felRegime(r, threshold)];
      return `<td class="omega-cell">
              <span class="omega-badge ${badge}">${formatOmegaInf(r.omega)}</span>
            </td>`;
    },
    csv: [{ header: 'dN/dS (Omega)', value: (r) => isFinite(r.omega) ? r.omega : 'Infinity' }],
    sort: (r) => r.omega
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
    id: 'regime',
    csv: [{ header: 'Regime', value: (r, { threshold }) => felRegime(r, threshold) }]
  },
  SUBSTITUTIONS_COLUMN,
  SITE_ACTIONS_COLUMN
];

// FEL filters also carry the empty state's wording: the unit that missed the
// cutoff, the evidence it lacked and the opposite direction to jump to
export const FEL_FILTERS = [
  {
    id: 'diversifying',
    key: 'div',
    glyph: '▲',
    label: 'Diversifying',
    unit: 'diversifying site',
    evidence: 'positive selection (β > α)',
    opposite: 'purifying',
    leans: (r) => r.beta > r.alpha,
    count: (stats) => stats.diversifying,
    test: felDiversifying,
    subtitle: (n, t) => `Showing ${n} diversifying sites (β > α, p ≤ ${formatThreshold(t)})`
  },
  {
    id: 'purifying',
    key: 'pur',
    glyph: '▼',
    label: 'Purifying',
    unit: 'purifying site',
    evidence: 'purifying selection (α > β)',
    opposite: 'diversifying',
    leans: (r) => r.alpha > r.beta,
    count: (stats) => stats.purifying,
    test: felPurifying,
    subtitle: (n, t) => `Showing ${n} purifying sites (α > β, p ≤ ${formatThreshold(t)})`
  },
  {
    id: 'significant',
    key: 'sig',
    label: 'All Significant',
    unit: 'site',
    count: (stats) => stats.significant,
    test: (r, t) => r.pValue <= t,
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

// Jump to the opposite direction when it has sites, else show this direction's
// closest calls by p
function felEmptyState(t, { filter, stats, rows }) {
  const title = `No ${filter.unit} reaches p ≤ ${formatThreshold(t)}`;
  const desc = filter.evidence && `No codon site showed evidence of ${filter.evidence} at the current threshold.`;
  const other = FEL_FILTERS.find(f => f.id === filter.opposite);
  const otherCount = other ? other.count(stats) : 0;
  if (otherCount > 0) {
    return { icon: filter.glyph, title, desc, jump: { filter: other.id, label: `View ${other.label} Sites (${other.glyph} ${otherCount})` } };
  }
  const leaning = filter.leans ? rows.filter(filter.leans) : rows;
  return { icon: filter.glyph, title, desc, fallback: closestCalls(leaning, r => r.pValue) };
}

export const SITE_SEARCH = {
  placeholder: 'Search codon or AA...',
  match: (r, q) => String(r.site).includes(q)
    || (r.history || '').toLowerCase().includes(q)
    || (r.rootAA || '').toLowerCase().includes(q),
  noMatch: { title: 'No Matching Sites', desc: 'No codon sites match your current filter and search query.' }
};

export const FEL_TABLE = {
  title: 'Candidate Sites',
  placeholder: 'Filtered by selective regime',
  rowKey: 'site',
  rows: (data) => data.siteResults,
  columns: FEL_COLUMNS,
  filters: FEL_FILTERS,
  defaultFilter: 'diversifying',
  defaultSort: { column: 'codon', dir: 'asc' },
  search: SITE_SEARCH,
  emptyState: felEmptyState,
  exportName: 'hyphy_fel_candidates'
};

export function siteResultsTable(parent, hooks, table) {
  return new ResultsTable(parent, Object.assign({}, table, {
    threshold: hooks.threshold,
    filterType: hooks.filterType,
    onSelect: (site, row, action) => hooks.onSiteSelect(site, action)
  }));
}

// Each factory gets the shell's hooks: { threshold, filterType, onThresholdChange(p),
// onSiteSelect(site, tab), onClose() }
export const FEL_CONFIG = {
  id: 'fel',
  Parser: FELData,
  defaultFilter: FEL_TABLE.defaultFilter,
  kpis: (parent, hooks) => new KpiBar(parent, {
    cards: FEL_KPI_CARDS,
    controller: SITE_P_THRESHOLD,
    threshold: hooks.threshold,
    onThresholdChange: hooks.onThresholdChange
  }),
  hero: (parent, hooks) => new EvidencePlot(parent, {
    threshold: hooks.threshold,
    onSiteSelect: (site) => hooks.onSiteSelect(site)
  }),
  table: (parent, hooks) => siteResultsTable(parent, hooks, FEL_TABLE),
  drawer: (container, hooks) => new SiteDrawer(container, {
    threshold: hooks.threshold,
    onClose: hooks.onClose
  }),
  diagnostics: felDiagnostics
};
