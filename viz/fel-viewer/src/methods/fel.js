// FEL method config: site-level pervasive selection, alpha vs beta per codon
// The site components below also serve MEME; they read data.isMEME for its wording

import { FELData } from '../data/fel-data.js';
import { KpiBar } from '../components/kpi-bar.js';
import { EvidencePlot } from '../components/evidence-plot.js';
import { CandidatesTable } from '../components/candidates-table.js';
import { SiteDrawer } from '../components/site-drawer.js';
import { formatOmegaInf, formatThreshold } from '../core/format.js';

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

// Each factory gets the shell's hooks: { threshold, filterType, onThresholdChange(p),
// onSiteSelect(site, tab), onClose() }
export const FEL_CONFIG = {
  id: 'fel',
  Parser: FELData,
  defaultFilter: 'diversifying',
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
  table: (parent, hooks) => new CandidatesTable(parent, {
    threshold: hooks.threshold,
    filterType: hooks.filterType,
    onSiteSelect: (site, siteResult, action) => hooks.onSiteSelect(site, action)
  }),
  drawer: (container, hooks) => new SiteDrawer(container, {
    threshold: hooks.threshold,
    onClose: hooks.onClose
  }),
  diagnostics: felDiagnostics
};
