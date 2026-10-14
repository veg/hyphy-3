// FEL method config: site-level pervasive selection, alpha vs beta per codon
// The site components below also serve MEME; they read data.isMEME for its wording

import { FELData } from '../data/fel-data.js';
import { KpiBar } from '../components/kpi-bar.js';
import { EvidencePlot } from '../components/evidence-plot.js';
import { CandidatesTable } from '../components/candidates-table.js';
import { SiteDrawer } from '../components/site-drawer.js';
import { formatOmegaInf } from '../core/format.js';

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

// Each factory gets the shell's hooks: { threshold, filterType, onThresholdChange(p),
// onSiteSelect(site, tab), onClose() }
export const FEL_CONFIG = {
  id: 'fel',
  Parser: FELData,
  defaultFilter: 'diversifying',
  kpis: (parent, hooks) => new KpiBar(parent, {
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
