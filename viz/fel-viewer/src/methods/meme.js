// MEME method config: episodic diversifying selection (beta+ vs alpha) per codon
// Shares FEL's site components and parser; the KPI cards and diagnostics differ here

import { FEL_CONFIG, mg94FitItems, SITE_P_THRESHOLD, FEL_KPI_CARDS, ALIGNMENT_KPI, sitePercent } from './fel.js';
import { KpiBar } from '../components/kpi-bar.js';
import { formatThreshold } from '../core/format.js';

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
    onThresholdChange: hooks.onThresholdChange
  }),
  diagnostics: memeDiagnostics
});
