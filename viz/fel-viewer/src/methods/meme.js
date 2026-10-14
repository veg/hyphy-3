// MEME method config: episodic diversifying selection (beta+ vs alpha) per codon
// Shares FEL's site components and parser; only the diagnostics summary differs here

import { FEL_CONFIG, mg94FitItems } from './fel.js';

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
  diagnostics: memeDiagnostics
});
