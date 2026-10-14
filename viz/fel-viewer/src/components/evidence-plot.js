// Evidence Plot: the FEL/MEME preset of the site track (Manhattan -log10 p, signed
// evidence and rate difference modes, ▲/▼ glyphs by direction at the live threshold)

import { formatPValue, formatOmegaInf } from '../core/format.js';
import { SiteTrackCanvas } from './site-track.js';
import { EvidenceBand } from './evidence-band.js';

// FEL calls a significant site by direction; MEME's test is one-sided, so every
// significant site is episodic diversifying
export function evidenceTone(r, track) {
  if (r.pValue > track.threshold) return null;
  if (track.data.isMEME || r.beta > r.alpha) return 'div';
  return r.alpha > r.beta ? 'pur' : null;
}

function evidenceTooltip(r, track) {
  const isMEME = track.data.isMEME;
  const tone = evidenceTone(r, track);

  let selBadge = '<span class="tip-badge tip-neutral">Neutral</span>';
  if (tone === 'div') {
    selBadge = isMEME
      ? '<span class="tip-badge tip-div">▲ Episodic Diversifying</span>'
      : '<span class="tip-badge tip-div">▲ Diversifying Selection</span>';
  } else if (tone === 'pur') {
    selBadge = '<span class="tip-badge tip-pur">▼ Purifying Selection</span>';
  }

  const pStr = formatPValue(r.pValue);
  const subStr = r.history ? `<div class="tip-row tip-sub"><span class="tip-lbl">Substitutions:</span> <span class="tip-val-sub">${r.history}</span></div>` : '';

  let statsRows = '';
  if (isMEME) {
    const betaPlus = r.beta_plus ?? r.beta;
    const pPlusPct = Math.round((r.p_plus ?? 0) * 100);
    const branches = r.branchesUnderSelection || (r.branchEBF ? Object.values(r.branchEBF).filter(v => v >= 100).length : 0);
    statsRows = `
        <div class="tip-row"><span class="tip-lbl">p-value:</span> <span class="tip-val font-mono">${pStr}</span></div>
        <div class="tip-row"><span class="tip-lbl">LRT:</span> <span class="tip-val font-mono">${r.lrt.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">α (dS):</span> <span class="tip-val font-mono">${r.alpha.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">β⁺ (dN):</span> <span class="tip-val font-mono">${betaPlus.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">p⁺ (% branches):</span> <span class="tip-val font-mono">${pPlusPct}%</span></div>
        <div class="tip-row"><span class="tip-lbl">Branches:</span> <span class="tip-val font-mono">${branches} (EBF≥100)</span></div>
      `;
  } else {
    const omegaStr = formatOmegaInf(r.omega);
    statsRows = `
        <div class="tip-row"><span class="tip-lbl">p-value:</span> <span class="tip-val font-mono">${pStr}</span></div>
        <div class="tip-row"><span class="tip-lbl">LRT:</span> <span class="tip-val font-mono">${r.lrt.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">α (dS):</span> <span class="tip-val font-mono">${r.alpha.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">β (dN):</span> <span class="tip-val font-mono">${r.beta.toFixed(2)}</span></div>
        <div class="tip-row"><span class="tip-lbl">dN/dS (ω):</span> <span class="tip-val font-mono">${omegaStr}</span></div>
      `;
  }

  return `
      <div class="tip-header">
        <span class="tip-title">Codon #${r.site}</span>
        ${selBadge}
      </div>
      <div class="tip-grid">
        ${statsRows}
      </div>
      ${subStr}
      <div class="tip-hint">Click point to inspect site tree &amp; branches</div>
    `;
}

// Manhattan: 0 up to 15% above the strongest site, never below 3
function manhattanDomain(rows) {
  let yMax = 3.0;
  for (const r of rows) {
    if (r.negLogP > yMax) yMax = r.negLogP;
  }
  return [0, Math.ceil(yMax * 1.15)];
}

function signedDomain(rows) {
  let maxAbs = 2.0;
  for (const r of rows) {
    const v = Math.abs(r.negLogP);
    if (v > maxAbs) maxAbs = v;
  }
  maxAbs = Math.ceil(maxAbs * 1.15);
  return [-maxAbs, maxAbs];
}

function rateDiffDomain(rows) {
  let minD = 0;
  let maxD = 1.0;
  for (const r of rows) {
    const diff = r.beta - r.alpha;
    if (diff < minD) minD = diff;
    if (diff > maxD) maxD = diff;
  }
  const margin = Math.max(1.0, (maxD - minD) * 0.1);
  return [Math.floor(minD - margin), Math.ceil(maxD + margin)];
}

// MEME plots its episodic rate ratio in signed mode
function signedEvidence(r, track) {
  if (track.data.isMEME) {
    const bp = r.beta_plus ?? r.beta;
    return r.alpha > 1e-4 ? (bp / r.alpha) : (bp > 0 ? 10 : 1);
  }
  const sign = r.beta >= r.alpha ? 1 : -1;
  return sign * r.negLogP;
}

function rateDiff(r, track) {
  const b = track.data.isMEME ? (r.beta_plus ?? r.beta) : r.beta;
  return b - r.alpha;
}

// Fresh bands per plot, since a band caches the points it last drew
export function evidencePlotModes() {
  const band = (options) => new EvidenceBand(Object.assign({ tone: evidenceTone, tooltip: evidenceTooltip }, options));
  return [
    {
      id: 'manhattan',
      label: 'Evidence (−log₁₀ p)',
      subtitle: (data) => data?.isMEME
        ? '−log₁₀(p) Evidence for Episodic Diversifying Selection'
        : '−log₁₀(p) vs Codon Position',
      band: band({ value: (r) => r.negLogP, domain: manhattanDomain })
    },
    {
      id: 'signed',
      label: 'Signed Evidence',
      subtitle: (data) => data?.isMEME
        ? 'Episodic Rate Ratio: β⁺ / α'
        : 'Signed Evidence: +Diversifying / −Purifying',
      band: band({ value: signedEvidence, domain: signedDomain, mirrorRule: true })
    },
    {
      id: 'diff',
      label: 'Rate Diff (β − α)',
      subtitle: (data) => data?.isMEME
        ? 'Episodic Rate Difference: β⁺ − α'
        : 'Absolute Rate Difference: β (dN) − α (dS)',
      band: band({ value: rateDiff, domain: rateDiffDomain, rule: null })
    }
  ];
}

export class EvidencePlot extends SiteTrackCanvas {
  constructor(container, options = {}) {
    super(container, Object.assign({
      title: 'Selection Evidence Landscape',
      subtitle: '−log₁₀(p) vs Codon Position',
      modes: evidencePlotModes(),
      exportName: 'hyphy_fel_evidence_plot'
    }, options));
  }
}
