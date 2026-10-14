// Evidence Plot: the site track preset for alpha/beta site methods (Manhattan -log10 p,
// signed evidence and rate difference modes, ▲/▼ glyphs by direction at the live threshold).
// The method config supplies the wording and the per-site values through a preset:
// { tone(r, track), badges: { div, pur }, rows(r), signed(r), diff(r), subtitles: { manhattan, signed, diff } }

import { formatPValue } from '../core/format.js';
import { SiteTrackCanvas } from './site-track.js';
import { EvidenceBand } from './evidence-band.js';

// Tooltip row for a label and its monospaced value
export function evidenceTipRow(label, value) {
  return `<div class="tip-row"><span class="tip-lbl">${label}:</span> <span class="tip-val font-mono">${value}</span></div>`;
}

// p-value, LRT and α lead every tooltip; the preset adds its own rates after them
function evidenceTooltip(preset, r, track) {
  const tone = preset.tone(r, track);
  const badge = tone ? preset.badges[tone] : null;
  const selBadge = badge
    ? `<span class="tip-badge tip-${tone}">${badge}</span>`
    : '<span class="tip-badge tip-neutral">Neutral</span>';

  const subStr = r.history ? `<div class="tip-row tip-sub"><span class="tip-lbl">Substitutions:</span> <span class="tip-val-sub">${r.history}</span></div>` : '';
  const rows = [
    evidenceTipRow('p-value', formatPValue(r.pValue)),
    evidenceTipRow('LRT', r.lrt.toFixed(2)),
    evidenceTipRow('α (dS)', r.alpha.toFixed(2)),
    ...preset.rows(r)
  ];
  const statsRows = rows.map(row => `
        ${row}`).join('') + `
      `;

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

// Fresh bands per plot, since a band caches the points it last drew
export function evidencePlotModes(preset) {
  const band = (options) => new EvidenceBand(Object.assign({
    tone: preset.tone,
    tooltip: (r, track) => evidenceTooltip(preset, r, track)
  }, options));
  return [
    {
      id: 'manhattan',
      label: 'Evidence (−log₁₀ p)',
      subtitle: () => preset.subtitles.manhattan,
      band: band({ value: (r) => r.negLogP, domain: manhattanDomain })
    },
    {
      id: 'signed',
      label: 'Signed Evidence',
      subtitle: () => preset.subtitles.signed,
      band: band({ value: (r) => preset.signed(r), domain: signedDomain, mirrorRule: true })
    },
    {
      id: 'diff',
      label: 'Rate Diff (β − α)',
      subtitle: () => preset.subtitles.diff,
      band: band({ value: (r) => preset.diff(r), domain: rateDiffDomain, rule: null })
    }
  ];
}

export class EvidencePlot extends SiteTrackCanvas {
  constructor(container, options = {}) {
    super(container, Object.assign({
      title: 'Selection Evidence Landscape',
      subtitle: '−log₁₀(p) vs Codon Position',
      exportName: 'hyphy_fel_evidence_plot'
    }, options, { modes: evidencePlotModes(options.preset) }));
  }
}
