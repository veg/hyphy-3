// Posterior Plot: the site track preset for posterior-probability site methods (FUBAR).
// Three modes share the call colouring at the live cutoff: P[β>α] on a fixed 0–1 axis
// with the cutoff rule at t; a signed mode, +P[β>α] or −P[α>β] for whichever side is
// larger, on a fixed −1..1 axis with rules at ±t (its ticks read as probabilities, the
// sign carried by the axis captions); and the posterior mean difference β − α.
// Rule labels sit at the plot's right end, over the dots.
// The method config supplies values and wording through a preset:
// { tone(r, track), pos(r), neg(r), diff(r), tooltip(r, track), ruleLabel(t),
//   labels: { pos, signed: [above, below], diff }, subtitles: { pos, signed, diff } }
// where each subtitle is (data) => text, so it can name the grid and the estimation method.

import { SiteTrackCanvas } from './site-track.js';
import { EvidenceBand } from './evidence-band.js';

// The signed mode's value: the larger posterior, negative for alpha > beta
export function signedPosterior(pos, neg) {
  return pos >= neg ? pos : -neg;
}

// A round tick step (1, 2 or 5 times a power of ten) giving about count intervals
export function niceStep(span, count = 5) {
  const raw = span / count;
  const power = 10 ** Math.floor(Math.log10(raw));
  return [1, 2, 5, 10].map(m => m * power).find(step => step >= raw);
}

// The β − α axis: the plotted values and zero, widened to whole tick steps
export function roundedDomain(values) {
  let lo = 0;
  let hi = 0;
  for (const v of values) {
    if (v < lo) lo = v;
    if (v > hi) hi = v;
  }
  const step = niceStep(hi - lo || 1);
  return [Math.floor(lo / step) * step, Math.ceil(hi / step) * step];
}

export function roundedTicks(min, max) {
  const step = niceStep(max - min);
  const ticks = [];
  for (let v = min; v <= max + step / 2; v += step) ticks.push(Number(v.toPrecision(12)));
  return ticks;
}

// Fresh bands per plot, since a band caches the points it last drew
export function posteriorPlotModes(preset) {
  const band = (options) => new EvidenceBand(Object.assign({
    tone: preset.tone,
    tooltip: preset.tooltip,
    shape: 'dot',
    rule: 'ge',
    ruleLabel: preset.ruleLabel,
    ruleLabelAt: 'end',
    tickFormat: (v) => v.toFixed(1)
  }, options));
  return [
    {
      id: 'pos',
      label: 'Posterior P[β>α]',
      subtitle: preset.subtitles.pos,
      // On the P[β>α] axis only the diversifying calls take a hue
      band: band({
        value: preset.pos,
        tone: (r, track) => (preset.tone(r, track) === 'div' ? 'div' : null),
        domain: [0, 1],
        ticks: [0, 0.5, 0.9, 1],
        label: preset.labels.pos
      })
    },
    {
      id: 'signed',
      label: 'Signed Posterior',
      subtitle: preset.subtitles.signed,
      band: band({
        value: (r) => signedPosterior(preset.pos(r), preset.neg(r)),
        domain: [-1, 1],
        ticks: [-0.9, -0.5, 0, 0.5, 0.9],
        tickFormat: (v) => Math.abs(v).toFixed(1),
        mirrorRule: true,
        label: preset.labels.signed
      })
    },
    {
      id: 'diff',
      label: 'Rate Diff (β − α)',
      subtitle: preset.subtitles.diff,
      band: band({
        value: preset.diff,
        domain: (rows) => roundedDomain(rows.map(preset.diff)),
        ticks: roundedTicks,
        tickFormat: (v) => String(v),
        rule: null,
        label: preset.labels.diff
      })
    }
  ];
}

export class PosteriorPlot extends SiteTrackCanvas {
  constructor(container, options = {}) {
    super(container, Object.assign({
      title: 'Posterior Evidence Landscape',
      exportName: 'hyphy_posterior_plot'
    }, options, { modes: posteriorPlotModes(options.preset) }));
  }

  // The subtitle names what the loaded result was run with
  setData(data) {
    this.subLabel.textContent = this.modes.find(m => m.id === this.mode).subtitle(data);
    super.setData(data);
  }
}
