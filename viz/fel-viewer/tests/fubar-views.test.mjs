import { test } from 'node:test';
import assert from 'node:assert/strict';
import { rateGridAxes, ordinalIndex, rateGridShade, rateGridTicks } from '../src/components/rate-grid-heatmap.js';
import { signedPosterior, posteriorPlotModes, niceStep, roundedDomain, roundedTicks } from '../src/components/posterior-plot.js';
import { EvidenceBand } from '../src/components/evidence-band.js';
import { FUBAR_EVIDENCE } from '../src/methods/fubar.js';

const FUBAR_AXIS = [0, 1 / 14, 2 / 14, 3 / 14, 4 / 14, 5 / 14, 6 / 14, 7 / 14, 8 / 14, 9 / 14, 10 / 14, 11 / 14, 12 / 14,
  1, 1.227, 2.815, 7.125, 15.52, 29.36, 50];

test('the grid axes are read alpha-major: alpha steps by row of n, beta within it', () => {
  const alpha = [];
  const beta = [];
  for (const a of [0, 1, 5]) for (const b of [0, 2, 9]) { alpha.push(a); beta.push(b); }
  assert.deepEqual(rateGridAxes({ alpha, beta, dims: 3 }), { alpha: [0, 1, 5], beta: [0, 2, 9] });
});

test('a posterior mean sits between its grid neighbours on the ordinal axis', () => {
  assert.equal(ordinalIndex(FUBAR_AXIS, 0), 0);
  assert.equal(ordinalIndex(FUBAR_AXIS, -1), 0);
  assert.equal(ordinalIndex(FUBAR_AXIS, 1), 13);
  assert.ok(Math.abs(ordinalIndex(FUBAR_AXIS, (1.227 + 2.815) / 2) - 14.5) < 1e-12);
  assert.equal(ordinalIndex(FUBAR_AXIS, 80), 19);
});

test('the fill is a sqrt greyscale from the panel to near-ink', () => {
  assert.equal(rateGridShade(0, 0.01), 'rgb(248, 248, 248)');
  assert.equal(rateGridShade(0.01, 0.01), 'rgb(34, 34, 34)');
  assert.equal(rateGridShade(0.0025, 0.01), 'rgb(141, 141, 141)');
});

test('grid tick labels thin out to fit and always include the last value', () => {
  assert.deepEqual(rateGridTicks(20, 16), [0, 2, 4, 6, 8, 10, 12, 14, 16, 19]);
  assert.deepEqual(rateGridTicks(20, 8), [0, 4, 8, 12, 19]);
  assert.deepEqual(rateGridTicks(5, 40), [0, 1, 2, 3, 4]);
});

test('the signed mode shows whichever posterior is larger, negative for α > β', () => {
  assert.equal(signedPosterior(0.91, 0.05), 0.91);
  assert.equal(signedPosterior(0.02, 0.95), -0.95);
  assert.equal(signedPosterior(0.4, 0.4), 0.4);
});

test('the β − α axis is widened to round steps', () => {
  assert.equal(niceStep(38.9), 10);
  assert.equal(niceStep(1), 0.2);
  assert.deepEqual(roundedDomain([-28.8, 10.1, 3]), [-30, 20]);
  assert.deepEqual(roundedTicks(-30, 20), [-30, -20, -10, 0, 10, 20]);
  assert.deepEqual(roundedDomain([]), [0, 0]);
});

test('posterior modes: fixed probability axes with the cutoff rule, only diversifying in hue on P[β>α]', () => {
  const [pos, signed, diff] = posteriorPlotModes(FUBAR_EVIDENCE);
  assert.deepEqual([pos.id, signed.id, diff.id], ['pos', 'signed', 'diff']);
  assert.deepEqual(pos.band.options.domain, [0, 1]);
  assert.deepEqual(signed.band.options.domain, [-1, 1]);
  assert.equal(pos.band.options.rule, 'ge');
  assert.equal(signed.band.options.mirrorRule, true);
  assert.equal(diff.band.options.rule, null);
  assert.equal(pos.band.ruleLabel(0.999), 'P ≥ 0.999');

  const track = { threshold: 0.9 };
  const purifying = { probPos: 0.02, probNeg: 0.95, diff: -3 };
  assert.equal(pos.band.options.tone(purifying, track), null);
  assert.equal(signed.band.options.tone(purifying, track), 'pur');
  assert.equal(diff.band.options.tone(purifying, track), 'pur');
  assert.equal(signed.band.options.value(purifying), -0.95);

  const data = { settings: { gridPoints: 20, method: 'Variational-Bayes' } };
  assert.equal(pos.subtitle(data), 'Posterior probability of β > α per codon (FUBAR, 20×20 grid, Variational-Bayes)');
});

test('an evidence band can derive its ticks from the domain', () => {
  const band = new EvidenceBand({ ticks: (min, max) => [min, (min + max) / 2, max] });
  band.yMin = -4;
  band.yMax = 2;
  assert.deepEqual(band.tickValues(), [-4, -1, 2]);
});
