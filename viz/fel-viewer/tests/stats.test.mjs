import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import {
  erfc, logGamma, gammaQ, chi2Survival, benjaminiHochberg, holm, enforceStepDown, poissonBinomialCI
} from '../src/core/stats.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';
const fixture = name => JSON.parse(readFileSync(new URL(name, FIXTURES), 'utf8'));

const close = (actual, expected, rel = 1e-6) =>
  assert.ok(Math.abs(actual - expected) <= rel * Math.abs(expected), `${actual} vs ${expected}`);
const allClose = (actual, expected) => {
  assert.equal(actual.length, expected.length);
  actual.forEach((a, i) => close(a, expected[i], 1e-12));
};

test('erfc matches reference values across both branches', () => {
  close(erfc(0), 1);
  close(erfc(1), 0.157299207050285);
  close(erfc(-1), 1.842700792949715);
  close(erfc(3), 2.209049699858544e-5);
  close(erfc(10), 2.088487583762545e-45, 1e-9);
  close(erfc(15), 7.212994172451206e-100, 1e-9);
  assert.equal(erfc(40), 0);
});

test('erfc recomputes the underflowed RELAX p-values', () => {
  // HyPhy 2.5 reports p = 0 at LR = 113.95; hyphy3 reports 3.96e-33 at LR = 143.78
  close(erfc(Math.sqrt(113.95 / 2)), 1.3364e-26, 1e-4);
  close(erfc(Math.sqrt(143.78 / 2)), 3.969e-33, 1e-3);
});

test('logGamma and gammaQ', () => {
  close(logGamma(5), Math.log(24), 1e-12);
  close(logGamma(0.5), Math.log(Math.sqrt(Math.PI)), 1e-12);
  close(gammaQ(1, 2), Math.exp(-2), 1e-12);
  assert.equal(gammaQ(3, 0), 1);
});

test('chi2Survival at the 5% critical values', () => {
  close(chi2Survival(3.841458820694124, 1), 0.05);
  close(chi2Survival(5.991464547107979, 2), 0.05, 1e-12);
  close(chi2Survival(7.814727903251178, 3), 0.05, 1e-12);
  close(chi2Survival(11.070497693516351, 5), 0.05, 1e-12);
  close(chi2Survival(100, 5), 5.285148360943e-20, 1e-9);
  assert.equal(chi2Survival(0, 3), 1);
  assert.equal(chi2Survival(-1), 1);
});

test('benjaminiHochberg is a monotone step-up', () => {
  // Without the running minimum, rank 2 would get 0.06 > rank 3's 0.04 * 4 / 3
  allClose(benjaminiHochberg([0.01, 0.04, 0.03, 0.2]), [0.04, 4 * 0.04 / 3, 4 * 0.04 / 3, 0.2]);
  allClose(benjaminiHochberg([0.01, 0.01]), [0.01, 0.01]);
  allClose(benjaminiHochberg([0.01, 0.5], 10), [0.1, 1]);
});

test('holm is a monotone step-down', () => {
  allClose(holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06]);
  allClose(holm([0.4, 0.3]), [0.6, 0.6]);
  allClose(holm([0.6, 0.01]), [0.6, 0.02]);
  allClose(holm([0.6, 0.01], 5), [1, 0.05]);
});

test('enforceStepDown takes a running maximum in raw-p order', () => {
  assert.deepEqual(enforceStepDown([0.56, 0.54, 0.9], [0.023, 0.024, 0.1]), [0.56, 0.56, 0.9]);
  assert.deepEqual(enforceStepDown([0.2, 0.1, 1.5], [0.2, 0.1, 0.3]), [0.2, 0.1, 1]);
});

test('PRIME v3 q-values and within-site Holm match the recomputation', { skip: noFixtures }, () => {
  const data = fixture('prime_cd2.json').site_results.data;
  assert.deepEqual(benjaminiHochberg(data.omnibus_p_value), data.omnibus_q_value);
  const properties = Object.keys(data).filter(k => k.startsWith('adj_p_')).map(k => k.slice(6));
  data.site.forEach((_, s) => {
    const adjusted = holm(properties.map(name => data[`p_${name}`][s]));
    properties.forEach((name, k) => close(adjusted[k], data[`adj_p_${name}`][s], 1e-12));
  });
});

test('HyPhy 2.5 aBSREL corrected p is made monotone', { skip: noFixtures }, () => {
  const branches = fixture('absrel_bglobin.hyphy25.json')['branch attributes']['0'];
  const tested = Object.keys(branches).filter(name => 'Corrected P-value' in branches[name]);
  const pHolm = enforceStepDown(
    tested.map(name => branches[name]['Corrected P-value']),
    tested.map(name => branches[name]['Uncorrected P-value'])
  );
  const changed = tested.filter((name, i) => pHolm[i] !== branches[name]['Corrected P-value']);
  assert.deepEqual(changed, ['HAMSTER']);
  assert.equal(pHolm[tested.indexOf('HAMSTER')], branches.Node26['Corrected P-value']);
});

test('poissonBinomialCI bounds the false-call count by the CDF crossings', () => {
  assert.deepEqual(poissonBinomialCI([]), [0, 0]);
  assert.deepEqual(poissonBinomialCI([0, 0, 0]), [0, 0]);
  assert.deepEqual(poissonBinomialCI([1, 1]), [2, 2]);
  // Binomial(10, 0.5): CDF(1) = 0.0107, CDF(2) = 0.0547, CDF(7) = 0.9453, CDF(8) = 0.9893
  assert.deepEqual(poissonBinomialCI(new Array(10).fill(0.5)), [2, 8]);
  // A wider level narrows the interval
  assert.deepEqual(poissonBinomialCI(new Array(10).fill(0.5), 0.5), [4, 6]);
});

test('poissonBinomialCI reproduces the FUBAR run summary', { skip: noFixtures }, () => {
  const json = fixture('fubar_cd2.json');
  const { prob_pos } = json.site_results.data;
  const errors = prob_pos.filter(p => p >= 0.9).map(p => 1 - p);
  assert.deepEqual(poissonBinomialCI(errors), json.summary.false_positives_ci_95);
});
