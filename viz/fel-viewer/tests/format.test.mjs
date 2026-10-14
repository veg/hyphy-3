import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  formatThreshold, formatPValue, formatOmegaInf, formatEBF,
  formatP, formatComputedP, formatOmega, formatPosterior, formatCount, formatPosteriorCutoff
} from '../src/core/format.js';

test('FEL formatters keep their existing output', () => {
  assert.equal(formatThreshold(0.1), '0.10');
  assert.equal(formatThreshold(0.005), '0.005');
  assert.equal(formatPValue(0.0123), '0.0123');
  assert.equal(formatPValue(3.96e-33), '3.96e-33');
  assert.equal(formatPValue(3.96e-33, 3), '3.960e-33');
  assert.equal(formatOmegaInf(0.4567), '0.46');
  assert.equal(formatOmegaInf(0.4567, 3), '0.457');
  assert.equal(formatOmegaInf(Infinity), '∞');
});

test('formatEBF caps and abbreviates', () => {
  assert.equal(formatEBF(null), '—');
  assert.equal(formatEBF(12.34), '12.3');
  assert.equal(formatEBF(123.4), '123');
  assert.equal(formatEBF(1234), '1.2k');
  assert.equal(formatEBF(20000), '≥10,000');
  assert.equal(formatEBF(20000, '≥10k'), '≥10k');
});

test('formatP flags an underflowed zero', () => {
  assert.equal(formatP(0), 'p < 1e-16');
  assert.equal(formatP(0.05), '0.0500');
});

test('formatComputedP recomputes the chi-square tail from LR', () => {
  assert.equal(formatComputedP(0), null);
  assert.equal(formatComputedP(undefined), null);
  assert.equal(formatComputedP(113.95), '1.34e-26 (computed)');
});

test('formatOmega clamps sentinels and the floor', () => {
  assert.equal(formatOmega(1e10), '>1000');
  assert.equal(formatOmega(1000), '>1000');
  assert.equal(formatOmega(0), '0');
  assert.equal(formatOmega(2.345), '2.35');
});

test('formatPosterior and formatCount', () => {
  assert.equal(formatPosterior(0.95123), '0.951');
  assert.equal(formatCount(341), '341');
  assert.equal(formatCount(193.75), '193.8');
});

test('posterior cutoffs keep the third decimal only when they need it', () => {
  assert.equal(formatPosteriorCutoff(0.9), '0.90');
  assert.equal(formatPosteriorCutoff(0.95), '0.95');
  assert.equal(formatPosteriorCutoff(0.5), '0.50');
  assert.equal(formatPosteriorCutoff(0.995), '0.995');
  assert.equal(formatPosteriorCutoff(0.999), '0.999');
  assert.equal(formatPosteriorCutoff(0.9 + 0.001 * 3), '0.903');
});
