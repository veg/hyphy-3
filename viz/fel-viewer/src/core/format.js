/**
 * House number formatters shared by every method viewer.
 */

import { erfc } from './stats.js';

// Threshold labels: three decimals below 0.01 so 0.005 does not read as 0.01
export function formatThreshold(t) {
  return t.toFixed(t < 0.01 ? 3 : 2);
}

// FEL/MEME site p-values: scientific notation below 1e-4
export function formatPValue(p, digits = 2) {
  return p < 0.0001 ? p.toExponential(digits) : p.toFixed(4);
}

// FEL site omega, where alpha = 0 makes beta/alpha unbounded
export function formatOmegaInf(omega, digits = 2) {
  return omega >= 100 ? '∞' : omega.toFixed(digits);
}

// Empirical Bayes factors; huge values are capped rather than spelled out
export function formatEBF(val, cap = '≥10,000') {
  if (val == null) return '—';
  if (val >= 10000) return cap;
  if (val >= 1000) return `${(val / 1000).toFixed(1)}k`;
  if (val >= 100) return val.toFixed(0);
  return val.toFixed(1);
}

// HyPhy 2.5 underflows tiny p-values to exactly 0; say so instead of printing 0
export function formatP(p, digits = 2) {
  return p === 0 ? 'p < 1e-16' : formatPValue(p, digits);
}

// The '(computed)' sub-value for an underflowed p: the 1-df chi-square tail
// erfc(sqrt(LR/2)), or null when there is no positive LR to recompute from
export function formatComputedP(lrt) {
  if (!(lrt > 0)) return null;
  return `${formatPValue(erfc(Math.sqrt(lrt / 2)))} (computed)`;
}

// Omega from rate-class fits: legacy sentinels (1e10) and large v3 values
// read as '>1000', and a class pinned at the floor reads as a plain '0'
export function formatOmega(omega, digits = 2) {
  if (omega >= 1000) return '>1000';
  if (omega <= 0) return '0';
  return omega.toFixed(digits);
}

export function formatPosterior(p) {
  return p.toFixed(3);
}

// Substitution counts are fractional when ties are split between paths
export function formatCount(n) {
  return Number.isInteger(n) ? String(n) : n.toFixed(1);
}

// Posterior cutoffs: two decimals unless the cutoff needs three (0.995, 0.999), so a
// cutoff near 1 never reads as 1.00
export function formatPosteriorCutoff(t) {
  const two = t.toFixed(2);
  return Math.abs(Number(two) - t) < 1e-9 ? two : t.toFixed(3);
}
