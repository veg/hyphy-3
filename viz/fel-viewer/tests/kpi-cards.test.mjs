import { test } from 'node:test';
import assert from 'node:assert/strict';
import { thresholdLabel } from '../src/components/kpi-bar.js';
import { FEL_KPI_CARDS, SITE_P_THRESHOLD } from '../src/methods/fel.js';
import { MEME_KPI_CARDS } from '../src/methods/meme.js';

const stats = { diversifying: 3, purifying: 12, neutral: 185, totalSites: 200 };
const view = { data: { numSites: 1200, numTaxa: 8, globalOmega: 0.45 }, threshold: 0.05 };
const render = (cards) => cards.map(c => [c.label, c.value(stats, view), c.sub(stats, view)]);

test('the threshold controller labels its cutoff in either direction', () => {
  assert.equal(thresholdLabel(SITE_P_THRESHOLD, 0.1), 'p ≤ 0.10');
  assert.equal(thresholdLabel(SITE_P_THRESHOLD, 0.005), 'p ≤ 0.005');
  assert.equal(thresholdLabel({ ...SITE_P_THRESHOLD, statistic: 'PP', direction: 'ge' }, 0.9), 'PP ≥ 0.90');
});

test('FEL cards report the alignment, omega and both selection directions', () => {
  assert.deepEqual(render(FEL_KPI_CARDS), [
    ['Alignment', `${(1200).toLocaleString()} Codons`, '8 Taxa'],
    ['Global Mean dN/dS (ω)', '0.45', 'Pervasive Purifying Selection'],
    ['Diversifying (Positive)', '▲ 3', '1.5% of sites (β > α)'],
    ['Purifying (Constraint)', '▼ 12', '6.0% of sites (α > β)']
  ]);
});

test('MEME cards stand in for FEL\'s blank cards one for one', () => {
  assert.deepEqual(MEME_KPI_CARDS.map(c => c.id), FEL_KPI_CARDS.map(c => c.id));
  assert.deepEqual(render(MEME_KPI_CARDS), [
    ['Alignment', `${(1200).toLocaleString()} Codons`, '8 Taxa'],
    ['Analysis Model', 'MEME', 'Episodic Diversifying'],
    ['Episodically Diversifying', '▲ 3', '1.5% of sites (p ≤ 0.05)'],
    ['Neutral / Pervasive', '185', 'Conserved or Neutral sites']
  ]);
});
