import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import { FUBARData, fubarSiteClass, fubarFalsePositives, fubarPriorNonPositive } from '../src/data/fubar-data.js';
import { configForResult } from '../src/methods/registry.js';
import { FUBAR_CONFIG, FUBAR_KPI_CARDS, FUBAR_FILTERS, FUBAR_POSTERIOR_CUTOFF, fubarEmptyState, fubarDiagnostics } from '../src/methods/fubar.js';
import { thresholdLabel, sliderAttributes, thresholdPasses } from '../src/components/kpi-bar.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';
const fixture = name => JSON.parse(readFileSync(new URL(name, FIXTURES), 'utf8'));

const close = (actual, expected, tol = 1e-9) =>
  assert.ok(Math.abs(actual - expected) <= tol, `${actual} vs ${expected}`);

// hyphy3 --legacy-json output, shaped as fubar.hpp to_legacy_json writes it: key-sorted
// object headers, one partition in site order, no tree and no branch attributes
function hyphy3Legacy(v3) {
  const d = v3.site_results.data;
  return {
    analysis: { info: 'Perform a Fast Unbiased AppRoximate Bayesian (FUBAR) analysis of a coding sequence alignment.', version: '2.2', reference: 'FUBAR (2013)' },
    input: { file: 'cd2.fna', sequences: 10, sites: d.site.length, 'partition count': 1 },
    fits: v3.fits,
    settings: v3.settings,
    grid: v3.grid,
    posterior: v3.posterior,
    MLE: {
      headers: { 'BayesFactor[alpha<beta]': '', 'Prob[alpha<beta]': '', 'Prob[alpha>beta]': '', alpha: '', beta: '', 'beta-alpha': '' },
      content: { 0: d.site.map((_, i) => [d.alpha[i], d.beta[i], d.beta_minus_alpha[i], d.prob_neg[i], d.prob_pos[i], d.bayes_factor[i], 0, 0]) }
    },
    timers: v3.timers,
    runtime: '3.0.0'
  };
}

test('every FUBAR fixture parses into the same 187 sites and calls', { skip: noFixtures }, () => {
  const v3 = new FUBARData(fixture('fubar_cd2.json'));
  const legacy = new FUBARData(fixture('fubar_cd2.hyphy25.json'));
  const h3legacy = new FUBARData(hyphy3Legacy(fixture('fubar_cd2.json')));
  assert.equal(v3.format, 'v3');
  assert.equal(legacy.format, 'hyphy25');
  assert.equal(h3legacy.format, 'hyphy3-legacy');
  assert.deepEqual([v3.formatLabel, legacy.formatLabel, h3legacy.formatLabel],
    ['hyphy3 v3', 'HyPhy 2.5 legacy', 'hyphy3 --legacy-json']);

  for (const data of [v3, legacy, h3legacy]) {
    assert.equal(data.numSites, 187);
    assert.equal(data.numTaxa, 10);
    assert.equal(data.siteResults.length, 187);
    assert.deepEqual(data.siteResults.map(r => r.site), Array.from({ length: 187 }, (_, i) => i + 1));
    assert.equal(data.runThreshold, 0.9);
    assert.equal(data.settings.method, 'Variational-Bayes');
    assert.equal(data.settings.gridPoints, 20);
    assert.equal(data.grid.dims, 20);
    assert.equal(data.grid.weight.length, 400);
    assert.equal(data.capabilities.sitePosterior, true);
    assert.equal(data.capabilities.multiPartition, false);
    const stats = data.getStats(0.9);
    assert.deepEqual([stats.diversifying, stats.purifying, stats.neutral], [13, 14, 160]);
    assert.deepEqual(stats.falseCI, [0, 3]);
    assert.equal(configForResult(data.raw), FUBAR_CONFIG);
  }
  // hyphy3 and HyPhy 2.5 differ by up to ~1e-3 per posterior; the calls agree
  assert.deepEqual(v3.siteResults.map(r => fubarSiteClass(r, 0.9)), legacy.siteResults.map(r => fubarSiteClass(r, 0.9)));
  close(legacy.getStats(0.9).expectedFalse, 0.94875, 1e-5);
  close(h3legacy.getStats(0.9).expectedFalse, v3.getStats(0.9).expectedFalse, 1e-12);
});

test('the v3 run summary is reproduced at the run cutoff', { skip: noFixtures }, () => {
  const json = fixture('fubar_cd2.json');
  const data = new FUBARData(json);
  const stats = data.getStats(data.runThreshold);
  const { summary } = json;
  assert.equal(stats.diversifying, summary.positive_sites_count);
  assert.equal(stats.purifying, summary.purifying_sites_count);
  close(stats.expectedFalse, summary.expected_false_positives, 1e-12);
  assert.deepEqual(stats.falseCI, summary.false_positives_ci_95);
  assert.deepEqual(data.siteResults.filter(r => fubarSiteClass(r, 0.9) === 'diversifying').map(r => r.site), summary.positive_sites);
  assert.deepEqual(data.siteResults.filter(r => fubarSiteClass(r, 0.9) === 'purifying').map(r => r.site), summary.purifying_sites);
  assert.deepEqual(data.runSummary, {
    threshold: 0.9, diversifying: 13, purifying: 14, expected: summary.expected_false_positives, ci: [0, 3]
  });
  // The classification column is frozen at the run cutoff and agrees with the recomputed calls
  assert.deepEqual(data.siteResults.map(r => fubarSiteClass(r, 0.9)), json.site_results.data.classification);
});

test('site posteriors reproduce P[β>α] and the EBF in both formats', { skip: noFixtures }, () => {
  for (const name of ['fubar_cd2.json', 'fubar_cd2.hyphy25.json']) {
    const data = new FUBARData(fixture(name));
    const w = data.priorNonPositive;
    for (const r of data.siteResults) {
      const post = data.sitePosterior(r.site);
      assert.equal(post.length, 400);
      let pos = 0;
      post.forEach((x, g) => { if (data.grid.beta[g] > data.grid.alpha[g]) pos += x; });
      close(pos, r.probPos, 1e-12);
      close(r.probEq, 1 - r.probPos - r.probNeg, 1e-15);
      if (r.probPos < 1) close((r.probPos / (1 - r.probPos)) / ((1 - w) / w), r.ebf, 1e-6 * r.ebf);
    }
  }
});

test('the prior mass on β ≤ α is read from v3 and computed from the legacy grid', { skip: noFixtures }, () => {
  const v3 = new FUBARData(fixture('fubar_cd2.json'));
  close(v3.priorNonPositive, 0.5077684277516545, 1e-15);
  close(fubarPriorNonPositive(v3.grid), v3.priorNonPositive, 1e-12);
  const legacy = new FUBARData(fixture('fubar_cd2.hyphy25.json'));
  close(legacy.priorNonPositive, 0.507805, 1e-6);
  close(legacy.priorPositive, 1 - legacy.priorNonPositive, 1e-15);

  const json = fixture('fubar_cd2.json');
  delete json.rate_grid.prior_weight_non_positive;
  close(new FUBARData(json).priorNonPositive, v3.priorNonPositive, 1e-12);
});

test('GTR fit, tree and run details come from schema paths, then the compat block', { skip: noFixtures }, () => {
  const v3 = new FUBARData(fixture('fubar_cd2.json'));
  const legacy = new FUBARData(fixture('fubar_cd2.hyphy25.json'));
  close(v3.gtr.treeLength, 1.693, 1e-3);
  close(legacy.gtr.treeLength, 1.693, 1e-3);
  close(v3.gtr.aicc, 7108.367, 1e-3);
  close(legacy.gtr.aicc, 7112.858, 1e-3);
  assert.equal(v3.gtr.params, 24);
  assert.equal(v3.gtr.rates.length, 6);
  assert.deepEqual(v3.gtr.rates.find(([k]) => k === 'A↔G'), ['A↔G', 1]);
  assert.ok(!v3.treeNewick.includes('{'));
  assert.ok(legacy.treeNewick.startsWith('(((('));
  assert.equal(v3.settings.chains, 1);
  assert.equal(v3.settings.treeScaling, 'grid');
  assert.equal(legacy.settings.treeScaling, undefined);
  assert.equal(v3.citation.doi, '10.1093/molbev/mst030');

  // Without the compat block the v3 schema paths still give the fit's 5 thetas plus A↔G
  const json = fixture('fubar_cd2.json');
  for (const key of ['MLE', 'grid', 'posterior', 'settings', 'fits', 'input', 'timers']) delete json[key];
  const bare = new FUBARData(json);
  assert.equal(bare.siteResults.length, 187);
  assert.equal(bare.capabilities.sitePosterior, false);
  assert.equal(bare.sitePosterior(7), null);
  assert.equal(bare.gtr.rates.length, 6);
  assert.equal(bare.settings.chains, undefined);
  assert.equal(bare.grid.dims, 20);
  assert.match(fubarDiagnostics(bare), /Variational-Bayes/);
  assert.doesNotMatch(fubarDiagnostics(bare), /Chains/);
});

test('HyPhy 2.5 partitions are placed by their coverage', () => {
  const row = (pos, neg) => [1, 2, 1, neg, pos, 3, 0, 0];
  const data = new FUBARData({
    runtime: '2.5.101',
    analysis: { info: 'Perform a Fast Unbiased AppRoximate Bayesian (FUBAR)\n analysis', version: '2.2' },
    input: { 'number of sites': 4, 'number of sequences': 3 },
    'data partitions': { 0: { coverage: [[0, 2]] }, 1: { coverage: [[1, 3]] } },
    grid: [[0, 0, 0.5], [0, 1, 0.25], [1, 0, 0.15], [1, 1, 0.1]],
    posterior: { 0: { 0: [[0, 1, 0, 0]], 1: [[0, 0, 1, 0]] }, 1: { 1: [[1, 0, 0, 0]] } },
    MLE: { headers: [], content: { 0: [row(0.95, 0.01), row(0.1, 0.92)], 1: [row(0.5, 0.4), row(0.3, 0.6)] } },
    settings: { posterior: 0.9, method: 'Collapsed-Gibbs', chains: 2, 'chain-length': 1000, 'burn-in': 100, samples: 50 }
  });
  assert.deepEqual(data.siteResults.map(r => [r.site, r.partition, r.probPos]), [[1, 1, 0.95], [2, 2, 0.5], [3, 1, 0.1], [4, 2, 0.3]]);
  assert.equal(data.capabilities.multiPartition, true);
  assert.deepEqual(data.sitePosterior(1), [0, 1, 0, 0]);
  assert.deepEqual(data.sitePosterior(3), [0, 0, 1, 0]);
  assert.deepEqual(data.sitePosterior(4), [1, 0, 0, 0]);
  assert.equal(data.sitePosterior(2), null);
  // beta <= alpha: cells (0,0), (1,0) and (1,1)
  close(data.priorNonPositive, 0.75);
  assert.equal(data.settings.gridPoints, 2);
  assert.match(fubarDiagnostics(data), /2 × 1,000 steps, burn-in 100, 50 samples/);
});

test('calls are diversifying first, and the false-positive count covers diversifying calls only', () => {
  const rows = [
    { probPos: 0.97, probNeg: 0.01 },
    { probPos: 0.91, probNeg: 0.05 },
    { probPos: 0.02, probNeg: 0.95 },
    { probPos: 0.6, probNeg: 0.3 }
  ];
  assert.deepEqual(rows.map(r => fubarSiteClass(r, 0.9)), ['diversifying', 'diversifying', 'purifying', 'neutral']);
  assert.equal(fubarSiteClass({ probPos: 0.9, probNeg: 0.9 }, 0.9), 'diversifying');
  const efp = fubarFalsePositives(rows, 0.9);
  assert.equal(efp.calls, 2);
  close(efp.expected, 0.12, 1e-12);
  assert.deepEqual(efp.ci, [0, 1]);
  assert.deepEqual(fubarFalsePositives(rows, 0.99), { calls: 0, expected: 0, ci: [0, 0] });
});

test('the posterior cutoff passes at or above t and labels three decimals near 1', () => {
  assert.equal(thresholdPasses(FUBAR_POSTERIOR_CUTOFF, 0.9, 0.9), true);
  assert.equal(thresholdPasses(FUBAR_POSTERIOR_CUTOFF, 0.89, 0.9), false);
  assert.equal(thresholdLabel(FUBAR_POSTERIOR_CUTOFF, 0.9), 'P ≥ 0.90');
  assert.equal(thresholdLabel(FUBAR_POSTERIOR_CUTOFF, 0.999), 'P ≥ 0.999');
  assert.deepEqual(sliderAttributes(FUBAR_POSTERIOR_CUTOFF, 0.95), { min: '0.50', max: '0.999', step: '0.001', value: '0.95' });
});

test('KPI cards follow the cutoff', { skip: noFixtures }, () => {
  const data = new FUBARData(fixture('fubar_cd2.json'));
  const read = (t) => FUBAR_KPI_CARDS.map(card => [card.value(data.getStats(t), { data, threshold: t }), card.sub(data.getStats(t), { data, threshold: t })]);
  assert.deepEqual(read(0.9), [
    ['187 Codons', '10 Taxa'],
    ['▲ 13', 'P[β>α] ≥ 0.90'],
    ['▼ 14', 'P[α>β] ≥ 0.90'],
    ['0.95', '95% CI 0–3 of 13 diversifying calls']
  ]);
  assert.deepEqual(read(0.99).slice(1), [
    ['▲ 0', 'P[β>α] ≥ 0.99'],
    ['▼ 2', 'P[α>β] ≥ 0.99'],
    ['0.00', '95% CI 0–0 of 0 diversifying calls']
  ]);
});

test('an empty class jumps to its mirror, else lists the closest calls', () => {
  const [div, pur] = FUBAR_FILTERS;
  const rows = [{ site: 1, probPos: 0.5, probNeg: 0.2 }, { site: 2, probPos: 0.8, probNeg: 0.1 }];
  const jump = fubarEmptyState(0.9, { filter: div, stats: { diversifying: 0, purifying: 3 }, rows });
  assert.deepEqual(jump.jump, { filter: 'purifying', label: 'View Purifying Sites (▼ 3)' });
  assert.match(jump.title, /^No diversifying site reaches <span class="case-literal">P\[β>α\]<\/span> ≥ 0\.90$/);
  const fallback = fubarEmptyState(0.9, { filter: pur, stats: { diversifying: 0, purifying: 0 }, rows });
  assert.deepEqual(fallback.fallback.map(r => r.site), [1, 2]);
  const closest = fubarEmptyState(0.9, { filter: div, stats: { diversifying: 0, purifying: 0 }, rows });
  assert.deepEqual(closest.fallback.map(r => r.site), [2, 1]);
});
