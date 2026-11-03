import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import {
  GardResults, gardBreakCall, gardSplits, gardRF, gardMovedTaxa, gardFormat, gardP, gardCutoff
} from '../src/data/gard-data.js';
import { detectMethod } from '../src/core/detect-method.js';
import { GARD_KPI_CARDS, GARD_ALPHA, gardTopologyValue, gardTopologySub } from '../src/methods/gard.js';
import { METHOD_CONFIGS } from '../src/methods/registry.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(new URL('gard_adh.json', FIXTURES)) && 'parity fixtures are not in this checkout';

const TREE_A = '((A:0.1,B:0.1)Node1:0.2,(C:0.1,D:0.1)Node2:0.2,E:0.3)Node0;';
const TREE_B = '((A:0.1,C:0.1)Node1:0.2,(B:0.1,D:0.1)Node2:0.2,E:0.3)Node0;';

const test1 = (overrides = {}) => Object.assign({
  breakpoint: 279,
  left_partition: [1, 279],
  right_partition: [280, 561],
  delta_log_l_left: 164.5,
  delta_log_l_right: 161.0,
  kh_p_value_left: 0,
  kh_p_value_right: 0,
  sh_p_value_left: 0,
  sh_p_value_right: 0,
  statistically_significant: true
}, overrides);

// A modern v3 result with one breakpoint after site 279 (0-based position 278)
function modernDoc(overrides = {}) {
  return Object.assign({
    analysis: { id: 'gard', name: 'GARD', version: '3.0.0', citation: 'Mol Biol Evol 23(10), 1891–1901' },
    software: { version: '3.0.0', git_commit: 'abc1234' },
    dataset: { file: '/home/user/data/recomb.fna', taxa_count: 5, sites_count: 561, variable_sites_count: 388, potential_breakpoints: 387, minimum_partition_size: 17 },
    baseline: { aicc: 7975.5, log_likelihood: -3963.5, tree: TREE_A },
    model_selection: {
      best_model_breakpoints: [278],
      best_model_aicc: 7071.5,
      best_model_log_likelihood: -3492.7,
      improvements: [{ breakpoints_count: 1, breakpoints: [278], delta_aicc: 904 }]
    },
    partitions: [
      { partition_index: 0, start_site: 1, end_site: 279, log_likelihood: -1720.9, tree: TREE_A },
      { partition_index: 1, start_site: 280, end_site: 561, log_likelihood: -1771.7, tree: TREE_B }
    ],
    site_breakpoint_support: { 100: 0.05, 278: 0.9977, 400: 1e-320 },
    topology_test: {
      single_tree_aicc: 7216.6,
      delta_aicc_topology: 145,
      topological_incongruence_supported: true,
      test_type: 'both',
      bootstrap_replicates: 1000,
      any_incongruence_significant: true,
      breakpoints: [test1()],
      partition_tree_matrix: [
        { partition_index: 0, interval: [1, 279], best_tree_index: 0, tree_log_likelihoods: [-1720.9, -1885.5], kh_p_values: [1, 0], sh_p_values: [1, 0] },
        { partition_index: 1, interval: [280, 561], best_tree_index: 1, tree_log_likelihoods: [-1932.7, -1771.7], kh_p_values: [0, 1], sh_p_values: [0, 1] }
      ]
    },
    ohmic_network: {
      bridges_detected: true, num_bridges: 2, ohmic_log_likelihood: -3960.5, ohmic_aicc: 7969.2, delta_aicc_vs_baseline: 5.9,
      bridge_candidate_breakpoints: [262, 17, 222], bridges: []
    },
    rate_variation: { mode: 'discrete_gamma', classes: 4, shape_alpha: 1.8 },
    tree_engine: { method: 'ohmic_spectral', complexity: 'O(N log N)' },
    total_models_evaluated: 780,
    time_elapsed: 1.5
  }, overrides);
}

function noBreakDoc() {
  return modernDoc({
    model_selection: { best_model_breakpoints: [], best_model_aicc: 7975.5, best_model_log_likelihood: -3963.5, improvements: [] },
    partitions: [{ partition_index: 0, start_site: 1, end_site: 561, log_likelihood: -3963.5, tree: TREE_A }],
    topology_test: { single_tree_aicc: 7975.5, delta_aicc_topology: 0, topological_incongruence_supported: false, test_type: 'both', bootstrap_replicates: 1000, any_incongruence_significant: false, breakpoints: [], partition_tree_matrix: [] }
  });
}

// An old hyphy3 writer: the AICc verdict only, no partitions, tests, matrix or Ohmic block
function oldV3Doc() {
  const doc = modernDoc();
  for (const key of ['ohmic_network', 'rate_variation', 'tree_engine', 'provenance']) delete doc[key];
  doc.topology_test = { single_tree_aicc: 7216.6, delta_aicc_topology: 145, topological_incongruence_supported: true };
  return doc;
}

function legacyDoc(bps, extra = {}) {
  return Object.assign({
    analysis: { info: 'GARD : Genetic Algorithms for Recombination Detection.', version: '0.2', citation: '**Automated Detection**, _Mol Biol Evol 23(10)' },
    baselineScore: 7975.5,
    bestModelAICc: 7071.5,
    singleTreeAICc: 7216.6,
    breakpointData: { 0: { bps: [bps[0]], tree: TREE_A }, 1: { bps: [bps[1]], tree: TREE_B } },
    improvements: { 0: { breakpoints: [278], deltaAICc: 904 } },
    input: { 'file name': '/x/recomb.fna', 'number of sequences': 5, 'number of sites': 561 },
    potentialBreakpoints: 387,
    siteBreakPointSupport: { 278: 0.9977 },
    timeElapsed: 2,
    totalModelCount: 780,
    trees: {}
  }, extra);
}

test('a recombinant result maps coordinates: support key 278 is the end of partition 1-279', () => {
  const data = new GardResults(modernDoc());
  assert.equal(data.off, 1);
  assert.deepEqual(data.segments.map(s => [s.key, s.start, s.end]), [['S0', 1, 279], ['S1', 280, 561]]);
  assert.deepEqual(data.support.map(s => s.site), [101, 279, 401]);
  assert.equal(data.support.find(s => s.site === 279).p, 0.9977);
  assert.equal(data.breaks.length, 1);
  assert.equal(data.breaks[0].site, 279);
  assert.equal(data.breaks[0].test.khLeft, 0);
  assert.equal(data.deltaAicc, 7975.5 - 7071.5);
  assert.deepEqual([data.numSites, data.numTaxa], [561, 5]);
  assert.deepEqual(data.capabilities, { topology: true, matrix: true, ohmic: true, legacy: false });
  assert.equal(data.runThreshold, 0.05);
  assert.equal(data.treeNewick, TREE_A);
  assert.deepEqual(data.steps, [{ k: 1, deltaAicc: 904, bps: [279] }]);
  assert.equal(data.matrix[1].bestTree, 1);
  assert.equal(data.getBreak('279'), data.breaks[0]);
  assert.deepEqual(data.getStats(0.05), { breakpoints: 1, incongruent: 1, branchOnly: 0, untested: 0, segments: 2, totalSites: 561 });
});

test('unsorted Ohmic candidates are sorted and shifted to display sites', () => {
  const data = new GardResults(modernDoc());
  assert.deepEqual(data.ohmic.candidates, [18, 223, 263]);
  assert.equal(data.ohmic.count, 2);
  assert.equal(data.ohmic.deltaAicc, 5.9);
});

test('no breakpoints is one segment, no calls, and the support and Ohmic blocks still read', () => {
  const data = new GardResults(noBreakDoc());
  assert.equal(data.segments.length, 1);
  assert.equal(data.off, 1);
  assert.deepEqual(data.breaks, []);
  assert.equal(data.deltaAicc, 0);
  assert.deepEqual(data.getStats(0.05), { breakpoints: 0, incongruent: 0, branchOnly: 0, untested: 0, segments: 1, totalSites: 561 });
  assert.equal(data.capabilities.topology, false);
  assert.equal(data.capabilities.ohmic, true);
  assert.equal(data.runThreshold, null);
  assert.equal(data.support.length, 3);
});

test('an old v3 file has no tests, matrix or Ohmic block: every break is untested, the verdict stays', () => {
  const data = new GardResults(oldV3Doc());
  assert.deepEqual(data.capabilities, { topology: false, matrix: false, ohmic: false, legacy: false });
  assert.equal(data.breaks[0].test, null);
  assert.equal(data.callOf(data.breaks[0], 0.05), 'untested');
  assert.equal(data.getStats(0.05).untested, 1);
  assert.deepEqual([data.topology.supported, data.topology.deltaAicc], [true, 145]);

  // Without partitions the baseline tree stands for the one segment
  const bare = modernDoc();
  delete bare.partitions;
  const single = new GardResults(bare);
  assert.deepEqual(single.segments.map(s => [s.start, s.end, s.tree]), [[1, 561, TREE_A]]);
});

test('hyphy3-legacy string bounds and HyPhy 2.5 numeric bounds give the same segments', () => {
  const strings = new GardResults(legacyDoc([['1', '279'], ['280', '561']]));
  const numbers = new GardResults(legacyDoc([[1, 279], [280, 561]], {
    analysis: { info: 'GARD : Genetic Algorithms\n    for Recombination Detection.', version: '0.2' }
  }));
  assert.equal(strings.format, 'hyphy3-legacy');
  assert.equal(numbers.format, 'hyphy25');
  assert.deepEqual(strings.segments, numbers.segments);
  assert.deepEqual(strings.segments.map(s => [s.start, s.end]), [[1, 279], [280, 561]]);
  assert.equal(strings.capabilities.legacy, true);
  assert.equal(strings.off, 1);
  assert.equal(strings.support.find(s => s.site === 279).p, 0.9977);
  assert.deepEqual(strings.steps, [{ k: 1, deltaAicc: 904, bps: [279] }]);
  assert.equal(strings.topology.deltaAicc, 7216.6 - 7071.5);
  assert.equal(strings.numTaxa, 5);
  // No topology tests in the file: untested, not branch-only
  assert.equal(strings.callOf(strings.breaks[0], 0.05), 'untested');
});

test('legacy topology tests read like the v3 ones; a legacy no-breakpoint file keeps its one tree', () => {
  const tests = { testType: 'kh', replicates: 100, anySignificant: true, breakpoints: [test1()] };
  const data = new GardResults(legacyDoc([['1', '279'], ['280', '561']], { topologyTests: tests }));
  assert.equal(data.capabilities.topology, true);
  assert.equal(data.testMode, 'kh');
  assert.equal(data.callOf(data.breaks[0], 0.05), 'incongruent');

  const none = new GardResults({
    analysis: { info: 'GARD : Genetic Algorithms' },
    baselineScore: 100, bestModelAICc: 100, singleTreeAICc: 100,
    breakpointData: {}, improvements: { 0: { breakpoints: null, deltaAICc: 0 } },
    input: { 'number of sequences': 5, 'number of sites': 300 },
    trees: { newickString: TREE_A }
  });
  assert.deepEqual(none.segments.map(s => [s.start, s.end, s.tree]), [[1, 300, TREE_A]]);
  assert.deepEqual(none.steps, []);
  assert.equal(none.deltaAicc, 0);
});

test('the live call is the writer\'s Bonferroni rule, per flank at alpha / 2B', () => {
  const t = (khL, khR, shL, shR) => ({ khLeft: khL, khRight: khR, shLeft: shL, shRight: shR });
  // B = 1: the cutoff at alpha 0.05 is 0.025, the writer's bonf_alpha
  assert.equal(gardBreakCall(t(0.02, 0.02, 0.5, 0.5), 0.05, 1, 'both'), 'incongruent');
  assert.equal(gardBreakCall(t(0.02, 0.02, 0.5, 0.5), 0.05, 1, 'kh'), 'incongruent');
  assert.equal(gardBreakCall(t(0.02, 0.02, 0.5, 0.5), 0.05, 1, 'sh'), 'branch');
  assert.equal(gardBreakCall(t(0.02, 0.03, 0.01, 0.01), 0.05, 1, 'kh'), 'branch');
  assert.equal(gardBreakCall(t(0.02, 0.03, 0.01, 0.01), 0.05, 1, 'both'), 'incongruent');
  // Both flanks must reject
  assert.equal(gardBreakCall(t(0.001, 0.2, 0.001, 0.2), 0.05, 1, 'both'), 'branch');
  // More breakpoints tighten the cutoff
  assert.equal(gardBreakCall(t(0.02, 0.02, 1, 1), 0.05, 2, 'both'), 'branch');
  assert.equal(gardBreakCall(t(0.02, 0.02, 1, 1), 0.2, 2, 'both'), 'incongruent');
  assert.equal(gardBreakCall(null, 0.05, 1), 'untested');
  assert.equal(gardBreakCall(t(null, 0, 0, 0), 0.05, 1, 'kh'), 'branch');

  // At alpha 0.05 it agrees with statistically_significant as the writer stored it
  for (const [kh, sh, mode, stored] of [[0.02, 0.5, 'both', true], [0.02, 0.5, 'sh', false], [0.2, 0.2, 'both', false], [0.5, 0.01, 'kh', false]]) {
    const doc = modernDoc();
    doc.topology_test.test_type = mode;
    Object.assign(doc.topology_test.breakpoints[0], {
      kh_p_value_left: kh, kh_p_value_right: kh, sh_p_value_left: sh, sh_p_value_right: sh, statistically_significant: stored
    });
    const data = new GardResults(doc);
    assert.equal(data.callOf(data.breaks[0], 0.05) === 'incongruent', stored, `${mode} kh ${kh} sh ${sh}`);
  }
});

test('RF and moved taxa read leaf-name splits, ignoring NodeN labels and the root', () => {
  assert.deepEqual(gardRF(TREE_A, TREE_A), { rf: 0, max: 4 });
  assert.deepEqual(gardRF(TREE_A, TREE_B), { rf: 4, max: 4 });
  assert.equal(gardSplits(TREE_A).splits.size, 2);
  // Renaming the internal nodes, or rerooting on E, changes nothing
  assert.equal(gardRF(TREE_A, '((A,B)X,(C,D)Y,E)Z;').rf, 0);
  assert.equal(gardRF(TREE_A, '(E,((A,B),(C,D)));').rf, 0);
  // One leaf moved: A leaves (A,B) for (C,D)
  const moved = '(B,(C,D,A)N,E);';
  assert.deepEqual(gardMovedTaxa(TREE_A, moved), ['A', 'B', 'C', 'D']);
  assert.deepEqual(gardMovedTaxa(TREE_A, TREE_A), []);
  assert.deepEqual(gardMovedTaxa(TREE_A, TREE_B), ['A', 'B', 'C', 'D']);
  // Different taxa sets compare on what they share
  assert.equal(gardRF(TREE_A, '((A,B),(C,D),F);').rf, 0);
  assert.deepEqual(gardMovedTaxa(TREE_A, '((A,B),(C,D),F);'), []);
});

test('numbers print as dashes when null, NaN or infinite', () => {
  assert.equal(gardFormat(null), '–');
  assert.equal(gardFormat(undefined), '–');
  assert.equal(gardFormat(Infinity), '–');
  assert.equal(gardFormat(NaN), '–');
  assert.equal(gardFormat(3.14159, 1), '3.1');
  assert.equal(gardP(null), '–');
  assert.equal(gardP(0), '0');
  assert.equal(gardP(0.00012), '1.2e-4');
  assert.equal(gardP(0.5), '0.5000');
  assert.equal(gardCutoff(0.05 / 2), '0.025');
  const data = new GardResults(modernDoc({ baseline: { aicc: null, tree: TREE_A } }));
  assert.equal(data.deltaAicc, null);
  assert.equal(GARD_KPI_CARDS[1].sub(data.getStats(0.05), { data }), 'ΔAICc – vs single tree');
});

test('the KPI cards read the result at the live alpha', () => {
  const view = (data, threshold = 0.05) => ({ data, threshold, stats: data.getStats(threshold) });
  const text = (card, v) => [card.value(v.stats, v), card.sub(v.stats, v)];
  const [alignment, breaks, topology, ohmic] = GARD_KPI_CARDS;
  const data = new GardResults(modernDoc());
  const v = view(data);
  assert.deepEqual(text(alignment, v), ['5 taxa · 561 sites', '388 variable · 387 candidate breakpoints']);
  assert.deepEqual(text(breaks, v), ['1 breakpoint', 'ΔAICc 904.0 vs single tree']);
  assert.deepEqual(text(topology, v), ['1 of 1 incongruent', 'ΔAICc_topology 145.0 (topology, not just branch lengths)']);
  assert.equal(topology.tone(v.stats, v), 'accent');
  assert.deepEqual(text(ohmic, v), ['2 bridges', '3 candidate breakpoints · ΔAICc vs baseline 5.9']);
  assert.equal(ohmic.available(data), true);

  // Branch-length heterogeneity: the AICc verdict is negative
  const rate = modernDoc();
  rate.topology_test.topological_incongruence_supported = false;
  Object.assign(rate.topology_test.breakpoints[0], { kh_p_value_left: 0.4, sh_p_value_left: 0.4, statistically_significant: false });
  const r = view(new GardResults(rate));
  assert.deepEqual(text(topology, r), ['0 of 1 incongruent', 'Branch-length heterogeneity only · ΔAICc_topology 145.0']);
  assert.equal(topology.tone(r.stats, r), 'ink');

  const none = view(new GardResults(noBreakDoc()));
  assert.deepEqual(text(breaks, none), ['No recombination detected', 'A single tree fits best (ΔAICc 0)']);
  assert.equal(breaks.tone(none.stats, none), 'muted');
  assert.deepEqual(text(topology, none), ['Not applicable', 'No breakpoints to test']);

  const old = view(new GardResults(oldV3Doc()));
  assert.deepEqual(text(topology, old), ['Supported', 'Topology tests not in this file · ΔAICc_topology 145.0']);
  assert.equal(topology.tone(old.stats, old), 'muted');
  assert.equal(ohmic.available(old.data), false);
  assert.equal(gardTopologyValue(old.stats, old), 'Supported');
  assert.match(gardTopologySub(old.stats, old), /not in this file/);
});

test('the alpha controller is a log slider that applies only where tests exist', () => {
  assert.equal(GARD_ALPHA.scale, 'log');
  assert.equal(GARD_ALPHA.seed, 0.05);
  const data = new GardResults(modernDoc());
  assert.equal(GARD_ALPHA.enabled(data), true);
  assert.equal(GARD_ALPHA.enabled(new GardResults(oldV3Doc())), false);
  assert.equal(GARD_ALPHA.enabled(new GardResults(noBreakDoc())), false);
  assert.equal(GARD_ALPHA.disabledCaption, 'No topology test');
  assert.equal(GARD_ALPHA.caption({ threshold: 0.05, data }), 'Bonferroni per test α/2B = 0.025, B = 1');
  assert.equal(GARD_ALPHA.caption({ threshold: 0.05, data: null }), 'Bonferroni per test α/2B');
});

test('GARD is detected in all three layouts and has a config', () => {
  assert.equal(detectMethod(modernDoc()).method, 'gard');
  assert.equal(detectMethod(legacyDoc([['1', '279'], ['280', '561']])).method, 'gard');
  assert.equal(detectMethod({ breakpointData: {} }).method, 'gard');
  assert.equal(METHOD_CONFIGS.gard.Parser, GardResults);
});

test('the parity fixtures parse', { skip: noFixtures }, () => {
  for (const [file, format] of [['gard_adh.json', 'v3'], ['gard_adh.legacy.json', 'hyphy3-legacy'], ['gard_adh.hyphy25.json', 'hyphy25']]) {
    const data = new GardResults(JSON.parse(readFileSync(new URL(file, FIXTURES), 'utf8')));
    assert.equal(data.format, format, file);
    assert.ok(data.segments.length >= 1, file);
    assert.equal(data.breaks.length, data.segments.length - 1, file);
    assert.ok(data.numSites > 0 && data.numTaxa > 0, file);
    assert.ok(data.support.length > 0, file);
    assert.equal(data.segments.at(-1).end, data.numSites, file);
  }
});
