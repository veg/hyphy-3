import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import { SLACData, slacSiteClass, consensusResidue } from '../src/data/slac-data.js';
import { configForResult } from '../src/methods/registry.js';
import {
  SLAC_CONFIG, SLAC_KPI_CARDS, SLAC_FILTERS, SLAC_TABLE, SLAC_COLUMNS, SLAC_P_THRESHOLD,
  SLAC_DRAWER, slacEmptyState, slacDiagnostics, slacCountingSentence, slacFitsTable, slacP, slacEffectBar,
  slacLeafLabel, slacShares, slacTipPValues, slacNewickText
} from '../src/methods/slac.js';
import { thresholdLabel, thresholdCaption, sliderAttributes } from '../src/components/kpi-bar.js';
import { formatThreshold } from '../src/core/format.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';
const fixture = name => JSON.parse(readFileSync(new URL(name, FIXTURES), 'utf8'));
const FILES = ['slac_cd2.json', 'slac_cd2.hyphy25.json', 'slac_cd2.legacy.json'];

const close = (actual, expected, tol = 1e-9) =>
  assert.ok(Math.abs(actual - expected) <= tol, `${actual} vs ${expected}`);

const PURIFYING_AT_01 = [47, 65, 78, 82, 87, 110, 116, 123, 130, 164, 166];
const NO_SUBSTITUTIONS = [2, 3, 35, 38, 39, 66, 73, 108, 114, 168, 181, 185, 186];

test('every SLAC fixture parses into the same 187 sites and calls', { skip: noFixtures }, () => {
  const [v3, hyphy25, legacy] = FILES.map(f => new SLACData(fixture(f)));
  assert.deepEqual([v3.format, hyphy25.format, legacy.format], ['v3', 'hyphy25', 'hyphy3-legacy']);
  assert.deepEqual([v3.formatLabel, hyphy25.formatLabel, legacy.formatLabel],
    ['hyphy3 v3', 'HyPhy 2.5 legacy', 'hyphy3 --legacy-json']);

  for (const data of [v3, hyphy25, legacy]) {
    assert.equal(data.numSites, 187);
    assert.equal(data.numTaxa, 10);
    assert.deepEqual(data.siteResults.map(r => r.site), Array.from({ length: 187 }, (_, i) => i + 1));
    const stats = data.getStats(0.1);
    assert.deepEqual([stats.diversifying, stats.purifying, stats.neutral, stats.invariant], [0, 11, 176, 13]);
    assert.deepEqual(data.siteResults.filter(r => slacSiteClass(r, 0.1) === 'purifying').map(r => r.site), PURIFYING_AT_01);
    assert.deepEqual(data.siteResults.filter(r => r.noSubs).map(r => r.site), NO_SUBSTITUTIONS);
    // dS is undefined (never 0) where ES = 0
    assert.deepEqual(data.siteResults.filter(r => r.dSUndefined).map(r => r.site), [35, 114]);
    assert.equal(data.getSite(35).dNmdS, null);
    assert.equal(data.capabilities.sampling, true);
    assert.equal(data.capabilities.multiPartition, false);
    assert.equal(data.capabilities.leafCodons, true);
    assert.equal(Object.keys(data.leafCodons).length, 10);
    assert.equal(data.testing.tested, 16);
    assert.equal(configForResult(data.raw), SLAC_CONFIG);
    // p_positive + p_negative >= 1 everywhere, so no site can pass both below 0.5
    assert.ok(data.siteResults.every(r => r.pPos + r.pNeg >= 1 - 1e-12));
  }
  assert.equal(v3.getStats(0.01).purifying, 1);
});

test('v3 and hyphy3 --legacy-json carry the same resolved per-site values', { skip: noFixtures }, () => {
  const v3 = new SLACData(fixture('slac_cd2.json'));
  const legacy = new SLACData(fixture('slac_cd2.legacy.json'));
  for (const key of ['ES', 'EN', 'S', 'N', 'P_S', 'dN', 'pPos', 'pNeg', 'tbl']) {
    v3.siteResults.forEach((r, i) => close(legacy.siteResults[i][key], r[key], 1e-8));
  }
  // Site 7: S 1.5, N 8.5; v3 and legacy both read the RESOLVED table (AVERAGED differs at 70 and 130)
  assert.deepEqual([v3.getSite(7).S, v3.getSite(7).N], [1.5, 8.5]);
  assert.equal(legacy.getSite(70).N, v3.getSite(70).N);
  close(legacy.getSite(130).pNeg, v3.getSite(130).pNeg, 1e-8);
});

test('totals, omega, settings and the run summary come from each format\'s own keys', { skip: noFixtures }, () => {
  const [v3, hyphy25, legacy] = FILES.map(f => new SLACData(fixture(f)));
  assert.equal(v3.totalS.toFixed(1), '193.8');
  assert.equal(v3.totalN.toFixed(1), '545.2');
  assert.equal(v3.globalOmega.toFixed(3), '0.993');
  assert.equal(hyphy25.globalOmega.toFixed(3), '1.017');
  assert.equal(legacy.globalOmega.toFixed(3), '0.993');

  assert.equal(v3.runThreshold, 0.1);
  assert.equal(legacy.runThreshold, 0.1);
  assert.equal(hyphy25.runThreshold, null);
  assert.deepEqual(v3.runSummary, { threshold: 0.1, diversifying: 0, purifying: 11, neutral: 176 });
  assert.equal(legacy.runSummary, null);
  // The v3 counts frozen at the run cutoff agree with the live recount there
  const stats = v3.getStats(v3.runThreshold);
  assert.deepEqual([stats.diversifying, stats.purifying, stats.neutral], [0, 11, 176]);
  v3.siteResults.forEach((r, i) => {
    const frozen = v3.raw.site_results.data.classification[i];
    const live = slacSiteClass(r, 0.1);
    assert.equal(live === 'invariant' ? 'Neutral' : live[0].toUpperCase() + live.slice(1), frozen);
  });

  assert.deepEqual(v3.settings, { samples: 100, testBranches: 'All', geneticCode: 'Universal' });
  assert.deepEqual(legacy.settings, { samples: 100, testBranches: 'All', geneticCode: 'Universal' });
  assert.deepEqual(hyphy25.settings, { samples: null, testBranches: null, geneticCode: null });

  // Fits keep each producer's spelling: hyphy3 --legacy-json writes 'Log-likelihood' and no count
  assert.deepEqual(v3.fits.map(f => [f.name, f.params]), [['Global MG94xREV', 22], ['Nucleotide GTR', 21]]);
  assert.deepEqual(hyphy25.fits.map(f => f.params), [31, 24]);
  assert.deepEqual(legacy.fits.map(f => f.params), [null, null]);
  close(legacy.fits[0].logL, -3466.999707318957);
});

test('trees take MG94 lengths from node metadata, with {tags} stripped', { skip: noFixtures }, () => {
  for (const file of FILES) {
    const data = new SLACData(fixture(file));
    assert.ok(!data.treeNewick.includes('{'), file);
    assert.equal(data.capabilities.tree, true);
    assert.ok(data.nodeMeta.Human.branch_lengths['Global MG94xREV'] >= 0);
  }
  const v3 = new SLACData(fixture('slac_cd2.json'));
  close(Object.values(v3.nodeMeta).reduce((s, m) => s + m.branch_lengths['Global MG94xREV'], 0), 1.723261353786068, 1e-9);
});

test('sampled ancestors give S, N, dS and dN ranges; p-value ranges are not read', { skip: noFixtures }, () => {
  for (const file of FILES) {
    const r = new SLACData(fixture(file)).getSite(7);
    assert.deepEqual(Object.keys(r.sampled), ['S', 'N', 'dS', 'dN']);
    assert.ok(r.sampled.N.lo <= r.sampled.N.median && r.sampled.N.median <= r.sampled.N.hi);
  }
  const v3 = new SLACData(fixture('slac_cd2.json'));
  assert.deepEqual(v3.getSite(7).sampled.S, { median: 1.5, lo: 1 / 6, hi: 4 });
  assert.equal(v3.getSite(35).sampled.dS, null);
});

test('where dS is undefined its sampled range is blank in every format, not --legacy-json\'s zeros', { skip: noFixtures }, () => {
  for (const file of FILES) {
    const data = new SLACData(fixture(file));
    for (const site of [35, 114]) {
      const r = data.getSite(site);
      assert.ok(r.dSUndefined, `${file} ${site}`);
      assert.equal(r.sampled.dS, null, `${file} ${site}`);
      assert.ok(r.sampled.dN, `${file} ${site}`);
    }
  }
});

test('the consensus residue ignores gaps and breaks ties alphabetically', { skip: noFixtures }, () => {
  assert.equal(consensusResidue(['GTC', 'ACC', 'ACC', 'GTC', '---', 'NNN']), 'T');
  assert.equal(consensusResidue(['---']), null);
  const consensus = FILES.map(f => new SLACData(fixture(f)).siteResults.map(r => r.consensusAA).join(''));
  assert.equal(consensus[0].slice(0, 12), 'KEIRNALETWGA');
  assert.deepEqual(consensus.slice(1), [consensus[0], consensus[0]]);
});

// A HyPhy 2.5-shaped file with two partitions, and a degenerate site: ES = 0 with N > 0
function twoPartitionLegacy() {
  const row = (ES, EN, S, N, pPos, pNeg) => [ES, EN, S, N, ES / (ES + EN), ES ? S / ES : null, N / EN, null, pPos, pNeg, 1];
  return {
    analysis: { info: 'SLAC (Single Likelihood Ancestor Counting)\n uses ...', version: '2.00' },
    runtime: '2.5.101',
    input: { 'number of sequences': 2, 'number of sites': 3, trees: { 0: '(A,B)', 1: '(A,B)' } },
    'data partitions': { 0: { coverage: [[0, 2]] }, 1: { coverage: [[1]] } },
    MLE: {
      headers: [['ES'], ['EN'], ['S'], ['N'], ['P[S]'], ['dS'], ['dN'], ['dN-dS'], ['P [dN/dS > 1]'], ['P [dN/dS < 1]'], ['Total branch length']],
      content: {
        0: { 'by-site': { RESOLVED: [row(1, 2, 0, 3, 0.04, 1), row(1, 2, 0, 0, 1, 1)] } },
        1: { 'by-site': { RESOLVED: [row(0, 3, 0, 2, 0, 1)] } }
      }
    },
    'branch attributes': {
      0: { A: { codon: [['AAA', 'GGG']] }, B: { codon: [['AAG', '---']] } },
      1: { A: { codon: [['TTT']] }, B: { codon: [['TTC']] } }
    },
    fits: {}
  };
}

test('legacy partitions are placed by coverage; a degenerate call is kept, not reclassified', () => {
  const data = new SLACData(twoPartitionLegacy());
  assert.equal(data.format, 'hyphy25');
  assert.deepEqual(data.siteResults.map(r => [r.site, r.partition]), [[1, 1], [2, 2], [3, 1]]);
  assert.equal(data.capabilities.multiPartition, true);
  assert.equal(data.capabilities.sampling, false);
  assert.deepEqual(data.siteCodons(2).map(c => c.codon), ['TTT', 'TTC']);
  assert.equal(data.getSite(3).consensusAA, 'G');
  const degenerate = data.getSite(2);
  assert.equal(degenerate.dSUndefined, true);
  assert.equal(slacSiteClass(degenerate, 0.1), 'diversifying');
  assert.equal(slacSiteClass(data.getSite(3), 0.1), 'invariant');
  assert.match(slacCountingSentence(degenerate, 0.1), /degenerate/);
  assert.equal(slacP(0), '0');
  assert.ok(SLAC_COLUMNS.find(c => c.id === 'class').cell(degenerate, { threshold: 0.1 }).includes('slac-flag'));
});

test('the controller is the shared p cutoff with the uncorrected caption', () => {
  assert.equal(SLAC_CONFIG.controller, SLAC_P_THRESHOLD);
  assert.equal(thresholdLabel(SLAC_P_THRESHOLD, 0.1), 'p ≤ 0.10');
  assert.deepEqual(SLAC_P_THRESHOLD.presets, [0.1, 0.05, 0.01]);
  assert.equal(sliderAttributes(SLAC_P_THRESHOLD, 0.1).max, '0.20');
  assert.equal(thresholdCaption(SLAC_P_THRESHOLD, 0.1, 'p'), 'Uncorrected per-site p-values');
});

// The browser only lets a range input land on min + k·step; each such stop must be the
// cutoff its label prints, and the presets and the maximum must be stops
test('every slider stop is the cutoff its label prints, and the presets and maximum are stops', () => {
  const { min, max, step, presets } = SLAC_P_THRESHOLD;
  const stops = [];
  for (let k = 0; min + k * step <= max + 1e-9; k++) stops.push(Number((min + k * step).toFixed(3)));
  for (const t of stops) assert.equal(Number(formatThreshold(t)), t, `stop ${t}`);
  assert.equal(stops[stops.length - 1], max);
  for (const p of presets) assert.ok(stops.includes(p), `preset ${p}`);
  assert.deepEqual(sliderAttributes(SLAC_P_THRESHOLD, 0.1), { min: '0.005', max: '0.20', step: '0.005', value: '0.10' });
  assert.equal(thresholdLabel(SLAC_P_THRESHOLD, 0.015), 'p ≤ 0.015');
});

test('p-values switch to exponential below 1e-3; shares add up to 100%', () => {
  assert.equal(slacP(0.0005), '5.00e-4');
  assert.equal(slacP(0.001), '0.0010');
  assert.equal(slacP(0.0086), '0.0086');
  assert.equal(slacShares(0.625), 'S 63% · N 37%');
  assert.equal(slacShares(0.2), 'S 20% · N 80%');
});

test('the tooltip puts both p-values on one row and colours the one that passes', () => {
  const r = { pPos: 1, pNeg: 0.0086, S: 3, N: 0, noSubs: false };
  assert.equal(slacTipPValues(r, 0.1), '1.0000 / <span class="tip-pass-pur">0.0086</span>');
  assert.equal(slacTipPValues(r, 0.005), '1.0000 / 0.0086');
  assert.equal(slacTipPValues({ pPos: 0.03, pNeg: 0.99 }, 0.05), '<span class="tip-pass-div">0.0300</span> / 0.9900');
});

test('the drawer says a failing p-value is not significant, and invariant sites get the short badge', { skip: noFixtures }, () => {
  const data = new SLACData(fixture('slac_cd2.json'));
  const render = (site, threshold) => {
    const panel = { innerHTML: '' };
    SLAC_DRAWER.tabs[0].render(panel, data.getSite(site), { data, threshold });
    return panel.innerHTML;
  };
  const purifying = render(87, 0.1);
  assert.doesNotMatch(purifying, /below the cutoff/);
  assert.match(purifying, /1\.0000<\/div>\s*<div class="dstat-sub">not significant at p ≤ 0\.10</);
  assert.match(purifying, /0\.0086<\/div>\s*<div class="dstat-sub">purifying at this cutoff</);
  assert.equal((render(1, 0.1).match(/not significant at p ≤ 0\.10/g) || []).length, 2);
  assert.match(render(47, 0.1), /S 63% · N 37%/);
  const header = (site) => SLAC_DRAWER.header(data.getSite(site), { data, threshold: 0.1 });
  assert.equal(header(2).badge, '○ Invariant');
  assert.equal(header(35).badge, '○ Invariant · ⚠ dS undefined');
  assert.equal(header(87).badge, '▼ Purifying selection');
});

test('codon and selection columns cannot be hidden', () => {
  const pinned = SLAC_COLUMNS.filter(c => c.alwaysShown).map(c => c.id);
  assert.deepEqual(pinned, ['codon', 'class']);
});

test('KPI cards: Alignment first with the substitution totals, the controller after the context card', { skip: noFixtures }, () => {
  const data = new SLACData(fixture('slac_cd2.json'));
  const view = { data, threshold: 0.1 };
  const stats = data.getStats(0.1);
  assert.deepEqual(SLAC_KPI_CARDS.map(c => c.id), ['alignment', 'div', 'pur', 'omega']);
  const [alignment, div, pur, omega] = SLAC_KPI_CARDS;
  assert.equal(alignment.value(stats, view), '187 Codons');
  assert.equal(alignment.sub(stats, view), '10 Taxa · ΣS 193.8 · ΣN 545.2');
  assert.equal(div.value(stats, view), '▲ 0');
  assert.equal(pur.value(stats, view), '▼ 11');
  assert.equal(pur.sub(stats, view), '5.9% of sites (dN < dS)');
  assert.equal(omega.value(stats, view), '0.993');
});

test('the empty diversifying tab explains itself, jumps to purifying and lists the closest calls', { skip: noFixtures }, () => {
  const data = new SLACData(fixture('slac_cd2.json'));
  const [div, pur, sig] = SLAC_FILTERS;
  const state = slacEmptyState(0.1, { filter: div, stats: data.getStats(0.1), rows: data.siteResults });
  assert.deepEqual(state.jump, { filter: 'purifying', label: 'View Purifying Sites (▼ 11)' });
  assert.match(state.title, /p\(dN&gt;dS\)/);
  assert.match(state.desc, /^No codon has significantly more non-synonymous substitutions/);
  assert.equal(state.fallback.length, 5);
  assert.ok(state.fallback.every(r => r.pPos > 0.1));
  const strict = slacEmptyState(0.001, { filter: sig, stats: data.getStats(0.001), rows: data.siteResults });
  assert.equal(strict.jump, undefined);
  assert.equal(strict.desc, undefined);
  assert.deepEqual(strict.fallback.map(r => r.site), [87, 166, 164, 65, 82]);
  assert.equal(pur.sort.column, 'pneg');
  assert.equal(SLAC_TABLE.defaultFilter, 'diversifying');
});

test('table cells: one-decimal counts, dashes for undefined rates, hue only when called', { skip: noFixtures }, () => {
  const data = new SLACData(fixture('slac_cd2.json'));
  const view = { threshold: 0.1, data };
  const cell = (id, site) => SLAC_COLUMNS.find(c => c.id === id).cell(data.getSite(site), view);
  assert.equal(cell('sn', 7), '<td class="font-mono">1.5 / 8.5</td>');
  assert.equal(cell('ds', 35), '<td class="font-mono">—</td>');
  assert.equal(cell('pneg', 87), '<td class="font-mono font-bold text-blue">0.0086</td>');
  assert.equal(cell('pneg', 7), '<td class="font-mono ">0.9110</td>');
  assert.match(cell('class', 2), /Invariant/);
  assert.ok(SLAC_COLUMNS.find(c => c.id === 'esen').hiddenByDefault);
  assert.match(slacEffectBar(data.getSite(87), 0.1, data.effectExtent), /effect-pur/);
  assert.doesNotMatch(slacEffectBar(data.getSite(87), 0.001, data.effectExtent), /effect-pur/);
  assert.equal(slacEffectBar(data.getSite(35), 0.1, data.effectExtent), '<span class="effect-bar"></span>');
  assert.ok(SLAC_TABLE.search.match(data.getSite(87), 'y'));
});

test('the counting sentence quotes the p of the direction the counts lean', { skip: noFixtures }, () => {
  const data = new SLACData(fixture('slac_cd2.json'));
  const purifying = slacCountingSentence(data.getSite(87), 0.1);
  assert.match(purifying, /^Of 3\.0 inferred substitutions, 3\.0 are synonymous \(100%\); 20% would be expected if neutral/);
  assert.match(purifying, /p\(dN&lt;dS\) = <strong class="text-blue">0\.0086<\/strong>/);
  assert.match(slacCountingSentence(data.getSite(7), 0.1), /8\.5 are non-synonymous \(85%\); 70% .* p\(dN&gt;dS\) = <strong>0\.2665<\/strong>/);
  assert.equal(slacCountingSentence(data.getSite(2), 0.1), 'No substitutions inferred at this site; SLAC cannot test it.');
});

test('diagnostics list the run as each producer wrote it', { skip: noFixtures }, () => {
  const [v3, hyphy25] = FILES.map(f => new SLACData(fixture(f)));
  const html = slacDiagnostics(v3);
  assert.match(html, /as written by <span class="case-literal">hyphy3<\/span>/);
  assert.match(html, /p ≤ 0\.10: 0 diversifying · 11 purifying · 176 neutral/);
  assert.match(html, /All · 16 tested · 0 background · tested tree length 1\.7233</);
  assert.match(html, /doi:10\.1093\/molbev\/msi105/);
  assert.match(html, /<div class="diag-item diag-item-wide">\s*<span class="diag-lbl">Citation<\/span>\s*<span class="diag-val diag-val-prose">/);
  assert.doesNotMatch(html, /branches: All/);
  const legacy = slacDiagnostics(hyphy25);
  assert.doesNotMatch(legacy, /Run cutoff/);
  assert.doesNotMatch(legacy, /Ancestor samples/);
  assert.match(slacFitsTable(hyphy25.fits, hyphy25.formatLabel), /<td class="font-mono">31<\/td>/);
  assert.match(slacFitsTable(new SLACData(fixture('slac_cd2.legacy.json')).fits, 'x'), /<td class="font-mono">—<\/td>/);
});

test('the Newick text drops v3\'s zero lengths and END root, and keeps real lengths', { skip: noFixtures }, () => {
  const v3 = new SLACData(fixture('slac_cd2.json'));
  const text = SLAC_CONFIG.newick(v3);
  assert.doesNotMatch(text, /:0|END/);
  assert.match(text, /^\(\(\(\(Pig,Cow\)Node3,Horse,Cat\)Node2,/);
  assert.match(text, /Rat,Mouse\);$/);
  assert.equal(slacNewickText('(A:0.12,B:0)C:0.0;'), '(A:0.12,B)C;');
});

test('leaf labels carry the codon and residue at the site', () => {
  assert.equal(slacLeafLabel('Human', { codon: 'GAG', aa: 'E' }), 'Human  GAG · E');
  assert.equal(slacLeafLabel('Node1', undefined), 'Node1');
});
