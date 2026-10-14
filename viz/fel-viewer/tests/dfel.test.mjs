import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import {
  DFELData, parseDfelSubstitutions, dfelChangeTally, dfelSiteClass, dfelOmega, dfelAtFloor
} from '../src/data/dfel-data.js';
import { configForResult } from '../src/methods/registry.js';
import {
  DFEL_CONFIG, DFEL_KPI_CARDS, DFEL_FILTERS, DFEL_COLUMNS, DFEL_STRIP, DFEL_P_THRESHOLD,
  dfelRate, dfelOmegaText, dfelModelVerdict, dfelDecidingTest, dfelRootPhrase, dfelRootCaption, dfelRootNotes,
  dfelRootToken, dfelRootWarning, dfelRootRanked, DFEL_EVIDENCE,
  dfelEmptyState, dfelChangeSummary, dfelChangesCell, dfelGcAtCell, dfelTallySentence, dfelDriftNote,
  dfelBranchStyle, dfelWeightsText, dfelDiagnostics, dfelFitsTable, dfelDiagnosticsPanel
} from '../src/methods/dfel.js';
import { sliderAttributes } from '../src/components/kpi-bar.js';
import { PALETTE } from '../src/data/palette.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';
const fixture = name => JSON.parse(readFileSync(new URL(name, FIXTURES), 'utf8'));

const DIVERSIFYING_AT_01 = [9, 34, 53, 55, 64, 81, 98, 102, 106, 113, 117, 149, 172];

// A hyphy3 --legacy-json dFEL document as to_legacy_json writes it: FEL-shaped, one row
// per site, nothing of the directional model
function legacyDoc(rows) {
  return {
    analysis: {
      authors: 'Sergei L Kosakovsky Pond and HyPhy Team',
      citation: 'Directional Fixed Effects Likelihood under Non-Stationary Composition Dynamics (2026)',
      info: 'Directional FEL detects site-level selection while accounting for lineage-specific compositional heterogeneity',
      version: '3.0.0'
    },
    MLE: {
      headers: [['alpha', ''], ['beta', ''], ['alpha=beta', ''], ['LRT', ''], ['p-value', ''], ['Total branch length', '']],
      content: { 0: rows }
    },
    input: { 'file name': '/data/cd2.fna', 'genetic code': 'Universal', 'number of sequences': 4, 'number of sites': rows.length, trees: { 0: '((a:1,b:1):1,c:1,d:1);' } }
  };
}

test('the substitution string parses every form the writer emits', () => {
  assert.deepEqual(parseDfelSubstitutions('Conserved'), { conserved: true, synOnly: false, events: [], synCount: 0, unparsed: [] });
  const synOnly = parseDfelSubstitutions('(+2 syn)');
  assert.deepEqual([synOnly.synOnly, synOnly.synCount, synOnly.events.length], [true, 2, 0]);
  const mixed = parseDfelSubstitutions('T → R, T → A [GC↑], K → E [AT↑]⚡, T → R (+1 syn)');
  assert.deepEqual(mixed.events, [
    { from: 'T', to: 'R', dir: null, surprising: false },
    { from: 'T', to: 'A', dir: 'GC', surprising: false },
    { from: 'K', to: 'E', dir: 'AT', surprising: true },
    { from: 'T', to: 'R', dir: null, surprising: false }
  ]);
  assert.equal(mixed.synCount, 1);
  assert.equal(mixed.synOnly, false);
  assert.deepEqual(dfelChangeTally(mixed), { nonsyn: 4, gc: 1, at: 1, flat: 2, surprising: 1, syn: 1 });
  assert.deepEqual(parseDfelSubstitutions('L → *⚡').events, [{ from: 'L', to: '*', dir: null, surprising: true }]);
  // Unreachable in the writer; read defensively as synonymous-only
  const defensive = parseDfelSubstitutions('Conserved (+syn)');
  assert.deepEqual([defensive.conserved, defensive.synOnly, defensive.events.length], [false, true, 0]);
  assert.deepEqual(parseDfelSubstitutions('garbage').unparsed, ['garbage']);
});

test('the cd2 result parses into 187 sites with live calls that match the stored column', { skip: noFixtures }, () => {
  const data = new DFELData(fixture('dfel_cd2.json'));
  assert.equal(configForResult(data.raw), DFEL_CONFIG);
  assert.deepEqual([data.format, data.formatLabel], ['v3', 'hyphy3 v3']);
  assert.deepEqual([data.numSites, data.numTaxa, data.patterns], [187, 10, 186]);
  assert.deepEqual(data.siteResults.map(r => r.site), Array.from({ length: 187 }, (_, i) => i + 1));

  const stats = data.getStats(0.1);
  assert.deepEqual([stats.diversifying, stats.purifying, stats.neutral, stats.invariant], [13, 27, 147, 13]);
  assert.deepEqual(data.siteResults.filter(r => dfelSiteClass(r, 0.1) === 'diversifying').map(r => r.site), DIVERSIFYING_AT_01);
  assert.deepEqual(data.runCheck(), { threshold: 0.1, written: { diversifying: 13, purifying: 27, neutral: 147 }, mismatches: 0 });
  assert.equal(data.runThreshold, 0.1);

  // Every token parses: 182 GC-ward, 152 AT-ward, 124 untagged, 24 surprising
  const sum = (key) => data.siteResults.reduce((n, r) => n + r.tally[key], 0);
  assert.deepEqual([sum('gc'), sum('at'), sum('flat'), sum('surprising')], [182, 152, 124, 24]);
  assert.ok(data.siteResults.every(r => r.subs.unparsed.length === 0));
  assert.equal(data.siteResults.filter(r => r.subs.conserved).length, 13);
  assert.equal(data.siteResults.filter(r => r.subs.synOnly).length, 25);

  // 44 sites have α exactly at its floor with β above it; 9 more sit there up to float noise
  const exact = data.siteResults.filter(r => r.alpha <= 1e-4 && r.beta > 1e-4).length;
  assert.equal(exact, 44);
  assert.equal(data.siteResults.filter(r => r.omega === Infinity).length, 53);
  assert.equal(data.getSite(9).omega, Infinity);
  assert.equal(dfelRate(data.getSite(144).alpha), '≥1000 (bound)');
  assert.equal(data.getSite(2).omega, null);
  assert.equal(data.getSite(87).omega, 0);

  assert.equal(data.K, 1);
  assert.deepEqual(data.regime(1).name, 'Neutral');
  assert.equal(data.regime(0), null);
  assert.equal(Object.keys(data.branchRegimes).length, 17);
  assert.equal(data.composition.heterogeneous, false);
  assert.equal(data.composition.taxa.length, 10);
  assert.equal(data.root.candidates.length, 16);
  assert.equal(data.root.candidates[0].branch, data.root.branch);
  assert.ok(Math.abs(data.root.candidates.reduce((s, c) => s + c.weight, 0) - 1) < 1e-9);
  assert.deepEqual(data.capabilities, { substitutions: true, directionalModel: true, tree: true });
  assert.equal(data.gcAtExtent, 6);
});

test('the KPI cards read the alignment, both calls and the model, never a global ω', { skip: noFixtures }, () => {
  const data = new DFELData(fixture('dfel_cd2.json'));
  const stats = data.getStats(0.1);
  const view = { data, threshold: 0.1 };
  assert.deepEqual(DFEL_KPI_CARDS.map(c => c.id), ['alignment', 'div', 'pur', 'model']);
  assert.deepEqual(DFEL_KPI_CARDS.map(c => c.value(stats, view)), ['187 Codons', '▲ 13', '▼ 27', 'Homogeneous screen\nK = 1']);
  assert.deepEqual(DFEL_KPI_CARDS.map(c => c.sub(stats, view)), [
    '10 Taxa · 186 patterns',
    '7.0% of sites (β > α)',
    '14.4% of sites (α > β)',
    'screen p = 0.061 (CP3) · root w = 8.1% (weak)'
  ]);
  assert.equal(DFEL_KPI_CARDS[3].available(data), true);
});

test('the screen verdict and the fitted K are separate facts; a K the screen would not pick is forced', () => {
  const comp = (heterogeneous) => ({ heterogeneous, overall: { p: 0.2 }, cp3: { p: 0.01 } });
  assert.deepEqual(dfelModelVerdict({ composition: comp(false), K: 1 }).forced, false);
  const forced = dfelModelVerdict({ composition: comp(false), K: 3 });
  // The KPI headline carries the qualifier: one line for the screen, one for the fit
  assert.deepEqual([forced.screen, forced.fit, forced.summary], ['Homogeneous screen', 'K = 3 (forced)', 'Homogeneous screen · K = 3 (forced)']);
  const view = { data: { composition: comp(false), K: 3, root: null } };
  assert.equal(DFEL_KPI_CARDS[3].value({}, view), 'Homogeneous screen\nK = 3 (forced)');
  assert.equal(DFEL_KPI_CARDS[3].sub({}, view), 'screen p = 0.010 (CP3)');
  assert.equal(dfelModelVerdict({ composition: comp(true), K: 1 }).summary, 'Heterogeneous screen · K = 1 (forced)');
  assert.equal(dfelModelVerdict({ composition: comp(true), K: 3 }).forced, false);
  assert.equal(dfelModelVerdict({ composition: null }), null);
  assert.equal(dfelDecidingTest(comp(true)), 'screen p = 0.010 (CP3)');
  assert.equal(dfelDecidingTest({ overall: { p: 0.0001 }, cp3: { p: 0.2 } }), 'screen p = 1.0e-4 (Overall)');
});

test('the root reads by strategy: ML edge, kept input, outgroup, or an outgroup that was not found', () => {
  const candidates = [{ branch: 'Node1 - Node2', weight: 0.6 }, { branch: 'Cow', weight: 0.4 }];
  const auto = { strategy: 'auto', branch: 'Node1 - Node2', weight: 0.081, candidates };
  assert.equal(dfelRootPhrase(auto), 'root Node1 – Node2, w = 8.1% (weak)');
  assert.equal(dfelRootPhrase(Object.assign({}, auto, { weight: 0.6 })), 'root Node1 – Node2, w = 60.0%');
  assert.equal(dfelRootCaption(auto), '◆ ML root on Node1 – Node2 (w = 8.1%)');
  assert.equal(dfelRootToken(auto), 'root w = 8.1% (weak)');
  assert.equal(dfelRootWarning(auto), null);
  assert.equal(dfelRootRanked(auto), true);
  // keep: the written weight 1.0 and log L 0 are placeholders and never shown
  const keep = { strategy: 'keep', branch: 'User Input Root (Preserved)', weight: 1, logL: 0, candidates: [] };
  assert.equal(dfelRootPhrase(keep), 'root: input (kept)');
  assert.equal(dfelRootToken(keep), 'root kept');
  assert.equal(dfelRootCaption(keep), '◆ input root, kept from the input tree (root_strategy = keep)');
  assert.equal(dfelRootRanked(keep), false);
  assert.doesNotMatch(dfelFitsTable({ fits: {}, root: keep }), /chosen root/);
  // Weak support is said for every strategy that weighs the root, not only auto
  const cow = { strategy: 'Cow', branch: 'Cow', weight: 0.4, candidates };
  assert.equal(dfelRootPhrase(cow), 'root: outgroup Cow, w = 40.0% (weak)');
  assert.equal(dfelRootToken(cow), 'outgroup root w = 40.0% (weak)');
  assert.equal(dfelRootNotes(cow)[0], 'Outgroup requested: Cow.');
  assert.equal(dfelRootWarning(cow), null);
  const zebra = { strategy: 'Zebra', branch: 'Node1 - Node2', weight: 0.081, candidates };
  assert.equal(dfelRootPhrase(zebra), 'ML root Node1 – Node2, w = 8.1% (weak)');
  assert.equal(dfelRootToken(zebra), '⚠ outgroup not found');
  assert.equal(dfelRootWarning(zebra), '⚠ Outgroup Zebra not found in the tree; the ML root is used.');
  assert.equal(dfelRootNotes(zebra)[0], 'Outgroup requested: Zebra (not found).');
  assert.equal(dfelRootNotes(auto)[1], '95% of the weight needs 2 of 2 placements.');
  assert.equal(dfelRootPhrase(Object.assign({}, cow, { weight: 0.7 })), 'root: outgroup Cow, w = 70.0%');
});

test('the hero tooltip prints α with its bounds, as the table and drawer do', () => {
  assert.equal(DFEL_EVIDENCE.alpha({ alpha: 1000 }), '≥1000 (bound)');
  assert.equal(DFEL_EVIDENCE.alpha({ alpha: 1e-4 }), '0');
});

test('rates and ω spell out the optimiser bounds', () => {
  assert.equal(dfelAtFloor(1.0000000000000186e-4), true);
  assert.equal(dfelAtFloor(1.01e-4), false);
  assert.deepEqual([dfelRate(0), dfelRate(1e-4), dfelRate(2.3159), dfelRate(1000)], ['0', '0', '2.32', '≥1000 (bound)']);
  assert.equal(dfelOmega(1e-4, 1.5), Infinity);
  assert.equal(dfelOmega(0, 0), null);
  assert.equal(dfelOmega(2, 1e-4), 0);
  assert.equal(dfelOmega(2, 1), 0.5);
  assert.deepEqual([Infinity, null, 0, 0.5, 2500].map(dfelOmegaText), ['∞', '—', '0', '0.50', '>1000']);
});

test('calls are recomputed from p, α and β; a conserved site is the invariant state of neutral', () => {
  const row = (alpha, beta, pValue, invariant = false) => ({ alpha, beta, pValue, invariant });
  assert.equal(dfelSiteClass(row(0.1, 2, 0.05), 0.1), 'diversifying');
  assert.equal(dfelSiteClass(row(0.1, 2, 0.05), 0.01), 'neutral');
  assert.equal(dfelSiteClass(row(2, 0.1, 0.05), 0.1), 'purifying');
  assert.equal(dfelSiteClass(row(1, 1, 0.05), 0.1), 'neutral');
  assert.equal(dfelSiteClass(row(0, 0, 1, true), 0.1), 'invariant');
});

test('a --legacy-json export maps onto the same site rows and nothing of the model', () => {
  const data = new DFELData(legacyDoc([
    [1e-4, 8.45, 2.94, 1.22, 0.04, 12.8],
    [0, 0, 0, 0, 1, 0],
    [3.1, 0.2, 1.1, 5.0, 0.02, 3.3]
  ]));
  assert.equal(configForResult(data.raw), DFEL_CONFIG);
  assert.deepEqual([data.format, data.numSites, data.numTaxa], ['hyphy3-legacy', 3, 4]);
  assert.deepEqual(data.siteResults.map(r => [r.site, r.alphaNull, r.omega, r.invariant]),
    [[1, 2.94, Infinity, false], [2, 0, null, true], [3, 1.1, 0.2 / 3.1, false]]);
  assert.deepEqual(data.getStats(0.05), { totalSites: 3, diversifying: 1, purifying: 1, significant: 2, neutral: 1, invariant: 1 });
  assert.deepEqual(data.capabilities, { substitutions: false, directionalModel: false, tree: true });
  assert.equal(data.K, null);
  assert.equal(data.runCheck(), null);
  assert.equal(data.dataset.file, '/data/cd2.fna');

  // Columns, the strip and the model KPI that need the directional model are dropped; the
  // format notes say why
  const shown = DFEL_COLUMNS.filter(c => !c.availableWhen || c.availableWhen(data.capabilities)).map(c => c.id);
  assert.deepEqual(shown, ['codon', 'alpha', 'beta', 'omega', 'lrt', 'p', 'detail', 'actions']);
  assert.equal(DFEL_STRIP.available(data), false);
  assert.deepEqual(DFEL_KPI_CARDS.filter(c => !c.available || c.available(data)).map(c => c.id), ['alignment', 'div', 'pur']);
  assert.match(dfelDiagnostics(data), /hyphy3 --legacy-json/);
});

test('the GC/AT strip starts shown only when directional regimes were fitted', () => {
  const caps = { substitutions: true };
  assert.equal(DFEL_STRIP.shown({ K: 1, capabilities: caps }), false);
  assert.equal(DFEL_STRIP.shown({ K: 3, capabilities: caps }), true);
  assert.match(DFEL_STRIP.note({ K: 1 }), /^Direction of reconstructed changes \(single Neutral regime fitted\)/);
  assert.doesNotMatch(DFEL_STRIP.note({ K: 3 }), /single Neutral/);
});

test('table cells: GC/AT counts, changes cut to whole tokens with the synonymous count apart', () => {
  const conserved = { subs: parseDfelSubstitutions('Conserved'), substitutions: 'Conserved' };
  assert.match(dfelChangesCell(conserved), /Invariant \(no substitutions\)/);
  const text = 'N → A [GC↑], N → D [GC↑], N → S [GC↑], N → T [GC↑], N → S [GC↑], N → K [AT↑] (+1 syn)';
  const cell = dfelChangesCell({ subs: parseDfelSubstitutions(text), substitutions: text });
  assert.match(cell, />N → A \[GC↑\], N → D \[GC↑\], N → S \[GC↑\], N → T \[GC↑\], … <span class="dfel-muted">\+1 syn<\/span>/);
  assert.match(cell, new RegExp(`title="${text.replace(/[()+[\]]/g, '\\$&')}"`));
  const synOnly = dfelChangesCell({ subs: parseDfelSubstitutions('(+2 syn)'), substitutions: '(+2 syn)' });
  assert.match(synOnly, /<div class="sub-text-wrapper"><span class="dfel-muted">\+2 syn<\/span><\/div>/);

  const tally = { gc: 3, at: 1 };
  assert.match(dfelGcAtCell({ tally }, 6), /<span>3 \/ 1<\/span>.*width: 50%.*width: 16\.6/);
  assert.equal(dfelGcAtCell({ tally: null }, 6), '<td></td>');
  // No GC- or AT-ward change: a muted dash, no empty track
  assert.equal(dfelGcAtCell({ tally: { gc: 0, at: 0 } }, 6), '<td class="font-mono"><span class="dfel-muted">—</span></td>');
});

test('the default filter is diversifying by p; an empty one jumps to purifying or shows closest calls', () => {
  assert.deepEqual(DFEL_FILTERS.map(f => f.id), ['diversifying', 'purifying', 'significant', 'all']);
  assert.deepEqual(DFEL_FILTERS[0].sort, { column: 'p', dir: 'asc' });
  const div = DFEL_FILTERS[0];
  const jump = dfelEmptyState(0.01, { filter: div, stats: { purifying: 4 }, rows: [] });
  assert.equal(jump.title, 'No diversifying site reaches p ≤ 0.01');
  assert.deepEqual(jump.jump, { filter: 'purifying', label: 'View Purifying Sites (▼ 4)' });
  const rows = [{ site: 1, alpha: 1, beta: 2, pValue: 0.3 }, { site: 2, alpha: 2, beta: 1, pValue: 0.02 }, { site: 3, alpha: 0, beta: 1, pValue: 0.2 }];
  const fallback = dfelEmptyState(0.01, { filter: div, stats: { purifying: 0 }, rows });
  assert.deepEqual(fallback.fallback.map(r => r.site), [3, 1]);
});

test('the drawer tallies changes and flags drift only for one-directional diversifying sites under K > 1', () => {
  assert.equal(dfelTallySentence({ nonsyn: 4, gc: 1, at: 2, flat: 1, surprising: 0 }),
    '4 non-synonymous: 1 GC-ward, 2 AT-ward, 1 no net GC change; 0 high-surprise (⚡ = −log₂ P ≥ 10 bits under the branch\'s transition matrix)');
  const site = (gc, at, flat) => ({ alpha: 0.1, beta: 2, pValue: 0.01, tally: { nonsyn: gc + at + flat, gc, at, flat } });
  assert.match(dfelDriftNote(site(3, 1, 0), 0.1, 3), /same base-composition direction/);
  assert.equal(dfelDriftNote(site(3, 1, 0), 0.1, 1), null);
  assert.equal(dfelDriftNote(site(2, 1, 1), 0.1, 3), null);
  assert.equal(dfelDriftNote(Object.assign(site(4, 0, 0), { pValue: 0.5 }), 0.1, 3), null);
  assert.equal(dfelChangeSummary({ subs: { conserved: false }, tally: { gc: 3, at: 1, surprising: 1, syn: 2 } }), 'nonsyn: 3 GC↑ / 1 AT↑ / 1 ⚡, +2 syn');
  assert.equal(dfelChangeSummary({ tally: null }), null);
});

test('regime branch styles: GC-drift solid, AT-drift dashed, Neutral thin grey; plain ink with one regime', () => {
  assert.deepEqual(dfelBranchStyle(2, 3, PALETTE), { stroke: PALETTE.ink, width: 2.5 });
  assert.deepEqual(dfelBranchStyle(0, 3, PALETTE), { stroke: PALETTE.ink, width: 2.5, dash: [5, 3] });
  assert.deepEqual(dfelBranchStyle(1, 3, PALETTE), { stroke: PALETTE.slateLight, width: 1 });
  assert.deepEqual(dfelBranchStyle(1, 1, PALETTE), {});
  assert.equal(dfelWeightsText([0.36, 0.58, 0.06]), 'AT 36% · Neutral 58% · GC 6%');
});

test('fits are listed without a ΔAICc; the cutoff controller spans 0.001-0.20 in steps of 0.001', { skip: noFixtures }, () => {
  const data = new DFELData(fixture('dfel_cd2.json'));
  const fits = dfelFitsTable(data);
  assert.match(fits, /Nucleotide GTR<\/td><td class="font-mono">-3532\.33<\/td><td class="font-mono">7110\.53/);
  assert.match(fits, /Directional MG94 \(codon\)<\/td><td class="font-mono">-3467\.00/);
  assert.match(fits, /Directional nucleotide model at chosen root<\/td><td class="font-mono">-3668\.72<\/td><td class="font-mono">—/);
  assert.doesNotMatch(fits, /Δ/);
  const diag = dfelDiagnostics(data);
  assert.match(diag, /<span title="[^"]*cd2\.fna">cd2\.fna<\/span>/);
  assert.match(diag, /13 diversifying · 27 purifying · 147 neutral; matches the live call/);
  assert.match(diag, /1\.7233 \(baseline units\)/);
  assert.deepEqual(sliderAttributes(DFEL_P_THRESHOLD, 0.1), { min: '0.001', max: '0.20', step: '0.001', value: '0.10' });
  // A run cutoff under 0.005 sits on a stop, so the thumb shows it and it can be restored
  assert.equal(sliderAttributes(DFEL_P_THRESHOLD, 0.003).value, '0.003');
  assert.equal(DFEL_P_THRESHOLD.caption, 'Uncorrected per-site p-values');
});

test('format notes speak only of what the file has', () => {
  const notes = (directionalModel) => {
    const el = { html: '', insertAdjacentHTML(where, html) { this.html += html; } };
    dfelDiagnosticsPanel(el, { fits: {}, root: null, formatLabel: 'hyphy3 --legacy-json', capabilities: { directionalModel } });
    return el.html;
  };
  assert.match(notes(false), /not reported in hyphy3 --legacy-json output/);
  assert.doesNotMatch(notes(false), /GTR|regimes 1 or 3/);
  assert.match(notes(true), /GTR is a nucleotide model/);
});
