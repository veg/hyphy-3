import { test } from 'node:test';
import assert from 'node:assert/strict';
import { installFakeDOM, FakeElement } from './fake-dom.mjs';
import { GardResults } from '../src/data/gard-data.js';
import { SegmentBand } from '../src/components/segment-band.js';
import { SiteTrackCanvas } from '../src/components/site-track.js';
import { resolveColumns } from '../src/components/results-table.js';
import { PALETTE } from '../src/data/palette.js';
import {
  GARD_BRUSH, GARD_COLUMNS, GARD_DRAWER, GARD_FILTERS, GardTickBand, gardBreakLabel, gardBreakTone, gardDiagnostics,
  gardEmptyState, gardHero, gardLogSupport, gardSegmentBreak, gardSupportTick, gardTestSentence, gardTrackModes, gardWindow
} from '../src/methods/gard.js';
import {
  GardTreeCard, GardTreePair, gardAddedLabel, gardAddedSites, gardDefaultPair, gardLeafLabel, gardMatrixMarkup, gardStepsSVG
} from '../src/components/gard-trees.js';

installFakeDOM();

function recordingContext() {
  const calls = [];
  return new Proxy({ calls }, {
    get(target, key) {
      if (key in target) return target[key];
      if (key === 'measureText') return (text) => ({ width: String(text).length * 6 });
      return (...args) => calls.push([key, ...args]);
    },
    set(target, key, value) {
      calls.push([`=${String(key)}`, value]);
      return true;
    }
  });
}

FakeElement.prototype.getContext = () => recordingContext();
FakeElement.prototype.append = function append(...els) { this.children.push(...els); };
FakeElement.prototype.getBoundingClientRect = () => ({ left: 0, top: 0, width: 800, height: 290 });
globalThis.ResizeObserver = class { observe() {} disconnect() {} };

const TREE_A = '((A:0.1,B:0.1)Node1:0.2,(C:0.1,D:0.1)Node2:0.2,E:0.3)Node0;';
const TREE_B = '((A:0.1,C:0.1)Node1:0.2,(B:0.1,D:0.1)Node2:0.2,E:0.3)Node0;';

// Three segments, two breakpoints (after 100 and 200), the first incongruent and the second not
function threeSegmentDoc({ tests = true, ohmic = true, sites = 300 } = {}) {
  const breakTest = (breakpoint, p, left, right) => ({
    breakpoint,
    left_partition: left,
    right_partition: right,
    delta_log_l_left: 12.5,
    delta_log_l_right: 9.25,
    kh_p_value_left: p,
    kh_p_value_right: p,
    sh_p_value_left: p,
    sh_p_value_right: p,
    statistically_significant: p < 0.0125
  });
  return {
    analysis: { id: 'gard' },
    dataset: { file: '/private/dir/aln.fna', taxa_count: 5, sites_count: sites, variable_sites_count: 80, potential_breakpoints: 70, minimum_partition_size: 10 },
    baseline: { aicc: 1000, log_likelihood: -480, tree: TREE_A },
    model_selection: {
      best_model_breakpoints: [99, 199],
      best_model_aicc: 900,
      improvements: [
        { breakpoints_count: 1, breakpoints: [99], delta_aicc: 70 },
        { breakpoints_count: 2, breakpoints: [99, 199], delta_aicc: 30 }
      ]
    },
    partitions: [
      { partition_index: 0, start_site: 1, end_site: 100, log_likelihood: -100, tree: TREE_A },
      { partition_index: 1, start_site: 101, end_site: 200, log_likelihood: -110, tree: TREE_B },
      { partition_index: 2, start_site: 201, end_site: sites, log_likelihood: -120, tree: TREE_A }
    ],
    site_breakpoint_support: { 99: 0.9, 150: 1e-320, 199: 0.4, 250: 0 },
    topology_test: {
      single_tree_aicc: 990,
      delta_aicc_topology: 90,
      topological_incongruence_supported: true,
      test_type: 'both',
      bootstrap_replicates: 1000,
      breakpoints: tests ? [breakTest(100, 0.001, [1, 100], [101, 200]), breakTest(200, 0.3, [101, 200], [201, sites])] : [],
      partition_tree_matrix: tests ? [
        { partition_index: 0, interval: [1, 100], best_tree_index: 0, tree_log_likelihoods: [-100, -150, -101], kh_p_values: [1, 0.001, 0.5], sh_p_values: [1, 0.002, 0.6] },
        { partition_index: 1, interval: [101, 200], best_tree_index: 1, tree_log_likelihoods: [-160, -110, -170], kh_p_values: [0.001, 1, 0.9], sh_p_values: [0.004, 1, 0.9] }
      ] : []
    },
    ohmic_network: ohmic ? { bridges_detected: true, num_bridges: 1, delta_aicc_vs_baseline: 4, bridge_candidate_breakpoints: [250, 40], bridges: [] } : undefined,
    total_models_evaluated: 50,
    time_elapsed: 0.5
  };
}

const view = (data, threshold = 0.05) => ({ data, threshold, statistic: 'α' });

test('a break\'s tone colours its rule and label on the canvas and in the SVG export', () => {
  const data = new GardResults(threeSegmentDoc());
  const band = new SegmentBand({ segments: (d) => d.segments, breakTone: gardBreakTone, breakLabel: gardBreakLabel });
  const numSites = data.numSites;
  const frame = { left: 52, top: 248, width: 400, height: 28, scaleX: (s) => 52 + ((s - 1) / (numSites - 1)) * 400 };
  const main = Object.assign({}, frame, { top: 24, height: 200 });
  const track = { rows: [], numSites, data, threshold: 0.05, frames: [main, frame], selectedRegion: null, hoveredRegion: null, colors: Object.assign({}, PALETTE, { ink: '#F8FAFC', stone: '#94A3B8' }) };
  const ctx = recordingContext();
  band.draw(ctx, frame, track);
  const strokes = ctx.calls.filter(c => c[0] === '=strokeStyle').map(c => c[1]);
  // The first break is incongruent (accent); the second is branch-length only (the themed ink)
  assert.deepEqual(strokes.slice(-2), [PALETTE.accent, '#F8FAFC']);
  assert.deepEqual(ctx.calls.filter(c => c[0] === 'fillText').map(c => c[1]).slice(-2), ['BP 1 ✱', 'BP 2']);

  const svg = band.svg(frame, track);
  assert.match(svg, new RegExp(`stroke="${PALETTE.accent}" stroke-width="1.5"`));
  assert.match(svg, new RegExp(`stroke="${PALETTE.ink}" stroke-width="1.5"`));
  assert.match(svg, />BP 1 ✱<\/text>/);

  // Untested breaks are muted; raising the cutoff past the data turns the first to branch-only
  assert.deepEqual([0.05, 0.001].map(t => gardBreakTone({ index: 1 }, { data, threshold: t })), ['accent', 'ink']);
  const untested = new GardResults(threeSegmentDoc({ tests: false }));
  assert.equal(gardBreakTone({ index: 1 }, { data: untested, threshold: 0.05 }), 'muted');
  // The default band keeps its single accent
  const plain = new SegmentBand({ segments: (d) => d.segments });
  const plainCtx = recordingContext();
  plain.draw(plainCtx, frame, track);
  assert.deepEqual(plainCtx.calls.filter(c => c[0] === '=strokeStyle').map(c => c[1]).slice(-2), [PALETTE.accent, PALETTE.accent]);
});

test('a clicked segment opens the breakpoint on its right, the last one the break on its left', () => {
  const data = new GardResults(threeSegmentDoc());
  assert.deepEqual([0, 1, 2].map(i => gardSegmentBreak(data, i).site), [100, 200, 200]);
  assert.equal(gardSegmentBreak(new GardResults({ ...threeSegmentDoc(), partitions: [] }), 0), null);
});

test('log support is clamped at 1e-6, so underflowed values sit on the axis floor; both modes are offered', () => {
  const [log, linear] = gardTrackModes();
  assert.deepEqual([log.id, linear.id], ['log', 'linear']);
  assert.equal(gardLogSupport(1e-320), -6);
  assert.equal(gardLogSupport(0), -6);
  assert.equal(gardLogSupport(0.1), -1);
  assert.deepEqual([-6, -3, -1, 0].map(gardSupportTick), ['≤1e-6', '1e-3', '0.1', '1']);

  const data = new GardResults(threeSegmentDoc());
  const frame = { left: 52, top: 24, width: 400, height: 200, scaleX: (s) => 52 + (s / 300) * 400 };
  const track = { rows: data.support, data, threshold: 0.05, selectedSite: null, hoveredSite: null, siteRange: () => [1, 300] };
  const ctx = recordingContext();
  log.band.drawBackground(ctx, frame, track);
  const points = log.band.drawMarks(ctx, frame, track);
  assert.deepEqual([log.band.yMin, log.band.yMax], [-6, 0]);
  const at = (site) => points.find(p => p.site === site);
  assert.equal(at(151).y, frame.top + frame.height);
  assert.equal(at(251).y, frame.top + frame.height);
  assert.equal(at(100).y, frame.scaleY(Math.log10(0.9)));
  assert.equal(at(100).tone, 'div');
  assert.equal(at(200).tone, null);
  // The ticks read as supports, from the floor to 1
  assert.deepEqual(ctx.calls.filter(c => c[0] === 'fillText').map(c => c[1]), ['≤1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '0.1', '1', 'support'].slice(0, 7).concat('support'));

  const flat = recordingContext();
  linear.band.drawBackground(flat, frame, track);
  const linearPoints = linear.band.drawMarks(flat, frame, track);
  assert.deepEqual([linear.band.yMin, linear.band.yMax], [0, 1]);
  assert.equal(linearPoints.find(p => p.site === 100).y, frame.scaleY(0.9));
  assert.match(log.band.tooltip(data.support[0], track), /Breakpoint after site 100[\s\S]*BP 1 · Incongruent[\s\S]*0\.900/);
  assert.doesNotMatch(log.band.tooltip(data.support.find(s => s.site === 151), track), /BP/);
});

test('the overview brush starts a tenth of a long alignment wide, at least 300 sites', () => {
  assert.equal(gardWindow(561), 300);
  assert.equal(gardWindow(3000), 300);
  assert.equal(gardWindow(15000), 1500);
  assert.equal(GARD_BRUSH.minSites, 1000);
  assert.equal(GARD_BRUSH.tone({ p: 0.5 }), 'accent');
  assert.equal(GARD_BRUSH.tone({ p: 0.49 }), null);
  assert.equal(GARD_BRUSH.value({ p: 1 }), 6);
  assert.equal(GARD_BRUSH.value({ p: 1e-300 }), 0);

  const hero = gardHero(new FakeElement(), { threshold: 0.05, onSiteSelect() {} });
  hero.setData(new GardResults(threeSegmentDoc({ sites: 15000 })));
  assert.equal(hero.zoomed, true);
  assert.equal(hero.span, 1500);
  hero.setData(new GardResults(threeSegmentDoc()));
  assert.equal(hero.zoomed, false);
});

test('the Ohmic candidates are hit regions under the segments, and the band takes no room without them', () => {
  const data = new GardResults(threeSegmentDoc());
  const band = new GardTickBand();
  const track = { data, siteRange: () => [1, 300], colors: PALETTE };
  band.measure(track);
  assert.equal(band.height, 20);
  const frame = { left: 52, top: 260, width: 400, height: 20, scaleX: (s) => 52 + ((s - 1) / 299) * 400 };
  const ctx = recordingContext();
  const regions = band.draw(ctx, frame, track);
  assert.deepEqual(regions.map(r => [r.kind, r.key]), [['candidate', 41], ['candidate', 251]]);
  assert.equal(regions[0].x0, frame.scaleX(41) - 5);
  assert.deepEqual([regions[0].y0, regions[0].y1], [260, 280]);
  assert.equal(ctx.calls.filter(c => c[0] === 'closePath').length, 2);
  assert.match(band.regionTooltip(regions[1]), /After site 251[\s\S]*Ohmic bridge candidate[\s\S]*seeded GA/);
  assert.equal((band.svg(frame, track).match(/<polygon/g) || []).length, 2);
  // A zoomed window keeps only the candidates in view
  assert.deepEqual(band.draw(recordingContext(), frame, Object.assign({}, track, { siteRange: () => [100, 300] })).map(r => r.key), [251]);

  const bare = { data: new GardResults(threeSegmentDoc({ ohmic: false })), siteRange: () => [1, 300] };
  band.measure(bare);
  assert.equal(band.height, 0);
  assert.equal(band.draw(recordingContext(), frame, bare), null);
  assert.equal(band.svg(frame, bare), '');
});

test('the stack lays out the Ohmic strip only when the result has one', () => {
  const track = new SiteTrackCanvas(new FakeElement(), {
    modes: gardTrackModes(),
    stack: [new SegmentBand({ segments: (d) => d.segments }), new GardTickBand()],
    sites: (d) => d.support
  });
  track.setData(new GardResults(threeSegmentDoc()));
  assert.deepEqual(track.layout().slice(1).map(f => f.height), [28, 20]);
  track.setData(new GardResults(threeSegmentDoc({ ohmic: false })));
  assert.deepEqual(track.layout().slice(1).map(f => f.height), [28, 0]);
});

test('the table drops the topology columns for a result without tests', () => {
  const ids = (data) => resolveColumns(GARD_COLUMNS, data).map(c => c.id);
  assert.deepEqual(ids(new GardResults(threeSegmentDoc())), ['bp', 'site', 'left', 'right', 'call', 'dlogl', 'kh', 'sh', 'rf', 'moved', 'actions']);
  assert.deepEqual(ids(new GardResults(threeSegmentDoc({ tests: false }))), ['bp', 'site', 'left', 'right', 'call', 'rf', 'moved', 'actions']);
});

test('table cells show the call at the live alpha, the flank pairs and the tree distance', () => {
  const data = new GardResults(threeSegmentDoc());
  const [first, second] = data.breaks;
  const cell = (id, r, threshold = 0.05) => GARD_COLUMNS.find(c => c.id === id).cell(r, view(data, threshold));
  assert.match(cell('call', first), /gard-call-incongruent[^>]*>✱ Incongruent/);
  assert.match(cell('call', second), /gard-call-branch[^>]*>≈ Branch lengths only/);
  assert.match(cell('call', first, 0.001), /gard-call-branch/);
  // B = 2: the cutoff is 0.0125 at alpha 0.05, so p = 0.001 passes and 0.3 does not
  assert.match(cell('kh', first), /font-bold text-orange">0\.0010 \/ 0\.0010/);
  assert.doesNotMatch(cell('kh', second), /font-bold/);
  assert.equal(cell('dlogl', first), '<td class="font-mono">12.5 / 9.3</td>');
  assert.equal(cell('left', second), '<td class="font-mono">101–200</td>');
  assert.match(cell('rf', first), /^<td class="font-mono">4 \/ 4<\/td>$/);
  assert.equal(cell('moved', first), '<td class="font-mono">4</td>');
  assert.match(cell('actions', first), /data-action="trees" data-site="100"/);

  const csv = GARD_COLUMNS.flatMap(c => c.csv || []);
  assert.deepEqual(csv.map(f => f.header).slice(0, 4), ['Breakpoint', 'After site', 'Left start', 'Left end']);
  assert.equal(csv.find(f => f.header === 'Call').value(first, view(data)), 'Incongruent');
  assert.equal(csv.find(f => f.header === 'KH p left').value(first, view(data)), 0.001);

  const bare = new GardResults(threeSegmentDoc({ tests: false }));
  const untested = GARD_COLUMNS.find(c => c.id === 'call').cell(bare.breaks[0], view(bare));
  assert.match(untested, /gard-call-untested[^>]*>– Untested/);
});

test('the filters count and select breakpoints by call; an empty table says why', () => {
  const data = new GardResults(threeSegmentDoc());
  const [incongruent, all] = GARD_FILTERS;
  assert.deepEqual([incongruent, all].map(f => f.count(data.getStats(0.05))), [1, 2]);
  assert.deepEqual(data.breaks.filter(b => incongruent.test(b, 0.05, view(data))).map(b => b.index), [1]);
  assert.deepEqual(data.breaks.filter(b => incongruent.test(b, 0.2, view(data, 0.2))).map(b => b.index), [1, 2].slice(0, 1));
  assert.equal(gardEmptyState(0.05, { filter: all, stats: { breakpoints: 0 } }).title, 'No breakpoints — a single tree fits best (ΔAICc 0)');
  const empty = gardEmptyState(0.01, { filter: incongruent, stats: { breakpoints: 2 } });
  assert.match(empty.title, /No incongruent breakpoint at α ≤ 0\.01/);
  assert.equal(empty.jump.filter, 'all');
});

test('the drawer headers a breakpoint by call and steps through its two tabs', () => {
  const data = new GardResults(threeSegmentDoc());
  assert.deepEqual(GARD_DRAWER.header(data.breaks[0], view(data)), { tag: 'Incongruent · after site 100', tone: 'div', badge: 'BP 1' });
  assert.deepEqual(GARD_DRAWER.header(data.breaks[1], view(data)), { tag: 'Branch lengths only · after site 200', tone: 'neutral', badge: 'BP 2' });
  assert.equal(GARD_DRAWER.row(data, 200), data.breaks[1]);
  assert.deepEqual(GARD_DRAWER.tabs.map(t => t.id), ['trees', 'test']);

  const [trees, test] = GARD_DRAWER.tabs;
  assert.equal(test.availableWhen(data), true);
  assert.equal(test.availableWhen(new GardResults(threeSegmentDoc({ tests: false }))), false);
  assert.equal(trees.availableWhen, undefined);

  const panel = new FakeElement();
  test.render(panel, data.breaks[0], { data, threshold: 0.05 });
  assert.match(panel.innerHTML, /Topology test at α ≤ 0\.05/);
  assert.match(panel.innerHTML, /KH p, left[\s\S]*text-orange">0\.0010/);
  assert.match(panel.innerHTML, /left flank 1–100/);
  assert.match(panel.innerHTML, /α \/ 2B = 0\.0125 for B = 2 breakpoints/);
  assert.match(gardTestSentence(data, data.breaks[0], 0.05), /recombination/);
  assert.match(gardTestSentence(data, data.breaks[1], 0.05), /rate heterogeneity/);
  test.render(panel, data.breaks[1], { data, threshold: 0.05 });
  assert.doesNotMatch(panel.innerHTML, /text-orange/);

  // The trees tab draws the two segments either side
  const host = new FakeElement();
  const mounted = trees.mount(host);
  trees.render(host, data.breaks[1], { data, view: mounted });
  assert.deepEqual(mounted.pair.indexes, [1, 2]);
  assert.equal(mounted.pair.sides[0].title.textContent, 'S2 · 101–200');
  assert.equal(mounted.pair.sides[1].title.textContent, 'S3 · 201–300');
});

test('the tree pair marks moved taxa by an accent stroke and a diamond in the label', () => {
  const data = new GardResults(threeSegmentDoc());
  const pair = new GardTreePair(new FakeElement());
  pair.setSegments(data, 0, 1);
  assert.deepEqual([...pair.moved].sort(), ['A', 'B', 'C', 'D']);
  const [left, right] = pair.trees;
  const leaf = left.layout.leaves.find(n => n.name === 'A');
  assert.deepEqual(left.options.style(leaf, undefined, PALETTE), { stroke: PALETTE.accent, width: 2 });
  assert.equal(left.options.style(left.layout.leaves.find(n => n.name === 'E'), undefined, PALETTE), null);
  assert.equal(left.options.style(left.layout.nodes.find(n => n.name === 'Node1'), undefined, PALETTE), null);
  assert.equal(left.leafLabel(leaf), '◆ A');
  assert.equal(right.leafLabel(right.layout.leaves.find(n => n.name === 'E')), 'E');
  assert.equal(gardLeafLabel('A', new Set(['A'])), '◆ A');
  assert.match(pair.summary.innerHTML, /Robinson–Foulds distance <strong>4<\/strong> of 4 · 4 moved: A, B, C, D/);
  assert.equal(left.legend.style.display, '');
  assert.equal(pair.sides[1].col.style.display, '');

  // The same tree twice: nothing moved, no legend
  pair.setSegments(data, 0, 2);
  assert.equal(pair.moved.size, 0);
  assert.match(pair.summary.innerHTML, /no taxon moved/);
  assert.equal(left.legend.style.display, 'none');

  // Hovering a leaf lights its taxon in the other tree, an internal node nothing
  const highlighted = [];
  right.highlight = (name) => highlighted.push(name);
  left.options.onHover('A', { children: [] });
  left.options.onHover('Node1', { children: [{}] });
  assert.deepEqual(highlighted, ['A', null]);

  // A single segment draws one tree
  pair.setSegments(new GardResults(Object.assign(threeSegmentDoc(), { partitions: [] })), 0, null);
  assert.equal(pair.sides[1].col.style.display, 'none');
  assert.match(pair.summary.innerHTML, /no breakpoints/);
});

test('trees past 60 taxa are laid out as cladograms', () => {
  const data = new GardResults(threeSegmentDoc());
  data.numTaxa = 266;
  const pair = new GardTreePair(new FakeElement());
  pair.setSegments(data, 0, 1);
  assert.equal(pair.trees[0].layout.cladogram, true);
  assert.match(pair.summary.innerHTML, /cladogram above 60 taxa/);
  data.numTaxa = 5;
  pair.setSegments(data, 0, 1);
  assert.equal(pair.trees[0].layout.cladogram, false);
});

test('the tree card starts on the pair around the selected breakpoint, else the first two segments', () => {
  const data = new GardResults(threeSegmentDoc());
  assert.deepEqual(gardDefaultPair(data, null), [0, 1]);
  assert.deepEqual(gardDefaultPair(data, 200), [1, 2]);
  assert.deepEqual(gardDefaultPair(new GardResults(Object.assign(threeSegmentDoc(), { partitions: [] })), null), [0, null]);

  const card = new GardTreeCard(new FakeElement());
  card.setData(data);
  assert.deepEqual(card.indexes, [0, 1]);
  assert.match(card.selects[0].innerHTML, /<option value="0">S1 · 1–100<\/option><option value="1">S2 · 101–200<\/option><option value="2">S3 · 201–300<\/option>/);
  assert.equal(card.selects[0].value, '0');
  assert.equal(card.selects[1].value, '1');
  card.selectSite(200);
  assert.deepEqual(card.pair.indexes, [1, 2]);
  assert.equal(card.selects[0].value, '1');
  assert.equal(card.selects[1].value, '2');
  card.selectSite(999);
  assert.deepEqual(card.pair.indexes, [1, 2]);
  card.selectSite(null);
  assert.deepEqual(card.pair.indexes, [1, 2]);
  card.selects[1].value = '0';
  card.selects[1].fire('change');
  assert.deepEqual(card.pair.indexes, [1, 0]);

  const single = new GardTreeCard(new FakeElement());
  single.setData(new GardResults(Object.assign(threeSegmentDoc(), { partitions: [] })));
  assert.equal(single.el.querySelector('#gard-selects').style.display, 'none');
  assert.equal(single.pair.sides[1].col.style.display, 'none');
});

test('the model-selection chart has a bar per step with the sites each step added', () => {
  const data = new GardResults(threeSegmentDoc());
  assert.deepEqual(data.steps, [{ k: 1, deltaAicc: 70, bps: [100] }, { k: 2, deltaAicc: 30, bps: [100, 200] }]);
  assert.deepEqual([0, 1].map(i => gardAddedSites(data.steps, i)), [[100], [200]]);
  assert.deepEqual(gardAddedSites([{ k: 1, bps: null }], 0), []);
  assert.deepEqual([[], [5], [1, 2, 3], [1, 2, 3, 4]].map(gardAddedLabel), ['', '+5', '+1, 2, 3', '+1, 2, 3 …']);
  const { svg, width } = gardStepsSVG(data.steps, { best: 2, width: 400, colors: PALETTE });
  assert.equal(width, 400);
  assert.equal((svg.match(/<rect/g) || []).length, 3);
  assert.match(svg, />1 BP<\/text>/);
  assert.match(svg, />\+200<\/text>/);
  assert.match(svg, new RegExp(`font-weight="bold" style="fill: ${PALETTE.ink}">2 BP</text>`));
  assert.match(svg, />70\.0<\/text>/);

  const card = new GardTreeCard(new FakeElement());
  card.setData(new GardResults(threeSegmentDoc()));
  assert.match(card.el.querySelector('#gard-steps').innerHTML, /<svg[^>]*class="gard-steps"/);
  card.setData(new GardResults(Object.assign(threeSegmentDoc(), { model_selection: { best_model_breakpoints: [], improvements: [] } })));
  assert.equal(card.el.querySelector('#gard-steps').innerHTML, '');
  assert.match(card.el.querySelector('#gard-steps-caption').textContent, /No breakpoint improved/);
});

test('the partition-by-tree matrix bolds the best tree and shades rejected ones; it hides without data', () => {
  const data = new GardResults(threeSegmentDoc());
  const cutoff = 0.05 / (2 * data.breaks.length);
  const html = gardMatrixMarkup(data.matrix, { mode: 'logL', cutoff });
  assert.match(html, /Tree 1<\/th>.*Tree 3<\/th>/);
  assert.match(html, /<th>P1<\/th><td class="font-mono">1–100<\/td><td class="font-mono gard-matrix-best">-100\.0<\/td><td class="font-mono gard-matrix-rejected">-150\.0<\/td><td class="font-mono">-101\.0<\/td>/);
  assert.match(html, /<th>P2<\/th>[^]*gard-matrix-rejected">-160\.0[^]*gard-matrix-best">-110\.0/);
  assert.match(gardMatrixMarkup(data.matrix, { mode: 'kh', cutoff }), /gard-matrix-best">1\.0000/);
  assert.match(gardMatrixMarkup(data.matrix, { mode: 'sh', cutoff }), /gard-matrix-rejected">0\.0020/);
  // No tree is shaded when no test reaches the cutoff
  assert.doesNotMatch(gardMatrixMarkup(data.matrix, { cutoff: 1e-9 }), /gard-matrix-rejected/);

  const card = new GardTreeCard(new FakeElement());
  card.setThreshold(0.05);
  card.setData(data);
  assert.equal(card.el.querySelector('#gard-matrix-panel').style.display, '');
  assert.match(card.el.querySelector('#gard-matrix').innerHTML, /gard-matrix-rejected/);
  assert.match(card.el.querySelector('#gard-matrix-caption').textContent, /α \/ 2B = 0\.0125/);
  card.setThreshold(0.001);
  assert.match(card.el.querySelector('#gard-matrix-caption').textContent, /α \/ 2B = 0\.00025/);
  card.setData(new GardResults(threeSegmentDoc({ tests: false })));
  assert.equal(card.el.querySelector('#gard-matrix-panel').style.display, 'none');
});

test('diagnostics name the run and show the input file by its name, never its path', () => {
  const html = gardDiagnostics(new GardResults(threeSegmentDoc()));
  assert.match(html, />aln\.fna</);
  assert.doesNotMatch(html, /private/);
  assert.match(html, /Taxa \/ Sequences<\/span>\s*<span class="diag-val ">5</);
  assert.match(html, /Topology test[\s\S]*BOTH, 1000 replicates/);
  assert.match(html, /Models evaluated/);
  assert.match(html, /hyphy3 v3/);
});
