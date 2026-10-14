import { test } from 'node:test';
import assert from 'node:assert/strict';
import { stackBandFrames, codonTickStep, TRACK_PADDING } from '../src/components/site-track.js';
import { EvidenceBand, autoEvidenceDomain, evidenceRuleValue } from '../src/components/evidence-band.js';
import { TickStripBand, tickStripCeiling } from '../src/components/tick-strip-band.js';
import { evidencePlotModes, evidenceTone } from '../src/components/evidence-plot.js';

// Canvas context that records the calls bands make
function recordingContext() {
  const calls = [];
  const ctx = new Proxy({ calls }, {
    get(target, key) {
      if (key in target) return target[key];
      return (...args) => calls.push([key, ...args]);
    },
    set(target, key, value) {
      calls.push([`=${String(key)}`, value]);
      return true;
    }
  });
  return ctx;
}

const frameFor = (numSites, top = 24, height = 200) => ({
  left: 52, top, width: 400, height,
  scaleX: (site) => 52 + ((site - 1) / Math.max(1, numSites - 1)) * 400
});

const track = (rows, extra = {}) => Object.assign({
  rows, threshold: 0.1, data: { isMEME: false }, selectedSite: null, hoveredSite: null
}, extra);

test('the main band keeps the full plot when nothing is stacked', () => {
  assert.deepEqual(stackBandFrames(800, 290, []), [{ left: 52, top: 24, width: 724, height: 228 }]);
});

test('stacked bands sit under the codon axis and shrink the main band', () => {
  const [main, strip] = stackBandFrames(800, 290, [{ height: 20 }]);
  assert.equal(main.height, 228 - 20 - TRACK_PADDING.gap);
  assert.equal(strip.top, main.top + main.height + TRACK_PADDING.axis);
  assert.equal(strip.height, 20);
  assert.equal(strip.top + strip.height, 290 - 20);
});

test('codon ticks widen with the alignment', () => {
  assert.deepEqual([codonTickStep(80), codonTickStep(200), codonTickStep(1200)], [10, 20, 50]);
});

test('the rule sits at -log10 t for p and at t for posteriors', () => {
  assert.equal(evidenceRuleValue('le', 0.01), 2);
  assert.equal(evidenceRuleValue('le', 0), 15);
  assert.equal(evidenceRuleValue('ge', 0.9), 0.9);
});

test('auto domains include zero and are symmetric for a two-sided band', () => {
  assert.deepEqual(autoEvidenceDomain([0.5, 4.2], false), [0, 5]);
  assert.deepEqual(autoEvidenceDomain([-3, 4], false), [-4, 5]);
  assert.deepEqual(autoEvidenceDomain([-3, 4], true), [-5, 5]);
  assert.deepEqual(autoEvidenceDomain([], false), [0, 1]);
});

test('a posterior band uses a fixed axis, explicit ticks and a ge rule', () => {
  const rows = [{ site: 1, pp: 0.2 }, { site: 2, pp: 0.95 }, { site: 3, pp: 0.5 }];
  const band = new EvidenceBand({
    value: (r) => r.pp,
    tone: (r, t) => (r.pp >= t.threshold ? 'div' : null),
    domain: [0, 1],
    ticks: [0, 0.5, 0.9, 1],
    tickFormat: (v) => v.toFixed(2),
    rule: 'ge',
    shape: 'dot'
  });
  const ctx = recordingContext();
  const frame = frameFor(3);
  const t = track(rows, { threshold: 0.9, selectedSite: 2 });
  band.drawBackground(ctx, frame, t);
  const points = band.drawMarks(ctx, frame, t);

  const labels = ctx.calls.filter(c => c[0] === 'fillText').map(c => c[1]);
  assert.deepEqual(labels, ['0.00', '0.50', '0.90', '1.00', 'P ≥ 0.90']);
  assert.equal(frame.scaleY(1), frame.top);
  assert.deepEqual(points.map(p => p.tone), [null, 'div', null]);
  // The selected site gets the ring after its own mark
  const arcs = ctx.calls.filter(c => c[0] === 'arc').map(c => c[3]);
  assert.deepEqual(arcs, [2.5, 2.5, 6, 11]);
});

test('a mirrored band plots the second series below zero with rules on both sides', () => {
  const rows = [
    { site: 1, pPos: 0.01, pNeg: 1 },
    { site: 2, pPos: 0.9, pNeg: 0.001 },
    { site: 3, pPos: 1, pNeg: null }
  ];
  const negLog = (p) => -Math.log10(Math.max(p, 1e-15));
  const band = new EvidenceBand({
    value: (r) => negLog(r.pPos),
    tone: (r, t) => (r.pPos <= t.threshold ? 'div' : null),
    mirror: {
      value: (r) => (r.pNeg == null ? null : negLog(r.pNeg)),
      tone: (r, t) => (r.pNeg <= t.threshold ? 'pur' : null)
    },
    label: ['more non-synonymous ▲', 'more synonymous ▼']
  });
  const ctx = recordingContext();
  const frame = frameFor(3);
  const t = track(rows);
  band.drawBackground(ctx, frame, t);
  const points = band.drawMarks(ctx, frame, t);

  assert.deepEqual([band.yMin, band.yMax], [-4, 4]);
  assert.deepEqual(points.map(p => [p.site, p.yVal === 0 ? 0 : p.yVal, p.tone]),
    [[1, 2, 'div'], [1, 0, null], [2, negLog(0.9), null], [2, -3, 'pur'], [3, 0, null]]);
  const dashed = ctx.calls.filter(c => c[0] === 'setLineDash').length;
  assert.equal(dashed, 2);
  const captions = ctx.calls.filter(c => c[0] === 'fillText').map(c => c[1]);
  assert.ok(captions.includes('more non-synonymous ▲') && captions.includes('more synonymous ▼'));
});

test('the tick strip clamps tick heights at the 95th percentile of counts', () => {
  const counts = Array.from({ length: 40 }, (_, i) => (i < 38 ? 1 : 50));
  assert.equal(tickStripCeiling(counts, 0.95), 1);
  assert.equal(tickStripCeiling([0, 0], 0.95), 1);
  assert.equal(tickStripCeiling([2, 4, 6, 8], 0.95), 8);

  const rows = [
    { site: 1, gc: 0, at: 0, odd: false },
    { site: 2, gc: 2, at: 1, odd: false },
    { site: 3, gc: 8, at: 0, odd: true },
    { site: 4, gc: 0, at: 0, odd: true }
  ];
  const strip = new TickStripBand({ up: r => r.gc, down: r => r.at, flag: r => r.odd, label: 'GC/AT' });
  const frame = frameFor(4, 260, 20);
  const marks = strip.marks(frame, track(rows));
  assert.deepEqual(marks.map(m => [m.site, m.up, m.down, m.flag]),
    [[2, 2.5, 1.25, false], [3, 10, 0, true], [4, 0, 0, true]]);
  const svg = strip.svg(frame, track(rows));
  assert.equal((svg.match(/<circle/g) || []).length, 2);
  assert.ok(svg.includes('>GC/AT</text>'));
});

test('the FEL preset calls sites by direction and MEME sites as episodic', () => {
  const fel = track([], { data: { isMEME: false } });
  const meme = track([], { data: { isMEME: true } });
  const site = (alpha, beta, pValue) => ({ alpha, beta, pValue });
  assert.equal(evidenceTone(site(0.1, 2, 0.05), fel), 'div');
  assert.equal(evidenceTone(site(2, 0.1, 0.05), fel), 'pur');
  assert.equal(evidenceTone(site(1, 1, 0.05), fel), null);
  assert.equal(evidenceTone(site(2, 0.1, 0.5), fel), null);
  assert.equal(evidenceTone(site(2, 0.1, 0.05), meme), 'div');
});

test('the FEL preset keeps its three modes and their domains', () => {
  const rows = [
    { site: 1, negLogP: 0.2, alpha: 1, beta: 0.5, pValue: 0.6 },
    { site: 2, negLogP: 5.0, alpha: 0.1, beta: 6, pValue: 1e-5 }
  ];
  const modes = evidencePlotModes();
  assert.deepEqual(modes.map(m => m.id), ['manhattan', 'signed', 'diff']);
  assert.equal(modes[0].subtitle({ isMEME: true }), '−log₁₀(p) Evidence for Episodic Diversifying Selection');
  const t = track(rows);
  assert.deepEqual(modes[0].band.domain(rows, t), [0, 6]);
  assert.deepEqual(modes[1].band.domain(rows, t), [-6, 6]);
  assert.deepEqual(modes[2].band.domain(rows, t), [-2, 7]);
  assert.equal(modes[2].band.options.rule, null);
  assert.equal(modes[1].band.options.mirrorRule, true);
});
