import { test } from 'node:test';
import assert from 'node:assert/strict';
import { installFakeDOM, FakeElement } from './fake-dom.mjs';
import { EvidenceBand, evidenceBarWidth, evidenceFlagPosition } from '../src/components/evidence-band.js';
import { SiteTrackCanvas, trackColors } from '../src/components/site-track.js';
import { PALETTE } from '../src/data/palette.js';
import { OverviewTrack, revealWindow } from '../src/components/overview-track.js';
import {
  slacTrackModes, slacEvidenceDomain, slacEvidenceTicks, slacPTick, slacMarker, slacNegLog10, SLAC_BRUSH
} from '../src/methods/slac.js';

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

const frameFor = (numSites, width = 400) => ({
  left: 52, top: 24, width, height: 200,
  scaleX: (site) => 52 + ((site - 1) / Math.max(1, numSites - 1)) * width
});

const site = (s, extra = {}) => Object.assign({
  site: s, S: 1, N: 1, pPos: 0.5, pNeg: 0.8, dNmdS: 0.2, noSubs: false, dSUndefined: false
}, extra);

test('the evidence axis is symmetric, at least p = 0.001 either way, and ticks read as p', () => {
  assert.deepEqual(slacEvidenceDomain([site(1, { pPos: 0.2, pNeg: 0.9 })]), [-3, 3]);
  assert.deepEqual(slacEvidenceDomain([site(1, { pPos: 1e-7, pNeg: 1 })]), [-8, 8]);
  assert.deepEqual(slacEvidenceTicks(-3, 3), [-3, -2, -1, 0, 1, 2, 3]);
  assert.deepEqual(slacEvidenceTicks(-8, 8), [-8, -6, -4, -2, 0, 2, 4, 6, 8]);
  assert.deepEqual([0, 1, -2, 3].map(v => slacPTick(v, 3)), ['1', '0.1', '0.01', '0.001']);
  // An axis past 0.001 reads in powers of ten throughout, never mixing the two
  assert.deepEqual([0, 2, -4, 6].map(v => slacPTick(v, 6)), ['1', '1e-2', '1e-4', '1e-6']);
  assert.equal(slacNegLog10(0), 15);
});

test('a degenerate p(dN>dS) = 0 does not set the scale; its flagged mark is held at the top edge', () => {
  const degenerate = site(2, { dSUndefined: true, N: 2, pPos: 0, pNeg: 1 });
  const rows = [site(1, { pPos: 1, pNeg: 0.0086 }), degenerate];
  assert.deepEqual(slacEvidenceDomain(rows), [-3, 3]);
  const [evidence] = slacTrackModes();
  const frame = frameFor(2);
  const ctx = recordingContext();
  const track = { rows, threshold: 0.1, selectedSite: null, hoveredSite: null };
  evidence.band.drawBackground(ctx, frame, track);
  const points = evidence.band.drawMarks(ctx, frame, track);
  const top = points.find(p => p.site === 2 && p.series === 0);
  assert.equal(top.y, frame.top);
  assert.ok(top.clipped);
  assert.equal(top.marker, 'flag');
  assert.ok(!points.find(p => p.site === 1 && p.series === 1).clipped);
  assert.ok(ctx.calls.some(c => c[0] === 'fillText' && c[1] === '⚠'));
});

test('the selection ring sits on the called mark, in its own hue', () => {
  const [evidence] = slacTrackModes();
  const rows = [site(1, { pPos: 1, pNeg: 0.0086 }), site(2)];
  const frame = frameFor(2);
  const ctx = recordingContext();
  const track = { rows, threshold: 0.1, selectedSite: 1, selectedSeries: null, hoveredSite: null };
  evidence.band.drawBackground(ctx, frame, track);
  evidence.band.drawMarks(ctx, frame, track);
  const ring = evidence.band.selectedPoint(track);
  assert.equal(ring.series, 1);
  assert.equal(evidence.band.ringColor(ring), PALETTE.purifying);
  assert.equal(evidence.band.ringColor({ tone: null }), PALETTE.ink);
  assert.equal(new EvidenceBand().ringColor({ tone: 'pur' }), PALETTE.amberVivid);
});

test('canvas colours fall back to PALETTE where the page defines no tokens', () => {
  assert.deepEqual(trackColors(null), PALETTE);
  const themed = recordingContext();
  const track = { rows: [site(1)], threshold: 0.1, selectedSite: null, hoveredSite: null, colors: Object.assign({}, PALETTE, { gridLine: '#2D2D2D', panel: '#1E1E1E' }) };
  const band = new EvidenceBand({ value: () => 0, marker: () => 'hollow', shape: 'dot' });
  const frame = frameFor(1);
  band.drawBackground(themed, frame, track);
  band.drawMarks(themed, frame, track);
  assert.ok(themed.calls.some(c => c[0] === '=strokeStyle' && c[1] === '#2D2D2D'));
  assert.ok(themed.calls.some(c => c[0] === '=fillStyle' && c[1] === '#1E1E1E'));
});

test('a pair caption is set level in the frame corners, in the two hues', () => {
  const [evidence] = slacTrackModes();
  const frame = frameFor(2);
  const ctx = recordingContext();
  evidence.band.drawBackground(ctx, frame, { rows: [site(1)], threshold: 0.1 });
  assert.ok(!ctx.calls.some(c => c[0] === 'rotate'));
  assert.deepEqual(evidence.band.cornerCaptions(frame).map(([text, color, y]) => [text, color, y]), [
    ['more non-synonymous ▲', PALETTE.accent, frame.top + 4],
    ['more synonymous ▼', PALETTE.purifying, frame.top + frame.height - 4]
  ]);
});

test('untestable sites are hollow rings; a degenerate diversifying call is flagged above', () => {
  assert.equal(slacMarker(site(1, { noSubs: true, S: 0, N: 0 }), {}, 0), 'hollow');
  const degenerate = site(2, { dSUndefined: true, N: 2, pPos: 0 });
  assert.equal(slacMarker(degenerate, {}, 0), 'flag');
  assert.equal(slacMarker(degenerate, {}, 1), null);
  assert.equal(slacMarker(site(3, { dSUndefined: true, N: 0, S: 0, noSubs: true }), {}, 0), 'hollow');
  assert.equal(slacMarker(site(4), {}, 0), null);
});

test('evidence mode: p(dN>dS) above zero in the accent, p(dN<dS) below in purifying, one point for an untestable site', () => {
  const [evidence] = slacTrackModes();
  const band = evidence.band;
  const rows = [site(1, { pPos: 0.05, pNeg: 1 }), site(2, { pPos: 1, pNeg: 0.01 }), site(3, { S: 0, N: 0, pPos: 1, pNeg: 1, noSubs: true })];
  const frame = frameFor(3);
  const ctx = recordingContext();
  const track = { rows, threshold: 0.1, selectedSite: null, hoveredSite: null };
  band.drawBackground(ctx, frame, track);
  const points = band.drawMarks(ctx, frame, track);
  assert.deepEqual(points.map(p => [p.site, p.series, p.tone, p.marker]), [
    [1, 0, 'div', null], [1, 1, null, null],
    [2, 0, null, null], [2, 1, 'pur', null],
    [3, 0, null, 'hollow']
  ]);
  assert.ok(Math.abs(points[3].yVal + 2) < 1e-12);
  assert.equal(band.ruleLabel(0.1), 'p ≤ 0.10');
  // The uncalled carpet is drawn faint
  assert.ok(ctx.calls.some(c => c[0] === '=globalAlpha' && c[1] === 0.45));
  // Rules at ±1 (p = 0.1)
  assert.equal(band.tickValues().length, 7);
});

test('effect mode draws bars from zero, coloured only when the site is called', () => {
  const [, effect] = slacTrackModes();
  const band = effect.band;
  assert.equal(band.options.shape, 'bar');
  assert.equal(band.options.rule, null);
  const rows = [site(1, { dNmdS: -2, pNeg: 0.01, pPos: 1 }), site(2, { dNmdS: 1.5 }), site(3, { dNmdS: null, dSUndefined: true })];
  const frame = frameFor(3);
  const ctx = recordingContext();
  const track = { rows, threshold: 0.1, selectedSite: null, hoveredSite: null };
  band.drawBackground(ctx, frame, track);
  const points = band.drawMarks(ctx, frame, track);
  assert.deepEqual(points.map(p => [p.site, p.tone]), [[1, 'pur'], [2, null]]);
  const bars = ctx.calls.filter(c => c[0] === 'fillRect');
  assert.equal(bars.length, 2);
  assert.equal(bars[0][3], evidenceBarWidth(frame));
  const svg = band.svg(frame, track);
  assert.equal((svg.match(/<rect /g) || []).length, 2);
});

test('bars take most of a site, capped; flags sit away from zero', () => {
  assert.equal(evidenceBarWidth(frameFor(3)), 6);
  assert.equal(evidenceBarWidth(frameFor(1000)), 1);
  assert.ok(Math.abs(evidenceBarWidth(frameFor(200, 1000)) - 0.6 * 1000 / 199) < 1e-9);
  assert.deepEqual(evidenceFlagPosition({ x: 5, y: 40 }, 100), { x: 5, y: 29 });
  assert.deepEqual(evidenceFlagPosition({ x: 5, y: 140 }, 100), { x: 5, y: 151 });
});

test('a hollow ring stands in for the mark, and a flag adds a glyph, in canvas and SVG', () => {
  const band = new EvidenceBand({
    value: (r) => r.v,
    tone: (r) => r.tone,
    marker: (r) => r.marker,
    shape: 'dot'
  });
  const rows = [{ site: 1, v: 0, tone: null, marker: 'hollow' }, { site: 2, v: 3, tone: 'div', marker: 'flag' }];
  const frame = frameFor(2);
  const ctx = recordingContext();
  const track = { rows, threshold: 0.1, selectedSite: null, hoveredSite: null };
  band.drawBackground(ctx, frame, track);
  band.drawMarks(ctx, frame, track);
  assert.ok(ctx.calls.some(c => c[0] === 'fillText' && c[1] === '⚠'));
  assert.equal(ctx.calls.filter(c => c[0] === 'stroke').length > 0, true);
  const svg = band.svg(frame, track);
  assert.match(svg, /fill="#FFFFFF" stroke="#94A3B8"/);
  assert.match(svg, />⚠<\/text>/);
  // A band that never sets neutralAlpha keeps its SVG free of opacity attributes
  assert.doesNotMatch(svg, /opacity=/);
});

// A track without a DOM: just the fields layout() and siteRange() read
function bareTrack(numSites, viewport = null) {
  const t = Object.create(SiteTrackCanvas.prototype);
  Object.assign(t, {
    viewport,
    stack: [],
    width: 800,
    height: 300,
    data: {},
    options: { siteCount: () => numSites, gutter: 52 }
  });
  return t;
}

test('a viewport narrows the codon scale to its sites; without one the scale spans the alignment', () => {
  const whole = bareTrack(500).layout()[0];
  assert.equal(whole.scaleX(1), whole.left);
  assert.equal(whole.scaleX(500), whole.left + whole.width);
  const zoomed = bareTrack(500, { start: 101, end: 250 });
  const [frame] = zoomed.layout();
  assert.deepEqual(zoomed.siteRange(), [101, 250]);
  assert.equal(frame.scaleX(101), frame.left);
  assert.equal(frame.scaleX(250), frame.left + frame.width);
  assert.ok(frame.scaleX(100) < frame.left);
  assert.ok(zoomed.windowPad(frame) > 0 && zoomed.windowPad(frame) <= 6);
});

test('revealing a site keeps the window when it holds the site, else centres on it within the alignment', () => {
  const view = { start: 0, end: 150 };
  assert.equal(revealWindow(view, 40, 150, 510), view);
  assert.deepEqual(revealWindow(view, 400, 150, 510), { start: 324, end: 474 });
  assert.deepEqual(revealWindow(view, 505, 150, 510), { start: 360, end: 510 });
});

// A stand-in site track that records what the overview asks of it
function stubTrack(numSites) {
  return {
    canvasContainer: new FakeElement(),
    threshold: 0.1,
    numSites,
    rows: Array.from({ length: numSites }, (_, i) => site(i + 1, { pNeg: i % 7 ? 0.9 : 0.01 })),
    viewports: [],
    selected: [],
    setData() {},
    setViewport(start, end) { this.viewports.push(start == null ? null : [start, end]); },
    selectSite(s) { this.selected.push(s); },
    setThreshold(t) { this.threshold = t; },
    render() {},
    destroy() {}
  };
}

function overview(numSites) {
  FakeElement.prototype.getContext = () => recordingContext();
  globalThis.ResizeObserver = class { observe() {} disconnect() {} };
  const picked = [];
  const track = stubTrack(numSites);
  const view = new OverviewTrack(new FakeElement(), { track: () => track, brush: SLAC_BRUSH, onSiteSelect: (s) => picked.push(s) });
  view.setData({});
  return { view, track, picked };
}

test('long alignments zoom the track to the brushed window; short ones hide the brush', () => {
  const short = overview(187);
  assert.equal(short.view.zoomed, false);
  assert.equal(short.view.brushHost.style.display, 'none');
  assert.deepEqual(short.track.viewports, [null]);

  const long = overview(510);
  assert.equal(long.view.zoomed, true);
  assert.equal(long.view.brushHost.style.display, '');
  assert.deepEqual(long.track.viewports, [[1, 150]]);
  assert.deepEqual(long.view.brush.viewport, { start: 0, end: 150 });
  // The brush counts the purifying calls in the purifying hue
  assert.ok(long.view.brush.significant.length > 0);

  // Focus from the table scrolls the window to the site; a pick on the brush reports it
  long.view.selectSite(400);
  assert.deepEqual(long.track.viewports.at(-1), [325, 474]);
  long.view.pick(30);
  assert.deepEqual(long.track.viewports.at(-1), [1, 150]);
  assert.deepEqual(long.picked, [30]);
  assert.deepEqual(long.track.selected, [400, 30]);
});
