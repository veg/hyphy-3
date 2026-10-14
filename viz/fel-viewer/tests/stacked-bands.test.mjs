import { test } from 'node:test';
import assert from 'node:assert/strict';
import { SiteTrackCanvas, resolveTrackHit, stackBandFrames } from '../src/components/site-track.js';
import { EvidenceBand } from '../src/components/evidence-band.js';
import { MatrixBand, matrixColumns, matrixExtent } from '../src/components/matrix-band.js';
import { SegmentBand, segmentBreaks } from '../src/components/segment-band.js';

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

const frameFor = (numSites, top = 24, height = 200, width = 400) => ({
  left: 52, top, width, height,
  scaleX: (site) => 52 + ((site - 1) / Math.max(1, numSites - 1)) * width
});

const track = (rows, extra = {}) => Object.assign({
  rows, numSites: rows.length, threshold: 0.1, data: {}, selectedSite: null, hoveredSite: null
}, extra);

test('a pointer inside a band with regions resolves to its region, never a main point', () => {
  const hit = {
    points: [{ site: 1, x: 60, y: 100 }, { site: 2, x: 100, y: 100 }],
    regions: [{ x0: 52, x1: 80, y0: 260, y1: 274, key: 1, kind: 'site' }, { x0: 0, x1: 50, y0: 260, y1: 274, key: 'pI', kind: 'row' }],
    hitFrames: [{ top: 260, height: 28 }]
  };
  assert.equal(resolveTrackHit(hit, 60, 265, 12).region.key, 1);
  assert.equal(resolveTrackHit(hit, 10, 265, 12).region.kind, 'row');
  // Inside the band's frame but on a blank cell: nothing, not the nearest main point
  assert.equal(resolveTrackHit(hit, 100, 280, 12), null);
  // Above the band: the main band's nearest point by x
  assert.equal(resolveTrackHit(hit, 98, 40, 12).point.site, 2);
  assert.equal(resolveTrackHit(hit, 160, 40, 12), null);
});

test('a band that returns no regions does not capture the pointer', () => {
  const hit = { points: [{ site: 3, x: 100, y: 50 }], regions: [], hitFrames: [] };
  assert.equal(resolveTrackHit(hit, 101, 270, 12).point.site, 3);
});

test('the two series of a mirrored band are picked by pointer y', () => {
  const hit = {
    points: [{ site: 1, x: 60, y: 50, series: 0 }, { site: 1, x: 60, y: 150, series: 1 }],
    regions: [],
    hitFrames: []
  };
  assert.equal(resolveTrackHit(hit, 61, 60, 12).point.series, 0);
  assert.equal(resolveTrackHit(hit, 61, 140, 12).point.series, 1);
});

test('a mirrored band rings only the selected series of the selected site', () => {
  const rows = [{ site: 1, up: 2, down: 1 }, { site: 2, up: 3, down: 4 }];
  const band = new EvidenceBand({
    value: r => r.up,
    tone: () => null,
    mirror: { value: r => r.down, tone: () => null }
  });
  const frame = frameFor(2);
  const ring = (selectedSeries) => {
    const ctx = recordingContext();
    const t = track(rows, { selectedSite: 2, selectedSeries });
    band.drawBackground(ctx, frame, t);
    band.drawMarks(ctx, frame, t);
    return ctx.calls.filter(c => c[0] === 'arc' && c[3] === 11).map(c => c[2]);
  };
  assert.deepEqual(ring(1), [frame.scaleY(-4)]);
  assert.deepEqual(ring(0), [frame.scaleY(3)]);
  assert.deepEqual(ring(null), [frame.scaleY(3)]);
});

// A SiteTrackCanvas without a DOM: just the state its click handler reads
function bareTrack(options, regions) {
  const t = Object.create(SiteTrackCanvas.prototype);
  t.options = Object.assign({ sites: d => d.rows, siteCount: d => d.rows.length }, options);
  t.data = { rows: [{ site: 1 }, { site: 2, p: 0.01 }] };
  t.canvas = { getBoundingClientRect: () => ({ left: 0, top: 0 }) };
  t.plotWidth = 400;
  t.points = [];
  t.regions = regions;
  t.hitFrames = [{ top: 260, height: 28 }];
  t.selectedSite = null;
  t.render = () => {};
  return t;
}

test('clicking a site region selects the site; clicking another region goes to its band', () => {
  const selected = [];
  const regionsSelected = [];
  const band = { select: (region) => selected.push(region.key) };
  const regions = [
    { x0: 52, x1: 80, y0: 260, y1: 274, key: 2, kind: 'site', band },
    { x0: 0, x1: 50, y0: 260, y1: 274, key: 'pI', kind: 'row', band }
  ];
  const sites = [];
  const t = bareTrack({
    onSiteSelect: (site, row) => sites.push([site, row]),
    onRegionSelect: region => regionsSelected.push(region.key)
  }, regions);

  t.handleClick({ clientX: 60, clientY: 265 });
  assert.deepEqual(sites, [[2, { site: 2, p: 0.01 }]]);
  assert.equal(t.selectedSite, 2);
  assert.deepEqual(selected, []);

  t.handleClick({ clientX: 10, clientY: 265 });
  assert.deepEqual(selected, ['pI']);
  assert.deepEqual(regionsSelected, ['pI']);
  assert.equal(t.selectedRegion.key, 'pI');
  assert.equal(t.selectedSite, 2);
});

test('a wider gutter moves every frame right by the same amount', () => {
  const [main, strip] = stackBandFrames(800, 290, [{ height: 20 }], 110);
  assert.equal(main.left, 110);
  assert.equal(strip.left, 110);
  assert.equal(main.width, 800 - 110 - 24);
});

const PROPS = [{ key: 'h', label: 'Hydro' }, { key: 'v', label: 'Volume' }, { key: 'bg', label: 'Background', muted: true }];

test('the matrix extent comes from the data, not a fixed domain', () => {
  const sites = [{ site: 1, h: 0.4, v: -1.9 }, { site: 2, h: 1.2, v: null }];
  const value = (s, row) => s[row.key] ?? null;
  assert.equal(matrixExtent(sites, PROPS, value), 1.9);
  assert.equal(matrixExtent([{ site: 1 }], PROPS, value), 1);
});

test('matrix columns bin to pixels once sites are narrower than a pixel', () => {
  const few = matrixColumns([{ site: 1 }, { site: 2 }, { site: 3 }], frameFor(3), 3);
  assert.deepEqual(few.map(c => [c.x0, c.x1]), [[52, 152], [152, 352], [352, 452]]);

  const many = Array.from({ length: 1000 }, (_, i) => ({ site: i + 1 }));
  const binned = matrixColumns(many, frameFor(1000), 1000);
  assert.ok(binned.length <= 400);
  assert.ok(binned.every(c => c.x1 - c.x0 === 1 && c.members.length >= 2));
  assert.equal(binned.reduce((n, c) => n + c.members.length, 0), 1000);
});

test('each binned cell shows its strongest member: tier first, then magnitude', () => {
  const many = Array.from({ length: 1000 }, (_, i) => ({ site: i + 1, h: 0.1, tier: 0 }));
  many[1] = { site: 2, h: -0.2, tier: 2 };
  many[2] = { site: 3, h: 0.9, tier: 0 };
  const band = new MatrixBand({
    rows: [PROPS[0]],
    value: (s, row) => s[row.key],
    emphasis: s => s.tier
  });
  const cells = band.cells(frameFor(1000), track(many));
  const first = cells.find(c => c.x0 === 52);
  assert.equal(first.site.site, 2);
  assert.equal(first.fill, band.options.negative);
});

test('matrix cells colour by sign, fade by tier and mute background rows', () => {
  const sites = [{ site: 1, h: 2, v: -1, bg: 2 }, { site: 2, h: 1, v: null, bg: 1 }];
  const band = new MatrixBand({
    rows: PROPS,
    value: (s, row) => s[row.key],
    emphasis: s => (s.site === 1 ? 2 : 0)
  });
  const cells = band.cells(frameFor(2), track(sites));
  const at = (site, key) => cells.find(c => c.site.site === site && c.row.key === key);
  assert.equal(at(1, 'h').fill, band.options.positive);
  assert.equal(at(1, 'h').alpha, 1);
  assert.equal(at(1, 'v').fill, band.options.negative);
  assert.equal(at(1, 'v').alpha, 0.5);
  assert.equal(at(2, 'h').alpha, 0.5 * 0.15);
  assert.equal(at(1, 'bg').alpha, 0.4);
  assert.equal(at(2, 'v'), undefined);
});

test('matrix rows resolve per load, return label and cell regions, and label clicks toggle focus', () => {
  const focus = [];
  const band = new MatrixBand({
    rows: data => data.props,
    value: (s, row) => s[row.key] ?? null,
    onFocus: key => focus.push(key)
  });
  const sites = [{ site: 1, h: 1, v: 0.5 }, { site: 2, h: -1 }];
  const t = track(sites, { data: { props: PROPS.slice(0, 2) } });
  band.measure(t);
  assert.equal(band.height, 28);

  const frame = frameFor(2, 260, band.height);
  const regions = band.draw(recordingContext(), frame, t);
  const rows = regions.filter(r => r.kind === 'row');
  assert.deepEqual(rows.map(r => [r.key, r.y0, r.y1, r.x1]), [['h', 260, 274, 50], ['v', 274, 288, 50]]);
  const cells = regions.filter(r => r.kind === 'site');
  assert.deepEqual(cells.map(r => [r.key, r.row.key]), [[1, 'h'], [2, 'h'], [1, 'v']]);
  assert.equal(cells[1].raw, sites[1]);

  band.select(rows[1]);
  assert.equal(band.focused, 'v');
  band.select(rows[1]);
  assert.equal(band.focused, null);
  band.select(cells[0]);
  assert.deepEqual(focus, ['v', null]);
});

const SEGMENTS = [{ key: 0, start: 1, end: 460 }, { key: 1, start: 461, end: 762 }];

test('breakpoints sit half a site past each left segment', () => {
  assert.deepEqual(segmentBreaks(SEGMENTS), [{ index: 1, at: 460.5 }]);
  assert.deepEqual(segmentBreaks([SEGMENTS[0]]), []);
});

test('segment blocks are clickable regions, with rules through the main band', () => {
  const band = new SegmentBand({ segments: data => data.segments });
  const main = frameFor(762, 24, 200, 700);
  const frame = frameFor(762, 248, 28, 700);
  const ctx = recordingContext();
  const t = track([], { numSites: 762, data: { segments: SEGMENTS }, frames: [main, frame] });
  const regions = band.draw(ctx, frame, t);

  assert.deepEqual(regions.map(r => [r.kind, r.key, r.x0, r.y0, r.y1]), [
    ['segment', 0, 52, 248, 276],
    ['segment', 1, frame.scaleX(460.5), 248, 276]
  ]);
  assert.equal(regions[1].x1, 752);
  const rule = ctx.calls.find(c => c[0] === 'moveTo' && c[1] === frame.scaleX(460.5));
  assert.equal(rule[2], main.top);
  const texts = ctx.calls.filter(c => c[0] === 'fillText').map(c => c[1]);
  assert.deepEqual(texts, ['Segment 1 · 1–460', 'Segment 2 · 461–762', 'BP 1']);
  assert.match(band.regionTooltip(regions[1], t), /Segment 2 · 461–762[\s\S]*302/);
});

test('a single segment reads as no recombination and draws no rule', () => {
  const band = new SegmentBand({ segments: () => [{ key: 0, start: 1, end: 762 }] });
  const main = frameFor(762, 24, 200, 700);
  const frame = frameFor(762, 248, 28, 700);
  const ctx = recordingContext();
  band.draw(ctx, frame, track([], { numSites: 762, frames: [main, frame] }));
  assert.deepEqual(ctx.calls.filter(c => c[0] === 'fillText').map(c => c[1]),
    ['No recombination detected — single tree fits best']);
  assert.equal(ctx.calls.filter(c => c[0] === 'moveTo').length, 0);
});
