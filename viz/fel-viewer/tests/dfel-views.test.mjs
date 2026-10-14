import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import { installFakeDOM, FakeElement } from './fake-dom.mjs';
import { SiteTrackCanvas } from '../src/components/site-track.js';
import { EvidenceBand } from '../src/components/evidence-band.js';
import { TickStripBand } from '../src/components/tick-strip-band.js';
import { BranchTree, diamondPoints } from '../src/components/branch-tree.js';
import {
  compositionTableMarkup, compositionP, gcOutliers, gcDomain, gcStripSVG, COMPOSITION_CSS_COLORS
} from '../src/components/composition-panel.js';
import { rootSupportSVG, rootCredibleCount, rootEdgeLabel } from '../src/components/root-support-bars.js';
import { DFELData } from '../src/data/dfel-data.js';
import { PALETTE } from '../src/data/palette.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';
const fixture = name => JSON.parse(readFileSync(new URL(name, FIXTURES), 'utf8'));

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

const SCREEN = {
  overall: { chi2: 33.33, df: 27, p: 0.186 },
  cp1: { chi2: 8.0, df: 27, p: 0.9999 },
  cp2: { chi2: 18.6, df: 27, p: 0.885 },
  cp3: { chi2: 39.2, df: 27, p: 0.0405 }
};

test('the screen table marks the deciding tests and bolds a p under the heterogeneity threshold', () => {
  const html = compositionTableMarkup(SCREEN, 0.05);
  assert.match(html, /<td>Overall •<\/td>/);
  assert.match(html, /<td>CP1<\/td>/);
  assert.match(html, /<td>CP3 •<\/td><td class="font-mono">39\.20<\/td><td class="font-mono">27<\/td><td class="font-mono font-bold">0\.041<\/td>/);
  assert.doesNotMatch(html, /font-bold">0\.186/);
  assert.deepEqual([compositionP(0.0004), compositionP(null)], ['4.0e-4', '—']);
});

test('the GC strip draws a filled GC3 and a hollow overall dot per taxon and labels the outliers', { skip: noFixtures }, () => {
  const taxa = new DFELData(fixture('dfel_cd2.json')).composition.taxa;
  assert.deepEqual([...gcOutliers(taxa)].sort(), ['Cow', 'Human', 'Mouse', 'Rat']);
  const [lo, hi] = gcDomain(taxa);
  assert.ok(Math.abs(lo - (0.3928 - 0.02)) < 1e-3 && Math.abs(hi - (0.6298 + 0.02)) < 1e-3);
  const { svg, height } = gcStripSVG(taxa, { width: 400 });
  assert.equal((svg.match(/<circle/g) || []).length, 20);
  // A hollow dot carries both paints in one style attribute, so its outline survives
  assert.equal((svg.match(/style="fill: var\(--panel\); stroke: var\(--slate-light\)"/g) || []).length, 10);
  assert.ok(!/<[^>]*style="[^"]*"[^>]*style="/.test(svg));
  assert.match(svg, /font-weight="bold" style="fill: var\(--ink\)">0\.63<\/text>/);
  assert.equal(height, 14 + 10 * 16 + 26);
  // Exports carry literal colours
  assert.ok(gcStripSVG(taxa, { colors: PALETTE }).svg.includes(`fill: ${PALETTE.ink}`));
  assert.ok(!gcStripSVG(taxa, { colors: PALETTE }).svg.includes('var(--'));
});

test('regime marks label the strip, a crowded label moving up a line', () => {
  const taxa = [{ taxon: 'a', gcCP3: 0.4, gcOverall: 0.38 }, { taxon: 'b', gcCP3: 0.6, gcOverall: 0.5 }];
  const marks = [{ value: 0.47, label: 'AT-drift' }, { value: 0.48, label: 'Neutral' }, { value: 0.6, label: 'GC-drift' }];
  const { svg } = gcStripSVG(taxa, { width: 400, marks, colors: COMPOSITION_CSS_COLORS });
  const y = (label) => Number(svg.match(new RegExp(`y="([\\d.]+)"[^>]*>${label}</text>`))[1]);
  assert.equal(y('AT-drift'), 36 - 13);
  assert.equal(y('Neutral'), 36 - 25);
  assert.equal(y('GC-drift'), 36 - 13);
});

test('root bars: the chosen edge by name in ink, the 95% rule, leaf glyphs and the top-N limit', () => {
  const candidates = [
    { branch: 'Node1 - Node2', weight: 0.5, deltaLogL: 0, terminal: false },
    { branch: 'Cow', weight: 0.46, deltaLogL: 0.08, terminal: true },
    { branch: 'Pig', weight: 0.04, deltaLogL: 2.5, terminal: true }
  ];
  assert.equal(rootCredibleCount(candidates), 2);
  assert.equal(rootEdgeLabel('Node1 - Node2'), 'Node1 – Node2');
  // An outgroup root need not be rank 1: the bar is matched by name
  const { svg } = rootSupportSVG(candidates, { chosen: 'Cow', colors: PALETTE });
  assert.match(svg, new RegExp(`font-weight="bold" style="fill: ${PALETTE.ink}">Cow</text>`));
  assert.match(svg, new RegExp(`style="fill: ${PALETTE.stone}">Node1 – Node2</text>`));
  assert.match(svg, />w = 46\.0% · ΔlnL 0\.08<\/text>/);
  assert.equal((svg.match(/stroke-dasharray/g) || []).length, 1);
  assert.equal((svg.match(/<circle/g) || []).length, 2);
  const top = rootSupportSVG(candidates, { chosen: 'Cow', limit: 1 });
  assert.equal((top.svg.match(/<rect/g) || []).length, 2);
  assert.equal(top.height, 6 + 20 + 6);
});

function stripTrack(toggle) {
  const band = new EvidenceBand({ value: (r) => r.y });
  const strip = new TickStripBand({ up: (r) => r.gc, down: (r) => r.at });
  return new SiteTrackCanvas(new FakeElement(), { band, stack: [strip], stackToggle: toggle });
}

const trackData = (K) => ({
  K,
  siteResults: [{ site: 1, y: 1, gc: 2, at: 0 }, { site: 2, y: 2, gc: 0, at: 1 }],
  numSites: 2
});

test('a stack toggle hides the stacked bands by default when shown(data) says so, and the main band takes the room', () => {
  const toggle = { label: 'GC/AT direction', shown: (data) => data.K > 1, note: (data) => `note K=${data.K}` };
  const track = stripTrack(toggle);
  assert.match(track.wrapper.innerHTML, /id="btn-stack-toggle"[^>]*>GC\/AT direction<\/button>/);
  assert.match(track.wrapper.innerHTML, /id="plot-stack-note"/);

  track.setData(trackData(1));
  assert.equal(track.stackShown, false);
  assert.deepEqual(track.activeStack, []);
  assert.equal(track.layout().length, 1);
  assert.equal(track.stackButton.getAttribute('aria-pressed'), 'false');
  assert.equal(track.wrapper.querySelector('#plot-stack-note').style.display, 'none');

  track.stackButton.fire('click');
  assert.equal(track.stackShown, true);
  assert.equal(track.layout().length, 2);
  assert.equal(track.wrapper.querySelector('#plot-stack-note').innerHTML, 'note K=1');

  track.setData(trackData(3));
  assert.equal(track.stackShown, true);
});

test('an unavailable stack hides its toggle and stays off; without a toggle the markup is unchanged', () => {
  const track = stripTrack({ label: 'GC/AT', shown: () => true, available: () => false });
  track.setData(trackData(3));
  assert.equal(track.stackShown, false);
  assert.equal(track.stackButton.style.display, 'none');

  const plain = stripTrack(null);
  assert.doesNotMatch(plain.wrapper.innerHTML, /btn-stack-toggle|plot-stack-note/);
  assert.match(plain.wrapper.innerHTML, /<\/div>\n    $/);
  plain.setData(trackData(1));
  assert.equal(plain.layout().length, 2);
});

test('the tick strip draws with the track\'s themed colours', () => {
  const strip = new TickStripBand({ up: (r) => r.gc, down: () => 0, flag: () => true, label: 'GC/AT' });
  const ctx = recordingContext();
  const frame = { left: 52, top: 0, width: 400, height: 20, scaleX: (s) => 52 + s * 10 };
  strip.draw(ctx, frame, { rows: [{ site: 1, gc: 1 }], colors: Object.assign({}, PALETTE, { ink: '#F8FAFC', gridLine: '#2D2D2D' }) });
  const fills = ctx.calls.filter(c => c[0] === '=fillStyle' || c[0] === '=strokeStyle').map(c => c[1]);
  assert.ok(fills.includes('#F8FAFC'));
  assert.ok(fills.includes('#2D2D2D'));
});

test('a placed root draws as a diamond; branch styles get the palette in use', () => {
  const styled = [];
  const container = new FakeElement();
  container.clientWidth = 300;
  globalThis.window = {};
  const tree = new BranchTree(container, {
    rootShape: 'diamond',
    style: (node, meta, colors) => { styled.push(colors); return {}; }
  });
  tree.setTree('((a:1,b:1)n1:1,c:2)Root;');
  const { dots } = tree.scene(PALETTE);
  assert.equal(dots[0].shape, 'diamond');
  assert.deepEqual(diamondPoints({ x: 10, y: 20, r: 5 }), [[10, 15], [15, 20], [10, 25], [5, 20]]);
  assert.match(tree.svg(), /<polygon points="[^"]+" fill="#1A1A18"\/>/);
  assert.ok(styled.length > 0 && styled.every(c => c.ink));
  globalThis.window = { addEventListener() {}, removeEventListener() {} };
});
