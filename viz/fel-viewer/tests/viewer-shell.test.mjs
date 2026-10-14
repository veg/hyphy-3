import { test } from 'node:test';
import assert from 'node:assert/strict';
import { ViewerShell, effectiveThreshold, thresholdAttribute } from '../src/shell/viewer-shell.js';
import { Drawer } from '../src/components/drawer.js';
import { thresholdLabel } from '../src/components/kpi-bar.js';
import { FELData, runPThreshold } from '../src/data/fel-data.js';
import { SITE_P_THRESHOLD } from '../src/methods/fel.js';

// Just enough DOM for the shell and the drawer to build their layout
class FakeElement {
  constructor() {
    this.children = [];
    this.events = [];
    this.style = {};
    this.classList = { add() {}, remove() {}, toggle() {} };
  }
  appendChild(el) { this.children.push(el); return el; }
  prepend(el) { this.children.unshift(el); }
  replaceChildren() { this.children = []; }
  remove() {}
  addEventListener() {}
  querySelector() { return new FakeElement(); }
  dispatchEvent(e) { this.events.push(e); }
}
globalThis.document = { createElement: () => new FakeElement() };
globalThis.window = { addEventListener() {}, removeEventListener() {} };

// A view that records what the shell fans out to it
class RecordingView {
  constructor(hooks, log) {
    this.log = log;
    this.hooks = hooks;
    log.push(['mount', hooks.threshold, hooks.statistic]);
  }
  setData() {}
  setThreshold(t) { this.log.push(['threshold', t]); }
  setStatistic(s) { this.log.push(['statistic', s]); }
}

// A branch-keyed method: its parser has no getSite, only the config's row lookup
function branchConfig(controller, log = []) {
  class Parser {
    constructor(json) {
      this.branches = json.branches;
      this.runThreshold = json.runThreshold ?? null;
    }
    getStats() { return {}; }
  }
  return {
    id: 'branchy',
    Parser,
    controller,
    row: (data, key) => data.branches.find(b => b.branch === key) || null,
    kpis: (parent, hooks) => new RecordingView(hooks, log),
    hero: () => ({ setData() {} }),
    table: () => ({ setData() {} }),
    drawer: () => ({ setData() {}, inspectSite() {}, close() {} }),
    diagnostics: () => ''
  };
}
const GE_CONTROLLER = { statistic: 'PP', direction: 'ge', seed: 0.9 };
const Q_CONTROLLER = { statistic: 'p', direction: 'le', seed: 0.05 };
const branches = [{ branch: 'Node3', p: 0.01 }, { branch: 'Pig', p: 0.4 }];

function shellWith(config, options = {}) {
  const shell = new ViewerShell(new FakeElement(), options);
  shell.mount(config);
  return shell;
}

test('the threshold comes from the caller, else the run, else the controller seed', () => {
  const config = { controller: GE_CONTROLLER };
  assert.equal(effectiveThreshold(0.2, { runThreshold: 0.5 }, config), 0.2);
  assert.equal(effectiveThreshold(null, { runThreshold: 0.5 }, config), 0.5);
  assert.equal(effectiveThreshold(undefined, { runThreshold: null }, config), 0.9);
  assert.equal(effectiveThreshold(null, null, config), 0.9);
  assert.equal(effectiveThreshold(null, null, null), null);
});

test('a threshold attribute is explicit only when it is a positive number', () => {
  assert.equal(thresholdAttribute('0.05'), 0.05);
  assert.equal(thresholdAttribute(null), null);
  assert.equal(thresholdAttribute(''), null);
  assert.equal(thresholdAttribute('abc'), null);
  assert.equal(thresholdAttribute('0'), null);
});

test('FEL and MEME files expose the cutoff the run used, v3 or legacy', () => {
  assert.equal(runPThreshold({ pvalue_threshold: 0.05 }), 0.05);
  assert.equal(runPThreshold({ pvalue: 0.1 }), 0.1);
  assert.equal(runPThreshold({ 'p-value': 0.01 }), 0.01);
  assert.equal(runPThreshold({}), null);
  assert.equal(runPThreshold(undefined), null);
  assert.equal(new FELData({ schema_version: '3', analysis: { settings: { pvalue_threshold: 0.05 } } }).runThreshold, 0.05);
  assert.equal(new FELData({ schema_version: '3', analysis: {} }).runThreshold, null);
});

test('views mount at the controller seed, not a p-value default', () => {
  const log = [];
  const shell = shellWith(branchConfig(GE_CONTROLLER, log));
  assert.equal(shell.threshold, 0.9);
  assert.equal(shell.statistic, 'PP');
  assert.deepEqual(log, [['mount', 0.9, 'PP']]);
});

test('each load reseeds from the run, and an explicit threshold wins', () => {
  const log = [];
  const shell = shellWith(branchConfig(GE_CONTROLLER, log));
  shell.load({ branches, runThreshold: 0.95 });
  assert.equal(shell.threshold, 0.95);
  shell.setThreshold(0.5);
  shell.load({ branches });
  assert.equal(shell.threshold, 0.9);
  assert.deepEqual(log.slice(1), [['threshold', 0.95], ['threshold', 0.5], ['threshold', 0.9]]);

  const pinned = shellWith(branchConfig(GE_CONTROLLER), { threshold: 0.8 });
  pinned.load({ branches, runThreshold: 0.95 });
  assert.equal(pinned.threshold, 0.8);
});

test('a remount for another method does not carry the previous threshold over', () => {
  const shell = shellWith(branchConfig(Q_CONTROLLER));
  shell.setThreshold(0.01);
  shell.mount(branchConfig(GE_CONTROLLER));
  assert.equal(shell.threshold, 0.9);
  shell.mount(branchConfig(Q_CONTROLLER), { runThreshold: 0.2 });
  assert.equal(shell.threshold, 0.2);
});

test('the statistic is shared state fanned out to the views that use it', () => {
  const log = [];
  const shell = shellWith(branchConfig(Q_CONTROLLER, log));
  assert.equal(shell.statistic, 'p');
  shell.setStatistic('q');
  assert.equal(shell.statistic, 'q');
  assert.deepEqual(log.at(-1), ['statistic', 'q']);
  assert.equal(shellWith(branchConfig(Q_CONTROLLER), { statistic: 'q' }).statistic, 'q');
});

test('a statistic switch in the KPI bar reaches the shell, its views and its host', () => {
  const log = [];
  const shell = shellWith(branchConfig(Q_CONTROLLER, log));
  shell.kpiBar.hooks.onStatisticChange('q');
  assert.equal(shell.statistic, 'q');
  assert.deepEqual(log.at(-1), ['statistic', 'q']);
  const event = shell.container.events.at(-1);
  assert.equal(event.type, 'statistic-change');
  assert.deepEqual(event.detail, { statistic: 'q' });
});

test('selection events carry the row from the config lookup, so branch keys work', () => {
  const picked = [];
  const shell = shellWith(branchConfig(Q_CONTROLLER), { onSiteSelect: (key, row) => picked.push([key, row]) });
  shell.load({ branches });
  shell.selectSite('Node3');
  const event = shell.container.events.find(e => e.type === 'site-select');
  assert.deepEqual(event.detail.raw, branches[0]);
  assert.deepEqual(picked, [['Node3', branches[0]]]);
});

test('the drawer redraws the open row with the statistic the shell switched to', () => {
  const seen = [];
  const drawer = new Drawer(new FakeElement(), {
    row: (data, key) => data.branches.find(b => b.branch === key),
    header: (row, ctx) => ({ tag: row.branch, tone: 'neutral', badge: ctx.statistic }),
    tabs: [{ id: 'summary', label: 'Summary', render: (panel, row, ctx) => seen.push([row.branch, ctx.statistic]) }],
    statistic: 'p'
  });
  drawer.setData({ branches });
  drawer.inspectSite('Node3');
  drawer.setStatistic('q');
  assert.deepEqual(seen, [['Node3', 'p'], ['Node3', 'q']]);
  assert.equal(drawer.regimeBadge.textContent, 'q');
});

test('the cutoff label names the statistic in force', () => {
  assert.equal(thresholdLabel(SITE_P_THRESHOLD, 0.05), 'p ≤ 0.05');
  assert.equal(thresholdLabel(SITE_P_THRESHOLD, 0.05, 'q'), 'q ≤ 0.05');
});

test('a config that names its format hands the label to the host on load, outside the layout', () => {
  const loaded = (shell) => shell.container.events.find(e => e.type === 'data-loaded').detail.formatLabel;
  const config = Object.assign(branchConfig(GE_CONTROLLER), { formatLabel: (data) => `${data.branches.length} branches` });
  const shell = shellWith(config);
  shell.load({ branches });
  assert.equal(loaded(shell), '2 branches');
  assert.ok(!shell.layout.children.some(el => /chip/.test(el.className)));
  const plain = shellWith(branchConfig(GE_CONTROLLER));
  plain.load({ branches });
  assert.equal(loaded(plain), null);
});

test('a drawer with one tab hides its tab strip; two tabs show it', () => {
  const tab = (id) => ({ id, label: id, render: () => {} });
  const strip = (tabs) => new Drawer(new FakeElement(), { tabs }).drawer.innerHTML.match(/<div class="segmented-control drawer-nav-tabs"[^>]*>/)[0];
  assert.equal(strip([tab('a')]), '<div class="segmented-control drawer-nav-tabs" style="display: none;">');
  assert.equal(strip([tab('a'), tab('b')]), '<div class="segmented-control drawer-nav-tabs">');
});

test('a diagnostics panel is drawn once per load, not on every threshold change', () => {
  const drawn = [];
  const config = Object.assign(branchConfig(GE_CONTROLLER), {
    diagnosticsTitle: 'Run details',
    diagnosticsPanel: (el, data) => drawn.push(data.branches.length)
  });
  const shell = shellWith(config);
  assert.match(shell.diagWrapper.innerHTML, /<span>Run details<\/span>/);
  assert.match(shell.diagWrapper.innerHTML, /<\/div>\n        <div class="diag-panel" id="diag-panel"><\/div>/);
  shell.load({ branches });
  shell.setThreshold(0.95);
  shell.load({ branches: branches.slice(1) });
  assert.deepEqual(drawn, [2, 1]);

  // Without one the card keeps its original markup
  const plain = shellWith(branchConfig(GE_CONTROLLER));
  assert.doesNotMatch(plain.diagWrapper.innerHTML, /diag-panel/);
  assert.match(plain.diagWrapper.innerHTML, /<span>Model Fits &amp; Phylogeny Diagnostics<\/span>/);
});
