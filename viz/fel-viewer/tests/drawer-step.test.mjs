import { test } from 'node:test';
import assert from 'node:assert/strict';
import { installFakeDOM, FakeElement } from './fake-dom.mjs';
import { Drawer, stepKey, stepPosition } from '../src/components/drawer.js';
import { ResultsTable } from '../src/components/results-table.js';
import { ViewerShell } from '../src/shell/viewer-shell.js';
import { FEL_DRAWER } from '../src/methods/fel.js';

installFakeDOM();

test('stepKey moves within the list and stops at either end', () => {
  const keys = ['Cow', 'Cat', 'Horse'];
  assert.equal(stepKey(keys, 'Cat', 1), 'Horse');
  assert.equal(stepKey(keys, 'Cat', -1), 'Cow');
  assert.equal(stepKey(keys, 'Horse', 1), null);
  assert.equal(stepKey(keys, 'Cow', -1), null);
  assert.equal(stepKey([3, 9, 12], '9', 1), 12);
});

test('a row the filter dropped steps to the first row forward and the last back', () => {
  const keys = ['Cow', 'Cat', 'Horse'];
  assert.equal(stepKey(keys, 'Pig', 1), 'Cow');
  assert.equal(stepKey(keys, 'Pig', -1), 'Horse');
  assert.equal(stepKey([], 'Pig', 1), null);
  assert.equal(stepPosition(keys, 'Cat'), '2 / 3');
  assert.equal(stepPosition(keys, 'Pig'), '– / 3');
});

const BRANCHES = { Cow: { branch: 'Cow' }, Cat: { branch: 'Cat' }, Horse: { branch: 'Horse' }, Pig: { branch: 'Pig' } };

function branchDrawer(siblings, steps) {
  const drawer = new Drawer(new FakeElement(), {
    row: (data, key) => data[key] || null,
    header: (r) => ({ tag: r.branch, tone: 'neutral', badge: '' }),
    tabs: [{ id: 'summary', label: 'Branch', render() {} }],
    siblings,
    onStep: (key) => steps.push(key)
  });
  drawer.setData(BRANCHES);
  return drawer;
}

test('the drawer steps through the siblings it is given and reports its place', () => {
  const steps = [];
  let keys = ['Cow', 'Cat', 'Horse'];
  const drawer = branchDrawer(() => keys, steps);
  assert.ok(drawer.drawer.innerHTML.includes('id="drawer-prev-btn"'));

  drawer.inspectSite('Cow');
  assert.equal(drawer.stepPos.textContent, '1 / 3');
  assert.equal(drawer.prevBtn.disabled, true);
  assert.equal(drawer.nextBtn.disabled, false);

  drawer.nextBtn.fire('click');
  drawer.prevBtn.fire('click');
  assert.deepEqual(steps, ['Cat']);

  keys = ['Cat', 'Horse'];
  drawer.setThreshold(0.01);
  assert.equal(drawer.stepPos.textContent, '– / 2');
  drawer.step(-1);
  assert.deepEqual(steps, ['Cat', 'Horse']);
});

test('arrow keys step an open drawer, but not from a text field', () => {
  const steps = [];
  const drawer = branchDrawer(() => ['Cow', 'Cat', 'Horse'], steps);
  const key = (k, tagName = 'BODY') => drawer.onKeyDown({ key: k, target: { tagName } });
  key('ArrowRight');
  assert.deepEqual(steps, []);
  drawer.inspectSite('Cat');
  key('ArrowRight');
  key('ArrowLeft');
  key('ArrowRight', 'INPUT');
  assert.deepEqual(steps, ['Horse', 'Cow']);
  key('Escape');
  assert.equal(drawer.isOpen, false);
});

test('without siblings the header keeps its FEL markup and has no step buttons', () => {
  const drawer = new Drawer(new FakeElement(), FEL_DRAWER);
  assert.ok(drawer.drawer.innerHTML.includes(`<div class="drawer-title-group">
          <span class="drawer-site-tag" id="drawer-site-tag">`));
  assert.ok(!drawer.drawer.innerHTML.includes('drawer-step'));
  drawer.step(1);
});

test('the shell hands the drawer the table\'s visible rows and steps through its focus path', () => {
  let hooks;
  const rows = Object.values(BRANCHES).map((r, i) => ({ ...r, p: i / 10 }));
  class Parser {
    constructor(json) { this.rows = json.rows; }
    getStats() { return {}; }
  }
  const config = {
    id: 'stepper',
    Parser,
    controller: { seed: 0.15, statistic: 'p' },
    row: (data, key) => data.rows.find(r => r.branch === key) || null,
    kpis: () => ({ setData() {} }),
    hero: () => ({ setData() {} }),
    table: (parent, h) => new ResultsTable(parent, {
      rowKey: 'branch',
      rows: (data) => data.rows,
      columns: [{ id: 'branch', header: 'Branch', cell: (r) => `<td>${r.branch}</td>` }],
      filters: [{ id: 'sig', key: 'sig', label: 'Sig', count: () => 0, test: (r, t) => r.p <= t, subtitle: () => '' }],
      threshold: h.threshold
    }),
    drawer: (container, h) => {
      hooks = h;
      return { setData() {}, inspectSite(key) { this.key = key; }, close() {} };
    },
    diagnostics: () => ''
  };
  const shell = new ViewerShell(new FakeElement(), { method: null });
  shell.mount(config);
  shell.load({ rows });
  assert.deepEqual(hooks.siblings(), ['Cow', 'Cat']);
  shell.setThreshold(0.25);
  assert.deepEqual(hooks.siblings(), ['Cow', 'Cat', 'Horse']);

  hooks.onStep('Horse');
  assert.equal(shell.selectedSite, 'Horse');
  assert.equal(shell.drawer.key, 'Horse');
  assert.deepEqual(shell.container.events.at(-1).detail.raw, rows[2]);
});
