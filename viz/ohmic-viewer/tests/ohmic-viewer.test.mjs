import { test } from 'node:test';
import assert from 'node:assert/strict';
import { OhmicDataModel, bridgePassesTau, countActiveBridges, bridgeCountText } from '../src/data/ohmic-data.js';
import { OhmicKpiBar } from '../src/components/ohmic-kpi-bar.js';

// Just enough DOM for the KPI bar: innerHTML is searched for ids, listeners are recorded
class FakeContainer {
  constructor() { this.els = new Map(); this._html = ''; }
  set innerHTML(html) {
    this._html = html;
    this.els.clear();
    for (const m of html.matchAll(/<(\w+)[^>]*\bid="([^"]+)"[^>]*>([^<]*)/g)) {
      this.els.set(m[2], { textContent: m[3], listeners: {}, addEventListener(t, fn) { this.listeners[t] = fn; } });
    }
  }
  get innerHTML() { return this._html; }
  querySelector(sel) { return this.els.get(sel.slice(1)) ?? null; }
}

const stresses = [0.0, 0.2, 0.5, 0.5, 0.9, 1.0];
function makeModel(extra = {}) {
  return new OhmicDataModel({
    schema_version: '3.0.0',
    dataset: { taxa_count: 8, sites_count: 120 },
    wheatstone_bridges: stresses.map((s, i) => ({ index: i, shunt_stress: s, shuntStress: s })),
    ...extra
  });
}
// The model's own bridge objects, whatever fields it normalises to
function bridgesWith(model) {
  return model.bridges.map((b, i) => ({ ...b, shuntStress: stresses[i] }));
}

test('bridgePassesTau is inclusive at the threshold', () => {
  assert.equal(bridgePassesTau({ shuntStress: 0.5 }, 0.5), true);
  assert.equal(bridgePassesTau({ shuntStress: 0.49 }, 0.5), false);
});

test('KPI card follows the slider, including the max, and matches the table filter', () => {
  const model = makeModel();
  model.bridges = bridgesWith(model);
  const container = new FakeContainer();
  const seen = [];
  const bar = new OhmicKpiBar(container, { onTauThresholdChange: (t) => seen.push(t) });
  bar.render(model, 0);
  assert.equal(container.querySelector('#ohmic-showing').textContent, 'Showing 6 of 6 bridges');

  for (const tau of [0.2, 0.5, 0.52, 0.9, 1.0]) {
    container.querySelector('#ohmic-tau-slider').listeners.input({ target: { value: String(tau) } });
    // The table's filter: bridges with shuntStress >= tau
    const tableRows = model.bridges.filter((b) => b.shuntStress >= tau).length;
    assert.equal(countActiveBridges(model.bridges, tau), tableRows);
    assert.equal(container.querySelector('#ohmic-showing').textContent, bridgeCountText(tableRows, 6));
    assert.equal(container.querySelector('#ohmic-active-count').textContent, String(tableRows));
    assert.equal(container.querySelector('#ohmic-tau-val').textContent, `≥ ${tau.toFixed(2)}`);
  }
  assert.deepEqual(seen, [0.2, 0.5, 0.52, 0.9, 1.0]);
});

test('card reports zero bridges when none reach the threshold', () => {
  const model = makeModel();
  model.bridges = model.bridges.map((b) => ({ ...b, shuntStress: 0.3 }));
  const container = new FakeContainer();
  const bar = new OhmicKpiBar(container);
  bar.render(model, 0);
  container.querySelector('#ohmic-tau-slider').listeners.input({ target: { value: '1' } });
  assert.equal(container.querySelector('#ohmic-showing').textContent, 'Showing 0 of 6 bridges');
});

test('format label names the hyphy3 producer', () => {
  assert.equal(makeModel().formatLabel, 'hyphy3 v3');
  assert.equal(makeModel({ input: { 'file name': 'x.fna', 'number of sequences': 8, 'number of sites': 120 } }).formatLabel,
    'hyphy3 --legacy-json');
  const m = makeModel();
  assert.equal(m.taxaCount, 8);
  assert.equal(m.sitesCount, 120);
});
