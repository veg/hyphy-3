import { test } from 'node:test';
import assert from 'node:assert/strict';
import { drawerTabButtons, drawerTabPanels } from '../src/components/drawer.js';
import { FEL_DRAWER } from '../src/methods/fel.js';
import { MEME_DRAWER } from '../src/methods/meme.js';

const site = (alpha, beta, pValue, extra = {}) => Object.assign({
  site: 9, alpha, beta, pValue, lrt: 4.19, omega: beta / alpha, history: 'V → I', events: []
}, extra);
const data = { eventsBySite: {}, sequences: { a: 'AAA'.repeat(8) + 'CCC', b: 'AAA'.repeat(8) + 'CCG' } };

test('the first tab starts active and the others start hidden', () => {
  assert.equal(drawerTabButtons(FEL_DRAWER.tabs),
    '<button type="button" class="seg-btn active" data-tab="summary" id="tab-btn-summary">Overview &amp; Events</button>\n'
    + '          <button type="button" class="seg-btn" data-tab="tree" id="tab-btn-tree">Site Tree 🌳</button>');
  assert.equal(drawerTabPanels(FEL_DRAWER.tabs),
    '<div id="drawer-summary-panel">\n          <!-- Rates and events table -->\n        </div>\n'
    + '        <div id="drawer-tree-panel" style="display: none;">\n          <!-- SiteTree component mounted here -->\n        </div>');
});

test('FEL badges a site by direction at the live threshold', () => {
  const badge = (r, threshold) => FEL_DRAWER.header(r, { data, threshold });
  assert.deepEqual(badge(site(0.1, 1.5, 0.04), 0.1), { tag: 'Codon #9', tone: 'div', badge: '▲ Diversifying Selection' });
  assert.deepEqual(badge(site(1.5, 0.1, 0.04), 0.1), { tag: 'Codon #9', tone: 'pur', badge: '▼ Purifying Selection' });
  assert.equal(badge(site(0.1, 1.5, 0.04), 0.01).tone, 'neutral');
});

test('MEME badges significant sites as episodic and never as purifying', () => {
  const badge = (r) => MEME_DRAWER.header(r, { data, threshold: 0.1 });
  assert.equal(badge(site(1.5, 0.1, 0.04)).badge, '▲ Episodic Diversifying Selection');
  assert.equal(badge(site(1.5, 0.1, 0.4)).badge, 'Neutral / Constrained');
});

test('the overview tab draws the method rates grid above the shared history', () => {
  const render = (drawer, r) => {
    const panel = {};
    drawer.tabs[0].render(panel, r, { data, threshold: 0.1 });
    return panel.innerHTML;
  };
  const fel = render(FEL_DRAWER, site(0.1, 1.5, 0.04));
  assert.ok(fel.includes('Maximum Likelihood Rates</div>'));
  assert.ok(fel.includes('<span class="extant-codon font-mono">CCC</span>'));
  assert.ok(!fel.includes('Branches Under Episodic Selection</div>'));
  const meme = render(MEME_DRAWER, site(0.1, 1.5, 0.04, { branchEBF: { Pig: 250 }, p_plus: 0.2 }));
  assert.ok(meme.includes('Maximum Likelihood Rates (MEME Mixture Model)'));
  assert.ok(meme.includes('<td class="font-mono font-bold">Pig</td>'));
});
