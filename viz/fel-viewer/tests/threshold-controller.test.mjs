import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  KpiBar, sliderPosition, sliderThreshold, sliderAttributes, cutoffAttribute, presetActive,
  thresholdPasses, thresholdRelation, thresholdCaption, thresholdLabel, thresholdControllerMarkup,
  cardsFitBlank
} from '../src/components/kpi-bar.js';
import { SITE_P_THRESHOLD } from '../src/methods/fel.js';
import { formatThreshold } from '../src/core/format.js';

const POSTERIOR = {
  title: 'Posterior Cutoff', statistic: 'PP', direction: 'ge',
  min: 0.5, max: 0.999, step: 0.001, presets: [0.9, 0.95, 0.99], seed: 0.9, format: formatThreshold
};
const LOG_P = {
  title: 'Significance Cutoff', statistic: 'p', direction: 'le', scale: 'log',
  min: 0.001, max: 1, step: 0.01, presets: [0.1, 0.05, 0.01], seed: 0.05, format: formatThreshold,
  caption: ({ statistic }) => `${statistic === 'q' ? 'FDR-corrected' : 'Uncorrected'} site test`,
  toggle: { statistics: ['p', 'q'] }
};

test('FEL slider attributes keep their label form when it reads back exactly', () => {
  assert.deepEqual(sliderAttributes(SITE_P_THRESHOLD, 0.1), { min: '0.001', max: '0.20', step: '0.005', value: '0.10' });
});

test('slider attributes carry the raw cutoff where the label would round it', () => {
  assert.equal(cutoffAttribute(SITE_P_THRESHOLD, 0.0501), '0.0501');
  assert.deepEqual(sliderAttributes(POSTERIOR, 0.995), { min: '0.50', max: '0.999', step: '0.001', value: '0.995' });
  assert.equal(cutoffAttribute(POSTERIOR, 0.999), '0.999');
});

test('a log scale maps slider position x to the cutoff 10^x', () => {
  assert.deepEqual(sliderAttributes(LOG_P, 0.01), { min: '-3', max: '0', step: '0.01', value: '-2' });
  assert.equal(sliderThreshold(LOG_P, -2), 0.01);
  assert.equal(sliderThreshold(LOG_P, -1.3), 0.0501);
  assert.equal(sliderPosition(LOG_P, 0.001), -3);
  assert.equal(sliderThreshold(SITE_P_THRESHOLD, 0.051), 0.051);
  assert.equal(sliderPosition(SITE_P_THRESHOLD, 0.051), 0.051);
});

test('presets match the cutoff in slider positions, so small log cutoffs stay distinct', () => {
  assert.ok(presetActive(LOG_P, 0.05, 0.0501));
  assert.ok(!presetActive(LOG_P, 0.001, 0.002));
  assert.ok(presetActive(SITE_P_THRESHOLD, 0.05, 0.051));
  assert.ok(!presetActive(SITE_P_THRESHOLD, 0.05, 0.046));
});

test('direction sets the comparison, the operator and the wording', () => {
  assert.ok(thresholdPasses(POSTERIOR, 0.95, 0.9));
  assert.ok(!thresholdPasses(POSTERIOR, 0.85, 0.9));
  assert.ok(thresholdPasses(SITE_P_THRESHOLD, 0.05, 0.05));
  assert.ok(!thresholdPasses(SITE_P_THRESHOLD, 0.2, 0.1));
  assert.equal(thresholdLabel(POSTERIOR, 0.95), 'PP ≥ 0.95');
  assert.equal(thresholdRelation(POSTERIOR), 'at least');
  assert.equal(thresholdRelation(SITE_P_THRESHOLD), 'at most');
});

test('a caption is a fixed line or follows the statistic and cutoff', () => {
  assert.equal(thresholdCaption(SITE_P_THRESHOLD, 0.1, 'p'), null);
  assert.equal(thresholdCaption({ caption: 'Does not change the c-AIC verdict' }, 0.1, 'p'), 'Does not change the c-AIC verdict');
  assert.equal(thresholdCaption(LOG_P, 0.05, 'q'), 'FDR-corrected site test');
});

test('the controller draws a caption line and a statistic toggle only when configured', () => {
  const plain = thresholdControllerMarkup(SITE_P_THRESHOLD, 0.1, 'p', 4);
  assert.ok(!plain.includes('kpi-control-caption'));
  assert.ok(!plain.includes('kpi-stat-toggle'));

  const full = thresholdControllerMarkup(LOG_P, 0.05, 'q', 2);
  assert.ok(full.includes('id="kpi-control-caption"'));
  assert.ok(full.includes('<button type="button" class="seg-btn" data-statistic="p">p</button>'));
  assert.ok(full.includes('<button type="button" class="seg-btn active" data-statistic="q">q</button>'));
  assert.ok(full.includes('data-p="0.05">q ≤ 0.05</button>'));
  assert.ok(full.includes(`value="${Math.log10(0.05)}"`));
});

test('posterior presets are drawn at their own cutoffs with ≥ labels', () => {
  const markup = thresholdControllerMarkup(POSTERIOR, 0.95, 'PP', 0);
  assert.ok(markup.includes('<button type="button" class="preset-pill active" data-p="0.95">PP ≥ 0.95</button>'));
  assert.ok(markup.includes('max="0.999"'));
});

test('cards fit the blank cards only with the same ids in the same order', () => {
  const ids = (...list) => list.map(id => ({ id }));
  assert.ok(cardsFitBlank(ids('a', 'b'), ids('a', 'b')));
  assert.ok(!cardsFitBlank(ids('a', 'b', 'c'), ids('a', 'b')));
  assert.ok(!cardsFitBlank(ids('b', 'a'), ids('a', 'b')));
});

// Enough of a card element for fillCards to relabel, remove or extend it
function fakeCard(id, widget = null) {
  const card = { id, removed: false, children: [], label: { textContent: '' }, value: { className: '' }, widget };
  card.root = { remove() { card.removed = true; }, appendChild(el) { card.children.push(el); } };
  return card;
}

function fakeBar(cards, blank) {
  const bar = Object.create(KpiBar.prototype);
  bar.cards = cards;
  bar.blank = blank;
  bar.inserted = [];
  bar.controlEl = { insertAdjacentHTML: (where, html) => bar.inserted.push([where, html]) };
  bar.cardEls = blank.map(card => fakeCard(card.id));
  bar.bindCards = list => list.map(card => fakeCard(card.id, card.widget ? { slot: card.id } : null));
  return bar;
}

test('cards whose ids differ from the blank cards replace them with their own markup', () => {
  const blank = [{ id: 'align', label: 'Alignment' }];
  const cards = [{ id: 'model', label: 'Model' }, { id: 'tested', label: 'Tested', widget() {} }];
  const bar = fakeBar(cards, blank);
  const old = bar.cardEls;
  bar.fillCards();
  assert.ok(old.every(c => c.removed));
  assert.equal(bar.inserted.length, 1);
  assert.equal(bar.inserted[0][0], 'beforebegin');
  assert.ok(bar.inserted[0][1].includes('id="kpi-model"'));
  assert.ok(bar.inserted[0][1].includes('id="kpi-tested"'));
  assert.ok(bar.inserted[0][1].includes('<div class="kpi-widget"></div>'));
  assert.equal(bar.cardEls.length, 2);
  assert.deepEqual(bar.cardEls.map(c => c.label.textContent), ['Model', 'Tested']);
  assert.equal(bar.blank, cards);
});

test('a card the result cannot fill is left out, and the rest replace the blank cards', () => {
  const cards = [{ id: 'align', label: 'Alignment' }, { id: 'model', label: 'Model', available: (data) => data.model }];
  const bar = fakeBar(cards, cards);
  bar.options = { cards };
  bar.updateStats = () => {};
  bar.setData({ model: false });
  assert.deepEqual(bar.cards.map(c => c.id), ['align']);
  assert.equal(bar.inserted.length, 1);
  assert.doesNotMatch(bar.inserted[0][1], /kpi-model/);
  bar.setData({ model: true });
  assert.deepEqual(bar.cards.map(c => c.id), ['align', 'model']);
});

test('a card with a widget gets a slot even when its blank twin had none', () => {
  globalThis.document = { createElement: () => ({ className: '' }) };
  const blank = [{ id: 'omega', label: 'ω' }];
  const cards = [{ id: 'omega', label: 'Mean ω', tone: 'accent', widget() {} }];
  const bar = fakeBar(cards, blank);
  bar.fillCards();
  assert.equal(bar.inserted.length, 0);
  const [els] = bar.cardEls;
  assert.equal(els.label.textContent, 'Mean ω');
  assert.equal(els.value.className, 'kpi-value kpi-text-orange');
  assert.equal(els.widget.className, 'kpi-widget');
  assert.deepEqual(els.children, [els.widget]);
});

test('a cutoff set from outside moves the slider to its position', () => {
  const bar = Object.create(KpiBar.prototype);
  Object.assign(bar, {
    controller: LOG_P, statistic: 'p', options: {}, data: null, caption: null,
    slider: { value: '-1' }, display: { textContent: '' }, presetBtns: []
  });
  bar.setThreshold(0.01);
  assert.equal(bar.slider.value, -2);
  assert.equal(bar.display.textContent, 'p ≤ 0.01');
});
