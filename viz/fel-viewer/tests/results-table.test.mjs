import { test } from 'node:test';
import assert from 'node:assert/strict';
import { closestCalls, resultsHeaderMarkup, resultsRowMarkup } from '../src/components/results-table.js';
import { FEL_COLUMNS, FEL_FILTERS, FEL_TABLE } from '../src/methods/fel.js';
import { MEME_COLUMNS, MEME_FILTERS, MEME_TABLE } from '../src/methods/meme.js';

const site = (n, alpha, beta, pValue, extra = {}) =>
  Object.assign({ site: n, alpha, beta, pValue, lrt: 4.2, omega: beta / alpha, history: 'I ↔ V' }, extra);
const rows = [site(1, 0.5, 2, 0.02), site(2, 2, 0.5, 0.01), site(3, 1, 1.2, 0.4), site(4, 1, 0.8, 0.3)];
const stats = { diversifying: 0, purifying: 1, significant: 1, totalSites: 4 };

test('headers carry their width and alignment, and skip export-only columns', () => {
  const header = resultsHeaderMarkup(FEL_COLUMNS, ' '.repeat(8));
  assert.ok(header.startsWith('\n        <tr>\n          <th style="width: 70px;">Codon</th>'));
  assert.ok(header.includes('\n          <th>Substitutions Reconstructed</th>'));
  assert.ok(header.endsWith('<th style="width: 140px; text-align: right;">Inspect / Tree</th>\n        </tr>\n      '));
  assert.equal(header.match(/<th/g).length, 8);
});

test('a row is keyed by rowKey and renders one cell per displayed column', () => {
  const html = resultsRowMarkup(rows[0], FEL_COLUMNS, { threshold: 0.1 }, 'site', true);
  assert.ok(html.startsWith('\n          <tr data-site="1" class="row-selected">\n            <td class="font-mono font-bold site-cell">'));
  assert.ok(html.includes('<td class="font-mono pval-cell pval-sig font-bold">0.0200</td>'));
  assert.ok(html.includes('<span class="omega-badge omega-div">4.00</span>'));
  assert.equal(html.match(/<td/g).length, 8);
});

test('FEL filters split significant sites by direction at the live threshold', () => {
  const pass = (f, t) => rows.filter(r => f.test(r, t)).map(r => r.site);
  const [div, pur, sig, all] = FEL_FILTERS;
  assert.deepEqual(pass(div, 0.05), [1]);
  assert.deepEqual(pass(pur, 0.05), [2]);
  assert.deepEqual(pass(sig, 0.05), [1, 2]);
  assert.deepEqual(pass(all, 0.05), [1, 2, 3, 4]);
  assert.deepEqual(pass(div, 0.015), []);
  assert.equal(div.subtitle(3, 0.05), 'Showing 3 diversifying sites (β > α, p ≤ 0.05)');
});

test('MEME keeps FEL\'s pills as blank slots and tests one direction', () => {
  assert.deepEqual(MEME_FILTERS.map(f => f.id), ['diversifying', 'all']);
  assert.deepEqual(MEME_TABLE.blank.filters, FEL_FILTERS);
  assert.deepEqual(rows.filter(r => MEME_FILTERS[0].test(r, 0.05)).map(r => r.site), [1, 2]);
  assert.equal(MEME_FILTERS[0].subtitle(2, 0.1), 'Showing 2 episodically diversifying sites (p ≤ 0.10)');
});

test('CSV fields export full precision and the live regime', () => {
  const fields = FEL_COLUMNS.flatMap(c => c.csv || []);
  assert.deepEqual(fields.map(f => f.header),
    ['Codon', 'Alpha (dS)', 'Beta (dN)', 'dN/dS (Omega)', 'LRT', 'p-value', 'Regime', 'Substitutions']);
  const r = site(9, 1 / 3, Infinity, 0.0123456789);
  assert.deepEqual(fields.map(f => f.value(r, { threshold: 0.05 })),
    [9, 1 / 3, Infinity, 'Infinity', 4.2, 0.0123456789, 'Diversifying', '"I ↔ V"']);

  const memeFields = MEME_COLUMNS.flatMap(c => c.csv || []);
  const m = site(5, 1, 0, 0.2, { beta_plus: 3, p_plus: 0.25, branchEBF: { a: 150, b: 20 } });
  assert.deepEqual(memeFields.map(f => [f.header, f.value(m, { threshold: 0.1 })]).slice(1, 5),
    [['Alpha (dS)', 1], ['Beta+ (dN+)', 3], ['Weight p+', 0.25], ['Branches (EBF>=100)', 1]]);
});

test('an empty FEL filter jumps to the other direction, or lists the closest calls', () => {
  const [div, pur] = FEL_FILTERS;
  const jump = FEL_TABLE.emptyState(0.05, { filter: div, stats, rows });
  assert.equal(jump.title, 'No diversifying site reaches p ≤ 0.05');
  assert.deepEqual(jump.jump, { filter: 'purifying', label: 'View Purifying Sites (▼ 1)' });

  const fallback = FEL_TABLE.emptyState(0.005, { filter: pur, stats: { ...stats, diversifying: 0 }, rows });
  assert.equal(fallback.jump, undefined);
  assert.deepEqual(fallback.fallback.map(r => r.site), [2, 4]);
});

test('MEME\'s empty state offers the closest calls by p', () => {
  const state = MEME_TABLE.emptyState(0.001, { filter: MEME_FILTERS[0], stats, rows });
  assert.equal(state.title, 'No episodically diversifying site reaches p ≤ 0.001');
  assert.deepEqual(state.fallback.map(r => r.site), [2, 1, 4, 3]);
  assert.deepEqual(closestCalls(rows, r => r.pValue, 2).map(r => r.site), [2, 1]);
});
