import { test } from 'node:test';
import assert from 'node:assert/strict';
import { installFakeDOM, FakeElement } from './fake-dom.mjs';
import { ResultsTable, resolveColumns, virtualWindow, columnChooserItems } from '../src/components/results-table.js';
import { FEL_TABLE } from '../src/methods/fel.js';

installFakeDOM();

// A branch-keyed result whose columns come from the file: one per declared statistic
const column = (id, extra = {}) => Object.assign({
  id,
  header: id.toUpperCase(),
  cell: (r) => `<td>${r[id]}</td>`,
  csv: [{ header: id, value: (r) => r[id] }],
  sort: (r) => r[id]
}, extra);

function branchData(n, declared = ['lrt', 'length'], capabilities = {}) {
  const branches = Array.from({ length: n }, (_, i) => ({ branch: `B${i}`, lrt: i, length: i / 10, ebf: 2 * i }));
  return { branches, declared, capabilities, getStats: () => ({ total: n }) };
}

const dynamicColumns = (data) => [
  column('branch'),
  ...data.declared.map(id => column(id, { hiddenByDefault: id === 'length' })),
  column('ebf', { availableWhen: (caps) => caps.ebf })
];

function branchTable(options = {}) {
  return new ResultsTable(new FakeElement(), Object.assign({
    rowKey: 'branch',
    rows: (data) => data.branches,
    columns: dynamicColumns,
    filters: [{ id: 'all', key: 'all', label: 'All', count: (s) => s.total, test: () => true, subtitle: (n) => `${n}` }]
  }, options));
}

const headers = (table) => [...table.thead.innerHTML.matchAll(/<th[^>]*>([^<]*)<\/th>/g)].map(m => m[1]);

test('columns(data) resolve per result and availableWhen reads its capabilities', () => {
  const data = branchData(2, ['lrt']);
  assert.deepEqual(resolveColumns(dynamicColumns, data).map(c => c.id), ['branch', 'lrt']);
  data.capabilities = { ebf: true };
  assert.deepEqual(resolveColumns(dynamicColumns, data).map(c => c.id), ['branch', 'lrt', 'ebf']);
  const noCaps = { declared: [] };
  assert.deepEqual(resolveColumns(dynamicColumns, noCaps).map(c => c.id), ['branch']);
  assert.deepEqual(resolveColumns(FEL_TABLE.columns, {}), FEL_TABLE.columns);
});

test('the table takes its columns from each load and hides hiddenByDefault ones', () => {
  const table = branchTable();
  assert.deepEqual(table.columns, []);
  table.setData(branchData(3));
  assert.deepEqual(headers(table), ['BRANCH', 'LRT']);
  assert.equal(table.tbody.innerHTML.match(/<tr /g).length, 3);
  assert.ok(!table.tbody.innerHTML.includes('<td>0.1</td>'));

  table.setData(branchData(3, ['lrt'], { ebf: true }));
  assert.deepEqual(headers(table), ['BRANCH', 'LRT', 'EBF']);
});

test('the chooser appears for hidden columns and its choices survive a reload', () => {
  const table = branchTable();
  table.setData(branchData(3));
  assert.ok(table.chooser, 'a hiddenByDefault column brings the chooser');
  assert.equal(table.exportBtn.inserted, table.chooser);
  const menu = table.chooser.querySelector('.column-chooser-menu').innerHTML;
  assert.ok(menu.includes('data-column="length" /> LENGTH'));
  assert.ok(menu.includes('data-column="lrt" checked /> LRT'));

  table.chooser.fire('change', { target: { dataset: { column: 'length' }, checked: true } });
  assert.deepEqual(headers(table), ['BRANCH', 'LRT', 'LENGTH']);
  table.chooser.fire('change', { target: { dataset: { column: 'lrt' }, checked: false } });
  assert.deepEqual(table.shownColumns().flatMap(c => c.csv).map(f => f.header), ['branch', 'length']);

  table.setData(branchData(4));
  assert.deepEqual(headers(table), ['BRANCH', 'LENGTH']);
  assert.ok(table.chooser.querySelector('.column-chooser-menu').innerHTML.includes('data-column="length" checked'));
});

test('a result without hidden columns drops the chooser unless the config asks for it', () => {
  const table = branchTable();
  table.setData(branchData(2));
  const chooser = table.chooser;
  table.setData(branchData(2, ['lrt']));
  assert.equal(table.chooser, null);
  assert.ok(chooser.removed);

  const forced = branchTable({ columnChooser: true });
  forced.setData(branchData(2, ['lrt']));
  assert.ok(forced.chooser);
  assert.equal(columnChooserItems([column('lrt')], () => false),
    '<label class="column-chooser-item"><input type="checkbox" data-column="lrt" /> LRT</label>');
});

test('virtualWindow covers the viewport plus overscan and sizes the spacers', () => {
  assert.deepEqual(virtualWindow(0, 380, 32, 5000, 10), { first: 0, last: 22, before: 0, after: (5000 - 22) * 32 });
  assert.deepEqual(virtualWindow(3200, 380, 32, 5000, 10), { first: 90, last: 122, before: 2880, after: (5000 - 122) * 32 });
  assert.deepEqual(virtualWindow(159800, 380, 32, 5000, 10), { first: 4983, last: 5000, before: 4983 * 32, after: 0 });
});

test('a virtual table above minRows draws only its window between spacers', () => {
  const table = branchTable({ virtualize: true });
  table.tableScroll.clientHeight = 320;
  table.setData(branchData(2500));
  assert.ok(table.virtual);
  const drawn = () => table.tbody.innerHTML.match(/<tr data-branch="B(\d+)"/g).map(s => Number(s.match(/\d+/)[0]));
  assert.deepEqual(drawn(), Array.from({ length: 22 }, (_, i) => i));
  assert.ok(table.tbody.innerHTML.endsWith(`<tr class="virtual-spacer" aria-hidden="true"><td colspan="2" style="height: ${(2500 - 22) * 32}px;"></td></tr>`));

  table.tableScroll.scrollTop = 32 * 1000;
  table.tableScroll.fire('scroll');
  assert.equal(drawn()[0], 988);
  assert.ok(table.tbody.innerHTML.startsWith(`<tr class="virtual-spacer" aria-hidden="true"><td colspan="2" style="height: ${988 * 32}px;"></td></tr>`));

  table.selectSite('B2000');
  assert.ok(drawn().includes(2000));
  assert.ok(table.tbody.innerHTML.includes('<tr data-branch="B2000" class="row-selected">'));
});

test('below minRows, or without virtualize, every row is drawn', () => {
  const small = branchTable({ virtualize: { minRows: 100 } });
  small.setData(branchData(99));
  assert.equal(small.virtual, false);
  assert.equal(small.tbody.innerHTML.match(/<tr /g).length, 99);

  const plain = branchTable();
  plain.setData(branchData(2500));
  assert.equal(plain.virtual, false);
  assert.equal(plain.tbody.innerHTML.match(/<tr /g).length, 2500);
  assert.equal(FEL_TABLE.virtualize, undefined);
});

test('a table without search draws no search box and ignores a stray query', () => {
  const table = branchTable();
  assert.ok(!table.wrapper.innerHTML.includes('search-input'));
  table.setData(branchData(3));
  table.searchQuery = 'B1';
  assert.equal(table.visibleRows().length, 3);

  const searched = branchTable({ search: { placeholder: 'Find branch', match: (r, q) => r.branch.toLowerCase() === q } });
  assert.ok(searched.wrapper.innerHTML.includes('placeholder="Find branch"'));
  searched.setData(branchData(3));
  searched.searchQuery = 'B1';
  assert.deepEqual(searched.visibleRows().map(r => r.branch), ['B1']);
});

test('a filter with its own sort orders its rows until the viewer picks a column', () => {
  const filters = [
    { id: 'all', key: 'all', label: 'All', count: (s) => s.total, test: () => true, subtitle: (n) => `${n}` },
    { id: 'rev', key: 'rev', label: 'Rev', sort: { column: 'lrt', dir: 'desc' }, count: (s) => s.total, test: () => true, subtitle: (n) => `${n}` }
  ];
  const table = branchTable({ filters, defaultSort: { column: 'lrt', dir: 'asc' } });
  table.setData(branchData(3));
  const order = () => table.visibleRows().map(r => r.branch);
  assert.deepEqual(order(), ['B0', 'B1', 'B2']);
  table.setFilterType('rev');
  assert.deepEqual(order(), ['B2', 'B1', 'B0']);
  assert.equal(headers(table)[1], 'LRT');
  table.sortBy(table.displayColumns()[1]);
  assert.deepEqual(order(), ['B0', 'B1', 'B2']);
  assert.equal(headers(table)[1], 'LRT ↑');
  table.setFilterType('all');
  assert.deepEqual(order(), ['B0', 'B1', 'B2']);
});
