import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readFileSync } from 'node:fs';
import { FELData } from '../src/data/fel-data.js';
import { FEL_CONFIG } from '../src/methods/fel.js';
import { MEME_CONFIG } from '../src/methods/meme.js';
import { ResultsTable } from '../src/components/results-table.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';
const fixture = name => JSON.parse(readFileSync(new URL(name, FIXTURES), 'utf8'));

test('FEL and MEME reports name their output format in the header chip', { skip: noFixtures }, () => {
  const cases = [
    ['fel_cd2.json', 'v3'], ['fel_cd2.legacy.json', 'hyphy3-legacy'], ['fel_cd2.hyphy25.json', 'hyphy25'],
    ['meme_cd2.json', 'v3'], ['meme_cd2.legacy.json', 'hyphy3-legacy'], ['meme_cd2.hyphy25.json', 'hyphy25']
  ];
  for (const [name, format] of cases) {
    const data = new FELData(fixture(name));
    assert.equal(data.format, format, name);
    assert.ok(FEL_CONFIG.formatLabel(data), name);
    assert.equal(MEME_CONFIG.formatLabel(data), data.formatLabel, name);
  }
});

test('the drawer steps through the closest calls when the filter leaves the table empty', () => {
  const rows = [{ site: 1 }, { site: 2 }];
  const table = { visibleRows: () => [], fallbackRows: [{ site: 7 }, { site: 3 }] };
  assert.deepEqual(ResultsTable.prototype.stepRows.call(table).map(r => r.site), [7, 3]);
  table.visibleRows = () => rows;
  assert.deepEqual(ResultsTable.prototype.stepRows.call(table), rows);
});
