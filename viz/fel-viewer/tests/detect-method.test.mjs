import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { detectMethod } from '../src/core/detect-method.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';

// Fixture names encode the answer: <method>_<dataset>[.hyphy25|.legacy].json
function expected(file) {
  const method = file.slice(0, file.indexOf('_'));
  if (file.endsWith('.hyphy25.json')) return { method, format: 'hyphy25' };
  if (file.endsWith('.legacy.json')) return { method, format: 'hyphy3-legacy' };
  return { method, format: 'v3' };
}

test('every parity fixture is detected', { skip: noFixtures }, () => {
  const files = readdirSync(FIXTURES).filter(f => f.endsWith('.json'));
  assert.ok(files.length > 0);
  for (const file of files) {
    const json = JSON.parse(readFileSync(new URL(file, FIXTURES), 'utf8'));
    assert.deepEqual(detectMethod(json), expected(file), file);
  }
});

test('legacy fixtures are recognized from their layout alone', { skip: noFixtures }, () => {
  const files = readdirSync(FIXTURES).filter(f => f.endsWith('.hyphy25.json') || f.endsWith('.legacy.json'));
  for (const file of files) {
    const json = JSON.parse(readFileSync(new URL(file, FIXTURES), 'utf8'));
    delete json.analysis.info;
    assert.equal(detectMethod(json).method, expected(file).method, file);
  }
});

test('v3 id wins and is lowercased', () => {
  assert.deepEqual(detectMethod({ analysis: { id: 'MEME', info: 'FEL (Fixed Effects Likelihood)' } }),
    { method: 'meme', format: 'v3' });
});

// tests/data/detect_cases.json is shared with the C++ detect_analysis test, so
// both detectors are held to the same answers
test('detection agrees with the shared test vectors', () => {
  const { cases } = JSON.parse(readFileSync(new URL('../../../tests/data/detect_cases.json', import.meta.url), 'utf8'));
  assert.ok(cases.length > 0);
  for (const c of cases) {
    assert.equal(detectMethod(c.json).method, c.method, c.name);
  }
});

test('legacy flavour comes from runtime, version or info line breaks', () => {
  const format = json => detectMethod(json).format;
  assert.equal(format({ runtime: '2.5.101', analysis: { info: 'FEL' } }), 'hyphy25');
  assert.equal(format({ runtime: '3.0.0', analysis: { info: 'SLAC' } }), 'hyphy3-legacy');
  assert.equal(format({ runtime: 12.5, analysis: { info: 'BUSTED' } }), 'hyphy3-legacy');
  assert.equal(format({ analysis: { info: 'FEL', version: '3.0-next' } }), 'hyphy3-legacy');
  assert.equal(format({ analysis: { info: 'GARD : one line', version: '0.2' } }), 'hyphy3-legacy');
  assert.equal(format({ analysis: { info: 'GARD : wrapped\n  text', version: '0.2' } }), 'hyphy25');
  assert.deepEqual(detectMethod(null), { method: null, format: 'hyphy25' });
  assert.deepEqual(detectMethod({ something: 'else' }), { method: null, format: 'hyphy25' });
});
