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

test('legacy info names the method', () => {
  const method = info => detectMethod({ analysis: { info } }).method;
  assert.equal(method('FEL (Fixed Effects Likelihood)\n estimates'), 'fel');
  assert.equal(method('aBSREL (Adaptive branch-site random effects likelihood)'), 'absrel');
  assert.equal(method('Contrast-FEL (Fixed Effects Likelihood) investigates'), 'contrast-fel');
  assert.equal(method('GARD : Genetic Algorithms for Recombination'), 'gard');
  assert.equal(method('Perform a Fast Unbiased AppRoximate Bayesian (FUBAR) analysis'), 'fubar');
  assert.equal(method('Directional FEL detects site-level selection'), 'dfel');
  assert.equal(method('Something else entirely'), null);
});

test('layout signatures that share columns resolve in order', () => {
  const headers = names => ({ MLE: { headers: names.map(n => [n, '']) } });
  assert.equal(detectMethod(headers(['alpha', 'beta (Group1)', 'P-value (overall)'])).method, 'contrast-fel');
  assert.equal(detectMethod(headers(['alpha', 'beta', 'alpha=beta'])).method, 'fel');
  assert.equal(detectMethod(headers(['&alpha;', '&beta;<sup>+</sup>'])).method, 'meme');
  assert.equal(detectMethod({ MLE: { headers: { '&beta;<sup>+</sup>': '' } } }).method, 'meme');
  assert.equal(detectMethod({ MLE: { headers: 'bogus', content: { 0: { 'by-site': {} } } } }).method, 'slac');
  assert.equal(detectMethod({ grid: [], posterior: {}, ...headers(['alpha=beta']) }).method, 'fubar');
  assert.equal(detectMethod({ 'test results': { 'p-value': 0.1, 'relaxation or intensification parameter': 1 } }).method, 'relax');
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
