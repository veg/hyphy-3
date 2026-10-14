import { test } from 'node:test';
import assert from 'node:assert/strict';
import { existsSync, readdirSync, readFileSync } from 'node:fs';
import { METHOD_CONFIGS, methodConfig, configForResult, resolveConfig } from '../src/methods/registry.js';

const FIXTURES = new URL('../../../parity/viewer/fixtures/', import.meta.url);
const noFixtures = !existsSync(FIXTURES) && 'parity fixtures are not in this checkout';

test('every config supplies what the shell mounts', () => {
  for (const [id, config] of Object.entries(METHOD_CONFIGS)) {
    assert.equal(config.id, id);
    assert.equal(typeof config.Parser, 'function');
    assert.equal(typeof config.defaultFilter, 'string');
    for (const part of ['kpis', 'hero', 'table', 'drawer', 'diagnostics', 'row']) {
      assert.equal(typeof config[part], 'function', `${id}.${part}`);
    }
    assert.equal(typeof config.controller.seed, 'number', `${id}.controller.seed`);
    // p-values and q-values pass at or below the cutoff, posteriors at or above it
    const { statistic, direction } = config.controller;
    assert.ok(direction === 'ge' ? statistic === 'Posterior' : ['p', 'q'].includes(statistic), `${id}.controller.statistic`);
  }
});

test('configs are looked up by method id, case-insensitively', () => {
  assert.equal(methodConfig('fel'), METHOD_CONFIGS.fel);
  assert.equal(methodConfig('MEME'), METHOD_CONFIGS.meme);
  assert.equal(methodConfig('busted'), null);
  assert.equal(methodConfig(null), null);
  assert.equal(methodConfig('constructor'), null);
});

test('FEL, MEME, FUBAR and SLAC results get their own config, v3 or legacy', { skip: noFixtures }, () => {
  const files = readdirSync(FIXTURES).filter(f => /^(fel|meme|fubar|slac)_.*\.json$/.test(f));
  assert.ok(files.length > 0);
  for (const file of files) {
    const json = JSON.parse(readFileSync(new URL(file, FIXTURES), 'utf8'));
    assert.equal(configForResult(json)?.id, file.slice(0, file.indexOf('_')), file);
  }
});

test('results of methods without a viewer get no config', () => {
  assert.equal(configForResult({ analysis: { id: 'busted' } }), null);
  assert.equal(configForResult({}), null);
});

test('only unrecognized results fall back; a recognized method without a viewer does not', () => {
  const fallback = METHOD_CONFIGS.meme;
  assert.deepEqual(resolveConfig({ analysis: { id: 'busted' } }, fallback), { method: 'busted', config: null });
  assert.deepEqual(resolveConfig({ analysis: { info: 'Directional FEL detects site-level selection' } }, fallback),
    { method: 'dfel', config: null });
  assert.deepEqual(resolveConfig({ breakpointData: {} }, fallback), { method: 'gard', config: null });
  assert.deepEqual(resolveConfig({ MLE: { content: { 0: { 'by-site': {} } } } }, fallback), { method: 'slac', config: METHOD_CONFIGS.slac });
  assert.deepEqual(resolveConfig({ grid: [], posterior: {} }, fallback), { method: 'fubar', config: METHOD_CONFIGS.fubar });
  assert.deepEqual(resolveConfig({ analysis: { id: 'fel' } }, fallback), { method: 'fel', config: METHOD_CONFIGS.fel });
  assert.deepEqual(resolveConfig({ something: 'else' }, fallback), { method: null, config: fallback });
});

test('diagnostics name the method and list the MG94 fit', () => {
  const data = {
    analysisType: 'fel', numSites: 10, numTaxa: 4, globalOmega: 0.5,
    raw: { model_fits: { codon_mg94: { log_likelihood: -123.456, aicc: 260.1 } } }
  };
  const fel = METHOD_CONFIGS.fel.diagnostics(data);
  assert.match(fel, /FEL \(Fixed Effects Likelihood\)/);
  assert.match(fel, /-123\.46/);
  assert.match(fel, /260\.10/);
  const meme = METHOD_CONFIGS.meme.diagnostics(data);
  assert.match(meme, /MEME \(Mixed Effects Model of Episodic Selection\)/);
  assert.doesNotMatch(meme, /Global Mean dN\/dS/);
});
