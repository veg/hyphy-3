// The :root tokens in viewer.css and PALETTE in palette.js must carry identical values
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { PALETTE, withAlpha } from '../src/data/palette.js';

const css = readFileSync(new URL('../src/styles/viewer.css', import.meta.url), 'utf8');
const rootBlock = css.match(/:root\s*\{([^}]*)\}/)[1];
const cssTokens = Object.fromEntries(
  [...rootBlock.matchAll(/--([\w-]+):\s*([^;]+);/g)].map(([, name, value]) => [name, value.trim()])
);
const camel = name => name.replace(/-(\w)/g, (_, c) => c.toUpperCase());

test('every solid CSS token matches PALETTE', () => {
  for (const [name, value] of Object.entries(cssTokens)) {
    if (!value.startsWith('#')) continue;
    assert.equal(PALETTE[camel(name)], value, `--${name}`);
  }
});

test('every CSS tint is a PALETTE colour at some alpha', () => {
  const tints = { 'accent-tint': 'accent', 'accent-tint-border': 'accent', 'purifying-tint': 'purifying',
                  'synonymous-tint': 'synonymous', 'tooltip-bg': 'ink' };
  for (const [name, base] of Object.entries(tints)) {
    const alpha = cssTokens[name].match(/([\d.]+)\)$/)[1];
    assert.equal(withAlpha(PALETTE[base], alpha), cssTokens[name], `--${name}`);
  }
});

test('withAlpha formats rgba() like the hand-written literals it replaced', () => {
  assert.equal(withAlpha('#EB5E28', 0.15), 'rgba(235, 94, 40, 0.15)');
  assert.equal(withAlpha('#1E40AF', 0.12), 'rgba(30, 64, 175, 0.12)');
});
