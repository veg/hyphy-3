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

// Rule blocks with their selectors, so token blocks can be told apart from component rules
const ruleBlocks = [...css.replace(/\/\*[\s\S]*?\*\//g, '').matchAll(/([^{}]+)\{([^}]*)\}/g)]
  .map(([, selector, body]) => ({ selector: selector.trim(), body }));
const darkBlocks = ruleBlocks.filter(b => b.selector.includes('hyphy-theme-dark'));

test('the dark theme is one override of the :root tokens', () => {
  assert.equal(darkBlocks.length, 1);
  const declarations = darkBlocks[0].body.split(';').map(d => d.trim()).filter(Boolean);
  for (const d of declarations) {
    const name = d.match(/^--([\w-]+):/)?.[1];
    assert.ok(name, `dark theme sets a plain property: ${d}`);
    assert.ok(name in cssTokens, `--${name} is not a :root token`);
  }
});

test('colour literals live only in the token blocks', () => {
  for (const { selector, body } of ruleBlocks) {
    if (selector === ':root' || selector.includes('hyphy-theme-dark')) continue;
    assert.doesNotMatch(body, /#[0-9A-Fa-f]{3,8}\b/, selector);
  }
  assert.doesNotMatch(css, /#3B82F6/i);
});

test('filled controls read the fill tokens, which default to the ink pair', () => {
  assert.equal(cssTokens['ink-fill'], 'var(--ink)');
  assert.equal(cssTokens['on-ink-fill'], 'var(--knockout)');
  const filled = ruleBlocks.filter(b => /background-color:\s*var\(--ink-fill\)/.test(b.body));
  assert.ok(filled.length >= 6);
  for (const { selector, body } of filled) {
    assert.doesNotMatch(body, /var\(--knockout\)/, selector);
  }
});
