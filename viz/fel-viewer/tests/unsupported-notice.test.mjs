import { test } from 'node:test';
import assert from 'node:assert/strict';
import { unsupportedNotice } from '../src/components/unsupported-notice.js';

test('the notice names the method and the methods that do have a viewer', () => {
  const html = unsupportedNotice('busted', ['fel', 'meme']);
  assert.match(html, /class="hyphy-unsupported-notice"/);
  assert.match(html, /<div class="unsupported-notice-title">No viewer for BUSTED results yet<\/div>/);
  assert.match(html, /The viewer can display: FEL, MEME\./);
});

test('a method id from the file cannot inject markup', () => {
  const html = unsupportedNotice('<img src=x onerror=alert(1)>', ['fel']);
  assert.ok(!html.includes('<img'));
  assert.match(html, /No viewer for &lt;IMG SRC=X ONERROR=ALERT\(1\)&gt; results yet/);
});
