// Method registry: one config per analysis the viewer can draw. A config supplies
// the parser, the KPI bar, the hero, the results table, the drawer and the
// diagnostics summary; ViewerShell lays them out and keeps them in sync.

import { detectMethod } from '../core/detect-method.js';
import { FEL_CONFIG } from './fel.js';
import { MEME_CONFIG } from './meme.js';
import { FUBAR_CONFIG } from './fubar.js';

export const METHOD_CONFIGS = {
  fel: FEL_CONFIG,
  meme: MEME_CONFIG,
  fubar: FUBAR_CONFIG
};

// Config for a method id ('fel', 'MEME', ...), or null when the viewer has none
export function methodConfig(id) {
  const key = String(id ?? '').toLowerCase();
  return Object.hasOwn(METHOD_CONFIGS, key) ? METHOD_CONFIGS[key] : null;
}

// Config for a result JSON, v3 or legacy, or null when its method has no viewer
export function configForResult(json) {
  return methodConfig(detectMethod(json).method);
}

// What the shell mounts for a result JSON: { method, config }. Only a result whose
// method goes unrecognized gets the fallback; a recognized method without a viewer
// gets a null config, so it is never drawn with another method's parser.
export function resolveConfig(json, fallback) {
  const method = detectMethod(json).method;
  return { method, config: method == null ? fallback : methodConfig(method) };
}
