/**
 * Identify the analysis and output flavour behind a result JSON. A port of
 * detect_analysis in include/hyphy/report/html_report.hpp: the v3 id first,
 * then the method name that leads legacy analysis.info. Legacy files without
 * a usable info fall back to structural signatures, ordered so that no method
 * is taken for one whose layout it shares.
 */

const LEGACY_INFO_METHODS = ['fel', 'meme', 'slac', 'fubar', 'busted', 'absrel',
  'relax', 'prime', 'gard', 'contrast-fel', 'bgm', 'fade'];

// Returns { method, format }: method is a lowercase id ('fel', 'contrast-fel', ...)
// or null when unknown; format is 'v3' (hyphy3 native output), 'hyphy3-legacy'
// (hyphy3 --legacy-json) or 'hyphy25' (HyPhy 2.5)
export function detectMethod(json) {
  const analysis = isPlainObject(json?.analysis) ? json.analysis : {};
  if (typeof analysis.id === 'string') {
    return { method: analysis.id.toLowerCase(), format: 'v3' };
  }
  return {
    method: methodFromInfo(analysis.info) ?? methodFromLayout(json ?? {}),
    format: legacyFormat(json ?? {}, analysis)
  };
}

function isPlainObject(x) {
  return x !== null && typeof x === 'object' && !Array.isArray(x);
}

// Legacy info text begins with the method name ('FEL (Fixed Effects ...', 'GARD : ...')
function methodFromInfo(info) {
  if (typeof info !== 'string') return null;
  const text = info.toLowerCase().trimStart();
  const lead = text.match(/^[a-z0-9-]*/)[0];
  if (LEGACY_INFO_METHODS.includes(lead)) return lead;
  if (text.includes('(fubar')) return 'fubar';
  if (text.startsWith('directional fel')) return 'dfel';
  return null;
}

// Site-table column names, from either the [[name, description], ...] array
// or the { name: description } object that hyphy3 MEME writes
function mleHeaderNames(json) {
  const headers = json.MLE?.headers;
  if (Array.isArray(headers)) {
    return headers.map(h => (Array.isArray(h) ? h[0] : h)).filter(h => typeof h === 'string');
  }
  return isPlainObject(headers) ? Object.keys(headers) : [];
}

function methodFromLayout(json) {
  const headers = mleHeaderNames(json);
  const tests = isPlainObject(json['test results']) ? json['test results'] : {};
  const fits = isPlainObject(json.fits) ? json.fits : {};
  const firstPartition = json.MLE?.content?.['0'];

  if (isPlainObject(firstPartition) && 'by-site' in firstPartition) return 'slac';
  if ('grid' in json && 'posterior' in json) return 'fubar';
  if (headers.includes('P-value (overall)')) return 'contrast-fel';
  if (headers.includes('PRIME LogL')) return 'prime';
  if ('Unconstrained model' in fits && 'p-value' in tests) return 'busted';
  if ('positive test results' in tests) return 'absrel';
  if ('relaxation or intensification parameter' in tests) return 'relax';
  if (isPlainObject(json.breakpointData)) return 'gard';

  const names = headers.map(h => h.toLowerCase());
  if (names.some(h => h.includes('beta') && (h.includes('+') || h.includes('plus')))) return 'meme';
  if (names.includes('alpha=beta')) return 'fel';
  return null;
}

// HyPhy 2.5 stamps its own version into runtime ('2.5.101'); hyphy3 writes
// '3.0.0' or elapsed seconds there, or tags analysis.version ('3.0-next',
// '4.7-next', '3.0.0'). GARD carries neither, but HyPhy 2.5 keeps the line
// breaks of its multi-line info text where hyphy3 writes a single line.
function legacyFormat(json, analysis) {
  const runtime = json.runtime;
  if (typeof runtime === 'string' && runtime.startsWith('2.')) return 'hyphy25';
  if (runtime != null) return 'hyphy3-legacy';
  const version = String(analysis.version ?? '');
  if (version.includes('next') || version.startsWith('3.')) return 'hyphy3-legacy';
  if (typeof analysis.info === 'string' && !analysis.info.includes('\n')) return 'hyphy3-legacy';
  return 'hyphy25';
}
