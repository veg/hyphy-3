// Directional FEL results: per-site alpha/beta fits under a model that lets base
// composition drift along lineages (K regimes: AT-drift, Neutral, GC-drift) on a
// likelihood-chosen root, each site tested by a 1-df LRT of alpha = beta. hyphy3 v3 output
// is read from its schema paths (site_results, dataset, compositional_heterogeneity,
// directional_model, root_selection, phylogeny, model_fits). hyphy3 --legacy-json writes
// a FEL-shaped site table (MLE.headers / MLE.content['0'], one row per site) and nothing
// else of the directional model; it maps onto the same site rows, and the composition
// screen, rooting, regimes and reconstructed substitutions stay absent.

import { detectMethod, FORMAT_LABELS } from '../core/detect-method.js';

// The optimiser's rate bounds: a rate at the floor is 'zero', one at the ceiling is
// pinned; floating-point noise above the floor (1.0000000000000186e-4) still counts
export const DFEL_RATE_FLOOR = 1e-4;
export const DFEL_RATE_CEILING = 999.5;

export function dfelAtFloor(rate) {
  return rate <= DFEL_RATE_FLOOR * (1 + 1e-9);
}

// Legacy site columns by header name, with the writer's position as the fallback
const DFEL_LEGACY_COLUMNS = {
  alpha: ['alpha', 0],
  beta: ['beta', 1],
  alphaNull: ['alpha=beta', 2],
  lrt: ['LRT', 3],
  pValue: ['p-value', 4]
};

// Regime ids as the writer numbers them; weights[] is indexed by the same id
export const DFEL_REGIMES = [
  { id: 0, name: 'AT-drift', short: 'AT' },
  { id: 1, name: 'Neutral', short: 'Neutral' },
  { id: 2, name: 'GC-drift', short: 'GC' }
];

const DFEL_TOKEN = /^([A-Z*]) → ([A-Z*])(?: \[(GC|AT)↑\])?(⚡)?$/;
const DFEL_SYN_SUFFIX = /^(.*?)\s*\(\+(\d+) syn\)$/;

// One site's reconstructed-substitution string: 'Conserved', 'K → E [GC↑], T → R⚡ (+1 syn)'
// or '(+2 syn)' alone. Returns { conserved, synOnly, events: [{ from, to, dir, surprising }],
// synCount, unparsed }. Repeated tokens are separate branch events and stay separate.
// 'Conserved (+syn)' is unreachable in the writer and reads as synonymous-only.
export function parseDfelSubstitutions(text) {
  const s = String(text ?? '').trim();
  const out = { conserved: false, synOnly: false, events: [], synCount: 0, unparsed: [] };
  if (s === 'Conserved') {
    out.conserved = true;
    return out;
  }
  let body = s;
  const syn = s.match(DFEL_SYN_SUFFIX);
  if (syn) {
    body = syn[1].trim();
    out.synCount = Number(syn[2]);
  } else if (/^Conserved\s*\(\+syn\)$/.test(s)) {
    body = '';
  }
  if (body === 'Conserved') body = '';
  for (const token of body ? body.split(', ') : []) {
    const m = token.trim().match(DFEL_TOKEN);
    if (m) out.events.push({ from: m[1], to: m[2], dir: m[3] ?? null, surprising: Boolean(m[4]) });
    else out.unparsed.push(token);
  }
  out.synOnly = out.events.length === 0 && out.unparsed.length === 0;
  return out;
}

// Counts the strip, table, drawer and CSV read: GC-ward, AT-ward, untagged, ⚡ and synonymous
export function dfelChangeTally(subs) {
  const tally = { nonsyn: subs.events.length, gc: 0, at: 0, flat: 0, surprising: 0, syn: subs.synCount };
  for (const e of subs.events) {
    if (e.dir === 'GC') tally.gc++;
    else if (e.dir === 'AT') tally.at++;
    else tally.flat++;
    if (e.surprising) tally.surprising++;
  }
  return tally;
}

// A site's call at cutoff t, recomputed from p, alpha and beta. A site with nothing to fit
// (alpha = beta = 0) is the 'invariant' display state of neutral.
export function dfelSiteClass(r, t) {
  if (r.pValue <= t && r.beta > r.alpha) return 'diversifying';
  if (r.pValue <= t && r.alpha > r.beta) return 'purifying';
  return r.invariant ? 'invariant' : 'neutral';
}

// omega = beta / alpha with the bounds spelled out: Infinity when alpha sits at the floor
// and beta does not, 0 when beta does and alpha does not, null when both do (nothing to
// compare)
export function dfelOmega(alpha, beta) {
  if (dfelAtFloor(alpha)) return dfelAtFloor(beta) ? null : Infinity;
  return dfelAtFloor(beta) ? 0 : beta / alpha;
}

function dfelNumber(x) {
  return typeof x === 'number' && Number.isFinite(x) ? x : null;
}

function dfelHeaderNames(headers) {
  if (Array.isArray(headers)) return headers.map(h => (Array.isArray(h) ? h[0] : h));
  return headers && typeof headers === 'object' ? Object.keys(headers) : [];
}

function dfelChi2(block) {
  if (!block) return null;
  return { chi2: dfelNumber(block.chi2), df: block.df ?? null, p: dfelNumber(block.p_value) };
}

export class DFELData {
  constructor(json) {
    this.raw = json;
    this.format = detectMethod(json).format;
    this.numSites = 0;
    this.numTaxa = 0;
    this.patterns = null;
    this.siteResults = [];     // one row per codon site, 1-based
    this.composition = null;   // { overall, cp1, cp2, cp3, heterogeneous, taxa[] }
    this.regimes = [];         // fitted regimes, as written: [{ id, name, psi, meanGC, branchCount }]
    this.branchRegimes = {};   // node -> { regimeId, weights[3] }; optional in v3
    this.root = null;          // { strategy, branch, weight, logL, candidates[] }
    this.fits = {};            // { gtr, directional } -> { logL, aicc }
    this.settings = {};        // { pThreshold, heterogeneityThreshold, rootStrategy, numRegimes }
    this.runThreshold = null;
    this.treeNewick = '';
    this.treeLength = null;
    this.dataset = {};         // { file, geneticCode }
    this.description = null;
    this.citation = null;
    this.version = {};
    this.provenance = {};
    this.statsCache = new Map();

    if (this.format === 'v3') this.parseV3(json);
    else this.parseLegacy(json);

    this.siteBySite = new Map(this.siteResults.map(r => [r.site, r]));
    // The largest GC-ward + AT-ward count, which the table's inline bars are drawn against
    this.gcAtExtent = this.siteResults.reduce((max, r) => Math.max(max, r.tally ? r.tally.gc + r.tally.at : 0), 0);
    this.capabilities = {
      substitutions: this.siteResults.some(r => r.subs),
      directionalModel: this.composition != null || this.regimes.length > 0,
      tree: Boolean(this.treeNewick)
    };
  }

  get formatLabel() {
    return FORMAT_LABELS[this.format];
  }

  // The number of regimes the model fitted: the realised count, which --regimes can force
  // away from what the composition screen alone would pick
  get K() {
    return this.regimes.length || this.settings.numRegimes || null;
  }

  // A fitted regime by its id (not its position: with K = 1 the only regime has id 1)
  regime(id) {
    return this.regimes.find(r => r.id === id) || null;
  }

  parseV3(json) {
    const d = json.site_results?.data || {};
    this.siteResults = (d.site || []).map((site, i) => this.siteRow(site, {
      alpha: d.alpha[i],
      beta: d.beta[i],
      alphaNull: d.alpha_null?.[i],
      lrt: d.lrt[i],
      pValue: d.p_value[i],
      classification: d.classification?.[i] ?? null,
      substitutions: d.substitutions?.[i]
    }));

    const dataset = json.dataset || {};
    this.numSites = dataset.codon_sites ?? this.siteResults.length;
    this.numTaxa = dataset.taxa_count ?? 0;
    this.patterns = dataset.patterns_count ?? null;
    this.dataset = { file: dataset.file ?? null, geneticCode: dataset.genetic_code ?? null };

    const ch = json.compositional_heterogeneity;
    if (ch) {
      this.composition = {
        overall: dfelChi2(ch.overall),
        cp1: dfelChi2(ch.codon_position_1),
        cp2: dfelChi2(ch.codon_position_2),
        cp3: dfelChi2(ch.codon_position_3),
        heterogeneous: Boolean(ch.is_heterogeneous),
        taxa: (ch.taxa || []).map(t => ({ taxon: t.taxon, gcCP3: dfelNumber(t.gc_cp3), gcOverall: dfelNumber(t.gc_overall) }))
      };
    }

    const model = json.directional_model || {};
    this.regimes = (model.regimes || []).map(r => ({
      id: r.id,
      name: r.name,
      psi: dfelNumber(r.psi),
      meanGC: dfelNumber(r.mean_gc),
      branchCount: dfelNumber(r.branch_count)
    }));
    for (const [name, b] of Object.entries(model.branch_regimes || {})) {
      this.branchRegimes[name] = { regimeId: b.regime_id, weights: b.weights || [] };
    }

    const s = json.analysis?.settings || {};
    this.settings = {
      pThreshold: dfelNumber(s.p_value_threshold),
      heterogeneityThreshold: dfelNumber(s.heterogeneity_threshold),
      rootStrategy: s.root_strategy ?? null,
      numRegimes: s.num_regimes ?? null
    };
    this.runThreshold = this.settings.pThreshold;

    const phylogeny = json.phylogeny || {};
    this.treeNewick = phylogeny.newick || '';
    this.treeLength = dfelNumber(phylogeny.tree_length);
    this.root = {
      strategy: this.settings.rootStrategy,
      branch: phylogeny.root_branch ?? null,
      weight: dfelNumber(phylogeny.root_weight),
      logL: dfelNumber(phylogeny.root_log_likelihood),
      candidates: (json.root_selection?.candidates || []).map(c => ({
        branch: c.branch,
        rank: c.rank,
        weight: dfelNumber(c.weight),
        deltaLogL: dfelNumber(c.delta_log_l),
        logL: dfelNumber(c.log_likelihood),
        terminal: Boolean(c.is_terminal)
      })).sort((a, b) => a.rank - b.rank)
    };

    for (const [key, name] of [['gtr', 'gtr'], ['directional', 'directional_mg94']]) {
      const fit = json.model_fits?.[name];
      if (fit) this.fits[key] = { logL: dfelNumber(fit.log_likelihood), aicc: dfelNumber(fit.aicc) };
    }

    const analysis = json.analysis || {};
    this.description = analysis.description ?? null;
    this.citation = json.software?.citation ?? null;
    this.version = { analysis: analysis.version ?? null, software: json.software?.version ?? null };
    this.provenance = { createdAt: json.provenance?.created_at ?? null };
  }

  // hyphy3 --legacy-json: one MLE row per site, columns found by header name
  parseLegacy(json) {
    const names = dfelHeaderNames(json.MLE?.headers);
    const column = Object.fromEntries(Object.entries(DFEL_LEGACY_COLUMNS).map(([key, [name, fallback]]) => {
      const at = names.indexOf(name);
      return [key, at < 0 ? fallback : at];
    }));
    const rows = json.MLE?.content?.['0'] || [];
    this.siteResults = rows.map((row, k) => this.siteRow(k + 1,
      Object.fromEntries(Object.entries(column).map(([key, at]) => [key, row[at]]))));

    const input = json.input || {};
    this.numSites = input['number of sites'] ?? this.siteResults.length;
    this.numTaxa = input['number of sequences'] ?? 0;
    this.dataset = { file: input['file name'] || null, geneticCode: input['genetic code'] ?? null };
    this.treeNewick = input.trees?.['0'] || '';

    const analysis = json.analysis || {};
    this.description = typeof analysis.info === 'string' ? analysis.info : null;
    this.citation = analysis.citation ?? null;
    this.version = { analysis: analysis.version ?? null, software: null };
  }

  siteRow(site, v) {
    const subs = v.substitutions == null ? null : parseDfelSubstitutions(v.substitutions);
    const alpha = v.alpha;
    const beta = v.beta;
    return {
      site,
      alpha,
      beta,
      alphaNull: dfelNumber(v.alphaNull),
      lrt: v.lrt,
      pValue: v.pValue,
      negLogP: -Math.log10(Math.max(v.pValue, 1e-15)),
      omega: dfelOmega(alpha, beta),
      classification: v.classification ?? null,
      substitutions: v.substitutions ?? null,
      subs,
      tally: subs ? dfelChangeTally(subs) : null,
      // Nothing to fit: the writer's 'Conserved' sites, or legacy rows with both rates 0
      invariant: subs ? subs.conserved : (alpha === 0 && beta === 0)
    };
  }

  getSite(site) {
    return this.siteBySite.get(Number(site)) || null;
  }

  getStats(t) {
    if (!this.statsCache.has(t)) {
      const counts = { diversifying: 0, purifying: 0, neutral: 0, invariant: 0 };
      for (const r of this.siteResults) counts[dfelSiteClass(r, t)]++;
      this.statsCache.set(t, {
        totalSites: this.siteResults.length,
        diversifying: counts.diversifying,
        purifying: counts.purifying,
        significant: counts.diversifying + counts.purifying,
        neutral: counts.neutral + counts.invariant,
        invariant: counts.invariant
      });
    }
    return this.statsCache.get(t);
  }

  // The stored classification column against the live call at the run's own cutoff: the
  // counts the file wrote, and how many sites disagree (a parse check)
  runCheck() {
    const t = this.runThreshold;
    const rows = this.siteResults.filter(r => r.classification != null);
    if (t == null || !rows.length) return null;
    const written = { diversifying: 0, purifying: 0, neutral: 0 };
    let mismatches = 0;
    for (const r of rows) {
      if (r.classification in written) written[r.classification]++;
      const live = dfelSiteClass(r, t);
      if ((live === 'invariant' ? 'neutral' : live) !== r.classification) mismatches++;
    }
    return { threshold: t, written, mismatches };
  }
}
