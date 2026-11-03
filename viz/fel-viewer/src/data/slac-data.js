// SLAC results: per-site synonymous and non-synonymous substitution counts over the
// maximum-likelihood ancestral reconstruction, each direction tested by its own one-sided
// binomial p-value. hyphy3 v3 output is read from its schema paths (site_results, dataset,
// model_fits, phylogeny, statistical_tests). HyPhy 2.5 and hyphy3 --legacy-json files map
// onto the same feature set: per-site values from the ambiguity-resolved table (RESOLVED,
// which is what v3 reports), leaf codons from the leaves' branch attributes, and the
// sampled-ancestor ranges of S, N, dS and dN. What only they carry (ancestral codons,
// per-branch counts, the by-branch table, the AVERAGED table, sampled p-value ranges) is
// not read.

import { detectMethod, FORMAT_LABELS } from '../core/detect-method.js';
import { parseNewick } from '../core/newick.js';
import { translateCodon } from './genetic-code.js';

// The branch lengths the reconstruction ran on, as tree views read them from nodeMeta
export const SLAC_LENGTH_KEY = 'Global MG94xREV';

// Legacy by-site columns by header name, with the HyPhy 2.5 position as the fallback
const SLAC_LEGACY_COLUMNS = {
  ES: ['ES', 0],
  EN: ['EN', 1],
  S: ['S', 2],
  N: ['N', 3],
  P_S: ['P[S]', 4],
  dS: ['dS', 5],
  dN: ['dN', 6],
  dNmdS: ['dN-dS', 7],
  pPos: ['P [dN/dS > 1]', 8],
  pNeg: ['P [dN/dS < 1]', 9],
  tbl: ['Total branch length', 10]
};

// The quantities sampled over ancestors in both formats
const SLAC_SAMPLED = ['S', 'N', 'dS', 'dN'];

// A site's call at cutoff t. p_positive + p_negative >= 1, so below t = 0.5 at most one
// direction passes; diversifying is checked first, as hyphy3 classifies. A site with no
// inferred substitutions has both p = 1 and is the 'invariant' display state of neutral.
export function slacSiteClass(r, t) {
  if (r.pPos <= t) return 'diversifying';
  if (r.pNeg <= t) return 'purifying';
  return r.noSubs ? 'invariant' : 'neutral';
}

// The most common amino acid among the leaves' codons, ignoring gaps and codons that do
// not translate to one residue; a tie goes to the first alphabetically, so the call does
// not depend on taxon order. null when no codon translates.
export function consensusResidue(codons) {
  const counts = new Map();
  for (const codon of codons) {
    const aa = translateCodon(codon);
    if (/^[A-Z*]$/.test(aa)) counts.set(aa, (counts.get(aa) || 0) + 1);
  }
  let best = null;
  for (const [aa, n] of counts) {
    if (!best || n > counts.get(best) || (n === counts.get(best) && aa < best)) best = aa;
  }
  return best;
}

// Column names of legacy MLE.headers: [[name, description], ...] or { name: description }
function slacHeaderNames(headers) {
  if (Array.isArray(headers)) return headers.map(h => (Array.isArray(h) ? h[0] : h));
  return headers && typeof headers === 'object' ? Object.keys(headers) : [];
}

function slacNumber(x) {
  return typeof x === 'number' && Number.isFinite(x) ? x : null;
}

// Rate-distribution values may be wrapped as [[value, weight]]
function slacUnwrap(x) {
  let v = x;
  while (Array.isArray(v)) v = v[0];
  return slacNumber(v);
}

function slacLeafNames(newick) {
  const leaves = [];
  const walk = (node) => {
    if (!node.children.length) leaves.push(node.name);
    node.children.forEach(walk);
  };
  const root = parseNewick(newick);
  if (root) walk(root);
  return leaves;
}

export class SLACData {
  constructor(json) {
    this.raw = json;
    this.format = detectMethod(json).format;
    this.numSites = 0;
    this.numTaxa = 0;
    this.numPartitions = 1;
    this.siteResults = [];     // one row per codon site, global 1-based site order
    this.leafCodons = {};      // taxon -> codons in site order
    this.nodeMeta = {};        // node -> { branch_lengths: { [SLAC_LENGTH_KEY]: length } }
    this.fits = [];            // [{ name, logL, aicc, params }] in file order
    this.globalOmega = null;
    this.settings = {};        // { samples, testBranches, geneticCode }
    this.testing = {};         // { tested, background, testedLength }
    this.runThreshold = null;  // the run's own p cutoff
    this.runSummary = null;    // v3 counts, frozen at the run's cutoff
    this.treeNewick = '';
    this.description = null;
    this.citation = null;
    this.version = {};
    this.provenance = {};
    this.statsCache = new Map();

    if (this.format === 'v3') this.parseV3(json);
    else this.parseLegacy(json);

    for (const r of this.siteResults) {
      r.consensusAA = consensusResidue(Object.values(this.leafCodons).map(codons => codons[r.site - 1]));
    }
    this.siteBySite = new Map(this.siteResults.map(r => [r.site, r]));
    this.totalS = this.siteResults.reduce((sum, r) => sum + r.S, 0);
    this.totalN = this.siteResults.reduce((sum, r) => sum + r.N, 0);
    // The largest |scaled dN - dS|, which inline effect bars are drawn against
    this.effectExtent = this.siteResults.reduce((max, r) => Math.max(max, Math.abs(r.dNmdS ?? 0)), 0);
    this.capabilities = {
      sampling: this.siteResults.some(r => r.sampled),
      multiPartition: this.numPartitions > 1,
      leafCodons: Object.keys(this.leafCodons).length > 0,
      tree: Boolean(this.treeNewick)
    };
  }

  get formatLabel() {
    return FORMAT_LABELS[this.format];
  }

  parseV3(json) {
    const d = json.site_results?.data || {};
    const partitions = json.dataset?.partitions || [];
    this.numPartitions = Math.max(1, partitions.length);
    const partitionOf = (site) => {
      const i = partitions.findIndex(p => (p.span || []).some(([from, to]) => site >= from && site <= to));
      return Math.max(0, i) + 1;
    };
    const sampled = (key, i) => (d[`${key}_median`] ? {
      median: slacNumber(d[`${key}_median`][i]),
      lo: slacNumber(d[`${key}_ci_lower`]?.[i]),
      hi: slacNumber(d[`${key}_ci_upper`]?.[i])
    } : null);

    this.siteResults = (d.site || []).map((site, i) => this.siteRow(site, partitionOf(site), {
      ES: d.ES[i], EN: d.EN[i], S: d.S[i], N: d.N[i], P_S: d.P_S?.[i],
      dS: d.dS[i], dN: d.dN[i], dNmdS: d.dN_minus_dS[i],
      pPos: d.p_positive[i], pNeg: d.p_negative[i], tbl: d.total_branch_length?.[i]
    }, Object.fromEntries(SLAC_SAMPLED.map(k => [k, sampled(k, i)]))));

    const dataset = json.dataset || {};
    this.numSites = dataset.codon_sites ?? this.siteResults.length;
    this.numTaxa = dataset.taxa_count ?? (dataset.taxa || []).length;
    for (const [taxon, seq] of Object.entries(dataset.sequences || {})) {
      this.leafCodons[taxon] = Array.from({ length: this.numSites }, (_, k) => seq.substr(k * 3, 3).toUpperCase());
    }

    for (const [name, fit] of Object.entries(json.model_fits || {})) {
      this.fits.push({ name, logL: slacNumber(fit.log_likelihood), aicc: slacNumber(fit.aicc), params: fit.parameters_count ?? null });
    }
    this.globalOmega = slacNumber(json.model_fits?.[SLAC_LENGTH_KEY]?.global_parameters?.omega);

    const s = json.analysis?.settings || {};
    const summary = json.statistical_tests?.site_level_summary || {};
    this.settings = { samples: s.samples ?? null, testBranches: s.test_branches ?? null, geneticCode: dataset.genetic_code?.name ?? null };
    this.testing = {
      tested: summary.tested_branches_count ?? null,
      background: summary.background_branches_count ?? null,
      testedLength: slacNumber(summary.total_tested_tree_length)
    };
    this.runThreshold = slacNumber(s.pvalue_threshold) ?? slacNumber(summary.threshold);
    if (summary.diversifying_count != null && slacNumber(summary.threshold) != null) {
      this.runSummary = {
        threshold: summary.threshold,
        diversifying: summary.diversifying_count,
        purifying: summary.purifying_count,
        neutral: summary.neutral_count
      };
    }

    const phylogeny = json.phylogeny || {};
    this.treeNewick = (phylogeny.newick || '').replace(/\{[^}]*\}/g, '');
    for (const [name, meta] of Object.entries(phylogeny.nodes || {})) {
      if (slacNumber(meta.length) != null) this.nodeMeta[name] = { branch_lengths: { [SLAC_LENGTH_KEY]: meta.length } };
    }

    const analysis = json.analysis || {};
    this.description = analysis.description ?? null;
    const citation = analysis.citations?.[0];
    if (citation) this.citation = { text: citation.citation, doi: citation.doi ?? null };
    this.version = {
      analysis: analysis.version,
      software: [json.software?.version, json.software?.git_commit].filter(Boolean).join(' · ')
    };
    const prov = json.provenance || {};
    this.provenance = {
      wallSeconds: slacNumber(prov.execution?.wall_time_seconds),
      command: prov.invocation?.cli_command ?? null,
      alignmentSha256: prov.inputs?.alignment?.sha256 ?? null
    };
  }

  // Legacy site tables are keyed by partition, rows partition-local and placed through
  // each partition's coverage
  parseLegacy(json) {
    const content = json.MLE?.content || {};
    const names = slacHeaderNames(json.MLE?.headers);
    const column = Object.fromEntries(Object.entries(SLAC_LEGACY_COLUMNS).map(([key, [name, fallback]]) => {
      const at = names.indexOf(name);
      return [key, at < 0 ? fallback : at];
    }));
    const partitions = json['data partitions'] || {};
    const keys = Object.keys(content).sort((a, b) => a - b);
    this.numPartitions = Math.max(1, keys.length);
    const attributes = json['branch attributes'] || {};
    const trees = json.input?.trees || {};
    const sample = (name, key, k) => json[name]?.[key]?.RESOLVED?.[k];

    for (const [p, key] of keys.entries()) {
      const coverage = partitions[key]?.coverage?.[0];
      const leaves = slacLeafNames(trees[key] ?? trees['0']);
      (content[key]['by-site']?.RESOLVED || []).forEach((row, k) => {
        const site = coverage ? coverage[k] + 1 : k + 1;
        const values = Object.fromEntries(Object.entries(column).map(([name, at]) => [name, row[at]]));
        const sampled = Object.fromEntries(SLAC_SAMPLED.map(name => {
          const [median, lo, hi] = ['sample-median', 'sample-2.5', 'sample-97.5'].map(s => sample(s, key, k));
          return [name, median ? {
            median: slacNumber(median[column[name]]),
            lo: slacNumber(lo?.[column[name]]),
            hi: slacNumber(hi?.[column[name]])
          } : null];
        }));
        this.siteResults.push(this.siteRow(site, p + 1, values, sampled));
        for (const leaf of leaves) {
          const codon = attributes[key]?.[leaf]?.codon?.[0]?.[k];
          if (typeof codon === 'string') (this.leafCodons[leaf] ||= [])[site - 1] = codon.toUpperCase();
        }
      });
    }
    this.siteResults.sort((a, b) => a.site - b.site);

    const input = json.input || {};
    this.numSites = input['number of sites'] ?? this.siteResults.length;
    this.numTaxa = input['number of sequences'] ?? 0;
    this.treeNewick = (trees['0'] || '').replace(/\{[^}]*\}/g, '');
    for (const [name, attrs] of Object.entries(attributes['0'] || {})) {
      const length = slacNumber(attrs?.[SLAC_LENGTH_KEY]);
      if (length != null) this.nodeMeta[name] = { branch_lengths: { [SLAC_LENGTH_KEY]: length } };
    }

    for (const [name, fit] of Object.entries(json.fits || {})) {
      this.fits.push({
        name,
        logL: slacNumber(fit['Log Likelihood'] ?? fit['Log-likelihood']),
        aicc: slacNumber(fit['AIC-c']),
        params: fit['estimated parameters'] ?? null
      });
    }
    const rates = json.fits?.[SLAC_LENGTH_KEY]?.['Rate Distributions'] || {};
    const omegaKey = Object.keys(rates).find(k => /non-synonymous\/synonymous rate ratio/.test(k));
    this.globalOmega = omegaKey ? slacUnwrap(rates[omegaKey]) : null;

    const analysis = json.analysis || {};
    const s = analysis.settings || {};
    this.settings = { samples: s.samples ?? null, testBranches: s.branches ?? null, geneticCode: input['genetic code'] ?? null };
    const tested = Object.values(json.tested?.['0'] || {});
    this.testing = {
      tested: tested.length ? tested.filter(v => v === 'test').length : null,
      background: tested.length ? tested.filter(v => v !== 'test').length : null,
      testedLength: null
    };
    this.runThreshold = slacNumber(s.pvalue);

    this.description = typeof analysis.info === 'string' ? analysis.info.replace(/\s+/g, ' ').trim() : null;
    if (analysis.citation) this.citation = { text: analysis.citation, doi: null };
    this.version = { analysis: analysis.version, software: json.runtime };
    const total = json.timers?.['Total time'];
    this.provenance = { wallSeconds: slacNumber(total?.timer ?? total), command: null, alignmentSha256: null };
  }

  // Where dS is undefined (ES = 0) its sampled range is too: v3 writes it as null,
  // --legacy-json as 0, and both read as blank
  siteRow(site, partition, v, sampled) {
    const S = v.S ?? 0;
    const N = v.N ?? 0;
    const dSUndefined = slacNumber(v.dS) == null;
    if (dSUndefined) sampled.dS = null;
    return {
      site,
      partition,
      ES: v.ES,
      EN: v.EN,
      S,
      N,
      P_S: slacNumber(v.P_S),
      dS: slacNumber(v.dS),
      dN: slacNumber(v.dN),
      dNmdS: slacNumber(v.dNmdS),
      pPos: v.pPos,
      pNeg: v.pNeg,
      tbl: slacNumber(v.tbl),
      sampled: SLAC_SAMPLED.some(k => sampled[k]) ? sampled : null,
      noSubs: S + N === 0,
      dSUndefined
    };
  }

  getSite(site) {
    return this.siteBySite.get(Number(site)) || null;
  }

  // The leaves' codons at a site, by taxon name: [{ taxon, codon, aa }]
  siteCodons(site) {
    return Object.entries(this.leafCodons)
      .map(([taxon, codons]) => ({ taxon, codon: codons[site - 1] }))
      .filter(c => typeof c.codon === 'string')
      .sort((a, b) => a.taxon.localeCompare(b.taxon))
      .map(c => Object.assign(c, { aa: translateCodon(c.codon) }));
  }

  getStats(t) {
    if (!this.statsCache.has(t)) {
      const counts = { diversifying: 0, purifying: 0, neutral: 0, invariant: 0 };
      for (const r of this.siteResults) counts[slacSiteClass(r, t)]++;
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
}
