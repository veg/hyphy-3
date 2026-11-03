// GARD results: the best recombination model found by a genetic-algorithm search over
// breakpoint placements, as a set of alignment segments, each with its own tree. A
// breakpoint sits between consecutive segments; KH and SH topology tests of each flank's
// tree against its neighbour's tell topological incongruence (recombination) from
// branch-length heterogeneity alone. hyphy3 v3 output is read from its schema paths
// (dataset, baseline, model_selection, partitions, site_breakpoint_support, topology_test,
// ohmic_network, rate_variation, tree_engine, provenance). HyPhy 2.5 and hyphy3
// --legacy-json files (baselineScore, breakpointData, improvements, siteBreakPointSupport,
// topologyTests) map onto the same segments, breaks and steps; they carry no tree-by-
// partition matrix and no Ohmic block, and HyPhy 2.5 no topology tests. Files from an
// older hyphy3 writer lack the per-breakpoint tests, the matrix and the Ohmic block too
// and read as untested.
//
// Coordinates: partitions, test intervals and breakpoints in topology_test are 1-based
// sites. Model breakpoints, support keys and Ohmic candidates are 0-based positions of
// the last site left of the break; display site = value + off, where off is read off the
// first partition so the two conventions meet (1 in every file seen).

import { detectMethod, FORMAT_LABELS } from '../core/detect-method.js';
import { parseNewick } from '../core/newick.js';

// The run's own family-wise alpha for the KH/SH tests (hyphy3 applies 0.05 / 2B per flank)
export const GARD_RUN_ALPHA = 0.05;

// Breakpoints as a flat list of numbers: HyPhy 2.5 writes each step's set as a one-row matrix
// ([[459]]) and hyphy3 as a plain list
function gardFlat(list) {
  return (Array.isArray(list) ? list.flat(Infinity) : []).map(Number).filter(Number.isFinite);
}

function gardNumber(x) {
  return typeof x === 'number' && Number.isFinite(x) ? x : null;
}

// A number as the report prints it, with a dash for null, NaN and infinities
export function gardFormat(v, digits = 2) {
  const n = gardNumber(v);
  return n == null ? '–' : n.toFixed(digits);
}

// A p-value: four decimals, exponential below 1e-3, '0' for the exact zeros the
// resampling tests report
export function gardP(p) {
  const n = gardNumber(p);
  if (n == null) return '–';
  if (n === 0) return '0';
  return n < 0.001 ? n.toExponential(1) : n.toFixed(4);
}

// A per-test cutoff (alpha / 2B) to three significant digits: 0.025, 0.0125, 2.5e-7
export function gardCutoff(c) {
  return String(Number(c.toPrecision(3)));
}

// A break's live call at family-wise alpha over B breakpoints: each flank must reject
// congruence at alpha / 2B (the writer's Bonferroni rule), by KH, SH or either ('both').
// 'untested' without a test, 'branch' when the tree differs only in branch lengths.
export function gardBreakCall(test, alpha, B, mode = 'both') {
  if (!test) return 'untested';
  const cutoff = alpha / (2 * Math.max(1, B));
  const passes = (left, right) => left != null && right != null && left <= cutoff && right <= cutoff;
  const kh = passes(test.khLeft, test.khRight);
  const sh = passes(test.shLeft, test.shRight);
  const significant = mode === 'kh' ? kh : (mode === 'sh' ? sh : kh || sh);
  return significant ? 'incongruent' : 'branch';
}

// Tree comparison. Splits are over leaf names only (internal NodeN labels are ignored),
// written by the side without the lexicographically first taxon, so they do not depend on
// where the tree is rooted.

function gardLeaves(node, into = []) {
  if (!node.children.length) into.push(node.name);
  node.children.forEach(child => gardLeaves(child, into));
  return into;
}

// Whether every leaf of a Newick tree has its own name
export function gardUniqueTaxa(newick) {
  const root = parseNewick(newick);
  if (!root) return false;
  const names = gardLeaves(root);
  return new Set(names).size === names.length;
}

const GARD_SEP = '\u001f';

// The side of a bipartition of taxa (sorted) that excludes the first taxon, as a key
function gardSideKey(side, taxa) {
  const inside = new Set(side);
  const keyed = inside.has(taxa[0]) ? taxa.filter(t => !inside.has(t)) : [...inside];
  return keyed.sort().join(GARD_SEP);
}

// { taxa (sorted), splits: Set of non-trivial split keys } for a Newick tree, over the
// taxa in keep when given
export function gardSplits(newick, keep = null) {
  const root = parseNewick(newick);
  const all = root ? gardLeaves(root).filter(t => !keep || keep.has(t)).sort() : [];
  const splits = new Set();
  const walk = (node) => {
    const leaves = node.children.length ? node.children.flatMap(walk) : (!keep || keep.has(node.name) ? [node.name] : []);
    if (node.children.length && leaves.length >= 2 && all.length - leaves.length >= 2) splits.add(gardSideKey(leaves, all));
    return leaves;
  };
  if (root) walk(root);
  return { taxa: all, splits };
}

const gardSplitCache = new Map();

function gardSplitsOf(tree, keep = null) {
  if (typeof tree !== 'string') return tree;
  if (keep) return gardSplits(tree, keep);
  if (!gardSplitCache.has(tree)) gardSplitCache.set(tree, gardSplits(tree));
  return gardSplitCache.get(tree);
}

// Robinson-Foulds distance between two trees (Newick text or gardSplits results) and its
// maximum, the number of non-trivial splits in both. Trees over different taxa are
// compared on the taxa they share.
export function gardRF(a, b) {
  let x = gardSplitsOf(a);
  let y = gardSplitsOf(b);
  if (x.taxa.length !== y.taxa.length || x.taxa.some((t, i) => t !== y.taxa[i])) {
    const shared = new Set(x.taxa.filter(t => y.taxa.includes(t)));
    x = gardSplits(typeof a === 'string' ? a : '', shared);
    y = gardSplits(typeof b === 'string' ? b : '', shared);
  }
  let common = 0;
  for (const split of x.splits) if (y.splits.has(split)) common++;
  return { rf: x.splits.size + y.splits.size - 2 * common, max: x.splits.size + y.splits.size };
}

// Each taxon's sister clade, as a key: the other leaves under its parent, as a bipartition
// of the remaining taxa so the key does not depend on the root
function gardSisters(newick, taxa) {
  const root = parseNewick(newick);
  const sisters = new Map();
  const keep = new Set(taxa);
  const walk = (node) => {
    if (!node.children.length) return keep.has(node.name) ? [node.name] : [];
    const below = node.children.map(walk);
    const leaves = below.flat();
    node.children.forEach((child, i) => {
      if (child.children.length || !below[i].length) return;
      const x = below[i][0];
      const rest = taxa.filter(t => t !== x);
      const sister = leaves.filter(t => t !== x);
      if (sister.length) sisters.set(x, gardSideKey(sister, rest));
    });
    return leaves;
  };
  if (root) walk(root);
  return sisters;
}

// The taxa whose sister clade differs between two trees, over the taxa they share, in
// alphabetical order: where the topologies disagree about who a taxon's nearest relatives are.
// A moved clade is flagged at its junction only (its members keep their sisters); finding the
// smallest set of taxa to prune for the trees to agree is a maximum-agreement-subtree search,
// far too slow at hundreds of taxa.
export function gardMovedTaxa(a, b) {
  const leavesOf = (tree) => { const r = parseNewick(tree); return r ? gardLeaves(r) : []; };
  const inB = new Set(leavesOf(b));
  const taxa = leavesOf(a).filter(t => inB.has(t)).sort();
  const sa = gardSisters(a, taxa);
  const sb = gardSisters(b, taxa);
  return taxa.filter(t => sa.get(t) !== sb.get(t));
}

// One test row as the writer keys it (snake_case, same keys in hyphy3 --legacy-json)
function gardTest(t) {
  return {
    site: t.breakpoint,
    leftInterval: t.left_partition ?? null,
    rightInterval: t.right_partition ?? null,
    dLogLLeft: gardNumber(t.delta_log_l_left),
    dLogLRight: gardNumber(t.delta_log_l_right),
    khLeft: gardNumber(t.kh_p_value_left),
    khRight: gardNumber(t.kh_p_value_right),
    shLeft: gardNumber(t.sh_p_value_left),
    shRight: gardNumber(t.sh_p_value_right),
    significant: Boolean(t.statistically_significant)
  };
}

function gardMode(testType) {
  const t = String(testType ?? 'both').toLowerCase();
  return t === 'kh' || t === 'sh' ? t : 'both';
}

export class GardResults {
  constructor(json) {
    this.raw = json;
    this.format = detectMethod(json).format;
    this.legacy = this.format !== 'v3';
    this.numSites = 0;
    this.numTaxa = 0;
    this.dataset = {};         // { file, variableSites, potentialBreakpoints, minPartition }
    this.baselineAicc = null;
    this.bestAicc = null;
    this.singleTreeAicc = null;
    this.baseline = {};        // { logL, tree, rates, frequencies }
    this.segments = [];        // [{ key: 'S0', start, end, tree, logL }], 1-based inclusive
    this.breaks = [];          // [{ index, site, left, right, test, diff }], site after which the tree changes
    this.support = [];         // [{ site, p }] sorted: Akaike weight of a breakpoint after the site
    this.steps = [];           // [{ k, deltaAicc, bps }]: the step-up procedure, bps as display sites
    this.matrix = [];          // [{ partition, interval, bestTree, logL, kh, sh }]
    this.topology = null;      // { supported, deltaAicc, testType, replicates, anySignificant }
    this.testMode = 'both';
    this.ohmic = null;         // { candidates, bridges, detected, count, logL, aicc, deltaAicc } | null
    this.rateVariation = null;
    this.treeEngine = null;
    this.modelsEvaluated = null;
    this.timeElapsed = null;
    this.software = {};
    this.provenance = {};
    this.description = null;
    this.citation = null;
    this.off = 1;
    this.statsCache = new Map();
    this.compareCache = new Map();

    if (this.format === 'v3') this.parseV3(json);
    else this.parseLegacy(json);

    this.deltaAicc = this.baselineAicc == null || this.bestAicc == null ? null : this.baselineAicc - this.bestAicc;
    this.treeNewick = this.segments[0]?.tree ?? '';
    this.capabilities = {
      topology: this.breaks.some(b => b.test),
      matrix: this.matrix.length > 0,
      ohmic: Boolean(this.ohmic),
      legacy: this.legacy
    };
    this.runThreshold = this.capabilities.topology ? GARD_RUN_ALPHA : null;
  }

  get formatLabel() {
    return FORMAT_LABELS[this.format];
  }

  // Display site of a 0-based position
  display(position) {
    return position + this.off;
  }

  parseV3(json) {
    const ds = json.dataset || {};
    const ms = json.model_selection || {};
    const tt = json.topology_test || {};
    this.numSites = ds.sites_count ?? 0;
    this.numTaxa = ds.taxa_count ?? 0;
    this.dataset = {
      file: ds.file ?? null,
      variableSites: ds.variable_sites_count ?? null,
      potentialBreakpoints: ds.potential_breakpoints ?? null,
      minPartition: ds.minimum_partition_size ?? null
    };
    const base = json.baseline || {};
    this.baselineAicc = gardNumber(base.aicc);
    this.baseline = { logL: gardNumber(base.log_likelihood), tree: base.tree ?? '', rates: base.rates ?? null, frequencies: base.frequencies ?? null };
    this.bestAicc = gardNumber(ms.best_model_aicc);
    this.singleTreeAicc = gardNumber(tt.single_tree_aicc);
    this.analysis = json.analysis || {};
    this.citation = json.analysis?.citation ?? null;
    this.software = json.software || {};
    this.modelsEvaluated = json.total_models_evaluated ?? null;
    this.timeElapsed = gardNumber(json.time_elapsed);
    this.rateVariation = json.rate_variation ?? null;
    this.treeEngine = json.tree_engine ?? null;
    const input = json.provenance?.inputs?.alignment;
    this.provenance = { alignmentSha256: input?.sha256 ?? null, alignmentFile: input?.path ?? null };

    const parts = (json.partitions || []).map(p => ({ start: p.start_site, end: p.end_site, tree: p.tree ?? '', logL: gardNumber(p.log_likelihood) }));
    // Files from the oldest writer have no partitions at all: one segment on the baseline tree
    this.setSegments(parts.length ? parts : [{ start: 1, end: this.numSites, tree: this.baseline.tree, logL: null }]);
    this.setOffset(ms.best_model_breakpoints);
    this.setBreaks((tt.breakpoints || []).map(gardTest));
    this.steps = (ms.improvements || []).map(s => ({
      k: s.breakpoints_count,
      deltaAicc: gardNumber(s.delta_aicc),
      bps: Array.isArray(s.breakpoints) ? gardFlat(s.breakpoints).map(b => this.display(b)) : null
    }));
    this.setSupport(json.site_breakpoint_support);
    this.topology = tt.single_tree_aicc === undefined && tt.topological_incongruence_supported === undefined ? null : {
      supported: tt.topological_incongruence_supported ?? null,
      deltaAicc: gardNumber(tt.delta_aicc_topology),
      testType: tt.test_type ?? null,
      replicates: tt.bootstrap_replicates ?? null,
      anySignificant: tt.any_incongruence_significant ?? null
    };
    this.testMode = gardMode(tt.test_type);
    this.matrix = (tt.partition_tree_matrix || []).map(m => ({
      partition: m.partition_index,
      interval: m.interval,
      bestTree: m.best_tree_index,
      logL: m.tree_log_likelihoods || [],
      kh: m.kh_p_values || [],
      sh: m.sh_p_values || []
    }));
    const oh = json.ohmic_network;
    if (oh) {
      this.ohmic = {
        candidates: (oh.bridge_candidate_breakpoints || []).map(b => this.display(b)).sort((a, b) => a - b),
        bridges: (oh.bridges || []).map(b => ({
          index: b.bridge_index, hub: b.hub_node, node: b.bridge_node, support: gardNumber(b.ohmic_support), stress: gardNumber(b.shunt_stress)
        })),
        detected: Boolean(oh.bridges_detected),
        count: oh.num_bridges ?? (oh.bridges || []).length,
        logL: gardNumber(oh.ohmic_log_likelihood),
        aicc: gardNumber(oh.ohmic_aicc),
        deltaAicc: gardNumber(oh.delta_aicc_vs_baseline)
      };
    }
  }

  parseLegacy(json) {
    const input = json.input || {};
    this.numSites = input['number of sites'] ?? 0;
    this.numTaxa = input['number of sequences'] ?? 0;
    this.dataset = { file: input['file name'] ?? null, variableSites: null, potentialBreakpoints: json.potentialBreakpoints ?? null, minPartition: null };
    this.baselineAicc = gardNumber(json.baselineScore);
    this.bestAicc = gardNumber(json.bestModelAICc);
    this.singleTreeAicc = gardNumber(json.singleTreeAICc);
    this.modelsEvaluated = json.totalModelCount ?? null;
    this.timeElapsed = gardNumber(json.timeElapsed);
    this.description = typeof json.analysis?.info === 'string' ? json.analysis.info.replace(/\s+/g, ' ').trim() : null;
    this.citation = json.analysis?.citation?.replace(/\*\*|_/g, '') ?? null;

    // hyphy3 writes the interval ends as strings, HyPhy 2.5 as numbers
    const byIndex = Object.entries(json.breakpointData || {}).sort((a, b) => Number(a[0]) - Number(b[0]));
    const parts = byIndex.map(([, d]) => ({ start: Number(d.bps?.[0]?.[0]), end: Number(d.bps?.[0]?.[1]), tree: d.tree ?? '', logL: null }));
    if (!parts.length) parts.push({ start: 1, end: this.numSites, tree: json.trees?.newickString ?? '', logL: null });
    this.setSegments(parts);
    this.baseline = { logL: null, tree: parts.length === 1 ? parts[0].tree : '', rates: null, frequencies: null };

    const steps = Object.entries(json.improvements || {}).sort((a, b) => Number(a[0]) - Number(b[0]))
      .filter(([, s]) => Array.isArray(s.breakpoints))
      .map(([, s]) => ({ bps: gardFlat(s.breakpoints), deltaAicc: gardNumber(s.deltaAICc) }));
    this.setOffset(steps.length ? steps[steps.length - 1].bps : []);
    this.steps = steps.map(s => ({ k: s.bps.length, deltaAicc: s.deltaAicc, bps: s.bps.map(b => this.display(b)) }));
    this.setSupport(json.siteBreakPointSupport);

    const delta = this.singleTreeAicc == null || this.bestAicc == null ? null : this.singleTreeAicc - this.bestAicc;
    const tt = json.topologyTests;
    this.topology = { supported: delta == null ? null : delta > 0, deltaAicc: delta, testType: tt?.testType ?? null, replicates: tt?.replicates ?? null, anySignificant: tt?.anySignificant ?? null };
    this.testMode = gardMode(tt?.testType);
    this.setBreaks((tt?.breakpoints || []).map(gardTest));
  }

  setSegments(parts) {
    this.segments = parts.map((p, i) => ({ key: `S${i}`, start: p.start, end: p.end, tree: p.tree, logL: p.logL }));
  }

  // The shift from 0-based breakpoints to display sites: the first partition ends at the
  // site after the first breakpoint's position
  // HyPhy 2.5 writes the 0-based break position itself as the left interval's end (1..459,
  // 460..762 for a break at 459, though its filter is sites 0..459): the offset reads 0, and the
  // intervals are moved onto 1-based sites (1..460, 461..762) so every producer agrees.
  setOffset(bestBreakpoints) {
    const flat = gardFlat(bestBreakpoints);
    const first = flat.length ? Math.min(...flat) : null;
    this.off = first != null && this.segments.length > 1 ? this.segments[0].end - first : 1;
    if (this.off === 0) {
      this.segments.forEach((s, i) => {
        if (i > 0) s.start += 1;
        if (i < this.segments.length - 1) s.end += 1;
      });
      this.off = 1;
    }
  }

  // diff is the trees' comparison across the break ({ rf, max, moved }), computed on first use
  setBreaks(tests) {
    const bySite = new Map(tests.map(t => [t.site, t]));
    const data = this;
    this.breaks = this.segments.slice(0, -1).map((left, i) => ({
      index: i + 1,
      site: left.end,
      left,
      right: this.segments[i + 1],
      test: bySite.get(left.end) ?? null,
      get diff() { return data.compare(i, i + 1); }
    }));
  }

  setSupport(support) {
    this.support = Object.entries(support || {})
      .map(([key, p]) => ({ site: Number(key) + this.off, p: gardNumber(p) }))
      .filter(s => s.p != null)
      .sort((a, b) => a.site - b.site);
  }

  getBreak(site) {
    return this.breaks.find(b => b.site === Number(site)) ?? null;
  }

  // A break's live call at family-wise alpha
  callOf(brk, alpha) {
    return gardBreakCall(brk.test, alpha, this.breaks.length, this.testMode);
  }

  // RF distance and moved taxa between two segments' trees, by index; cached. Null when a tree
  // is missing or repeats a taxon name (splits over names would then be meaningless).
  compare(i, j) {
    const key = `${i}-${j}`;
    if (!this.compareCache.has(key)) {
      const a = this.segments[i]?.tree;
      const b = this.segments[j]?.tree;
      const ok = a && b && gardUniqueTaxa(a) && gardUniqueTaxa(b);
      this.compareCache.set(key, ok ? Object.assign(gardRF(a, b), { moved: gardMovedTaxa(a, b) }) : null);
    }
    return this.compareCache.get(key);
  }

  // The breaks' counts at family-wise alpha, as the KPI bar and table filters read them
  getStats(alpha) {
    if (!this.statsCache.has(alpha)) {
      const calls = this.breaks.map(b => this.callOf(b, alpha));
      const count = (call) => calls.filter(c => c === call).length;
      this.statsCache.set(alpha, {
        breakpoints: this.breaks.length,
        incongruent: count('incongruent'),
        branchOnly: count('branch'),
        untested: count('untested'),
        segments: this.segments.length,
        totalSites: this.numSites
      });
    }
    return this.statsCache.get(alpha);
  }
}
