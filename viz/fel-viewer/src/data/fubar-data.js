// FUBAR results: per-site posterior summaries over the alpha x beta rate grid.
// hyphy3 v3 output is read from its schema paths first (site_results, rate_grid,
// model_fits, summary, dataset), with its undeclared legacy-compatibility block (MLE,
// grid, posterior, settings, fits, timers) as the fallback and the only source of the
// per-site posterior grids and chain settings. HyPhy 2.5 and hyphy3 --legacy-json files
// are mapped onto the same feature set; what only they carry (PSRF/Neff columns,
// equilibrium frequencies) is not read.

import { detectMethod, FORMAT_LABELS } from '../core/detect-method.js';
import { poissonBinomialCI } from '../core/stats.js';

// Positions in a legacy MLE.content row. hyphy3 writes MLE.headers as a key-sorted
// object, so the header order cannot name the columns.
const FUBAR_ROW = { alpha: 0, beta: 1, diff: 2, probNeg: 3, probPos: 4, ebf: 5 };

// A site's call at cutoff t: diversifying first, as hyphy3 classifies
export function fubarSiteClass(r, t) {
  if (r.probPos >= t) return 'diversifying';
  return r.probNeg >= t ? 'purifying' : 'neutral';
}

// Expected false diversifying calls at cutoff t: the sum of 1 - P[beta>alpha] over the
// called sites, with its exact 95% Poisson-binomial interval
export function fubarFalsePositives(rows, t) {
  const errors = rows.filter(r => r.probPos >= t).map(r => 1 - r.probPos);
  return {
    calls: errors.length,
    expected: errors.reduce((sum, e) => sum + e, 0),
    ci: poissonBinomialCI(errors)
  };
}

// Prior mass on beta <= alpha (the diagonal included) from the grid weights
export function fubarPriorNonPositive(grid) {
  let mass = 0;
  grid.weight.forEach((w, g) => {
    if (grid.beta[g] <= grid.alpha[g]) mass += w;
  });
  return mass;
}

// 'Substitution rate from nucleotide A to nucleotide G' -> 'A↔G'
function fubarRateLabel(key) {
  const m = key.match(/nucleotide (\w) to nucleotide (\w)/);
  return m ? `${m[1]}↔${m[2]}` : key;
}

function fubarFinite(x) {
  return typeof x === 'number' && Number.isFinite(x) ? x : null;
}

function fubarSum(values) {
  return values.reduce((sum, v) => sum + v, 0);
}

export class FUBARData {
  constructor(json) {
    this.raw = json;
    this.format = detectMethod(json).format;
    this.numSites = 0;
    this.numTaxa = 0;
    this.numPartitions = 1;
    this.siteResults = [];        // one row per codon site, global 1-based site order
    this.grid = null;             // { alpha[], beta[], weight[], dims }, alpha-major
    this.priorNonPositive = null; // grid mass with beta <= alpha
    this.posteriors = new Map();  // site -> its posterior weights over the grid
    this.settings = {};
    this.runThreshold = null;     // the run's own posterior cutoff
    this.runSummary = null;       // v3 summary, frozen at the run's cutoff
    this.gtr = null;
    this.treeNewick = '';
    this.timers = {};
    this.citation = null;
    this.version = {};
    this.statsCache = new Map();

    if (this.format === 'v3') this.parseV3(json);
    else this.parseLegacy(json);
    this.parseCompat(json);

    this.siteBySite = new Map(this.siteResults.map(r => [r.site, r]));
    this.capabilities = {
      sitePosterior: this.posteriors.size > 0,
      multiPartition: this.numPartitions > 1
    };
  }

  get formatLabel() {
    return FORMAT_LABELS[this.format];
  }

  parseV3(json) {
    const data = json.site_results?.data || {};
    const partitions = json.dataset?.partitions || [];
    const partitionOf = new Map();
    partitions.forEach((p, i) => (p.coverage || []).forEach(k => partitionOf.set(k + 1, i)));
    this.numPartitions = Math.max(1, partitions.length);

    this.siteResults = (data.site || []).map((site, i) => this.siteRow(site, partitionOf.get(site) ?? 0, {
      alpha: data.alpha[i],
      beta: data.beta[i],
      diff: data.beta_minus_alpha?.[i] ?? data.beta[i] - data.alpha[i],
      probNeg: data.prob_neg[i],
      probPos: data.prob_pos[i],
      ebf: data.bayes_factor?.[i]
    }));

    const dataset = json.dataset || {};
    this.numSites = dataset.codon_sites ?? this.siteResults.length;
    this.numTaxa = dataset.taxa_count ?? (dataset.taxa || []).length;

    const rg = json.rate_grid;
    if (rg?.data?.posterior_weight) {
      this.grid = {
        alpha: rg.data.alpha,
        beta: rg.data.beta,
        weight: rg.data.posterior_weight,
        dims: rg.dimensions?.[0] ?? Math.round(Math.sqrt(rg.data.posterior_weight.length))
      };
      this.priorNonPositive = fubarFinite(rg.prior_weight_non_positive);
    }

    const s = json.analysis?.settings || {};
    Object.assign(this.settings, {
      method: s.method,
      gridPoints: s.grid_points ?? rg?.dimensions?.[0],
      totalGridPoints: s.total_grid_points ?? rg?.total_points,
      concentration: s.concentration,
      nonZero: s.non_zero,
      treeScaling: s.tree_scaling,
      bestScalerAlpha: rg?.best_scaler_alpha,
      bestScalerBeta: rg?.best_scaler_beta
    });

    const summary = json.summary;
    this.runThreshold = fubarFinite(summary?.posterior_threshold) ?? fubarFinite(s.posterior_threshold);
    if (summary && this.runThreshold != null) {
      this.runSummary = {
        threshold: this.runThreshold,
        diversifying: summary.positive_sites_count,
        purifying: summary.purifying_sites_count,
        expected: summary.expected_false_positives,
        ci: summary.false_positives_ci_95
      };
    }

    const fit = json.model_fits?.nucleotide_gtr;
    if (fit) {
      this.gtr = {
        logL: fit.log_likelihood,
        aicc: fit.aicc,
        rates: fit.rates
          ? [['A↔G', 1], ...Object.entries(fit.rates).map(([k, v]) => [`${k.slice(-2, -1)}↔${k.slice(-1)}`, v])]
          : null,
        treeLength: fit.branch_lengths ? fubarSum(Object.values(fit.branch_lengths)) : null
      };
    }

    const tree = json.phylogeny?.tree;
    if (tree) this.treeNewick = tree.replace(/\{[^}]*\}/g, '');

    const citation = json.analysis?.citations?.[0];
    if (citation) this.citation = { text: citation.citation, doi: citation.doi };
    this.version = {
      analysis: json.analysis?.version,
      software: [json.software?.version, json.software?.git_commit].filter(Boolean).join(' · ')
    };
    const prov = json.provenance || {};
    this.timers.wall = prov.execution?.wall_time_seconds;
    this.timers.command = prov.invocation?.cli_command;
  }

  // Legacy site tables are keyed by partition. HyPhy 2.5 rows are partition-local and
  // placed through each partition's coverage; hyphy3 writes one partition in site order.
  parseLegacy(json) {
    const content = json.MLE?.content || {};
    const partitions = json['data partitions'] || {};
    const keys = Object.keys(content).sort((a, b) => a - b);
    this.numPartitions = Math.max(1, keys.length);
    this.legacySites = {};

    for (const [p, key] of keys.entries()) {
      const coverage = partitions[key]?.coverage?.[0];
      this.legacySites[key] = content[key].map((row, k) => {
        const site = coverage ? coverage[k] + 1 : k + 1;
        this.siteResults.push(this.siteRow(site, p, {
          alpha: row[FUBAR_ROW.alpha],
          beta: row[FUBAR_ROW.beta],
          diff: row[FUBAR_ROW.diff],
          probNeg: row[FUBAR_ROW.probNeg],
          probPos: row[FUBAR_ROW.probPos],
          ebf: row[FUBAR_ROW.ebf]
        }));
        return site;
      });
    }
    this.siteResults.sort((a, b) => a.site - b.site);

    const input = json.input || {};
    this.numSites = input['number of sites'] ?? input.sites ?? this.siteResults.length;
    this.numTaxa = input['number of sequences'] ?? input.sequences ?? 0;
    this.treeNewick = input.trees?.['0'] || '';

    const lengths = json['branch attributes']?.['0'];
    if (lengths) {
      this.gtrTreeLength = fubarSum(Object.values(lengths)
        .map(b => b['Nucleotide GTR'])
        .filter(v => typeof v === 'number'));
    }

    const citation = json.analysis?.citation ?? json.analysis?.reference;
    if (citation) this.citation = { text: citation, doi: null };
    this.version = { analysis: json.analysis?.version, software: json.runtime };
  }

  // The compat block hyphy3 v3 also writes, and the whole of a legacy file: grid,
  // settings, the GTR fit, timers and the per-site posterior grids
  parseCompat(json) {
    if (!this.grid && Array.isArray(json.grid) && json.grid.length) {
      this.grid = {
        alpha: json.grid.map(g => g[0]),
        beta: json.grid.map(g => g[1]),
        weight: json.grid.map(g => g[2]),
        dims: Math.round(Math.sqrt(json.grid.length))
      };
    }
    if (this.grid && this.priorNonPositive == null) {
      this.priorNonPositive = fubarPriorNonPositive(this.grid);
    }

    const s = json.settings || {};
    const settings = this.settings;
    settings.method ??= s.method;
    settings.gridPoints ??= s['grid size'] ?? this.grid?.dims;
    settings.totalGridPoints ??= settings.gridPoints && settings.gridPoints ** 2;
    settings.concentration ??= s.concentration;
    if (settings.nonZero == null && s['non-zero'] != null) settings.nonZero = Boolean(s['non-zero']);
    Object.assign(settings, {
      chains: s.chains,
      chainLength: s['chain-length'],
      burnIn: s['burn-in'],
      samples: s.samples
    });
    this.runThreshold ??= fubarFinite(s.posterior);

    const fit = json.fits?.['Nucleotide GTR'];
    if (fit) {
      const rates = fit['Rate Distributions'];
      this.gtr = Object.assign({ logL: fit['Log Likelihood'], aicc: fit['AIC-c'], rates: null, treeLength: null }, this.gtr, {
        params: fit['estimated parameters']
      });
      if (rates) this.gtr.rates = Object.entries(rates).map(([k, v]) => [fubarRateLabel(k), v]);
    }
    if (this.gtrTreeLength != null) {
      this.gtr ??= { rates: null };
      this.gtr.treeLength = this.gtrTreeLength;
    }

    this.timers.overall = json.timers?.Overall?.timer;

    // hyphy3 indexes posterior['0'] by global site; HyPhy 2.5 by partition-local row
    const posterior = json.posterior || {};
    for (const [key, sites] of Object.entries(posterior)) {
      for (const [k, entry] of Object.entries(sites || {})) {
        const weights = Array.isArray(entry?.[0]) ? entry[0] : entry;
        const site = this.legacySites?.[key]?.[Number(k)] ?? (key === '0' ? Number(k) + 1 : null);
        if (site != null && Array.isArray(weights) && weights.length) this.posteriors.set(site, weights);
      }
    }
  }

  siteRow(site, partition, v) {
    return {
      site,
      partition: partition + 1,
      alpha: v.alpha,
      beta: v.beta,
      diff: v.diff,
      probPos: v.probPos,
      probNeg: v.probNeg,
      probEq: Math.max(0, 1 - v.probPos - v.probNeg),
      ebf: v.ebf ?? null
    };
  }

  getSite(site) {
    return this.siteBySite.get(Number(site)) || null;
  }

  // The site's posterior weights over the grid (alpha-major), or null when the file has none
  sitePosterior(site) {
    return this.posteriors.get(Number(site)) || null;
  }

  // Prior P[beta > alpha]: the grid mass strictly above the diagonal
  get priorPositive() {
    return this.priorNonPositive == null ? null : 1 - this.priorNonPositive;
  }

  getStats(t) {
    if (!this.statsCache.has(t)) {
      let diversifying = 0;
      let purifying = 0;
      for (const r of this.siteResults) {
        const call = fubarSiteClass(r, t);
        if (call === 'diversifying') diversifying++;
        else if (call === 'purifying') purifying++;
      }
      const efp = fubarFalsePositives(this.siteResults, t);
      this.statsCache.set(t, {
        totalSites: this.siteResults.length,
        diversifying,
        purifying,
        neutral: this.siteResults.length - diversifying - purifying,
        expectedFalse: efp.expected,
        falseCI: efp.ci
      });
    }
    return this.statsCache.get(t);
  }
}
