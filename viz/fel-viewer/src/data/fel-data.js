import { translateCodon, getRateSparkGlyphs } from './genetic-code.js';

// HyPhy 2.5 column order of each legacy site table. hyphy3 --legacy-json writes MEME's
// MLE.headers as an object, whose keys come out alphabetized while the content rows keep this
// order, so object-form headers are mapped back through the matching layout.
const LEGACY_COLUMN_LAYOUTS = [
  ['alpha', 'beta', 'alpha=beta', 'LRT', 'p-value', 'Total branch length'],
  ['&alpha;', '&beta;<sup>1</sup>', 'p<sup>1</sup>', '&beta;<sup>+</sup>', 'p<sup>+</sup>', 'LRT', 'p-value',
   '# branches under selection', 'Total branch length', 'MEME LogL', 'FEL LogL', 'LRT MEME vs FEL',
   'FEL &alpha;', 'FEL &beta;'],
];

function orderedLegacyHeaders(headers) {
  if (Array.isArray(headers)) return headers.map(h => (Array.isArray(h) ? h[0] : h) || '');
  if (!headers || typeof headers !== 'object') return [];
  const keys = Object.keys(headers);
  const layout = LEGACY_COLUMN_LAYOUTS.find(l => l.length === keys.length && l.every(name => name in headers));
  if (!layout) throw new Error('Unrecognized legacy site-table headers; column order cannot be determined');
  return layout;
}

// The p-value cutoff the run itself used, or null when the file does not say: v3
// writes analysis.settings.pvalue_threshold, legacy FEL 'pvalue' and legacy MEME 'p-value'
export function runPThreshold(settings) {
  const t = settings?.pvalue_threshold ?? settings?.pvalue ?? settings?.['p-value'];
  return Number.isFinite(t) && t > 0 ? t : null;
}

export class FELData {
  constructor(json) {
    this.raw = json;
    this.taxa = [];
    this.numSites = 0;
    this.numTaxa = 0;
    this.sequences = {};          // taxon -> nucleotide string
    this.translated = {};         // taxon -> array of amino acids
    this.consensusAA = [];        // site -> most frequent AA
    this.treeNewick = '';
    this.treeNodes = {};          // node name -> metadata
    this.siteResults = [];        // site 1..N
    this.eventsBySite = {};       // 1-based site -> events[]
    this.eventsByBranch = {};     // branch name -> events[]
    this.substitutionEvents = []; // raw substitution_events from modern JSON
    this.ancestralSequences = {}; // internal node name -> reconstructed nucleotide sequence
    this.ancestralTranslated = {};// internal node name -> array of amino acids
    this.ancestralTaxa = [];      // list of internal node names in tree preorder
    this.rootCodons = [];         // site -> root codon
    this.rootAA = [];             // site -> root amino acid
    this.globalOmega = 1.0;
    this.hasBackground = false;
    this.numTested = 0;
    this.numBackground = 0;
    this.provenance = null;
    this.analysis = null;
    this.runThreshold = null;     // the run's own p-value cutoff, the shell's starting threshold
    this.analysisType = 'fel';
    this.isMEME = false;

    this.parse(json);
  }

  parse(json) {
    // 1. Detect analysis method (FEL vs MEME vs Contrast-FEL)
    const aId = (json.analysis?.id || '').toLowerCase();
    const aInfo = (json.analysis?.info || '').toLowerCase();
    const isContrastData = aId.includes('contrast') || !!json.site_results?.data?.penalized_lrt;
    const isMemeData = !isContrastData && (
      aId === 'meme' ||
      aInfo.includes('meme') ||
      !!json.site_results?.data?.beta_plus ||
      orderedLegacyHeaders(json.MLE?.headers).some(h => {
        const name = h.toLowerCase();
        return name.includes('beta') && (name.includes('+') || name.includes('2') || name.includes('plus'));
      })
    );
    
    this.analysisType = isContrastData ? 'contrast-fel' : (isMemeData ? 'meme' : 'fel');
    this.isMEME = isMemeData;
    this.isContrast = isContrastData;

    // 2. Detect format (Modern v3 vs Legacy 2.5)
    if (json.schema_version || json.site_results?.columns || json.dataset) {
      this.parseModern(json);
    } else if (json.MLE?.content || json.fits) {
      this.parseLegacy(json);
    } else {
      throw new Error('Unrecognized HyPhy JSON format: missing site_results or MLE.');
    }
    this.runThreshold = runPThreshold(this.analysis.settings);

    this.postProcess();
  }

  parseModern(j) {
    this.analysis = j.analysis || {};
    this.provenance = j.provenance || {};

    // Dataset & Taxa
    const ds = j.dataset || {};
    this.taxa = ds.taxa || [];
    this.numTaxa = this.taxa.length;
    this.numSites = ds.codon_sites || 0;
    this.sequences = ds.sequences || {};

    // Tree & Phylogeny
    const phy = j.phylogeny || {};
    this.treeNewick = phy.newick || '';
    this.treeNodes = phy.nodes || {};

    // Check branch testing partitions
    for (const [name, node] of Object.entries(this.treeNodes)) {
      if (node.tested === true) this.numTested++;
      else if (node.tested === false && node.type === 'leaf') this.numBackground++;
      if (node.tested === false) this.hasBackground = true;
    }

    // Global omega
    const codonFit = j.model_fits?.codon_mg94 || {};
    this.globalOmega = codonFit.global_parameters?.omega ?? 1.0;

    // Site results (Columnar format)
    const sr = j.site_results || {};
    const cols = sr.columns || {};
    const data = sr.data || {};

    // Contrast-FEL groups
    this.contrastGroups = [];
    if (this.isContrast) {
      for (const colName of Object.keys(data)) {
        if (colName.startsWith('beta_') && colName !== 'beta_null' && colName !== 'beta_background') {
          this.contrastGroups.push(colName.substring(5));
        }
      }
    }

    const sitesCount = data.site?.length || this.numSites;
    this.siteResults = [];

    for (let i = 0; i < sitesCount; ++i) {
      const alpha = data.alpha?.[i] ?? 0;
      const beta = data.beta?.[i] ?? (data.beta_plus?.[i] ?? 0);
      const beta1 = data.beta1?.[i] ?? data.beta?.[i] ?? 0;
      const p1 = data.p1?.[i] ?? (this.isMEME ? 0.0 : 1.0);
      const beta_plus = data.beta_plus?.[i] ?? 0;
      const p_plus = data.p_plus?.[i] ?? 0;
      const branches_under_selection = data.branches_under_selection?.[i] ?? 0;
      const branch_ebf = data.branch_ebf?.[i] || {};
      const log_l_meme = data.log_l_meme?.[i] ?? null;
      const log_l_fel = data.log_l_fel?.[i] ?? null;

      const alpha_null = data.alpha_null?.[i] ?? 0;
      const lrt = data.penalized_lrt?.[i] ?? (data.lrt?.[i] ?? 0);
      const pval = data.p_value?.[i] ?? 1.0;
      const negLogP = pval > 0 ? -Math.log10(Math.max(pval, 1e-15)) : 15.0;
      const classification = data.classification?.[i] || (pval <= 0.05 ? 'different' : 'neutral');
      const history = data.substitution_history?.[i] || (data.substitutions?.[i] || '');
      const surprise = data.mean_substitution_surprise?.[i] ?? 0;
      const beta_bg = data.beta_background?.[i] ?? null;

      const groupBetas = {};
      if (this.isContrast) {
        for (const grp of this.contrastGroups) {
          groupBetas[grp] = data['beta_' + grp]?.[i] ?? 0;
        }
      }
      const contrastTag = data.contrast_tag?.[i] || '';
      const qVal = data.q_value?.[i] ?? 1.0;

      const omega = this.isMEME 
        ? (alpha > 1e-6 ? (beta_plus / alpha) : (beta_plus > 0 ? Infinity : 1.0))
        : (alpha > 1e-6 ? (beta / alpha) : (beta > 0 ? Infinity : 1.0));
      const spark = this.isMEME
        ? getRateSparkGlyphs(alpha, beta_plus)
        : getRateSparkGlyphs(alpha, beta);

      this.siteResults.push({
        site: data.site?.[i] ?? (i + 1),
        index: i,
        alpha,
        beta,
        beta1,
        p1,
        beta_plus,
        p_plus,
        branchesUnderSelection: branches_under_selection,
        branchEBF: branch_ebf,
        logL_meme: log_l_meme,
        logL_fel: log_l_fel,
        beta_bg,
        alpha_null,
        lrt,
        pValue: pval,
        qValue: qVal,
        contrastTag,
        groupBetas,
        negLogP,
        omega,
        classification,
        history,
        meanSurprise: surprise,
        spark,
        events: []
      });
    }

    // Substitution events
    const subEvents = j.substitution_events || [];
    this.substitutionEvents = subEvents;
    for (const item of subEvents) {
      const siteIdx = item.site;
      const events = item.events || [];
      this.eventsBySite[siteIdx] = events;

      if (this.siteResults[siteIdx - 1]) {
        this.siteResults[siteIdx - 1].events = events;
        this.siteResults[siteIdx - 1].rootCodon = (item.root_codon || '---').toUpperCase();
        this.siteResults[siteIdx - 1].rootAA = item.root_aa || translateCodon(item.root_codon);
      }

      for (const ev of events) {
        ev.site = siteIdx;
        if (!this.eventsByBranch[ev.branch]) {
          this.eventsByBranch[ev.branch] = [];
        }
        this.eventsByBranch[ev.branch].push(ev);
      }
    }
  }

  // Reconstruct full ancestral codon and amino acid sequences for all internal nodes
  reconstructAncestralSequences(treeRoot) {
    if (!treeRoot) return;

    // Helper: link parent references in tree
    function linkParents(node, parent = null) {
      node.parent = parent;
      for (const c of node.children || []) linkParents(c, node);
    }
    linkParents(treeRoot);

    const eventsBySite = {};
    for (const se of this.substitutionEvents || []) {
      eventsBySite[se.site] = se;
    }

    const self = this;
    function getAncestralState(node, site) {
      const se = eventsBySite[site];
      if (!se) {
        // 100% conserved site: identical to consensus
        const cAA = self.consensusAA[site - 1] || '-';
        const firstTaxon = Object.keys(self.sequences)[0];
        const codon = self.sequences[firstTaxon]?.substr((site - 1) * 3, 3) || '---';
        return { codon, aa: cAA };
      }
      // Check if substitution event occurred on branch to node
      const ev = (se.events || []).find(e => e.branch === node.name);
      if (ev) {
        return { codon: ev.child_codon, aa: ev.child_aa };
      }
      if (!node.parent) {
        return { codon: se.root_codon, aa: se.root_aa };
      }
      return getAncestralState(node.parent, site);
    }

    this.ancestralSequences = {};
    this.ancestralTranslated = {};
    this.ancestralTaxa = [];
    this.rootCodons = [];
    this.rootAA = [];

    // Precompute root sequence across all sites
    for (let s = 1; s <= this.numSites; ++s) {
      const st = getAncestralState(treeRoot, s);
      this.rootCodons.push(st.codon);
      this.rootAA.push(st.aa);
    }

    // Traverse all internal nodes to construct full sequences
    const traverse = (node) => {
      if (node.children && node.children.length > 0) {
        const name = node.name || 'Node';
        this.ancestralTaxa.push(name);
        const codonList = [];
        const aaList = [];
        for (let s = 1; s <= this.numSites; ++s) {
          const st = getAncestralState(node, s);
          codonList.push(st.codon);
          aaList.push(st.aa);
        }
        this.ancestralSequences[name] = codonList.join('');
        this.ancestralTranslated[name] = aaList;
        for (const c of node.children) traverse(c);
      }
    };
    traverse(treeRoot);
  }

  parseLegacy(j) {
    this.analysis = j.analysis || {};
    const input = j.input || {};
    this.numSites = input['number of sites'] || 0;
    this.numTaxa = input['number of sequences'] || 0;
    this.treeNewick = input.trees?.[0] || '';

    // Legacy site content
    const headers = orderedLegacyHeaders(j.MLE?.headers);
    let content = j.MLE?.content?.['by-site']?.['0'] || j.MLE?.content?.['0'] || [];
    if (!content.length && j.MLE?.content) {
      const firstKey = Object.keys(j.MLE.content)[0];
      if (firstKey && Array.isArray(j.MLE.content[firstKey])) {
        content = j.MLE.content[firstKey];
      }
    }

    const colIdx = {};
    headers.forEach((name, idx) => {
      const raw = name.toLowerCase();
      colIdx[raw] = idx;
      if (raw.includes('alpha') && !raw.includes('fel') && colIdx['alpha'] === undefined) colIdx['alpha'] = idx;
      if ((raw === 'lrt' || (raw.includes('lrt') && !raw.includes('fel'))) && colIdx['lrt'] === undefined) colIdx['lrt'] = idx;
      if (raw.includes('p-value') && colIdx['p-value'] === undefined) colIdx['p-value'] = idx;
      if (raw.includes('beta') && (raw.includes('1') || raw.includes('-')) && colIdx['beta1'] === undefined) colIdx['beta1'] = idx;
      if (raw.startsWith('p') && (raw.includes('1') || raw.includes('-')) && colIdx['p1'] === undefined) colIdx['p1'] = idx;
      if (raw.includes('beta') && (raw.includes('+') || raw.includes('2') || raw.includes('plus')) && colIdx['beta_plus'] === undefined) colIdx['beta_plus'] = idx;
      if (raw.startsWith('p') && (raw.includes('+') || raw.includes('2') || raw.includes('plus')) && colIdx['p_plus'] === undefined) colIdx['p_plus'] = idx;
      if ((raw.includes('branches under selection') || raw.includes('# branches')) && colIdx['branches'] === undefined) colIdx['branches'] = idx;
      if (raw.includes('meme logl')) colIdx['meme_logl'] = idx;
      if (raw.includes('fel logl')) colIdx['fel_logl'] = idx;
    });

    const aIdx = colIdx['alpha'] ?? 0;
    const bpIdx = colIdx['beta_plus'] ?? -1;
    const bIdx = this.isMEME ? (bpIdx >= 0 ? bpIdx : 1) : (colIdx['beta'] ?? 1);
    const b1Idx = colIdx['beta1'] ?? bIdx;
    const p1Idx = colIdx['p1'] ?? -1;
    const ppIdx = colIdx['p_plus'] ?? -1;
    const lrtIdx = colIdx['lrt'] ?? 3;
    const pIdx = colIdx['p-value'] ?? 4;
    const brIdx = colIdx['branches'] ?? -1;
    const mllIdx = colIdx['meme_logl'] ?? -1;
    const fllIdx = colIdx['fel_logl'] ?? -1;

    // Check for branch attributes posterior probabilities in legacy MEME
    const branchAttrs = j['branch attributes']?.['0'] || {};

    this.siteResults = content.map((row, i) => {
      const alpha = row[aIdx] ?? 0;
      const beta = this.isMEME ? (bpIdx >= 0 ? row[bpIdx] : 0) : (row[bIdx] ?? 0);
      const beta1 = b1Idx >= 0 ? row[b1Idx] : 0;
      const p1 = p1Idx >= 0 ? row[p1Idx] : 1.0;
      const beta_plus = bpIdx >= 0 ? row[bpIdx] : 0;
      const p_plus = ppIdx >= 0 ? row[ppIdx] : 0;
      const lrt = row[lrtIdx] ?? 0;
      const pval = row[pIdx] ?? 1.0;
      const negLogP = pval > 0 ? -Math.log10(Math.max(pval, 1e-15)) : 15.0;
      const branches = brIdx >= 0 ? row[brIdx] : 0;
      const logL_meme = mllIdx >= 0 ? row[mllIdx] : null;
      const logL_fel = fllIdx >= 0 ? row[fllIdx] : null;

      let classification = 'neutral';
      if (pval <= 0.1) {
        if (this.isMEME) {
          classification = 'diversifying';
        } else if (beta > alpha) {
          classification = 'diversifying';
        } else if (alpha > beta) {
          classification = 'purifying';
        }
      }

      // Compute branch EBF from legacy posterior array if present
      const branchEBF = {};
      if (this.isMEME && p_plus > 0 && p_plus < 1) {
        for (const [bName, bData] of Object.entries(branchAttrs)) {
          const postMatrix = bData['Posterior prob omega class by site'];
          if (postMatrix && postMatrix[1]) {
            const q = postMatrix[1][i];
            if (q !== undefined && q > 0) {
              if (q >= 0.9999) {
                branchEBF[bName] = 10000.0;
              } else if (q < 1.0) {
                const ebf = (q * (1 - p_plus)) / ((1 - q) * p_plus);
                if (ebf >= 1.0) branchEBF[bName] = ebf;
              }
            }
          }
        }
      }

      const omega = this.isMEME
        ? (alpha > 1e-6 ? (beta_plus / alpha) : (beta_plus > 0 ? Infinity : 1.0))
        : (alpha > 1e-6 ? (beta / alpha) : (beta > 0 ? Infinity : 1.0));
      const spark = this.isMEME
        ? getRateSparkGlyphs(alpha, beta_plus)
        : getRateSparkGlyphs(alpha, beta);

      const legSub = this.recoverLegacySubstitutions(j.substitutions?.['0'], i);

      return {
        site: i + 1,
        index: i,
        alpha,
        beta,
        beta1,
        p1,
        beta_plus,
        p_plus,
        branchesUnderSelection: branches,
        branchEBF,
        logL_meme,
        logL_fel,
        alpha_null: alpha,
        lrt,
        pValue: pval,
        negLogP,
        omega,
        classification,
        history: legSub.history,
        meanSurprise: 0,
        spark,
        events: legSub.events,
        rootCodon: legSub.rootCodon,
        rootAA: legSub.rootAA
      };
    });

    // Populate eventsBySite for legacy data
    for (const res of this.siteResults) {
      if (res.events && res.events.length > 0) {
        this.eventsBySite[res.site] = res.events;
      }
    }

    // Try extracting legacy global omega from fits if available
    const gFit = j.fits?.['Global MG94xREV'] || j.fits?.['Nucleotide GTR'];
    if (gFit) {
      const omegaDist = gFit['Rate Distributions']?.['non-synonymous/synonymous rate ratio'];
      if (omegaDist && omegaDist[0] && typeof omegaDist[0][0] === 'number') {
        this.globalOmega = omegaDist[0][0];
      }
    }
    if (this.globalOmega === 1.0 && this.siteResults.length > 0) {
      // Compute robust median/mean omega from site rates
      const validOmegas = this.siteResults.map(r => r.omega).filter(w => isFinite(w) && w > 0);
      if (validOmegas.length > 0) {
        validOmegas.sort((a, b) => a - b);
        this.globalOmega = validOmegas[Math.floor(validOmegas.length / 2)];
      }
    }

    this.numSites = this.siteResults.length;
  }

  // Recover substitutions from HyPhy 2.5 legacy substitutions map
  recoverLegacySubstitutions(substitutions0, siteIndex) {
    if (!substitutions0) return { history: '', events: [], rootCodon: '---', rootAA: '-' };
    const siteDict = substitutions0[String(siteIndex)] || substitutions0[siteIndex];
    if (!siteDict || typeof siteDict !== 'object') {
      return { history: '', events: [], rootCodon: '---', rootAA: '-' };
    }
    const root = (siteDict.root || '---').toUpperCase();
    const rootAA = translateCodon(root);
    const events = [];
    const changes = {};
    let synCount = 0;
    let multiHitCount = 0;

    for (const [node, rawCodon] of Object.entries(siteDict)) {
      if (node === 'root' || !rawCodon || rawCodon === '---') continue;
      const codon = rawCodon.toUpperCase();
      const childAA = translateCodon(codon);
      let dist = 0;
      for (let c = 0; c < 3; ++c) {
        if (root[c] && codon[c] && root[c] !== codon[c]) dist++;
      }
      const isSyn = (childAA === rootAA);
      if (dist > 1) multiHitCount++;
      if (isSyn) {
        synCount++;
      } else {
        const pair = dist > 1 ? `${rootAA} ↔ ${childAA} [${dist}nt]⚡` : `${rootAA} ↔ ${childAA}`;
        changes[pair] = (changes[pair] || 0) + 1;
      }
      events.push({
        branch: node,
        parent_codon: root,
        parent_aa: rootAA,
        child_codon: codon,
        child_aa: childAA,
        is_synonymous: isSyn,
        nt_distance: dist
      });
    }

    const parts = [];
    for (const [pair, count] of Object.entries(changes)) {
      parts.push(count > 1 ? `${pair} (${count}×)` : pair);
    }
    if (synCount > 0) parts.push(`(+${synCount} syn)`);
    const history = parts.length > 0 ? parts.join(', ') : (events.length > 0 ? 'Conserved' : '');
    return { history, events, rootCodon: root, rootAA, multiHitCount };
  }

  // Get dynamic selection statistics at any p-value cutoff
  getStats(threshold = 0.1) {
    let diversifying = 0;
    let purifying = 0;
    let neutral = 0;

    for (const r of this.siteResults) {
      if (r.pValue <= threshold) {
        if (this.isMEME) {
          diversifying++;
        } else if (r.beta > r.alpha) {
          diversifying++;
        } else if (r.alpha > r.beta) {
          purifying++;
        } else {
          neutral++;
        }
      } else {
        neutral++;
      }
    }

    return {
      diversifying,
      purifying,
      significant: diversifying + purifying,
      neutral,
      totalSites: this.siteResults.length,
      globalOmega: this.globalOmega
    };
  }

  // Get specific site result
  getSite(siteNumber) {
    return this.siteResults[siteNumber - 1] || null;
  }

  postProcess() {
    // 1. Translate nucleotide sequences to amino acids
    this.translated = {};
    for (const [taxon, seq] of Object.entries(this.sequences)) {
      const aaArr = [];
      const codonCount = Math.floor(seq.length / 3);
      for (let c = 0; c < codonCount; ++c) {
        const codon = seq.substr(c * 3, 3);
        aaArr.push(translateCodon(codon));
      }
      this.translated[taxon] = aaArr;
    }

    // 2. Compute consensus Codon & Amino Acid per site
    this.consensusCodon = [];
    this.consensusAA = [];
    const taxonList = Object.keys(this.translated);
    for (let c = 0; c < this.numSites; ++c) {
      const countsCodon = {};
      const countsAA = {};
      let maxCountCodon = 0, topCodon = '---';
      let maxCountAA = 0, topAA = '-';

      for (const t of taxonList) {
        const codon = this.sequences[t]?.substr(c * 3, 3).toUpperCase() || '---';
        const aa = this.translated[t]?.[c] || '-';

        if (codon && codon.length === 3 && codon !== '---') {
          countsCodon[codon] = (countsCodon[codon] || 0) + 1;
          if (countsCodon[codon] > maxCountCodon) {
            maxCountCodon = countsCodon[codon];
            topCodon = codon;
          }
        }

        if (aa !== '-' && aa !== '?') {
          countsAA[aa] = (countsAA[aa] || 0) + 1;
          if (countsAA[aa] > maxCountAA) {
            maxCountAA = countsAA[aa];
            topAA = aa;
          }
        }
      }
      this.consensusCodon.push(topCodon);
      this.consensusAA.push(topAA);
    }

    // 3. Ensure rootCodon and rootAA are populated for 100% of sites (including conserved sites)
    this.rootCodons = [];
    this.rootAA = [];
    for (let c = 0; c < this.numSites; ++c) {
      const r = this.siteResults[c];
      if (r) {
        if (!r.rootCodon || r.rootCodon === '---' || r.rootCodon === '') {
          if (this.consensusCodon[c] && this.consensusCodon[c] !== '---') {
            r.rootCodon = this.consensusCodon[c];
            r.rootAA = translateCodon(r.rootCodon);
          }
        }
        if (!r.rootAA || r.rootAA === '-' || r.rootAA === '?') {
          r.rootAA = translateCodon(r.rootCodon);
        }
        this.rootCodons.push(r.rootCodon);
        this.rootAA.push(r.rootAA);
      }
    }
  }

  // Load alignment from raw FASTA or NEXUS string if not embedded in JSON
  loadExternalAlignment(text) {
    const lines = text.split(/\r?\n/);
    const seqs = {};
    let currentTaxon = null;

    // Try FASTA first
    for (const line of lines) {
      const trimmed = line.trim();
      if (!trimmed) continue;
      if (trimmed.startsWith('>')) {
        currentTaxon = trimmed.substring(1).split(/\s+/)[0];
        seqs[currentTaxon] = '';
      } else if (currentTaxon) {
        seqs[currentTaxon] += trimmed.replace(/\s+/g, '').toUpperCase();
      }
    }

    if (Object.keys(seqs).length > 0) {
      this.sequences = seqs;
      this.taxa = Object.keys(seqs);
      this.numTaxa = this.taxa.length;
      this.postProcess();
      return true;
    }

    return false;
  }
}
