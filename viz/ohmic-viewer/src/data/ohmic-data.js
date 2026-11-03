/**
 * Data parser and model normalizer for HyPhy OHMIC (v3 modern & legacy JSON formats).
 * Includes Newick parser with Wheatstone bridge hub/node mapping.
 */

export function parseNewick(newickStr) {
  if (!newickStr || typeof newickStr !== 'string') {
    return { root: null, nodes: [], leaves: [], nodeByName: new Map() };
  }
  const s = newickStr.trim().replace(/;$/, '');
  let pos = 0;
  let nextId = 0;
  const nodes = [];
  const leaves = [];
  const nodeByName = new Map();

  function parseSubtree(parent = null) {
    const node = {
      id: nextId++,
      name: '',
      branchLength: 0,
      parent,
      children: [],
      isLeaf: false,
      depth: parent ? parent.depth + 1 : 0,
      x: 0,
      y: 0
    };
    nodes.push(node);

    if (s[pos] === '(') {
      pos++; // consume '('
      while (pos < s.length) {
        const child = parseSubtree(node);
        node.children.push(child);
        while (pos < s.length && /\s/.test(s[pos])) pos++;
        if (s[pos] === ',') {
          pos++;
        } else if (s[pos] === ')') {
          pos++;
          break;
        } else {
          break;
        }
      }
    } else {
      node.isLeaf = true;
    }

    // Parse node name and optional :branchLength
    let label = '';
    while (pos < s.length && !/[,():;]/.test(s[pos])) {
      label += s[pos++];
    }
    label = label.trim();
    if (label) {
      node.name = label;
    } else {
      node.name = node.isLeaf ? `Taxon_${node.id}` : `Node_${node.id}`;
    }

    if (s[pos] === ':') {
      pos++;
      let lenStr = '';
      while (pos < s.length && !/[,();]/.test(s[pos])) {
        lenStr += s[pos++];
      }
      const val = parseFloat(lenStr.trim());
      node.branchLength = Number.isFinite(val) ? Math.max(0, val) : 0;
    }

    if (node.isLeaf) {
      leaves.push(node);
    }
    if (node.name) {
      nodeByName.set(node.name, node);
    }
    return node;
  }

  const root = parseSubtree(null);

  // Assign leaf clade sets for each node
  function collectLeaves(n) {
    if (n.isLeaf) {
      n.leafSet = new Set([n.name]);
      return n.leafSet;
    }
    const set = new Set();
    for (const c of n.children) {
      for (const leafName of collectLeaves(c)) {
        set.add(leafName);
      }
    }
    n.leafSet = set;
    return set;
  }
  if (root) collectLeaves(root);

  return { root, nodes, leaves, nodeByName };
}

export class OhmicDataModel {
  constructor(rawJson) {
    this.raw = rawJson || {};
    this._parse();
  }

  _parse() {
    const raw = this.raw;
    const isModern = Boolean(raw.electrical_network || raw.wheatstone_bridges || raw.dataset);

    // 1. Dataset metadata
    const ds = raw.dataset || raw.input || {};
    this.fileName = ds.file || ds['file name'] || 'alignment.fna';
    this.datatype = ds.datatype || 'codon';
    this.taxaCount = ds.taxa_count ?? ds['number of sequences'] ?? 0;
    this.sitesCount = ds.sites_count ?? ds['number of sites'] ?? 0;
    this.patternsCount = ds.patterns_count ?? this.sitesCount;
    this.superNodesCount = ds.super_nodes_count ?? this.taxaCount;
    this.taxa = Array.isArray(ds.taxa) ? ds.taxa : [];

    // 2. Electrical Network summary
    const en = raw.electrical_network || raw['electrical network'] || {};
    this.solver = en.solver || 'ExactDense (Full Laplacian Pseudoinverse)';
    this.usedHyphaeon = Boolean(en.used_hyphaeon_embeddings);
    this.embeddingDim = en.embedding_dim || (this.usedHyphaeon ? 384 : 0);
    this.fosterSum = Number(en.foster_sum ?? (this.superNodesCount - 1));
    this.expectedFosterSum = Math.max(0, this.superNodesCount - 1);
    this.fosterResidual = Math.abs(this.fosterSum - this.expectedFosterSum);
    this.numBridges = Number(en.num_bridges ?? 0);
    this.numCollapsedPolytomies = Number(en.num_collapsed_polytomies ?? 0);
    this.backboneNewick = en.backbone_tree || '';
    this.collapsedNewick = en.collapsed_tree || this.backboneNewick;

    // Parse trees
    this.backboneTree = parseNewick(this.backboneNewick);
    this.collapsedTree = parseNewick(this.collapsedNewick);

    if (this.taxa.length === 0 && this.backboneTree.leaves.length > 0) {
      this.taxa = this.backboneTree.leaves.map((l) => l.name);
    }

    // Pairwise matrices (N x N)
    this.pairwiseOmega = Array.isArray(en.pairwise_effective_resistances)
      ? en.pairwise_effective_resistances
      : null;
    this.pairwiseFoster = Array.isArray(en.pairwise_foster_probabilities)
      ? en.pairwise_foster_probabilities
      : null;
    this.pairwiseDistances = Array.isArray(en.pairwise_tn93_distances)
      ? en.pairwise_tn93_distances
      : null;

    // Compute pairwise conductance c_ij = P_ij / Omega_ij if available
    if (this.pairwiseFoster && this.pairwiseOmega) {
      const N = this.pairwiseFoster.length;
      this.pairwiseConductance = Array.from({ length: N }, (_, i) =>
        Array.from({ length: N }, (__, j) => {
          if (i === j) return 0;
          const om = this.pairwiseOmega[i]?.[j] || 0;
          const fp = this.pairwiseFoster[i]?.[j] || 0;
          return om > 1e-12 ? fp / om : 0;
        })
      );
    } else {
      this.pairwiseConductance = null;
    }

    // 3. Model fits
    const fits = raw.fits || {};
    if (isModern && fits.single_tree && fits.ohmic_network) {
      this.singleTreeLogL = Number(fits.single_tree.log_likelihood ?? 0);
      this.singleTreeAICc = Number(fits.single_tree.aicc ?? 0);
      this.ohmicLogL = Number(fits.ohmic_network.log_likelihood ?? 0);
      this.ohmicAICc = Number(fits.ohmic_network.aicc ?? 0);
      this.deltaLogL = Number(fits.improvement?.delta_log_l ?? (this.ohmicLogL - this.singleTreeLogL));
      this.deltaAICc = Number(fits.improvement?.delta_aicc ?? (this.singleTreeAICc - this.ohmicAICc));
      this.mg94Params = fits.mg94_parameters || null;
      this.gtrParams = fits.gtr_parameters || null;
    } else {
      const st = fits['Single Tree'] || {};
      const on = fits['Ohmic Network'] || {};
      this.singleTreeLogL = Number(st['Log Likelihood'] ?? 0);
      this.singleTreeAICc = Number(st['AIC-c'] ?? 0);
      this.ohmicLogL = Number(on['Log Likelihood'] ?? 0);
      this.ohmicAICc = Number(on['AIC-c'] ?? 0);
      this.deltaLogL = this.ohmicLogL - this.singleTreeLogL;
      this.deltaAICc = this.singleTreeAICc - this.ohmicAICc;
      this.mg94Params = null;
      this.gtrParams = null;
    }

    // 4. Wheatstone bridges
    const rawBridges = raw.wheatstone_bridges || raw['wheatstone bridges'] || [];
    this.bridges = rawBridges.map((b, idx) => {
      const children = Array.isArray(b.children) ? b.children : ['A', 'B', 'C'];
      const initW = Array.isArray(b.initial_weights) ? b.initial_weights : [0.7, 0.3, 0.0];
      const fitW = Array.isArray(b.fitted_weights) ? b.fitted_weights : initW;
      const fitL = Array.isArray(b.fitted_branch_lengths) ? b.fitted_branch_lengths : [0.01, 0.01, 0.0];

      // Determine quartet leaf sets [A, B, C, D]
      let quartetClades = Array.isArray(b.quartet_clades) && b.quartet_clades.length === 4
        ? b.quartet_clades
        : null;

      if (!quartetClades) {
        const getLeaves = (nodeName) => {
          const n = this.backboneTree.nodeByName.get(nodeName) || this.collapsedTree.nodeByName.get(nodeName);
          return n && n.leafSet ? Array.from(n.leafSet) : [nodeName];
        };
        const cladeA = getLeaves(children[0]);
        const cladeB = getLeaves(children[1]);
        const cladeC = getLeaves(children[2]);
        const used = new Set([...cladeA, ...cladeB, ...cladeC]);
        const cladeD = this.taxa.filter((t) => !used.has(t));
        quartetClades = [cladeA, cladeB, cladeC, cladeD];
      }

      return {
        index: b.bridge_index ?? idx,
        hubNode: b.hub_node || `Hub_${idx}`,
        bridgeNode: b.bridge_node || `Bridge_${idx}`,
        children,
        quartetClades,
        shuntStress: Number(b.shunt_stress ?? 0),
        ohmicSupport: Number(b.ohmic_support ?? initW[0] ?? 0),
        initialWeights: initW.map(Number),
        fittedWeights: fitW.map(Number),
        fittedBranchLengths: fitL.map(Number),
        activeShuntWeight: Number(fitW[1] ?? 0) + Number(fitW[2] ?? 0),
        initialShuntWeight: Number(initW[1] ?? 0) + Number(initW[2] ?? 0)
      };
    });

    // Map bridges by hubNode and bridgeNode name for quick tree overlay lookup
    this.bridgeByHubName = new Map();
    this.bridgeByBridgeNodeName = new Map();
    for (const br of this.bridges) {
      this.bridgeByHubName.set(br.hubNode, br);
      this.bridgeByBridgeNodeName.set(br.bridgeNode, br);
    }

    // 5. Site-level log-likelihood profile
    const sr = raw.site_results || {};
    const deltaArr = Array.isArray(sr.delta_log_l) ? sr.delta_log_l : [];
    const singleArr = Array.isArray(sr.single_tree_log_l) ? sr.single_tree_log_l : [];
    const ohmicArr = Array.isArray(sr.ohmic_log_l) ? sr.ohmic_log_l : [];

    this.siteResults = deltaArr.map((d, i) => ({
      site: i + 1,
      deltaLogL: Number(d),
      singleTreeLogL: Number(singleArr[i] ?? 0),
      ohmicLogL: Number(ohmicArr[i] ?? 0)
    }));

    this.networkFavoringSitesCount = this.siteResults.filter((s) => s.deltaLogL > 0.05).length;
    this.treeFavoringSitesCount = this.siteResults.filter((s) => s.deltaLogL < -0.05).length;
    this.maxSiteDeltaLogL = this.siteResults.reduce(
      (max, s) => Math.max(max, Math.abs(s.deltaLogL)),
      0.01
    );

    this.runtimeSeconds = Number(raw.timers?.total_seconds ?? 0);
  }
}
