// Detail Inspector Panel for selected codon site and ancestral substitution events

import { PALETTE, withAlpha } from '../data/palette.js';

export class Inspector {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      onBranchHighlight: null
    }, options);

    this.siteData = null;
    this.el = document.createElement('div');
    this.el.className = 'hyphy-inspector';
    this.container.appendChild(this.el);

    this.renderEmpty();
  }

  setSite(siteResult, consensusAA, isMEME = false) {
    this.siteData = siteResult;
    this.consensusAA = consensusAA;
    this.isMEME = isMEME;
    this.render();
  }

  renderEmpty() {
    this.el.innerHTML = `
      <div class="hyphy-inspector-empty">
        <span class="hyphy-inspector-icon">🔍</span>
        <span>Click any codon site or track above to inspect site-level rates and ancestral substitution events.</span>
      </div>
    `;
  }

  render() {
    if (!this.siteData) {
      this.renderEmpty();
      return;
    }

    const r = this.siteData;
    const isDiversifying = r.classification === 'diversifying';
    const isPurifying = r.classification === 'purifying';

    let badgeClass = 'badge-neutral';
    let badgeText = 'Neutral Evolution';
    const isContrast = !!r.contrastTag || (r.groupBetas && Object.keys(r.groupBetas).length > 0);

    if (isContrast) {
      if (r.pValue <= 0.05) {
        badgeClass = 'badge-diversifying';
        badgeText = r.contrastTag || 'Differential Selection';
      } else {
        badgeClass = 'badge-neutral';
        badgeText = 'Equal Selection Across Groups';
      }
    } else if (this.isMEME) {
      if (isDiversifying) {
        badgeClass = 'badge-diversifying';
        const nBr = r.branchesUnderSelection || 0;
        badgeText = nBr > 0 
          ? `▲ Episodic Diversifying (${nBr} branch${nBr > 1 ? 'es' : ''}, EBF ≥ 100)`
          : `▲ Episodic Diversifying Selection`;
      }
    } else {
      if (isDiversifying) {
        badgeClass = 'badge-diversifying';
        badgeText = '▲ Diversifying Selection';
      } else if (isPurifying) {
        badgeClass = 'badge-purifying';
        badgeText = '▼ Purifying Selection';
      }
    }

    const omegaStr = r.omega === Infinity ? 'Inf' : r.omega.toFixed(3);
    const pValStr = r.pValue < 0.0001 ? '<0.0001' : r.pValue.toFixed(4);

    // Contrast-FEL vs MEME vs FEL quick stats
    let statsHtml = '';
    if (isContrast) {
      let groupBetasHtml = '';
      for (const [grp, val] of Object.entries(r.groupBetas || {})) {
        groupBetasHtml += `<span class="stat-pill"><span class="stat-label">β (${grp}):</span> <strong>${val.toFixed(2)}</strong></span>`;
      }
      statsHtml = `
        <span class="stat-pill"><span class="stat-label">p-value:</span> <strong>${pValStr}</strong></span>
        <span class="stat-pill"><span class="stat-label">q-value (FDR):</span> <strong>${r.qValue != null ? (r.qValue < 0.0001 ? '<0.0001' : r.qValue.toFixed(4)) : 'N/A'}</strong></span>
        <span class="stat-pill"><span class="stat-label">PLRT:</span> <strong>${r.lrt.toFixed(2)}</strong></span>
        <span class="stat-pill"><span class="stat-label">α (Syn):</span> <strong>${r.alpha.toFixed(2)}</strong></span>
        ${groupBetasHtml}
      `;
    } else if (this.isMEME) {
      const deltaLogL = (r.logL_meme != null && r.logL_fel != null) 
        ? (r.logL_meme - r.logL_fel).toFixed(2) 
        : null;
      statsHtml = `
        <span class="stat-pill"><span class="stat-label">p-value:</span> <strong>${pValStr}</strong></span>
        <span class="stat-pill"><span class="stat-label">LRT:</span> <strong>${r.lrt.toFixed(2)}</strong></span>
        <span class="stat-pill"><span class="stat-label">α (Syn):</span> <strong>${r.alpha.toFixed(2)}</strong></span>
        <span class="stat-pill"><span class="stat-label">β1 (p1):</span> <strong>${r.beta1.toFixed(2)} (${(r.p1 * 100).toFixed(0)}%)</strong></span>
        <span class="stat-pill"><span class="stat-label">β+ (p+):</span> <strong style="color:${PALETTE.accent};">${r.beta_plus.toFixed(2)} (${(r.p_plus * 100).toFixed(0)}%)</strong></span>
        <span class="stat-pill"><span class="stat-label">Rates (α:β+):</span> <strong class="spark-code">${r.spark.text}</strong></span>
        ${deltaLogL ? `<span class="stat-pill"><span class="stat-label">ΔLogL (MEME–FEL):</span> <strong>+${deltaLogL}</strong></span>` : ''}
      `;
    } else {
      statsHtml = `
        <span class="stat-pill"><span class="stat-label">p-value:</span> <strong>${pValStr}</strong></span>
        <span class="stat-pill"><span class="stat-label">LRT:</span> <strong>${r.lrt.toFixed(2)}</strong></span>
        <span class="stat-pill"><span class="stat-label">dN/dS (ω):</span> <strong>${omegaStr}</strong></span>
        <span class="stat-pill"><span class="stat-label">Rates (α:β):</span> <strong class="spark-code">${r.alpha.toFixed(2)} : ${r.beta.toFixed(2)} [${r.spark.text}]</strong></span>
      `;
    }

    // Candidate Episodic Branches section (for MEME)
    let episodicBranchesHtml = '';
    if (this.isMEME) {
      const ebfEntries = Object.entries(r.branchEBF || {}).sort((a, b) => b[1] - a[1]);
      if (ebfEntries.length > 0) {
        const branchRows = ebfEntries.map(([bName, ebfVal]) => {
          const isHigh = ebfVal >= 100.0;
          const isMod = ebfVal >= 10.0 && ebfVal < 100.0;
          const badge = isHigh
            ? `<span class="event-badge badge-nonsyn">EBF ≥ 100 (Strong Selection)</span>`
            : (isMod 
                ? `<span class="event-badge" style="background-color:${withAlpha(PALETTE.amber, 0.12)};color:${PALETTE.amberOchre};border:1px solid ${withAlpha(PALETTE.amber, 0.3)};">EBF ≥ 10 (Moderate)</span>`
                : `<span class="event-badge badge-syn">EBF ≥ 1.0 (Weak)</span>`);
          const ebfFmt = ebfVal >= 1000 ? ebfVal.toFixed(0) : ebfVal.toFixed(1);
          return `
            <tr class="event-row" data-branch="${bName}">
              <td class="cell-branch"><strong>${bName}</strong></td>
              <td class="cell-prob"><strong style="${isHigh ? `color:${PALETTE.accent};` : ''}">${ebfFmt}</strong></td>
              <td class="cell-type">${badge}</td>
            </tr>
          `;
        }).join('');

        episodicBranchesHtml = `
          <div class="inspector-section" style="margin-bottom: 12px;">
            <div class="inspector-section-title">
              <span>Lineages Subject to Episodic Diversifying Selection (${ebfEntries.length})</span>
              <span class="title-hint">Hover row to highlight branch on tree</span>
            </div>
            <div class="events-table-wrapper">
              <table class="events-table">
                <thead>
                  <tr>
                    <th>Branch / Node</th>
                    <th>Empirical Bayes Factor (EBF)</th>
                    <th>Evidence Support</th>
                  </tr>
                </thead>
                <tbody>
                  ${branchRows}
                </tbody>
              </table>
            </div>
          </div>
        `;
      }
    }

    // Substitution events rows
    const events = r.events || [];
    let eventsHtml = '';

    if (events.length === 0) {
      eventsHtml = `
        <div class="events-empty">
          <span>No ancestral substitution events detected on this tree (Completely conserved site).</span>
        </div>
      `;
    } else {
      const rows = events.map(ev => {
        const isSyn = ev.is_synonymous;
        const typeBadge = isSyn 
          ? `<span class="event-badge badge-syn">Synonymous</span>` 
          : `<span class="event-badge badge-nonsyn">Non-synonymous</span>`;

        const hitBadge = ev.nt_distance > 1 
          ? `<span class="hit-badge">[${ev.nt_distance}nt hit]</span>` 
          : '';

        const surpriseHigh = ev.surprise_bits >= 10.0;
        const surpriseClass = surpriseHigh ? 'surprise-high' : '';
        const surpriseIcon = surpriseHigh ? '⚡ ' : '';

        const pProbStr = ev.transition_probability < 0.001 
          ? ev.transition_probability.toExponential(2) 
          : ev.transition_probability.toFixed(4);

        return `
          <tr class="event-row" data-branch="${ev.branch}">
            <td class="cell-branch"><strong>${ev.branch}</strong></td>
            <td class="cell-trans">${ev.parent_aa} ↔ ${ev.child_aa}</td>
            <td class="cell-codon"><code>${ev.parent_codon}</code> ↔ <code>${ev.child_codon}</code> ${hitBadge}</td>
            <td class="cell-type">${typeBadge}</td>
            <td class="cell-prob">${pProbStr}</td>
            <td class="cell-surprise ${surpriseClass}">${surpriseIcon}${ev.surprise_bits.toFixed(1)} bits</td>
            <td class="cell-tested">${ev.tested ? '<span class="tag-tested">Tested</span>' : '<span class="tag-bg">Background</span>'}</td>
          </tr>
        `;
      }).join('');

      eventsHtml = `
        <div class="events-table-wrapper">
          <table class="events-table">
            <thead>
              <tr>
                <th>Branch</th>
                <th>Residue</th>
                <th>Codon Transition</th>
                <th>Substitution Type</th>
                <th>Transition P</th>
                <th>Shannon Surprise</th>
                <th>Branch Set</th>
              </tr>
            </thead>
            <tbody>
              ${rows}
            </tbody>
          </table>
        </div>
      `;
    }

    this.el.innerHTML = `
      <div class="hyphy-inspector-header">
        <div class="inspector-title-group">
          <span class="inspector-site-tag">Site ${r.site}</span>
          <span class="inspector-consensus">Consensus: <strong>${this.consensusAA || '-'}</strong></span>
          <span class="inspector-badge ${badgeClass}">${badgeText}</span>
        </div>
        <div class="inspector-quick-stats">
          ${statsHtml}
        </div>
      </div>
      <div class="hyphy-inspector-body">
        ${episodicBranchesHtml}
        <div class="inspector-section-title">
          <span>Ancestral Substitution Events (${events.length})</span>
          <span class="title-hint">Hover over a branch row to highlight it on the phylogenetic tree</span>
        </div>
        ${eventsHtml}
      </div>
    `;

    // Setup hover interactions for table rows
    const tableRows = this.el.querySelectorAll('.event-row');
    tableRows.forEach(row => {
      const branch = row.getAttribute('data-branch');
      row.addEventListener('mouseenter', () => {
        if (this.options.onBranchHighlight) this.options.onBranchHighlight(branch);
      });
      row.addEventListener('mouseleave', () => {
        if (this.options.onBranchHighlight) this.options.onBranchHighlight(null);
      });
    });
  }
}
