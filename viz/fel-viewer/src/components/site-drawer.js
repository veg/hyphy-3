// Site Inspector Drawer: Progressive Disclosure of Granular Branch Events & Site Phylogeny Tree
// Slides open on demand from the right to reveal exact substitution events, rates, and dedicated site tree

import { SiteTree } from './site-tree.js';
import { translateCodon } from '../data/genetic-code.js';
import { PALETTE, withAlpha } from '../data/palette.js';
import { formatEBF, formatOmegaInf, formatPValue } from '../core/format.js';

export class SiteDrawer {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      threshold: 0.1,
      onClose: null
    }, options);

    this.threshold = this.options.threshold;
    this.data = null;
    this.currentSite = null;
    this.currentTab = 'summary'; // 'summary' | 'tree'
    this.isOpen = false;

    this.initDOM();
    this.initKeyboardListener();
  }

  initDOM() {
    // Backdrop
    this.backdrop = document.createElement('div');
    this.backdrop.className = 'hyphy-drawer-backdrop';
    this.backdrop.addEventListener('click', () => this.close());
    this.container.appendChild(this.backdrop);

    // Drawer container
    this.drawer = document.createElement('div');
    this.drawer.className = 'hyphy-site-drawer';
    this.drawer.innerHTML = `
      <div class="drawer-header">
        <div class="drawer-title-group">
          <span class="drawer-site-tag" id="drawer-site-tag">Codon #—</span>
          <span class="drawer-regime-badge" id="drawer-regime-badge">—</span>
        </div>
        <div class="segmented-control drawer-nav-tabs">
          <button type="button" class="seg-btn active" data-tab="summary" id="tab-btn-summary">Overview &amp; Events</button>
          <button type="button" class="seg-btn" data-tab="tree" id="tab-btn-tree">Site Tree 🌳</button>
        </div>
        <button type="button" class="drawer-close-btn" id="drawer-close-btn" title="Close Inspector (Esc)">✕</button>
      </div>

      <div class="drawer-body" id="drawer-body">
        <div id="drawer-summary-panel">
          <!-- Rates and events table -->
        </div>
        <div id="drawer-tree-panel" style="display: none;">
          <!-- SiteTree component mounted here -->
        </div>
      </div>
    `;

    this.container.appendChild(this.drawer);

    this.drawer.querySelector('#drawer-close-btn').addEventListener('click', () => this.close());
    this.siteTag = this.drawer.querySelector('#drawer-site-tag');
    this.regimeBadge = this.drawer.querySelector('#drawer-regime-badge');
    this.summaryPanel = this.drawer.querySelector('#drawer-summary-panel');
    this.treePanel = this.drawer.querySelector('#drawer-tree-panel');

    // Tab buttons
    this.tabBtns = this.drawer.querySelectorAll('.drawer-nav-tabs .seg-btn');
    this.tabBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        this.setTab(btn.dataset.tab);
      });
    });

    // Initialize SiteTree inside tree panel
    this.siteTree = new SiteTree(this.treePanel);
  }

  setTab(tab) {
    this.currentTab = tab;
    this.tabBtns.forEach(b => {
      if (b.dataset.tab === tab) b.classList.add('active');
      else b.classList.remove('active');
    });

    if (tab === 'summary') {
      this.summaryPanel.style.display = 'block';
      this.treePanel.style.display = 'none';
    } else {
      this.summaryPanel.style.display = 'none';
      this.treePanel.style.display = 'block';
      if (this.currentSite && this.data) {
        this.siteTree.setSite(this.currentSite, this.data.getSite(this.currentSite));
      }
    }
  }

  initKeyboardListener() {
    this.onKeyDown = (e) => {
      if (e.key === 'Escape' && this.isOpen) {
        this.close();
      }
    };
    window.addEventListener('keydown', this.onKeyDown);
  }

  // Detach from the page when the viewer swaps in another method's components
  destroy() {
    window.removeEventListener('keydown', this.onKeyDown);
    this.backdrop.remove();
    this.drawer.remove();
  }

  setData(data) {
    this.data = data;
    if (this.siteTree && data) {
      this.siteTree.setTree(data.treeNewick);
      this.siteTree.setSequences(data.sequences);
    }
  }

  setThreshold(p) {
    this.threshold = p;
    if (this.isOpen && this.currentSite) {
      this.inspectSite(this.currentSite, this.currentTab);
    }
  }

  open() {
    this.isOpen = true;
    this.backdrop.classList.add('active');
    this.drawer.classList.add('active');
  }

  close() {
    this.isOpen = false;
    this.backdrop.classList.remove('active');
    this.drawer.classList.remove('active');
    if (typeof this.options.onClose === 'function') {
      this.options.onClose();
    }
  }

  inspectSite(siteNumber, activeTab = null) {
    if (!this.data) return;
    const r = this.data.getSite(siteNumber);
    if (!r) return;

    this.currentSite = siteNumber;
    if (activeTab) {
      this.setTab(activeTab);
    }
    this.open();

    const isSig = r.pValue <= this.threshold;
    const isDiv = this.data.isMEME ? isSig : (isSig && r.beta > r.alpha);
    const isPur = !this.data.isMEME && isSig && r.alpha > r.beta;

    this.siteTag.textContent = `Codon #${r.site}`;
    if (isDiv) {
      this.regimeBadge.className = 'drawer-regime-badge tag-div';
      this.regimeBadge.textContent = this.data.isMEME ? '▲ Episodic Diversifying Selection' : '▲ Diversifying Selection';
    } else if (isPur) {
      this.regimeBadge.className = 'drawer-regime-badge tag-pur';
      this.regimeBadge.textContent = '▼ Purifying Selection';
    } else {
      this.regimeBadge.className = 'drawer-regime-badge tag-neutral';
      this.regimeBadge.textContent = this.data.isMEME ? 'Neutral / Constrained' : 'Neutral Evolution';
    }

    const omegaStr = formatOmegaInf(r.omega);
    const pStr = formatPValue(r.pValue, 3);

    // Get substitution events
    const events = r.events || this.data.eventsBySite[r.site] || [];
    const rootCodon = r.rootCodon || (this.data.rootCodons?.[r.site - 1]) || '---';
    const rootAA = r.rootAA || (this.data.rootAA?.[r.site - 1]) || translateCodon(rootCodon);

    // Update SiteTree
    if (this.siteTree) {
      this.siteTree.setSite(siteNumber, r);
    }

    // Extant codon counts
    const extantCodons = {};
    const extantAAs = {};
    if (this.data.sequences && Object.keys(this.data.sequences).length > 0) {
      for (const [taxon, seq] of Object.entries(this.data.sequences)) {
        const codon = seq.substr((r.site - 1) * 3, 3).toUpperCase();
        if (codon && codon.length === 3 && codon !== '---') {
          extantCodons[codon] = (extantCodons[codon] || 0) + 1;
        }
      }
      for (const [taxon, aaArr] of Object.entries(this.data.translated)) {
        const aa = aaArr[r.site - 1];
        if (aa && aa !== '-' && aa !== '?') {
          extantAAs[aa] = (extantAAs[aa] || 0) + 1;
        }
      }
    }

    let ratesGridHtml = '';
    if (this.data.isMEME) {
      const betaPlus = (r.beta_plus != null) ? r.beta_plus : r.beta;
      const pPlusPct = Math.round((r.p_plus ?? 0) * 100);
      const p1Pct = Math.round((r.p1 ?? (1 - (r.p_plus ?? 0))) * 100);
      const brCount = r.branchesUnderSelection || (r.branchEBF ? Object.values(r.branchEBF).filter(v => v >= 100).length : 0);

      ratesGridHtml = `
        <div class="drawer-section">
          <div class="drawer-section-title">Maximum Likelihood Rates (MEME Mixture Model)</div>
          <div class="drawer-stats-grid">
            <div class="drawer-stat-card">
              <div class="dstat-label">Synonymous Rate (α)</div>
              <div class="dstat-value font-mono font-bold">${r.alpha.toFixed(2)}</div>
              <div class="dstat-sub">dS estimate</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Constrained Rate (β₁)</div>
              <div class="dstat-value font-mono font-bold">${(r.beta1 ?? 0).toFixed(2)}</div>
              <div class="dstat-sub">p₁ = ${p1Pct}% of branches</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Episodic Rate (β⁺)</div>
              <div class="dstat-value font-mono font-bold text-orange">${betaPlus.toFixed(2)}</div>
              <div class="dstat-sub">p⁺ = ${pPlusPct}% under selection</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Ratio (β⁺ / α)</div>
              <div class="dstat-value font-mono font-bold">${omegaStr}</div>
              <div class="dstat-sub">${isDiv ? 'Episodic Diversification' : 'Unconstrained'}</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">LRT Statistic</div>
              <div class="dstat-value font-mono font-bold">${r.lrt.toFixed(2)}</div>
              <div class="dstat-sub">2ΔlnL (vs β⁺ ≤ α)</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">p-value</div>
              <div class="dstat-value font-mono font-bold ${isSig ? 'text-orange' : ''}">${pStr}</div>
              <div class="dstat-sub">${isSig ? 'Episodically Diversifying' : 'Not Significant'}</div>
            </div>
          </div>
        </div>
      `;
    } else {
      ratesGridHtml = `
        <div class="drawer-section">
          <div class="drawer-section-title">Maximum Likelihood Rates</div>
          <div class="drawer-stats-grid">
            <div class="drawer-stat-card">
              <div class="dstat-label">Synonymous Rate (α)</div>
              <div class="dstat-value font-mono font-bold">${r.alpha.toFixed(2)}</div>
              <div class="dstat-sub">dS estimate</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Non-synonymous Rate (β)</div>
              <div class="dstat-value font-mono font-bold">${r.beta.toFixed(2)}</div>
              <div class="dstat-sub">dN estimate</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Ratio (ω = dN/dS)</div>
              <div class="dstat-value font-mono font-bold">${omegaStr}</div>
              <div class="dstat-sub">${isDiv ? 'Positive Selection' : (isPur ? 'Constraint' : 'Neutral')}</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">LRT Statistic</div>
              <div class="dstat-value font-mono font-bold">${r.lrt.toFixed(2)}</div>
              <div class="dstat-sub">2ΔlnL (vs neutral)</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">p-value</div>
              <div class="dstat-value font-mono font-bold ${isSig ? 'text-orange' : ''}">${pStr}</div>
              <div class="dstat-sub">${isSig ? 'Statistically Significant' : 'Not Significant'}</div>
            </div>
            <div class="drawer-stat-card">
              <div class="dstat-label">Neutral Rate (α=β)</div>
              <div class="dstat-value font-mono font-bold">${(r.alpha_null ?? 0).toFixed(2)}</div>
              <div class="dstat-sub">Null model rate</div>
            </div>
          </div>
        </div>
      `;
    }

    let ebfSectionHtml = '';
    if (this.data.isMEME) {
      const supportedBranches = Object.entries(r.branchEBF || {})
        .filter(([_, ebf]) => ebf >= 100)
        .sort((a, b) => b[1] - a[1]);

      if (supportedBranches.length > 0) {
        ebfSectionHtml = `
          <div class="drawer-section">
            <div class="drawer-section-title">
              <span>Branches Under Episodic Selection (EBF ≥ 100)</span>
              <span class="root-badge" style="background: ${withAlpha(PALETTE.accent, 0.15)}; color: ${PALETTE.accent}; border-color: ${withAlpha(PALETTE.accent, 0.3)};">
                <strong>${supportedBranches.length}</strong> branch(es) supported
              </span>
            </div>
            <div class="events-table-wrapper">
              <table class="events-table">
                <thead>
                  <tr>
                    <th>Branch</th>
                    <th>Empirical Bayes Factor (EBF)</th>
                    <th>Reconstructed Substitution</th>
                    <th>Type</th>
                  </tr>
                </thead>
                <tbody>
                  ${supportedBranches.map(([bName, ebf]) => {
                    const ev = events.find(e => e.branch === bName);
                    const subStr = ev
                      ? `<strong>${ev.parent_codon || rootCodon}</strong> (${ev.parent_aa || rootAA}) → <strong>${ev.child_codon}</strong> (${ev.child_aa})`
                      : `<span style="color: ${PALETTE.neutral}; font-style: italic;">No discrete mutation reconstructed</span>`;
                    const typeHtml = ev
                      ? `<span class="ev-tag ${ev.is_synonymous ? 'ev-syn' : 'ev-nonsyn'}">${ev.is_synonymous ? 'Synonymous' : 'Non-synonymous'}</span>`
                      : '—';
                    return `
                      <tr>
                        <td class="font-mono font-bold">${bName}</td>
                        <td class="font-mono"><span class="ebf-branch-badge ebf-active">⚡ EBF ${formatEBF(ebf)}</span></td>
                        <td class="font-mono">${subStr}</td>
                        <td>${typeHtml}</td>
                      </tr>
                    `;
                  }).join('')}
                </tbody>
              </table>
            </div>
          </div>
        `;
      } else {
        ebfSectionHtml = `
          <div class="drawer-section">
            <div class="drawer-section-title">Branches Under Episodic Selection</div>
            <div class="empty-events-notice">
              No individual branches reached the strict EBF ≥ 100 evidentiary cutoff at this codon site.
              ${r.p_plus > 0 ? `The episodic signal is distributed across branches (overall p⁺ = ${(r.p_plus * 100).toFixed(1)}%).` : ''}
            </div>
          </div>
        `;
      }
    }

    let html = `
      <!-- 1. Rates Grid -->
      ${ratesGridHtml}

      <!-- 2. Branches Under Episodic Selection (MEME Only) -->
      ${ebfSectionHtml}

      <!-- 3. Reconstructed Substitutions -->
      <div class="drawer-section">
        <div class="drawer-section-title">
          <span>Reconstructed Substitutions</span>
          <span class="root-badge">Root: <strong>${rootCodon}</strong> (${rootAA})</span>
        </div>
        <div class="sub-history-callout">
          <strong>Summary:</strong> ${r.history || 'No substitutions detected (fully conserved site)'}
        </div>
    `;

    if (events.length > 0) {
      html += `
        <div class="events-table-wrapper">
          <table class="events-table">
            <thead>
              <tr>
                <th>Branch</th>
                <th>Transition</th>
                <th>Type</th>
                <th>Distance</th>
                <th style="text-align: right;">Surprise</th>
              </tr>
            </thead>
            <tbody>
      `;

      for (const ev of events) {
        const isMultiHit = (ev.nt_distance > 1);
        const typeClass = ev.is_synonymous ? 'ev-syn' : 'ev-nonsyn';
        const typeLabel = ev.is_synonymous ? 'Synonymous' : 'Non-synonymous';
        const distLabel = `${ev.nt_distance || 1}nt${isMultiHit ? ' ⚡' : ''}`;
        const surpriseStr = ev.surprise_bits != null ? `${ev.surprise_bits.toFixed(1)} bits` : '—';
        const ebf = r.branchEBF?.[ev.branch];
        const ebfTag = (ebf != null && ebf >= 100)
          ? `<span class="ebf-branch-badge ebf-active" style="margin-left: 6px; font-size: 8.5px; padding: 1px 4px;">⚡ EBF ${formatEBF(ebf)}</span>`
          : '';

        html += `
          <tr>
            <td class="font-mono">${ev.branch}${ebfTag}</td>
            <td class="font-mono">
              <strong>${ev.parent_codon || rootCodon}</strong> (${ev.parent_aa || rootAA}) → 
              <strong>${ev.child_codon}</strong> (${ev.child_aa})
            </td>
            <td><span class="ev-tag ${typeClass}">${typeLabel}</span></td>
            <td class="font-mono ${isMultiHit ? 'text-orange font-bold' : ''}">${distLabel}</td>
            <td style="text-align: right;" class="font-mono">${surpriseStr}</td>
          </tr>
        `;
      }

      html += `
            </tbody>
          </table>
        </div>
      `;
    } else {
      html += `
        <div class="empty-events-notice">
          No branch-specific substitution events reconstructed for this codon site.
        </div>
      `;
    }

    html += `</div>`;

    // 3. Extant Residues (if present)
    if (Object.keys(extantCodons).length > 0) {
      html += `
        <div class="drawer-section">
          <div class="drawer-section-title">Observed Extant Residues</div>
          <div class="extant-chips-grid">
      `;
      for (const [codon, count] of Object.entries(extantCodons)) {
        html += `
          <div class="extant-chip">
            <span class="extant-codon font-mono">${codon}</span>
            <span class="extant-count font-mono">×${count}</span>
          </div>
        `;
      }
      html += `
          </div>
        </div>
      `;
    }

    this.summaryPanel.innerHTML = html;
  }
}
