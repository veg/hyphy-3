// Candidate Sites Table: Diversifying Only Default
// Dieter Rams minimalist data table with category pills, rate sparkbars, live search, and CSV export

import { PALETTE } from '../data/palette.js';
import { formatPValue, formatThreshold, formatOmegaInf } from '../core/format.js';
import { downloadCSV, csvQuote } from '../core/export.js';

export class CandidatesTable {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      threshold: 0.1,
      filterType: 'diversifying', // Diversifying only by default!
      onSiteSelect: null
    }, options);

    this.threshold = this.options.threshold;
    this.filterType = this.options.filterType;
    this.data = null;
    this.selectedSite = null;
    this.searchQuery = '';

    this.initDOM();
  }

  initDOM() {
    this.wrapper = document.createElement('div');
    this.wrapper.className = 'hyphy-candidates-card';
    this.wrapper.innerHTML = `
      <div class="candidates-header">
        <div class="candidates-title-group">
          <span class="candidates-title">Candidate Sites</span>
          <span class="candidates-subtitle" id="candidates-subtitle">Filtered by selective regime</span>
        </div>
        <div class="candidates-controls">
          <!-- Category Filter Pills: Diversifying Default -->
          <div class="segmented-control category-filter-control">
            <button type="button" class="seg-btn active" data-filter="diversifying" id="btn-filter-div">
              ▲ Diversifying <span class="badge-count" id="count-div">0</span>
            </button>
            <button type="button" class="seg-btn" data-filter="purifying" id="btn-filter-pur">
              ▼ Purifying <span class="badge-count" id="count-pur">0</span>
            </button>
            <button type="button" class="seg-btn" data-filter="significant" id="btn-filter-sig">
              All Significant <span class="badge-count" id="count-sig">0</span>
            </button>
            <button type="button" class="seg-btn" data-filter="all" id="btn-filter-all">
              All Sites <span class="badge-count" id="count-all">0</span>
            </button>
          </div>

          <!-- Search Input -->
          <div class="candidates-search-box">
            <input type="text" class="candidates-search-input" id="search-input" placeholder="Search codon or AA..." />
          </div>

          <!-- CSV Export -->
          <button type="button" class="plot-tool-btn" id="btn-export-csv" title="Export Candidate Sites to CSV">Export CSV</button>
        </div>
      </div>

      <!-- Table Body Container -->
      <div class="candidates-table-scroll" id="table-scroll">
        <table class="candidates-table" id="candidates-table">
          <thead>
            <tr>
              <th style="width: 70px;">Codon</th>
              <th style="width: 130px;">Rates (α : β)</th>
              <th style="width: 140px;">Balance (dS vs dN)</th>
              <th style="width: 90px;">ω (dN/dS)</th>
              <th style="width: 80px;">LRT</th>
              <th style="width: 90px;">p-value</th>
              <th>Substitutions Reconstructed</th>
              <th style="width: 140px; text-align: right;">Inspect / Tree</th>
            </tr>
          </thead>
          <tbody id="candidates-tbody">
            <!-- Dynamic rows -->
          </tbody>
        </table>
        <div class="candidates-empty-state" id="empty-state" style="display: none;">
          <!-- Dynamic empty message -->
        </div>
      </div>
    `;

    this.container.appendChild(this.wrapper);

    this.tbody = this.wrapper.querySelector('#candidates-tbody');
    this.thead = this.wrapper.querySelector('#candidates-table thead');
    this.emptyState = this.wrapper.querySelector('#empty-state');
    this.tableScroll = this.wrapper.querySelector('#table-scroll');
    this.tableEl = this.wrapper.querySelector('#candidates-table');
    this.subtitle = this.wrapper.querySelector('#candidates-subtitle');
    this.btnDiv = this.wrapper.querySelector('#btn-filter-div');
    this.btnPur = this.wrapper.querySelector('#btn-filter-pur');
    this.btnSig = this.wrapper.querySelector('#btn-filter-sig');
    this.btnAll = this.wrapper.querySelector('#btn-filter-all');

    // Filter Buttons
    const filterBtns = this.wrapper.querySelectorAll('.category-filter-control .seg-btn');
    filterBtns.forEach(btn => {
      btn.addEventListener('click', () => {
        filterBtns.forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        this.setFilterType(btn.dataset.filter);
      });
    });

    // Search input
    const searchInput = this.wrapper.querySelector('#search-input');
    searchInput.addEventListener('input', (e) => {
      this.searchQuery = e.target.value;
      this.render();
    });

    // Export CSV
    this.wrapper.querySelector('#btn-export-csv').addEventListener('click', () => this.exportCSV());
  }

  setData(data) {
    this.data = data;
    this.updateCounts();
    this.render();
  }

  setThreshold(p) {
    this.threshold = p;
    this.updateCounts();
    this.render();
  }

  setFilterType(filterType) {
    this.filterType = filterType;
    this.render();
  }

  selectSite(siteNumber) {
    this.selectedSite = siteNumber;
    // Update active row highlighting
    const rows = this.tbody.querySelectorAll('tr');
    rows.forEach(r => {
      if (parseInt(r.dataset.site, 10) === siteNumber) {
        r.classList.add('row-selected');
        r.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      } else {
        r.classList.remove('row-selected');
      }
    });
  }

  updateCounts() {
    if (!this.data) return;
    const stats = this.data.getStats(this.threshold);
    this.wrapper.querySelector('#count-div').textContent = stats.diversifying;
    this.wrapper.querySelector('#count-pur').textContent = stats.purifying;
    this.wrapper.querySelector('#count-sig').textContent = stats.significant;
    this.wrapper.querySelector('#count-all').textContent = stats.totalSites;

    if (this.data.isMEME) {
      this.btnPur.style.display = 'none';
      this.btnSig.style.display = 'none';
      this.btnDiv.childNodes[0].nodeValue = '▲ Episodic Diversifying ';
      if (this.filterType === 'purifying' || this.filterType === 'significant') {
        this.filterType = 'diversifying';
        const filterBtns = this.wrapper.querySelectorAll('.category-filter-control .seg-btn');
        filterBtns.forEach(b => b.classList.remove('active'));
        this.btnDiv.classList.add('active');
      }
    } else {
      this.btnPur.style.display = '';
      this.btnSig.style.display = '';
      this.btnDiv.childNodes[0].nodeValue = '▲ Diversifying ';
    }
  }

  updateHeader() {
    if (this.data && this.data.isMEME) {
      this.thead.innerHTML = `
        <tr>
          <th style="width: 70px;">Codon</th>
          <th style="width: 140px;">Rates (α : β⁺)</th>
          <th style="width: 140px;">Weight (p⁺)</th>
          <th style="width: 110px;">Branches (EBF≥100)</th>
          <th style="width: 80px;">LRT</th>
          <th style="width: 90px;">p-value</th>
          <th>Substitutions Reconstructed</th>
          <th style="width: 140px; text-align: right;">Inspect / Tree</th>
        </tr>
      `;
    } else {
      this.thead.innerHTML = `
        <tr>
          <th style="width: 70px;">Codon</th>
          <th style="width: 130px;">Rates (α : β)</th>
          <th style="width: 140px;">Balance (dS vs dN)</th>
          <th style="width: 90px;">ω (dN/dS)</th>
          <th style="width: 80px;">LRT</th>
          <th style="width: 90px;">p-value</th>
          <th>Substitutions Reconstructed</th>
          <th style="width: 140px; text-align: right;">Inspect / Tree</th>
        </tr>
      `;
    }
  }

  render() {
    if (!this.data) return;
    this.updateCounts();
    this.updateHeader();

    const filtered = this.data.getFilteredSites(this.threshold, this.filterType, this.searchQuery);
    const stats = this.data.getStats(this.threshold);
    const threshStr = formatThreshold(this.threshold);

    // Update Subtitle
    if (this.data.isMEME) {
      if (this.filterType === 'diversifying') {
        this.subtitle.textContent = `Showing ${filtered.length} episodically diversifying sites (p ≤ ${threshStr})`;
      } else {
        this.subtitle.textContent = `Showing all ${filtered.length} codon sites`;
      }
    } else {
      if (this.filterType === 'diversifying') {
        this.subtitle.textContent = `Showing ${filtered.length} diversifying sites (β > α, p ≤ ${threshStr})`;
      } else if (this.filterType === 'purifying') {
        this.subtitle.textContent = `Showing ${filtered.length} purifying sites (α > β, p ≤ ${threshStr})`;
      } else if (this.filterType === 'significant') {
        this.subtitle.textContent = `Showing ${filtered.length} significant sites (p ≤ ${threshStr})`;
      } else {
        this.subtitle.textContent = `Showing all ${filtered.length} codon sites`;
      }
    }

    if (filtered.length === 0) {
      this.tableEl.style.display = 'none';
      this.emptyState.style.display = 'block';

      if (this.data.isMEME) {
        this.emptyState.innerHTML = `
          <div class="empty-icon">▲</div>
          <div class="empty-title">No Episodically Diversifying Sites at p ≤ ${threshStr}</div>
          <div class="empty-desc">
            No codon sites showed evidence of episodic diversifying selection at the current threshold.<br/>
            Adjust the significance cutoff slider in the KPI bar to explore higher p-values.
          </div>
        `;
      } else if (this.filterType === 'diversifying') {
        this.emptyState.innerHTML = `
          <div class="empty-icon">▲</div>
          <div class="empty-title">No Diversifying Sites at p ≤ ${threshStr}</div>
          <div class="empty-desc">
            No codon sites showed evidence of positive selection (β > α) at the current threshold.<br/>
            ${stats.purifying > 0 ? `<strong>${stats.purifying} purifying sites</strong> were detected.` : ''}
          </div>
          ${stats.purifying > 0 ? `<button type="button" class="btn-jump" id="btn-switch-pur">View Purifying Sites (▼ ${stats.purifying})</button>` : ''}
        `;
        const switchBtn = this.emptyState.querySelector('#btn-switch-pur');
        if (switchBtn) {
          switchBtn.addEventListener('click', () => {
            const purBtn = this.wrapper.querySelector('#btn-filter-pur');
            purBtn.click();
          });
        }
      } else {
        this.emptyState.innerHTML = `
          <div class="empty-icon">🔍</div>
          <div class="empty-title">No Matching Sites</div>
          <div class="empty-desc">No codon sites match your current filter and search query.</div>
        `;
      }
      return;
    }

    this.tableEl.style.display = 'table';
    this.emptyState.style.display = 'none';

    // Render Rows
    let html = '';
    for (const r of filtered) {
      const isSelected = (this.selectedSite === r.site);
      const isSig = r.pValue <= this.threshold;
      const isDiv = this.data.isMEME ? isSig : (isSig && r.beta > r.alpha);
      const isPur = !this.data.isMEME && isSig && r.alpha > r.beta;

      const pStr = formatPValue(r.pValue);
      const historyStr = r.history || (r.events?.length ? `${r.events.length} events` : '—');

      if (this.data.isMEME) {
        const betaPlus = (r.beta_plus != null) ? r.beta_plus : r.beta;
        const pPlusPct = Math.round((r.p_plus ?? 0) * 100);
        const brCount = r.branchesUnderSelection || (r.branchEBF ? Object.values(r.branchEBF).filter(v => v >= 100).length : 0);
        const brBadge = brCount > 0
          ? `<span class="ebf-branch-badge ebf-active" title="${brCount} branch(es) with Empirical Bayes Factor ≥ 100">⚡ ${brCount}</span>`
          : `<span class="font-mono" style="color: ${PALETTE.neutral}; font-size: 11px;">0</span>`;

        html += `
          <tr data-site="${r.site}" class="${isSelected ? 'row-selected' : ''}">
            <td class="font-mono font-bold site-cell">
              <span class="site-number">#${r.site}</span>
            </td>
            <td class="font-mono rates-cell">
              <span class="rate-alpha font-mono" title="Synonymous rate (α)">α=${r.alpha.toFixed(2)}</span>
              <span class="rate-sep">|</span>
              <span class="rate-beta font-mono" title="Episodic non-synonymous rate (β⁺)">β⁺=${betaPlus.toFixed(2)}</span>
            </td>
            <td class="balance-cell">
              <div class="meme-weight-bar-wrapper" title="Weight allocated to episodic selection: ${pPlusPct}%">
                <div class="micro-balance-bar">
                  <div class="bar-meme-weight" style="width: ${pPlusPct}%;"></div>
                </div>
                <span class="meme-weight-label">${pPlusPct}%</span>
              </div>
            </td>
            <td class="branches-cell">
              ${brBadge}
            </td>
            <td class="font-mono lrt-cell">${r.lrt.toFixed(2)}</td>
            <td class="font-mono pval-cell ${isSig ? 'pval-sig font-bold' : ''}" style="${isSig ? `color: ${PALETTE.accent};` : ''}">${pStr}</td>
            <td class="sub-cell" title="${historyStr}">
              <div class="sub-text-wrapper">${historyStr}</div>
            </td>
            <td style="text-align: right; white-space: nowrap;">
              <button type="button" class="btn-inspect" data-action="summary" data-site="${r.site}">Inspect</button>
              <button type="button" class="btn-tree" data-action="tree" data-site="${r.site}" title="View Site Phylogeny Tree">🌳 Tree</button>
            </td>
          </tr>
        `;
      } else {
        // Micro balance bar for FEL
        const maxRate = Math.max(0.1, r.alpha + r.beta);
        const alphaPct = Math.round((r.alpha / maxRate) * 100);
        const betaPct = Math.round((r.beta / maxRate) * 100);
        const omegaStr = formatOmegaInf(r.omega);

        html += `
          <tr data-site="${r.site}" class="${isSelected ? 'row-selected' : ''}">
            <td class="font-mono font-bold site-cell">
              <span class="site-number">#${r.site}</span>
            </td>
            <td class="font-mono rates-cell">
              <span class="rate-alpha font-mono" title="Synonymous rate (α)">α=${r.alpha.toFixed(2)}</span>
              <span class="rate-sep">|</span>
              <span class="rate-beta font-mono" title="Non-synonymous rate (β)">β=${r.beta.toFixed(2)}</span>
            </td>
            <td class="balance-cell">
              <div class="micro-balance-bar" title="dS: ${r.alpha.toFixed(2)} vs dN: ${r.beta.toFixed(2)}">
                <div class="bar-alpha" style="width: ${alphaPct}%;"></div>
                <div class="bar-beta" style="width: ${betaPct}%;"></div>
              </div>
            </td>
            <td class="omega-cell">
              <span class="omega-badge ${isDiv ? 'omega-div' : (isPur ? 'omega-pur' : '')}">${omegaStr}</span>
            </td>
            <td class="font-mono lrt-cell">${r.lrt.toFixed(2)}</td>
            <td class="font-mono pval-cell ${isSig ? 'pval-sig font-bold' : ''}">${pStr}</td>
            <td class="sub-cell" title="${historyStr}">
              <div class="sub-text-wrapper">${historyStr}</div>
            </td>
            <td style="text-align: right; white-space: nowrap;">
              <button type="button" class="btn-inspect" data-action="summary" data-site="${r.site}">Inspect</button>
              <button type="button" class="btn-tree" data-action="tree" data-site="${r.site}" title="View Site Phylogeny Tree">🌳 Tree</button>
            </td>
          </tr>
        `;
      }
    }

    this.tbody.innerHTML = html;

    // Row click listeners
    const rows = this.tbody.querySelectorAll('tr');
    rows.forEach(row => {
      row.addEventListener('click', (e) => {
        const site = parseInt(row.dataset.site, 10);
        this.selectedSite = site;
        rows.forEach(r => r.classList.remove('row-selected'));
        row.classList.add('row-selected');

        const btn = e.target.closest('button');
        const action = btn?.dataset?.action || 'summary';

        if (typeof this.options.onSiteSelect === 'function') {
          this.options.onSiteSelect(site, this.data.getSite(site), action);
        }
      });
    });
  }

  exportCSV() {
    if (!this.data) return;
    const filtered = this.data.getFilteredSites(this.threshold, this.filterType, this.searchQuery);
    const isMEME = this.data.isMEME;

    const headers = isMEME
      ? ['Codon', 'Alpha (dS)', 'Beta+ (dN+)', 'Weight p+', 'Branches (EBF>=100)', 'LRT', 'p-value', 'Regime', 'Substitutions']
      : ['Codon', 'Alpha (dS)', 'Beta (dN)', 'dN/dS (Omega)', 'LRT', 'p-value', 'Regime', 'Substitutions'];

    const rows = filtered.map(r => {
      const isSig = r.pValue <= this.threshold;
      if (isMEME) {
        const betaPlus = r.beta_plus ?? r.beta;
        const pPlus = r.p_plus ?? 0;
        const brCount = r.branchesUnderSelection || (r.branchEBF ? Object.values(r.branchEBF).filter(v => v >= 100).length : 0);
        return [
          r.site,
          r.alpha.toFixed(4),
          betaPlus.toFixed(4),
          pPlus.toFixed(4),
          brCount,
          r.lrt.toFixed(4),
          r.pValue,
          isSig ? 'Episodic Diversifying' : 'Neutral',
          csvQuote(r.history)
        ];
      } else {
        const isDiv = isSig && r.beta > r.alpha;
        const isPur = isSig && r.alpha > r.beta;
        const regime = isDiv ? 'Diversifying' : (isPur ? 'Purifying' : 'Neutral');
        return [
          r.site,
          r.alpha.toFixed(4),
          r.beta.toFixed(4),
          isFinite(r.omega) ? r.omega.toFixed(4) : 'Infinity',
          r.lrt.toFixed(4),
          r.pValue,
          regime,
          csvQuote(r.history)
        ];
      }
    });

    downloadCSV(headers, rows, `hyphy_${isMEME ? 'meme' : 'fel'}_candidates_${this.filterType}.csv`);
  }
}
