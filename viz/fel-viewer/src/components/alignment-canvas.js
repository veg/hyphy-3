// Virtualized 2D canvas sequence alignment matrix renderer (Codons & Amino Acids)

import { ZAPPO_COLORS, NUCLEOTIDE_COLORS, SELECTION_COLORS } from '../data/genetic-code.js';
import { PALETTE, withAlpha } from '../data/palette.js';

export class AlignmentCanvas {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      rowHeight: 22,
      aaColWidth: 24,
      codonColWidth: 46,
      font: '12px -apple-system, BlinkMacSystemFont, "SF Mono", Monaco, "Cascadia Code", Consolas, monospace',
      bgColor: PALETTE.panel,
      cellBorderColor: PALETTE.gridLine,
      displayMode: 'amino_acid',     // 'amino_acid' | 'codon'
      colorScheme: 'difference',      // 'difference' | 'selection' | 'chemistry'
      onSiteSelect: null,
      onCellHover: null,
      onScroll: null
    }, options);

    this.canvas = document.createElement('canvas');
    this.canvas.className = 'hyphy-alignment-canvas';
    this.ctx = this.canvas.getContext('2d');
    this.container.appendChild(this.canvas);

    this.data = null;
    this.taxaOrder = [];            // Taxa in vertical tree display order
    this.scrollLeft = 0;
    this.scrollTop = 0;
    this.selectedSite = 1;          // 1-based site index
    this.hoveredCol = -1;
    this.hoveredRow = -1;

    this.setupEvents();
  }

  setData(felData, taxaOrder = []) {
    this.data = felData;
    this.taxaOrder = taxaOrder.length > 0 ? taxaOrder : felData.taxa;
    this.render();
  }

  setTaxaOrder(taxaOrder) {
    this.taxaOrder = taxaOrder;
    this.render();
  }

  setDisplayMode(mode) {
    this.options.displayMode = mode;
    this.render();
  }

  setColorScheme(scheme) {
    this.options.colorScheme = scheme;
    this.render();
  }

  setSelectedSite(siteIndex) {
    this.selectedSite = siteIndex;
    this.render();
  }

  getColWidth() {
    return this.options.displayMode === 'codon' 
      ? this.options.codonColWidth 
      : this.options.aaColWidth;
  }

  getTotalWidth() {
    if (!this.data) return 0;
    return this.data.numSites * this.getColWidth();
  }

  getTotalHeight() {
    if (!this.data) return 0;
    return this.taxaOrder.length * this.options.rowHeight;
  }

  getVisibleRange() {
    const colW = this.getColWidth();
    const rowH = this.options.rowHeight;
    const startCol = Math.max(0, Math.floor(this.scrollLeft / colW));
    const endCol = Math.min(this.data?.numSites || 0, Math.ceil((this.scrollLeft + (this.width || 0)) / colW));
    const startRow = Math.max(0, Math.floor(this.scrollTop / rowH));
    const endRow = Math.min(this.taxaOrder.length, Math.ceil((this.scrollTop + (this.height || 0)) / rowH));
    return { startCol, endCol, startRow, endRow };
  }

  setScroll(sl, st) {
    const maxSl = Math.max(0, this.getTotalWidth() - (this.width || 0));
    const maxSt = Math.max(0, this.getTotalHeight() - (this.height || 0));
    this.scrollLeft = Math.max(0, Math.min(sl, maxSl));
    this.scrollTop = Math.max(0, Math.min(st, maxSt));
    this.render();
    if (this.options.onScroll) {
      this.options.onScroll(this.scrollLeft, this.scrollTop, this.getVisibleRange());
    }
  }

  resize(width, height) {
    const dpr = window.devicePixelRatio || 1;
    this.width = width;
    this.height = height;
    this.canvas.width = width * dpr;
    this.canvas.height = height * dpr;
    this.canvas.style.width = width + 'px';
    this.canvas.style.height = height + 'px';
    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    this.ctx.scale(dpr, dpr);
    this.render();
    if (this.options.onScroll) {
      this.options.onScroll(this.scrollLeft, this.scrollTop, this.getVisibleRange());
    }
  }

  render() {
    if (!this.data || !this.ctx || !this.width) return;
    const ctx = this.ctx;
    const w = this.width;
    const h = this.height;
    const colW = this.getColWidth();
    const rowH = this.options.rowHeight;
    const isCodon = this.options.displayMode === 'codon';
    const scheme = this.options.colorScheme;

    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = this.options.bgColor;
    ctx.fillRect(0, 0, w, h);

    // Compute visible window
    const startCol = Math.max(0, Math.floor(this.scrollLeft / colW));
    const endCol = Math.min(this.data.numSites, Math.ceil((this.scrollLeft + w) / colW) + 1);
    const startRow = Math.max(0, Math.floor(this.scrollTop / rowH));
    const endRow = Math.min(this.taxaOrder.length, Math.ceil((this.scrollTop + h) / rowH) + 1);

    ctx.save();
    ctx.translate(-this.scrollLeft, -this.scrollTop);

    // 1. Draw Column Selection Highlight Bar (Iconic Braun Signal Orange precision focus)
    if (this.selectedSite >= 1 && this.selectedSite <= this.data.numSites) {
      const sx = (this.selectedSite - 1) * colW;
      ctx.fillStyle = withAlpha(PALETTE.accent, 0.08);
      ctx.fillRect(sx, this.scrollTop, colW, h);
      ctx.strokeStyle = PALETTE.accent;
      ctx.lineWidth = 1.5;
      ctx.strokeRect(sx + 0.5, this.scrollTop, colW - 1, h);
    }

    // 2. Draw Visible Cells
    ctx.font = this.options.font;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    const lightZappoResidues = new Set(['A', 'I', 'L', 'M', 'V']);

    for (let r = startRow; r < endRow; ++r) {
      const taxon = this.taxaOrder[r];
      const isAncestor = this.data.ancestralSequences?.[taxon] !== undefined;
      const rawSeq = isAncestor 
        ? (this.data.ancestralSequences[taxon] || '') 
        : (this.data.sequences[taxon] || '');
      const aaSeq = isAncestor 
        ? (this.data.ancestralTranslated[taxon] || []) 
        : (this.data.translated[taxon] || []);
      const y = r * rowH;

      // Ancestral row subtle background wash
      if (isAncestor) {
        ctx.fillStyle = PALETTE.cellBg;
        ctx.fillRect(this.scrollLeft, y, w, rowH);
      }

      // Draw subtle row hairline separator
      ctx.strokeStyle = this.options.cellBorderColor;
      ctx.lineWidth = 0.5;
      ctx.beginPath();
      ctx.moveTo(this.scrollLeft, y + rowH);
      ctx.lineTo(this.scrollLeft + w, y + rowH);
      ctx.stroke();

      for (let c = startCol; c < endCol; ++c) {
        const x = c * colW;
        const siteResult = this.data.siteResults[c];
        const refAA = scheme === 'difference_root'
          ? (this.data.rootAA?.[c] || this.data.consensusAA[c] || '-')
          : (this.data.consensusAA[c] || '-');
        const aa = aaSeq[c] || '-';
        const codon = rawSeq.substr(c * 3, 3) || '---';

        const isMutated = aa !== refAA && aa !== '-';

        // Background cell fill
        let cellFill = null;
        let textColor = isAncestor ? PALETTE.ancestorInk : PALETTE.ink;

        if (scheme === 'difference' || scheme === 'difference_root') {
          if (isMutated) {
            cellFill = ZAPPO_COLORS[aa] || PALETTE.stone;
            textColor = lightZappoResidues.has(aa) ? PALETTE.ink : PALETTE.knockout;
          } else {
            textColor = isAncestor ? PALETTE.dotMuted : PALETTE.dotFaint; // Restrained Dieter Rams muted dot
          }
        } else if (scheme === 'chemistry') {
          cellFill = ZAPPO_COLORS[aa] || PALETTE.gridLine;
          textColor = lightZappoResidues.has(aa) ? PALETTE.ink : PALETTE.knockout;
        } else if (scheme === 'selection') {
          if (siteResult) {
            if (siteResult.classification === 'diversifying') {
              cellFill = isMutated ? PALETTE.accent : withAlpha(PALETTE.accent, 0.12);
              textColor = isMutated ? PALETTE.knockout : PALETTE.accent;
            } else if (siteResult.classification === 'purifying') {
              cellFill = isMutated ? PALETTE.purifying : withAlpha(PALETTE.purifying, 0.10);
              textColor = isMutated ? PALETTE.knockout : PALETTE.purifying;
            } else {
              if (isMutated) {
                cellFill = PALETTE.gridLine;
                textColor = PALETTE.ink;
              } else {
                textColor = PALETTE.dotFaint;
              }
            }
          }
        }

        // Draw machined tactile tile
        if (cellFill) {
          ctx.fillStyle = cellFill;
          ctx.beginPath();
          ctx.roundRect(x + 1.5, y + 2, colW - 3, rowH - 4, 2);
          ctx.fill();
        }

        // Text character rendering
        ctx.fillStyle = textColor;
        if (isCodon) {
          if (scheme === 'difference' && !isMutated && aa !== '-') {
            ctx.fillStyle = PALETTE.branchGrey;
          }
          ctx.fillText(codon, x + colW / 2, y + rowH / 2);
        } else {
          // If difference mode and not mutated, display a subtle dot
          if (scheme === 'difference' && !isMutated && aa !== '-') {
            ctx.fillText('·', x + colW / 2, y + rowH / 2);
          } else {
            ctx.fillText(aa, x + colW / 2, y + rowH / 2);
          }
        }
      }
    }

    ctx.restore();
  }

  setupEvents() {
    this.canvas.addEventListener('mousemove', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left + this.scrollLeft;
      const my = e.clientY - rect.top + this.scrollTop;
      const colW = this.getColWidth();
      const rowH = this.options.rowHeight;

      const col = Math.floor(mx / colW);
      const row = Math.floor(my / rowH);

      if (col >= 0 && col < this.data?.numSites && row >= 0 && row < this.taxaOrder.length) {
        if (this.hoveredCol !== col || this.hoveredRow !== row) {
          this.hoveredCol = col;
          this.hoveredRow = row;
          if (this.options.onCellHover) {
            const taxon = this.taxaOrder[row];
            const siteResult = this.data.siteResults[col];
            const rawSeq = this.data.sequences[taxon] || '';
            const codon = rawSeq.substr(col * 3, 3);
            const aa = this.data.translated[taxon]?.[col] || '-';
            this.options.onCellHover({
              site: col + 1,
              taxon,
              codon,
              aa,
              siteResult,
              clientX: e.clientX,
              clientY: e.clientY
            });
          }
        }
      } else if (this.hoveredCol !== -1) {
        this.hoveredCol = -1;
        this.hoveredRow = -1;
        if (this.options.onCellHover) this.options.onCellHover(null);
      }
    });

    this.canvas.addEventListener('mouseleave', () => {
      this.hoveredCol = -1;
      this.hoveredRow = -1;
      if (this.options.onCellHover) this.options.onCellHover(null);
    });

    this.canvas.addEventListener('click', (e) => {
      const rect = this.canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left + this.scrollLeft;
      const colW = this.getColWidth();
      const site = Math.floor(mx / colW) + 1;

      if (site >= 1 && site <= this.data?.numSites) {
        this.setSelectedSite(site);
        if (this.options.onSiteSelect) {
          this.options.onSiteSelect(site, this.data.siteResults[site - 1]);
        }
      }
    });

    // Native mouse wheel virtual scroll
    this.canvas.addEventListener('wheel', (e) => {
      e.preventDefault();
      const maxScrollLeft = Math.max(0, this.getTotalWidth() - this.width);
      const maxScrollTop = Math.max(0, this.getTotalHeight() - this.height);

      let newSl = this.scrollLeft + e.deltaX;
      let newSt = this.scrollTop + e.deltaY;

      // Handle Shift+Wheel for horizontal scroll
      if (e.shiftKey && e.deltaX === 0) {
        newSl = this.scrollLeft + e.deltaY;
        newSt = this.scrollTop;
      }

      newSl = Math.max(0, Math.min(newSl, maxScrollLeft));
      newSt = Math.max(0, Math.min(newSt, maxScrollTop));

      this.scrollLeft = newSl;
      this.scrollTop = newSt;
      this.render();

      if (this.options.onScroll) {
        this.options.onScroll(this.scrollLeft, this.scrollTop, this.getVisibleRange());
      }
    }, { passive: false });
  }
}
