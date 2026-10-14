// Sticky multi-track selection header aligned with codon columns

import { SELECTION_COLORS } from '../data/genetic-code.js';
import { PALETTE, withAlpha } from '../data/palette.js';

export class TrackHeader {
  constructor(container, options = {}) {
    this.container = container;
    this.options = Object.assign({
      height: 105,
      aaColWidth: 24,
      codonColWidth: 46,
      displayMode: 'amino_acid',
      pThreshold: 0.1,
      bgColor: PALETTE.trackBg,
      onSiteSelect: null
    }, options);

    this.canvas = document.createElement('canvas');
    this.canvas.className = 'hyphy-track-header-canvas';
    this.ctx = this.canvas.getContext('2d');
    this.container.appendChild(this.canvas);

    this.data = null;
    this.scrollLeft = 0;
    this.selectedSite = 1;
    this.hoveredCol = -1;

    this.setupEvents();
  }

  setData(felData) {
    this.data = felData;
    this.render();
  }

  setDisplayMode(mode) {
    this.options.displayMode = mode;
    this.render();
  }

  setPThreshold(th) {
    this.options.pThreshold = th;
    this.render();
  }

  setSelectedSite(siteIndex) {
    this.selectedSite = siteIndex;
    this.render();
  }

  setScrollLeft(sl) {
    const totalW = this.data ? this.data.numSites * this.getColWidth() : 0;
    const maxSl = Math.max(0, totalW - (this.width || 0));
    this.scrollLeft = Math.max(0, Math.min(sl, maxSl));
    this.render();
  }

  getColWidth() {
    return this.options.displayMode === 'codon' 
      ? this.options.codonColWidth 
      : this.options.aaColWidth;
  }

  resize(width) {
    const dpr = window.devicePixelRatio || 1;
    this.width = width;
    this.height = this.options.height;
    this.canvas.width = width * dpr;
    this.canvas.height = this.height * dpr;
    this.canvas.style.width = width + 'px';
    this.canvas.style.height = this.height + 'px';
    this.ctx.setTransform(1, 0, 0, 1, 0, 0);
    this.ctx.scale(dpr, dpr);
    this.render();
  }

  render() {
    if (!this.data || !this.ctx || !this.width) return;
    const ctx = this.ctx;
    const w = this.width;
    const h = this.height;
    const colW = this.getColWidth();

    ctx.clearRect(0, 0, w, h);
    ctx.fillStyle = this.options.bgColor;
    ctx.fillRect(0, 0, w, h);

    const startCol = Math.max(0, Math.floor(this.scrollLeft / colW));
    const endCol = Math.min(this.data.numSites, Math.ceil((this.scrollLeft + w) / colW) + 1);

    ctx.save();
    ctx.translate(-this.scrollLeft, 0);

    // Selected Column Highlight (Braun Signal Orange)
    if (this.selectedSite >= 1 && this.selectedSite <= this.data.numSites) {
      const sx = (this.selectedSite - 1) * colW;
      ctx.fillStyle = withAlpha(PALETTE.accent, 0.08);
      ctx.fillRect(sx, 0, colW, h);
      ctx.strokeStyle = PALETTE.accent;
      ctx.lineWidth = 1.5;
      ctx.strokeRect(sx + 0.5, 0, colW - 1, h);
    }

    // Track 1: Site Position Index (Y: 0..20)
    ctx.font = '10px -apple-system, BlinkMacSystemFont, "SF Mono", Monaco, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    for (let c = startCol; c < endCol; ++c) {
      const siteNum = c + 1;
      const x = c * colW;
      const midX = x + colW / 2;

      // Draw tick & label for every 10th site, or every 5th if wide
      if (siteNum === 1 || siteNum % 10 === 0 || (colW >= 36 && siteNum % 5 === 0)) {
        ctx.fillStyle = PALETTE.ink;
        ctx.fillText(siteNum.toString(), midX, 10);
        ctx.strokeStyle = PALETTE.hairline;
        ctx.lineWidth = 1;
        ctx.beginPath();
        ctx.moveTo(midX, 16);
        ctx.lineTo(midX, 20);
        ctx.stroke();
      }
    }

    // Track 2: Selection Status Pills (Y: 22..40)
    for (let c = startCol; c < endCol; ++c) {
      const r = this.data.siteResults[c];
      if (!r) continue;
      const x = c * colW;
      const midX = x + colW / 2;

      const isSig = r.pValue <= this.options.pThreshold;
      if (isSig) {
        if (this.data.isMEME) {
          ctx.fillStyle = PALETTE.accent; // Braun Signal Orange
          const nBr = r.branchesUnderSelection || 0;
          if (colW >= 34 && nBr > 0) {
            ctx.font = 'bold 9px sans-serif';
            ctx.fillText(`▲${nBr}`, midX, 31);
          } else {
            ctx.font = 'bold 11px sans-serif';
            ctx.fillText('▲', midX, 31);
          }
        } else if (r.beta > r.alpha) {
          ctx.fillStyle = PALETTE.accent; // Braun Signal Orange
          ctx.font = 'bold 11px sans-serif';
          ctx.fillText('▲', midX, 31);
        } else if (r.alpha > r.beta) {
          ctx.fillStyle = PALETTE.purifying; // Technical Cobalt Blue
          ctx.font = 'bold 11px sans-serif';
          ctx.fillText('▼', midX, 31);
        }
      } else {
        ctx.fillStyle = PALETTE.dotMuted;
        ctx.fillText('·', midX, 31);
      }
    }

    // Track 3: Significance -log10(p) Bar Chart (Y: 42..72, height 30px)
    const track3Top = 42;
    const track3H = 30;
    const maxNegLogP = 6.0; // Scale up to p = 1e-6

    // Threshold baseline line (Precision dashed hairline in Braun Orange)
    const threshNegLog = -Math.log10(this.options.pThreshold);
    const threshY = track3Top + track3H - (threshNegLog / maxNegLogP) * track3H;
    ctx.strokeStyle = withAlpha(PALETTE.accent, 0.7);
    ctx.lineWidth = 1;
    ctx.setLineDash([2, 2]);
    ctx.beginPath();
    ctx.moveTo(this.scrollLeft, threshY);
    ctx.lineTo(this.scrollLeft + w, threshY);
    ctx.stroke();
    ctx.setLineDash([]);

    for (let c = startCol; c < endCol; ++c) {
      const r = this.data.siteResults[c];
      if (!r) continue;
      const x = c * colW;

      const val = Math.min(r.negLogP, maxNegLogP);
      const barH = Math.max(1, (val / maxNegLogP) * track3H);
      const barY = track3Top + track3H - barH;

      const isSig = r.pValue <= this.options.pThreshold;
      if (isSig) {
        if (this.data.isMEME) {
          ctx.fillStyle = PALETTE.accent; // All episodic selection peaks in Signal Orange
        } else {
          ctx.fillStyle = r.beta > r.alpha ? PALETTE.accent : PALETTE.purifying;
        }
      } else {
        ctx.fillStyle = PALETTE.hairline;
      }

      const barW = Math.max(3, colW - 6);
      ctx.fillRect(x + (colW - barW) / 2, barY, barW, barH);
    }

    // Track 4: Paradigm 2 Rate Sparkline / Balance Gauge (Y: 74..90)
    for (let c = startCol; c < endCol; ++c) {
      const r = this.data.siteResults[c];
      if (!r) continue;
      const x = c * colW;
      const midX = x + colW / 2;
      const isSig = r.pValue <= this.options.pThreshold;

      const spk = r.spark;
      if (colW >= 38) {
        // Draw textual glyphs
        ctx.font = '9px monospace';
        ctx.fillStyle = isSig && r.alpha > r.beta ? PALETTE.purifying : PALETTE.stone; // alpha
        ctx.fillText(spk.sa, midX - 8, 82);
        ctx.fillStyle = PALETTE.dotFaint;
        ctx.fillText(':', midX, 82);
        ctx.fillStyle = isSig && r.beta > r.alpha ? PALETTE.accent : PALETTE.stone; // beta
        ctx.fillText(spk.sb, midX + 8, 82);
      } else {
        // Draw proportional micro-segments
        const segW = Math.max(1.5, (colW - 8) / 7);
        let curX = x + 4;
        const alphaCol = isSig && r.alpha > r.beta ? PALETTE.purifying : PALETTE.stone;
        const betaCol = isSig && r.beta > r.alpha ? PALETTE.accent : PALETTE.stone;

        for (let i = 0; i < 3; ++i) {
          ctx.fillStyle = i < spk.na ? alphaCol : PALETTE.gridLine;
          ctx.fillRect(curX, 78, segW, 6);
          curX += segW + 0.5;
        }
        curX += 2;
        for (let i = 0; i < 3; ++i) {
          ctx.fillStyle = i < spk.nb ? betaCol : PALETTE.gridLine;
          ctx.fillRect(curX, 78, segW, 6);
          curX += segW + 0.5;
        }
      }
    }

    // Track 5: Shannon Perplexity / Surprise Flag ⚡ (Y: 92..104)
    for (let c = startCol; c < endCol; ++c) {
      const r = this.data.siteResults[c];
      if (!r) continue;
      const x = c * colW;
      const midX = x + colW / 2;

      let hasHighSurprise = false;
      if (r.events) {
        hasHighSurprise = r.events.some(ev => ev.surprise_bits >= 10.0);
      }

      if (hasHighSurprise) {
        ctx.fillStyle = PALETTE.amberOchre; // Amber Ochre
        ctx.font = '11px sans-serif';
        ctx.fillText('⚡', midX, 98);
      }
    }

    // Draw track baseline border
    ctx.strokeStyle = PALETTE.hairline;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(this.scrollLeft, h);
    ctx.lineTo(this.scrollLeft + w, h);
    ctx.stroke();

    ctx.restore();
  }

  setupEvents() {
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
  }
}
