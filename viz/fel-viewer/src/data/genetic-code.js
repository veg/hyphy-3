// Genetic code tables, biophysical chemistry palettes, and selection colour schemes

import { PALETTE, withAlpha } from './palette.js';

export const CODON_TABLE = {
  // Phenylalanine
  'TTT': 'F', 'TTC': 'F',
  // Leucine
  'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
  // Isoleucine
  'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
  // Methionine (Start)
  'ATG': 'M',
  // Valine
  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
  // Serine
  'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
  // Proline
  'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
  // Threonine
  'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
  // Alanine
  'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
  // Tyrosine
  'TAT': 'Y', 'TAC': 'Y',
  // Stop
  'TAA': '*', 'TAG': '*', 'TGA': '*',
  // Histidine
  'CAT': 'H', 'CAC': 'H',
  // Glutamine
  'CAA': 'Q', 'CAG': 'Q',
  // Asparagine
  'AAT': 'N', 'AAC': 'N',
  // Lysine
  'AAA': 'K', 'AAG': 'K',
  // Aspartate
  'GAT': 'D', 'GAC': 'D',
  // Glutamate
  'GAA': 'E', 'GAG': 'E',
  // Cysteine
  'TGT': 'C', 'TGC': 'C',
  // Tryptophan
  'TGG': 'W',
  // Arginine
  'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
  // Glycine
  'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
};

const IUPAC = {
  'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
  'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
  'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T']
};

export function translateCodon(codon) {
  if (!codon || codon.length !== 3) return '-';
  const clean = codon.toUpperCase().replace(/U/g, 'T');
  if (CODON_TABLE[clean]) return CODON_TABLE[clean];

  // Resolve IUPAC ambiguities
  const pos = [
    IUPAC[clean[0]] || [clean[0]],
    IUPAC[clean[1]] || [clean[1]],
    IUPAC[clean[2]] || [clean[2]]
  ];
  const possibleAAs = new Set();
  for (const b0 of pos[0]) {
    for (const b1 of pos[1]) {
      for (const b2 of pos[2]) {
        const tr = CODON_TABLE[b0 + b1 + b2];
        if (tr) possibleAAs.add(tr);
      }
    }
  }
  if (possibleAAs.size === 1) return Array.from(possibleAAs)[0];
  if (possibleAAs.size > 1) return Array.from(possibleAAs).join('/');
  return '?';
}

// Calibrated Matte Zappo Palette (Dieter Rams aesthetic: desaturated, tactile industrial tones)
export const ZAPPO_COLORS = {
  // Aliphatic / hydrophobic (matte stone chalk)
  'A': '#E2DDD9', 'I': '#D6CFCA', 'L': '#D6CFCA', 'M': '#D6CFCA', 'V': '#D6CFCA',
  // Aromatic (warm raw sienna)
  'F': '#D97706', 'W': '#B45309', 'Y': '#D97706',
  // Positive / basic (Braun technical cobalt)
  'K': '#1E40AF', 'R': '#1E40AF', 'H': '#3B82F6',
  // Negative / acidic (Braun matte vermilion)
  'D': '#DC2626', 'E': '#DC2626',
  // Conformationally special (warm burnt umber)
  'P': '#9A3412', 'G': '#9A3412',
  // Cysteine (matte brass)
  'C': '#CA8A04',
  // Hydrophilic / polar (slate graphite)
  'S': '#4B5563', 'T': '#4B5563', 'N': '#4B5563', 'Q': '#4B5563',
  // Stop / Gap / Unknown
  '*': '#6B7280', '-': '#E5E7EB', '?': '#D1D5DB'
};

export const NUCLEOTIDE_COLORS = {
  'A': '#1E40AF', // Cobalt
  'C': '#047857', // Sage green
  'G': '#EB5E28', // Braun Orange
  'T': '#DC2626', // Vermilion
  'U': '#DC2626',
  '-': '#D1D5DB',
  'N': '#9CA3AF'
};

// Selection Palette (Disciplined Braun signal colors on neutral ground)
export const SELECTION_COLORS = {
  'diversifying': {
    fill: PALETTE.accent,
    fillLight: withAlpha(PALETTE.accent, 0.12),
    text: PALETTE.accent,
    border: PALETTE.accentBorder,
    label: '▲ Diversifying'
  },
  'purifying': {
    fill: PALETTE.purifying,
    fillLight: withAlpha(PALETTE.purifying, 0.10),
    text: PALETTE.purifying,
    border: PALETTE.purifyingDeep,
    label: '▼ Purifying'
  },
  'neutral': {
    fill: PALETTE.neutral,
    fillLight: withAlpha(PALETTE.neutral, 0.06),
    text: PALETTE.neutralText,
    border: PALETTE.neutralBorder,
    label: '-'
  }
};

// Paradigm 2 Rate Sparkline: ▰▰▰:▰▰▱ (Machined gauge blocks)
export function getRateSparkGlyphs(alpha, beta) {
  function getSegs(v) {
    if (v == null || isNaN(v) || v <= 1e-4) return 0;
    if (v < 0.5) return 1;
    if (v < 3.0) return 2;
    return 3;
  }
  const na = getSegs(alpha);
  const nb = getSegs(beta);
  const sa = '■'.repeat(na) + '□'.repeat(3 - na);
  const sb = '■'.repeat(nb) + '□'.repeat(3 - nb);
  return { sa, sb, na, nb, text: `${sa}:${sb}` };
}
