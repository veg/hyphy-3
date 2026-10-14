// Design tokens for code that draws on canvas or writes SVG/HTML strings.
// Values mirror the :root custom properties in styles/viewer.css; keep the two in sync.

export const PALETTE = {
  // Signal hues: accent marks the finding (diversifying), purifying is its mirror
  accent: '#EB5E28',          // --accent: Braun Signal Orange
  accentDeep: '#9A3412',      // Burnt outline for accent glyphs
  accentBorder: '#D9480F',
  accentPale: '#FED7AA',
  purifying: '#1E40AF',       // --purifying: Technical Cobalt Blue
  purifyingDeep: '#1E3A8A',   // Outline for purifying glyphs
  synonymous: '#047857',      // --synonymous: Deep green
  amber: '#D97706',           // Moderate EBF
  amberOchre: '#B45309',
  amberVivid: '#F59E0B',

  // Surfaces
  chassis: '#F2F2EF',         // --chassis: Warm matte off-white
  panel: '#FFFFFF',           // --panel: Crisp flat white
  panelHeader: '#FAF8F5',     // --panel-header
  inset: '#EDEDE9',           // --inset: Recessed wells and table heads
  rowHover: '#F8F8F5',        // --row-hover
  rowSelected: '#FEF3C7',     // --row-selected: Soft warm amber
  cellBg: '#FAF9F6',          // Canvas cell and label backdrop
  trackBg: '#F8F8F6',
  overviewBg: '#F4F4F1',

  // Ink and lines
  ink: '#1A1A18',             // --ink: Rich charcoal
  knockout: '#FFFFFF',        // --knockout: Text on filled ink or signal backgrounds
  ancestorInk: '#4A4A45',
  stone: '#767672',           // --stone: Secondary warm label
  branchGrey: '#8E8E8A',      // Untested branches
  dotFaint: '#A8A8A2',
  dotMuted: '#B8B8B2',        // Muted consensus dots
  hairline: '#D1D1CB',        // --hairline: Precision 1px structural line
  overviewEmpty: '#D6D6D0',
  gridLine: '#EAEAE6',        // --grid-line: Subtle cell grid

  // Cool greys
  slate: '#475569',           // --slate
  slateLight: '#94A3B8',      // --slate-light
  treeLine: '#CBD5E1',
  slatePale: '#E2E8F0',       // --slate-pale
  slateFaint: '#F8FAFC',      // --slate-faint
  neutral: '#9CA3AF',
  neutralText: '#6B7280',
  neutralBorder: '#D1D5DB'
};

// '#RRGGBB' at the given opacity, as an rgba() string
export function withAlpha(hex, alpha) {
  const n = parseInt(hex.slice(1), 16);
  return `rgba(${n >> 16}, ${(n >> 8) & 255}, ${n & 255}, ${alpha})`;
}
