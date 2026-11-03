/**
 * Pure numeric helpers: tail probabilities and multiple-testing corrections.
 */

// Complementary error function. Below x = 10 this is Numerical Recipes erfcc,
// whose relative error stays under 1.2e-7 where 1 - erf(x) would round to 0.
// Beyond that the asymptotic series is more accurate: its first omitted term
// is under 2e-10 of the value at x = 10 and shrinks from there.
export function erfc(x) {
  const z = Math.abs(x);
  let r;
  if (z >= 10) {
    const inv2z2 = 1 / (2 * z * z);
    let term = 1, sum = 1;
    for (let k = 1; k <= 5; k++) {
      term *= -(2 * k - 1) * inv2z2;
      sum += term;
    }
    r = Math.exp(-z * z) / (z * Math.sqrt(Math.PI)) * sum;
  } else {
    const t = 1 / (1 + 0.5 * z);
    r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
      t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
      t * (-0.82215223 + t * 0.17087277)))))))));
  }
  return x >= 0 ? r : 2 - r;
}

// ln Γ(x) for x > 0 (Lanczos, g = 7, nine coefficients)
const LANCZOS_COEFFICIENTS = [
  0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
  -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
  1.5056327351493116e-7
];

export function logGamma(x) {
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1 - x);
  const y = x - 1;
  let a = LANCZOS_COEFFICIENTS[0];
  const t = y + 7.5;
  for (let i = 1; i < 9; i++) a += LANCZOS_COEFFICIENTS[i] / (y + i);
  return 0.5 * Math.log(2 * Math.PI) + (y + 0.5) * Math.log(t) - t + Math.log(a);
}

// Regularized upper incomplete gamma Q(a, x) = Γ(a, x) / Γ(a): the series for P
// below x = a + 1, the Lentz continued fraction for Q above (Numerical Recipes)
export function gammaQ(a, x) {
  if (x <= 0) return 1;
  const lead = -x + a * Math.log(x) - logGamma(a);
  if (x < a + 1) {
    let ap = a, del = 1 / a, sum = del;
    for (let n = 0; n < 500 && Math.abs(del) > Math.abs(sum) * 1e-15; n++) {
      ap += 1;
      del *= x / ap;
      sum += del;
    }
    return Math.max(0, 1 - sum * Math.exp(lead));
  }
  const tiny = 1e-300;
  let b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
  for (let i = 1; i < 500; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < tiny) d = tiny;
    c = b + an / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-15) break;
  }
  return Math.exp(lead) * h;
}

// P(χ²_df > x). df = 1 goes through erfc, which keeps full relative precision
// far into the tail; other df use the regularized gamma.
export function chi2Survival(x, df = 1) {
  if (!(x > 0)) return 1;
  if (df === 1) return erfc(Math.sqrt(x / 2));
  return gammaQ(df / 2, x / 2);
}

// Indices of p sorted ascending; ties keep their input order
function ascendingOrder(p) {
  return p.map((_, i) => i).sort((a, b) => p[a] - p[b] || a - b);
}

// Benjamini-Hochberg q-values in input order. Step-up: q at rank i is the
// minimum of p·m/j over ranks j ≥ i, so q never decreases with p and tied p
// share one q. m defaults to the number of p-values but can be larger, e.g.
// all codon sites when only some were tested.
export function benjaminiHochberg(p, m = p.length) {
  const order = ascendingOrder(p);
  const q = new Array(p.length);
  let running = 1;
  for (let rank = order.length; rank >= 1; rank--) {
    const i = order[rank - 1];
    running = Math.min(running, p[i] * m / rank);
    q[i] = running;
  }
  return q;
}

// Holm-Bonferroni adjusted p-values in input order: (m - rank + 1)·p, then a
// running maximum down the ranks so adjusted p never decreases with raw p
export function holm(p, m = p.length) {
  const order = ascendingOrder(p);
  const scaled = new Array(p.length);
  order.forEach((i, k) => { scaled[i] = (m - k) * p[i]; });
  return runningMaxInOrder(scaled, order);
}

// Running maximum of already-adjusted p-values taken in ascending raw-p order,
// capped at 1. This makes stored corrections that lost step-down monotonicity
// (HyPhy 2.5 aBSREL) behave like a true step-down at every threshold; on
// monotone input it changes nothing.
export function enforceStepDown(adjusted, raw) {
  return runningMaxInOrder(adjusted, ascendingOrder(raw));
}

function runningMaxInOrder(values, order) {
  const out = new Array(values.length);
  let running = 0;
  for (const i of order) {
    running = Math.max(running, Math.min(1, values[i]));
    out[i] = running;
  }
  return out;
}

// Exact 95% (by default) interval for the number of false calls among calls with
// the given error probabilities: the Poisson-binomial distribution by convolution,
// bounded by the first counts whose CDF reaches alpha/2 and 1 - alpha/2. Matches
// hyphy3 compute_enfp_ci and HyPhy 2.5 fubar.ComputeENFP_CI; no calls gives [0, 0].
export function poissonBinomialCI(errProbs, alpha = 0.05) {
  let pdf = [1];
  for (const p of errProbs) {
    const next = new Array(pdf.length + 1).fill(0);
    pdf.forEach((mass, k) => {
      next[k] += mass * (1 - p);
      next[k + 1] += mass * p;
    });
    pdf = next;
  }
  const firstReaching = (level) => {
    let cdf = 0;
    for (let k = 0; k < pdf.length; k++) {
      cdf += pdf[k];
      if (cdf >= level) return k;
    }
    return pdf.length - 1;
  };
  return [firstReaching(alpha / 2), firstReaching(1 - alpha / 2)];
}
