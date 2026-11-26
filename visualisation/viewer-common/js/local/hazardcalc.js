// hazardcalc.js — the hazard variables and statistics over a run, for runs
// read in the browser (local files). A line-by-line port of
// analysis/hazard.jl (the server's and FloodA5Analysis's implementation):
// keep the two in step — test/test_vis_local.jl checks this port against the
// server on a generated run (viewer-common/test/fixtures/local_run.json).
// Pure: unit-tested with `node --test`.
//
// Arithmetic follows the Julia code: inputs are Float32 (as the server reads
// frames), computed in double, results rounded to Float32. Results agree with
// the server's to the last bit except where Math.pow / Math.sin differ from
// Julia's by an ulp (the Chen et al. hazard degree).

/** Debris factor by land use, rows d ≤ 0.25; 0.25 < d ≤ 0.75 (v ≤ 2); d > 0.75 or v > 2 (FD2321 Table 3.1). */
export const DEBRIS_TABLES = {
  fd2320: [0.5, 1.0, 1.0],
  pasture: [0.0, 0.0, 0.5],
  woodland: [0.0, 0.5, 1.0],
  urban: [0.0, 1.0, 1.0],
  none: [0.0, 0.0, 0.0],
};
export const LAND_USES = ['fd2320', 'pasture', 'woodland', 'urban', 'none'];

/** HazardParams defaults (analysis/hazard.jl). */
export const HAZARD_DEFAULTS = Object.freeze({
  landuse: 'fd2320', mass: 75.1, height: 1.77, theta: 75.0, alpha: 3.064, beta: 0.303, flood_depth: 0.1,
});

// Julia's tryparse(Float64, s): nothing for an empty or malformed string
function tryparseFloat(s) {
  if (s === undefined || s === null) return null;
  const t = String(s);
  if (!/^\s*[+-]?(\d+\.?\d*([eE][+-]?\d+)?|\.\d+([eE][+-]?\d+)?|Inf|NaN|Infinity)\s*$/.test(t)) return null;
  const u = t.trim().replace(/^([+-]?)Inf(inity)?$/, '$1Infinity');
  return Number(u);
}
const num = (q, k, dflt, lo, hi) => {
  const x = tryparseFloat(q[k]);
  return x === null || !Number.isFinite(x) ? dflt : Math.min(Math.max(x, lo), hi);
};

/**
 * hazard_params(q): parameters from a query (an object of strings, e.g.
 * Object.fromEntries(URLSearchParams)); missing or invalid values take the
 * defaults, numbers are clamped (as the server).
 */
export function hazardParams(q = {}) {
  const d = HAZARD_DEFAULTS;
  const lu = q.landuse === undefined ? d.landuse : String(q.landuse);
  return {
    landuse: LAND_USES.includes(lu) ? lu : d.landuse,
    mass: num(q, 'mass', d.mass, 1.0, 300.0),
    height: num(q, 'height', d.height, 0.3, 2.5),
    theta: num(q, 'theta', d.theta, 1.0, 90.0),
    alpha: num(q, 'alpha', d.alpha, 0.01, 100.0),
    beta: num(q, 'beta', d.beta, 0.0, 5.0),
    flood_depth: num(q, 'flood_depth', d.flood_depth, 0.0, 100.0),
  };
}

/** A string key of parameters `p` (caches). */
export const hazardKey = p => [p.landuse, p.mass, p.height, p.theta, p.alpha, p.beta, p.flood_depth].join('|');

export const depthVelocity = (d, v) => (d > 0 ? d * Math.max(v, 0) : 0);

export function debrisFactor(d, v, t) {
  if (!(d > 0)) return 0.0;
  if (d <= 0.25) return t[0];
  return (d > 0.75 || v > 2) ? t[2] : t[1];
}

export const hazardRating = (d, v, t) => (d > 0 ? d * (Math.max(v, 0) + 0.5) + debrisFactor(d, v, t) : 0.0);

export function arrHazardClass(d, v) {
  if (!(d > 0)) return 0;
  v = Math.max(v, 0);
  const dv = d * v;
  if (dv <= 0.3 && d <= 0.3 && v <= 2) return 1;
  if (dv <= 0.6 && d <= 0.5 && v <= 2) return 2;
  if (dv <= 0.6 && d <= 1.2 && v <= 2) return 3;
  if (dv <= 1.0 && d <= 2.0 && v <= 2) return 4;
  if (dv <= 4.0 && d <= 4.0 && v <= 4) return 5;
  return 6;
}

// sind/cosd: exact at multiples of 90° (as Julia's sincosd)
function sincosd(deg) {
  if (deg === 90) return [1, 0];
  if (deg === 0) return [0, 1];
  const r = deg * Math.PI / 180;
  return [Math.sin(r), Math.cos(r)];
}

/** Chen et al. (2019) eq. (10): incipient velocity for toppling, U_c (m/s), in water of depth h. */
export function chenCriticalVelocity(h, p) {
  const { mass, height, theta, alpha, beta } = p;
  const rho = 1000.0, a1 = 0.735, b1 = 0.265, a2 = 1.015e-3, b2 = -4.927e-3, gamma = 10.0;
  if (!(h > 0)) return Infinity;
  if (h >= height) return 0.0;
  const r = h / height;
  const [s, c] = sincosd(theta);
  const rad = mass / (rho * (h * h)) * (gamma * c + s) -
              (a1 * (r * r) + b1 * r) * (a2 * mass + b2) * (gamma * c / (h * height) + s / (h * h));
  if (!(rad > 0)) return 0.0;
  return alpha * r ** beta * Math.sqrt(rad);
}

/** Chen et al. (2019) eq. (14): HD = min(1, v / U_c). */
export function hazardDegree(d, v, p) {
  if (!(d > 0)) return 0.0;
  if (d >= p.height) return 1.0;
  const uc = chenCriticalVelocity(d, p);
  if (!(uc > 0)) return 1.0;
  return Math.min(1.0, Math.max(v, 0) / uc);
}

/** The hazard variables, by the names the server serves: name → f(d, v, p, debrisTable). */
export const HAZARD_VARS = {
  depth_velocity: (d, v) => depthVelocity(d, v),
  hazard_rating: (d, v, p, t) => hazardRating(d, v, t),
  hazard_degree: (d, v, p) => hazardDegree(d, v, p),
  hazard_class_arr: (d, v) => arrHazardClass(d, v),
};
export const HAZARD_NAMES = Object.keys(HAZARD_VARS).sort();
export const isDerived = name => Object.prototype.hasOwnProperty.call(HAZARD_VARS, name);

/** The inputs of the hazard variables. */
export const DERIVED_INPUTS = ['water_depth', 'velocity'];

/**
 * derive_hazard: hazard variable `name` elementwise over `depth` and
 * `velocity` (equal lengths) → Float32Array; NaN where an input is not finite.
 */
export function deriveHazard(name, depth, velocity, p = HAZARD_DEFAULTS) {
  const f = HAZARD_VARS[name];
  if (!f) throw new Error(`not a hazard variable: ${name}`);
  if (depth.length !== velocity.length) throw new Error('depth and velocity lengths differ');
  const t = DEBRIS_TABLES[p.landuse];
  const out = new Float32Array(depth.length);
  for (let i = 0; i < depth.length; i++) {
    const d = depth[i], v = velocity[i];
    out[i] = Number.isFinite(d) && Number.isFinite(v) ? f(d, v, p, t) : NaN;
  }
  return out;
}

/** The statistics runStatistics returns (STAT_VARS). */
export const STAT_VARS = ['max_water_depth', 'max_velocity', 'max_depth_velocity', 'max_hazard_rating',
                          'max_hazard_degree', 'max_hazard_class_arr', 'time_max_water_depth',
                          'time_first_flooded', 'flood_duration'];

/**
 * run_statistics: per-cell statistics over frames at `times` (s), in one
 * pass. `frame(k)` (k 0-based) returns [depth, velocity] (length n;
 * velocity may be null: then the velocity and hazard statistics are left
 * out). → { name: Float32Array }.
 */
export function runStatistics(times, frame, n, p = HAZARD_DEFAULTS) {
  const nf = times.length;
  const mx = { max_water_depth: new Float32Array(n).fill(NaN) };
  const tfirst = new Float32Array(n).fill(NaN), dur = new Float32Array(n), tmax = new Float32Array(n).fill(NaN);
  let haveV = false;
  const t = DEBRIS_TABLES[p.landuse];
  const fd = p.flood_depth;
  const maxin = (a, i, x) => {
    const y = Math.fround(x);
    if (Number.isNaN(a[i]) || y > a[i]) a[i] = y;
  };
  for (let k = 0; k < nf; k++) {
    const [d, v] = frame(k);
    if (d.length !== n) throw new Error(`frame ${k + 1} has ${d.length} values, expected ${n}`);
    const w = nf === 1 ? 0.0 : (times[Math.min(k + 1, nf - 1)] - times[Math.max(k - 1, 0)]) / 2;
    const md = mx.max_water_depth;
    for (let i = 0; i < n; i++) {
      const x = d[i];
      if (!Number.isFinite(x)) continue;
      if (Number.isNaN(md[i]) || x > md[i]) {
        md[i] = x;
        if (x > 0) tmax[i] = times[k];
      }
      if (x >= fd) {
        if (Number.isNaN(tfirst[i])) tfirst[i] = times[k];
        dur[i] += w;
      }
    }
    if (!v) continue;
    if (!haveV) {
      haveV = true;
      for (const s of ['max_velocity', 'max_depth_velocity', 'max_hazard_rating', 'max_hazard_degree', 'max_hazard_class_arr']) {
        mx[s] = new Float32Array(n).fill(NaN);
      }
    }
    const { max_velocity: mv, max_depth_velocity: mdv, max_hazard_rating: mhr,
            max_hazard_degree: mhd, max_hazard_class_arr: mar } = mx;
    for (let i = 0; i < n; i++) {
      const x = d[i], u = v[i];
      if (!(Number.isFinite(x) && Number.isFinite(u))) continue;
      maxin(mv, i, u);
      maxin(mdv, i, depthVelocity(x, u));
      maxin(mhr, i, hazardRating(x, u, t));
      maxin(mhd, i, hazardDegree(x, u, p));
      maxin(mar, i, arrHazardClass(x, u));
    }
  }
  mx.time_max_water_depth = tmax;
  mx.time_first_flooded = tfirst;
  mx.flood_duration = dur;
  return mx;
}
