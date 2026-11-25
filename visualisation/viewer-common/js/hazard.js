// hazard.js — the hazard parameters of the viewer (kept in the browser, sent
// to the server in the query string of hazard variables, ranges and
// statistics). The hazard variables themselves are computed by the server
// (analysis/hazard.jl). Pure: unit-tested with `node --test`.

/** Defaults: analysis/hazard.jl HazardParams. */
export const DEFAULT_HAZARD = Object.freeze({
  landuse: 'fd2320', person: 'adult', mass: 75.1, height: 1.77,
  theta: 75, alpha: 3.064, beta: 0.303, flood_depth: 0.1,
});

/** Land uses for the Defra debris factor (one for the whole domain). */
export const LAND_USES = [
  ['fd2320', 'any (FD2320)'],
  ['pasture', 'pasture / arable'],
  ['woodland', 'woodland'],
  ['urban', 'urban'],
  ['none', 'no debris (DF 0)'],
];

/** Chen et al. (2019) §4: the adult and child of their tests. */
export const PERSON_PRESETS = {
  adult: { mass: 75.1, height: 1.77, label: 'adult (75.1 kg, 1.77 m)' },
  child: { mass: 25.5, height: 1.26, label: 'child (25.5 kg, 1.26 m)' },
};

// plausible ranges (the server clamps to the same: hazard_params)
const LIMITS = {
  mass: [1, 300], height: [0.3, 2.5], theta: [1, 90], alpha: [0.01, 100], beta: [0, 5], flood_depth: [0, 100],
};

/**
 * Complete, valid parameters from a partial object (e.g. parsed from
 * storage): unknown land uses and invalid numbers take the defaults; numbers
 * are clamped; a preset person's mass and height are the preset's.
 */
export function normalizeHazard(p = {}) {
  const out = { ...DEFAULT_HAZARD };
  if (p && LAND_USES.some(([k]) => k === p.landuse)) out.landuse = p.landuse;
  for (const [k, [lo, hi]] of Object.entries(LIMITS)) {
    const x = Number(p?.[k]);
    if (p && p[k] !== undefined && p[k] !== null && p[k] !== '' && Number.isFinite(x)) out[k] = Math.min(Math.max(x, lo), hi);
  }
  out.person = p && (p.person === 'custom' || PERSON_PRESETS[p.person]) ? p.person : DEFAULT_HAZARD.person;
  if (PERSON_PRESETS[out.person]) {
    out.mass = PERSON_PRESETS[out.person].mass;
    out.height = PERSON_PRESETS[out.person].height;
  }
  return out;
}

/** The query string of parameters `p` (normalised), in a fixed order: "landuse=…&mass=…". */
export function hazardQuery(p) {
  const n = normalizeHazard(p);
  return ['landuse', 'mass', 'height', 'theta', 'alpha', 'beta', 'flood_depth']
    .map(k => `${k}=${encodeURIComponent(String(n[k]))}`).join('&');
}

/** One-line description of the parameters (legend, log). */
export function describeHazard(p) {
  const n = normalizeHazard(p);
  const lu = (LAND_USES.find(([k]) => k === n.landuse) || [])[1] || n.landuse;
  return `land use: ${lu}; person: ${n.mass} kg, ${n.height} m (θ ${n.theta}°, α ${n.alpha}, β ${n.beta}); ` +
         `flooded if depth ≥ ${n.flood_depth} m`;
}
