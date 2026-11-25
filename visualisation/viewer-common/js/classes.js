// classes.js — classed (stepped) colour scales: breaks, class labels and
// colours, and the classification of values. Pure (no DOM): unit-tested with
// `node --test`. A classed variable is drawn by giving the surface each
// cell's class index (0 … k−1) with a colour table of k equal steps
// (`colormaps.js`, names 'cls:<set>:<k>') and the range −0.5 … k−0.5.

/** Colours of the class ramps (k colours, low to high hazard / depth). */
const RAMPS = {
  blues: [[198, 219, 239], [158, 202, 225], [107, 174, 214], [49, 130, 189], [8, 81, 156]],
  ylorrd: [[255, 237, 160], [254, 178, 76], [253, 141, 60], [240, 59, 32], [189, 0, 38]],
  hazard4: [[255, 237, 160], [254, 178, 76], [240, 59, 32], [128, 0, 38]],
  arr: [[26, 152, 80], [166, 217, 106], [255, 255, 153], [254, 174, 97], [215, 48, 39], [128, 0, 38]],
};

/**
 * Class sets. `breaks`: the default class boundaries (k − 1 for k classes;
 * a value v is in class i when breaks[i−1] ≤ v < breaks[i]); `labels`: names
 * for the default number of classes (otherwise the ranges are shown);
 * `editable`: whether the user may change the breaks.
 */
export const CLASS_SETS = {
  depth: {
    title: 'Depth bands', unit: 'm', ramp: 'blues', editable: true,
    breaks: [0.15, 0.3, 0.6, 1.2],
  },
  depth_velocity: {
    title: 'D·V (ARR hazard-curve limits)', unit: 'm²/s', ramp: 'ylorrd', editable: true,
    breaks: [0.3, 0.6, 1.0, 4.0],
  },
  hazard_rating: {
    title: 'Hazard to people (FD2320)', unit: '', ramp: 'hazard4', editable: true,
    breaks: [0.75, 1.25, 2.0],
    labels: ['Low – caution', 'Danger for some', 'Danger for most', 'Danger for all'],
  },
  hazard_class_arr: {
    title: 'ARR combined hazard curves', unit: '', ramp: 'arr', editable: false,
    breaks: [1.5, 2.5, 3.5, 4.5, 5.5],
    labels: ['H1 generally safe', 'H2 unsafe for small vehicles', 'H3 unsafe for vehicles, children, elderly',
             'H4 unsafe for vehicles and people', 'H5 + buildings: structural damage', 'H6 + buildings: failure'],
  },
};

/** k colours ([r, g, b], 0–255) of ramp `name`, interpolated when k differs from its length. */
export function classColors(name, k) {
  const pts = RAMPS[name] || RAMPS.ylorrd;
  if (k === pts.length) return pts.map(c => [...c]);
  const out = [];
  for (let i = 0; i < k; i++) {
    const t = k === 1 ? 0.5 : (i / (k - 1)) * (pts.length - 1);
    const j = Math.min(Math.floor(t), pts.length - 2), u = t - j;
    out.push(pts[j].map((c, m) => Math.round(c + (pts[j + 1][m] - c) * u)));
  }
  return out;
}

/** Name of the stepped colour table of class set `set` with k classes (colormaps.js). */
export const classedName = (set, k) => `cls:${set}:${k}`;

/** {set, k} from a 'cls:<set>:<k>' name, or null. */
export function parseClassedName(name) {
  const m = /^cls:([a-z_]+):(\d+)$/.exec(String(name));
  if (!m || !CLASS_SETS[m[1]]) return null;
  const k = +m[2];
  return k >= 1 && k <= 16 ? { set: m[1], k } : null;
}

/** 256-entry RGBA8 table of k equal steps (the shader samples index i at (i + 0.5)/k). */
export function classTable(set, k, n = 256) {
  const cols = classColors(CLASS_SETS[set].ramp, k), out = new Uint8Array(4 * n);
  for (let i = 0; i < n; i++) {
    const c = cols[Math.min(Math.floor((i / n) * k), k - 1)];
    out[4 * i] = c[0]; out[4 * i + 1] = c[1]; out[4 * i + 2] = c[2]; out[4 * i + 3] = 255;
  }
  return out;
}

/**
 * Class index of each value: the number of breaks ≤ v (0 … breaks.length);
 * NaN stays NaN. `out`: an optional Float32Array to fill.
 */
export function classify(values, breaks, out = new Float32Array(values.length)) {
  const b = breaks, nb = b.length;
  for (let i = 0; i < values.length; i++) {
    const v = values[i];
    if (!(v === v)) { out[i] = NaN; continue; }          // NaN
    let lo = 0, hi = nb;                                   // first break > v
    while (lo < hi) { const mid = (lo + hi) >> 1; if (b[mid] <= v) lo = mid + 1; else hi = mid; }
    out[i] = lo;
  }
  return out;
}

/** Class index of one value (−1 for NaN). */
export function classOf(v, breaks) {
  if (!Number.isFinite(v)) return -1;
  let i = 0;
  while (i < breaks.length && breaks[i] <= v) i++;
  return i;
}

const fmtNum = x => (Math.abs(x) >= 100 ? x.toFixed(0) : String(+x.toFixed(3)));

/**
 * Labels of the classes of `set` with `breaks`: the set's own names when the
 * breaks are its defaults' count, otherwise ranges ("< 0.3", "0.3 – 0.6", "≥ 4").
 */
export function classLabels(set, breaks) {
  const cs = CLASS_SETS[set], u = cs.unit ? ` ${cs.unit}` : '';
  const ranges = [];
  for (let i = 0; i <= breaks.length; i++) {
    if (i === 0) ranges.push(`< ${fmtNum(breaks[0])}${u}`);
    else if (i === breaks.length) ranges.push(`≥ ${fmtNum(breaks[i - 1])}${u}`);
    else ranges.push(`${fmtNum(breaks[i - 1])} – ${fmtNum(breaks[i])}${u}`);
  }
  if (cs.labels && cs.labels.length === breaks.length + 1) {
    return cs.editable ? cs.labels.map((l, i) => `${l} (${ranges[i]})`) : [...cs.labels];
  }
  return ranges;
}

/**
 * Breaks typed by the user ("0.75, 1.25, 2"): finite, strictly increasing,
 * 1 to 15 of them — or null.
 */
export function parseBreaks(text) {
  const parts = String(text ?? '').split(/[\s,;]+/).filter(Boolean);
  if (parts.length < 1 || parts.length > 15) return null;
  const b = parts.map(Number);
  if (!b.every(Number.isFinite)) return null;
  for (let i = 1; i < b.length; i++) if (!(b[i] > b[i - 1])) return null;
  return b;
}

/** CSS gradient with hard steps for a legend bar. */
export function classGradient(set, k) {
  const cols = classColors(CLASS_SETS[set].ramp, k), parts = [];
  cols.forEach((c, i) => {
    const a = (100 * i) / k, b = (100 * (i + 1)) / k;
    parts.push(`rgb(${c.join(',')}) ${a.toFixed(2)}%`, `rgb(${c.join(',')}) ${b.toFixed(2)}%`);
  });
  return `linear-gradient(90deg, ${parts.join(', ')})`;
}
