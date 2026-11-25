// colormaps.js — colour maps for the globe viewer, as 256-entry RGBA8 lookup
// tables (uploaded as a texture) and CSS gradients for the legend. Pure.
// Names 'cls:<set>:<k>' are classed (stepped) tables: classes.js.

import { parseClassedName, classTable, classGradient } from './classes.js';

// Control points [t, r, g, b] (0–1 / 0–255), interpolated linearly.
const MAPS = {
  // water depth: pale cyan → deep blue
  depth: [[0, 198, 235, 255], [0.25, 110, 190, 240], [0.5, 40, 130, 220], [0.75, 20, 70, 170], [1, 10, 25, 100]],
  // brown flood water: a narrow range around the colour of open brown floodwater in
  // aerial photographs (≈ rgb(112, 100, 92)), lighter where shallow
  brown: [[0, 152, 138, 120], [0.5, 114, 101, 90], [1, 80, 68, 60]],
  // perceptually uniform, readable with colour-vision deficiencies
  cividis: [[0, 0, 34, 78], [0.25, 61, 78, 110], [0.5, 124, 123, 120], [0.75, 188, 175, 111], [1, 255, 234, 70]],
  // perceptually uniform, for water level / general use
  viridis: [[0, 68, 1, 84], [0.25, 59, 82, 139], [0.5, 33, 145, 140], [0.75, 94, 201, 98], [1, 253, 231, 37]],
  // speed: dark → hot
  magma: [[0, 0, 0, 4], [0.25, 81, 18, 124], [0.5, 183, 55, 121], [0.75, 252, 137, 97], [1, 252, 253, 191]],
  // turbo-like rainbow (kept for continuity with the earlier viewer)
  turbo: [[0, 48, 18, 59], [0.15, 70, 118, 238], [0.3, 27, 208, 213], [0.45, 97, 252, 108],
          [0.6, 210, 233, 53], [0.75, 254, 155, 45], [0.9, 218, 56, 9], [1, 122, 4, 3]],
  // bed elevation
  terrain: [[0, 40, 80, 60], [0.2, 80, 130, 70], [0.45, 160, 170, 100], [0.7, 175, 140, 100],
            [0.9, 210, 200, 190], [1, 250, 250, 250]],
  // fraction 0–1
  blues: [[0, 247, 251, 255], [0.5, 107, 174, 214], [1, 8, 48, 107]],
  greys: [[0, 245, 245, 245], [1, 40, 40, 40]],
};

export const COLORMAP_NAMES = Object.keys(MAPS);

/** Names shown in the palette menu. */
export const COLORMAP_LABELS = {
  depth: 'Water (blue)', brown: 'Brown flood water', viridis: 'Viridis', cividis: 'Cividis',
  magma: 'Magma', turbo: 'Turbo', blues: 'Blues', greys: 'Greys', terrain: 'Terrain',
};

/** Colour (r, g, b in 0–255) of map `name` at t ∈ [0, 1]. */
export function colorAt(name, t) {
  const pts = MAPS[name] || MAPS.viridis;
  t = Math.min(Math.max(Number.isFinite(t) ? t : 0, 0), 1);
  let i = 0;
  while (i < pts.length - 2 && t > pts[i + 1][0]) i++;
  const [t0, r0, g0, b0] = pts[i], [t1, r1, g1, b1] = pts[i + 1];
  const u = t1 > t0 ? (t - t0) / (t1 - t0) : 0;
  return [r0 + (r1 - r0) * u, g0 + (g1 - g0) * u, b0 + (b1 - b0) * u].map(Math.round);
}

/** 256 × 1 RGBA8 lookup table for map `name`. */
export function colormapTable(name, n = 256) {
  const cls = parseClassedName(name);
  if (cls) return classTable(cls.set, cls.k, n);
  const out = new Uint8Array(4 * n);
  for (let i = 0; i < n; i++) {
    const [r, g, b] = colorAt(name, i / (n - 1));
    out[4 * i] = r; out[4 * i + 1] = g; out[4 * i + 2] = b; out[4 * i + 3] = 255;
  }
  return out;
}

/** CSS linear-gradient for a legend bar. */
export function cssGradient(name, stops = 9) {
  const cls = parseClassedName(name);
  if (cls) return classGradient(cls.set, cls.k);
  const parts = [];
  for (let i = 0; i < stops; i++) {
    const t = i / (stops - 1);
    const [r, g, b] = colorAt(name, t);
    parts.push(`rgb(${r},${g},${b}) ${(100 * t).toFixed(0)}%`);
  }
  return `linear-gradient(90deg, ${parts.join(', ')})`;
}
