// palettes.js — colour palettes beyond the built-in ones: the scientific
// colour maps shipped as fixed 256-entry tables (colormaps.json: cmocean,
// Crameri), user-defined palettes (colour stops interpolated in OKLab) and the
// strict validation of palette files (imported by the user, or the server's
// shared palettes.json). Pure (no DOM): unit-tested with `node --test`.
//
// Palette file (docs/DATA_FORMATS.md §6):
//   {"flooda5_palettes": 1,
//    "palettes": [{"name": "My depth",
//                  "stops": [{"pos": 0, "color": "#c6ebff"}, {"pos": 1, "color": "#0a1964"}]}]}

/** Limits of a palette file. */
export const LIMITS = { minStops: 2, maxStops: 32, maxName: 40, maxPalettes: 64, maxBytes: 256 * 1024 };
export const FORMAT_KEY = 'flooda5_palettes';
export const FORMAT_VERSION = 1;

// ---------------------------------------------------------------------------
// Colours: hex, sRGB ↔ OKLab (Ottosson 2020, https://bottosson.github.io/posts/oklab/)
// ---------------------------------------------------------------------------

const HEX = /^#[0-9a-fA-F]{6}$/;

/** '#rrggbb' → [r, g, b] (0–255), or null. */
export function hexToRgb(hex) {
  if (typeof hex !== 'string' || !HEX.test(hex)) return null;
  const v = parseInt(hex.slice(1), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}

/** [r, g, b] (0–255, rounded and clamped) → '#rrggbb'. */
export function rgbToHex(rgb) {
  return '#' + rgb.map(c => Math.min(Math.max(Math.round(c), 0), 255).toString(16).padStart(2, '0')).join('');
}

const toLinear = c => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
const toGamma = c => (c <= 0.0031308 ? 12.92 * c : 1.055 * c ** (1 / 2.4) - 0.055);

/** sRGB [r, g, b] (0–255) → OKLab [L, a, b]. */
export function srgbToOklab([R, G, B]) {
  const r = toLinear(R / 255), g = toLinear(G / 255), b = toLinear(B / 255);
  const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b);
  const m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b);
  const s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b);
  return [
    0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
    1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
    0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s,
  ];
}

/** OKLab [L, a, b] → sRGB [r, g, b] (0–255, unrounded, clamped to the gamut). */
export function oklabToSrgb([L, A, B]) {
  const l = (L + 0.3963377774 * A + 0.2158037573 * B) ** 3;
  const m = (L - 0.1055613458 * A - 0.0638541728 * B) ** 3;
  const s = (L - 0.0894841775 * A - 1.2914855480 * B) ** 3;
  const lin = [
    4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
    -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
    -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s,
  ];
  return lin.map(c => 255 * toGamma(Math.min(Math.max(c, 0), 1)));
}

// ---------------------------------------------------------------------------
// Stops → table
// ---------------------------------------------------------------------------

/**
 * Colour [r, g, b] (0–255, rounded) of the stops at t ∈ [0, 1], interpolated
 * in OKLab. Stops are sorted by `pos`; before the first / after the last the
 * end colours hold; equal positions make a hard step.
 */
export function stopsColorAt(stops, t, labs = stops.map(s => srgbToOklab(hexToRgb(s.color)))) {
  t = Math.min(Math.max(Number.isFinite(t) ? t : 0, 0), 1);
  const n = stops.length;
  if (t <= stops[0].pos) return hexToRgb(stops[0].color);
  if (t >= stops[n - 1].pos) return hexToRgb(stops[n - 1].color);
  let i = 0;
  while (i < n - 2 && t >= stops[i + 1].pos) i++;
  const p0 = stops[i].pos, p1 = stops[i + 1].pos;
  const u = p1 > p0 ? (t - p0) / (p1 - p0) : 1;
  const a = labs[i], b = labs[i + 1];
  return oklabToSrgb([0, 1, 2].map(k => a[k] + (b[k] - a[k]) * u)).map(Math.round);
}

/** n × 1 RGBA8 table (n = 256: the viewer's colour texture) of the stops. */
export function stopsTable(stops, n = 256) {
  const labs = stops.map(s => srgbToOklab(hexToRgb(s.color)));
  const out = new Uint8Array(4 * n);
  for (let i = 0; i < n; i++) {
    const c = stopsColorAt(stops, n === 1 ? 0 : i / (n - 1), labs);
    out[4 * i] = c[0]; out[4 * i + 1] = c[1]; out[4 * i + 2] = c[2]; out[4 * i + 3] = 255;
  }
  return out;
}

/** The stops reversed (positions mirrored: pos → 1 − pos). */
export function reverseStops(stops) {
  return stops.map(s => ({ pos: +(1 - s.pos).toFixed(6), color: s.color })).reverse();
}

/** Sorted copy of the stops (stable for equal positions). */
export const sortStops = stops => stops.map((s, i) => [s, i])
  .sort((a, b) => a[0].pos - b[0].pos || a[1] - b[1]).map(([s]) => ({ pos: s.pos, color: s.color }));

/** k stops sampled evenly from an RGBA8 table (e.g. to start a palette from another). */
export function stopsFromTable(table, k = 5) {
  const n = table.length / 4, out = [];
  for (let j = 0; j < k; j++) {
    const i = Math.round((j / (k - 1)) * (n - 1));
    out.push({ pos: +(j / (k - 1)).toFixed(4), color: rgbToHex([table[4 * i], table[4 * i + 1], table[4 * i + 2]]) });
  }
  return out;
}

// ---------------------------------------------------------------------------
// Fixed tables (colormaps.json)
// ---------------------------------------------------------------------------

/** A table stored as 256 concatenated 'rrggbb' → RGBA8 Uint8Array, or null. */
export function decodeHexTable(hex, n = 256) {
  if (typeof hex !== 'string' || hex.length !== 6 * n || !/^[0-9a-f]+$/i.test(hex)) return null;
  const out = new Uint8Array(4 * n);
  for (let i = 0; i < n; i++) {
    const v = parseInt(hex.slice(6 * i, 6 * i + 6), 16);
    out[4 * i] = (v >> 16) & 255; out[4 * i + 1] = (v >> 8) & 255; out[4 * i + 2] = v & 255; out[4 * i + 3] = 255;
  }
  return out;
}

/**
 * The colour-map families of colormaps.json → [{id, label, group, table}]:
 * ids '<family>:<map>' (e.g. 'cmocean:deep', 'crameri:batlow'); maps whose
 * table does not decode are skipped.
 */
export function parseColormapFile(json) {
  const out = [];
  const fams = json && typeof json === 'object' ? json.families : null;
  if (!fams || typeof fams !== 'object') return out;
  for (const [fid, fam] of Object.entries(fams)) {
    if (!/^[a-z0-9_]+$/.test(fid) || !fam || typeof fam.maps !== 'object') continue;
    for (const [mid, hex] of Object.entries(fam.maps)) {
      const table = /^[a-z0-9_]+$/.test(mid) ? decodeHexTable(hex) : null;
      if (table) out.push({ id: `${fid}:${mid}`, label: mid[0].toUpperCase() + mid.slice(1), group: fid, table });
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Validation of palette files
// ---------------------------------------------------------------------------

const keysOnly = (o, allowed) => Object.keys(o).every(k => allowed.includes(k));
const isObj = o => o !== null && typeof o === 'object' && !Array.isArray(o);

/** A palette name: a string of 1–40 characters (trimmed), no control characters; or null. */
export function cleanName(name) {
  if (typeof name !== 'string') return null;
  const s = name.trim();
  // eslint-disable-next-line no-control-regex
  if (s.length < 1 || s.length > LIMITS.maxName || /[\u0000-\u001f\u007f-\u009f]/.test(s)) return null;
  return s;
}

/**
 * Check one palette {name, stops: [{pos, color}]} strictly: no other keys;
 * 2–32 stops; positions numbers in [0, 1], sorted (non-decreasing); colours
 * '#rrggbb'. → {ok: true, palette} (name trimmed, colours lower case) or
 * {ok: false, error}.
 */
export function validatePalette(p) {
  if (!isObj(p)) return { ok: false, error: 'a palette must be an object' };
  if (!keysOnly(p, ['name', 'stops'])) return { ok: false, error: 'a palette has only "name" and "stops"' };
  const name = cleanName(p.name);
  if (name === null) return { ok: false, error: `a name must be text of 1–${LIMITS.maxName} characters` };
  const st = p.stops;
  if (!Array.isArray(st) || st.length < LIMITS.minStops || st.length > LIMITS.maxStops) {
    return { ok: false, error: `"${name}": ${LIMITS.minStops}–${LIMITS.maxStops} stops are needed` };
  }
  const stops = [];
  for (const s of st) {
    if (!isObj(s) || !keysOnly(s, ['pos', 'color'])) return { ok: false, error: `"${name}": a stop is {"pos", "color"}` };
    if (typeof s.pos !== 'number' || !Number.isFinite(s.pos) || s.pos < 0 || s.pos > 1) {
      return { ok: false, error: `"${name}": positions are numbers from 0 to 1` };
    }
    if (typeof s.color !== 'string' || !HEX.test(s.color)) return { ok: false, error: `"${name}": colours are "#rrggbb"` };
    if (stops.length && s.pos < stops[stops.length - 1].pos) return { ok: false, error: `"${name}": positions must be sorted` };
    stops.push({ pos: s.pos, color: s.color.toLowerCase() });
  }
  return { ok: true, palette: { name, stops } };
}

/**
 * Check a palette file (an object, or its JSON text): {"flooda5_palettes": 1,
 * "palettes": [1–64 palettes]}, no other keys, names unique (ignoring case),
 * text at most 256 KiB. → {ok: true, palettes} or {ok: false, error}.
 */
export function validatePaletteFile(input) {
  let doc = input;
  if (typeof input === 'string') {
    if (input.length > LIMITS.maxBytes) return { ok: false, error: 'the file is larger than 256 KiB' };
    try { doc = JSON.parse(input); } catch { return { ok: false, error: 'not JSON' }; }
  }
  if (!isObj(doc) || !keysOnly(doc, [FORMAT_KEY, 'palettes']) || doc[FORMAT_KEY] !== FORMAT_VERSION) {
    return { ok: false, error: `not a palette file ({"${FORMAT_KEY}": ${FORMAT_VERSION}, "palettes": […]})` };
  }
  const ps = doc.palettes;
  if (!Array.isArray(ps) || ps.length < 1 || ps.length > LIMITS.maxPalettes) {
    return { ok: false, error: `1–${LIMITS.maxPalettes} palettes are needed` };
  }
  const out = [], seen = new Set();
  for (const p of ps) {
    const r = validatePalette(p);
    if (!r.ok) return r;
    const key = r.palette.name.toLowerCase();
    if (seen.has(key)) return { ok: false, error: `"${r.palette.name}" appears twice` };
    seen.add(key);
    out.push(r.palette);
  }
  return { ok: true, palettes: out };
}

/** A palette file (pretty JSON text) of palettes [{name, stops}]. */
export function paletteFileText(palettes) {
  const doc = { [FORMAT_KEY]: FORMAT_VERSION,
                palettes: palettes.map(p => ({ name: p.name, stops: p.stops.map(s => ({ pos: s.pos, color: s.color })) })) };
  return JSON.stringify(doc, null, 2) + '\n';
}

/** `name`, or "name (2)", "name (3)", … — the first not in `taken` (ignoring case), ≤ 40 characters. */
export function uniqueName(name, taken) {
  const has = n => [...taken].some(t => t.toLowerCase() === n.toLowerCase());
  if (!has(name)) return name;
  for (let k = 2; ; k++) {
    const suf = ` (${k})`, n = name.slice(0, LIMITS.maxName - suf.length) + suf;
    if (!has(n)) return n;
  }
}
