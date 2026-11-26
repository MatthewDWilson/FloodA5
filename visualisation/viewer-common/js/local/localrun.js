// localrun.js — the viewer's /api endpoints for a FloodA5 output file read in
// the browser (a local file: nothing is uploaded). The client-side twin of
// the server's replay path (visualisation/web_replay.jl +
// VisualisationServer.jl): the same responses, computed from the file
// through a `reader` (h5reader.js wraps h5wasm; the tests use a fixture).
// Keep the two in step — test/test_vis_local.jl compares them on a
// generated run. Pure (no DOM, no worker globals).
//
// Reader interface (synchronous; paths absolute, e.g. "/results/t"):
//   has(path) → bool               isDataset(path) → bool
//   keys(path) → [names]           attrs(path) → {name: value}
//   shape(path) → [dims]           read(path) → typed array | array | string | number
//   row(path, k) → frame k (0-based) of an (n_frames, n_cells) dataset
//   col(path, c) → column c of it (all frames)

import {
  hazardParams, hazardKey, isDerived, HAZARD_NAMES, DERIVED_INPUTS, deriveHazard, STAT_VARS, runStatistics,
} from './hazardcalc.js';

export const PROTOCOL = 2;

/** The frame variables the viewer serves (web_replay.jl _WEB_REPLAY_VARS), in that order. */
export const WEB_VARS = ['water_depth', 'velocity', 'vel_u', 'vel_v', 'saturation'];

// /mesh datasets that are not per-cell static variables (FA.mesh_geometry)
const RESERVED = ['cell_ids', 'center_lons', 'center_lats', 'elevations', 'cell_area', 'polygon_offsets',
                  'polygon_lon', 'polygon_lat', 'neighbours', 'edge_cells', 'source'];

const FRAME_CACHE = 48;          // frames × variables kept (as the server's cache_size)

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** A Float32 value as the shortest decimal that reads back as the same Float32 (as JSON3 writes it). */
export function f32short(x) {
  if (!Number.isFinite(x)) return x;
  const f = Math.fround(x);
  for (let p = 1; p <= 9; p++) {
    const y = Number(f.toPrecision(p));
    if (Math.fround(y) === f) return y;
  }
  return f;
}

const fl = x => (Number.isFinite(x) ? x : null);              // JSON has no NaN
const fl32 = x => (Number.isFinite(x) ? f32short(x) : null);

/**
 * JSON as JSON3 writes it with allow_inf (bare NaN, Infinity, -Infinity),
 * non-finite numbers given as the strings run_record uses ("NaN", "Inf", "-Inf").
 */
export function parseJsonInf(text) {
  const s = String(text).replace(/"(?:[^"\\]|\\.)*"|-?Infinity|NaN/g,
    m => (m[0] === '"' ? m : m === 'NaN' ? '"NaN"' : m[0] === '-' ? '"-Inf"' : '"Inf"'));
  return JSON.parse(s);
}

/** Julia's round(Int, x) on 1 + num/den (ties to even), exactly. */
function roundRational(num, den) {
  let q = Math.floor(num / den);
  const r = num - q * den;
  if (2 * r > den || (2 * r === den && (1 + q) % 2 === 1)) q += 1;   // 1 + q + r/den, ties to the even integer
  return 1 + q;
}

/** Frame ids (1-based) of up to `m` frames spread evenly over n (Julia: unique(round.(Int, range(1, n; length = m)))). */
export function pickFrames(n, m) {
  if (n <= m) return Array.from({ length: n }, (_, i) => i + 1);
  const out = [];
  for (let j = 0; j < m; j++) {
    const k = m === 1 ? 1 : roundRational(j * (n - 1), m - 1);
    if (out[out.length - 1] !== k && !out.includes(k)) out.push(k);
  }
  return out;
}

/** q-quantile of sorted values, as the server's _sorted_quantile (NaN if empty). */
export function sortedQuantile(a, q) {
  if (!a.length) return NaN;
  const h = (a.length - 1) * q + 1;
  const lo = Math.floor(h), hi = Math.min(lo + 1, a.length);
  return a[lo - 1] + (h - lo) * (a[hi - 1] - a[lo - 1]);
}

/** JSON-safe plain data as FO._plain: non-finite numbers → "Inf" / "-Inf" / "NaN". */
function plain(x) {
  if (typeof x === 'number') return Number.isFinite(x) ? x : Number.isNaN(x) ? 'NaN' : x > 0 ? 'Inf' : '-Inf';
  if (typeof x === 'bigint') return Number(x);
  if (Array.isArray(x) || ArrayBuffer.isView(x)) return Array.from(x, plain);
  if (x && typeof x === 'object') {
    const o = {};
    for (const [k, v] of Object.entries(x)) o[k] = plain(v);
    return o;
  }
  return x;
}

const toF32 = a => (a instanceof Float32Array ? a : Float32Array.from(a));
const toNum = a => (ArrayBuffer.isView(a) && !(a instanceof BigInt64Array) && !(a instanceof BigUint64Array)
  ? a : Array.from(a, Number));

// ---------------------------------------------------------------------------
// Mesh (VisualisationServer.encode_mesh, depth_datum)
// ---------------------------------------------------------------------------

/**
 * The binary mesh (protocol 2, `/api/mesh`) from CSR polygons: cell i's
 * vertices are lon/lat[offsets[i] … offsets[i+1]); a closing repeat of the
 * first vertex is dropped. Byte-identical to the server's encode_mesh.
 */
export function encodeMesh(offsets, plon, plat, elevation, base) {
  const n = offsets.length - 1;
  if (n <= 0) throw new Error('encode_mesh: empty mesh');
  if (elevation.length !== n || base.length !== n) throw new Error('encode_mesh: lengths differ');
  const start = new Array(n), len = new Array(n);
  let nv = 0, loMin = Infinity, loMax = -Infinity, laMin = Infinity, laMax = -Infinity;
  for (let i = 0; i < n; i++) {
    const a = Number(offsets[i]), b = Number(offsets[i + 1]);
    let m = b - a;
    if (m > 1 && plon[a] === plon[b - 1] && plat[a] === plat[b - 1]) m -= 1;
    if (m < 3) throw new Error('encode_mesh: a cell has fewer than 3 vertices');
    start[i] = a; len[i] = m; nv += m;
    for (let j = a; j < a + m; j++) {
      const x = plon[j], y = plat[j];
      if (x < loMin) loMin = x; if (x > loMax) loMax = x;
      if (y < laMin) laMin = y; if (y > laMax) laMax = y;
    }
  }
  const olon = (loMin + loMax) / 2, olat = (laMin + laMax) / 2;
  const buf = new ArrayBuffer(32 + 4 * (n + 1) + 8 * nv + 8 * n);
  const dv = new DataView(buf);
  [70, 65, 53, 77].forEach((c, i) => dv.setUint8(i, c));               // "FA5M"
  dv.setUint32(4, 1, true); dv.setUint32(8, n, true); dv.setUint32(12, nv, true);
  dv.setFloat64(16, olon, true); dv.setFloat64(24, olat, true);
  let o = 32, off = 0;
  dv.setUint32(o, 0, true); o += 4;
  for (let i = 0; i < n; i++) { off += len[i]; dv.setUint32(o, off, true); o += 4; }
  for (let i = 0; i < n; i++) for (let j = start[i]; j < start[i] + len[i]; j++) { dv.setFloat32(o, plon[j] - olon, true); o += 4; }
  for (let i = 0; i < n; i++) for (let j = start[i]; j < start[i] + len[i]; j++) { dv.setFloat32(o, plat[j] - olat, true); o += 4; }
  for (let i = 0; i < n; i++) { dv.setFloat32(o, elevation[i], true); o += 4; }
  for (let i = 0; i < n; i++) { dv.setFloat32(o, base[i], true); o += 4; }
  return buf;
}

/** depth_datum: SGS depth is above sgs_z_min (where finite), the other models' above the cell elevation. */
export function depthDatum(flowModel, elevation, stat) {
  const z = String(flowModel) === 'sgs' ? stat.sgs_z_min : null;
  return Array.from(elevation, (e, i) => (z && Number.isFinite(z[i]) ? z[i] : e));
}

// ---------------------------------------------------------------------------
// The run
// ---------------------------------------------------------------------------

export class LocalRun {
  /**
   * @param reader  the file (see the header)
   * @param name    the file's name (label, source)
   * @param viewer  the viewer page ("cesium" / "giro3d"), as the server's info
   */
  constructor(reader, { name = 'run.h5', viewer = 'cesium' } = {}) {
    this.r = reader;
    this.name = name;
    this.viewer = viewer;
    const r = reader;
    if (!r.has('/results')) {
      if (r.has('/frames')) throw new Error(`${name} is a format-1 file (FloodA5 ≤ 0.3.x): it has no cell polygons — open it with web_replay.jl --mesh`);
      throw new Error(`${name} is not a FloodA5 output file`);
    }
    this.rootAttrs = r.attrs('/');
    this.formatVersion = Number(this.rootAttrs.format_version ?? 2);
    this.t = Array.from(r.read('/results/t'), Number);
    if (!this.t.length) throw new Error(`${name} has no frames yet`);
    const res = new Set(r.keys('/results'));
    this.vars = WEB_VARS.filter(v => res.has(v));
    if (!this.vars.includes('water_depth')) throw new Error(`${name} has no water_depth frames`);
    if (!r.has('/mesh/polygon_offsets')) throw new Error(`${name} has no cell polygons`);
    // the mesh (FA.mesh_geometry, format 2)
    this.ids = Array.from(r.read('/mesh/cell_ids'), String);
    const n = this.nCells = this.ids.length;
    this.lon = toNum(r.read('/mesh/center_lons'));
    this.lat = toNum(r.read('/mesh/center_lats'));
    const elev = toNum(r.read('/mesh/elevations'));
    this.static = {};
    for (const k of r.keys('/mesh')) {
      if (RESERVED.includes(k) || !r.isDataset(`/mesh/${k}`)) continue;
      const sh = r.shape(`/mesh/${k}`);
      if (sh.length !== 1 || sh[0] !== n) continue;
      const v = r.read(`/mesh/${k}`);
      if (ArrayBuffer.isView(v) || (Array.isArray(v) && typeof v[0] === 'number')) this.static[k] = toNum(v);
    }
    const runAttrs = r.has('/run') ? r.attrs('/run') : {};
    let fm = runAttrs.flow_model;
    if (fm === undefined && r.has('/run/options')) fm = parseJsonInf(r.read('/run/options')).flow_model;
    this.flowModel = fm === undefined ? '' : String(fm);
    this.elevation = Float32Array.from(elev);
    this.base = Float32Array.from(depthDatum(this.flowModel, elev, this.static));
    this.meshBytes = encodeMesh(toNum(r.read('/mesh/polygon_offsets')), toNum(r.read('/mesh/polygon_lon')),
                                toNum(r.read('/mesh/polygon_lat')), this.elevation, this.base);
    this.idScheme = r.attrs('/mesh').id_scheme ?? '';
    this.cache = new Map();
    this.stats = new Map();
    this.ranges = new Map();
    // the model's statistics over the run (/statistics, 0.8.13), when the file has them
    this.fileStats = null;
    if (r.has('/statistics')) {
      const st = {};
      for (const k of r.keys('/statistics')) {
        if (r.isDataset(`/statistics/${k}`)) st[k] = toF32(r.read(`/statistics/${k}`));
      }
      this.fileStats = { stats: st, attrs: r.attrs('/statistics') };
    }
  }

  /** file_stat: a statistic from the file when it answers the request (see VisualisationServer.file_stat). */
  fileStat(name, p) {
    const fs = this.fileStats;
    const v = fs?.stats[name];
    if (!v || v.length !== this.nCells) return null;
    if ((name === 'time_first_flooded' || name === 'flood_duration') && Number(fs.attrs.flood_depth) !== p.flood_depth) return null;
    return v;
  }

  get nFrames() { return this.t.length; }
  get hasDerived() { return DERIVED_INPUTS.every(v => this.vars.includes(v)); }
  derivedNames() { return this.hasDerived ? HAZARD_NAMES.slice() : []; }
  statNames() {
    return this.hasDerived ? STAT_VARS.slice()
      : this.vars.includes('water_depth') ? ['max_water_depth', 'time_max_water_depth', 'time_first_flooded', 'flood_duration'] : [];
  }

  /** Frame `id` (1-based) of frame variable `v` as Float32 (cached), or null. */
  read(id, v) {
    if (!(id >= 1 && id <= this.t.length) || !this.vars.includes(v)) return null;
    const key = `${id}/${v}`;
    let a = this.cache.get(key);
    if (a) { this.cache.delete(key); this.cache.set(key, a); return a; }
    a = toF32(this.r.row(`/results/${v}`, id - 1));
    this.cache.set(key, a);
    while (this.cache.size > FRAME_CACHE) this.cache.delete(this.cache.keys().next().value);
    return a;
  }

  /** _frame_data: a frame variable or a hazard variable, Float32, or null. */
  frameData(id, v, p) {
    if (isDerived(v)) {
      const d = this.read(id, 'water_depth'), u = this.read(id, 'velocity');
      return d && u ? deriveHazard(v, d, u, p) : null;
    }
    return this.read(id, v);
  }

  /** series_data: one cell's (0-based) values over the frames, Float32, or null. */
  series(cell, v, p) {
    if (isDerived(v)) {
      const d = this.series(cell, 'water_depth', p), u = this.series(cell, 'velocity', p);
      return d && u ? deriveHazard(v, d, u, p) : null;
    }
    if (!this.vars.includes(v)) return null;
    return toF32(this.r.col(`/results/${v}`, cell));
  }

  info() {
    return {
      label: this.name, flow_model: this.flowModel, source: this.name, id_scheme: this.idScheme,
      flooda5_version: String(this.rootAttrs.flooda5_version ?? ''), viewer: this.viewer,
      protocol: PROTOCOL, mode: 'replay', complete: true, frames: this.t.length, vars: this.vars.slice(),
      derived: this.derivedNames(), stats: this.statNames(), mesh_ready: true, n_cells: this.nCells,
      local: true,
    };
  }

  frameList() {
    return { complete: true, frames: this.t.map((t, i) => ({ id: i + 1, t })) };
  }

  /** run_ranges (replay): 2nd/98th percentiles over wet cells of up to 24 frames. */
  runRanges(wet = 0.05, p = hazardParams({}), { maxFrames = 24, maxValues = 4_000_000 } = {}) {
    const key = `${wet}|${hazardKey(p)}`;
    if (this.ranges.has(key)) return this.ranges.get(key);
    const n = this.t.length, nc = this.nCells;
    const picks = pickFrames(n, maxFrames);
    const names = this.vars, derived = this.derivedNames();
    const vals = {};
    for (const v of [...names, ...derived, 'level']) vals[v] = [];
    const stride = Math.max(1, Math.ceil(nc * Math.max(picks.length, 1) / maxValues));
    const base = this.base;
    for (const k of picks) {
      const d = this.read(k, 'water_depth');
      if (!d) continue;
      const idx = [];
      for (let i = 0; i < Math.min(d.length, nc); i += stride) if (d[i] >= wet) idx.push(i);
      let vel = null;
      for (const v of names) {
        const x = v === 'water_depth' ? d : this.read(k, v);
        if (v === 'velocity') vel = x;
        if (x) for (const i of idx) if (Number.isFinite(x[i])) vals[v].push(x[i]);
      }
      if (derived.length && vel) {
        const dw = Float32Array.from(idx, i => d[i]), vw = Float32Array.from(idx, i => vel[i]);
        for (const v of derived) for (const y of deriveHazard(v, dw, vw, p)) if (Number.isFinite(y)) vals[v].push(y);
      }
      for (const i of idx) { const y = Math.fround(base[i] + d[i]); if (Number.isFinite(y)) vals.level.push(y); }
    }
    const out = {};
    for (const [v, list] of Object.entries(vals)) {
      const a = Float32Array.from(list).sort();
      out[v] = { p02: fl(sortedQuantile(a, 0.02)), p98: fl(sortedQuantile(a, 0.98)),
                 min: a.length ? fl(a[0]) : null, max: a.length ? fl(a[a.length - 1]) : null, n: a.length };
    }
    const r = { wet, frames: picks.length, of: n, complete: true, vars: out };
    this.ranges.set(key, r);
    return r;
  }

  /** run_stats: statistics over every frame (cached per parameter set). */
  runStats(p) {
    const key = hazardKey(p);
    if (this.stats.has(key)) return this.stats.get(key);
    const hv = this.vars.includes('velocity');
    // read each frame once, bypassing the cache (a run can be longer than it)
    const rd = (k, v) => toF32(this.r.row(`/results/${v}`, k));
    const st = runStatistics(this.t, k => [rd(k, 'water_depth'), hv ? rd(k, 'velocity') : null], this.nCells, p);
    if (this.stats.size >= 4) this.stats.clear();
    this.stats.set(key, st);
    return st;
  }

  cell(i, q) {
    if (!(Number.isInteger(i) && i >= 0 && i < this.nCells)) return null;
    const vals = {};
    const fid = parseIntStrict(q.frame);
    if (fid !== null) {
      for (const v of this.vars) {
        const d = this.read(fid, v);
        if (d) vals[v] = fl32(d[i]);
      }
      const p = hazardParams(q);
      for (const v of this.derivedNames()) {
        const d = this.read(fid, 'water_depth'), u = this.read(fid, 'velocity');
        if (d && u) vals[v] = fl32(deriveHazard(v, d.subarray(i, i + 1), u.subarray(i, i + 1), p)[0]);
      }
    }
    return { index: i, id: this.ids[i], lon: this.lon[i], lat: this.lat[i], elevation: fl32(this.elevation[i]),
             base: fl32(this.base[i]), frame: fid, values: vals };
  }

  /** FO.run_record(f) plus the file and label (web_replay.jl). */
  runRecord() {
    const r = this.r;
    const info = { format_version: this.formatVersion };
    Object.assign(info, this.rootAttrs);
    if (r.has('/run')) {
      Object.assign(info, r.attrs('/run'));
      info.argv = r.has('/run/argv') ? parseJsonInf(r.read('/run/argv')).map(String) : [];
      info.options = r.has('/run/options') ? parseJsonInf(r.read('/run/options')) : {};
      for (const k of ['git_diff', 'project', 'manifest']) info[k] = r.has(`/run/${k}`) ? String(r.read(`/run/${k}`)) : '';
      info.inputs = r.has('/run/inputs') ? r.keys('/run/inputs').slice().sort().map(k => r.attrs(`/run/inputs/${k}`)) : [];
    }
    info.mesh_source = r.has('/mesh/source') ? r.attrs('/mesh/source') : {};
    info.summary = r.has('/summary') ? r.attrs('/summary') : {};
    if (r.has('/statistics')) info.statistics = r.attrs('/statistics');
    // run_record(info; mesh, n_cells, cell_area, times)
    const rec = {};
    for (const [k, v] of Object.entries(info)) if (!['project', 'manifest', 'git_diff'].includes(k)) rec[k] = v;
    const gd = info.git_diff;
    if (typeof gd === 'string' && gd.length) rec.git_diff_lines = (gd.match(/\n/g) || []).length;
    const mesh = { ...r.attrs('/mesh'), n_cells: this.nCells };
    if (r.has('/mesh/cell_area')) {
      const a = Array.from(r.read('/mesh/cell_area'), Number).filter(Number.isFinite);
      if (a.length) {
        let s = 0, lo = Infinity, hi = -Infinity;
        for (const x of a) { s += x; if (x < lo) lo = x; if (x > hi) hi = x; }
        mesh.cell_area_min_m2 = lo; mesh.cell_area_max_m2 = hi; mesh.cell_area_mean_m2 = s / a.length;
      }
    }
    rec.mesh = mesh;
    rec.frames = { n: this.t.length, t_first: this.t[0] ?? null, t_last: this.t[this.t.length - 1] ?? null };
    const out = plain(rec);
    out.file = this.name;
    out.label = this.name;
    return out;
  }

  /**
   * One request, as the server answers it: `path` with its query
   * ("/api/frames/3/water_depth?landuse=urban"). → {status, json} |
   * {status, bytes: ArrayBuffer} | {status, text}.
   */
  handle(path) {
    const u = new URL(path, 'http://local/');
    const parts = u.pathname.split('/').filter(Boolean).map(decodeURIComponent);
    const q = Object.fromEntries(u.searchParams);
    const json = x => ({ status: 200, json: x });
    const text = (status, t) => ({ status, text: t });
    const bin = a => {
      const b = a instanceof Float32Array ? a : Float32Array.from(a);
      return { status: 200, bytes: b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) };
    };
    if (parts[0] !== 'api' || parts.length < 2) return text(404, 'Not found');
    const a = parts[1];
    if (a === 'info' && parts.length === 2) return json(this.info());
    if (a === 'run' && parts.length === 2) return json(this.runRecord());
    if (a === 'mesh' && parts.length === 2) return { status: 200, bytes: this.meshBytes.slice(0) };
    if (a === 'frames') {
      if (parts.length === 2) return json(this.frameList());
      if (parts.length === 4) {
        const id = parseIntStrict(parts[2]);
        if (id === null) return text(400, 'Bad frame id');
        const v = this.frameData(id, parts[3], hazardParams(q));
        return v ? bin(v) : text(404, `No frame ${id} variable ${parts[3]}`);
      }
    }
    if (a === 'ranges' && parts.length === 2) {
      let w = q.wet === undefined || q.wet === '' ? null : Number(q.wet);
      if (w === null || Number.isNaN(w)) w = 0.05;
      w = Number.isFinite(w) ? Math.max(w, 0) : 0.05;
      return json(this.runRanges(w, hazardParams(q)));
    }
    if (a === 'stats' && parts.length === 3) {
      const name = parts[2];
      if (!STAT_VARS.includes(name)) return text(404, `No statistic ${name}`);
      const p = hazardParams(q);
      const v = this.fileStat(name, p) || this.runStats(p)[name];
      return v ? bin(v) : text(404, `No statistic ${name} for this run`);
    }
    if (a === 'cell' && parts.length === 3) {
      const c = this.cell(parseIntStrict(parts[2]), q);
      return c ? json(c) : text(404, 'No such cell');
    }
    if (a === 'series' && parts.length === 3) {
      const i = parseIntStrict(parts[2]);
      if (!(i !== null && i >= 0 && i < this.nCells)) return text(404, 'No such cell');
      const v = q.var ?? 'water_depth';
      const s = this.series(i, v, hazardParams(q));
      if (!s) return text(404, `No variable ${v}`);
      return json({ index: i, id: this.ids[i], var: v, frames: this.t.map((_, k) => k + 1), t: this.t.slice(),
                    values: Array.from(s, fl32) });
    }
    return text(404, 'Not found');
  }
}

// Julia's tryparse(Int, s): an optional sign and digits only
function parseIntStrict(s) {
  if (s === undefined || s === null || !/^[+-]?\d+$/.test(String(s))) return null;
  return Number(s);
}
