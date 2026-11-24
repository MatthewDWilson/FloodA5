// protocol.js — the viewer side of the FloodA5 visualisation protocol (version 2).
// Pure functions (no Cesium, no DOM): unit-tested with `node --test`.
// Layouts: visualisation/VisualisationServer.jl and docs/DATA_FORMATS.md §6.

export const PROTOCOL = 2;
const MAGIC = 'FA5M';

/**
 * Decode the binary mesh (`GET /api/mesh`).
 * @param {ArrayBuffer} buffer
 * @returns {{nCells:number, nVertices:number, originLon:number, originLat:number,
 *            offsets:Uint32Array, dlon:Float32Array, dlat:Float32Array,
 *            elevation:Float32Array, base:Float32Array}}
 */
export function parseMesh(buffer) {
  if (!(buffer instanceof ArrayBuffer)) throw new Error('parseMesh: expected an ArrayBuffer');
  if (buffer.byteLength < 32) throw new Error('mesh payload too short');
  const dv = new DataView(buffer);
  const magic = String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3));
  if (magic !== MAGIC) throw new Error(`not a FloodA5 mesh (magic "${magic}")`);
  const version = dv.getUint32(4, true);
  if (version !== 1) throw new Error(`unsupported mesh version ${version}`);
  const nCells = dv.getUint32(8, true);
  const nVertices = dv.getUint32(12, true);
  const originLon = dv.getFloat64(16, true);
  const originLat = dv.getFloat64(24, true);
  let o = 32;
  const need = 32 + 4 * (nCells + 1) + 8 * nVertices + 8 * nCells;
  if (buffer.byteLength !== need) {
    throw new Error(`mesh payload is ${buffer.byteLength} bytes, expected ${need}`);
  }
  const offsets = u32(buffer, o, nCells + 1); o += 4 * (nCells + 1);
  const dlon = f32(buffer, o, nVertices); o += 4 * nVertices;
  const dlat = f32(buffer, o, nVertices); o += 4 * nVertices;
  const elevation = f32(buffer, o, nCells); o += 4 * nCells;
  const base = f32(buffer, o, nCells);
  if (offsets[0] !== 0 || offsets[nCells] !== nVertices) throw new Error('mesh offsets are inconsistent');
  for (let i = 0; i < nCells; i++) {
    if (offsets[i + 1] - offsets[i] < 3) throw new Error(`mesh cell ${i} has fewer than 3 vertices`);
  }
  return { nCells, nVertices, originLon, originLat, offsets, dlon, dlat, elevation, base };
}

const LITTLE = new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;

function u32(buffer, offset, n) {
  if (LITTLE) return new Uint32Array(buffer, offset, n);
  const dv = new DataView(buffer, offset), a = new Uint32Array(n);
  for (let i = 0; i < n; i++) a[i] = dv.getUint32(4 * i, true);
  return a;
}

function f32(buffer, offset, n) {
  if (LITTLE) return new Float32Array(buffer, offset, n);
  return decodeFloat32(buffer.slice(offset, offset + 4 * n));
}

/** Decode a frame variable (`GET /api/frames/{id}/{var}`): little-endian Float32. */
export function decodeFloat32(buffer, nCells) {
  if (buffer.byteLength % 4 !== 0) throw new Error('frame payload is not a Float32 array');
  let a;
  if (LITTLE) a = new Float32Array(buffer);
  else {
    const dv = new DataView(buffer); a = new Float32Array(buffer.byteLength / 4);
    for (let i = 0; i < a.length; i++) a[i] = dv.getFloat32(4 * i, true);
  }
  if (nCells !== undefined && a.length !== nCells) {
    throw new Error(`frame has ${a.length} values, the mesh ${nCells} cells`);
  }
  return a;
}

/**
 * The frames known to the viewer: [{id, t}] in time order, with lookups.
 */
export class FrameList {
  constructor(frames = []) { this.set(frames); }
  set(frames) {
    this.frames = frames.map(f => ({ id: f.id, t: f.t })).sort((a, b) => a.t - b.t || a.id - b.id);
    this._pos = new Map(this.frames.map((f, i) => [f.id, i]));
  }
  get length() { return this.frames.length; }
  at(i) { return this.frames[i]; }
  last() { return this.frames[this.frames.length - 1]; }
  indexOfId(id) { const i = this._pos.get(id); return i === undefined ? -1 : i; }
  /** Append a frame announced live (ignored if already known). */
  add(id, t) {
    if (this._pos.has(id)) return false;
    const last = this.last();
    if (!last || t >= last.t) {
      this.frames.push({ id, t }); this._pos.set(id, this.frames.length - 1);
    } else {
      this.set([...this.frames, { id, t }]);
    }
    return true;
  }
  /** Index of the last frame with t <= time (0 if before the first). */
  indexAtTime(time) {
    const f = this.frames;
    if (f.length === 0) return -1;
    let lo = 0, hi = f.length - 1;
    if (time < f[0].t) return 0;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (f[mid].t <= time) lo = mid; else hi = mid - 1;
    }
    return lo;
  }
}

/** Small LRU cache of frame arrays keyed by `${id}/${var}`. */
export class FrameCache {
  constructor(capacity = 24) { this.capacity = capacity; this.map = new Map(); }
  get(key) {
    const v = this.map.get(key);
    if (v !== undefined) { this.map.delete(key); this.map.set(key, v); }
    return v;
  }
  set(key, value) {
    this.map.delete(key); this.map.set(key, value);
    while (this.map.size > this.capacity) this.map.delete(this.map.keys().next().value);
  }
  has(key) { return this.map.has(key); }
  clear() { this.map.clear(); }
}

/** Format a simulation time in seconds for display. */
export function formatTime(s) {
  if (!Number.isFinite(s)) return '–';
  const sign = s < 0 ? '-' : '';
  s = Math.abs(s);
  const h = Math.floor(s / 3600), m = Math.floor((s % 3600) / 60), sec = s % 60;
  if (h > 0) return `${sign}${h}:${String(m).padStart(2, '0')}:${String(Math.floor(sec)).padStart(2, '0')}`;
  if (m > 0) return `${sign}${m}:${String(Math.floor(sec)).padStart(2, '0')}`;
  return `${sign}${sec.toFixed(sec < 10 ? 1 : 0)} s`;
}
