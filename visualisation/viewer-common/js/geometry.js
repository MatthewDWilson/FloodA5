// geometry.js — mesh geometry for the globe viewer. Pure (no Cesium, no DOM):
// unit-tested with `node --test`.
//
// The water (and bed) surface is one triangle mesh for the whole domain, each
// cell with its own vertices so that it can sit at its own level. Vertices are
// stored relative to an ECEF centre (Float32 offsets keep ~1 mm precision over
// tens of kilometres); the per-frame level moves them along the ellipsoid
// normal in the vertex shader (water.js).

export const WGS84_A = 6378137.0;
export const WGS84_F = 1 / 298.257223563;
export const WGS84_E2 = WGS84_F * (2 - WGS84_F);
export const WGS84_B = WGS84_A * (1 - WGS84_F);
const DEG = Math.PI / 180;

/** Geodetic (degrees, metres above the ellipsoid) → ECEF [x, y, z] (metres). */
export function geodeticToEcef(lon, lat, h = 0) {
  const lo = lon * DEG, la = lat * DEG;
  const s = Math.sin(la), c = Math.cos(la);
  const N = WGS84_A / Math.sqrt(1 - WGS84_E2 * s * s);
  return [(N + h) * c * Math.cos(lo), (N + h) * c * Math.sin(lo), (N * (1 - WGS84_E2) + h) * s];
}

/**
 * ECEF (metres) → geodetic {lon, lat (degrees), h (metres above the ellipsoid)}.
 * Bowring's method (with the second eccentricity in the latitude numerator),
 * one step: better than 1 mm for heights within ±10 km of the surface.
 * (Giro3D 2.0.7's Ellipsoid.toGeodetic uses the first eccentricity there and
 * is ~130 m off in height at mid-latitudes: use this instead.)
 */
export function ecefToGeodetic(x, y, z) {
  const a = WGS84_A, b = WGS84_B, e2 = WGS84_E2, ep2 = (a * a - b * b) / (b * b);
  const p = Math.hypot(x, y);
  const th = Math.atan2(z * a, p * b);
  const st = Math.sin(th), ct = Math.cos(th);
  const lat = Math.atan2(z + ep2 * b * st * st * st, p - e2 * a * ct * ct * ct);
  const sl = Math.sin(lat), cl = Math.cos(lat);
  const N = a / Math.sqrt(1 - e2 * sl * sl);
  const h = Math.abs(cl) > 1e-10 ? p / cl - N : Math.abs(z) / Math.abs(sl) - N * (1 - e2);
  return { lon: Math.atan2(y, x) / DEG, lat: lat / DEG, h };
}

/** Absolute vertex longitudes/latitudes (Float64) of a parsed mesh. */
export function vertexLonLat(mesh) {
  const n = mesh.nVertices, lon = new Float64Array(n), lat = new Float64Array(n);
  for (let v = 0; v < n; v++) {
    lon[v] = mesh.originLon + mesh.dlon[v];
    lat[v] = mesh.originLat + mesh.dlat[v];
  }
  return { lon, lat };
}

/** Cell centres (vertex means), Float64. */
export function cellCentres(mesh) {
  const { lon, lat } = vertexLonLat(mesh);
  const cl = new Float64Array(mesh.nCells), ca = new Float64Array(mesh.nCells);
  for (let i = 0; i < mesh.nCells; i++) {
    const a = mesh.offsets[i], b = mesh.offsets[i + 1];
    let sx = 0, sy = 0;
    for (let v = a; v < b; v++) { sx += lon[v]; sy += lat[v]; }
    cl[i] = sx / (b - a); ca[i] = sy / (b - a);
  }
  return { lon: cl, lat: ca };
}

/** Median of the finite values (NaN if none). */
export function finiteMedian(values) {
  const a = Array.from(values).filter(Number.isFinite).sort((x, y) => x - y);
  if (a.length === 0) return NaN;
  const m = a.length >> 1;
  return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
}

/** q-quantile (0–1) of the finite values (NaN if none). */
export function finiteQuantile(values, q) {
  const a = Array.from(values).filter(Number.isFinite).sort((x, y) => x - y);
  if (a.length === 0) return NaN;
  const p = Math.min(Math.max(q, 0), 1) * (a.length - 1);
  const lo = Math.floor(p), hi = Math.ceil(p);
  return a[lo] + (a[hi] - a[lo]) * (p - lo);
}

/**
 * Triangulate the mesh for the surface layers.
 * Cells are convex (A5 pentagons, squares, hexagons): a fan from the first
 * vertex. Vertices are placed on the ellipsoid at height h0 (choose it near
 * the water level so the bounding volume is right) relative to `centre`.
 * @returns {{positions:Float32Array, cell:Float32Array, indices:Uint32Array,
 *            centre:number[], h0:number, min:number[], max:number[]}}
 */
export function buildSurface(mesh, h0 = 0) {
  const centre = geodeticToEcef(mesh.originLon, mesh.originLat, h0);
  const nv = mesh.nVertices;
  const positions = new Float32Array(3 * nv);
  const cell = new Float32Array(nv);
  let ntri = 0;
  for (let i = 0; i < mesh.nCells; i++) ntri += mesh.offsets[i + 1] - mesh.offsets[i] - 2;
  const indices = new Uint32Array(3 * ntri);
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  let t = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    const a = mesh.offsets[i], b = mesh.offsets[i + 1];
    for (let v = a; v < b; v++) {
      const p = geodeticToEcef(mesh.originLon + mesh.dlon[v], mesh.originLat + mesh.dlat[v], h0);
      for (let k = 0; k < 3; k++) {
        const x = Math.fround(p[k] - centre[k]);
        positions[3 * v + k] = x;
        if (x < min[k]) min[k] = x;
        if (x > max[k]) max[k] = x;
      }
      cell[v] = i;
    }
    for (let v = a + 1; v < b - 1; v++) {
      indices[t++] = a; indices[t++] = v; indices[t++] = v + 1;
    }
  }
  return { positions, cell, indices, centre, h0, min, max };
}

/**
 * A binary glTF 2.0 (GLB) of the surface: POSITION, a custom per-vertex
 * attribute `_CELL` (the cell index, float — exact below 2^24 cells),
 * UNSIGNED_INT indices and one unlit, double-sided material. The glTF is
 * Z-up (load it with upAxis Z, forwardAxis X so no axis correction applies).
 * @returns {Uint8Array}
 */
export function buildGlb(surface) {
  const { positions, cell, indices, min, max } = surface;
  const pad4 = n => (n + 3) & ~3;
  const posBytes = positions.byteLength, cellBytes = cell.byteLength, idxBytes = indices.byteLength;
  const posOff = 0, cellOff = pad4(posOff + posBytes), idxOff = pad4(cellOff + cellBytes);
  const binLength = pad4(idxOff + idxBytes);
  const nv = positions.length / 3;
  const json = {
    asset: { version: '2.0', generator: 'FloodA5 viewer' },
    extensionsUsed: ['KHR_materials_unlit'],
    scene: 0,
    scenes: [{ nodes: [0] }],
    nodes: [{ mesh: 0 }],
    meshes: [{ primitives: [{ attributes: { POSITION: 0, _CELL: 1 }, indices: 2, material: 0, mode: 4 }] }],
    materials: [{
      pbrMetallicRoughness: { baseColorFactor: [1, 1, 1, 1], metallicFactor: 0, roughnessFactor: 1 },
      doubleSided: true,
      alphaMode: 'BLEND',
      extensions: { KHR_materials_unlit: {} },
    }],
    buffers: [{ byteLength: binLength }],
    bufferViews: [
      { buffer: 0, byteOffset: posOff, byteLength: posBytes, target: 34962 },
      { buffer: 0, byteOffset: cellOff, byteLength: cellBytes, target: 34962 },
      { buffer: 0, byteOffset: idxOff, byteLength: idxBytes, target: 34963 },
    ],
    accessors: [
      { bufferView: 0, componentType: 5126, count: nv, type: 'VEC3', min, max },
      { bufferView: 1, componentType: 5126, count: nv, type: 'SCALAR' },
      { bufferView: 2, componentType: 5125, count: indices.length, type: 'SCALAR' },
    ],
  };
  let jsonText = JSON.stringify(json);
  while (jsonText.length % 4) jsonText += ' ';
  const jsonBytes = new TextEncoder().encode(jsonText);
  const total = 12 + 8 + jsonBytes.length + 8 + binLength;
  const out = new Uint8Array(total);
  const dv = new DataView(out.buffer);
  dv.setUint32(0, 0x46546c67, true);           // "glTF"
  dv.setUint32(4, 2, true);
  dv.setUint32(8, total, true);
  dv.setUint32(12, jsonBytes.length, true);
  dv.setUint32(16, 0x4e4f534a, true);          // "JSON"
  out.set(jsonBytes, 20);
  const b0 = 20 + jsonBytes.length;
  dv.setUint32(b0, binLength, true);
  dv.setUint32(b0 + 4, 0x004e4942, true);      // "BIN\0"
  const bin = b0 + 8;
  out.set(new Uint8Array(positions.buffer, positions.byteOffset, posBytes), bin + posOff);
  out.set(new Uint8Array(cell.buffer, cell.byteOffset, cellBytes), bin + cellOff);
  out.set(new Uint8Array(indices.buffer, indices.byteOffset, idxBytes), bin + idxOff);
  return out;
}

/** Texture size for n values packed four per RGBA texel (≤ maxWidth texels wide). */
export function packedSize(n, perTexel = 4, maxWidth = 2048) {
  const texels = Math.max(1, Math.ceil(n / perTexel));
  const width = Math.min(maxWidth, texels);
  return { width, height: Math.ceil(texels / width) };
}

/**
 * Pack per-cell values four to an RGBA32F texel (cell i → texel ⌊i/4⌋,
 * component i mod 4), the layout the shader reads. Non-finite values → 0.
 */
export function packCellValues(values, n, out) {
  const { width, height } = packedSize(n);
  const len = 4 * width * height;
  const data = out && out.length === len ? out : new Float32Array(len);
  data.fill(0);
  for (let i = 0; i < n; i++) {
    const v = values[i];
    data[i] = Number.isFinite(v) ? v : 0;
  }
  return { data, width, height };
}

/**
 * Per-cell static values, one cell per RGBA32F texel: (bed, base, 0, 0).
 * Missing values are filled (bed ← base ← bed ← median) so that no vertex
 * gets a NaN height.
 */
export function packStatic(elevation, base, n) {
  const { width, height } = packedSize(n, 1);
  const data = new Float32Array(4 * width * height);
  const fill = finiteMedian(Number.isFinite(finiteMedian(base)) ? base : elevation);
  const dflt = Number.isFinite(fill) ? fill : 0;
  for (let i = 0; i < n; i++) {
    let e = elevation[i], b = base[i];
    if (!Number.isFinite(b)) b = Number.isFinite(e) ? e : dflt;
    if (!Number.isFinite(e)) e = b;
    data[4 * i] = e; data[4 * i + 1] = b;
  }
  return { data, width, height };
}

/**
 * Locate the cell containing a point: a uniform grid over the cell centres,
 * then point-in-polygon tests on the nearby cells.
 */
export class CellLocator {
  constructor(mesh) {
    this.mesh = mesh;
    const { lon, lat } = cellCentres(mesh);
    this.cLon = lon; this.cLat = lat;
    const [x0, y0, x1, y1] = meshBounds(mesh);   // vertex extents: never degenerate
    // Square buckets of about two cells each (also for long, thin domains):
    // the ±2-bucket search then always reaches the containing cell's centre.
    const w = Math.max(x1 - x0, 1e-9), h = Math.max(y1 - y0, 1e-9);
    const size = Math.max(Math.sqrt((w * h) / Math.max(mesh.nCells / 2, 1)), 1e-9);
    this.x0 = x0; this.y0 = y0;
    this.dx = Math.min(size, w); this.dy = Math.min(size, h);
    this.nx = Math.ceil(w / this.dx) + 1; this.ny = Math.ceil(h / this.dy) + 1;
    this.buckets = new Map();
    for (let i = 0; i < mesh.nCells; i++) {
      const k = this._key(this._bx(lon[i]), this._by(lat[i]));
      let b = this.buckets.get(k);
      if (!b) { b = []; this.buckets.set(k, b); }
      b.push(i);
    }
  }
  _bx(x) { return Math.min(this.nx - 1, Math.max(0, Math.floor((x - this.x0) / this.dx))); }
  _by(y) { return Math.min(this.ny - 1, Math.max(0, Math.floor((y - this.y0) / this.dy))); }
  _key(i, j) { return j * this.nx + i; }

  /** Index of the cell containing (lon, lat), or -1. */
  find(lon, lat) {
    const bi = this._bx(lon), bj = this._by(lat);
    for (let j = bj - 2; j <= bj + 2; j++) {
      for (let i = bi - 2; i <= bi + 2; i++) {
        if (i < 0 || j < 0 || i >= this.nx || j >= this.ny) continue;
        const b = this.buckets.get(this._key(i, j));
        if (!b) continue;
        for (const c of b) if (this.contains(c, lon, lat)) return c;
      }
    }
    return -1;
  }

  /** Whether cell c's polygon contains (lon, lat) (even-odd rule). */
  contains(c, lon, lat) {
    const m = this.mesh, a = m.offsets[c], b = m.offsets[c + 1];
    const x = lon - m.originLon, y = lat - m.originLat;
    let inside = false;
    for (let v = a, w = b - 1; v < b; w = v++) {
      const xv = m.dlon[v], yv = m.dlat[v], xw = m.dlon[w], yw = m.dlat[w];
      if ((yv > y) !== (yw > y) && x < ((xw - xv) * (y - yv)) / (yw - yv) + xv) inside = !inside;
    }
    return inside;
  }
}

/** Bounding box of the mesh vertices: [west, south, east, north] (degrees). */
export function meshBounds(mesh) {
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (let v = 0; v < mesh.nVertices; v++) {
    x0 = Math.min(x0, mesh.dlon[v]); x1 = Math.max(x1, mesh.dlon[v]);
    y0 = Math.min(y0, mesh.dlat[v]); y1 = Math.max(y1, mesh.dlat[v]);
  }
  return [mesh.originLon + x0, mesh.originLat + y0, mesh.originLon + x1, mesh.originLat + y1];
}
