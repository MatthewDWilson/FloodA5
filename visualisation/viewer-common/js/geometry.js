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
 * Each cell is a fan from its centre (the vertex mean; the cells are convex
 * or star-shaped about it: A5 pentagons, squares, hexagons, Voronoi cells)
 * on its own vertices: vertex `offsets[i] + i` is cell i's centre, followed by
 * its ring. `edge` is 1 at the centre and 0 on the ring, so in the fragment
 * shader it falls to 0 exactly on the cell's boundary (the cell-border lines).
 * Vertices are placed on the ellipsoid at height h0 (choose it near the water
 * level so the bounding volume is right) relative to `centre`.
 * @returns {{positions:Float32Array, cell:Float32Array, edge:Float32Array,
 *            indices:Uint32Array, centre:number[], h0:number, min:number[], max:number[]}}
 */
export function buildSurface(mesh, h0 = 0) {
  const centre = geodeticToEcef(mesh.originLon, mesh.originLat, h0);
  const nv = mesh.nVertices + mesh.nCells;
  const positions = new Float32Array(3 * nv);
  const cell = new Float32Array(nv);
  const edge = new Float32Array(nv);
  const indices = new Uint32Array(3 * mesh.nVertices);     // one triangle per cell side
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  const put = (o, lon, lat) => {
    const p = geodeticToEcef(lon, lat, h0);
    for (let k = 0; k < 3; k++) {
      const x = Math.fround(p[k] - centre[k]);
      positions[3 * o + k] = x;
      if (x < min[k]) min[k] = x;
      if (x > max[k]) max[k] = x;
    }
  };
  let t = 0;
  for (let i = 0; i < mesh.nCells; i++) {
    const a = mesh.offsets[i], b = mesh.offsets[i + 1], c = a + i;   // c: the centre's vertex
    let sx = 0, sy = 0;
    for (let v = a; v < b; v++) { sx += mesh.dlon[v]; sy += mesh.dlat[v]; }
    put(c, mesh.originLon + sx / (b - a), mesh.originLat + sy / (b - a));
    cell[c] = i; edge[c] = 1;
    for (let v = a; v < b; v++) {
      const o = v + i + 1;
      put(o, mesh.originLon + mesh.dlon[v], mesh.originLat + mesh.dlat[v]);
      cell[o] = i;
      indices[t++] = c; indices[t++] = o; indices[t++] = v + 1 < b ? o + 1 : c + 1;
    }
  }
  return { positions, cell, edge, indices, centre, h0, min, max };
}

/**
 * A binary glTF 2.0 (GLB) of the surface: POSITION, custom per-vertex
 * attributes `_CELL` (the cell index, float — exact below 2^24 cells) and
 * `_FA5EDGE` (buildSurface's `edge`; a name of our own, since CesiumJS
 * already declares a `v_edge` of its own), UNSIGNED_INT indices and one unlit,
 * double-sided material. The glTF is Z-up (load it with upAxis Z, forwardAxis
 * X so no axis correction applies).
 * @returns {Uint8Array}
 */
export function buildGlb(surface, corner = null) {
  const { positions, cell, edge, indices, min, max } = surface;
  const pad4 = n => (n + 3) & ~3;
  const posBytes = positions.byteLength, cellBytes = cell.byteLength, edgeBytes = edge.byteLength,
        idxBytes = indices.byteLength;
  const cornerBytes = corner ? corner.byteLength : 0;
  const posOff = 0, cellOff = pad4(posOff + posBytes), edgeOff = pad4(cellOff + cellBytes),
        idxOff = pad4(edgeOff + edgeBytes), cornerOff = pad4(idxOff + idxBytes);
  const binLength = pad4(cornerOff + cornerBytes);
  const nv = positions.length / 3;
  const json = {
    asset: { version: '2.0', generator: 'FloodA5 viewer' },
    extensionsUsed: ['KHR_materials_unlit'],
    scene: 0,
    scenes: [{ nodes: [0] }],
    nodes: [{ mesh: 0 }],
    meshes: [{ primitives: [{ attributes: { POSITION: 0, _CELL: 1, _FA5EDGE: 2 }, indices: 3, material: 0, mode: 4 }] }],
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
      { buffer: 0, byteOffset: edgeOff, byteLength: edgeBytes, target: 34962 },
      { buffer: 0, byteOffset: idxOff, byteLength: idxBytes, target: 34963 },
    ],
    accessors: [
      { bufferView: 0, componentType: 5126, count: nv, type: 'VEC3', min, max },
      { bufferView: 1, componentType: 5126, count: nv, type: 'SCALAR' },
      { bufferView: 2, componentType: 5126, count: nv, type: 'SCALAR' },
      { bufferView: 3, componentType: 5125, count: indices.length, type: 'SCALAR' },
    ],
  };
  if (corner) {                       // the smooth surface (0.8.14): each vertex's corner (cornerIndex)
    json.meshes[0].primitives[0].attributes._FA5CORNER = 4;
    json.bufferViews.push({ buffer: 0, byteOffset: cornerOff, byteLength: cornerBytes, target: 34962 });
    json.accessors.push({ bufferView: 4, componentType: 5126, count: nv, type: 'SCALAR' });
  }
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
  out.set(new Uint8Array(edge.buffer, edge.byteOffset, edgeBytes), bin + edgeOff);
  out.set(new Uint8Array(indices.buffer, indices.byteOffset, idxBytes), bin + idxOff);
  if (corner) out.set(new Uint8Array(corner.buffer, corner.byteOffset, cornerBytes), bin + cornerOff);
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

// ---------------------------------------------------------------------------
// Smooth surface (0.8.14): the cells' shared corners
// ---------------------------------------------------------------------------

/**
 * The corners the cells share. Ring vertices at the same place (to `quantum`
 * degrees, ~1 cm by default) are one corner; found by sorting, not hashing,
 * so meshes of millions of cells stay affordable. Built when the smooth
 * surface is first asked for.
 * @returns {{vertexCorner: Float32Array, nCorners: number, offsets: Uint32Array, cells: Uint32Array}}
 *   vertexCorner: for each surface vertex (buildSurface's order) its corner,
 *   −1 at the cell centres; corner c's cells are cells[offsets[c] … offsets[c+1]).
 */
export function cornerIndex(mesh, quantum = 1e-7) {
  const nv = mesh.nVertices, nc = mesh.nCells;
  let ext = 0;
  for (let v = 0; v < nv; v++) ext = Math.max(ext, Math.abs(mesh.dlon[v]), Math.abs(mesh.dlat[v]));
  const q = Math.max(quantum, ext / 2 ** 24);            // |key parts| < 2^25: the key is exact in a double
  const OFF = 2 ** 25, MUL = 2 ** 26;
  const keys = new Float64Array(nv);
  for (let v = 0; v < nv; v++) {
    keys[v] = (Math.round(mesh.dlon[v] / q) + OFF) * MUL + (Math.round(mesh.dlat[v] / q) + OFF);
  }
  const sorted = Float64Array.from(keys).sort();
  let nCorners = 0;
  for (let k = 0; k < nv; k++) if (k === 0 || sorted[k] !== sorted[k - 1]) sorted[nCorners++] = sorted[k];
  const uniq = sorted.subarray(0, nCorners);
  const ringCorner = new Uint32Array(nv);
  for (let v = 0; v < nv; v++) {
    let lo = 0, hi = nCorners - 1;
    const key = keys[v];
    while (lo < hi) { const mid = (lo + hi) >> 1; if (uniq[mid] < key) lo = mid + 1; else hi = mid; }
    ringCorner[v] = lo;
  }
  // number the corners in the order the cells first use them (0.9.2): a cell's corners then
  // sit near each other in memory, which more than halves the per-frame corner work
  const renum = new Int32Array(nCorners).fill(-1);
  let next = 0;
  for (let v = 0; v < nv; v++) {
    if (renum[ringCorner[v]] < 0) renum[ringCorner[v]] = next++;
    ringCorner[v] = renum[ringCorner[v]];
  }
  // corner → cells (CSR); vertex → corner in the surface's vertex order
  const offsets = new Uint32Array(nCorners + 1);
  const vertexCorner = new Float32Array(nv + nc).fill(-1);
  for (let i = 0; i < nc; i++) {
    for (let v = mesh.offsets[i]; v < mesh.offsets[i + 1]; v++) {
      offsets[ringCorner[v] + 1]++;
      vertexCorner[v + i + 1] = ringCorner[v];
    }
  }
  for (let c = 0; c < nCorners; c++) offsets[c + 1] += offsets[c];
  const fillAt = offsets.slice(0, nCorners), cells = new Uint32Array(nv);
  for (let i = 0; i < nc; i++) {
    for (let v = mesh.offsets[i]; v < mesh.offsets[i + 1]; v++) cells[fillAt[ringCorner[v]]++] = i;
  }
  return { vertexCorner, nCorners, offsets, cells };
}

/**
 * The water at each corner: the means over the corner's wet cells (depth ≥
 * `minDepth`) of (level = base + depth, bed, depth, colour value), one corner
 * per RGBA32F texel (packedSize(nCorners, 1)). A corner with no wet cell gets
 * zeros (its cells are not drawn). `cvals` (optional): the per-cell colour
 * values; non-finite ones are left out of that mean.
 */
export function cornerValues(ci, depth, base, elevation, cvals, minDepth, out) {
  const { width, height } = packedSize(ci.nCorners, 1);
  const len = 4 * width * height;
  const data = out && out.length === len ? out : new Float32Array(len);
  data.fill(0);
  const { offsets, cells } = ci;
  for (let c = 0; c < ci.nCorners; c++) {
    let n = 0, sl = 0, se = 0, sd = 0, nv = 0, sv = 0;
    for (let k = offsets[c]; k < offsets[c + 1]; k++) {
      const i = cells[k], d = depth[i];
      if (!(d >= minDepth)) continue;
      const b = base[i], e = elevation[i];
      if (!Number.isFinite(b)) continue;
      n++; sl += b + d; sd += d; se += Number.isFinite(e) ? e : b;
      if (cvals) { const x = cvals[i]; if (Number.isFinite(x)) { nv++; sv += x; } }
    }
    if (!n) continue;
    data[4 * c] = sl / n; data[4 * c + 1] = se / n; data[4 * c + 2] = sd / n;
    data[4 * c + 3] = nv ? sv / nv : 0;
  }
  return { data, width, height };
}

/** The bed at each corner: the mean elevation of its cells (the base where it is missing), as cornerValues' second component. */
export function bedCornerValues(ci, elevation, base) {
  const { width, height } = packedSize(ci.nCorners, 1);
  const data = new Float32Array(4 * width * height);
  const { offsets, cells } = ci;
  for (let c = 0; c < ci.nCorners; c++) {
    let n = 0, s = 0;
    for (let k = offsets[c]; k < offsets[c + 1]; k++) {
      const i = cells[k];
      const e = Number.isFinite(elevation[i]) ? elevation[i] : base[i];
      if (Number.isFinite(e)) { n++; s += e; }
    }
    if (n) data[4 * c + 1] = s / n;
  }
  return { data, width, height };
}

// ---------------------------------------------------------------------------
// Smooth-surface shading (0.9.2): slopes per vertex, not per triangle
// ---------------------------------------------------------------------------

/**
 * Where the smooth surface's slopes are kept: one slot per corner (slot c)
 * and per cell centre (slot nCorners + i), with their positions in metres
 * east (x) and north (y) of the mesh origin (a local plane: the slopes are
 * for shading only). Also, per ring vertex (mesh vertex order), its corner
 * and its offset from its cell's centre, so that surfaceGradients reads them
 * in order. Built with the corners, once.
 */
export function gradientSlots(mesh, ci) {
  const nC = ci.nCorners, n = nC + mesh.nCells, nv = mesh.nVertices;
  const ky = WGS84_A * DEG, kx = ky * Math.cos(mesh.originLat * DEG);
  const x = new Float64Array(n), y = new Float64Array(n);
  const ringCorner = new Int32Array(nv), rx = new Float32Array(nv), ry = new Float32Array(nv);
  for (let i = 0; i < mesh.nCells; i++) {
    const a = mesh.offsets[i], b = mesh.offsets[i + 1];
    let sx = 0, sy = 0;
    for (let v = a; v < b; v++) {
      const c = ci.vertexCorner[v + i + 1];
      ringCorner[v] = c;
      x[c] = mesh.dlon[v] * kx; y[c] = mesh.dlat[v] * ky;
      sx += mesh.dlon[v]; sy += mesh.dlat[v];
    }
    const cx = (sx / (b - a)) * kx, cy = (sy / (b - a)) * ky;    // (the centre as buildSurface places it)
    x[nC + i] = cx; y[nC + i] = cy;
    for (let v = a; v < b; v++) { rx[v] = x[ringCorner[v]] - cx; ry[v] = y[ringCorner[v]] - cy; }
  }
  return { nCorners: nC, n, x, y, ringCorner, rx, ry, acc: null, touched: null, seen: null };
}

/**
 * The smooth surface's slope at each vertex, for shading without facets: the
 * sum of the normals of the fan triangles that meet there (each weighted by
 * its area) gives the area-weighted mean slope of those triangles. The
 * triangles are those buildSurface draws: each cell's centre at its own
 * height, its ring at the corners' heights (cornerValues: the means over the
 * wet cells). Only wet cells' triangles count (depth ≥ minDepth); depth null
 * (the bed surface): every cell's, with the bed as the level. A slope is
 * linear in the heights, so the shader exaggerates it and forms the
 * depth × exaggeration surface's from these two.
 *   corner   cornerValues (or bedCornerValues) data for the same frame
 * @returns {{data:Float32Array, width:number, height:number}} one slot per
 *   RGBA32F texel (packedSize(n, 1)): (∂level/∂x, ∂level/∂y, ∂bed/∂x, ∂bed/∂y),
 *   x east, y north; zero where no drawn triangle meets.
 */
export function surfaceGradients(mesh, ci, slots, corner, depth, base, elevation, minDepth, out) {
  const nC = slots.nCorners, n = slots.n, { ringCorner, rx, ry } = slots;
  if (!slots.acc || slots.acc.length !== 5 * nC) {     // (corner sums, kept zero between calls: only the corners touched are visited)
    slots.acc = new Float64Array(5 * nC); slots.touched = new Int32Array(nC); slots.seen = new Uint8Array(nC);
  }
  const { acc, touched, seen } = slots;
  const { width, height } = packedSize(n, 1);
  const len = 4 * width * height;
  const data = out && out.length === len ? out : new Float32Array(len);
  data.fill(0);
  // per ring vertex of the cell: heights, and the sums of its two triangles
  let hl = new Float64Array(16), hb = new Float64Array(16), sv = new Float64Array(80);
  let nt = 0;
  const bedOnly = !depth;
  for (let i = 0; i < mesh.nCells; i++) {
    const bi = base[i], ei = elevation[i];
    let hc, bc;
    if (bedOnly) {
      bc = Number.isFinite(ei) ? ei : bi;
      if (!Number.isFinite(bc)) continue;
      hc = bc;
    } else {
      const d = depth[i];
      if (!(d >= minDepth) || !Number.isFinite(bi)) continue;
      hc = bi + d; bc = Number.isFinite(ei) ? ei : bi;
    }
    const a = mesh.offsets[i], m = mesh.offsets[i + 1] - a;
    if (m > hl.length) { hl = new Float64Array(m); hb = new Float64Array(m); sv = new Float64Array(5 * m); }
    for (let k = 0; k < m; k++) {
      const c = ringCorner[a + k];
      hb[k] = corner[4 * c + 1] - bc;
      hl[k] = bedOnly ? hb[k] : corner[4 * c] - hc;
    }
    for (let q = 0; q < 5 * m; q++) sv[q] = 0;
    let cl = 0, cm = 0, cb = 0, cn = 0, cz = 0;          // the centre's sums
    for (let k = 0; k < m; k++) {
      const k2 = k + 1 < m ? k + 1 : 0;
      const e1x = rx[a + k], e1y = ry[a + k], e2x = rx[a + k2], e2y = ry[a + k2];
      let nz = e1x * e2y - e1y * e2x;
      if (nz === 0) continue;
      const s = nz < 0 ? -1 : 1;                       // (rings either way round: the normal points up)
      nz *= s;
      const h1 = hl[k], h2 = hl[k2], b1 = hb[k], b2 = hb[k2];
      const lx = s * (e1y * h2 - h1 * e2y), ly = s * (h1 * e2x - e1x * h2);
      const bx = s * (e1y * b2 - b1 * e2y), by = s * (b1 * e2x - e1x * b2);
      cl += lx; cm += ly; cb += bx; cn += by; cz += nz;
      let o = 5 * k;
      sv[o] += lx; sv[o + 1] += ly; sv[o + 2] += bx; sv[o + 3] += by; sv[o + 4] += nz;
      o = 5 * k2;
      sv[o] += lx; sv[o + 1] += ly; sv[o + 2] += bx; sv[o + 3] += by; sv[o + 4] += nz;
    }
    if (cz > 0) {                                      // the centre: its own cell's triangles only
      const o = 4 * (nC + i);
      data[o] = -cl / cz; data[o + 1] = -cm / cz; data[o + 2] = -cb / cz; data[o + 3] = -cn / cz;
    }
    for (let k = 0; k < m; k++) {                      // the corners: one write each
      const c = ringCorner[a + k], o = 5 * c, q = 5 * k;
      if (!seen[c]) { seen[c] = 1; touched[nt++] = c; }
      acc[o] += sv[q]; acc[o + 1] += sv[q + 1]; acc[o + 2] += sv[q + 2]; acc[o + 3] += sv[q + 3]; acc[o + 4] += sv[q + 4];
    }
  }
  for (let j = 0; j < nt; j++) {
    const c = touched[j], o = 5 * c, z = acc[o + 4];
    if (z > 0) {
      data[4 * c] = -acc[o] / z; data[4 * c + 1] = -acc[o + 1] / z;
      data[4 * c + 2] = -acc[o + 2] / z; data[4 * c + 3] = -acc[o + 3] / z;
    }
    acc[o] = acc[o + 1] = acc[o + 2] = acc[o + 3] = acc[o + 4] = 0;
    seen[c] = 0;
  }
  return { data, width, height };
}
