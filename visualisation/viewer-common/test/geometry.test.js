import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { parseMesh } from '../js/protocol.js';
import {
  geodeticToEcef, ecefToGeodetic, buildSurface, buildGlb, packCellValues, packStatic, packedSize,
  CellLocator, cellCentres, meshBounds, finiteMedian, finiteQuantile, WGS84_A, WGS84_B,
  cornerIndex, cornerValues, bedCornerValues, gradientSlots, surfaceGradients,
} from '../js/geometry.js';

const fixtureMesh = () => {
  const b = readFileSync(fileURLToPath(new URL('./fixtures/mesh_small.bin', import.meta.url)));
  return parseMesh(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength));
};

// A parsed-mesh-like grid of nx × ny square-ish cells (0.001° spacing)
function gridMesh(nx, ny, lon0 = -2.95, lat0 = 54.88, d = 0.001) {
  const nCells = nx * ny, nVertices = 4 * nCells;
  const originLon = lon0 + nx * d / 2, originLat = lat0 + ny * d / 2;
  const offsets = new Uint32Array(nCells + 1), dlon = new Float32Array(nVertices), dlat = new Float32Array(nVertices);
  let v = 0;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    const c = j * nx + i;
    offsets[c] = v;
    for (const [a, b] of [[0, 0], [1, 0], [1, 1], [0, 1]]) {
      dlon[v] = lon0 + (i + a) * d - originLon; dlat[v] = lat0 + (j + b) * d - originLat; v++;
    }
  }
  offsets[nCells] = v;
  return { nCells, nVertices, originLon, originLat, offsets, dlon, dlat,
           elevation: new Float32Array(nCells).fill(10), base: new Float32Array(nCells).fill(9) };
}

test('geodeticToEcef', () => {
  const close = (a, b, tol = 1e-6) => a.forEach((x, k) => assert.ok(Math.abs(x - b[k]) < tol, `${x} vs ${b[k]}`));
  close(geodeticToEcef(0, 0, 0), [WGS84_A, 0, 0]);
  close(geodeticToEcef(90, 0, 100), [0, WGS84_A + 100, 0], 1e-6);
  close(geodeticToEcef(0, 90, 0), [0, 0, WGS84_B], 1e-6);
  // height is along the ellipsoid normal
  const p0 = geodeticToEcef(-2.9, 54.9, 0), p1 = geodeticToEcef(-2.9, 54.9, 50);
  assert.ok(Math.abs(Math.hypot(...p1.map((x, k) => x - p0[k])) - 50) < 1e-6);
});

test('ecefToGeodetic inverts geodeticToEcef', () => {
  for (const [lon, lat, h] of [[-2.93, 54.893, 3000], [-2.93, 54.893, -60], [172.6, -43.5, 15], [0, 0, 0],
                               [179.9, 89.9, 100], [-120, -89.99, 2000], [45, 30, 8848], [10, 60, -400]]) {
    const g = ecefToGeodetic(...geodeticToEcef(lon, lat, h));
    assert.ok(Math.abs(g.lon - lon) < 1e-9, `lon ${g.lon} vs ${lon}`);
    assert.ok(Math.abs(g.lat - lat) < 1e-9, `lat ${g.lat} vs ${lat}`);
    assert.ok(Math.abs(g.h - h) < 1e-3, `h ${g.h} vs ${h}`);
  }
  const pole = ecefToGeodetic(0, 0, WGS84_B + 50);
  assert.ok(Math.abs(pole.lat - 90) < 1e-9 && Math.abs(pole.h - 50) < 1e-3);
});

test('buildSurface: a fan from each cell centre, on the cell\'s own vertices', () => {
  const m = fixtureMesh();
  const s = buildSurface(m, 20);
  assert.equal(s.positions.length, 3 * (12 + 3));          // ring vertices + one centre per cell
  assert.equal(s.indices.length, 3 * 12);                  // one triangle per side
  assert.deepEqual(Array.from(s.cell), [0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]);
  assert.deepEqual(Array.from(s.edge), [1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]);
  for (const i of s.indices) assert.ok(i < 15);
  // every triangle: the cell's centre, then two consecutive ring vertices (closing the ring)
  const tri = k => Array.from(s.indices.subarray(3 * k, 3 * k + 3));
  assert.deepEqual([tri(0), tri(3), tri(4), tri(6)], [[0, 1, 2], [0, 4, 1], [5, 6, 7], [5, 8, 6]]);
  for (let k = 0; k < 12; k++) {
    const [c, a, b] = tri(k);
    assert.equal(s.edge[c], 1); assert.equal(s.edge[a], 0); assert.equal(s.edge[b], 0);
    assert.ok(s.cell[c] === s.cell[a] && s.cell[a] === s.cell[b]);
  }
  // the centre is the mean of the ring
  for (let d = 0; d < 3; d++) {
    const mean = [1, 2, 3, 4].reduce((acc, v) => acc + s.positions[3 * v + d], 0) / 4;
    assert.ok(Math.abs(s.positions[d] - mean) < 0.05);
  }
  // relative to the centre: small numbers, and the square's side is 0.001° of longitude
  for (const x of s.positions) assert.ok(Math.abs(x) < 1000);
  const side = Math.hypot(s.positions[6] - s.positions[3], s.positions[7] - s.positions[4], s.positions[8] - s.positions[5]);
  const expect = (Math.PI / 180) * 0.001 * WGS84_A * Math.cos(43 * Math.PI / 180) / Math.sqrt(1 - 0.00669438 * Math.sin(43 * Math.PI / 180) ** 2);
  assert.ok(Math.abs(side - expect) < 0.01, `${side} vs ${expect}`);
  for (let k = 0; k < 3; k++) assert.ok(s.min[k] <= s.max[k]);
});

test('buildGlb writes a well-formed GLB', () => {
  const s = buildSurface(fixtureMesh(), 0);
  const glb = buildGlb(s);
  const dv = new DataView(glb.buffer, glb.byteOffset, glb.byteLength);
  assert.equal(dv.getUint32(0, true), 0x46546c67);
  assert.equal(dv.getUint32(4, true), 2);
  assert.equal(dv.getUint32(8, true), glb.byteLength);
  const jlen = dv.getUint32(12, true);
  assert.equal(jlen % 4, 0);
  assert.equal(dv.getUint32(16, true), 0x4e4f534a);
  const json = JSON.parse(new TextDecoder().decode(glb.subarray(20, 20 + jlen)));
  assert.equal(json.asset.version, '2.0');
  assert.deepEqual(json.meshes[0].primitives[0].attributes, { POSITION: 0, _CELL: 1, _FA5EDGE: 2 });
  assert.equal(json.accessors[0].count, 15);
  assert.equal(json.accessors[2].count, 15);
  assert.equal(json.accessors[3].count, 36);
  assert.equal(json.materials[0].doubleSided, true);
  const binLen = dv.getUint32(20 + jlen, true);
  assert.equal(dv.getUint32(24 + jlen, true), 0x004e4942);
  assert.equal(binLen, json.buffers[0].byteLength);
  assert.equal(28 + jlen + binLen, glb.byteLength);
  for (const bv of json.bufferViews) {
    assert.equal(bv.byteOffset % 4, 0);
    assert.ok(bv.byteOffset + bv.byteLength <= binLen);
  }
  // the indices round-trip
  const bin = 28 + jlen, iv = json.bufferViews[3];
  const idx = new Uint32Array(glb.buffer.slice(glb.byteOffset + bin + iv.byteOffset, glb.byteOffset + bin + iv.byteOffset + iv.byteLength));
  assert.deepEqual(Array.from(idx), Array.from(s.indices));
});

test('packCellValues: four cells per texel, non-finite → 0', () => {
  assert.deepEqual(packedSize(1), { width: 1, height: 1 });
  assert.deepEqual(packedSize(4 * 2048 + 1), { width: 2048, height: 2 });
  const v = new Float32Array([1, 2, NaN, 4, 5]);
  const p = packCellValues(v, 5);
  assert.equal(p.width, 2); assert.equal(p.height, 1);
  assert.deepEqual(Array.from(p.data), [1, 2, 0, 4, 5, 0, 0, 0]);
  const p2 = packCellValues(new Float32Array([7, 7, 7, 7, 7]), 5, p.data);
  assert.equal(p2.data, p.data);                     // buffer reused
  assert.deepEqual(Array.from(p2.data), [7, 7, 7, 7, 7, 0, 0, 0]);
});

test('packStatic fills missing values', () => {
  const p = packStatic(new Float32Array([10, NaN, NaN]), new Float32Array([9, 11, NaN]), 3);
  assert.equal(p.width, 3);
  assert.deepEqual(Array.from(p.data), [10, 9, 0, 0, 11, 11, 0, 0, 10, 10, 0, 0]);
});

test('CellLocator finds the containing cell', () => {
  const m = gridMesh(40, 7);
  const loc = new CellLocator(m);
  const { lon, lat } = cellCentres(m);
  for (let c = 0; c < m.nCells; c++) assert.equal(loc.find(lon[c], lat[c]), c);
  // points near corners
  let seed = 1;
  const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  for (let k = 0; k < 500; k++) {
    const x = -2.95 + 40 * 0.001 * rnd(), y = 54.88 + 7 * 0.001 * rnd();
    const c = loc.find(x, y);
    const i = Math.floor((x + 2.95) / 0.001), j = Math.floor((y - 54.88) / 0.001);
    assert.equal(c, j * 40 + i, `point ${x},${y}`);
  }
  assert.equal(loc.find(-3.5, 54.9), -1);
  // the fixture's triangle and pentagon
  const f = fixtureMesh(), lf = new CellLocator(f);
  assert.equal(lf.find(172.0015, -42.9997), 1);
  assert.equal(lf.find(172.0025, -42.9996), 2);
  // a single row of cells
  const row = gridMesh(30, 1), lr = new CellLocator(row), cr = cellCentres(row);
  for (let c = 0; c < row.nCells; c++) assert.equal(lr.find(cr.lon[c], cr.lat[c]), c);
});

test('meshBounds and quantiles', () => {
  const m = gridMesh(3, 2, 10, 20);
  const b = meshBounds(m);
  [10, 20, 10.003, 20.002].forEach((x, k) => assert.ok(Math.abs(b[k] - x) < 1e-6));
  assert.equal(finiteMedian([3, NaN, 1, 2]), 2);
  assert.equal(finiteMedian([4, 1, 2, 3]), 2.5);
  assert.ok(Number.isNaN(finiteMedian([NaN])));
  assert.equal(finiteQuantile([0, 10, 20, 30, 40], 0.25), 10);
  assert.equal(finiteQuantile([0, 10], 0.5), 5);
});

// ---------------------------------------------------------------------------
// Smooth surface (0.8.14): shared corners and corner values
// ---------------------------------------------------------------------------

// 2 × 2 squares of 0.001°, rings counter-clockwise: 9 distinct corners, the
// middle one shared by all four cells
function quadMesh() {
  const sq = (c, r) => [[c, r], [c + 1, r], [c + 1, r + 1], [c, r + 1]];
  const rings = [sq(0, 0), sq(1, 0), sq(0, 1), sq(1, 1)];
  const offsets = new Uint32Array(5), dlon = [], dlat = [];
  rings.forEach((ring, i) => {
    for (const [x, y] of ring) { dlon.push(0.001 * x - 0.001); dlat.push(0.001 * y - 0.001); }
    offsets[i + 1] = offsets[i] + ring.length;
  });
  return { nCells: 4, nVertices: 16, originLon: 172.6, originLat: -43.5, offsets,
           dlon: Float32Array.from(dlon), dlat: Float32Array.from(dlat),
           elevation: Float32Array.from([10, 11, 12, 13]), base: Float32Array.from([9, 10, 11, 12]) };
}

test('cornerIndex: shared corners, centres −1, corner → cells', () => {
  const m = quadMesh(), ci = cornerIndex(m);
  assert.equal(ci.nCorners, 9);
  assert.equal(ci.vertexCorner.length, m.nVertices + m.nCells);
  const surf = buildSurface(m);
  for (let v = 0; v < ci.vertexCorner.length; v++) assert.equal(ci.vertexCorner[v] === -1, surf.edge[v] === 1);
  // the middle corner: cell 0's third ring vertex (surface vertex 0 + 0 + 1 + 2)
  const mid = ci.vertexCorner[3];
  assert.deepEqual(Array.from(ci.cells.subarray(ci.offsets[mid], ci.offsets[mid + 1])).sort(), [0, 1, 2, 3]);
  // every corner's cells hold it; 16 ring vertices in all
  assert.equal(ci.offsets[ci.nCorners], 16);
  // vertices within the quantum weld; a different place does not
  const m2 = quadMesh(); m2.dlon = Float32Array.from(m2.dlon); m2.dlon[6] += 2e-8;   // (cell 1, its third vertex)
  assert.equal(cornerIndex(m2).nCorners, 9);
});

test('cornerValues: means over the wet cells; dry corners zero', () => {
  const m = quadMesh(), ci = cornerIndex(m);
  const depth = Float32Array.from([0.5, 0, 0.2, 1.0]);
  const cv = Float32Array.from([1, 2, NaN, 4]);
  const p = cornerValues(ci, depth, m.base, m.elevation, cv, 0.1);
  const at = c => Array.from(p.data.subarray(4 * c, 4 * c + 4));
  const mid = ci.vertexCorner[3];
  const [lv, bd, dp, cval] = at(mid);           // wet: cells 0, 2, 3
  assert.ok(Math.abs(lv - (9.5 + 11.2 + 13) / 3) < 1e-5);
  assert.ok(Math.abs(bd - (10 + 12 + 13) / 3) < 1e-5);
  assert.ok(Math.abs(dp - 1.7 / 3) < 1e-6);
  assert.ok(Math.abs(cval - (1 + 4) / 2) < 1e-6);          // (NaN left out)
  // cell 1's outer corner: only cell 1, dry → zeros
  const c1 = ci.vertexCorner[1 + 5 + 1];                    // cell 1's second ring vertex
  assert.deepEqual(at(c1), [0, 0, 0, 0]);
  const b = bedCornerValues(ci, Float32Array.from([10, NaN, 12, 13]), m.base);
  assert.ok(Math.abs(b.data[4 * mid + 1] - (10 + 10 + 12 + 13) / 4) < 1e-5);   // (missing elevation: the base)
});

test('buildGlb: the corner attribute when asked for', () => {
  const m = quadMesh(), surf = buildSurface(m), ci = cornerIndex(m);
  const json = glb => JSON.parse(new TextDecoder().decode(glb.subarray(20, 20 + new DataView(glb.buffer).getUint32(12, true))));
  assert.equal(json(buildGlb(surf)).meshes[0].primitives[0].attributes._FA5CORNER, undefined);
  const j = json(buildGlb(surf, ci.vertexCorner));
  assert.equal(j.meshes[0].primitives[0].attributes._FA5CORNER, 4);
  assert.equal(j.accessors[4].count, surf.cell.length);
});

// ---------------------------------------------------------------------------
// Smooth-surface shading (0.9.2): slopes at the vertices
// ---------------------------------------------------------------------------

/**
 * Heights from planes (level L, bed B, of x east and y north in metres) at
 * every corner and cell centre: corner data as cornerValues lays it out, and
 * per-cell base, elevation and depth giving the same centre heights.
 */
function planar(m, ci, slots, L, B) {
  const corner = new Float32Array(4 * ci.nCorners);
  for (let c = 0; c < ci.nCorners; c++) {
    corner[4 * c] = L(slots.x[c], slots.y[c]); corner[4 * c + 1] = B(slots.x[c], slots.y[c]);
  }
  const n = m.nCells, base = new Float32Array(n), elevation = new Float32Array(n), depth = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const x = slots.x[ci.nCorners + i], y = slots.y[ci.nCorners + i];
    elevation[i] = base[i] = B(x, y); depth[i] = L(x, y) - base[i];
  }
  return { corner, base, elevation, depth };
}
const L0 = (x, y) => 20 + 0.01 * x - 0.02 * y, B0 = (x, y) => 3 + 0.03 * x + 0.005 * y;
const slopeAt = (g, k) => Array.from(g.data.subarray(4 * k, 4 * k + 4));
const near = (a, b, eps) => a.every((v, k) => Math.abs(v - b[k]) <= eps);

test('gradientSlots: corners then cell centres, metres east and north of the origin', () => {
  const m = quadMesh(), ci = cornerIndex(m), s = gradientSlots(m, ci);
  assert.equal(s.n, ci.nCorners + m.nCells);
  const mid = ci.vertexCorner[3];                           // at the origin (dlon = dlat = 0)
  assert.ok(Math.abs(s.x[mid]) < 1e-6 && Math.abs(s.y[mid]) < 1e-6);
  // cell 0's centre: 0.0005° west and south of the origin
  const ky = WGS84_A * Math.PI / 180, kx = ky * Math.cos(-43.5 * Math.PI / 180);
  assert.ok(Math.abs(s.x[ci.nCorners] + 0.0005 * kx) < 1e-3);
  assert.ok(Math.abs(s.y[ci.nCorners] + 0.0005 * ky) < 1e-3);
});

test('surfaceGradients: a plane gives its slope at every vertex (squares and A5 pentagons)', () => {
  for (const m of [quadMesh(), fixtureMesh()]) {
    const ci = cornerIndex(m), slots = gradientSlots(m, ci);
    const { corner, base, elevation, depth } = planar(m, ci, slots, L0, B0);
    const g = surfaceGradients(m, ci, slots, corner, depth, base, elevation, 0.05);
    assert.deepEqual([g.width, g.height], [packedSize(slots.n, 1).width, packedSize(slots.n, 1).height]);
    for (let k = 0; k < slots.n; k++) assert.ok(near(slopeAt(g, k), [0.01, -0.02, 0.03, 0.005], 1e-5), `slot ${k}: ${slopeAt(g, k)}`);
    // the bed surface (depth null): its own slope in both pairs
    const b = surfaceGradients(m, ci, slots, corner, null, base, elevation, 0);
    for (let k = 0; k < slots.n; k++) assert.ok(near(slopeAt(b, k), [0.03, 0.005, 0.03, 0.005], 1e-5));
  }
});

test('surfaceGradients: rings either way round give the same slopes', () => {
  const m = quadMesh(), ci = cornerIndex(m), slots = gradientSlots(m, ci);
  const bump = (x, y) => 10 + Math.exp(-(x * x + y * y) / 3000);
  const { corner, base, elevation, depth } = planar(m, ci, slots, bump, B0);
  const g = surfaceGradients(m, ci, slots, corner, depth, base, elevation, 0.05);
  // reverse every ring (clockwise): same corners, same slopes
  const r = quadMesh();
  for (let i = 0; i < r.nCells; i++) {
    const a = r.offsets[i], b = r.offsets[i + 1];
    r.dlon.subarray(a, b).reverse(); r.dlat.subarray(a, b).reverse();
  }
  const cr = cornerIndex(r), sr = gradientSlots(r, cr);
  assert.equal(cr.nCorners, ci.nCorners);
  const pr = planar(r, cr, sr, bump, B0);
  const gr = surfaceGradients(r, cr, sr, pr.corner, pr.depth, pr.base, pr.elevation, 0.05);
  // match the slots by position (corners are numbered in the order the cells use them)
  const at = (sl, k) => `${Math.round(sl.x[k] * 1000)},${Math.round(sl.y[k] * 1000)}`;
  const byPos = new Map(Array.from({ length: slots.n }, (_, k) => [at(slots, k), k]));
  for (let k = 0; k < sr.n; k++) {
    const j = byPos.get(at(sr, k));
    assert.ok(j !== undefined && near(slopeAt(gr, k), slopeAt(g, j), 1e-7), `slot ${k}`);
  }
});

test('surfaceGradients: a peak at a corner is flat there — the vertex mean, not one triangle\'s slope', () => {
  const m = quadMesh(), ci = cornerIndex(m), slots = gradientSlots(m, ci);
  const peak = (x, y) => 12 - 0.01 * Math.hypot(x, y);           // a cone on the middle corner
  const { corner, base, elevation, depth } = planar(m, ci, slots, peak, B0);
  const g = surfaceGradients(m, ci, slots, corner, depth, base, elevation, 0.05);
  const mid = ci.vertexCorner[3];
  const [gx, gy] = slopeAt(g, mid);
  assert.ok(Math.hypot(gx, gy) < 1e-6, `${gx}, ${gy}`);          // (each triangle there slopes ~0.01)
  // cell 0's centre (south-west of the peak) slopes up towards it: north-east
  const [cx, cy] = slopeAt(g, ci.nCorners);
  assert.ok(cx > 0.003 && cy > 0.003, `${cx}, ${cy}`);
});

test('surfaceGradients: dry cells do not count; no triangle, no slope', () => {
  const m = quadMesh(), ci = cornerIndex(m), slots = gradientSlots(m, ci);
  const { corner, base, elevation, depth } = planar(m, ci, slots, L0, B0);
  depth[1] = 0.01;                                              // cell 1 dry (threshold 0.05)
  const tilted = Float32Array.from(corner);
  const c1 = ci.vertexCorner[1 + 5 + 1];                        // cell 1's outer corner (its alone)
  tilted[4 * c1] += 5;                                           // would tilt cell 1's triangles
  const g = surfaceGradients(m, ci, slots, tilted, depth, base, elevation, 0.05);
  assert.deepEqual(slopeAt(g, ci.nCorners + 1), [0, 0, 0, 0]);  // cell 1's centre
  assert.deepEqual(slopeAt(g, c1), [0, 0, 0, 0]);
  for (const k of [ci.nCorners, ci.nCorners + 2, ci.nCorners + 3, ci.vertexCorner[3]]) {
    assert.ok(near(slopeAt(g, k), [0.01, -0.02, 0.03, 0.005], 1e-5), `slot ${k}`);
  }
  // reuses its output buffer
  const again = surfaceGradients(m, ci, slots, corner, depth, base, elevation, 0.05, g.data);
  assert.equal(again.data, g.data);
});
