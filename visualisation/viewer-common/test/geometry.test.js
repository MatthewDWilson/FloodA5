import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { parseMesh } from '../js/protocol.js';
import {
  geodeticToEcef, ecefToGeodetic, buildSurface, buildGlb, packCellValues, packStatic, packedSize,
  CellLocator, cellCentres, meshBounds, finiteMedian, finiteQuantile, WGS84_A, WGS84_B,
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

test('buildSurface triangulates each cell on its own vertices', () => {
  const m = fixtureMesh();
  const s = buildSurface(m, 20);
  assert.equal(s.positions.length, 3 * 12);
  assert.equal(s.indices.length, 3 * (2 + 1 + 3));
  assert.deepEqual(Array.from(s.cell), [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2]);
  for (const i of s.indices) assert.ok(i < 12);
  // relative to the centre: small numbers, and the square's side is 0.001° of longitude
  for (const x of s.positions) assert.ok(Math.abs(x) < 1000);
  const side = Math.hypot(s.positions[3] - s.positions[0], s.positions[4] - s.positions[1], s.positions[5] - s.positions[2]);
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
  assert.deepEqual(json.meshes[0].primitives[0].attributes, { POSITION: 0, _CELL: 1 });
  assert.equal(json.accessors[0].count, 12);
  assert.equal(json.accessors[2].count, 18);
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
  const bin = 28 + jlen, iv = json.bufferViews[2];
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
