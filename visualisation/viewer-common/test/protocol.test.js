// node --test visualisation/cesium/test/
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { parseMesh, decodeFloat32, FrameList, FrameCache, formatTime, PROTOCOL } from '../js/protocol.js';

const fixture = () => {
  // written by test/test_vis_server.jl (FLOODA5_UPDATE_FIXTURES=1): three cells
  const b = readFileSync(fileURLToPath(new URL('./fixtures/mesh_small.bin', import.meta.url)));
  return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
};

test('protocol version', () => assert.equal(PROTOCOL, 2));

test('parseMesh decodes the Julia encoder output', () => {
  const m = parseMesh(fixture());
  assert.equal(m.nCells, 3);
  assert.equal(m.nVertices, 12);
  assert.deepEqual(Array.from(m.offsets), [0, 4, 7, 12]);
  assert.ok(Math.abs(m.originLon - (172.0 + 172.0035) / 2) < 1e-12);
  assert.ok(Math.abs(m.originLat - (-43.0 - 42.999) / 2) < 1e-12);
  assert.ok(Math.abs(m.originLon + m.dlon[0] - 172.0) < 1e-7);
  assert.ok(Math.abs(m.originLat + m.dlat[6] - -42.999) < 1e-7);   // triangle apex
  assert.equal(m.elevation[0], 10);
  assert.ok(Number.isNaN(m.elevation[2]));
  assert.deepEqual(Array.from(m.base), [9.5, 11.5, 12]);
});

test('parseMesh rejects bad payloads', () => {
  const buf = fixture();
  const bad = new Uint8Array(buf.slice(0)); bad[0] = 0x58;
  assert.throws(() => parseMesh(bad.buffer), /not a FloodA5 mesh/);
  assert.throws(() => parseMesh(buf.slice(0, buf.byteLength - 4)), /expected/);
  const v2 = new Uint8Array(buf.slice(0)); new DataView(v2.buffer).setUint32(4, 2, true);
  assert.throws(() => parseMesh(v2.buffer), /unsupported mesh version/);
  assert.throws(() => parseMesh(new ArrayBuffer(8)), /too short/);
  assert.throws(() => parseMesh(new Uint8Array(4)), /ArrayBuffer/);
});

test('decodeFloat32', () => {
  const a = new Float32Array([1, 2.5, -3]);
  assert.deepEqual(Array.from(decodeFloat32(a.buffer, 3)), [1, 2.5, -3]);
  assert.throws(() => decodeFloat32(a.buffer, 4), /3 values/);
  assert.throws(() => decodeFloat32(new ArrayBuffer(5)), /Float32/);
});

test('FrameList lookups', () => {
  const f = new FrameList([{ id: 2, t: 60 }, { id: 1, t: 0 }, { id: 3, t: 120 }]);
  assert.deepEqual(f.frames.map(x => x.id), [1, 2, 3]);
  assert.equal(f.indexOfId(3), 2);
  assert.equal(f.indexOfId(9), -1);
  assert.equal(f.indexAtTime(-5), 0);
  assert.equal(f.indexAtTime(0), 0);
  assert.equal(f.indexAtTime(59.9), 0);
  assert.equal(f.indexAtTime(60), 1);
  assert.equal(f.indexAtTime(1e9), 2);
  assert.equal(f.add(3, 120), false);
  assert.equal(f.add(4, 180), true);
  assert.equal(f.last().id, 4);
  assert.equal(f.add(5, 90), true);              // out of order: re-sorted
  assert.deepEqual(f.frames.map(x => x.id), [1, 2, 5, 3, 4]);
  assert.equal(f.indexOfId(4), 4);
  assert.equal(new FrameList().indexAtTime(0), -1);
});

test('FrameCache is LRU', () => {
  const c = new FrameCache(2);
  c.set('a', 1); c.set('b', 2);
  assert.equal(c.get('a'), 1);                    // a now most recent
  c.set('c', 3);
  assert.equal(c.has('b'), false);
  assert.equal(c.has('a'), true);
  assert.equal(c.has('c'), true);
});

test('formatTime', () => {
  assert.equal(formatTime(5.25), '5.3 s');
  assert.equal(formatTime(42), '42 s');
  assert.equal(formatTime(125), '2:05');
  assert.equal(formatTime(3 * 3600 + 7 * 60 + 9), '3:07:09');
  assert.equal(formatTime(NaN), '–');
});
