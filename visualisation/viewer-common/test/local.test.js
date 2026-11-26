// Local files (js/local/): the /api answered in the browser must agree with
// the server. fixtures/local_run.json (test/vis_local_fixture.jl) holds two
// small output files' contents and the server's answers to a list of
// requests; LocalRun answers the same requests from the same contents.
// FLOODA5_LOCAL_FIXTURE=path checks another fixture (test/test_vis_local.jl
// passes a freshly generated one).

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { LocalRun, pickFrames, sortedQuantile, parseJsonInf, f32short, encodeMesh, depthDatum } from '../js/local/localrun.js';
import { hazardParams, deriveHazard, runStatistics, chenCriticalVelocity, HAZARD_DEFAULTS } from '../js/local/hazardcalc.js';
import { parseMesh } from '../js/protocol.js';

const path = process.env.FLOODA5_LOCAL_FIXTURE || new URL('./fixtures/local_run.json', import.meta.url);
const FIX = JSON.parse(readFileSync(path, 'utf8'));

// ---------------------------------------------------------------------------
// The fixture as a reader (localrun.js header)
// ---------------------------------------------------------------------------

const num = x => (typeof x === 'string' ? (x === 'NaN' ? NaN : x === 'Inf' ? Infinity : x === '-Inf' ? -Infinity : Number(x)) : x);
const attrVal = v => (v && typeof v === 'object' && !Array.isArray(v) && 'float' in v ? num(v.float)
  : Array.isArray(v) ? v.map(attrVal) : v);

function fixtureReader(h5) {
  const obj = p => h5[p];
  const isNumeric = o => /^(Float|Int|UInt)/.test(o.dtype);
  const data = o => {
    if (!o.shape.length) return isNumeric(o) ? num(o.data) : o.data;
    return isNumeric(o) ? Float64Array.from(o.data, num) : o.data.slice();
  };
  return {
    has: p => p in h5,
    isDataset: p => obj(p)?.kind === 'dataset',
    keys: p => (obj(p)?.kind === 'group' ? obj(p).keys.slice() : []),
    attrs: p => Object.fromEntries(Object.entries(obj(p)?.attrs || {}).map(([k, v]) => [k, attrVal(v)])),
    shape: p => obj(p).shape.slice(),
    read: p => data(obj(p)),
    row: (p, k) => { const o = obj(p), n = o.shape[1]; return Float64Array.from(o.data.slice(k * n, (k + 1) * n), num); },
    col: (p, c) => { const o = obj(p), [nf, n] = o.shape; return Float64Array.from({ length: nf }, (_, k) => num(o.data[k * n + c])); },
  };
}

// ---------------------------------------------------------------------------
// Comparison
// ---------------------------------------------------------------------------

const SKIP = new Set(['file', 'source', 'local']);       // where the file is (server: a path; here: its name)
const close = (a, b) => a === b || (Number.isNaN(a) && Number.isNaN(b)) ||
  Math.abs(a - b) <= 1e-6 * Math.max(Math.abs(a), Math.abs(b));

function sameJson(a, b, at) {
  if (typeof a === 'number' && typeof b === 'number') return assert.ok(close(a, b), `${at}: ${a} vs ${b}`);
  if (Array.isArray(a)) {
    assert.ok(Array.isArray(b), `${at}: array expected`);
    assert.equal(a.length, b.length, `${at}: length`);
    return a.forEach((x, i) => sameJson(x, b[i], `${at}[${i}]`));
  }
  if (a && typeof a === 'object') {
    assert.ok(b && typeof b === 'object', `${at}: object expected`);
    const ka = Object.keys(a).filter(k => !SKIP.has(k)).sort(), kb = Object.keys(b).filter(k => !SKIP.has(k)).sort();
    assert.deepEqual(kb, ka, `${at}: keys`);
    return ka.forEach(k => sameJson(a[k], b[k], `${at}.${k}`));
  }
  assert.equal(b, a, at);
}

for (const run of FIX.runs) {
  test(`local file ${run.name}: every answer as the server's`, () => {
    const lr = new LocalRun(fixtureReader(run.h5), { name: run.name, viewer: 'cesium' });
    let exact = 0;
    for (const q of run.requests) {
      const r = lr.handle(q.path);
      assert.equal(r.status, q.status, `${q.path}: status`);
      if (q.status !== 200) continue;
      if (q.json !== undefined) {
        assert.ok(r.json !== undefined, `${q.path}: JSON expected`);
        sameJson(q.json, r.json, q.path);
      } else {
        const want = Buffer.from(q.base64, 'base64');
        const got = Buffer.from(r.bytes);
        if (q.path === '/api/mesh') { assert.ok(got.equals(want), 'mesh bytes'); continue; }
        assert.equal(got.length, want.length, `${q.path}: length`);
        const a = new Float32Array(want.buffer, want.byteOffset, want.length / 4);
        const b = new Float32Array(got.buffer, got.byteOffset, got.length / 4);
        // equal to the last bit; the hazard degree within 1e-6 (Math.pow vs Julia's ^)
        const tol = /hazard_degree/.test(q.path);
        a.forEach((x, i) => assert.ok(tol ? close(x, b[i]) : Object.is(x, b[i]) || (Number.isNaN(x) && Number.isNaN(b[i])),
                                      `${q.path}[${i}]: ${x} vs ${b[i]}`));
        if (got.equals(want)) exact++;
      }
    }
    assert.ok(exact > 0);
  });
}

test('local file: opening errors', () => {
  const h5 = FIX.runs[1].h5;
  const without = drop => fixtureReader(Object.fromEntries(Object.entries(h5).filter(([k]) => !drop(k))));
  assert.throws(() => new LocalRun(without(k => k.startsWith('/results'))), /not a FloodA5 output/);
  assert.throws(() => new LocalRun(without(k => k === '/mesh/polygon_offsets')), /no cell polygons/);
  const v1 = { '/': { kind: 'group', keys: ['frames'], attrs: {} }, '/frames': { kind: 'group', keys: [], attrs: {} } };
  assert.throws(() => new LocalRun(fixtureReader(v1)), /format-1/);
});

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

test('pickFrames is Julia\'s unique(round.(Int, range(1, n; length = m)))', () => {
  for (const { n, m, picks } of FIX.picks) assert.deepEqual(pickFrames(n, m), picks, `n=${n} m=${m}`);
});

test('sortedQuantile interpolates between order statistics', () => {
  assert.ok(Number.isNaN(sortedQuantile([], 0.5)));
  assert.equal(sortedQuantile([4], 0.98), 4);
  assert.equal(sortedQuantile([0, 10], 0.25), 2.5);
  assert.equal(sortedQuantile([1, 2, 3, 4, 5], 0.5), 3);
});

test('parseJsonInf reads JSON3 allow_inf output', () => {
  const o = parseJsonInf('{"a":Infinity,"b":-Infinity,"c":NaN,"s":"x Infinity \\"NaN\\"","v":[1,NaN]}');
  assert.deepEqual(o, { a: 'Inf', b: '-Inf', c: 'NaN', s: 'x Infinity "NaN"', v: [1, 'NaN'] });
});

test('f32short is the shortest decimal of the Float32', () => {
  assert.equal(f32short(Math.fround(0.1)), 0.1);
  assert.equal(f32short(Math.fround(15.379024)), 15.379024);
  assert.equal(Math.fround(f32short(Math.fround(1 / 3))), Math.fround(1 / 3));
  assert.equal(f32short(0), 0);
  assert.ok(Number.isNaN(f32short(NaN)));
});

test('encodeMesh round-trips through parseMesh, dropping closing vertices', () => {
  const offs = [0, 5, 8], lon = [1, 2, 2, 1, 1, 3, 4, 3.5], lat = [1, 1, 2, 2, 1, 1, 1, 2];
  const m = parseMesh(encodeMesh(offs, lon, lat, [5, NaN], depthDatum('sgs', [5, NaN], { sgs_z_min: [4, NaN] })));
  assert.equal(m.nCells, 2);
  assert.equal(m.nVertices, 7);                      // 4 + 3: the closing repeat dropped
  assert.deepEqual(Array.from(m.offsets), [0, 4, 7]);
  assert.equal(m.originLon, 2.5); assert.equal(m.originLat, 1.5);
  assert.equal(m.base[0], 4); assert.ok(Number.isNaN(m.base[1]));
  assert.throws(() => encodeMesh([0, 2], [0, 1], [0, 1], [0], [0]), /fewer than 3/);
});

test('hazardParams: defaults, clamping, invalid values', () => {
  assert.deepEqual(hazardParams({}), { ...HAZARD_DEFAULTS });
  const p = hazardParams({ landuse: 'urban', mass: '500', height: '0.1', theta: 'abc', alpha: '', beta: 'Inf' });
  assert.equal(p.landuse, 'urban'); assert.equal(p.mass, 300); assert.equal(p.height, 0.3);
  assert.equal(p.theta, 75); assert.equal(p.alpha, 3.064); assert.equal(p.beta, 0.303);
  assert.equal(hazardParams({ landuse: 'mars' }).landuse, 'fd2320');
});

test('hazard variables: dry, NaN, floating, still water', () => {
  const p = { ...HAZARD_DEFAULTS };
  const hr = deriveHazard('hazard_rating', [0, 0.2, 1.0, NaN], [1, 1, 3, 1], p);
  assert.deepEqual(Array.from(hr).slice(0, 3), [0, Math.fround(0.2 * 1.5 + 0.5), Math.fround(1.0 * 3.5 + 1)]);
  assert.ok(Number.isNaN(hr[3]));
  const hd = deriveHazard('hazard_degree', [0, 2.0, 0.5, 0.5], [1, 0, 0, 100], p);
  assert.deepEqual(Array.from(hd), [0, 1, 0, 1]);    // dry; floating; still; fast
  assert.equal(chenCriticalVelocity(0, p), Infinity);
  assert.deepEqual(Array.from(deriveHazard('hazard_class_arr', [0.2, 0.4, 5], [1, 1, 5], p)), [1, 2, 6]);
});

test('runStatistics: duration weights and the first flooded time', () => {
  const t = [0, 10, 30];
  const d = [[0, 0.2], [0.2, 0.2], [0, 0.05]];
  const st = runStatistics(t, k => [d[k], null], 2, { ...HAZARD_DEFAULTS, flood_depth: 0.1 });
  assert.deepEqual(Array.from(st.flood_duration), [15, 20]);       // frame 2 stands for (30 − 0)/2; cell 2: 10/2 + 15
  assert.deepEqual(Array.from(st.time_first_flooded), [10, 0]);
  assert.deepEqual(Array.from(st.time_max_water_depth), [10, 0]);     // the first frame with the maximum
  assert.equal(st.max_velocity, undefined);
});
