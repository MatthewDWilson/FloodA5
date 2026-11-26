// Files on a web server read in pieces (js/local/rangefile.js).

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { RangeCache, urlFileName, resolveRunUrl } from '../js/local/rangefile.js';

// a 1000-byte "file" and a fetcher that records its requests
function source(size = 1000) {
  const file = Uint8Array.from({ length: size }, (_, i) => (i * 7 + 3) & 255);
  const calls = [];
  const fetchRange = (from, to) => { calls.push([from, to]); return file.slice(from, to + 1); };
  return { file, calls, fetchRange };
}

test('reads equal the file, across blocks and at its end', () => {
  const { file, fetchRange } = source();
  const rc = new RangeCache(fetchRange, file.length, { block: 64 });
  for (const [pos, len] of [[0, 10], [60, 10], [63, 130], [990, 50], [0, 1000], [999, 1]]) {
    const out = new Uint8Array(len);
    const n = rc.read(out, 0, len, pos);
    assert.equal(n, Math.min(len, file.length - pos));
    assert.deepEqual(out.subarray(0, n), file.subarray(pos, pos + n));
  }
  assert.equal(rc.read(new Uint8Array(4), 0, 4, 1000), 0);              // past the end
  const into = new Int8Array(8);                                          // (the HEAP is an Int8Array)
  rc.read(into, 2, 4, 100);
  assert.deepEqual(new Uint8Array(into.buffer).subarray(2, 6), file.subarray(100, 104));
});

test('adjacent missing blocks are one request; cached blocks are not fetched again', () => {
  const { file, calls, fetchRange } = source();
  const rc = new RangeCache(fetchRange, file.length, { block: 100 });
  rc.read(new Uint8Array(250), 0, 250, 120);           // blocks 1–3: one request
  assert.deepEqual(calls, [[100, 399]]);
  rc.read(new Uint8Array(50), 0, 50, 150);             // cached
  assert.equal(calls.length, 1);
  rc.read(new Uint8Array(500), 0, 500, 0);             // block 0 and block 4: two requests
  assert.deepEqual(calls.slice(1), [[0, 99], [400, 499]]);
  assert.equal(rc.requests, 3);
  assert.equal(rc.fetched, 500);
  // one request is at most maxRequest
  const s2 = source(); const r2 = new RangeCache(s2.fetchRange, 1000, { block: 100, maxRequest: 200 });
  r2.read(new Uint8Array(1000), 0, 1000, 0);
  assert.deepEqual(s2.calls, [[0, 199], [200, 399], [400, 599], [600, 799], [800, 999]]);
});

test('the cache stays within its limit, least recently used first', () => {
  const { file, calls, fetchRange } = source();
  const rc = new RangeCache(fetchRange, file.length, { block: 100, maxBytes: 300 });
  for (let b = 0; b < 10; b++) rc.read(new Uint8Array(10), 0, 10, b * 100);
  assert.ok(rc.bytes <= 300);
  assert.equal(rc.cache.size, 3);
  rc.read(new Uint8Array(10), 0, 10, 950);              // the latest: still there
  assert.equal(calls.length, 10);
  rc.read(new Uint8Array(10), 0, 10, 0);                // long gone: fetched again
  assert.equal(calls.length, 11);
  // a read larger than the limit keeps its own blocks
  const big = new Uint8Array(1000);
  assert.equal(rc.read(big, 0, 1000, 0), 1000);
  assert.deepEqual(big, file);
});

test('a short answer or a failure is an error', () => {
  const short = new RangeCache(() => new Uint8Array(3), 1000, { block: 100 });
  assert.throws(() => short.read(new Uint8Array(10), 0, 10, 0), /got 3 bytes, expected 100/);
  const fail = new RangeCache(() => { throw new Error('HTTP 500'); }, 1000);
  assert.throws(() => fail.read(new Uint8Array(10), 0, 10, 0), /HTTP 500/);
  assert.throws(() => new RangeCache(() => null, NaN), /size/);
});

test('URLs: the file name, links the viewer may open', () => {
  assert.equal(urlFileName('https://example.org/data/runs/carlisle%20a5.h5?x=1'), 'carlisle a5.h5');
  assert.equal(urlFileName('https://example.org/'), 'example.org');
  const base = 'https://host.example/viz/giro3d/index.html?file=x';
  assert.equal(resolveRunUrl('../../data/run.h5', base), 'https://host.example/data/run.h5');
  assert.equal(resolveRunUrl('http://other.example/r.h5', base), 'http://other.example/r.h5');
  for (const bad of ['javascript:alert(1)', 'data:text/plain,x', 'file:///c:/r.h5', '', '   ']) {
    assert.equal(resolveRunUrl(bad, base), null, bad);
  }
});
