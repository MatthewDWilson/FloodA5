import { test } from 'node:test';
import assert from 'node:assert/strict';
import { SERIES_COLORS, MAX_SERIES, freeSlot, isClick, seriesValues, nearestIndex, seriesCsv, timeTicks, timeTick } from '../js/plot.js';

test('colour slots: lowest free, a removed cell frees its slot, others keep theirs', () => {
  assert.equal(MAX_SERIES, SERIES_COLORS.length);
  assert.equal(new Set(SERIES_COLORS).size, SERIES_COLORS.length);
  assert.equal(freeSlot([]), 0);
  assert.equal(freeSlot([0, 1, 2]), 3);
  assert.equal(freeSlot([0, 2, 3]), 1);
  assert.equal(freeSlot([...Array(MAX_SERIES).keys()]), -1);
});

test('a right click is short and still; a right drag is not', () => {
  const down = { x: 100, y: 100, t: 0 };
  assert.ok(isClick(down, { x: 102, y: 103, t: 150 }));
  assert.ok(!isClick(down, { x: 120, y: 100, t: 150 }));       // dragged
  assert.ok(!isClick(down, { x: 100, y: 100, t: 900 }));       // held
  assert.ok(!isClick(null, { x: 100, y: 100, t: 10 }));        // pressed elsewhere
  assert.ok(isClick(down, { x: 108, y: 100, t: 100 }, { px: 10 }));
});

test('series values: gaps for missing values, level only where wet', () => {
  assert.deepEqual(seriesValues('depth', [0, 0.5, null, NaN]), [0, 0.5, null, null]);
  assert.deepEqual(seriesValues('level', [0, 0.02, 0.5, null], { base: 10, wet: 0.05 }), [null, null, 10.5, null]);
  assert.deepEqual(seriesValues('level', [0.05], { base: 1, wet: 0.05 }), [1.05]);
});

test('nearest frame time', () => {
  const ts = [0, 60, 120, 300];
  assert.equal(nearestIndex(ts, -5), 0);
  assert.equal(nearestIndex(ts, 29), 0);
  assert.equal(nearestIndex(ts, 31), 1);
  assert.equal(nearestIndex(ts, 200), 2);
  assert.equal(nearestIndex(ts, 211), 3);
  assert.equal(nearestIndex(ts, 1e9), 3);
  assert.equal(nearestIndex([], 5), -1);
  assert.equal(nearestIndex([7], 5), 0);
});

test('CSV over the union of times, empty where missing, quoted labels', () => {
  const csv = seriesCsv([
    { label: 'cell 3 (a5,x)', t: [0, 60], y: [0.1, null] },
    { label: 'cell 9', t: [60, 120], y: [0.25, 0.5] },
  ], { variable: 'depth', unit: 'm' });
  assert.equal(csv, 'time_s,"depth (m) cell 3 (a5,x)",depth (m) cell 9\n0,0.1,\n60,,0.25\n120,,0.5\n');
  assert.equal(seriesCsv([], {}), 'time_s\n');
});

test('time axis: round steps and their labels', () => {
  assert.deepEqual(timeTicks(0, 3600, 6), { step: 600, ticks: [0, 600, 1200, 1800, 2400, 3000, 3600] });
  assert.deepEqual(timeTicks(100, 200000, 6).step, 43200);
  assert.deepEqual(timeTicks(100, 200000, 6).ticks.slice(0, 2), [43200, 86400]);
  assert.equal(timeTicks(0, 1e8, 4).step, 25056000);            // beyond the table: whole 10-day steps
  assert.deepEqual(timeTicks(5, 5), { step: 1, ticks: [5] });
  assert.equal(timeTick(86400, 43200), '24 h');
  assert.equal(timeTick(172800, 86400), '2 d');
  assert.equal(timeTick(1200, 600), '20 min');
  assert.equal(timeTick(5400, 1800), '1:30:00');
  assert.equal(timeTick(30, 10), '30 s');
});
