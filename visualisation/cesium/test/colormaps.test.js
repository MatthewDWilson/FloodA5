import { test } from 'node:test';
import assert from 'node:assert/strict';
import { colorAt, colormapTable, cssGradient, COLORMAP_NAMES } from '../js/colormaps.js';

test('colour maps', () => {
  for (const name of COLORMAP_NAMES) {
    const t = colormapTable(name);
    assert.equal(t.length, 4 * 256);
    for (let i = 0; i < 256; i++) assert.equal(t[4 * i + 3], 255);
    assert.deepEqual(colorAt(name, -1), colorAt(name, 0));     // clamped
    assert.deepEqual(colorAt(name, 2), colorAt(name, 1));
    assert.deepEqual(colorAt(name, NaN), colorAt(name, 0));
    assert.match(cssGradient(name), /^linear-gradient\(90deg, rgb\(/);
  }
  assert.deepEqual(colorAt('depth', 0), [198, 235, 255]);
  assert.deepEqual(colorAt('depth', 1), [10, 25, 100]);
  assert.deepEqual(colorAt('nonesuch', 0), colorAt('viridis', 0));
});
