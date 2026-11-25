import { test } from 'node:test';
import assert from 'node:assert/strict';
import { colorAt, colormapTable, cssGradient, COLORMAP_NAMES, COLORMAP_LABELS } from '../js/colormaps.js';

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

test('palette menu: every map has a label; sediment runs light to dark brown', () => {
  for (const name of COLORMAP_NAMES) assert.ok(COLORMAP_LABELS[name], name);
  const [r0, g0, b0] = colorAt('sediment', 0), [r1, g1, b1] = colorAt('sediment', 1);
  assert.ok(r0 > g0 && g0 > b0 && r1 > g1 && g1 > b1);          // brown hues
  assert.ok(r0 + g0 + b0 > 2 * (r1 + g1 + b1));                  // light → dark
});
