import { test } from 'node:test';
import assert from 'node:assert/strict';
import { COLOR_MODE, SURFACE_MODE, defaultStyle, glslHelpers } from '../js/surface.js';

test('surface style defaults', () => {
  const w = defaultStyle('water'), b = defaultStyle('bed');
  assert.equal(w.mode, SURFACE_MODE.water);
  assert.equal(w.cmode, COLOR_MODE.depth);
  assert.equal(b.mode, SURFACE_MODE.bed);
  assert.equal(b.cmode, COLOR_MODE.bed);
  assert.equal(b.alpha, 1);
  assert.ok(w.alpha < 1 && w.minDepth > 0);
  assert.notEqual(defaultStyle('water'), w);         // a fresh object each time
});

test('shared GLSL helpers', () => {
  const g = glslHelpers();
  for (const f of ['fa5Packed', 'fa5Static', 'fa5Level', 'fa5ColorValue', 'fa5Normal', 'fa5MapCoord']) {
    assert.match(g, new RegExp(`\\b${f}\\(`), f);
  }
  assert.match(g, /texelFetch/);
  // the mode numbers the shaders compare against
  assert.deepEqual([SURFACE_MODE.water, SURFACE_MODE.exaggerated, SURFACE_MODE.bed], [0, 1, 2]);
  assert.deepEqual([COLOR_MODE.depth, COLOR_MODE.level, COLOR_MODE.texture, COLOR_MODE.bed], [0, 1, 2, 3]);
  // braces balance (a cheap syntax check)
  assert.equal((g.match(/\{/g) || []).length, (g.match(/\}/g) || []).length);
});
