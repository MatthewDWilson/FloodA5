import { test } from 'node:test';
import assert from 'node:assert/strict';
import { COLOR_MODE, SURFACE_MODE, SHADE, defaultStyle, glslHelpers, glslFragmentHelpers,
         glslSurfaceVertex } from '../js/surface.js';

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

test('shared vertex body: two frames, the engine\'s uniform names', () => {
  const names = ['depth', 'depth1', 'cvar', 'cvar1', 'static', 'texW', 'staticW', 'corner', 'corner1', 'cornerW',
                 'grad', 'grad1', 'gradW', 'nCorners', 'w', 'minDepth', 'mode', 'exag', 'cmode', 'blend', 'smooth', 'shade'];
  const used = new Set();
  const body = glslSurfaceVertex(k => { used.add(k); return `u_${k}`; }, 'a_cell', 'a_corner');
  for (const k of names) assert.ok(used.has(k), k);
  for (const k of used) assert.ok(names.includes(k), `undocumented uniform ${k}`);
  assert.match(body, /float d = mix\(dA, dB, w\);/);                    // depths blend: the wet threshold is a crossing time
  assert.match(body, /w < 0\.5 \? cA : cB/);                           // class numbers: the nearer frame
  assert.match(body, /a_cell/); assert.match(body, /a_corner/);
  assert.equal((body.match(/\{/g) || []).length, (body.match(/\}/g) || []).length);
  // without a corner attribute: no smooth surface code
  const flat = glslSurfaceVertex(k => `u_${k}`, 'a_cell', null);
  assert.doesNotMatch(flat, /u_corner|u_grad/);
  for (const v of ['d', 'level', 'cval', 'grad']) assert.match(flat, new RegExp(`\\b(float|vec2) ${v} =`), v);
  // helpers it calls
  const g = glslHelpers();
  for (const f of ['fa5Pick', 'fa5GradMode', 'fa5Corner', 'fa5LevelC']) assert.match(g, new RegExp(`\\b${f}\\(`), f);
  const fr = glslFragmentHelpers();
  for (const f of ['fa5Shade', 'fa5ShadeSlope', 'fa5Shading']) assert.match(fr, new RegExp(`\\b${f}\\(`), f);
  assert.deepEqual([SHADE.off, SHADE.smooth, SHADE.facets], [0, 1, 2]);
  assert.equal(defaultStyle('water').shade, SHADE.smooth);
});
