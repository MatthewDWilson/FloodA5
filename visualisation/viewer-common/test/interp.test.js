import { test } from 'node:test';
import assert from 'node:assert/strict';
import { parseSpeed, timeAtPosition, positionAtTime, advancePosition, framePair, slotOf, slotWeight,
         frameWallMs, quantize, ema, planSteps, Pacer } from '../js/interp.js';

// frames as FrameList gives them: {length, at(i) → {id, t}}; uneven spacing
const list = ts => ({ length: ts.length, at: i => ({ id: i + 1, t: ts[i] }) });
const F = list([0, 60, 120, 300, 600]);
const close = (a, b, eps = 1e-9) => assert.ok(Math.abs(a - b) <= eps, `${a} vs ${b}`);

test('parseSpeed: frames/s or simulation time per second', () => {
  assert.deepEqual(parseSpeed('fps:4'), { kind: 'fps', val: 4 });
  assert.deepEqual(parseSpeed('x:600'), { kind: 'x', val: 600 });
  assert.deepEqual(parseSpeed('x:1'), { kind: 'x', val: 1 });
  assert.deepEqual(parseSpeed('nonsense'), { kind: 'fps', val: 4 });
  assert.deepEqual(parseSpeed(undefined), { kind: 'fps', val: 4 });
});

test('time ↔ position: linear between frames, clamped at the ends', () => {
  assert.equal(timeAtPosition(F, 0), 0);
  assert.equal(timeAtPosition(F, 2.5), 210);
  assert.equal(timeAtPosition(F, 3.5), 450);
  assert.equal(timeAtPosition(F, -1), 0);
  assert.equal(timeAtPosition(F, 9), 600);
  assert.ok(Number.isNaN(timeAtPosition(list([]), 0)));
  for (const p of [0, 0.25, 1, 1.5, 2.9, 3.1, 4]) close(positionAtTime(F, timeAtPosition(F, p)), p);
  assert.equal(positionAtTime(F, -5), 0);
  assert.equal(positionAtTime(F, 1e9), 4);
  // repeated times: no division by zero
  assert.equal(positionAtTime(list([0, 10, 10, 20]), 10), 2);
});

test('advancePosition: frames per second, or simulation seconds per second', () => {
  close(advancePosition(F, 0, 0.5, { kind: 'fps', val: 4 }), 2);
  close(advancePosition(F, 3.9, 1, { kind: 'fps', val: 4 }), 4);           // clamped to the last frame
  close(advancePosition(F, 0, 1, { kind: 'x', val: 30 }), 0.5);            // 30 s of a 60 s gap
  close(advancePosition(F, 2, 1, { kind: 'x', val: 90 }), 2.5);            // 90 s of a 180 s gap
  close(advancePosition(F, 1.5, 1, { kind: 'x', val: 120 }), 2 + 90 / 180); // across a frame: t 90 → 210
  assert.equal(advancePosition(list([5]), 0, 1, { kind: 'fps', val: 4 }), 0);
});

test('framePair and the texture slots', () => {
  assert.deepEqual(framePair(2.25, 5), { i0: 2, i1: 3, w: 0.25 });
  assert.deepEqual(framePair(4, 5), { i0: 4, i1: 4, w: 0 });               // the last frame exactly
  assert.deepEqual(framePair(7, 5), { i0: 4, i1: 4, w: 0 });
  assert.deepEqual(framePair(-1, 5), { i0: 0, i1: 1, w: 0 });
  assert.deepEqual(framePair(0, 0), { i0: -1, i1: -1, w: 0 });
  // frame k in slot k mod 2; the uniform is slot 1's weight
  assert.deepEqual([0, 1, 2, 3].map(slotOf), [0, 1, 0, 1]);
  assert.equal(slotWeight(2, 0.25), 0.25);   // frames 2 (slot 0) and 3 (slot 1)
  assert.equal(slotWeight(3, 0.25), 0.75);   // frames 3 (slot 1) and 4 (slot 0)
  assert.equal(slotWeight(3, 0), 1);         // frame 3 exactly
  assert.equal(slotWeight(2, 0), 0);
  // moving on from pair (2, 3) at w → 1 to (3, 4) at w = 0 keeps frame 3's weight: no jump
  close(slotWeight(2, 1), slotWeight(3, 0));
});

test('frameWallMs and quantize', () => {
  assert.equal(frameWallMs(F, 0, { kind: 'fps', val: 4 }), 250);
  assert.equal(frameWallMs(F, 2, { kind: 'x', val: 60 }), 3000);           // 180 s at 60 s/s
  assert.equal(frameWallMs(F, 4, { kind: 'x', val: 60 }), Infinity);       // no next frame
  assert.equal(quantize(0.37, 10), 0.3);
  assert.equal(quantize(0.3, 10), 0.3);                                    // (exact steps stay)
  assert.equal(quantize(0.999, 4), 0.75);
  assert.equal(quantize(0.5, 1), 0);                                       // one step: the frame itself
  assert.equal(quantize(0.5, 0), 0);
});

test('planSteps: interpolate while frames keep up, step otherwise', () => {
  // 4 frames/s on a 60 fps display: 15 blends per frame interval
  assert.deepEqual(planSteps({ frameWall: 250, tickMs: 16.7, readyMs: 20 }), { mode: 'interp', steps: 14, why: '' });
  assert.equal(planSteps({ frameWall: 250 }).steps, 15);                   // not yet measured: 60 fps assumed
  // a display faster than 60 fps is capped at 60 blends a second
  assert.equal(planSteps({ frameWall: 1000, tickMs: 7 }).steps, 60);
  // 40 fps: one blend per animation frame
  assert.equal(planSteps({ frameWall: 1000, tickMs: 25 }).steps, 40);
  // below 30 fps (heavy meshes): a blend every other animation frame
  assert.equal(planSteps({ frameWall: 1000, tickMs: 50 }).steps, 10);
  // frames slow to arrive: plain stepping, with hysteresis (0.9 to leave, 0.6 to return)
  assert.equal(planSteps({ frameWall: 250, readyMs: 240 }).mode, 'step');
  assert.equal(planSteps({ frameWall: 250, readyMs: 200 }).mode, 'interp');
  assert.equal(planSteps({ frameWall: 250, readyMs: 200, prev: 'step' }).mode, 'step');
  assert.equal(planSteps({ frameWall: 250, readyMs: 140, prev: 'step' }).mode, 'interp');
  // fewer than two steps per interval (fast playback, slow display): stepping
  assert.equal(planSteps({ frameWall: 50, tickMs: 40 }).mode, 'step');
  assert.equal(planSteps({ frameWall: Infinity }).mode, 'step');
  assert.equal(planSteps({ frameWall: 0 }).mode, 'step');
});

test('Pacer: moving averages and the plan it gives', () => {
  assert.equal(ema(null, 10), 10);
  assert.equal(ema(10, 20, 0.5), 15);
  const p = new Pacer();
  for (let k = 0; k < 50; k++) p.tick(16, true);
  close(p.tickMs, 16, 1e-6);
  assert.ok(Math.abs(p.fps - 62.5) < 0.5);
  p.tick(5000, true);                         // a background tab: ignored
  close(p.tickMs, 16, 1e-6);
  p.tick(40, false);                          // a frame that drew no blend: fps only
  close(p.tickMs, 16, 1e-6);
  assert.equal(p.plan(250).mode, 'interp');
  for (let k = 0; k < 20; k++) p.ready(400);  // frames take longer to arrive than a frame interval
  assert.equal(p.plan(250).mode, 'step');
  assert.equal(p.why, 'frames slow to arrive');
  for (let k = 0; k < 20; k++) p.ready(10);
  assert.equal(p.plan(250).mode, 'interp');
  p.reset();
  assert.equal(p.tickMs, null);
  assert.ok(p.readyMs < 20);                  // (kept: it belongs to the run, not the playback)
});
