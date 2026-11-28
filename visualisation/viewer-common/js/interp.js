// interp.js — smooth playback (0.9.2): the pure parts of temporal
// interpolation between output frames. No DOM, no engine: unit-tested with
// `node --test` (test/interp.test.js).
//
// Playback keeps a continuous position p in frames (frame i at p = i). Each
// animation frame advances p (by frames/s, or by simulation time per second
// through the frames' times) and shows the pair of frames either side of it,
// i0 = ⌊p⌋ and i1 = i0 + 1, blended with weight w = p − i0 in the vertex
// shader (surface.js glslSurfaceVertex). Frame k lives in texture slot k mod 2,
// so moving on to the next pair replaces one slot only (the one whose weight
// is then about 0), and the shader's uniform `w` is the weight of slot 1.
//
// The pacer decides, from what it measures, how finely to step: one blend
// per animation frame while the renderer keeps up (w quantised to the steps
// that fit in a frame interval at the current speed), and plain frame
// stepping (w = 0, as before 0.9.2) when frames take too long to arrive or the
// display cannot draw at least a few steps between frames.

/** Playback speed from the speed menu's value: 'fps:4' (frames per second) or 'x:600' (simulation seconds per second). */
export function parseSpeed(value) {
  const [kind, v] = String(value || '').split(':');
  const val = +v;
  if (!(val > 0)) return { kind: 'fps', val: 4 };
  return { kind: kind === 'x' ? 'x' : 'fps', val };
}

/** Simulation time at position p (frames: {length, at(i) → {t}}), linear between frames. */
export function timeAtPosition(frames, p) {
  const n = frames.length;
  if (n === 0) return NaN;
  if (!(p > 0)) return frames.at(0).t;
  if (p >= n - 1) return frames.at(n - 1).t;
  const i = Math.floor(p), w = p - i;
  const t0 = frames.at(i).t, t1 = frames.at(i + 1).t;
  return t0 + (t1 - t0) * w;
}

/** Position of simulation time t (the inverse of timeAtPosition; clamped to [0, n − 1]). */
export function positionAtTime(frames, t) {
  const n = frames.length;
  if (n === 0) return 0;
  if (!(t > frames.at(0).t)) return 0;
  if (t >= frames.at(n - 1).t) return n - 1;
  let lo = 0, hi = n - 1;                      // frames.at(lo).t <= t < frames.at(hi).t
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (frames.at(mid).t <= t) lo = mid; else hi = mid;
  }
  const t0 = frames.at(lo).t, t1 = frames.at(hi).t;
  return t1 > t0 ? lo + (t - t0) / (t1 - t0) : lo;
}

/** Position after dt seconds of playback at `speed` (parseSpeed), clamped to the last frame. */
export function advancePosition(frames, p, dt, speed) {
  const n = frames.length;
  if (n < 2) return 0;
  const q = speed.kind === 'x' ? positionAtTime(frames, timeAtPosition(frames, p) + dt * speed.val)
                               : p + dt * speed.val;
  return Math.min(Math.max(q, 0), n - 1);
}

/** The frames either side of position p among n: {i0, i1, w} (w: the weight of i1; 0 at the last frame). */
export function framePair(p, n) {
  if (n < 1) return { i0: -1, i1: -1, w: 0 };
  const q = Math.min(Math.max(p, 0), n - 1);
  const i0 = Math.min(Math.floor(q), n - 1);
  const i1 = Math.min(i0 + 1, n - 1);
  return { i0, i1, w: i1 > i0 ? q - i0 : 0 };
}

/** Texture slot of frame k. */
export const slotOf = k => k & 1;

/** The shader's blend uniform (the weight of slot 1) for frames i0, i0 + 1 at weight w (w = 0: frame i0 exactly). */
export function slotWeight(i0, w) {
  return slotOf(i0) ? 1 - w : w;
}

/** Wall-clock time (ms) that playback spends between frames i0 and i0 + 1 at `speed`. */
export function frameWallMs(frames, i0, speed) {
  if (speed.kind !== 'x') return 1000 / speed.val;
  if (i0 < 0 || i0 + 1 >= frames.length) return Infinity;
  return (1000 * (frames.at(i0 + 1).t - frames.at(i0).t)) / speed.val;
}

/** w rounded down to one of `steps` equal steps (fewer than 2 steps: 0, the frame itself). */
export function quantize(w, steps) {
  if (!(steps >= 2)) return 0;
  return Math.min(Math.floor(w * steps + 1e-9), steps - 1) / steps;
}

/** Exponential moving average (the first sample is taken as it is). */
export function ema(prev, x, a = 0.2) {
  return prev === null || !Number.isFinite(prev) ? x : prev + a * (x - prev);
}

/**
 * The step plan for one frame interval:
 *   frameWall  — wall time between the two frames at the current speed (ms)
 *   tickMs     — measured time of an animation frame that drew a blend (ms; null: not yet measured)
 *   readyMs    — measured time to make a frame ready: fetch (if not prefetched), decode, upload (ms)
 *   prev       — the previous plan's mode (hysteresis)
 * Returns {mode: 'interp' | 'step', steps, why}. Interpolation needs the frames
 * to arrive in well under a frame interval (< 0.6 of it to switch on, ≥ 0.9
 * to switch off) and at least 2 steps per interval. A step is one blend per
 * animation frame as long as the display keeps `minFps` doing that (capped at
 * `maxFps`); a display slower than that gets a blend every other animation
 * frame, so that blending takes at most about half the drawing time.
 */
export function planSteps({ frameWall, tickMs = null, readyMs = null, prev = 'interp', minFps = 30, maxFps = 60, maxSteps = 240 }) {
  if (!(frameWall > 0) || !Number.isFinite(frameWall)) return { mode: 'step', steps: 1, why: 'no next frame' };
  const ready = Number.isFinite(readyMs) ? readyMs : 0;
  const limit = (prev === 'step' ? 0.6 : 0.9) * frameWall;
  if (ready > limit) return { mode: 'step', steps: 1, why: 'frames slow to arrive' };
  const tick = Number.isFinite(tickMs) ? tickMs : 1000 / maxFps;
  const interval = tick <= 1000 / minFps ? Math.max(tick, 1000 / maxFps) : 2 * tick;
  const steps = Math.min(Math.floor(frameWall / interval + 1e-6), maxSteps);
  if (steps < 2) return { mode: 'step', steps: 1, why: 'frames too close together' };
  return { mode: 'interp', steps, why: '' };
}

/**
 * The pacer: the measurements behind planSteps, kept as moving averages.
 *   tick(dtMs, drew)   each animation frame while playing (drew: the frame before it set a new blend)
 *   ready(ms)          when a frame has been made ready (fetched, decoded, uploaded)
 *   plan(frameWall)    the step plan
 */
export class Pacer {
  constructor(opts = {}) {
    this.opts = opts;
    this.tickMs = null;      // animation frames that drew a blend
    this.readyMs = null;     // making a frame ready
    this.fps = null;         // animation frames per second while playing (any)
    this.mode = 'interp';
    this.why = '';
    this.steps = 0;
  }
  reset() { this.tickMs = null; this.fps = null; this.mode = 'interp'; }
  tick(dtMs, drew) {
    if (!(dtMs > 0) || dtMs > 1000) return;          // (a background tab, a debugger stop)
    this.fps = ema(this.fps, 1000 / dtMs, 0.1);
    if (drew) this.tickMs = ema(this.tickMs, dtMs);
  }
  ready(ms) { if (ms >= 0) this.readyMs = ema(this.readyMs, ms, 0.3); }
  plan(frameWall) {
    const p = planSteps({ ...this.opts, frameWall, tickMs: this.tickMs, readyMs: this.readyMs, prev: this.mode });
    this.mode = p.mode; this.why = p.why; this.steps = p.steps;
    return p;
  }
}
