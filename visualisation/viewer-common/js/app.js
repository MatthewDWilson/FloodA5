// app.js — the FloodA5 globe viewer, independent of the rendering engine:
// live view of a running simulation or replay of a recorded one, from
// visualisation/VisualisationServer.jl (protocol 2, docs/DATA_FORMATS.md §6).
//
// An engine (cesium/js/engine.js, giro3d/js/engine.js) provides:
//   name                                   'CesiumJS 1.147', …
//   init(container, cfg, log)              create the 3D view
//   basemaps(cfg) → [{value, label, available}]
//   defaultBasemap(cfg) → name
//   setBasemap(name) → name actually set   (falls back when unavailable)
//   photoreal(name) → bool                 basemap heights include buildings/trees
//   heightDatum(name) → 'ellipsoid' | 'geoid' | 'flat'   what the basemap's heights are
//   createSurface(mesh, {role, h0}) → surface {setDepth, setColorValues,
//                                      setStyle, show, destroy}
//   flyToBounds([w, s, e, n], height, duration)
//   flyTo(lon, lat, height, headingDeg, pitchDeg)
//   pickLonLat(x, y) → {lon, lat} | null   (canvas pixels)
//   sampleBasemapHeights([{lon, lat}], surfaces) → Float64Array (NaN: none)
//   tilesLoaded() → bool
//   credits() → html (or '' when the engine shows its own)
//   webgl2 → bool
//
// window.__flooda5 exposes the viewer's state and window.__flooda5api a small
// scripting interface (viewer-common/tools/smoke.mjs, the browser console).

import { parseMesh, decodeFloat32, FrameList, FrameCache, formatTime, PROTOCOL } from './protocol.js';
import { CellLocator, cellCentres, meshBounds, finiteMedian, finiteQuantile } from './geometry.js';
import { cssGradient, COLORMAP_NAMES, COLORMAP_LABELS } from './colormaps.js';
import { COLOR_MODE } from './surface.js';
import { UI_HTML } from './ui.js';

// source: the frame variable read; stat: its name in /api/ranges; colormap: the default palette
const VARS = {
  depth:      { label: 'Water depth',     unit: 'm',   source: 'water_depth', stat: 'water_depth', cmode: COLOR_MODE.depth,   colormap: 'depth',   fixed: null },
  level:      { label: 'Water level',     unit: 'm',   source: 'water_depth', stat: 'level',       cmode: COLOR_MODE.level,   colormap: 'viridis', fixed: null },
  speed:      { label: 'Speed',           unit: 'm/s', source: 'velocity',    stat: 'velocity',    cmode: COLOR_MODE.texture, colormap: 'magma',   fixed: null },
  saturation: { label: 'Wetted fraction', unit: '',    source: 'saturation',  stat: 'saturation',  cmode: COLOR_MODE.texture, colormap: 'blues',   fixed: [0, 1] },
};

/** Palettes offered for the water (the bed keeps 'terrain'). */
const PALETTES = COLORMAP_NAMES.filter(n => n !== 'terrain');

/** localStorage, where available (private windows, blocked storage: no-ops). */
const store = {
  get(k) { try { return localStorage.getItem(k); } catch { return null; } },
  set(k, v) { try { localStorage.setItem(k, v); } catch { /* unavailable */ } },
};

const $ = id => document.getElementById(id);

/** Start the viewer with `engine` drawing into the element #viewer. */
export function startViewer(engine) {
  const status = (window.__flooda5 = {
    state: 'starting', errors: [], meshCells: 0, frames: 0, shownFrame: null, renders: 0,
    lastRenderMs: null, basemap: null, offset: 0, complete: false, mode: null, webgl2: null,
    calibration: null, engine: engine.name, offsetSource: null, ranges: null,
  });
  const app = {
    cfg: null, info: null, mesh: null, centres: null, locator: null,
    water: null, bed: null, basemap: null, bedUser: null,
    frames: new FrameList(), cache: new FrameCache(32), pending: new Map(),
    current: -1, target: -1, live: true, complete: false, mode: null, wsDown: false,
    playing: false, playAcc: 0, simTime: 0, scrubbing: false, autoplay: false,
    varName: 'depth', range: { mode: 'run', min: 0, max: 1 }, runRanges: null, rangesAt: -Infinity,
    offset: 0, offsetJob: 0,
    depth: null, cvals: null, render: { busy: false, next: null, promise: null }, cellIds: new Map(),
  };

  document.body.insertAdjacentHTML('beforeend', UI_HTML);

  // -------------------------------------------------------------------------
  // Logging and status
  // -------------------------------------------------------------------------

  function log(msg, level = 'info') {
    const el = document.createElement('div');
    if (level !== 'info') el.className = level;
    el.textContent = `${new Date().toLocaleTimeString([], { hour12: false })}  ${msg}`;
    $('log').appendChild(el);
    while ($('log').children.length > 60) $('log').firstChild.remove();
    $('log').scrollTop = $('log').scrollHeight;
    (level === 'err' ? console.error : level === 'warn' ? console.warn : console.info)(`[FloodA5] ${msg}`);
    if (level === 'err') status.errors.push(msg);
  }

  function fatal(msg) {
    status.state = 'error';
    log(msg, 'err');
    $('fatal').textContent = msg;
    $('fatal').hidden = false;
    setStatus('error', 'error');
  }

  function setStatus(cls, label) {
    $('status-dot').className = `dot ${cls}`;
    $('status-label').textContent = label;
  }

  function refreshStatus() {
    if (status.state === 'error') return;
    if (app.wsDown) return setStatus('disconnected', app.mode === 'replay' || app.complete ? 'offline (cached)' : 'disconnected');
    if (app.mode === 'replay') return setStatus('replay', 'replay');
    if (app.complete) return setStatus('complete', 'run complete');
    return setStatus('live', 'running');
  }

  // -------------------------------------------------------------------------
  // Server access
  // -------------------------------------------------------------------------

  /**
   * fetch with retries on network errors (not HTTP errors): the browser can
   * drop a request to the local server, e.g. on a connection reused after the
   * page was busy for seconds building a large mesh.
   */
  async function fetchRetry(path, tries = 4) {
    for (let k = 0; ; k++) {
      try {
        return await fetch(path, { cache: 'no-store' });
      } catch (e) {
        if (k + 1 >= tries) throw e;
        await new Promise(r => setTimeout(r, 250 * 2 ** k));
      }
    }
  }

  async function fetchJson(path) {
    const r = await fetchRetry(path);
    if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
    return r.json();
  }

  /** Frame variable data (cached; concurrent requests share one fetch). */
  async function getVar(id, name) {
    const key = `${id}/${name}`;
    const hit = app.cache.get(key);
    if (hit) return hit;
    if (app.pending.has(key)) return app.pending.get(key);
    const p = (async () => {
      const r = await fetchRetry(`/api/frames/${id}/${name}`);
      if (!r.ok) {
        if (r.status === 404) return null;
        throw new Error(`frame ${id} ${name}: HTTP ${r.status}`);
      }
      const a = decodeFloat32(await r.arrayBuffer(), app.mesh.nCells);
      app.cache.set(key, a);
      return a;
    })().finally(() => app.pending.delete(key));
    app.pending.set(key, p);
    return p;
  }

  // -------------------------------------------------------------------------
  // Basemaps
  // -------------------------------------------------------------------------

  function populateBasemaps() {
    const sel = $('basemap');
    sel.innerHTML = '';
    for (const b of engine.basemaps(app.cfg)) {
      const o = document.createElement('option');
      o.value = b.value;
      o.textContent = b.label + (b.available ? '' : ' (needs a token)');
      o.disabled = !b.available;
      sel.appendChild(o);
    }
  }

  async function setBasemap(name) {
    let got;
    try {
      got = await engine.setBasemap(name);
    } catch (e) {
      log(`Basemap ${name} failed: ${e.message || e}`, 'err');
      got = name === 'osm' ? null : await engine.setBasemap('osm').catch(() => null);
    }
    if (!got) return;
    if (got !== name) log(`${name} is not available — using ${got}.`, 'warn');
    app.basemap = got; status.basemap = got;
    $('basemap').value = got;
    // the bed surface stands in for terrain on the flat basemaps
    if (app.bed && app.bedUser === null) setBed(got === 'osm' || got === 'none', false);
    if (app.mesh) applyOffset();
    log(`Basemap: ${$('basemap').selectedOptions[0]?.textContent || got}`);
    updateCredits();
  }

  function updateCredits() {
    $('credits').innerHTML = engine.credits ? engine.credits() : '';
  }

  // -------------------------------------------------------------------------
  // Mesh and surfaces
  // -------------------------------------------------------------------------

  async function loadMesh() {
    const r = await fetchRetry('/api/mesh');
    if (r.status === 204) { log('Waiting for the mesh…'); return false; }
    if (!r.ok) throw new Error(`/api/mesh: HTTP ${r.status}`);
    const t0 = performance.now();
    const mesh = parseMesh(await r.arrayBuffer());
    app.water?.destroy(); app.bed?.destroy();
    app.water = app.bed = null;
    app.cache.clear(); app.cellIds.clear();
    app.mesh = mesh;
    app.centres = cellCentres(mesh);
    app.locator = new CellLocator(mesh);
    app.offset = resolveOffset().value;
    const h0 = (finiteMedian(mesh.base) || 0) + app.offset;
    app.water = await engine.createSurface(mesh, { role: 'water', h0 });
    app.bed = await engine.createSurface(mesh, { role: 'bed', h0 });
    const zmin = finiteQuantile(mesh.elevation, 0.01), zmax = finiteQuantile(mesh.elevation, 0.99);
    app.bed.setStyle({ cmin: zmin, cmax: zmax, offset: app.offset });
    app.water.setStyle({ offset: app.offset, alpha: +$('opacity').value,
                         minDepth: Math.max(+$('min-depth').value || 0, 0),
                         mode: +$('display-mode').value, exag: +$('exag').value });
    const borders = $('show-borders').checked ? 1 : 0;
    app.water.setStyle({ borders }); app.bed.setStyle({ borders });
    setBed(app.bedUser ?? (app.basemap === 'osm' || app.basemap === 'none'), false);
    status.meshCells = mesh.nCells;
    status.meshLoadMs = performance.now() - t0;          // download excluded: decode + surfaces
    if (performance.memory) status.heapMB = Math.round(performance.memory.usedJSHeapSize / 2 ** 20);
    $('stat-cells').textContent = mesh.nCells.toLocaleString();
    applyVar();
    log(`Mesh: ${mesh.nCells.toLocaleString()} cells (${((performance.now() - t0) / 1000).toFixed(1)} s)`);
    flyHome(0);
    applyOffset();
    return true;
  }

  function flyHome(duration = 1.5) {
    if (!app.mesh) return;
    const h = (finiteMedian(app.mesh.base) || 0) + app.offset;
    engine.flyToBounds(meshBounds(app.mesh), h, duration);
  }

  function setBed(on, user = true) {
    if (user) app.bedUser = on;
    $('show-bed').checked = on;
    if (app.bed) app.bed.show = on;
  }

  // -------------------------------------------------------------------------
  // Vertical offset
  // -------------------------------------------------------------------------

  // Each basemap has its own offset: their heights need not share a datum
  // (Google tiles and Cesium World Terrain are ellipsoidal, Terrarium is above
  // mean sea level, a flat map is at height 0). The offset used is, in order:
  //   1. --offset (the server's vertical_offset), for an ellipsoidal basemap;
  //   2. the one saved in this browser for this basemap and place (set by
  //      hand or estimated earlier);
  //   3. a flat basemap: the lowest bed (2nd percentile) put at the map;
  //   4. otherwise estimated from the basemap's heights, now (then saved).

  function heightDatum(name = app.basemap) {
    return engine.heightDatum ? engine.heightDatum(name) : (name === 'osm' || name === 'none' ? 'flat' : 'ellipsoid');
  }

  function offsetKey(name = app.basemap) {
    const m = app.mesh;
    return `flooda5.offset.${name || ''}.${m.originLon.toFixed(4)},${m.originLat.toFixed(4)}`;
  }

  function flatOffset() {
    const z = finiteQuantile(app.mesh.elevation, 0.02);
    return Number.isFinite(z) ? -Math.round(z * 10) / 10 : 0;
  }

  /** The offset for the current basemap: {value, source, estimate}. */
  function resolveOffset() {
    const datum = heightDatum();
    const cli = app.cfg.vertical_offset;
    if (datum === 'ellipsoid' && cli !== undefined && cli !== null && Number.isFinite(+cli)) {
      return { value: +cli, source: 'from the command line (--offset)', estimate: false };
    }
    const saved = store.get(offsetKey());
    if (saved !== null && saved !== '' && Number.isFinite(+saved)) {
      return { value: +saved, source: 'saved in this browser', estimate: false };
    }
    if (datum === 'flat') return { value: flatOffset(), source: 'flat basemap: lowest bed at the map', estimate: false };
    return { value: 0, source: 'estimating…', estimate: true };
  }

  /** Apply the current basemap's offset; estimate it if there is none yet. */
  function applyOffset() {
    const r = resolveOffset();
    setOffset(r.value, { save: false, source: r.source });
    if (r.estimate) calibrateOffset({ auto: true });
  }

  function setOffset(v, { save = true, source = 'set by hand' } = {}) {
    if (!Number.isFinite(v)) return;
    app.offset = v; status.offset = v; status.offsetSource = source;
    $('offset').value = Math.round(v * 100) / 100;
    $('offset-source').textContent = `Now: ${source}.`;
    app.water?.setStyle({ offset: v });
    app.bed?.setStyle({ offset: v });
    if (save && app.mesh) store.set(offsetKey(), String(v));
  }

  /**
   * Estimate the offset: basemap height − cell elevation at sample cells (on a
   * flat basemap: the lowest bed at the map). The result is used and saved
   * for this basemap unless the basemap changed meanwhile.
   */
  async function calibrateOffset({ auto = false } = {}) {
    if (!app.mesh) return null;
    const basemap = app.basemap, job = ++app.offsetJob;
    if (heightDatum(basemap) === 'flat') {
      setOffset(flatOffset(), { source: 'flat basemap: lowest bed at the map' });
      return null;
    }
    const m = app.mesh, n = m.nCells, step = Math.max(1, Math.floor(n / 300));
    const cells = [];
    for (let i = 0; i < n; i += step) if (Number.isFinite(m.elevation[i])) cells.push(i);
    const pts = cells.map(i => ({ lon: app.centres.lon[i], lat: app.centres.lat[i] }));
    $('calibrate').disabled = true;
    log(`Estimating the vertical offset from ${cells.length} cells…`);
    try {
      const h = await engine.sampleBasemapHeights(pts, [app.water, app.bed].filter(Boolean));
      if (app.basemap !== basemap || job !== app.offsetJob) return null;     // superseded
      const d = [];
      for (let k = 0; k < cells.length; k++) if (Number.isFinite(h[k])) d.push(h[k] - m.elevation[cells[k]]);
      if (d.length < 10) {
        log(`Too few basemap heights (${d.length}) to estimate the offset.`, 'warn');
        if (auto) setOffset(app.offset, { save: false, source: 'not estimated (no basemap heights): Estimate to retry' });
        return null;
      }
      // Buildings and trees only raise photogrammetry above the bare-earth DEM:
      // a lower quantile is a better ground estimate than the median.
      const q = engine.photoreal(basemap) ? 0.3 : 0.5;
      const est = finiteQuantile(d, q);
      status.calibration = { n: d.length, of: cells.length, q, estimate: est, p10: finiteQuantile(d, 0.1),
                             p50: finiteQuantile(d, 0.5), p90: finiteQuantile(d, 0.9) };
      log(`Offset ≈ ${est.toFixed(2)} m from ${d.length} cells ` +
          `(basemap − bed: p10 ${status.calibration.p10.toFixed(1)}, p50 ${status.calibration.p50.toFixed(1)}, ` +
          `p90 ${status.calibration.p90.toFixed(1)} m)`);
      setOffset(Math.round(est * 10) / 10, { source: `estimated from the basemap (${d.length} cells)` });
      return status.calibration;
    } catch (e) {
      log(`Offset estimate failed: ${e.message || e}`, 'err');
      return null;
    } finally {
      if (job === app.offsetJob) $('calibrate').disabled = false;
    }
  }

  // -------------------------------------------------------------------------
  // Frames
  // -------------------------------------------------------------------------

  async function refreshFrames() {
    const r = await fetchJson('/api/frames');
    app.frames.set(r.frames);
    app.complete = r.complete; status.complete = r.complete;
    status.frames = app.frames.length;
    updateTimeline();
    refreshStatus();
    refreshRanges();
    if (!app.mesh || app.frames.length === 0) return;
    if (app.live) return showFrame(app.frames.length - 1);
    // keep the shown frame; a replay (or a finished run) opens at the start
    const keep = app.current >= 0 ? app.frames.indexOfId(status.shownFrame) : -1;
    const n = app.frames.length;
    return showFrame(keep >= 0 ? keep : app.target >= 0 ? Math.min(app.target, n - 1) : 0);
  }

  /** Show frame index i (latest request wins while a render is in flight). */
  function showFrame(i) {
    if (i < 0 || i >= app.frames.length) return Promise.resolve();
    app.target = i;
    app.render.next = i;
    if (app.render.busy) return app.render.promise;
    app.render.busy = true;
    app.render.promise = (async () => {
      try {
        while (app.render.next !== null) {
          const k = app.render.next; app.render.next = null;
          await renderFrame(k);
        }
      } catch (e) {
        log(`Render failed: ${e.message || e}`, 'err');
      } finally {
        app.render.busy = false;
      }
    })();
    return app.render.promise;
  }

  async function renderFrame(i) {
    const f = app.frames.at(i);
    if (!f || !app.water) return;
    const v = VARS[app.varName];
    const t0 = performance.now();
    const [depth, cv] = await Promise.all([
      getVar(f.id, 'water_depth'),
      v.cmode === COLOR_MODE.texture ? getVar(f.id, v.source) : Promise.resolve(null),
    ]);
    if (!depth) { log(`Frame ${f.id} has no water depth`, 'warn'); return; }
    app.water.setDepth(depth);
    if (cv) app.water.setColorValues(cv);
    app.depth = depth; app.cvals = cv;
    app.current = i;
    status.shownFrame = f.id; status.renders++; status.lastRenderMs = performance.now() - t0;
    updateRange();
    updateStats();
    updateTimeline();
    if (app.playing) prefetch(i + 1);
  }

  function prefetch(i) {
    const f = app.frames.at(i);
    if (!f) return;
    getVar(f.id, 'water_depth').catch(() => {});
    const v = VARS[app.varName];
    if (v.cmode === COLOR_MODE.texture) getVar(f.id, v.source).catch(() => {});
  }

  // -------------------------------------------------------------------------
  // Variable, colour range, statistics
  // -------------------------------------------------------------------------

  function availableVars() {
    const have = new Set(app.info?.vars || []);
    return Object.keys(VARS).filter(k => have.size === 0 || have.has(VARS[k].source));
  }

  function populateVars() {
    const sel = $('var'), cur = sel.value || app.varName;
    sel.innerHTML = '';
    for (const k of availableVars()) {
      const o = document.createElement('option');
      o.value = k; o.textContent = `${VARS[k].label}${VARS[k].unit ? ` (${VARS[k].unit})` : ''}`;
      sel.appendChild(o);
    }
    sel.value = availableVars().includes(cur) ? cur : 'depth';
    app.varName = sel.value;
  }

  /** The palette of variable `name` (chosen in the palette menu, else its default). */
  function paletteOf(name = app.varName) {
    const p = store.get(`flooda5.palette.${name}`);
    return PALETTES.includes(p) ? p : VARS[name].colormap;
  }

  function applyVar() {
    const v = VARS[app.varName], cmap = paletteOf();
    $('legend-bar').style.background = cssGradient(cmap);
    $('legend-bar').title = `Colour palette: ${COLORMAP_LABELS[cmap] || cmap} (click to change)`;
    $('legend-unit').textContent = v.unit;
    app.water?.setStyle({ cmode: v.cmode, colormap: cmap });
    status.palette = cmap;
  }

  function setPalette(name) {
    if (!PALETTES.includes(name)) return;
    store.set(`flooda5.palette.${app.varName}`, name);
    applyVar();
    closePaletteMenu();
  }

  function openPaletteMenu() {
    const menu = $('palette-menu'), cur = paletteOf();
    menu.innerHTML = '';
    for (const name of PALETTES) {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'palette'; b.setAttribute('role', 'option');
      b.setAttribute('aria-selected', String(name === cur));
      b.dataset.palette = name;
      const sw = document.createElement('span'); sw.className = 'swatch'; sw.style.background = cssGradient(name);
      const lab = document.createElement('span');
      lab.textContent = (COLORMAP_LABELS[name] || name) + (name === VARS[app.varName].colormap ? ' (default)' : '');
      b.append(sw, lab);
      b.onclick = () => setPalette(name);
      menu.appendChild(b);
    }
    menu.hidden = false;
    $('legend-bar').setAttribute('aria-expanded', 'true');
    menu.querySelector('[aria-selected="true"]')?.focus();
  }

  function closePaletteMenu() {
    $('palette-menu').hidden = true;
    $('legend-bar').setAttribute('aria-expanded', 'false');
  }

  function sampleWet(fn) {
    const d = app.depth, n = d.length, out = [], step = Math.max(1, Math.floor(n / 60000));
    const th = +$('min-depth').value || 0;
    for (let i = 0; i < n; i += step) if (d[i] >= th) { const x = fn(i); if (Number.isFinite(x)) out.push(x); }
    return out;
  }

  /**
   * The run's colour ranges (2nd–98th percentiles over wet cells and frames,
   * /api/ranges). A finished run's are fetched once; a live run's are
   * refreshed at most every 10 s as frames arrive.
   */
  async function refreshRanges(force = false) {
    const now = performance.now();
    if (!force && app.runRanges?.complete) return;
    if (!force && now - app.rangesAt < 10000) return;
    app.rangesAt = now;
    try {
      const th = Math.max(+$('min-depth').value || 0, 0);
      const r = await fetchJson(`/api/ranges?wet=${th}`);
      app.runRanges = r; status.ranges = r;
      if (app.range.mode === 'run') updateRange();
    } catch (e) {
      app.runRanges = null;                     // (an older server: frame ranges instead)
    }
  }

  function updateRange() {
    const v = VARS[app.varName];
    let lo, hi;
    const run = app.runRanges?.vars?.[v.stat];
    if (v.fixed) [lo, hi] = v.fixed;
    else if (app.range.mode === 'fixed') { lo = app.range.min; hi = app.range.max; }
    else if (app.range.mode === 'run' && run && Number.isFinite(run.p02) && Number.isFinite(run.p98)) {
      lo = run.p02; hi = run.p98;
      if (!(hi > lo)) hi = lo + (app.varName === 'level' ? 1 : 0.1);
      app.range.min = lo; app.range.max = hi;
    } else if (app.depth) {
      const base = app.mesh.base, d = app.depth, cv = app.cvals;
      const vals = app.varName === 'level' ? sampleWet(i => base[i] + d[i])
                 : app.varName === 'depth' ? sampleWet(i => d[i])
                 : cv ? sampleWet(i => cv[i]) : [];
      if (vals.length) {
        lo = app.varName === 'level' ? finiteQuantile(vals, 0.01) : 0;
        hi = finiteQuantile(vals, 0.99);
      } else { lo = 0; hi = 1; }
      if (!(hi > lo)) hi = lo + (app.varName === 'level' ? 1 : 0.1);
      app.range.min = lo; app.range.max = hi;
    } else return;
    app.water.setStyle({ cmin: lo, cmax: hi });
    const fmt = x => (Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : x.toFixed(2));
    $('legend-min').textContent = fmt(lo);
    $('legend-max').textContent = fmt(hi);
  }

  function updateStats() {
    const d = app.depth;
    if (!d) return;
    const th = +$('min-depth').value || 0;
    let wet = 0, max = 0;
    for (let i = 0; i < d.length; i++) { const x = d[i]; if (x >= th) wet++; if (x > max) max = x; }
    $('stat-wet').textContent = wet.toLocaleString();
    $('stat-max').textContent = `${max.toFixed(2)} m`;
    const f = app.frames.at(app.current);
    $('stat-time').textContent = f ? formatTime(f.t) : '–';
  }

  function updateTimeline() {
    const n = app.frames.length, tl = $('timeline');
    tl.max = Math.max(n - 1, 0);
    // the frame asked for (renders lag behind it); not while the slider is
    // being dragged, or it would jump back under the pointer
    const at = app.target >= 0 ? app.target : app.current;
    if (at >= 0 && !app.scrubbing) tl.value = String(at);
    $('stat-frames').textContent = n.toLocaleString();
    const f = app.frames.at(app.current);
    $('time-label').textContent = f ? `t = ${formatTime(f.t)}   frame ${app.current + 1} / ${n}` : `${n} frames`;
    const done = app.complete || app.mode === 'replay';
    $('btn-live').disabled = done;
    $('btn-live').classList.toggle('active', app.live && !done);
    $('btn-live').textContent = done ? '■ END' : '● LIVE';
  }

  // -------------------------------------------------------------------------
  // Live updates and playback
  // -------------------------------------------------------------------------

  function connectWs() {
    const ws = new WebSocket(`${location.protocol === 'https:' ? 'wss' : 'ws'}://${location.host}/api/live`);
    ws.onopen = () => { app.wsDown = false; refreshStatus(); };
    ws.onclose = () => {
      if (!app.wsDown) log('Connection to the server closed — retrying', 'warn');
      app.wsDown = true; refreshStatus();
      setTimeout(connectWs, 3000);
    };
    ws.onmessage = ev => {
      let msg;
      try { msg = JSON.parse(ev.data); } catch { return; }
      handleMessage(msg).catch(e => log(`${msg.type}: ${e.message || e}`, 'err'));
    };
  }

  async function handleMessage(msg) {
    switch (msg.type) {
      case 'hello':
        if (msg.protocol !== PROTOCOL) log(`Server protocol ${msg.protocol}, viewer ${PROTOCOL}`, 'warn');
        app.mode = msg.mode; status.mode = msg.mode;
        app.complete = msg.complete;
        if (msg.mode === 'replay' || msg.complete) app.live = false;
        await refreshInfo();
        if (msg.mesh_ready && !app.mesh) await loadMesh();
        await refreshFrames();
        status.state = 'ready';
        status.readyMs ??= Math.round(performance.now());   // page start → first frame shown
        if (app.autoplay && !app.live && !app.playing) { app.autoplay = false; startPlay(); }
        break;
      case 'mesh':
        app.runRanges = null;
        await refreshInfo();
        await loadMesh();
        await refreshFrames();
        refreshRanges(true);
        break;
      case 'frame':
        if (app.frames.add(msg.id, msg.t)) {
          status.frames = app.frames.length;
          if (app.info && !(app.info.vars || []).length) await refreshInfo();
          updateTimeline();
          if (app.live && app.mesh) showFrame(app.frames.length - 1);
          refreshRanges();
        }
        break;
      case 'reindex':
        await refreshFrames();
        break;
      case 'complete':
        app.complete = true; status.complete = true;
        log(`Run complete: ${msg.frames} frames`);
        await refreshFrames();
        refreshRanges(true);
        break;
    }
    refreshStatus();
  }

  async function refreshInfo() {
    app.info = await fetchJson('/api/info');
    const i = app.info;
    $('run-label').textContent = [i.label, i.flow_model].filter(Boolean).join(' · ');
    $('run-label').title = i.label || '';
    populateVars();
  }

  function setLive(on) {
    if (on && (app.complete || app.mode === 'replay')) on = false;
    app.live = on;
    if (on) { stopPlay(); showFrame(app.frames.length - 1); }
    updateTimeline();
  }

  function startPlay() {
    if (app.frames.length < 2) return;
    setLive(false);
    if (app.target >= app.frames.length - 1) showFrame(0);
    app.playing = true; app.playAcc = 0;
    app.simTime = app.frames.at(Math.max(app.target, 0)).t;
    $('btn-play').textContent = '⏸';
  }

  function stopPlay() {
    app.playing = false;
    $('btn-play').textContent = '▶';
  }

  function step(delta) {
    setLive(false);
    const i = Math.min(Math.max((app.target < 0 ? 0 : app.target) + delta, 0), app.frames.length - 1);
    showFrame(i);
    return i;
  }

  let lastTick = performance.now();
  function tick(now) {
    const dt = Math.min((now - lastTick) / 1000, 0.5);
    lastTick = now;
    if (app.playing && app.frames.length > 1 && !app.render.busy) {
      const [kind, val] = $('speed').value.split(':');
      if (kind === 'fps') {
        app.playAcc += dt * +val;
        if (app.playAcc >= 1) {
          app.playAcc = 0;
          if (app.target >= app.frames.length - 1) stopPlay(); else step(1);
        }
      } else {
        app.simTime += dt * +val;
        const i = app.frames.indexAtTime(app.simTime);
        if (i !== app.target) showFrame(i);
        if (app.simTime >= app.frames.last().t) stopPlay();
      }
    }
    requestAnimationFrame(tick);
  }

  // -------------------------------------------------------------------------
  // Hover read-out
  // -------------------------------------------------------------------------

  function setupPicking() {
    const canvas = engine.canvas;
    let last = 0;
    canvas.addEventListener('mousemove', ev => {
      const now = performance.now();
      if (now - last < 80) return;
      last = now;
      const r = canvas.getBoundingClientRect();
      pick(ev.clientX - r.left, ev.clientY - r.top, ev.clientX, ev.clientY);
    });
    canvas.addEventListener('mouseleave', () => { $('tooltip').hidden = true; });
  }

  function pick(x, y, cx, cy) {
    const tip = $('tooltip');
    if (!app.mesh) { tip.hidden = true; return; }
    let ll = null;
    try { ll = engine.pickLonLat(x, y); } catch { ll = null; }
    const i = ll ? app.locator.find(ll.lon, ll.lat) : -1;
    if (i < 0) { tip.hidden = true; return; }
    const m = app.mesh, d = app.depth ? app.depth[i] : NaN;
    const f2 = v => (Number.isFinite(v) ? v.toFixed(2) : '–');
    const rows = [
      ['Cell', `${i}${app.cellIds.has(i) ? ` (${app.cellIds.get(i)})` : ''}`],
      ['Depth', `${f2(d)} m`],
      ['Level', `${f2(m.base[i] + d)} m`],
      ['Bed', `${f2(m.elevation[i])} m`],
    ];
    const v = VARS[app.varName];
    if (v.cmode === COLOR_MODE.texture && app.cvals) rows.push([v.label, `${f2(app.cvals[i])} ${v.unit}`]);
    tip.innerHTML = rows.map(([k, val]) => `<div><span class="k">${k}</span> ${val}</div>`).join('');
    tip.style.left = `${cx + 14}px`; tip.style.top = `${cy + 10}px`;
    tip.hidden = false;
    if (!app.cellIds.has(i)) {
      app.cellIds.set(i, '…');
      fetchJson(`/api/cell/${i}`).then(r => app.cellIds.set(i, r.id)).catch(() => app.cellIds.delete(i));
    }
  }

  // -------------------------------------------------------------------------
  // Controls
  // -------------------------------------------------------------------------

  function wireControls() {
    $('panel-toggle').onclick = () => $('panel').classList.toggle('collapsed');
    $('var').onchange = () => {
      app.varName = $('var').value;
      applyVar();
      if (app.current >= 0) showFrame(app.current);
    };
    $('range-mode').onchange = () => {
      app.range.mode = $('range-mode').value;
      $('range-fixed').hidden = app.range.mode !== 'fixed';
      $('range-min').value = +app.range.min.toFixed(3); $('range-max').value = +app.range.max.toFixed(3);
      updateRange();
    };
    const fixedRange = () => {
      const a = +$('range-min').value, b = +$('range-max').value;
      if (Number.isFinite(a) && Number.isFinite(b) && b > a) { app.range.min = a; app.range.max = b; updateRange(); }
    };
    $('range-min').onchange = fixedRange; $('range-max').onchange = fixedRange;
    $('display-mode').onchange = () => {
      const mode = +$('display-mode').value;
      $('exag-row').hidden = mode !== 1;
      app.water?.setStyle({ mode });
    };
    $('exag').oninput = () => {
      $('exag-val').textContent = `${$('exag').value}×`;
      app.water?.setStyle({ exag: +$('exag').value });
    };
    $('opacity').oninput = () => {
      $('opacity-val').textContent = (+$('opacity').value).toFixed(2);
      app.water?.setStyle({ alpha: +$('opacity').value });
    };
    $('min-depth').onchange = () => {
      app.water?.setStyle({ minDepth: Math.max(+$('min-depth').value || 0, 0) });
      updateStats(); updateRange();
      refreshRanges(true);
    };
    $('show-bed').onchange = () => setBed($('show-bed').checked);
    $('show-borders').onchange = () => {
      const borders = $('show-borders').checked ? 1 : 0;
      app.water?.setStyle({ borders }); app.bed?.setStyle({ borders });
      store.set('flooda5.borders', String(borders));
    };
    $('legend-bar').onclick = () => ($('palette-menu').hidden ? openPaletteMenu() : closePaletteMenu());
    document.addEventListener('pointerdown', e => {
      if (!$('palette-menu').hidden && !e.target.closest('#palette-menu, #legend-bar')) closePaletteMenu();
    });
    $('palette-menu').addEventListener('keydown', e => {
      if (e.key === 'Escape') { closePaletteMenu(); $('legend-bar').focus(); }
    });
    $('basemap').onchange = () => setBasemap($('basemap').value);
    $('offset').onchange = () => setOffset(+$('offset').value);
    $('calibrate').onclick = () => calibrateOffset();
    $('home').onclick = () => flyHome();
    const tl = $('timeline');
    tl.addEventListener('pointerdown', () => { app.scrubbing = true; });
    const endScrub = () => { if (app.scrubbing) { app.scrubbing = false; updateTimeline(); } };
    tl.addEventListener('pointerup', endScrub);
    tl.addEventListener('pointercancel', endScrub);
    tl.addEventListener('change', endScrub);
    tl.oninput = () => { stopPlay(); setLive(false); showFrame(+tl.value); };
    $('btn-start').onclick = () => { stopPlay(); setLive(false); showFrame(0); };
    $('btn-end').onclick = () => { stopPlay(); setLive(false); showFrame(app.frames.length - 1); };
    $('btn-back').onclick = () => { stopPlay(); step(-1); };
    $('btn-fwd').onclick = () => { stopPlay(); step(1); };
    $('btn-play').onclick = () => (app.playing ? stopPlay() : startPlay());
    $('btn-live').onclick = () => setLive(!app.live);
    window.addEventListener('keydown', e => {
      if (e.target instanceof HTMLInputElement || e.target instanceof HTMLSelectElement ||
          e.target.closest?.('#palette-menu, #legend-bar')) return;
      if (e.key === ' ') { e.preventDefault(); app.playing ? stopPlay() : startPlay(); }
      else if (e.key === 'ArrowRight') { stopPlay(); step(1); }
      else if (e.key === 'ArrowLeft') { stopPlay(); step(-1); }
    });
  }

  // -------------------------------------------------------------------------
  // Scripting interface and boot
  // -------------------------------------------------------------------------

  window.__flooda5api = {
    engine,                    // for debugging from the console
    flyTo: (lon, lat, height, headingDeg = 0, pitchDeg = -35) => engine.flyTo(lon, lat, height, headingDeg, pitchDeg),
    home: () => flyHome(0),
    showFrame: i => { setLive(false); return showFrame(i < 0 ? app.frames.length + i : i); },
    setBasemap: name => setBasemap(name),
    setOffset: v => setOffset(v),
    calibrate: () => calibrateOffset(),
    setVar: name => { $('var').value = name; $('var').onchange(); },
    setDisplay: (mode, exag) => {
      $('display-mode').value = String(mode); $('display-mode').onchange();
      if (exag) { $('exag').value = String(exag); $('exag').oninput(); }
    },
    setBed: on => setBed(on),
    setBorders: on => { $('show-borders').checked = !!on; $('show-borders').onchange(); },
    setPalette: name => setPalette(name),
    setRangeMode: mode => { $('range-mode').value = mode; $('range-mode').onchange(); },
    setDepthBias: f => app.water?.setStyle({ depthBias: f }),
    play: () => startPlay(),
    tilesLoaded: () => engine.tilesLoaded(),
    wetCell: () => {           // the deepest cell, for pointing the camera
      if (!app.depth) return null;
      let best = -1, bd = -1;
      for (let i = 0; i < app.depth.length; i++) if (app.depth[i] > bd) { bd = app.depth[i]; best = i; }
      return best < 0 ? null : { index: best, depth: bd, lon: app.centres.lon[best], lat: app.centres.lat[best],
                                 level: app.mesh.base[best] + bd, bed: app.mesh.elevation[best] };
    },
  };

  window.addEventListener('error', e => status.errors.push(String(e.message || e)));
  window.addEventListener('unhandledrejection', e => status.errors.push(String(e.reason?.message || e.reason)));

  (async () => {
    wireControls();
    $('show-borders').checked = store.get('flooda5.borders') === '1';
    try {
      app.cfg = await fetchJson('/api/config');
      const q = new URLSearchParams(location.search);
      app.autoplay = q.has('autoplay') ? !['0', 'false', 'no'].includes(q.get('autoplay')) : !!app.cfg.autoplay;
      if (app.cfg.protocol !== PROTOCOL) log(`Server protocol ${app.cfg.protocol}, viewer ${PROTOCOL}`, 'warn');
      await engine.init($('viewer'), app.cfg, log);
      status.webgl2 = engine.webgl2;
      status.engine = engine.name;          // with its version, once loaded
      $('engine-label').textContent = engine.name;
      log(engine.name);
      setupPicking();
      populateBasemaps();
      const want = q.get('basemap') || app.cfg.default_basemap || engine.defaultBasemap(app.cfg);
      await setBasemap(want);
      connectWs();
      requestAnimationFrame(tick);
      setInterval(updateCredits, 2000);
    } catch (e) {
      console.error(e);
      fatal(`The viewer could not start: ${e.message || e}`);
    }
  })();
}
