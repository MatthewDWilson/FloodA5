// app.js — the FloodA5 globe viewer, independent of the rendering engine:
// live view of a running simulation or replay of a recorded one, from
// visualisation/VisualisationServer.jl (protocol 2, docs/DATA_FORMATS.md §6),
// or of an output file on this computer read in the browser (local/: the
// same /api answered by a worker, nothing uploaded). Without a server (the
// pages on a static web host) it starts with no data and offers to open a file.
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
//   project(lon, lat, h) → {x, y} | null   canvas pixels; null behind the globe or camera
//   sampleBasemapHeights([{lon, lat}], surfaces) → Float64Array (NaN: none)
//   tilesLoaded() → bool
//   credits() → html (or '' when the engine shows its own)
//   webgl2 → bool
//
// window.__flooda5 exposes the viewer's state and window.__flooda5api a small
// scripting interface (viewer-common/tools/smoke.mjs, the browser console).

import { parseMesh, decodeFloat32, FrameList, FrameCache, formatTime, PROTOCOL } from './protocol.js';
import { CellLocator, cellCentres, meshBounds, finiteMedian, finiteQuantile } from './geometry.js';
import { cssGradient, colormapKey, colormapLabel } from './colormaps.js';
import { initPalettes, isPalette, fillPaletteMenu, openPaletteEditor, paletteApi } from './palette-editor.js';
import { COLOR_MODE } from './surface.js';
import { UI_HTML } from './ui.js';
import { renderRun, renderAbout, shellLine, esc } from './info.js';
import { TimePlot, SERIES_COLORS, MAX_SERIES, freeSlot, isClick, seriesValues, seriesCsv, nearestIndex } from './plot.js';
import { CLASS_SETS, classify, classOf, classLabels, classedName, classColors, parseBreaks } from './classes.js';
import { DEFAULT_HAZARD, LAND_USES, PERSON_PRESETS, normalizeHazard, hazardQuery, describeHazard } from './hazard.js';
import { LocalSource, ABOUT_URL, DEFAULTS_URL } from './local/source.js';

// source: the variable read from the server; stat: its name in /api/ranges;
// colormap: the default palette; kind: 'frame' (a frame variable), 'derived'
// (a hazard variable the server derives from depth and speed, with the hazard
// parameters) or 'static' (a statistic over the run, /api/stats: no time);
// classes: its class set (classes.js), shown classed by default if
// classedDefault, always if classedOnly; time: values are times (s).
const T = COLOR_MODE.texture;
const VARS = {
  depth:      { label: 'Water depth',     unit: 'm',   source: 'water_depth', stat: 'water_depth', cmode: COLOR_MODE.depth, colormap: 'depth',   fixed: null, kind: 'frame', classes: 'depth' },
  level:      { label: 'Water level',     unit: 'm',   source: 'water_depth', stat: 'level',       cmode: COLOR_MODE.level, colormap: 'viridis', fixed: null, kind: 'frame' },
  speed:      { label: 'Speed',           unit: 'm/s', source: 'velocity',    stat: 'velocity',    cmode: T, colormap: 'magma',   fixed: null, kind: 'frame' },
  saturation: { label: 'Wetted fraction', unit: '',    source: 'saturation',  stat: 'saturation',  cmode: T, colormap: 'blues',   fixed: [0, 1], kind: 'frame' },
  dv:         { label: 'Depth × velocity', unit: 'm²/s', source: 'depth_velocity', stat: 'depth_velocity', cmode: T, colormap: 'magma', fixed: null, kind: 'derived', classes: 'depth_velocity' },
  hr:         { label: 'Hazard rating (Defra)', unit: '', source: 'hazard_rating', stat: 'hazard_rating', cmode: T, colormap: 'turbo', fixed: null, kind: 'derived', classes: 'hazard_rating', classedDefault: true },
  arr:        { label: 'Hazard class (ARR H1–H6)', unit: '', source: 'hazard_class_arr', stat: 'hazard_class_arr', cmode: T, colormap: 'turbo', fixed: null, kind: 'derived', classes: 'hazard_class_arr', classedOnly: true },
  hd:         { label: 'Hazard degree (Chen et al.)', unit: '', source: 'hazard_degree', stat: 'hazard_degree', cmode: T, colormap: 'magma', fixed: [0, 1], kind: 'derived' },
  max_depth:  { label: 'Maximum depth',   unit: 'm',   source: 'max_water_depth', cmode: T, colormap: 'depth', fixed: null, kind: 'static', classes: 'depth' },
  max_speed:  { label: 'Maximum speed',   unit: 'm/s', source: 'max_velocity', cmode: T, colormap: 'magma', fixed: null, kind: 'static' },
  max_dv:     { label: 'Maximum D·V',     unit: 'm²/s', source: 'max_depth_velocity', cmode: T, colormap: 'magma', fixed: null, kind: 'static', classes: 'depth_velocity' },
  max_hr:     { label: 'Maximum hazard rating (Defra)', unit: '', source: 'max_hazard_rating', cmode: T, colormap: 'turbo', fixed: null, kind: 'static', classes: 'hazard_rating', classedDefault: true },
  max_arr:    { label: 'Maximum hazard class (ARR)', unit: '', source: 'max_hazard_class_arr', cmode: T, colormap: 'turbo', fixed: null, kind: 'static', classes: 'hazard_class_arr', classedOnly: true },
  max_hd:     { label: 'Maximum hazard degree (Chen et al.)', unit: '', source: 'max_hazard_degree', cmode: T, colormap: 'magma', fixed: [0, 1], kind: 'static' },
  t_first:    { label: 'Time of first flooding', unit: 's', source: 'time_first_flooded', cmode: T, colormap: 'viridis', fixed: null, kind: 'static', time: true },
  duration:   { label: 'Duration of flooding', unit: 's', source: 'flood_duration', cmode: T, colormap: 'viridis', fixed: null, kind: 'static', time: true },
};
const VAR_GROUPS = [['frame', 'Each frame'], ['derived', 'Hazard (each frame)'], ['static', 'Over the run (no time)']];
const isStatic = name => VARS[name]?.kind === 'static';

/** The bed's default palette ('terrain': for the bed only; the palettes: palette-editor.js). */
const BED_DEFAULT_PALETTE = 'terrain';

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
    offset: 0, offsetJob: 0, bedRange: { mode: 'auto', min: 0, max: 1 }, tab: 'layers', about: null, run: null,
    depth: null, cvals: null, render: { busy: false, next: null, promise: null }, cellIds: new Map(),
    plot: { series: [], varName: null, dock: 'bottom', chart: null, timer: null, ctxCell: -1 }, chartLib: null,
    hazard: loadHazard(), stats: new Map(), statsAt: -Infinity, classBuf: null,
    source: null, ws: null, standalone: false, gen: 0,   // local file (LocalSource); no server; run generation
  };

  /** The hazard parameters saved in this browser (defaults otherwise). */
  function loadHazard() {
    try { return normalizeHazard(JSON.parse(store.get('flooda5.hazard') || '{}')); } catch { return normalizeHazard({}); }
  }

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
    if (app.source) return setStatus('replay', 'local file');
    if (app.standalone || app.cfg?.no_data) return setStatus('', 'no data');
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
    // a local file answers the run's requests itself; without a server, About is a static file
    if (app.source && path.startsWith('/api/') && !path.startsWith('/api/config')) return app.source.fetch(path);
    if (app.standalone && path.startsWith('/api/about')) return fetch(ABOUT_URL, { cache: 'no-cache' });
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

  /** Is `name` (a server variable) a hazard variable, which needs the hazard parameters? */
  const isDerivedSource = name => Object.values(VARS).some(v => v.kind === 'derived' && v.source === name);

  /** Frame variable data (cached; concurrent requests share one fetch). */
  async function getVar(id, name) {
    const q = isDerivedSource(name) ? `?${hazardQuery(app.hazard)}` : '';
    const key = `${id}/${name}${q}`;
    const hit = app.cache.get(key);
    if (hit) return hit;
    if (app.pending.has(key)) return app.pending.get(key);
    const gen = app.gen;
    const p = (async () => {
      const r = await fetchRetry(`/api/frames/${id}/${name}${q}`);
      if (!r.ok) {
        if (r.status === 404) return null;
        throw new Error(`frame ${id} ${name}: HTTP ${r.status}`);
      }
      const a = decodeFloat32(await r.arrayBuffer(), app.mesh.nCells);
      if (gen === app.gen) app.cache.set(key, a);       // (not a frame of the run before a file was opened)
      return a;
    })().finally(() => app.pending.delete(key));
    app.pending.set(key, p);
    return p;
  }

  /**
   * A statistic over the run (/api/stats/{name}, with the hazard parameters):
   * cached for the parameters and, live, until frames arrive.
   */
  async function getStat(name) {
    const key = `${name}?${hazardQuery(app.hazard)}`;
    if (app.stats.has(key)) return app.stats.get(key);
    const p = (async () => {
      const r = await fetchRetry(`/api/stats/${key}`);
      if (!r.ok) {
        if (r.status === 404) return null;
        throw new Error(`statistic ${name}: HTTP ${r.status}`);
      }
      return decodeFloat32(await r.arrayBuffer(), app.mesh.nCells);
    })();
    app.stats.set(key, p);
    p.catch(() => app.stats.delete(key));
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
    app.bed.setStyle({ offset: app.offset });
    applyBedStyle();
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

  /** The bed's palette (chosen in its menu, else 'terrain'). */
  function bedPalette() {
    const p = store.get('flooda5.palette.bed');
    return isPalette(p, true) ? p : BED_DEFAULT_PALETTE;
  }

  /** Bed colours: palette, range (auto: 1st–99th percentiles of the elevations, or fixed), opacity. */
  function applyBedStyle() {
    const cmap = bedPalette(), r = app.bedRange;
    if (r.mode !== 'fixed' && app.mesh) {
      let lo = finiteQuantile(app.mesh.elevation, 0.01), hi = finiteQuantile(app.mesh.elevation, 0.99);
      if (!Number.isFinite(lo)) lo = 0;
      if (!(hi > lo)) hi = lo + 1;
      r.min = lo; r.max = hi;
    }
    const alpha = Math.min(Math.max(+$('bed-opacity').value || 1, 0.1), 1);
    app.bed?.setStyle({ colormap: colormapKey(cmap), cmin: r.min, cmax: r.max, alpha });
    $('bed-legend').style.background = cssGradient(cmap);
    $('bed-legend').title = `Bed colour palette: ${colormapLabel(cmap)} (click to change)`;
    const fmt = x => (Math.abs(x) >= 100 ? x.toFixed(0) : x.toFixed(1));
    $('bed-min').textContent = fmt(r.min); $('bed-max').textContent = fmt(r.max);
    $('bed-opacity-val').textContent = alpha.toFixed(2);
    status.bed = { palette: cmap, min: r.min, max: r.max, alpha, mode: r.mode };
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

  /** The values the surface is coloured by: class indices when classed. */
  function colorValues(raw) {
    if (!raw || !classedNow()) return raw;
    const b = breaksOf(VARS[app.varName].classes);
    if (!app.classBuf || app.classBuf.length !== raw.length) app.classBuf = new Float32Array(raw.length);
    return classify(raw, b, app.classBuf);
  }

  async function renderFrame(i) {
    const f = app.frames.at(i);
    if (!f || !app.water) return;
    const v = VARS[app.varName];
    const t0 = performance.now();
    if (v.kind === 'static') return renderStatic(i, t0);
    const [depth, cv0] = await Promise.all([
      getVar(f.id, 'water_depth'),
      v.cmode === COLOR_MODE.texture ? getVar(f.id, v.source) : Promise.resolve(null),
    ]);
    if (!depth) { log(`Frame ${f.id} has no water depth`, 'warn'); return; }
    if (VARS[app.varName] !== v) return;                     // the variable changed meanwhile
    const cv = cv0 ?? (classedNow() ? depth : null);        // classed depth: coloured by class
    app.water.setDepth(depth);
    if (cv) app.water.setColorValues(colorValues(cv));
    app.depth = depth; app.cvals = cv;
    app.current = i;
    status.shownFrame = f.id; status.renders++; status.lastRenderMs = performance.now() - t0;
    updateRange();
    updateStats();
    updateTimeline();
    app.plot.chart?.setCursor(f.t);
    if (app.playing) prefetch(i + 1);
  }

  /**
   * A statistic over the run: the surface at each cell's maximum depth (cells
   * without a value, e.g. never flooded, are left out), coloured by the
   * statistic. Frame index i is kept for when a frame variable is chosen again.
   */
  async function renderStatic(i, t0) {
    const v = VARS[app.varName];
    const [hmax, cv] = await Promise.all([getStat('max_water_depth'), getStat(v.source)]);
    if (VARS[app.varName] !== v) return;
    if (!hmax || !cv) { log(`${v.label} is not available for this run`, 'warn'); return; }
    const depth = new Float32Array(hmax.length);
    for (let k = 0; k < depth.length; k++) depth[k] = Number.isFinite(cv[k]) ? hmax[k] : 0;
    app.water.setDepth(depth);
    app.water.setColorValues(colorValues(cv));
    app.depth = depth; app.cvals = cv;
    if (i >= 0) app.current = i;
    status.shownFrame = null; status.static = v.source; status.renders++;
    status.lastRenderMs = performance.now() - t0;
    updateRange();
    updateStats();
    updateTimeline();
  }

  function prefetch(i) {
    const f = app.frames.at(i);
    if (!f || isStatic(app.varName)) return;
    getVar(f.id, 'water_depth').catch(() => {});
    const v = VARS[app.varName];
    if (v.cmode === COLOR_MODE.texture) getVar(f.id, v.source).catch(() => {});
  }

  // -------------------------------------------------------------------------
  // Variable, colour range, statistics
  // -------------------------------------------------------------------------

  function availableVars() {
    const i = app.info || {};
    const have = new Set([...(i.vars || []), ...(i.derived || []), ...(i.stats || [])]);
    const known = (i.vars || []).length > 0;
    return Object.keys(VARS).filter(k => (VARS[k].kind === 'frame' && !known) || have.has(VARS[k].source));
  }

  function populateVars() {
    const sel = $('var'), cur = sel.value || app.varName;
    const avail = availableVars();
    sel.innerHTML = '';
    for (const [kind, title] of VAR_GROUPS) {
      const ks = avail.filter(k => VARS[k].kind === kind);
      if (!ks.length) continue;
      const g = document.createElement('optgroup');
      g.label = title;
      for (const k of ks) {
        const o = document.createElement('option');
        o.value = k; o.textContent = `${VARS[k].label}${VARS[k].unit && !VARS[k].time ? ` (${VARS[k].unit})` : ''}`;
        g.appendChild(o);
      }
      sel.appendChild(g);
    }
    sel.value = avail.includes(cur) ? cur : 'depth';
    app.varName = sel.value;
    $('hazard-section').hidden = !avail.some(k => VARS[k].kind !== 'frame');
  }

  // -------------------------------------------------------------------------
  // Classes and hazard parameters
  // -------------------------------------------------------------------------

  /** Whether the variable shown is drawn in classes. */
  function classedNow(name = app.varName) {
    const v = VARS[name];
    if (!v?.classes) return false;
    if (v.classedOnly) return true;
    const s = store.get(`flooda5.classed.${name}`);
    return s === null ? !!v.classedDefault : s === '1';
  }

  /** The breaks of class set `set` (the user's, else its defaults). */
  function breaksOf(set) {
    const cs = CLASS_SETS[set];
    if (!cs.editable) return cs.breaks;
    return parseBreaks(store.get(`flooda5.breaks.${set}`)) || cs.breaks;
  }

  /** Classed legend: one swatch and label per class. */
  function renderClassLegend() {
    const v = VARS[app.varName], on = classedNow();
    $('legend-classes').hidden = !on;
    $('class-row').hidden = !v.classes || !!v.classedOnly;
    $('show-classes').checked = on;
    const editable = on && CLASS_SETS[v.classes].editable;
    $('breaks-row').hidden = !editable;
    if (!on) return;
    const b = breaksOf(v.classes), labels = classLabels(v.classes, b), k = b.length + 1;
    const cols = classColors(CLASS_SETS[v.classes].ramp, k);
    const el = $('legend-classes');
    el.innerHTML = '';
    for (let i = k - 1; i >= 0; i--) {             // highest class first
      const row = document.createElement('div');
      row.className = 'class';
      const sw = document.createElement('span'); sw.className = 'swatch';
      sw.style.background = `rgb(${cols[i].join(',')})`;
      const lab = document.createElement('span'); lab.textContent = labels[i];
      row.append(sw, lab);
      el.appendChild(row);
    }
    if (editable && document.activeElement !== $('breaks')) $('breaks').value = b.join(', ');
  }

  function setClassed(on) {
    const v = VARS[app.varName];
    if (!v.classes || v.classedOnly) return;
    store.set(`flooda5.classed.${app.varName}`, on ? '1' : '0');
    applyVar();
    rerender();
  }

  function setBreaks(text) {
    const v = VARS[app.varName];
    if (!v.classes || !CLASS_SETS[v.classes].editable) return false;
    const b = text === null ? CLASS_SETS[v.classes].breaks : parseBreaks(text);
    if (!b) { log(`Breaks "${text}": give increasing numbers, e.g. 0.75, 1.25, 2`, 'warn'); return false; }
    store.set(`flooda5.breaks.${v.classes}`, b.join(', '));
    $('breaks').value = b.join(', ');
    applyVar();
    rerender();
    return true;
  }

  /** Show the shown frame (or statistic) again, e.g. after a style or parameter change. */
  function rerender() {
    if (app.current >= 0 || isStatic(app.varName)) showFrame(Math.max(app.current, 0));
  }

  function populateHazard() {
    const lu = $('hz-landuse');
    lu.innerHTML = LAND_USES.map(([k, l]) => `<option value="${k}">${l}</option>`).join('');
    $('hz-person').innerHTML = [...Object.entries(PERSON_PRESETS).map(([k, p]) => `<option value="${k}">${p.label}</option>`),
                                '<option value="custom">your own</option>'].join('');
    showHazard();
  }

  function showHazard() {
    const h = app.hazard;
    $('hz-landuse').value = h.landuse; $('hz-person').value = h.person;
    for (const k of ['mass', 'height', 'theta', 'alpha', 'beta', 'flood_depth']) $(`hz-${k}`).value = String(h[k]);
    $('hz-mass').disabled = $('hz-height').disabled = h.person !== 'custom';
    status.hazard = { ...h };
  }

  /** Set hazard parameters (a partial object), save them and redraw what depends on them. */
  function setHazard(p) {
    const next = normalizeHazard({ ...app.hazard, ...p });
    if (JSON.stringify(next) === JSON.stringify(app.hazard)) { showHazard(); return; }
    app.hazard = next;
    store.set('flooda5.hazard', JSON.stringify(next));
    showHazard();
    log(`Hazard parameters: ${describeHazard(next)}`);
    app.stats.clear();
    if (VARS[app.varName].kind !== 'frame') { refreshRanges(true); rerender(); }
    else app.runRanges && refreshRanges(true);
    if (app.plot.series.length && VARS[plotVar()].kind === 'derived') reloadSeries();
  }

  /** The palette of variable `name` (chosen in the palette menu, else its default). */
  function paletteOf(name = app.varName) {
    const p = store.get(`flooda5.palette.${name}`);
    return isPalette(p) ? p : VARS[name].colormap;
  }

  function applyVar() {
    const v = VARS[app.varName], classed = classedNow();
    const cmap = classed ? classedName(v.classes, breaksOf(v.classes).length + 1) : paletteOf();
    $('legend-bar').style.background = cssGradient(cmap);
    $('legend-bar').title = classed ? `Classes: ${CLASS_SETS[v.classes].title}`
                                    : `Colour palette: ${colormapLabel(cmap)} (click to change)`;
    $('legend-bar').disabled = classed;
    $('legend-unit').textContent = v.time ? 'time' : v.unit;
    $('legend-labels').hidden = classed;
    $('range-row').hidden = classed;
    $('range-fixed').hidden = classed || app.range.mode !== 'fixed';
    app.water?.setStyle({ cmode: classed ? COLOR_MODE.texture : v.cmode, colormap: colormapKey(cmap) });
    status.palette = cmap; status.classed = classed;
    if (v.kind !== 'static') status.static = null;
    renderClassLegend();
    updateTimeline();
  }

  function setPalette(name) {
    if (!isPalette(name) || classedNow()) return;
    store.set(`flooda5.palette.${app.varName}`, name);
    applyVar();
    closePaletteMenu();
  }

  function setBedPalette(name) {
    if (!isPalette(name, true)) return;
    store.set('flooda5.palette.bed', name);
    applyBedStyle();
    closePaletteMenu('bed');
  }

  const MENUS = {
    water: { menu: 'palette-menu', button: 'legend-bar' },
    bed: { menu: 'bed-palette-menu', button: 'bed-legend' },
  };

  function openPaletteMenu(target = 'water') {
    const bed = target === 'bed', ids = MENUS[target];
    const menu = $(ids.menu), cur = bed ? bedPalette() : paletteOf();
    const dflt = bed ? BED_DEFAULT_PALETTE : VARS[app.varName].colormap;
    fillPaletteMenu(menu, {
      bed, current: cur, dflt,
      onPick: name => (bed ? setBedPalette(name) : setPalette(name)),
      onEdit: () => { closePaletteMenu(target); openPaletteEditor({ target, current: cur }); },
    });
    menu.hidden = false;
    $(ids.button).setAttribute('aria-expanded', 'true');
    menu.querySelector('[aria-selected="true"]')?.focus();
  }

  function closePaletteMenu(target = 'water') {
    const ids = MENUS[target];
    $(ids.menu).hidden = true;
    $(ids.button).setAttribute('aria-expanded', 'false');
  }

  // -------------------------------------------------------------------------
  // Panel tabs: Layers, Run, About
  // -------------------------------------------------------------------------

  const TABS = ['layers', 'run', 'about'];

  async function showTab(name) {
    if (!TABS.includes(name)) return;
    app.tab = name; status.tab = name;
    for (const t of TABS) {
      $(`tab-${t}`).hidden = t !== name;
      $(`tab-btn-${t}`).setAttribute('aria-selected', String(t === name));
    }
    $('panel').classList.toggle('wide', name !== 'layers');
    try {
      if (name === 'run') {
        const rec = await fetchJson('/api/run');
        app.run = rec; status.run = rec;
        $('tab-run').innerHTML = renderRun(rec);
      } else if (name === 'about') {
        app.about ??= await fetchJson('/api/about');
        $('tab-about').innerHTML = renderAbout(app.about, { engine: engine.name });
      }
    } catch (e) {
      $(`tab-${name}`).innerHTML = '<p class="muted">Not available from this server.</p>';
      log(`${name}: ${e.message || e}`, 'warn');
    }
  }

  async function copyArgv() {
    const text = `julia FloodModel.jl ${shellLine(app.run?.argv)}`;
    try {
      await navigator.clipboard.writeText(text);
      log('Command line copied.');
    } catch {
      log('Could not copy (clipboard not available): select the text instead.', 'warn');
    }
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
      const r = await fetchJson(`/api/ranges?wet=${th}&${hazardQuery(app.hazard)}`);
      app.runRanges = r; status.ranges = r;
      if (app.range.mode === 'run') updateRange();
    } catch (e) {
      app.runRanges = null;                     // (an older server: frame ranges instead)
    }
  }

  function updateRange() {
    const v = VARS[app.varName];
    if (classedNow()) {                     // class indices 0 … k−1
      const k = breaksOf(v.classes).length + 1;
      app.water?.setStyle({ cmin: -0.5, cmax: k - 0.5 });
      return;
    }
    let lo, hi;
    const run = v.stat ? app.runRanges?.vars?.[v.stat] : null;    // (a statistic: its own map is the run)
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
        lo = app.varName === 'level' || v.time ? finiteQuantile(vals, 0.01) : 0;
        hi = finiteQuantile(vals, 0.99);
      } else { lo = 0; hi = 1; }
      if (!(hi > lo)) hi = lo + (app.varName === 'level' ? 1 : v.time ? 60 : 0.1);
      app.range.min = lo; app.range.max = hi;
    } else return;
    app.water.setStyle({ cmin: lo, cmax: hi });
    const fmt = v.time ? formatTime
      : x => (Math.abs(x) >= 100 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : x.toFixed(2));
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
    $('stat-time').textContent = isStatic(app.varName) ? 'whole run' : f ? formatTime(f.t) : '–';
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
    const stat = isStatic(app.varName);
    $('time-label').textContent = stat
      ? `over the run: ${n} frames${n ? `, ${formatTime(app.frames.at(0).t)} – ${formatTime(app.frames.last().t)}` : ''}`
      : f ? `t = ${formatTime(f.t)}   frame ${app.current + 1} / ${n}` : `${n} frames`;
    const done = app.complete || app.mode === 'replay';
    // a statistic has no time: no slider, no playback
    for (const id of ['timeline', 'btn-start', 'btn-back', 'btn-play', 'btn-fwd', 'btn-end', 'speed']) $(id).disabled = stat;
    $('playbar').classList.toggle('static', stat);
    $('btn-live').disabled = done || stat;
    $('btn-live').classList.toggle('active', app.live && !done);
    $('btn-live').textContent = done ? '■ END' : '● LIVE';
  }

  // -------------------------------------------------------------------------
  // Live updates and playback
  // -------------------------------------------------------------------------

  function connectWs() {
    const ws = app.ws = new WebSocket(`${location.protocol === 'https:' ? 'wss' : 'ws'}://${location.host}/api/live`);
    ws.onopen = () => { app.wsDown = false; refreshStatus(); };
    ws.onclose = () => {
      if (app.source || ws !== app.ws) return;       // a local file replaced the server's run
      if (!app.wsDown) log('Connection to the server closed — retrying', 'warn');
      app.wsDown = true; refreshStatus();
      setTimeout(connectWs, 3000);
    };
    ws.onmessage = ev => {
      if (app.source) return;
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
        closePlot();                       // cell indices refer to the old mesh
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
          if (isStatic(app.varName)) refreshStatic();
          else app.stats.clear();
          if (!isStatic(app.varName) && app.live && app.mesh) showFrame(app.frames.length - 1);
          refreshRanges();
          schedulePlotRefresh();
        }
        break;
      case 'reindex':
        app.stats.clear();
        await refreshFrames();
        schedulePlotRefresh(0);
        break;
      case 'complete':
        app.complete = true; status.complete = true;
        log(`Run complete: ${msg.frames} frames`);
        app.stats.clear();
        await refreshFrames();
        refreshRanges(true);
        schedulePlotRefresh(0);
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
    if (!$('plot').hidden) populatePlotVars();
  }

  // -------------------------------------------------------------------------
  // Local files: an output file on this computer, read in the browser
  // -------------------------------------------------------------------------

  const viewerName = () => (/giro3d/i.test(engine.name) ? 'giro3d' : 'cesium');

  /**
   * Open `file` (a File) in this page: it replaces the server's run here
   * (other pages are not affected; reload to return to the server's run).
   */
  async function openLocalFile(file) {
    if (!file) return false;
    if (!app.cfg?.h5wasm_url) { log('Cannot open files: no h5wasm_url in the configuration', 'err'); return false; }
    const src = app.source || new LocalSource({ h5wasmUrl: app.cfg.h5wasm_url, viewer: viewerName() });
    const mb = file.size / 2 ** 20;
    log(`Opening ${file.name} (${mb >= 1 ? `${mb.toFixed(0)} MB` : `${(file.size / 1024).toFixed(0)} kB`}) — read in this browser, not uploaded…`);
    status.local = { file: file.name, state: 'opening' };
    $('open-file').disabled = true;
    const t0 = performance.now();
    let fileMs;
    try {
      await src.open(file);
      fileMs = Math.round(performance.now() - t0);
    } catch (e) {
      if (!app.source) src.close();               // (else the file open before stays open)
      status.local = { file: file.name, state: 'error', error: String(e.message || e) };
      log(`Cannot open ${file.name}: ${e.message || e}`, 'err');
      return false;
    } finally {
      $('open-file').disabled = false;
    }
    const first = !app.source;
    app.source = src;
    if (first && app.ws) { app.ws.onclose = null; app.ws.close(); app.ws = null; }
    // forget the previous run
    stopPlay(); closePlot();
    app.gen++;
    app.cache.clear(); app.pending.clear(); app.stats.clear(); app.cellIds.clear();
    app.runRanges = null; app.run = null; app.current = -1; app.target = -1;
    status.shownFrame = null;
    app.mode = 'replay'; status.mode = 'replay'; app.complete = true; status.complete = true;
    app.live = false; app.wsDown = false;
    $('open-prompt').hidden = true;
    await refreshInfo();
    await loadMesh();
    await refreshFrames();
    refreshRanges(true);
    if (app.tab !== 'layers') showTab(app.tab);
    status.state = 'ready';
    status.local = { file: file.name, state: 'open', fileMs, openMs: Math.round(performance.now() - t0) };
    refreshStatus();
    log(`Local file: ${file.name}, ${app.frames.length} frames`);
    return true;
  }

  function chooseFile() {
    const inp = $('file-input');
    inp.value = '';
    inp.click();
  }

  /** Drop a file anywhere on the page to open it. */
  function setupDrop() {
    let depth = 0;
    const hasFiles = e => Array.from(e.dataTransfer?.types || []).includes('Files');
    window.addEventListener('dragenter', e => { if (!hasFiles(e)) return; e.preventDefault(); depth++; $('drop-hint').hidden = false; });
    window.addEventListener('dragover', e => { if (hasFiles(e)) e.preventDefault(); });
    window.addEventListener('dragleave', e => { if (!hasFiles(e)) return; if (--depth <= 0) { depth = 0; $('drop-hint').hidden = true; } });
    window.addEventListener('drop', e => {
      if (!hasFiles(e)) return;
      e.preventDefault(); depth = 0; $('drop-hint').hidden = true;
      const f = e.dataTransfer.files?.[0];
      if (f) openLocalFile(f);
    });
  }

  /** The configuration without a server: viewer-common/defaults.json. */
  async function standaloneConfig() {
    let d = {};
    try { d = await (await fetch(DEFAULTS_URL, { cache: 'no-cache' })).json(); } catch { /* defaults below */ }
    delete d._comment;
    return { protocol: PROTOCOL, ion_token: '', google_key: '', ...d, standalone: true };
  }

  /** A live run's statistics are recomputed as frames arrive: at most every 10 s. */
  function refreshStatic() {
    const now = performance.now();
    if (now - app.statsAt < 10000) return;
    app.statsAt = now;
    app.stats.clear();
    rerender();
  }

  function setLive(on) {
    if (on && (app.complete || app.mode === 'replay')) on = false;
    app.live = on;
    if (on) { stopPlay(); showFrame(app.frames.length - 1); }
    updateTimeline();
  }

  function startPlay() {
    if (app.frames.length < 2 || isStatic(app.varName)) return;
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
    if (app.plot.series.length) updateMarkers();
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
    if ((v.cmode === COLOR_MODE.texture || classedNow()) && app.cvals) {
      const x = app.cvals[i];
      let txt = v.time ? (Number.isFinite(x) ? formatTime(x) : 'never') : `${f2(x)} ${v.unit}`;
      if (classedNow()) {
        const c = classOf(x, breaksOf(v.classes));
        if (c >= 0) txt += ` · ${classLabels(v.classes, breaksOf(v.classes))[c]}`;
      }
      rows.push([v.label, txt]);
    }
    // (escaped: cell ids come from the file, which may be anyone's when opened locally)
    tip.innerHTML = rows.map(([k, val]) => `<div><span class="k">${esc(k)}</span> ${esc(val)}</div>`).join('');
    tip.style.left = `${cx + 14}px`; tip.style.top = `${cy + 10}px`;
    tip.hidden = false;
    if (!app.cellIds.has(i)) {
      app.cellIds.set(i, '…');
      fetchJson(`/api/cell/${i}`).then(r => app.cellIds.set(i, r.id)).catch(() => app.cellIds.delete(i));
    }
  }

  // -------------------------------------------------------------------------
  // Time series at cells (plot.js): right-click → Add to plot
  // -------------------------------------------------------------------------

  /** Chart.js (an ES module from cfg.chartjs_url), loaded when the first plot opens. */
  function loadChart() {
    app.chartLib ??= import(app.cfg.chartjs_url).then(m => {
      const Chart = m.Chart || m.default?.Chart || m.default || window.Chart;   // ES module, or the UMD build (sets window.Chart)
      if (m.registerables) Chart.register(...m.registerables);
      return Chart;
    }).catch(e => { app.chartLib = null; throw e; });
    return app.chartLib;
  }

  function plotVar() {
    const v = app.plot.varName || app.varName;
    return isStatic(v) ? 'depth' : v;            // statistics have no time series
  }

  /** The plot's variables: those with a value per frame (frame and hazard variables, not statistics). */
  function populatePlotVars() {
    const sel = $('plot-var'), cur = plotVar();
    const avail = availableVars().filter(k => !isStatic(k));
    sel.innerHTML = '';
    for (const [kind, title] of VAR_GROUPS) {
      const ks = avail.filter(k => VARS[k].kind === kind);
      if (!ks.length) continue;
      const g = document.createElement('optgroup');
      g.label = title;
      for (const k of ks) {
        const o = document.createElement('option');
        o.value = k; o.textContent = `${VARS[k].label}${VARS[k].unit ? ` (${VARS[k].unit})` : ''}`;
        g.appendChild(o);
      }
      sel.appendChild(g);
    }
    sel.value = avail.includes(cur) ? cur : 'depth';
    app.plot.varName = sel.value;
  }

  const seriesLabel = s => `cell ${s.cell}${s.id ? ` (${s.id})` : ''}`;

  /** The cell under canvas pixel (x, y), or -1. */
  function cellAt(x, y) {
    if (!app.mesh) return -1;
    let ll = null;
    try { ll = engine.pickLonLat(x, y); } catch { ll = null; }
    return ll ? app.locator.find(ll.lon, ll.lat) : -1;
  }

  async function addToPlot(i) {
    const p = app.plot;
    if (!app.mesh || !(i >= 0 && i < app.mesh.nCells)) return false;
    if (p.series.some(s => s.cell === i)) { log(`Cell ${i} is already in the plot.`); return false; }
    const slot = freeSlot(p.series.map(s => s.slot));
    if (slot < 0) { log(`At most ${MAX_SERIES} cells in the plot: remove one first.`, 'warn'); return false; }
    const s = { cell: i, slot, id: null, t: [], raw: [], source: null };
    p.series.push(s);
    renderMarkers();
    await openPlot();
    try {
      await loadSeries(s);
    } catch (e) {
      log(`Series of cell ${i}: ${e.message || e}`, 'err');
    }
    drawPlot();
    return true;
  }

  function removeFromPlot(i) {
    const p = app.plot;
    p.series = p.series.filter(s => s.cell !== i);
    if (!p.series.length) return closePlot();
    drawPlot();
  }

  function setPlotDock(dock) {
    const d = dock === 'right' ? 'right' : 'bottom';
    app.plot.dock = d;
    store.set('flooda5.plot.dock', d);
    $('plot').className = `dock-${d}`;
    const open = !$('plot').hidden;
    document.body.classList.toggle('plot-bottom', open && d === 'bottom');
    document.body.classList.toggle('plot-right', open && d === 'right');
  }

  async function openPlot() {
    $('plot').hidden = false;
    setPlotDock(app.plot.dock);
    if (app.plot.chart) return;
    populatePlotVars();
    try {
      const Chart = await loadChart();
      app.plot.chart ??= new TimePlot($('plot-canvas'), Chart, { onPick: t => goToTime(t) });
    } catch (e) {
      log(`Chart.js did not load (${app.cfg.chartjs_url}): ${e.message || e}`, 'err');
    }
  }

  function closePlot() {
    const p = app.plot;
    p.series = [];
    clearTimeout(p.timer); p.timer = null;
    $('plot').hidden = true;
    setPlotDock(p.dock);
    p.chart?.destroy(); p.chart = null;
    renderMarkers();
    status.plot = { cells: [], var: plotVar(), points: 0 };
  }

  /** Show the frame nearest to time t (a click on the plot). */
  function goToTime(t) {
    const ts = app.frames.frames.map(f => f.t), i = nearestIndex(ts, t);
    if (i < 0) return;
    stopPlay(); setLive(false); showFrame(i);
  }

  async function loadSeries(s) {
    const src = VARS[plotVar()].source;
    const q = isDerivedSource(src) ? `&${hazardQuery(app.hazard)}` : '';
    const r = await fetchJson(`/api/series/${s.cell}?var=${encodeURIComponent(src)}${q}`);
    s.id = r.id; s.t = r.t; s.raw = r.values; s.source = src;
  }

  async function reloadSeries() {
    await Promise.all(app.plot.series.map(s => loadSeries(s).catch(e => log(`Series of cell ${s.cell}: ${e.message || e}`, 'warn'))));
    drawPlot();
  }

  /** Re-read the plotted series (live frames arrive, frames reindexed): at most every 2 s. */
  function schedulePlotRefresh(ms = 2000) {
    const p = app.plot;
    if (!p.series.length || p.timer) return;
    p.timer = setTimeout(() => { p.timer = null; reloadSeries(); }, ms);
  }

  function plotLines() {
    const name = plotVar(), wet = Math.max(+$('min-depth').value || 0, 0);
    return app.plot.series.map(s => ({
      key: s.cell, label: `${s.slot + 1} ${seriesLabel(s)}`, color: SERIES_COLORS[s.slot], t: s.t,
      y: seriesValues(name, s.raw, { base: app.mesh?.base[s.cell] ?? 0, wet }),
    }));
  }

  function drawPlot() {
    const p = app.plot, name = plotVar(), v = VARS[name], lines = plotLines();
    p.chart?.setSeries(lines, { label: v.label, unit: v.unit });
    const f = app.frames.at(app.current);
    p.chart?.setCursor(f ? f.t : null);
    const leg = $('plot-legend');
    leg.innerHTML = '';
    for (const s of [...p.series].sort((a, b) => a.slot - b.slot)) {
      const e = document.createElement('span');
      e.className = 'entry';
      const line = document.createElement('span'); line.className = 'line'; line.style.background = SERIES_COLORS[s.slot];
      const num = document.createElement('span'); num.className = 'num'; num.textContent = String(s.slot + 1);
      const lab = document.createElement('span'); lab.textContent = `cell ${s.cell}`; lab.title = seriesLabel(s);
      const x = document.createElement('button');
      x.type = 'button'; x.textContent = '✕'; x.title = `Remove cell ${s.cell}`; x.dataset.cell = String(s.cell);
      e.append(line, num, lab, x);
      leg.appendChild(e);
    }
    renderMarkers();
    status.plot = { cells: p.series.map(s => s.cell), var: name, dock: p.dock,
                    points: lines.reduce((n, l) => n + l.y.filter(y => y !== null).length, 0) };
  }

  function plotCsv() {
    const lines = plotLines();
    return seriesCsv(lines.map((l, k) => ({ label: seriesLabel(app.plot.series[k]), t: l.t, y: l.y })),
                     { variable: plotVar(), unit: VARS[plotVar()].unit });
  }

  function downloadCsv() {
    if (!app.plot.series.length) return;
    const name = `${(app.info?.label || 'flooda5').replace(/[^\w.-]+/g, '_')}_${plotVar()}.csv`;
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([plotCsv()], { type: 'text/csv' }));
    a.download = name;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  // map markers: numbered circles at the plotted cells, placed every animation frame
  function renderMarkers() {
    const box = $('markers');
    box.innerHTML = '';
    for (const s of app.plot.series) {
      const m = document.createElement('div');
      m.className = 'marker'; m.textContent = String(s.slot + 1);
      m.style.background = SERIES_COLORS[s.slot]; m.dataset.cell = String(s.cell);
      m.hidden = true;
      box.appendChild(m);
    }
    if (app.plot.series.length) updateMarkers();
  }

  function updateMarkers() {
    if (!app.mesh || !engine.project) return;
    const m = app.mesh, d = app.depth, th = Math.max(+$('min-depth').value || 0, 0);
    for (const el of $('markers').children) {
      const i = +el.dataset.cell;
      const wet = d && d[i] >= th;
      const z = (wet ? m.base[i] + d[i] : m.elevation[i]) + app.offset;
      let q = null;
      try { q = engine.project(app.centres.lon[i], app.centres.lat[i], Number.isFinite(z) ? z : app.offset); } catch { q = null; }
      el.hidden = !q;
      if (q) { el.style.left = `${q.x}px`; el.style.top = `${q.y}px`; }
    }
  }

  // the right-click menu
  function openCtxMenu(cx, cy, x, y) {
    const i = cellAt(x, y), menu = $('ctx-menu');
    if (i < 0) return closeCtxMenu();
    app.plot.ctxCell = i;
    const has = app.plot.series.some(s => s.cell === i);
    $('ctx-add').textContent = `Add cell ${i} to plot`;
    $('ctx-add').disabled = has;
    $('ctx-remove').hidden = !has;
    menu.hidden = false;
    menu.style.left = `${Math.min(cx, window.innerWidth - menu.offsetWidth - 4)}px`;
    menu.style.top = `${Math.min(cy, window.innerHeight - menu.offsetHeight - 4)}px`;
    $(has ? 'ctx-remove' : 'ctx-add').focus();
  }

  function closeCtxMenu() { $('ctx-menu').hidden = true; app.plot.ctxCell = -1; }

  function setupContextMenu() {
    const canvas = engine.canvas;
    let down = null;
    // right-drag turns or zooms the view in both engines: only a right *click* opens the menu
    canvas.addEventListener('contextmenu', e => e.preventDefault());
    canvas.addEventListener('pointerdown', e => {
      down = e.button === 2 ? { x: e.clientX, y: e.clientY, t: performance.now() } : null;
      if (!$('ctx-menu').hidden) closeCtxMenu();
    }, true);
    canvas.addEventListener('pointerup', e => {
      if (e.button !== 2) return;
      const up = { x: e.clientX, y: e.clientY, t: performance.now() };
      if (isClick(down, up)) {
        const r = canvas.getBoundingClientRect();
        openCtxMenu(e.clientX, e.clientY, e.clientX - r.left, e.clientY - r.top);
      }
      down = null;
    }, true);
    $('ctx-add').onclick = () => { const i = app.plot.ctxCell; closeCtxMenu(); addToPlot(i); };
    $('ctx-remove').onclick = () => { const i = app.plot.ctxCell; closeCtxMenu(); removeFromPlot(i); };
    $('ctx-menu').addEventListener('keydown', e => { if (e.key === 'Escape') closeCtxMenu(); });
    document.addEventListener('pointerdown', e => { if (!$('ctx-menu').hidden && !e.target.closest('#ctx-menu')) closeCtxMenu(); });
  }

  // -------------------------------------------------------------------------
  // Controls
  // -------------------------------------------------------------------------

  function wireControls() {
    $('open-file').onclick = chooseFile;
    $('open-prompt-btn').onclick = chooseFile;
    $('file-input').onchange = () => { const f = $('file-input').files?.[0]; if (f) openLocalFile(f); };
    setupDrop();
    $('panel-toggle').onclick = () => $('panel').classList.toggle('collapsed');
    $('var').onchange = () => {
      const was = isStatic(app.varName);
      app.varName = $('var').value;
      if (isStatic(app.varName)) { stopPlay(); setLive(false); }
      applyVar();
      if (app.current >= 0 || isStatic(app.varName) || was) rerender();
    };
    $('show-classes').onchange = () => setClassed($('show-classes').checked);
    $('breaks').onchange = () => setBreaks($('breaks').value);
    $('breaks-reset').onclick = () => setBreaks(null);
    $('hz-landuse').onchange = () => setHazard({ landuse: $('hz-landuse').value });
    $('hz-person').onchange = () => setHazard({ person: $('hz-person').value });
    for (const k of ['mass', 'height', 'theta', 'alpha', 'beta', 'flood_depth']) {
      $(`hz-${k}`).onchange = () => setHazard({ [k]: +$(`hz-${k}`).value });
    }
    $('hz-reset').onclick = () => setHazard({ theta: DEFAULT_HAZARD.theta, alpha: DEFAULT_HAZARD.alpha, beta: DEFAULT_HAZARD.beta });
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
      if (app.plot.series.length && plotVar() === 'level') drawPlot();
    };
    $('plot-var').onchange = () => { app.plot.varName = $('plot-var').value; return reloadSeries(); };
    $('plot-dock').onclick = () => setPlotDock(app.plot.dock === 'bottom' ? 'right' : 'bottom');
    $('plot-close').onclick = () => closePlot();
    $('plot-csv').onclick = () => downloadCsv();
    $('plot-legend').addEventListener('click', e => {
      const b = e.target.closest('button[data-cell]');
      if (b) removeFromPlot(+b.dataset.cell);
    });
    $('show-bed').onchange = () => setBed($('show-bed').checked);
    $('bed-legend').onclick = () => ($('bed-palette-menu').hidden ? openPaletteMenu('bed') : closePaletteMenu('bed'));
    $('bed-range-mode').onchange = () => {
      app.bedRange.mode = $('bed-range-mode').value;
      $('bed-range-fixed').hidden = app.bedRange.mode !== 'fixed';
      $('bed-range-min').value = +app.bedRange.min.toFixed(2); $('bed-range-max').value = +app.bedRange.max.toFixed(2);
      applyBedStyle();
    };
    const bedFixed = () => {
      const a = +$('bed-range-min').value, b = +$('bed-range-max').value;
      if (Number.isFinite(a) && Number.isFinite(b) && b > a) { app.bedRange.min = a; app.bedRange.max = b; applyBedStyle(); }
    };
    $('bed-range-min').onchange = bedFixed; $('bed-range-max').onchange = bedFixed;
    $('bed-opacity').oninput = () => { applyBedStyle(); store.set('flooda5.bed.opacity', $('bed-opacity').value); };
    for (const b of document.querySelectorAll('#panel .tab')) b.onclick = () => showTab(b.dataset.tab);
    $('tab-run').addEventListener('click', e => { if (e.target.closest('[data-copy="argv"]')) copyArgv(); });
    $('show-borders').onchange = () => {
      const borders = $('show-borders').checked ? 1 : 0;
      app.water?.setStyle({ borders }); app.bed?.setStyle({ borders });
      store.set('flooda5.borders', String(borders));
    };
    $('legend-bar').onclick = () => ($('palette-menu').hidden ? openPaletteMenu() : closePaletteMenu());
    document.addEventListener('pointerdown', e => {
      for (const [target, ids] of Object.entries(MENUS)) {
        if (!$(ids.menu).hidden && !e.target.closest(`#${ids.menu}, #${ids.button}`)) closePaletteMenu(target);
      }
    });
    for (const [target, ids] of Object.entries(MENUS)) {
      $(ids.menu).addEventListener('keydown', e => {
        if (e.key === 'Escape') { closePaletteMenu(target); $(ids.button).focus(); }
      });
    }
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
          e.target.closest?.('.palette-menu, .palette-editor, .legend, .tabs, #plot, .ctx-menu')) return;
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
    setBedPalette: name => setBedPalette(name),
    palettes: paletteApi,      // user palettes and the editor (palette-editor.js)
    setBedStyle: ({ mode, min, max, alpha } = {}) => {
      if (alpha !== undefined) $('bed-opacity').value = String(alpha);
      if (mode === 'fixed' && Number.isFinite(min) && Number.isFinite(max) && max > min) {
        app.bedRange = { mode, min, max }; $('bed-range-mode').value = 'fixed';
      } else if (mode === 'auto') { app.bedRange.mode = 'auto'; $('bed-range-mode').value = 'auto'; }
      applyBedStyle();
    },
    showTab: name => showTab(name),
    addToPlot: i => addToPlot(i),                       // cell index (0-based, mesh order)
    addToPlotAt: (lon, lat) => addToPlot(app.locator ? app.locator.find(lon, lat) : -1),
    removeFromPlot: i => removeFromPlot(i),
    closePlot: () => closePlot(),
    setPlotVar: name => { $('plot-var').value = name; return $('plot-var').onchange(); },
    setPlotDock: dock => setPlotDock(dock),
    plotCsv: () => plotCsv(),
    plotPick: t => goToTime(t),
    markers: () => Array.from($('markers').children, el => ({ cell: +el.dataset.cell, shown: !el.hidden,
                                                               x: parseFloat(el.style.left), y: parseFloat(el.style.top) })),
    rightClick: (x, y) => { const r = engine.canvas.getBoundingClientRect(); openCtxMenu(r.left + x, r.top + y, x, y); },
    setRangeMode: mode => { $('range-mode').value = mode; $('range-mode').onchange(); },
    setHazard: p => setHazard(p || {}),
    openFile: file => openLocalFile(file),               // a File (tests: from an <input type=file>)
    localFetch: path => (app.source ? app.source.fetch(path).then(async r => ({ status: r.status,
                                       body: r.headers.get('Content-Type')?.startsWith('application/json') ? await r.json()
                                           : Array.from(new Uint8Array(await r.arrayBuffer())) })) : null),
    setClassed: on => setClassed(!!on),
    setBreaks: text => setBreaks(text),
    vars: () => availableVars(),
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
    populateHazard();
    $('show-borders').checked = store.get('flooda5.borders') === '1';
    app.plot.dock = store.get('flooda5.plot.dock') === 'right' ? 'right' : 'bottom';
    const bo = +store.get('flooda5.bed.opacity');
    if (bo >= 0.1 && bo <= 1) $('bed-opacity').value = String(bo);
    try {
      try {
        app.cfg = await fetchJson('/api/config');
      } catch (e) {
        // no server (a static web host): local files only
        app.cfg = await standaloneConfig();
        app.standalone = true; status.standalone = true;
      }
      const q = new URLSearchParams(location.search);
      app.autoplay = q.has('autoplay') ? !['0', 'false', 'no'].includes(q.get('autoplay')) : !!app.cfg.autoplay;
      if (app.cfg.protocol !== PROTOCOL) log(`Server protocol ${app.cfg.protocol}, viewer ${PROTOCOL}`, 'warn');
      await initPalettes({ store, fetchJson, log, onChange: () => { applyVar(); applyBedStyle(); },
                           onPick: (target, id) => (target === 'bed' ? setBedPalette(id) : setPalette(id)) });
      await engine.init($('viewer'), app.cfg, log);
      status.webgl2 = engine.webgl2;
      status.engine = engine.name;          // with its version, once loaded
      $('engine-label').textContent = engine.name;
      log(engine.name);
      setupPicking();
      setupContextMenu();
      populateBasemaps();
      const want = q.get('basemap') || app.cfg.default_basemap || engine.defaultBasemap(app.cfg);
      await setBasemap(want);
      if (app.standalone || app.cfg.no_data) {
        status.state = 'ready'; refreshStatus();
        $('run-label').textContent = 'no run open';
        $('open-prompt').hidden = false;
        log(app.standalone ? 'No server: open a FloodA5 output file from this computer.'
                           : 'Started without a run: open a FloodA5 output file from this computer.');
      } else connectWs();
      requestAnimationFrame(tick);
      setInterval(updateCredits, 2000);
    } catch (e) {
      console.error(e);
      fatal(`The viewer could not start: ${e.message || e}`);
    }
  })();
}
