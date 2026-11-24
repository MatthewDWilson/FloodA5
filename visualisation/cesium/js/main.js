// main.js — the FloodA5 globe viewer (CesiumJS): live view of a running
// simulation or replay of a recorded one, from visualisation/VisualisationServer.jl.
//
// window.__flooda5 exposes the viewer's state for scripted checks
// (tools/smoke.mjs).

import { parseMesh, decodeFloat32, FrameList, FrameCache, formatTime, PROTOCOL } from './protocol.js';
import { CellLocator, cellCentres, meshBounds, finiteMedian, finiteQuantile } from './geometry.js';
import { cssGradient } from './colormaps.js';
import { SurfaceLayer, COLOR_MODE } from './water.js';

const status = (window.__flooda5 = {
  state: 'starting', errors: [], meshCells: 0, frames: 0, shownFrame: null, renders: 0,
  lastRenderMs: null, basemap: null, offset: 0, complete: false, mode: null, webgl2: null,
  calibration: null,
});

const VARS = {
  depth:      { label: 'Water depth',     unit: 'm',   source: 'water_depth', cmode: COLOR_MODE.depth,   colormap: 'depth',   fixed: null },
  level:      { label: 'Water level',     unit: 'm',   source: 'water_depth', cmode: COLOR_MODE.level,   colormap: 'viridis', fixed: null },
  speed:      { label: 'Speed',           unit: 'm/s', source: 'velocity',    cmode: COLOR_MODE.texture, colormap: 'magma',   fixed: null },
  saturation: { label: 'Wetted fraction', unit: '',    source: 'saturation',  cmode: COLOR_MODE.texture, colormap: 'blues',   fixed: [0, 1] },
};

const $ = id => document.getElementById(id);
const app = {
  cfg: null, info: null, viewer: null, mesh: null, centres: null, locator: null,
  water: null, bed: null, tileset: null, basemap: null, bedUser: null,
  frames: new FrameList(), cache: new FrameCache(32), pending: new Map(),
  current: -1, target: -1, live: true, complete: false, playing: false, playAcc: 0, simTime: 0,
  varName: 'depth', range: { mode: 'frame', min: 0, max: 1 }, offset: 0,
  depth: null, cvals: null, render: { busy: false, next: null }, cellIds: new Map(),
};

// ---------------------------------------------------------------------------
// Logging and status
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// Server access
// ---------------------------------------------------------------------------

async function fetchJson(path) {
  const r = await fetch(path, { cache: 'no-store' });
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
    const r = await fetch(`/api/frames/${id}/${name}`, { cache: 'no-store' });
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

// ---------------------------------------------------------------------------
// CesiumJS
// ---------------------------------------------------------------------------

function loadCesium(base) {
  window.CESIUM_BASE_URL = base;
  const css = document.createElement('link');
  css.rel = 'stylesheet'; css.href = `${base}Widgets/widgets.css`;
  document.head.appendChild(css);
  return new Promise((resolve, reject) => {
    const s = document.createElement('script');
    s.src = `${base}Cesium.js`;
    s.onload = () => (window.Cesium ? resolve(window.Cesium) : reject(new Error('Cesium did not load')));
    s.onerror = () => reject(new Error(`cannot load CesiumJS from ${base}`));
    document.head.appendChild(s);
  });
}

function initViewer() {
  const C = window.Cesium;
  C.Ion.defaultAccessToken = app.cfg.ion_token || '';
  const viewer = new C.Viewer('cesiumContainer', {
    baseLayer: false,
    terrainProvider: new C.EllipsoidTerrainProvider(),
    baseLayerPicker: false, geocoder: false, homeButton: false, sceneModePicker: false,
    navigationHelpButton: false, animation: false, timeline: false, fullscreenButton: false,
    infoBox: false, selectionIndicator: false,
    contextOptions: { requestWebgl1: false },
  });
  const scene = viewer.scene;
  status.webgl2 = !!scene.context.webgl2;
  if (!scene.context.webgl2) throw new Error('This viewer needs WebGL 2.');
  scene.globe.baseColor = C.Color.fromCssColorString('#1b2731');
  scene.pickTranslucentDepth = true;
  scene.screenSpaceCameraController.enableCollisionDetection = true;
  scene.debugShowFramesPerSecond = new URLSearchParams(location.search).has('fps');
  app.viewer = viewer;
  setupPicking();
  return viewer;
}

function hasGoogleAccess() { return !!(app.cfg.ion_token || app.cfg.google_key); }

async function setBasemap(name) {
  const C = window.Cesium, viewer = app.viewer, scene = viewer.scene;
  if ((name === 'google3d' && !hasGoogleAccess()) || (name === 'terrain' && !app.cfg.ion_token)) {
    log(`${name === 'google3d' ? 'Google 3D tiles' : 'Cesium World Terrain'} need a Cesium ion token ` +
        '(FLOODA5_CESIUM_ION_TOKEN) — using OpenStreetMap.', 'warn');
    name = 'osm';
  }
  if (app.tileset) { scene.primitives.remove(app.tileset); app.tileset = null; }
  viewer.imageryLayers.removeAll();
  viewer.terrainProvider = new C.EllipsoidTerrainProvider();
  scene.globe.show = true;
  scene.globe.depthTestAgainstTerrain = false;
  try {
    if (name === 'google3d') {
      app.tileset = await C.createGooglePhotorealistic3DTileset(
        app.cfg.google_key ? { key: app.cfg.google_key, onlyUsingWithGoogleGeocoder: true }
                           : { onlyUsingWithGoogleGeocoder: true },
        { enableCollision: true });
      scene.primitives.add(app.tileset);
      scene.globe.show = false;
    } else if (name === 'terrain') {
      viewer.terrainProvider = await C.createWorldTerrainAsync();
      viewer.imageryLayers.add(osmLayer());
      scene.globe.depthTestAgainstTerrain = true;
    } else if (name === 'osm') {
      viewer.imageryLayers.add(osmLayer());
    }
  } catch (e) {
    log(`Basemap ${name} failed: ${e.message || e}`, 'err');
    if (name !== 'osm') return setBasemap('osm');
  }
  app.basemap = name; status.basemap = name;
  $('basemap').value = name;
  // the bed surface stands in for terrain on the flat basemaps
  if (app.bed && app.bedUser === null) setBed(name === 'osm' || name === 'none', false);
  log(`Basemap: ${$('basemap').selectedOptions[0].textContent}`);
}

function osmLayer() {
  const C = window.Cesium;
  return new C.ImageryLayer(new C.OpenStreetMapImageryProvider({ url: 'https://tile.openstreetmap.org/' }));
}

// ---------------------------------------------------------------------------
// Mesh and surfaces
// ---------------------------------------------------------------------------

async function loadMesh() {
  const r = await fetch('/api/mesh', { cache: 'no-store' });
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
  app.offset = loadOffset();
  const h0 = (finiteMedian(mesh.base) || 0) + app.offset;
  app.water = await SurfaceLayer.create(app.viewer, mesh, { role: 'water', h0 });
  app.bed = await SurfaceLayer.create(app.viewer, mesh, { role: 'bed', h0 });
  const zmin = finiteQuantile(mesh.elevation, 0.01), zmax = finiteQuantile(mesh.elevation, 0.99);
  app.bed.setStyle({ cmin: zmin, cmax: zmax, offset: app.offset });
  app.water.setStyle({ offset: app.offset, alpha: +$('opacity').value, minDepth: +$('min-depth').value });
  setBed(app.bedUser ?? (app.basemap === 'osm' || app.basemap === 'none'), false);
  status.meshCells = mesh.nCells;
  $('stat-cells').textContent = mesh.nCells.toLocaleString();
  $('offset').value = app.offset;
  status.offset = app.offset;
  applyVar();
  log(`Mesh: ${mesh.nCells.toLocaleString()} cells (${((performance.now() - t0) / 1000).toFixed(1)} s)`);
  flyHome(0);
  return true;
}

function flyHome(duration = 1.5) {
  const C = window.Cesium;
  if (!app.mesh) return;
  const [w, s, e, n] = meshBounds(app.mesh);
  const h = (finiteMedian(app.mesh.base) || 0) + app.offset;
  const pts = [[w, s], [e, s], [e, n], [w, n], [(w + e) / 2, (s + n) / 2]].map(([x, y]) => C.Cartesian3.fromDegrees(x, y, h));
  const sphere = C.BoundingSphere.fromPoints(pts);
  app.viewer.camera.flyToBoundingSphere(sphere, {
    offset: new C.HeadingPitchRange(0, C.Math.toRadians(-40), Math.max(sphere.radius * 2.2, 500)),
    duration,
  });
}

function setBed(on, user = true) {
  if (user) app.bedUser = on;
  $('show-bed').checked = on;
  if (app.bed) app.bed.show = on;
}

// ---------------------------------------------------------------------------
// Vertical offset
// ---------------------------------------------------------------------------

function offsetKey() {
  const m = app.mesh;
  return `flooda5.offset.${m.originLon.toFixed(4)},${m.originLat.toFixed(4)}`;
}

function loadOffset() {
  try {
    const v = localStorage.getItem(offsetKey());
    if (v !== null && Number.isFinite(+v)) return +v;
  } catch { /* storage unavailable */ }
  const d = Number(app.cfg.vertical_offset);
  return Number.isFinite(d) ? d : 0;
}

function setOffset(v) {
  if (!Number.isFinite(v)) return;
  app.offset = v; status.offset = v;
  $('offset').value = Math.round(v * 100) / 100;
  app.water?.setStyle({ offset: v });
  app.bed?.setStyle({ offset: v });
  try { localStorage.setItem(offsetKey(), String(v)); } catch { /* ignore */ }
}

/** Estimate the offset: basemap height − cell elevation at sample cells. */
async function calibrateOffset() {
  const C = window.Cesium;
  if (!app.mesh) return;
  if (app.basemap !== 'google3d' && app.basemap !== 'terrain') {
    log('Estimating the offset needs the Google 3D or World Terrain basemap.', 'warn');
    return;
  }
  const m = app.mesh, n = m.nCells, step = Math.max(1, Math.floor(n / 300));
  const cells = [];
  for (let i = 0; i < n; i += step) if (Number.isFinite(m.elevation[i])) cells.push(i);
  const cartos = cells.map(i => C.Cartographic.fromDegrees(app.centres.lon[i], app.centres.lat[i]));
  $('calibrate').disabled = true;
  log(`Estimating the vertical offset from ${cells.length} cells…`);
  try {
    const res = app.basemap === 'google3d'
      ? await app.viewer.scene.sampleHeightMostDetailed(cartos, [app.water.model, app.bed.model].filter(Boolean))
      : await C.sampleTerrainMostDetailed(app.viewer.terrainProvider, cartos);
    const d = res.map((c, k) => (c && Number.isFinite(c.height) ? c.height - m.elevation[cells[k]] : NaN))
                 .filter(Number.isFinite);
    if (d.length < 10) { log('Too few basemap heights to estimate the offset.', 'warn'); return; }
    // Buildings and trees only raise photogrammetry above the bare-earth DEM:
    // a lower quantile is a better ground estimate than the median.
    const q = app.basemap === 'google3d' ? 0.3 : 0.5;
    const est = finiteQuantile(d, q);
    status.calibration = { n: d.length, q, estimate: est, p10: finiteQuantile(d, 0.1),
                           p50: finiteQuantile(d, 0.5), p90: finiteQuantile(d, 0.9) };
    log(`Offset ≈ ${est.toFixed(2)} m from ${d.length} cells ` +
        `(basemap − bed: p10 ${status.calibration.p10.toFixed(1)}, p50 ${status.calibration.p50.toFixed(1)}, ` +
        `p90 ${status.calibration.p90.toFixed(1)} m)`);
    setOffset(Math.round(est * 10) / 10);
  } catch (e) {
    log(`Offset estimate failed: ${e.message || e}`, 'err');
  } finally {
    $('calibrate').disabled = false;
  }
}

// ---------------------------------------------------------------------------
// Frames
// ---------------------------------------------------------------------------

async function refreshFrames() {
  const r = await fetchJson('/api/frames');
  app.frames.set(r.frames);
  app.complete = r.complete; status.complete = r.complete;
  // cached arrays of frames that no longer exist are simply never asked for again
  status.frames = app.frames.length;
  updateTimeline();
  refreshStatus();
  if (!app.mesh || app.frames.length === 0) return;
  if (app.live) return showFrame(app.frames.length - 1);
  // keep the shown frame; a fresh replay opens at the end of the run
  const keep = app.current >= 0 ? app.frames.indexOfId(status.shownFrame) : -1;
  const n = app.frames.length;
  return showFrame(keep >= 0 ? keep : app.target >= 0 ? Math.min(app.target, n - 1) : n - 1);
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

// ---------------------------------------------------------------------------
// Variable, colour range, statistics
// ---------------------------------------------------------------------------

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

function applyVar() {
  const v = VARS[app.varName];
  $('legend-bar').style.background = cssGradient(v.colormap);
  $('legend-unit').textContent = v.unit;
  app.water?.setStyle({ cmode: v.cmode, colormap: v.colormap });
  if (app.range.mode === 'fixed' && !v.fixed) {
    $('range-min').value = app.range.min; $('range-max').value = app.range.max;
  }
}

function sampleWet(fn) {
  const d = app.depth, n = d.length, out = [], step = Math.max(1, Math.floor(n / 60000));
  const th = +$('min-depth').value || 0;
  for (let i = 0; i < n; i += step) if (d[i] >= th) { const x = fn(i); if (Number.isFinite(x)) out.push(x); }
  return out;
}

function updateRange() {
  const v = VARS[app.varName];
  let lo, hi;
  if (v.fixed) [lo, hi] = v.fixed;
  else if (app.range.mode === 'fixed') { lo = app.range.min; hi = app.range.max; }
  else if (app.depth) {
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
  if (app.current >= 0) tl.value = app.current;
  $('stat-frames').textContent = n.toLocaleString();
  const f = app.frames.at(app.current);
  $('time-label').textContent = f ? `t = ${formatTime(f.t)}   frame ${app.current + 1} / ${n}` : `${n} frames`;
  const done = app.complete || app.mode === 'replay';
  $('btn-live').disabled = done;
  $('btn-live').classList.toggle('active', app.live && !done);
  $('btn-live').textContent = done ? '■ END' : '● LIVE';
}

// ---------------------------------------------------------------------------
// Live updates and playback
// ---------------------------------------------------------------------------

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
      if (!app.live && app.current < 0 && app.frames.length) await showFrame(app.frames.length - 1);
      status.state = 'ready';
      break;
    case 'mesh':
      await refreshInfo();
      await loadMesh();
      await refreshFrames();
      break;
    case 'frame':
      if (app.frames.add(msg.id, msg.t)) {
        status.frames = app.frames.length;
        if (app.info && !(app.info.vars || []).length) { await refreshInfo(); populateVars(); }
        updateTimeline();
        if (app.live && app.mesh) showFrame(app.frames.length - 1);
      }
      break;
    case 'reindex':
      await refreshFrames();
      break;
    case 'complete':
      app.complete = true; status.complete = true;
      log(`Run complete: ${msg.frames} frames`);
      await refreshFrames();
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

// ---------------------------------------------------------------------------
// Hover read-out
// ---------------------------------------------------------------------------

function setupPicking() {
  const C = window.Cesium, viewer = app.viewer;
  const handler = new C.ScreenSpaceEventHandler(viewer.scene.canvas);
  let lastPick = 0, pendingPos = null;
  handler.setInputAction(mv => {
    pendingPos = mv.endPosition;
    const now = performance.now();
    if (now - lastPick < 80) return;
    lastPick = now;
    pick(pendingPos);
  }, C.ScreenSpaceEventType.MOUSE_MOVE);
  viewer.scene.canvas.addEventListener('mouseleave', () => { $('tooltip').hidden = true; });
}

function pick(pos) {
  const C = window.Cesium, tip = $('tooltip');
  if (!app.mesh || !pos) { tip.hidden = true; return; }
  const scene = app.viewer.scene;
  let p = scene.pickPositionSupported ? scene.pickPosition(pos) : undefined;
  if (!C.defined(p)) p = app.viewer.camera.pickEllipsoid(pos);
  if (!C.defined(p)) { tip.hidden = true; return; }
  const c = C.Cartographic.fromCartesian(p);
  const i = app.locator.find(C.Math.toDegrees(c.longitude), C.Math.toDegrees(c.latitude));
  if (i < 0) { tip.hidden = true; return; }
  const m = app.mesh, d = app.depth ? app.depth[i] : NaN;
  const f2 = x => (Number.isFinite(x) ? x.toFixed(2) : '–');
  const rows = [
    ['Cell', `${i}${app.cellIds.has(i) ? ` (${app.cellIds.get(i)})` : ''}`],
    ['Depth', `${f2(d)} m`],
    ['Level', `${f2(m.base[i] + d)} m`],
    ['Bed', `${f2(m.elevation[i])} m`],
  ];
  const v = VARS[app.varName];
  if (v.cmode === COLOR_MODE.texture && app.cvals) rows.push([v.label, `${f2(app.cvals[i])} ${v.unit}`]);
  tip.innerHTML = rows.map(([k, x]) => `<div><span class="k">${k}</span> ${x}</div>`).join('');
  tip.style.left = `${pos.x + 14}px`; tip.style.top = `${pos.y + 10}px`;
  tip.hidden = false;
  if (!app.cellIds.has(i)) {
    app.cellIds.set(i, '…');
    fetchJson(`/api/cell/${i}`).then(r => app.cellIds.set(i, r.id)).catch(() => app.cellIds.delete(i));
  }
}

// ---------------------------------------------------------------------------
// Controls
// ---------------------------------------------------------------------------

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
  };
  $('show-bed').onchange = () => setBed($('show-bed').checked);
  $('basemap').onchange = () => setBasemap($('basemap').value);
  $('offset').onchange = () => setOffset(+$('offset').value);
  $('calibrate').onclick = () => calibrateOffset();
  $('home').onclick = () => flyHome();
  $('timeline').oninput = () => { stopPlay(); setLive(false); showFrame(+$('timeline').value); };
  $('btn-start').onclick = () => { stopPlay(); setLive(false); showFrame(0); };
  $('btn-end').onclick = () => { stopPlay(); setLive(false); showFrame(app.frames.length - 1); };
  $('btn-back').onclick = () => { stopPlay(); step(-1); };
  $('btn-fwd').onclick = () => { stopPlay(); step(1); };
  $('btn-play').onclick = () => (app.playing ? stopPlay() : startPlay());
  $('btn-live').onclick = () => setLive(!app.live);
  window.addEventListener('keydown', e => {
    if (e.target instanceof HTMLInputElement || e.target instanceof HTMLSelectElement) return;
    if (e.key === ' ') { e.preventDefault(); app.playing ? stopPlay() : startPlay(); }
    else if (e.key === 'ArrowRight') { stopPlay(); step(1); }
    else if (e.key === 'ArrowLeft') { stopPlay(); step(-1); }
  });
}

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------

async function boot() {
  wireControls();
  try {
    app.cfg = await fetchJson('/api/config');
    if (app.cfg.protocol !== PROTOCOL) log(`Server protocol ${app.cfg.protocol}, viewer ${PROTOCOL}`, 'warn');
    await loadCesium(app.cfg.cesium_url);
    log(`CesiumJS ${window.Cesium.VERSION}`);
    initViewer();
    const want = new URLSearchParams(location.search).get('basemap') || app.cfg.default_basemap ||
                 (hasGoogleAccess() ? 'google3d' : 'osm');
    await setBasemap(want);
    if (!hasGoogleAccess()) log('No Cesium ion token: Google 3D tiles and World Terrain are unavailable.', 'warn');
    connectWs();
    requestAnimationFrame(tick);
  } catch (e) {
    console.error(e);
    fatal(`The viewer could not start: ${e.message || e}`);
  }
}

// A small scripting interface (tools/smoke.mjs, the browser console).
window.__flooda5api = {
  flyTo(lon, lat, height, headingDeg = 0, pitchDeg = -35, duration = 0) {
    const C = window.Cesium;
    app.viewer.camera.flyTo({
      destination: C.Cartesian3.fromDegrees(lon, lat, height),
      orientation: { heading: C.Math.toRadians(headingDeg), pitch: C.Math.toRadians(pitchDeg), roll: 0 },
      duration,
    });
  },
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
  tilesLoaded: () => (app.tileset ? app.tileset.tilesLoaded : app.viewer.scene.globe.tilesLoaded),
  wetCell: () => {           // a wet cell (the deepest), for pointing the camera
    if (!app.depth) return null;
    let best = -1, bd = -1;
    for (let i = 0; i < app.depth.length; i++) if (app.depth[i] > bd) { bd = app.depth[i]; best = i; }
    return best < 0 ? null : { index: best, depth: bd, lon: app.centres.lon[best], lat: app.centres.lat[best],
                               level: app.mesh.base[best] + bd, bed: app.mesh.elevation[best] };
  },
};

window.addEventListener('error', e => status.errors.push(String(e.message || e)));
window.addEventListener('unhandledrejection', e => status.errors.push(String(e.reason?.message || e.reason)));
boot();
