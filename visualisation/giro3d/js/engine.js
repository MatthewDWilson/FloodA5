// engine.js — the Giro3D engine for the globe viewer (the interface:
// viewer-common/js/app.js). Giro3D 2.0.7 on three.js r180, loaded without a
// build step through the import map in ../index.html.
//
// Basemaps: Google Photorealistic 3D Tiles (a Tiles3D entity with the Cesium
// ion or Google auth plugin of 3d-tiles-renderer), Terrarium terrain + OSM on
// a Giro3D Globe, OSM, none.

import * as THREE from 'three';
import Instance from '@giro3d/giro3d/core/Instance.js';
import CoordinateSystem from '@giro3d/giro3d/core/geographic/CoordinateSystem.js';
import Ellipsoid from '@giro3d/giro3d/core/geographic/Ellipsoid.js';
import Globe from '@giro3d/giro3d/entities/Globe.js';
import Atmosphere from '@giro3d/giro3d/entities/Atmosphere.js';
import Tiles3D from '@giro3d/giro3d/entities/Tiles3D.js';
import ColorLayer from '@giro3d/giro3d/core/layer/ColorLayer.js';
import ElevationLayer from '@giro3d/giro3d/core/layer/ElevationLayer.js';
import TiledImageSource from '@giro3d/giro3d/sources/TiledImageSource.js';
import GlobeControls from '@giro3d/giro3d/controls/GlobeControls.js';
import OSM from 'ol/source/OSM.js';
import XYZ from 'ol/source/XYZ.js';
import { CesiumIonAuthPlugin, GoogleCloudAuthPlugin } from '3d-tiles-renderer/core/plugins';
import { TerrariumFormat, TERRARIUM_URL, TERRARIUM_MAX_ZOOM, sampleTerrarium } from './terrarium.js';
import { SurfaceMesh } from './water.js';
import { ecefToGeodetic, WGS84_A, WGS84_B } from '../../viewer-common/js/geometry.js';

/** The official Google Maps logo with a dark outline (Google's attribution assets: viewer-common/img). */
const GOOGLE_LOGO = new URL('../../viewer-common/img/GoogleMaps_Logo_WithDarkOutline.svg', import.meta.url).href;

const ELL_A2 = WGS84_A * WGS84_A, ELL_B2 = WGS84_B * WGS84_B;   // (the geodetic normal)

export const GIRO3D_VERSION = '2.0.7';         // keep in step with ../index.html
const GOOGLE_ION_ASSET = '2275207';
const ELL = Ellipsoid.WGS84;
const DEG = Math.PI / 180;

export class Giro3DEngine {
  constructor() {
    this.name = `Giro3D ${GIRO3D_VERSION} (three.js r${THREE.REVISION})`;
    this.instance = null;
    this.globe = null;
    this.tiles = null;
    this.layers = [];
    this.basemap = null;
    this.webgl2 = null;
  }

  async init(container, cfg, log) {
    this.cfg = cfg; this.log = log;
    const instance = new Instance({ target: container, crs: CoordinateSystem.epsg4978, backgroundColor: '#0b1116' });
    this.instance = instance;
    this.canvas = instance.renderer.domElement;
    this.webgl2 = instance.renderer.capabilities.isWebGL2 !== false &&
                  typeof WebGL2RenderingContext !== 'undefined' &&
                  instance.renderer.getContext() instanceof WebGL2RenderingContext;
    if (!this.webgl2) throw new Error('This viewer needs WebGL 2.');
    this.globe = new Globe({ backgroundColor: '#1b2731' });
    await instance.add(this.globe);
    // sky (0.8.15): the atmosphere's halo seen from afar, the sky dome near the ground
    this.atmosphere = new Atmosphere();
    this.atmosphere.inner.visible = false;     // (its veil over the globe washes the light basemaps out: the halo only)
    await instance.add(this.atmosphere);
    this.sky = skyGradient();
    await instance.add(this.sky);
    this.setSun(0, 30);
    const camera = instance.view.camera;
    camera.up.set(0, 0, 1);
    camera.position.set(2.5e7, 0, 0);
    camera.lookAt(0, 0, 0);
    this.controls = new GlobeControls({
      scene: instance.scene, camera, domElement: instance.domElement, enableDamping: true,
    });
    this.controls.addEventListener('change', () => instance.notifyChange(camera));
    const loop = () => {
      this.controls.update();
      const alt = ecefToGeodetic(camera.position.x, camera.position.y, camera.position.z).h;
      instance.view.minNearPlane = Math.min(Math.max(alt / 200, 0.5), 500);
      // the sky: the local vertical at the camera; gone above ~150 km (space: the atmosphere's halo)
      const u = this.sky.material.uniforms;
      u.uUp.value.set(camera.position.x / ELL_A2, camera.position.y / ELL_A2, camera.position.z / ELL_B2).normalize();
      u.uAlpha.value = 1 - THREE.MathUtils.smoothstep(alt, 30e3, 150e3);
      this.sky.visible = u.uAlpha.value > 0;
      // the atmosphere's halo is for views from space: inside its shell it would whiten the sky
      const space = alt > 120e3;
      if (this.atmosphere.visible !== space) { this.atmosphere.visible = space; instance.notifyChange(); }
      if (this.tiles) instance.notifyChange(this.tiles);    // keep tiles streaming
      requestAnimationFrame(loop);
    };
    requestAnimationFrame(loop);
  }

  basemaps(cfg) {
    const google = !!(cfg.ion_token || cfg.google_key);
    return [
      { value: 'google3d', label: 'Google Photorealistic 3D', available: google },
      { value: 'terrain', label: 'Terrarium terrain + OSM', available: true },
      { value: 'osm', label: 'OpenStreetMap (flat)', available: true },
      { value: 'none', label: 'None', available: true },
    ];
  }

  defaultBasemap(cfg) { return cfg.ion_token || cfg.google_key ? 'google3d' : 'osm'; }

  photoreal(name) { return name === 'google3d'; }

  /**
   * The height reference of a basemap: 'ellipsoid' (Google tiles: heights
   * above WGS84), 'geoid' (Terrarium: heights above mean sea level, drawn as
   * if above the ellipsoid) or 'flat'.
   */
  heightDatum(name) { return name === 'google3d' ? 'ellipsoid' : name === 'terrain' ? 'geoid' : 'flat'; }

  async setBasemap(name) {
    const inst = this.instance, cfg = this.cfg;
    if (name === 'google3d' && !(cfg.ion_token || cfg.google_key)) name = 'osm';
    if (this.tiles) { inst.remove(this.tiles); this.tiles = null; }
    for (const l of this.layers) this.globe.removeLayer(l, { disposeLayer: true });
    this.layers = [];
    this.globe.visible = name !== 'google3d';
    if (name === 'google3d') {
      const tiles = new Tiles3D({ enableFetchPlugin: false });
      tiles.tiles.registerPlugin(cfg.google_key
        ? new GoogleCloudAuthPlugin({ apiToken: cfg.google_key, autoRefreshToken: true })
        : new CesiumIonAuthPlugin({ apiToken: cfg.ion_token, assetId: GOOGLE_ION_ASSET, autoRefreshToken: true }));
      await inst.add(tiles);
      this.tiles = tiles;
    } else if (name === 'terrain') {
      await this._addLayer(new ElevationLayer({
        // full resolution: coarser terrain smooths river channels away, and
        // the depth test then hides the water in them
        name: 'terrarium', resolutionFactor: 1,
        source: new TiledImageSource({
          retries: 1, format: new TerrariumFormat(),
          source: new XYZ({ url: TERRARIUM_URL, projection: 'EPSG:3857', crossOrigin: 'anonymous',
                            maxZoom: TERRARIUM_MAX_ZOOM }),
        }),
      }));
      await this._addLayer(this._osm());
    } else if (name === 'osm') {
      await this._addLayer(this._osm());
    }
    this.basemap = name;
    inst.notifyChange();
    return name;
  }

  async _addLayer(layer) {
    await this.globe.addLayer(layer);
    this.layers.push(layer);
  }

  _osm() {
    return new ColorLayer({ name: 'osm', source: new TiledImageSource({ source: new OSM() }) });
  }

  credits() {
    const parts = [];
    if (this.tiles) {
      // Google's terms (Map Tiles API policies): the Google Maps logo — the official
      // outlined one, on a map — with the data providers of the displayed tiles
      // beside it, kept apart from other credits.
      const data = new Set();
      try {
        for (const a of this.tiles.tiles.getAttributions([])) {
          if (a.type === 'string') String(a.value).split(/[;,]/).map(s => s.trim()).filter(Boolean).forEach(s => data.add(s));
        }
      } catch { /* not loaded yet */ }
      data.delete('Google');
      parts.push(`<img class="gmaps" src="${GOOGLE_LOGO}" alt="Google Maps" aria-label="Google Maps">` +
                 (data.size ? ` ${escapeHtml(Array.from(data).sort().join(', '))}` : ''));
      if (!this.cfg.google_key) parts.push('Data via <a href="https://cesium.com/ion" target="_blank" rel="noopener">Cesium ion</a>');
    }
    if (this.basemap === 'osm' || this.basemap === 'terrain') {
      parts.push('© <a href="https://www.openstreetmap.org/copyright" target="_blank" rel="noopener">OpenStreetMap</a> contributors');
    }
    if (this.basemap === 'terrain') parts.push('Terrain: Mapzen / Tilezen (AWS Open Data)');
    return parts.join(' · ');
  }

  createSurface(mesh, opts) { return SurfaceMesh.create(this.instance, mesh, opts); }

  /** Local east, north, up unit vectors at (lon, lat). */
  _enu(lon, lat) {
    const up = new THREE.Vector3(Math.cos(lat * DEG) * Math.cos(lon * DEG), Math.cos(lat * DEG) * Math.sin(lon * DEG),
                                 Math.sin(lat * DEG));
    const east = new THREE.Vector3(-Math.sin(lon * DEG), Math.cos(lon * DEG), 0);
    const north = new THREE.Vector3().crossVectors(up, east);
    return { east, north, up };
  }

  _place(pos, target, up) {
    const cam = this.instance.view.camera;
    cam.position.copy(pos);
    cam.up.copy(up);
    cam.lookAt(target);
    cam.updateMatrixWorld(true);
    this.controls.resetState?.();
    this.instance.notifyChange(cam);
  }

  /**
   * The sun of the sky and the atmosphere: over (lon, lat), 45° up in the
   * south-west (north-west in the southern hemisphere) — the light of
   * hill-shaded maps, not the sun of a date (that is on the wish list).
   */
  setSun(lon, lat) {
    const { east, north, up } = this._enu(lon, lat);
    const toward = north.clone().multiplyScalar(lat >= 0 ? -1 : 1);       // the equator side
    const dir = up.clone().multiplyScalar(Math.SQRT1_2)
      .addScaledVector(toward, 0.5).addScaledVector(east, -0.5).normalize();
    const sun = ELL.toCartesian(lat, lon, 0).addScaledVector(dir, 1.5e11);
    this.atmosphere?.setSunPosition(sun);
    this.instance.notifyChange();
  }

  flyToBounds([w, s, e, n], h) {
    const lon = (w + e) / 2, lat = (s + n) / 2;
    this.setSun(lon, lat);
    const dx = (e - w) * DEG * 6371000 * Math.cos(lat * DEG), dy = (n - s) * DEG * 6371000;
    const r = Math.hypot(dx, dy) / 2, dist = Math.max(r * 2.2, 500), pitch = 40 * DEG;
    const target = ELL.toCartesian(lat, lon, h);
    const { north, up } = this._enu(lon, lat);
    const pos = target.clone().addScaledVector(north, -dist * Math.cos(pitch)).addScaledVector(up, dist * Math.sin(pitch));
    this._place(pos, target, up);
  }

  flyTo(lon, lat, height, headingDeg = 0, pitchDeg = -35) {
    const pos = ELL.toCartesian(lat, lon, height);
    const { east, north, up } = this._enu(lon, lat);
    const hd = headingDeg * DEG, p = pitchDeg * DEG;
    const dir = east.clone().multiplyScalar(Math.sin(hd)).addScaledVector(north, Math.cos(hd))
      .multiplyScalar(Math.cos(p)).addScaledVector(up, Math.sin(p));
    this._place(pos, pos.clone().addScaledVector(dir, 1000), up);
  }

  // (not ELL.toGeodetic: wrong by ~130 m in height in Giro3D 2.0.7 — geometry.js)
  _toLonLat(p) { return ecefToGeodetic(p.x, p.y, p.z); }

  pickLonLat(x, y) {
    const where = [this.tiles || this.globe];
    const hits = this.instance.pickObjectsAt(new THREE.Vector2(x, y), { where, limit: 1 });
    if (hits.length && hits[0].point) return this._toLonLat(hits[0].point);
    // nothing under the cursor: intersect the ellipsoid
    const cam = this.instance.view.camera, rect = this.canvas.getBoundingClientRect();
    const ndc = new THREE.Vector2((x / rect.width) * 2 - 1, -(y / rect.height) * 2 + 1);
    const ray = new THREE.Raycaster(); ray.setFromCamera(ndc, cam);
    const q = rayEllipsoid(ray.ray);
    return q ? this._toLonLat(q) : null;
  }

  /** Canvas pixels of (lon, lat, h), or null when behind the globe or the camera. */
  project(lon, lat, h) {
    const cam = this.instance.view.camera, rect = this.canvas.getBoundingClientRect();
    const p = ELL.toCartesian(lat, lon, h);
    // below the horizon (a sphere through the point): hidden by the globe
    if (cam.position.clone().sub(p).dot(p) < 0) return null;
    const v = p.clone().project(cam);
    if (!(v.z > -1 && v.z < 1)) return null;
    return { x: (v.x + 1) / 2 * rect.width, y: (1 - v.y) / 2 * rect.height };
  }

  /**
   * Basemap heights at points. Terrain: read from the Terrarium tiles. 3D
   * tiles: vertical rays onto the loaded tiles — so first look straight down
   * over the points and wait until the tiles for that view have loaded (coarse
   * tiles sag below the surface: a 70 km tile by about 90 m), then put the
   * camera back.
   */
  async sampleBasemapHeights(points) {
    if (this.basemap === 'terrain') return sampleTerrarium(points);
    const out = new Float64Array(points.length).fill(NaN);
    if (!this.tiles || !points.length) return out;
    const cam = this.instance.view.camera;
    const saved = { pos: cam.position.clone(), quat: cam.quaternion.clone(), up: cam.up.clone() };
    let w = Infinity, e = -Infinity, s = Infinity, n = -Infinity;
    for (const p of points) { w = Math.min(w, p.lon); e = Math.max(e, p.lon); s = Math.min(s, p.lat); n = Math.max(n, p.lat); }
    const lon = (w + e) / 2, lat = (s + n) / 2;
    const span = Math.max((e - w) * DEG * 6371000 * Math.cos(lat * DEG), (n - s) * DEG * 6371000);
    const { north, up } = this._enu(lon, lat);
    const target = ELL.toCartesian(lat, lon, 0);
    this._place(target.clone().addScaledVector(up, Math.max(span * 1.1, 1500)), target, north);
    try {
      await this._waitTiles(60000);
      const group = this.tiles.tiles.group, rc = new THREE.Raycaster();
      rc.firstHitOnly = true;
      points.forEach((p, k) => {
        const { up: u } = this._enu(p.lon, p.lat);
        rc.set(ELL.toCartesian(p.lat, p.lon, 3000), u.clone().negate());
        rc.far = 6000;
        const hit = rc.intersectObject(group, true)[0];
        if (hit) out[k] = this._toLonLat(hit.point).h;
      });
    } finally {
      cam.position.copy(saved.pos); cam.quaternion.copy(saved.quat); cam.up.copy(saved.up);
      cam.updateMatrixWorld(true);
      this.controls.resetState?.();
      this.instance.notifyChange(cam);
    }
    return out;
  }

  /** Resolve when the basemap has stopped loading (at least 1.5 s, at most ms). */
  async _waitTiles(ms) {
    const t0 = performance.now();
    let quiet = 0;
    while (performance.now() - t0 < ms) {
      await new Promise(r => setTimeout(r, 250));
      this.instance.notifyChange(this.tiles || this.globe);
      quiet = this.tilesLoaded() ? quiet + 250 : 0;
      if (quiet >= 1500) return true;
    }
    return false;
  }

  tilesLoaded() {
    const e = this.tiles || this.globe;
    return e ? !e.loading : true;
  }
}

/** First intersection of a ray with the WGS84 ellipsoid (or null). */
function rayEllipsoid(ray) {
  const a = 6378137.0, b = 6356752.314245;
  const o = ray.origin, d = ray.direction;
  const ox = o.x / a, oy = o.y / a, oz = o.z / b, dx = d.x / a, dy = d.y / a, dz = d.z / b;
  const A = dx * dx + dy * dy + dz * dz, B = 2 * (ox * dx + oy * dy + oz * dz), C = ox * ox + oy * oy + oz * oz - 1;
  const disc = B * B - 4 * A * C;
  if (disc < 0) return null;
  const t = (-B - Math.sqrt(disc)) / (2 * A);
  return t > 0 ? o.clone().addScaledVector(d, t) : null;
}

function escapeHtml(s) {
  return s.replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

/**
 * A sky for views near the ground (0.8.15): pale blue at the horizon, deeper
 * blue overhead — a sphere around the camera, drawn first, the elevation of
 * each direction from the local vertical (uniform uUp, set every frame).
 * Giro3D's physical SkyDome renders almost white at the low angles an
 * oblique map view shows.
 */
function skyGradient() {
  const mesh = new THREE.Mesh(new THREE.SphereGeometry(1e5, 32, 16), new THREE.ShaderMaterial({
    side: THREE.BackSide, depthTest: false, depthWrite: false,      // opaque: drawn first (renderOrder)
    uniforms: { uUp: { value: new THREE.Vector3(0, 0, 1) }, uAlpha: { value: 1 } },
    vertexShader: `varying vec3 vDir;
void main() {
  vDir = position;
  gl_Position = projectionMatrix * viewMatrix * vec4(cameraPosition + position, 1.0);
}`,
    fragmentShader: `uniform vec3 uUp;
uniform float uAlpha;
varying vec3 vDir;
void main() {
  float e = dot(normalize(vDir), uUp);                    // sine of the elevation
  vec3 horizon = vec3(0.74, 0.84, 0.93), zenith = vec3(0.20, 0.42, 0.74);
  vec3 c = mix(horizon, zenith, pow(clamp(e, 0.0, 1.0), 0.5));
  gl_FragColor = vec4(mix(vec3(0.043, 0.067, 0.086), c, uAlpha), 1.0);   // (toward the background in space)
}`,
  }));
  mesh.name = 'flooda5-sky';
  mesh.frustumCulled = false;
  mesh.renderOrder = -10000;
  return mesh;
}
