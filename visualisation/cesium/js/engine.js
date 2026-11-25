// engine.js — the CesiumJS engine for the globe viewer (the interface:
// viewer-common/js/app.js). CesiumJS is loaded at run time from the URL in
// /api/config (default: the pinned 1.147 release on cesium.com).

import { SurfaceLayer } from './water.js';

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

export class CesiumEngine {
  constructor() {
    this.name = 'CesiumJS';
    this.viewer = null;
    this.tileset = null;
    this.cfg = null;
    this.webgl2 = null;
  }

  async init(container, cfg, log) {
    this.cfg = cfg; this.log = log;
    const C = await loadCesium(cfg.cesium_url);
    this.name = `CesiumJS ${C.VERSION}`;
    C.Ion.defaultAccessToken = cfg.ion_token || '';
    const viewer = new C.Viewer(container, {
      baseLayer: false,
      terrainProvider: new C.EllipsoidTerrainProvider(),
      baseLayerPicker: false, geocoder: false, homeButton: false, sceneModePicker: false,
      navigationHelpButton: false, animation: false, timeline: false, fullscreenButton: false,
      infoBox: false, selectionIndicator: false,
      contextOptions: { requestWebgl1: false },
    });
    const scene = viewer.scene;
    this.webgl2 = !!scene.context.webgl2;
    if (!this.webgl2) throw new Error('This viewer needs WebGL 2.');
    scene.globe.baseColor = C.Color.fromCssColorString('#1b2731');
    scene.pickTranslucentDepth = true;
    scene.screenSpaceCameraController.enableCollisionDetection = true;
    scene.debugShowFramesPerSecond = new URLSearchParams(location.search).has('fps');
    this.viewer = viewer;
    this.canvas = scene.canvas;
  }

  _google() { return !!(this.cfg.ion_token || this.cfg.google_key); }

  basemaps(cfg) {
    const ion = !!cfg.ion_token, google = !!(cfg.ion_token || cfg.google_key);
    return [
      { value: 'google3d', label: 'Google Photorealistic 3D', available: google },
      { value: 'terrain', label: 'Cesium World Terrain + OSM', available: ion },
      { value: 'osm', label: 'OpenStreetMap (flat)', available: true },
      { value: 'none', label: 'None', available: true },
    ];
  }

  defaultBasemap(cfg) { return cfg.ion_token || cfg.google_key ? 'google3d' : 'osm'; }

  photoreal(name) { return name === 'google3d'; }

  /** The height reference of a basemap: 'ellipsoid' (heights above WGS84) or 'flat'. */
  heightDatum(name) { return name === 'google3d' || name === 'terrain' ? 'ellipsoid' : 'flat'; }

  /**
   * Basemap layers are created once and then shown or hidden, never removed:
   * destroying an imagery layer while CesiumJS still has work queued for it
   * (reprojection of its tiles) fails inside the render loop (seen switching
   * from OpenStreetMap to the Google tiles: "e._va is undefined").
   */
  async setBasemap(name) {
    const C = window.Cesium, viewer = this.viewer, scene = viewer.scene;
    if ((name === 'google3d' && !this._google()) || (name === 'terrain' && !this.cfg.ion_token)) name = 'osm';
    // load what the new basemap needs before changing what is shown
    let tileset = null, terrain = null;
    if (name === 'google3d') {
      this._google3d ??= C.createGooglePhotorealistic3DTileset(
        this.cfg.google_key ? { key: this.cfg.google_key, onlyUsingWithGoogleGeocoder: true }
                            : { onlyUsingWithGoogleGeocoder: true },
        { enableCollision: true }).then(t => { t.show = false; scene.primitives.add(t); return t; })
        .catch(e => { this._google3d = null; throw e; });
      tileset = await this._google3d;
    } else if (name === 'terrain') {
      this._worldTerrain ??= C.createWorldTerrainAsync().catch(e => { this._worldTerrain = null; throw e; });
      terrain = await this._worldTerrain;
    }
    this._ellipsoid ??= new C.EllipsoidTerrainProvider();
    this._osmLayer ??= viewer.imageryLayers.addImageryProvider(
      new C.OpenStreetMapImageryProvider({ url: 'https://tile.openstreetmap.org/' }));
    const wantTerrain = terrain || this._ellipsoid;
    if (viewer.terrainProvider !== wantTerrain) viewer.terrainProvider = wantTerrain;
    this._osmLayer.show = name === 'osm' || name === 'terrain';
    scene.globe.show = name !== 'google3d';
    scene.globe.depthTestAgainstTerrain = name === 'terrain';
    this._google3d?.then(t => { t.show = name === 'google3d'; }).catch(() => {});
    this.tileset = tileset;
    return name;
  }

  credits() { return ''; }      // CesiumJS shows its own credits

  createSurface(mesh, opts) { return SurfaceLayer.create(this.viewer, mesh, opts); }

  flyToBounds([w, s, e, n], h, duration = 1.5) {
    const C = window.Cesium;
    const pts = [[w, s], [e, s], [e, n], [w, n], [(w + e) / 2, (s + n) / 2]].map(([x, y]) => C.Cartesian3.fromDegrees(x, y, h));
    const sphere = C.BoundingSphere.fromPoints(pts);
    this.viewer.camera.flyToBoundingSphere(sphere, {
      offset: new C.HeadingPitchRange(0, C.Math.toRadians(-40), Math.max(sphere.radius * 2.2, 500)),
      duration,
    });
  }

  flyTo(lon, lat, height, headingDeg = 0, pitchDeg = -35) {
    const C = window.Cesium;
    this.viewer.camera.flyTo({
      destination: C.Cartesian3.fromDegrees(lon, lat, height),
      orientation: { heading: C.Math.toRadians(headingDeg), pitch: C.Math.toRadians(pitchDeg), roll: 0 },
      duration: 0,
    });
  }

  pickLonLat(x, y) {
    const C = window.Cesium, scene = this.viewer.scene;
    const pos = new C.Cartesian2(x, y);
    let p = scene.pickPositionSupported ? scene.pickPosition(pos) : undefined;
    if (!C.defined(p)) p = this.viewer.camera.pickEllipsoid(pos);
    if (!C.defined(p)) return null;
    const c = C.Cartographic.fromCartesian(p);
    return { lon: C.Math.toDegrees(c.longitude), lat: C.Math.toDegrees(c.latitude) };
  }

  async sampleBasemapHeights(points, surfaces) {
    const C = window.Cesium;
    const cartos = points.map(p => C.Cartographic.fromDegrees(p.lon, p.lat));
    const res = this.tileset
      ? await this.viewer.scene.sampleHeightMostDetailed(cartos, surfaces.map(s => s.model).filter(Boolean))
      : await C.sampleTerrainMostDetailed(this.viewer.terrainProvider, cartos);
    return Float64Array.from(res, c => (c && Number.isFinite(c.height) ? c.height : NaN));
  }

  tilesLoaded() {
    return this.tileset ? this.tileset.tilesLoaded : this.viewer.scene.globe.tilesLoaded;
  }
}
