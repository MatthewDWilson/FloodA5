// water.js — the water and bed surfaces as CesiumJS models with a custom shader.
//
// One glTF model holds every cell (its own vertices, so each cell is flat at
// its own level). Per-cell values live in float textures: the frame's depth
// (and the colouring variable) four cells per texel, the static bed and depth
// datum one cell per texel. The vertex shader moves each vertex along the
// ellipsoid normal to
//     water surface:  level = datum + depth            (mode 0)
//     exaggerated:    level = bed + depth × exag       (mode 1)
//     bed:            level = bed                      (mode 2)
// plus the vertical offset (model heights → ellipsoidal heights). A new
// frame is one texture upload: CesiumJS keeps the old texture until the new
// one is created, so nothing flashes. Water is translucent and depth-tested,
// so 3D Tiles buildings stand out of it.
//
// Requires a WebGL 2 context (texelFetch): CesiumJS's default.

import { buildSurface, buildGlb, packCellValues, packStatic, finiteMedian } from './geometry.js';
import { colormapTable } from './colormaps.js';

const VERTEX = `
float fa5Packed(sampler2D tex, float w, int i) {
  int wi = int(w);
  int ti = i / 4;
  int c = i - ti * 4;
  vec4 v = texelFetch(tex, ivec2(ti - (ti / wi) * wi, ti / wi), 0);
  return c == 0 ? v.x : (c == 1 ? v.y : (c == 2 ? v.z : v.w));
}
void vertexMain(VertexInput vsInput, inout czm_modelVertexOutput vsOutput) {
  int i = int(vsInput.attributes.cell + 0.5);
  int sw = int(u_staticW);
  vec4 st = texelFetch(u_static, ivec2(i - (i / sw) * sw, i / sw), 0);
  float d = max(fa5Packed(u_depth, u_texW, i), 0.0);
  float level;
  if (u_mode > 1.5) level = st.x;
  else if (u_mode > 0.5) level = st.x + d * u_exag;
  else level = st.y + d;
  vec3 p = vsInput.attributes.positionMC + u_centre;
  vec3 n = normalize(vec3(p.x / 40680631590769.0, p.y / 40680631590769.0, p.z / 40408299984661.445));
  vsOutput.positionMC = vsInput.attributes.positionMC + n * (level + u_offset - u_h0);
  v_depth = d;
  float cv;
  if (u_cmode < 0.5) cv = d;
  else if (u_cmode < 1.5) cv = st.y + d;
  else if (u_cmode < 2.5) cv = fa5Packed(u_cvar, u_texW, i);
  else cv = st.x;
  v_cval = cv;
}
`;

const FRAGMENT = `
void fragmentMain(FragmentInput fsInput, inout czm_modelMaterial material) {
  if (u_mode < 1.5 && v_depth < u_minDepth) discard;
  float t = clamp((v_cval - u_cmin) / max(u_cmax - u_cmin, 1e-9), 0.0, 1.0);
  material.diffuse = texture(u_cmap, vec2(t * 0.99609375 + 0.001953125, 0.5)).rgb;
  material.alpha = u_alpha;
}
`;

/** Colour modes: what the surface is coloured by. */
export const COLOR_MODE = { depth: 0, level: 1, texture: 2, bed: 3 };

export class SurfaceLayer {
  /**
   * @param {Cesium.Viewer} viewer
   * @param {object} mesh       parsed mesh (protocol.js)
   * @param {object} opts       {role: 'water'|'bed', h0}
   */
  static async create(viewer, mesh, opts = {}) {
    const layer = new SurfaceLayer(viewer, mesh, opts);
    await layer._load();
    return layer;
  }

  constructor(viewer, mesh, { role = 'water', h0 } = {}) {
    this.viewer = viewer;
    this.mesh = mesh;
    this.role = role;
    this.h0 = Number.isFinite(h0) ? h0 : (finiteMedian(mesh.base) || 0);
    this.n = mesh.nCells;
    this.style = { mode: role === 'bed' ? 2 : 0, cmode: role === 'bed' ? COLOR_MODE.bed : COLOR_MODE.depth,
                   cmin: 0, cmax: 1, alpha: role === 'bed' ? 1 : 0.8, minDepth: 0.05, exag: 10,
                   offset: 0, colormap: role === 'bed' ? 'terrain' : 'depth' };
    this._depthBuf = null;
    this._cvarBuf = null;
  }

  async _load() {
    const C = window.Cesium;
    const surface = buildSurface(this.mesh, this.h0);
    this.centre = surface.centre;
    const glb = buildGlb(surface);
    const url = URL.createObjectURL(new Blob([glb], { type: 'model/gltf-binary' }));
    const st = packStatic(this.mesh.elevation, this.mesh.base, this.n);
    const zero = packCellValues(new Float32Array(0), this.n);
    this._texW = zero.width;
    const U = C.UniformType;
    this.shader = new C.CustomShader({
      uniforms: {
        u_depth: { type: U.SAMPLER_2D, value: floatTexture(zero) },
        u_cvar: { type: U.SAMPLER_2D, value: floatTexture(zero) },
        u_static: { type: U.SAMPLER_2D, value: floatTexture(st) },
        u_cmap: { type: U.SAMPLER_2D, value: colormapTexture(this.style.colormap) },
        u_texW: { type: U.FLOAT, value: zero.width },
        u_staticW: { type: U.FLOAT, value: st.width },
        u_centre: { type: U.VEC3, value: new C.Cartesian3(...surface.centre) },
        u_h0: { type: U.FLOAT, value: this.h0 },
        u_offset: { type: U.FLOAT, value: this.style.offset },
        u_mode: { type: U.FLOAT, value: this.style.mode },
        u_cmode: { type: U.FLOAT, value: this.style.cmode },
        u_cmin: { type: U.FLOAT, value: this.style.cmin },
        u_cmax: { type: U.FLOAT, value: this.style.cmax },
        u_alpha: { type: U.FLOAT, value: this.style.alpha },
        u_minDepth: { type: U.FLOAT, value: this.style.minDepth },
        u_exag: { type: U.FLOAT, value: this.style.exag },
      },
      varyings: { v_depth: C.VaryingType.FLOAT, v_cval: C.VaryingType.FLOAT },
      lightingModel: C.LightingModel.UNLIT,
      translucencyMode: this.role === 'bed' ? C.CustomShaderTranslucencyMode.OPAQUE
                                            : C.CustomShaderTranslucencyMode.TRANSLUCENT,
      vertexShaderText: VERTEX,
      fragmentShaderText: FRAGMENT,
    });
    try {
      this.model = await C.Model.fromGltfAsync({
        url,
        modelMatrix: C.Matrix4.fromTranslation(new C.Cartesian3(...surface.centre)),
        upAxis: C.Axis.Z,
        forwardAxis: C.Axis.X,
        customShader: this.shader,
        shadows: C.ShadowMode.DISABLED,
        backFaceCulling: false,
        allowPicking: false,
      });
    } finally {
      // the model keeps its own copy once loaded
      setTimeout(() => URL.revokeObjectURL(url), 60000);
    }
    this.viewer.scene.primitives.add(this.model);
    await new Promise((resolve, reject) => {
      const t0 = performance.now();
      const check = () => {
        if (this.model.ready) return resolve();
        if (performance.now() - t0 > 120000) return reject(new Error('surface model did not load'));
        requestAnimationFrame(check);
      };
      check();
    });
  }

  /** Frame depths (Float32Array, mesh order). */
  setDepth(depth) {
    const p = packCellValues(depth, this.n, this._depthBuf);
    this._depthBuf = p.data;
    this.shader.setUniform('u_depth', floatTexture(p));
  }

  /** Values for COLOR_MODE.texture (Float32Array, mesh order). */
  setColorValues(values) {
    const p = packCellValues(values, this.n, this._cvarBuf);
    this._cvarBuf = p.data;
    this.shader.setUniform('u_cvar', floatTexture(p));
  }

  /** Update any of: mode, cmode, cmin, cmax, alpha, minDepth, exag, offset, colormap. */
  setStyle(s) {
    const old = this.style;
    this.style = { ...old, ...s };
    const set = (name, key) => {
      if (s[key] !== undefined && s[key] !== old[key]) this.shader.setUniform(name, this.style[key]);
    };
    set('u_mode', 'mode'); set('u_cmode', 'cmode'); set('u_cmin', 'cmin'); set('u_cmax', 'cmax');
    set('u_alpha', 'alpha'); set('u_minDepth', 'minDepth'); set('u_exag', 'exag'); set('u_offset', 'offset');
    if (s.colormap && s.colormap !== old.colormap) this.shader.setUniform('u_cmap', colormapTexture(s.colormap));
  }

  set show(v) { if (this.model) this.model.show = v; }
  get show() { return !!(this.model && this.model.show); }

  destroy() {
    if (this.model) {
      this.viewer.scene.primitives.remove(this.model);   // destroys it
      this.model = null;
    }
  }
}

function floatTexture({ data, width, height }) {
  const C = window.Cesium;
  return new C.TextureUniform({
    typedArray: data, width, height,
    pixelFormat: C.PixelFormat.RGBA,
    pixelDatatype: C.PixelDatatype.FLOAT,
    minificationFilter: C.TextureMinificationFilter.NEAREST,
    magnificationFilter: C.TextureMagnificationFilter.NEAREST,
  });
}

const cmapCache = new Map();
function colormapTexture(name) {
  const C = window.Cesium;
  if (!cmapCache.has(name)) cmapCache.set(name, colormapTable(name));
  return new C.TextureUniform({
    typedArray: cmapCache.get(name), width: 256, height: 1,
    pixelFormat: C.PixelFormat.RGBA,
    pixelDatatype: C.PixelDatatype.UNSIGNED_BYTE,
    minificationFilter: C.TextureMinificationFilter.NEAREST,
    magnificationFilter: C.TextureMagnificationFilter.NEAREST,
  });
}
