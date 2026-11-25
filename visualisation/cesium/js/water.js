// water.js — the water and bed surfaces as CesiumJS models with a custom shader
// (the model and its uniforms: viewer-common/js/surface.js).
//
// One glTF model per surface, built in the browser; a new frame is one texture
// upload (CesiumJS keeps the old texture until the new one is created, so
// nothing flashes). Water is translucent and depth-tested; the bed is opaque
// unless its opacity is below 1 (CesiumJS fixes a custom shader's
// translucency when it is created: the shader is rebuilt when that changes).
// Needs WebGL 2.

import { buildSurface, buildGlb, packCellValues, packStatic, finiteMedian } from '../../viewer-common/js/geometry.js';
import { colormapTable } from '../../viewer-common/js/colormaps.js';
import { defaultStyle, glslHelpers, glslFragmentHelpers } from '../../viewer-common/js/surface.js';

const VERTEX = `${glslHelpers()}
void vertexMain(VertexInput vsInput, inout czm_modelVertexOutput vsOutput) {
  int i = int(vsInput.attributes.cell + 0.5);
  vec4 st = fa5Static(u_static, u_staticW, i);
  float d = max(fa5Packed(u_depth, u_texW, i), 0.0);
  float level = fa5Level(u_mode, u_exag, d, st);
  vec3 n = fa5Normal(vsInput.attributes.positionMC + u_centre);
  vec3 p = vsInput.attributes.positionMC + n * (level + u_offset - u_h0);
  vsOutput.positionMC = fa5DepthBias(p, (czm_inverseModelView * vec4(0.0, 0.0, 0.0, 1.0)).xyz, u_depthBias);
  v_depth = d;
  v_fa5Edge = vsInput.attributes.fa5edge;
  float cv = (u_cmode > 1.5 && u_cmode < 2.5) ? fa5Packed(u_cvar, u_texW, i) : 0.0;
  v_cval = fa5ColorValue(u_cmode, d, st, cv);
}
`;

// (CesiumJS compiles the custom vertex and fragment texts separately: both need the helpers)
const FRAGMENT = `${glslHelpers()}${glslFragmentHelpers()}
void fragmentMain(FragmentInput fsInput, inout czm_modelMaterial material) {
  float line = fa5Border(v_fa5Edge, u_borders);
  if (u_mode < 1.5 && v_depth < u_minDepth) discard;
  material.diffuse = mix(texture(u_cmap, fa5MapCoord(v_cval, u_cmin, u_cmax)).rgb, vec3(1.0), line);
  material.alpha = max(u_alpha, line);
}
`;

/** Style keys and their uniforms (surface.js defaultStyle). */
const STYLE_UNIFORMS = [
  ['u_offset', 'offset'], ['u_mode', 'mode'], ['u_cmode', 'cmode'], ['u_cmin', 'cmin'], ['u_cmax', 'cmax'],
  ['u_alpha', 'alpha'], ['u_minDepth', 'minDepth'], ['u_exag', 'exag'], ['u_borders', 'borders'],
  ['u_depthBias', 'depthBias'],
];

export class SurfaceLayer {
  /** @param {object} opts {role: 'water'|'bed', h0} */
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
    this.style = defaultStyle(role);
    this._depthBuf = null;
    this._cvarBuf = null;
  }

  async _load() {
    const C = window.Cesium;
    const surface = buildSurface(this.mesh, this.h0);
    this.centre = surface.centre;
    const url = URL.createObjectURL(new Blob([buildGlb(surface)], { type: 'model/gltf-binary' }));
    const st = packStatic(this.mesh.elevation, this.mesh.base, this.n);
    const zero = packCellValues(new Float32Array(0), this.n);
    const U = C.UniformType, s = this.style;
    // current uniform values (type, value): kept to rebuild the shader
    this._uni = {
      u_depth: [U.SAMPLER_2D, floatTexture(zero)],
      u_cvar: [U.SAMPLER_2D, floatTexture(zero)],
      u_static: [U.SAMPLER_2D, floatTexture(st)],
      u_cmap: [U.SAMPLER_2D, colormapTexture(s.colormap)],
      u_texW: [U.FLOAT, zero.width],
      u_staticW: [U.FLOAT, st.width],
      u_centre: [U.VEC3, new C.Cartesian3(...surface.centre)],
      u_h0: [U.FLOAT, this.h0],
    };
    for (const [name, key] of STYLE_UNIFORMS) this._uni[name] = [U.FLOAT, s[key]];
    this.translucent = this._wantsTranslucent();
    this.shader = this._makeShader();
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
      setTimeout(() => URL.revokeObjectURL(url), 60000);   // the model keeps its own copy
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

  _wantsTranslucent() { return this.role !== 'bed' || this.style.alpha < 1; }

  _makeShader() {
    const C = window.Cesium, uniforms = {};
    for (const [name, [type, value]] of Object.entries(this._uni)) uniforms[name] = { type, value };
    return new C.CustomShader({
      uniforms,
      varyings: { v_depth: C.VaryingType.FLOAT, v_cval: C.VaryingType.FLOAT, v_fa5Edge: C.VaryingType.FLOAT },
      lightingModel: C.LightingModel.UNLIT,
      translucencyMode: this.translucent ? C.CustomShaderTranslucencyMode.TRANSLUCENT
                                         : C.CustomShaderTranslucencyMode.OPAQUE,
      vertexShaderText: VERTEX,
      fragmentShaderText: FRAGMENT,
    });
  }

  _set(name, value) {
    this._uni[name][1] = value;
    this.shader.setUniform(name, value);
  }

  /** Frame depths (Float32Array, mesh order). */
  setDepth(depth) {
    const p = packCellValues(depth, this.n, this._depthBuf);
    this._depthBuf = p.data;
    this._set('u_depth', floatTexture(p));
  }

  /** Values for COLOR_MODE.texture (Float32Array, mesh order). */
  setColorValues(values) {
    const p = packCellValues(values, this.n, this._cvarBuf);
    this._cvarBuf = p.data;
    this._set('u_cvar', floatTexture(p));
  }

  /** Update any of: mode, cmode, cmin, cmax, alpha, minDepth, exag, offset, borders, colormap. */
  setStyle(s) {
    const old = this.style;
    this.style = { ...old, ...s };
    for (const [name, key] of STYLE_UNIFORMS) {
      if (s[key] !== undefined && s[key] !== old[key]) this._set(name, this.style[key]);
    }
    if (s.colormap && s.colormap !== old.colormap) this._set('u_cmap', colormapTexture(s.colormap));
    const tr = this._wantsTranslucent();
    if (tr !== this.translucent && this.model) {
      this.translucent = tr;
      this.shader = this._makeShader();
      this.model.customShader = this.shader;
    }
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
