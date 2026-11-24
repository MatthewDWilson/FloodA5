// water.js — the water and bed surfaces as CesiumJS models with a custom shader
// (the model and its uniforms: viewer-common/js/surface.js).
//
// One glTF model per surface, built in the browser; a new frame is one texture
// upload (CesiumJS keeps the old texture until the new one is created, so
// nothing flashes). Water is translucent and depth-tested. Needs WebGL 2.

import { buildSurface, buildGlb, packCellValues, packStatic, finiteMedian } from '../../viewer-common/js/geometry.js';
import { colormapTable } from '../../viewer-common/js/colormaps.js';
import { defaultStyle, glslHelpers } from '../../viewer-common/js/surface.js';

const VERTEX = `${glslHelpers()}
void vertexMain(VertexInput vsInput, inout czm_modelVertexOutput vsOutput) {
  int i = int(vsInput.attributes.cell + 0.5);
  vec4 st = fa5Static(u_static, u_staticW, i);
  float d = max(fa5Packed(u_depth, u_texW, i), 0.0);
  float level = fa5Level(u_mode, u_exag, d, st);
  vec3 n = fa5Normal(vsInput.attributes.positionMC + u_centre);
  vsOutput.positionMC = vsInput.attributes.positionMC + n * (level + u_offset - u_h0);
  v_depth = d;
  float cv = (u_cmode > 1.5 && u_cmode < 2.5) ? fa5Packed(u_cvar, u_texW, i) : 0.0;
  v_cval = fa5ColorValue(u_cmode, d, st, cv);
}
`;

// (CesiumJS compiles the custom vertex and fragment texts separately: both need the helpers)
const FRAGMENT = `${glslHelpers()}
void fragmentMain(FragmentInput fsInput, inout czm_modelMaterial material) {
  if (u_mode < 1.5 && v_depth < u_minDepth) discard;
  material.diffuse = texture(u_cmap, fa5MapCoord(v_cval, u_cmin, u_cmax)).rgb;
  material.alpha = u_alpha;
}
`;

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
    this.shader = new C.CustomShader({
      uniforms: {
        u_depth: { type: U.SAMPLER_2D, value: floatTexture(zero) },
        u_cvar: { type: U.SAMPLER_2D, value: floatTexture(zero) },
        u_static: { type: U.SAMPLER_2D, value: floatTexture(st) },
        u_cmap: { type: U.SAMPLER_2D, value: colormapTexture(s.colormap) },
        u_texW: { type: U.FLOAT, value: zero.width },
        u_staticW: { type: U.FLOAT, value: st.width },
        u_centre: { type: U.VEC3, value: new C.Cartesian3(...surface.centre) },
        u_h0: { type: U.FLOAT, value: this.h0 },
        u_offset: { type: U.FLOAT, value: s.offset },
        u_mode: { type: U.FLOAT, value: s.mode },
        u_cmode: { type: U.FLOAT, value: s.cmode },
        u_cmin: { type: U.FLOAT, value: s.cmin },
        u_cmax: { type: U.FLOAT, value: s.cmax },
        u_alpha: { type: U.FLOAT, value: s.alpha },
        u_minDepth: { type: U.FLOAT, value: s.minDepth },
        u_exag: { type: U.FLOAT, value: s.exag },
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
