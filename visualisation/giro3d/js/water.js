// water.js — the water and bed surfaces as three.js meshes for Giro3D
// (the model and its uniforms: viewer-common/js/surface.js).
//
// One mesh per surface, positions relative to an ECEF centre (Object3D
// position, float64 on the CPU: three forms the model-view matrix in double
// precision, so street-level views do not jitter). A new frame rewrites one
// float DataTexture. The shader includes three's logarithmic-depth chunks so
// that depth testing against Giro3D's tiles is right whichever depth buffer
// the renderer uses. Needs WebGL 2 (texelFetch).

import * as THREE from 'three';
import { buildSurface, packCellValues, packStatic, finiteMedian } from '../../viewer-common/js/geometry.js';
import { colormapTable } from '../../viewer-common/js/colormaps.js';
import { defaultStyle, glslHelpers, glslFragmentHelpers } from '../../viewer-common/js/surface.js';

const VERTEX = `
#include <common>
#include <logdepthbuf_pars_vertex>
${glslHelpers()}
attribute float cell;
attribute float edge;
attribute float corner;          // smooth surface: the vertex's corner (−1 at the centre; unset until smoothing is first on)
uniform sampler2D uDepth;
uniform sampler2D uCvar;
uniform sampler2D uStatic;
uniform sampler2D uCorner;
uniform float uTexW, uStaticW, uCornerW, uH0, uOffset, uMode, uCmode, uExag, uDepthBias, uBiasNear, uBiasFar, uSmooth, uBlend;
uniform vec3 uCentre;
varying float vDepth;
varying float vCval;
varying float vEdge;
varying vec3 vPos;
varying vec3 vUp;
void main() {
  int i = int(cell + 0.5);
  vec4 st = fa5Static(uStatic, uStaticW, i);
  float d = max(fa5Packed(uDepth, uTexW, i), 0.0);
  float level = fa5Level(uMode, uExag, d, st);
  float cv = (uCmode > 1.5 && uCmode < 2.5) ? fa5Packed(uCvar, uTexW, i) : 0.0;
  float cval = fa5ColorValue(uCmode, d, st, cv);
  if (uSmooth > 0.5 && corner > -0.5) {
    vec4 cc = fa5Corner(uCorner, uCornerW, int(corner + 0.5));
    level = fa5LevelC(uMode, uExag, cc);
    if (uBlend > 0.5) cval = fa5ColorValueC(uCmode, cc);
  }
  vec3 n = fa5Normal(position + uCentre);
  // (the mesh's matrix is a translation to uCentre: model = world − uCentre)
  vec3 p0 = position + n * (level + uOffset - uH0);
  vec3 p = fa5DepthBias(p0, cameraPosition - uCentre, uDepthBias, vec2(uBiasNear, uBiasFar));
  vDepth = d;
  vEdge = edge;
  vCval = cval;
  vPos = p0; vUp = n;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(p, 1.0);
  #include <logdepthbuf_vertex>
}
`;

const FRAGMENT = `
#include <common>
#include <logdepthbuf_pars_fragment>
${glslHelpers()}${glslFragmentHelpers()}
uniform sampler2D uCmap;
uniform float uCmin, uCmax, uAlpha, uMinDepth, uMode, uBorders, uSmooth, uShadeExag;
varying float vDepth;
varying float vCval;
varying float vEdge;
varying vec3 vPos;
varying vec3 vUp;
// three.js declares no fragment output for an explicit GLSL3 ShaderMaterial
layout(location = 0) out highp vec4 fa5FragColor;
void main() {
  float line = fa5Border(vEdge, uBorders);
  if (uMode < 1.5 && vDepth < uMinDepth) discard;
  #include <logdepthbuf_fragment>
  vec3 c = texture(uCmap, fa5MapCoord(vCval, uCmin, uCmax)).rgb;
  if (uSmooth > 0.5) c = clamp(c * fa5Shade(vPos, normalize(vUp), uShadeExag), 0.0, 1.0);
  c = mix(c, vec3(1.0), line);
  fa5FragColor = vec4(c, max(uAlpha, line));
}
`;

const STYLE_UNIFORMS = { mode: 'uMode', cmode: 'uCmode', cmin: 'uCmin', cmax: 'uCmax', alpha: 'uAlpha',
                         minDepth: 'uMinDepth', exag: 'uExag', offset: 'uOffset', borders: 'uBorders',
                         depthBias: 'uDepthBias', biasNear: 'uBiasNear', biasFar: 'uBiasFar', smooth: 'uSmooth', blend: 'uBlend', shadeExag: 'uShadeExag' };

export class SurfaceMesh {
  /**
   * @param {object} instance  Giro3D Instance
   * @param {object} mesh      parsed mesh (protocol.js)
   * @param {object} opts      {role: 'water'|'bed', h0}
   */
  static async create(instance, mesh, opts = {}) {
    const s = new SurfaceMesh(instance, mesh, opts);
    await instance.add(s.object3d);
    instance.notifyChange();
    return s;
  }

  constructor(instance, mesh, { role = 'water', h0 } = {}) {
    this.instance = instance;
    this.role = role;
    this.n = mesh.nCells;
    this.h0 = Number.isFinite(h0) ? h0 : (finiteMedian(mesh.base) || 0);
    this.style = defaultStyle(role);
    const surface = buildSurface(mesh, this.h0);
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(surface.positions, 3));
    g.setAttribute('cell', new THREE.BufferAttribute(surface.cell, 1));
    g.setAttribute('edge', new THREE.BufferAttribute(surface.edge, 1));
    g.setIndex(new THREE.BufferAttribute(surface.indices, 1));
    g.computeBoundingSphere();
    this.depthTex = floatTexture(packCellValues(new Float32Array(0), this.n));
    this.cvarTex = floatTexture(packCellValues(new Float32Array(0), this.n));
    const st = packStatic(mesh.elevation, mesh.base, this.n);
    const s = this.style;
    this.material = new THREE.ShaderMaterial({
      glslVersion: THREE.GLSL3,
      vertexShader: VERTEX,
      fragmentShader: FRAGMENT,
      uniforms: {
        uDepth: { value: this.depthTex },
        uCvar: { value: this.cvarTex },
        uStatic: { value: floatTexture(st) },
        uCmap: { value: colormapTexture(s.colormap) },
        uTexW: { value: this.depthTex.image.width },
        uStaticW: { value: st.width },
        uCentre: { value: new THREE.Vector3(...surface.centre) },
        uH0: { value: this.h0 },
        uOffset: { value: s.offset }, uMode: { value: s.mode }, uCmode: { value: s.cmode },
        uCmin: { value: s.cmin }, uCmax: { value: s.cmax }, uAlpha: { value: s.alpha },
        uMinDepth: { value: s.minDepth }, uExag: { value: s.exag },
        uBorders: { value: s.borders }, uDepthBias: { value: s.depthBias },
        uBiasNear: { value: s.biasNear }, uBiasFar: { value: s.biasFar },
        uCorner: { value: floatTexture({ data: new Float32Array(4), width: 1, height: 1 }) }, uCornerW: { value: 1 },
        uSmooth: { value: s.smooth }, uBlend: { value: s.blend }, uShadeExag: { value: s.shadeExag },
      },
      transparent: role !== 'bed' || s.alpha < 1,     // the bed: opaque unless its opacity is < 1
      depthWrite: role === 'bed' && s.alpha >= 1,
      side: THREE.DoubleSide,
    });
    this.object3d = new THREE.Mesh(g, this.material);
    this.object3d.name = `flooda5-${role}`;
    this.object3d.position.set(...surface.centre);
    this.object3d.frustumCulled = false;            // vertices move in the shader
    // Not raycastable: the globe controls raycast the scene every frame to find
    // the ground, and a brute-force pass over millions of vertices (whose
    // CPU positions are not where the shader draws them anyway) starves tile
    // loading. Picking uses the basemap (engine.pickLonLat).
    this.object3d.raycast = () => {};
    this.object3d.renderOrder = role === 'bed' ? 1 : 2;
    this.object3d.updateMatrixWorld(true);
    // Giro3D turns three.js's shader error checks off: check this program once
    // it has been compiled, so a broken shader is reported, not silently blank.
    this.object3d.onAfterRender = renderer => {
      this.object3d.onAfterRender = () => {};
      checkProgram(renderer, this.material, role);
    };
  }

  setDepth(depth) {
    packCellValues(depth, this.n, this.depthTex.image.data);
    this.depthTex.needsUpdate = true;
    this.instance.notifyChange();
  }

  /** Smooth surface: the corners (geometry.js cornerIndex), attached once. */
  async setCorners(ci) {
    const g = this.object3d.geometry;
    if (ci && !g.getAttribute('corner')) g.setAttribute('corner', new THREE.BufferAttribute(ci.vertexCorner, 1));
  }

  /** Smooth surface: the corner values (geometry.js cornerValues / bedCornerValues, packed). */
  setCornerValues(p) {
    const u = this.material.uniforms.uCorner;
    if (u.value.image.width === p.width && u.value.image.height === p.height) {
      u.value.image.data.set(p.data);
      u.value.needsUpdate = true;
    } else {
      u.value.dispose();
      u.value = floatTexture({ data: Float32Array.from(p.data), width: p.width, height: p.height });
      this.material.uniforms.uCornerW.value = p.width;
    }
    this.instance.notifyChange();
  }

  setColorValues(values) {
    packCellValues(values, this.n, this.cvarTex.image.data);
    this.cvarTex.needsUpdate = true;
    this.instance.notifyChange();
  }

  setStyle(s) {
    const old = this.style;
    this.style = { ...old, ...s };
    for (const [k, u] of Object.entries(STYLE_UNIFORMS)) {
      if (s[k] !== undefined) this.material.uniforms[u].value = this.style[k];
    }
    if (s.colormap && s.colormap !== old.colormap) {
      this.material.uniforms.uCmap.value.dispose();
      this.material.uniforms.uCmap.value = colormapTexture(s.colormap);
    }
    if (this.role === 'bed' && s.alpha !== undefined) {
      const tr = this.style.alpha < 1;
      if (tr !== this.material.transparent) {
        this.material.transparent = tr;
        this.material.depthWrite = !tr;
        this.material.needsUpdate = true;
      }
    }
    this.instance.notifyChange();
  }

  set show(v) { this.object3d.visible = v; this.instance.notifyChange(); }
  get show() { return this.object3d.visible; }

  destroy() {
    this.instance.remove(this.object3d);
    this.object3d.geometry.dispose();
    for (const u of Object.values(this.material.uniforms)) if (u.value && u.value.isTexture) u.value.dispose();
    this.material.dispose();
    this.instance.notifyChange();
  }
}

function checkProgram(renderer, material, role) {
  try {
    const gl = renderer.getContext(), p = renderer.properties.get(material).currentProgram;
    if (!p || gl.getProgramParameter(p.program, gl.LINK_STATUS)) return;
    const msg = `FloodA5 ${role} shader failed: ${gl.getShaderInfoLog(p.vertexShader) || ''} ` +
                `${gl.getShaderInfoLog(p.fragmentShader) || ''} ${gl.getProgramInfoLog(p.program) || ''}`;
    console.error(msg.trim());
    window.__flooda5?.errors.push(msg.trim());
  } catch { /* diagnostics only */ }
}

function floatTexture({ data, width, height }) {
  const t = new THREE.DataTexture(data, width, height, THREE.RGBAFormat, THREE.FloatType);
  t.minFilter = t.magFilter = THREE.NearestFilter;
  t.generateMipmaps = false;
  t.needsUpdate = true;
  return t;
}

function colormapTexture(name) {
  const t = new THREE.DataTexture(colormapTable(name), 256, 1, THREE.RGBAFormat, THREE.UnsignedByteType);
  t.minFilter = t.magFilter = THREE.NearestFilter;
  t.generateMipmaps = false;
  t.needsUpdate = true;
  return t;
}
