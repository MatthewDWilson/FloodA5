// surface.js — what the water and bed surfaces share across rendering engines:
// the style model, the colour modes and the GLSL that reads per-cell values.
//
// Both engines draw one mesh per surface, each cell on its own vertices (a fan
// from its centre) with float attributes `cell` (its index) and `edge` (1 at
// the centre, 0 on the boundary: the cell-border lines), and three float
// textures:
//   depth  — the frame's water depth, four cells per RGBA texel
//   cvar   — the colouring variable (speed, …), same packing
//   static — per cell (bed elevation, depth datum, 0, 0), one cell per texel
// The vertex shader moves each vertex along the ellipsoid normal to the level
//   mode 0 (water surface):  datum + depth
//   mode 1 (exaggerated):    bed + depth × exag
//   mode 2 (bed):            bed
// plus the vertical offset; the fragment shader drops dry cells, colours
// through a 256-entry colour map and, if `borders` is on, draws the cell
// borders as white hairlines.

/** What the surface is coloured by (uniform `cmode`). */
export const COLOR_MODE = { depth: 0, level: 1, texture: 2, bed: 3 };

/** Surface geometry (uniform `mode`). */
export const SURFACE_MODE = { water: 0, exaggerated: 1, bed: 2 };

/** Water drawn this fraction of its distance nearer the camera when zoomed out (fa5DepthBias). */
export const DEPTH_BIAS = 0.015;

/** Default style of a surface. */
export function defaultStyle(role) {
  const bed = role === 'bed';
  return {
    mode: bed ? SURFACE_MODE.bed : SURFACE_MODE.water,
    cmode: bed ? COLOR_MODE.bed : COLOR_MODE.depth,
    cmin: 0, cmax: 1, alpha: bed ? 1 : 0.8, minDepth: 0.05, exag: 10, offset: 0, borders: 0,
    depthBias: bed ? 0 : DEPTH_BIAS,
    colormap: bed ? 'terrain' : 'depth',
  };
}

/** WGS84 a² and b², for the ellipsoid normal in the shaders. */
export const GLSL_ELLIPSOID = `
const float FA5_A2 = 40680631590769.0;
const float FA5_B2 = 40408299984661.445;
`;

/**
 * GLSL ES 3.0 helpers (need texelFetch: WebGL 2). Uniform names are passed in
 * so each engine can follow its own naming.
 */
export function glslHelpers() {
  return `${GLSL_ELLIPSOID}
float fa5Packed(sampler2D tex, float w, int i) {
  int wi = int(w);
  int ti = i / 4;
  int c = i - ti * 4;
  vec4 v = texelFetch(tex, ivec2(ti - (ti / wi) * wi, ti / wi), 0);
  return c == 0 ? v.x : (c == 1 ? v.y : (c == 2 ? v.z : v.w));
}
vec4 fa5Static(sampler2D tex, float w, int i) {
  int wi = int(w);
  return texelFetch(tex, ivec2(i - (i / wi) * wi, i / wi), 0);
}
// level (model heights, before the offset) for a cell: depth d, static st
float fa5Level(float mode, float exag, float d, vec4 st) {
  if (mode > 1.5) return st.x;
  if (mode > 0.5) return st.x + d * exag;
  return st.y + d;
}
// the value the surface is coloured by
float fa5ColorValue(float cmode, float d, vec4 st, float cv) {
  if (cmode < 0.5) return d;
  if (cmode < 1.5) return st.y + d;
  if (cmode < 2.5) return cv;
  return st.x;
}
// geodetic (ellipsoid) normal at an ECEF position
vec3 fa5Normal(vec3 p) {
  return normalize(vec3(p.x / FA5_A2, p.y / FA5_A2, p.z / FA5_B2));
}
// colour-map lookup coordinate for value v in [cmin, cmax]
vec2 fa5MapCoord(float v, float cmin, float cmax) {
  float t = clamp((v - cmin) / max(cmax - cmin, 1e-9), 0.0, 1.0);
  return vec2(t * 0.99609375 + 0.001953125, 0.5);
}
// Move a vertex (model coordinates, metres) toward the camera by a fraction
// of its distance: its place on screen is unchanged, only its depth. Distant
// 3D tiles are coarse and stand metres above the ground they simplify, so
// without this they hide shallow water when zoomed out. 'bias' is the
// fraction reached beyond 4 km; it is 0 within 800 m (close views keep exact
// occlusion: semi-submerged buildings).
vec3 fa5DepthBias(vec3 p, vec3 cam, float bias) {
  vec3 v = cam - p;
  float d = length(v);
  return p + v * (bias * smoothstep(800.0, 4000.0, d));
}
`;
}

/** GLSL ES 3.0 helpers for fragment shaders only (screen-space derivatives). */
export function glslFragmentHelpers() {
  return `
// Cell-border line weight (0–1) from the interpolated 'edge' (0 on the cell's
// boundary): about one pixel wide, faded out where cells are under ~8 pixels
// across (they would turn the surface white). Call before any discard.
float fa5Border(float edge, float on) {
  float w = max(fwidth(edge), 1e-6);
  if (on < 0.5) return 0.0;
  return (1.0 - smoothstep(0.5 * w, 1.5 * w, edge)) * (1.0 - smoothstep(0.15, 0.3, w));
}
`;
}
