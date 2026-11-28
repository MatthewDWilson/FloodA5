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
// Smooth playback (0.9.2): depth and cvar come twice, for two frames (texture
// slots 0 and 1: frame k in slot k mod 2, interp.js), blended with the weight
// `w` of slot 1 (0 or 1: one frame exactly). Depths blend linearly, so a
// cell turns wet when the blended depth reaches the wet threshold — its
// crossing time under that linear depth — and the front advances cell by cell
// instead of a thin film fading in ahead of it. Class numbers never blend
// (`blend` 0: the nearer frame's).
// The vertex shader moves each vertex along the ellipsoid normal to the level
//   mode 0 (water surface):  datum + depth
//   mode 1 (exaggerated):    bed + depth × exag
//   mode 2 (bed):            bed
// plus the vertical offset; the fragment shader drops dry cells, colours
// through a 256-entry colour map and, if `borders` is on, draws the cell
// borders as white hairlines.
//
// Smooth surface (0.8.14, style `smooth`): each ring vertex takes the mean of
// the wet cells at its corner (a fourth texture, one corner per texel: level,
// bed, depth, colour value — geometry.js cornerValues), the centre keeps its
// cell's own value; colours blend too (`blend`, off for classed palettes),
// and the surface is shaded by its slope (exaggerated `shadeExag` times for
// the shading only), lit from the north-west. Shading (0.9.2, style `shade`):
// 1 from the slopes at the vertices (geometry.js surfaceGradients, a fifth
// texture with a slot per corner and per cell centre), interpolated across
// each triangle — no facets; 2 from each triangle's own slope (screen-space
// derivatives, the 0.8.14 shading: faceted; kept for comparison); 0 none.
// With two frames, corner values and slopes blend where the corner (or the
// centre) is wet in both frames and are the wet one's otherwise.

/** What the surface is coloured by (uniform `cmode`). */
export const COLOR_MODE = { depth: 0, level: 1, texture: 2, bed: 3 };

/** Surface geometry (uniform `mode`). */
export const SURFACE_MODE = { water: 0, exaggerated: 1, bed: 2 };

/** Water drawn this fraction of its distance nearer the camera when zoomed out (fa5DepthBias). */
export const DEPTH_BIAS = 0.015;

/**
 * Distances (m) over which that bias comes in, by basemap: 3D tiles stay
 * exact within 800 m (semi-submerged buildings); terrain tiles are coarse at
 * middle distances too — at a few hundred metres they stand metres above the
 * ground and would hide the water (0.8.15).
 */
export const DEPTH_BIAS_RANGE = { tiles: [800, 4000], terrain: [30, 600] };

/** Slopes are exaggerated this many times for the smooth surface's shading: water slopes are small, the bed's are not. */
export const SHADE_EXAG = { water: 20, bed: 3 };

/** Shading of the smooth surface (style `shade`). */
export const SHADE = { off: 0, smooth: 1, facets: 2 };

/** Default style of a surface. */
export function defaultStyle(role) {
  const bed = role === 'bed';
  return {
    mode: bed ? SURFACE_MODE.bed : SURFACE_MODE.water,
    cmode: bed ? COLOR_MODE.bed : COLOR_MODE.depth,
    cmin: 0, cmax: 1, alpha: bed ? 1 : 0.8, minDepth: 0.05, exag: 10, offset: 0, borders: 0,
    depthBias: bed ? 0 : DEPTH_BIAS, biasNear: DEPTH_BIAS_RANGE.tiles[0], biasFar: DEPTH_BIAS_RANGE.tiles[1],
    colormap: bed ? 'terrain' : 'depth',
    smooth: 0, blend: 1, shadeExag: bed ? SHADE_EXAG.bed : SHADE_EXAG.water, shade: SHADE.smooth,
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
// smooth surface: corner c's values (level, bed, depth, colour value)
vec4 fa5Corner(sampler2D tex, float w, int c) {
  int wi = int(w);
  return texelFetch(tex, ivec2(c - (c / wi) * wi, c / wi), 0);
}
float fa5LevelC(float mode, float exag, vec4 cv) {
  if (mode > 1.5) return cv.y;
  if (mode > 0.5) return cv.y + cv.z * exag;
  return cv.x;
}
float fa5ColorValueC(float cmode, vec4 cv) {
  if (cmode < 0.5) return cv.z;
  if (cmode < 1.5) return cv.x;
  if (cmode < 2.5) return cv.w;
  return cv.y;
}
// two frames' values at a vertex (slot 1 weighted w): blended where both are
// wet (or both dry: not drawn), the wet one's otherwise; w = 0 or 1: that frame
vec4 fa5Pick(vec4 a, vec4 b, bool wa, bool wb, float w) {
  if (w <= 0.0) return a;
  if (w >= 1.0) return b;
  if (wa == wb) return mix(a, b, w);
  return wa ? a : b;
}
// smooth shading: the slope (east, north) of the surface drawn in mode 'mode'
// from a slope texel (level slope, bed slope); depth × exag is bed + exag ×
// (level − bed): exact where the depth datum is the bed (not with SGS)
vec2 fa5GradMode(float mode, float exag, vec4 g) {
  if (mode > 1.5) return g.zw;
  if (mode > 0.5) return mix(g.zw, g.xy, exag);
  return g.xy;
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
// 3D tiles and terrain are coarse and stand metres above the ground they
// simplify, so without this they hide shallow water when zoomed out. 'bias'
// is the fraction reached beyond range.y metres; it is 0 within range.x
// (3D tiles: 800 m – 4 km, so close views keep exact occlusion of
// semi-submerged buildings; terrain: 30 – 600 m, DEPTH_BIAS_RANGE).
vec3 fa5DepthBias(vec3 p, vec3 cam, float bias, vec2 range) {
  vec3 v = cam - p;
  float d = length(v);
  return p + v * (bias * smoothstep(range.x, range.y, d));
}
`;
}

/**
 * GLSL ES 3.0 vertex-shader body shared by the engines: from the cell index
 * `cell` (and the smooth surface's `corner`, −1 at a centre; null where the
 * engine's model has no corner attribute) declares and sets
 *   float dA, dB, d (blended depth), level (model heights, before the offset),
 *   cval (colour value); vec2 grad (the shading slope, east and north)
 * `u(name)` gives the engine's uniform name for: depth, depth1, cvar, cvar1,
 * static, texW, staticW, corner, corner1, cornerW, grad, grad1, gradW,
 * nCorners, w, minDepth, mode, exag, cmode, blend, smooth, shade.
 */
export function glslSurfaceVertex(u, cellExpr, cornerExpr) {
  const smooth = cornerExpr ? `
  float fa5c = ${cornerExpr};
  if (${u('smooth')} > 0.5) {
    bool atCorner = fa5c > -0.5;
    bool wA = dA >= ${u('minDepth')}, wB = dB >= ${u('minDepth')};
    if (atCorner) {
      int c = int(fa5c + 0.5);
      vec4 ccA = fa5Corner(${u('corner')}, ${u('cornerW')}, c), ccB = fa5Corner(${u('corner1')}, ${u('cornerW')}, c);
      wA = ccA != vec4(0.0); wB = ccB != vec4(0.0);       // (a corner with no wet cell is all zeros)
      vec4 cc = fa5Pick(ccA, ccB, wA, wB, w);
      level = fa5LevelC(${u('mode')}, ${u('exag')}, cc);
      if (${u('blend')} > 0.5) cval = fa5ColorValueC(${u('cmode')}, cc);
    }
    if (${u('shade')} > 0.5 && ${u('shade')} < 1.5) {
      int g = atCorner ? int(fa5c + 0.5) : int(${u('nCorners')} + 0.5) + i;
      vec4 gA = fa5Corner(${u('grad')}, ${u('gradW')}, g), gB = fa5Corner(${u('grad1')}, ${u('gradW')}, g);
      grad = fa5GradMode(${u('mode')}, ${u('exag')}, fa5Pick(gA, gB, wA, wB, w));
    }
  }` : '';
  return `
  int i = int(${cellExpr} + 0.5);
  vec4 st = fa5Static(${u('static')}, ${u('staticW')}, i);
  float w = ${u('w')};
  float dA = max(fa5Packed(${u('depth')}, ${u('texW')}, i), 0.0);
  float dB = max(fa5Packed(${u('depth1')}, ${u('texW')}, i), 0.0);
  float d = mix(dA, dB, w);
  float level = fa5Level(${u('mode')}, ${u('exag')}, d, st);
  float cv = 0.0;
  if (${u('cmode')} > 1.5 && ${u('cmode')} < 2.5) {
    float cA = fa5Packed(${u('cvar')}, ${u('texW')}, i), cB = fa5Packed(${u('cvar1')}, ${u('texW')}, i);
    cv = ${u('blend')} > 0.5 ? mix(cA, cB, w) : (w < 0.5 ? cA : cB);
  }
  float cval = fa5ColorValue(${u('cmode')}, d, st, cv);
  vec2 grad = vec2(0.0);${smooth}
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
// Smooth surface: brightness from the slope (flat ground = 1), from the
// screen-space derivatives of the position p (metres) with up the local
// vertical; heights exaggerated k times for the shading; light from the
// north-west, 45° up, as on hill-shaded maps.
float fa5Shade(vec3 p, vec3 up, float k) {
  vec3 t1 = dFdx(p), t2 = dFdy(p);
  t1 += up * (dot(t1, up) * (k - 1.0));
  t2 += up * (dot(t2, up) * (k - 1.0));
  vec3 n = cross(t1, t2);
  float len = length(n);
  if (len < 1e-12) return 1.0;
  n /= len;
  if (dot(n, up) < 0.0) n = -n;
  vec3 east = normalize(cross(vec3(0.0, 0.0, 1.0), up));
  vec3 north = cross(up, east);
  vec3 l = normalize(up + 0.7071 * (north - east));
  return clamp(dot(n, l) / dot(up, l), 0.55, 1.35);
}
// Smooth surface, shading without facets: brightness from the slope g (east,
// north; geometry.js surfaceGradients, interpolated from the vertices), the
// same light and exaggeration as fa5Shade
float fa5ShadeSlope(vec2 g, float k) {
  vec3 n = normalize(vec3(-k * g, 1.0));
  vec3 l = vec3(-0.5, 0.5, 0.70710678);       // normalize(−east + north + √2 up)
  return clamp(dot(n, l) / l.z, 0.55, 1.35);
}
// the smooth surface's shading factor: shade 1 from the slope g, 2 per triangle
float fa5Shading(float shade, vec2 g, vec3 p, vec3 up, float k) {
  if (shade < 0.5) return 1.0;
  return shade < 1.5 ? fa5ShadeSlope(g, k) : fa5Shade(p, up, k);
}
`;
}
