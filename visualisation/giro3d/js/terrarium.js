// terrarium.js — Terrarium elevation tiles (AWS Open Data, Mapzen/Tilezen;
// no key): a Giro3D image format, and point sampling for the offset estimate.
//
//   https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png
//   elevation (m) = R × 256 + G + B / 256 − 32768
//
// Ported from the Open Geo Platform prototype (3DGeo, TerrariumFormat.ts).

import { DataTexture, FloatType, LinearFilter, RGFormat } from 'three';
import ImageFormat from '@giro3d/giro3d/formats/ImageFormat.js';

export const TERRARIUM_URL = 'https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png';
export const TERRARIUM_MAX_ZOOM = 15;

async function pixels(blob) {
  const image = await createImageBitmap(blob);
  const canvas = new OffscreenCanvas(image.width, image.height);
  const ctx = canvas.getContext('2d', { willReadFrequently: true });
  ctx.drawImage(image, 0, 0);
  return { data: ctx.getImageData(0, 0, image.width, image.height).data, width: image.width, height: image.height };
}

const decode = (d, k) => d[k] * 256 + d[k + 1] + d[k + 2] / 256 - 32768;

/** Giro3D image format for Terrarium tiles (as its MapboxTerrainFormat). */
export class TerrariumFormat extends ImageFormat {
  constructor() {
    super(true, FloatType);
    this.isTerrariumFormat = true;
    this.type = 'TerrariumFormat';
  }

  async decode(blob, options) {
    const { data, width, height } = await pixels(blob);
    const n = width * height, out = new Float32Array(2 * n), noData = options?.noDataValue;
    let min = Infinity, max = -Infinity;
    for (let i = 0; i < n; i++) {
      const e = decode(data, 4 * i);
      out[2 * i] = e;
      if (noData != null && noData === e) out[2 * i + 1] = 0;
      else { out[2 * i + 1] = 1; if (e < min) min = e; if (e > max) max = e; }
    }
    const texture = new DataTexture(out, width, height, RGFormat, FloatType);
    texture.needsUpdate = true;
    texture.generateMipmaps = false;
    texture.magFilter = texture.minFilter = LinearFilter;
    return { texture, min: Number.isFinite(min) ? min : 0, max: Number.isFinite(max) ? max : 0 };
  }
}

/** Elevations (m) at points [{lon, lat}] from zoom-z Terrarium tiles (NaN on failure). */
export async function sampleTerrarium(points, z = 14) {
  const tiles = new Map(), out = new Float64Array(points.length).fill(NaN);
  const n = 2 ** z;
  const where = points.map(({ lon, lat }) => {
    const x = ((lon + 180) / 360) * n;
    const la = (lat * Math.PI) / 180;
    const y = ((1 - Math.log(Math.tan(la) + 1 / Math.cos(la)) / Math.PI) / 2) * n;
    return { tx: Math.floor(x), ty: Math.floor(y), fx: x - Math.floor(x), fy: y - Math.floor(y) };
  });
  for (const w of where) {
    const key = `${w.tx}/${w.ty}`;
    if (!tiles.has(key)) {
      const url = TERRARIUM_URL.replace('{z}', z).replace('{x}', w.tx).replace('{y}', w.ty);
      tiles.set(key, fetch(url).then(r => (r.ok ? r.blob() : null)).then(b => (b ? pixels(b) : null)).catch(() => null));
    }
  }
  for (let k = 0; k < points.length; k++) {
    const w = where[k], px = await tiles.get(`${w.tx}/${w.ty}`);
    if (!px) continue;
    const i = Math.min(px.width - 1, Math.floor(w.fx * px.width));
    const j = Math.min(px.height - 1, Math.floor(w.fy * px.height));
    out[k] = decode(px.data, 4 * (j * px.width + i));
  }
  return out;
}
