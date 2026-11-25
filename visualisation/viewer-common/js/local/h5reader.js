// h5reader.js — the `reader` of localrun.js over an h5wasm File (HDF5
// compiled to WebAssembly). Values come back as plain JavaScript: numeric
// datasets as typed arrays (64-bit integers as Float64Array), attributes as
// numbers, strings, booleans (h5py-style FALSE/TRUE enums) or arrays.

const BITFIELD = 4;            // H5T_class_t.H5T_BITFIELD

const plainAttr = v => {
  if (typeof v === 'bigint') return Number(v);
  if (ArrayBuffer.isView(v)) return Array.from(v, x => (typeof x === 'bigint' ? Number(x) : x));
  if (Array.isArray(v)) return v.map(plainAttr);
  return v;
};

const plainData = v => {
  if (v instanceof BigInt64Array || v instanceof BigUint64Array) return Float64Array.from(v, Number);
  if (typeof v === 'bigint') return Number(v);
  return v;
};

/**
 * @param h5wasm  the h5wasm module (its exports: File, Dataset, Group, …)
 * @param file    an open h5wasm.File
 */
export function h5Reader(h5wasm, file) {
  const get = path => (path === '/' ? file : file.get(path));
  const node = path => {
    try { return get(path); } catch { return null; }
  };
  const isDataset = obj => obj instanceof h5wasm.Dataset;
  const isGroup = obj => obj instanceof h5wasm.Group;
  const dataset = path => {
    const d = node(path);
    if (!isDataset(d)) throw new Error(`${path}: no such dataset`);
    return d;
  };
  return {
    has: path => { const o = node(path); return isDataset(o) || isGroup(o); },
    isDataset: path => isDataset(node(path)),
    keys: path => { const g = node(path); return isGroup(g) ? g.keys() : []; },
    attrs: path => {
      const o = node(path);
      if (!o) return {};
      const out = {};
      for (const [k, a] of Object.entries(o.attrs)) {
        const m = a.metadata, v = a.json_value;
        // HDF5.jl writes Bool as an 8-bit bitfield (H5T_NATIVE_B8), read by Julia as Bool
        if (m && m.type === BITFIELD && m.size === 1 && ArrayBuffer.isView(v)) {
          out[k] = m.shape && m.shape.length ? Array.from(v, x => !!x) : !!v[0];
        } else out[k] = plainAttr(v);
      }
      return out;
    },
    shape: path => (dataset(path).shape || []).map(Number),
    read: path => plainData(dataset(path).value),
    row: (path, k) => plainData(dataset(path).slice([[k, k + 1], []])),
    col: (path, c) => plainData(dataset(path).slice([[], [c, c + 1]])),
  };
}
