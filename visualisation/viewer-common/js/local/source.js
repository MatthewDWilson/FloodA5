// source.js — a FloodA5 output file on this computer as the viewer's data
// source: the main-thread side of worker.js. `fetch(path)` answers the
// viewer's /api requests as the server would (a real Response), so the rest
// of the viewer does not know the difference. The file is read in the
// browser (h5wasm in a worker) and never uploaded.

/** Variable descriptions without a server (generated from io/variables.toml: tools/make_about.jl). */
export const ABOUT_URL = new URL('../../about.json', import.meta.url).href;

/** Defaults of the viewer configuration without a server (also read by VisualisationServer.jl). */
export const DEFAULTS_URL = new URL('../../defaults.json', import.meta.url).href;

/** Whether `name` looks like an HDF5 file. */
export const isHdf5Name = name => /\.(h5|hdf5|he5)$/i.test(String(name || ''));

const response = r => {
  if (r.bytes) return new Response(r.bytes, { status: r.status, headers: { 'Content-Type': 'application/octet-stream' } });
  if (r.json !== undefined) return new Response(JSON.stringify(r.json), { status: r.status, headers: { 'Content-Type': 'application/json' } });
  return new Response(r.text ?? '', { status: r.status, headers: { 'Content-Type': 'text/plain; charset=utf-8' } });
};

export class LocalSource {
  /** @param h5wasmUrl  the h5wasm ES module (dist/esm/hdf5_hl.js); @param viewer  'cesium' | 'giro3d' */
  constructor({ h5wasmUrl, viewer = 'cesium' } = {}) {
    if (!h5wasmUrl) throw new Error('no h5wasm URL in the configuration');
    this.h5wasmUrl = h5wasmUrl;
    this.viewer = viewer;
    this.worker = null;
    this.pending = new Map();
    this.next = 1;
    this.opening = null;
    this.name = null;
  }

  _start() {
    if (this.worker) return;
    this.worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module', name: 'flooda5-local' });
    this.worker.onmessage = ev => {
      const m = ev.data || {};
      if (m.op === 'open') {
        const o = this.opening; this.opening = null;
        if (o) m.ok ? o.resolve(m.info) : o.reject(new Error(m.error));
        return;
      }
      const p = this.pending.get(m.id);
      if (p) { this.pending.delete(m.id); p.resolve(response(m)); }
    };
    this.worker.onerror = ev => {
      const err = new Error(ev.message || 'the file reader stopped');
      this.opening?.reject(err); this.opening = null;
      for (const p of this.pending.values()) p.reject(err);
      this.pending.clear();
    };
  }

  /** Open `file` (a File from a file input or a drop) → the run's /api/info. Replaces the file open before. */
  open(file) {
    this._start();
    if (this.opening) this.opening.reject(new Error('superseded by another file'));
    return new Promise((resolve, reject) => {
      this.opening = { resolve, reject };
      this.worker.postMessage({ op: 'open', file, h5wasmUrl: this.h5wasmUrl, viewer: this.viewer });
    }).then(info => { this.name = file.name; return info; });
  }

  /** The viewer's request `path` ("/api/…?…") → Promise<Response>. */
  fetch(path) {
    if (path.startsWith('/api/about')) return fetch(ABOUT_URL, { cache: 'no-cache' });
    if (!this.worker) return Promise.resolve(response({ status: 503, text: 'No file open' }));
    return new Promise((resolve, reject) => {
      const id = this.next++;
      this.pending.set(id, { resolve, reject });
      this.worker.postMessage({ op: 'fetch', id, path });
    });
  }

  close() {
    this.worker?.terminate();
    this.worker = null;
    for (const p of this.pending.values()) p.reject(new Error('closed'));
    this.pending.clear();
  }
}
