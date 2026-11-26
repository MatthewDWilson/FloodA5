// worker.js — reads a FloodA5 output file from this computer in a Web Worker
// (a module worker; started by source.js). h5wasm (HDF5 in WebAssembly) opens
// the File through Emscripten's WORKERFS: the file is read where it lies, in
// the pieces needed (a frame, a cell's series), never copied whole into
// memory and never uploaded. localrun.js answers the viewer's /api requests.
//
// Messages in:  {op: 'open', file: File, h5wasmUrl, viewer}  → {op: 'open', ok, info} | {op: 'open', ok: false, error}
//               {op: 'fetch', id, path}                       → {id, status, json | bytes | text}
//               {op: 'close'}
// A file that fails to open leaves the one open before it in place.

import { h5Reader } from './h5reader.js';
import { LocalRun } from './localrun.js';

let h5wasm = null, FS = null, cur = null, n = 0;     // cur: {dir, file, run}

async function open({ file: f, h5wasmUrl, viewer }) {
  if (!h5wasm) {
    const m = await import(h5wasmUrl);
    h5wasm = m.default || m.h5wasm || m;
  }
  ({ FS } = await h5wasm.ready);
  const dir = `/work${++n}`;            // each file its own mount: the previous stays usable until this one opens
  FS.mkdir(dir);
  FS.mount(FS.filesystems.WORKERFS, { files: [f] }, dir);
  let file = null;
  try {
    file = new h5wasm.File(`${dir}/${f.name}`, 'r');
    const run = new LocalRun(h5Reader(h5wasm, file), { name: f.name, viewer });
    release(cur);
    cur = { dir, file, run };
    return run.info();
  } catch (e) {
    release({ dir, file });
    throw e;
  }
}

function release(c) {
  if (!c) return;
  if (c.file) { try { c.file.close(); } catch { /* already closed */ } }
  try { FS.unmount(c.dir); FS.rmdir(c.dir); } catch { /* not mounted */ }
}

self.onmessage = async ev => {
  const m = ev.data || {};
  if (m.op === 'open') {
    try {
      self.postMessage({ op: 'open', ok: true, info: await open(m) });
    } catch (e) {
      self.postMessage({ op: 'open', ok: false, error: String(e?.message || e) });
    }
  } else if (m.op === 'fetch') {
    let r;
    try {
      r = cur ? cur.run.handle(m.path) : { status: 503, text: 'No file open' };
    } catch (e) {
      r = { status: 500, text: String(e?.message || e) };
    }
    self.postMessage({ id: m.id, ...r }, r.bytes ? [r.bytes] : []);
  } else if (m.op === 'close') {
    release(cur); cur = null;
  }
};
