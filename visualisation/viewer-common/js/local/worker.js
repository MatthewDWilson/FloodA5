// worker.js — reads a FloodA5 output file in a Web Worker (a module worker;
// started by source.js), for the viewer's /api requests (localrun.js), with
// h5wasm (HDF5 in WebAssembly):
//   - a file on this computer (a File): through Emscripten's WORKERFS, read
//     where it lies, in the pieces needed, never uploaded;
//   - a file on a web server (a URL, 0.8.15): read in pieces with HTTP range
//     requests (rangefile.js) when the server allows them, otherwise
//     downloaded whole into memory.
//
// Messages in:  {op: 'open', file: File | url: string, h5wasmUrl, viewer}
//                 → {op: 'open', ok, info} | {op: 'open', ok: false, error}
//               {op: 'fetch', id, path}  → {id, status, json | bytes | text}
//               {op: 'stats', id}        → {id, stats}   (a URL's requests and bytes so far)
//               {op: 'close'}
// Messages out also: {op: 'progress', loaded, total} while a whole file downloads.
// A file that fails to open leaves the one open before it in place.

import { h5Reader } from './h5reader.js';
import { LocalRun } from './localrun.js';
import { RangeCache, urlFileName } from './rangefile.js';

let h5wasm = null, FS = null, cur = null, n = 0;     // cur: {dir, path, kind, file, run, rc}
let lastError = '';

async function loadH5wasm(h5wasmUrl) {
  if (!h5wasm) {
    const m = await import(h5wasmUrl);
    h5wasm = m.default || m.h5wasm || m;
  }
  ({ FS } = await h5wasm.ready);
}

async function open({ file: f, url, h5wasmUrl, viewer }) {
  await loadH5wasm(h5wasmUrl);
  const dir = `/work${++n}`;            // each file its own directory: the previous stays usable until this one opens
  const c = { dir, kind: null, file: null };
  let name;
  try {
    if (f) {
      name = f.name;
      FS.mkdir(dir);
      FS.mount(FS.filesystems.WORKERFS, { files: [f] }, dir);
      c.kind = 'workerfs';
    } else {
      name = urlFileName(url);
      const how = await probe(url);
      FS.mkdir(dir);
      if (how.ranged) {
        c.rc = new RangeCache(syncRange(url), how.size);
        lazyFile(dir, name, c.rc);
      } else {
        FS.writeFile(`${dir}/${name}`, how.data);
      }
      c.kind = 'memfs';
      c.remote = { url, ranged: how.ranged, size: how.size };
    }
    c.path = `${dir}/${name}`;
    lastError = '';
    try {
      c.file = new h5wasm.File(c.path, 'r');
    } catch (e) {
      throw new Error(`${name}: not readable as HDF5 (${lastError || e?.message || e})`);
    }
    c.run = new LocalRun(h5Reader(h5wasm, c.file), { name, viewer, source: url || name });
    release(cur);
    cur = c;
    const info = c.run.info();
    if (c.remote) info.remote = { ...c.remote };
    return info;
  } catch (e) {
    release(c);
    throw e;
  }
}

/**
 * How to read `url`: a range request for its first bytes. 206: the server
 * serves ranges (the size from Content-Range — exposed to cross-origin pages
 * only if the server says so — else from a HEAD request's Content-Length);
 * 200: it does not, and the answer is the whole file.
 */
async function probe(url) {
  let r;
  try {
    r = await fetch(url, { headers: { Range: 'bytes=0-7' }, cache: 'no-store' });
  } catch (e) {
    throw new Error(`cannot fetch ${url} (offline, or the server does not let other sites read it: CORS)`);
  }
  if (r.status === 206) {
    let size = Number((r.headers.get('Content-Range') || '').split('/')[1]);
    await r.arrayBuffer().catch(() => {});
    if (!(size > 0)) {
      const h = await fetch(url, { method: 'HEAD', cache: 'no-store' }).catch(() => null);
      size = h && h.ok ? Number(h.headers.get('Content-Length')) : NaN;
    }
    if (size > 0) return { ranged: true, size };
    r = await fetch(url, { cache: 'no-store' });          // size unknown: the whole file
  }
  if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
  return { ranged: false, ...(await download(r)) };
}

/** The whole body of response `r`, reporting progress. */
async function download(r) {
  const total = Number(r.headers.get('Content-Length')) || 0;
  if (!r.body) { const data = new Uint8Array(await r.arrayBuffer()); return { data, size: data.length }; }
  const parts = []; let loaded = 0, last = 0;
  const reader = r.body.getReader();
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    parts.push(value); loaded += value.length;
    if (loaded - last > 8 * 2 ** 20) { last = loaded; self.postMessage({ op: 'progress', loaded, total }); }
  }
  const data = new Uint8Array(loaded);
  let o = 0;
  for (const p of parts) { data.set(p, o); o += p.length; }
  return { data, size: loaded };
}

/** A synchronous range request (workers may block): bytes from…to of `url`. */
const syncRange = url => (from, to) => {
  const x = new XMLHttpRequest();
  x.open('GET', url, false);
  x.setRequestHeader('Range', `bytes=${from}-${to}`);
  x.responseType = 'arraybuffer';
  x.send();
  if (x.status !== 206) throw new Error(`range request ${from}–${to}: HTTP ${x.status}`);
  return new Uint8Array(x.response);
};

/** A read-only file `dir/name` whose reads come from `rc` (a RangeCache). */
function lazyFile(dir, name, rc) {
  const node = FS.createFile(dir, name, {}, true, false);
  node.contents = null;
  Object.defineProperty(node, 'usedBytes', { get: () => rc.size, configurable: true });
  const ops = { ...node.stream_ops };
  ops.read = (stream, buffer, offset, length, position) => {
    try {
      return rc.read(buffer, offset, length, position);
    } catch (e) {
      lastError = String(e?.message || e);
      throw new FS.ErrnoError(29);                         // EIO
    }
  };
  ops.write = () => { throw new FS.ErrnoError(29); };
  ops.mmap = () => { throw new FS.ErrnoError(43); };       // ENODEV
  node.stream_ops = ops;
}

function release(c) {
  if (!c) return;
  if (c.file) { try { c.file.close(); } catch { /* already closed */ } }
  if (c.kind === 'workerfs') {
    try { FS.unmount(c.dir); } catch { /* not mounted */ }
  } else if (c.kind === 'memfs' && c.path) {
    try { FS.unlink(c.path); } catch { /* gone */ }
  }
  try { FS.rmdir(c.dir); } catch { /* not there */ }
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
    lastError = '';
    try {
      r = cur ? cur.run.handle(m.path) : { status: 503, text: 'No file open' };
    } catch (e) {
      r = { status: 500, text: `${e?.message || e}${lastError ? ` (${lastError})` : ''}` };
    }
    self.postMessage({ id: m.id, ...r }, r.bytes ? [r.bytes] : []);
  } else if (m.op === 'stats') {
    const rc = cur?.rc;
    self.postMessage({ id: m.id, stats: cur?.remote ? { ...cur.remote, requests: rc ? rc.requests : 1,
      fetched: rc ? rc.fetched : cur.remote.size, cached: rc ? rc.bytes : cur.remote.size } : null });
  } else if (m.op === 'close') {
    release(cur); cur = null;
  }
};
