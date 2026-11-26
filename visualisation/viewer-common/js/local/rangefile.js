// rangefile.js — a file on a web server read in pieces (HTTP range requests),
// for h5wasm in the worker (worker.js): HDF5 reads a few kilobytes here and
// there (metadata, one compressed chunk of a frame), so only those parts are
// downloaded. Pure (the fetching is passed in): unit-tested with `node --test`.

/**
 * Reads served from fixed-size blocks of the file, fetched on demand and
 * kept in a least-recently-used cache. Adjacent missing blocks of one read
 * are fetched with one request.
 */
export class RangeCache {
  /**
   * @param fetchRange (from, to) → Uint8Array of bytes from…to (inclusive); throws on failure
   * @param size       the file's length (bytes)
   * @param block      block size (bytes)
   * @param maxBytes   cache limit (bytes; at least the blocks of one read are kept)
   * @param maxRequest largest single request (bytes)
   */
  constructor(fetchRange, size, { block = 64 * 1024, maxBytes = 256 * 2 ** 20, maxRequest = 16 * 2 ** 20 } = {}) {
    if (!(size >= 0)) throw new Error('RangeCache: the size must be known');
    this.fetchRange = fetchRange;
    this.size = size;
    this.block = block;
    this.maxBytes = maxBytes;
    this.maxBlocks = Math.max(1, Math.floor(maxRequest / block));
    this.cache = new Map();          // block index → Uint8Array (insertion order = recency)
    this.bytes = 0;
    this.requests = 0;
    this.fetched = 0;
  }

  _get(b) {
    const v = this.cache.get(b);
    if (v) { this.cache.delete(b); this.cache.set(b, v); }
    return v;
  }

  _put(b, v) {
    this.cache.set(b, v);
    this.bytes += v.length;
  }

  _evict(keep) {
    for (const [b, v] of this.cache) {
      if (this.bytes <= this.maxBytes) break;
      if (keep.has(b)) continue;
      this.cache.delete(b);
      this.bytes -= v.length;
    }
  }

  /** Fetch blocks first…last (inclusive) with one request. */
  _fetch(first, last) {
    const from = first * this.block, to = Math.min((last + 1) * this.block, this.size) - 1;
    const data = this.fetchRange(from, to);
    if (!data || data.length !== to - from + 1) {
      throw new Error(`range ${from}–${to}: got ${data ? data.length : 0} bytes, expected ${to - from + 1}`);
    }
    this.requests++; this.fetched += data.length;
    for (let b = first; b <= last; b++) {
      const o = (b - first) * this.block;
      this._put(b, data.subarray(o, Math.min(o + this.block, data.length)));
    }
  }

  /**
   * Copy up to `length` bytes at `position` into `out` (any typed array of
   * bytes) at `offset`. → the number of bytes copied (fewer at the end of
   * the file).
   */
  read(out, offset, length, position) {
    if (position >= this.size || length <= 0) return 0;
    const n = Math.min(length, this.size - position);
    const B = this.block, first = Math.floor(position / B), last = Math.floor((position + n - 1) / B);
    // fetch the missing blocks, adjacent ones together
    let run = -1;
    for (let b = first; b <= last + 1; b++) {
      const missing = b <= last && !this.cache.has(b);
      if (missing && run < 0) run = b;
      if (run >= 0 && (!missing || b - run >= this.maxBlocks)) {
        this._fetch(run, b - 1);
        run = missing ? b : -1;
      }
    }
    let done = 0;
    for (let b = first; b <= last; b++) {
      const v = this._get(b);
      const start = b === first ? position - b * B : 0;
      const end = Math.min(v.length, position + n - b * B);
      out.set(v.subarray(start, end), offset + done);
      done += end - start;
    }
    const keep = new Set();
    for (let b = first; b <= last; b++) keep.add(b);
    this._evict(keep);
    return done;
  }
}

/** The file's name from its URL (the last path segment), for labels. */
export function urlFileName(url) {
  try {
    const u = new URL(url);
    const seg = u.pathname.split('/').filter(Boolean).pop();
    return seg ? decodeURIComponent(seg) : u.hostname;
  } catch {
    return 'remote.h5';
  }
}

/**
 * A link the viewer may open: an absolute http(s) URL, or a path relative to
 * the page (resolved against `base`). → the absolute URL, or null.
 */
export function resolveRunUrl(text, base) {
  const s = String(text || '').trim();
  if (!s) return null;
  let u;
  try { u = new URL(s, base); } catch { return null; }
  return u.protocol === 'http:' || u.protocol === 'https:' ? u.href : null;
}
