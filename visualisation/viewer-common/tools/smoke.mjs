#!/usr/bin/env node
// tools/smoke.mjs — headless browser check of the globe viewer.
//
//   node visualisation/cesium/tools/smoke.mjs --url http://localhost:8091/ --out shot.png
//        [--browser PATH] [--timeout 120] [--settle 8] [--size 1400x900]
//        [--eval "JS" ...]          (each run after the viewer is ready, then a
//                                    screenshot shot_1.png, shot_2.png, …)
//        [--file RUN.h5 ...]        open local files in the page, in turn (as "Open file…"):
//                                   for pages started without a run, or to replace it
//        [--files-last]             open them after the --eval steps (default: before)
//
// Starts Edge/Chrome headless with a DevTools port, opens the viewer, waits
// until it reports a rendered frame (window.__flooda5), takes screenshots and
// prints the viewer state and any console errors as JSON. Exit code 1 if the
// viewer reported an error or never showed a frame. Needs Node ≥ 22 (built-in
// WebSocket); no npm packages. Not part of CI (needs a browser and network).

import { spawn, execFileSync } from 'node:child_process';
import { mkdtempSync, writeFileSync, existsSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const args = process.argv.slice(2);
const opt = (name, dflt) => { const i = args.indexOf(name); return i >= 0 ? args[i + 1] : dflt; };
const evals = []; args.forEach((a, i) => { if (a === '--eval') evals.push(args[i + 1]); });
const url = opt('--url', 'http://localhost:8080/');
const out = opt('--out', 'smoke.png');
const timeout = +opt('--timeout', '120') * 1000;
const settle = +opt('--settle', '8') * 1000;
const localFiles = []; args.forEach((a, i) => { if (a === '--file') localFiles.push(args[i + 1]); });
const localFile = localFiles[0] || null;
const [W, H] = opt('--size', '1400x900').split('x').map(Number);
const browser = opt('--browser', [
  'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
  'C:/Program Files/Google/Chrome/Application/chrome.exe',
  '/usr/bin/google-chrome', '/usr/bin/chromium', '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
].find(p => existsSync(p)));
if (!browser) { console.error('no browser found: use --browser PATH'); process.exit(2); }

const sleep = ms => new Promise(r => setTimeout(r, ms));
const port = 9300 + Math.floor(Math.random() * 500);
const profile = mkdtempSync(join(tmpdir(), 'fa5-smoke-'));
const proc = spawn(browser, [
  '--headless=new', `--remote-debugging-port=${port}`, `--user-data-dir=${profile}`,
  `--window-size=${W},${H}`, '--no-first-run', '--no-default-browser-check',
  '--ignore-gpu-blocklist', '--enable-webgl', 'about:blank',
], { stdio: 'ignore' });

function killBrowser() {
  try {
    if (process.platform === 'win32') execFileSync('taskkill', ['/PID', String(proc.pid), '/T', '/F'], { stdio: 'ignore' });
    else proc.kill('SIGKILL');
  } catch { /* already gone */ }
  setTimeout(() => { try { rmSync(profile, { recursive: true, force: true }); } catch { /* locked */ } }, 500);
}

async function cdpTarget() {
  for (let k = 0; k < 100; k++) {
    try {
      const list = await (await fetch(`http://127.0.0.1:${port}/json/list`)).json();
      const page = list.find(t => t.type === 'page');
      if (page) return page.webSocketDebuggerUrl;
    } catch { /* not up yet */ }
    await sleep(200);
  }
  throw new Error('browser DevTools endpoint did not come up');
}

const consoleErrors = [];
const requests = new Map(), failed = [];      // network failures (query secrets redacted)
const netlog = args.includes('--netlog'), apiLog = [], tStart = Date.now();   // --netlog: the viewer's API requests
let ws, nextId = 1;
const waiting = new Map();
function send(method, params = {}) {
  const id = nextId++;
  ws.send(JSON.stringify({ id, method, params }));
  return new Promise((resolve, reject) => waiting.set(id, { resolve, reject }));
}
async function evaluate(expr) {
  const r = await send('Runtime.evaluate', { expression: expr, awaitPromise: true, returnByValue: true });
  if (r.exceptionDetails) throw new Error(r.exceptionDetails.exception?.description || r.exceptionDetails.text);
  return r.result.value;
}
async function shot(path) {
  const r = await send('Page.captureScreenshot', { format: 'png' });
  writeFileSync(path, Buffer.from(r.data, 'base64'));
  console.error(`screenshot: ${path}`);
}

let code = 0;
try {
  ws = new WebSocket(await cdpTarget());
  await new Promise((res, rej) => { ws.onopen = res; ws.onerror = rej; });
  ws.onmessage = ev => {
    const m = JSON.parse(ev.data);
    if (m.id && waiting.has(m.id)) {
      const w = waiting.get(m.id); waiting.delete(m.id);
      m.error ? w.reject(new Error(m.error.message)) : w.resolve(m.result);
    } else if (m.method === 'Network.requestWillBeSent') {
      requests.set(m.params.requestId, m.params.request.url);
      if (netlog && /\/api\//.test(m.params.request.url)) apiLog.push(`${Date.now() - tStart} send ${m.params.request.url.replace(/^https?:\/\/[^/]+/, '')}`);
    } else if (netlog && (m.method === 'Network.loadingFinished' || m.method === 'Network.responseReceived') &&
               /\/api\//.test(requests.get(m.params.requestId) || '')) {
      apiLog.push(`${Date.now() - tStart} ${m.method.slice(8)} ${requests.get(m.params.requestId).replace(/^https?:\/\/[^/]+/, '')}`);
    } else if (m.method === 'Network.loadingFailed') {
      const u = (requests.get(m.params.requestId) || '?').replace(/([?&](key|session|access_token)=)[^&]*/g, '$1…');
      failed.push(`${m.params.errorText}${m.params.canceled ? ' (canceled)' : ''}${m.params.blockedReason ? ` (${m.params.blockedReason})` : ''}: ${u.slice(0, 160)}`);
      if (netlog && /\/api\//.test(u)) apiLog.push(`${Date.now() - tStart} FAILED ${m.params.errorText} ${u}`);
    } else if (m.method === 'Runtime.exceptionThrown') {
      consoleErrors.push(m.params.exceptionDetails.exception?.description || m.params.exceptionDetails.text);
    } else if (m.method === 'Runtime.consoleAPICalled' && (m.params.type === 'error' || m.params.type === 'warning')) {
      consoleErrors.push(`${m.params.type}: ${m.params.args.map(a => a.value ?? a.description ?? '').join(' ')}`);
    }
  };
  await send('Runtime.enable');
  await send('Page.enable');
  await send('Network.enable');
  await send('Emulation.setDeviceMetricsOverride', { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  await send('Page.navigate', { url });
  const t0 = Date.now();
  let st = null;
  while (Date.now() - t0 < timeout) {
    await sleep(500);
    try { st = await evaluate('JSON.stringify(window.__flooda5 || null)').then(JSON.parse); } catch { continue; }
    if (st && (st.state === 'error' || (st.shownFrame !== null && st.renders > 0))) break;
    if (localFile && st && st.state === 'ready') break;           // (no run yet: open the file)
  }
  const filesLast = args.includes('--files-last');
  const openFiles = async () => { for (const f of (st && st.state !== 'error' ? localFiles : [])) {   // one after another, as a user would
    const { root } = await send('DOM.getDocument');
    const { nodeId } = await send('DOM.querySelector', { nodeId: root.nodeId, selector: '#file-input' });
    await evaluate("window.__flooda5.local && (window.__flooda5.local.state = 'chosen'), 0");
    await send('DOM.setFileInputFiles', { nodeId, files: [f] });
    const t2 = Date.now();
    while (Date.now() - t2 < timeout) {
      await sleep(500);
      try { st = await evaluate('JSON.stringify(window.__flooda5 || null)').then(JSON.parse); } catch { continue; }
      const l = st.local;
      if (st.state === 'error' || l?.state === 'error' || (l?.state === 'open' && st.shownFrame !== null)) break;
    }
    console.error(`local file: ${f} → ${JSON.stringify(st.local)}`);
  } };
  if (!filesLast) await openFiles();
  await sleep(settle);
  await shot(out);
  for (let k = 0; k < evals.length; k++) {
    console.error(`eval: ${evals[k]}`);
    try { const v = await evaluate(evals[k]); if (v !== undefined) console.error(`  → ${JSON.stringify(v)}`); }
    catch (e) { console.error(`  ✗ ${e.message}`); code = 1; }
    const t1 = Date.now();
    await sleep(1000);
    while (Date.now() - t1 < settle * 3) {      // let tiles stream in
      try { if (await evaluate('window.__flooda5api && window.__flooda5api.tilesLoaded()')) break; } catch { break; }
      await sleep(500);
    }
    await sleep(settle);
    await shot(out.replace(/\.png$/i, `_${k + 1}.png`));
  }
  if (filesLast) { await openFiles(); await sleep(settle); await shot(out.replace(/\.png$/i, '_files.png')); }
  const profSecs = +opt('--profile', '0');        // --profile N: CPU profile, top self time
  let profile;
  if (profSecs > 0) {
    await send('Profiler.enable');
    await send('Profiler.start');
    await sleep(profSecs * 1000);
    const { profile: p } = await send('Profiler.stop');
    const dt = (p.endTime - p.startTime) / p.samples.length / 1000;    // ms per sample
    const self = new Map(), byId = new Map(p.nodes.map(nd => [nd.id, nd]));
    for (const id of p.samples) {
      const cf = byId.get(id).callFrame;
      const key = `${cf.functionName || '(anon)'} ${cf.url.split('/').slice(-2).join('/')}:${cf.lineNumber}`;
      self.set(key, (self.get(key) || 0) + dt);
    }
    profile = [...self].sort((a, b) => b[1] - a[1]).slice(0, 15).map(([k, v]) => `${v.toFixed(0).padStart(6)} ms  ${k.slice(0, 140)}`);
  }
  st = await evaluate('JSON.stringify(window.__flooda5 || null)').then(JSON.parse);
  console.log(JSON.stringify({ status: st, ...(profile ? { profile } : {}), consoleErrors, failedRequests: failed.slice(0, 20),
                               failedCount: failed.length, ...(netlog ? { apiLog } : {}) }, null, 2));
  if (!st || st.state === 'error' || st.shownFrame === null || st.errors.length) code = 1;
} catch (e) {
  console.error(`smoke test failed: ${e.message}`);
  code = 1;
} finally {
  try { ws && ws.close(); } catch { /* ignore */ }
  killBrowser();
}
process.exit(code);
