#!/usr/bin/env node
// tools/smoke.mjs — headless browser check of the globe viewer.
//
//   node visualisation/cesium/tools/smoke.mjs --url http://localhost:8091/ --out shot.png
//        [--browser PATH] [--timeout 120] [--settle 8] [--size 1400x900]
//        [--eval "JS" ...]          (each run after the viewer is ready, then a
//                                    screenshot shot_1.png, shot_2.png, …)
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
    } else if (m.method === 'Runtime.exceptionThrown') {
      consoleErrors.push(m.params.exceptionDetails.exception?.description || m.params.exceptionDetails.text);
    } else if (m.method === 'Runtime.consoleAPICalled' && (m.params.type === 'error' || m.params.type === 'warning')) {
      consoleErrors.push(`${m.params.type}: ${m.params.args.map(a => a.value ?? a.description ?? '').join(' ')}`);
    }
  };
  await send('Runtime.enable');
  await send('Page.enable');
  await send('Emulation.setDeviceMetricsOverride', { width: W, height: H, deviceScaleFactor: 1, mobile: false });
  await send('Page.navigate', { url });
  const t0 = Date.now();
  let st = null;
  while (Date.now() - t0 < timeout) {
    await sleep(500);
    try { st = await evaluate('JSON.stringify(window.__flooda5 || null)').then(JSON.parse); } catch { continue; }
    if (st && (st.state === 'error' || (st.shownFrame !== null && st.renders > 0))) break;
  }
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
  st = await evaluate('JSON.stringify(window.__flooda5 || null)').then(JSON.parse);
  console.log(JSON.stringify({ status: st, consoleErrors }, null, 2));
  if (!st || st.state === 'error' || st.shownFrame === null || st.errors.length) code = 1;
} catch (e) {
  console.error(`smoke test failed: ${e.message}`);
  code = 1;
} finally {
  try { ws && ws.close(); } catch { /* ignore */ }
  killBrowser();
}
process.exit(code);
