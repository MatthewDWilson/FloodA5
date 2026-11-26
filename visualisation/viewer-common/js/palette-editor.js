// palette-editor.js — the palettes of the globe viewer beyond the built-in
// ones, and their DOM: loading the scientific maps (colormaps.json), the
// server's shared palettes (/api/palettes, a palettes.json beside config.json)
// and the user's palettes (this browser's storage); the grouped palette menu;
// the editor of user palettes (colour stops, OKLab, reverse, export/import as
// JSON). The tables and checks are in palettes.js (pure, unit-tested).
//
// Names (ids) given to colormaps.js: 'cmocean:<map>', 'crameri:<map>',
// 'shared:<name>', 'user:<n>'; built-ins keep their names ('depth', …), so a
// choice stored by an earlier version still applies.

import { COLORMAP_NAMES, colormapTable, cssGradient, colormapLabel, hasColormap,
         registerColormap, unregisterColormap, registeredColormaps } from './colormaps.js';
import { LIMITS, stopsTable, reverseStops, sortStops, stopsFromTable, rgbToHex, stopsColorAt,
         parseColormapFile, validatePalette, validatePaletteFile, paletteFileText, uniqueName } from './palettes.js';

const USER_KEY = 'flooda5.palettes.user';
const FAMILY_TITLES = { cmocean: 'cmocean', crameri: 'Crameri (Scientific colour maps)' };

const st = {
  store: null, log: () => {}, onChange: () => {}, onPick: () => {},
  families: {},              // colormaps.json family → {title, tip}
  user: [],                  // [{id, name, stops}]
  nextUser: 1,
  shared: 0,
  editor: null,              // the editor's element, once built
  draft: null,               // {id | null, name, stops}
  target: 'water',
};

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

/**
 * Load the scientific maps, the shared palettes and the user's palettes and
 * register them (colormaps.js). `store`: app.js's localStorage wrapper;
 * `fetchJson(path)`; `onChange()`: palettes were edited (re-apply the styles);
 * `onPick(target, id)`: use palette id for 'water' or 'bed'. Failures are
 * logged and leave the built-in palettes.
 */
export async function initPalettes({ store, fetchJson, log, onChange, onPick }) {
  Object.assign(st, { store, log: log || st.log, onChange: onChange || st.onChange, onPick: onPick || st.onPick });
  loadCss();
  loadUser();
  const [maps, shared] = await Promise.allSettled([
    fetchJson(new URL('../colormaps.json', import.meta.url).pathname),
    fetchJson('/api/palettes'),
  ]);
  if (maps.status === 'fulfilled') {
    for (const [fid, f] of Object.entries(maps.value.families || {})) {
      const lic = `${f.licence || ''}${f.copyright ? `, ${f.copyright}` : ''}`;
      st.families[fid] = { title: FAMILY_TITLES[fid] || String(f.title || fid),
                           tip: `${f.title || fid} ${f.version ? `(${f.version}) ` : ''}— ${lic}; ${f.url || ''}` };
    }
    for (const m of parseColormapFile(maps.value)) registerColormap(m.id, m.table, { label: m.label, group: m.group });
  } else st.log(`Scientific colour maps not loaded: ${maps.reason?.message || maps.reason}`, 'warn');
  if (shared.status === 'fulfilled' && shared.value) {
    const s = shared.value;
    if (s.error) st.log(`Shared palettes: ${s.error}`, 'warn');
    if (Array.isArray(s.palettes) && s.palettes.length) {
      const r = validatePaletteFile({ flooda5_palettes: 1, palettes: s.palettes });
      if (r.ok) {
        for (const p of r.palettes) registerColormap(`shared:${p.name}`, stopsTable(p.stops), { label: p.name, group: 'shared' });
        st.shared = r.palettes.length;
      } else st.log(`Shared palettes ignored: ${r.error}`, 'warn');
    }
  }
  registerUser();
}

/** The user's palettes from the browser's storage (entries that fail the checks are dropped). */
function loadUser() {
  st.user = []; st.nextUser = 1;
  let d = null;
  try { d = JSON.parse(st.store?.get(USER_KEY) || 'null'); } catch { d = null; }
  if (!d || !Array.isArray(d.palettes)) return;
  for (const e of d.palettes) {
    const id = e && typeof e.id === 'string' && /^user:\d{1,9}$/.test(e.id) ? e.id : null;
    const r = id ? validatePalette({ name: e.name, stops: e.stops }) : { ok: false };
    if (r.ok && !st.user.some(u => u.id === id)) st.user.push({ id, ...r.palette });
  }
  const n = Number.isInteger(d.next) ? d.next : 1;
  st.nextUser = Math.max(n, 1 + Math.max(0, ...st.user.map(u => +u.id.slice(5))));
}

function saveUser() {
  st.store?.set(USER_KEY, JSON.stringify({ next: st.nextUser, palettes: st.user }));
}

function registerUser() {
  for (const e of registeredColormaps()) if (e.group === 'user') unregisterColormap(e.id);
  for (const u of st.user) registerColormap(u.id, stopsTable(u.stops), { label: u.name, group: 'user' });
}

// ---------------------------------------------------------------------------
// The palette menu
// ---------------------------------------------------------------------------

/** Groups of the palette menu: [{title, tip, ids}] (built-in first; the bed adds 'terrain'). */
export function paletteGroups(bed = false) {
  const reg = registeredColormaps(), of = g => reg.filter(e => e.group === g).map(e => e.id);
  const groups = [{ title: 'Built-in', ids: COLORMAP_NAMES.filter(n => bed || n !== 'terrain') }];
  for (const [fid, f] of Object.entries(st.families)) groups.push({ title: f.title, tip: f.tip, ids: of(fid) });
  groups.push({ title: 'Shared', tip: 'palettes.json beside config.json (read-only)', ids: of('shared') });
  groups.push({ title: 'Your palettes', tip: 'Kept in this browser', ids: of('user') });
  return groups.filter(g => g.ids.length);
}

/** Is `id` a palette the water (bed = false) or the bed may use? */
export const isPalette = (id, bed = false) => typeof id === 'string' && hasColormap(id) && (bed || id !== 'terrain');

/**
 * Fill the palette menu element: a heading per group, a button per palette
 * (swatch and name, as text), then "Edit palettes…".
 */
export function fillPaletteMenu(menu, { bed = false, current, dflt, onPick, onEdit }) {
  menu.replaceChildren();
  for (const g of paletteGroups(bed)) {
    const h = document.createElement('div');
    h.className = 'palette-group'; h.textContent = g.title; h.setAttribute('role', 'presentation');
    if (g.tip) h.title = g.tip;
    menu.appendChild(h);
    for (const id of g.ids) {
      const b = document.createElement('button');
      b.type = 'button'; b.className = 'palette'; b.setAttribute('role', 'option');
      b.setAttribute('aria-selected', String(id === current));
      b.dataset.palette = id;
      const sw = document.createElement('span'); sw.className = 'swatch'; sw.style.background = cssGradient(id);
      const lab = document.createElement('span');
      lab.textContent = colormapLabel(id) + (id === dflt ? ' (default)' : '');
      b.append(sw, lab);
      b.onclick = () => onPick(id);
      menu.appendChild(b);
    }
  }
  const e = document.createElement('button');
  e.type = 'button'; e.className = 'palette palette-edit'; e.textContent = 'Edit palettes…';
  e.onclick = () => onEdit();
  menu.appendChild(e);
}

// ---------------------------------------------------------------------------
// The editor
// ---------------------------------------------------------------------------

const DEFAULT_STOPS = [{ pos: 0, color: '#e8f4f8' }, { pos: 0.5, color: '#4a90c2' }, { pos: 1, color: '#14285a' }];

const EDITOR_HTML = `
<div class="pe-head">
  <h2 id="pe-title">Palettes</h2>
  <button class="icon-btn" type="button" data-act="close" title="Close (Esc)">✕</button>
</div>
<div class="pe-row">
  <select id="pe-list" aria-label="Your palettes"></select>
  <button class="btn" type="button" data-act="new" title="A new palette">New</button>
  <button class="btn" type="button" data-act="delete" title="Delete this palette">Delete</button>
</div>
<div class="pe-row">
  <label for="pe-name">Name</label>
  <input id="pe-name" type="text" maxlength="${LIMITS.maxName}" autocomplete="off">
</div>
<canvas id="pe-preview" width="256" height="1" aria-label="Preview"></canvas>
<div id="pe-stops" class="pe-stops"></div>
<div class="pe-row">
  <button class="btn" type="button" data-act="add" title="Add a stop in the widest gap">Add stop</button>
  <button class="btn" type="button" data-act="reverse">Reverse</button>
  <button class="btn" type="button" data-act="even" title="Space the stops evenly from 0 to 1">Even</button>
</div>
<div class="pe-row">
  <label for="pe-from">Start from</label>
  <select id="pe-from" aria-label="Start from a palette"></select>
  <button class="btn" type="button" data-act="from" title="Replace the stops by 7 samples of that palette">Copy</button>
</div>
<div class="pe-row">
  <button class="btn" type="button" data-act="save">Save</button>
  <button class="btn" type="button" data-act="water" title="Save, and colour the water with it">Use for water</button>
  <button class="btn" type="button" data-act="bed" title="Save, and colour the bed with it">Use for bed</button>
</div>
<div class="pe-row">
  <button class="btn" type="button" data-act="export" title="Download your palettes as a JSON file">Export…</button>
  <button class="btn" type="button" data-act="import" title="Add palettes from a JSON file (stays in this browser)">Import…</button>
  <input id="pe-file" type="file" accept=".json,application/json" hidden>
</div>
<p id="pe-msg" class="help" role="status"></p>
<p class="help">Colours are interpolated in OKLab. Kept in this browser; Export
  to keep a copy or share them (docs/DATA_FORMATS.md).</p>
`;

const q = id => st.editor.querySelector(`#${id}`);

function loadCss() {
  if (document.getElementById('pe-css')) return;
  const l = document.createElement('link');
  l.id = 'pe-css'; l.rel = 'stylesheet'; l.href = new URL('../css/palettes.css', import.meta.url).href;
  document.head.appendChild(l);
}

function buildEditor() {
  loadCss();
  const el = document.createElement('div');
  el.id = 'palette-editor'; el.className = 'palette-editor';
  el.setAttribute('role', 'dialog'); el.setAttribute('aria-labelledby', 'pe-title');
  el.hidden = true;
  el.innerHTML = EDITOR_HTML;              // fixed markup: names are only ever set as text
  document.body.appendChild(el);
  st.editor = el;
  el.addEventListener('click', e => {
    const b = e.target.closest('button[data-act]');
    if (b) act(b.dataset.act);
    const rm = e.target.closest('button[data-rm]');
    if (rm) removeStop(+rm.dataset.rm);
  });
  el.addEventListener('keydown', e => { if (e.key === 'Escape') closePaletteEditor(); });
  q('pe-list').onchange = () => selectDraft(q('pe-list').value);
  q('pe-name').oninput = () => { st.draft.name = q('pe-name').value; };
  q('pe-stops').addEventListener('input', e => {
    const i = +e.target.dataset.i;
    if (e.target.type === 'color') { st.draft.stops[i].color = e.target.value; preview(); }
  });
  q('pe-stops').addEventListener('change', e => {
    const i = +e.target.dataset.i;
    if (e.target.type === 'number') {
      const v = +e.target.value;
      if (Number.isFinite(v)) st.draft.stops[i].pos = Math.min(Math.max(v, 0), 1);
      st.draft.stops = sortStops(st.draft.stops);
      renderStops();
    }
  });
  q('pe-file').onchange = async () => {
    const f = q('pe-file').files[0];
    q('pe-file').value = '';
    if (!f) return;
    if (f.size > LIMITS.maxBytes) return msg(`${f.name}: larger than 256 KiB`, true);
    msg(importPalettes(await f.text(), f.name).message, false);
  };
}

function msg(text, bad = false) {
  const m = q('pe-msg');
  m.textContent = text; m.classList.toggle('bad', !!bad);
}

/** Open the editor for colour target 'water' or 'bed', on palette `current` if it is the user's. */
export function openPaletteEditor({ target = 'water', current = null } = {}) {
  if (!st.editor) buildEditor();
  st.target = target;
  st.editor.hidden = false;
  fillFrom();
  selectDraft(st.user.some(u => u.id === current) ? current : st.draft?.id ?? st.user[0]?.id ?? '');
  q('pe-name').focus();
}

export function closePaletteEditor() {
  if (st.editor) st.editor.hidden = true;
}

function fillList(sel) {
  const list = q('pe-list');
  list.replaceChildren();
  for (const u of st.user) {
    const o = document.createElement('option'); o.value = u.id; o.textContent = u.name; list.appendChild(o);
  }
  const o = document.createElement('option'); o.value = ''; o.textContent = '(new palette)'; list.appendChild(o);
  list.value = sel ?? '';
}

function fillFrom() {
  const from = q('pe-from');
  from.replaceChildren();
  for (const g of paletteGroups(true)) {
    const og = document.createElement('optgroup'); og.label = g.title;
    for (const id of g.ids) {
      const o = document.createElement('option'); o.value = id; o.textContent = colormapLabel(id); og.appendChild(o);
    }
    from.appendChild(og);
  }
}

function selectDraft(id) {
  const u = st.user.find(x => x.id === id);
  st.draft = u ? { id: u.id, name: u.name, stops: u.stops.map(s => ({ ...s })) }
               : { id: null, name: '', stops: DEFAULT_STOPS.map(s => ({ ...s })) };
  fillList(st.draft.id ?? '');
  q('pe-name').value = st.draft.name;
  renderStops();
  msg('');
}

function renderStops() {
  const box = q('pe-stops'), n = st.draft.stops.length;
  box.replaceChildren();
  st.draft.stops.forEach((s, i) => {
    const row = document.createElement('div'); row.className = 'pe-stop';
    const c = document.createElement('input');
    c.type = 'color'; c.value = s.color; c.dataset.i = String(i); c.setAttribute('aria-label', `Stop ${i + 1} colour`);
    const p = document.createElement('input');
    p.type = 'number'; p.min = '0'; p.max = '1'; p.step = '0.01'; p.value = String(s.pos); p.dataset.i = String(i);
    p.setAttribute('aria-label', `Stop ${i + 1} position (0–1)`);
    const rm = document.createElement('button');
    rm.type = 'button'; rm.className = 'icon-btn'; rm.textContent = '✕'; rm.dataset.rm = String(i);
    rm.title = 'Remove this stop'; rm.disabled = n <= LIMITS.minStops;
    row.append(c, p, rm);
    box.appendChild(row);
  });
  st.editor.querySelector('[data-act="add"]').disabled = n >= LIMITS.maxStops;
  preview();
}

function preview() {
  const cv = q('pe-preview'), ctx = cv.getContext('2d');
  if (!ctx) return;
  const img = ctx.createImageData(256, 1);
  img.data.set(stopsTable(sortStops(st.draft.stops)));
  ctx.putImageData(img, 0, 0);
}

function removeStop(i) {
  if (st.draft.stops.length <= LIMITS.minStops) return;
  st.draft.stops.splice(i, 1);
  renderStops();
}

/** A stop in the middle of the widest gap, with the colour there. */
function addStop() {
  const s = st.draft.stops;
  if (s.length >= LIMITS.maxStops) return;
  let best = 0, gap = -1;
  for (let i = 0; i + 1 < s.length; i++) if (s[i + 1].pos - s[i].pos > gap) { gap = s[i + 1].pos - s[i].pos; best = i; }
  const pos = +((s[best].pos + s[best + 1].pos) / 2).toFixed(4);
  s.splice(best + 1, 0, { pos, color: rgbToHex(stopsColorAt(s, pos)) });
  renderStops();
}

/** Save the draft as a user palette → its id, or null (the message says why). */
function saveDraft() {
  const d = st.draft;
  const r = validatePalette({ name: d.name, stops: sortStops(d.stops) });
  if (!r.ok) { msg(r.error, true); return null; }
  const p = r.palette;
  if (st.user.some(u => u.id !== d.id && u.name.toLowerCase() === p.name.toLowerCase())) {
    msg(`You already have a palette called "${p.name}"`, true); return null;
  }
  if (!d.id && st.user.length >= LIMITS.maxPalettes) { msg(`At most ${LIMITS.maxPalettes} palettes`, true); return null; }
  if (!d.id) d.id = `user:${st.nextUser++}`;
  const i = st.user.findIndex(u => u.id === d.id), entry = { id: d.id, ...p };
  if (i >= 0) st.user[i] = entry; else st.user.push(entry);
  saveUser(); registerUser();
  d.name = p.name; d.stops = p.stops.map(s => ({ ...s }));
  fillList(d.id); fillFrom();
  msg(`Saved "${p.name}"`);
  st.onChange();
  return d.id;
}

function deleteDraft() {
  const id = st.draft.id;
  if (!id) return selectDraft('');
  const name = st.user.find(u => u.id === id)?.name;
  st.user = st.user.filter(u => u.id !== id);
  saveUser(); registerUser();
  selectDraft(st.user[0]?.id ?? '');
  fillFrom();
  msg(`Deleted "${name}"`);
  st.onChange();
}

function act(a) {
  const d = st.draft;
  if (a === 'close') closePaletteEditor();
  else if (a === 'new') selectDraft('');
  else if (a === 'delete') deleteDraft();
  else if (a === 'add') addStop();
  else if (a === 'reverse') { d.stops = reverseStops(sortStops(d.stops)); renderStops(); }
  else if (a === 'even') { d.stops = d.stops.map((s, i, all) => ({ ...s, pos: +(i / (all.length - 1)).toFixed(4) })); renderStops(); }
  else if (a === 'from') {
    const id = q('pe-from').value;
    if (!id) return;
    d.stops = stopsFromTable(colormapTable(id), 7);
    if (!d.name.trim()) { d.name = uniqueName(`${colormapLabel(id)} (edited)`.slice(0, LIMITS.maxName), st.user.map(u => u.name)); q('pe-name').value = d.name; }
    renderStops();
  } else if (a === 'save') saveDraft();
  else if (a === 'water' || a === 'bed') { const id = saveDraft(); if (id) st.onPick(a, id); }
  else if (a === 'export') exportPalettes();
  else if (a === 'import') q('pe-file').click();
}

/** The user's palettes as a palette file (JSON text), or null when there are none. */
export function exportText() {
  return st.user.length ? paletteFileText(st.user) : null;
}

function exportPalettes() {
  const text = exportText();
  if (!text) return msg('No palettes to export: save one first', true);
  const url = URL.createObjectURL(new Blob([text], { type: 'application/json' }));
  const a = document.createElement('a');
  a.href = url; a.download = 'flooda5-palettes.json';
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  msg(`Exported ${st.user.length} palette${st.user.length === 1 ? '' : 's'}`);
}

/**
 * Add the palettes of a palette file (JSON text) to the user's; a name in use
 * gets " (2)", … → {ok, message, ids}. Nothing is added when any check fails.
 */
export function importPalettes(text, label = 'file') {
  const r = validatePaletteFile(text);
  if (!r.ok) return { ok: false, message: `${label}: ${r.error}`, ids: [] };
  if (st.user.length + r.palettes.length > LIMITS.maxPalettes) {
    return { ok: false, message: `${label}: at most ${LIMITS.maxPalettes} palettes in all`, ids: [] };
  }
  const ids = [];
  for (const p of r.palettes) {
    const id = `user:${st.nextUser++}`;
    st.user.push({ id, name: uniqueName(p.name, st.user.map(u => u.name)), stops: p.stops });
    ids.push(id);
  }
  saveUser(); registerUser();
  if (st.editor && !st.editor.hidden) { fillFrom(); selectDraft(ids[0]); }
  st.onChange();
  return { ok: true, message: `Imported ${ids.length} palette${ids.length === 1 ? '' : 's'} from ${label}`, ids };
}

/** A scripting interface (window.__flooda5api.palettes; tools/smoke.mjs). */
export const paletteApi = {
  groups: (bed = false) => paletteGroups(bed).map(g => ({ title: g.title, ids: [...g.ids] })),
  user: () => st.user.map(u => ({ id: u.id, name: u.name, stops: u.stops.map(s => ({ ...s })) })),
  open: (target = 'water', current = null) => openPaletteEditor({ target, current }),
  close: () => closePaletteEditor(),
  /** Set the editor's draft {name, stops} (opens the editor), as typed by the user. */
  edit: ({ name, stops } = {}) => {
    if (!st.editor || st.editor.hidden) openPaletteEditor({ target: st.target });
    if (typeof name === 'string') { st.draft.name = name; q('pe-name').value = name; }
    if (Array.isArray(stops)) st.draft.stops = stops.map(s => ({ pos: +s.pos, color: String(s.color) }));
    renderStops();
    return { name: st.draft.name, stops: st.draft.stops.map(s => ({ ...s })) };
  },
  /** Press an editor button: new, delete, add, reverse, even, from, save, water, bed. */
  press: a => { if (!st.editor || st.editor.hidden) openPaletteEditor({ target: st.target }); act(a); return q('pe-msg').textContent; },
  select: id => { openPaletteEditor({ target: st.target, current: id }); return st.draft.id; },
  draft: () => st.draft && { id: st.draft.id, name: st.draft.name, stops: st.draft.stops.map(s => ({ ...s })) },
  exportText: () => exportText(),
  importText: text => importPalettes(text, 'text'),
  table: id => Array.from(colormapTable(id)),
};
