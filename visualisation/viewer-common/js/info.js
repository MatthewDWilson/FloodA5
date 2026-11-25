// info.js — the side panel's Run and About tabs: HTML from the run record
// (/api/run, FloodA5Output.run_record) and the variable descriptions
// (/api/about, io/variables.toml). Pure: no DOM access, testable in node.

/** HTML-escape a value (records hold paths and command lines). */
export function esc(x) {
  return String(x ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
}

/** The command line, quoted for a POSIX shell where needed. */
export function shellLine(argv) {
  return (argv || []).map(a => (/^[\w@%+=:,./-]+$/.test(a) ? a : `'${String(a).replace(/'/g, `'\\''`)}'`)).join(' ');
}

/** Seconds as "1 h 05 min", "12 min 30 s", "45 s". */
export function duration(s) {
  if (!Number.isFinite(s)) return '–';
  s = Math.round(s);
  const h = Math.floor(s / 3600), m = Math.floor((s % 3600) / 60), r = s % 60;
  if (h) return `${h} h ${String(m).padStart(2, '0')} min`;
  if (m) return `${m} min${r ? ` ${r} s` : ''}`;
  return `${r} s`;
}

/** An area in m², in m² or km². */
export function area(a) {
  if (!Number.isFinite(a)) return '–';
  return a >= 1e5 ? `${(a / 1e6).toPrecision(3)} km²` : `${a.toPrecision(3)} m²`;
}

/** A byte count. */
export function bytes(n) {
  if (!Number.isFinite(n) || n < 0) return '–';
  const u = ['B', 'kB', 'MB', 'GB', 'TB'];
  let k = 0;
  while (n >= 1024 && k < u.length - 1) { n /= 1024; k++; }
  return `${k ? n.toFixed(1) : n} ${u[k]}`;
}

function value(v) {
  if (v === null || v === undefined || v === '') return '–';
  if (typeof v === 'boolean') return v ? 'yes' : 'no';
  if (Array.isArray(v)) return v.length ? v.map(value).join(', ') : '–';
  if (typeof v === 'object') return JSON.stringify(v);
  if (typeof v === 'number' && !Number.isInteger(v)) return String(+v.toPrecision(6));
  return String(v);
}

function dl(rows) {
  const r = rows.filter(([, v]) => v !== undefined && v !== null && v !== '');
  if (!r.length) return '';
  return `<dl>${r.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${esc(v)}</dd>`).join('')}</dl>`;
}

const MESH_LABELS = {
  mesh_family: 'Family', id_scheme: 'Cell ids', dggs_library: 'DGGS library', dggs_system: 'DGGS',
  dggs_level: 'Level', dggs_refinement_ratio: 'Refinement ratio', dggs_equal_area: 'Equal area',
  dggs_interior_valence: 'Interior valence', target_dx_m: 'Target spacing (m)',
};

/** The Run tab. `rec`: /api/run ({} when the server has none). */
export function renderRun(rec) {
  if (!rec || !Object.keys(rec).length) return '<p class="muted">The server has no run record (yet).</p>';
  const out = [];
  const commit = rec.git_commit && rec.git_commit !== 'unknown' ? rec.git_commit.slice(0, 10) : rec.git_commit;
  const dirty = rec.git_dirty ? ` + uncommitted changes${rec.git_diff_lines ? ` (${rec.git_diff_lines} diff lines)` : ''}` : '';
  out.push('<h3>Run</h3>', dl([
    ['Label', rec.label],
    ['File', rec.file || rec.output],
    ['FloodA5', [rec.flooda5_version, commit && `git ${commit}${dirty}`].filter(Boolean).join(' · ')],
    ['Flow model', rec.flow_model],
    ['Backend', [rec.backend, rec.precision, rec.threads && `${rec.threads} threads`].filter(Boolean).join(' · ')],
    ['Cell order', rec.cell_order],
    ['Started', rec.started_utc],
    ['Finished', rec.summary?.finished_utc],
    ['Host', rec.hostname],
    ['Julia', rec.julia_version],
    ['Output format', rec.format_version],
  ]));
  if (rec.argv && rec.argv.length) {
    out.push('<h3>Command line <button class="btn copy" type="button" data-copy="argv" title="Copy to the clipboard">copy</button></h3>',
             `<pre id="run-argv">julia FloodModel.jl ${esc(shellLine(rec.argv))}</pre>`);
  }
  const m = rec.mesh || {};
  const meshRows = [['Cells', Number.isFinite(m.n_cells) ? m.n_cells.toLocaleString() : undefined]];
  for (const [k, lab] of Object.entries(MESH_LABELS)) if (m[k] !== undefined) meshRows.push([lab, value(m[k])]);
  if (Number.isFinite(m.cell_area_mean_m2)) {
    meshRows.push(['Cell area', `${area(m.cell_area_mean_m2)} mean (${area(m.cell_area_min_m2)} – ${area(m.cell_area_max_m2)})`]);
  }
  const ms = rec.mesh_source || {};
  if (ms.original_path) meshRows.push(['Mesh file', ms.original_path + (ms.embedded ? ` (${ms.embedded} in the output)` : '')]);
  out.push('<h3>Mesh</h3>', dl(meshRows));
  const f = rec.frames;
  if (f) out.push('<h3>Frames</h3>', dl([['Frames', f.n], ['Time', f.n ? `${duration(f.t_first)} – ${duration(f.t_last)}` : '–']]));
  const inputs = rec.inputs || [];
  if (inputs.length) {
    out.push('<h3>Inputs</h3>', dl(inputs.map(i => [i.role || i.flag, [
      i.original_path || i.argv_value,
      Number.isFinite(i.size) && i.size >= 0 ? bytes(i.size) : null,
      i.sha256 ? `sha256 ${String(i.sha256).slice(0, 12)}…` : null,
      i.embedded ? 'embedded' : null,
    ].filter(Boolean).join(' · ')])));
  }
  const o = rec.options || {};
  const keys = Object.keys(o).filter(k => !(o[k] === null || (Array.isArray(o[k]) && !o[k].length))).sort();
  if (keys.length) {
    out.push(`<details><summary>Options (${keys.length})</summary>${dl(keys.map(k => [k, value(o[k])]))}</details>`);
  }
  const sm = rec.summary || {};
  if (Object.keys(sm).length) {
    out.push(`<details><summary>Summary at the end of the run</summary>${dl(Object.keys(sm).sort().map(k => [k, value(sm[k])]))}</details>`);
  }
  return out.join('\n');
}

/** Order of the variables in the About tab (others follow). */
const ABOUT_ORDER = ['water_depth', 'level', 'velocity', 'saturation', 'elevation'];

/** The About tab. `docs`: /api/about; `ctx`: {engine}. (The run's FloodA5 version is in the Run tab.) */
export function renderAbout(docs, ctx = {}) {
  const out = [];
  out.push(`<p>FloodA5 globe viewer${ctx.engine ? ` · ${esc(ctx.engine)}` : ''}.</p>`);
  const vars = docs?.variables || {};
  const names = [...ABOUT_ORDER.filter(k => vars[k]), ...Object.keys(vars).filter(k => !ABOUT_ORDER.includes(k)).sort()];
  if (names.length) out.push('<h3>Variables</h3>');
  for (const k of names) {
    const v = vars[k];
    out.push(`<div class="var"><b>${esc(v.label || k)}</b>${v.unit ? ` <span class="unit">(${esc(v.unit)})</span>` : ''}`,
             `<p>${esc(v.definition)}</p>`);
    const notes = [...(v.notes || []), ...(v.viewer ? [v.viewer] : [])];
    if (notes.length) out.push(`<ul>${notes.map(n => `<li>${esc(n)}</li>`).join('')}</ul>`);
    if (v.references?.length) out.push(`<ul>${v.references.map(r => `<li>${esc(r)}</li>`).join('')}</ul>`);
    out.push('</div>');
  }
  const topics = docs?.topics || {};
  for (const t of Object.values(topics)) out.push(`<h3>${esc(t.title)}</h3><p>${esc(t.text)}</p>`);
  if (!names.length && !Object.keys(topics).length) out.push('<p class="muted">No descriptions from the server.</p>');
  return out.join('\n');
}
