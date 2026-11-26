import { test } from 'node:test';
import assert from 'node:assert/strict';
import { esc, shellLine, duration, area, bytes, renderRun, renderAbout } from '../js/info.js';

test('escaping and formatting', () => {
  assert.equal(esc(`<a href="x">'&'</a>`), '&lt;a href=&quot;x&quot;&gt;&#39;&amp;&#39;&lt;/a&gt;');
  assert.equal(esc(null), '');
  assert.equal(shellLine(['--dx', '25', '--output', 'out dir/r.h5', "it's"]),
               `--dx 25 --output 'out dir/r.h5' 'it'\\''s'`);
  assert.equal(shellLine(undefined), '');
  assert.equal(duration(45), '45 s');
  assert.equal(duration(750), '12 min 30 s');
  assert.equal(duration(3900), '1 h 05 min');
  assert.equal(duration(NaN), '–');
  assert.equal(area(25), '25.0 m²');
  assert.equal(area(2.5e6), '2.50 km²');
  assert.equal(bytes(512), '512 B');
  assert.equal(bytes(2048), '2.0 kB');
  assert.equal(bytes(-1), '–');
});

test('Run tab: the record, escaped; empty record', () => {
  assert.match(renderRun({}), /no run record/);
  const html = renderRun({
    label: 'r<1>.h5', started_utc: '2026-11-25T01:00:00Z', flooda5_version: '0.8.7', git_commit: '0123456789abcdef', git_dirty: true, git_diff_lines: 12,
    flow_model: 'standard', backend: 'cpu', precision: 'float64', threads: 8,
    argv: ['--mesh', 'a b.parquet', '--sim-duration', '600'],
    mesh: { n_cells: 12345, mesh_family: 'a5', dggs_level: 14, cell_area_mean_m2: 30, cell_area_min_m2: 28, cell_area_max_m2: 33 },
    frames: { n: 11, t_first: 0, t_last: 600 },
    inputs: [{ role: 'inflow hydrograph', original_path: '/x/q.csv', size: 2048, sha256: 'abcdef0123456789', embedded: true }],
    options: { manning_n: 0.035, dt_max: 'Inf', injection_points: [], dem: null },
    summary: { completed: true, n_steps: 1234, finished_utc: '2026-11-25T02:00:00Z' },
  });
  assert.match(html, /r&lt;1&gt;\.h5/);                     // escaped
  assert.ok(!html.includes('<1>'));
  assert.match(html, /0\.8\.7 · git 0123456789 \+ uncommitted changes \(12 diff lines\)/);
  assert.match(html, /cpu · float64 · 8 threads/);
  assert.match(html, /Started<\/dt><dd>2026-11-25T01:00:00Z<\/dd><dt>Finished<\/dt><dd>2026-11-25T02:00:00Z/);
  assert.match(html, /julia FloodModel\.jl --mesh &#39;a b\.parquet&#39; --sim-duration 600/);
  assert.match(html, /12,345|12 345|12 345/);
  assert.match(html, /30\.0 m² mean \(28\.0 m² – 33\.0 m²\)/);
  assert.match(html, /0 s – 10 min/);
  assert.match(html, /2\.0 kB · sha256 abcdef012345… · embedded/);
  assert.match(html, /Options \(2\)/);                       // empty list and null left out
  assert.match(html, /manning_n<\/dt><dd>0\.035/);
  assert.match(html, /n_steps<\/dt><dd>1234/);
  assert.ok(!html.includes('Statistics over the run'));
  // the model's statistics (0.8.13)
  const h2 = renderRun({ label: 'x', statistics: { cadence: 'step', interval_s: 0, flood_depth: 0.05, n_states: 12345 } });
  assert.match(h2, /by the model, every time step \((12,345|12 345|12 345) states\)/);
  assert.match(h2, /depth ≥ 0\.05 m/);
  assert.match(renderRun({ label: 'x', statistics: { cadence: 'interval', interval_s: 30, flood_depth: 0.1, n_states: 5 } }), /every 30 s/);
});

test('About tab: variables in order, notes, topics', () => {
  const docs = {
    variables: {
      velocity: { label: 'Speed', unit: 'm/s', definition: 'Magnitude.', notes: ['Zero below 1 cm.'] },
      water_depth: { label: 'Water depth', unit: 'm', definition: 'Depth <above> datum.', viewer: 'Not drawn when dry.' },
      zeta: { label: 'Zeta', definition: 'Last.' },
    },
    topics: { colour_range: { title: 'Colour range', text: '2nd to 98th percentiles.' } },
  };
  const html = renderAbout(docs, { engine: 'CesiumJS 1.147' });
  assert.ok(html.indexOf('Water depth') < html.indexOf('Speed') && html.indexOf('Speed') < html.indexOf('Zeta'));
  assert.match(html, /Depth &lt;above&gt; datum\./);
  assert.match(html, /<li>Not drawn when dry\.<\/li>/);
  assert.match(html, /<h3>Colour range<\/h3><p>2nd to 98th percentiles\.<\/p>/);
  assert.match(html, /FloodA5 globe viewer · CesiumJS 1\.147/);
  assert.match(renderAbout(null), /No descriptions/);
});
