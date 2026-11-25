// plot.js — time series at cells (right-click on the map → Add to plot), drawn
// with Chart.js (loaded by app.js when the first plot opens). The pure helpers
// at the top are unit-tested (test/plot.test.js); TimePlot needs a DOM.

import { formatTime } from './protocol.js';

/**
 * Line colours, one per plotted cell: the dataviz reference categorical palette,
 * dark-surface steps (validated on the panel, #0d1318: CVD ΔE ≥ 8.4 between
 * neighbours, all ≥ 3:1). A cell keeps its colour (slot) until removed, and
 * map markers carry the slot number, so identity is never colour alone.
 */
export const SERIES_COLORS = ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#0ca30c', '#9085e9', '#e66767'];
export const MAX_SERIES = SERIES_COLORS.length;

/** The lowest free colour slot (0-based) given the slots in use, or -1 if all are taken. */
export function freeSlot(used) {
  const u = new Set(used);
  for (let k = 0; k < MAX_SERIES; k++) if (!u.has(k)) return k;
  return -1;
}

/**
 * A right *click*, not a right drag (which rotates or zooms the view): the
 * button released within `px` pixels and `ms` milliseconds of being pressed.
 */
export function isClick(down, up, { px = 5, ms = 300 } = {}) {
  if (!down || !up) return false;
  return Math.hypot(up.x - down.x, up.y - down.y) <= px && up.t - down.t <= ms;
}

/**
 * A plotted variable's values from the server's series ({t, values}, values
 * null where missing). `level` is derived from the depth: bed datum + depth
 * where the depth is at least `wet`, else null (dry: no water level).
 */
export function seriesValues(varName, values, { base = 0, wet = 0 } = {}) {
  if (varName !== 'level') return values.map(v => (v === null || !Number.isFinite(v) ? null : v));
  return values.map(d => (d === null || !Number.isFinite(d) || d < wet ? null : base + d));
}

/** The index of the time in sorted `ts` nearest to `t` (-1 if `ts` is empty). */
export function nearestIndex(ts, t) {
  if (!ts.length) return -1;
  let lo = 0, hi = ts.length - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (ts[mid] <= t) lo = mid; else hi = mid;
  }
  return Math.abs(ts[hi] - t) < Math.abs(ts[lo] - t) ? hi : lo;
}

function csvField(s) {
  s = String(s);
  return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

/**
 * CSV of the plotted series: a time column (seconds) and one column per cell,
 * over the union of their times (empty where a series has no value).
 * series: [{label, t: [...], y: [...]}]; header names the variable and unit.
 */
export function seriesCsv(series, { variable = 'value', unit = '' } = {}) {
  const times = Array.from(new Set(series.flatMap(s => s.t))).sort((a, b) => a - b);
  const maps = series.map(s => new Map(s.t.map((t, k) => [t, s.y[k]])));
  const u = unit ? ` (${unit})` : '';
  const head = ['time_s', ...series.map(s => `${variable}${u} ${s.label}`)].map(csvField).join(',');
  const rows = times.map(t => [t, ...maps.map(m => {
    const v = m.get(t);
    return v === null || v === undefined || !Number.isFinite(v) ? '' : v;
  })].join(','));
  return [head, ...rows].join('\n') + '\n';
}

// tick steps of a time axis (seconds): round seconds, minutes, hours, days
const TIME_STEPS = [1, 2, 5, 10, 15, 30, 60, 120, 300, 600, 900, 1800, 3600, 7200, 10800, 21600, 43200,
                    86400, 172800, 432000, 864000];

/** About n ticks at a round time step over [min, max] (s): {step, ticks}. */
export function timeTicks(min, max, n = 6) {
  if (!(max > min)) return { step: 1, ticks: Number.isFinite(min) ? [min] : [] };
  const want = (max - min) / n;
  const step = TIME_STEPS.find(x => x >= want) ?? Math.ceil(want / 864000) * 864000;
  const ticks = [];
  for (let k = Math.ceil(min / step - 1e-9); k * step <= max + 1e-9; k++) ticks.push(k * step + 0);   // (+ 0: no -0)
  return { step, ticks };
}

/** Label of tick t (s) on an axis with the given step: whole days or hours as such, else h:mm:ss. */
export function timeTick(t, step = 1) {
  if (step % 86400 === 0) return `${t / 86400} d`;
  if (step % 3600 === 0) return `${t / 3600} h`;
  if (step % 60 === 0 && Math.abs(t) < 3600) return `${t / 60} min`;
  return formatTime(t);
}

// ---------------------------------------------------------------------------
// The plot panel
// ---------------------------------------------------------------------------

/**
 * A Chart.js line chart of cell series against time, with a dashed line at the
 * frame shown and a click handler that reports the time clicked.
 *   new TimePlot(canvas, Chart, {onPick: t => …})
 *   plot.setSeries([{key, label, color, t, y}], {label, unit})
 *   plot.setCursor(t)
 */
export class TimePlot {
  constructor(canvas, Chart, { onPick = null, text = '#d2e3ec', dim = '#6b8293', grid = '#1f2d39' } = {}) {
    this.cursor = null;
    const plot = this;
    const cursorPlugin = {
      id: 'flooda5Cursor',
      afterDatasetsDraw(chart) {
        const t = plot.cursor, x = chart.scales.x;
        if (t === null || !x || t < x.min || t > x.max) return;
        const px = x.getPixelForValue(t), { top, bottom } = chart.chartArea, ctx = chart.ctx;
        ctx.save();
        ctx.strokeStyle = text; ctx.globalAlpha = 0.6; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
        ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, bottom); ctx.stroke();
        ctx.restore();
      },
    };
    this.chart = new Chart(canvas, {
      type: 'line',
      data: { datasets: [] },
      plugins: [cursorPlugin],
      options: {
        animation: false, responsive: true, maintainAspectRatio: false, parsing: false, normalized: true,
        spanGaps: false,
        interaction: { mode: 'nearest', axis: 'x', intersect: false },
        elements: { point: { radius: 0, hoverRadius: 4, hitRadius: 6 }, line: { borderWidth: 2, tension: 0 } },
        scales: {
          x: { type: 'linear', bounds: 'data', grid: { color: grid }, border: { color: grid },
               afterBuildTicks: axis => {
                 const { step, ticks } = timeTicks(axis.min, axis.max, Math.max(2, Math.floor((axis.width || 600) / 90)));
                 axis.flooda5Step = step;
                 axis.ticks = ticks.map(value => ({ value }));
               },
               ticks: { color: dim, autoSkip: false, maxRotation: 0,
                        callback(v) { return timeTick(v, this.flooda5Step); } } },
          y: { ticks: { color: dim, maxTicksLimit: 6 }, grid: { color: grid }, border: { color: grid },
               title: { display: true, color: dim, text: '' } },
        },
        plugins: {
          legend: { display: false },          // the panel's own legend (removable entries)
          decimation: { enabled: true, algorithm: 'min-max' },
          tooltip: {
            backgroundColor: 'rgba(13,19,24,.96)', borderColor: grid, borderWidth: 1,
            titleColor: text, bodyColor: text, titleFont: { family: 'monospace' }, bodyFont: { family: 'monospace' },
            callbacks: {
              title: items => (items.length ? `t = ${formatTime(items[0].parsed.x)}` : ''),
              label: item => `${item.dataset.label}: ${Number.isFinite(item.parsed.y) ? item.parsed.y.toFixed(3) : '–'}`,
            },
          },
        },
        onClick: (ev, _els, chart) => {
          if (!onPick) return;
          const x = chart.scales.x.getValueForPixel(ev.x);
          if (Number.isFinite(x)) onPick(x);
        },
      },
    });
  }

  /** Replace the lines: [{key, label, color, t, y}] (y null: gap). */
  setSeries(series, { label = '', unit = '' } = {}) {
    const c = this.chart;
    c.data.datasets = series.map(s => ({
      label: s.label, borderColor: s.color, backgroundColor: s.color,
      data: s.t.map((t, k) => ({ x: t, y: s.y[k] })),
    }));
    c.options.scales.y.title.text = unit ? `${label} (${unit})` : label;
    c.update('none');
  }

  /** Move the dashed line to time t (null: none). */
  setCursor(t) {
    if (t === this.cursor) return;
    this.cursor = t;
    this.chart.draw();
  }

  resize() { this.chart.resize(); }
  destroy() { this.chart.destroy(); }
}
