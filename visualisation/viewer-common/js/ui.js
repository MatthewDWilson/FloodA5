// ui.js — the viewer's controls (the same for every rendering engine).
// app.js inserts this markup; the engine draws into #viewer.

export const UI_HTML = `
<header id="topbar">
  <div class="logo">Flood<span>A5</span></div>
  <div id="status-dot" class="dot"></div>
  <div id="status-label" class="status">starting</div>
  <div id="run-label" class="run" title="Run"></div>
  <div class="spacer"></div>
  <div class="stat"><div class="k">Cells</div><div class="v" id="stat-cells">–</div></div>
  <div class="stat"><div class="k">Frames</div><div class="v" id="stat-frames">0</div></div>
  <div class="stat"><div class="k">Wet cells</div><div class="v accent" id="stat-wet">–</div></div>
  <div class="stat"><div class="k">Max depth</div><div class="v warn" id="stat-max">–</div></div>
  <div class="stat"><div class="k">Sim time</div><div class="v" id="stat-time">–</div></div>
</header>

<aside id="panel">
  <button id="panel-toggle" class="icon-btn" title="Hide/show the panel">☰</button>
  <div id="panel-body">
    <nav class="tabs" role="tablist" aria-label="Panel">
      <button id="tab-btn-layers" class="tab" role="tab" data-tab="layers" aria-selected="true" aria-controls="tab-layers">Layers</button>
      <button id="tab-btn-run" class="tab" role="tab" data-tab="run" aria-selected="false" aria-controls="tab-run">Run</button>
      <button id="tab-btn-about" class="tab" role="tab" data-tab="about" aria-selected="false" aria-controls="tab-about">About</button>
    </nav>
    <div id="tab-layers" class="tab-pane" role="tabpanel" aria-labelledby="tab-btn-layers">
    <section>
      <h2>Variable</h2>
      <select id="var"></select>
      <button id="legend-bar" class="legend" type="button" title="Colour palette: click to change"
              aria-haspopup="listbox" aria-expanded="false"></button>
      <div id="palette-menu" class="palette-menu" role="listbox" aria-label="Colour palette" hidden></div>
      <div class="legend-labels"><span id="legend-min">0</span><span id="legend-unit"></span><span id="legend-max">1</span></div>
      <div class="row">
        <label for="range-mode">Range</label>
        <select id="range-mode" title="run: the 2nd–98th percentiles of the whole run (wet cells); frame: of the frame shown">
          <option value="run">run (2–98 %)</option>
          <option value="frame">frame (1–99 %)</option>
          <option value="fixed">fixed</option>
        </select>
      </div>
      <div class="row" id="range-fixed" hidden>
        <input id="range-min" type="number" step="any" aria-label="Range minimum">
        <span>to</span>
        <input id="range-max" type="number" step="any" aria-label="Range maximum">
      </div>
    </section>

    <section>
      <h2>Display</h2>
      <div class="row">
        <label for="display-mode">Surface</label>
        <select id="display-mode">
          <option value="0">water level</option>
          <option value="1">depth × exaggeration</option>
        </select>
      </div>
      <div class="row" id="exag-row" hidden>
        <label for="exag">Exaggeration</label>
        <input id="exag" type="range" min="1" max="100" step="1" value="10">
        <span id="exag-val" class="val">10×</span>
      </div>
      <div class="row">
        <label for="opacity">Opacity</label>
        <input id="opacity" type="range" min="0.1" max="1" step="0.05" value="0.8">
        <span id="opacity-val" class="val">0.80</span>
      </div>
      <div class="row">
        <label for="min-depth">Wet if depth ≥</label>
        <input id="min-depth" type="number" min="0" step="0.01" value="0.05">
        <span class="val">m</span>
      </div>
      <div class="row">
        <label for="show-borders">Cell borders</label>
        <input id="show-borders" type="checkbox">
      </div>
    </section>

    <section>
      <h2>Bed surface</h2>
      <div class="row">
        <label for="show-bed">Show</label>
        <input id="show-bed" type="checkbox">
      </div>
      <button id="bed-legend" class="legend" type="button" title="Bed colour palette: click to change"
              aria-haspopup="listbox" aria-expanded="false"></button>
      <div id="bed-palette-menu" class="palette-menu" role="listbox" aria-label="Bed colour palette" hidden></div>
      <div class="legend-labels"><span id="bed-min">–</span><span>m</span><span id="bed-max">–</span></div>
      <div class="row">
        <label for="bed-range-mode">Range</label>
        <select id="bed-range-mode" title="auto: the 1st–99th percentiles of the bed elevations">
          <option value="auto">auto (1–99 %)</option>
          <option value="fixed">fixed</option>
        </select>
      </div>
      <div class="row" id="bed-range-fixed" hidden>
        <input id="bed-range-min" type="number" step="any" aria-label="Bed range minimum">
        <span>to</span>
        <input id="bed-range-max" type="number" step="any" aria-label="Bed range maximum">
      </div>
      <div class="row">
        <label for="bed-opacity">Opacity</label>
        <input id="bed-opacity" type="range" min="0.1" max="1" step="0.05" value="1">
        <span id="bed-opacity-val" class="val">1.00</span>
      </div>
    </section>

    <section>
      <h2>Basemap</h2>
      <select id="basemap"></select>
    </section>

    <section>
      <h2>Vertical offset</h2>
      <div class="row">
        <input id="offset" type="number" step="0.1" value="0" aria-label="Vertical offset in metres">
        <span class="val">m</span>
        <button id="calibrate" class="btn" title="Estimate from the basemap's heights under the mesh">Estimate</button>
      </div>
      <p class="help">Added to the model's heights to give the basemap's. For
        a DEM on a mean-sea-level datum and a basemap with ellipsoidal heights
        (3D tiles), it is the geoid height there: between about −107 and +86 m
        worldwide. Each basemap keeps its own; it is estimated from the
        basemap when not given (<code>--offset</code>) — Estimate redoes it.</p>
      <p id="offset-source" class="help"></p>
    </section>

    <section class="buttons">
      <button id="home" class="btn">Home view</button>
      <span id="engine-label" class="engine"></span>
    </section>
    </div>
    <div id="tab-run" class="tab-pane info" role="tabpanel" aria-labelledby="tab-btn-run" hidden></div>
    <div id="tab-about" class="tab-pane info" role="tabpanel" aria-labelledby="tab-btn-about" hidden></div>
  </div>
</aside>

<footer id="playbar">
  <input id="timeline" type="range" min="0" max="0" value="0" step="1" aria-label="Frame">
  <div class="controls">
    <button id="btn-start" class="btn" title="First frame">⏮</button>
    <button id="btn-back" class="btn" title="Previous frame (←)">◀</button>
    <button id="btn-play" class="btn" title="Play/pause (space)">▶</button>
    <button id="btn-fwd" class="btn" title="Next frame (→)">▶|</button>
    <button id="btn-end" class="btn" title="Last frame">⏭</button>
    <select id="speed" title="Playback speed">
      <option value="fps:2">2 frames/s</option>
      <option value="fps:4" selected>4 frames/s</option>
      <option value="fps:10">10 frames/s</option>
      <option value="x:60">1 min/s</option>
      <option value="x:600">10 min/s</option>
      <option value="x:3600">1 h/s</option>
    </select>
    <button id="btn-live" class="btn active" title="Follow the newest frame">● LIVE</button>
    <div id="time-label" class="time">–</div>
  </div>
</footer>

<div id="credits"></div>
<div id="tooltip" hidden></div>
<div id="log"></div>
<div id="fatal" hidden></div>
`;
