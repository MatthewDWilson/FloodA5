import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import {
  hexToRgb, rgbToHex, srgbToOklab, oklabToSrgb, stopsColorAt, stopsTable, reverseStops, sortStops,
  stopsFromTable, decodeHexTable, parseColormapFile, validatePalette, validatePaletteFile, paletteFileText,
  uniqueName, cleanName, LIMITS,
} from '../js/palettes.js';
import {
  registerColormap, unregisterColormap, hasColormap, colormapKey, colormapLabel, colormapTable,
  colorAt, cssGradient, registeredColormaps,
} from '../js/colormaps.js';

const COLORMAPS = JSON.parse(readFileSync(new URL('../colormaps.json', import.meta.url), 'utf8'));

test('hex colours', () => {
  assert.deepEqual(hexToRgb('#ff8000'), [255, 128, 0]);
  assert.deepEqual(hexToRgb('#FF8000'), [255, 128, 0]);
  for (const bad of ['ff8000', '#ff800', '#ff80000', '#gg8000', 12, null, 'red']) assert.equal(hexToRgb(bad), null);
  assert.equal(rgbToHex([255, 128, 0]), '#ff8000');
  assert.equal(rgbToHex([300, -4, 15.6]), '#ff0010');
});

test('OKLab: reference values and round trip', () => {
  const near = (a, b, tol) => a.forEach((x, k) => assert.ok(Math.abs(x - b[k]) < tol, `${a} vs ${b}`));
  near(srgbToOklab([255, 255, 255]), [1, 0, 0], 1e-4);
  near(srgbToOklab([0, 0, 0]), [0, 0, 0], 1e-6);
  near(srgbToOklab([255, 0, 0]), [0.62796, 0.22486, 0.12585], 1e-3);      // Ottosson's sRGB red
  near(srgbToOklab([0, 0, 255]), [0.45201, -0.03246, -0.31153], 1e-3);
  for (let r = 0; r < 256; r += 15) for (let g = 0; g < 256; g += 17) for (let b = 0; b < 256; b += 51) {
    const back = oklabToSrgb(srgbToOklab([r, g, b]));
    near(back, [r, g, b], 1e-6 * 255 + 0.01);
  }
});

test('stops: interpolation in OKLab, ends, steps', () => {
  const stops = [{ pos: 0, color: '#000000' }, { pos: 1, color: '#ffffff' }];
  assert.deepEqual(stopsColorAt(stops, 0), [0, 0, 0]);
  assert.deepEqual(stopsColorAt(stops, 1), [255, 255, 255]);
  assert.deepEqual(stopsColorAt(stops, -3), [0, 0, 0]);
  assert.deepEqual(stopsColorAt(stops, NaN), [0, 0, 0]);
  // OKLab midpoint of black and white: L = 0.5 → sRGB ≈ 99 (RGB would give 128)
  const mid = stopsColorAt(stops, 0.5);
  assert.ok(mid[0] === mid[1] && mid[1] === mid[2] && Math.abs(mid[0] - 99) <= 1, `${mid}`);
  // inner stops: ends hold before the first / after the last
  const inner = [{ pos: 0.2, color: '#ff0000' }, { pos: 0.8, color: '#0000ff' }];
  assert.deepEqual(stopsColorAt(inner, 0.1), [255, 0, 0]);
  assert.deepEqual(stopsColorAt(inner, 0.9), [0, 0, 255]);
  // equal positions: a hard step
  const step = [{ pos: 0, color: '#ff0000' }, { pos: 0.5, color: '#ff0000' }, { pos: 0.5, color: '#00ff00' }, { pos: 1, color: '#00ff00' }];
  assert.deepEqual(stopsColorAt(step, 0.49), [255, 0, 0]);
  assert.deepEqual(stopsColorAt(step, 0.51), [0, 255, 0]);
  const t = stopsTable(stops);
  assert.equal(t.length, 4 * 256);
  assert.deepEqual([...t.slice(0, 4)], [0, 0, 0, 255]);
  assert.deepEqual([...t.slice(1020)], [255, 255, 255, 255]);
  for (let i = 1; i < 256; i++) assert.ok(t[4 * i] >= t[4 * i - 4]);     // monotone grey ramp
});

test('stops: reverse, sort, sample', () => {
  const s = [{ pos: 0, color: '#112233' }, { pos: 0.25, color: '#445566' }, { pos: 1, color: '#778899' }];
  assert.deepEqual(reverseStops(s), [{ pos: 0, color: '#778899' }, { pos: 0.75, color: '#445566' }, { pos: 1, color: '#112233' }]);
  assert.deepEqual(reverseStops(reverseStops(s)), s);
  const t1 = stopsTable(s), t2 = stopsTable(reverseStops(s));
  for (let i = 0; i < 256; i++) for (let k = 0; k < 3; k++) assert.ok(Math.abs(t1[4 * i + k] - t2[4 * (255 - i) + k]) <= 1);
  assert.deepEqual(sortStops([{ pos: 0.5, color: '#000001' }, { pos: 0, color: '#000002' }, { pos: 0.5, color: '#000003' }]).map(x => x.color),
                   ['#000002', '#000001', '#000003']);
  const k5 = stopsFromTable(t1, 5);
  assert.equal(k5.length, 5);
  assert.deepEqual(k5[0], { pos: 0, color: '#112233' });
  assert.deepEqual(k5[4], { pos: 1, color: '#778899' });
});

test('colormaps.json: the shipped scientific maps', () => {
  const fams = COLORMAPS.families;
  assert.deepEqual(Object.keys(fams.cmocean.maps),
                   ['deep', 'dense', 'turbid', 'haline', 'speed', 'amp', 'tempo', 'thermal', 'ice', 'rain']);
  assert.deepEqual(Object.keys(fams.crameri.maps),
                   ['batlow', 'oslo', 'devon', 'lajolla', 'davos', 'imola', 'lapaz', 'hawaii', 'roma', 'vik']);
  for (const f of Object.values(fams)) {
    assert.equal(f.licence, 'MIT');
    assert.match(f.licence_text, /Permission is hereby granted, free of charge/);
    assert.match(f.copyright, /^Copyright \(c\)/);
  }
  const maps = parseColormapFile(COLORMAPS);
  assert.equal(maps.length, 20);
  for (const m of maps) assert.equal(m.table.length, 4 * 256);
  // known entries (ColorSchemes.jl data/scicolor.jl, data/cmocean.jl): batlow
  // starts at RGB(0.0052, 0.0982, 0.3498); deep runs from RGB(0.9928, 0.9944, 0.8001)
  // to RGB(0.1561, 0.1026, 0.1727)
  const batlow = maps.find(m => m.id === 'crameri:batlow'), deep = maps.find(m => m.id === 'cmocean:deep');
  assert.deepEqual([...batlow.table.slice(0, 4)], [1, 25, 89, 255]);
  assert.deepEqual([...deep.table.slice(0, 4)], [253, 254, 204, 255]);
  assert.deepEqual([...deep.table.slice(1020)], [40, 26, 44, 255]);
  assert.equal(deep.group, 'cmocean');
  assert.equal(deep.label, 'Deep');
  assert.equal(decodeHexTable('00ff'), null);
  assert.equal(decodeHexTable('zz'.repeat(768)), null);
  assert.deepEqual(parseColormapFile(null), []);
  assert.deepEqual(parseColormapFile({ families: { bad: { maps: { x: '0' } } } }), []);
});

test('registered palettes: tables, keys, labels', () => {
  const tab = stopsTable([{ pos: 0, color: '#ff0000' }, { pos: 1, color: '#0000ff' }]);
  assert.ok(registerColormap('user:t1', tab, { label: 'Test', group: 'user' }));
  assert.ok(hasColormap('user:t1') && hasColormap('viridis') && !hasColormap('nonesuch') && !hasColormap('toString'));
  assert.equal(colormapLabel('user:t1'), 'Test');
  assert.equal(colormapLabel('depth'), 'Water (blue)');
  const k1 = colormapKey('user:t1');
  assert.notEqual(k1, 'user:t1');
  assert.deepEqual(colormapTable(k1), tab);
  assert.deepEqual(colormapTable('user:t1'), tab);
  assert.deepEqual(colorAt('user:t1', 0), [255, 0, 0]);
  assert.deepEqual(colorAt('user:t1', 1), [0, 0, 255]);
  assert.match(cssGradient('user:t1'), /^linear-gradient\(90deg, rgb\(255,0,0\) 0%/);
  // re-registering (an edited palette) gives a new key, so a surface reloads its table
  const tab2 = stopsTable([{ pos: 0, color: '#00ff00' }, { pos: 1, color: '#0000ff' }]);
  registerColormap('user:t1', tab2, { label: 'Test', group: 'user' });
  assert.notEqual(colormapKey('user:t1'), k1);
  assert.deepEqual(colormapTable(colormapKey('user:t1')), tab2);
  assert.ok(registeredColormaps().some(e => e.id === 'user:t1' && e.group === 'user'));
  assert.ok(!registerColormap('viridis', tab));                  // built-ins are not replaced
  assert.ok(!registerColormap('user:bad', new Uint8Array(12)));
  unregisterColormap('user:t1');
  assert.ok(!hasColormap('user:t1'));
  assert.equal(colormapKey('depth'), 'depth');
});

const good = { name: 'My depth', stops: [{ pos: 0, color: '#C6EBFF' }, { pos: 0.5, color: '#2882dc' }, { pos: 1, color: '#0a1964' }] };
const file = ps => ({ flooda5_palettes: 1, palettes: ps });

test('validation: accepts', () => {
  const r = validatePalette(good);
  assert.ok(r.ok);
  assert.equal(r.palette.stops[0].color, '#c6ebff');                  // lower case
  assert.equal(validatePalette({ ...good, name: '  padded  ' }).palette.name, 'padded');
  const f = validatePaletteFile(JSON.stringify(file([good, { ...good, name: 'Other' }])));
  assert.ok(f.ok && f.palettes.length === 2);
  // equal positions (a hard step) are sorted
  assert.ok(validatePalette({ name: 's', stops: [{ pos: 0, color: '#000000' }, { pos: 0.5, color: '#000000' }, { pos: 0.5, color: '#ffffff' }] }).ok);
  // a written file reads back
  const back = validatePaletteFile(paletteFileText([good]));
  assert.ok(back.ok);
  assert.deepEqual(back.palettes[0], validatePalette(good).palette);
});

test('validation: rejects', () => {
  const stops = good.stops;
  const bad = [
    null, [], 'x', { name: 'x' }, { stops },
    { ...good, extra: 1 },                                           // unknown key
    { ...good, name: '' }, { ...good, name: '   ' }, { ...good, name: 'x'.repeat(LIMITS.maxName + 1) },
    { ...good, name: 'a\u0007b' }, { ...good, name: 42 },
    { name: 'x', stops: [stops[0]] },                                // < 2 stops
    { name: 'x', stops: Array.from({ length: 33 }, (_, i) => ({ pos: i / 32, color: '#000000' })) },
    { name: 'x', stops: [stops[2], stops[0]] },                      // unsorted
    { name: 'x', stops: [{ pos: -0.1, color: '#000000' }, stops[2]] },
    { name: 'x', stops: [stops[0], { pos: 1.5, color: '#000000' }] },
    { name: 'x', stops: [stops[0], { pos: '1', color: '#000000' }] },
    { name: 'x', stops: [stops[0], { pos: NaN, color: '#000000' }] },
    { name: 'x', stops: [stops[0], { pos: 1, color: 'blue' }] },
    { name: 'x', stops: [stops[0], { pos: 1, color: '#fff' }] },
    { name: 'x', stops: [stops[0], { pos: 1, color: '#000000', alpha: 1 }] },
    { name: 'x', stops: [stops[0], [1, '#000000']] },
    { name: '<img src=x onerror=alert(1)>', stops: 'nope' },
  ];
  for (const p of bad) assert.equal(validatePalette(p).ok, false, JSON.stringify(p));
  const badFiles = [
    'not json', '[]', '{}', JSON.stringify({ palettes: [good] }),               // no format key
    JSON.stringify({ flooda5_palettes: 2, palettes: [good] }),                 // other version
    JSON.stringify({ flooda5_palettes: 1, palettes: [good], x: 1 }),           // unknown key
    JSON.stringify(file([])),                                                  // no palettes
    JSON.stringify(file(Array.from({ length: 65 }, (_, i) => ({ ...good, name: `p${i}` })))),
    JSON.stringify(file([good, { ...good, name: 'MY DEPTH' }])),               // duplicate name
    JSON.stringify(file([good, { name: 'x', stops: [] }])),                    // one bad palette rejects the file
    ' '.repeat(LIMITS.maxBytes + 1),
  ];
  for (const f of badFiles) assert.equal(validatePaletteFile(f).ok, false, f.slice(0, 80));
  const r = validatePaletteFile(JSON.stringify(file([{ name: 'x', stops: [] }])));
  assert.match(r.error, /stops/);
});

test('names', () => {
  assert.equal(cleanName(' a '), 'a');
  assert.equal(cleanName(''), null);
  assert.equal(uniqueName('Depth', ['other']), 'Depth');
  assert.equal(uniqueName('Depth', ['depth']), 'Depth (2)');
  assert.equal(uniqueName('Depth', ['Depth', 'Depth (2)']), 'Depth (3)');
  const long = 'x'.repeat(LIMITS.maxName);
  const u = uniqueName(long, [long]);
  assert.ok(u.length <= LIMITS.maxName && u.endsWith(' (2)'));
});
