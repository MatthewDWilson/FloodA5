// classes.js and hazard.js: node --test test/*.test.js
import test from 'node:test';
import assert from 'node:assert/strict';
import { CLASS_SETS, classify, classOf, classLabels, classColors, classTable, classedName, parseClassedName,
         parseBreaks, classGradient } from '../js/classes.js';
import { colormapTable, cssGradient } from '../js/colormaps.js';
import { DEFAULT_HAZARD, normalizeHazard, hazardQuery, describeHazard, PERSON_PRESETS } from '../js/hazard.js';

test('classify: breaks are lower bounds; NaN stays NaN', () => {
  const b = [0.75, 1.25, 2.0];
  const c = classify(new Float32Array([0, 0.74, 0.75, 1.3, 2.0, 9, NaN]), b);
  assert.deepEqual([...c].slice(0, 6), [0, 0, 1, 2, 3, 3]);
  assert.ok(Number.isNaN(c[6]));
  assert.equal(classOf(1.0, b), 1);
  assert.equal(classOf(NaN, b), -1);
  // ARR classes 1–6 (and 0 = dry) onto H1–H6
  assert.deepEqual([...classify([1, 2, 3, 4, 5, 6], CLASS_SETS.hazard_class_arr.breaks)], [0, 1, 2, 3, 4, 5]);
});

test('class labels: named for the default count, ranges otherwise', () => {
  const hr = classLabels('hazard_rating', CLASS_SETS.hazard_rating.breaks);
  assert.equal(hr.length, 4);
  assert.match(hr[0], /^Low – caution \(< 0\.75\)$/);
  assert.match(hr[3], /Danger for all \(≥ 2\)/);
  assert.deepEqual(classLabels('depth', [0.15, 0.3]), ['< 0.15 m', '0.15 – 0.3 m', '≥ 0.3 m']);
  assert.equal(classLabels('hazard_class_arr', CLASS_SETS.hazard_class_arr.breaks)[0].slice(0, 2), 'H1');
});

test('class colours and the stepped table', () => {
  assert.equal(classColors('hazard4', 4).length, 4);
  const c7 = classColors('ylorrd', 7);
  assert.equal(c7.length, 7);
  assert.deepEqual(c7[0], [255, 237, 160]);
  assert.deepEqual(c7[6], [189, 0, 38]);
  for (const k of [2, 4, 6, 9]) {
    const t = classTable('depth', k), cols = classColors('blues', k);
    // the shader samples class i at texel round(255 (i + 0.5) / k)
    for (let i = 0; i < k; i++) {
      const x = Math.round(255 * (i + 0.5) / k);
      assert.deepEqual([t[4 * x], t[4 * x + 1], t[4 * x + 2]], cols[i], `k=${k} class ${i}`);
    }
  }
  assert.deepEqual(parseClassedName(classedName('hazard_rating', 4)), { set: 'hazard_rating', k: 4 });
  assert.equal(parseClassedName('cls:nope:3'), null);
  assert.equal(parseClassedName('viridis'), null);
  assert.deepEqual(colormapTable('cls:hazard_rating:4'), classTable('hazard_rating', 4));
  assert.equal(cssGradient('cls:depth:5'), classGradient('depth', 5));
});

test('parseBreaks', () => {
  assert.deepEqual(parseBreaks('0.75, 1.25, 2'), [0.75, 1.25, 2]);
  assert.deepEqual(parseBreaks(' 0.1 0.2;0.4 '), [0.1, 0.2, 0.4]);
  assert.equal(parseBreaks('1, 1'), null);
  assert.equal(parseBreaks('2, 1'), null);
  assert.equal(parseBreaks('a'), null);
  assert.equal(parseBreaks(''), null);
});

test('hazard parameters: defaults, presets, clamping, query', () => {
  assert.deepEqual(normalizeHazard({}), { ...DEFAULT_HAZARD });
  const child = normalizeHazard({ person: 'child', mass: 99 });
  assert.equal(child.mass, PERSON_PRESETS.child.mass);
  assert.equal(child.height, PERSON_PRESETS.child.height);
  const own = normalizeHazard({ person: 'custom', mass: 50, height: 1.6, theta: 200, landuse: 'moon' });
  assert.deepEqual([own.mass, own.height, own.theta, own.landuse], [50, 1.6, 90, 'fd2320']);
  assert.equal(normalizeHazard({ mass: 'x' }).mass, 75.1);
  assert.equal(hazardQuery({ landuse: 'urban' }),
               'landuse=urban&mass=75.1&height=1.77&theta=75&alpha=3.064&beta=0.303&flood_depth=0.1');
  assert.match(describeHazard({ landuse: 'woodland' }), /woodland/);
});
