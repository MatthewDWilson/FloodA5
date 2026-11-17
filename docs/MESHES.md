# FloodA5 — Meshes: the mesh contract and the mesh families

FloodA5's solvers depend on the mesh only through per-cell and per-face
geometry (areas, centres, face widths, normals, centre distances, vertex
stencils) and the neighbour lists. Any polygon mesh that satisfies the
contract below runs under every flow model and option. This document
defines that contract, lists the families that produce it, and records
what is known about each family's geometry. A5-specific behaviour (ID
schemes, the dual sublattice, pya5 versions) is in `A5_QUIRKS.md`.

Plan and status of the mesh-family work: `docs/notes/MeshFamilies_DGGS_Plan_2026-11.md`.

## 1. The mesh contract

A FloodA5 mesh is a GeoParquet file (one row per cell, CRS EPSG:4326).
It is read by `A5Grid.load_mesh_geoparquet` (the model, `--meshload`) and by
`FA.load_mesh_parquet` (analysis, pure Julia).

| Column | Type | Meaning |
|---|---|---|
| `cell_id` | string | Unique cell ID. A 64-bit integer as 16 hex digits wherever the family has one (FloodA5 normalises hex IDs to 16 digits, `_norm_cell_id`). |
| `resolution` | int | The family's level (A5 resolution, DGGS level); 0 for Cartesian and the planar families. |
| `center_lon`, `center_lat` | float | The cell centre the solver uses (degrees, WGS84). |
| `geometry` | polygon | The cell: exterior ring of its **corner vertices only**, one straight side per face. |
| `neighbours` | list of string | IDs of the cells that share a **side** (not only a corner). Required for every non-A5 mesh: without it FloodA5 falls back to pya5's `grid_disk`, which is A5-only. |
| `cell_area` | float | Informational. The model recomputes areas from the polygons (`A5Grid._polygon_area_m2`). |
| `elevation`, SGS columns, … | | Added by DEM sampling (`--dem`) and SGS builds; see `DATA_FORMATS.md`. |

Schema metadata `flooda5` (JSON) holds `a5_id_scheme` (the ID scheme, despite
the name: `a5-v0`, `a5-v1`, `cartesian`, `polygon`, `h3`, `igeo7`, `isea3h`,
`healpix`, `rhealpix`), plus `mesh_family` and family-specific keys (for DGGS
meshes: `dggs_library`, `dggs_level`, `dggs_refinement_ratio`,
`dggs_equal_area`, `dggs_system`, `dggs_interior_valence`). Only `a5-*`
schemes are checked against the installed pya5.

Geometric rules:

1. A vertex shared by several cells is **bitwise identical** in each of them.
   FloodA5 finds faces by matching vertices after rounding; copies that differ
   by round-off can straddle a rounding boundary and drop the face. Generators
   snap vertices to one canonical copy (`mesh_families.write_mesh`,
   `dggs_families.VertexSnapper`).
2. Two neighbouring cells share exactly one side (two vertices). No hanging
   nodes (T-junctions): a cell side split by a neighbour's vertex is not a face.
3. Rings are simple polygons. Orientation is free (FloodA5 orients faces
   itself).
4. No AOI across the antimeridian or containing a pole.

What FloodA5 does with a mesh on load (any family): faces from shared sides;
geometry per face in an equirectangular projection about the face's
mid-latitude (`A5Grid._edge_geometry`); boundary cells are those with an
unshared polygon side, which get ghost edges for open/closed boundaries.
One consequence on a sphere: Σ w·n̂ over a cell closes only to
O(L/R · tan φ) (L the cell size, φ the latitude) — about 1e-5 for 40 m cells
at 55° N, and exact at the equator.

## 2. Families

| Family | Generator | Cell | Interior valence | Equal-area | Refinement | ID scheme |
|---|---|---|---|---|---|---|
| A5 | `FloodModel.jl --meshgen AOI --meshres N` (pya5) | pentagon | 5 | yes | ×4 | `a5-v1` (`a5-v0`) |
| Cartesian | `--meshgen-cartesian AOI --dx M`, or a raster (`dem.tif`) | square | 4 | no (lon/lat grid) / yes (projected raster) | — | `cartesian` |
| hex, rotated Cartesian, Voronoi, import | `mesh/mesh_families.py` | planar polygon | 6 / 4 / mixed | locally | — | `polygon` |
| H3 | `mesh/dggs_families.py h3` (h3) | hexagon | 6 | no | ×7 | `h3` |
| IGEO7 (ISEA7H, Z7 index) | `mesh/dggs_families.py igeo7` (DGGAL) | hexagon | 6 | yes | ×7 | `igeo7` |
| ISEA3H | `mesh/dggs_families.py isea3h` (DGGAL) | hexagon | 6 | yes | ×3 | `isea3h` |
| HEALPix (nested) | `mesh/dggs_families.py healpix` (DGGAL) | quadrilateral | 4 | yes | ×4 | `healpix` |
| rHEALPix | `mesh/dggs_families.py rhealpix` (DGGAL) | quadrilateral | 4 | yes | ×9 | `rhealpix` |

Generating and running a DGGS mesh (the generator runs in FloodA5's Python,
Julia's Conda.jl Python with its `Library\bin` on PATH; from Julia,
`A5Grid.PYTHON_CMD[]`):

```bash
python mesh/dggs_families.py igeo7 --aoi AOI.geojson --list-levels        # cell size per level
python mesh/dggs_families.py igeo7 --aoi AOI.geojson --target-dx 25 --out m.parquet
julia --project=. --threads auto FloodModel.jl --meshload m.parquet --dem dem.tif ...
```

`--target-dx M` picks the level whose cell area at the AOI centre is closest
to M² (in log terms); `--res N` sets it directly. Cells are those whose
polygon intersects the AOI (A5's rule). The generator checks its shared-side
adjacency against the library's own neighbour list and refuses a mesh where
they differ, or where two neighbours share other than two vertices.

### Levels near flood-model cell sizes (at 54.9° N, −2.9° E)

Cell size = √(cell area) from the library, at the point:

| Family | Levels (cell size) |
|---|---|
| A5 | 17 (44 m), **18 (22 m)**, 19 (11 m) |
| H3 | 10 (115 m), **11 (43 m)**, **12 (16 m)**, 13 (6 m) |
| IGEO7 | 12 (61 m), **13 (23 m)**, 14 (9 m) |
| ISEA3H | 21 (70 m), 22 (40 m), **23 (23 m)**, 24 (13 m) |
| HEALPix | 16 (100 m), 17 (50 m), **18 (25 m)**, 19 (12 m) |
| rHEALPix | 11 (52 m), **12 (17 m)** |

Families with large refinement ratios (H3, IGEO7 ×7; rHEALPix ×9) cannot
match a target size closely: a comparison "at matched resolution" across
families has to state the actual mean areas.

## 3. Geometry of each family at Carlisle

`FA.mesh_quality` (`analysis/mesh_quality.jl`) on meshes over the same
≈ 640 × 670 m box (`--bbox=-2.905,54.897,-2.895,54.903 --target-dx 30`), and
on the whole Carlisle A5 res-18 mesh (`test/carlisle/carlisle_mesh18_standard.parquet`).
Interior cells and their faces; non-orthogonality = angle between the face
normal and the centre-to-centre line, as the model computes it; skew = where
that line crosses the face, from its midpoint, / face length; centre offset =
declared centre to polygon centroid / √area; vertex valence = cells
meeting at a vertex, over vertices whose cells are all interior; m_k = |Σ w e^{ikφ}| / Σ w over
face-normal azimuths φ (m4 = 1: four face directions; m6 = 1: six).

| Mesh | Mean area (m²) | Area CV | Non-orth. median / max (°) | Skew median / max | Centre offset max | m2 / m4 / m6 | Vertex valence |
|---|---|---|---|---|---|---|---|
| A5 res 18 | 493 | 1.9e-6 | 20.0 / 37.2 | 0.09 / 0.24 | 0.014 | 0.38 / 0.19 / 0.56 | 3 and 4 |
| planar hex (`mesh_families.py`, 15 m side) | 585 | 4.2e-5 | 0.0 / 0.0 | 0 / 0 | 0.000 | 0.00 / 0.00 / 1.00 | 3 |
| H3 level 11 | 1885 | 4.6e-5 | 4.7 / 5.7 | 0 / 0 | 0.000 | 0.09 / 0.08 / 0.98 | 3 |
| IGEO7 level 13 | 524 | 4.0e-7 | 9.6 / 15.7 | 0 / 0 | 0.000 | 0.21 / 0.20 / 0.84 | 3 |
| ISEA3H level 23 | 539 | 4.3e-7 | 13.5 / 14.0 | 0 / 0 | 0.000 | 0.20 / 0.23 / 0.84 | 3 |
| HEALPix level 18 | 616 | 4.1e-7 | 33.3 / 33.3 | 0 / 0 | 0.000 | 0.60 / 0.47 / 0.99 | 4 |
| rHEALPix level 12 | 300 | 3.9e-7 | 7.1 / 7.1 | 0 / 0 | 0.000 | 0.22 / 0.97 / 0.40 | 4 |

Reading the table:

- The equal-area DGGS stay uniform in FloodA5's own area measure (CV ≈ 4e-7);
  H3 and the planar meshes vary by 4e-5 over the box (H3 cells are not
  equal-area; the planar families are equal-area in their local plane, not
  on the sphere). FloodA5's areas are 0.3–0.5% below the libraries'
  ellipsoidal areas at this latitude, the same offset as against British
  National Grid pixels (CHANGELOG, WP3 part 1): FloodA5 measures on a 6371 km
  sphere.
- None of the DGGS is orthogonal in FloodA5's metric at 55° N, though every
  one has zero skew and exact centres (the faces are crossed at their
  midpoints by the centre line). Projection distortion makes the cells
  irregular: hexagons at 5–16°, HEALPix rhombi at 33°. HEALPix is in its
  polar region above 41.8° (authalic) latitude, where its cells are skewed
  rhombi; near the equator they are close to squares. Non-orthogonality is
  what the diamond gradient (`--face-flux-method diamond`) corrects; A5's 20–37° is larger still.
- The face-direction moments show each family's preferred directions: H3
  and the planar hexagons six (m6 ≈ 1), rHEALPix four (m4 ≈ 1), HEALPix
  two pairs of skewed directions (m2 = 0.6). A5's are spread (m6 = 0.56),
  but A5 is not direction-free either (HYDRAULICS §15.12).

## 4. Family notes

**H3** (`h3` 4.x). Not equal-area: cell area varies by up to about 2× over the
globe, smoothly. Aperture 7, so each level is rotated about 19° against the
one above, and children do not tile their parent exactly. Twelve pentagons
on the globe (at icosahedron vertices, mostly in oceans); a mesh containing
one has a 5-valent cell. Near icosahedron edges some H3 cells carry extra
"distortion" vertices; the generator refuses a mesh where such a vertex splits
a face (choose another AOI or level). IDs are H3's 15-hex-digit index,
zero-padded to 16.

**IGEO7** (DGGAL `ISEA7H_Z7`). Equal-area (ISEA projection on the WGS84
authalic sphere) hexagons, aperture 7, indexed with Z7. Shapes vary with
position on the icosahedron face (non-orthogonality 9.6–15.7° within the
640 m box). DGGAL also has `ISEA7H` (other indexing), and IVEA/RTEA
variants of each grid (not exposed).

**ISEA3H** (DGGAL). Equal-area hexagons, aperture 3: the finest choice of
levels (×3 area per level); orientation alternates between levels.

**HEALPix** (DGGAL `HEALPix`, nested). Equal-area quadrilaterals on the WGS84
authalic sphere: identical to `healpix-geo` with `ellipsoid="WGS84"` (corners
agree to < 1e-9°, test X); `healpix-geo`'s default `"sphere"` places cells
differently. Only the four side-sharing neighbours are faces; HEALPix's
corner-only neighbours are not. Polar region above 41.8°: skewed rhombi
(see §3).

**rHEALPix** (DGGAL). Equal-area quadrilaterals, aperture 9 (few usable
levels: 52 m and 17 m around Carlisle), close to a north-aligned square
lattice at mid-latitudes (m4 = 0.97).

**DGGAL on Windows.** `dggal` 0.0.7's Windows wheel crashes in
`pydggal_setup` (access violation): its Python extension loads the eC runtime
as `ecrt-<hash>.dll`, while `dggal.dll` loads `ecrt.dll`, so two copies of
the runtime end up in the process. 0.0.6 (pinned in `setup.jl` and
`ci/python-requirements.txt`) works; the generator refuses 0.0.7 on Windows.

## 5. Tests

- `test/test_mesh_families.jl` — planar families, at the equator.
- `test/test_dggs_meshes.jl` — DGGS families near Carlisle: generator (level
  choice, vertex snapping, reproducibility), geometry and closure, lake at
  rest under four solver options, mass conservation, HEALPix against
  healpix-geo. A family whose library is missing is skipped, except in CI
  (`FLOODA5_REQUIRE_DGGS=1`).
- `test/analysis/test_mesh_quality.jl` — `FA.mesh_quality` on meshes with
  known answers (squares, regular hexagons, 60° rhombi).
