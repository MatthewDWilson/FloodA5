# Mesh families and DGGS meshes: review and plan (WP3 part 3)

2026-11-17; updated 2026-11-23. Branch `wp3c-dggs-meshes`, built on 0.6.0 and
rebased onto 0.7.0 after the SGS conveyance merge; **Phase A released as 0.7.1**.
Continues WP3 of `FloodA5_Review_and_ForwardPlan_2026-10.md` ("Remaining:
H3/IGEO7 meshes") and Phase 2 of `FloodA5_Generalisation_BenchmarkingPlan.md`
(multi-DGGS). The same text is in the FloodA5 repository
(`docs/notes/MeshFamilies_DGGS_Plan_2026-11.md`) and in `projectinstructions/`
(`FloodA5_MeshFamilies_DGGS_Plan_2026-11.md`).

**Status (2026-11-23).** Phase A done (§3). Phase B waits for the fixes
flagged during the SGS work that affect it, first the DEM registration fix
(§4, "Sequencing").

## 1. Where things stood

- **The solver already works on any mesh.** It depends on the mesh only
  through per-face and per-cell geometry and the neighbour lists;
  `max_neighbours` comes from the mesh.
- **Three ways to get a mesh, all writing one GeoParquet layout:**
  - A5: `--meshgen`, via `mesh/a5_bridge.py` and pya5;
  - Cartesian: `--meshgen-cartesian`, AOI or raster-aligned;
  - planar polygon families: `mesh/mesh_families.py` (hex, rotated
    Cartesian, Voronoi, import), loaded with `--meshload`.
- **Everything downstream works on any polygon mesh:** DEM sampling, SGS
  tables, output, analysis, replay. Point sources use a nearest-centre
  search, which is not A5-specific.
- **Remaining A5 assumptions** (to fix in Phase B):
  - a mesh without an `id_scheme` is guessed to be `cartesian` (resolution 0)
    or `a5-v0`: `A5Grid._load_from_geoparquet` and
    `FA._parquet_id_scheme` (which mirrors it);
  - pya5's `grid_disk` is the fallback for adjacency when a mesh has no
    `neighbours` column (loader and `build_sgs_tables!`), which is A5-only;
  - `--meshgen`/`--meshres` generate only A5, so a family cannot be chosen
    with one flag;
  - the type names `A5Mesh`/`A5Cell` are cosmetic (add a `PolyMesh` alias,
    don't rename).
- **No DGGS other than A5 had been built.** The library survey of
  2026-09-14 (DGGAL for IGEO7 and HEALPix, healpix-geo as the alternative
  for HEALPix; xdggs and Vgrid not adopted) had not been tested.

## 2. Decisions (with the user, 2026-11-17)

- Families: H3, IGEO7, ISEA3H, HEALPix, rHEALPix. Not S2.
- Phase A: build now on a separate branch; merge after the SGS conveyance
  work (`wps-sgs-conveyance`, 0.7.0) has merged. (Done: 0.7.1.)
- 2026-11-23: the fixes flagged during the SGS work (DEM registration first)
  come before Phase B, which they affect.
- The plan is kept in the branch (this file) and in `projectinstructions/`.

## 3. Phase A: done on this branch (no files the SGS branch changes, except
one line each in `analysis/FloodA5Analysis.jl`, `ci/run_tests.sh` and the
`CHANGELOG.md` `[Unreleased]` heading)

1. **Library spike.**
   - h3 4.5.0, dggal and healpix-geo 0.5.0 install from wheels into the
     Conda.jl Python (Python 3.13) without changing other packages. Linux
     wheels exist for CI's Python 3.12.
   - **dggal 0.0.7's Windows wheel crashes** in `pydggal_setup`. Its
     extension loads the eC runtime as `ecrt-<hash>.dll` while `dggal.dll`
     loads `ecrt.dll`, so two copies of the runtime are in the process.
     0.0.6 works and is pinned; the generator refuses 0.0.7 on Windows.
     Worth reporting upstream.
   - **DGGAL's HEALPix is identical to healpix-geo with
     `ellipsoid="WGS84"`** (authalic latitude; corners agree to < 1e-9°).
     healpix-geo's default `"sphere"` differs.
   - **Decision: DGGAL for all four DGGAL grids.** One dependency;
     healpix-geo is kept as an independent cross-check in the tests.
   - DGGAL returns corner vertices only and side-sharing neighbours only
     (HEALPix: 4). Shared vertices differ by up to 3e-14° between cells, so
     snapping is needed.
   - Speed: about 150 µs per cell for geometry and neighbours, so roughly
     20 s for a Carlisle-sized mesh.
2. **`mesh/dggs_families.py`** — `h3`, `igeo7`, `isea3h`, `healpix`,
   `rhealpix`:
   - `--res N` or `--target-dx M`; `--list-levels`;
   - cells intersecting the AOI (A5's rule);
   - corner vertices snapped to one copy;
   - shared-side adjacency checked against the library's neighbour lists;
   - IDs are the 64-bit index as 16 hex digits;
   - metadata records family, library and version, level, refinement ratio
     and whether the grid is equal-area.
3. **`docs/MESHES.md`** — the mesh contract, families, levels near 20–50 m,
   geometry of each family at Carlisle, family notes.
4. **`analysis/mesh_quality.jl`** (`FA.face_geometry`,
   `FA.orientation_moments`, `FA.mesh_quality`): valence, area uniformity,
   non-orthogonality (as the model computes it), skew, centre offset,
   face-orientation moments, vertex valence. Test
   `test/analysis/test_mesh_quality.jl` (squares, hexagons, 60° rhombi,
   hexagons at 55° N).
5. **`test/test_dggs_meshes.jl`** (in CI; `FLOODA5_REQUIRE_DGGS=1` there):
   generator, geometry and closure, lake at rest (4 options), mass
   conservation, HEALPix against healpix-geo. Libraries added to `setup.jl`
   and `ci/python-requirements.txt`.
6. `visualisation/replay_data.jl`: only A5 meshes get an A5 level in the
   replay diagnostics panel. `FA._parquet_id_scheme` keeps mirroring the
   model's guess until Phase B changes both.

Findings worth carrying forward (detail in `docs/MESHES.md` §3):

- **No DGGS is orthogonal at 55° N in FloodA5's metric:**

  | Family | Non-orthogonality |
  |---|---|
  | H3 | 4.7° |
  | rHEALPix | 7.1° |
  | IGEO7 | 9.6–15.7° |
  | ISEA3H | 13.5° |
  | HEALPix (polar region) | 33° |
  | A5 | 20° median, 37° max |

  All DGGS have zero skew and exact centres. So the "near-orthogonal
  hexagons" of the forward plan §6.1 holds for H3 only roughly, and HEALPix
  at UK latitudes is a strongly non-orthogonal mesh: a useful test for the
  diamond gradient.
- **Equal-area grids stay uniform in FloodA5's area measure** (CV 4e-7).
  FloodA5's areas are 0.3–0.5% below the libraries' ellipsoidal ones
  (6371 km sphere).
- **On a sphere, cells close only approximately.** FloodA5's per-face
  projection gives Σ w·n̂ = O(L/R · tan φ) per cell (1e-5 at 40 m, 55° N),
  not 0; the old family tests ran at the equator and did not see this.
- **A centroid bug was caught and fixed during development.** The shoelace
  centroid in absolute degrees cancels catastrophically at 1e-4° cells
  (0.9 m error at 55° N); `mesh_quality` computes it relative to the first
  vertex. The same formula elsewhere should be checked before use at real
  latitudes.
- `mesh/__pycache__/a5_bridge.cpython-313.pyc` is tracked in git and is
  rewritten whenever a generator imports the bridge. It should be untracked
  (Phase B, or any maintenance branch).

## 4. Phase B: after the SGS merge (and after the DEM registration fix)

7. **CLI.**
   - Add `--meshgen-dggs FAMILY AOI (--meshres N | --target-dx M)`, through
     a thin `mesh/DGGSGrid.jl` that calls the generator, as tests do.
   - Mesh-generation flags go into `FO._DROP_FLAGS`; update `--help` and
     `USER_GUIDE.md`.
   - With this, the mesh is one flag in any run or case.
8. **`A5Grid.jl`:**
   - no `id_scheme` guess for files that record one; for files that don't,
     keep the old guess but warn;
   - the `grid_disk` fallback becomes an error for non-A5 schemes (loader
     and `build_sgs_tables!`);
   - add a `const PolyMesh = A5Mesh` alias;
   - write the DGGS metadata into the output's `/mesh` attributes
     (`io/FloodA5Output.jl`).
9. **SGS on DGGS meshes.** Extend `test/test_sgs_general.jl` (now hex) to
   H3, IGEO7 and HEALPix with the new face tables.
10. **Cases** (`cases/registry/mesh_families.jl`, a new file), each on
    Cartesian, rotated Cartesian, planar hex, centroidal and jittered
    Voronoi, A5, H3, IGEO7, ISEA3H, HEALPix and rHEALPix at matched mean
    cell area:
    - point spread (isotropy);
    - planar slope (asymmetry);
    - a subset of the channel orientation series;
    - Carlisle peak at about 22 m against the 5 m `carlisle_peak_fine`
      references.
11. Docs: ARCHITECTURE (mesh layer), DATA_FORMATS (metadata keys),
    TEST_CASES, README (roadmap line "Further mesh families"), CHANGELOG.

Sequencing: items 10 and 11 come after the **DEM registration fix**
(`_affine_to_pixel` half-pixel offset; first on the SGS thread's "next"
list). It moves every DEM sample, so baselines recorded before it would be
recorded twice.

## 5. Phase C: notebook `notebooks/06_mesh_families.ipynb`

1. Why swappable meshes; the mesh contract (`docs/MESHES.md`).
2. Gallery: one AOI at matched mean cell area, one panel per family.
3. Geometry: the `FA.mesh_quality` table and histograms (non-orthogonality,
   face-orientation roses), and hierarchy properties: refinement ratio,
   whether children tile the parent, equal-area.
4. Hydraulics from the Phase B cases: lake at rest and mass balance,
   point-spread isotropy, planar asymmetry, channel orientation, Carlisle
   against the fine reference.
5. Cost: cells, faces per cell, step time.
6. Advantages and disadvantages per family, from the evidence above.
   Starting hypotheses to test:
   - Cartesian: raster- and LISFLOOD-FP-compatible, but four face directions.
   - A5: equal-area, fewest faces, but non-orthogonal and chiral.
   - H3: the largest ecosystem, but not equal-area, and children don't tile
     the parent.
   - IGEO7: equal-area hexagons, but coarse level steps.
   - ISEA3H: fine level steps.
   - HEALPix: equal-area and its children tile the parent, but skewed
     outside ±41.8°.
   - rHEALPix: equal-area near-squares, but ×9 level steps.
   - Voronoi: can follow features, but has no global index.

## 6. Later (out of scope)

- A5 through the same generator interface.
- Parent IDs for multi-resolution meshes (AMR).
- Triangular meshes (comparison with Sridharan et al. 2021).
- An IVEA/RTEA variant of a grid (DGGAL has them).
