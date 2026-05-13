# FloodA5 — Testing Reference

_Test plan, acceptance criteria, and analytical benchmarks for the FloodA5 hydraulic model. Intended as project context for testing and validation conversations._

---

## 1. Philosophy

Testing proceeds from simple to complex:

1. **Unit tests** — individual functions (edge geometry, hypsometric lookup, sill computation)
2. **Flat-terrain point-source** — zero-DEM mesh, single injection point; tests routing without topographic complexity
3. **Analytical benchmarks** — Thacker planar surface, circular dam-break; compare to known solutions
4. **Real-domain validation** — Christchurch/Kaiapoi event data

No test should require external data files (GeoTIFF, historical records) at the unit or benchmark level. All data is generated procedurally or from small embedded fixtures.

---

## 2. Existing Test Files

| File | Status | Coverage |
|------|--------|----------|
| `test_a5grid.jl` | Exists — basic mesh/cell API tests | `lonlat_to_cell`, `cell_to_lonlat`, mesh generation |
| `test_edge_geometry.jl` | Exists — edge geometry functions | `_shared_edge`, `_edge_cos_theta`, `_edge_length_m` |

---

## 3. Flat-Terrain Point-Source Test

### Purpose
The simplest possible routing test: a flat mesh (all elevations zero), a single cell receiving constant input volume, no DEM. Expected behaviour is symmetric radial spread from the source cell outward across the domain.

### Setup
```bash
# Generate mesh (or use saved mesh — no DEM needed)
julia --threads auto FloodModel.jl \
    --meshgen example_aoi.geojson --meshres 9 --meshout flat_mesh.parquet \
    --flow-model standard --mesh-only

# Run with single rainfall point (the new --rainpoint flag)
julia --threads auto FloodModel.jl \
    --meshload flat_mesh.parquet \
    --flow-model standard \
    --rainpoint -43.531,172.636,0.05 \
    --sim-duration 3600 \
    --dt-max 30 \
    --output flat_test.h5 \
    --output-interval 300 \
    --vis makie
```

### Acceptance criteria
1. **Water reaches all cells** adjacent to the source cell within a physically reasonable time
2. **No uphill flow** — on a flat mesh all elevations are zero; WSE should propagate outward symmetrically
3. **Mass conservation** — `mb_err` logged every 50 steps should be < 0.1% of domain volume
4. **No NaN/Inf** in `water_depth` or `volume` at any timestep
5. **Velocity is non-zero** in wet cells (after velocity computation is implemented)

### Diagnostic checks (Python)
```python
import h5py, numpy as np

with h5py.File("flat_test.h5") as f:
    # Check mass balance at final frame
    frames = sorted(f["frames"].keys())
    last = f[f"frames/{frames[-1]}"]
    vol = last["volume"][:]
    print(f"Total domain volume: {vol.sum():.2f} m³")
    print(f"Wet cells: {(vol > 0.001).sum()} / {len(vol)}")
    print(f"Max depth: {last['water_depth'][:].max():.4f} m")
    print(f"NaN check: {np.isnan(vol).any()}")
```

---

## 4. Thacker Planar Surface Benchmark

### Description
The Thacker (1981) planar surface test provides an analytical solution for flow in a parabolic bowl. On a flat tilted surface the free surface oscillates as a standing wave. The simpler 1D version is a translation of a water slab across a sloped surface.

For a pentagonal grid the most accessible variant is a **uniform slope** test: place a rectangular water block on one end of a sloped domain and observe it spreading downslope. The diffusive wave approximation does not conserve momentum exactly, so we test:
- Correct flow direction (downslope, not upslope)
- Approximate spreading rate consistent with Manning's equation
- No negative volumes

### Setup
Since FloodA5 does not currently ingest analytical elevation fields directly, this test uses a GeoTIFF DEM generated from a linear elevation function. A small Python script generates the synthetic DEM:

```python
# generate_slope_dem.py — creates a 0.1% east-west slope GeoTIFF
import numpy as np
import rasterio
from rasterio.transform import from_bounds

# 1° × 1° domain, 10m resolution, 0.001 m/m east-west slope
lon0, lat0, lon1, lat1 = -1.35, 51.65, -1.10, 51.80
nx, ny = 2500, 1500
lons = np.linspace(lon0, lon1, nx)
lats = np.linspace(lat0, lat1, ny)  # note: north-up, so reverse for array
elev = np.outer(np.ones(ny), (lons - lon0) / (lon1 - lon0))  # 0→1 m west→east
transform = from_bounds(lon0, lat0, lon1, lat1, nx, ny)
with rasterio.open("slope_dem.tif", "w", driver="GTiff",
                   height=ny, width=nx, count=1, dtype="float32",
                   crs="EPSG:4326", transform=transform) as dst:
    dst.write(elev.astype("float32"), 1)
```

### Acceptance criteria
1. Water moves in the correct direction (eastward = downslope)
2. No reverse flow under quiescent conditions
3. Volume conservation within 1% over 1-hour simulation

---

## 5. Circular Dam-Break Benchmark

### Description
A circular region of elevated water surrounded by a dry flat bed. On release the wave front should propagate radially outward at a speed consistent with the shallow-water wave speed `c = √(g·h)`. The diffusive wave approximation damps this front relative to the full Saint-Venant solution, so the test checks:
- Correct radial symmetry (no directional bias — A5 grid key advantage)
- Positive spreading (no reverse flow)
- Approximate front speed consistent with diffusive wave theory

### Comparison with rectangular grid
A key scientific result for the A5 grid paper: repeat the circular dam-break on a rectangular grid at the same resolution and measure the angular variance of the flood front radius over time. The A5 grid should show lower directional bias (more circular front).

### Setup
```bash
# Initial condition: inject large volume into a cluster of cells at domain centre
# Use --injection-point with a high rate and short duration (or a modified initial condition)
julia --threads auto FloodModel.jl \
    --meshload flat_mesh.parquet \
    --flow-model standard \
    --injection-point -43.531,172.636,1000.0 \
    --sim-duration 60 \
    --dt-max 5 \
    --output dambreak.h5 \
    --output-interval 5
```

For a proper dam-break initial condition, a small code addition is needed: `--initial-depth` or `--initial-volume` to set a non-zero starting volume in cells within a specified radius. This is a planned test infrastructure item.

---

## 6. Mass Conservation Test

### Description
For any closed-domain run (no outflow BC), the total domain volume should equal the cumulative input from all sources.

```
domain_vol(t) = Σ rainfall_rate × Δt × cell_area   +   Σ injection_rate × Δt
```

The model already logs `mb_err` every 50 steps. This test formalises the criterion.

### Acceptance criterion
`|mb_err| / input_vol < 0.5%` at any logged step for a 1-hour simulation with `--dt-max 60`.

Volume loss of up to ~0.5% is expected due to the 50% volume limiter clipping during wetting-front advance. Losses substantially larger than this indicate a bug in the flux loop or the volume update.

---

## 7. Dry-Cell Stability Test

### Description
Run a simulation where cells alternate between wet and dry (intermittent rainfall or a source that runs for part of the simulation). Check that cells that become dry (volume → 0) do not produce NaN or negative volumes, and that rewetting works correctly.

### Setup
```bash
# Inject for first 30 minutes only — then let domain drain (no outflow = water stays)
# Use very small domain and short duration
julia --threads auto FloodModel.jl \
    --meshload flat_mesh.parquet \
    --flow-model standard \
    --injection-point -43.531,172.636,0.1 \
    --sim-duration 7200 \
    --dt-max 30 \
    --output drywet.h5 --output-interval 60
```

Modify `run_simulation!` to stop injection after 1800s for this test (or implement `--source-duration`, which is a useful general feature).

---

## 8. Sign Convention Regression Test

### Description
Directly test the `_bates_flux` function to confirm that water flows downhill, not uphill.

```julia
# test_sign_convention.jl
include("FloodModel.jl")

# Cell i (index 1) is higher: wse_i = 2.0 m, wse_j = 1.0 m
# Expect: Q < 0 (flow from i to j, so ci loses volume)
Q = FloodA5._bates_flux(0.0, 2.0, 1.0, 0.0, 100.0, 1000.0, 1.0, 0.03, 30.0)
@assert Q < 0 "Flow should be from i to j (Q < 0) when WSE_i > WSE_j"

# With Q < 0, dV[ci] += Q*dt < 0 → ci loses volume ✓
# With Q < 0, dV[cj] -= Q*dt > 0 → cj gains volume ✓
println("Sign convention test passed: Q = $Q")
```

---

## 9. Edge Geometry Tests (`test_edge_geometry.jl`)

Already exists. Ensure it covers:
- `_shared_edge` returns `nothing` for non-adjacent cells
- `_edge_cos_theta` returns 1.0 when `_shared_edge` returns `nothing`
- `_edge_cos_theta` returns a value in [0, 1] for all real A5 cell pairs
- `_edge_length_m` is positive and consistent with cell resolution
- `cos_theta` values for typical resolution-14 cells are logged for inspection

---

## 10. Velocity Computation Test (post-implementation)

Once velocity is implemented (see HYDRAULICS.md §5):

```julia
# After one step with non-zero flux:
# velocity[i] > 0 for all wet cells adjacent to the source
# velocity[i] == 0 for all dry cells
# velocity[i] is finite (no NaN, no Inf)
@assert all(isfinite, state.velocity)
@assert all(v -> v >= 0.0, state.velocity)   # scalar magnitude
wet_cells = state.water_depth .> 1e-4
@assert any(state.velocity[wet_cells] .> 0)  # at least some non-zero velocity
```

---

## 11. Test Infrastructure Notes

- All test files guard against being run by `include()` from a harness by checking `abspath(PROGRAM_FILE) == @__FILE__` (same pattern as FloodModel.jl).
- Tests that require mesh generation should save to `test/` subdirectory and check for existence before regenerating (mesh generation is slow).
- HDF5 outputs from tests should go to `test/` and be gitignored.
- Python diagnostic scripts that read HDF5 output live in `test/` alongside the Julia tests.
