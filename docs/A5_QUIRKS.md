# FloodA5 — A5 Grid Quirks Reference

_Consolidated reference for A5 DGGS-specific behaviour that has caused bugs or requires non-obvious handling. Read this before touching mesh generation, adjacency, or any code that processes cell IDs. Intended as project context for mesh and grid-related development conversations._

---

## 1. Dual Sublattice — The Most Important Quirk

**Root cause of Bug 35 — the most significant mesh bug in the project.**

At every A5 resolution level, pentagons tile in **two interleaved sublattices**. At resolution 14, these correspond to cells whose hex IDs begin with `8e…` and `9e…` (the leading nibble encodes the sublattice). The two sublattices are geometrically complementary — every `9e…` cell is surrounded exclusively by `8e…` neighbours, and vice versa.

**`fill_polygon + uncompact` returns only one sublattice.** If you call `pya5.fill_polygon(aoi, resolution)` followed by `pya5.uncompact(cells, resolution)`, you get only the `9e…` cells (or only `8e…` — which sublattice depends on the AOI and resolution). Their edge-neighbours belong entirely to the other sublattice, which `fill_polygon` never returns.

**Consequence:** A mesh built from `fill_polygon` alone has zero edges — no cell has any neighbour in the mesh. The EdgeList is empty and the model produces no flux.

**Fix (in `a5_bridge.py`):**
After `fill_polygon + uncompact`, collect all `grid_disk(cell, 1)` neighbours of the primary cells, `uncompact` them to the target resolution, filter to those whose centres lie inside the AOI, and add them:

```python
primary = set(pya5.uncompact(pya5.fill_polygon(aoi_wkt, resolution), resolution))
neighbours = set()
for cell in primary:
    disk = pya5.uncompact(pya5.grid_disk(cell, 1), resolution)
    neighbours.update(disk)
# Add neighbours whose centres are inside the AOI
inside = {c for c in neighbours - primary if point_in_polygon(cell_centre(c), aoi)}
all_cells = primary | inside
```

**Any mesh file generated before this fix (Bug 35, 2026-03-31) is incorrect and must be regenerated with `--meshgen`.**

---

## 2. Cell ID Padding — `u64_to_hex` vs `_to_hex`

**Root cause of Bug 32.**

pya5's `u64_to_hex` function may omit leading zeros from cell IDs:
- `u64_to_hex(0x08a2a1072b59ffff)` → `"8a2a1072b59ffff"` (15 chars)

Julia's internal `A5Grid._to_hex` always pads to 16 characters:
- `_to_hex(0x08a2a1072b59ffff)` → `"08a2a1072b59ffff"` (16 chars)

If you compare IDs from parquet (15-char) with IDs from `grid_disk_neighbours` (16-char), no matches are found and the adjacency dict is empty — zero edges, no flux.

**Fix:** Normalise every cell ID to 16-character zero-padded hex at the earliest opportunity:

```julia
_norm_id(id::String) = A5Grid._to_hex(parse(UInt64, id, base=16))
```

This normalisation is applied:
- At the start of `initialise_flow_model` for all IDs from `mesh.cells`
- Inside `_build_adjacency_grid_disk`, `_build_adjacency_matrix!`, `_build_edge_list`
- In Python: `a5_bridge.py` writes IDs via `cell.id.hex()` padded to 16 chars

**Rule:** Never compare cell IDs without normalising both sides first.

---

## 3. `grid_disk` Returns Compact Mixed-Resolution Results

**Root cause of Bug 34.**

`pya5.grid_disk(cell, 1)` returns a **compact** result — a mix of cells at different resolution levels. For a cell at resolution 14, the disk may include cells at resolutions 13, 15, or even coarser. When you filter `filter(id -> id in all_ids, disk)`, the mismatched resolutions give zero matches.

**Fix:** Always `uncompact` after `grid_disk`:

```python
disk = pya5.uncompact(pya5.grid_disk(cell, 1), target_resolution)
```

This expands the compact mixed-resolution disk to a uniform set of cells at `target_resolution`.

The same applies in Julia's `grid_disk_neighbours`:
```julia
function grid_disk_neighbours(cell_id::String, resolution::Int)
    # ... calls pya5.grid_disk via Python bridge
    # Must call uncompact before returning
end
```

---

## 4. Longitude Normalisation

pya5 returns cell centre longitudes in the range [-180, 180]. However, some coordinate paths can produce longitudes outside this range (e.g., cells near the antimeridian).

**Always normalise:**
```python
lon = ((lon + 180) % 360) - 180
```

In Julia:
```julia
_norm_lon(lon::Float64) = mod(lon + 180.0, 360.0) - 180.0
```

This is applied in `a5_bridge.py` during GeoParquet write. If implementing new coordinate paths, ensure normalisation is applied before storing or comparing longitudes.

---

## 5. Resolution Level Conventions

A5 resolution levels increase with finer cell size (higher number = smaller cells, more cells). The relationship between level and cell area is approximately:

| Level | Approx. cell area | Cell spacing (≈√area) | Typical use |
|-------|-------------------|----------------------|-------------|
| 5  | ~33,100 km²  | ~182 km    | Continental |
| 6  | ~8,295 km²   | ~91 km     | Sub-continental |
| 7  | ~2,071 km²   | ~45.5 km   | Regional |
| 8  | ~518 km²     | ~22.8 km   | Large catchment |
| 9  | ~129 km²     | ~11.4 km   | Medium-large catchment |
| 10 | ~32.4 km²    | ~5.7 km    | Medium catchment |
| 11 | ~8.09 km²    | ~2.84 km   | Small-medium catchment |
| 12 | ~2.02 km²    | ~1.42 km   | Small catchment |
| 13 | ~50.6 ha     | ~711 m     | Sub-catchment |
| 14 | ~12.6 ha     | ~356 m     | Urban / detailed |
| 15 | ~3.16 ha     | ~178 m     | Fine urban |
| 16 | ~7,902 m²    | ~89 m      | High-res urban |
| 17 | ~1,976 m²    | ~44 m      | LiDAR-scale |
| 18 | ~494 m²      | ~22 m      | High-res channel |
| 19 | ~124 m²      | ~11 m      | Very high-res |
| 20 | ~30.9 m²     | ~5.6 m     | Ultra high-res |

A5 cell areas are exactly equal for all cells at a given resolution — this equal-area property is a fundamental design characteristic of the A5 DGGS and holds globally regardless of latitude. The values above were measured at lon=172.636°, lat=−43.531° (Christchurch, NZ) but apply everywhere. Use `_polygon_area_m2(cell.boundary)` for the exact geodetic area of each cell as computed from its polygon boundary.

---

## 6. Adjacency: 5 Neighbours for Interior Cells, Fewer for Boundary Cells

Interior A5 cells have exactly 5 edge-sharing neighbours. Cells near the boundary of the AOI may have 3 or 4 neighbours (some neighbours fall outside the AOI and are excluded from the mesh).

The `adj_matrix` is sized `(max_nb=5, n_cells)`. Empty slots are `0` (zero-indexed as "no neighbour"). The EdgeList is built only from non-zero slots.

**In the flux loop:** Never assume 5 neighbours per cell. The EdgeList only contains edges for pairs that are both in the mesh, so boundary effects are handled automatically.

---

## 7. Parent–Child Relationship and Phase 3 AMR

At adjacent resolution levels, A5 cells obey a strict parent–child relationship:
- One parent cell at level L contains exactly 5 children at level L+1
- Children cover the same geographic area as the parent without gaps or overlaps

This property is exploited in Phase 3 (AMR) for conservative refinement: a coarse cell can be split into 5 fine cells, conserving volume, and reassembled without re-meshing.

**Implication for current code:** `_edge_cos_theta` and `_edge_length_m` are already designed to accept cell pairs at different resolution levels (they take raw boundary arrays and centre coordinates). No interface change is needed for Phase 3.

---

## 8. A5 Cell ID Structure

A5 cell IDs are 64-bit unsigned integers, represented as 16-character hex strings. The bit structure encodes:
- Resolution level (high bits)
- Hierarchical position (remaining bits)

The leading nibble (`8` vs `9`) encodes the sublattice for a given resolution (see §1). Do not rely on specific bit patterns — use `pya5` functions for all topological operations (`grid_disk`, `uncompact`, `get_resolution`, `cell_to_boundary`).

**Serialisation:** Always use `_to_hex` / `u64_to_hex` with 16-char padding. Never use Python's `hex()` without padding — it omits leading zeros.

---

## 9. Coordinate Reference System

All coordinates are in **EPSG:4326** (WGS84 geographic, lon/lat in degrees). There is no CRS reprojection in the A5 grid layer — pya5 works natively in geographic coordinates.

Area and distance computations use haversine/geodetic formulas:
- `_haversine_m(lon1, lat1, lon2, lat2)` — great-circle distance
- `_polygon_area_m2(boundary)` — Shoelace formula on an equirectangular projection centred on the polygon centroid (< 0.1% error at A5 res 14)
- `_edge_cos_theta` — local equirectangular projection centred on the edge midpoint (< 0.05% error at res 14)

This avoids the need for projected CRS selection, which would introduce distortion over large or high-latitude domains.

---

## 10. Phase 3 Multi-Resolution Considerations

Several design decisions in the current single-resolution code anticipate Phase 3:

| Feature | Current state | Phase 3 relevance |
|---------|--------------|-------------------|
| `EdgeList.cell_i < cell_j` ordering | Enforced | Works for mixed-resolution cells as long as indices are consistent |
| `_edge_cos_theta` accepts raw boundaries | Already takes raw arrays | Can handle coarse/fine boundary pairs |
| `adj_matrix` retained alongside EdgeList | Yes | Needed for refinement trigger queries (wet/dry front detection) |
| `sgs_cell_area` stored in parquet | Yes | Cell area needed for conservative volume redistribution at refinement |
| `FlowState.sgs_tables` is `Vector{Any}` | Placeholder | Will become `Vector{Union{SGSTable, Nothing}}` in Phase 3 |

The primary Phase 3 data structure addition is `MultiResMesh`: a multi-level cell store with parent/child pointers and a coarse/fine interface flux handler.
