"""
A5Grid.jl
---------
Julia interface to the A5 pentagonal DGGS.

Architecture
------------
Mesh generation is delegated entirely to `a5_bridge.py` via subprocess.
Python handles all pya5 calls, coordinate conventions, WKB encoding, and
GeoParquet/GeoJSON writing — eliminating all PyCall interop issues with
UInt64 cell IDs, pya5's internal coordinate conventions, and WKB geometry.

Individual cell queries (lonlat_to_cell, cell_to_boundary, etc.) use PyCall
directly for use in the flow model, with a GIL lock for thread safety.

GPU/CPU parallelism is used for the point-in-polygon sampling step, which
is the bottleneck at high resolution. The PIP kernel runs in parallel; the
pya5 cell indexing runs serially in Python (GIL-safe).

Static variables
----------------
Each cell can carry static (time-invariant) variables: elevation, friction,
porosity, etc. These are stored in `A5Mesh.static_vars` as a
`Dict{String, Vector{Float64}}` keyed by variable name. The vector indices
correspond to cells in `A5Mesh.cells` order.

  mesh.static_vars["elevation"]   # Vector{Float64}, NaN where unsampled
  mesh.static_vars["friction"]    # added the same way when available

Static variables are persisted as extra columns in the GeoParquet file so
that `--meshload` restores them without re-sampling. The `A5Cell.elevation`
field is a convenience accessor (kept for backwards compatibility with code
that references it directly); it mirrors `static_vars["elevation"]`.

DEM ingestion
--------------
  sample_dem_centroid!(mesh, dem_path; strict=false)
      Sample DEM at each cell centre. Fast. Alias: sample_dem!
  sample_dem_mean!(mesh, dem_path; n_samples=256, halton_seed=0, strict=false)
      Sample DEM at n_samples Halton quasi-random points per cell polygon
      and store the arithmetic mean. Better physics, ~n_samples× slower.
  check_dem_coverage(mesh, dem_path)
      Coverage report without sampling.

DEM operations use ArchGDAL.jl (required — Pkg.add("ArchGDAL")).
Query points (EPSG:4326) are reprojected into the DEM native CRS for
sampling — faster and more accurate than warping the raster.

DEM sources are abstracted behind the `DEMSource` type hierarchy so that
future global providers can be added without touching callers:

  Abstract type: DEMSource
    FileDEM(path)           — user-supplied local GeoTIFF   [implemented]
    CopernicusDEM()         — GLO-30, tiles auto-downloaded  [future]
    SRTMDEM()               — SRTM v3, tiles auto-downloaded [future]

  sample_dem!(mesh, source::DEMSource; strict=false)

Public API
----------
    mesh_for_aoi(geojson_path_or_str, resolution)           → A5Mesh
    save_mesh(mesh, path; format=:geoparquet)               # :geoparquet | :geojson
    save_mesh_geoparquet(mesh, path)
    save_mesh_geojson(mesh, path)
    load_mesh_geoparquet(path)                              → A5Mesh
    mesh_to_geojson_string(mesh)                            → String
    sample_dem!(mesh, source; strict=false)                 # mutates mesh.static_vars
    check_dem_coverage(mesh, dem_path)                      → NamedTuple
    cell_to_boundary(cell_id_hex)                           → Vector{Vector{Float64}}
    cell_to_lonlat(cell_id_hex)                             → Tuple{Float64,Float64}
    lonlat_to_cell(lon, lat, resolution)                    → String
    cell_to_children(cell_id_hex, resolution)               → Vector{String}
    cell_to_parent(cell_id_hex)                             → String
    get_resolution(cell_id_hex)                             → Int
    load_aoi(path_or_str)                                   → String
    mesh_summary(mesh)                                      → String
"""
module A5Grid

export A5Mesh, A5Cell, PolyMesh, PolyCell
export pya5_version, a5_id_scheme, check_id_scheme, migrate_cell_id
export DEMSource, FileDEM
export mesh_for_aoi, save_mesh, save_mesh_geoparquet, save_mesh_geojson, load_mesh_geoparquet
export mesh_to_geojson_string
export sample_dem!, sample_dem_centroid!, sample_dem_mean!, check_dem_coverage
export build_sgs_tables!, SGSTable
export wse_from_volume, wetted_area_from_wse, wet_centroid, sgs_table
export face_apk, face_crit, face_table_from_profile, SGS_FACE_LEVELS, boundary_side_slots
export cell_to_boundary, cell_to_lonlat, lonlat_to_cell
export cell_to_children, cell_to_parent, get_resolution, grid_disk_neighbours
export load_aoi, mesh_summary
export _edge_geometry

using PyCall
using JSON3
using ArchGDAL
using Statistics: mean

# ---------------------------------------------------------------------------
# Optional GPU support — loaded at runtime, no compile-time CUDA dependency
# ---------------------------------------------------------------------------

const _cuda_available   = Ref{Bool}(false)
const _cuda_tried       = Ref{Bool}(false)
const _module_init_done = Ref{Bool}(false)

"""
    _gpu_available() → Bool

Load CUDA.jl on first use (WP4, 2026-10: it was loaded in `__init__`, ~5 s
of every start-up, although only point-in-polygon sampling during mesh
generation can use it). `FLOODA5_NO_CUDA=1` skips CUDA entirely. Callers
must reach the GPU path through `Base.invokelatest`, since `_pip_mask_gpu`
is defined at run time.
"""
function _gpu_available()::Bool
    if !_cuda_tried[]
        _cuda_tried[] = true
        if get(ENV, "FLOODA5_NO_CUDA", "") == "1"
            @info "A5Grid: FLOODA5_NO_CUDA=1 — using CPU threads"
        else
            _try_load_cuda()
        end
    end
    return _cuda_available[]
end

function _try_load_cuda()
    try
        @eval A5Grid using CUDA
        _cuda_available[] = @eval A5Grid CUDA.functional()
        if _cuda_available[]
            @info "A5Grid: CUDA GPU detected — PIP sampling will run on GPU"
            @eval A5Grid _define_gpu_path()
        else
            @info "A5Grid: CUDA loaded but no functional GPU — using CPU threads"
        end
    catch
        @info "A5Grid: CUDA.jl not available — using CPU threads"
    end
end

function _define_gpu_path()
    @eval A5Grid begin
        function _pip_mask_gpu(lons::Vector{Float64}, lats::Vector{Float64},
                                ring_xs::Vector{Float64}, ring_ys::Vector{Float64})::Vector{Bool}
            n_points = length(lons)
            n_verts  = length(ring_xs)
            d_lons  = CUDA.CuArray(lons)
            d_lats  = CUDA.CuArray(lats)
            d_rxs   = CUDA.CuArray(ring_xs)
            d_rys   = CUDA.CuArray(ring_ys)
            d_mask  = CUDA.zeros(Bool, n_points)

            function _pip_kernel(lons, lats, rxs, rys, mask, np, nv)
                i = (CUDA.blockIdx().x - 1) * CUDA.blockDim().x + CUDA.threadIdx().x
                i > np && return
                lon, lat = lons[i], lats[i]
                inside = false
                j = nv
                for k in 1:nv
                    xi, yi = rxs[k], rys[k]
                    xj, yj = rxs[j], rys[j]
                    if ((yi > lat) != (yj > lat)) &&
                       (lon < (xj - xi) * (lat - yi) / (yj - yi + 1f-15) + xi)
                        inside = !inside
                    end
                    j = k
                end
                mask[i] = inside
                return
            end

            threads = 256
            blocks  = cld(n_points, threads)
            CUDA.@cuda threads=threads blocks=blocks _pip_kernel(
                d_lons, d_lats, d_rxs, d_rys, d_mask, n_points, n_verts)
            return Array(d_mask)
        end
    end
end

# ---------------------------------------------------------------------------
# Python bridge helpers
# ---------------------------------------------------------------------------

const BRIDGE_DIR    = @__DIR__
# Suffix for temporary files in BRIDGE_DIR (process id + random): several
# mesh generations may run in one checkout at once (the test-case runner).
_tmp_tag() = string(getpid(), "_", string(rand(UInt32); base = 16))
const BRIDGE_SCRIPT = "a5_bridge.py"
const PYTHON_CMD    = Ref{String}("python")   # safe fallback; overwritten by _find_python()

function _find_python()::Bool
    # Already resolved to a non-default path in this session
    !isempty(PYTHON_CMD[]) && PYTHON_CMD[] != "python" && return true

    # Priority 1: use PyCall's configured Python — this is the same interpreter
    # that setup.jl installed pya5 into, so it's guaranteed to have the right
    # packages and avoids the Windows App Execution Alias problem entirely.
    try
        py = PyCall.python
        if !isempty(py) && isfile(py)
            PYTHON_CMD[] = py
            return true
        end
    catch end

    # Priority 2: search common fixed paths and PATH entries.
    # On Windows, "python" and "python3" may be intercepted by App Execution
    # Aliases which launch the Microsoft Store instead of a real interpreter.
    # We detect this by checking that the output contains "Python X.Y.Z" and
    # that the command is not the Windows stub (which prints a Store redirect
    # message and exits with code 9009 or 0 depending on the system).
    candidates = if Sys.iswindows()
        [raw"C:\Python312\python.exe", raw"C:\Python311\python.exe",
         raw"C:\Python310\python.exe", raw"C:\Python39\python.exe",
         "python", "python3"]
    else
        ["python3", "python"]
    end

    for py in candidates
        try
            buf = IOBuffer()
            # Capture both stdout and stderr: the real Python prints to stderr,
            # the Windows Store alias prints to stdout.
            result = run(pipeline(ignorestatus(Cmd([py, "--version"])),
                                  stdout=buf, stderr=buf))
            output = String(take!(buf))
            # Require "Python X.Y" in the output AND a sane exit code.
            # The Windows App Execution Alias exits 9009 or prints no version.
            if result.exitcode == 0 && occursin(r"Python \d+\.\d+", output)
                PYTHON_CMD[] = py
                return true
            end
        catch end
    end

    # Last resort: warn with actionable advice
    @warn "Python not found via PyCall or PATH. " *
          "If you are on Windows, check Settings → Apps → Advanced app settings → " *
          "App execution aliases and disable the Python alias. " *
          "Then re-run setup.jl to verify the bridge."
    return false
end

"""Ensure Python is resolved before any inline subprocess call."""
macro ensure_python()
    return :((_find_python() || error("Python not found. Install Python 3.9+ and ensure it is on PATH.")))
end

"""
    _run_bridge(args...) → Dict

Run a5_bridge.py with the given arguments and return the parsed JSON result.
"""
function _run_bridge(args::String...)
    _find_python() || error("Python not found — run setup.jl to configure the environment.")
    bridge_path = joinpath(BRIDGE_DIR, BRIDGE_SCRIPT)
    isfile(bridge_path) || error("Bridge script not found: $bridge_path")

    cmd = Cmd(Cmd([PYTHON_CMD[], BRIDGE_SCRIPT, args...]); dir=BRIDGE_DIR)
    stdout_buf = IOBuffer()
    stderr_buf = IOBuffer()
    try
        run(pipeline(cmd, stdout=stdout_buf, stderr=stderr_buf))
    catch e
        py_err = String(take!(stderr_buf))
        py_out = String(take!(stdout_buf))
        error("Python bridge failed.\n" *
              "Python: $(PYTHON_CMD[])\n" *
              "Bridge: $bridge_path\n" *
              "STDERR:\n$(isempty(py_err) ? "(empty)" : py_err)\n" *
              "STDOUT:\n$(isempty(py_out) ? "(empty)" : py_out)\n" *
              "Hint: run setup.jl to verify the Python environment and bridge.")
    end
    out   = String(take!(stdout_buf))
    isempty(strip(out)) && error("Python bridge returned empty output")
    # The bridge may print intermediate status lines before the final result.
    # Use the last non-empty JSON line.
    lines = filter(l -> !isempty(strip(l)), split(out, "\n"))
    isempty(lines) && error("Python bridge returned no output lines")
    result = JSON3.read(lines[end], Dict)
    haskey(result, "error") && error("Python bridge error: $(result["error"])")
    return result
end

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

"""
A single A5 pentagon cell.
`elevation` is a convenience field mirroring `mesh.static_vars["elevation"]`;
it is NaN when no DEM has been sampled onto this mesh.
"""
struct A5Cell
    id::String
    resolution::Int
    center_lon::Float64
    center_lat::Float64
    boundary::Vector{Vector{Float64}}
    elevation::Float64   # metres above datum; NaN if not yet sampled
end

# Constructor without elevation (backwards compatibility)
A5Cell(id, res, lon, lat, boundary) = A5Cell(id, res, lon, lat, boundary, NaN)

"""
A5 pentagon mesh covering an AOI.

Fields
------
  resolution   — A5 resolution level
  aoi_geojson  — original AOI GeoJSON string (empty when loaded from parquet)
  cells        — ordered vector of A5Cell
  static_vars  — scalar static variables per cell, keyed by variable name.
                 Index i in each vector corresponds to cells[i].
                 Standard keys: "elevation", "friction", "porosity".
                 NaN = variable not yet assigned for that cell.

                 Example:
                   mesh.static_vars["elevation"][42]  # metres for cell 42
"""
mutable struct A5Mesh
    resolution::Int
    aoi_geojson::String
    cells::Vector{A5Cell}
    static_vars::Dict{String, Vector{Float64}}
    array_vars::Dict{String, Matrix{Float64}}
    adjacency::Dict{String, Vector{String}}   # cell_id → neighbour hex IDs
    # Cell-ID scheme (WP-M1, 2026-10): "a5-v0" (pya5 0.x IDs), "a5-v1"
    # (pya5 1.x IDs) or "cartesian" (CartesianGrid's own row/column IDs).
    # The same cell has different 64-bit IDs under pya5 0.x and 1.x, and an
    # old ID is still a *valid* 1.x ID for a different cell, so meshes record
    # their scheme and load_mesh_geoparquet refuses a mismatch with the
    # installed pya5 — see check_id_scheme and mesh/migrate_mesh_ids.py.
    id_scheme::String
    # The file's other `flooda5` schema metadata (0.8.1): for DGGS meshes the
    # family, library, level, refinement ratio and equal-area flag
    # (mesh/dggs_families.py, docs/MESHES.md); for planar families the family.
    # Kept through re-saves (save_mesh_geoparquet) and written to the output's
    # /mesh attributes. Empty for meshes generated in this process by pya5 or
    # CartesianGrid.
    meta::Dict{String,Any}
end

# Backwards-compatible constructors. A mesh built in this process without an
# explicit scheme uses the installed pya5's scheme.
A5Mesh(res, aoi, cells, sv, av, adj, scheme) =
    A5Mesh(res, aoi, cells, sv, av, adj, scheme, Dict{String,Any}())
A5Mesh(res, aoi, cells) =
    A5Mesh(res, aoi, cells,
           Dict{String,Vector{Float64}}(),
           Dict{String,Matrix{Float64}}(),
           Dict{String,Vector{String}}())
A5Mesh(res, aoi, cells, sv) =
    A5Mesh(res, aoi, cells, sv,
           Dict{String,Matrix{Float64}}(),
           Dict{String,Vector{String}}())
A5Mesh(res, aoi, cells, sv, av) =
    A5Mesh(res, aoi, cells, sv, av, Dict{String,Vector{String}}())
A5Mesh(res, aoi, cells, sv, av, adj) =
    A5Mesh(res, aoi, cells, sv, av, adj, a5_id_scheme())

# The mesh types hold any polygon mesh (A5, Cartesian, planar families, DGGS);
# the A5 names are historical. New code may use the neutral aliases.
const PolyMesh = A5Mesh
const PolyCell = A5Cell

Base.length(m::A5Mesh) = length(m.cells)
Base.iterate(m::A5Mesh, s=1) = s > length(m.cells) ? nothing : (m.cells[s], s+1)
Base.show(io::IO, m::A5Mesh) =
    print(io, "A5Mesh(resolution=$(m.resolution), cells=$(length(m.cells)), " *
              "id_scheme=$(m.id_scheme), " *
              "static_vars=[$(join(keys(m.static_vars), ", "))])")

# ---------------------------------------------------------------------------
# DEM source abstraction
# ---------------------------------------------------------------------------

"""
Abstract supertype for DEM data sources.

Subtypes
--------
  FileDEM(path)      — Local GeoTIFF supplied by the user [implemented]
  CopernicusDEM()    — GLO-30 global 30m DEM              [future]
  SRTMDEM()          — SRTM v3 global 30m DEM             [future]

Implementing a new provider
---------------------------
1. Define `struct MyDEM <: DEMSource` with whatever config fields you need.
2. Implement `_resolve_dem_path(source::MyDEM, mesh::A5Mesh) → String`
   which returns a local GeoTIFF path, downloading/merging tiles as needed.
3. The `sample_dem!` and `check_dem_coverage` functions call `_resolve_dem_path`
   and then delegate to the Python bridge — no other changes needed.
"""
abstract type DEMSource end

"""
    FileDEM(path)

A local GeoTIFF DEM file supplied directly by the user.
The file is validated for existence and basic readability before sampling.
"""
struct FileDEM <: DEMSource
    path::String
end

# ---------------------------------------------------------------------------
# TODO: Future global DEM providers
# ---------------------------------------------------------------------------
# struct CopernicusDEM <: DEMSource
#     resolution_arcsec :: Int   # 30 (default) or 90
# end
# function _resolve_dem_path(src::CopernicusDEM, mesh::A5Mesh)::String
#     # 1. Compute bbox from mesh cell centres
#     # 2. Download Copernicus GLO-30 tiles from the AWS open-data bucket
#     #    (s3://copernicus-dem-30m/) or ESA OpenTopography API
#     # 3. Merge tiles using ArchGDAL.gdalwarp / GDAL.gdalwarp_util
#     # 4. Return merged GeoTIFF path (then sample_dem! handles the rest)
#     error("CopernicusDEM not yet implemented")
# end
#
# struct SRTMDEM <: DEMSource end
# function _resolve_dem_path(src::SRTMDEM, mesh::A5Mesh)::String
#     # Same pattern: download → merge with ArchGDAL → return GeoTIFF path
#     error("SRTMDEM not yet implemented")
# end

"""
    _resolve_dem_path(source::FileDEM, mesh::A5Mesh) → String

Validate and return the path to the user-supplied GeoTIFF.
ArchGDAL can open any GDAL-supported format, but GeoTIFF is expected here.
"""
function _resolve_dem_path(source::FileDEM, ::A5Mesh)::String
    isfile(source.path) ||
        error("DEM file not found: $(source.path)")
    ext = lowercase(splitext(source.path)[2])
    ext in (".tif", ".tiff", ".geotiff") ||
        @warn "DEM file extension '$ext' is unusual — expected .tif/.tiff. " *
              "Proceeding anyway; ArchGDAL will attempt to open it."
    return source.path
end

# ---------------------------------------------------------------------------
# DEM coverage check — pure Julia via ArchGDAL
# ---------------------------------------------------------------------------

"""
    _crs_gis(epsg_or_wkt) → ISpatialRef

Return an ArchGDAL spatial reference with axis mapping forced to traditional
GIS order (x = Easting/Longitude, y = Northing/Latitude).

Without this, GDAL ≥ 3.0 honours the official CRS axis order from the
EPSG registry — NZTM2000 (EPSG:2193) is officially (Northing, Easting),
which means GDAL treats x-coordinates as northings and y-coordinates as
eastings.  Our coordinate arrays are always (lon/easting, lat/northing), so
we must override the axis mapping on every CRS object we create.
"""
function _crs_gis(epsg::Int)
    crs = ArchGDAL.importEPSG(epsg)
    ArchGDAL.GDAL.osrsetaxismappingstrategy(crs, ArchGDAL.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    crs
end

function _crs_gis(wkt::String)
    crs = ArchGDAL.importWKT(wkt)
    ArchGDAL.GDAL.osrsetaxismappingstrategy(crs, ArchGDAL.GDAL.OAMS_TRADITIONAL_GIS_ORDER)
    crs
end

"""
    check_dem_coverage(mesh, dem_source) → NamedTuple

Check how well a DEM covers the mesh without sampling it.
Opens the file with ArchGDAL, reads the geotransform + CRS, reprojects the
bounding box to EPSG:4326, then does a fast bounding-box test against all
cell centres.

Returns a named tuple with:
  - total_cells    :: Int
  - cells_covered  :: Int
  - cells_outside  :: Int
  - coverage_pct   :: Float64
  - dem_bounds     :: NamedTuple(west, south, east, north)  [EPSG:4326]
  - dem_crs        :: String
  - dem_path       :: String
"""
function check_dem_coverage(mesh::A5Mesh, dem_source::DEMSource)
    dem_path = _resolve_dem_path(dem_source, mesh)

    west, south, east, north, crs_str = ArchGDAL.read(dem_path) do ds
        gt  = ArchGDAL.getgeotransform(ds)
        nx  = ArchGDAL.width(ds)
        ny  = ArchGDAL.height(ds)
        # Affine geotransform: [ulx, xres, xskew, uly, yskew, yres]
        ulx = gt[1];  uly = gt[4]
        lrx = gt[1] + nx * gt[2] + ny * gt[3]
        lry = gt[4] + nx * gt[5] + ny * gt[6]
        native_west  = min(ulx, lrx);  native_east  = max(ulx, lrx)
        native_south = min(uly, lry);  native_north = max(uly, lry)

        proj    = ArchGDAL.getproj(ds)
        crs_str = isempty(proj) ? "unknown" : proj

        # Reproject the four corners to EPSG:4326 to get lat/lon bounds
        w, s, e, n = if isempty(proj)
            @warn "DEM has no CRS — assuming EPSG:4326 for coverage check."
            native_west, native_south, native_east, native_north
        else
            src_crs = _crs_gis(proj)
            wgs84   = _crs_gis(4326)
            xs = [native_west, native_east, native_east, native_west]
            ys = [native_south, native_south, native_north, native_north]
            zs = [0.0, 0.0, 0.0, 0.0]
            ArchGDAL.createcoordtrans(src_crs, wgs84) do xform
                ArchGDAL.transform!(xs, ys, zs, xform)
            end
            minimum(xs), minimum(ys), maximum(xs), maximum(ys)
        end

        w, s, e, n, crs_str
    end

    n_total   = length(mesh.cells)
    n_outside = count(c -> c.center_lon < west || c.center_lon > east ||
                            c.center_lat < south || c.center_lat > north,
                      mesh.cells)
    n_covered = n_total - n_outside

    return (
        total_cells   = n_total,
        cells_covered = n_covered,
        cells_outside = n_outside,
        coverage_pct  = 100.0 * n_covered / max(1, n_total),
        dem_bounds    = (west=west, south=south, east=east, north=north),
        dem_crs       = crs_str,
        dem_path      = dem_path,
    )
end

# Convenience overload accepting a path string directly
check_dem_coverage(mesh::A5Mesh, dem_path::String) =
    check_dem_coverage(mesh, FileDEM(dem_path))

# ---------------------------------------------------------------------------
# DEM sampling — pure Julia via ArchGDAL
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Shared raster loading helper
# ---------------------------------------------------------------------------

"""
    _load_dem_raster(dem_path) → NamedTuple

Open a GeoTIFF with ArchGDAL and read all metadata plus the full band-1
data array into a Julia Matrix{Float64}.  Loading the entire raster once
and working from the in-memory array is vastly faster than making one
ArchGDAL.read window call per sample point.

Returns:
  (band, gt, nx, ny, nodata_val, proj)
  band      — Matrix{Float64}, size (nx, ny), column-major (col, row)
  gt        — geotransform vector (length 6)
  nx, ny    — raster width and height in pixels
  nodata_val — nodata sentinel as Float64, or NaN if none defined
  proj      — WKT projection string (may be empty)
"""
function _load_dem_raster(dem_path::String)
    ArchGDAL.read(dem_path) do ds
        gt      = ArchGDAL.getgeotransform(ds)
        nx      = Int(ArchGDAL.width(ds))
        ny      = Int(ArchGDAL.height(ds))
        raw_nd  = ArchGDAL.getnodatavalue(ArchGDAL.getband(ds, 1))
        nodata_val = raw_nd === nothing ? NaN : Float64(raw_nd)
        proj    = ArchGDAL.getproj(ds)
        # Read full band into a (nx × ny) matrix — ArchGDAL returns (nx, ny)
        band = Float64.(ArchGDAL.read(ds, 1))
        (band=band, gt=gt, nx=nx, ny=ny, nodata_val=nodata_val, proj=proj)
    end
end

"""
    _reproject_to_dem(lons, lats, proj) → (native_xs, native_ys)

Reproject lon/lat (EPSG:4326) coordinate vectors into the DEM's native CRS
described by `proj` (WKT string).  Returns two Float64 vectors of equal
length.  If `proj` is empty, returns the inputs unchanged (assumes WGS84).
"""
function _reproject_to_dem(lons::Vector{Float64}, lats::Vector{Float64},
                            proj::String)
    if isempty(proj)
        return copy(lons), copy(lats)
    end
    xs = copy(lons)
    ys = copy(lats)
    zs = zeros(length(xs))
    wgs84   = _crs_gis(4326)
    dem_crs = _crs_gis(proj)
    ArchGDAL.createcoordtrans(wgs84, dem_crs) do xform
        ArchGDAL.transform!(xs, ys, zs, xform)
    end
    xs, ys
end

"""
    _affine_to_pixel(px, py, gt) → (col_f, row_f)

Convert a native-CRS coordinate (px, py) to fractional 0-based pixel
indices using the inverse of the raster's affine geotransform.
"""
@inline function _affine_to_pixel(px::Float64, py::Float64, gt::Vector{Float64})
    # gt = [ulx, xres, xskew, uly, yskew, yres]
    if gt[3] ≈ 0.0 && gt[5] ≈ 0.0
        # Common case — no rotation/skew
        return (px - gt[1]) / gt[2], (py - gt[4]) / gt[6]
    else
        det   = gt[2] * gt[6] - gt[3] * gt[5]
        dx    = px - gt[1];  dy = py - gt[4]
        return (dx * gt[6] - dy * gt[3]) / det,
               (dy * gt[2] - dx * gt[5]) / det
    end
end

"""
    _bilinear(band, col_f, row_f, nx, ny, nodata_val) → Float64

Bilinear interpolation from an in-memory raster array at fractional pixel
(col_f, row_f) (0-based, *corner*-based as `_affine_to_pixel` returns them:
pixel (c, r) covers [c, c+1) × [r, r+1)). A pixel's value is its value at
its centre (c + ½, r + ½) — GDAL's pixel-is-area convention; GDAL also
reports AREA_OR_POINT = Point rasters with a corner-based geotransform — so
the interpolation runs between pixel centres, and within half a pixel of the
raster's edge the value is constant across the edge (nearest centre).

Returns NaN outside the raster. If any of the 2×2 centres is nodata, the
value of the pixel containing the point (NaN if that is nodata).

0.8.0 (⚠ moves every DEM-sampled elevation): until 0.7.1 the corner-based
coordinates were used as centre-based ones, so every sample was taken half a
pixel east and south of its point (north-up rasters).
"""
@inline function _bilinear(band::Matrix{Float64},
                            col_f::Float64, row_f::Float64,
                            nx::Int, ny::Int,
                            nodata_val::Float64)::Float64
    (col_f < 0 || col_f >= nx || row_f < 0 || row_f >= ny) && return NaN

    nd = nodata_val
    isnd(v) = !isnan(nd) && v == nd
    if nx < 2 || ny < 2
        v = band[floor(Int, col_f) + 1, floor(Int, row_f) + 1]
        return isnd(v) ? NaN : v
    end

    u = col_f - 0.5; w = row_f - 0.5            # centre-based coordinates
    c0 = clamp(floor(Int, u), 0, nx - 2)
    r0 = clamp(floor(Int, w), 0, ny - 2)
    fc = clamp(u - c0, 0.0, 1.0)
    fr = clamp(w - r0, 0.0, 1.0)

    # ArchGDAL band array: band[col+1, row+1] (1-based, column-major)
    v00 = band[c0+1, r0+1];  v10 = band[c0+2, r0+1]
    v01 = band[c0+1, r0+2];  v11 = band[c0+2, r0+2]

    if isnd(v00) || isnd(v10) || isnd(v01) || isnd(v11)
        # The pixel containing the point — avoids propagating nodata into
        # an otherwise valid interpolation neighbourhood
        v = band[floor(Int, col_f) + 1, floor(Int, row_f) + 1]
        return isnd(v) ? NaN : v
    end

    return (1.0 - fr) * ((1.0 - fc) * v00 + fc * v10) +
                   fr  * ((1.0 - fc) * v01 + fc * v11)
end

"""
    _store_elevation!(mesh, elevations, var_name)

Write an elevation vector into `mesh.static_vars[var_name]` and, when
`var_name == "elevation"`, rebuild the cells vector with updated fields.
"""
function _store_elevation!(mesh::A5Mesh, elevations::Vector{Float64},
                            var_name::String)
    mesh.static_vars[var_name] = elevations
    if var_name == "elevation"
        # 0.8.0: elevations sampled at pixel centres (see `_bilinear`); meshes
        # without this column were sampled half a pixel off if from a DEM
        mesh.static_vars["dem_registration"] = fill(1.0, length(elevations))
        mesh.cells = [A5Cell(c.id, c.resolution, c.center_lon, c.center_lat,
                              c.boundary, elevations[i])
                      for (i, c) in enumerate(mesh.cells)]
    end
end

# ---------------------------------------------------------------------------
# Halton sequence generator
# ---------------------------------------------------------------------------

"""
    _halton(base, n, offset=0) → Vector{Float64}

Generate `n` values of the 1-D Halton low-discrepancy sequence in base
`base`, starting at position `offset` (0-indexed).  Values are in [0, 1).

Using bases 2 and 3 gives a 2-D Halton sequence with good space-filling
properties and no axis-aligned aliasing (unlike regular grids).

The `offset` parameter shifts the starting position in the sequence, enabling:
  - Reproducible runs (same offset → same points)
  - Monte Carlo uncertainty estimation (different offsets → independent samples)
"""
function _halton(base::Int, n::Int, offset::Int=0)::Vector{Float64}
    out = Vector{Float64}(undef, n)
    for i in 1:n
        idx  = i + offset
        f    = 1.0
        r    = 0.0
        while idx > 0
            f   = f / base
            r   = r + f * (idx % base)
            idx = idx ÷ base
        end
        out[i] = r
    end
    return out
end

# ---------------------------------------------------------------------------
# DEM sampling — centroid method
# ---------------------------------------------------------------------------

"""
    sample_dem_centroid!(mesh, dem_source; strict=false, var_name="elevation")

Sample a DEM at each cell's centre point and store in `mesh.static_vars`.

This is the fast method (~milliseconds for a typical mesh).  It uses bilinear
interpolation from an in-memory raster array — one GDAL read for the whole
raster, then pure Julia arithmetic per cell.

For a physically-motivated mean elevation (better for storage-cell flood
models), use `sample_dem_mean!` instead.

See `sample_dem_mean!` for argument documentation.
"""
function sample_dem_centroid!(mesh::A5Mesh, dem_source::DEMSource;
                               strict::Bool=false, var_name::String="elevation")
    dem_path = _resolve_dem_path(dem_source, mesh)
    n = length(mesh.cells)
    n == 0 && error("Cannot sample DEM onto an empty mesh.")

    @info "DEM sampling (centroid): $(basename(dem_path)) → $n cells..."
    t0 = time()

    raster = _load_dem_raster(dem_path)

    # Reproject all cell centres in one batch
    lons = [c.center_lon for c in mesh.cells]
    lats = [c.center_lat for c in mesh.cells]
    native_xs, native_ys = _reproject_to_dem(lons, lats, raster.proj)

    elevs  = Vector{Float64}(undef, n)
    n_oob  = 0; n_nd = 0
    for i in 1:n
        col_f, row_f = _affine_to_pixel(native_xs[i], native_ys[i], raster.gt)
        v = _bilinear(raster.band, col_f, row_f, raster.nx, raster.ny,
                      raster.nodata_val)
        if col_f < 0 || col_f >= raster.nx || row_f < 0 || row_f >= raster.ny
            n_oob += 1
        elseif isnan(v)
            n_nd += 1
        end
        elevs[i] = v
    end

    if n_oob > 0
        msg = "$n_oob cell centre(s) outside DEM extent"
        strict ? error("$msg. Re-run without --dem-strict to assign NaN instead.") :
                 @warn "$msg — NaN assigned for those cells."
    end

    n_valid = count(isfinite, elevs)
    @info "DEM centroid sampling complete in $(round(time()-t0, digits=1))s: " *
          "$n_valid / $n cells with valid elevation" *
          (n_oob > 0 ? "  ⚠ $n_oob outside extent" : "") *
          (n_nd  > 0 ? "  ⚠ $n_nd nodata"           : "")

    _store_elevation!(mesh, elevs, var_name)
    return mesh
end

sample_dem_centroid!(mesh::A5Mesh, dem_path::String; kwargs...) =
    sample_dem_centroid!(mesh, FileDEM(dem_path); kwargs...)

# Keep the old name as a deprecated alias so existing code doesn't break
const sample_dem! = sample_dem_centroid!

# ---------------------------------------------------------------------------
# DEM sampling — cell-mean method (Halton quasi-random sampling)
# ---------------------------------------------------------------------------

"""
    sample_dem_mean!(mesh, dem_source;
                     n_samples=256, halton_seed=0,
                     strict=false, var_name="elevation")

Sample a DEM at `n_samples` quasi-random points within each cell's polygon
footprint and store the mean of valid samples in `mesh.static_vars[var_name]`.

This gives a physically better representative elevation for storage-cell
flood models than the centroid, particularly for sloped terrain or cells
that straddle ridges and channels.

Method
------
1. For each cell, generate `n_samples` candidate points in the cell's
   lon/lat bounding box using a 2-D Halton sequence (bases 2 and 3).
2. Filter to points inside the cell polygon using the existing PIP
   infrastructure (`_pip_mask_gpu` if CUDA available, else `_pip_mask_cpu`
   with all threads).
3. Reproject surviving points from EPSG:4326 to the DEM's native CRS in
   one vectorised batch (ArchGDAL `createcoordtrans`).
4. For each point, apply the inverse affine geotransform to get fractional
   pixel indices, then bilinear-interpolate from the in-memory raster array.
5. Average the finite (non-NaN, non-nodata) values for the cell.
   Cells with zero valid samples are assigned NaN.

Arguments
---------
  dem_source   — a `DEMSource` (e.g. `FileDEM("/path/to/dem.tif")`) or path
  n_samples    — Halton candidates per cell before PIP filtering (default 256)
  halton_seed  — offset into the Halton sequence (default 0).
                 Changing this gives independent sample sets for Monte Carlo
                 uncertainty estimation.
  strict       — error if any cell has zero valid samples (default false)
  var_name     — key in `mesh.static_vars` (default "elevation")

Performance
-----------
The raster is loaded into memory once.  PIP tests use the GPU kernel if
CUDA is available and total candidate points exceed CUDA_PIP_THRESHOLD,
otherwise Julia CPU threads.  All bilinear sampling is pure Julia with no
GDAL calls in the inner loop.
"""
function sample_dem_mean!(mesh::A5Mesh, dem_source::DEMSource;
                           n_samples::Int=256,
                           halton_seed::Int=0,
                           strict::Bool=false,
                           var_name::String="elevation")
    dem_path = _resolve_dem_path(dem_source, mesh)
    n = length(mesh.cells)
    n == 0 && error("Cannot sample DEM onto an empty mesh.")

    @info "DEM sampling (mean, $(n_samples) Halton pts/cell, seed=$(halton_seed)): " *
          "$(basename(dem_path)) → $n cells..."
    t0 = time()

    raster = _load_dem_raster(dem_path)

    # Pre-generate Halton sequences once (reused for every cell via bbox scaling)
    h2 = _halton(2, n_samples, halton_seed)   # X offsets in [0,1)
    h3 = _halton(3, n_samples, halton_seed)   # Y offsets in [0,1)

    elevs   = Vector{Float64}(undef, n)
    n_oob   = 0
    n_nd    = 0

    # We process cells in parallel on CPU threads.
    # The raster array is read-only so no locking needed.
    Threads.@threads for ci in 1:n
        cell = mesh.cells[ci]
        bnd  = cell.boundary   # Vector{Vector{Float64}} — lon/lat ring

        # 1. Bounding box of the cell polygon in lon/lat
        bb_lon_min = minimum(Float64(v[1]) for v in bnd)
        bb_lon_max = maximum(Float64(v[1]) for v in bnd)
        bb_lat_min = minimum(Float64(v[2]) for v in bnd)
        bb_lat_max = maximum(Float64(v[2]) for v in bnd)
        dlon = bb_lon_max - bb_lon_min
        dlat = bb_lat_max - bb_lat_min

        # 2. Generate candidate points in the bbox
        cand_lons = bb_lon_min .+ h2 .* dlon
        cand_lats = bb_lat_min .+ h3 .* dlat

        # 3. PIP filter — keep only points inside the cell polygon.
        #    We call _point_in_ring directly (cheaper than the full
        #    _pip_mask dispatcher which has GPU overhead for small n).
        #    For large meshes the outer Threads.@threads loop already
        #    saturates all CPU cores; GPU PIP is not beneficial here
        #    because each cell individually has too few points.
        inside = [_point_in_ring(cand_lons[j], cand_lats[j], bnd)
                  for j in 1:n_samples]
        in_lons = cand_lons[inside]
        in_lats = cand_lats[inside]
        n_in = length(in_lons)

        if n_in == 0
            # Extremely unlikely for a valid cell — fall back to centroid
            elevs[ci] = NaN
            n_nd += 1
            continue
        end

        # 4. Reproject surviving points to DEM native CRS
        native_xs, native_ys = if isempty(raster.proj)
            copy(in_lons), copy(in_lats)
        else
            xs = copy(in_lons);  ys = copy(in_lats);  zs = zeros(n_in)
            wgs84   = _crs_gis(4326)
            dem_crs = _crs_gis(raster.proj)
            ArchGDAL.createcoordtrans(wgs84, dem_crs) do xform
                ArchGDAL.transform!(xs, ys, zs, xform)
            end
            xs, ys
        end

        # 5. Sample raster and accumulate mean
        s    = 0.0
        nval = 0
        for j in 1:n_in
            col_f, row_f = _affine_to_pixel(native_xs[j], native_ys[j], raster.gt)
            v = _bilinear(raster.band, col_f, row_f, raster.nx, raster.ny,
                          raster.nodata_val)
            if isfinite(v)
                s    += v
                nval += 1
            end
        end

        if nval == 0
            elevs[ci] = NaN
            n_nd += 1
        else
            elevs[ci] = s / nval
        end
    end

    if n_nd > 0
        msg = "$n_nd cell(s) with no valid DEM samples"
        strict ? error("$msg.") : @warn "$msg — NaN assigned."
    end

    n_valid = count(isfinite, elevs)
    @info "DEM mean sampling complete in $(round(time()-t0, digits=1))s: " *
          "$n_valid / $n cells with valid elevation" *
          (n_oob > 0 ? "  ⚠ $n_oob outside extent" : "") *
          (n_nd  > 0 ? "  ⚠ $n_nd no valid samples" : "")

    _store_elevation!(mesh, elevs, var_name)
    return mesh
end

sample_dem_mean!(mesh::A5Mesh, dem_path::String; kwargs...) =
    sample_dem_mean!(mesh, FileDEM(dem_path); kwargs...)


# ---------------------------------------------------------------------------
# PyCall setup — pya5 for individual cell queries
# ---------------------------------------------------------------------------

const _a5 = PyNULL()
const _PYA5_VERSION = Ref{String}("")
const _A5_MAJOR     = Ref{Int}(0)

function __init__()
    _module_init_done[] && return
    _module_init_done[] = true
    try
        copy!(_a5, pyimport("a5"))
    catch e
        error("Failed to import pya5: $e\nRun: pip install pya5")
    end
    # pya5 version → one code path per major version (WP-M1). `a5` has no
    # __version__; `a5.migrate` (new in 1.0) is the fallback signal.
    try
        _PYA5_VERSION[] = String(pyimport("importlib.metadata").version("pya5"))
    catch
        _PYA5_VERSION[] = "unknown"
    end
    m = tryparse(Int, first(split(_PYA5_VERSION[], '.')))
    _A5_MAJOR[] = m !== nothing ? m : (PyCall.hasproperty(_a5, "migrate") ? 1 : 0)
    # CUDA is loaded on first use (_gpu_available), not here (WP4).
end

"""    pya5_version() → String   (e.g. "0.7.3", "1.0.0") """
pya5_version() = _PYA5_VERSION[]

"""
    a5_id_scheme() → String

Cell-ID scheme of the installed pya5: "a5-v1" for 1.x, "a5-v0" for 0.x. See
`A5Mesh.id_scheme` and mesh/a5_bridge.py ("pya5 versions").
"""
a5_id_scheme() = _A5_MAJOR[] >= 1 ? "a5-v1" : "a5-v0"

"""
    check_id_scheme(mesh_scheme, source) → Nothing

Refuse to use a mesh whose A5 cell-ID scheme differs from the installed
pya5's: an "a5-v0" ID passed to pya5 1.x (or the reverse) is still a valid
cell ID, just of a different cell, so mixing them fails silently. Non-A5
schemes ("cartesian") are not checked. Set ENV FLOODA5_ALLOW_ID_SCHEME_MISMATCH=1
to downgrade the error to a warning (e.g. to inspect an old mesh); a
simulation does not call pya5 with the cell IDs, but mesh generation,
re-sampling and grid_disk adjacency do.
"""
function check_id_scheme(mesh_scheme::AbstractString, source::AbstractString)
    startswith(mesh_scheme, "a5-") || return nothing
    installed = a5_id_scheme()
    mesh_scheme == installed && return nothing
    hint = mesh_scheme == "a5-v0" ?
        "Migrate it with:  python mesh/migrate_mesh_ids.py <in.parquet> <out.parquet>  (needs pya5 >= 1.0)" :
        "It was written with pya5 >= 1.0; upgrade pya5 to use it."
    msg = "$source uses A5 cell-ID scheme '$mesh_scheme' but the installed pya5 " *
          "$(pya5_version()) uses '$installed'. $hint"
    if get(ENV, "FLOODA5_ALLOW_ID_SCHEME_MISMATCH", "") == "1"
        @warn msg * " (continuing: FLOODA5_ALLOW_ID_SCHEME_MISMATCH=1)"
    else
        error(msg)
    end
    return nothing
end

# GIL lock — all PyCall calls into pya5 must hold this lock.
# @py is a module-local macro (NOT from PyCall) that wraps expressions in the lock.
const _pycall_lock = ReentrantLock()
macro py(expr)
    return :(lock(_pycall_lock) do; $(esc(expr)); end)
end

# ---------------------------------------------------------------------------
# Individual cell queries via PyCall (for flow model use)
# ---------------------------------------------------------------------------

# pya5 uses signed int64 internally (struct.pack '>q').
# Pass UInt64 as reinterpreted Int64; receive back as Int64 or BigInt from PyCall.
@inline _to_pyint(u::UInt64)  = PyObject(reinterpret(Int64, u))
@inline _hex_to_pyint(h::String) = _to_pyint(parse(UInt64, h, base=16))

# BigInt mask handles both Int64 (possibly negative) and BigInt returns from PyCall
@inline function _to_hex(cell_id)::String
    u = UInt64(BigInt(cell_id) & BigInt(0xffffffffffffffff))
    return string(u, base=16, pad=16)
end

# Like _to_hex, but for contexts (mesh re-save/adjacency preservation) that
# receive an already-string cell ID which may or may not be an A5 hex ID —
# e.g. a Cartesian or other future backend's own ID scheme. Falls back to
# the raw string rather than throwing, so non-A5 backends round-trip through
# save_mesh_geoparquet/load_mesh_geoparquet without needing to imitate A5's
# hex-ID convention. Mirrors FloodModel.jl's `_norm_cell_id` (same rationale,
# duplicated here rather than shared because A5Grid is a self-contained
# module included before FloodModel.jl's own definitions exist).
@inline function _to_hex_safe(id::AbstractString)::String
    try
        return _to_hex(parse(UInt64, id, base=16))
    catch
        return String(id)
    end
end

"""    lonlat_to_cell(lon, lat, resolution) → String """
function lonlat_to_cell(lon::Real, lat::Real, resolution::Int)::String
    cell_id = @py _a5.lonlat_to_cell((Float64(lon), Float64(lat)), resolution)
    return _to_hex(cell_id)
end

"""    cell_to_lonlat(cell_id_hex) → Tuple{Float64,Float64} """
function cell_to_lonlat(cell_id_hex::String)::Tuple{Float64,Float64}
    r = @py _a5.cell_to_lonlat(_hex_to_pyint(cell_id_hex))
    return (Float64(r[1]), Float64(r[2]))
end

"""    cell_to_boundary(cell_id_hex) → Vector{Vector{Float64}} """
function cell_to_boundary(cell_id_hex::String)::Vector{Vector{Float64}}
    raw = @py _a5.cell_to_boundary(_hex_to_pyint(cell_id_hex))
    return [[Float64(v[1]), Float64(v[2])] for v in raw]
end

"""    cell_to_children(cell_id_hex, resolution) → Vector{String} """
function cell_to_children(cell_id_hex::String, resolution::Int)::Vector{String}
    ch = @py _a5.cell_to_children(_hex_to_pyint(cell_id_hex), resolution)
    return [_to_hex(c) for c in ch]
end

"""    cell_to_parent(cell_id_hex) → String """
function cell_to_parent(cell_id_hex::String)::String
    return _to_hex(@py _a5.cell_to_parent(_hex_to_pyint(cell_id_hex)))
end

"""    get_resolution(cell_id_hex) → Int """
function get_resolution(cell_id_hex::String)::Int
    return Int(@py _a5.get_resolution(_hex_to_pyint(cell_id_hex)))
end

"""
    grid_disk_neighbours(cell_id_hex) → Vector{String}

Return the topological neighbours of an A5 pentagon cell at the same
resolution level.

`pya5.grid_disk(cell, 1)` returns a *compact* set — cells may be at
coarser resolution levels than the input.  We call `pya5.uncompact` to
expand the disk to the target resolution, then exclude the centre cell.
This guarantees all returned neighbours are at the same resolution as
the input cell, giving true same-level edge-sharing adjacency.
"""
function grid_disk_neighbours(cell_id_hex::String)::Vector{String}
    py_id      = _hex_to_pyint(cell_id_hex)
    resolution = Int(@py _a5.get_resolution(py_id))
    disk       = @py _a5.grid_disk(py_id, 1)
    expanded = if _A5_MAJOR[] >= 1
        # pya5 1.x: uncompact() has no target-resolution argument (a compact
        # set carries a compaction marker instead), so drop the markers and
        # keep the cells at this resolution.
        [c for c in (@py _a5.uncompact(disk))
         if !(@py _a5.is_compaction_marker(c)) && Int(@py _a5.get_resolution(c)) == resolution]
    else
        @py _a5.uncompact(disk, resolution)
    end
    return [_to_hex(c) for c in expanded if _to_hex(c) != cell_id_hex]
end

"""
    migrate_cell_id(cell_id_hex) → String

pya5 0.x ("a5-v0") cell ID → the same cell's pya5 1.x ID (`a5.migrate`).
Requires pya5 >= 1.0. For whole meshes use mesh/migrate_mesh_ids.py.
"""
function migrate_cell_id(cell_id_hex::String)::String
    _A5_MAJOR[] >= 1 || error("migrate_cell_id needs pya5 >= 1.0 (installed: $(pya5_version()))")
    return _to_hex(@py _a5.migrate(_hex_to_pyint(cell_id_hex)))
end

"""
    grid_disk_neighbours_batch(cell_ids) → Dict{String, Vector{String}}

Batch version of grid_disk_neighbours. Calls pya5 once per cell but
releases the GIL lock between cells so Julia threads can overlap.
Returns a Dict mapping each cell_id to its neighbour list.
"""
function grid_disk_neighbours_batch(cell_ids::Vector{String})::Dict{String,Vector{String}}
    adj = Dict{String,Vector{String}}()
    for id in cell_ids
        adj[id] = grid_disk_neighbours(id)
    end
    return adj
end

# ---------------------------------------------------------------------------
# AOI helpers & geometry (pure Julia — used for GPU/CPU PIP sampling)
# ---------------------------------------------------------------------------

function load_aoi(path_or_str::String)::String
    isfile(path_or_str) && return read(path_or_str, String)
    startswith(strip(path_or_str), "{") && return path_or_str
    error("AOI must be a GeoJSON file path or GeoJSON string.")
end

@inline function _point_in_ring(lon::Float64, lat::Float64, ring)::Bool
    inside = false; n = length(ring); j = n
    for i in 1:n
        xi, yi = Float64(ring[i][1]), Float64(ring[i][2])
        xj, yj = Float64(ring[j][1]), Float64(ring[j][2])
        if ((yi > lat) != (yj > lat)) &&
           (lon < (xj - xi) * (lat - yi) / (yj - yi + 1e-15) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

function _point_in_geometry(lon::Float64, lat::Float64, geometry::Dict)::Bool
    t = geometry["type"]
    t == "Polygon" && return _point_in_ring(lon, lat, geometry["coordinates"][1])
    t == "MultiPolygon" && return any(
        _point_in_ring(lon, lat, poly[1]) for poly in geometry["coordinates"])
    error("Unsupported geometry type: $t")
end

function _geometry_bbox(geometry::Dict)
    t = geometry["type"]
    rings = t == "Polygon" ? geometry["coordinates"] :
            [r for poly in geometry["coordinates"] for r in poly]
    coords = [c for ring in rings for c in ring]
    return (minimum(Float64(c[1]) for c in coords), minimum(Float64(c[2]) for c in coords),
            maximum(Float64(c[1]) for c in coords), maximum(Float64(c[2]) for c in coords))
end

function _aoi_geometry(geojson::Dict)::Dict
    t = geojson["type"]
    t == "Feature"      && return Dict(geojson["geometry"])
    t == "Polygon"      && return geojson
    t == "MultiPolygon" && return geojson
    if t == "FeatureCollection"
        feats = geojson["features"]
        length(feats) == 1 && return Dict(feats[1]["geometry"])
        polys = []
        for f in feats
            g = f["geometry"]
            g["type"] == "Polygon"      && push!(polys, g["coordinates"])
            g["type"] == "MultiPolygon" && append!(polys, g["coordinates"])
        end
        return Dict("type" => "MultiPolygon", "coordinates" => polys)
    end
    error("Cannot extract polygon geometry from GeoJSON type: $t")
end

# ---------------------------------------------------------------------------
# Mesh generation — PIP in Julia (GPU or CPU), cell indexing in Python
# ---------------------------------------------------------------------------

function _build_sample_grid(min_lon, min_lat, max_lon, max_lat, step)
    lons = Float64[]; lats = Float64[]
    lat = min_lat
    while lat <= max_lat + step
        lon = min_lon
        while lon <= max_lon + step
            push!(lons, clamp(lon, -179.9, 179.9))
            push!(lats, min(lat, 89.9))
            lon += step
        end
        lat += step
    end
    return lons, lats
end

function _pip_mask_cpu(lons, lats, geometry)::Vector{Bool}
    mask = Vector{Bool}(undef, length(lons))
    Threads.@threads for i in eachindex(lons)
        mask[i] = _point_in_geometry(lons[i], lats[i], geometry)
    end
    return mask
end

# Minimum point count before GPU is worth the PCIe transfer overhead.
const CUDA_PIP_THRESHOLD = Ref{Int}(500_000)

function _pip_mask(lons, lats, geometry)::Vector{Bool}
    use_gpu = length(lons) >= CUDA_PIP_THRESHOLD[] && _gpu_available()
    if use_gpu
        t    = geometry["type"]
        ring = t == "Polygon" ? geometry["coordinates"][1] :
                                geometry["coordinates"][1][1]
        rxs  = Float64[Float64(v[1]) for v in ring]
        rys  = Float64[Float64(v[2]) for v in ring]
        return Base.invokelatest(_pip_mask_gpu, lons, lats, rxs, rys)
    else
        return _pip_mask_cpu(lons, lats, geometry)
    end
end

"""
    mesh_for_aoi(geojson_path_or_str, resolution) → A5Mesh

Generate a pentagonal A5 mesh covering the AOI.
"""
function mesh_for_aoi(geojson_path_or_str::String, resolution::Int;
                      output_path::Union{String,Nothing}=nothing,
                      format::Symbol=:geoparquet)::A5Mesh

    geojson_str = load_aoi(geojson_path_or_str)
    geojson     = JSON3.read(geojson_str, Dict)
    geometry    = _aoi_geometry(geojson)

    min_lon, min_lat, max_lon, max_lat = _geometry_bbox(geometry)
    approx_cell_deg = max(0.0001, 90.0 / (2.0^(resolution * 1.16)))
    sample_step     = approx_cell_deg * 0.45

    @info "Building sample grid (resolution=$resolution, step=$(round(sample_step,sigdigits=4))°)..."
    lons, lats = _build_sample_grid(min_lon, min_lat, max_lon, max_lat, sample_step)
    use_gpu = length(lons) >= CUDA_PIP_THRESHOLD[] && _gpu_available()
    backend = use_gpu ? "GPU" : "CPU ($(Threads.nthreads()) threads)"
    @info "  Sample points: $(length(lons))  |  Backend: $backend (threshold=$(CUDA_PIP_THRESHOLD[]) pts)"

    mask      = _pip_mask(lons, lats, geometry)
    kept_lons = lons[mask]
    kept_lats = lats[mask]
    @info "  Points inside AOI: $(length(kept_lons)) — handing off to Python..."

    # unique names: several mesh generations may run in one checkout at once
    tag     = _tmp_tag()
    ext     = format == :geoparquet ? ".parquet" : ".geojson"
    aoi_nm  = "_tmp_aoi_$tag.geojson"
    out_nm  = "_tmp_mesh_$tag" * ext
    tmp_aoi = joinpath(BRIDGE_DIR, aoi_nm)
    tmp_out = joinpath(BRIDGE_DIR, out_nm)
    fmt_str = format == :geoparquet ? "geoparquet" : "geojson"

    try
        write(tmp_aoi, geojson_str)
        result = _run_bridge("mesh_for_aoi", aoi_nm, string(resolution), out_nm, fmt_str)
        @info "  Python bridge: $(result["cell_count"]) cells written to $(result["output"])"

        mesh = load_mesh_from_file(tmp_out, geojson_str, resolution)

        if output_path !== nothing
            cp(tmp_out, output_path; force=true)
            endswith(output_path, ".parquet") && _write_prj_sidecar(output_path)
        end

        return mesh
    finally
        isfile(tmp_aoi) && rm(tmp_aoi)
        isfile(tmp_out) && rm(tmp_out)
    end
end

# ---------------------------------------------------------------------------
# Serialisation — GeoJSON string (for visualisation server)
# ---------------------------------------------------------------------------

"""
    mesh_to_geojson_string(mesh) → String

Serialise an A5Mesh to a GeoJSON FeatureCollection string.
Includes elevation in properties when present (for Cesium depth visualisation).
"""
function mesh_to_geojson_string(mesh::A5Mesh)::String
    has_elev = haskey(mesh.static_vars, "elevation")
    elevations = has_elev ? mesh.static_vars["elevation"] : nothing

    features = [Dict(
        "type"     => "Feature",
        "id"       => cell.id,
        "geometry" => Dict(
            "type"        => "Polygon",
            "coordinates" => [cell.boundary]
        ),
        "properties" => merge(
            Dict(
                "cell_id"    => cell.id,
                "resolution" => cell.resolution,
                "center_lon" => cell.center_lon,
                "center_lat" => cell.center_lat,
            ),
            has_elev ? Dict("elevation" => (isfinite(elevations[i]) ? elevations[i] : nothing)) : Dict()
        )
    ) for (i, cell) in enumerate(mesh.cells)]

    fc = Dict(
        "type"     => "FeatureCollection",
        "features" => features,
        "properties" => Dict(
            "resolution"   => mesh.resolution,
            "cell_count"   => length(mesh.cells),
            "has_elevation" => has_elev,
            "static_vars"  => collect(keys(mesh.static_vars)),
        )
    )
    return JSON3.write(fc)
end

# ---------------------------------------------------------------------------
# Read mesh from GeoParquet or GeoJSON
# ---------------------------------------------------------------------------

function load_mesh_from_file(path::String, aoi_geojson::String, resolution::Int)::A5Mesh
    if endswith(path, ".parquet")
        return _load_from_geoparquet(path, aoi_geojson, resolution)
    else
        return _load_from_geojson(path, aoi_geojson, resolution)
    end
end

function _load_from_geojson(path::String, aoi_geojson::String, resolution::Int)::A5Mesh
    raw = JSON3.read(read(path, String), Dict)
    cells = A5Cell[]
    for feat in raw["features"]
        props = feat["properties"]
        coords = feat["geometry"]["coordinates"][1]
        boundary = [[Float64(v[1]), Float64(v[2])] for v in coords]
        push!(cells, A5Cell(
            String(props["cell_id"]), Int(props["resolution"]),
            Float64(props["center_lon"]), Float64(props["center_lat"]),
            boundary))
    end
    return A5Mesh(resolution, aoi_geojson, cells)
end

"""
    _load_from_geoparquet(path, aoi_geojson, resolution) → A5Mesh

Load mesh from GeoParquet. Reads all extra columns beyond the base set as
static variables, so elevation (and future fields like friction, porosity)
are automatically restored from a previously saved mesh.
"""
function _load_from_geoparquet(path::String, aoi_geojson::String, resolution::Int)::A5Mesh
    @ensure_python
    # Resolve to absolute path before handing to the Python subprocess.
    # The bridge runs with dir=BRIDGE_DIR (mesh/), so relative paths would be
    # resolved from there rather than from the Julia working directory.
    #
    # WP4 (2026-10): the columns come back as raw Float64 files plus a small
    # JSON manifest, read with vectorised numpy/shapely — previously every
    # cell (≈1,300 numbers per SGS cell) went through iterrows and JSON text,
    # ~75 s for Carlisle res-18 SGS. The values are identical: array columns
    # map non-finite values to NaN, scalar columns map NaN/None to NaN (±Inf
    # kept), as before.
    path_fwd = replace(abspath(path), "\\" => "/")
    out_dir  = mktempdir()
    out_fwd  = replace(out_dir, "\\" => "/")
    py_code = """
import geopandas as gpd, json, math, os
import numpy as np
import shapely
gdf = gpd.read_parquet('$(path_fwd)')
out = '$(out_fwd)'

import pyarrow.parquet as _pq
_raw_md = (_pq.read_schema('$(path_fwd)').metadata or {}).get(b'flooda5')
mesh_meta = json.loads(_raw_md.decode()) if _raw_md else {}
id_scheme = mesh_meta.get('a5_id_scheme')

BASE_COLS = {'cell_id', 'resolution', 'center_lon', 'center_lat', 'geometry'}
extra_cols = [c for c in gdf.columns if c not in BASE_COLS]
scalar_cols, array_cols, strlist_cols = [], [], []
for col in extra_cols:
    sample = gdf[col].iloc[0] if len(gdf) > 0 else None
    if hasattr(sample, '__len__') and not isinstance(sample, str):
        try:
            first_elem = list(sample)[0] if len(sample) > 0 else None
            if isinstance(first_elem, str):
                strlist_cols.append(col)
            else:
                array_cols.append(col)
        except Exception:
            array_cols.append(col)
    else:
        scalar_cols.append(col)

n = len(gdf)
def wr(name, a):
    np.ascontiguousarray(a, dtype=np.float64).tofile(os.path.join(out, name))

wr('center.bin', np.column_stack([gdf['center_lon'].astype(float).to_numpy(),
                                  gdf['center_lat'].astype(float).to_numpy()]))
coords, idx = shapely.get_coordinates(gdf.geometry.exterior.values, return_index=True)
wr('coords.bin', coords)
counts = np.bincount(idx, minlength=n).astype(np.float64)
wr('counts.bin', counts)

for k, col in enumerate(scalar_cols):
    v = gdf[col].to_numpy(dtype=object)
    a = np.array([np.nan if x is None else float(x) for x in v], dtype=np.float64)
    wr('s%d.bin' % k, a)

array_lens = []
for k, col in enumerate(array_cols):
    v = gdf[col].to_numpy(dtype=object)
    first = v[0] if n > 0 else None
    L = len(first) if (first is not None and len(first) > 0) else 0
    array_lens.append(L)
    if L == 0:
        continue
    m = np.full((n, L), np.nan)
    for i, x in enumerate(v):
        if x is None:
            continue
        a = np.asarray(x, dtype=np.float64)
        m[i, :len(a)] = a
    m[~np.isfinite(m)] = np.nan
    wr('a%d.bin' % k, m)

adjacency = {}
if 'neighbours' in strlist_cols:
    for cid, v in zip(gdf['cell_id'], gdf['neighbours']):
        if v is not None:
            adjacency[cid] = [str(x) for x in v]

with open(os.path.join(out, 'manifest.json'), 'w') as f:
    json.dump({'n': n, 'cell_id': [str(c) for c in gdf['cell_id']],
               'resolution': [int(r) for r in gdf['resolution']],
               'scalar_cols': scalar_cols, 'array_cols': array_cols,
               'array_lens': array_lens, 'adjacency': adjacency,
               'n_coords': int(len(coords)), 'id_scheme': id_scheme,
               'meta': mesh_meta}, f)
"""
    tmp_script = tempname() * ".py"
    write(tmp_script, py_code)
    errbuf = IOBuffer()
    raw = try
        run(pipeline(Cmd(Cmd([PYTHON_CMD[], tmp_script]); dir=BRIDGE_DIR),
                     stdout=devnull, stderr=errbuf))
        JSON3.read(read(joinpath(out_dir, "manifest.json"), String), Dict)
    catch e
        rm(out_dir; recursive = true, force = true)
        error("Failed to read GeoParquet: $(String(take!(errbuf)))")
    finally
        isfile(tmp_script) && rm(tmp_script)
    end

    _rd(name, dims...) = read!(joinpath(out_dir, name), Array{Float64}(undef, dims...))
    nc          = Int(raw["n"])
    scalar_cols = [String(k) for k in raw["scalar_cols"]]
    array_cols  = [String(k) for k in raw["array_cols"]]
    array_lens  = Int.(raw["array_lens"])
    raw_adj     = get(raw, "adjacency", Dict())
    ids         = raw["cell_id"]; ress = raw["resolution"]

    center = _rd("center.bin", 2, nc)
    coords = _rd("coords.bin", 2, Int(raw["n_coords"]))
    counts = round.(Int, _rd("counts.bin", nc))

    static_vars = Dict{String, Vector{Float64}}(
        col => _rd("s$(k-1).bin", nc) for (k, col) in enumerate(scalar_cols))
    array_vars = Dict{String, Matrix{Float64}}(
        col => _rd("a$(k-1).bin", array_lens[k], nc)
        for (k, col) in enumerate(array_cols) if array_lens[k] > 0)
    rm(out_dir; recursive = true, force = true)

    elev = get(static_vars, "elevation", nothing)
    cells = Vector{A5Cell}(undef, nc)
    off = 0
    for i in 1:nc
        bnd = [[coords[1, off + k], coords[2, off + k]] for k in 1:counts[i]]
        off += counts[i]
        cells[i] = A5Cell(String(ids[i]), Int(ress[i]), center[1, i], center[2, i], bnd,
                          elev === nothing ? NaN : elev[i])
    end

    res  = resolution > 0 ? resolution : (isempty(cells) ? 0 : cells[1].resolution)

    # Build adjacency dict from the neighbours column (if present in parquet)
    adjacency = Dict{String,Vector{String}}()
    for (cell_id, nbrs) in raw_adj
        adjacency[String(cell_id)] = [String(nb) for nb in nbrs]
    end
    if !isempty(adjacency)
        @info "Loaded adjacency for $(length(adjacency)) cells from parquet neighbours column"
    end

    # Cell-ID scheme: recorded in the file since WP-M1 (2026-10) by every
    # generator. A file without one predates that — a pya5 0.7.x A5 mesh
    # ("a5-v0") or a Cartesian mesh (resolution 0) — and is guessed, with a
    # warning: a mesh written by other software would be guessed wrongly.
    sch = get(raw, "id_scheme", nothing)
    id_scheme = if sch !== nothing
        String(sch)
    else
        guess = (!isempty(cells) && all(c -> c.resolution == 0, cells)) ? "cartesian" : "a5-v0"
        @warn "$path records no cell-ID scheme (written before 2026-10, or not by FloodA5): " *
              "assuming \"$guess\". Re-save it with FloodA5 to record one."
        guess
    end
    check_id_scheme(id_scheme, path)

    # The rest of the file's flooda5 metadata (DGGS family, level, …). The ID
    # scheme is the id_scheme field; pya5_version is rewritten on every save.
    meta = Dict{String,Any}(String(k) => (v isa AbstractString ? String(v) : v)
                            for (k, v) in get(raw, "meta", Dict())
                            if !(String(k) in ("a5_id_scheme", "pya5_version")))

    mesh = A5Mesh(res, aoi_geojson, cells, static_vars, array_vars, adjacency, id_scheme, meta)

    parts = String[]
    !isempty(scalar_cols) && push!(parts, "scalar: [$(join(scalar_cols, ", "))]")
    !isempty(array_cols)  && push!(parts, "arrays: [$(join(array_cols, ", "))]")
    !isempty(adjacency)   && push!(parts, "adjacency: $(length(adjacency)) cells")
    isempty(parts) || @info "Loaded mesh — $(join(parts, ", "))"
    return mesh
end

"""
    load_mesh_geoparquet(path) → A5Mesh

Load a previously saved mesh from a GeoParquet file.
Static variables (elevation, etc.) are restored automatically.
"""
function load_mesh_geoparquet(path::String)::A5Mesh
    return _load_from_geoparquet(path, "", 0)
end

# ---------------------------------------------------------------------------
# Save helpers
# ---------------------------------------------------------------------------

function save_mesh(mesh::A5Mesh, path::String; format::Symbol=:geoparquet)
    fmt = if format == :geoparquet || endswith(path, ".parquet")
        :geoparquet
    else
        :geojson
    end
    fmt == :geoparquet ? save_mesh_geoparquet(mesh, path) : save_mesh_geojson(mesh, path)
end

"""
    save_mesh_geoparquet(mesh, path)

Save mesh to GeoParquet. Scalar static variables (elevation, friction, etc.)
are written as Float64 columns.  Array-valued variables (SGS hypsometric
curves, edge sills) are written as Arrow list columns — each cell stores a
variable-length Float64 array in a single column, which GeoParquet/Arrow
supports natively.  Both types are automatically restored on load.
"""
function save_mesh_geoparquet(mesh::A5Mesh, path::String)
    @ensure_python
    mkpath(dirname(path) == "" ? "." : dirname(path))

    # Build per-cell dict list including scalar static_vars and array_vars
    # Build a cell_id → neighbours lookup for fast per-cell access
    adj_by_id = mesh.adjacency   # Dict{String,Vector{String}}, may be empty

    cells_data = []
    for (i, c) in enumerate(mesh.cells)
        d = Dict{String,Any}(
            "cell_id"    => c.id,
            "resolution" => c.resolution,
            "center_lon" => c.center_lon,
            "center_lat" => c.center_lat,
            "boundary"   => c.boundary,
        )
        for (varname, vals) in mesh.static_vars
            d[varname] = isfinite(vals[i]) ? vals[i] : nothing
        end
        # Array vars: column i of the matrix → list for this cell.
        # Replace NaN/Inf with nothing so JSON3 can serialise (null in JSON,
        # read back as NaN by the Python loader).
        for (varname, mat) in mesh.array_vars
            col = mat[:, i]
            d[varname] = [isfinite(v) ? v : nothing for v in col]
        end
        # Neighbours: preserve adjacency through re-saves (e.g. after DEM sampling).
        # Look up by both raw ID and normalised ID to handle padding differences.
        norm_id = _to_hex_safe(c.id)
        nbrs = get(adj_by_id, norm_id, get(adj_by_id, c.id, nothing))
        if nbrs !== nothing
            d["neighbours"] = nbrs
        end
        push!(cells_data, d)
    end

    tag        = _tmp_tag()       # unique: concurrent saves in one checkout
    tmp_cells  = joinpath(BRIDGE_DIR, "_tmp_save_cells_$tag.json")
    tmp_out    = joinpath(BRIDGE_DIR, "_tmp_save_$tag.parquet")
    tmp_meta   = joinpath(BRIDGE_DIR, "_tmp_save_meta_$tag.json")
    path_fwd   = replace(path, "\\" => "/")
    tmp_fwd    = replace(tmp_out, "\\" => "/")

    try
        open(tmp_cells, "w") do io
            JSON3.write(io, cells_data)
        end
        write(tmp_meta, JSON3.write(mesh.meta))

        tmp_cells_fwd = replace(tmp_cells, "\\" => "/")
        tmp_meta_fwd  = replace(tmp_meta, "\\" => "/")
        py_code = """
import sys, json
sys.argv = ['a5_bridge.py']
import importlib.util, os
spec = importlib.util.spec_from_file_location('a5_bridge', os.path.join(os.path.dirname('$(tmp_cells_fwd)'), 'a5_bridge.py'))
bridge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)

with open('$(tmp_cells_fwd)', 'r') as f:
    cells = json.load(f)
with open('$(tmp_meta_fwd)', 'r') as f:
    meta = json.load(f)
bridge.mesh_to_geoparquet(cells, '$(tmp_fwd)', id_scheme='$(mesh.id_scheme)', **meta)
print('ok')
"""
        tmp_script = joinpath(BRIDGE_DIR, "_tmp_save_$tag.py")
        write(tmp_script, py_code)
        buf = IOBuffer(); errbuf = IOBuffer()
        try
            run(pipeline(Cmd(Cmd([PYTHON_CMD[], basename(tmp_script)]); dir=BRIDGE_DIR),
                         stdout=buf, stderr=errbuf))
        catch e
            error("GeoParquet save failed: $(String(take!(errbuf)))")
        finally
            isfile(tmp_script) && rm(tmp_script)
        end

        cp(tmp_out, path; force=true)
        endswith(path, ".parquet") && _write_prj_sidecar(path)
        info_parts = String[]
        !isempty(mesh.static_vars) &&
            push!(info_parts, "static: [$(join(keys(mesh.static_vars), ", "))]")
        !isempty(mesh.array_vars) &&
            push!(info_parts, "arrays: [$(join(keys(mesh.array_vars), ", "))]")
        @info "Saved $(length(mesh.cells)) cells to $path (GeoParquet)" *
              (isempty(info_parts) ? "" : " — " * join(info_parts, ", "))
    finally
        isfile(tmp_cells) && rm(tmp_cells)
        isfile(tmp_meta)  && rm(tmp_meta)
        isfile(tmp_out)   && rm(tmp_out)
    end
end

function save_mesh_geojson(mesh::A5Mesh, path::String)
    mkpath(dirname(path) == "" ? "." : dirname(path))
    open(path, "w") do io
        JSON3.write(io, JSON3.read(mesh_to_geojson_string(mesh)))
    end
    @info "Saved $(length(mesh.cells)) cells to $path (GeoJSON)"
end

"""Write a WKT CRS sidecar .prj file alongside a GeoParquet mesh.
Needed because GDAL < 3.6 silently drops the PROJJSON CRS metadata
when opening .parquet files in QGIS."""
function _write_prj_sidecar(parquet_path::String)
    prj_path = parquet_path[1:findlast('.', parquet_path)-1] * ".prj"
    wkt = """GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433],AUTHORITY["EPSG","4326"]]"""
    write(prj_path, wkt)
    @info "Wrote CRS sidecar: $prj_path"
end


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

function diagnose()
    println("=== A5Grid diagnostics ===")
    println("BRIDGE_DIR    : ", BRIDGE_DIR)
    println("BRIDGE_SCRIPT : ", BRIDGE_SCRIPT)
    _find_python()
    println("PYTHON_CMD    : ", PYTHON_CMD[])
    println("CUDA_AVAILABLE: ", _gpu_available())
    println()
    println("Running bridge check...")
    try
        result = _run_bridge("check")
        println("  Python executable : ", get(result, "python", "unknown"))
        println("  pya5 available    : ", get(result, "pya5", false))
        println("  geopandas         : ", get(result, "geopandas", false))
        println("  geopandas path    : ", get(result, "geopandas_location", "not found"))
        println()
        if !get(result, "geopandas", false)
            println("  ⚠ geopandas not found. Fix:")
            println("""    run(`$(PYTHON_CMD[]) -m pip install geopandas pyarrow shapely`)""")
        else
            println("  ✓ Python dependencies found")
        end
        println()
        println("Julia DEM packages:")
        println("  ArchGDAL : ✓ loaded (hard dependency)")
    catch e
        println("  Bridge check failed: $e")
    end
    println("==========================")
end

function mesh_summary(mesh::A5Mesh)::String
    isempty(mesh.cells) && return "A5Mesh: empty (0 cells) at resolution $(mesh.resolution)"
    lons  = [c.center_lon for c in mesh.cells]
    lats  = [c.center_lat for c in mesh.cells]
    elvs  = get(mesh.static_vars, "elevation", nothing)
    elev_summary = if elvs !== nothing
        n_valid = count(isfinite, elvs)
        n_nan   = length(elvs) - n_valid
        valid   = filter(isfinite, elvs)
        "$(round(minimum(valid),digits=1)) – $(round(maximum(valid),digits=1)) m " *
        "  ($(n_valid) valid, $(n_nan) NaN)"
    else
        "not sampled"
    end
    static_keys = isempty(mesh.static_vars) ? "none" :
                  join(keys(mesh.static_vars), ", ")
    array_keys  = isempty(mesh.array_vars)  ? "none" :
                  join(("$k($(size(v,1))×$(size(v,2)))" for (k,v) in mesh.array_vars), ", ")
    has_sgs = haskey(mesh.array_vars, "sgs_elev_bins")
    return """
A5 Pentagon Mesh Summary
  Resolution   : $(mesh.resolution)
  Cell count   : $(length(mesh.cells))
  Lon range    : $(round(minimum(lons),digits=4)) → $(round(maximum(lons),digits=4))
  Lat range    : $(round(minimum(lats),digits=4)) → $(round(maximum(lats),digits=4))
  Elevation    : $elev_summary
  SGS tables   : $(has_sgs ? "present ($(size(mesh.array_vars["sgs_elev_bins"],1)) bins)" : "not built")
  Static vars  : $static_keys
  Array vars   : $array_keys
  Backend      : $(_cuda_available[] ? "GPU (CUDA, threshold=$(CUDA_PIP_THRESHOLD[]))" : "CPU ($(Threads.nthreads()) threads)")
"""
end

# ---------------------------------------------------------------------------
# Sub-Grid Sampling (SGS) — hypsometric curve pre-processing
# ---------------------------------------------------------------------------

"""
    SGSTable

Pre-computed hypsometric (area-elevation) relationship for a single A5 cell.

The table is built from the high-resolution DEM by binning all sample-point
elevations within the cell polygon into `n_bins` quantile levels.  From the
elevation histogram we derive cumulative volume and wetted-area curves as
functions of water surface elevation (WSE).

Fields
------
  elev_bins  — sorted elevation knot points (length n_bins), metres
  vol_curve  — cumulative storage volume at each WSE knot (m³)
  area_curve — cumulative wetted plan area at each WSE knot (m²)
  cell_area  — total plan area of the cell (m²), used for depth estimation
  z_min      — minimum terrain elevation in the cell (m) — the "thalweg"
  z_max      — maximum terrain elevation in the cell (m)

Usage
-----
Given a current stored volume V, look up the corresponding WSE via the
inverse of vol_curve (linear interpolation).  The wetted area A follows
directly.  Both lookups use Interpolations.jl for O(log n) binary search.
"""
struct SGSTable
    elev_bins        :: Vector{Float64}   # WSE knots (n_bins)
    vol_curve        :: Vector{Float64}   # cumulative volume  (n_bins) m³
    area_curve       :: Vector{Float64}   # cumulative wet area (n_bins) m²
    cell_area        :: Float64           # total cell plan area m²
    z_min            :: Float64           # minimum terrain elevation m
    z_max            :: Float64           # maximum terrain elevation m
    # Face tables (layout 3, 2026-11): the conveyance geometry of each face,
    # from the DEM profile sampled along the face line. Column s = adjacency
    # slot s (the s-th entry of the cell's neighbour list; unused slots hold
    # NaN / zeros). The two cells of a face hold the same table. Row k is the
    # level  face_z0[s] + face_zr[s]·((k−1)/(NF−1))²  (denser near the sill);
    # the levels reach above the highest point of the profile (see
    # `face_table_from_profile`), so the top level has the whole face wet:
    #   face_A[k, s]  flow area            Σ dx·max(η − z_p, 0)         (m²)
    #   face_P[k, s]  wetted width         Σ dx·[z_p ≤ η]               (m)
    #   face_K[k, s]  conveyance integral  Σ dx·max(η − z_p, 0)^(5/3)   (m^(8/3))
    #   face_C[k, s]  critical integral    Σ dx·max(η − z_p, 0)^(3/2)   (m^(5/2))
    #                 (layout 4, 0.8.0; `face_crit`)
    # Above the top level `face_apk` / `face_crit` continue them analytically
    # with the full width face_P[NF, s].
    face_z0 :: Vector{Float64}   # (n_slots) lowest point of the face = its sill (m)
    face_zr :: Vector{Float64}   # (n_slots) height of the tabulated range (m)
    face_A  :: Matrix{Float64}   # (NF × n_slots)
    face_P  :: Matrix{Float64}   # (NF × n_slots)
    face_K  :: Matrix{Float64}   # (NF × n_slots)
    # Wetted-area centroid at each level of elev_bins (m east / north of the
    # cell centre), of the ground connected to the cell's main basin: where
    # the water of a partly wet cell is (`wet_centroid`). Zeros: the centre.
    wet_cx  :: Vector{Float64}   # (n_bins)
    wet_cy  :: Vector{Float64}   # (n_bins)
    face_C  :: Matrix{Float64}   # (NF × n_slots), layout 4 (see above)
end
SGSTable(elev_bins, vol_curve, area_curve, cell_area, z_min, z_max, face_z0, face_zr, face_A, face_P, face_K) =
    SGSTable(elev_bins, vol_curve, area_curve, cell_area, z_min, z_max, face_z0, face_zr, face_A, face_P, face_K,
             zeros(length(elev_bins)), zeros(length(elev_bins)))
"""
Layout-3 form (no `face_C`): the critical integral is approximated from A and
P as if the wet part of the face had a uniform depth, C ≈ P·(A/P)^(3/2) — a
lower bound of the profile's own Σ h^(3/2) (fixtures; real tables are built
with `face_table_from_profile`).
"""
SGSTable(elev_bins, vol_curve, area_curve, cell_area, z_min, z_max, face_z0, face_zr, face_A, face_P, face_K,
         wet_cx, wet_cy) =
    SGSTable(elev_bins, vol_curve, area_curve, cell_area, z_min, z_max, face_z0, face_zr, face_A, face_P, face_K,
             wet_cx, wet_cy, _face_C_uniform(face_A, face_P))

"Critical integral of a uniform-depth section with flow area A and wetted width P (element-wise)."
_face_C_uniform(A::AbstractMatrix, P::AbstractMatrix) =
    [P[k] > 1e-12 ? Float64(P[k]) * (max(Float64(A[k]), 0.0) / P[k])^1.5 : 0.0 for k in CartesianIndices(A)]

"Levels per face table (see `SGSTable`)."
const SGS_FACE_LEVELS = 64

"Relative height (0–1) of face-table level k of nf."
@inline _face_level(k::Int, nf::Int) = ((k - 1) / (nf - 1))^2

"""
    face_table_from_profile(z, dx; nf = SGS_FACE_LEVELS) → (z0, zr, A, P, K, C)

The face table of a profile: elevations `z` (m) of points that each stand for
a width `dx` (m) of the face. `A`, `P`, `K`, `C` are the flow area, wetted
width, conveyance integral Σ dx·h^(5/3) and critical integral Σ dx·h^(3/2)
(0.8.0) at the `nf` levels
z0 + zr·((k−1)/(nf−1))². z0 is the lowest point; the range zr is the relief
of the profile plus max(2·relief, 0.5 m), so that the tabulated part covers
the levels at which the shape of the profile still matters to K (a channel in
a floodplain: K is not a function of the mean depth until the floodplain is
several channel depths under water) and `face_apk`'s one-term continuation
above it is accurate (≈ 1% for a channel a quarter of the face wide).
With `act` (same length as `z`) a point conveys only at levels ≥ act[p]
(≥ z[p]): its activation level (connectivity, Phase 3); z0 and the relief are
then taken over the activation levels.
"""
function face_table_from_profile(z::AbstractVector{<:Real}, dx::Real;
                                 act::AbstractVector{<:Real} = z, nf::Int = SGS_FACE_LEVELS)
    A = zeros(nf); P = zeros(nf); K = zeros(nf); C = zeros(nf)
    isempty(z) && return NaN, 0.0, A, P, K, C
    z0 = Float64(minimum(act)); relief = Float64(maximum(act)) - z0
    zr = relief + max(2.0 * relief, 0.5)
    for k in 1:nf
        η = k == nf ? z0 + zr : z0 + zr * _face_level(k, nf)
        a = 0.0; w = 0.0; c = 0.0; q = 0.0
        for i in eachindex(z)
            act[i] <= η || continue
            h = η - Float64(z[i])
            w += 1.0
            h > 0.0 && (a += h; c += h^(5 / 3); q += h * sqrt(h))
        end
        A[k] = a * dx; P[k] = w * dx; K[k] = c * dx; C[k] = q * dx
    end
    return z0, zr, A, P, K, C
end

"""
    face_apk(table, slot, wse) → (A, P, K)

Flow area (m²), wetted width (m) and conveyance integral Σ dx·h^(5/3) of the
face in adjacency `slot` at water level `wse`, from the face table: linear
between levels; above the top level A grows with the full width W = P_top,
P = W, and K = W·(d + c)^(5/3) with d the height above the top level and c
set by K at the top level (exact there and in the deep limit). Zeros at or
below the sill and for a slot without a table.
"""
@inline function face_apk(t, slot::Int, wse::Float64)::Tuple{Float64,Float64,Float64}   # SGSTable or a view
    @inbounds z0 = t.face_z0[slot]
    d = wse - z0
    d > 0.0 || return 0.0, 0.0, 0.0           # also a NaN sill (no table)
    fA = t.face_A; fP = t.face_P; fK = t.face_K
    nf = size(fA, 1)
    @inbounds zr = t.face_zr[slot]
    @inbounds Pt = fP[nf, slot]
    Pt > 0.0 || return 0.0, 0.0, 0.0
    if d >= zr
        dd = d - zr
        @inbounds At = fA[nf, slot]
        @inbounds Kt = fK[nf, slot]
        c = (Kt / Pt)^0.6
        return At + Pt * dd, Pt, Pt * (dd + c)^(5 / 3)
    end
    u  = sqrt(d / zr) * (nf - 1)
    lo = min(unsafe_trunc(Int, u), nf - 2)    # level index − 1 below wse
    r0 = (lo / (nf - 1))^2; r1 = ((lo + 1) / (nf - 1))^2
    f  = (d / zr - r0) / (r1 - r0)
    f  = f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f)
    @inbounds begin
        A = fA[lo + 1, slot] + f * (fA[lo + 2, slot] - fA[lo + 1, slot])
        P = fP[lo + 1, slot] + f * (fP[lo + 2, slot] - fP[lo + 1, slot])
        K = fK[lo + 1, slot] + f * (fK[lo + 2, slot] - fK[lo + 1, slot])
    end
    return A, P, K
end

"""
    face_crit(table, slot, wse) → C

Critical integral Σ dx·max(wse − a_p, 0)^(3/2) (m^(5/2)) of the face in
adjacency `slot` at level `wse` (a_p the points' activation levels), from the
face table as `face_apk` reads A, P, K: linear between levels; above the top
level C = W·(d + c)^(3/2) with W = P_top and c set by C at the top level.
`√g·(2/3)^(3/2)·C` is the discharge of the face at critical flow with the
energy level `wse`, each strip of the face at its own critical depth
(`_sgs_critical_q`, 0.8.0). Zero at or below the sill and for a slot without
a table.
"""
@inline function face_crit(t, slot::Int, wse::Float64)::Float64
    @inbounds z0 = t.face_z0[slot]
    d = wse - z0
    d > 0.0 || return 0.0
    fC = t.face_C; fP = t.face_P
    nf = size(fC, 1)
    @inbounds zr = t.face_zr[slot]
    @inbounds Pt = fP[nf, slot]
    Pt > 0.0 || return 0.0
    if d >= zr
        @inbounds Ct = fC[nf, slot]
        c = cbrt(Ct / Pt)^2
        x = d - zr + c
        return Pt * x * sqrt(x)
    end
    u  = sqrt(d / zr) * (nf - 1)
    lo = min(unsafe_trunc(Int, u), nf - 2)
    r0 = (lo / (nf - 1))^2; r1 = ((lo + 1) / (nf - 1))^2
    f  = (d / zr - r0) / (r1 - r0)
    f  = f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f)
    @inbounds return fC[lo + 1, slot] + f * (fC[lo + 2, slot] - fC[lo + 1, slot])
end

"""
    SGSTable(elev_bins, vol_curve, area_curve, cell_area, z_min, z_max,
             edge_area_curves, edge_perim_curves)

Test-fixture form (the pre-layout-3 fields): per-slot flow area and wetted
width tabulated on the cell's own levels (n_bins × n_slots). They are
resampled to face tables; with no profile to integrate, K is the
single-radius conveyance A^(5/3)/P^(2/3) and C the uniform-depth critical
integral P·(A/P)^(3/2). The sill of a slot is the highest
level at which its area is still zero.
"""
function SGSTable(elev_bins::Vector{Float64}, vol_curve::Vector{Float64}, area_curve::Vector{Float64},
                  cell_area::Real, z_min::Real, z_max::Real,
                  edge_area_curves::AbstractMatrix{<:Real}, edge_perim_curves::AbstractMatrix{<:Real})
    nb, ns = size(edge_area_curves); nf = SGS_FACE_LEVELS
    z0 = fill(NaN, ns); zr = zeros(ns)
    A = zeros(nf, ns); P = zeros(nf, ns); K = zeros(nf, ns)
    for s in 1:ns
        edge_area_curves[nb, s] > 0 || continue
        k0 = something(findlast(k -> edge_area_curves[k, s] <= 0, 1:nb), 1)
        z0[s] = elev_bins[k0]; zr[s] = elev_bins[nb] - z0[s]
        for k in 1:nf
            η = k == nf ? elev_bins[nb] : z0[s] + zr[s] * _face_level(k, nf)
            lo = clamp(searchsortedlast(elev_bins, η), 1, nb - 1)
            f = clamp((η - elev_bins[lo]) / (elev_bins[lo+1] - elev_bins[lo] + 1e-15), 0.0, 1.0)
            a = edge_area_curves[lo, s] + f * (edge_area_curves[lo+1, s] - edge_area_curves[lo, s])
            w = edge_perim_curves[lo, s] + f * (edge_perim_curves[lo+1, s] - edge_perim_curves[lo, s])
            A[k, s] = a; P[k, s] = w
            K[k, s] = w > 1e-9 ? a^(5 / 3) / w^(2 / 3) : 0.0
        end
    end
    return SGSTable(elev_bins, vol_curve, area_curve, Float64(cell_area), Float64(z_min), Float64(z_max),
                    z0, zr, A, P, K)
end

"""
    wse_from_volume(table, V) → Float64

Given a stored volume V (m³), return the water surface elevation (m) by
inverse-interpolating the SGS volume curve.

For V between 0 and vol_curve[end] the hypsometric curve is used.
For V > vol_curve[end] (cell overfull — water above the terrain ceiling) the WSE
is extrapolated linearly above z_max using the full cell_area as the wet area:
    WSE = z_max + (V - vol_curve[end]) / cell_area
This preserves a positive driving head between neighbouring overfull cells and
allows the SGS solver to propagate volume across a flat upstream basin even when
individual cells have filled their hypsometric range.
"""
@inline function wse_from_volume(t, V::Float64)::Float64      # SGSTable or a view of one (WP5)
    V <= 0.0 && return t.z_min
    V >= t.vol_curve[end] && return t.z_max + (V - t.vol_curve[end]) / t.cell_area
    # Binary search for bracketing interval
    lo, hi = 1, length(t.vol_curve)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        t.vol_curve[mid] <= V ? lo = mid : hi = mid
    end
    frac = (V - t.vol_curve[lo]) / (t.vol_curve[hi] - t.vol_curve[lo] + 1e-15)
    return t.elev_bins[lo] + frac * (t.elev_bins[hi] - t.elev_bins[lo])
end

"""
    wet_centroid(table, wse) → (cx, cy)

Centroid of the wetted area of the cell at level `wse` (m east / north of the
cell centre), interpolated from `wet_cx`/`wet_cy`: the lowest point of the
main basin for a dry cell, the polygon's centroid for a submerged one.
"""
@inline function wet_centroid(t, wse::Float64)::Tuple{Float64,Float64}   # SGSTable or a view
    eb = t.elev_bins; n = length(eb)
    @inbounds begin
        wse <= eb[1] && return t.wet_cx[1], t.wet_cy[1]
        wse >= eb[n] && return t.wet_cx[n], t.wet_cy[n]
        lo, hi = 1, n
        while hi - lo > 1
            mid = (lo + hi) >>> 1
            eb[mid] <= wse ? lo = mid : hi = mid
        end
        f = (wse - eb[lo]) / (eb[hi] - eb[lo] + 1e-15)
        return t.wet_cx[lo] + f * (t.wet_cx[hi] - t.wet_cx[lo]), t.wet_cy[lo] + f * (t.wet_cy[hi] - t.wet_cy[lo])
    end
end

"""
    wetted_area_from_wse(table, wse) → Float64

Return the wetted plan area (m²) at a given WSE by interpolating area_curve.
"""
@inline function wetted_area_from_wse(t, wse::Float64)::Float64   # SGSTable or a view (WP5)
    wse <= t.z_min && return 0.0
    wse >= t.z_max && return t.cell_area
    lo, hi = 1, length(t.elev_bins)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        t.elev_bins[mid] <= wse ? lo = mid : hi = mid
    end
    frac = (wse - t.elev_bins[lo]) / (t.elev_bins[hi] - t.elev_bins[lo] + 1e-15)
    return t.area_curve[lo] + frac * (t.area_curve[hi] - t.area_curve[lo])
end

# ---------------------------------------------------------------------------
# Great-circle arc geometry helpers (used for edge sill and cell area)
# ---------------------------------------------------------------------------

const _EARTH_R = 6_371_000.0   # mean Earth radius in metres

"""Haversine distance in metres between two lon/lat points."""
@inline function _haversine_m(lon1::Float64, lat1::Float64,
                               lon2::Float64, lat2::Float64)::Float64
    φ1, φ2 = deg2rad(lat1), deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)
    a  = sin(Δφ/2)^2 + cos(φ1)*cos(φ2)*sin(Δλ/2)^2
    return 2 * _EARTH_R * asin(min(1.0, sqrt(a)))
end

"""
    _arc_sample_lons_lats(lon1, lat1, lon2, lat2, spacing_m) → (lons, lats)

Generate sample points along the great-circle arc from (lon1,lat1) to
(lon2,lat2) at approximately `spacing_m` metre intervals.  Used to sample
the DEM along a shared cell edge to find the sill elevation.
"""
function _arc_sample_lons_lats(lon1::Float64, lat1::Float64,
                                lon2::Float64, lat2::Float64,
                                spacing_m::Float64)
    d   = _haversine_m(lon1, lat1, lon2, lat2)
    d < 1.0 && return [lon1], [lat1]
    n   = max(2, round(Int, d / spacing_m) + 1)
    ts  = range(0.0, 1.0, length=n)
    # Interpolate on a sphere via slerp approximation (flat for short arcs)
    lons = lon1 .+ ts .* (lon2 - lon1)
    lats = lat1 .+ ts .* (lat2 - lat1)
    return collect(lons), collect(lats)
end

"""
    _polygon_area_m2(boundary) → Float64

Approximate planar area of a lon/lat polygon in m², using the shoelace
formula on an equirectangular projection centred on the polygon centroid.
Accurate to ~0.1% for A5 cells at resolution 14 (~km scale).
"""
function _polygon_area_m2(boundary::Vector{Vector{Float64}})::Float64
    n = length(boundary)
    n < 3 && return 0.0
    lat0 = mean(v[2] for v in boundary)
    cos_lat = cosd(lat0)
    # Shoelace on projected coordinates (metres)
    area = 0.0
    j = n
    for i in 1:n
        xi = deg2rad(boundary[i][1]) * _EARTH_R * cos_lat
        yi = deg2rad(boundary[i][2]) * _EARTH_R
        xj = deg2rad(boundary[j][1]) * _EARTH_R * cos_lat
        yj = deg2rad(boundary[j][2]) * _EARTH_R
        area += (xj + xi) * (yj - yi)
        j = i
    end
    return abs(area) / 2.0
end

"""
    _shared_edge(bnd_i, bnd_j) → (lon1, lat1, lon2, lat2) or nothing

Find the shared edge between two A5 pentagon boundaries.  Pentagon vertices
are stored as [lon, lat] pairs.  Returns the two shared vertices if found,
or `nothing` if no shared edge exists.

We match vertices with a tolerance of 1e-9 degrees (~0.1 mm) to handle
floating-point rounding in pya5's boundary computation.
"""
function _shared_edge(bnd_i::Vector{Vector{Float64}},
                      bnd_j::Vector{Vector{Float64}})
    # Thread-safe implementation: linear scan with no heap-allocated Set or Dict.
    # The original Set{Tuple} construction triggered PyCall GC finalisation on
    # non-main threads (EXCEPTION_ACCESS_VIOLATION in PyObject_ClearWeakRefs),
    # because Julia's GC can run on any thread during Dict rehash allocation.
    # A pentagon boundary has at most 6 vertices (5 + closing repeat), so the
    # O(n²) scan over ≤36 pairs is faster than a Set for this input size.
    tol = 1e-9
    shared = Vector{Tuple{Float64,Float64}}()
    for vi in bnd_i
        ri = (round(vi[1], digits=9), round(vi[2], digits=9))
        for vj in bnd_j
            rj = (round(vj[1], digits=9), round(vj[2], digits=9))
            if ri == rj
                push!(shared, (vi[1], vi[2]))
                length(shared) == 2 && @goto done
                break   # each vertex in bnd_i matches at most one in bnd_j
            end
        end
    end
    @label done
    length(shared) < 2 && return nothing
    return shared[1][1], shared[1][2], shared[2][1], shared[2][2]
end

"""
    _edge_geometry(bnd_i, bnd_j, lon_i, lat_i, lon_j, lat_j)
        → (cos_theta::Float64, corr_x::Float64, corr_y::Float64)

Combined non-orthogonality geometry for one A5 cell-pair edge: the
unsigned magnitude `cos_theta` (used by the legacy/uncorrected flux
kernels), and the WLSQ gradient-correction direction vector `(corr_x,
corr_y)` (used by the corrected kernels).

⚠️  REVISION NOTE (2026-06-24): this function previously returned a
    `skew_vec` (face-midpoint-to-line-intersection offset) for the second
    pair of values, intended for use as `dWSE_n = (wse_i-wse_j) +
    ∇WSE_f · skew_vec`. That formula was found to be geometrically
    incorrect during T-NOC5b investigation (it used a *positional* offset
    vector where the rigorous correction requires a *directional* one —
    they coincide only in the trivial orthogonal limit). This function now
    returns the correct vector, re-derived and independently verified
    against five distinct synthetic geometries (see implementation notes
    below and `FloodA5_NonOrthogonal_Correction_Plan.md` §10.6 for the full
    derivation record). The field names `skew_x`/`skew_y` in `EdgeList`
    were intentionally NOT renamed to `corr_x`/`corr_y` to keep this diff
    minimal — they now hold the corrected vector's components, not a
    skewness offset. See `EdgeList`'s field comment for the updated
    meaning.

Returns:

  cos_theta — |d̂ · n̂|, the unsigned magnitude of alignment between the
              centre-to-centre unit vector d̂ (cell i → cell j) and the
              edge face normal. 1.0 = perfectly orthogonal. Used by the
              legacy kernels' `L_eff = L × cos θ` magnitude-only
              projection. Unsigned because the legacy kernels have no
              concept of a consistently-oriented normal — they only ever
              use this as a scalar distance-scaling factor.

  corr_x,
  corr_y    — components (m⁻¹, dimensionless direction vector — NOT a
              positional offset) of `V̂ = n̂ − c·d̂`, where `n̂` is the face
              normal **consistently oriented so that d̂·n̂ ≥ 0** (i.e.
              re-oriented to point from cell i toward cell j; this
              orientation step is essential — the raw candidate normal
              from rotating the edge vector has an orientation that
              depends only on shared-vertex ordering, unrelated to which
              cell is i and which is j) and `c = d̂·n̂` (always ≥ 0 after
              orientation). `(corr_x, corr_y)` is used by the WLSQ-corrected
              kernels to build the corrected driving head:
                  dWSE_n = c·(wse_i − wse_j) − L·∇WSE_f·(corr_x, corr_y)
              `(corr_x, corr_y) = (0, 0)` exactly at perfect orthogonality
              (c = ±1), recovering the legacy term exactly with no
              correction needed — verified analytically and numerically.
              See FloodA5_NonOrthogonal_Correction_Plan.md §5.1/§10.6 for
              the full derivation, including why the orientation step is
              necessary and why the value returned here is `c`-free (the
              caller multiplies by the now-signed, always-nonnegative `c`
              separately — NOT by `cos_theta` from this same return tuple,
              since `cos_theta` is for the *legacy* kernel and is
              positionally identical to `c` in magnitude but the caller
              must use the oriented version consistently within the
              corrected-kernel code path).

Returns `(1.0, 0.0, 0.0)` if no shared edge is found (non-adjacent cells)
— an orthogonal, uncorrected fallback rather than NaN. A NaN cos_theta
would propagate through the NaN guard in the step functions and suppress
all flux on that edge, which is far worse than a small error in the slope
projection.

Returns `(NaN, 0.0, 0.0)` for degenerate edge/distance geometry (near-zero
edge length or near-zero centre separation) — callers already guard NaN
cos_theta by falling back to 1.0, and a zero correction vector is the safe
choice when the geometry itself is degenerate.

Implementation notes
---------------------
All geometry is performed in a local equirectangular projection centred on
the edge midpoint. At A5 resolution 14 (~1.4 km cell spacing) the planar
approximation introduces < 0.05% error. The same projection is used at finer
resolutions (including multi-resolution meshes in Phase 3), where the error
is even smaller.

The function accepts raw boundary arrays and cell-centre coordinates so it
can be called for cell pairs at *different* A5 resolution levels (coarse/fine
boundaries in a multi-resolution mesh) without any change to its interface.
"""
function _edge_geometry(bnd_i  :: Vector{Vector{Float64}},
                         bnd_j  :: Vector{Vector{Float64}},
                         lon_i  :: Float64, lat_i  :: Float64,
                         lon_j  :: Float64, lat_j  :: Float64
                         )::Tuple{Float64,Float64,Float64}

    # ── 1. Find the shared edge vertices ─────────────────────────────
    edge = _shared_edge(bnd_i, bnd_j)
    edge === nothing && return (1.0, 0.0, 0.0)

    elon1, elat1, elon2, elat2 = edge

    # ── 2. Local equirectangular projection ────────────────────────
    lat0     = (elat1 + elat2) * 0.5
    cos_lat0 = cosd(lat0)
    R        = _EARTH_R

    to_xy(lon, lat) = (deg2rad(lon) * R * cos_lat0,
                       deg2rad(lat) * R)

    # ── 3. Edge vector and face normal candidate ──────────────────
    x1, y1 = to_xy(elon1, elat1)
    x2, y2 = to_xy(elon2, elat2)
    ex, ey  = x2 - x1, y2 - y1
    e_len   = sqrt(ex^2 + ey^2)
    # Degenerate-length guards: 1 mm (WP3, 2026-10; was 1 m, which is a
    # normal face length on Voronoi meshes and on A5 at resolution >= 22).
    e_len < 1e-3 && return (NaN, 0.0, 0.0)   # degenerate edge guard

    # Rotate the edge vector 90° to get a normal candidate. This has NO
    # guaranteed orientation relative to i/j — it depends only on the
    # (elon1,elat1)/(elon2,elat2) vertex ordering from _shared_edge, which
    # is unrelated to which cell is i and which is j. The orientation step
    # below (step 5) fixes this; do not skip it.
    nx, ny  = -ey / e_len, ex / e_len

    # ── 4. Centre-to-centre unit vector (i → j) ────────────────────
    ci_x, ci_y = to_xy(lon_i, lat_i)
    cj_x, cj_y = to_xy(lon_j, lat_j)
    dx, dy      = cj_x - ci_x, cj_y - ci_y
    d_len       = sqrt(dx^2 + dy^2)
    d_len < 1e-3 && return (NaN, 0.0, 0.0)   # degenerate distance guard

    dx_hat, dy_hat = dx / d_len, dy / d_len

    # ── 5. Orient n̂ so that d̂·n̂ ≥ 0 (points from i toward j) ─────────
    # This is the critical step that was missing before the 2026-06-24
    # revision. Without it, the signed correction vector below would have
    # an arbitrary, geometry-vertex-ordering-dependent sign, unrelated to
    # the physically meaningful i→j direction.
    c_raw = dx_hat * nx + dy_hat * ny
    if c_raw < 0.0
        nx, ny = -nx, -ny
        c_raw = -c_raw
    end
    c = clamp(c_raw, 0.0, 1.0)   # = cos_theta, but kept separate for clarity
                                  # at the call site; guard float rounding.

    cos_theta = c   # unsigned magnitude — identical to c now that n̂ is
                     # oriented (c ≥ 0 by construction), kept as a separate
                     # named return value for the legacy kernels' clarity.

    # ── 6. WLSQ correction direction: V̂ = n̂ − c·d̂ ─────────────────────
    # The tangential remainder of the (now consistently i→j-oriented) face
    # normal after removing its component along the centre-to-centre
    # direction. Verified to vanish exactly at orthogonality (c = 1) and to
    # exactly recover the true face-normal gradient for an arbitrary linear
    # field when used in the formula documented above. See the docstring's
    # revision note for the verification record.
    corr_x = nx - c * dx_hat
    corr_y = ny - c * dy_hat

    return (cos_theta, corr_x, corr_y)
end

# ---------------------------------------------------------------------------
# Boundary sides of a cell and their table slots (0.7.0)
# ---------------------------------------------------------------------------

"The sides of a cell polygon as vertex pairs, in polygon order (a closing vertex is dropped)."
function _polygon_sides(boundary)
    nv = length(boundary)
    n  = (nv >= 2 && boundary[nv][1] ≈ boundary[1][1] && boundary[nv][2] ≈ boundary[1][2]) ? nv - 1 : nv
    return [(boundary[k], boundary[mod1(k + 1, n)]) for k in 1:n]
end

"Whether the side va–vb is one of the `shared` edges (lon1, lat1, lon2, lat2), in either direction."
function _side_is_shared(va, vb, shared)
    near(a, b) = abs(a - b) < 1e-8
    for sp in shared
        (near(sp[1], va[1]) && near(sp[2], va[2]) && near(sp[3], vb[1]) && near(sp[4], vb[2])) && return true
        (near(sp[1], vb[1]) && near(sp[2], vb[2]) && near(sp[3], va[1]) && near(sp[4], va[2])) && return true
    end
    return false
end

"""
    boundary_side_slots(n_slots, in_mesh) → Vector{Int}

The table slots of a cell that hold no face: slots of neighbour-list entries
that are not in the mesh (`in_mesh[s] == false`) and slots beyond the list.
The cell's boundary sides (polygon sides shared with no mesh neighbour) take
them in polygon-side order: the b-th boundary side has slot
`boundary_side_slots(…)[b]`, in `build_sgs_tables!` and in the ghost-edge
builder alike. A boundary side beyond the free slots has no table.
"""
boundary_side_slots(n_slots::Int, in_mesh::AbstractVector{Bool}) =
    [s for s in 1:n_slots if s > length(in_mesh) || !in_mesh[s]]

# ---------------------------------------------------------------------------
# SGS pre-processing — build_sgs_tables!
# ---------------------------------------------------------------------------

"""
    _cell_connectivity(z, idx, sidemask, need) → (seed, conn)

Connectivity of a cell's working grid (0.7.0). `z[p]` are the elevations of
the grid points inside the cell, `idx[i, j]` their index on the regular grid
(0: outside the cell or no data) and `sidemask[p]` the polygon sides the
point lies next to (bit k: side k).

The cell's **main basin** is the first region of ground below a rising water
level (4-connected, so a wall one point thick and diagonal to the grid still
holds) that reaches `need` different sides of the cell — the lowest route
*through* the cell. A pit that touches no side, or a hollow cut by one side,
does not qualify until it has spilled into something that does. `seed` is
the lowest point of that region.

`conn[p]` is the lowest water level at which point p is connected to the seed
through ground below that level (≥ z[p]; = z[p] for the main basin itself).
"""
function _cell_connectivity(z::Vector{Float64}, idx::Matrix{Int32}, gi::Vector{Int32},
                            gj::Vector{Int32}, sidemask::Vector{UInt32}, need::Int)
    np = length(z)
    np == 0 && return 0, Float64[]
    order  = sortperm(z)
    parent = collect(1:np)
    find(i) = (@inbounds while parent[i] != i; parent[i] = parent[parent[i]]; i = parent[i]; end; i)
    nxg, nyg = size(idx)
    @inline function nbr(p, k)
        i = Int(gi[p]); j = Int(gj[p])
        k == 1 ? (i > 1   ? Int(idx[i-1, j]) : 0) :
        k == 2 ? (i < nxg ? Int(idx[i+1, j]) : 0) :
        k == 3 ? (j > 1   ? Int(idx[i, j-1]) : 0) :
                 (j < nyg ? Int(idx[i, j+1]) : 0)
    end

    # pass 1: the first region to touch `need` sides
    mask  = copy(sidemask); low = collect(1:np); added = falses(np)
    seed  = 0
    for p in order
        added[p] = true
        for k in 1:4
            q = nbr(p, k); (q == 0 || !added[q]) && continue
            a = find(p); b = find(q); a == b && continue
            parent[b] = a
            mask[a] |= mask[b]
            z[low[b]] < z[low[a]] && (low[a] = low[b])
        end
        r = find(p)
        if count_ones(mask[r]) >= need
            seed = low[r]; break
        end
    end
    seed == 0 && (seed = order[1])

    # pass 2: the level at which each point joins the seed's region
    conn = copy(z)
    parent .= 1:np; fill!(added, false)
    next = zeros(Int, np); tail = collect(1:np)       # member lists per root
    for p in order
        added[p] = true
        for k in 1:4
            q = nbr(p, k); (q == 0 || !added[q]) && continue
            a = find(p); b = find(q); a == b && continue
            if added[seed]
                s = find(seed)
                if s == a || s == b
                    o = s == a ? b : a                # the region joining the main basin now
                    m = o
                    while m != 0
                        conn[m] = max(conn[m], z[p]); m = next[m]
                    end
                end
            end
            next[tail[a]] = b; tail[a] = tail[b]      # append b's members to a's
            parent[b] = a
        end
    end
    return seed, conn
end

"Distance from (x, y) to the segment (x1, y1)–(x2, y2)."
@inline function _seg_dist(x, y, x1, y1, x2, y2)
    dx = x2 - x1; dy = y2 - y1; l2 = dx * dx + dy * dy
    t = l2 > 0 ? clamp(((x - x1) * dx + (y - y1) * dy) / l2, 0.0, 1.0) : 0.0
    return hypot(x - (x1 + t * dx), y - (y1 + t * dy))
end

"""
    _pixel_blocks(mesh, raster) → NamedTuple

For `build_sgs_tables!(…; sampling = :pixels)`: the block of the raster's
pixels that each cell covers. Every cell must be a quadrilateral whose
corners, transformed to the raster's CRS, fall on pixel corners (within
1e-3 pixel) — a raster-aligned Cartesian grid whose cells are whole blocks of
the DEM's pixels; anything else is an error. Returns `block[ci] = (c0, c1,
r0, r1)` (0-based pixel columns c0 … c1−1, rows r0 … r1−1), `corner[ci]` the
pixel corners (c, r) of the cell's polygon vertices in side order, the
pixel centres' lon/lat (`cen_lon`, `cen_lat`, cell ci's block row-major from
`cen_off[ci] + 1`), and `max_pixels`.
"""
function _pixel_blocks(mesh::A5Mesh, raster)
    n = length(mesh.cells)
    sides = [_polygon_sides(c.boundary) for c in mesh.cells]
    nv = [length(s) for s in sides]
    all(==(4), nv) || error("SGS pixel sampling: every cell must be a quadrilateral " *
                            "(a raster-aligned Cartesian grid); found cells with $(sort(unique(nv))) sides")
    lons = Float64[]; lats = Float64[]
    for s in sides, (va, _) in s
        push!(lons, Float64(va[1])); push!(lats, Float64(va[2]))
    end
    xs, ys = _reproject_to_dem(lons, lats, raster.proj)
    block = Vector{NTuple{4,Int}}(undef, n)
    corner = Vector{Vector{NTuple{2,Int}}}(undef, n)
    worst = 0.0
    for ci in 1:n
        cs = NTuple{2,Int}[]
        for k in 1:4
            cf, rf = _affine_to_pixel(xs[4(ci-1)+k], ys[4(ci-1)+k], raster.gt)
            c = round(Int, cf); r = round(Int, rf)
            worst = max(worst, abs(cf - c), abs(rf - r))
            push!(cs, (c, r))
        end
        c0, c1 = extrema(first.(cs)); r0, r1 = extrema(last.(cs))
        (c1 > c0 && r1 > r0 && all(p -> p[1] in (c0, c1) && p[2] in (r0, r1), cs) &&
         length(unique(cs)) == 4) ||
            error("SGS pixel sampling: cell $(mesh.cells[ci].id) is not a block of the DEM's pixels " *
                  "(its corners fall on pixel corners $cs)")
        block[ci] = (c0, c1, r0, r1); corner[ci] = cs
    end
    worst > 1e-3 && error("SGS pixel sampling: the cells' corners are up to $(round(worst; sigdigits = 3)) " *
                          "pixel from the DEM's pixel corners — the mesh is not aligned with this DEM " *
                          "(make it with --meshgen-cartesian from a block average of it, cases/scripts/block_dem.jl)")
    # pixel centres, raster CRS → lon/lat
    cen_off = zeros(Int, n); tot = 0
    for ci in 1:n
        c0, c1, r0, r1 = block[ci]; cen_off[ci] = tot; tot += (c1 - c0) * (r1 - r0)
    end
    cx = Vector{Float64}(undef, tot); cy = similar(cx)
    gt = raster.gt
    for ci in 1:n
        c0, c1, r0, r1 = block[ci]; q = cen_off[ci]
        for r in r0:(r1 - 1), c in c0:(c1 - 1)
            q += 1
            cx[q] = gt[1] + (c + 0.5) * gt[2] + (r + 0.5) * gt[3]
            cy[q] = gt[4] + (c + 0.5) * gt[5] + (r + 0.5) * gt[6]
        end
    end
    if !isempty(raster.proj)
        zz = zeros(tot)
        ArchGDAL.createcoordtrans(_crs_gis(raster.proj), _crs_gis(4326)) do xform
            ArchGDAL.transform!(cx, cy, zz, xform)
        end
    end
    max_pixels = maximum((b[2] - b[1]) * (b[4] - b[3]) for b in block)
    @info "  Pixel sampling: $n cells, blocks of up to $max_pixels pixels " *
          "(corners within $(round(worst; sigdigits = 2)) pixel of the DEM's)"
    return (; block, corner, cen_lon = cx, cen_lat = cy, cen_off, max_pixels)
end

"""
    _sgs_pixel_cells!(mesh, raster, pix, connectivity, …)

Step 1 of `build_sgs_tables!` with `sampling = :pixels` (see there): per
cell, the hypsometric curve, connectivity, wetted-area centroids and side
profiles from the pixels of its block (`_pixel_blocks`).
"""
function _sgs_pixel_cells!(mesh, raster, pix, connectivity, elev_bins_mat, vol_curve_mat, area_curve_mat,
                           wet_cx_mat, wet_cy_mat, cell_areas, z_mins, z_maxs, n_grid,
                           side_z, side_act, side_W, side_zin, side_actin)
    n = length(mesh.cells); n_bins = size(elev_bins_mat, 1)
    d2m = _EARTH_R * π / 180.0
    nd = raster.nodata_val
    function pixel(c, r)
        (0 <= c < raster.nx && 0 <= r < raster.ny) || return NaN
        v = raster.band[c + 1, r + 1]
        return (isfinite(v) && !(isfinite(nd) && v == nd)) ? v : NaN
    end
    Threads.@threads for ci in 1:n
        cell = mesh.cells[ci]
        bnd  = cell.boundary
        lonc = Float64(cell.center_lon); latc = Float64(cell.center_lat)
        mx = d2m * cosd(latc); my = d2m
        A_cell = _polygon_area_m2(bnd); cell_areas[ci] = A_cell
        sides = _polygon_sides(bnd)
        c0, c1, r0, r1 = pix.block[ci]
        kc = c1 - c0; kr = r1 - r0; ntot = kc * kr
        dA = A_cell / ntot; hg = sqrt(dA)

        # the block's pixels (the working grid)
        idx = zeros(Int32, kc, kr)
        zg = Float64[]; gi = Int32[]; gj = Int32[]; gx = Float64[]; gy = Float64[]
        q = pix.cen_off[ci]
        for r in r0:(r1 - 1), c in c0:(c1 - 1)
            q += 1
            v = pixel(c, r); isfinite(v) || continue
            push!(zg, v); push!(gi, Int32(c - c0 + 1)); push!(gj, Int32(r - r0 + 1))
            push!(gx, (pix.cen_lon[q] - lonc) * mx); push!(gy, (pix.cen_lat[q] - latc) * my)
            idx[c - c0 + 1, r - r0 + 1] = Int32(length(zg))
        end
        np = length(zg); n_grid[ci] = np

        # sides: per pixel along the side (from vertex a to b, reversed into the
        # canonical direction), the pixel inside and the one across the side
        ns = length(sides)
        zin = Vector{Vector{Float64}}(undef, ns); zmx = similar(zin)
        pin = Vector{Vector{Int}}(undef, ns); sW = zeros(ns)
        for (k, (va, vb)) in enumerate(sides)
            (ca, ra) = pix.corner[ci][k]; (cb, rb) = pix.corner[ci][mod1(k + 1, 4)]
            pairs = NTuple{4,Int}[]                     # (c_in, r_in, c_out, r_out)
            if ra == rb                                  # along a row boundary
                rin = ra == r0 ? ra : ra - 1; rout = ra == r0 ? ra - 1 : ra
                for c in (ca < cb ? (ca:(cb - 1)) : ((ca - 1):-1:cb))
                    push!(pairs, (c, rin, c, rout))
                end
            else                                         # along a column boundary
                cin = ca == c0 ? ca : ca - 1; cout = ca == c0 ? ca - 1 : ca
                for r in (ra < rb ? (ra:(rb - 1)) : ((ra - 1):-1:rb))
                    push!(pairs, (cin, r, cout, r))
                end
            end
            a1, a2, b1, b2 = Float64(va[1]), Float64(va[2]), Float64(vb[1]), Float64(vb[2])
            (a1, a2) > (b1, b2) && reverse!(pairs)       # the canonical direction (as bilinear sampling)
            sW[k] = _haversine_m(a1, a2, b1, b2)
            zi = [pixel(p[1], p[2]) for p in pairs]
            zo = [pixel(p[3], p[4]) for p in pairs]
            zin[k] = zi
            zmx[k] = [isfinite(zo[t]) && isfinite(zi[t]) ? max(zi[t], zo[t]) : zi[t] for t in eachindex(zi)]
            pin[k] = [Int(idx[p[1] - c0 + 1, p[2] - r0 + 1]) for p in pairs]
        end
        side_W[ci] = sW; side_z[ci] = zmx; side_zin[ci] = zin
        if np == 0
            elev_bins_mat[:, ci] .= NaN; vol_curve_mat[:, ci] .= 0.0; area_curve_mat[:, ci] .= 0.0
            side_act[ci] = [copy(v) for v in zmx]; side_actin[ci] = [copy(v) for v in zin]
            continue
        end

        z_min = minimum(zg); z_max = maximum(zg)
        for v in (zin, zmx), s in v, zz in s
            isfinite(zz) || continue
            zz < z_min && (z_min = zz); zz > z_max && (z_max = zz)
        end
        z_mins[ci] = z_min; z_maxs[ci] = z_max

        # connectivity on the pixel grid
        seed, conn = if connectivity
            smask = zeros(UInt32, np); touched = UInt32(0)
            for p in 1:np, (k, (va, vb)) in enumerate(sides)
                d = _seg_dist(gx[p], gy[p], (Float64(va[1]) - lonc) * mx, (Float64(va[2]) - latc) * my,
                              (Float64(vb[1]) - lonc) * mx, (Float64(vb[2]) - latc) * my)
                if d <= 0.75 * hg
                    smask[p] |= UInt32(1) << (k - 1); touched |= UInt32(1) << (k - 1)
                end
            end
            _cell_connectivity(zg, idx, gi, gj, smask, min(2, count_ones(touched)))
        else
            argmin(zg), zg
        end
        # activation: a side point conveys from the level at which the pixel
        # inside it is wet and connected to the main basin
        actf(zs, k) = [isfinite(zs[t]) ? (pin[k][t] > 0 ? max(zs[t], conn[pin[k][t]]) : zs[t]) : NaN
                       for t in eachindex(zs)]
        side_act[ci] = [actf(zmx[k], k) for k in 1:ns]
        side_actin[ci] = [actf(zin[k], k) for k in 1:ns]

        # bins: every pixel elevation (the volume is linear between them),
        # filled out to n_bins with a uniform range and then the largest gaps' midpoints
        zu = sort(unique(zg))
        b = sort(unique(vcat(collect(range(z_min, z_max; length = max(2, n_bins - length(zu)))), zu)))
        if length(b) == 1                               # a flat block: one level, repeated
            b = fill(b[1], n_bins)
        end
        while length(b) < n_bins
            g = argmax(diff(b)); insert!(b, g + 1, 0.5 * (b[g] + b[g + 1]))
        end
        # exact level-pool volume, wetted area and connected centroid at each bin
        zs = sort(zg); oc = sortperm(conn)
        cum_area = 0.0; cum_vol = 0.0; ptr = 1; n_in = 0; prev = b[1]
        cptr = 1; sx = 0.0; sy = 0.0; nc = 0
        for k in 1:n_bins
            w = b[k]
            cum_vol += n_in * dA * (w - prev)
            while ptr <= np && zs[ptr] <= w
                cum_area += dA; cum_vol += dA * (w - zs[ptr]); n_in += 1; ptr += 1
            end
            while cptr <= np && conn[oc[cptr]] <= w
                p = oc[cptr]; sx += gx[p]; sy += gy[p]; nc += 1; cptr += 1
            end
            elev_bins_mat[k, ci] = w; area_curve_mat[k, ci] = cum_area; vol_curve_mat[k, ci] = cum_vol
            wet_cx_mat[k, ci] = nc > 0 ? sx / nc : gx[seed]
            wet_cy_mat[k, ci] = nc > 0 ? sy / nc : gy[seed]
            prev = w
        end
    end
    return nothing
end

"""
    build_sgs_tables!(mesh, dem_source;
                      n_bins=100, n_samples=1600, edge_spacing_m=nothing,
                      connectivity=true)

Pre-compute the sub-grid (SGS) tables of every cell and face and store them
in `mesh.array_vars` / `mesh.static_vars` (layout 4, 0.8.0: the DEM sampled
at pixel centres, and the faces' critical integrals `sgs_face_C`).

**Working grid.** Each cell is covered by a regular grid in a local metric
frame, spaced min(DEM pixel, √(bounding box / `n_samples`)), and the DEM is
sampled (bilinear) at the grid points inside the polygon. The grid is used
for everything below, so the tables of a cell and of its faces describe the
same surface. (Until 0.7.0: 512 Halton candidates per cell whatever its size
or the DEM's resolution, about two thirds of them inside a pentagon.)

**Pixel sampling** (`sampling = :pixels`, 0.9.1; `--sgs-sampling pixels`).
For a Cartesian mesh whose cells are whole blocks of the DEM's pixels (a
raster-aligned grid from a block average of the DEM,
`cases/scripts/block_dem.jl`), the tables are built from the pixels
themselves rather than a bilinear surface: the cell is the k × k pixels of
its block (the working grid), its hypsometric curve is exactly
Σ_p (A/k²)·max(η − z_p, 0) — the bins include every pixel elevation, between
which the volume is linear, so the table holds the volume the DEM's pixels
store at every level — and a face's profile is, per pixel pair across it,
max(z_in, z_out): the bed a cell-mean model on the pixels puts on that face
(a boundary side takes the pixel inside). The SGS model then stores, at any
level, the water a model on the DEM's own grid stores at that level, so the
two can be compared exactly. `n_bins` is raised to the block's pixel count
+ 2 if smaller. Any other mesh is an error.

**Per cell**
1. Hypsometric curve: stored volume and wetted area at `n_bins` levels from
   the cell's lowest to its highest ground, from four samples per grid cell
   (so that the grid does not alias with the DEM's lattice); the extremes are
   taken over them, the grid *and* the profiles of the cell's sides (so that
   no face point lies outside the cell's range).
2. Connectivity (`connectivity = true`): the level at which each grid point
   connects to the cell's main basin, `_cell_connectivity`.
3. Wetted-area centroid at each level (of the ground connected to the main
   basin), relative to the cell centre, for the flow-path length and
   direction of a partly wet cell (`wet_centroid`).

**Per face** (and per boundary side): the DEM profile along the side at
`edge_spacing_m` (default half a pixel) and, per point, its activation
level: the lowest level at which it is wet *and* connected to the main basins
of both cells. One face table (`face_table_from_profile`) per face, held by
both cells; the sill of a face is its lowest activation level. A barrier
inside a cell therefore closes the faces on the side cut off from the cell's
main basin up to the level at which the water gets over or round it — a gap
through the barrier is found by the connectivity search and keeps its own
sill. `sgs_face_raise` records by how much each face's sill was raised.

Output arrays in `mesh.array_vars`:
  "sgs_elev_bins", "sgs_vol_curve", "sgs_area_curve" — (n_bins × n_cells)
  "sgs_wet_cx", "sgs_wet_cy" — (n_bins × n_cells) wetted-area centroid (m east /
                      north of the cell centre)
  "sgs_edge_sills"  — (n_slots × n_cells) sill per slot = lowest level of the face table
  "sgs_face_zr"     — (n_slots × n_cells) height of the tabulated range (m)
  "sgs_face_A", "sgs_face_P", "sgs_face_K", "sgs_face_C" — (NF·n_slots × n_cells) face tables
                      (see `SGSTable`), column-major (level fastest)
  "sgs_face_raise"  — (n_slots × n_cells) sill − lowest point of the profile (m)

Slot k of a cell is the k-th entry of its neighbour list `mesh.adjacency[id]`
(the order `initialise_flow_model` uses for `adj_matrix`); boundary sides
take the slots without an in-mesh neighbour (`boundary_side_slots`).
`n_slots` is the longest neighbour list, so any valence is supported. If the
mesh has no adjacency, it is computed with `grid_disk_neighbours` (A5 only;
for any other mesh that is an error).

Output scalars in `mesh.static_vars`: "sgs_cell_area", "sgs_z_min",
"sgs_z_max", "sgs_n_bins", "sgs_layout" (4), "sgs_connectivity" (1 or 0),
"sgs_n_grid" (grid points used per cell).
"""
function build_sgs_tables!(mesh::A5Mesh, dem_source::DEMSource;
                            n_bins::Int=100,
                            n_samples::Int=1600,
                            halton_seed::Int=0,        # unused since 0.7.0 (kept for callers)
                            edge_spacing_m::Union{Float64,Nothing}=nothing,
                            connectivity::Bool=true,
                            sampling::Symbol=:bilinear)
    n = length(mesh.cells)
    n == 0 && error("Cannot build SGS tables on an empty mesh.")
    sampling in (:bilinear, :pixels) ||
        throw(ArgumentError("build_sgs_tables!: sampling must be :bilinear or :pixels (got :$sampling)"))
    # pya5's grid_disk (the fallback below) understands A5 cell IDs only; every
    # other generator writes a neighbours column (the mesh contract,
    # docs/MESHES.md §1).
    isempty(mesh.adjacency) && !startswith(mesh.id_scheme, "a5-") &&
        error("build_sgs_tables!: the mesh (id scheme \"$(mesh.id_scheme)\") has no " *
              "neighbour lists, and only A5 meshes can compute them (pya5 grid_disk). " *
              "Load it from a GeoParquet file with a `neighbours` column, or regenerate it.")
    dem_path = _resolve_dem_path(dem_source, mesh)

    @info "SGS pre-processing: $(basename(dem_path)), $n cells, $n_bins bins, " *
          "working grid ≈ $n_samples points/cell, connectivity $(connectivity ? "on" : "off")..."
    t0 = time()

    raster = _load_dem_raster(dem_path)
    lat_ref = sum(c.center_lat for c in mesh.cells) / n

    # DEM pixel size in metres. A geographic DEM (pixel < 1: degrees) is
    # converted with the mesh's latitude (0.7.0: the east–west size was taken
    # without cos φ, so faces were sampled at 0.8 rather than 0.5 pixels at 51°).
    px_size_m = if abs(raster.gt[2]) < 1.0
        d2m = _EARTH_R * π / 180.0
        min(abs(raster.gt[2]) * d2m * cosd(lat_ref), abs(raster.gt[6]) * d2m)
    else
        min(abs(raster.gt[2]), abs(raster.gt[6]))
    end
    spacing = something(edge_spacing_m, 0.5 * px_size_m)
    @info "  DEM pixel ≈ $(round(px_size_m, digits=2)) m, face sample spacing = $(round(spacing, digits=2)) m"

    # DEM elevations at lon/lat points (NaN: no data)
    function _sample(lons::Vector{Float64}, lats::Vector{Float64})
        np = length(lons)
        xs, ys = if isempty(raster.proj)
            lons, lats
        else
            x = copy(lons); y = copy(lats); zz = zeros(np)
            ArchGDAL.createcoordtrans(_crs_gis(4326), _crs_gis(raster.proj)) do xform
                ArchGDAL.transform!(x, y, zz, xform)
            end
            x, y
        end
        z = Vector{Float64}(undef, np)
        for j in 1:np
            col_f, row_f = _affine_to_pixel(xs[j], ys[j], raster.gt)
            v = _bilinear(raster.band, col_f, row_f, raster.nx, raster.ny, raster.nodata_val)
            z[j] = isfinite(v) ? v : NaN
        end
        return z
    end

    # Neighbour lists define the slots
    if isempty(mesh.adjacency)      # A5 only (checked above)
        @info "  No adjacency on the mesh — computing it via grid_disk..."
        gd = grid_disk_neighbours_batch([c.id for c in mesh.cells])
        in_mesh_ids = Set(c.id for c in mesh.cells)
        for c in mesh.cells
            mesh.adjacency[c.id] = [nb for nb in get(gd, c.id, String[]) if nb in in_mesh_ids]
        end
    end
    _nbrs(c) = get(mesh.adjacency, c.id, get(mesh.adjacency, _to_hex_safe(c.id), String[]))
    n_slots = max(1, maximum(length(_nbrs(c)) for c in mesh.cells))
    nf = SGS_FACE_LEVELS

    # pixel sampling: each cell's block of the DEM's pixels (an error unless
    # every cell is one); the bins must hold every pixel elevation of a block
    pix = sampling === :pixels ? _pixel_blocks(mesh, raster) : nothing
    if pix !== nothing && n_bins < pix.max_pixels + 2
        @info "  Pixel sampling: $(n_bins) bins raised to $(pix.max_pixels + 2) (pixels per block + 2)"
        n_bins = pix.max_pixels + 2
    end

    elev_bins_mat  = Matrix{Float64}(undef, n_bins, n)
    vol_curve_mat  = Matrix{Float64}(undef, n_bins, n)
    area_curve_mat = Matrix{Float64}(undef, n_bins, n)
    wet_cx_mat     = zeros(Float64, n_bins, n)
    wet_cy_mat     = zeros(Float64, n_bins, n)
    edge_sills_mat = fill(NaN, n_slots, n)
    face_zr_mat    = zeros(Float64, n_slots, n)
    face_raise_mat = zeros(Float64, n_slots, n)
    face_A_mat     = zeros(Float64, nf, n_slots, n)
    face_P_mat     = zeros(Float64, nf, n_slots, n)
    face_K_mat     = zeros(Float64, nf, n_slots, n)
    face_C_mat     = zeros(Float64, nf, n_slots, n)
    cell_areas     = Vector{Float64}(undef, n)
    z_mins         = fill(NaN, n)
    z_maxs         = fill(NaN, n)
    n_grid         = zeros(Float64, n)
    # per cell, per polygon side: profile elevations, activation levels in
    # this cell, width (the sides in a canonical direction, so that the two
    # cells of a face sample the same points)
    side_z   = Vector{Vector{Vector{Float64}}}(undef, n)
    side_act = Vector{Vector{Vector{Float64}}}(undef, n)
    side_W   = Vector{Vector{Float64}}(undef, n)
    # pixel sampling: a boundary side takes the pixels inside (side_z holds
    # max(inside, outside), the bed of the face between two cells)
    side_zin   = Vector{Vector{Vector{Float64}}}(undef, n)
    side_actin = Vector{Vector{Vector{Float64}}}(undef, n)

    # ── Step 1: per cell — working grid, hypsometry, connectivity, sides ──
    @info "  Cells: working grid, hypsometric curves, connectivity..."
    d2m = _EARTH_R * π / 180.0
    pix === nothing ||
        _sgs_pixel_cells!(mesh, raster, pix, connectivity, elev_bins_mat, vol_curve_mat, area_curve_mat,
                          wet_cx_mat, wet_cy_mat, cell_areas, z_mins, z_maxs, n_grid,
                          side_z, side_act, side_W, side_zin, side_actin)
    pix === nothing && Threads.@threads for ci in 1:n      # bilinear sampling
        cell = mesh.cells[ci]
        bnd  = cell.boundary
        lonc = Float64(cell.center_lon); latc = Float64(cell.center_lat)
        mx = d2m * cosd(latc); my = d2m
        cell_areas[ci] = _polygon_area_m2(bnd)
        sides = _polygon_sides(bnd)
        ns = length(sides)

        # sides: profiles (midpoint samples, canonical direction)
        sz = Vector{Vector{Float64}}(undef, ns); sW = zeros(ns)
        slon = Vector{Vector{Float64}}(undef, ns); slat = Vector{Vector{Float64}}(undef, ns)
        for (k, (va, vb)) in enumerate(sides)
            a1, a2, b1, b2 = Float64(va[1]), Float64(va[2]), Float64(vb[1]), Float64(vb[2])
            (a1, a2) > (b1, b2) && ((a1, a2, b1, b2) = (b1, b2, a1, a2))
            W = _haversine_m(a1, a2, b1, b2); sW[k] = W
            m = max(2, ceil(Int, W / spacing))
            f = [(q - 0.5) / m for q in 1:m]
            slon[k] = a1 .+ f .* (b1 - a1); slat[k] = a2 .+ f .* (b2 - a2)
            sz[k] = _sample(slon[k], slat[k])
        end
        side_z[ci] = sz; side_W[ci] = sW
        side_act[ci] = [copy(v) for v in sz]

        # working grid over the bounding box, local metric frame at the cell centre
        px = [(Float64(v[1]) - lonc) * mx for v in bnd]; py = [(Float64(v[2]) - latc) * my for v in bnd]
        xmin, xmax = extrema(px); ymin, ymax = extrema(py)
        wx = max(xmax - xmin, 1e-6); wy = max(ymax - ymin, 1e-6)
        h = min(px_size_m, sqrt(wx * wy / n_samples))
        h = max(h, sqrt(wx * wy / 250_000))                    # at most 250 000 points
        nxg = max(1, ceil(Int, wx / h)); nyg = max(1, ceil(Int, wy / h))
        dxg = wx / nxg; dyg = wy / nyg
        glon = Float64[]; glat = Float64[]; gi = Int32[]; gj = Int32[]
        for j in 1:nyg, i in 1:nxg
            lo = lonc + (xmin + (i - 0.5) * dxg) / mx; la = latc + (ymin + (j - 0.5) * dyg) / my
            _point_in_ring(lo, la, bnd) || continue
            push!(glon, lo); push!(glat, la); push!(gi, i); push!(gj, j)
        end
        gz = _sample(glon, glat)
        keep = findall(isfinite, gz)
        zg = gz[keep]; gi = gi[keep]; gj = gj[keep]
        gx = [xmin + (Int(gi[p]) - 0.5) * dxg for p in eachindex(zg)]
        gy = [ymin + (Int(gj[p]) - 0.5) * dyg for p in eachindex(zg)]
        np = length(zg); n_grid[ci] = np
        if np == 0
            elev_bins_mat[:, ci] .= NaN; vol_curve_mat[:, ci] .= 0.0; area_curve_mat[:, ci] .= 0.0
            continue
        end

        # storage samples: four per grid cell (quarter points). The grid itself
        # can line up with the DEM's lattice (a 0.5 m grid on a 0.5 m DEM sees
        # the foot and the top of a one-pixel bank and nothing between), which
        # biases the volume at low levels; connectivity keeps the grid.
        hlon = Vector{Float64}(undef, 4np); hlat = Vector{Float64}(undef, 4np)
        for p in 1:np, (q, (sx_, sy_)) in enumerate(((-0.25, -0.25), (0.25, -0.25), (-0.25, 0.25), (0.25, 0.25)))
            hlon[4(p-1)+q] = lonc + (gx[p] + sx_ * dxg) / mx
            hlat[4(p-1)+q] = latc + (gy[p] + sy_ * dyg) / my
        end
        zh = filter(isfinite, _sample(hlon, hlat))
        isempty(zh) && (zh = copy(zg))

        # shared extremes: the grid, the storage samples and the side profiles
        z_min = min(minimum(zg), minimum(zh)); z_max = max(maximum(zg), maximum(zh))
        for v in sz, zz in v
            isfinite(zz) || continue
            zz < z_min && (z_min = zz); zz > z_max && (z_max = zz)
        end
        z_mins[ci] = z_min; z_maxs[ci] = z_max

        # connectivity to the main basin
        conn = zg
        seed = argmin(zg)
        if connectivity
            idx = zeros(Int32, nxg, nyg)
            for p in 1:np; idx[gi[p], gj[p]] = p; end
            hg = max(dxg, dyg)
            smask = zeros(UInt32, np); touched = UInt32(0)
            for p in 1:np, (k, (va, vb)) in enumerate(sides)
                k > 32 && break
                d = _seg_dist(gx[p], gy[p], (Float64(va[1]) - lonc) * mx, (Float64(va[2]) - latc) * my,
                              (Float64(vb[1]) - lonc) * mx, (Float64(vb[2]) - latc) * my)
                if d <= 0.75 * hg
                    smask[p] |= UInt32(1) << (k - 1); touched |= UInt32(1) << (k - 1)
                end
            end
            seed, conn = _cell_connectivity(zg, idx, gi, gj, smask, min(2, count_ones(touched)))
            # activation level of each side point in this cell: the easiest
            # connection among the grid points next to it
            for k in 1:ns
                act = side_act[ci][k]
                for q in eachindex(act)
                    zp = sz[k][q]; isfinite(zp) || continue
                    x = (slon[k][q] - lonc) * mx; y = (slat[k][q] - latc) * my
                    i0 = clamp(floor(Int, (x - xmin) / dxg) + 1, 1, nxg)
                    j0 = clamp(floor(Int, (y - ymin) / dyg) + 1, 1, nyg)
                    best = Inf
                    for jj in max(1, j0 - 1):min(nyg, j0 + 1), ii in max(1, i0 - 1):min(nxg, i0 + 1)
                        p = Int(idx[ii, jj]); p == 0 && continue
                        hypot(gx[p] - x, gy[p] - y) <= 1.5 * hg || continue
                        c = max(zp, conn[p]); c < best && (best = c)
                    end
                    act[q] = isfinite(best) ? best : zp
                end
            end
        end

        # hypsometric curve (level pool over the whole cell) and wetted-area centroid
        bins = range(z_min, z_max, length = n_bins)
        A_cell = cell_areas[ci]; nh = length(zh); dA = A_cell / nh
        zs = sort(zh)
        oc = sortperm(conn)
        cum_area = 0.0; cum_vol = 0.0; ptr = 1; n_in = 0; prev = z_min
        cptr = 1; sx = 0.0; sy = 0.0; nc = 0
        for k in 1:n_bins
            w = Float64(bins[k])
            cum_vol += n_in * dA * (w - prev)
            while ptr <= nh && zs[ptr] <= w
                cum_area += dA; cum_vol += dA * (w - zs[ptr]); n_in += 1; ptr += 1
            end
            while cptr <= np && conn[oc[cptr]] <= w
                p = oc[cptr]; sx += gx[p]; sy += gy[p]; nc += 1; cptr += 1
            end
            elev_bins_mat[k, ci] = w; area_curve_mat[k, ci] = cum_area; vol_curve_mat[k, ci] = cum_vol
            wet_cx_mat[k, ci] = nc > 0 ? sx / nc : gx[seed]
            wet_cy_mat[k, ci] = nc > 0 ? sy / nc : gy[seed]
            prev = w
        end
    end

    # ── Step 2: faces and boundary sides ──────────────────────────────────
    @info "  Faces: activation levels and face tables..."
    id_to_idx = Dict{String,Int}(cell.id => i for (i, cell) in enumerate(mesh.cells))
    for (i, cell) in enumerate(mesh.cells)
        id_to_idx[_to_hex_safe(cell.id)] = i
    end
    # index of the polygon side of cell c that is the edge (lon1, lat1, lon2, lat2)
    _side_of(c, edge) = findfirst(s -> _side_is_shared(s[1], s[2], (edge,)), _polygon_sides(mesh.cells[c].boundary))
    function _store!(slot, ci, z, act, W)
        ok = [q for q in eachindex(z) if isfinite(z[q])]
        length(ok) >= 2 || return false
        z0, zr, fa, fp, fk, fc = face_table_from_profile(z[ok], W / length(ok); act = act[ok], nf = nf)
        edge_sills_mat[slot, ci] = z0
        face_zr_mat[slot, ci]    = zr
        face_raise_mat[slot, ci] = z0 - minimum(z[ok])
        face_A_mat[:, slot, ci] .= fa; face_P_mat[:, slot, ci] .= fp; face_K_mat[:, slot, ci] .= fk
        face_C_mat[:, slot, ci] .= fc
        return true
    end
    n_raised = 0; n_faces = 0; max_raise = 0.0
    for ci in 1:n
        cell_i = mesh.cells[ci]
        nbrs   = _nbrs(cell_i)
        in_mesh = falses(length(nbrs))
        used = falses(length(side_z[ci]))                  # sides that are faces
        for (slot, nb_id) in enumerate(nbrs)               # slot = position in the neighbour list
            cj = get(id_to_idx, nb_id, get(id_to_idx, _to_hex_safe(nb_id), 0))
            cj == 0 && continue
            in_mesh[slot] = true
            edge = _shared_edge(cell_i.boundary, mesh.cells[cj].boundary)
            ki = edge === nothing ? nothing : _side_of(ci, edge)
            kj = edge === nothing ? nothing : _side_of(cj, edge)
            if ki === nothing || kj === nothing
                zi = z_mins[ci]; zj = z_mins[cj]           # no common side found: no table
                edge_sills_mat[slot, ci] = isnan(zi) || isnan(zj) ? NaN : min(zi, zj)
                continue
            end
            used[ki] = true
            slot_nb = findfirst(id -> get(id_to_idx, id, get(id_to_idx, _to_hex_safe(id), 0)) == ci,
                                _nbrs(mesh.cells[cj]))
            if ci > cj && slot_nb !== nothing && face_P_mat[nf, slot_nb, cj] > 0.0
                # one table per face: computed from the lower-index cell, copied here
                edge_sills_mat[slot, ci] = edge_sills_mat[slot_nb, cj]
                face_zr_mat[slot, ci]    = face_zr_mat[slot_nb, cj]
                face_raise_mat[slot, ci] = face_raise_mat[slot_nb, cj]
                face_A_mat[:, slot, ci] .= view(face_A_mat, :, slot_nb, cj)
                face_P_mat[:, slot, ci] .= view(face_P_mat, :, slot_nb, cj)
                face_K_mat[:, slot, ci] .= view(face_K_mat, :, slot_nb, cj)
                face_C_mat[:, slot, ci] .= view(face_C_mat, :, slot_nb, cj)
                continue
            end
            z  = side_z[ci][ki]; ai = side_act[ci][ki]; aj = side_act[cj][kj]
            m  = length(z)
            act = Vector{Float64}(undef, m)
            for q in 1:m
                qj = length(aj) == m ? q : clamp(round(Int, (q - 0.5) / m * length(aj) + 0.5), 1, length(aj))
                a = max(ai[q], isfinite(aj[qj]) ? aj[qj] : ai[q])
                act[q] = isfinite(z[q]) ? max(z[q], a) : NaN
            end
            if _store!(slot, ci, z, act, side_W[ci][ki])
                n_faces += 1
                r = face_raise_mat[slot, ci]
                r > 0.05 && (n_raised += 1); r > max_raise && (max_raise = r)
            end
        end
        # boundary sides: the slots without an in-mesh neighbour, in polygon-side order
        free = boundary_side_slots(n_slots, in_mesh)
        b = 0
        for k in eachindex(used)
            used[k] && continue
            b += 1
            b > length(free) && break
            if pix === nothing
                _store!(free[b], ci, side_z[ci][k], side_act[ci][k], side_W[ci][k])
            else
                _store!(free[b], ci, side_zin[ci][k], side_actin[ci][k], side_W[ci][k])
            end
        end
    end

    # ── Store results ─────────────────────────────────────────────────────
    mesh.array_vars["sgs_elev_bins"]  = elev_bins_mat
    mesh.array_vars["sgs_vol_curve"]  = vol_curve_mat
    mesh.array_vars["sgs_area_curve"] = area_curve_mat
    mesh.array_vars["sgs_wet_cx"]     = wet_cx_mat
    mesh.array_vars["sgs_wet_cy"]     = wet_cy_mat
    mesh.array_vars["sgs_edge_sills"] = edge_sills_mat
    mesh.array_vars["sgs_face_zr"]    = face_zr_mat
    mesh.array_vars["sgs_face_raise"] = face_raise_mat
    mesh.array_vars["sgs_face_A"]     = reshape(face_A_mat, nf*n_slots, n)
    mesh.array_vars["sgs_face_P"]     = reshape(face_P_mat, nf*n_slots, n)
    mesh.array_vars["sgs_face_K"]     = reshape(face_K_mat, nf*n_slots, n)
    mesh.array_vars["sgs_face_C"]     = reshape(face_C_mat, nf*n_slots, n)

    mesh.static_vars["sgs_cell_area"]    = cell_areas
    mesh.static_vars["sgs_z_min"]        = z_mins
    mesh.static_vars["sgs_z_max"]        = z_maxs
    mesh.static_vars["sgs_n_bins"]       = fill(Float64(n_bins), n)
    mesh.static_vars["sgs_layout"]       = fill(4.0, n)   # 4 (0.8.0): + critical integral, pixel-centre DEM; 3: face tables; 2: per-cell edge curves
    mesh.static_vars["sgs_connectivity"] = fill(connectivity ? 1.0 : 0.0, n)
    mesh.static_vars["sgs_n_grid"]       = n_grid
    mesh.static_vars["sgs_sampling"]     = fill(pix === nothing ? 0.0 : 1.0, n)   # 1: pixels (0.9.1)

    n_ok = count(isfinite, z_mins)
    @info "SGS tables built in $(round(time()-t0, digits=1))s: $n_ok / $n cells with valid curves; " *
          "grid points per cell: median $(round(Int, sort(n_grid)[max(1, n ÷ 2)]))"
    connectivity && @info "  Connectivity: $n_raised of $n_faces faces have their sill raised by more than " *
          "5 cm (largest $(round(max_raise, digits=2)) m) — see sgs_face_raise"
    return mesh
end

build_sgs_tables!(mesh::A5Mesh, dem_path::String; kwargs...) =
    build_sgs_tables!(mesh, FileDEM(dem_path); kwargs...)

"""
    sgs_table(mesh, i) → SGSTable

Reconstruct an `SGSTable` for cell index `i` from the stored array_vars.
Used by the flow solver; called once per cell during initialisation to
build a cached vector of SGSTable objects.
"""
function sgs_table(mesh::A5Mesh, i::Int)::SGSTable
    haskey(mesh.array_vars, "sgs_elev_bins") ||
        error("SGS tables not built — call build_sgs_tables! first.")
    n_bins = size(mesh.array_vars["sgs_elev_bins"], 1)

    haskey(mesh.array_vars, "sgs_face_A") ||
        error("This mesh's SGS tables are older than layout 3 (no face tables). " *
              "Rebuild them: --meshgen … --flow-model sgs, or build_sgs_tables!.")
    n_slots = size(mesh.array_vars["sgs_edge_sills"], 1)
    nf      = size(mesh.array_vars["sgs_face_A"], 1) ÷ n_slots
    fcol(key) = reshape(mesh.array_vars[key][:, i], nf, n_slots)

    # WP0 fix (2026-10): the previous `get(dict, key, fill(NaN, n))[i]` form
    # evaluated the `fill` default EAGERLY on every call, i.e. three length-n
    # allocations per cell → O(n²) work and memory during SGS initialisation
    # (~21 GB of allocation at 30k cells). Look the scalar up directly instead.
    _sv(key) = haskey(mesh.static_vars, key) ? Float64(mesh.static_vars[key][i]) : NaN

    SGSTable(
        mesh.array_vars["sgs_elev_bins"][:, i],
        mesh.array_vars["sgs_vol_curve"][:, i],
        mesh.array_vars["sgs_area_curve"][:, i],
        _sv("sgs_cell_area"),
        _sv("sgs_z_min"),
        _sv("sgs_z_max"),
        mesh.array_vars["sgs_edge_sills"][:, i],
        mesh.array_vars["sgs_face_zr"][:, i],
        fcol("sgs_face_A"), fcol("sgs_face_P"), fcol("sgs_face_K"),
        haskey(mesh.array_vars, "sgs_wet_cx") ? mesh.array_vars["sgs_wet_cx"][:, i] : zeros(n_bins),
        haskey(mesh.array_vars, "sgs_wet_cy") ? mesh.array_vars["sgs_wet_cy"][:, i] : zeros(n_bins),
        haskey(mesh.array_vars, "sgs_face_C") ? fcol("sgs_face_C") :
            _face_C_uniform(fcol("sgs_face_A"), fcol("sgs_face_P")),
    )
end

end # module A5Grid
