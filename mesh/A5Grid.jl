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

export A5Mesh, A5Cell
export DEMSource, FileDEM
export mesh_for_aoi, save_mesh, save_mesh_geoparquet, save_mesh_geojson, load_mesh_geoparquet
export mesh_to_geojson_string
export sample_dem!, sample_dem_centroid!, sample_dem_mean!, check_dem_coverage
export build_sgs_tables!, SGSTable
export wse_from_volume, wetted_area_from_wse, sgs_table
export flow_area_from_wse, wetted_perim_from_wse, hydraulic_radius_from_wse
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
const _module_init_done = Ref{Bool}(false)

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
end

# Backwards-compatible constructors
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

Base.length(m::A5Mesh) = length(m.cells)
Base.iterate(m::A5Mesh, s=1) = s > length(m.cells) ? nothing : (m.cells[s], s+1)
Base.show(io::IO, m::A5Mesh) =
    print(io, "A5Mesh(resolution=$(m.resolution), cells=$(length(m.cells)), " *
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
(col_f, row_f) (0-based).  Returns NaN for out-of-bounds or nodata corners;
falls back to nearest-neighbour if any of the 2×2 window corners are nodata.
"""
@inline function _bilinear(band::Matrix{Float64},
                            col_f::Float64, row_f::Float64,
                            nx::Int, ny::Int,
                            nodata_val::Float64)::Float64
    (col_f < 0 || col_f >= nx || row_f < 0 || row_f >= ny) && return NaN

    c0 = clamp(floor(Int, col_f), 0, nx - 2)
    r0 = clamp(floor(Int, row_f), 0, ny - 2)
    fc = col_f - c0
    fr = row_f - r0

    # ArchGDAL band array: band[col+1, row+1] (1-based, column-major)
    v00 = band[c0+1, r0+1];  v10 = band[c0+2, r0+1]
    v01 = band[c0+1, r0+2];  v11 = band[c0+2, r0+2]

    nd = nodata_val
    isnd(v) = !isnan(nd) && v == nd

    if isnd(v00) || isnd(v10) || isnd(v01) || isnd(v11)
        # Fall back to nearest neighbour — avoids propagating nodata into
        # an otherwise valid interpolation neighbourhood
        v = band[c0 + round(Int, fc) + 1, r0 + round(Int, fr) + 1]
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

function __init__()
    _module_init_done[] && return
    _module_init_done[] = true
    try
        copy!(_a5, pyimport("a5"))
    catch e
        error("Failed to import pya5: $e\nRun: pip install pya5")
    end
    _try_load_cuda()
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
    expanded   = @py _a5.uncompact(disk, resolution)
    return [_to_hex(c) for c in expanded if _to_hex(c) != cell_id_hex]
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
    use_gpu = _cuda_available[] && length(lons) >= CUDA_PIP_THRESHOLD[]
    if use_gpu
        t    = geometry["type"]
        ring = t == "Polygon" ? geometry["coordinates"][1] :
                                geometry["coordinates"][1][1]
        rxs  = Float64[Float64(v[1]) for v in ring]
        rys  = Float64[Float64(v[2]) for v in ring]
        return _pip_mask_gpu(lons, lats, rxs, rys)
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
    use_gpu = _cuda_available[] && length(lons) >= CUDA_PIP_THRESHOLD[]
    backend = use_gpu ? "GPU" : "CPU ($(Threads.nthreads()) threads)"
    @info "  Sample points: $(length(lons))  |  Backend: $backend (threshold=$(CUDA_PIP_THRESHOLD[]) pts)"

    mask      = _pip_mask(lons, lats, geometry)
    kept_lons = lons[mask]
    kept_lats = lats[mask]
    @info "  Points inside AOI: $(length(kept_lons)) — handing off to Python..."

    tmp_aoi = joinpath(BRIDGE_DIR, "_tmp_aoi.geojson")
    ext     = format == :geoparquet ? ".parquet" : ".geojson"
    tmp_out = joinpath(BRIDGE_DIR, "_tmp_mesh" * ext)
    fmt_str = format == :geoparquet ? "geoparquet" : "geojson"

    try
        write(tmp_aoi, geojson_str)
        result = _run_bridge("mesh_for_aoi", "_tmp_aoi.geojson", string(resolution), "_tmp_mesh" * ext, fmt_str)
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
    path_fwd = replace(abspath(path), "\\" => "/")
    py_code = """
import geopandas as gpd, json, math
import numpy as np
gdf = gpd.read_parquet('$(path_fwd)')

BASE_COLS = {'cell_id', 'resolution', 'center_lon', 'center_lat', 'geometry'}
extra_cols = [c for c in gdf.columns if c not in BASE_COLS]

# Detect column types: string-list (neighbours), numeric list, or scalar
scalar_cols  = []
array_cols   = []
strlist_cols = []  # e.g. 'neighbours'
for col in extra_cols:
    sample = gdf[col].iloc[0] if len(gdf) > 0 else None
    if hasattr(sample, '__len__') and not isinstance(sample, str):
        # List column — check element type
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

cells = []
adjacency = {}  # cell_id -> list of neighbour hex IDs
for _, row in gdf.iterrows():
    geom   = row.geometry
    coords = list(geom.exterior.coords)
    cell_id = row['cell_id']
    cell   = {
        'cell_id':    cell_id,
        'resolution': int(row['resolution']),
        'center_lon': float(row['center_lon']),
        'center_lat': float(row['center_lat']),
        'boundary':   [[c[0], c[1]] for c in coords],
    }
    for col in scalar_cols:
        v = row[col]
        cell[col] = None if (v is None or (isinstance(v, float) and math.isnan(v))) else float(v)
    for col in array_cols:
        v = row[col]
        if v is None:
            cell[col] = None
        else:
            arr = np.asarray(v, dtype=np.float64)
            cell[col] = [None if (x != x or x == float('inf') or x == float('-inf')) else x
                         for x in arr.tolist()]
    for col in strlist_cols:
        v = row[col]
        if col == 'neighbours' and v is not None:
            adjacency[cell_id] = [str(x) for x in v]
        # other string-list cols ignored for now
    cells.append(cell)

print(json.dumps({'cells': cells, 'scalar_cols': scalar_cols,
                  'array_cols': array_cols, 'adjacency': adjacency}))
"""
    tmp_script = tempname() * ".py"
    write(tmp_script, py_code)
    buf = IOBuffer()
    errbuf = IOBuffer()
    try
        run(pipeline(Cmd(Cmd([PYTHON_CMD[], tmp_script]); dir=BRIDGE_DIR),
                     stdout=buf, stderr=errbuf))
    catch e
        error("Failed to read GeoParquet: $(String(take!(errbuf)))")
    finally
        isfile(tmp_script) && rm(tmp_script)
    end

    raw        = JSON3.read(String(take!(buf)), Dict)
    raw_cells  = raw["cells"]
    scalar_cols = [String(k) for k in raw["scalar_cols"]]
    array_cols  = [String(k) for k in raw["array_cols"]]
    raw_adj    = get(raw, "adjacency", Dict())
    nc = length(raw_cells)

    # Build static_vars (scalar) accumulator
    static_vars = Dict{String, Vector{Float64}}(
        col => Vector{Float64}(undef, nc) for col in scalar_cols
    )

    # We need to know each array column's length to allocate the matrix.
    # Peek at the first cell.
    array_lengths = Dict{String,Int}()
    if nc > 0
        for col in array_cols
            v = get(raw_cells[1], col, nothing)
            if v !== nothing && !isempty(v)
                array_lengths[col] = length(v)
            end
        end
    end
    array_vars = Dict{String, Matrix{Float64}}(
        col => Matrix{Float64}(undef, array_lengths[col], nc)
        for col in array_cols if haskey(array_lengths, col)
    )

    cells = A5Cell[]
    for (i, c) in enumerate(raw_cells)
        elev = haskey(c, "elevation") ?
               (c["elevation"] === nothing ? NaN : Float64(c["elevation"])) : NaN
        push!(cells, A5Cell(
            String(c["cell_id"]), Int(c["resolution"]),
            Float64(c["center_lon"]), Float64(c["center_lat"]),
            [[Float64(v[1]), Float64(v[2])] for v in c["boundary"]],
            elev))
        for col in scalar_cols
            v = get(c, col, nothing)
            static_vars[col][i] = (v === nothing) ? NaN : Float64(v)
        end
        for col in array_cols
            haskey(array_vars, col) || continue
            v = get(c, col, nothing)
            if v === nothing
                array_vars[col][:, i] .= NaN
            else
                for (k, x) in enumerate(v)
                    array_vars[col][k, i] = (x === nothing) ? NaN : Float64(x)
                end
            end
        end
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

    mesh = A5Mesh(res, aoi_geojson, cells, static_vars, array_vars, adjacency)

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
        norm_id = _to_hex(parse(UInt64, c.id, base=16))
        nbrs = get(adj_by_id, norm_id, get(adj_by_id, c.id, nothing))
        if nbrs !== nothing
            d["neighbours"] = nbrs
        end
        push!(cells_data, d)
    end

    tmp_cells  = joinpath(BRIDGE_DIR, "_tmp_save_cells.json")
    tmp_out    = joinpath(BRIDGE_DIR, "_tmp_save.parquet")
    path_fwd   = replace(path, "\\" => "/")
    tmp_fwd    = replace(tmp_out, "\\" => "/")

    try
        open(tmp_cells, "w") do io
            JSON3.write(io, cells_data)
        end

        tmp_cells_fwd = replace(tmp_cells, "\\" => "/")
        py_code = """
import sys, json
sys.argv = ['a5_bridge.py']
import importlib.util, os
spec = importlib.util.spec_from_file_location('a5_bridge', os.path.join(os.path.dirname('$(tmp_cells_fwd)'), 'a5_bridge.py'))
bridge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)

with open('$(tmp_cells_fwd)', 'r') as f:
    cells = json.load(f)
bridge.mesh_to_geoparquet(cells, '$(tmp_fwd)')
print('ok')
"""
        tmp_script = joinpath(BRIDGE_DIR, "_tmp_save.py")
        write(tmp_script, py_code)
        buf = IOBuffer(); errbuf = IOBuffer()
        try
            run(pipeline(Cmd(Cmd([PYTHON_CMD[], "_tmp_save.py"]); dir=BRIDGE_DIR),
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
    println("CUDA_AVAILABLE: ", _cuda_available[])
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
    # Per-edge hydraulic geometry for the R-A flux kernel.
    # Rows = n_bins elevation knots (same as elev_bins).
    # Cols = adjacency slots 1..5 (unused slots hold zeros).
    # edge_area_curves[k, s]  = cross-sectional flow area (m²) through slot-s edge
    #                           at WSE = elev_bins[k].
    # edge_perim_curves[k, s] = wetted perimeter (m) at same WSE.
    # Both are zero when the mesh was built without R-A tables (legacy parquet).
    # The R-A flux kernel (_manning_flux_ra) requires these; falls back to Bates
    # form if they are all zero.
    edge_area_curves :: Matrix{Float64}   # (n_bins × 5)
    edge_perim_curves:: Matrix{Float64}   # (n_bins × 5)
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
@inline function wse_from_volume(t::SGSTable, V::Float64)::Float64
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
    wetted_area_from_wse(table, wse) → Float64

Return the wetted plan area (m²) at a given WSE by interpolating area_curve.
"""
@inline function wetted_area_from_wse(t::SGSTable, wse::Float64)::Float64
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

"""
    flow_area_from_wse(table, slot, wse) → Float64

Cross-sectional flow area (m²) at the shared edge in adjacency `slot` (1–5) at
water surface elevation `wse`, interpolated from the pre-computed edge_area_curves.

Returns 0.0 for dry edges (wse at or below the edge sill elevation, detected as
zero area in the curve) or for legacy meshes without R-A tables (all zeros).
"""
@inline function flow_area_from_wse(t::SGSTable, slot::Int, wse::Float64)::Float64
    col = view(t.edge_area_curves, :, slot)
    # If tables not built (legacy mesh), col is all zeros
    col[end] <= 0.0 && return 0.0
    wse <= t.elev_bins[1]   && return 0.0
    wse >= t.elev_bins[end] && return col[end]
    lo, hi = 1, length(t.elev_bins)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        t.elev_bins[mid] <= wse ? lo = mid : hi = mid
    end
    frac = (wse - t.elev_bins[lo]) / (t.elev_bins[hi] - t.elev_bins[lo] + 1e-15)
    return col[lo] + frac * (col[hi] - col[lo])
end

"""
    wetted_perim_from_wse(table, slot, wse) → Float64

Wetted perimeter (m) at the shared edge in adjacency `slot` (1–5) at WSE `wse`,
interpolated from the pre-computed edge_perim_curves.
"""
@inline function wetted_perim_from_wse(t::SGSTable, slot::Int, wse::Float64)::Float64
    col = view(t.edge_perim_curves, :, slot)
    col[end] <= 0.0 && return 0.0
    wse <= t.elev_bins[1]   && return 0.0
    wse >= t.elev_bins[end] && return col[end]
    lo, hi = 1, length(t.elev_bins)
    while hi - lo > 1
        mid = (lo + hi) >>> 1
        t.elev_bins[mid] <= wse ? lo = mid : hi = mid
    end
    frac = (wse - t.elev_bins[lo]) / (t.elev_bins[hi] - t.elev_bins[lo] + 1e-15)
    return col[lo] + frac * (col[hi] - col[lo])
end

"""
    hydraulic_radius_from_wse(table, slot, wse) → Float64

Hydraulic radius R = A/P (m) at the shared edge in adjacency `slot` at WSE `wse`.
Returns 0.0 when the edge is dry (P < 1e-6 m).
"""
@inline function hydraulic_radius_from_wse(t::SGSTable, slot::Int,
                                            wse::Float64)::Float64
    A = flow_area_from_wse(t, slot, wse)
    P = wetted_perim_from_wse(t, slot, wse)
    P < 1e-6 && return 0.0
    return A / P
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
        → (cos_theta::Float64, skew_x::Float64, skew_y::Float64)

Combined non-orthogonality and skewness geometry for one A5 cell-pair edge.

Supersedes the former standalone `_edge_cos_theta` by computing both the
non-orthogonality cosine *and* the skewness vector from a single local
equirectangular projection of the shared edge and cell centres — avoiding
running the projection setup twice per edge at mesh-build time, and keeping
all edge-geometry maths in one place for future modification.

Returns:

  cos_theta  — |d̂ · n̂|, the magnitude of alignment between the centre-to-
               centre unit vector d̂ (cell i → cell j) and the outward face
               normal n̂ of their shared edge.  1.0 = perfectly orthogonal.
               Used to project the centre-to-centre distance onto the
               face-normal direction in the *uncorrected* flux kernels
               (`_bates_flux`, `_manning_flux_ra`):  L_eff = L × cos θ.

  skew_x,
  skew_y     — components (m, local equirectangular frame) of the vector
               from the shared-edge face midpoint to the point where the
               centre-to-centre line intersects the face.  (0.0, 0.0) for
               an orthogonal, unskewed cell pair; non-zero on skewed A5
               edges.  Used by the WLSQ-corrected flux kernels
               (`_bates_flux_corrected`, `_manning_flux_ra_corrected`) to
               build a skewness-corrected driving head:
                   dWSE_n = (WSE_j − WSE_i) + ∇WSE_f · (skew_x, skew_y)
               See FloodA5_NonOrthogonal_Correction_Plan.md §2, §5.1.

Returns `(1.0, 0.0, 0.0)` if no shared edge is found (non-adjacent cells) —
an orthogonal, unskewed fallback rather than NaN.  A NaN cos_theta would
propagate through the NaN guard in the step functions and suppress all flux
on that edge, which is far worse than a small error in the slope projection.

Returns `(NaN, 0.0, 0.0)` for degenerate edge/distance geometry (near-zero
edge length or near-zero centre separation) — callers already guard NaN
cos_theta by falling back to 1.0, and a zero skewness correction is the
safe choice when the geometry itself is degenerate.

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

Skewness is computed by intersecting the centre-to-centre line with the
shared-edge line (2D line/line intersection via Cramer's rule). If the two
lines are parallel, or the intersection point falls outside the shared edge
segment (a genuinely degenerate skewed case), the skewness correction is
set to zero rather than propagating an unbounded or ill-defined value.
"""
function _edge_geometry(bnd_i  :: Vector{Vector{Float64}},
                         bnd_j  :: Vector{Vector{Float64}},
                         lon_i  :: Float64, lat_i  :: Float64,
                         lon_j  :: Float64, lat_j  :: Float64
                         )::Tuple{Float64,Float64,Float64}

    # ── 1. Find the shared edge vertices ─────────────────────────────────
    edge = _shared_edge(bnd_i, bnd_j)
    # If no shared edge is found, fall back to cos θ = 1.0 (orthogonal
    # assumption) with zero skewness, rather than NaN.  A NaN would
    # propagate through the NaN guard in the step functions and suppress
    # all flux on that edge, which is far worse than a small error in the
    # slope projection.
    edge === nothing && return (1.0, 0.0, 0.0)

    elon1, elat1, elon2, elat2 = edge

    # ── 2. Local equirectangular projection ───────────────────────────────
    # Centre the projection on the edge midpoint to minimise distortion.
    lat0     = (elat1 + elat2) * 0.5
    cos_lat0 = cosd(lat0)
    R        = _EARTH_R

    # Project a lon/lat point to local (x, y) metres
    to_xy(lon, lat) = (deg2rad(lon) * R * cos_lat0,
                       deg2rad(lat) * R)

    # ── 3. Edge vector and face normal ────────────────────────────────────
    x1, y1 = to_xy(elon1, elat1)
    x2, y2 = to_xy(elon2, elat2)
    ex, ey  = x2 - x1, y2 - y1
    e_len   = sqrt(ex^2 + ey^2)
    e_len < 1.0 && return (NaN, 0.0, 0.0)   # degenerate edge guard

    # Face normal: rotate edge vector 90° (two choices; we take the one
    # pointing from i toward j, so the sign of the dot product is positive)
    nx, ny  = -ey / e_len, ex / e_len   # unit normal candidate

    # ── 4. Centre-to-centre unit vector ───────────────────────────────────
    ci_x, ci_y = to_xy(lon_i, lat_i)
    cj_x, cj_y = to_xy(lon_j, lat_j)
    dx, dy      = cj_x - ci_x, cj_y - ci_y
    d_len       = sqrt(dx^2 + dy^2)
    d_len < 1.0 && return (NaN, 0.0, 0.0)   # degenerate distance guard

    dx_hat, dy_hat = dx / d_len, dy / d_len

    # ── 5. cos θ = |d̂ · n̂| ───────────────────────────────────────────────
    # Absolute value: we want the magnitude of the alignment regardless of
    # which side of the edge we are on.
    cos_theta = abs(dx_hat * nx + dy_hat * ny)
    cos_theta = clamp(cos_theta, 0.0, 1.0)   # guard floating-point rounding above 1.0

    # ── 6. Skewness: face midpoint vs. centre-to-centre / face intersection ─
    # The face midpoint is the average of the two shared-edge vertices.
    f_mid_x, f_mid_y = 0.5 * (x1 + x2), 0.5 * (y1 + y2)

    # Solve for the intersection of the centre-to-centre line
    # (ci_x, ci_y) + s·(dx, dy)  with the edge line  (x1, y1) + t·(ex, ey):
    #     ci_x + s·dx = x1 + t·ex
    #     ci_y + s·dy = y1 + t·ey
    # via Cramer's rule on the 2×2 system in (s, t).
    denom = dx * (-ey) - dy * (-ex)   # = -dx·ey + dy·ex
    if abs(denom) < 1e-9
        # Lines parallel / degenerate — no well-defined intersection.
        # Treat as zero skewness rather than propagating NaN/Inf.
        return (cos_theta, 0.0, 0.0)
    end
    rhs_x = x1 - ci_x
    rhs_y = y1 - ci_y
    t = (dx * rhs_y - dy * rhs_x) / denom

    if abs(t) > 1.0
        # Intersection falls outside the shared edge segment (using the
        # edge vector's own parametrisation as the unit scale) — a
        # genuinely degenerate skewed case.  Fall back to zero skewness
        # correction for this edge rather than extrapolating.
        return (cos_theta, 0.0, 0.0)
    end

    d_int_x = x1 + t * ex
    d_int_y = y1 + t * ey

    skew_x = f_mid_x - d_int_x
    skew_y = f_mid_y - d_int_y

    return (cos_theta, skew_x, skew_y)
end

# ---------------------------------------------------------------------------
# SGS pre-processing — build_sgs_tables!
# ---------------------------------------------------------------------------

"""
    build_sgs_tables!(mesh, dem_source;
                      n_bins=100, n_samples=512,
                      halton_seed=0,
                      edge_spacing_m=nothing,
                      var_prefix="sgs")

Pre-compute Sub-Grid Sampling (SGS) hypsometric tables for every cell in
the mesh and store them in `mesh.array_vars`.

For each cell:
1. Generate `n_samples` Halton quasi-random points within the cell polygon.
2. Sample the DEM at all points (bilinear, from in-memory raster).
3. Build an elevation histogram with `n_bins` quantile-spaced bins.
4. Derive cumulative volume and wetted-area curves from the histogram.

For each pair of topologically adjacent cells:
5. Sample the DEM along the shared boundary edge at ~`edge_spacing_m`
   intervals (default: 0.5 × DEM pixel size) to find the sill elevation.
   The sill is stored in `mesh.static_vars["sgs_sill_<i>_<j>"]` — but
   since sills are keyed by cell-pair, they are instead stored in a
   flat `mesh.static_vars["sgs_sill"]` vector for each cell, holding the
   minimum sill elevation to any neighbour (conservative for flow).
   Per-edge sill values are kept in `mesh.array_vars["sgs_edge_sills"]`
   as a (max_neighbours × n_cells) matrix.

Output arrays in `mesh.array_vars`:
  "sgs_elev_bins"   — (n_bins × n_cells)  elevation knot points (m)
  "sgs_vol_curve"   — (n_bins × n_cells)  cumulative volume (m³)
  "sgs_area_curve"  — (n_bins × n_cells)  cumulative wetted area (m²)
  "sgs_edge_sills"  — (5 × n_cells)       sill elevation per neighbour slot

Output scalars in `mesh.static_vars`:
  "sgs_cell_area"   — total plan area of each cell (m²)
  "sgs_z_min"       — minimum DEM elevation in cell
  "sgs_z_max"       — maximum DEM elevation in cell
  "sgs_n_bins"      — n_bins (stored as Float64 for parquet compatibility)

The arrays are automatically persisted when you call `save_mesh_geoparquet`.
"""
function build_sgs_tables!(mesh::A5Mesh, dem_source::DEMSource;
                            n_bins::Int=100,
                            n_samples::Int=512,
                            halton_seed::Int=0,
                            edge_spacing_m::Union{Float64,Nothing}=nothing)
    dem_path = _resolve_dem_path(dem_source, mesh)
    n = length(mesh.cells)
    n == 0 && error("Cannot build SGS tables on an empty mesh.")

    @info "SGS pre-processing: $(basename(dem_path)), $n cells, " *
          "$n_bins bins, $n_samples Halton pts/cell..."
    t0 = time()

    raster = _load_dem_raster(dem_path)

    # Default edge spacing: half the DEM pixel size in metres.
    # gt[2] is the pixel width in the raster's native CRS units.
    # For a projected CRS (e.g. NZTM2000) it is already in metres.
    # For a geographic CRS (degrees) we convert via the Earth radius.
    # Heuristic: if gt[2] < 1.0 it's almost certainly degrees.
    px_size_m = if abs(raster.gt[2]) < 1.0
        abs(raster.gt[2]) * _EARTH_R * π / 180.0
    else
        abs(raster.gt[2])   # already metres
    end
    spacing   = something(edge_spacing_m, 0.5 * px_size_m)
    @info "  DEM pixel ≈ $(round(px_size_m, digits=1)) m, " *
          "edge sample spacing = $(round(spacing, digits=1)) m"

    # Pre-generate Halton sequences (shared across all cells, scaled per bbox)
    h2 = _halton(2, n_samples, halton_seed)
    h3 = _halton(3, n_samples, halton_seed)

    # Output matrices
    elev_bins_mat  = Matrix{Float64}(undef, n_bins, n)
    vol_curve_mat  = Matrix{Float64}(undef, n_bins, n)
    area_curve_mat = Matrix{Float64}(undef, n_bins, n)
    edge_sills_mat = fill(NaN, 5, n)   # up to 5 neighbours per pentagon
    # R-A flux tables: cross-sectional flow area and wetted perimeter per edge slot
    # Shape: (n_bins, 5, n_cells) — stored flat as (n_bins×5, n_cells) in parquet
    edge_area_mat  = zeros(Float64, n_bins, 5, n)   # m²
    edge_perim_mat = zeros(Float64, n_bins, 5, n)   # m
    cell_areas     = Vector{Float64}(undef, n)
    z_mins         = Vector{Float64}(undef, n)
    z_maxs         = Vector{Float64}(undef, n)

    # ── Step 1: Build hypsometric curves for each cell ────────────────────
    @info "  Building hypsometric curves..."
    Threads.@threads for ci in 1:n
        cell = mesh.cells[ci]
        bnd  = cell.boundary

        # Bounding box
        bb_lon_min = minimum(Float64(v[1]) for v in bnd)
        bb_lon_max = maximum(Float64(v[1]) for v in bnd)
        bb_lat_min = minimum(Float64(v[2]) for v in bnd)
        bb_lat_max = maximum(Float64(v[2]) for v in bnd)
        dlon = bb_lon_max - bb_lon_min
        dlat = bb_lat_max - bb_lat_min

        # Halton candidates → PIP filter
        cand_lons = bb_lon_min .+ h2 .* dlon
        cand_lats = bb_lat_min .+ h3 .* dlat
        inside    = [_point_in_ring(cand_lons[j], cand_lats[j], bnd)
                     for j in 1:n_samples]
        in_lons   = cand_lons[inside]
        in_lats   = cand_lats[inside]
        n_in      = length(in_lons)

        cell_areas[ci] = _polygon_area_m2(bnd)

        if n_in == 0
            # Degenerate cell — fill with NaN
            elev_bins_mat[:, ci]  .= NaN
            vol_curve_mat[:, ci]  .= 0.0
            area_curve_mat[:, ci] .= 0.0
            z_mins[ci] = NaN;  z_maxs[ci] = NaN
            continue
        end

        # Reproject to DEM CRS
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

        # Sample elevations at all interior points
        sample_elevs = Vector{Float64}(undef, n_in)
        for j in 1:n_in
            col_f, row_f = _affine_to_pixel(native_xs[j], native_ys[j], raster.gt)
            sample_elevs[j] = _bilinear(raster.band, col_f, row_f,
                                         raster.nx, raster.ny, raster.nodata_val)
        end
        valid_elevs = filter(isfinite, sample_elevs)
        if isempty(valid_elevs)
            elev_bins_mat[:, ci]  .= NaN
            vol_curve_mat[:, ci]  .= 0.0
            area_curve_mat[:, ci] .= 0.0
            z_mins[ci] = NaN;  z_maxs[ci] = NaN
            continue
        end

        sort!(valid_elevs)
        z_min = valid_elevs[1]
        z_max = valid_elevs[end]
        z_mins[ci] = z_min
        z_maxs[ci] = z_max

        # Build n_bins quantile-spaced elevation knots spanning [z_min, z_max]
        bins = range(z_min, z_max, length=n_bins)
        A_cell = cell_areas[ci]
        n_valid = length(valid_elevs)
        dA = A_cell / n_valid   # area element per sample point (uniform weights)

        # We scan through sorted valid_elevs once per bin.
        # For each WSE knot, volume = sum over all inundated samples of
        #   dA × (wse_k - z_sample)
        # This can be maintained incrementally:
        #   when a new sample s crosses below wse_k, add dA*(wse_k - z_s)
        #   for already-inundated samples, add dA*(wse_{k} - wse_{k-1}) to cum_vol
        #
        #   cum_vol[k] = cum_vol[k-1]
        #              + n_newly_inundated * dA * (wse_k - their_elevation)   [new samples]
        #              + n_already_inundated_prev * dA * (wse_k - wse_{k-1}) [depth increase]

        cum_area   = 0.0
        cum_vol    = 0.0
        bin_idx    = 1        # next unprocessed sample index in sorted valid_elevs
        n_inundated = 0       # count of samples already below WSE
        prev_wse   = z_min

        for k in 1:n_bins
            wse_k = Float64(bins[k])

            # Volume increase for samples already inundated in previous bins:
            # each gains (wse_k - prev_wse) × dA additional depth
            cum_vol += n_inundated * dA * (wse_k - prev_wse)

            # Bring in any new samples that fall below wse_k
            while bin_idx <= n_valid && valid_elevs[bin_idx] <= wse_k
                z_s      = valid_elevs[bin_idx]
                cum_area += dA
                cum_vol  += dA * (wse_k - z_s)
                n_inundated += 1
                bin_idx  += 1
            end

            elev_bins_mat[k,  ci] = wse_k
            area_curve_mat[k, ci] = cum_area
            vol_curve_mat[k,  ci] = cum_vol
            prev_wse = wse_k
        end
    end

    # ── Step 2: Edge sill elevations ──────────────────────────────────────
    @info "  Sampling edge sill elevations..."

    # Build a cell_id → index map
    id_to_idx = Dict{String,Int}(cell.id => i for (i, cell) in enumerate(mesh.cells))

    # We need true topological neighbours; build adjacency via grid_disk
    @info "  Building exact adjacency via grid_disk..."
    adj = grid_disk_neighbours_batch([c.id for c in mesh.cells])

    # For each cell, sample the DEM along each shared edge.
    # NOTE: must be serial (not Threads.@threads). The ArchGDAL coordinate transform
    # inside this loop uses PyCall, whose GC finaliser (pydecref) is not thread-safe
    # with Python 3.13 on Windows — a GC cycle on a non-main thread triggers
    # EXCEPTION_ACCESS_VIOLATION in PyObject_ClearWeakRefs (Bug 57, same root cause
    # as Bug 50 in _shared_edge).  The hypsometric curve build (Step 1) dominates
    # wall time; the edge sill loop is fast even serially at res 18 (29,902 cells).
    for ci in 1:n
        cell_i = mesh.cells[ci]
        nbrs   = get(adj, cell_i.id, String[])
        slot   = 0
        for nb_id in nbrs
            slot += 1
            slot > 5 && break
            ci_nb = get(id_to_idx, nb_id, 0)
            ci_nb == 0 && continue

            cell_j  = mesh.cells[ci_nb]
            edge    = _shared_edge(cell_i.boundary, cell_j.boundary)
            if edge === nothing
                # Fallback: use min of the two cell z_mins
                zmin_i = z_mins[ci];  zmin_j = z_mins[ci_nb]
                edge_sills_mat[slot, ci] = isnan(zmin_i) || isnan(zmin_j) ?
                    NaN : min(zmin_i, zmin_j)
                continue
            end

            lon1, lat1, lon2, lat2 = edge

            # Sample DEM along the edge arc
            e_lons, e_lats = _arc_sample_lons_lats(lon1, lat1, lon2, lat2, spacing)
            n_ep = length(e_lons)

            native_xs, native_ys = if isempty(raster.proj)
                e_lons, e_lats
            else
                xs = copy(e_lons);  ys = copy(e_lats);  zs = zeros(n_ep)
                wgs84   = _crs_gis(4326)
                dem_crs = _crs_gis(raster.proj)
                ArchGDAL.createcoordtrans(wgs84, dem_crs) do xform
                    ArchGDAL.transform!(xs, ys, zs, xform)
                end
                xs, ys
            end

            sill = Inf
            valid_edge_elevs = Float64[]
            for j in 1:n_ep
                col_f, row_f = _affine_to_pixel(native_xs[j], native_ys[j], raster.gt)
                v = _bilinear(raster.band, col_f, row_f,
                              raster.nx, raster.ny, raster.nodata_val)
                if isfinite(v)
                    v < sill && (sill = v)
                    push!(valid_edge_elevs, v)
                end
            end
            edge_sills_mat[slot, ci] = isfinite(sill) ? sill : NaN

            # ── R-A edge hydraulic curves ──────────────────────────────────
            # Build cumulative flow area A(wse) and wetted perimeter P(wse)
            # from the sorted edge elevation profile.  Uses the same n_bins
            # elevation knots as the cell hypsometric curve (elev_bins_mat[:, ci]).
            # This keeps A/P lookups at runtime cheap (same binary-search grid).
            if length(valid_edge_elevs) >= 2 && !isnan(z_mins[ci])
                sort!(valid_edge_elevs)
                # Physical edge length and width element per sample
                W  = _haversine_m(lon1, lat1, lon2, lat2)
                dx = W / length(valid_edge_elevs)

                bins     = view(elev_bins_mat, :, ci)
                cum_A    = 0.0
                cum_P    = 0.0
                ptr      = 1
                n_wet    = 0
                prev_wse = valid_edge_elevs[1]   # sill elevation
                ne_pts   = length(valid_edge_elevs)

                for k in 1:n_bins
                    wse_k = bins[k]
                    # depth increase for already-wet arc segments
                    cum_A += n_wet * dx * (wse_k - prev_wse)
                    # bring in newly wet arc segments
                    while ptr <= ne_pts && valid_edge_elevs[ptr] <= wse_k
                        cum_A += dx * (wse_k - valid_edge_elevs[ptr])
                        cum_P += dx
                        n_wet += 1
                        ptr   += 1
                    end
                    edge_area_mat[k,  slot, ci] = cum_A
                    edge_perim_mat[k, slot, ci] = cum_P
                    prev_wse = wse_k
                end
            end
            # If valid_edge_elevs is empty or z_min is NaN, edge_area_mat and
            # edge_perim_mat remain zero for this slot (dry / no-data edge).
        end
    end

    # ── Store results ─────────────────────────────────────────────────────
    mesh.array_vars["sgs_elev_bins"]       = elev_bins_mat
    mesh.array_vars["sgs_vol_curve"]       = vol_curve_mat
    mesh.array_vars["sgs_area_curve"]      = area_curve_mat
    mesh.array_vars["sgs_edge_sills"]      = edge_sills_mat

    # R-A flux tables: reshape (n_bins, 5, n) → (n_bins×5, n) for parquet storage
    mesh.array_vars["sgs_edge_area_curve"]  = reshape(edge_area_mat,  n_bins*5, n)
    mesh.array_vars["sgs_edge_perim_curve"] = reshape(edge_perim_mat, n_bins*5, n)

    mesh.static_vars["sgs_cell_area"] = cell_areas
    mesh.static_vars["sgs_z_min"]     = z_mins
    mesh.static_vars["sgs_z_max"]     = z_maxs
    mesh.static_vars["sgs_n_bins"]    = fill(Float64(n_bins), n)

    n_ok = count(isfinite, z_mins)
    @info "SGS tables built in $(round(time()-t0, digits=1))s: " *
          "$n_ok / $n cells with valid curves"
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

    # R-A edge hydraulic curves: stored as (n_bins×5, n_cells) in parquet,
    # reshaped to (n_bins, 5) per cell.
    # Backward compatibility: meshes built before the R-A tables were added
    # (pre sgs_ra_flux branch) lack these columns.  Fill with zeros and warn
    # once — the R-A kernel will detect all-zero curves and fall back to Bates.
    if haskey(mesh.array_vars, "sgs_edge_area_curve")
        edge_area  = reshape(mesh.array_vars["sgs_edge_area_curve"][:,  i], n_bins, 5)
        edge_perim = reshape(mesh.array_vars["sgs_edge_perim_curve"][:, i], n_bins, 5)
    else
        edge_area  = zeros(Float64, n_bins, 5)
        edge_perim = zeros(Float64, n_bins, 5)
    end

    SGSTable(
        mesh.array_vars["sgs_elev_bins"][:, i],
        mesh.array_vars["sgs_vol_curve"][:, i],
        mesh.array_vars["sgs_area_curve"][:, i],
        get(mesh.static_vars, "sgs_cell_area", fill(NaN, length(mesh.cells)))[i],
        get(mesh.static_vars, "sgs_z_min",     fill(NaN, length(mesh.cells)))[i],
        get(mesh.static_vars, "sgs_z_max",     fill(NaN, length(mesh.cells)))[i],
        edge_area,
        edge_perim,
    )
end

end # module A5Grid
