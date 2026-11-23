# mesh/DGGSGrid.jl
# ----------------
# Discrete-global-grid meshes from the command line (WP3 part 3, Phase B;
# 0.8.1): `--meshgen-dggs AOI --dggs FAMILY (--meshres N | --dx M)`.
#
# A thin wrapper. The meshes are built by mesh/dggs_families.py (H3 with h3;
# IGEO7, ISEA3H, HEALPix and rHEALPix with DGGAL), which writes a GeoParquet
# file in the mesh contract (docs/MESHES.md §1); this file runs it with
# FloodA5's Python (PyCall's, as the A5 bridge) and loads the result. `a5` is
# accepted as a family too and goes to the A5 generator (A5Grid.mesh_for_aoi),
# so one flag selects any global grid.
#
# Plain include (as mesh/CartesianGrid.jl): defines functions in Main.

"Discrete global grids `--dggs` accepts (`a5` through pya5, the rest through mesh/dggs_families.py)."
const DGGS_FAMILIES = ("a5", "h3", "igeo7", "isea3h", "healpix", "rhealpix")

const _DGGS_SCRIPT = joinpath(@__DIR__, "dggs_families.py")

"""
    a5_level_for_dx(dx_m) → Int

The A5 resolution whose mean cell area is nearest dx² in log terms (the rule
mesh/dggs_families.py uses for the other grids): resolution r ≥ 1 has
60·4^(r−1) cells over the sphere (pya5's authalic radius, 6 371 007.2 m).
Resolution 18 has 22.2 m cells.
"""
function a5_level_for_dx(dx_m::Real)::Int
    dx_m > 0 || error("--dx must be positive (got $dx_m)")
    R = 6_371_007.2
    area(r) = 4π * R^2 / (60 * 4.0^(r - 1))
    return argmin(r -> abs(log(area(r)) - 2log(dx_m)), 1:30)
end

"""
    dggs_generate(family, aoi_path; level = nothing, target_dx = nothing,
                  output_path = nothing) → A5Mesh

A mesh of discrete-global-grid `family` ("h3", "igeo7", "isea3h", "healpix",
"rhealpix") over the GeoJSON AOI: every cell that intersects it, at `level`,
or at the level whose cell area is nearest `target_dx`² (exactly one of the
two). Written to `output_path` (GeoParquet) if given, else to a temporary file,
and loaded. The file records the family, library and version, level,
refinement ratio and whether the grid is equal-area (`mesh.meta`).
"""
function dggs_generate(family::AbstractString, aoi_path::AbstractString;
                       level::Union{Integer,Nothing} = nothing,
                       target_dx::Union{Real,Nothing} = nothing,
                       output_path::Union{AbstractString,Nothing} = nothing)
    family in DGGS_FAMILIES && family != "a5" ||
        error("unknown DGGS family \"$family\"; expected one of $(join(DGGS_FAMILIES[2:end], ", "))")
    (level === nothing) == (target_dx === nothing) &&
        error("give exactly one of a level (--meshres) or a target cell size (--dx)")
    isfile(aoi_path) || error("AOI file not found: $aoi_path")
    if output_path !== nothing
        endswith(lowercase(output_path), ".parquet") ||
            error("--meshout for a DGGS mesh must be a .parquet file (got $output_path)")
    end
    A5Grid._find_python() || error("Python not found — run setup.jl to configure the environment.")
    out = output_path === nothing ? joinpath(mktempdir(), "$(family).parquet") : abspath(output_path)
    mkpath(dirname(out))
    args = [_DGGS_SCRIPT, family, "--aoi", abspath(aoi_path), "--out", out]
    level === nothing ? append!(args, ["--target-dx", string(Float64(target_dx))]) :
                        append!(args, ["--res", string(level)])
    obuf = IOBuffer(); ebuf = IOBuffer()
    ok = success(pipeline(ignorestatus(Cmd(Cmd([A5Grid.PYTHON_CMD[], args...]); dir = dirname(_DGGS_SCRIPT))),
                          stdout = obuf, stderr = ebuf))
    lines = filter(l -> !isempty(strip(l)), split(String(take!(obuf)), '\n'))
    res = isempty(lines) ? nothing : (try JSON3.read(lines[end], Dict{String,Any}) catch; nothing end)
    if !ok || res === nothing || haskey(res, "error")
        msg = res !== nothing && haskey(res, "error") ? res["error"] : String(take!(ebuf))
        error("DGGS mesh generation ($family) failed:\n$msg")
    end
    @info "DGGS mesh: $family level $(res["level"]) ($(get(res, "library", "?"))), " *
          "$(res["cells"]) cells, mean cell $(round(sqrt(res["mean_area_m2"]); sigdigits = 4)) m"
    output_path === nothing || A5Grid._write_prj_sidecar(out)
    return load_mesh_geoparquet(out)
end
