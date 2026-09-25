#!/usr/bin/env julia
# check_boundary_proximity.jl
# ----------------------------
# Answers a specific question raised in the cartesian-formulation branch's
# testing-strategy review: "the boundaries were closed — could there still
# be effects if water reaches the boundary?" This reports, for a given
# frame, how close the wetted region's extent is to the mesh's own domain
# boundary in each cardinal direction. If the wetted region is closing in
# on (or has reached) a boundary, any symmetry/asymmetry metric measured
# at that frame is potentially contaminated by the closed-boundary
# reflection/ponding regime, not measuring free-spreading behaviour alone
# — see the precedent for this exact caution in
# FloodA5_GradientCorrection_PentagonChirality_Handover.md §4 ("Caveat on
# the late-run 'recovery'").
#
# Usage:
#   julia --project=. test/check_boundary_proximity.jl --input FILE.h5 --frame N
#   julia --project=. test/check_boundary_proximity.jl --input FILE.h5 --sweep 10
#     (prints the wet-region extent at every 10th frame, so you can see
#      whether/when it starts approaching a boundary over the course of
#      the run, without needing to guess a single frame up front)

using HDF5
using Printf

function _argval(args::Vector{String}, flag::String, default=nothing)
    idx = findfirst(==(flag), args)
    idx === nothing && return default
    idx == length(args) && error("$flag requires a value")
    return args[idx + 1]
end

const ARGS_V = String.(ARGS)

function _report_frame(f, fname::String, cell_lons, cell_lats,
                        lon_min, lon_max, lat_min, lat_max)
    grp = f["frames/$fname"]
    sat = read(grp["saturation"])
    t   = haskey(grp, "t") ? read(grp["t"]) : NaN
    wet = sat .>= 1.0
    n_wet = count(wet)
    if n_wet == 0
        @printf("frame %-6s  t=%9.1fs  n_wet=%6d  (no wet cells)\n", fname, t, n_wet)
        return
    end
    wlon = cell_lons[wet]; wlat = cell_lats[wet]
    wlon_min, wlon_max = extrema(wlon)
    wlat_min, wlat_max = extrema(wlat)

    # Distance from the wetted extent to each domain edge, in metres
    # (equirectangular approx at the domain's own mid-latitude — adequate
    # for a proximity check, not a precision measurement).
    R = 6_371_000.0
    cos_lat0 = cosd(0.5 * (lat_min + lat_max))
    d_west  = (wlon_min - lon_min) * cos_lat0 * (π/180) * R
    d_east  = (lon_max - wlon_max) * cos_lat0 * (π/180) * R
    d_south = (wlat_min - lat_min) * (π/180) * R
    d_north = (lat_max - wlat_max) * (π/180) * R

    # Report ALL FOUR distances, not just the overall closest — a fixed
    # point source near one edge (e.g. west, for this project's planar-slope
    # tests) pins that one direction's distance for the entire run,
    # which would otherwise permanently mask whatever is happening at the
    # OTHER boundaries (north/south, the ones actually in question for a
    # north/south symmetry investigation). Printing one number was the bug
    # in the previous version of this script — every frame reported the
    # same ~150m "closest" distance because west was always the minimum,
    # regardless of what north/south were doing.
    flag(d, name) = d < 0 ? "  ⚠ AT/PAST $name BOUNDARY" :
                    d < 100 ? "  ⚠ within 100m of $name" : ""

    @printf("frame %-6s  t=%9.1fs  n_wet=%6d  N=%8.1fm%s  S=%8.1fm%s  E=%8.1fm%s  W=%8.1fm%s\n",
            fname, t, n_wet,
            d_north, flag(d_north, "north"),
            d_south, flag(d_south, "south"),
            d_east,  flag(d_east,  "east"),
            d_west,  flag(d_west,  "west"))
end

function main()
    input_path = _argval(ARGS_V, "--input")
    input_path === nothing && error("--input FILE.h5 is required")

    h5open(input_path, "r") do f
        cell_lons = read(f["mesh/center_lons"])
        cell_lats = read(f["mesh/center_lats"])
        lon_min, lon_max = extrema(cell_lons)
        lat_min, lat_max = extrema(cell_lats)
        println("Domain bounds: lon [$lon_min, $lon_max]  lat [$lat_min, $lat_max]")
        println("(bounds are cell CENTRES, not the outer edge of the boundary cells —")
        println(" true domain edge is up to half a cell-width further out than this)")
        println()

        names = sort(collect(keys(f["frames"])))

        sweep_str = _argval(ARGS_V, "--sweep")
        if sweep_str !== nothing
            stride = parse(Int, sweep_str)
            for (i, nm) in enumerate(names)
                (i - 1) % stride == 0 &&
                    _report_frame(f, nm, cell_lons, cell_lats, lon_min, lon_max, lat_min, lat_max)
            end
        else
            frame_str = _argval(ARGS_V, "--frame")
            frame_str === nothing && error("Provide --frame N or --sweep STRIDE")
            fname = lpad(frame_str, 6, '0')
            _report_frame(f, fname, cell_lons, cell_lats, lon_min, lon_max, lat_min, lat_max)
        end
    end
end

main()
