#!/usr/bin/env julia
# test_planar_symmetry.jl
# ------------------------
# Planar slope north/south VOLUME SYMMETRY benchmark.
#
# This script does not run the FloodA5 model. It analyses existing HDF5
# simulation outputs.
#
# In addition to the normal human-readable output, --csv FILE can be used
# to write machine-readable results for batch benchmarking.
#
# Usage:
#   julia --project=. test/test_planar_symmetry.jl `
#       --baseline FILE.h5 `
#       --corrected FILE.h5 `
#       --source-lat 51.0001 `
#       --source-lon -0.0434 `
#       --frame 241 `
#       --csv result.csv
#
#   julia --project=. test/test_planar_symmetry.jl `
#       --baseline FILE.h5 `
#       --list-frames

using HDF5
using Statistics
using Printf

const M_PER_DEG_LAT = 111_320.0

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------

function _argval(args::Vector{String}, flag::String, default=nothing)
    idx = findfirst(==(flag), args)
    idx === nothing && return default
    idx == length(args) && error("$flag requires a value")
    return args[idx + 1]
end

const ARGS_V = String.(ARGS)
const DO_LIST = "--list-frames" in ARGS_V

# ---------------------------------------------------------------------------
# HDF5 access
# ---------------------------------------------------------------------------

function _frame_names(f)
    return sort(collect(keys(f["frames"])))
end

function _list_frames(h5_path::String)
    h5open(h5_path, "r") do f
        names = _frame_names(f)
        println("Frames in $h5_path:")
        @printf("%8s  %10s  %8s  %10s\n",
                "frame#", "t (s)", "n_wet", "max_depth")

        for nm in names
            grp = f["frames/$nm"]
            t = haskey(grp, "t") ? read(grp["t"]) : NaN
            sat = read(grp["saturation"])
            n_wet = count(>=(1.0), sat)
            depth = read(grp["water_depth"])
            md = isempty(depth) ? NaN : maximum(depth)

            @printf("%8s  %10.1f  %8d  %10.4f\n",
                    nm, t, n_wet, md)
        end
    end
end

"""
    _read_frame_volume(h5_path, frame_idx) → (center_lats, volume, t)
"""
function _read_frame_volume(h5_path::String, frame_idx::Int)
    h5open(h5_path, "r") do f
        center_lats = read(f["mesh/center_lats"])
        names = _frame_names(f)

        fname = lpad(string(frame_idx), 6, '0')

        haskey(f["frames"], fname) ||
            error("Frame $fname not found in $h5_path. Available: " *
                  "$(first(names, 3))…$(last(names, 3)) ($(length(names)) total). " *
                  "Run with --list-frames to inspect.")

        grp = f["frames/$fname"]

        volume = read(grp["volume"])
        t = haskey(grp, "t") ? read(grp["t"]) : NaN

        return center_lats, volume, t
    end
end

# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

function analyse_symmetry(
    h5_path::String,
    frame_idx::Int;
    source_lat::Float64,
    lat_tolerance_m::Float64=50.0
)
    center_lats, volume, t = _read_frame_volume(h5_path, frame_idx)

    dlat_m = (center_lats .- source_lat) .* M_PER_DEG_LAT

    north_mask   = dlat_m .>  lat_tolerance_m
    south_mask   = dlat_m .< -lat_tolerance_m
    symline_mask = .!north_mask .& .!south_mask

    n_excluded = count(symline_mask)

    V_north   = sum(volume[north_mask])
    V_south   = sum(volume[south_mask])
    V_symline = sum(volume[symline_mask])

    V_ns = V_north + V_south

    asymmetry =
        V_ns > 1e-9 ?
        (V_north - V_south) / V_ns :
        0.0

    return (
        t=t,
        n_north=count(north_mask),
        n_south=count(south_mask),
        n_excluded=n_excluded,
        V_north=V_north,
        V_south=V_south,
        V_symline=V_symline,
        V_ns=V_ns,
        asymmetry=asymmetry
    )
end

# ---------------------------------------------------------------------------
# Human-readable output
# ---------------------------------------------------------------------------

function _print_result(label::String, r)
    println(label, ":")

    @printf("  t = %.1f s\n", r.t)

    @printf(
        "  cells:  north=%d  south=%d  excluded (within tolerance band)=%d\n",
        r.n_north,
        r.n_south,
        r.n_excluded
    )

    pct_n = r.V_ns > 1e-9 ? 100 * r.V_north / r.V_ns : NaN
    pct_s = r.V_ns > 1e-9 ? 100 * r.V_south / r.V_ns : NaN

    @printf(
        "  volume: north=%.1f m³ (%.1f%%)   south=%.1f m³ (%.1f%%)   north+south=%.1f m³\n",
        r.V_north,
        pct_n,
        r.V_south,
        pct_s,
        r.V_ns
    )

    @printf(
        "  (symline band volume, excluded from the ratio above: %.1f m³)\n",
        r.V_symline
    )

    @printf(
        "  asymmetry = %+.4f   (0 = perfectly symmetric, i.e. 50%%/50%%; sign indicates which side is favoured)\n",
        r.asymmetry
    )

    println()
end

# ---------------------------------------------------------------------------
# CSV support
# ---------------------------------------------------------------------------

function _csv_escape(value)
    s = string(value)

    if occursin('"', s) || occursin(',', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end

    return s
end

function _csv_row(values)
    println(join(_csv_escape.(values), ","))
end

function _write_csv(
    csv_path::String,
    frame_idx::Int,
    source_lat::Float64,
    lat_tolerance_m::Float64,
    baseline_result,
    corrected_result
)
    open(csv_path, "w") do io
        redirect_stdout(io) do

            _csv_row([
                "case",
                "frame",
                "t",
                "source_lat",
                "lat_tolerance_m",
                "n_north",
                "n_south",
                "n_excluded",
                "V_north",
                "V_south",
                "V_symline",
                "V_ns",
                "pct_north",
                "pct_south",
                "asymmetry"
            ])

            function write_result(case_name, r)
                pct_n =
                    r.V_ns > 1e-9 ?
                    100 * r.V_north / r.V_ns :
                    NaN

                pct_s =
                    r.V_ns > 1e-9 ?
                    100 * r.V_south / r.V_ns :
                    NaN

                _csv_row([
                    case_name,
                    frame_idx,
                    r.t,
                    source_lat,
                    lat_tolerance_m,
                    r.n_north,
                    r.n_south,
                    r.n_excluded,
                    r.V_north,
                    r.V_south,
                    r.V_symline,
                    r.V_ns,
                    pct_n,
                    pct_s,
                    r.asymmetry
                ])
            end

            baseline_result !== nothing &&
                write_result("baseline", baseline_result)

            corrected_result !== nothing &&
                write_result("corrected", corrected_result)
        end
    end
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()

    if DO_LIST
        for flag in ("--baseline", "--corrected")
            path = _argval(ARGS_V, flag)
            path !== nothing && isfile(path) && _list_frames(path)
        end
        return
    end

    baseline_path  = _argval(ARGS_V, "--baseline")
    corrected_path = _argval(ARGS_V, "--corrected")
    frame_str      = _argval(ARGS_V, "--frame")
    sweep_str      = _argval(ARGS_V, "--sweep")
    source_lat_str = _argval(ARGS_V, "--source-lat")
    tol_str        = _argval(ARGS_V, "--lat-tolerance-m", "50.0")
    csv_path       = _argval(ARGS_V, "--csv")

    if baseline_path === nothing && corrected_path === nothing
        println("""
        Usage:
          julia test_planar_symmetry.jl \\
              --baseline FILE.h5 \\
              --corrected FILE.h5 \\
              --source-lat LAT \\
              --frame N \\
              [--lat-tolerance-m T] \\
              [--csv FILE]

          julia test_planar_symmetry.jl \\
              --baseline FILE.h5 \\
              --list-frames
        """)
        return
    end

    source_lat_str === nothing &&
        error("--source-lat LAT is required")

    source_lat = parse(Float64, source_lat_str)
    lat_tolerance_m = parse(Float64, tol_str)

    if sweep_str !== nothing
        (baseline_path === nothing || corrected_path === nothing) &&
            error("--sweep requires both --baseline and --corrected")

        run_sweep(
            baseline_path,
            corrected_path,
            parse(Int, sweep_str);
            source_lat=source_lat,
            lat_tolerance_m=lat_tolerance_m
        )

        return
    end

    frame_str === nothing &&
        error("--frame N is required")

    frame_idx = parse(Int, frame_str)

    baseline_result = nothing
    corrected_result = nothing

    if baseline_path !== nothing
        baseline_result = analyse_symmetry(
            baseline_path,
            frame_idx;
            source_lat=source_lat,
            lat_tolerance_m=lat_tolerance_m
        )

        _print_result(
            "Baseline (--gradient-correction off)",
            baseline_result
        )
    end

    if corrected_path !== nothing
        corrected_result = analyse_symmetry(
            corrected_path,
            frame_idx;
            source_lat=source_lat,
            lat_tolerance_m=lat_tolerance_m
        )

        _print_result(
            "Corrected (--gradient-correction on)",
            corrected_result
        )
    end

    if baseline_result !== nothing && corrected_result !== nothing

        a = baseline_result
        b = corrected_result

        println("─" ^ 60)
        println("Comparison (corrected vs baseline):")

        @printf(
            "  |asymmetry|:  %.4f → %.4f   (%s)\n",
            abs(a.asymmetry),
            abs(b.asymmetry),
            abs(b.asymmetry) < abs(a.asymmetry) ?
                "IMPROVED ✓" :
                "NOT improved ✗"
        )

        println()

        println(
            "Acceptance criterion (adapted for §10 Stage 3): the corrected run's " *
            "|asymmetry| should be markedly smaller than baseline's, since a pure " *
            "east-west slope has no physical north/south preference — any imbalance " *
            "is a numerical artifact of mesh non-orthogonality. No mass-balance/" *
            "stability regression (check mb_err separately from the model's own run logs)."
        )
    end

    # -----------------------------------------------------------------------
    # Machine-readable output
    # -----------------------------------------------------------------------

    if csv_path !== nothing
        _write_csv(
            csv_path,
            frame_idx,
            source_lat,
            lat_tolerance_m,
            baseline_result,
            corrected_result
        )

        println()
        println("CSV written to: $csv_path")
    end
end

# ---------------------------------------------------------------------------
# Sweep mode
# ---------------------------------------------------------------------------

function run_sweep(
    baseline_path::String,
    corrected_path::String,
    step::Int;
    source_lat::Float64,
    lat_tolerance_m::Float64
)
    n_frames = h5open(baseline_path, "r") do f
        length(_frame_names(f))
    end

    println("Sweep — every $step frames, $n_frames total frames available")

    @printf(
        "%6s  %8s | %10s %10s %9s | %10s %10s %9s | %9s\n",
        "frame",
        "t(s)",
        "base Vn",
        "base Vs",
        "base asym",
        "corr Vn",
        "corr Vs",
        "corr asym",
        "|asym| gain"
    )

    for fidx in 1:step:n_frames
        local a, b

        try
            a = analyse_symmetry(
                baseline_path,
                fidx;
                source_lat=source_lat,
                lat_tolerance_m=lat_tolerance_m
            )

            b = analyse_symmetry(
                corrected_path,
                fidx;
                source_lat=source_lat,
                lat_tolerance_m=lat_tolerance_m
            )
        catch e
            @printf(
                "%6d  (skipped — %s)\n",
                fidx,
                sprint(showerror, e)
            )
            continue
        end

        gain = abs(a.asymmetry) - abs(b.asymmetry)

        @printf(
            "%6d  %8.0f | %10.1f %10.1f %9.4f | %10.1f %10.1f %9.4f | %+9.4f\n",
            fidx,
            a.t,
            a.V_north,
            a.V_south,
            a.asymmetry,
            b.V_north,
            b.V_south,
            b.asymmetry,
            gain
        )
    end

    println()

    println(
        "'|asym| gain' = |baseline asymmetry| − |corrected asymmetry|; positive means"
    )

    println(
        "the corrected run is closer to the physically-expected symmetric ideal at that"
    )

    println(
        "frame."
    )
end

main()
