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
#
#   --frame N is the 1-based frame index (the same number as the old
#   /frames/NNNNNN group name); output-format 1 and 2 files both work
#   (io/FloodA5Output.jl).
#
# The metric lives in analysis/shape.jl (FloodA5Analysis: `FA.ns_split`);
# this script is a CLI over it. Its convention, passed explicitly: the
# north/south offset of a cell is Δlat · 111_320 m/° (FA.M_PER_DEG_LAT) from
# the nominal --source-lat (wp1_metrics.jl planar uses the equirectangular
# R = 6_371_008.8 m about the source cell's centre instead); non-finite
# volumes are summed as they are.

using Printf
include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
const FA = FloodA5Analysis

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
# Core analysis
# ---------------------------------------------------------------------------

function analyse_symmetry(
    r::FA.Run,
    frame_idx::Int;
    source_lat::Float64,
    lat_tolerance_m::Float64=50.0
)
    FA.check_frame(r, frame_idx; hint = " Run with --list-frames to inspect.")
    _, dlat_m = FA.shape_xy(r.mesh, 0.0, source_lat; snap = false,
                            m_per_deg = FA.M_PER_DEG_LAT)
    s = FA.ns_split(FA.volume(r, frame_idx), dlat_m; band = lat_tolerance_m)
    return (
        t=r.t[frame_idx],
        n_north=s.n_north,
        n_south=s.n_south,
        n_excluded=s.n_band,
        V_north=s.V_north,
        V_south=s.V_south,
        V_symline=s.V_band,
        V_ns=s.V_ns,
        asymmetry=s.asym
    )
end
analyse_symmetry(h5_path::String, frame_idx::Int; kw...) =
    analyse_symmetry(FA.load_run(h5_path), frame_idx; kw...)

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
            path !== nothing && isfile(path) && FA.print_frames(FA.load_run(path))
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
    rb, rc = FA.load_run(baseline_path), FA.load_run(corrected_path)
    n_frames = FA.nframes(rb)

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
                rb,
                fidx;
                source_lat=source_lat,
                lat_tolerance_m=lat_tolerance_m
            )

            b = analyse_symmetry(
                rc,
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
