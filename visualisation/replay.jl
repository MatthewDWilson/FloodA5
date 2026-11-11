# visualisation/replay.jl — replay an output file in the GLMakie viewer
#
#   julia --project=. --threads auto visualisation/replay.jl RUN.h5 [options]
#   julia --project=. --threads auto FloodModel.jl --replay RUN.h5 [options]
#
# Opens the live viewer's window (`--vis makie`: map, colorbar, controls,
# diagnostics, time-series strip) on a finished — or still running — run's
# output file, whether or not it was visualised live. Adds a frame slider and
# a playback-speed slider; Pause/Play drives playback. Format-2 files
# (FloodA5 0.4.0+) carry the mesh; format-1 files need it (--mesh).
#
# Options:
#   --mesh FILE.parquet  the run's GeoParquet mesh (format-1 files only)
#   --fps N              playback speed, output frames per second (default 4)
#   --var V              depth | saturation | volume | velocity | ring
#                        (initial display variable; default depth)
#   --from T             start at the frame nearest simulation time T (s)
#   --paused             start paused
#   --png FILE           render the window at --from (default: last frame)
#                        to FILE and exit, without opening a window
#
# What the strip shows: domain volume and wet cells are read from the file
# for every frame; input volume is rebuilt from the run record and checked
# against the recorded total; removed volume is stored only at the end of
# the run (see visualisation/replay_data.jl).

isdefined(Main, :FloodA5Output) || include(joinpath(@__DIR__, "..", "io", "FloodA5Output.jl"))
isdefined(Main, :FloodA5Analysis) || include(joinpath(@__DIR__, "..", "analysis", "FloodA5Analysis.jl"))
isdefined(Main, :FloodA5Replay) || include(joinpath(@__DIR__, "replay_data.jl"))
isdefined(Main, :MakieVisualiser) || include(joinpath(@__DIR__, "MakieVisualiser.jl"))

const _REPLAY_VARS = Dict("depth" => 1, "saturation" => 2, "volume" => 3, "velocity" => 4,
                          "ring" => "Ring index")

function _replay_usage(code::Int)
    println("""
Usage: julia --project=. --threads auto visualisation/replay.jl RUN.h5 [options]

  --mesh FILE.parquet  the run's GeoParquet mesh (format-1 files only)
  --fps N              playback speed, output frames per second (default 4)
  --var V              depth | saturation | volume | velocity | ring (default depth)
  --from T             start at the frame nearest simulation time T (s)
  --paused             start paused
  --png FILE           save the window at --from (default: last frame) to FILE, no window
""")
    code == 0 || exit(code)
    return nothing
end

"""
    replay_main(args) → MakieVis or nothing

The command line of visualisation/replay.jl (also `FloodModel.jl --replay`).
"""
function replay_main(args::Vector{String})
    ("--help" in args || "-h" in args) && return _replay_usage(0)
    isdefined(MakieVisualiser, :replay) ||
        error("the GLMakie viewer is not loaded (FLOODA5_NO_VIS=1?): unset it to replay")
    args = copy(args)
    opt(flag) = (i = findfirst(==(flag), args);
                 i === nothing ? nothing :
                 (i == length(args) ? error("$flag needs a value") :
                  (v = args[i+1]; deleteat!(args, [i, i+1]); v)))
    flag(f) = (i = findfirst(==(f), args); i === nothing ? false : (deleteat!(args, i); true))
    mesh  = opt("--mesh")
    fps   = something(tryparse(Float64, something(opt("--fps"), "4")), 4.0)
    var   = something(opt("--var"), "depth")
    from  = opt("--from")
    png   = opt("--png")
    paused = flag("--paused")
    haskey(_REPLAY_VARS, var) ||
        error("--var $var: one of $(join(sort(collect(keys(_REPLAY_VARS))), ", "))")
    files = filter(a -> !startswith(a, "--"), args)
    length(files) == 1 || (println("ERROR: give one output file"); return _replay_usage(1))
    length(args) == 1 || error("unknown option(s): $(join(filter(startswith("--"), args), " "))")
    h5 = files[1]

    rd = FloodA5Replay.open_replay(h5; parquet = mesh)
    @info "Replay: $(rd.label): $(length(rd.t)) frames, t = $(rd.t[1])…$(rd.t[end]) s, " *
          "$(length(rd.mesh)) cells"
    k0 = from === nothing ? (png === nothing ? 1 : length(rd.t)) :
         FloodA5Analysis.frame_at(rd.run, parse(Float64, from))
    try
        if png !== nothing
            vis = MakieVisualiser.replay(rd; fps, var = _REPLAY_VARS[var], from = k0,
                                         show_window = false)
            MakieVisualiser.GLMakie.save(png, vis.figure)
            @info "Replay: frame $k0 (t = $(rd.t[k0]) s) saved to $png"
            return vis
        end
        return MakieVisualiser.replay(rd; fps, var = _REPLAY_VARS[var], from = k0,
                                      playing = !paused)
    finally
        close(rd)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    replay_main(copy(ARGS))
end
