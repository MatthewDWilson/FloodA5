# surfacewater/limiter_stats.jl
# -----------------------------
# D1b (damping work, 2026-10): optional instrumentation of the limiters and
# switches in the edge kernels — `--limiter-stats`. Results are unchanged
# (the kernels only record flags); with stats on, Phase A always takes the
# full pass so every edge is recorded.
#
# Per edge and step (bit flags):
#   WET     the edge carries flow (passed the dry-edge exit)
#   FROUDE  the Froude cap bound            (standard kernels)
#   FIXB    the per-edge volume limiter bound (standard kernels; q ≤ h_d·w/(N·Δt))
#   DONOR   the per-edge donor cap reduced the edge volume (V/(2N) per step)
#   POS     positivity scaling reduced it    (--limiter positivity|hybrid)
#   WD      diamond-wetdry fell back to the two-point gradient
# A clamp returns its bound exactly, so binding is detected by comparing the
# kernel's output with the bound (an exact test for FROUDE with component
# friction and for FIXB; 1e-12 relative for the vector-friction Froude cap).
# Aggregated per step by flow-depth band (h_flow < 0.3, 0.3–1, 1–3, ≥ 3 m): the
# fraction of wet edge-steps on which each flag is set, and "toggles" — wet
# edges whose flags changed since the previous step, a direct measure of how
# much discrete switching the dynamics undergo.

const LIMITER_STATS = Ref(false)

# D1b experiment (`--boundary-limiter positivity`): boundary (ghost-edge)
# outflow without the per-step limits — no volume limiter (Fix B) in the ghost
# kernel and no V/(2N) donor cap, only "no more than the cell holds". Both of
# those limits scale with 1/Δt. Off (legacy limits) by default.
const BOUNDARY_POSITIVITY = Ref(false)

const LF_WET    = 0x20
const LF_FROUDE = 0x01
const LF_FIXB   = 0x02
const LF_DONOR  = 0x04
const LF_POS    = 0x08
const LF_WD     = 0x10
const _LF_NAMES = ((LF_FROUDE, "Froude cap"), (LF_FIXB, "volume limiter"), (LF_DONOR, "donor cap"),
                   (LF_POS, "positivity"), (LF_WD, "wet/dry fallback"))
const _LF_BANDS = ("< 0.3 m", "0.3–1 m", "1–3 m", "≥ 3 m")

mutable struct LimiterStats
    flag   :: Vector{UInt8}       # this step, per edge (written by the kernels)
    prev   :: Vector{UInt8}       # previous step
    hf     :: Vector{Float64}     # flow depth per edge this step
    steps  :: Int
    wet    :: Vector{Float64}     # wet edge-steps per band
    count  :: Matrix{Float64}     # (flag, band) edge-steps
    toggle :: Vector{Float64}     # toggles per band
    # hourly window (all bands)
    w_t0   :: Float64
    w_steps :: Int
    w_wet  :: Float64
    w_count :: Vector{Float64}
    w_toggle :: Float64
end
LimiterStats() = LimiterStats(UInt8[], UInt8[], Float64[], 0, zeros(4), zeros(5, 4), zeros(4),
                              0.0, 0, 0.0, zeros(5), 0.0)
const _LSTATS = Ref(LimiterStats())

_lf_band(h) = h < 0.3 ? 1 : h < 1.0 ? 2 : h < 3.0 ? 3 : 4

function _lstats_reset!(ne::Int)
    s = LimiterStats()
    s.flag = zeros(UInt8, ne); s.prev = zeros(UInt8, ne); s.hf = zeros(ne)
    _LSTATS[] = s
    return s
end

# arrays for the kernels' NamedTuples (empty when stats are off)
_lf_flag() = LIMITER_STATS[] ? _LSTATS[].flag : UInt8[]
_lf_hf()   = LIMITER_STATS[] ? _LSTATS[].hf : Float64[]

"""
    _lstats_step!(t)

Fold this step's per-edge flags into the totals (serial; stats mode only),
and print an hourly line.
"""
function _lstats_step!(t::Float64)
    s = _LSTATS[]
    ne = length(s.flag)
    s.steps += 1; s.w_steps += 1
    @inbounds for e in 1:ne
        f = s.flag[e]; p = s.prev[e]
        (f & LF_WET) == 0 && (p & LF_WET) == 0 && continue
        b = _lf_band(s.hf[e])
        if (f & LF_WET) != 0
            s.wet[b] += 1; s.w_wet += 1
            for (k, (bit, _)) in enumerate(_LF_NAMES)
                if (f & bit) != 0
                    s.count[k, b] += 1; s.w_count[k] += 1
                end
            end
        end
        if (f & 0x3f) != (p & 0x3f)
            s.toggle[b] += 1; s.w_toggle += 1
        end
    end
    copyto!(s.prev, s.flag)
    if t - s.w_t0 >= 3600.0
        w = max(s.w_wet, 1.0)
        @info "limiter stats t=$(round(t, digits = 0)) s: wet edges/step $(round(s.w_wet / s.w_steps, digits = 0)); " *
              "binding per wet edge-step: " *
              join(["$(n) $(round(100 * s.w_count[k] / w, sigdigits = 3))%" for (k, (_, n)) in enumerate(_LF_NAMES)], ", ") *
              "; toggles $(round(100 * s.w_toggle / w, sigdigits = 3))%"
        s.w_t0 = t; s.w_steps = 0; s.w_wet = 0.0; fill!(s.w_count, 0.0); s.w_toggle = 0.0
    end
    return nothing
end

"""
    _lstats_report()

End-of-run table: per depth band, the share of wet edge-steps on which each
limiter/switch was active, and the toggle rate.
"""
function _lstats_report()
    s = _LSTATS[]
    io = IOBuffer()
    println(io, "Limiter / switch activity over $(s.steps) steps (% of wet edge-steps per flow-depth band)")
    println(io, rpad("", 18), join(lpad.(_LF_BANDS, 10)), lpad("all", 10))
    tot = sum(s.wet)
    println(io, rpad("wet edges/step", 18), join(lpad.(string.(round.(Int, s.wet ./ max(s.steps, 1))), 10)),
            lpad(string(round(Int, tot / max(s.steps, 1))), 10))
    for (k, (_, n)) in enumerate(_LF_NAMES)
        row = [s.wet[b] > 0 ? 100 * s.count[k, b] / s.wet[b] : 0.0 for b in 1:4]
        all_ = tot > 0 ? 100 * sum(s.count[k, :]) / tot : 0.0
        println(io, rpad(n, 18), join(lpad.(string.(round.(row, sigdigits = 3)), 10)),
                lpad(string(round(all_, sigdigits = 3)), 10))
    end
    row = [s.wet[b] > 0 ? 100 * s.toggle[b] / s.wet[b] : 0.0 for b in 1:4]
    println(io, rpad("toggles", 18), join(lpad.(string.(round.(row, sigdigits = 3)), 10)),
            lpad(string(round(tot > 0 ? 100 * sum(s.toggle) / tot : 0.0, sigdigits = 3)), 10))
    @info String(take!(io))
end
