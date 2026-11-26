# analysis/hazard.jl — flood hazard to people, and statistics over a run.
#
# Pure functions of water depth d (m) and velocity magnitude v (m/s), scalar
# (broadcast them over a frame's cells or a cell's series), plus a
# one-pass accumulator of per-cell statistics over frames. No dependencies:
# included by FloodA5Analysis (FA) and by the globe viewers' server
# (visualisation/VisualisationServer.jl), so the viewer, the static plots and
# the notebooks share one implementation. Descriptions: io/variables.toml.
#
# Methods
#   depth_velocity     D·V (m²/s)
#   hazard_rating      Defra/EA flood hazard rating HR = d (v + 0.5) + DF
#                      (FD2321/TR1 §3.5; debris factor Table 3.1; classes
#                      FD2320/TR2 Table 13.1 as in the EA/HR Wallingford
#                      supplementary note of May 2008)
#   arr_hazard_class   combined hazard curves H1–H6 (AIDR Guideline 7-3,
#                      2017, Table 2; Smith et al. 2014; ARR 2019 Book 6
#                      ch. 7)
#   hazard_degree      Chen, Xia, Falconer & Guo (2019), J. Flood Risk
#                      Management 12(3) e12486, eq. (10) and (14): toppling
#                      of a person, HD = min(1, v / U_c)
#   run_statistics     per-cell maxima, time of first flooding and duration
#                      of flooding over the frames of a run

# (exports: in FloodA5Analysis.jl — the server includes this file into its own
# module and must not export these names into Main)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

"""
Debris factor DF by land use, for the three depth/velocity rows of FD2321/TR1
Table 3.1 (also FD2321/TR2 Table 3.1 and the 2008 supplementary note):
d ≤ 0.25 m; 0.25 < d ≤ 0.75 m (v ≤ 2 m/s); d > 0.75 m and/or v > 2 m/s.

- `fd2320` (default): the simplified, precautionary DF of FD2320/TR2 Table
  13.1 for any land use — 0.5 for d ≤ 0.25 m, 1.0 above.
- `pasture`, `woodland`, `urban`: FD2321 Table 3.1 by dominant land use.
- `none`: DF = 0 (FD2321/TR1 Figure 3.2).

Reading of the table's "and/or v > 2": the third row applies to depths above
0.25 m with v > 2 m/s as well as to d > 0.75 m; depths ≤ 0.25 m stay in the
first row (all its factors are 0, except FD2320's 0.5).
"""
const DEBRIS_TABLES = Dict(
    "fd2320"   => (0.5, 1.0, 1.0),
    "pasture"  => (0.0, 0.0, 0.5),
    "woodland" => (0.0, 0.5, 1.0),
    "urban"    => (0.0, 1.0, 1.0),
    "none"     => (0.0, 0.0, 0.0),
)
const LAND_USES = ("fd2320", "pasture", "woodland", "urban", "none")
const LAND_USE_LABELS = Dict(
    "fd2320" => "any (FD2320 default: DF 0.5 / 1)", "pasture" => "pasture / arable",
    "woodland" => "woodland", "urban" => "urban", "none" => "no debris (DF 0)")

"Body mass (kg) and height (m) of the people in Chen et al. (2019) §4."
const PERSON_PRESETS = Dict("adult" => (mass = 75.1, height = 1.77),
                            "child" => (mass = 25.5, height = 1.26))

"""
    HazardParams(; landuse = "fd2320", mass = 75.1, height = 1.77, theta = 75.0,
                   alpha = 3.064, beta = 0.303, flood_depth = 0.1)

Parameters of the hazard variables and statistics:
- `landuse`: one land use for the whole domain (`LAND_USES`), for the Defra
  debris factor (a land-use raster is a later development);
- `mass` (kg), `height` (m): the person of Chen et al. (2019) (defaults: their
  adult); `theta` (°, body tilt), `alpha` (m^0.5/s), `beta`: their
  recommended calibration (θ = 75°, unified α and β of Table 2);
- `flood_depth` (m): a cell counts as flooded for the time of first flooding
  and the duration when its depth is at least this.
"""
Base.@kwdef struct HazardParams
    landuse::String = "fd2320"
    mass::Float64 = 75.1
    height::Float64 = 1.77
    theta::Float64 = 75.0
    alpha::Float64 = 3.064
    beta::Float64 = 0.303
    flood_depth::Float64 = 0.1
end

"A hashable key of `p` (caches)."
hazard_key(p::HazardParams) = (p.landuse, p.mass, p.height, p.theta, p.alpha, p.beta, p.flood_depth)

_num(q, k, dflt, lo, hi) = (x = tryparse(Float64, string(get(q, k, ""))); x === nothing || !isfinite(x) ? dflt : clamp(x, lo, hi))

"""
    hazard_params(q) → HazardParams

Parameters from a string dictionary (an HTTP query): `landuse`, `mass`,
`height`, `theta`, `alpha`, `beta`, `flood_depth`. Missing or invalid values
take the defaults; numbers are clamped to plausible ranges.
"""
function hazard_params(q::AbstractDict)
    d = HazardParams()
    lu = string(get(q, "landuse", d.landuse))
    return HazardParams(landuse = lu in LAND_USES ? lu : d.landuse,
                        mass = _num(q, "mass", d.mass, 1.0, 300.0),
                        height = _num(q, "height", d.height, 0.3, 2.5),
                        theta = _num(q, "theta", d.theta, 1.0, 90.0),
                        alpha = _num(q, "alpha", d.alpha, 0.01, 100.0),
                        beta = _num(q, "beta", d.beta, 0.0, 5.0),
                        flood_depth = _num(q, "flood_depth", d.flood_depth, 0.0, 100.0))
end

# ---------------------------------------------------------------------------
# Hazard variables (scalar; broadcast them)
# ---------------------------------------------------------------------------

"Depth × velocity (m²/s); 0 where dry."
depth_velocity(d::Real, v::Real) = d > 0 ? d * max(v, zero(v)) : zero(d * v)

"The debris factors of land use `landuse` (see `DEBRIS_TABLES`)."
debris_table(landuse::AbstractString) = DEBRIS_TABLES[landuse]

"""
    debris_factor(d, v, table) → DF

FD2321 Table 3.1 row for depth `d` and velocity `v`, from a `debris_table`
(0 where dry).
"""
function debris_factor(d::Real, v::Real, t::NTuple{3,Float64})
    d > 0 || return 0.0
    d <= 0.25 && return t[1]
    return (d > 0.75 || v > 2) ? t[3] : t[2]
end
debris_factor(d::Real, v::Real, landuse::AbstractString) = debris_factor(d, v, debris_table(landuse))

"""
    hazard_rating(d, v, table) → HR

Defra/EA flood hazard rating HR = d (v + 0.5) + DF (FD2321/TR1 §3.5), DF from
`debris_factor` (`table` = `debris_table(landuse)`, or a land-use name).
0 where dry. Classes (FD2320/TR2 Table 13.1): < 0.75 low (caution), 0.75–1.25
danger for some, 1.25–2.0 danger for most, > 2.0 danger for all (FD2321/TR1:
2.5 for the last).
"""
hazard_rating(d::Real, v::Real, t) = d > 0 ? d * (max(v, zero(v)) + 0.5) + debris_factor(d, v, t) : 0.0

"""
    arr_hazard_class(d, v) → 0–6

Combined flood hazard curves (AIDR Guideline 7-3, 2017, Table 2; ARR 2019):
the lowest class whose D·V, depth and velocity limits all hold —
H1 D·V ≤ 0.3, d ≤ 0.3, v ≤ 2; H2 0.6, 0.5, 2; H3 0.6, 1.2, 2; H4 1.0, 2.0, 2;
H5 4.0, 4.0, 4.0; H6 otherwise. 0 where dry.
"""
function arr_hazard_class(d::Real, v::Real)
    d > 0 || return 0
    v = max(v, zero(v)); dv = d * v
    (dv <= 0.3 && d <= 0.3 && v <= 2) && return 1
    (dv <= 0.6 && d <= 0.5 && v <= 2) && return 2
    (dv <= 0.6 && d <= 1.2 && v <= 2) && return 3
    (dv <= 1.0 && d <= 2.0 && v <= 2) && return 4
    (dv <= 4.0 && d <= 4.0 && v <= 4) && return 5
    return 6
end

"""
    chen_critical_velocity(h; mass = 75.1, height = 1.77, theta = 75, alpha = 3.064,
                           beta = 0.303, rho = 1000, a1 = 0.735, b1 = 0.265,
                           a2 = 1.015e-3, b2 = -4.927e-3, gamma = 10) → U_c (m/s)

Incipient velocity for toppling of a person (Chen et al. 2019, eq. 10) in
water of depth `h` (m):

    U_c = α (h/h_p)^β sqrt( m_p/(ρ h²) (γ cosθ + sinθ)
                            − [a1 (h/h_p)² + b1 (h/h_p)] (a2 m_p + b2) (γ cosθ/(h h_p) + sinθ/h²) )

with body mass m_p = `mass`, height h_p = `height`, tilt θ = `theta` (degrees).
Valid for h < h_p: at h ≥ h_p (floating) and where the radicand is negative,
U_c = 0; at h ≤ 0, U_c = Inf. θ = 90° gives eq. (11) (Xia et al. 2014).
"""
function chen_critical_velocity(h::Real; mass::Real = 75.1, height::Real = 1.77, theta::Real = 75.0,
                                alpha::Real = 3.064, beta::Real = 0.303, rho::Real = 1000.0,
                                a1::Real = 0.735, b1::Real = 0.265, a2::Real = 1.015e-3,
                                b2::Real = -4.927e-3, gamma::Real = 10.0)
    h > 0 || return Inf
    h >= height && return 0.0
    r = h / height
    s, c = sincosd(theta)
    rad = mass / (rho * h^2) * (gamma * c + s) -
          (a1 * r^2 + b1 * r) * (a2 * mass + b2) * (gamma * c / (h * height) + s / h^2)
    rad > 0 || return 0.0
    return alpha * r^beta * sqrt(rad)
end
chen_critical_velocity(h::Real, p::HazardParams) =
    chen_critical_velocity(h; mass = p.mass, height = p.height, theta = p.theta, alpha = p.alpha, beta = p.beta)

"""
    hazard_degree(d, v, p::HazardParams) → HD ∈ [0, 1]

Chen et al. (2019) eq. (14): HD = min(1, v / U_c) with U_c from
`chen_critical_velocity`; 1 where the person floats (d ≥ height) or U_c = 0;
0 where dry or still.
"""
function hazard_degree(d::Real, v::Real, p::HazardParams)
    d > 0 || return 0.0
    d >= p.height && return 1.0
    uc = chen_critical_velocity(d, p)
    uc > 0 || return 1.0
    return min(1.0, max(v, zero(v)) / uc)
end
hazard_degree(d::Real, v::Real; kw...) = hazard_degree(d, v, HazardParams(; kw...))

# ---------------------------------------------------------------------------
# Derived variables of a frame (or of a cell's series)
# ---------------------------------------------------------------------------

"""
The hazard variables derived from water depth and velocity, by name (the
names the viewer's server serves): name → f(d, v, p::HazardParams).
"""
const HAZARD_VARS = Dict{String,Function}(
    "depth_velocity"   => (d, v, p) -> depth_velocity(d, v),
    "hazard_rating"    => (d, v, p) -> hazard_rating(d, v, debris_table(p.landuse)),
    "hazard_degree"    => (d, v, p) -> hazard_degree(d, v, p),
    "hazard_class_arr" => (d, v, p) -> arr_hazard_class(d, v),
)

"""
    derive_hazard(var, depth, velocity, p = HazardParams()) → Vector{Float32}

Hazard variable `var` (a key of `HAZARD_VARS`) elementwise over `depth` and
`velocity` (equal lengths: a frame's cells, or one cell's frames). NaN where
either input is not finite.
"""
function derive_hazard(var::AbstractString, depth::AbstractVector, velocity::AbstractVector,
                       p::HazardParams = HazardParams())
    f = HAZARD_VARS[var]
    length(depth) == length(velocity) || throw(DimensionMismatch("depth and velocity lengths differ"))
    out = Vector{Float32}(undef, length(depth))
    t = debris_table(p.landuse)
    @inbounds for i in eachindex(depth, velocity, out)
        d = Float64(depth[i]); v = Float64(velocity[i])
        out[i] = (isfinite(d) && isfinite(v)) ? Float32(f(d, v, p)) : NaN32
    end
    return out
end

# ---------------------------------------------------------------------------
# Statistics over a run (per cell)
# ---------------------------------------------------------------------------

"""
The statistics `run_statistics` returns, by name, with the frame variable or
hazard each one summarises.
"""
const STAT_VARS = ("max_water_depth", "max_velocity", "max_depth_velocity", "max_hazard_rating",
                   "max_hazard_degree", "max_hazard_class_arr", "time_max_water_depth",
                   "time_first_flooded", "flood_duration")

"""
    run_statistics(times, frame, n, p = HazardParams()) → Dict{String,Vector{Float32}}

Per-cell statistics over frames at `times` (s), in one pass: `frame(k)`
returns `(depth, velocity)` of frame k (length-`n` vectors; `velocity` may be
`nothing`, and then the velocity and hazard statistics are left out).

- `max_water_depth`, `max_velocity`, `max_<hazard>` (`HAZARD_VARS`): maxima over
  the frames — at frame resolution, so peaks between output frames are missed;
- `time_max_water_depth` (s): the time of the (first) frame with the
  maximum depth (NaN where the depth never exceeds 0);
- `time_first_flooded` (s): the time of the first frame with depth ≥
  `p.flood_depth` (NaN if never);
- `flood_duration` (s): the time with depth ≥ `p.flood_depth`, each frame
  standing for the interval between the midpoints to its neighbours (the
  first and last frames for half an interval).

Non-finite values are skipped; a cell with none gets NaN maxima. The frames
go through a `StatsAccumulator`, which the model also uses during a run
(`--stats-interval`), so both give the same values from the same states.
"""
function run_statistics(times::AbstractVector, frame::Function, n::Integer, p::HazardParams = HazardParams())
    a = StatsAccumulator(n, p)
    for k in eachindex(times)
        d, v = frame(k)
        length(d) == n || throw(DimensionMismatch("frame $k has $(length(d)) values, expected $n"))
        accumulate_frame!(a, times[k], d, v)
    end
    return accumulated_statistics(a)
end

"""
    StatsAccumulator(n, p = HazardParams(); hazards = true)

The state of `run_statistics` between frames, for statistics gathered while
a run goes on (the model, `--stats-interval`): add each state with
`accumulate_frame!(a, t, depth, velocity)` (velocity may be `nothing`) and
read the statistics so far with `accumulated_statistics(a)`. `hazards =
false` leaves out the maxima that depend on the viewer's hazard parameters
(hazard rating: land use; hazard degree: the person) — D·V and the ARR
class depend on none. `p.flood_depth` sets the flooded threshold.

A frame's duration weight (half the interval between its neighbours) is
known only when the next frame comes; it is added then, in frame order, so
the arithmetic is that of a pass over all frames.
"""
mutable struct StatsAccumulator
    n       :: Int
    p       :: HazardParams
    table   :: NTuple{3,Float64}
    hazards :: Bool
    nf      :: Int                  # frames added
    t1      :: Float64              # time of the frame before the latest
    t2      :: Float64              # time of the latest frame
    pending :: BitVector            # cells flooded at the latest frame (weight not yet known)
    md      :: Vector{Float32}
    tmax    :: Vector{Float32}
    tfirst  :: Vector{Float32}
    dur     :: Vector{Float32}
    have_v  :: Bool
    mv      :: Vector{Float32}
    mdv     :: Vector{Float32}
    mar     :: Vector{Float32}
    mhr     :: Vector{Float32}
    mhd     :: Vector{Float32}
end

StatsAccumulator(n::Integer, p::HazardParams = HazardParams(); hazards::Bool = true) =
    StatsAccumulator(Int(n), p, debris_table(p.landuse), hazards, 0, NaN, NaN, falses(n),
                     fill(NaN32, n), fill(NaN32, n), fill(NaN32, n), zeros(Float32, n), false,
                     Float32[], Float32[], Float32[], Float32[], Float32[])

"Add the state at time `t` (depth, velocity or `nothing`; length-n vectors, any real type)."
function accumulate_frame!(a::StatsAccumulator, t::Real, d::AbstractVector, v = nothing)
    n = a.n
    length(d) == n || throw(DimensionMismatch("depth has $(length(d)) values, expected $n"))
    v === nothing || length(v) == n || throw(DimensionMismatch("velocity has $(length(v)) values, expected $n"))
    # the previous frame's duration weight is known now
    if a.nf >= 1
        w = (Float64(t) - (a.nf >= 2 ? a.t1 : a.t2)) / 2
        _add_pending!(a, w)
    end
    fd = a.p.flood_depth
    md = a.md; tmax = a.tmax; tfirst = a.tfirst; pend = a.pending
    @inbounds for i in 1:n
        x = Float64(d[i])
        if !isfinite(x)
            pend[i] = false
            continue
        end
        if isnan(md[i]) || x > md[i]
            md[i] = x
            x > 0 && (tmax[i] = t)
        end
        f = x >= fd
        pend[i] = f
        f && isnan(tfirst[i]) && (tfirst[i] = t)
    end
    if v !== nothing
        if !a.have_v
            a.have_v = true
            a.mv = fill(NaN32, n); a.mdv = fill(NaN32, n); a.mar = fill(NaN32, n)
            if a.hazards
                a.mhr = fill(NaN32, n); a.mhd = fill(NaN32, n)
            end
        end
        mv = a.mv; mdv = a.mdv; mar = a.mar; mhr = a.mhr; mhd = a.mhd
        tb = a.table; p = a.p; hz = a.hazards
        @inbounds for i in 1:n
            x = Float64(d[i]); u = Float64(v[i])
            (isfinite(x) && isfinite(u)) || continue
            _maxin!(mv, i, u)
            _maxin!(mdv, i, depth_velocity(x, u))
            if hz
                _maxin!(mhr, i, hazard_rating(x, u, tb))
                _maxin!(mhd, i, hazard_degree(x, u, p))
            end
            _maxin!(mar, i, arr_hazard_class(x, u))
        end
    end
    a.t1 = a.t2; a.t2 = Float64(t); a.nf += 1
    return a
end

function _add_pending!(a::StatsAccumulator, w::Float64)
    dur = a.dur; pend = a.pending
    @inbounds for i in 1:a.n
        pend[i] && (dur[i] += w)
    end
end

"""
    accumulated_statistics(a) → Dict{String,Vector{Float32}}

The statistics of the frames added so far (as `run_statistics`; copies — the
accumulator can go on). The latest frame counts for half the interval since
the one before it.
"""
function accumulated_statistics(a::StatsAccumulator)
    out = Dict{String,Vector{Float32}}("max_water_depth" => copy(a.md), "time_max_water_depth" => copy(a.tmax),
                                       "time_first_flooded" => copy(a.tfirst))
    dur = copy(a.dur)
    if a.nf >= 1
        w = a.nf >= 2 ? (a.t2 - a.t1) / 2 : 0.0
        @inbounds for i in 1:a.n
            a.pending[i] && (dur[i] += w)
        end
    end
    out["flood_duration"] = dur
    if a.have_v
        out["max_velocity"] = copy(a.mv); out["max_depth_velocity"] = copy(a.mdv)
        out["max_hazard_class_arr"] = copy(a.mar)
        if a.hazards
            out["max_hazard_rating"] = copy(a.mhr); out["max_hazard_degree"] = copy(a.mhd)
        end
    end
    return out
end

@inline function _maxin!(a::Vector{Float32}, i::Int, x::Real)
    y = Float32(x)
    (isnan(a[i]) || y > a[i]) && (a[i] = y)
    return nothing
end
