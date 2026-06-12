# boundaryinputs/sources.jl
# -------------------------
# Abstract source hierarchy and concrete source types for all water-adding
# boundary inputs.
#
# Design
# ------
# All volume additions to the domain are modelled as subtypes of AbstractSource.
# The simulation loop calls apply_source!(state, src, t, dt) for each source
# each timestep, then routes flux.  This pattern allows new source types
# (SpatialRainfall, Infiltration, etc.) to be added without changing the loop.
#
# Existing types InjectionPoint and RainPoint are promoted to AbstractSource
# subtypes.  Their structs are unchanged; only the supertype is added.
#
# Volume additions happen BEFORE the flux routing step so that injected water
# participates in that timestep's flux computation.  Multiple sources targeting
# the same cell are summed — this is by design (e.g. a .bci inflow on the same
# cell as a --injection-point).
#
# Mass balance
# ------------
# Each source implements cumulative_volume(src, t) → Float64 (m³) for the
# mass balance logger.  The logger sums across all sources.
#
# Future additions
# ----------------
# SpatialRainfall <: AbstractSource — gridded rainfall field applied per-cell
# UniformRainfallSource <: AbstractSource — wrapper for the scalar rainfall_rate
#   (currently handled as a special case in the loop; could be unified here
#    in a future cleanup pass)

# ---------------------------------------------------------------------------
# Abstract type
# ---------------------------------------------------------------------------

"""
    AbstractSource

Supertype for all water-adding boundary inputs.

Concrete subtypes must implement:
  apply_source!(state, source, t, dt)  — add dV to state.volume
  cumulative_volume(source, t)         — total volume injected up to time t (m³)
  source_label(source)                 — short string for logging
"""
abstract type AbstractSource end

# ---------------------------------------------------------------------------
# InflowPoint — time-varying hydrograph at a point
# ---------------------------------------------------------------------------

"""
    InflowPoint <: AbstractSource

Time-varying volumetric inflow at a single mesh cell, driven by a hydrograph.

Rate is linearly interpolated from (t_s, Q_m3s) at each timestep.
Flat extrapolation beyond the first and last knot (not clipped to zero —
a hydrograph record starting mid-storm carries its first value before t_offset).

`t_s` is seconds from simulation start (relative time).
For absolute-time hydrographs, the caller subtracts the epoch before
constructing the InflowPoint (see timeseries_io.jl `t_offset` kwarg).
"""
struct InflowPoint <: AbstractSource
    cell_index :: Int             # index into state.cell_ids
    cell_id    :: String          # hex cell ID (for logging)
    lon        :: Float64         # requested longitude (degrees)
    lat        :: Float64         # requested latitude (degrees)
    t_s        :: Vector{Float64} # knot times (seconds from sim start), monotonically ↑
    Q_m3s      :: Vector{Float64} # discharge at each knot (m³/s), same length as t_s
    label      :: String          # gauge ID / name for logging and mass balance output
end

source_label(src::InflowPoint) = src.label

"""
    _interp_hydrograph(t_s, Q_m3s, t) → Float64

Linear interpolation of a (time, discharge) hydrograph at simulation time `t`.
Flat extrapolation before the first knot and after the last knot.
Uses binary search (searchsortedlast) for O(log n) bracket lookup.
"""
@inline function _interp_hydrograph(t_s   :: Vector{Float64},
                                     Q_m3s :: Vector{Float64},
                                     t     :: Float64)::Float64
    isempty(t_s)    && return 0.0
    t <= t_s[1]     && return Q_m3s[1]
    t >= t_s[end]   && return Q_m3s[end]
    lo = searchsortedlast(t_s, t)
    frac = (t - t_s[lo]) / (t_s[lo+1] - t_s[lo])
    return Q_m3s[lo] + frac * (Q_m3s[lo+1] - Q_m3s[lo])
end

function apply_source!(state::FlowState, src::InflowPoint, t::Float64, dt::Float64)
    Q = _interp_hydrograph(src.t_s, src.Q_m3s, t)
    Q > 0.0 && (state.volume[src.cell_index] += Q * dt)
end

"""
    cumulative_volume(src::InflowPoint, t) → Float64

Cumulative volume injected from t=0 to simulation time `t`, computed by
trapezoidal integration over the stored hydrograph knots with a final
partial-interval interpolation.
"""
function cumulative_volume(src::InflowPoint, t::Float64)::Float64
    isempty(src.t_s) && return 0.0
    t <= 0.0         && return 0.0
    # Before first knot: flat extrapolation at Q[1]
    if t <= src.t_s[1]
        return max(0.0, src.Q_m3s[1]) * t
    end
    vol = 0.0
    # Pre-first-knot segment (t_s[1] > 0 case — rare but handled)
    if src.t_s[1] > 0.0
        vol += max(0.0, src.Q_m3s[1]) * src.t_s[1]
    end
    for k in 1:length(src.t_s)-1
        t1, t2 = src.t_s[k], src.t_s[k+1]
        Q1 = max(0.0, src.Q_m3s[k])
        Q2 = max(0.0, src.Q_m3s[k+1])
        if t2 <= t
            vol += 0.5 * (Q1 + Q2) * (t2 - t1)
        else
            frac = (t - t1) / (t2 - t1)
            Q_t  = Q1 + frac * (Q2 - Q1)
            vol += 0.5 * (Q1 + Q_t) * (t - t1)
            return vol
        end
    end
    # After last knot: flat extrapolation at Q[end]
    vol += max(0.0, src.Q_m3s[end]) * (t - src.t_s[end])
    return vol
end

# ---------------------------------------------------------------------------
# AbstractSource dispatch for existing types (InjectionPoint, RainPoint)
# ---------------------------------------------------------------------------
# These types are defined in FloodModel.jl (before this file is included).
# We add the dispatch methods here so that the unified all_sources loop works.

function apply_source!(state::FlowState, src::InjectionPoint, t::Float64, dt::Float64)
    state.volume[src.cell_index] += src.rate_m3s * dt
end

function apply_source!(state::FlowState, src::RainPoint, t::Float64, dt::Float64)
    state.volume[src.cell_index] += src.rate_m3s * dt
end

cumulative_volume(src::InjectionPoint, t::Float64) = src.rate_m3s * t
cumulative_volume(src::RainPoint,      t::Float64) = src.rate_m3s * t

source_label(src::InjectionPoint) = "injection @ $(src.cell_id)"
source_label(src::RainPoint)      = "rainpoint @ $(src.cell_id)"

# ---------------------------------------------------------------------------
# Aggregate helpers used by the simulation loop
# ---------------------------------------------------------------------------

"""
    total_cumulative_input(sources, rainfall_rate, cell_areas, t) → Float64

Total volume input to the domain up to simulation time `t` from all sources
plus uniform rainfall.  Used by the mass balance logger.
"""
function total_cumulative_input(sources       :: Vector{<:AbstractSource},
                                 rainfall_rate  :: Float64,
                                 cell_areas     :: Vector{Float64},
                                 t              :: Float64)::Float64
    vol = 0.0
    for src in sources
        vol += cumulative_volume(src, t)
    end
    if rainfall_rate > 0.0
        vol += rainfall_rate * t * sum(a for a in cell_areas if a >= 1.0; init=0.0)
    end
    return vol
end
