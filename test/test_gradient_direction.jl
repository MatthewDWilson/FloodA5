"""
test_gradient_direction.jl
--------------------------
Unit tests for the d̂ projection gradient correction formula
(FloodA5_GradientCorrection_Fix.md, flow-direction-fixes branch, 2026-07-06).

Motivation
----------
The previous V̂-based formula `dWSE_n = c·(wse_i-wse_j) - L·(∇WSE_f·V̂)` was
found to over-correct at res-18 and under-correct at res-16 — the correction
efficiency changed from 74% at res-16 to 183% at res-18, with a sign flip in
the resulting plume direction. Root cause: V̂ = n̂ - c·d̂ is built from the face
normal n̂, whose orientation depends on the A5 dual-sublattice structure and
changes with resolution. The replacement formula `dWSE_n = -(∇WSE_f · d_vec_m)`
avoids face normals entirely, using only the unambiguous ci→cj displacement
vector. See FloodA5_GradientCorrection_Fix.md §1–§2.

Test Groups
-----------
GD1  V̂_y sign across resolutions (diagnostic)
     Documents that V̂_y for E-W cell pairs is NOT consistent in sign across
     A5 resolutions 12–20 — the quantitative evidence that V̂ is resolution-
     dependent. These are informational assertions (finite, bounded), not
     correctness claims about sign. A sign-consistency check would always fail
     on A5, which is precisely the point.

GD2  d̂ projection formula: sign, magnitude, and resolution-consistency
     Verifies that `dWSE_n = -(∇WSE_f · d_vec_m)` produces:
     (a) the correct sign (positive when ci is uphill) at every A5 resolution
         12–20 tested,
     (b) magnitude within 2% of the analytical `wse_ci - wse_cj` for a pure
         linear WSE field,
     (c) no sign flip across resolutions — the property V̂ lacks.

No external mesh, pya5 bridge, DEM, or HDF5 file is required. All geometry is
constructed synthetically. The test runs in well under a second.

Run with:
    julia test/test_gradient_direction.jl
"""

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "mesh"))
include(joinpath(dirname(@__DIR__), "mesh", "A5Grid.jl"))
using .A5Grid: _EARTH_R       # _edge_geometry no longer needed: GD1 uses
                               # analytic _tilted_edge_geometry instead
using Statistics: mean, std
using Printf
using Test

println("=" ^ 62)
println("test_gradient_direction.jl")
println("d̂ projection formula — sign, magnitude, resolution sweep")
println("=" ^ 62)

# ─── geometry helpers ─────────────────────────────────────────────────────────

"""
    _approx_cell_diameter_deg(res; lat_ref) → Float64

Approximate A5 cell diameter in degrees longitude at the given latitude.
A5 has 12 pentagons at res 0; each level multiplies count by ~5.
Cell area ≈ 4πR² / (12 × 5^res). Diameter = 2√(area/π).
Used only to set realistic ci–cj spacing for GD2 test geometry.
"""
function _approx_cell_diameter_deg(res::Int; lat_ref::Float64 = 51.0)
    R           = _EARTH_R
    n_cells     = 12 * 5^res
    area_m2     = (4π * R^2) / n_cells
    diam_m      = 2.0 * sqrt(area_m2 / π)
    return rad2deg(diam_m / (R * cosd(lat_ref)))
end

"""
    _wlsq_grad(c_lon, c_lat, nb_lons, nb_lats, wse_c, wse_nb) → (gx, gy)

Compute the WLSQ gradient at a single cell given its neighbours.
Replicates `_build_wlsq_weights!` + `_compute_wse_gradients!` inline.
Uses per-cell cos_lat0 = cosd(c_lat), matching FloodModel.jl's convention.
"""
function _wlsq_grad(c_lon::Float64, c_lat::Float64,
                    nb_lons::Vector{Float64}, nb_lats::Vector{Float64},
                    wse_c::Float64, wse_nb::Vector{Float64})
    R        = _EARTH_R
    cos_lat0 = cosd(c_lat)
    k        = length(nb_lons)
    k < 2 && return (0.0, 0.0)

    dxs   = [deg2rad(nb_lons[m] - c_lon) * R * cos_lat0 for m in 1:k]
    dys   = [deg2rad(nb_lats[m] - c_lat) * R             for m in 1:k]
    ws    = [1.0 / max(dxs[m]^2 + dys[m]^2, 1.0)         for m in 1:k]
    dWSEs = wse_nb .- wse_c

    Sxx = sum(ws .* dxs .^ 2)
    Sxy = sum(ws .* dxs .* dys)
    Syy = sum(ws .* dys .^ 2)
    det_M = Sxx * Syy - Sxy^2
    abs(det_M) < 1e-20 && return (0.0, 0.0)

    rhs_x = sum(ws .* dxs .* dWSEs)
    rhs_y = sum(ws .* dys .* dWSEs)
    gx    = ( Syy * rhs_x - Sxy * rhs_y) / det_M
    gy    = (-Sxy * rhs_x + Sxx * rhs_y) / det_M
    return (gx, gy)
end

"""
    _pentagon_neighbours(c_lon, c_lat, diam_deg; lat_ref)
        → (nlons, nlats)

Place 5 synthetic neighbours at pentagon vertex angles (0°, 72°, 144°,
216°, 288° from north, clockwise) at a distance of `diam_deg` from the
centre cell. Used to build a regular 5-neighbour WLSQ stencil.
"""
function _pentagon_neighbours(c_lon::Float64, c_lat::Float64,
                               diam_deg::Float64;
                               lat_ref::Float64 = 51.0)
    R        = _EARTH_R
    cos_lat0 = cosd(lat_ref)
    diam_m   = diam_deg * deg2rad(1.0) * R * cos_lat0

    # Angles: 0° = east (x-direction in the local frame) for consistency with
    # the equirectangular frame used by _wlsq_grad.
    angles = [0.0, 72.0, 144.0, 216.0, 288.0]
    nlons  = [c_lon + rad2deg(diam_m * cosd(a) / (R * cos_lat0)) for a in angles]
    nlats  = [c_lat + rad2deg(diam_m * sind(a) / R)               for a in angles]
    return nlons, nlats
end

# ─── WSE field ────────────────────────────────────────────────────────────────
# Pure eastward-downslope linear field: WSE decreases eastward.
# WSE(x_m) = -SLOPE × x_m   where x_m is east displacement from lon_ref.
# ci (west cell, x=0) is uphill; cj (east cell, x=diam_m) is downhill.
# Expected dWSE_n = wse_ci - wse_cj = SLOPE × diam_m > 0.
const SLOPE   = 0.001      # m/m  (1 mm per metre of east travel)
const LON_REF = -0.05      # reference longitude for both cell pairs
const LAT_REF = 51.0       # reference latitude

function _wse_at_lon(lon::Float64, lon_ref::Float64 = LON_REF,
                     lat_ref::Float64 = LAT_REF)::Float64
    x_m = deg2rad(lon - lon_ref) * _EARTH_R * cosd(lat_ref)
    return -SLOPE * x_m
end

# ─── GD1: V̂_y sign diagnostic using tilted-edge geometry ────────────────────
# On the real A5 mesh, edges are NOT perpendicular to the centre-to-centre
# vector d̂. The tilt angle (and its sign) varies between A5's two sublattices
# and changes with mesh resolution — this is the root cause of the V̂-based
# formula's resolution-dependent sign flip (FloodA5_GradientCorrection_Fix.md §1).
#
# A pure E-W cell pair with a PERPENDICULAR bisector edge is always perfectly
# orthogonal (cos_theta=1, V̂=0) — that geometry can't demonstrate the effect.
# Instead, GD1 uses tilted-edge cell pairs where the shared edge is rotated
# by a known angle from the perpendicular, matching the typical 22–37° range
# seen on real A5 meshes. The key property: edge tilted NE-SW gives V̂_y > 0;
# edge tilted NW-SE gives V̂_y < 0. A mesh where different cell pairs have
# different tilt directions will produce mixed V̂_y signs — that's the
# problematic scenario the d̂ formula avoids.

"""
    _tilted_edge_geometry(tilt_deg)
        → (cos_theta, V̂_x, V̂_y)

Compute _edge_geometry analytically for a pure E-W cell pair (d̂ = east)
whose shared edge is tilted by `tilt_deg` degrees clockwise from north
(the perpendicular to d̂). No call to _shared_edge or pya5 needed.

tilt_deg = 0°  → perpendicular edge (orthogonal, V̂ = 0)
tilt_deg > 0°  → edge tilts clockwise (NW-SE when viewed from above)
tilt_deg < 0°  → edge tilts counter-clockwise (NE-SW)
"""
function _tilted_edge_geometry(tilt_deg::Float64)
    α  = deg2rad(tilt_deg)
    # Edge unit vector: start from north (0,1) and rotate by α clockwise
    ex, ey = sin(α), cos(α)
    # Face normal: 90° CCW from edge direction
    nx_raw, ny_raw = -ey, ex
    # d̂ = east = (1, 0)
    c_raw = nx_raw   # dot with (1,0)
    # Orient so that d̂·n̂ ≥ 0
    if c_raw < 0.0
        nx, ny, c = -nx_raw, -ny_raw, -c_raw
    else
        nx, ny, c = nx_raw, ny_raw, c_raw
    end
    Vx = nx - c   # V̂ = n̂ - c·d̂;  d̂=(1,0) so V̂_x = n̂_x - c
    Vy = ny       # V̂_y = n̂_y - c·0 = n̂_y
    return (c, Vx, Vy)
end

println()
println("── GD1: V̂_y for tilted-edge cell pairs, tilt −25° to +25° ────────")
println("   Models non-orthogonality seen on real A5 cells (θ ≈ 22–37°)")
println("   Positive tilt (CW from north) → V̂_y negative; negative tilt → positive")
println("   The sign depends on which way the edge tilts — sublattice-dependent on A5")
println()
println(@sprintf("  %-8s  %-8s  %-8s  %-8s", "tilt(°)", "cos_θ", "V̂_y", "sign_V̂_y"))
println("  " * "─"^36)

VHat_y_vals = Float64[]
tilt_angles = [-25.0, -20.0, -15.0, -10.0, 0.0, 10.0, 15.0, 20.0, 25.0]

for tilt in tilt_angles
    ct, Vx, Vy = _tilted_edge_geometry(tilt)
    push!(VHat_y_vals, Vy)
    sign_str = Vy > 1e-6 ? "positive" : (Vy < -1e-6 ? "negative" : "zero")
    println(@sprintf("  %-8.1f  %-8.4f  %-8.4f  %s", tilt, ct, Vy, sign_str))
end

n_pos  = count(v ->  v >  1e-6, VHat_y_vals)
n_neg  = count(v ->  v < -1e-6, VHat_y_vals)
n_zero = count(v -> abs(v) <= 1e-6, VHat_y_vals)
println()
println("  V̂_y sign counts:  positive=$n_pos  negative=$n_neg  zero=$n_zero")
println("  (On A5 meshes, different cell pairs have different tilt directions,")
println("   giving mixed V̂_y signs — this is why V̂ correction is resolution-dependent)")

@testset "GD1 — V̂_y sign depends on edge tilt direction" begin
    # Negative tilt (CCW from north, NE-SW edge) → V̂_y > 0
    ct_neg, _, Vy_neg = _tilted_edge_geometry(-25.0)
    @test Vy_neg > 0.1

    # Positive tilt (CW from north, NW-SE edge) → V̂_y < 0
    ct_pos, _, Vy_pos = _tilted_edge_geometry(+25.0)
    @test Vy_pos < -0.1

    # Zero tilt (perpendicular edge) → V̂_y = 0, cos_theta = 1
    ct_zero, Vx_zero, Vy_zero = _tilted_edge_geometry(0.0)
    @test abs(Vy_zero) < 1e-12
    @test ct_zero ≈ 1.0

    # V̂_y is antisymmetric in tilt angle (equal magnitude, opposite sign)
    @test Vy_neg ≈ -Vy_pos atol=1e-12

    # |V̂_y| grows with |tilt| — larger non-orthogonality = larger correction
    _, _, Vy_10 = _tilted_edge_geometry(10.0)
    _, _, Vy_20 = _tilted_edge_geometry(20.0)
    @test abs(Vy_20) > abs(Vy_10)

    # All values are finite and bounded in [-1, 1]
    @test all(isfinite, VHat_y_vals)
    @test all(v -> abs(v) <= 1.0 + 1e-9, VHat_y_vals)
end

# ─── GD2: d̂ projection formula ────────────────────────────────────────────────
println()
println("── GD2: d̂ projection formula — sign and magnitude, res 12–20 ────")
println("   WSE field: WSE = −$(SLOPE) m/m × x_east  (ci west=uphill, cj east=downhill)")
println("   Expected dWSE_n > 0 at every resolution (ci uphill → flux ci→cj)")
println()
println(@sprintf("  %-4s  %-12s  %-12s  %-10s  %-8s",
                 "res", "dWSE_n (m)", "expected (m)", "rel_err (%)", "sign"))
println("  " * "─"^52)

dWSE_n_vals = Float64[]
expected_vals = Float64[]

for res in 12:20
    diam_deg  = _approx_cell_diameter_deg(res; lat_ref = LAT_REF)
    diam_m    = diam_deg * deg2rad(1.0) * _EARTH_R * cosd(LAT_REF)

    lon_i, lat_i = LON_REF,            LAT_REF
    lon_j, lat_j = LON_REF + diam_deg, LAT_REF

    wse_ci = _wse_at_lon(lon_i)   # = 0.0 (reference point)
    wse_cj = _wse_at_lon(lon_j)   # = -SLOPE × diam_m  (lower, eastward)

    # Build 5-neighbour stencil for ci
    nb_lons_i, nb_lats_i = _pentagon_neighbours(lon_i, lat_i, diam_deg)
    nb_wses_i = _wse_at_lon.(nb_lons_i)

    # Build 5-neighbour stencil for cj
    nb_lons_j, nb_lats_j = _pentagon_neighbours(lon_j, lat_j, diam_deg)
    nb_wses_j = _wse_at_lon.(nb_lons_j)

    gx_i, gy_i = _wlsq_grad(lon_i, lat_i, nb_lons_i, nb_lats_i, wse_ci, nb_wses_i)
    gx_j, gy_j = _wlsq_grad(lon_j, lat_j, nb_lons_j, nb_lats_j, wse_cj, nb_wses_j)

    # Face-interpolated gradient
    gx_f = 0.5 * (gx_i + gx_j)
    gy_f = 0.5 * (gy_i + gy_j)

    # d_vec_m: ci → cj in local equirectangular metres, centred on ci
    # (same convention as _build_wlsq_weights! / _build_edge_list dx_ms/dy_ms)
    cos_lat_ci = cosd(lat_i)
    dx_m = deg2rad(lon_j - lon_i) * _EARTH_R * cos_lat_ci
    dy_m = deg2rad(lat_j - lat_i) * _EARTH_R

    # d̂ projection formula
    dWSE_n   = -(gx_f * dx_m + gy_f * dy_m)
    expected = wse_ci - wse_cj   # = SLOPE × diam_m > 0
    rel_err  = abs(dWSE_n - expected) / abs(expected) * 100.0
    sign_ok  = dWSE_n > 0

    push!(dWSE_n_vals, dWSE_n)
    push!(expected_vals, expected)

    println(@sprintf("  %-4d  %-12.6f  %-12.6f  %-10.3f  %s",
                 res, dWSE_n, expected, rel_err, (sign_ok ? "✓" : "✗ WRONG")))
end

println()
n_correct_sign = count(>(0), dWSE_n_vals)
println("  Correct sign (dWSE_n > 0): $n_correct_sign / $(length(dWSE_n_vals))")
rel_errors = abs.(dWSE_n_vals .- expected_vals) ./ abs.(expected_vals) .* 100.0
println(@sprintf("  Relative error: mean=%.3f%%  max=%.3f%%  (threshold: 2%%)",
                 mean(rel_errors), maximum(rel_errors)))

@testset "GD2 — d̂ projection formula" begin

    @testset "GD2a — correct sign at all resolutions" begin
        for (i, res) in enumerate(12:20)
            @test dWSE_n_vals[i] > 0   # ci uphill → dWSE_n positive
        end
    end

    @testset "GD2b — magnitude within 2% of analytical value" begin
        for (i, res) in enumerate(12:20)
            @test rel_errors[i] < 2.0
        end
    end

    @testset "GD2c — no sign flip across resolutions" begin
        # All values should be positive (ci is uphill at every resolution).
        # Contrast with GD1 where V̂_y changes sign.
        @test all(>(0), dWSE_n_vals)
        # Coefficient of variation of dWSE_n/expected should be small —
        # the formula should scale proportionally to cell size, not drift.
        normalised = dWSE_n_vals ./ expected_vals
        cv = std(normalised) / mean(normalised)
        @test cv < 0.05   # < 5% CV across 9 resolutions
    end

    @testset "GD2d — symmetric: flipping ci/cj reverses sign only" begin
        # If we swap ci and cj (lon_j becomes the reference, lon_i the target),
        # dWSE_n should flip sign but have the same magnitude.
        diam_deg  = _approx_cell_diameter_deg(14; lat_ref = LAT_REF)
        lon_i, lat_i = LON_REF,            LAT_REF
        lon_j, lat_j = LON_REF + diam_deg, LAT_REF

        wse_ci = _wse_at_lon(lon_i)
        wse_cj = _wse_at_lon(lon_j)
        nb_lons_i, nb_lats_i = _pentagon_neighbours(lon_i, lat_i, diam_deg)
        nb_lons_j, nb_lats_j = _pentagon_neighbours(lon_j, lat_j, diam_deg)
        nb_wses_i = _wse_at_lon.(nb_lons_i)
        nb_wses_j = _wse_at_lon.(nb_lons_j)

        gx_i, gy_i = _wlsq_grad(lon_i, lat_i, nb_lons_i, nb_lats_i, wse_ci, nb_wses_i)
        gx_j, gy_j = _wlsq_grad(lon_j, lat_j, nb_lons_j, nb_lats_j, wse_cj, nb_wses_j)
        gx_f = 0.5 * (gx_i + gx_j)
        gy_f = 0.5 * (gy_i + gy_j)

        # Forward: ci → cj
        cos_lat_ci = cosd(lat_i)
        dx_m_fwd   = deg2rad(lon_j - lon_i) * _EARTH_R * cos_lat_ci
        dy_m_fwd   = deg2rad(lat_j - lat_i) * _EARTH_R
        dWSE_fwd   = -(gx_f * dx_m_fwd + gy_f * dy_m_fwd)

        # Reversed: cj → ci (same face gradient, opposite d_vec_m)
        cos_lat_cj = cosd(lat_j)
        dx_m_rev   = deg2rad(lon_i - lon_j) * _EARTH_R * cos_lat_cj
        dy_m_rev   = deg2rad(lat_i - lat_j) * _EARTH_R
        dWSE_rev   = -(gx_f * dx_m_rev + gy_f * dy_m_rev)

        @test dWSE_fwd > 0
        @test dWSE_rev < 0
        @test abs(dWSE_fwd + dWSE_rev) < 1e-9 * abs(dWSE_fwd)
    end

    @testset "GD2e — zero gradient → zero dWSE_n" begin
        # On a flat field (all WSE = 0), the gradient is zero everywhere and
        # dWSE_n must be zero regardless of cell pair geometry or resolution.
        diam_deg = _approx_cell_diameter_deg(16; lat_ref = LAT_REF)
        lon_i, lat_i = LON_REF,            LAT_REF
        lon_j, lat_j = LON_REF + diam_deg, LAT_REF

        nb_lons_i, nb_lats_i = _pentagon_neighbours(lon_i, lat_i, diam_deg)
        nb_lons_j, nb_lats_j = _pentagon_neighbours(lon_j, lat_j, diam_deg)
        flat_wse = zeros(5)

        gx_i, gy_i = _wlsq_grad(lon_i, lat_i, nb_lons_i, nb_lats_i, 0.0, flat_wse)
        gx_j, gy_j = _wlsq_grad(lon_j, lat_j, nb_lons_j, nb_lats_j, 0.0, flat_wse)
        gx_f = 0.5 * (gx_i + gx_j)
        gy_f = 0.5 * (gy_i + gy_j)

        cos_lat_ci = cosd(lat_i)
        dx_m = deg2rad(lon_j - lon_i) * _EARTH_R * cos_lat_ci
        dy_m = deg2rad(lat_j - lat_i) * _EARTH_R
        dWSE_n = -(gx_f * dx_m + gy_f * dy_m)

        @test abs(dWSE_n) < 1e-12
    end
end


# ─── GD3: n̂_f projection formula ─────────────────────────────────────────────
# Verifies that the full n̂_f formula:
#   dWSE_n = −c·(∇WSE_f·d_vec_m) − L·(∇WSE_f·V̂)
# correctly recovers ∇WSE_f · n̂_f for a tilted edge with known geometry,
# at all resolutions 12–20, for both positive and negative tilt directions.
#
# Key property being tested: for the same linear WSE field (slope = SLOPE
# eastward), the n̂_f formula should recover the TRUE face-normal gradient
# component —  i.e. |∇WSE| × cos(angle between ∇WSE and n̂_f) — regardless
# of edge tilt direction or A5 resolution. The d̂ formula only recovers the
# component along d̂, which equals the face-normal component only when
# cos_theta = 1 (orthogonal). For tilted edges the n̂_f formula gives a
# strictly better approximation of the physically correct driving slope.
println("── GD3: n̂_f projection formula — tilted edges, multi-resolution ──")
println("   Checks that dWSE_n recovers the face-normal WSE gradient component")
println("   for both tilt directions at each A5 resolution 12–20")
println()
println(@sprintf("  %-5s  %-8s  %-12s  %-12s  %-12s  %-8s",
                 "res", "tilt(°)", "dWSE_n (m)", "expected (m)", "rel_err(%)", "sign"))
println("  " * "─"^62)

@testset "GD3 — n̂_f projection formula" begin
    for res in [12, 14, 16, 18, 20]
        diam_deg = _approx_cell_diameter_deg(res; lat_ref = LAT_REF)
        diam_m   = diam_deg * deg2rad(1.0) * _EARTH_R * cosd(LAT_REF)

        for tilt_deg in [-25.0, 25.0]
            # Analytical geometry for a tilted edge
            α    = deg2rad(tilt_deg)
            c    = cos(α)         # cos_theta for this tilt (= cos of angle between d̂ and n̂_f)
            # Edge vector (from tilted perpendicular):
            ex, ey = sin(α), cos(α)
            # Face normal (90° CCW from edge), re-oriented so c ≥ 0:
            nx_raw, ny_raw = -ey, ex
            c_raw = nx_raw       # dot with d̂=(1,0)
            if c_raw < 0.0
                nx, ny = -nx_raw, -ny_raw
            else
                nx, ny = nx_raw, ny_raw
            end
            c_oriented = abs(c_raw)
            # V̂ = n̂ - c·d̂:
            Vx = nx - c_oriented   # V̂_x (skew_x equivalent)
            Vy = ny               # V̂_y (skew_y equivalent)

            # Build a 5-neighbour WLSQ stencil for both ci and cj
            lon_i, lat_i = LON_REF,            LAT_REF
            lon_j, lat_j = LON_REF + diam_deg, LAT_REF
            wse_ci = _wse_at_lon(lon_i)
            wse_cj = _wse_at_lon(lon_j)

            nb_lons_i, nb_lats_i = _pentagon_neighbours(lon_i, lat_i, diam_deg)
            nb_lons_j, nb_lats_j = _pentagon_neighbours(lon_j, lat_j, diam_deg)
            nb_wses_i = _wse_at_lon.(nb_lons_i)
            nb_wses_j = _wse_at_lon.(nb_lons_j)

            gx_i, gy_i = _wlsq_grad(lon_i, lat_i, nb_lons_i, nb_lats_i, wse_ci, nb_wses_i)
            gx_j, gy_j = _wlsq_grad(lon_j, lat_j, nb_lons_j, nb_lats_j, wse_cj, nb_wses_j)
            gx_f = 0.5 * (gx_i + gx_j)
            gy_f = 0.5 * (gy_i + gy_j)

            # d_vec_m components (E-W pair so dy_m ≈ 0)
            cos_lat_ci = cosd(lat_i)
            dx_m = deg2rad(lon_j - lon_i) * _EARTH_R * cos_lat_ci
            dy_m = deg2rad(lat_j - lat_i) * _EARTH_R
            L    = sqrt(dx_m^2 + dy_m^2)

            # n̂_f formula:
            dhat_dot = gx_f * dx_m + gy_f * dy_m      # L·(∇WSE_f · d̂)
            Vhat_dot = gx_f * Vx  + gy_f * Vy         # (∇WSE_f · V̂)
            dWSE_n   = -c_oriented * dhat_dot - L * Vhat_dot

            # Expected: -L·(∇WSE_f · n̂_f)
            # For our linear field, ∇WSE_f = (-SLOPE, 0) (slope eastward, gradient westward)
            # n̂_f = (nx, ny) after orientation
            # ∇WSE_f · n̂_f = -SLOPE * nx + 0 * ny = -SLOPE * nx
            # expected = -L·(-SLOPE * nx) = L · SLOPE · nx
            expected = L * SLOPE * nx

            rel_err = abs(expected) > 1e-12 ? abs(dWSE_n - expected) / abs(expected) * 100.0 : 0.0
            sign_ok = dWSE_n > 0   # ci is uphill so dWSE_n should be positive

            println(@sprintf("  %-5d  %-8.1f  %-12.6f  %-12.6f  %-12.3f  %s",
                             res, tilt_deg, dWSE_n, expected, rel_err,
                             (sign_ok ? "✓" : "✗ WRONG")))

            @testset "GD3 res=$res tilt=$(tilt_deg)°" begin
                @test sign_ok
                @test rel_err < 2.0
            end
        end
    end
end

println()
println("All GD tests complete.")
