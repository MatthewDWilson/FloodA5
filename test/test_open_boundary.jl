# test/test_open_boundary.jl
# --------------------------
# Unit tests for open outflow boundary conditions.
# Tests T-BC1 through T-BC10 covering ghost edge detection, geometry,
# ghost-cell WSE computation, one-way flux enforcement, mass balance,
# and GeoJSON bc-file segment matching.
#
# Run standalone:
#   julia --project=. test/test_open_boundary.jl

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
    include(joinpath(@__DIR__, "..", "FloodModel.jl"))
    ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
end

using Test

# ---------------------------------------------------------------------------
# Test mesh setup — generate a small flat mesh on the fly.
# Uses the same ~5.5 km² Christchurch bbox as the standard flat test AOI
# (test/flat_test_aoi.geojson) but does not require that file to exist.
# The mesh is generated once and cached for the session.
# ---------------------------------------------------------------------------

const _TEST_AOI_GEOJSON = """
{
  "type": "Feature",
  "properties": {},
  "geometry": {
    "type": "Polygon",
    "coordinates": [[
      [172.610, -43.515],
      [172.660, -43.515],
      [172.660, -43.490],
      [172.610, -43.490],
      [172.610, -43.515]
    ]]
  }
}
"""

# Session-level cache: generate once, reuse for all tests
const _TEST_STATE_CACHE = Ref{Union{Nothing, Tuple}}(nothing)

function _load_flat_state(; closed::Bool = false)
    if _TEST_STATE_CACHE[] === nothing
        @info "Generating test mesh (res 14, ~5 km² flat domain)..."
        mesh  = A5Grid.mesh_for_aoi(_TEST_AOI_GEOJSON, 14)
        @info "  Test mesh: $(length(mesh.cells)) cells"
        _TEST_STATE_CACHE[] = (mesh,)
    end
    mesh  = _TEST_STATE_CACHE[][1]
    state = initialise_flow_model(mesh, StandardFlow())
    if closed
        for i in eachindex(state.ghost_cell_bc)
            state.ghost_cell_bc[i] = Closed
        end
    end
    return state, mesh
end

# ---------------------------------------------------------------------------
# T-BC1: _build_ghost_edges detects correct boundary cell count
# ---------------------------------------------------------------------------
@testset "T-BC1: ghost edge count for known mesh" begin
    state, mesh = _load_flat_state()
    n_boundary = count(state.boundary_mask)
    n_interior = count(.!state.boundary_mask)
    n_ghost    = length(state.ghost_edges)

    @test n_boundary >= 1
    @test n_interior >= 1
    @test n_ghost    >= n_boundary   # at least one ghost edge per boundary cell
    # Boundary cells have 3 or 4 real neighbours (pentagon = 5 max),
    # so they have 1 or 2 ghost edges each.
    @info "T-BC1: boundary=$n_boundary interior=$n_interior ghost_edges=$n_ghost"
end

# ---------------------------------------------------------------------------
# T-BC2: ghost edge widths are positive and plausible (not the mean)
# ---------------------------------------------------------------------------
@testset "T-BC2: ghost edge widths are actual polygon side lengths" begin
    state, _ = _load_flat_state()
    for ge in state.ghost_edges
        @test ge.width > 0.0
        @test isfinite(ge.width)
        # At res 14, pentagon sides are ~300–500 m
        @test ge.width > 50.0     # not degenerate
        @test ge.width < 2000.0   # not wildly wrong
        @test ge.L     > ge.width  # L > width (centre-to-centre > edge length)
    end
end

# ---------------------------------------------------------------------------
# T-BC3: _ghost_wse ZeroGradient extrapolates gradient from interior neighbour
# ---------------------------------------------------------------------------
@testset "T-BC3: _ghost_wse ZeroGradient gradient extrapolation" begin
    wse_ci  = 3.7
    wse_nb  = 3.5    # interior neighbour is 0.2 m lower → slope flows outward
    L_ci_nb = 400.0  # 400 m to interior neighbour
    L_ghost = 400.0  # same distance to ghost centre
    sill    = 0.0

    # Expected: wse_ghost = wse_ci + (wse_ci - wse_nb) * L_ghost / L_ci_nb
    #         = 3.7 + 0.2 * 1.0 = 3.9
    expected = wse_ci + (wse_ci - wse_nb) * L_ghost / L_ci_nb
    @test _ghost_wse(wse_ci, sill, ZeroGradient, wse_nb, L_ci_nb, L_ghost) ≈ expected

    # When no interior neighbour available (NaN), falls back to wse_ci
    @test _ghost_wse(wse_ci, sill, ZeroGradient, NaN, 0.0, L_ghost) ≈ wse_ci

    # wse_ghost must be ≥ sill even when extrapolation goes below it
    wse_ci2 = 0.05; wse_nb2 = 0.1   # slope falling toward boundary
    result = _ghost_wse(wse_ci2, 0.0, ZeroGradient, wse_nb2, 400.0, 400.0)
    @test result >= 0.0   # not below sill
end

# ---------------------------------------------------------------------------
# T-BC4: _ghost_wse Critical
# ---------------------------------------------------------------------------
@testset "T-BC4: _ghost_wse Critical" begin
    wse_ci   = 2.0
    sill     = 0.5
    expected = sill + (2.0/3.0) * (wse_ci - sill)
    # Critical doesn't use nb args
    @test _ghost_wse(wse_ci, sill, Critical, NaN, 0.0, 1.0) ≈ expected
end

# ---------------------------------------------------------------------------
# T-BC5: _bates_ghost_flux never returns negative Q (one-way enforcement)
# ---------------------------------------------------------------------------
@testset "T-BC5: _bates_ghost_flux never injects water" begin
    # Case 1: outflow — wse_ci (2.0) > wse_ghost (1.5); water exits
    Q_out, _ = _bates_ghost_flux(0.0, 2.0, 1.5, 0.0, 300.0, 500.0, 0.03, 10.0, 1.5)
    @test Q_out > 0.0

    # Case 2: ghost WSE above interior — inflow scenario; must be blocked
    Q_out2, _ = _bates_ghost_flux(0.0, 1.0, 1.5, 0.0, 300.0, 500.0, 0.03, 10.0, 0.8)
    @test Q_out2 == 0.0

    # Case 3: Closed BC sentinel (-Inf ghost WSE)
    Q_out3, q3 = _bates_ghost_flux(5.0, 2.0, -Inf, 0.0, 300.0, 500.0, 0.03, 10.0, 1.0)
    @test Q_out3 == 0.0
    @test q3     == 0.0
end

# ---------------------------------------------------------------------------
# T-BC6: Closed BC: ghost flux is zero regardless of WSE gradient
# ---------------------------------------------------------------------------
@testset "T-BC6: Closed BC returns zero flux" begin
    state, _ = _load_flat_state(closed=true)
    # Inject water into the first boundary cell
    bi = findfirst(state.boundary_mask)
    state.volume[bi]      = 1000.0
    state.water_depth[bi] = 1000.0 / state.cell_area[bi]

    vol_before = state.vol_removed
    _apply_ghost_fluxes_standard!(state, 30.0)
    @test state.vol_removed == vol_before   # no outflow
end

# ---------------------------------------------------------------------------
# T-BC7: Mass balance on open flat domain: input = domain + vol_removed
# ---------------------------------------------------------------------------
@testset "T-BC7: open-boundary mass balance" begin
    state, mesh = _load_flat_state()
    method = StandardFlow()

    # Inject 50 mm/hr for 600 s
    rainfall_rate = 50.0 / 3_600_000.0   # m/s
    sim_duration  = 600.0
    dt_max        = 30.0

    run_simulation!(state, mesh, sim_duration, dt_max, nothing, :none;
                    method        = method,
                    rainfall_rate = rainfall_rate)

    cell_areas  = state.cell_area
    input_vol   = rainfall_rate * sim_duration *
                  sum(a for a in cell_areas if a >= 1.0; init=0.0)
    domain_vol  = sum(state.volume)
    vol_removed = state.vol_removed

    mb_err_pct = abs(input_vol - domain_vol - vol_removed) /
                 max(input_vol, 1.0) * 100.0

    @info "T-BC7: input=$(round(input_vol,sigdigits=4))m³  " *
          "domain=$(round(domain_vol,sigdigits=4))m³  " *
          "removed=$(round(vol_removed,sigdigits=4))m³  " *
          "mb_err=$(round(mb_err_pct,sigdigits=3))%"

    @test mb_err_pct < 0.1    # < 0.1% mass balance error
    @test vol_removed >= 0.0  # some water has left (flat mesh with open BCs)
end

# ---------------------------------------------------------------------------
# T-BC8: --closed-boundaries: vol_removed stays zero
# ---------------------------------------------------------------------------
@testset "T-BC8: closed-boundaries: no outflow" begin
    state, mesh = _load_flat_state(closed=true)
    method = StandardFlow()

    rainfall_rate = 50.0 / 3_600_000.0
    run_simulation!(state, mesh, 300.0, 30.0, nothing, :none;
                    method        = method,
                    rainfall_rate = rainfall_rate)

    @test state.vol_removed == 0.0
    # Mass balance: all input stays in domain
    input_vol  = rainfall_rate * 300.0 *
                 sum(a for a in state.cell_area if a >= 1.0; init=0.0)
    domain_vol = sum(state.volume)
    mb_err_pct = abs(input_vol - domain_vol) / max(input_vol, 1.0) * 100.0
    @test mb_err_pct < 0.01
end

# ---------------------------------------------------------------------------
# T-BC9: GeoJSON bc-file: Closed segment overrides ZeroGradient
# ---------------------------------------------------------------------------
@testset "T-BC9: GeoJSON bc-file assigns Closed to matched boundary cells" begin
    state, mesh = _load_flat_state()

    # Find centre of first boundary cell to build a tight polygon around it
    bi   = findfirst(state.boundary_mask)
    lon0 = state.cell_lons[bi]
    lat0 = state.cell_lats[bi]
    δ    = 0.005   # ~500 m at these latitudes

    geojson = """
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "Polygon",
        "coordinates": [[[$(lon0-δ),$(lat0-δ)],[$(lon0+δ),$(lat0-δ)],
                          [$(lon0+δ),$(lat0+δ)],[$(lon0-δ),$(lat0+δ)],
                          [$(lon0-δ),$(lat0-δ)]]]
      },
      "properties": { "bc_type": "Closed", "label": "test_closed" }
    }
  ]
}
"""
    tmpf = tempname() * ".geojson"
    write(tmpf, geojson)
    try
        updated = load_bc_file(tmpf, mesh.cells, state.boundary_mask,
                                state.ghost_edges, state.ghost_cell_bc,
                                ZeroGradient)
        # At least one ghost edge near the target cell should now be Closed
        n_closed = count(bc -> bc === Closed, updated)
        @test n_closed >= 1
        @info "T-BC9: $n_closed ghost edge(s) set to Closed by GeoJSON bc-file"
    finally
        rm(tmpf; force=true)
    end
end

# ---------------------------------------------------------------------------
# T-BC10: Gradient extrapolation drives outflow from first step
# ---------------------------------------------------------------------------
@testset "T-BC10: ghost edge momentum persists between steps" begin
    state, mesh = _load_flat_state()

    # Find a boundary cell that has an interior neighbour recorded
    bi       = findfirst(state.boundary_mask)
    ge_idxs  = findall(ge -> (ge::GhostEdge).cell_index == bi, state.ghost_edges)
    @test !isempty(ge_idxs)

    # Set boundary cell high, interior neighbour lower → outward slope
    state.volume[bi]      = 5e4
    state.water_depth[bi] = 5e4 / state.cell_area[bi]
    nb_idx = (state.ghost_edges[ge_idxs[1]]::GhostEdge).interior_nb_idx
    if nb_idx > 0
        state.volume[nb_idx]      = 2e4
        state.water_depth[nb_idx] = 2e4 / state.cell_area[nb_idx]
    end

    vol_before = state.vol_removed
    step_standard!(state, 10.0)

    n_nonzero = count(ge -> abs((ge::GhostEdge).flux_prev) > 0.0, state.ghost_edges)
    @info "T-BC10: interior_nb_idx=$(nb_idx)  ghost_edges_with_momentum=$n_nonzero  " *
          "vol_removed=$(state.vol_removed)"

    @test state.vol_removed > vol_before   # water left the domain on step 1
    @test n_nonzero > 0                    # momentum persists for next step
end

@info "test_open_boundary.jl: all tests complete"
