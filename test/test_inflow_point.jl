# test/test_inflow_point.jl
# -------------------------
# Unit tests for dynamic fluvial inflow sources.
# Tests T-IP1 through T-IP12 covering hydrograph interpolation,
# cumulative volume, BDY/BCI file parsing, and source application.
#
# Run standalone:
#   julia --project=. test/test_inflow_point.jl
#
# The file guards against accidental execution when included by
# a harness that has already loaded the FloodModel symbols.

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
    include(joinpath(@__DIR__, "..", "FloodModel.jl"))
    ENV["FLOODMODEL_INCLUDE_ONLY"] = ""
end

using Test

# ---------------------------------------------------------------------------
# T-IP1: interpolation at knot points returns exact values
# ---------------------------------------------------------------------------
@testset "T-IP1: _interp_hydrograph at knot points" begin
    t_s   = [0.0, 3600.0, 7200.0]
    Q_m3s = [0.0, 10.0,   5.0]
    @test _interp_hydrograph(t_s, Q_m3s, 0.0)    ≈ 0.0
    @test _interp_hydrograph(t_s, Q_m3s, 3600.0) ≈ 10.0
    @test _interp_hydrograph(t_s, Q_m3s, 7200.0) ≈ 5.0
end

# ---------------------------------------------------------------------------
# T-IP2: linear interpolation between knots
# ---------------------------------------------------------------------------
@testset "T-IP2: linear interpolation between knots" begin
    t_s   = [0.0, 100.0]
    Q_m3s = [0.0, 10.0]
    @test _interp_hydrograph(t_s, Q_m3s, 50.0) ≈ 5.0
    @test _interp_hydrograph(t_s, Q_m3s, 25.0) ≈ 2.5
    @test _interp_hydrograph(t_s, Q_m3s, 75.0) ≈ 7.5
end

# ---------------------------------------------------------------------------
# T-IP3: flat extrapolation before first knot
# ---------------------------------------------------------------------------
@testset "T-IP3: flat extrapolation before first knot" begin
    t_s   = [100.0, 200.0]
    Q_m3s = [5.0,   10.0]
    @test _interp_hydrograph(t_s, Q_m3s, 0.0)  ≈ 5.0
    @test _interp_hydrograph(t_s, Q_m3s, 50.0) ≈ 5.0
end

# ---------------------------------------------------------------------------
# T-IP4: flat extrapolation after last knot
# ---------------------------------------------------------------------------
@testset "T-IP4: flat extrapolation after last knot" begin
    t_s   = [0.0, 100.0]
    Q_m3s = [0.0, 10.0]
    @test _interp_hydrograph(t_s, Q_m3s, 200.0) ≈ 10.0
    @test _interp_hydrograph(t_s, Q_m3s, 1e6)   ≈ 10.0
end

# ---------------------------------------------------------------------------
# T-IP5: apply_source! adds Q(t)×dt to correct cell
# ---------------------------------------------------------------------------
@testset "T-IP5: apply_source! adds correct volume" begin
    # Build a minimal stub FlowState just for the volume vector
    # (apply_source! only accesses state.volume[cell_index])
    n    = 5
    vol  = zeros(Float64, n)
    # Manually construct enough of FlowState to run apply_source!
    # We use setfield! to bypass the constructor for testing
    state = FlowState(
        fill("", n), zeros(n), vol, zeros(n), zeros(n), zeros(n),
        zeros(n), fill(0.03, n), fill(1e4, n), zeros(n), zeros(n),
        Dict{String,Vector{String}}(), zeros(Int,5,n),
        EdgeList(0,Int[],Int[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Int[],Int[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[]),
        Any[],
        falses(n), Any[], Any[], 0.0,
        zeros(Float64, 2, n), zeros(Float64, 10, n), 0.9, false,
        1.0,                                            # gradient_correction_alpha
        zeros(Float64, n), zeros(Float64, n),          # qvec_u, qvec_v
        zeros(Int, N_SIDES, n), zeros(Float64, 10, n), # cell_edge_index, mom_weights
        :edge,
        :legacy, nothing)   # face_flux_method, diamond_table (Phase C)

    t_s   = [0.0, 3600.0]
    Q_m3s = [5.0, 5.0]   # constant 5 m³/s
    src   = InflowPoint(3, "test", 0.0, 0.0, t_s, Q_m3s, "test")

    dt = 10.0
    apply_source!(state, src, 0.0, dt)
    @test state.volume[3] ≈ 5.0 * dt   # 50 m³
    @test all(state.volume[i] == 0.0 for i in [1,2,4,5])
end

# ---------------------------------------------------------------------------
# T-IP6: cumulative_volume matches trapezoidal integral
# ---------------------------------------------------------------------------
@testset "T-IP6: cumulative_volume trapezoidal integration" begin
    # Triangle hydrograph: 0 at t=0, peak 10 at t=3600, 0 at t=7200
    t_s   = [0.0,   3600.0, 7200.0]
    Q_m3s = [0.0,   10.0,   0.0]
    src   = InflowPoint(1, "test", 0.0, 0.0, t_s, Q_m3s, "test")

    # Full triangle area = 0.5 × base × height = 0.5 × 7200 × 10 = 36000 m³
    @test cumulative_volume(src, 7200.0) ≈ 36000.0 atol=1e-6

    # At t=3600 (peak): area under rising limb = 0.5 × 3600 × 10 = 18000 m³
    # (right triangle from origin to peak)
    v_half = cumulative_volume(src, 3600.0)
    @test v_half ≈ 18000.0 atol=1e-6

    # At t=1800 (quarter): area = 0.5 × 1800 × 5 = 4500 m³
    # (Q at t=1800 interpolates to 5.0 m³/s; triangle from 0 to 1800)
    @test cumulative_volume(src, 1800.0) ≈ 0.5 * 1800.0 * 5.0 atol=1e-6

    @test cumulative_volume(src, 0.0) == 0.0
end

# ---------------------------------------------------------------------------
# T-IP7: BDY file with 'hours' time unit converts correctly
# ---------------------------------------------------------------------------
# NOTE: .bdy data rows are "Value  Time" (discharge first, time second) per
# the LISFLOOD-FP manual §3.2.5 — NOT "Time  Value". This fixture previously
# had the columns the wrong way round (a docstring error in timeseries_io.jl
# was copied into this test), which would have made the *test* wrong while
# the production parser was actually correct all along. Confirmed against
# the official LISFLOOD-FP manual before fixing — see timeseries_io.jl.
@testset "T-IP7: LisfloodBDYReader hours→seconds" begin
    bdy_content = """
Test BDY header
GAUGE_MAIN
3	hours
0.0	0.0
5.5	1.0
3.2	2.0
"""
    tmpfile = tempname() * ".bdy"
    write(tmpfile, bdy_content)
    try
        reader  = LisfloodBDYReader(tmpfile)
        series  = read_timeseries(reader)
        @test haskey(series, "GAUGE_MAIN")
        t_s, Q  = series["GAUGE_MAIN"]
        @test length(t_s) == 3
        @test t_s[1] ≈ 0.0
        @test t_s[2] ≈ 3600.0
        @test t_s[3] ≈ 7200.0
        @test Q[2]   ≈ 5.5
    finally
        rm(tmpfile; force=true)
    end
end

# ---------------------------------------------------------------------------
# T-IP8: BDY file with two series parses both correctly
# ---------------------------------------------------------------------------
# NOTE: .bdy columns are "Value  Time" (discharge first) — see T-IP7 note
# above. Fixture corrected to preserve the intended scenario (discharge
# rising gently from 1.0 to 4.0 m³/s over 300s on UPSTREAM; 0.5 to 2.0 m³/s
# on DOWNSTREAM), which the original time-first ordering did not represent.
@testset "T-IP8: LisfloodBDYReader two series" begin
    bdy_content = """
Test BDY header
UPSTREAM
2	seconds
1.0	0.0
4.0	300.0
DOWNSTREAM
2	seconds
0.5	0.0
2.0	300.0
"""
    tmpfile = tempname() * ".bdy"
    write(tmpfile, bdy_content)
    try
        series = read_timeseries(LisfloodBDYReader(tmpfile))
        @test haskey(series, "UPSTREAM")
        @test haskey(series, "DOWNSTREAM")
        @test series["UPSTREAM"][2][2]   ≈ 4.0
        @test series["DOWNSTREAM"][2][2] ≈ 2.0
    finally
        rm(tmpfile; force=true)
    end
end

# ---------------------------------------------------------------------------
# T-IP9: BCI QVAR entry parses and links to series
# ---------------------------------------------------------------------------
# NOTE: .bdy columns are "Value  Time" — see T-IP7 note above. Fixture
# corrected to preserve the intended scenario (discharge rising 2.0 → 8.0
# m³/s over 600s).
@testset "T-IP9: parse_bci_file QVAR entry" begin
    bdy_content = """
Test BDY header
INFLOW_A
2	seconds
2.0	0.0
8.0	600.0
"""
    bci_content = "P  172.648  -43.386  QVAR  INFLOW_A\n"
    tmpbdy = tempname() * ".bdy"
    tmpbci = replace(tmpbdy, ".bdy" => ".bci")
    write(tmpbdy, bdy_content)
    write(tmpbci, bci_content)
    try
        entries, series = parse_bci_file(tmpbci; bdy_path=tmpbdy)
        @test length(entries) == 1
        @test entries[1].boundary_type == 'P'
        @test entries[1].bc_code       == "QVAR"
        @test entries[1].bc_value      == "INFLOW_A"
        @test entries[1].x1            ≈ 172.648
        @test entries[1].y1            ≈ -43.386
        @test haskey(series, "INFLOW_A")
        @test series["INFLOW_A"][2][2] ≈ 8.0
    finally
        rm(tmpbdy; force=true)
        rm(tmpbci; force=true)
    end
end

# ---------------------------------------------------------------------------
# T-IP10: BCI QFIX entry parses with correct rate
# ---------------------------------------------------------------------------
@testset "T-IP10: parse_bci_file QFIX entry" begin
    bci_content = "P  172.700  -43.400  QFIX  3.75\n"
    tmpbci = tempname() * ".bci"
    write(tmpbci, bci_content)
    try
        entries, _ = parse_bci_file(tmpbci)
        @test length(entries) == 1
        @test entries[1].bc_code  == "QFIX"
        @test entries[1].bc_value == "3.75"
        rate = tryparse(Float64, entries[1].bc_value)
        @test rate ≈ 3.75
    finally
        rm(tmpbci; force=true)
    end
end

# ---------------------------------------------------------------------------
# T-IP11: N/E/S/W entries produce unsupported warning (not an error)
# ---------------------------------------------------------------------------
@testset "T-IP11: BCI N/E/S/W entries log warning, not error" begin
    bci_content = "N  172.55  172.75  FREE\n"
    tmpbci = tempname() * ".bci"
    write(tmpbci, bci_content)
    try
        # Should not throw
        entries, _ = parse_bci_file(tmpbci)
        @test length(entries) == 1
        @test entries[1].boundary_type == 'N'
        # apply_bci_free_entries! should warn but not error
        # (tested indirectly; warning is checked via log capture in integration tests)
    finally
        rm(tmpbci; force=true)
    end
end

# ---------------------------------------------------------------------------
# T-IP12: Multiple InflowPoints on same cell sum correctly
# ---------------------------------------------------------------------------
@testset "T-IP12: two InflowPoints on same cell sum volumes" begin
    n   = 3
    vol = zeros(Float64, n)
    state = FlowState(
        fill("", n), zeros(n), vol, zeros(n), zeros(n), zeros(n),
        zeros(n), fill(0.03, n), fill(1e4, n), zeros(n), zeros(n),
        Dict{String,Vector{String}}(), zeros(Int,5,n),
        EdgeList(0,Int[],Int[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[],Int[],Int[],Float64[],Float64[],Float64[],Float64[],Float64[],Float64[]),
        Any[],
        falses(n), Any[], Any[], 0.0,
        zeros(Float64, 2, n), zeros(Float64, 10, n), 0.9, false,
        1.0,                                            # gradient_correction_alpha
        zeros(Float64, n), zeros(Float64, n),          # qvec_u, qvec_v
        zeros(Int, N_SIDES, n), zeros(Float64, 10, n), # cell_edge_index, mom_weights
        :edge,
        :legacy, nothing)   # face_flux_method, diamond_table (Phase C)

    src1 = InflowPoint(2, "c1", 0.0, 0.0, [0.0, 100.0], [3.0, 3.0], "A")
    src2 = InflowPoint(2, "c1", 0.0, 0.0, [0.0, 100.0], [7.0, 7.0], "B")
    dt   = 10.0

    apply_source!(state, src1, 0.0, dt)
    apply_source!(state, src2, 0.0, dt)

    @test state.volume[2] ≈ (3.0 + 7.0) * dt   # = 100 m³
end

@info "test_inflow_point.jl: all tests complete"
