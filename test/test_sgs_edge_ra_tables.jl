"""
test_sgs_edge_ra_tables.jl — Stage 1 unit tests for SGS R-A hydraulic tables.

Tests the pre-computed per-edge flow area A(wse) and wetted perimeter P(wse)
curves introduced in the sgs_ra_flux branch, and the derived lookup functions
flow_area_from_wse, wetted_perim_from_wse, hydraulic_radius_from_wse.

Tests are purely synthetic — no DEM files required.  They construct SGSTable
objects directly from known elevation profiles and verify the curves against
analytical cross-section geometry.

Run:
    julia --project=. test/test_sgs_edge_ra_tables.jl
"""

# Top-level imports — must be before the if block so macros are available
# at parse time (Julia expands macros before executing any using statements).
using Test, Printf

if abspath(PROGRAM_FILE) == @__FILE__

    # Load FloodModel.jl in include-only mode (suppresses main() execution)
    # to get access to A5Grid.SGSTable and the new lookup functions.
    project_dir      = dirname(dirname(abspath(@__FILE__)))
    flood_model_path = joinpath(project_dir, "FloodModel.jl")
    isfile(flood_model_path) || error("FloodModel.jl not found at $flood_model_path")

    ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
    include(flood_model_path)
    delete!(ENV, "FLOODMODEL_INCLUDE_ONLY")

    println("=" ^ 62)
    println("SGS edge R-A hydraulic table unit tests")
    println("=" ^ 62)

    # ── Helper: build a synthetic SGSTable from a known edge profile ────────
    #
    # We construct the edge_area_curve and edge_perim_curve by hand using the
    # same incremental algorithm as build_sgs_tables! Step 2.  This lets us
    # verify the build logic against analytical solutions independently of
    # the mesh/DEM pipeline.
    #
    # edge_profile: sorted vector of DEM elevations along the edge arc
    # W:            total edge width (m)
    # elev_bins:    elevation knots for the curves (n_bins vector)
    # slot:         which adjacency slot (1–5) to fill (others stay zero)

    function _make_edge_curves(edge_profile::Vector{Float64},
                                W::Float64,
                                elev_bins::Vector{Float64},
                                slot::Int)
        n_bins = length(elev_bins)
        n_pts  = length(edge_profile)
        dx     = W / n_pts

        edge_area  = zeros(n_bins, 5)
        edge_perim = zeros(n_bins, 5)

        sorted_e = sort(edge_profile)
        cum_A    = 0.0;  cum_P = 0.0
        ptr      = 1;    n_wet = 0
        prev_wse = sorted_e[1]

        for k in 1:n_bins
            wse_k  = elev_bins[k]
            cum_A += n_wet * dx * (wse_k - prev_wse)
            while ptr <= n_pts && sorted_e[ptr] <= wse_k
                cum_A += dx * (wse_k - sorted_e[ptr])
                cum_P += dx
                n_wet += 1
                ptr   += 1
            end
            edge_area[k,  slot] = cum_A
            edge_perim[k, slot] = cum_P
            prev_wse = wse_k
        end
        return edge_area, edge_perim
    end

    function _make_table(elev_bins, edge_area, edge_perim;
                          z_min=elev_bins[1], z_max=elev_bins[end])
        n = length(elev_bins)
        A5Grid.SGSTable(
            elev_bins,
            zeros(n),        # vol_curve  (not used in these tests)
            zeros(n),        # area_curve (not used in these tests)
            1000.0,          # cell_area
            z_min,
            z_max,
            edge_area,
            edge_perim,
        )
    end

    # ────────────────────────────────────────────────────────────────────────
    # T-EA1: Rectangular flat-bottomed trench
    #
    # Profile: 20 sample points all at z = 5.0 m (flat trench bottom)
    # Edge width W = 10.0 m → dx = 0.5 m per sample
    # At WSE = 5.0 + d:
    #   A_expected = d × W = d × 10.0
    #   P_expected = W = 10.0  (all segments wet, walls vertical = dx-approximation)
    # ────────────────────────────────────────────────────────────────────────
    @testset "T-EA1 — Rectangular trench (flat bottom)" begin
        z_sill = 5.0
        W      = 10.0
        n_pts  = 20
        profile = fill(z_sill, n_pts)   # all samples at sill elevation

        # Elevation knots: sill to sill+3m in 100 steps
        n_bins    = 100
        elev_bins = collect(range(z_sill, z_sill + 3.0, length=n_bins))
        slot      = 1
        ea, ep    = _make_edge_curves(profile, W, elev_bins, slot)
        tbl       = _make_table(elev_bins, ea, ep; z_min=z_sill, z_max=z_sill+3.0)

        # Test at several depths
        for d in [0.1, 0.5, 1.0, 2.0, 2.9]
            wse = z_sill + d
            A_got = A5Grid.flow_area_from_wse(tbl, slot, wse)
            P_got = A5Grid.wetted_perim_from_wse(tbl, slot, wse)
            A_exp = d * W
            P_exp = W   # flat bottom: perimeter = full width

            @test isapprox(A_got, A_exp, rtol=0.02) broken=false
            @test isapprox(P_got, P_exp, rtol=0.02) broken=false
        end
    end

    # ────────────────────────────────────────────────────────────────────────
    # T-EA2: Dry edge — WSE at or below sill returns zero
    # ────────────────────────────────────────────────────────────────────────
    @testset "T-EA2 — Dry edge: A = P = 0 below sill" begin
        z_sill = 3.0
        W      = 8.0
        n_pts  = 16
        profile   = fill(z_sill, n_pts)
        elev_bins = collect(range(z_sill, z_sill + 2.0, length=50))
        slot      = 2
        ea, ep    = _make_edge_curves(profile, W, elev_bins, slot)
        tbl       = _make_table(elev_bins, ea, ep; z_min=z_sill, z_max=z_sill+2.0)

        @test A5Grid.flow_area_from_wse(tbl, slot, z_sill - 0.1)  == 0.0
        @test A5Grid.flow_area_from_wse(tbl, slot, z_sill - 1.0)  == 0.0
        @test A5Grid.wetted_perim_from_wse(tbl, slot, z_sill)     == 0.0
        @test A5Grid.hydraulic_radius_from_wse(tbl, slot, z_sill - 0.1) == 0.0
    end

    # ────────────────────────────────────────────────────────────────────────
    # T-EA3: Monotonically non-decreasing curves
    # ────────────────────────────────────────────────────────────────────────
    @testset "T-EA3 — Monotonicity" begin
        # Irregular profile (V-shaped channel)
        z_sill = 2.0
        W      = 12.0
        n_pts  = 24
        profile = [z_sill + abs(i - n_pts/2) * 0.1 for i in 1:n_pts]
        sort!(profile)

        elev_bins = collect(range(z_sill, z_sill + 3.0, length=100))
        slot      = 3
        ea, ep    = _make_edge_curves(profile, W, elev_bins, slot)
        tbl       = _make_table(elev_bins, ea, ep; z_min=z_sill, z_max=z_sill+3.0)

        A_prev = 0.0
        P_prev = 0.0
        for k in 1:length(elev_bins)
            wse = elev_bins[k]
            A   = A5Grid.flow_area_from_wse(tbl, slot, wse)
            P   = A5Grid.wetted_perim_from_wse(tbl, slot, wse)
            @test A >= A_prev - 1e-12   # non-decreasing
            @test P >= P_prev - 1e-12   # non-decreasing
            A_prev = A
            P_prev = P
        end
    end

    # ────────────────────────────────────────────────────────────────────────
    # T-EA4: Hydraulic radius bounds
    # ────────────────────────────────────────────────────────────────────────
    @testset "T-EA4 — Hydraulic radius: 0 ≤ R ≤ max_depth/2" begin
        z_sill    = 1.0
        W         = 5.0
        n_pts     = 10
        profile   = fill(z_sill, n_pts)
        max_depth = 4.0
        elev_bins = collect(range(z_sill, z_sill + max_depth, length=100))
        slot      = 1
        ea, ep    = _make_edge_curves(profile, W, elev_bins, slot)
        tbl       = _make_table(elev_bins, ea, ep; z_min=z_sill, z_max=z_sill+max_depth)

        for d in [0.01, 0.1, 0.5, 1.0, 2.0, 3.9]
            R = A5Grid.hydraulic_radius_from_wse(tbl, slot, z_sill + d)
            @test R >= 0.0
            # For a wide flat channel R = d × W / W = d — upper bound is d itself
            @test R <= d + 1e-6
        end
    end

    # ────────────────────────────────────────────────────────────────────────
    # T-EA5: Backward compatibility — legacy table (zeros) returns 0.0
    # ────────────────────────────────────────────────────────────────────────
    @testset "T-EA5 — Legacy table (zeros): no crash, returns 0.0" begin
        n_bins    = 50
        elev_bins = collect(range(0.0, 5.0, length=n_bins))
        zero_ea   = zeros(Float64, n_bins, 5)
        zero_ep   = zeros(Float64, n_bins, 5)
        tbl       = _make_table(elev_bins, zero_ea, zero_ep)

        for slot in 1:5
            @test A5Grid.flow_area_from_wse(tbl, slot, 2.5)         == 0.0
            @test A5Grid.wetted_perim_from_wse(tbl, slot, 2.5)      == 0.0
            @test A5Grid.hydraulic_radius_from_wse(tbl, slot, 2.5)  == 0.0
        end
    end

    # ────────────────────────────────────────────────────────────────────────
    # T-EA6: Unused adjacency slots stay zero (boundary cells)
    # ────────────────────────────────────────────────────────────────────────
    @testset "T-EA6 — Unused adjacency slots remain zero" begin
        z_sill    = 0.0
        W         = 10.0
        n_pts     = 20
        profile   = fill(z_sill, n_pts)
        elev_bins = collect(range(z_sill, z_sill + 3.0, length=100))
        slot      = 2   # only slot 2 is built
        ea, ep    = _make_edge_curves(profile, W, elev_bins, slot)
        tbl       = _make_table(elev_bins, ea, ep; z_min=z_sill, z_max=z_sill+3.0)

        # Slots 1, 3, 4, 5 were not built — should return zero
        for s in [1, 3, 4, 5]
            @test A5Grid.flow_area_from_wse(tbl, s, z_sill + 1.0)  == 0.0
            @test A5Grid.wetted_perim_from_wse(tbl, s, z_sill + 1.0) == 0.0
        end
        # Slot 2 should be non-zero above sill
        @test A5Grid.flow_area_from_wse(tbl, 2, z_sill + 1.0)  > 0.0
        @test A5Grid.wetted_perim_from_wse(tbl, 2, z_sill + 1.0) > 0.0
    end

    println()
    println("=" ^ 62)
    println("All SGS edge R-A table unit tests passed.")
    println("=" ^ 62)
end
