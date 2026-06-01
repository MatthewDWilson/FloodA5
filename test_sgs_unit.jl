# test_sgs_unit.jl — Minimal SGS physics unit test
#
# Tests the SGS solver with 5 synthetic cells arranged as a linear chain:
#
#   [1] source → [2] → [3] → [4] → [5]   (1-based Julia indices)
#    z=5m         z=6m   z=7m   z=8m   z=9m
#
# Each cell has a hand-crafted SGSTable using a flat-bed (bathtub) hypsometric
# model so that wse_from_volume is analytically exact.  Edge sills are set 0.5m
# below the lower cell's z_min so that flow is possible once WSE exceeds the
# sill, but dry cells (volume=0) should NOT drive flux.
#
# Tests:
#   1. wse_from_volume / wetted_area_from_wse roundtrip accuracy
#   2. Dry cells drive no flux (Bug 48 regression)
#   3. Mass conservation over 100 steps with constant injection
#   4. Water flows downhill — source cell accumulates most volume
#   5. No sloshing — domain volume never exceeds injected volume

# Top-level imports — must be before the if block so macros are available
# when Julia parses the body (even though execution is guarded).
using Test, Printf

if abspath(PROGRAM_FILE) == @__FILE__

    project_dir      = dirname(abspath(@__FILE__))
    flood_model_path = joinpath(project_dir, "FloodModel.jl")
    isfile(flood_model_path) || error("FloodModel.jl not found at $flood_model_path")

    ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
    include(flood_model_path)
    delete!(ENV, "FLOODMODEL_INCLUDE_ONLY")
    _step_debug_count[] = typemax(Int)   # suppress first-step debug logging

    println("\n", "="^60)
    println("SGS unit test — 5-cell linear chain")
    println("="^60)

    # ── 1. Build synthetic SGSTables ──────────────────────────────────────────
    # Flat-bed (bathtub) model: every point within the cell is at z_min,
    # so area is constant and vol = cell_area * (wse - z_min).
    # This makes wse_from_volume analytically invertible: wse = z_min + V/area.

    function make_flat_sgs_table(z_min::Float64, z_max::Float64;
                                  cell_area::Float64 = 500.0, n_bins::Int = 20)
        bins       = collect(range(z_min, z_max; length = n_bins))
        vol_curve  = [cell_area * (b - z_min)          for b in bins]
        area_curve = [b > z_min ? cell_area : 0.0       for b in bins]
        return A5Grid.SGSTable(bins, vol_curve, area_curve, cell_area, z_min, z_max)
    end

    # z_min increases away from source — water naturally pools at cell 1
    z_mins     = [5.0, 6.0, 7.0, 8.0, 9.0]
    z_maxs     = z_mins .+ 3.0
    cell_area  = 500.0          # m² per cell
    n_cells    = 5
    tbls       = [make_flat_sgs_table(z_mins[i], z_maxs[i]; cell_area) for i in 1:n_cells]

    # ── 2. Verify SGSTable lookups ────────────────────────────────────────────
    @testset "SGSTable lookups" begin
        for (i, tbl) in enumerate(tbls)
            # Dry cell → z_min
            @test A5Grid.wse_from_volume(tbl, 0.0) == tbl.z_min

            # Full cell → z_max
            v_full = cell_area * (tbl.z_max - tbl.z_min)
            @test A5Grid.wse_from_volume(tbl, v_full) ≈ tbl.z_max atol=1e-6

            # Analytical roundtrip: wse = z_min + V/area  →  V = area*(wse-z_min)
            for frac in [0.1, 0.25, 0.5, 0.75, 0.9]
                V_test   = frac * cell_area * (tbl.z_max - tbl.z_min)
                wse_got  = A5Grid.wse_from_volume(tbl, V_test)
                wse_exp  = tbl.z_min + V_test / cell_area
                @test wse_got ≈ wse_exp atol=0.02   # tolerance = one bin width
            end

            # vol_curve is monotone non-decreasing
            diffs = diff(tbl.vol_curve)
            @test all(d -> d >= -1e-10, diffs)
        end
        println("  SGSTable lookups: OK")
    end

    # ── 3. Build FlowState manually ───────────────────────────────────────────
    ids = ["0000000000000$(lpad(i,3,'0'))" for i in 1:n_cells]

    # Linear chain adjacency: cell i neighbours cell i-1 and i+1
    adj_dict   = Dict{String,Vector{String}}()
    for i in 1:n_cells
        nbrs = String[]
        i > 1          && push!(nbrs, ids[i-1])
        i < n_cells    && push!(nbrs, ids[i+1])
        adj_dict[ids[i]] = nbrs
    end

    adj_matrix = zeros(Int, 5, n_cells)   # (max_nb=5) × n_cells
    for i in 1:n_cells
        for (slot, nb_id) in enumerate(adj_dict[ids[i]])
            j = findfirst(==(nb_id), ids)
            j !== nothing && (adj_matrix[slot, i] = j)
        end
    end

    # Edges: one per adjacent pair, cell_i < cell_j
    n_edges  = n_cells - 1
    e_ci     = collect(1:n_edges)
    e_cj     = collect(2:n_cells)
    e_width  = fill(10.0,  n_edges)   # 10 m wide shared edge
    e_L      = fill(30.0,  n_edges)   # 30 m centre-to-centre
    e_ct     = fill(1.0,   n_edges)   # perfectly orthogonal
    # Sill = min(z_min_i, z_min_j) - 0.5 m  (channel below both cell thalwegs)
    e_sill   = [min(z_mins[i], z_mins[i+1]) - 0.5 for i in 1:n_edges]
    e_flux   = zeros(Float64, n_edges)

    edges = EdgeList(n_edges, e_ci, e_cj, e_width, e_L, e_ct, e_sill, e_flux)

    # Lons/lats: ~30m spacing along a longitude line near (0,0)
    cell_lons = Float64[i * 0.00027 for i in 0:(n_cells-1)]   # ~30m per step
    cell_lats = zeros(Float64, n_cells)

    function fresh_state()
        FlowState(
            copy(ids),
            zeros(n_cells),     # water_depth
            zeros(n_cells),     # volume
            zeros(n_cells),     # velocity
            zeros(n_cells),     # vel_u
            zeros(n_cells),     # vel_v
            copy(z_mins),       # elevation = z_min
            fill(0.03, n_cells),# manning_n
            fill(cell_area, n_cells),
            copy(cell_lons),
            copy(cell_lats),
            deepcopy(adj_dict),
            copy(adj_matrix),
            EdgeList(n_edges, copy(e_ci), copy(e_cj), copy(e_width),
                     copy(e_L), copy(e_ct), copy(e_sill), copy(e_flux)),
            deepcopy(tbls),
        )
    end

    # ── 4. Bug 48 regression: dry cells drive no spurious flux ────────────────
    @testset "Bug 48 — dry cells drive no flux" begin
        state = fresh_state()
        # Inject 1 m³ into source cell only; all others dry
        state.volume[1] = 1.0

        vol_before = sum(state.volume)
        step_sgs!(state, 60.0)

        # Non-adjacent cells (3,4,5) must remain completely dry
        @test state.volume[3] ≈ 0.0 atol=1e-10
        @test state.volume[4] ≈ 0.0 atol=1e-10
        @test state.volume[5] ≈ 0.0 atol=1e-10

        # Mass must be conserved
        @test sum(state.volume) ≈ vol_before atol=1e-8

        # Cell 1 must not GAIN volume from its dry neighbour (cell 2 is higher, dry)
        # so volume[1] should be ≤ starting value (it can only lose to cell 2)
        @test state.volume[1] <= 1.0 + 1e-8

        println("  Bug 48 regression: dry cells drive no spurious flux — OK")
        println("    vol after 1 step: ", round.(state.volume; digits=5))
    end

    # ── 5. Mass conservation over 100 steps ───────────────────────────────────
    @testset "Mass conservation — 100 steps, constant injection" begin
        state          = fresh_state()
        rate_m3s       = 0.5       # m³/s injection into cell 1
        dt             = 30.0      # s
        n_steps        = 100
        total_injected = 0.0

        for _ in 1:n_steps
            added           = rate_m3s * dt
            state.volume[1] += added
            total_injected  += added
            step_sgs!(state, dt)
        end

        domain_vol = sum(state.volume)
        mb_err     = abs(domain_vol - total_injected)
        mb_pct     = 100.0 * mb_err / total_injected

        @test mb_pct < 0.01
        println("  Mass conservation:")
        println("    injected = $(round(total_injected; digits=2)) m³")
        println("    domain   = $(round(domain_vol;    digits=6)) m³")
        println("    error    = $(round(mb_err;        digits=8)) m³  " *
                "($(round(mb_pct; digits=6))%)")
    end

    # ── 6. Water flows downhill — correct volume gradient ────────────────────
    @testset "Downhill flow — source retains most volume" begin
        state = fresh_state()
        for _ in 1:300
            state.volume[1] += 0.5 * 30.0
            step_sgs!(state, 30.0)
        end

        vols = state.volume
        println("  Volume profile after 300 steps:")
        for i in 1:n_cells
            @printf("    cell %d (z_min=%.1f m): %.2f m³\n", i, z_mins[i], vols[i])
        end

        # Source (cell 1, lowest) should hold at least as much as the highest cell
        @test vols[1] >= vols[n_cells] - 1.0
        println("  Downhill flow: source volume ≥ far cell — OK")
    end

    # ── 7. No sloshing — free drainage conserves and decreases total ─────────
    @testset "No sloshing — free drainage" begin
        state = fresh_state()
        # Give source a pulse; no further injection
        state.volume[1] = 200.0
        initial_total   = sum(state.volume)

        max_any_cell = initial_total   # no cell should ever hold more than started
        for _ in 1:200
            step_sgs!(state, 30.0)
            for v in state.volume
                v > max_any_cell + 1e-6 && (max_any_cell = v)
            end
        end

        final_total = sum(state.volume)

        # Total should be conserved (no outflow BC, closed domain)
        @test final_total ≈ initial_total atol=1e-6

        # No individual cell should ever exceed the starting total (sloshing check)
        @test max_any_cell <= initial_total + 1e-6

        println("  No sloshing: max cell volume = $(round(max_any_cell; digits=4)) m³ " *
                "(started with $(initial_total) m³ total) — OK")
    end

    println("\n", "="^60)
    println("All SGS unit tests passed.")
    println("="^60, "\n")

end  # if abspath(PROGRAM_FILE) == @__FILE__
