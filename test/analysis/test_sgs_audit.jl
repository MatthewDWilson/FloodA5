# SGS audit functions (analysis/sgs_audit.jl): hypsometry, sub-basins, Courant.
@testset "A-sgs-audit — hypsometry, basins, Courant" begin
    # A plane z = x on a 100 × 100 grid over a 10 m × 10 m cell: V(η) = 10·η²/2 for 0 ≤ η ≤ 10.
    xs = [(i - 0.5) / 10 for i in 1:100]
    z = [x for x in xs, _ in 1:100]
    c = FA.HypsCurve(vec(z), 100.0)
    @test FA.hyps_volume(c, 4.0) ≈ 80.0 rtol = 1e-3
    @test FA.hyps_volume(c, -1.0) == 0.0
    @test FA.hyps_volume(c, 12.0) ≈ 100.0 * (12.0 - 5.0) rtol = 1e-9      # above the top: full area
    @test FA.hyps_wet_area(c, 4.0) ≈ 40.0
    @test FA.hyps_wse(c, 80.0) ≈ 4.0 atol = 0.01
    @test FA.hyps_wse(c, FA.hyps_volume(c, 3.21)) ≈ 3.21 atol = 1e-9      # inverse of the forward curve
    @test FA.hyps_wse(c, 0.0) == minimum(z)
    @test FA.hyps_wse(c, FA.hyps_volume(c, 12.0)) ≈ 12.0

    # The model's table from the same samples: exact at its levels, chord (≥) between them.
    eb, vc, ac = FA.uniform_level_table(c.z, 100.0; n_bins = 100)
    @test eb[1] == c.z[1] && eb[end] == c.z[end] && length(eb) == 100
    @test all(k -> isapprox(vc[k], FA.hyps_volume(c, eb[k]); atol = 1e-9), 1:100)
    @test ac[end] ≈ 100.0
    @test vc[1] == 0.0
    ηm = (eb[40] + eb[41]) / 2
    @test FA.sgs_table_volume(eb, vc, 100.0, ηm) >= FA.hyps_volume(c, ηm) - 1e-9  # convex curve, linear interpolation
    @test FA.sgs_table_volume(eb, vc, 100.0, eb[1] - 1) == 0.0
    @test FA.sgs_table_volume(eb, vc, 100.0, eb[end] + 2) ≈ vc[end] + 200.0
    @test FA.sgs_table_wse(eb, vc, 100.0, 0.0) == eb[1]
    @test FA.sgs_table_wse(eb, vc, 100.0, vc[end] + 200.0) ≈ eb[end] + 2
    for η in (0.5, 3.3, 9.0)
        @test FA.sgs_table_wse(eb, vc, 100.0, FA.sgs_table_volume(eb, vc, 100.0, η)) ≈ η atol = 1e-9
    end
    @test FA.sgs_table_wet_area(eb, ac, 100.0, 5.0) ≈ 50.0 atol = 1.0
    @test FA.sgs_table_wet_area(eb, ac, 100.0, 20.0) == 100.0
    # few levels on a high cell: the level for a thin layer is too low (the chord overstates V)
    e5, v5, _ = FA.uniform_level_table(c.z, 100.0; n_bins = 5)
    @test FA.sgs_table_wse(e5, v5, 100.0, 1.0) < FA.hyps_wse(c, 1.0) - 0.1

    # Sub-basins. A plane has one; a ridge across the cell makes two.
    t1 = FA.basin_merge_tree(z)
    @test t1.n == 10_000
    @test FA.significant_basins(t1; area = 0.01, depth = 0.01) == 0
    @test FA.lowest_point_basin(t1) === nothing
    @test t1.basins[t1.root].depth ≈ maximum(z) - minimum(z)
    # ridge at x = 6 (crest 3 m): west pool floor 0 m (60 % of the cell), east pool floor 1 m (40 %)
    zr = [x < 6 ? 3.0 * x / 6 : 1.0 + 2.0 * (10 - x) / 4 for x in xs, _ in 1:100]
    t2 = FA.basin_merge_tree(zr)
    side = [b for b in t2.basins if b.parent == t2.root && b.depth > 0.5]
    @test length(side) == 1
    @test side[1].area_frac ≈ 0.4 atol = 0.02
    @test side[1].depth ≈ 2.0 atol = 0.1                    # crest − its floor
    @test side[1].main_depth ≈ 3.0 atol = 0.1
    @test FA.significant_basins(t2; area = 0.2, depth = 0.5) == 1
    @test FA.significant_basins(t2; area = 0.5, depth = 0.5) == 0
    @test FA.significant_basins(t2; area = 0.2, depth = 2.5) == 0
    @test FA.significant_basins(t2; area = 0.2, depth = 1.5, two_sided = true) == 1
    @test FA.lowest_point_basin(t2) === nothing             # the lowest point is in the larger pool
    # the deeper pool is the smaller one: the lowest point is then outside the main basin
    zs = [x < 6 ? 1.0 + 2.0 * x / 6 : 3.0 * (10 - x) / 4 for x in xs, _ in 1:100]
    t3 = FA.basin_merge_tree(zs)
    lb = FA.lowest_point_basin(t3)
    @test lb !== nothing && isapprox(lb.area_frac, 0.4; atol = 0.02) && isapprox(lb.depth, 3.0; atol = 0.1)
    # a gap in the ridge near the pool floor joins the pools low down: no significant second basin
    zg = copy(zr); zg[:, 50:51] .= min.(zg[:, 50:51], 1.2)
    @test FA.significant_basins(FA.basin_merge_tree(zg); area = 0.2, depth = 0.5) == 0
    # a diagonal wall one point thick holds (4-connectivity)
    zd = [i == j ? 5.0 : 0.001 * (i + j) for i in 1:40, j in 1:40]
    @test FA.significant_basins(FA.basin_merge_tree(zd); area = 0.2, depth = 1.0) == 1
    # a flat crest over most of the cell beside sloping ground: a plateau is not a basin
    zp = [x < 4 ? 5.0 - 0.1 * x : 7.0 for x in xs, _ in 1:100]
    tp = FA.basin_merge_tree(zp)
    @test FA.significant_basins(tp; area = 0.01, depth = 0.1) == 0
    @test FA.lowest_point_basin(tp) === nothing
    # mask: only the west pool
    m = [x < 6 for x in xs, _ in 1:100]
    tm = FA.basin_merge_tree(zr, m)
    @test tm.n == count(m) && FA.significant_basins(tm; area = 0.01, depth = 0.01) == 0

    # Courant: dt √(g h)/ℓ
    e = FA.courant_exposure([1.0, 4.0, 0.0, NaN], [20.0, 20.0, 20.0, 20.0], 2.0)
    @test e.n == 2
    @test e.max ≈ 2.0 * sqrt(9.81 * 4.0) / 20.0
    @test e.share_above_0p7 == 0.0 && e.share_above_1 == 0.0
    @test e.dt_required ≈ 0.7 * 20.0 / sqrt(9.81 * 4.0)
    @test FA.courant_exposure([0.0], [20.0], 1.0).n == 0
    d = FA.dist_summary([1.0, 2.0, 3.0, NaN])
    @test d["n"] == 3 && d["median"] == 2.0 && d["max"] == 3.0
    @test FA.dist_summary(Float64[])["n"] == 0
end
