# test/analysis/test_hazard.jl — analysis/hazard.jl (run by test/test_analysis.jl)

@testset "hazard: D·V, Defra HR, ARR classes" begin
    @test FA.depth_velocity(0.5, 2.0) == 1.0
    @test FA.depth_velocity(0.0, 3.0) == 0.0
    @test FA.depth_velocity(0.5, -1.0) == 0.0
    # FD2321/TR1 Table 3.1 rows: d ≤ 0.25; 0.25 < d ≤ 0.75 (v ≤ 2); d > 0.75 and/or v > 2
    for (lu, row) in (("pasture", (0.0, 0.0, 0.5)), ("woodland", (0.0, 0.5, 1.0)),
                      ("urban", (0.0, 1.0, 1.0)), ("fd2320", (0.5, 1.0, 1.0)), ("none", (0.0, 0.0, 0.0)))
        @test FA.debris_factor(0.2, 0.5, lu) == row[1]
        @test FA.debris_factor(0.5, 0.5, lu) == row[2]
        @test FA.debris_factor(1.0, 0.5, lu) == row[3]
        @test FA.debris_factor(0.5, 2.5, lu) == row[3]       # v > 2
        @test FA.debris_factor(0.2, 2.5, lu) == row[1]       # shallow stays in the first row
        @test FA.debris_factor(0.0, 1.0, lu) == 0.0          # dry
    end
    # HR = d (v + 0.5) + DF; FD2321 Figure 3.2 (DF = 0): d = 1 m, v = 1 m/s → 1.5
    @test FA.hazard_rating(1.0, 1.0, "none") ≈ 1.5
    @test FA.hazard_rating(1.0, 1.0, "urban") ≈ 2.5
    @test FA.hazard_rating(0.25, 0.0, "fd2320") ≈ 0.625                 # 0.125 + 0.5
    @test FA.hazard_rating(0.0, 1.0, "fd2320") == 0.0
    # ARR / AIDR Guideline 7-3 Table 2
    @test FA.arr_hazard_class(0.0, 1.0) == 0
    @test FA.arr_hazard_class(0.2, 1.0) == 1
    @test FA.arr_hazard_class(0.2, 2.5) == 5        # v > 2
    @test FA.arr_hazard_class(0.4, 1.0) == 2
    @test FA.arr_hazard_class(1.0, 0.5) == 3
    @test FA.arr_hazard_class(1.0, 0.9) == 4        # D·V 0.9
    @test FA.arr_hazard_class(1.5, 0.5) == 4
    @test FA.arr_hazard_class(3.0, 1.0) == 5
    @test FA.arr_hazard_class(5.0, 0.1) == 6        # d > 4
    @test FA.arr_hazard_class(1.5, 3.0) == 6        # D·V 4.5
end

@testset "hazard: Chen et al. (2019) toppling" begin
    uc(h; kw...) = FA.chen_critical_velocity(h; kw...)
    # Fig. 6 check values with the recommended θ = 75°, α = 3.064, β = 0.303
    # (projectinstructions/FloodA5_GlobeViewer_Plan_2026-11.md, appendix B)
    for (h, u) in ((0.4, 2.47), (0.7, 1.60), (1.0, 1.15), (1.2, 0.93))
        @test isapprox(uc(h), u; atol = 0.006)
    end
    for (h, u) in ((0.4, 1.57), (0.7, 0.98), (1.0, 0.63))
        @test isapprox(uc(h; mass = 25.5, height = 1.26), u; atol = 0.006)
    end
    @test isapprox(uc(0.7; mass = 50.0, height = 1.60), 1.33; atol = 0.006)
    # θ = 90°: eq. (11), Xia et al. (2014)
    let m = 60.0, hp = 1.7, a = 3.0, b = 0.3, h = 0.8
        eq11 = a * (h / hp)^b * sqrt(m / (1000 * h^2) - (0.735 / hp^2 + 0.265 / (h * hp)) * (1.015e-3 * m - 4.927e-3))
        @test uc(h; mass = m, height = hp, theta = 90.0, alpha = a, beta = b) ≈ eq11
    end
    # edge cases: dry, floating, negative radicand
    @test uc(0.0) == Inf
    @test uc(1.8) == 0.0
    @test uc(1.769; mass = 1000.0, height = 1.77) == 0.0     # (a2 m + b2) ρ > m: buoyancy wins
    p = FA.HazardParams()
    @test FA.hazard_degree(0.0, 2.0, p) == 0.0
    @test FA.hazard_degree(2.0, 0.0, p) == 1.0
    @test FA.hazard_degree(0.7, 0.8, p) ≈ 0.8 / uc(0.7)
    @test FA.hazard_degree(0.7, 5.0, p) == 1.0
    @test FA.hazard_degree(0.7, 0.8; mass = 25.5, height = 1.26) ≈ 0.8 / uc(0.7; mass = 25.5, height = 1.26)
end

@testset "hazard: parameters and derived variables" begin
    p = FA.hazard_params(Dict("landuse" => "urban", "mass" => "50", "height" => "abc", "theta" => "200"))
    @test p.landuse == "urban" && p.mass == 50.0 && p.height == 1.77 && p.theta == 90.0
    @test FA.hazard_params(Dict("landuse" => "moon")).landuse == "fd2320"
    @test FA.hazard_key(FA.HazardParams()) == FA.hazard_key(FA.hazard_params(Dict{String,String}()))
    d = [0.0, 0.3, 1.0, NaN]; v = [1.0, 1.0, 2.5, 1.0]
    dv = FA.derive_hazard("depth_velocity", d, v)
    @test eltype(dv) == Float32 && dv[1:3] ≈ Float32[0, 0.3, 2.5] && isnan(dv[4])
    hr = FA.derive_hazard("hazard_rating", d, v, FA.HazardParams(landuse = "urban"))
    @test hr[2] ≈ 0.3 * 1.5 + 1.0 && hr[3] ≈ 1.0 * 3.0 + 1.0
    @test FA.derive_hazard("hazard_class_arr", d, v)[1:3] == Float32[0, 1, 5]      # H1; v > 2: H5
    @test_throws DimensionMismatch FA.derive_hazard("hazard_degree", d, v[1:2])
    @test Set(keys(FA.HAZARD_VARS)) == Set(["depth_velocity", "hazard_rating", "hazard_degree", "hazard_class_arr"])
end

@testset "hazard: statistics over frames" begin
    t = [0.0, 60.0, 120.0, 180.0]
    D = [0.0 0.2 0.5 0.1;          # cell 1 over the frames
         0.0 0.0 0.0 0.0;          # cell 2: never wet
         1.0 1.0 1.0 1.0]          # cell 3: wet throughout
    V = [0.0 1.0 0.5 0.2; 0.0 0.0 0.0 0.0; 0.1 0.1 3.0 0.1]
    s = FA.run_statistics(t, k -> (D[:, k], V[:, k]), 3, FA.HazardParams(flood_depth = 0.15))
    @test Set(keys(s)) == Set(FA.STAT_VARS)
    @test s["max_water_depth"] == Float32[0.5, 0.0, 1.0]
    @test s["max_velocity"] == Float32[1.0, 0.0, 3.0]
    @test s["max_depth_velocity"] ≈ Float32[0.25, 0.0, 3.0]
    @test s["max_hazard_class_arr"] == Float32[2, 0, 5]      # d 1, v 3: H5
    # time of the maximum: frame 3 for cell 1; never wet: NaN; constant: the first frame
    @test s["time_max_water_depth"][1] == 120 && isnan(s["time_max_water_depth"][2]) && s["time_max_water_depth"][3] == 0
    @test s["time_first_flooded"][1] == 60 && isnan(s["time_first_flooded"][2]) && s["time_first_flooded"][3] == 0
    # midpoint weights 30, 60, 60, 30: cell 1 flooded at frames 2 and 3
    @test s["flood_duration"] == Float32[120, 0, 180]
    @test s["max_hazard_degree"][3] == 1.0f0
    # without velocity: depth statistics only
    s2 = FA.run_statistics(t, k -> (D[:, k], nothing), 3)
    @test !haskey(s2, "max_velocity") && s2["max_water_depth"] == Float32[0.5, 0.0, 1.0]
    @test FA.run_statistics([5.0], k -> ([0.3], [0.0]), 1)["flood_duration"] == Float32[0]
end

@testset "hazard: on a run (statistics_maps, hazard_field)" begin
    r = FA.load_run(FIX_H5)
    s = FA.statistics_maps(r)
    n = length(r.mesh)
    @test all(length(s[k]) == n for k in keys(s))
    hmax = reduce((a, k) -> max.(a, FA.depth(r, k)), 1:FA.nframes(r); init = fill(-Inf, n))
    @test s["max_water_depth"] ≈ Float32.(hmax)
    k = FA.nframes(r)
    @test FA.hazard_field(r, "depth_velocity", k) ≈ Float32.(FA.depth(r, k) .* FA.speed(r, k))
    @test maximum(s["max_hazard_degree"]) <= 1
    @test count(isfinite, s["time_first_flooded"]) == count(>=(0.1), hmax)
end
