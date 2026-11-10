# Shape and direction metrics (analysis/shape.jl): synthetic fields with known
# answers, then a smoke test on the format-2 fixture run.

# Points on rings r = 5, 10, …, rmax at n angles (offset half a step from 0°),
# stretched by (a, b) along axes rotated by θ degrees.
function _rings(; rmax = 50.0, n = 360, a = 1.0, b = 1.0, θ = 0.0)
    xs, ys = Float64[], Float64[]
    for r in 5.0:5.0:rmax, k in 0:n-1
        φ = (k + 0.5) * 2π / n
        u, v = a * r * cos(φ), b * r * sin(φ)
        push!(xs, u * cosd(θ) - v * sind(θ)); push!(ys, u * sind(θ) + v * cosd(θ))
    end
    return xs, ys
end

@testset "A-shape — isotropic disc" begin
    x, y = _rings()
    h = ones(length(x))
    @test abs(FA.harmonic(h, x, y, 4)) < 1e-12
    @test abs(FA.harmonic(h, x, y, 6)) < 1e-12
    @test abs(FA.harmonic(h, x, y, 2; sinus = true)) < 1e-12
    @test isnan(FA.harmonic(zeros(length(x)), x, y, 4))
    @test FA.sector_ratio(h, x, y, 0.5) ≈ 1.0 atol = 1e-12
    @test isnan(FA.sector_ratio(h[1:10], x[1:10], y[1:10], 0.5))      # < 20 wet
    hc = FA.hull_circularity(x, y)
    @test 0.9999 < hc.circularity <= 1.0
    @test hc.area ≈ π * 50^2 rtol = 1e-3
    for centred in (true, false), n_dirs in (8, 16)
        b, rad, cnt = FA.front_radius_by_sector(x, y; n_dirs, centred)
        @test length(b) == n_dirs && all(cnt .> 0)
        @test all(r -> isapprox(r, 50.0; rtol = 1e-3), rad)
        cv, n_empty = FA.radius_cv(rad)
        @test cv < 1e-3 && n_empty == 0
    end
end

@testset "A-shape — elongation sign" begin
    # four-fold star r(φ) = 50 (1 + ε cos 4φ): ε > 0 reaches further along the
    # axes, ε < 0 along the diagonals (ratios 1.2, 1/1.2 up to the 1° sampling)
    star(ε) = (φs = [(k + 0.5) * 2π / 360 for k in 0:359];
               pts = [(s * 50 * (1 + ε * cos(4φ)), φ) for φ in φs for s in 0.1:0.1:1.0];
               ([r * cos(φ) for (r, φ) in pts], [r * sin(φ) for (r, φ) in pts]))
    xa, ya = star(0.2); xd, yd = star(-0.2)
    h = ones(length(xa))
    @test FA.harmonic(h, xa, ya, 4) > 0.05
    @test FA.harmonic(h, xd, yd, 4) < -0.05
    @test FA.sector_ratio(h, xa, ya, 0.5) ≈ 1.2 rtol = 0.02
    @test FA.sector_ratio(h, xd, yd, 0.5) ≈ 1 / 1.2 rtol = 0.02
    b, rad, _ = FA.front_radius_by_sector(xa, ya; centred = true)
    cd = FA.cardinal_diagonal_means(b, rad)
    @test cd.n_cardinal == 4 && cd.n_diagonal == 4
    @test cd.ratio ≈ 1.2 rtol = 0.02
    b, rad, _ = FA.front_radius_by_sector(xd, yd; centred = true)
    @test FA.cardinal_diagonal_means(b, rad).ratio ≈ 1 / 1.2 rtol = 0.02
    @test FA.radius_cv(rad)[1] > 0.05
    # ellipse along the NE–SW diagonal: lateral skew A2 > 0; NW–SE: A2 < 0
    he = ones(length(_rings()[1]))
    xe, ye = _rings(a = 1.5, θ = 45.0)
    @test FA.harmonic(he, xe, ye, 2; sinus = true) > 0.1
    xm, ym = _rings(a = 1.5, θ = -45.0)
    @test FA.harmonic(he, xm, ym, 2; sinus = true) < -0.1
    @test FA.hull_circularity(xe, ye).circularity < 0.97
end

@testset "A-shape — sectors, bearings, sampling" begin
    # one point per compass bearing at radius 10, plus one at bearing 30°
    bs = [0.0, 45, 90, 135, 180, 225, 270, 315]
    x = 10 .* sind.(bs); y = 10 .* cosd.(bs)
    b, rad, cnt = FA.front_radius_by_sector(x, y; centred = true)
    @test b == bs && cnt == ones(Int, 8) && all(rad .≈ 10)
    b2, _, cnt2 = FA.front_radius_by_sector([20sind(30.0)], [20cosd(30.0)]; centred = false)
    @test b2 == collect(22.5:45:337.5) && cnt2[1] == 1               # [0°, 45°)
    _, _, cnt3 = FA.front_radius_by_sector([20sind(30.0)], [20cosd(30.0)]; centred = true)
    @test cnt3[2] == 1                                               # NE, centred on 45°
    _, rad0, cnt0 = FA.front_radius_by_sector([0.0], [0.0])          # r < 1e-6 skipped
    @test all(iszero, cnt0) && isnan(FA.radius_cv(rad0)[1]) && FA.radius_cv(rad0)[2] == 8
    @test FA.cardinal_mask(bs) == [true, false, true, false, true, false, true, false]
    cd = FA.cardinal_diagonal_means(bs, [2.0, 1, 2, 1, 2, 1, 2, 0])
    @test cd.cardinal == 2 && cd.diagonal == 1 && cd.n_diagonal == 3 && cd.ratio == 2
    # sampling on a 1 m grid
    gx = repeat(-20.0:20.0, outer = 41); gy = repeat(-20.0:20.0, inner = 41)
    s = FA.sample_at_bearings(gx, gy, 10.0)
    @test [p.bearing for p in s] == bs
    @test (gx[s[3].cell], gy[s[3].cell]) == (10.0, 0.0)              # east
    @test (gx[s[1].cell], gy[s[1].cell]) == (0.0, 10.0)              # north
    @test all(p -> p.dist <= sqrt(0.5) + 1e-12, s)
end

@testset "A-shape — hull of a square" begin
    gx = repeat(0.0:10.0, outer = 11); gy = repeat(0.0:10.0, inner = 11)
    hull = FA.convex_hull(gx, gy)
    @test sort(hull) == [(0.0, 0.0), (0.0, 10.0), (10.0, 0.0), (10.0, 10.0)]
    @test FA.polygon_area_perimeter(hull) == (100.0, 40.0)
    @test FA.polsby_popper(100.0, 40.0) ≈ π / 4
    @test FA.polsby_popper(1.0, 0.0) == 0.0
    @test FA.polygon_area_perimeter(FA.convex_hull([0.0, 1.0], [0.0, 1.0])) == (0.0, 0.0)
end

@testset "A-shape — planar plume" begin
    x = [100.0, 100.0, 100.0, 200.0]
    y = [100.0, -100.0, 10.0, 0.0]
    V = [3.0, 1.0, 5.0, 0.0]
    s = FA.ns_split(V, y)
    @test (s.V_north, s.V_south, s.V_band, s.V_ns) == (3.0, 1.0, 5.0, 4.0)
    @test (s.n_north, s.n_south, s.n_band) == (1, 1, 2)
    @test s.asym ≈ 0.5
    @test FA.ns_split(V, y; band = 5.0).asym ≈ (8 - 1) / 9
    @test FA.ns_split(zeros(4), y).asym == 0.0
    @test FA.plume_deflection(V, x, y) ≈ atand(3 * 100 - 100 + 50, 900)
    @test FA.plume_deflection([1.0], [100.0], [100.0]) ≈ 45.0
end

@testset "A-shape — flow direction" begin
    n = 50
    u = fill(cosd(30.0), n); v = fill(sind(30.0), n); h = ones(n)
    h[1] = 0.0                                      # dry
    u[2] = v[2] = 0.0                               # wet but not moving
    ang, spd, n_wet, n_moving = FA.flow_directions(u, v, h)
    @test n_wet == n - 1 && n_moving == n - 2 && length(ang) == length(spd) == n - 2
    @test all(a -> a ≈ 30.0, ang) && all(s -> s ≈ 1.0, spd)
    m, sd, R = FA.circular_mean_std(ang)
    @test m ≈ 30.0 && R ≈ 1.0 && sd < 1e-5
    @test FA.circular_mean_std(ang; weights = spd)[1] ≈ 30.0
    # wrap-around: ±179° average to 180°, not 0°
    @test abs(FA.circular_mean_std([179.0, -179.0])[1]) ≈ 180.0
    # bimodal ±60°: mean 0° but R = cos 60° and only half within ±30°
    st = FA.direction_stats([60.0, -60.0, 60.0, -60.0], ones(4))
    @test abs(st.mean_dir) < 1e-9 && st.R ≈ 0.5 && st.R_w ≈ 0.5
    @test st.circ_std ≈ rad2deg(sqrt(-2log(0.5)))
    @test st.mean_abs_dev ≈ 60.0 && st.frac_within_30 == 0.0
    # speed weighting picks the fast cells' direction
    @test FA.circular_mean_std([0.0, 90.0]; weights = [3.0, 0.0])[1] ≈ 0.0
    @test FA.circular_mean_std([0.0, 180.0])[2] == Inf
    @test all(isnan, FA.circular_mean_std(Float64[]))
end

@testset "A-shape — shape_xy" begin
    lon = repeat(0.0:0.001:0.004, outer = 5); lat = repeat(51.0:0.001:51.004, inner = 5)
    n = length(lon)
    m = FA.MeshGeom(string.(1:n), lon, lat, zeros(n), Float64[], Matrix{Float64}[],
                    zeros(Int32, 0, n), Dict{String,Vector{Float64}}(), "")
    x, y = FA.shape_xy(m, 0.0021, 51.0019)                    # snaps to (0.002, 51.002)
    i = FA.source_cell(m, 0.0021, 51.0019)
    @test x[i] == 0.0 && y[i] == 0.0 && (lon[i], lat[i]) == (0.002, 51.002)
    xn, yn = FA.shape_xy(m, 0.0021, 51.0019; snap = false)
    @test xn[i] < 0 && yn[i] > 0
    @test (xn, yn) == FA.local_xy(lon, lat, 0.0021, 51.0019)
    _, y2 = FA.shape_xy(m, 0.0, 51.0; R = 6_371_008.8, snap = false)
    @test y2[end] ≈ deg2rad(0.004) * 6_371_008.8
    _, y3 = FA.shape_xy(m, 0.0, 51.0; snap = false, m_per_deg = FA.M_PER_DEG_LAT)
    @test y3 == (lat .- 51.0) .* 111_320.0
end

@testset "A-shape — Run level (format-2 fixture)" begin
    r = FA.load_run(FIX_H5)
    k = FA.nframes(r)
    lon0, lat0 = FIX_SRC
    sm = FA.spread_metrics(r, k; lon0, lat0)
    @test sm.t == r.t[k] && sm.n_wet > 20
    @test isfinite(sm.A4) && abs(sm.A4) < 0.3 && isfinite(sm.A6)
    @test length(sm.ratio) == 3 && 0.5 < sm.ratio[1] < 2
    ps = FA.point_spread_shape(r, k; lon0, lat0)
    @test ps.n_wet > 20 && 0.6 < ps.circularity <= 1 && ps.area > 0
    @test length(ps.radii) == 8 && ps.n_empty_sectors == 0 && ps.cv_radius < 0.5
    pm = FA.planar_metrics(r, k; lon0, lat0)
    @test pm.t == r.t[k] && -1 <= pm.asym <= 1 && isfinite(pm.A2) && isfinite(pm.deflection)
    @test pm.V_north + pm.V_south + pm.V_band ≈ FA.stored_volume(r, k) rtol = 1e-9
    if FA.has_field(r, "vel_u")
        fd = FA.flow_direction_metrics(r, k)
        @test fd.n_moving <= fd.n_wet && length(fd.angles) == fd.n_moving
        @test fd.n_moving == 0 || 0 <= fd.R <= 1
    end
    io = IOBuffer(); FA.print_frames(io, r)
    out = String(take!(io))
    @test startswith(out, "Frames in $(r.path):") && count('\n', out) == k + 2
    @test FA.check_frame(r, 1) === nothing
    @test_throws ErrorException FA.check_frame(r, 0)
    @test_throws ErrorException FA.check_frame(r, k + 1)
end
