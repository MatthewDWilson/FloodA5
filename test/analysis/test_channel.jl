# Channel-series metrics (analysis/channel.jl): synthetic meshes, logs, manifest.
@testset "A-channel — channel-series metrics" begin
    C = FA.CHANNEL
    HN = FA.channel_normal_depth()
    @test HN ≈ (C.Q * C.NMAN / (C.W * sqrt(C.S)))^(3 / 5)

    # Cells placed at chosen (s, n) for a channel at angle α, offset off.
    function chan_mesh(sn, α, off; area = 400.0)
        lon = Float64[]; lat = Float64[]
        for (s, n) in sn
            x = s * cosd(α) - (n + off) * sind(α)
            y = s * sind(α) + (n + off) * cosd(α)
            push!(lon, C.LON0 + rad2deg(x / (FA.R_EARTH * cosd(C.LAT0))))
            push!(lat, C.LAT0 + rad2deg(y / FA.R_EARTH))
        end
        nc = length(sn)
        FA.MeshGeom(string.(1:nc), lon, lat, zeros(nc), fill(area, nc), Matrix{Float64}[],
                    zeros(Int32, 0, nc), Dict{String,Vector{Float64}}(), "")
    end
    cell = 20.0                         # chan: |n| ≤ 20, flood: |n| > 30
    sn = [(-500.0, 0.0), (-200.0, 15.0), (50.0, -5.0),      # channel, in the reach
          (-300.0, 25.0),                                   # neither (bank cell)
          (-400.0, 50.0), (0.0, -80.0), (-100.0, 40.0),     # floodplain, in the reach
          (-700.0, 0.0), (150.0, 0.0), (-650.0, 60.0)]      # outside the reach
    for (α, off) in ((0.0, 0.0), (37.5, 11.1))
        mesh = chan_mesh(sn, α, off)
        s, n = FA.channel_coords(mesh; angle = α, offset = off)
        @test maximum(abs.(s .- first.(sn))) < 1e-6 && maximum(abs.(n .- last.(sn))) < 1e-6
        zone, chan, flood = FA.channel_masks(s, n; cell)
        @test zone == [true, true, true, true, true, true, true, false, false, false]
        @test findall(chan) == [1, 2, 3]
        @test findall(flood) == [5, 6, 7]

        # bed on the normal-flow profile, depth h_n + δ in the channel
        z = [C.Z0 - C.S * s[i] for i in eachindex(s)]
        δ = [0.10, 0.20, 0.30]
        h1 = zeros(length(sn)); h1[1:3] .= HN .+ δ
        h1[5] = 0.06; h1[6] = 0.04; h1[7] = 0.5; h1[8] = 3.0; h1[10] = 1.0   # 5 and 7 spill
        h0 = copy(h1); h0[2] -= 0.003; h0[8] -= 1.0                          # 8 not in reach
        m = FA.channel_metrics(mesh, h0, h1; angle = α, offset = off, cell, z)
        @test m.deta ≈ 0.2 atol = 1e-9
        @test m.spill ≈ 2 * 400.0 / (0.35 * C.L)
        @test m.unsteady ≈ 0.003 atol = 1e-12
        # a dry channel: deta and unsteady NaN; nothing wet on the floodplain: spill 0
        md = FA.channel_metrics(mesh, zeros(10), zeros(10); angle = α, offset = off, cell, z)
        @test isnan(md.deta) && isnan(md.unsteady) && md.spill == 0.0
        # constants as keyword overrides: a 1 m higher bed reference lowers deta by 1 m
        @test FA.channel_metrics(mesh, h0, h1; angle = α, offset = off, cell, z,
                                 Z0 = C.Z0 + 1.0).deta ≈ -0.8 atol = 1e-9
    end
    mesh0 = chan_mesh(sn, 0.0, 0.0; area = 0.0)
    @test_throws ErrorException FA.channel_metrics(FA.MeshGeom(mesh0.ids, mesh0.lon, mesh0.lat,
        mesh0.elevation, Float64[], mesh0.polygons, mesh0.neighbours, mesh0.static, ""),
        zeros(10), zeros(10); angle = 0.0, offset = 0.0, cell)        # no areas

    # Log parser
    logp = joinpath(FIX_DIR, "chan_synth.log")
    write(logp, """
    [ Info: starting
    [ Info:   step=   50  t=500.0s  dt=10.0s  wet=11  max_depth=0.751m  domain_vol=8507.7m³  vol_removed=1492.3m³  mb_err=0.0m³
    [ Info:   step=  100  t=999.3s  dt=9.88s  wet=18  max_depth=0.956m  domain_vol=14418.3m³  vol_removed=5567.7m³  mb_err=-0.0m³
    [ Info: an unrelated line t=5s
    [ Info:   step=  150  t=1500.0s  dt=9.73s  wet=20  max_depth=0.983m  domain_vol=18817.5m³  vol_removed=1.55677e4m³  mb_err=-0.0m³
    """)
    e = FA.parse_log_volumes(logp)
    @test [x.t for x in e] == [500.0, 999.3, 1500.0]
    @test [x.vol_removed for x in e] == [1492.3, 5567.7, 15567.7]
    @test FA.log_volume_at(e, 1000.0) == (999.3, 5567.7)
    @test FA.log_volume_at(logp, 0.0) == (500.0, 1492.3)
    @test all(isnan, FA.log_volume_at(empty(e), 10.0))
    @test FA.log_outflow(logp, 1000.0, 1500.0) ≈ (15567.7 - 5567.7) / (C.Q * (1500.0 - 999.3))
    @test FA.log_outflow(logp, 1000.0, 1500.0; Q = 2C.Q) ≈ (15567.7 - 5567.7) / (2C.Q * (1500.0 - 999.3))
    @test isnan(FA.log_outflow(logp, 1500.0, 1500.0))

    # Run level and the metrics.csv writer, on a manifest pointing at the fixture
    # (format 2: no mesh file needed). The fixture is no channel: check plumbing.
    r = FA.load_run(FIX_H5)
    nf = FA.nframes(r)
    rm1 = FA.channel_run_metrics(r; angle = 0.0, offset = 0.0, cell = 25.0)
    @test isnan(rm1.outflow) && (rm1.t0, rm1.t1) == (r.t[nf-1], r.t[nf])
    @test isequal(rm1.spill, FA.channel_metrics(r.mesh, FA.depth(r, nf - 1), FA.depth(r, nf);
                                                angle = 0.0, offset = 0.0, cell = 25.0).spill)
    md, rd = mkpath(joinpath(FIX_DIR, "chm")), mkpath(joinpath(FIX_DIR, "chr"))
    write(joinpath(md, "manifest.csv"),
          "case,mesh,res,cell_m,angle,offset_m,n_cells,inflow_lat,inflow_lon,mesh_file\n" *
          "c1,cart,0,25.000,0.0,0.000,144,0,0,/nonexistent.parquet\n")
    cp(FIX_H5, joinpath(rd, "chan_c1_li.h5"))
    cp(logp, joinpath(rd, "chan_c1_li.log"))
    buf = IOBuffer()
    rows = with_logger(NullLogger()) do                     # "missing …_swe.h5" warning
        FA.channel_series_metrics(md, rd; io = buf)
    end
    @test length(rows) == 1 && rows[1].solver == "li" && rows[1].case == "c1"
    lines = readlines(joinpath(rd, "metrics.csv"))
    @test lines[1] == "case,mesh,res,cell_m,angle,offset_m,solver,deta,spill,unsteady,outflow"
    @test length(lines) == 2 && startswith(lines[2], "c1,cart,0,25.000,0.0,0.000,li,")
    @test isequal(rows[1].outflow, FA.log_outflow(logp, r.t[nf-1], r.t[nf]))
    @test occursin("c1                 li  Δη", String(take!(buf)))
    @test_logs (:warn, r"missing") match_mode = :any FA.channel_series_metrics(md, rd; io = nothing, csv = nothing)
end
