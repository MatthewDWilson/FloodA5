# test/test_cases.jl — the test-case registry (cases/FloodA5Cases.jl, FC)
#
#   FLOODA5_NO_VIS=1 julia --project=. test/test_cases.jl
#
# C1  every registered case is well formed: placeholders resolve, outputs
#     lie in the case directory, commands display, `--output` is the runner's
# C2  the smoke tier runs end to end through the runner (mesh generation,
#     run, analysis, physical checks, expected values) in a temp directory
# C3  a completed run is reused; a changed command line or a missing output
#     is run again; --dry-run runs nothing
# C4  expected failures (Check xfail): reported XFAIL/XPASS, never failing a
#     case; an ordinary failing check still fails it
# C5  the planar_wall case (a wall with no gap; SGS hops it) is registered,
#     well formed and carries its expected-failure marks — without running it

ENV["FLOODA5_NO_VIS"] = "1"
include(joinpath(@__DIR__, "..", "cases", "FloodA5Cases.jl"))
const FC = FloodA5Cases
using Test

const ROOT = mktempdir()
ENV["FLOODA5_CASES_DIR"] = ROOT

@testset "C1 — registered cases are well formed" begin
    cs = FC.cases()
    @test !isempty(cs) && !isempty(FC.cases(; tier = :smoke))
    for c in cs
        @test c.tier in (:smoke, :standard, :full)
        @test !isempty(c.jobs)
        for j in c.jobs
            argv = FC.job_argv(c, j)
            @test !any(a -> occursin('{', a), argv)             # placeholders resolved
            @test count(==("--output"), argv) == 1              # only the runner's
            @test startswith(argv[end], FC.case_dir(c))
        end
        for p in c.prep
            @test p.kind in (:floodmodel, :julia, :python, :write)
            @test all(o -> startswith(FC._subst(o, c), FC.case_dir(c)), p.outputs)
        end
        isempty(c.doc) || @test isfile(joinpath(FC.REPO, first(split(c.doc))))
        @test occursin(c.name, sprint(io -> FC.show_case(c; io)))
    end
    # the public names the notebooks use (notebooks/*.ipynb)
    for s in (:python_command, :python_env, :planar_params, :planar_point, :PLANAR_SLOPE_BAND,
              :EMBANKMENT_WET_DEPTH, :embankment_xrange, :log_mass_balance, :channel_series_table,
              :sgs_strip_params, :sgs_strip_analytic, :sgs_strip_zf, :sgs_strip_zrel, :sgs_strip_conveyance,
              :sgs_strip_storage, :sgs_strip_profile, :sgs_strip_tables, :SGS_STRIP_REACH, :SGS_STRIP_N,
              :MESH_FAMILIES, :MESH_FAMILIES_DX, :MESH_FAMILIES_CARLISLE, :CHECKERBOARD_EXTREMA,
              :CARLISLE_RISE, :CARLISLE_PEAK_WINDOW)
        @test isdefined(FC, s)
    end
    @test FC.planar_point() == FC._planar_point()
end

@testset "C1b — the SGS strip's analytic answers (no runs)" begin
    # cases/scripts/sgs_strip_geometry.jl on the case's own terrain parameters
    p = FC._SS_DEM
    a = FC.sgs_strip_analytic(p)
    @test isapprox(a["delta_2d_q3"], -0.267; atol = 1e-3)      # in bank: 0.733 m in the channel
    @test isapprox(a["delta_2d_q20"], 0.372; atol = 1e-3)      # overbank
    @test isapprox(a["delta_cells20_q3"], 0.045; atol = 1e-3)  # the 20 m cell-mean grid, in bank
    # full rows: the exact-table (lumped) and cell-mean conveyances coincide
    @test isapprox(a["delta_cells20_q20"], a["delta_lumped20_q20"]; atol = 1e-6)
    # above the floodplain the volume above the cell-mean beds is the water volume
    @test isapprox(a["store_cells20_q20"], FC.sgs_strip_storage(p, a["delta_cells20_q20"]); rtol = 1e-3)
    @test FC.sgs_strip_zrel(p, 0.0) ≈ -p.D && FC.sgs_strip_zrel(p, 10.0) ≈ 0.0
end

@testset "C2 — smoke tier end to end" begin
    for c in FC.cases(; tier = :smoke)
        io = IOBuffer()
        ok = FC.run_case(c; io)
        out = String(take!(io))
        ok || println(out)
        @test ok
        @test !occursin("FAIL", out)
        @test isfile(FC.expected_path(c))               # smoke cases carry expected values
        @test occursin("metrics match expected values", out)
    end
end

@testset "C3 — reuse, changed commands, dry run" begin
    c = first(FC.cases(; tier = :smoke))
    j = first(c.jobs)
    @test FC.job_status(c, j) === :done
    io = IOBuffer()
    @test FC.run_jobs(c; io)
    @test occursin("up to date", String(take!(io)))
    j2 = FC.Job(j.name, vcat(j.args, ["--manning-n", "0.031"]), j.threads, j.tags)
    @test FC.job_status(c, j2) === :changed
    c2 = FC.Case(name = c.name, title = c.title, tier = c.tier, prep = c.prep, jobs = [j2])
    io = IOBuffer()
    @test FC.run_jobs(c2; dry_run = true, io)
    s = String(take!(io))
    @test occursin("--manning-n 0.031", s) && FC.job_status(c, j) === :done   # nothing ran
    mv(FC.job_output(c, j), FC.job_output(c, j) * ".away")
    @test FC.job_status(c, j) === :missing
    mv(FC.job_output(c, j) * ".away", FC.job_output(c, j))

    # A run made in another checkout (worktree) of the repository is the same
    # job: its {repo} paths are rebased onto this one; case-root paths are not.
    sep = Sys.iswindows() ? "\\" : "/"
    other = joinpath(dirname(FC.REPO), "FloodA5-otherworktree")
    rec = [joinpath(other, "test", "carlisle", "x.bci"), "--sgs-bins", joinpath(FC.case_root(), "c", "m.parquet")]
    @test FC._rebase_argv(rec, other * sep) ==
          [joinpath(FC.REPO, "test", "carlisle", "x.bci"), "--sgs-bins", rec[3]]
    @test FC._rebase_argv(rec, FC.REPO) == rec                     # same checkout: unchanged
    @test FC._rebase_argv(rec, "") == rec                          # cwd not recorded
    # a case root inside the recording checkout stays put
    withenv("FLOODA5_CASES_DIR" => joinpath(other, "test", "case_runs")) do
        inroot = joinpath(other, "test", "case_runs", "c", "m.parquet")
        @test FC._rebase_argv([inroot], other) == [inroot]
    end
    @test FC._same_argv(FC._rebase_argv([other * sep * "a.tif"], other), [joinpath(FC.REPO, "a.tif")])
    # Windows paths are case-insensitive (d:\… is D:\…); elsewhere case matters
    Sys.iswindows() &&
        @test FC._same_argv(FC._rebase_argv([lowercase(other) * sep * "a.tif"], other), [joinpath(FC.REPO, "a.tif")])
end

@testset "C4 — expected failures (xfail)" begin
    m = Dict("run" => Dict("x" => 1.0))
    mk(checks) = FC.Case(name = "xfail_test_case", title = "t", jobs = FC.Job[], checks = checks)
    io = IOBuffer()
    @test FC.check(mk([FC.Check("known defect", m -> false; xfail = "documented")]), m; io)
    s = String(take!(io)); @test occursin("XFAIL", s) && occursin("documented", s)
    @test FC.check(mk([FC.Check("fixed defect", m -> true; xfail = "documented")]), m; io)
    @test occursin("XPASS", String(take!(io)))
    @test !FC.check(mk([FC.Check("ordinary", m -> false), FC.Check("known", m -> false; xfail = "d")]), m; io)
    # sgs_strip had eight expected failures until 0.7.0 (the SGS conveyance
    # defects) and planar_wall two (the hop over a wall inside a cell); the
    # checks stay, unweakened, and none is marked now
    for name in ("sgs_strip", "planar_wall")
        @test all(ch -> isempty(ch.xfail), FC.case(name).checks)
        @test count(ch -> occursin("SGS", ch.name), FC.case(name).checks) >= 4
    end
end

@testset "C5 — planar_wall: registered, well formed, checks unweakened (no runs)" begin
    c = FC.case("planar_wall")
    @test c.tier === :standard
    # 2 walls × 2 meshes × SGS/standard, the flow model always explicit, closed, no shared inputs
    @test length(c.jobs) == 8
    names = [j.name for j in c.jobs]
    @test Set(names) == Set("$(m)_$(f)_$(w)" for m in ("a5", "cart"), f in ("sgs", "std"), w in ("wide", "thin"))
    for j in c.jobs
        argv = FC.job_argv(c, j)
        w = j.tags[:wall]
        i = findfirst(==("--flow-model"), argv)
        @test i !== nothing && argv[i+1] == j.tags[:flow_model] && j.tags[:flow_model] in ("sgs", "standard")
        @test "--closed-boundaries" in argv
        mesh = argv[findfirst(==("--meshload"), argv)+1]
        @test dirname(mesh) == joinpath(FC.case_dir(c), w)        # the wall's own mesh
        @test occursin(j.tags[:flow_model], basename(mesh))
    end
    # one DEM per wall, in its own directory, without a notch; SGS meshes built as SGS
    dems = [p for p in c.prep if p.kind === :python]
    @test length(dems) == 2
    for p in dems
        i = findfirst(==("--notch-width"), p.args)
        @test i !== nothing && p.args[i+1] == "0"
    end
    @test length(unique(first(p.outputs) for p in dems)) == 2
    meshes = [p for p in c.prep if p.kind === :floodmodel]
    @test length(meshes) == 8 && all(p -> "--flow-model" in p.args && "--mesh-only" in p.args, meshes)
    # the design: the level pool stays below the crest with margin, for both walls
    for (w, geom) in FC._PW_WALLS
        g = FC._pw_design(geom)
        @test isapprox(g.crest, w == "wide" ? 8.0 : 8.05; atol = 1e-9) && g.freeboard >= FC._PW_MARGIN
        @test isapprox(g.wse_end - g.z_foot, sqrt(213000 / 5e4); atol = 1e-9)   # the wedge formula
    end
    # the thin wall is narrower than a cell and inside one Cartesian 20 m column (faces at
    # 981.1 and 1001.1 m in the generator's frame, where the 20 m are 20.02 m), not on a face
    thin = Dict(FC._PW_WALLS)["thin"]
    @test thin.width < 20 && thin.x - thin.width / 2 > 981.1 + 5 && thin.x + thin.width / 2 < 1001.1 - 5
    # checks: "no water downstream" for every run but standard flow on the thin wall (cell-mean
    # beds: checked against the ridge they see), stating the correct behaviour; expected
    nodown = [ch for ch in c.checks if occursin("no water downstream of the wall", ch.name)]
    @test length(nodown) == 6
    @test count(ch -> startswith(ch.name, "standard") && occursin("thin wall", ch.name) &&
                      occursin("cell-mean", ch.name), c.checks) == 2
    # no expected failures since 0.7.0 (connectivity): the two SGS thin-wall
    # checks are ordinary checks again
    @test all(ch -> isempty(ch.xfail), c.checks)
    @test count(ch -> startswith(ch.name, "SGS") && occursin("thin wall", ch.name), nodown) == 2
    # a run with water downstream fails its check; a dry one passes (the check is not weakened)
    job = "a5_sgs_thin"
    ch = only(ch for ch in nodown if startswith(ch.name, "SGS, A5") && occursin("thin wall", ch.name))
    @test ch.test(Dict(job => Dict("v_down_end" => 0.0)))
    @test !ch.test(Dict(job => Dict("v_down_end" => 5.0)))
    @test !ch.test(Dict{String,Dict{String,Float64}}())          # no metrics: not a pass
end

println("\nTest-case registry tests complete.")
