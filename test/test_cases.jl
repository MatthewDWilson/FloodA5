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
              :sgs_strip_storage, :sgs_strip_profile, :sgs_strip_tables, :SGS_STRIP_REACH, :SGS_STRIP_N)
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
    # the registered SGS strip case marks only its SGS-behaviour checks
    xs = [ch for ch in FC.case("sgs_strip").checks if !isempty(ch.xfail)]
    @test !isempty(xs) && all(ch -> occursin("SGS", ch.name), xs)
end

println("\nTest-case registry tests complete.")
