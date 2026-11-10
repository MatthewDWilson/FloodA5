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

println("\nTest-case registry tests complete.")
