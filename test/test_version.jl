# The version string lives in several places (VERSIONING.md "Practical
# guidance"): VERSION, FloodModel.jl's FLOODA5_VERSION, the README badge
# (alt text and image URL) and a CHANGELOG entry. Check they agree, so a
# bump that misses one fails CI. The README badge sat at 0.8.0 through
# 0.8.1–0.8.16, and the 0.9.0 fix first changed only its alt text.
# Plain text checks: nothing is included.
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
rd(f) = read(joinpath(ROOT, f), String)

@testset "version strings agree" begin
    v = strip(rd("VERSION"))
    @test occursin(r"^\d+\.\d+\.\d+$", v)

    m = match(r"const\s+FLOODA5_VERSION\s*=\s*\"([^\"]+)\"", rd("FloodModel.jl"))
    @test m !== nothing && m[1] == v

    readme = rd("README.md")
    b = match(r"!\[Version ([^\]]+)\]\(https://img\.shields\.io/badge/Version-([^-)]+)-", readme)
    @test b !== nothing
    if b !== nothing
        @test b[1] == v          # alt text
        @test b[2] == v          # the image actually shown
    end

    @test occursin("## [" * v * "]", rd("CHANGELOG.md"))
end
