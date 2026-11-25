# visualisation/viewer-common/tools/make_about.jl — write viewer-common/about.json,
# the variable descriptions (io/variables.toml) for the viewer without a server
# (local files on a static web host; the server serves /api/about from the TOML).
#
#   julia --project=. visualisation/viewer-common/tools/make_about.jl
#
# test/test_vis_server.jl checks that about.json matches io/variables.toml.

import TOML, JSON3

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TOML_FILE = joinpath(ROOT, "io", "variables.toml")
const JSON_FILE = joinpath(@__DIR__, "..", "about.json")

# keys in sorted order, so the file only changes when the descriptions do
_sorted(x::AbstractDict) = [k => _sorted(x[k]) for k in sort!(collect(keys(x)))]
_sorted(x::AbstractVector) = Any[_sorted(v) for v in x]
_sorted(x) = x
_json(x::AbstractVector{<:Pair}) = "{" * join([JSON3.write(String(k)) * ":" * _json(v) for (k, v) in x], ",") * "}"
_json(x::AbstractVector) = "[" * join(_json.(x), ",") * "]"
_json(x) = JSON3.write(x)

function make_about(out = JSON_FILE)
    d = TOML.parsefile(TOML_FILE)
    io = IOBuffer()
    JSON3.pretty(io, _json(_sorted(d)), JSON3.AlignmentContext(indent = 1))
    write(out, String(take!(io)) * "\n")
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("wrote ", make_about())
end
