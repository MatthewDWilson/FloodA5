# test/test_notebooks.jl — the Jupyter notebooks in notebooks/ (committed
# without outputs) and notebooks/nbtools.jl.
#
#   julia --project=. test/test_notebooks.jl
#
# N1  nbtools writes notebooks exactly as Jupyter does (a notebook Jupyter
#     wrote, in test/fixtures, round-trips byte for byte), and strip clears
#     outputs, counts and execution metadata and is idempotent
# N2  every committed notebook: nbformat 4, the FloodA5 kernel, no outputs,
#     every code cell parses as Julia, and the first cell is a title
# (Executing the notebooks takes the case runs — hours for Carlisle — so it is
# not done here: `julia --project=. notebooks/nbtools.jl execute NB.ipynb`.)

include(joinpath(@__DIR__, "..", "notebooks", "nbtools.jl"))
using Test

const FIXTURE = joinpath(@__DIR__, "fixtures", "notebook_jupyter_written.ipynb")
_lf(s) = replace(s, "\r\n" => "\n")

@testset "N1 — nbtools writes notebooks as Jupyter does" begin
    io = IOBuffer(); NBTools._write(io, NBTools.read_nb(FIXTURE), 0); print(io, "\n")
    @test String(take!(io)) == _lf(read(FIXTURE, String))
    @test NBTools.has_outputs(NBTools.read_nb(FIXTURE))
    tmp = joinpath(mktempdir(), "nb.ipynb"); cp(FIXTURE, tmp)
    @test NBTools.strip!(tmp)
    nb = NBTools.read_nb(tmp)
    @test !NBTools.has_outputs(nb)
    code = [c for c in nb["cells"] if c["cell_type"] == "code"]
    @test all(c -> isempty(c["outputs"]) && c["execution_count"] === nothing &&
                   !haskey(c["metadata"], "execution"), code)
    s1 = read(tmp, String)
    @test !NBTools.strip!(tmp) && read(tmp, String) == s1       # idempotent
    @test [c["source"] for c in nb["cells"]] ==
          [c["source"] for c in NBTools.read_nb(FIXTURE)["cells"]]   # sources untouched
end

@testset "N2 — committed notebooks: kernel, no outputs, code parses" begin
    nbs = NBTools.notebooks()
    for p in nbs
        nb = NBTools.read_nb(p)
        name = basename(p)
        @test nb["nbformat"] == 4
        @test nb["metadata"]["kernelspec"]["name"] == NBTools.KERNEL
        @test !NBTools.has_outputs(nb)
        @test nb["cells"][1]["cell_type"] == "markdown" && startswith(join(nb["cells"][1]["source"]), "# ")
        for (i, c) in enumerate(nb["cells"])
            c["cell_type"] == "code" || continue
            ex = Meta.parseall(join(c["source"]))
            bad = any(a -> a isa Expr && a.head in (:error, :incomplete), ex.args)
            bad && @warn "$name cell $i does not parse"
            @test !bad
        end
    end
    println("checked ", length(nbs), " notebook(s)")
end

println("\nNotebook tests complete.")
