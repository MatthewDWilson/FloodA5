# notebooks/nbtools.jl — notebook housekeeping without Python tools.
#
#   julia --project=. notebooks/nbtools.jl strip [NOTEBOOK.ipynb …]      # clear outputs in place (default: all)
#   julia --project=. notebooks/nbtools.jl check [NOTEBOOK.ipynb …]      # error if any has outputs (what CI checks)
#   julia --project=. notebooks/nbtools.jl execute NOTEBOOK.ipynb [OUT.ipynb]
#       # run it headless on the FloodA5 kernel (needs Jupyter: see notebooks/README.md);
#       # the executed copy goes to OUT (default test/notebook_runs/<name>), the
#       # notebook itself is not modified
#
# Notebooks are committed without outputs. `strip` writes the JSON the way
# Jupyter does (sorted keys, one-space indent), so stripping a notebook that
# Jupyter or VS Code saved changes only the outputs.

module NBTools

using JSON3

const NB_DIR = @__DIR__
const KERNEL = "julia-flooda5"

notebooks(dir = NB_DIR) = sort([joinpath(dir, f) for f in readdir(dir) if endswith(f, ".ipynb")])

read_nb(path) = JSON3.read(read(path, String), Dict{String,Any})

# ── JSON output as nbformat writes it: json.dumps(nb, sort_keys=True, indent=1, ensure_ascii=False)
function _esc(io, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"'; print(io, "\\\"")
        elseif c == '\\'; print(io, "\\\\")
        elseif c == '\n'; print(io, "\\n")
        elseif c == '\r'; print(io, "\\r")
        elseif c == '\t'; print(io, "\\t")
        elseif c == '\b'; print(io, "\\b")
        elseif c == '\f'; print(io, "\\f")
        elseif c < ' '; print(io, "\\u", lpad(string(UInt16(c); base = 16), 4, '0'))
        else print(io, c)
        end
    end
    print(io, '"')
end
function _write(io, x, ind)
    pad(n) = print(io, '\n', ' '^n)
    if x isa AbstractDict
        isempty(x) && return print(io, "{}")
        print(io, '{')
        ks = sort!(collect(keys(x)))
        for (i, k) in enumerate(ks)
            pad(ind + 1); _esc(io, string(k)); print(io, ": "); _write(io, x[k], ind + 1)
            i < length(ks) && print(io, ',')
        end
        pad(ind); print(io, '}')
    elseif x isa AbstractVector
        isempty(x) && return print(io, "[]")
        print(io, '[')
        for (i, v) in enumerate(x)
            pad(ind + 1); _write(io, v, ind + 1)
            i < length(x) && print(io, ',')
        end
        pad(ind); print(io, ']')
    elseif x isa AbstractString; _esc(io, x)
    elseif x === nothing; print(io, "null")
    elseif x isa Bool; print(io, x ? "true" : "false")
    elseif x isa Integer; print(io, x)
    elseif x isa AbstractFloat; isinteger(x) && abs(x) < 1e15 ? print(io, Int(x), ".0") : print(io, x)
    else error("unsupported JSON value $(typeof(x))")
    end
end
write_nb(path, nb) = open(io -> (_write(io, nb, 0); print(io, '\n')), path, "w")

"Does the notebook carry outputs, execution counts or execution metadata?"
function has_outputs(nb)
    for c in nb["cells"]
        c["cell_type"] == "code" || continue
        isempty(get(c, "outputs", [])) || return true
        get(c, "execution_count", nothing) === nothing || return true
        haskey(get(c, "metadata", Dict()), "execution") && return true
    end
    return false
end

"Clear outputs, execution counts and execution metadata (in place); true if it changed."
function strip!(path)
    nb = read_nb(path)
    changed = has_outputs(nb)
    for c in nb["cells"]
        c["cell_type"] == "code" || continue
        c["outputs"] = Any[]
        c["execution_count"] = nothing
        md = get!(c, "metadata", Dict{String,Any}())
        delete!(md, "execution")
    end
    before = read(path, String)
    write_nb(path, nb)
    return changed || read(path, String) != before
end

"Run a notebook headless (copy → nbconvert --execute) and return the executed copy's path."
function execute(path; out = joinpath(dirname(NB_DIR), "test", "notebook_runs", basename(path)),
                 timeout = 7200)
    jup = _jupyter()
    mkpath(dirname(out))
    d = dirname(dirname(jup))                        # …/envs/jupyter
    sep = Sys.iswindows() ? ";" : ":"
    env = copy(ENV)
    env["PATH"] = join([d, joinpath(d, "Library", "bin"), dirname(jup), get(ENV, "PATH", "")], sep)
    cmd = `$jup nbconvert --to notebook --execute --ExecutePreprocessor.timeout=$timeout
           --ExecutePreprocessor.kernel_name=$KERNEL --output $(abspath(out)) $(abspath(path))`
    run(Cmd(setenv(cmd, env); dir = dirname(NB_DIR)))
    return out
end

function _jupyter()
    p = get(ENV, "FLOODA5_JUPYTER", "")
    isempty(p) || return p
    for c in (joinpath(homedir(), ".julia", "conda", "3", "x86_64", "envs", "jupyter", "Scripts", "jupyter.exe"),
              joinpath(homedir(), ".julia", "conda", "3", "envs", "jupyter", "bin", "jupyter"))
        isfile(c) && return c
    end
    j = Sys.which("jupyter")
    j === nothing && error("Jupyter not found: see notebooks/README.md, or set FLOODA5_JUPYTER")
    return j
end

function main(args)
    isempty(args) && (println(read(@__FILE__, String)[1:findfirst("\n\nmodule", read(@__FILE__, String))[1]]); return 0)
    cmd = args[1]
    files = length(args) > 1 ? args[2:end] : notebooks()
    if cmd == "strip"
        for f in files
            println(strip!(f) ? "stripped  " : "unchanged ", f)
        end
    elseif cmd == "check"
        bad = filter(f -> has_outputs(read_nb(f)), files)
        isempty(bad) || (println("notebooks with outputs (run `julia --project=. notebooks/nbtools.jl strip`):\n  ",
                                 join(bad, "\n  ")); return 1)
        println("ok: ", length(files), " notebook(s) without outputs")
    elseif cmd == "execute"
        out = length(args) > 2 ? execute(args[2]; out = args[3]) : execute(args[2])
        println("executed copy: ", out)
    else
        error("unknown command \"$cmd\" (strip, check, execute)")
    end
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    exit(NBTools.main(ARGS))
end
