"""
    FloodA5Output

The FloodA5 HDF5 output format: writer and reader in one place. Format
version 2 (FloodA5 0.4.0) is self-describing, so that a run can be
reproduced exactly from its output file (`FloodModel.jl --rerun out.h5`) and
replayed or analysed without the original mesh file.

Layout (format_version 2)
-------------------------
    /                     attrs: format, format_version, flooda5_version, created_utc
    /run                  attrs: git_commit, git_dirty, julia_version, hostname, cwd,
                                 backend, threads, cell_order, flow_model, started_utc
      argv                JSON array: the command line (empty if not run from the CLI)
      options             JSON object: every option after the --scheme defaults
      git_diff            the uncommitted patch (only if the working tree was dirty)
      project, manifest   Project.toml / Manifest.toml text of the active environment
      inputs/NN           one group per input file: attrs role, flag, argv_value,
                          original_path, sha256, size, embedded; dataset `bytes`
    /mesh                 static data, written once
      cell_ids, center_lons, center_lats, elevations, <static vars>   (as in format 1)
      cell_area           plan area (m²)
      polygon_offsets     CSR offsets (n_cells+1, 0-based) into polygon_lon/lat
      polygon_lon, polygon_lat   cell polygon vertices, as stored in the mesh file
      neighbours          (n_cells × max_nb in C order) neighbour indices, 0-based, -1 = none
      edge_cells          (n_edges × 2) file-order cells of each edge (only with flux_Q)
      source/             the mesh file: attrs embedded ("full"|"geometry"), kind,
                          original_path, sha256, size; datasets `parquet` (bytes, full
                          only) and `prj`
    /results              time series, one dataset per variable
      t                   (n_frames,) simulation time (s)
      <var>               (n_frames × n_cells in C order; n_cells × n_frames in Julia)
                          water_depth, volume, saturation, velocity, vel_u, vel_v;
                          flux_Q (n_edges) only for SGS R-A runs
    /summary              attrs, written when the run ends: completed, n_steps, …

Chunks are one frame × up to 4096 cells (deflate 4, as format 1): a frame is a
contiguous read, a cell's time series one small read per frame.

Format 1 (FloodA5 ≤ 0.3.x: one `/frames/NNNNNN` group per output time, no
`/run`, no geometry) is still readable through the same functions.

Indices returned by the reader are 1-based Julia indices; the file stores
0-based indices (`start_index` attribute) so that Python readers need no
conversion.
"""
module FloodA5Output

using HDF5
using JSON3
using SHA
using Dates

const FORMAT_NAME    = "FloodA5 output"
const FORMAT_VERSION = 2
const CELL_CHUNK     = 4096
const EMBED_MAX_BYTES = 64 * 1024^2     # larger input files are referenced, not embedded

_filters() = [HDF5.Filters.Deflate(4)]   # as format 1; shuffle made Carlisle outputs 11% larger
_nowutc()  = Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") * "Z"
sha256_file(path::AbstractString) = bytes2hex(open(SHA.sha256, path))

# ═══════════════════════════════════════════════════════════════════════════
# Provenance
# ═══════════════════════════════════════════════════════════════════════════

function _git(repo::AbstractString, args::Vector{String})
    try
        return read(pipeline(Cmd(vcat(["git", "-C", repo], args)); stderr = devnull), String)
    catch
        return nothing
    end
end

"""
    provenance(repo_dir) → Dict

Code and environment identity: git commit, whether the working tree had
uncommitted changes (and the patch), Julia version, host, and the active
environment's Project.toml/Manifest.toml.
"""
function provenance(repo_dir::AbstractString)
    d = Dict{String,Any}()
    commit = _git(repo_dir, ["rev-parse", "HEAD"])
    d["git_commit"] = commit === nothing ? "unknown" : strip(commit)
    status = _git(repo_dir, ["status", "--porcelain", "--untracked-files=no", "--ignore-submodules"])
    d["git_dirty"] = status === nothing ? false : !isempty(strip(status))
    if d["git_dirty"]
        diff = _git(repo_dir, ["diff", "HEAD", "--ignore-submodules"])
        diff !== nothing && (d["git_diff"] = diff)
    end
    d["julia_version"] = string(VERSION)
    d["hostname"]      = gethostname()
    proj = Base.active_project()
    if proj !== nothing && isfile(proj)
        d["project"] = read(proj, String)
        man = joinpath(dirname(proj), "Manifest.toml")
        isfile(man) && (d["manifest"] = read(man, String))
    end
    return d
end

# ═══════════════════════════════════════════════════════════════════════════
# Writer
# ═══════════════════════════════════════════════════════════════════════════

"""
    create_output(path; flooda5_version, cells, run, inputs, mesh_source)

Create a format-2 output file (overwriting `path`) with its static parts.

- `cells`: NamedTuple with `ids`, `center_lons`, `center_lats`, `elevations`,
  `static_vars` (Dict name → vector), `id_scheme`, and optionally `cell_area`,
  `polygons` (per cell a vector of `[lon, lat]` vertices) and `neighbours`
  (max_nb × n_cells, 1-based, 0 = none) — all in the mesh file's cell order.
- `run`: Dict of run metadata (see the module docstring); `argv` (vector of
  strings) and `options` (Dict) are stored as JSON, `git_diff`, `project`,
  `manifest` as text, everything else scalar as attributes.
- `inputs`: vector of NamedTuples `(role, flag, argv_value, path)`.
- `mesh_source`: `nothing` or NamedTuple `(embed, path, original_path, kind)`;
  `embed` is `:full` (store the file's bytes) or `:geometry` (checksum only).
"""
function create_output(path::AbstractString; flooda5_version::AbstractString = "unknown",
                       cells, run::AbstractDict = Dict{String,Any}(),
                       inputs = NamedTuple[], mesh_source = nothing,
                       embed_max_bytes::Integer = EMBED_MAX_BYTES)
    h5open(path, "w") do f
        a = attributes(f)
        a["format"] = FORMAT_NAME
        a["format_version"] = FORMAT_VERSION
        a["flooda5_version"] = String(flooda5_version)
        a["created_utc"] = _nowutc()
        _write_run!(f, run, inputs, embed_max_bytes)
        _write_mesh!(f, cells, mesh_source)
        r = create_group(f, "results")
        create_dataset(r, "t", Float64, dataspace((0,); max_dims = (-1,)); chunk = (1024,))
    end
    return path
end

function _write_run!(f, run::AbstractDict, inputs, embed_max_bytes)
    g = create_group(f, "run")
    g["argv"]    = JSON3.write(get(run, "argv", String[]); allow_inf = true)
    g["options"] = JSON3.write(get(run, "options", Dict{String,Any}()); allow_inf = true)
    for k in ("git_diff", "project", "manifest")
        v = get(run, k, "")
        isempty(v) || (g[k] = String(v))
    end
    for (k, v) in run
        k in ("argv", "options", "git_diff", "project", "manifest") && continue
        attributes(g)[String(k)] = v isa Union{Symbol,AbstractString} ? String(v) : v
    end
    gi = create_group(g, "inputs")
    for (i, inp) in enumerate(inputs)
        gk = create_group(gi, lpad(string(i), 2, '0'))
        a = attributes(gk)
        a["role"] = String(inp.role); a["flag"] = String(inp.flag)
        a["argv_value"] = String(inp.argv_value)
        a["original_path"] = abspath(inp.path)
        if isfile(inp.path)
            sz = filesize(inp.path)
            a["size"] = sz
            a["sha256"] = sha256_file(inp.path)
            emb = sz <= embed_max_bytes
            a["embedded"] = emb
            emb && sz > 0 && (gk["bytes"] = read(inp.path))
        else
            a["size"] = -1; a["sha256"] = ""; a["embedded"] = false
        end
    end
end

function _write_mesh!(f, c, mesh_source)
    g = create_group(f, "mesh")
    n = length(c.ids)
    g["cell_ids"] = c.ids
    attributes(g)["id_scheme"] = c.id_scheme
    g["center_lons"] = c.center_lons
    g["center_lats"] = c.center_lats
    g["elevations"]  = c.elevations
    reserved = ("cell_ids", "center_lons", "center_lats", "elevations", "cell_area",
                "polygon_offsets", "polygon_lon", "polygon_lat", "neighbours",
                "edge_cells", "source")
    for (name, vals) in c.static_vars
        (name == "elevation" || name in reserved) && continue
        g[name] = vals
    end
    if haskey(c, :cell_area) && c.cell_area !== nothing
        g["cell_area"] = c.cell_area
    end
    if haskey(c, :polygons) && c.polygons !== nothing
        offs = Vector{Int64}(undef, n + 1); offs[1] = 0
        for i in 1:n
            offs[i+1] = offs[i] + length(c.polygons[i])
        end
        nv = offs[end]
        lon = Vector{Float64}(undef, nv); lat = Vector{Float64}(undef, nv)
        for i in 1:n, (j, v) in enumerate(c.polygons[i])
            lon[offs[i] + j] = v[1]; lat[offs[i] + j] = v[2]
        end
        g["polygon_offsets"] = offs
        attributes(g["polygon_offsets"])["start_index"] = 0
        if nv > 0
            ch = (min(nv, 4 * CELL_CHUNK),)
            write(create_dataset(g, "polygon_lon", Float64, (nv,); chunk = ch, filters = _filters()), lon)
            write(create_dataset(g, "polygon_lat", Float64, (nv,); chunk = ch, filters = _filters()), lat)
        end
    end
    if haskey(c, :neighbours) && c.neighbours !== nothing
        nb = Int32.(c.neighbours) .- Int32(1)        # 0 → -1 (none), k → k-1
        g["neighbours"] = nb
        attributes(g["neighbours"])["start_index"] = 0
        attributes(g["neighbours"])["fill_value"] = -1
    end
    mesh_source === nothing && return
    s = create_group(g, "source")
    a = attributes(s)
    a["embedded"] = String(mesh_source.embed)
    a["kind"] = String(mesh_source.kind)
    a["original_path"] = mesh_source.original_path === nothing ? "" : abspath(mesh_source.original_path)
    p = mesh_source.path
    if p !== nothing && isfile(p)
        a["sha256"] = sha256_file(p)
        a["size"] = filesize(p)
        mesh_source.embed === :full && (s["parquet"] = read(p))
        prj = splitext(p)[1] * ".prj"
        isfile(prj) && (s["prj"] = read(prj, String))
    else
        a["sha256"] = ""; a["size"] = -1
    end
end

function _create_series(r, name::AbstractString, n::Integer, k0::Integer)
    create_dataset(r, name, Float64, dataspace((n, k0); max_dims = (n, -1));
                   chunk = (min(max(n, 1), CELL_CHUNK), 1), filters = _filters())
end

"""
    append_frame!(path, t, vars; edge_vars = (), edge_cells = nothing)

Append one output time: `vars` and `edge_vars` are iterables of
`name => vector` (per cell in file order, per edge). An edge variable
(flux_Q) is created at its first appearance; earlier frames read as 0.
`edge_cells` (n_edges × 2 file-order cells, 1-based) is written once, with
the first edge variable.
"""
function append_frame!(path::AbstractString, t::Real, vars; edge_vars = (), edge_cells = nothing)
    h5open(path, "r+") do f
        r = f["results"]
        dt = r["t"]
        k = length(dt) + 1
        HDF5.set_extent_dims(dt, (k,))
        dt[k:k] = [Float64(t)]
        for (name, data) in Iterators.flatten((vars, edge_vars))
            n = length(data)
            ds = haskey(r, name) ? r[name] : _create_series(r, name, n, k - 1)
            HDF5.set_extent_dims(ds, (n, k))
            ds[:, k] = data
        end
        if !isempty(edge_vars) && edge_cells !== nothing && !haskey(f["mesh"], "edge_cells")
            f["mesh"]["edge_cells"] = Int32.(permutedims(edge_cells)) .- Int32(1)
            attributes(f["mesh"]["edge_cells"])["start_index"] = 0
        end
    end
    return nothing
end

"""
    write_summary!(path, summary::AbstractDict)

Write the end-of-run summary (`/summary` attributes). A file without
`/summary` belongs to a run that did not finish.
"""
function write_summary!(path::AbstractString, summary::AbstractDict)
    h5open(path, "r+") do f
        haskey(f, "summary") && delete_object(f, "summary")
        g = create_group(f, "summary")
        attributes(g)["finished_utc"] = _nowutc()
        for (k, v) in summary
            attributes(g)[String(k)] = v isa Symbol ? String(v) : v
        end
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════════════════
# Reader (format 2 and format 1)
# ═══════════════════════════════════════════════════════════════════════════

_reader(fn, path::AbstractString) = h5open(fn, path, "r")

"""`format_version(f)` → 2 (this format) or 1 (FloodA5 ≤ 0.3.x)."""
format_version(f::HDF5.File) =
    haskey(f, "results") ? Int(read(attributes(f)["format_version"])) :
    haskey(f, "frames")  ? 1 : error("$(HDF5.filename(f)) is not a FloodA5 output file")

_v1_keys(f) = sort(collect(keys(f["frames"])), by = x -> parse(Int, x))

"""Number of output frames."""
n_frames(f::HDF5.File) = format_version(f) >= 2 ? length(f["results/t"]) : length(keys(f["frames"]))

"""Simulation time (s) of every frame."""
frame_times(f::HDF5.File) = format_version(f) >= 2 ? read(f["results/t"]) :
    Float64[read(f["frames"][k]["t"]) for k in _v1_keys(f)]

"""Index of the frame closest to simulation time `t` (s)."""
frame_at(f::HDF5.File, t::Real) = argmin(abs.(frame_times(f) .- t))

"""Names of the per-frame variables (`t` excluded)."""
function frame_vars(f::HDF5.File)
    if format_version(f) >= 2
        return sort([k for k in keys(f["results"]) if k != "t"])
    end
    ks = _v1_keys(f)
    isempty(ks) && return String[]
    return sort([k for k in keys(f["frames"][ks[end]]) if k != "t"])
end

has_var(f::HDF5.File, var::AbstractString) = var in frame_vars(f)

"""
    read_frame(f, var, k = n_frames(f)) → Vector{Float64}

Variable `var` at frame `k` (1-based; default the last frame), in the mesh
file's cell order (edge order for flux_Q). A format-1 frame without flux_Q
reads as zeros.
"""
function read_frame(f::HDF5.File, var::AbstractString, k::Integer = n_frames(f))
    nf = n_frames(f)
    1 <= k <= nf || throw(BoundsError("frame $k of $nf"))
    if format_version(f) >= 2
        return f["results"][var][:, k]
    end
    g = f["frames"][_v1_keys(f)[k]]
    haskey(g, var) && return read(g[var])
    var == "flux_Q" && return zeros(Float64, size(f["mesh/edge_cells"], 2))
    error("frame $k has no variable '$var'")
end

"""
    read_var(f, var) → Matrix{Float64}   (n_cells × n_frames)

All frames of `var` at once.
"""
function read_var(f::HDF5.File, var::AbstractString)
    format_version(f) >= 2 && return read(f["results"][var])
    ks = _v1_keys(f)
    return reduce(hcat, (read(f["frames"][k][var]) for k in ks))
end

"""
    read_series(f, var, cells) → Matrix{Float64}   (n_frames × length(cells))

Time series of `var` at the given cells (1-based, file order).
"""
function read_series(f::HDF5.File, var::AbstractString, cells::AbstractVector{<:Integer})
    if format_version(f) >= 2
        ds = f["results"][var]
        return reduce(hcat, (vec(ds[c, :]) for c in cells))
    end
    return permutedims(read_var(f, var)[cells, :])
end
read_series(f::HDF5.File, var::AbstractString, cell::Integer) = vec(read_series(f, var, [cell]))

"""`cell_table(f)` → (ids, lon, lat, elevation), file order."""
cell_table(f::HDF5.File) = (ids = read(f["mesh/cell_ids"]), lon = read(f["mesh/center_lons"]),
                            lat = read(f["mesh/center_lats"]), elevation = read(f["mesh/elevations"]))

"""`cell_area(f)` → plan area (m²) per cell (format 2)."""
cell_area(f::HDF5.File) = (_need_v2(f, "cell areas"); read(f["mesh/cell_area"]))

"""
    cell_polygons(f) → Vector{Matrix{Float64}}

Polygon of every cell as a 2 × n_vertices matrix (row 1 lon, row 2 lat),
vertices as stored in the mesh file (format 2).
"""
function cell_polygons(f::HDF5.File)
    _need_v2(f, "cell polygons")
    offs = read(f["mesh/polygon_offsets"])
    lon = read(f["mesh/polygon_lon"]); lat = read(f["mesh/polygon_lat"])
    return [permutedims(hcat(lon[offs[i]+1:offs[i+1]], lat[offs[i]+1:offs[i+1]]))
            for i in 1:length(offs)-1]
end

"""`cell_neighbours(f)` → max_nb × n_cells, 1-based, 0 = none (format 2)."""
function cell_neighbours(f::HDF5.File)
    _need_v2(f, "neighbour lists")
    return Int.(read(f["mesh/neighbours"])) .+ 1
end

_need_v2(f, what) = format_version(f) >= 2 ||
    error("$(HDF5.filename(f)) is a format-1 file (FloodA5 ≤ 0.3.x): it has no $what — re-run it")

"""
    run_info(f) → Dict

Everything the file records about its run: version, git, environment,
command line (`argv`), resolved `options`, `inputs`, the mesh `source`, and
the `summary` (empty if the run did not finish). Format 1 files give only
`format_version`.
"""
function run_info(f::HDF5.File)
    d = Dict{String,Any}("format_version" => format_version(f))
    d["format_version"] >= 2 || return d
    for (k, v) in attrs_dict(f) ; d[k] = v; end
    g = f["run"]
    for (k, v) in attrs_dict(g); d[k] = v; end
    d["argv"] = String.(JSON3.read(read(g["argv"])))
    d["options"] = Dict{String,Any}(String(k) => v for (k, v) in JSON3.read(read(g["options"]); allow_inf = true))
    for k in ("git_diff", "project", "manifest")
        d[k] = haskey(g, k) ? read(g[k]) : ""
    end
    d["inputs"] = [attrs_dict(g["inputs"][k]) for k in sort(collect(keys(g["inputs"])))]
    d["mesh_source"] = haskey(f["mesh"], "source") ? attrs_dict(f["mesh/source"]) : Dict{String,Any}()
    d["summary"] = haskey(f, "summary") ? attrs_dict(f["summary"]) : Dict{String,Any}()
    return d
end

attrs_dict(obj) = Dict{String,Any}(k => read(attributes(obj)[k]) for k in keys(attributes(obj)))

for fn in (:format_version, :n_frames, :frame_times, :frame_at, :frame_vars, :has_var,
           :read_frame, :read_var, :read_series, :cell_table, :cell_area, :cell_polygons,
           :cell_neighbours, :run_info)
    @eval $fn(path::AbstractString, args...) = _reader(f -> $fn(f, args...), path)
end

# ═══════════════════════════════════════════════════════════════════════════
# Reproduction
# ═══════════════════════════════════════════════════════════════════════════

"""
    extract_mesh(f, dir) → path

Write the run's mesh file into `dir` (with its .prj) and return its path.
For a `geometry`-only file, return the original mesh path after checking
that it still has the recorded checksum.
"""
function extract_mesh(f::HDF5.File, dir::AbstractString)
    _need_v2(f, "embedded mesh")
    haskey(f["mesh"], "source") || error("$(HDF5.filename(f)) records no mesh source")
    s = f["mesh/source"]; a = attrs_dict(s)
    if haskey(s, "parquet")
        mkpath(dir)
        p = joinpath(dir, "mesh.parquet")
        write(p, read(s["parquet"]))
        haskey(s, "prj") && write(joinpath(dir, "mesh.prj"), read(s["prj"]))
        sha256_file(p) == a["sha256"] || error("embedded mesh fails its checksum")
        return p
    end
    p = a["original_path"]
    isempty(a["sha256"]) &&
        error("the mesh was not embedded (--embed-mesh geometry) and was never saved to a file")
    isfile(p) || error("the mesh was not embedded (--embed-mesh geometry) and $p no longer exists")
    sha256_file(p) == a["sha256"] ||
        error("the mesh was not embedded (--embed-mesh geometry) and $p has changed since the run")
    return p
end
extract_mesh(path::AbstractString, dir::AbstractString) = _reader(f -> extract_mesh(f, dir), path)

# CLI flags (FloodModel.jl main) that rerun_argv must know about.
const _BOOL_FLAGS = Set(["--closed-boundaries", "--dem-strict", "--limiter-stats", "--mesh-only"])
# Mesh-generation and output flags: dropped (the rerun loads the embedded mesh).
const _DROP_FLAGS = Set(["--meshgen", "--meshres", "--meshgen-cartesian", "--dx",
                         "--cartesian-lat-ref", "--meshout", "--meshload", "--dem",
                         "--dem-method", "--dem-samples", "--dem-seed", "--dem-strict",
                         "--sgs-bins", "--sgs-samples", "--mesh-only", "--output",
                         "--vis-port", "--embed-mesh"])
const _VIS_MODES = Set(["cesium", "makie"])

# Split argv into (flag, values) items; `--vis` takes an optional mode.
function _argv_items(argv::Vector{String})
    items = Tuple{String,Vector{String}}[]
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--vis"
            nxt = get(argv, i + 1, "")
            if nxt in _VIS_MODES
                push!(items, (a, [nxt])); i += 2
            else
                push!(items, (a, String[])); i += 1
            end
        elseif a in _BOOL_FLAGS || !startswith(a, "--")
            push!(items, (a, String[])); i += 1
        else
            push!(items, (a, i < length(argv) ? [argv[i+1]] : String[])); i += 2
        end
    end
    return items
end

"""
    rerun_argv(h5path; output, workdir, overrides = String[]) → Vector{String}

Command line that reproduces the run recorded in `h5path`: its stored
arguments with the mesh replaced by the embedded mesh (`--meshload`), every
input file replaced by its embedded copy (extracted into `workdir`), no
viewer, and `--output output`. `overrides` (further CLI arguments) are
appended and replace the stored value of any flag they repeat.
"""
function rerun_argv(h5path::AbstractString; output::AbstractString,
                    workdir::AbstractString = mktempdir(), overrides::Vector{String} = String[])
    h5open(h5path, "r") do f
        _need_v2(f, "run record")
        info = run_info(f)
        argv = info["argv"]
        isempty(argv) && error("$h5path was not written by a command-line run: no argv to rerun")
        mesh = extract_mesh(f, joinpath(workdir, "mesh"))
        # Input files: original argv value → path to use now
        remap = Dict{String,String}()
        extra = String[]
        gi = f["run/inputs"]
        for k in sort(collect(keys(gi)))
            g = gi[k]; a = attrs_dict(g)
            new = if a["embedded"]
                d = joinpath(workdir, "inputs", k); mkpath(d)
                p = joinpath(d, basename(a["original_path"]))
                write(p, haskey(g, "bytes") ? read(g["bytes"]) : UInt8[])
                p
            else
                p = a["original_path"]
                isfile(p) || error("input $(a["role"]) was not embedded and $p no longer exists")
                sha256_file(p) == a["sha256"] ||
                    @warn "input $(a["role"]) was not embedded and $p has changed since the run"
                p
            end
            if isempty(a["argv_value"])      # used implicitly (e.g. the .bdy next to a .bci)
                a["flag"] == "--inflow-bdy" && append!(extra, ["--inflow-bdy", new])
            else
                remap[a["argv_value"]] = new
            end
        end
        over = Set(first.(_argv_items(overrides)))
        out = String[]
        for (flag, vals) in _argv_items(argv)
            (flag in _DROP_FLAGS || flag == "--vis" || flag in over) && continue
            if flag == "--inflow-point" && !isempty(vals)
                parts = split(vals[1], ","; limit = 4)
                length(parts) >= 3 && (parts[3] = get(remap, strip(parts[3]), parts[3]))
                vals = [join(parts, ",")]
            elseif !isempty(vals)
                vals = [get(remap, vals[1], vals[1])]
            end
            push!(out, flag); append!(out, vals)
        end
        # Defaults change between versions (0.5.0: --flow-model sgs → standard):
        # pin the recorded flow model and scheme when the command line left them
        # to the default, so that the rerun is the recorded run.
        given = Set(first.(_argv_items(argv)))
        opts = get(info, "options", Dict{String,Any}())
        for (flag, val) in (("--flow-model", get(info, "flow_model", get(opts, "flow_model", nothing))),
                            ("--scheme", get(opts, "scheme", nothing)))
            (val === nothing || isempty(string(val)) || flag in given || flag in over) && continue
            append!(out, [flag, string(val)])
        end
        append!(out, extra)
        append!(out, ["--meshload", mesh])
        "--output" in over || append!(out, ["--output", output])
        append!(out, overrides)
        return out
    end
end

end # module
