"""
setup.jl
--------
One-time setup script for FloodA5.

Run this from the project root before launching FloodModel.jl for the first time,
or after a fresh Julia / CUDA install:

    julia setup.jl

What it does
------------
1. Activates the project environment and instantiates all Julia packages.
2. Installs required Python packages into PyCall's bundled Conda Python
   (avoids system Python conflicts).
3. Verifies the a5_bridge.py subprocess pathway.
4. Checks for a functional CUDA GPU; warns and confirms CPU fallback if absent.
5. Applies the Julia 1.12 Printf compatibility fix in FloodModel.jl if needed
   (replaces @sprintf with string interpolation; removes unused Printf import).
"""

# ---------------------------------------------------------------------------
# 0. Activate project environment
# ---------------------------------------------------------------------------

using Pkg

println("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
println("  FloodA5 — setup")
println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

project_dir = @__DIR__
println("[ 0/4 ] Activating project at: $project_dir")
Pkg.activate(project_dir)

# ---------------------------------------------------------------------------
# 1. Julia packages
# ---------------------------------------------------------------------------

println("\n[ 1/4 ] Installing / updating Julia packages …")

# Required
required = [
    "PyCall", "JSON3", "ArchGDAL", "HDF5",
    "Oxygen", "HTTP",
]

# Optional (warn but don't abort if they fail)
optional = Dict(
    "CUDA"           => "GPU-accelerated PIP sampling",
    "GLMakie"        => "--vis makie visualisation backend",
    "BenchmarkTools" => "benchmarking utilities",
)

# Instantiate from Project.toml/Manifest.toml when present
if isfile(joinpath(project_dir, "Project.toml"))
    println("  Found Project.toml — running Pkg.instantiate()")
    Pkg.instantiate()
else
    println("  No Project.toml found — adding packages individually")
    Pkg.add(required)
end

# Attempt optional packages
for (pkg, purpose) in optional
    try
        Pkg.add(pkg)
        println("  ✓ $pkg  ($purpose)")
    catch e
        @warn "  Could not install $pkg ($purpose): $e"
    end
end

println("  ✓ Julia packages ready")

# ---------------------------------------------------------------------------
# 2. Python packages (into PyCall's Conda Python — no system conflicts)
# ---------------------------------------------------------------------------

println("\n[ 2/4 ] Installing Python packages into PyCall's Python …")

using PyCall

py_packages = ["pya5", "geopandas", "pyarrow", "shapely", "numpy"]

python_exe = PyCall.python
println("  PyCall Python: $python_exe")

for pkg in py_packages
    print("  Installing $pkg … ")
    result = run(ignorestatus(Cmd([python_exe, "-m", "pip", "install", "--quiet", pkg])))
    if result.exitcode == 0
        println("✓")
    else
        println()
        @warn "  pip install $pkg returned exit code $(result.exitcode) — check output above"
    end
end

# Verify pya5 is importable
print("  Verifying pya5 import … ")
try
    pyimport("a5")
    println("✓")
catch e
    @error "pya5 could not be imported after installation: $e"
    println("""
  Try manually:
      $(python_exe) -m pip install pya5
  then re-run this script.
""")
end

# ---------------------------------------------------------------------------
# 3. Python bridge smoke-test
# ---------------------------------------------------------------------------

println("\n[ 3/4 ] Verifying Python bridge (a5_bridge.py) …")

bridge = joinpath(project_dir, "mesh", "a5_bridge.py")

if !isfile(bridge)
    @error "a5_bridge.py not found at $bridge — expected at mesh/a5_bridge.py"
else
    buf_out = IOBuffer()
    buf_err = IOBuffer()
    result  = run(pipeline(
        ignorestatus(setenv(Cmd([python_exe, bridge, "check"]), dir=joinpath(project_dir, "mesh"))),
        stdout=buf_out, stderr=buf_err
    ))
    out = String(take!(buf_out))
    err = String(take!(buf_err))

    if result.exitcode == 0
        println("  ✓ Bridge check passed")
        for line in filter(!isempty, split(out, "\n"))
            println("    $line")
        end
    else
        @warn "  Bridge check exited with code $(result.exitcode)"
        isempty(out) || println("  STDOUT:\n$out")
        isempty(err) || println("  STDERR:\n$err")
        println("""
  If pya5 installed correctly above, this may be a PATH issue.
  Manual check: $(python_exe) $(bridge) check
""")
    end
end

# ---------------------------------------------------------------------------
# 4. CUDA check
# ---------------------------------------------------------------------------

println("\n[ 4/4 ] Checking CUDA / GPU …")

# Use `global` so assignments inside try/catch are visible outside the block.
global cuda_ok = false
try
    @eval using CUDA
    # Query functional() — may warn about toolkit/driver version mismatch but
    # still return true when the GPU works despite the version difference.
    global cuda_ok = try
        @eval CUDA.functional()
    catch cuda_err
        @warn "  CUDA.jl initialisation error: $cuda_err"
        false
    end
    if cuda_ok
        # Wrap device query — can fail if runtime/driver are slightly out of sync.
        gpu_name = try
            dev = @eval CUDA.device()
            @eval CUDA.name($dev)
        catch
            "unknown GPU"
        end
        println("  ✓ CUDA functional — GPU: $gpu_name")
        println("    PIP sampling will run on GPU (fast path).")
    else
        @warn "  CUDA.jl loaded but no functional GPU detected."
        println("    FloodA5 will fall back to multi-threaded CPU for PIP sampling.")
        println("    If you have a CUDA GPU, check your driver with: nvidia-smi")
    end
catch e
    @warn "  CUDA.jl not available or failed to load: $e"
    println("    FloodA5 will use multi-threaded CPU for PIP sampling (normal on systems without a CUDA GPU).")
    println("    To enable GPU support: Pkg.add(\"CUDA\") and ensure CUDA toolkit drivers are installed.")
end

# ---------------------------------------------------------------------------
# 5. @sprintf / Printf compatibility fix (Julia 1.12+)
# ---------------------------------------------------------------------------

println("\n[ * ] Checking FloodModel.jl for Julia 1.12 Printf compatibility …")

flood_model = joinpath(project_dir, "FloodModel.jl")

if !isfile(flood_model)
    @warn "FloodModel.jl not found — skipping Printf check"
else
    src = read(flood_model, String)

    # Julia 1.12 changed @sprintf to require a compile-time string literal as its
    # first argument. When DocStringExtensions is active, macro expansion during
    # docstring lowering causes even correctly-scoped @sprintf calls to fail.
    # The robust fix is to replace all @sprintf/@printf with Julia string
    # interpolation, and remove the now-unused `using Printf` import.

    # Each new_pat is a raw string to be written into FloodModel.jl.
    # Dollar signs must be escaped (\$) so Julia does not interpolate them
    # here in setup.jl — they should only be evaluated inside FloodModel.jl.
    patches = [
        # (old_pattern, new_pattern, description)
        (
            """    @info @sprintf(\"Adjacency built: %d cells, %d undirected edges (shared-vertex method)\",
                   n, n_edges)""",
            """    @info \"Adjacency built: \$(n) cells, \$(n_edges) undirected edges (shared-vertex method)\" """,
            "Adjacency built message"
        ),
        (
            """    @info @sprintf(\"Edge list built: %d edges for %d cells (%.2f edges/cell)\",
                   n_edges, n, n_edges / max(n, 1))""",
            """    @info \"Edge list built: \$(n_edges) edges for \$(n) cells (\$(round(n_edges/max(n,1),digits=2)) edges/cell)\" """,
            "Edge list built message"
        ),
        (
            """        @info @sprintf(\"Edge non-orthogonality (cos θ):  min=%.3f  mean=%.3f  max=%.3f\",
                       minimum(valid_ct), mean(valid_ct), maximum(valid_ct))""",
            """        @info \"Edge non-orthogonality (cos θ):  min=\$(round(minimum(valid_ct),digits=3))  mean=\$(round(mean(valid_ct),digits=3))  max=\$(round(maximum(valid_ct),digits=3))\" """,
            "Edge non-orthogonality message"
        ),
        (
            """            @info @sprintf(
                \"  step=%5d  t=%.1fs  dt=%.2fs  wet=%d  max_depth=%.3fm  \" *
                \"domain_vol=%.1fm³  mb_err=%.1fm³\",
                step, t, dt, n_wet, max_depth, domain_vol, mb_err)""",
            """            @info \"  step=\$(lpad(step,5))  t=\$(round(t,digits=1))s  dt=\$(round(dt,digits=2))s  wet=\$(n_wet)  max_depth=\$(round(max_depth,digits=3))m  domain_vol=\$(round(domain_vol,digits=1))m³  mb_err=\$(round(mb_err,digits=1))m³\" """,
            "Simulation step progress message"
        ),
        (
            """        @info @sprintf(\"Injection point: (%.5f, %.5f) → cell %s  \" *
                       \"(dist=%.0fm)  rate=%.4f m³/s\",
                       lon, lat, cid, dist_m, rate)""",
            """        @info \"Injection point: (\$(round(lon,digits=5)), \$(round(lat,digits=5))) → cell \$(cid)  (dist=\$(round(dist_m,digits=0))m)  rate=\$(round(rate,digits=4)) m³/s\" """,
            "Injection point message"
        ),
    ]

    local patched = src
    local n_fixes = 0

    for (old_pat, new_pat, desc) in patches
        if occursin(old_pat, patched)
            patched = replace(patched, old_pat => new_pat)
            n_fixes += 1
            println("    fixed: $desc")
        end
    end

    # Remove `using Printf` if no Printf macros remain
    if !occursin("@sprintf", patched) && !occursin("@printf", patched) && !occursin("Printf.", patched)
        patched = replace(patched, "using Printf\n" => "")
        println("    removed: unused `using Printf`")
    end

    if n_fixes > 0
        backup = flood_model * ".pre_setup_backup"
        cp(flood_model, backup; force=true)
        write(flood_model, patched)
        println("  ✓ Applied $n_fixes Printf fix(es) in FloodModel.jl")
        println("    Original backed up to: $(basename(backup))")
    else
        println("  ✓ No Printf fix needed (already patched or pattern not found)")
    end
end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

println("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Setup complete.

  Quick-start:
    julia --threads auto FloodModel.jl \\
        --meshgen examples/example_aoi.geojson --meshres 14 --vis

  GPU acceleration: $(cuda_ok ? "✓ enabled" : "✗ not available — CPU fallback active")
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")
