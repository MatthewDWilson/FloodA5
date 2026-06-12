# boundaryinputs/timeseries_io.jl
# --------------------------------
# Time-series file I/O for boundary condition inputs.
#
# Supported formats
# -----------------
#   LisfloodBDYReader  — LISFLOOD-FP .bdy file (one or more named series)
#   TwoColumnCSVReader — plain two-column CSV (t_s, Q_m3s); auto-detects header
#
# Reader abstraction
# ------------------
# All readers implement:
#   read_timeseries(reader) → Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
# where the Dict maps series name → (t_seconds, Q_m3s).
# For single-series files (CSV, single-series .bdy), the dict has one entry
# keyed by the filename stem.
#
# BCI file parsing
# ----------------
# parse_bci_file(path; epsg=nothing) → Vector{BCIEntry}
# Parses a LISFLOOD-FP .bci file and returns structured entries.
# The companion .bdy file (same directory, referenced by QVAR entries) is
# located and read automatically.
#
# Future additions
# ----------------
# WaterML2Reader <: AbstractTimeSeriesReader
# CFNetCDFGaugeReader <: AbstractTimeSeriesReader
# HydroMLReader <: AbstractTimeSeriesReader
# All follow the same read_timeseries interface.

abstract type AbstractTimeSeriesReader end

# ---------------------------------------------------------------------------
# .bdy file reader
# ---------------------------------------------------------------------------

"""
    LisfloodBDYReader

Reader for LISFLOOD-FP .bdy boundary data files.

File format (one or more named series per file):

    SERIES_NAME
    N                     # number of time steps
    seconds               # time unit: seconds | hours | days
    0.0    0.0
    3600   5.2
    ...
    NEXT_SERIES_NAME
    ...

Comments (lines beginning # or !) and blank lines are ignored.
Multiple series in a single file are all parsed and returned.
"""
struct LisfloodBDYReader <: AbstractTimeSeriesReader
    path :: String
end

function read_timeseries(r::LisfloodBDYReader
        )::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}

    isfile(r.path) || error("BDY file not found: $(r.path)")

    result = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()

    # Split on whitespace (handles both spaces and tabs).
    # Filter blank and comment lines.
    lines = filter(l -> !isempty(strip(l)) &&
                        !startswith(strip(l), '#') &&
                        !startswith(strip(l), '!'),
                   readlines(r.path))

    valid_units = Set([
        "second", "seconds",
        "minute", "minutes", "mins",
        "hour",   "hours",
        "day",    "days"
    ])

    unit_mult = Dict(
        "second"  => 1.0,  "seconds" => 1.0,
        "minute"  => 60.0, "minutes" => 60.0, "mins" => 60.0,
        "hour"    => 3600.0, "hours" => 3600.0,
        "day"     => 86400.0, "days" => 86400.0,
    )

    # Is this line purely a time-unit keyword?
    _is_unit_line(l::String) = lowercase(strip(l)) in valid_units

    # Is this line a plain integer (step count)?
    _is_int_line(l::String) = tryparse(Int, strip(l)) !== nothing

    # Is this line a two-column data row (both numeric)?
    _is_data_row(l::String) = begin
        p = split(strip(l))
        length(p) >= 2 &&
            tryparse(Float64, p[1]) !== nothing &&
            tryparse(Float64, p[2]) !== nothing
    end

    # Skip a file-level header line: a line that is not a data row,
    # not a unit keyword, and not an integer (i.e. plaintext descriptor).
    i = 1
    if i <= length(lines)
        l = strip(lines[i])
        if !_is_data_row(lines[i]) && !_is_unit_line(lines[i]) && !_is_int_line(lines[i])
            # This looks like a text descriptor — skip it if it contains no digits
            # at the start (i.e. it's not a series name that happens to start with text).
            # We always skip the first non-numeric, non-unit line as a file header.
            @info "BDY parse: file header: '$l'"
            i += 1
        end
    end

    while i <= length(lines)

        # ── Series name ───────────────────────────────────────────────────
        # The series name is a plain text label — not a data row, not an int,
        # not a unit keyword.
        name = strip(lines[i])
        i += 1

        if i > length(lines)
            @warn "BDY parse: unexpected end-of-file after series name '$name'."
            break
        end

        # ── Step count and time unit — on the same line, whitespace-delimited ─
        # Format: "241   hours"  or  "241\thours"
        step_parts = split(strip(lines[i]))
        if length(step_parts) < 2
            @warn "BDY parse: expected '<n_steps> <unit>' after series name '$name', " *
                  "got '$(strip(lines[i]))'. Skipping series."
            i += 1
            continue
        end

        n_steps = tryparse(Int, step_parts[1])
        if n_steps === nothing
            @warn "BDY parse: invalid step count '$(step_parts[1])' for series '$name'. " *
                  "Skipping series."
            i += 1
            continue
        end

        unit_str = lowercase(step_parts[2])
        if !(unit_str in valid_units)
            @warn "BDY parse: unrecognised time unit '$unit_str' for series '$name'. " *
                  "Defaulting to seconds."
        end
        t_mult = get(unit_mult, unit_str, 1.0)
        i += 1

        if i > length(lines)
            @warn "BDY parse: unexpected end-of-file before data for series '$name'."
            break
        end

        # ── Data rows ─────────────────────────────────────────────────────
        t_vec  = Vector{Float64}(undef, n_steps)
        Q_vec  = Vector{Float64}(undef, n_steps)
        n_read = 0

        while n_read < n_steps && i <= length(lines)
            parts = split(strip(lines[i]))

            if length(parts) >= 2
                q_val = tryparse(Float64, parts[1])
                t_val = tryparse(Float64, parts[2])
                
                if t_val !== nothing && q_val !== nothing
                    n_read += 1
                    t_vec[n_read] = t_val * t_mult
                    Q_vec[n_read] = q_val
                    length(parts) > 2 &&
                        @debug "BDY parse: extra columns ignored in series '$name' row $n_read."
                else
                    # Non-numeric line — likely start of next series name; stop.
                    break
                end
            else
                # Single-column or empty — likely next series name; stop.
                break
            end
            i += 1
        end

        # ── Handle truncated series ────────────────────────────────────────
        if n_read < n_steps
            @warn "BDY parse: series '$name' expected $n_steps rows, got $n_read. Truncating."
            resize!(t_vec, n_read)
            resize!(Q_vec, n_read)
        end

        n_read == 0 && continue

        # ── Validate monotonicity ─────────────────────────────────────────
        # Use a tolerance to absorb floating-point rounding from t_mult multiplication
        tol = t_mult * 1e-9
        for k in 2:n_read
            if !isfinite(t_vec[k])
                @warn "BDY parse: series '$name' contains non-finite time at step $k."
                break
            end
            if t_vec[k] < t_vec[k-1] - tol
                @warn "BDY parse: series '$name' non-monotonic time at step $k " *
                      "(t=$(t_vec[k]) ≤ t=$(t_vec[k-1])). Results may be incorrect."
                break
            end
        end

        # ── Duplicate name protection ─────────────────────────────────────
        haskey(result, name) &&
            @warn "BDY parse: duplicate series name '$name'. Overwriting previous series."

        result[name] = (t_vec, Q_vec)
    end

    isempty(result) && error("BDY file '$(r.path)' contained no valid series.")
    return result
end


# ---------------------------------------------------------------------------
# Two-column CSV reader
# ---------------------------------------------------------------------------

"""
    TwoColumnCSVReader

Reader for plain two-column CSV files (t_s, Q_m3s).

Auto-detects a header row (skipped if the first data item is non-numeric).
Returns a single series keyed by the filename stem.
Also handles single-series LISFLOOD-FP .bdy files transparently — the series
name from the .bdy is used as the key.
"""
struct TwoColumnCSVReader <: AbstractTimeSeriesReader
    path  :: String
    t_col :: Int   # column index for time (default 1)
    Q_col :: Int   # column index for discharge (default 2)
end
TwoColumnCSVReader(path::String) = TwoColumnCSVReader(path, 1, 2)

function read_timeseries(r::TwoColumnCSVReader
        )::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
    isfile(r.path) || error("File not found: $(r.path)")

    # Delegate .bdy files to the BDY reader; take first (or only) series
    if lowercase(last(splitext(r.path))) == ".bdy"
        all = read_timeseries(LisfloodBDYReader(r.path))
        isempty(all) && error("No series in BDY file: $(r.path)")
        if length(all) > 1
            @warn "CSV reader applied to multi-series .bdy file '$(r.path)': " *
                  "using first series '$(first(keys(all)))'. " *
                  "Use LisfloodBDYReader to access all series."
        end
        key = first(keys(all))
        return Dict(key => all[key])
    end

    label = splitext(basename(r.path))[1]
    lines = filter(l -> !isempty(strip(l)) &&
                        !startswith(strip(l), '#') &&
                        !startswith(strip(l), '!'),
                   readlines(r.path))

    isempty(lines) && error("File is empty: $(r.path)")

    t_vec = Float64[]
    Q_vec = Float64[]

    for (li, line) in enumerate(lines)
        parts = split(strip(line), r"[,\s]+")
        length(parts) < max(r.t_col, r.Q_col) && continue
        t_val = tryparse(Float64, parts[r.t_col])
        Q_val = tryparse(Float64, parts[r.Q_col])
        if t_val === nothing || Q_val === nothing
            li == 1 && continue   # skip header row silently
            @warn "CSV parse: skipping non-numeric row $li in '$(r.path)': $line"
            continue
        end
        push!(t_vec, t_val)
        push!(Q_vec, Q_val)
    end

    isempty(t_vec) && error("No numeric data rows in '$(r.path)'.")

    return Dict(label => (t_vec, Q_vec))
end

# ---------------------------------------------------------------------------
# Auto-detect reader
# ---------------------------------------------------------------------------

"""
    auto_timeseries_reader(path) → AbstractTimeSeriesReader

Return the appropriate reader for `path` based on file extension:
  .bdy  → LisfloodBDYReader
  other → TwoColumnCSVReader (handles .csv, .txt, .dat, etc.)
"""
function auto_timeseries_reader(path::String)::AbstractTimeSeriesReader
    ext = lowercase(last(splitext(path)))
    ext == ".bdy" ? LisfloodBDYReader(path) : TwoColumnCSVReader(path)
end

# ---------------------------------------------------------------------------
# BCI file parser
# ---------------------------------------------------------------------------
# BCIEntry is defined in boundary_conditions.jl (included before this file)
# so that boundary_conditions.jl can use it without a forward dependency.

"""
    parse_bci_file(path; bdy_path=nothing) → (Vector{BCIEntry}, Dict series)

Parse a LISFLOOD-FP .bci file and return:
  1. A vector of BCIEntry structs (one per valid non-comment line)
  2. A Dict of time series loaded from the companion .bdy file (for QVAR entries)

`bdy_path` overrides the automatic .bdy search (replaces .bci extension with .bdy
in the same directory).  Pass `nothing` to use automatic search.

Coordinate system: coordinates are returned as-is from the file.  The caller
is responsible for CRS conversion (see `_convert_bci_coords` in FloodModel.jl).

Boundary types N, E, S, W and F are returned in the entry vector but flagged
as unsupported — the caller logs a helpful message and skips them.
See PROJECT_STATE.md for planned implementation.
"""
function parse_bci_file(path     :: String;
                         bdy_path :: Union{String,Nothing} = nothing
                         )::Tuple{Vector{BCIEntry},
                                  Dict{String,Tuple{Vector{Float64},Vector{Float64}}}}
    isfile(path) || error("BCI file not found: $path")

    # Locate companion .bdy file — only required if QVAR entries are present.
    # Defer the error to the QVAR-handling block below so QFIX-only .bci files
    # work without a companion .bdy.
    auto_bdy = replace(path, r"\.[bB][cC][iI]$" => ".bdy")
    if !isfile(auto_bdy)
        auto_bdy = replace(path, r"\.[bB][cC][iI]$" => ".BDY")
    end
    bdy_file = bdy_path !== nothing ? bdy_path :
               isfile(auto_bdy)     ? auto_bdy : nothing

    entries = BCIEntry[]
    for (lineno, raw) in enumerate(readlines(path))
        line = strip(raw)
        isempty(line)           && continue
        startswith(line, '#')   && continue
        startswith(line, '!')   && continue
        startswith(line, "//")  && continue

        parts = split(line)
        length(parts) < 4 && begin
            @warn "BCI line $lineno: expected ≥4 columns, got $(length(parts)). Skipping."
            continue
        end

        btype_str = uppercase(strip(parts[1]))
        length(btype_str) != 1 && begin
            @warn "BCI line $lineno: boundary type must be a single character, " *
                  "got '$btype_str'. Skipping."
            continue
        end
        btype = btype_str[1]
        btype ∉ ('N','E','S','W','P','F') && begin
            @warn "BCI line $lineno: unknown boundary type '$btype'. Skipping."
            continue
        end

        x1 = tryparse(Float64, parts[2])
        y1 = tryparse(Float64, parts[3])
        (x1 === nothing || y1 === nothing) && begin
            @warn "BCI line $lineno: non-numeric coordinates '$(parts[2])' '$(parts[3])'. Skipping."
            continue
        end

        bc_code  = length(parts) >= 4 ? uppercase(strip(parts[4])) : "FREE"
        bc_value = length(parts) >= 5 ? strip(parts[5]) : ""

        # For P type: x1=lon, y1=lat, x2=x1, y2=y1 (point, not segment)
        push!(entries, BCIEntry(btype, x1, y1, x1, y1, bc_code, bc_value))
    end

    isempty(entries) && @warn "BCI file '$path' contained no valid entries."

    # Load .bdy for QVAR entries
    series_dict = Dict{String, Tuple{Vector{Float64},Vector{Float64}}}()
    qvar_names  = Set{String}(e.bc_value for e in entries if e.bc_code == "QVAR" && !isempty(e.bc_value))

    if !isempty(qvar_names)
        bdy_file === nothing && error(
            "BCI file references QVAR series $(collect(qvar_names)) but no " *
            "companion .bdy file was found alongside '$path'. " *
            "Pass bdy_path= explicitly or place the .bdy file in the same directory.")
        all_series = read_timeseries(LisfloodBDYReader(bdy_file))
        for name in qvar_names
            haskey(all_series, name) || @warn "BCI references series '$name' not found in '$bdy_file'."
        end
        merge!(series_dict, all_series)
    end

    return entries, series_dict
end
