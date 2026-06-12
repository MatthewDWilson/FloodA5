# boundaryinputs/boundary_conditions.jl
# ---------------------------------------
# Open/closed outflow boundary conditions for domain-edge cells.
#
# Design
# ------
# Boundary cells are mesh cells with fewer than N_SIDES (5) edge-sharing
# neighbours.  Their "missing" edges face the exterior of the domain.
# Currently these are implicit closed walls.  This module makes the BC
# type explicit per ghost edge, defaulting to ZeroGradient (open outflow).
#
# Ghost edges are pre-computed at initialise_flow_model time and stored in
# FlowState.  Each step, Phase E (boundary outflow) applies a flux kernel
# across each non-Closed ghost edge and accumulates vol_removed.
#
# BCType hierarchy
# ----------------
#   Closed       — no flux (legacy implicit, enabled by --closed-boundaries)
#   ZeroGradient — ghost WSE = boundary cell WSE; transmissive/non-reflective
#   Critical     — Q = A × √(g × h³); free outfall (future, stub ready)
#   FixedWSE     — fixed water surface (tide); requires time series (future)
#   FixedQ       — prescribed outflow discharge (future)
#
# Unsupported .bci boundary types
# --------------------------------
# N/E/S/W (rectangular edge boundaries): these require an axis-aligned
#   rectangular domain.  FloodA5 supports arbitrary polygon domains, so
#   cardinal-direction edge BCs have no natural mapping.  If encountered in
#   a .bci file, a helpful message is logged directing the user to specify
#   boundaries using a GeoJSON --bc-file instead.  See PROJECT_STATE.md.
#
# F (sub-grid channel internal free boundary): LISFLOOD-FP uses this for
#   internal channel boundaries in the SGC solver.  FloodA5 has no direct
#   equivalent in the current architecture.  Parsed without error; logged as
#   unsupported.  See PROJECT_STATE.md for planned reconsideration.
#
# GeoJSON boundary file
# ---------------------
# --bc-file GEOJSON  specifies BC types per geographic segment (LineString or
# Polygon).  Each feature's "bc_type" property overrides the default for
# boundary cells whose centres fall within 1.5× the cell diameter of the
# feature geometry.

# ---------------------------------------------------------------------------
# BCType enum
# ---------------------------------------------------------------------------

"""
    BCType

Boundary condition type for domain-edge (ghost) edges.

  Closed       — no flux through this boundary edge.
                 Legacy implicit behaviour; enabled by --closed-boundaries.
  ZeroGradient — ghost cell WSE = boundary cell WSE.
                 Transmissive / non-reflective outflow.
                 New default.
  Critical     — critical-depth outflow: Q = A × √(g h³).
                 Planned; enum value reserved.
  FixedWSE     — prescribed water surface elevation (tide / downstream stage).
                 Planned (tide session); requires a time series input.
  FixedQ       — prescribed outflow discharge.
                 Planned.
"""
# ---------------------------------------------------------------------------
# BCType and GhostEdge are defined in FloodModel.jl (before FlowState) so
# that FlowState can reference them in its field declarations.
# This file uses them but does not redefine them.
# ---------------------------------------------------------------------------

function BCType_from_string(s::AbstractString)::BCType
    u = uppercase(strip(s))
    u == "CLOSED"        && return Closed
    u == "ZEROGRADIENT"  && return ZeroGradient
    u == "ZERO_GRADIENT" && return ZeroGradient
    u == "FREE"          && return ZeroGradient   # .bci FREE → ZeroGradient
    u == "CRITICAL"      && return Critical
    u == "FIXEDWSE"      && return FixedWSE
    u == "FIXED_WSE"     && return FixedWSE
    u == "HFIX"          && return FixedWSE       # LISFLOOD-FP code
    u == "FIXEDQ"        && return FixedQ
    u == "FIXED_Q"       && return FixedQ
    error("Unknown BCType string: '$s'. " *
          "Valid values: Closed, ZeroGradient, Critical, FixedWSE, FixedQ")
end

# ---------------------------------------------------------------------------
# GhostEdge — one virtual edge per missing neighbour slot on a boundary cell
# ---------------------------------------------------------------------------

# GhostEdge is defined in FloodModel.jl — see preamble above EdgeList.

# ---------------------------------------------------------------------------
# BoundarySegment — associates a BCType with a set of boundary cells
# ---------------------------------------------------------------------------

"""
    BoundarySegment

Associates a BC type with a set of boundary cell indices, typically loaded
from a GeoJSON --bc-file or from .bci FREE entries.

`wse_series` is used for FixedWSE (tide) BCs; it holds (t_seconds, wse_m).
It is `nothing` for all currently implemented BC types.
"""
struct BoundarySegment
    cell_indices :: Vector{Int}
    bc_type      :: BCType
    label        :: String
    wse_series   :: Union{Tuple{Vector{Float64},Vector{Float64}}, Nothing}
end
BoundarySegment(idxs, bc, lbl) = BoundarySegment(idxs, bc, lbl, nothing)

# ---------------------------------------------------------------------------
# Ghost-cell WSE selector
# ---------------------------------------------------------------------------

"""
    _ghost_wse(wse_ci, sill, bc, wse_nb, L_ci_nb, L_ghost) → Float64

Water surface elevation to assign to the virtual ghost cell.

  ZeroGradient — linear extrapolation of the interior WSE gradient:
                   wse_ghost = wse_ci + (wse_ci - wse_nb) * L_ghost / L_ci_nb
                 This gives a non-zero dWSE across the ghost edge whenever the
                 interior is flowing toward the boundary, producing immediate
                 outflow from the first timestep without needing accumulated
                 momentum.  Clamped so wse_ghost ≥ sill (no negative depth).
                 Falls back to wse_ci (true zero-gradient) when no interior
                 neighbour is available (wse_nb = NaN or L_ci_nb ≤ 0).

  Critical      — critical-depth: wse_ghost = sill + (2/3)·max(0, wse_ci - sill)

  Closed        — returns -Inf (sentinel; caller skips the edge)
"""
@inline function _ghost_wse(wse_ci   :: Float64,
                              sill     :: Float64,
                              bc       :: BCType,
                              wse_nb   :: Float64 = NaN,
                              L_ci_nb  :: Float64 = 0.0,
                              L_ghost  :: Float64 = 0.0)::Float64
    bc === Closed   && return -Inf
    bc === Critical && return sill + (2.0/3.0) * max(0.0, wse_ci - sill)

    if bc === ZeroGradient
        # Gradient extrapolation when an interior neighbour is available
        if !isnan(wse_nb) && L_ci_nb > 1.0 && L_ghost > 0.0
            slope      = (wse_ci - wse_nb) / L_ci_nb
            wse_extrap = wse_ci + slope * L_ghost
            # Clamp: ghost WSE must be ≥ sill (no below-bed ghost cell)
            # but we do NOT clamp wse_extrap < wse_ci — if the interior slope
            # is falling (wse_nb > wse_ci), the extrapolated ghost WSE will be
            # below wse_ci, which correctly reduces the flux for decelerating flow.
            return max(wse_extrap, sill)
        else
            # Fallback: true zero-gradient (flat extrapolation)
            return wse_ci
        end
    end

    # FixedWSE / FixedQ: caller must inject the correct value.
    return wse_ci   # safe fallback
end

# ---------------------------------------------------------------------------
# Ghost edge geometry: build from cell polygon
# ---------------------------------------------------------------------------

"""
    _ghost_edge_sides(boundary) → Vector{Tuple{Int,Int}}

Return the (v1_idx, v2_idx) index pairs for each side of a pentagon polygon.
Pentagon boundaries are stored as 5 vertices + closing repeat (6 entries).
Returns 5 side pairs (closing side connects vertex 5 back to vertex 1).
"""
function _ghost_edge_sides(boundary::Vector{Vector{Float64}})
    nv = length(boundary)
    # Remove closing vertex if present (same as first)
    n  = (nv >= 2 &&
          boundary[nv][1] ≈ boundary[1][1] &&
          boundary[nv][2] ≈ boundary[1][2]) ? nv - 1 : nv
    return [(k, mod1(k+1, n)) for k in 1:n]
end

"""
    _build_ghost_edges(cells, adj, id_idx, edges, elevations, sgs_tables,
                        n_sides, default_bc) → (BitVector, Vector{GhostEdge}, Vector{BCType})

Identify all domain-edge (boundary) cells and pre-compute a GhostEdge for
each missing neighbour slot.

Returns:
  boundary_mask  :: BitVector         — true for each domain-edge cell
  ghost_edges    :: Vector{GhostEdge} — flat list of all ghost edges
  ghost_cell_bc  :: Vector{BCType}    — BC type per ghost edge (default_bc for all)

Ghost edge geometry:
  width — haversine length of the actual missing polygon side, identified by
          checking which cell polygon sides are not shared with any mesh neighbour.
          Falls back to mean of existing edge widths if polygon analysis fails.
  L     — haversine distance from cell centre to the midpoint of the ghost edge,
          doubled (to approximate ghost cell centre-to-cell-centre distance).
  sill  — elevation[ci] for standard flow; sgs_tables[ci].z_min for SGS.

Called once from initialise_flow_model after _build_edge_list.
"""
function _build_ghost_edges(cells       :: Vector{A5Grid.A5Cell},
                              adj         :: Dict{String,Vector{String}},
                              id_idx      :: Dict{String,Int},
                              edges_list  :: EdgeList,
                              elevations  :: Vector{Float64},
                              sgs_tables  :: AbstractVector,
                              n_sides     :: Int,
                              default_bc  :: BCType
                              )::Tuple{BitVector, Vector{GhostEdge}, Vector{BCType}}
    n             = length(cells)
    use_sgs       = !isempty(sgs_tables)
    _norm(id)     = A5Grid._to_hex(parse(UInt64, id, base=16))

    # Count actual neighbours per cell from EdgeList
    n_neighbours  = zeros(Int, n)
    for e in 1:edges_list.n_edges
        n_neighbours[edges_list.cell_i[e]] += 1
        n_neighbours[edges_list.cell_j[e]] += 1
    end

    boundary_mask = BitVector(n_neighbours[i] < n_sides for i in 1:n)
    n_boundary    = count(boundary_mask)
    n_boundary == 0 && @info "No boundary cells found (fully enclosed domain)."
    n_boundary  > 0 && @info "Boundary cells: $n_boundary / $n ($(round(100n_boundary/n, digits=1))%)"

    ghost_edges   = GhostEdge[]
    ghost_cell_bc = BCType[]

    # For each boundary cell, find which polygon sides are not shared with
    # any mesh neighbour and create a GhostEdge for each such side.
    for i in 1:n
        boundary_mask[i] || continue
        cell  = cells[i]
        bnd   = cell.boundary

        # Collect shared vertex pairs from all actual neighbours
        nid     = _norm(cell.id)
        nbr_ids = get(adj, nid, String[])
        shared_pairs = Set{Tuple{Float64,Float64,Float64,Float64}}()
        for nb_id in nbr_ids
            j = get(id_idx, nb_id, 0)
            j == 0 && continue
            edge = A5Grid._shared_edge(bnd, cells[j].boundary)
            edge === nothing && continue
            lon1, lat1, lon2, lat2 = edge
            # Store canonically ordered pair (smaller lon first) for comparison
            key = lon1 <= lon2 ? (lon1, lat1, lon2, lat2) : (lon2, lat2, lon1, lat1)
            push!(shared_pairs, key)
        end

        # Pentagon sides: check each side against shared_pairs
        sides      = _ghost_edge_sides(bnd)
        nv         = length(bnd)
        n_actual   = length(bnd) - (nv >= 2 &&
                                    bnd[nv][1] ≈ bnd[1][1] &&
                                    bnd[nv][2] ≈ bnd[1][2] ? 1 : 0)

        # Mean edge width fallback (in case polygon analysis misidentifies sides)
        existing_widths = Float64[]
        for e in 1:edges_list.n_edges
            if edges_list.cell_i[e] == i || edges_list.cell_j[e] == i
                push!(existing_widths, edges_list.width[e])
            end
        end
        mean_width = isempty(existing_widths) ? sqrt(A5Grid._polygon_area_m2(bnd)) :
                                                sum(existing_widths) / length(existing_widths)

        cell_lon = cell.center_lon
        cell_lat = cell.center_lat

        for (v1_idx, v2_idx) in sides
            v1 = bnd[v1_idx]
            v2 = bnd[v2_idx]

            # Check whether this side is shared with a neighbour
            key_fwd = v1[1] <= v2[1] ? (v1[1], v1[2], v2[1], v2[2]) :
                                        (v2[1], v2[2], v1[1], v1[2])
            is_shared = false
            for sp in shared_pairs
                if (abs(sp[1] - key_fwd[1]) < 1e-8 && abs(sp[2] - key_fwd[2]) < 1e-8 &&
                    abs(sp[3] - key_fwd[3]) < 1e-8 && abs(sp[4] - key_fwd[4]) < 1e-8)
                    is_shared = true
                    break
                end
            end
            is_shared && continue

            # This side is a ghost edge — compute geometry
            width = A5Grid._haversine_m(v1[1], v1[2], v2[1], v2[2])
            isnan(width) || width < 1.0 && (width = mean_width)

            # Ghost cell centre approximation: 2× distance from cell centre to edge midpoint
            mid_lon = 0.5 * (v1[1] + v2[1])
            mid_lat = 0.5 * (v1[2] + v2[2])
            d_to_mid = A5Grid._haversine_m(cell_lon, cell_lat, mid_lon, mid_lat)
            L = 2.0 * max(d_to_mid, width)

            # Sill elevation
            elev_ci = isnan(elevations[i]) ? 0.0 : elevations[i]
            sill = if use_sgs && i <= length(sgs_tables) && sgs_tables[i] !== nothing
                tbl = sgs_tables[i]
                isnan(tbl.z_min) ? elev_ci : tbl.z_min
            else
                elev_ci
            end

            # Find the interior neighbour most aligned with the outward normal of
            # this ghost edge.  The outward normal points from cell_centre toward
            # edge midpoint.  We pick the real neighbour whose direction from ci
            # is most *opposite* to the outward normal (i.e. most inward-facing),
            # which is the cell most directly "upstream" of this boundary edge.
            # This neighbour's WSE is used for gradient extrapolation at runtime.
            #
            # outward unit vector (equirectangular, sufficient for small distances):
            cos_lat  = cosd(cell_lat)
            out_dx   = (mid_lon - cell_lon) * cos_lat
            out_dy   = mid_lat - cell_lat
            out_len  = sqrt(out_dx^2 + out_dy^2)
            out_dx  /= max(out_len, 1e-12)
            out_dy  /= max(out_len, 1e-12)

            best_nb_idx = 0
            best_nb_L   = L       # fallback: use ghost L as interior distance
            best_dot    = -Inf    # most negative dot = most opposite to outward = most inward
            for nb_id in nbr_ids
                j = get(id_idx, nb_id, 0)
                j == 0 && continue
                nb  = cells[j]
                ndx = (nb.center_lon - cell_lon) * cos_lat
                ndy = nb.center_lat - cell_lat
                nlen = sqrt(ndx^2 + ndy^2)
                nlen < 1e-12 && continue
                dot = (ndx * out_dx + ndy * out_dy) / nlen   # cosine of angle with outward
                if dot > best_dot
                    best_dot    = dot
                    best_nb_idx = j
                    best_nb_L   = A5Grid._haversine_m(cell_lon, cell_lat,
                                                       nb.center_lon, nb.center_lat)
                end
            end
            # If best_dot < 0, the best neighbour is on the inward side — good.
            # If all neighbours are on the outward side (rare corner cell), fall back to 0.
            best_dot > 0.5 && (best_nb_idx = 0)   # no clearly inward neighbour

            push!(ghost_edges,   GhostEdge(i, width, L, sill, 0.0, use_sgs,
                                            best_nb_idx, best_nb_L))
            push!(ghost_cell_bc, default_bc)
        end
    end

    n_ghost = length(ghost_edges)
    @info "Ghost edges: $n_ghost created for $n_boundary boundary cells " *
          "(default BC: $default_bc)"

    return boundary_mask, ghost_edges, ghost_cell_bc
end

# ---------------------------------------------------------------------------
# GeoJSON boundary file parser
# ---------------------------------------------------------------------------

"""
    load_bc_file(path, cells, boundary_mask, ghost_edges, ghost_cell_bc,
                  default_bc) → Vector{BCType}

Parse a GeoJSON boundary-condition file and return an updated ghost_cell_bc
vector.  Each GeoJSON feature with a "bc_type" property overrides the BC type
for boundary cells whose centres fall within 1.5× the cell diameter of the
feature geometry (LineString or Polygon perimeter).

Valid `bc_type` property values: "Closed", "ZeroGradient", "Critical".
"FixedWSE" and "FixedQ" are parsed but log a not-yet-implemented warning.
Cells not matched by any feature retain `default_bc`.
"""
function load_bc_file(path          :: String,
                       cells         :: Vector{A5Grid.A5Cell},
                       boundary_mask :: BitVector,
                       ghost_edges   :: AbstractVector,
                       ghost_cell_bc :: AbstractVector,
                       default_bc    :: BCType)::Vector{BCType}
    isfile(path) || error("BC GeoJSON file not found: $path")

    updated_bc = copy(ghost_cell_bc)

    # Parse GeoJSON manually (no additional deps required; we only need
    # feature geometries and properties)
    raw = JSON3.read(read(path, String))
    features = if haskey(raw, :features)
        raw.features
    elseif haskey(raw, :type) && raw[:type] == "Feature"
        [raw]
    else
        error("BC file '$path' is not a GeoJSON Feature or FeatureCollection.")
    end

    for feat in features
        props = get(feat, :properties, nothing)
        props === nothing && continue
        bc_str = get(props, :bc_type, nothing)
        bc_str === nothing && (bc_str = get(props, :BCType, nothing))
        bc_str === nothing && continue

        # Warn on unimplemented types but continue (don't crash)
        bc_type = try
            BCType_from_string(string(bc_str))
        catch e
            @warn "BC file: unrecognised bc_type '$bc_str'. Skipping feature."
            continue
        end
        if bc_type in (FixedWSE, FixedQ)
            label = get(props, :label, "unnamed")
            @warn "BC type '$bc_type' (feature: '$label') is not yet implemented. " *
                  "Skipping. It will be available in a future session."
            continue
        end

        label = string(get(props, :label, "segment"))
        geom  = get(feat, :geometry, nothing)
        geom  === nothing && continue

        # Extract coordinate ring(s) for distance testing
        coord_rings = _extract_geojson_coords(geom)
        isempty(coord_rings) && continue

        # Estimate a per-cell tolerance: 1.5 × cell diameter
        # Cell diameter ≈ 2 × sqrt(area/π)
        n_matched = 0
        for (ge_idx, ge) in enumerate(ghost_edges)
            ci     = ge.cell_index
            c      = cells[ci]
            diam   = 2.0 * sqrt(A5Grid._polygon_area_m2(c.boundary) / π)
            tol    = 1.5 * diam

            # Distance from cell centre to nearest point on any coord ring
            min_d = Inf
            for ring in coord_rings
                for k in 1:length(ring)-1
                    lon1, lat1 = ring[k][1], ring[k][2]
                    lon2, lat2 = ring[k+1][1], ring[k+1][2]
                    d = _point_to_segment_dist(c.center_lon, c.center_lat,
                                               lon1, lat1, lon2, lat2)
                    d < min_d && (min_d = d)
                end
                # Also check last-to-first closing segment for polygons
                if length(ring) >= 2
                    lon1, lat1 = ring[end][1], ring[end][2]
                    lon2, lat2 = ring[1][1], ring[1][2]
                    d = _point_to_segment_dist(c.center_lon, c.center_lat,
                                               lon1, lat1, lon2, lat2)
                    d < min_d && (min_d = d)
                end
            end

            if min_d <= tol
                updated_bc[ge_idx] = bc_type
                n_matched += 1
            end
        end
        @info "BC file: feature '$label' ($bc_type) matched $n_matched ghost edges"
    end

    return updated_bc
end

# Helper: extract coordinate rings from a GeoJSON geometry object
function _extract_geojson_coords(geom)::Vector{Vector{Vector{Float64}}}
    gtype = string(get(geom, :type, ""))
    coords = get(geom, :coordinates, nothing)
    coords === nothing && return Vector{Vector{Vector{Float64}}}()

    if gtype == "LineString"
        ring = [[Float64(pt[1]), Float64(pt[2])] for pt in coords]
        return [ring]
    elseif gtype == "Polygon"
        # Outer ring only (index 1); inner rings (holes) ignored
        ring = [[Float64(pt[1]), Float64(pt[2])] for pt in coords[1]]
        return [ring]
    elseif gtype == "MultiLineString"
        return [[[Float64(pt[1]), Float64(pt[2])] for pt in seg] for seg in coords]
    elseif gtype == "MultiPolygon"
        return [[[Float64(pt[1]), Float64(pt[2])] for pt in poly[1]] for poly in coords]
    else
        @warn "BC file: unsupported geometry type '$gtype'. Only LineString and Polygon are supported."
        return Vector{Vector{Vector{Float64}}}()
    end
end

# Helper: approximate haversine distance from point (plon, plat) to segment
# (lon1,lat1)→(lon2,lat2), in metres.
function _point_to_segment_dist(plon :: Float64, plat :: Float64,
                                  lon1 :: Float64, lat1 :: Float64,
                                  lon2 :: Float64, lat2 :: Float64)::Float64
    # Project onto equirectangular plane centred at segment midpoint
    mlon = 0.5*(lon1+lon2)
    mlat = 0.5*(lat1+lat2)
    cos_lat = cosd(mlat)

    dx(a, b) = (a - b) * cos_lat * 111_320.0
    dy(a, b) = (a - b) * 111_320.0

    ax, ay = dx(lon1, mlon), dy(lat1, mlat)
    bx, by = dx(lon2, mlon), dy(lat2, mlat)
    px, py = dx(plon, mlon), dy(plat, mlat)

    abx, aby = bx - ax, by - ay
    len2 = abx^2 + aby^2
    if len2 < 1e-6
        return sqrt((px-ax)^2 + (py-ay)^2)
    end
    t = clamp(((px-ax)*abx + (py-ay)*aby) / len2, 0.0, 1.0)
    cx, cy = ax + t*abx, ay + t*aby
    return sqrt((px-cx)^2 + (py-cy)^2)
end

# ---------------------------------------------------------------------------
# BCIEntry — parsed record from a LISFLOOD-FP .bci file
# ---------------------------------------------------------------------------
# Defined here (rather than timeseries_io.jl) so that boundary_conditions.jl
# can reference it in apply_bci_free_entries! without a forward-include
# dependency.  timeseries_io.jl uses BCIEntry but does not define it.

"""
    BCIEntry

Parsed record from a LISFLOOD-FP .bci file.

Fields match the .bci column specification:
  boundary_type : one of 'N','E','S','W','P','F'
  x1, y1        : start of segment (lon/easting) or point location (for P)
  x2, y2        : end of segment — equal to x1,y1 for P entries
  bc_code       : "QVAR", "QFIX", "FREE", "HFIX", etc.
  bc_value      : series name (QVAR), numeric rate (QFIX/HFIX), or "" (FREE)
"""
struct BCIEntry
    boundary_type :: Char      # N, E, S, W, P, F
    x1            :: Float64   # start lon / easting
    y1            :: Float64   # start lat / northing
    x2            :: Float64   # end lon (= x1 for P entries)
    y2            :: Float64   # end lat (= y1 for P entries)
    bc_code       :: String    # QVAR, QFIX, FREE, HFIX, HVAR, ...
    bc_value      :: String    # series name, numeric string, or ""
end

# ---------------------------------------------------------------------------
# Apply .bci FREE entries to ghost_cell_bc
# ---------------------------------------------------------------------------

"""
    apply_bci_free_entries!(ghost_edges, ghost_cell_bc, cells, bci_entries)

For each FREE entry in a parsed .bci file, match to boundary ghost edges and
set their BC type to ZeroGradient.  Called after _build_ghost_edges when a
--inflow-bci file is loaded.

N/E/S/W entries produce a helpful unsupported warning rather than an error.
F (internal free boundary) entries are logged as unsupported.
"""
function apply_bci_free_entries!(ghost_edges   :: AbstractVector,
                                  ghost_cell_bc :: AbstractVector,
                                  cells         :: Vector{A5Grid.A5Cell},
                                  bci_entries   :: Vector{BCIEntry})
    warned_cardinal = false
    warned_F        = false

    for entry in bci_entries
        if entry.boundary_type in ('N', 'E', 'S', 'W')
            if !warned_cardinal
                @warn "BCI file contains cardinal-direction boundary entries " *
                      "(N/E/S/W). These are not supported in FloodA5 because the " *
                      "domain can be any polygon shape — cardinal directions have no " *
                      "natural mapping. To control boundary conditions by location, " *
                      "use a GeoJSON boundary file with --bc-file. " *
                      "See the FloodA5 documentation for the GeoJSON format."
                warned_cardinal = true
            end
            continue
        end

        if entry.boundary_type == 'F'
            if !warned_F
                @warn "BCI file contains 'F' (internal free boundary) entries. " *
                      "This boundary type is not yet supported in FloodA5 " *
                      "(no direct equivalent to the LISFLOOD-FP SGC internal boundary). " *
                      "These entries will be ignored. See PROJECT_STATE.md for " *
                      "planned reconsideration."
                warned_F = true
            end
            continue
        end

        # P entries: only FREE bc_code applies here; QVAR/QFIX handled elsewhere
        entry.bc_code != "FREE" && continue

        # Find the nearest boundary ghost edge to this point
        best_ge_idx = 0
        best_d      = Inf
        for (ge_idx, ge) in enumerate(ghost_edges)
            ci = ge.cell_index
            d  = A5Grid._haversine_m(entry.x1, entry.y1,
                                      cells[ci].center_lon, cells[ci].center_lat)
            d < best_d && (best_d = d; best_ge_idx = ge_idx)
        end

        best_ge_idx == 0 && continue
        best_d > 5000.0 && @warn "BCI FREE entry at ($(entry.x1), $(entry.y1)) " *
                                  "nearest boundary cell is $(round(best_d/1000, digits=1)) km away."

        ghost_cell_bc[best_ge_idx] = ZeroGradient
    end
end
