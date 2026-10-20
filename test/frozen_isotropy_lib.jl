# test/frozen_isotropy_lib.jl
# ---------------------------
# Frozen-coefficient friction-limit isotropy test (WP1, 2026-10): the Julia
# port of projectinstructions/review_scripts/friction_isotropy.py, run through
# the real solver code path rather than a re-implementation.
#
# Set a uniform depth h on a planar bed of slope S at bearing φ (so the WSE is
# planar too), then iterate step_standard! with the volumes reset after every
# step. The face fluxes relax to the friction-dominated steady state of the
# scheme for that fixed WSE field. The cell discharge vectors are then
# reconstructed (Phase F, WLSQ) and compared with the downslope direction −∇η.
#
#   deflection   angle(q_i) − angle(−∇η), per interior cell, degrees
#   conductance  |q_i| / √S  (equal for every bearing on an isotropic scheme)
#
# Theory (review §3.1): component friction gives q_n ∝ √|g_n|, a deflection of
# up to 12.8° on a Cartesian grid and a diagonal/cardinal conductance ratio of
# 2^¼; vector friction with an exact face-normal gradient gives none.
#
# Requires FloodModel.jl to be included already (FLOODMODEL_INCLUDE_ONLY).

using Statistics

"""
    frozen_isotropy(mesh; angles_deg, S, h, dt, iters, kwargs...) → NamedTuple

Run the frozen-coefficient test on `mesh` for each bearing in `angles_deg`.
`kwargs` are passed to `initialise_flow_model` (friction_model,
momentum_model, gradient_correction, face_flux_method, …). Interior cells are
those that are not on the domain edge and have no domain-edge neighbour.

Returns `(angles, mean_defl, rms_defl, cond, max_abs_mean_defl, cond_ratio,
n_interior)`; `mean_defl`/`rms_defl` in degrees, one entry per bearing.
"""
function frozen_isotropy(mesh::A5Mesh;
                         angles_deg = 0:5:85,
                         S::Float64  = 1e-4,
                         h::Float64  = 0.5,
                         dt::Float64 = 20.0,
                         iters::Int  = 300,
                         kwargs...)
    R    = A5Grid._EARTH_R
    lons = [c.center_lon for c in mesh.cells]
    lats = [c.center_lat for c in mesh.cells]
    lon0 = mean(lons); lat0 = mean(lats)
    x = deg2rad.(lons .- lon0) .* R .* cosd(lat0)
    y = deg2rad.(lats .- lat0) .* R

    mean_defl = Float64[]; rms_defl = Float64[]; cond = Float64[]
    n_interior = 0
    for φ in angles_deg
        gx, gy = S * cosd(φ), S * sind(φ)
        mesh.static_vars["elevation"] = gx .* x .+ gy .* y
        state = initialise_flow_model(mesh, StandardFlow(); kwargs...)
        for k in eachindex(state.ghost_cell_bc); state.ghost_cell_bc[k] = Closed; end
        n = length(state.cell_ids)
        vol0 = [h * state.cell_area[i] for i in 1:n]
        for _ in 1:iters
            copyto!(state.volume, vol0); fill!(state.water_depth, h)
            step_standard!(state, dt)
        end
        _reconstruct_qvec!(state, false)

        interior = [!state.boundary_mask[i] &&
                    all(j -> j == 0 || !state.boundary_mask[j], view(state.adj_matrix, :, i))
                    for i in 1:n]
        idx = findall(interior)
        n_interior = length(idx)
        ψ0 = atand(-gy, -gx)
        d  = [mod(atand(state.qvec_v[i], state.qvec_u[i]) - ψ0 + 180.0, 360.0) - 180.0 for i in idx]
        push!(mean_defl, mean(d))
        push!(rms_defl,  sqrt(mean(d .^ 2)))
        push!(cond, mean(hypot(state.qvec_u[i], state.qvec_v[i]) for i in idx) / sqrt(S))
    end
    return (angles = collect(angles_deg), mean_defl = mean_defl, rms_defl = rms_defl,
            cond = cond, max_abs_mean_defl = maximum(abs, mean_defl),
            cond_ratio = maximum(cond) / minimum(cond), n_interior = n_interior)
end
