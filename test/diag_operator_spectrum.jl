# test/diag_operator_spectrum.jl
# ------------------------------
# Damping diagnostics (D1, 2026-10): the spatial free-surface operator on its
# own. K maps a cell WSE field to cell volume tendencies per unit area,
# (Kη)_i = Σ_f w_f·dWSE_n,f(η)/(L_f·A_i), for
#   two-point   dWSE_n = cosθ·(η_i − η_j)              (symmetric in the area inner product)
#   diamond     the diamond face-normal gradient (falls back to two-point where invalid)
#   blend ω     (1−ω)·diamond + ω·two-point
#   Ks          the symmetric part of the diamond operator in the area inner
#               product, A⁻¹·(AK + (AK)ᵀ)/2 — a diagnostic only (§ recommendations 1)
#
# For each: number of complex eigenvalues, max |Im λ|·ℓ², max Re λ·ℓ², the
# frictionless forward–backward growth factor per step at Courant 0.7 and
# 0.35 (uniform depth, local CFL length ℓ), and the linear-field residual —
# how far K·η is from zero for a planar η = a·x + b·y on interior cells (0 for
# a consistent operator; the diamond operator is exact there).
#
#   FLOODA5_NO_VIS=1 julia --project=. -t 4 test/diag_operator_spectrum.jl [MESH.parquet]
# (default: the planar res-16 mesh test/planar_embankment/planar_mesh16_std.parquet)

ENV["FLOODMODEL_INCLUDE_ONLY"] = "1"
include(joinpath(@__DIR__, "..", "FloodModel.jl"))
ENV["FLOODMODEL_INCLUDE_ONLY"] = ""

using LinearAlgebra, Logging, Printf, Statistics

p = isempty(ARGS) ? joinpath(@__DIR__, "planar_embankment", "planar_mesh16_std.parquet") : ARGS[1]
st = with_logger(SimpleLogger(stderr, Logging.Error)) do
    initialise_flow_model(load_mesh_geoparquet(p), StandardFlow(); cfl_mode = :local,
                          gradient_correction = true, face_flux_method = :diamond)
end
E = st.edges; n = length(st.cell_ids); A = st.cell_area; T = st.diamond_table
ℓmin = minimum(st.cell_length)

function Kmat(ω)                    # ω = 1: two-point; ω = 0: diamond
    K = zeros(n, n)
    for c in 1:n
        x = zeros(n); x[c] = 1.0
        for e in 1:E.n_edges
            i, j = E.cell_i[e], E.cell_j[e]
            tp = E.cos_theta[e] * (x[i] - x[j])
            d = ω == 1.0 ? NaN : diamond_dWSE_n(T, e, x[i], x[j], x, E.L[e])
            d = isnan(d) ? tp : (1 - ω) * d + ω * tp
            f = E.width[e] * d / E.L[e]; isfinite(f) || continue
            K[i, c] += f / A[i]; K[j, c] -= f / A[j]
        end
    end
    return K
end

# growth per step of the frictionless forward–backward update, worst mode
function growth(λ, C)
    maximum(λ) do l
        κ = C^2 * ℓmin^2 * l
        maximum(abs, ((2 - κ) .+ [1, -1] .* sqrt(complex((2 - κ)^2 - 4))) ./ 2)
    end
end

# planar field residual on interior cells
lon = st.cell_lons; lat = st.cell_lats
η_lin = 1e-3 .* 111_000 .* (lon .- mean(lon)) .+ 0.7e-3 .* 111_000 .* (lat .- mean(lat))
interior = [!st.boundary_mask[i] && all(j -> j == 0 || !st.boundary_mask[j], view(st.adj_matrix, :, i)) for i in 1:n]
resid(K) = maximum(abs, (K * η_lin)[interior]) * ℓmin^2

Kd = Kmat(0.0); Kt = Kmat(1.0)
Da = Diagonal(A); Ks = Diagonal(1 ./ A) * ((Da * Kd + (Da * Kd)') / 2)
ops = [("two-point", Kt), ("diamond", Kd), ("blend ω=0.25", 0.75 .* Kd .+ 0.25 .* Kt),
       ("blend ω=0.5", 0.5 .* Kd .+ 0.5 .* Kt), ("Ks (sym. part of diamond)", Ks)]
println("Operator spectrum: ", basename(p), " ($n cells), min ℓ = ", round(ℓmin, digits = 2), " m")
@printf("%-26s %9s %13s %11s %11s %11s %13s\n", "operator", "#complex", "max|Im|·ℓ²", "max Re·ℓ²",
        "growth 0.7", "growth 0.35", "lin. residual")
for (lab, K) in ops
    λ = eigvals(K)
    nc = count(z -> abs(imag(z)) > 1e-12 * maximum(abs, λ), λ)
    @printf("%-26s %9d %13.3e %11.3f %11.5f %11.5f %13.3e\n", lab, nc,
            maximum(z -> abs(imag(z)), λ) * ℓmin^2, maximum(real, λ) * ℓmin^2,
            growth(λ, 0.7), growth(λ, 0.35), resid(K))
end
