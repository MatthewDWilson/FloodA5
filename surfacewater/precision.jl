# surfacewater/precision.jl
# -------------------------
# 0.6.0 (2026-11): working precision of the device solvers (`--precision`).
#
# Float64 is the default and the reference. `--precision float32` (device
# backends: gpu, ka-cpu; standard flow and SWE) runs the fluxes and the
# momentum state in Float32, following the plan of record §7.3:
#
#   - cell volume, and the outflow / clamp totals, stay in Float64: every
#     face's exchange is added to one cell and subtracted from the other as
#     the same number, summed in Float64, so mass is conserved as in a
#     Float64 run;
#   - no water level z + h is formed in absolute terms and then differenced
#     (at z = 500 m a Float32 level resolves 3·10⁻⁵ m). Levels are taken
#     relative to a local datum — the face's sill (standard flow) or cell i's
#     bed (SWE) — from offsets that are differenced in Float64 once, when the
#     arrays are built: `RelLevels`. `AbsLevels` is the original arithmetic
#     (bit-identical to earlier versions) and is what Float64 runs use.
#     `RelLevels` in Float64 is available to tests: it separates the change
#     of formulation (round-off) from the change of precision.
#
# The kernel bodies are written once. `@fp T expr` rewrites the Float64
# literals and the named constants below in `expr` to `convert(T, ·)`, which
# is the identity for T = Float64 (the code generated is unchanged) and keeps
# a Float32 kernel free of Float64 arithmetic. A kernel takes T from one of
# its arrays.

"Water levels as z + h (the original arithmetic; Float64 runs)."
struct AbsLevels end
"Water levels relative to a local datum, from precomputed bed offsets (Float32-safe)."
struct RelLevels end

# the kernels' argument tuple may carry `lev`; the CPU step's tuples do not
@inline _lev(a::NamedTuple{names}) where {names} = Base.sym_in(:lev, names) ? a.lev : AbsLevels()

# Float64 constants that `@fp` converts along with the literals
const _FP_NAMES = (:_G, :HFLOW_THRESHOLD, :FROUDE_LIMIT, :SWE_HDRY, :SWE_HVEL, :MOM_DAMP_HW,
                   :ETA_FILTER_HW, :NaN, :Inf)

_fp_rewrite(T, x) = x
_fp_rewrite(T, x::Float64) = :(convert($T, $x))
_fp_rewrite(T, x::Symbol) = x in _FP_NAMES ? :(convert($T, $x)) : x
function _fp_rewrite(T, ex::Expr)
    ex.head in (:quote, :inert) && return ex
    if ex.head === :kw || ex.head === :(=)              # keep the name / target, rewrite the value
        return Expr(ex.head, ex.args[1] isa Expr ? _fp_rewrite(T, ex.args[1]) : ex.args[1],
                    (_fp_rewrite(T, a) for a in ex.args[2:end])...)
    elseif ex.head === :. && length(ex.args) == 2 && ex.args[2] isa QuoteNode
        return Expr(:., _fp_rewrite(T, ex.args[1]), ex.args[2])   # a.NaN is a field, not the constant
    end
    return Expr(ex.head, (_fp_rewrite(T, a) for a in ex.args)...)
end

"""
    @fp T expr

`expr` with every Float64 literal, and every constant in `_FP_NAMES`,
converted to the float type `T` (see the file header).
"""
macro fp(T, ex)
    return esc(_fp_rewrite(T, ex))
end

# "the ghost level is above the interior level" (_bates_ghost_flux: no inflow
# through an open boundary). Float64: above by 1e-10 m, as always. Float32
# cannot resolve that margin; there it is a few units of round-off of the
# level (which is relative to the sill), so that an extrapolated level equal
# to the interior one up to round-off still counts as equal.
@inline _ghost_above(wse_ghost::Float64, wse_ci::Float64) = wse_ghost >= wse_ci + 1e-10
@inline _ghost_above(wse_ghost::Float32, wse_ci::Float32) =
    wse_ghost >= wse_ci + 4 * eps(Float32) * max(abs(wse_ci), 1f0)
# Tikhonov factor of the η-filter fit (_eta_residual): only exact degeneracy
@inline _fit_reg(::Type{Float64}) = 1e-10
@inline _fit_reg(::Type{Float32}) = 1f-5

# Device arrays that stay Float64 at any working precision (mass accounting)
const _KEEP_F64 = (:volume, :V, :V0, :clamp, :removed, :g_dv)

"""
    _to_precision(T, a) → NamedTuple

The argument tuple `a` with its Float64 arrays and scalars converted to `T`
(nested tuples too), except the arrays in `_KEEP_F64`. For `T = Float64`
this is `a` itself: the same arrays, not copies.
"""
_to_precision(::Type{Float64}, a::NamedTuple) = a
_to_precision(::Type{T}, a::NamedTuple) where {T} =
    NamedTuple{keys(a)}(map((k, v) -> k in _KEEP_F64 ? v : _to_precision(T, v), keys(a), values(a)))
_to_precision(::Type{T}, x::AbstractArray{Float64}) where {T} = convert(Array{T}, x)
_to_precision(::Type{T}, x::Float64) where {T} = convert(T, x)
_to_precision(::Type{T}, x) where {T} = x
