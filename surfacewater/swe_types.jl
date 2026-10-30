# surfacewater/swe_types.jl
# -------------------------
# R1 (2026-10): data layout of the full shallow-water reference solver
# (`--flow-model swe`, surfacewater/swe_reference.jl). Included before
# FlowState so that FlowState.swe can be concretely typed.
#
# Struct-of-arrays throughout, and every per-step pass is either a per-face
# kernel writing only its own face slots or a per-cell gather over the fixed
# cell → face index (StepWorkspace.cptr/cedge for interior faces, gptr/gface
# here for boundary faces). No atomics, no shared writes, no Dicts or
# Symbols in the hot loops, boundary types as Int8 codes: each pass maps onto
# one GPU kernel (KernelAbstractions) with the same body. GPU code itself
# belongs on its own branch.

const SWE_BC_WALL         = Int8(0)   # Closed: reflecting wall
const SWE_BC_TRANSMISSIVE = Int8(1)   # ZeroGradient (and FixedWSE/FixedQ until implemented)
const SWE_BC_CRITICAL     = Int8(2)   # Critical: outflow at least at critical speed

"""
    SWEWork

Geometry and work arrays of the R1 shallow-water solver, built once by
`_build_swe_work` in `initialise_flow_model` (`SWEFlow`). The conserved
state is not here: cell volume is `state.volume` (so sources and rainfall
act on it unchanged) and the cell discharge per unit width (hu, hv) is
`state.qvec_u` / `state.qvec_v` (m²/s, local east/north frame).
"""
mutable struct SWEWork
    order    :: Int                # 1 (first order, Euler) or 2 (MUSCL + SSP-RK2)
    visc_ce  :: Float64            # C_E of the eddy viscosity ν_t = C_E·h·u* (0 = off)
    open_inflow :: Bool            # ZeroGradient faces also admit inflow (benchmarks); default false
    # interior faces (EdgeList order): unit normal cell_i → cell_j, and the
    # face-midpoint offsets from the two cell centres (m, local frames)
    nx  :: Vector{Float64}
    ny  :: Vector{Float64}
    rxi :: Vector{Float64}
    ryi :: Vector{Float64}
    rxj :: Vector{Float64}
    ryj :: Vector{Float64}
    # boundary faces (ghost_edges order): cell, outward normal, length, bed,
    # midpoint offset, BC code; and the cell → boundary-face CSR
    gcell :: Vector{Int}
    gnx   :: Vector{Float64}
    gny   :: Vector{Float64}
    gw    :: Vector{Float64}
    gz    :: Vector{Float64}
    grx   :: Vector{Float64}
    gry   :: Vector{Float64}
    gbc   :: Vector{Int8}
    gptr  :: Vector{Int}
    gface :: Vector{Int}
    perim :: Vector{Float64}       # cell perimeter incl. boundary faces (m)
    # per-stage scratch: primitives and limited gradients (n)
    h   :: Vector{Float64}
    u   :: Vector{Float64}
    v   :: Vector{Float64}
    eta :: Vector{Float64}
    dex :: Vector{Float64}
    dey :: Vector{Float64}
    dux :: Vector{Float64}
    duy :: Vector{Float64}
    dvx :: Vector{Float64}
    dvy :: Vector{Float64}
    dhx :: Vector{Float64}         # depth gradient (order 2: face bed z_f = η_f − h_f)
    dhy :: Vector{Float64}
    # interior-face fluxes per unit length: mass, x/y momentum, the two
    # hydrostatic-reconstruction pressure corrections, viscous stress (ne)
    F1 :: Vector{Float64}
    F2 :: Vector{Float64}
    F3 :: Vector{Float64}
    PL :: Vector{Float64}
    PR :: Vector{Float64}
    T2 :: Vector{Float64}
    T3 :: Vector{Float64}
    # boundary-face fluxes out of the cell, per unit length (ng)
    G1 :: Vector{Float64}
    G2 :: Vector{Float64}
    G3 :: Vector{Float64}
    # RK stage storage and the mass removed by clamping negative volumes (n)
    V0  :: Vector{Float64}
    qx0 :: Vector{Float64}
    qy0 :: Vector{Float64}
    clamp_vol :: Vector{Float64}
end
