# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — Z-axis (longitudinal) thermodynamic-limit
# observables.
#
# All quantities here are in the thermodynamic limit (`Infinite()`) and so
# carry no `N` dependence.  The closed forms come from the Pfeuty (1970)
# free-fermion + Bogoliubov solution at T = 0:
#
#   m_z = (1 - (h/J)²)^{1/8}    (h < J, ordered)
#       = 0                     (h ≥ J, disordered)
#
#   ξ   = 1 / (2|h - J|)         (T = 0 correlation length, gapped phase)
#       = ∞                     (h = J, critical)
#
# `SusceptibilityZZ` at `Infinite` is a per-site uniform longitudinal
# susceptibility.  A closed-form has been worked out in the
# Calabrese–Mussardo / Sachdev tradition (form-factor expansion +
# integration over ω), but the result is intricate and prone to typos.
# QAtlas defines `χ_zz` at `Infinite` as the large-N OBC limit of the
# fluctuation-dissipation expression
#
#   χ_zz(β) = (β/N) · Var(M_z) = (β/N) · Σ_{ij} ⟨σᶻ_i σᶻ_j⟩_β  (since ⟨σᶻ⟩ = 0
#                                                                in OBC by Z₂)
#
# computed via the existing `_zz_uniform_susceptibility(N, …)` Pfaffian
# routine at a large N (see kwargs).  The Infinite() return value is just
# the OBC value at the chosen N — this is a deliberate practical compromise
# that lets `SusceptibilityZZ, Infinite()` cohabit with the other Z-axis
# Infinite quantities, while the user is told via the docstring how to
# tighten the bound by raising `N_proxy`.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Spontaneous magnetisation / MagnetizationZ at T = 0 (Pfeuty)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_pfeuty_mz(J, h) -> Float64

Pfeuty (1970) closed-form spontaneous longitudinal magnetisation per
site for the infinite TFIM at T = 0 (Pauli convention, eigenvalues
±1):

    m_z = (1 - (h/J)²)^{1/8}    (h < J, ordered phase)
        = 0                     (h ≥ J, disordered phase)

The exponent `β = 1/8` is the Onsager–Yang exponent of the 2D Ising
universality class.

Returns the *positive* branch of the spontaneously-broken Z₂ doublet
(the negative branch is `-m_z`).  `m_z = 0` exactly at the critical
point `h = J`.
"""
function _tfim_pfeuty_mz(J::Real, h::Real)
    if abs(h) >= abs(J)
        return 0.0
    else
        return (1 - (h / J)^2)^(1 / 8)
    end
end

"""
    fetch(model::TFIM, ::MagnetizationZ, ::Infinite; kwargs...) -> Float64

Spontaneous longitudinal magnetisation per site of the infinite TFIM
at T = 0, `m_z = (1 - (h/J)²)^{1/8}` for `h < J`, else 0
(Pfeuty 1970).  Returns the positive branch of the Z₂-broken doublet.

The result is the *T = 0* order parameter; the function does not take
a `beta` kwarg because `m_z(T > 0) = 0` for any finite chain in the
thermodynamic limit (Mermin-Wagner is irrelevant in 1D, but the broken
phase requires explicit symmetry breaking; m_z(T,h) ≠ 0 only at
T = 0).  Pass `beta = Inf` if you want to be explicit; it is ignored.
"""
function fetch(model::TFIM, ::MagnetizationZ, ::Infinite; kwargs...)
    return _tfim_pfeuty_mz(model.J, model.h)
end

"""
    fetch(model::TFIM, ::SpontaneousMagnetization, ::Infinite; kwargs...) -> Float64

Same value as `fetch(::TFIM, ::MagnetizationZ, ::Infinite)`, exposed
under the order-parameter name commonly used in the Pfeuty / 2D-Ising
universality literature.  The struct `SpontaneousMagnetization` is
shared with `IsingSquare` (defined in
`src/models/classical/IsingSquare/IsingSquare.jl`); this method adds
the TFIM-at-`Infinite` branch.
"""
function fetch(model::TFIM, ::SpontaneousMagnetization, ::Infinite; kwargs...)
    return _tfim_pfeuty_mz(model.J, model.h)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Correlation length at T = 0
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::TFIM, ::CorrelationLength, ::Infinite; kwargs...) -> Float64

T = 0 correlation length of the infinite TFIM in the **relativistic
continuum convention** (inverse mass gap),

    ξ = 1 / (2|h - J|)        (gapped phase)
    ξ = Inf                   (critical point h = J)

set by the lattice mass gap `Δ = 2|h - J|` via the universal IR
relation `ξ = 1/Δ` (with `v_F = 1` implicit; lattice units). Tracks
[`MassGap`](@ref) at `Infinite`.

# Convention note

Three legitimate conventions exist for the TFIM correlation length on
the lattice; QAtlas exposes the first by default for consistency with
`MassGap`:

| Convention                        | Formula                       | Origin                                 |
|-----------------------------------|-------------------------------|----------------------------------------|
| **Inverse mass gap** (this fetch) | `1 / (2|h - J|)`              | relativistic continuum / Sachdev IR    |
| Pfeuty 1970 longitudinal          | `1 / log(max(J,h) / min(J,h))` | lattice JW-fermion <σᶻσᶻ> decay exact  |
| Sachdev lattice relativistic      | `min(J,h) / |h - J|`           | E²(k) = Δ² + v²k² with v_F = 2·min     |

The three agree to leading order near criticality (|h - J| << max). For
exact lattice decay of the longitudinal correlator at any (J, h), use
the Pfeuty form externally.

In QAtlas convention ξ is dimensionless (in units of the lattice
spacing).
"""
function fetch(model::TFIM, ::CorrelationLength, ::Infinite; kwargs...)
    Δ = 2 * abs(model.h - model.J)
    return Δ ≤ 0 ? Inf : 1 / Δ
end

# ═══════════════════════════════════════════════════════════════════════════════
# Longitudinal susceptibility at Infinite — OBC large-N proxy
# ═══════════════════════════════════════════════════════════════════════════════

# NOTE: The (ω-omitted) static SusceptibilityZZ, Infinite was originally
# defined here.  It is now folded into the dynamic-aware router inside
# TFIM_infinite_dynamics.jl (same pattern as ZZStructureFactor) so the
# dynamic-ω branch can be added without method-overwrite conflicts.
# The static behaviour is preserved bit-for-bit: the router falls
# through to _zz_uniform_susceptibility(N_proxy, J, h, β) when
# ω === nothing.

# ═══════════════════════════════════════════════════════════════════════════════
# Infinite-N static σᶻσᶻ structure factor
# ═══════════════════════════════════════════════════════════════════════════════

# NOTE: The static (ω-omitted) `ZZStructureFactor, Infinite` was originally
# defined here.  In v0.18 it was unified with the dynamic variant inside
# `TFIM_infinite_dynamics.jl` to avoid a method-overwrite conflict — the
# new method is a router that branches on whether `ω` is supplied.  The
# static behaviour is preserved bit-for-bit (the router falls through to
# `_zz_static_structure_factor(N_proxy, ...)` when `ω === nothing`).
