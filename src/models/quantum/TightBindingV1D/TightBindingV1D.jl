# ─────────────────────────────────────────────────────────────────────────────
# TightBindingV1D — 1D spinless fermion chain with nearest-neighbour interaction
# (the "t-V model"; Jordan-Wigner-equivalent to the spin-1/2 XXZ chain).
#
# Hamiltonian:
#
#     H = -t  Σ_i (c†_i c_{i+1} + h.c.)
#         + V  Σ_i n_i n_{i+1}
#         - μ  Σ_i n_i ,
#
#     t > 0,   V ∈ ℝ,   μ ∈ ℝ.
#
# At V = 0 the chain reduces to the free-fermion tight-binding model with
# closed-form spectrum
#
#     ε(k) = -2t cos(k) - μ ,
#
# from which
#
#     MassGap         = max(0, |μ| - 2t)              (insulating for |μ| ≥ 2t)
#     FermiVelocity   = 2t · sin(arccos(-μ/(2t)))     (for |μ| < 2t)
#
# follow immediately.  At V ≠ 0 the model is, via the Jordan-Wigner
# transformation, equivalent to the XXZ chain with V/t setting the
# anisotropy Δ_XXZ; the gap and Fermi velocity then require the
# Yang-Yang (1966) Bethe-ansatz result.  Phase 1 of this file exposes
# only the V = 0 free-fermion point; V ≠ 0 raises `DomainError` and
# will be addressed in Phase 2 via the JW ↔ XXZ1D mapping.
#
# NOTE: We deliberately do NOT delegate to `TightBinding1D` (PR #316,
# not yet merged on `main`); the free-fermion closed forms are
# implemented inline here.
#
# References:
#   - C. N. Yang, C. P. Yang, "One-dimensional chain of anisotropic
#     spin-spin interactions", Phys. Rev. 150, 321 (1966).
#   - N. W. Ashcroft, N. D. Mermin, *Solid State Physics* (1976),
#     Chapter 9 (tight-binding bands).
# ─────────────────────────────────────────────────────────────────────────────

"""
    TightBindingV1D(t::Real, V::Real, μ::Real)
    TightBindingV1D(; t=1.0, V=0.0, μ=0.0)

1D spinless-fermion t-V chain (NN hopping `t > 0`, NN density-density
interaction `V`, chemical potential `μ`).  At `V = 0` the model is the
free-fermion tight-binding chain; at `V ≠ 0` it is Jordan-Wigner
equivalent to the spin-1/2 XXZ chain (Phase 2).

The default keyword constructor lands at the closed-form free-fermion
point `(t, V, μ) = (1, 0, 0)`.
"""
struct TightBindingV1D <: AbstractQAtlasModel
    t::Float64
    V::Float64
    μ::Float64
    function TightBindingV1D(t::Real, V::Real, μ::Real)
        t > 0 || throw(DomainError(t, "TightBindingV1D requires t > 0; got t = $t."))
        return new(Float64(t), Float64(V), Float64(μ))
    end
end
TightBindingV1D(; t::Real=1.0, V::Real=0.0, μ::Real=0.0) = TightBindingV1D(t, V, μ)

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — V = 0 free-fermion closed form
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::MassGap, ::Infinite;
          t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Single-particle gap of the tight-binding t-V chain at `V = 0`:

    Δ = max(0, |μ| - 2t) .

Gapless inside the band (|μ| < 2t), zero at the band edges (|μ| = 2t),
linear-in-|μ| insulator gap outside (|μ| > 2t).

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::MassGap,
    ::Infinite;
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(DomainError(t, "TightBindingV1D MassGap requires t > 0; got t = $t."))
    if !isapprox(V, 0.0; atol=1e-12)
        throw(
            DomainError(
                V,
                "TightBindingV1D MassGap: V ≠ 0 (Jordan-Wigner-equivalent to interacting XXZ, " *
                "Yang-Yang 1966) deferred to Phase 2. Got V = $V.",
            ),
        )
    end
    return max(0.0, abs(μ) - 2t)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Fermi velocity — V = 0 free-fermion closed form
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::TightBindingV1D, ::FermiVelocity, ::Infinite;
          t=m.t, V=m.V, μ=m.μ, kwargs...) -> Float64

Fermi velocity of the V = 0 tight-binding t-V chain,

    v_F = 2t · sin(k_F),    k_F = arccos(-μ/(2t)) ,

valid in the metallic regime `|μ| < 2t`.  For `|μ| ≥ 2t` there is no
Fermi surface (band insulator); a `DomainError` is raised.

`V ≠ 0` raises `DomainError` (Phase 2 via JW ↔ XXZ1D).
"""
function fetch(
    m::TightBindingV1D,
    ::FermiVelocity,
    ::Infinite;
    t::Real=m.t,
    V::Real=m.V,
    μ::Real=m.μ,
    kwargs...,
)
    t > 0 || throw(
        DomainError(t, "TightBindingV1D FermiVelocity requires t > 0; got t = $t."),
    )
    if !isapprox(V, 0.0; atol=1e-12)
        throw(
            DomainError(
                V,
                "TightBindingV1D FermiVelocity: V ≠ 0 deferred to Phase 2. Got V = $V.",
            ),
        )
    end
    abs(μ) < 2t || throw(
        DomainError(
            μ,
            "TightBindingV1D FermiVelocity: no Fermi surface for |μ| ≥ 2t (insulating regime). " *
            "Got μ = $μ, t = $t.",
        ),
    )
    return 2t * sin(acos(-μ / (2t)))
end
