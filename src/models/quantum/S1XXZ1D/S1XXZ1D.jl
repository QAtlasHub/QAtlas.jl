# ─────────────────────────────────────────────────────────────────────────────
# Spin-1 XXZ chain (S1XXZ1D)
#
# Hamiltonian:
#   H = J Σᵢ [Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Δ Sᶻᵢ Sᶻᵢ₊₁]
#
# with S = 1, J > 0 antiferromagnetic, Δ ∈ ℝ the easy-axis anisotropy.
#
# At Δ = 1 this reduces to the spin-1 Heisenberg chain (Haldane phase,
# gap Δ_∞ ≈ 0.41048 J — White-Huse 1993 DMRG).  For Δ ≠ 1 the model
# realises a rich phase diagram:
#
#   Δ → +∞   : large-Δ (Néel) gapped phase, c1 = c2 = 1
#   Δ ∈ (0,1): XY1 phase, gapless TLL-like
#   Δ < 0    : ferromagnetic / non-trivial regimes
#
# Phase 1 (this PR) supports MassGap only at Δ = 1 by delegation to
# `S1Heisenberg1D`.  Δ ≠ 1 throws `DomainError` and is deferred to
# Phase 2 (DMRG / TLL).
#
# References:
#   F. D. M. Haldane, [Haldane1983](@cite).
#   S. R. White, D. A. Huse, [WhiteHuse1993](@cite).
#   H. J. Schulz, [Schulz1986](@cite).
#   Y.-C. Tzeng, H. Yang, M.-F. Yang, Phys. Rev. B 96, 064419 (2017).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S         (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    S1XXZ1D(; J::Real = 1.0, Δ::Real = 1.0) <: AbstractQAtlasModel

Spin-1 antiferromagnetic XXZ chain,

    H = J Σᵢ [Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Δ Sᶻᵢ Sᶻᵢ₊₁],   S = 1, J > 0.

The default `Δ = 1.0` places the model at the SU(2)-symmetric Haldane
point, where Phase-1 closed-form `MassGap` is available by delegation
to [`S1Heisenberg1D`](@ref).
"""
struct S1XXZ1D <: AbstractQAtlasModel
    J::Float64
    Δ::Float64
    function S1XXZ1D(J::Real, Δ::Real)
        J > 0 || throw(DomainError(J, "S1XXZ1D requires J > 0; got J = $J."))
        return new(Float64(J), Float64(Δ))
    end
end
S1XXZ1D(; J::Real=1.0, Δ::Real=1.0) = S1XXZ1D(J, Δ)

"""
    fetch(model::S1XXZ1D, ::MassGap, ::Infinite; J, Δ, kwargs...) -> Float64

Bulk mass gap of the spin-1 XXZ chain.

Phase 1: supported only at the Heisenberg point `Δ = 1`, where the
Haldane gap is delegated to `S1Heisenberg1D`,

    Δ_∞ ≈ 0.41048 J   (White-Huse 1993 DMRG).

For `Δ ≠ 1` the result is not a closed-form literature constant
(XY1 / large-Δ Néel phase diagram, Schulz 1986; Tzeng-Yang-Hsu 2017)
and a `DomainError` is raised — Phase 2 will plug in DMRG / TLL.
"""
function fetch(m::S1XXZ1D, ::MassGap, ::Infinite; J::Real=m.J, Δ::Real=m.Δ, kwargs...)
    J > 0 || throw(DomainError(J, "S1XXZ1D MassGap requires J > 0; got J = $J."))
    if !isone(Δ)
        throw(
            DomainError(
                Δ,
                "S1XXZ1D MassGap: closed-form Haldane gap supported only at Δ = 1 " *
                "(spin-1 Heisenberg, White-Huse 1993 DMRG). Δ ≠ 1 traverses XY1/large-Δ Néel " *
                "phase diagram (Schulz 1986; Tzeng-Yang-Hsu 2017) — deferred to Phase 2. " *
                "Got Δ = $Δ.",
            ),
        )
    end
    return QAtlas.fetch(S1Heisenberg1D(; J=J), MassGap(), Infinite())
end

"""
    fetch(model::S1XXZ1D, ::Energy{:per_site}, ::Infinite; J, Δ, kwargs...) -> Float64

Ground-state energy per site of the spin-1 XXZ chain.

Phase 1: supported only at the Heisenberg point `Δ = 1`, where the
result is delegated to `S1Heisenberg1D`,

    e₀ ≈ -1.40148403897 J   (White-Huse 1993 DMRG).

For `Δ ≠ 1` no closed-form literature constant is available (XY1 /
large-Δ Néel phase diagram, Schulz 1986; Tzeng-Yang-Hsu 2017); a
`DomainError` is raised — Phase 2 will plug in DMRG / TLL.
"""
function fetch(
    m::S1XXZ1D, ::Energy{:per_site}, ::Infinite; J::Real=m.J, Δ::Real=m.Δ, kwargs...
)
    J > 0 || throw(DomainError(J, "S1XXZ1D Energy{:per_site} requires J > 0; got J = $J."))
    if !isone(Δ)
        throw(
            DomainError(
                Δ,
                "S1XXZ1D Energy{:per_site}: closed-form energy density supported only at Δ = 1 " *
                "(spin-1 Heisenberg, White-Huse 1993 DMRG). Δ ≠ 1 traverses XY1/large-Δ Néel " *
                "phase diagram (Schulz 1986; Tzeng-Yang-Hsu 2017) — deferred to Phase 2. " *
                "Got Δ = $Δ.",
            ),
        )
    end
    return QAtlas.fetch(S1Heisenberg1D(; J=J), Energy(:per_site), Infinite())
end
