# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: XX (Δ = 0) free-fermion quench observables for
# `XXZ1D` at `Infinite()`.
#
# Validates `fetch(::XXZ1D, ::LoschmidtRateFunction, ::Infinite;
#                  initial::XXZ1D, t::Real)` against the closed-form
# analytical result derived in `src/models/quantum/XXZ/XXZ_xx_quench.jl`:
#
#   λ(t) ≡ 0   for every (J_initial, J_final, t) at Δ = 0.
#
# Reason: |ψ₀⟩ = |GS(J_initial)⟩ is a number eigenstate of H_XX(J_final)
# in the shared plane-wave basis (no Bogoliubov pairing rotation
# between the two XX Hamiltonians), so e^{-iH_f t}|ψ₀⟩ is a pure phase
# and |L(t)| = 1.  This holds at every sgn-J combination including the
# sign-flip case (the static overlap |⟨GS(J₀)|GS(J_f)⟩| does vanish at
# sign-flip — Anderson orthogonality — but that is a different quantity
# from the Loschmidt autocorrelation).
#
# A `Δ ≠ 0` final or initial model raises `DomainError` — this is the
# phase-1 contract; the general-Δ Loschmidt route lives behind issue
# #108 / #143.
#
# Conventions match PR #157 (`feat/issue-95-xx-finite-t`):
#   * spin convention `Sᵅ = σᵅ/2`
#   * single-particle dispersion `ε_J(k) = J cos(k)`
#   * GS energy density `-J/π` (XX point)
#
# Test running:
#
#   julia --project=test test/standalone/test_xxz_xx_quench.jl
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

const _XX = XXZ1D(; J=1.0, Δ=0.0)
const _LE_RATE = LoschmidtRateFunction()

# ─────────────────────────────────────────────────────────────────────────────

@testset "XXZ1D Δ=0 LoschmidtEcho — smoke / public API" begin
    m_f = XXZ1D(; J=1.0, Δ=0.0)
    m_0 = XXZ1D(; J=1.0, Δ=0.0)
    λ = QAtlas.fetch(m_f, _LE_RATE, Infinite(); initial=m_0, t=1.0)
    @test λ isa Float64
end

@testset "XXZ1D Δ=0 LoschmidtEcho — no-quench gives λ = 0" begin
    # H_initial = H_final exactly — λ(t) = 0 for every t.
    for t in (0.0, 0.5, 1.0, 5.0, 17.3)
        λ = QAtlas.fetch(_XX, _LE_RATE, Infinite(); initial=_XX, t=t)
        @test λ == 0.0
    end
end

@testset "XXZ1D Δ=0 LoschmidtEcho — t = 0 always gives λ = 0" begin
    # Initial and final differ but t = 0: |⟨ψ₀|ψ₀⟩|² = 1 trivially.
    for (J0, Jf) in ((1.0, 0.5), (0.5, 1.0), (2.0, 0.3))
        m_0 = XXZ1D(; J=J0, Δ=0.0)
        m_f = XXZ1D(; J=Jf, Δ=0.0)
        λ = QAtlas.fetch(m_f, _LE_RATE, Infinite(); initial=m_0, t=0.0)
        @test λ == 0.0
    end
end

@testset "XXZ1D Δ=0 LoschmidtEcho — same-sign J quench: λ = 0" begin
    # Same-sign Fermi sea ⇒ amplitude is a pure phase ⇒ λ ≡ 0.
    for (J0, Jf) in ((0.5, 1.0), (1.0, 0.5), (1.0, 2.0), (2.5, 0.7))
        m_0 = XXZ1D(; J=J0, Δ=0.0)
        m_f = XXZ1D(; J=Jf, Δ=0.0)
        for t in (0.5, 1.0, 3.7)
            λ = QAtlas.fetch(m_f, _LE_RATE, Infinite(); initial=m_0, t=t)
            @test λ == 0.0
        end
    end
end

@testset "XXZ1D Δ=0 LoschmidtEcho — pinned value (J0=1, Jf=1, t=1.0)" begin
    # Brief asks for atol=1e-8 at (J=1, h_0=0.5, h_f=0.0, t=1.0).  The
    # XXZ1D struct has no `h` field in this release; we pin the
    # closest meaningful no-quench evaluation instead — λ is exactly
    # 0 by the (♣)/(♠) derivation.  This freezes the contract that
    # any future generalisation (magnetic-field XX or full Bogoliubov
    # XY) must reproduce λ = 0 at the no-quench fixed point.
    λ = QAtlas.fetch(
        XXZ1D(; J=1.0, Δ=0.0), _LE_RATE, Infinite(); initial=XXZ1D(; J=1.0, Δ=0.0), t=1.0
    )
    @test isapprox(λ, 0.0; atol=1e-8)
end

@testset "XXZ1D Δ=0 LoschmidtEcho — λ(t) ≥ 0 and finite (same-sign)" begin
    # Same-sign quench yields finite (= 0) λ.
    for (J0, Jf) in ((1.0, 0.5), (1.0, 2.0), (3.0, 0.1))
        m_0 = XXZ1D(; J=J0, Δ=0.0)
        m_f = XXZ1D(; J=Jf, Δ=0.0)
        for t in (0.0, 0.1, 1.0, 10.0)
            λ = QAtlas.fetch(m_f, _LE_RATE, Infinite(); initial=m_0, t=t)
            @test λ ≥ 0.0
            @test isfinite(λ)
        end
    end
end

@testset "XXZ1D Δ=0 LoschmidtEcho — sign-flip quench: λ ≡ 0" begin
    # sgn J flips → the two Fermi seas are complementary, but |ψ₀⟩ is
    # still a number eigenstate of H_f in the shared plane-wave basis,
    # so |L(t)| = 1 and λ ≡ 0.  Anderson orthogonality of the static
    # |⟨GS(J₀)|GS(J_f)⟩| vanishing is a different quantity from the
    # Loschmidt autocorrelation and does not enter λ(t).  Verified
    # numerically by Slater-determinant evolution det(D₀† U(t) D₀) at
    # finite N up to N = 20: |L|² = 1.0 to machine precision for all
    # (J₀, J_f, t) combinations including sgn-flip.
    for (J0, Jf) in ((1.0, -1.0), (-1.0, 1.0), (0.5, -2.0), (-0.7, 1.3))
        m_0 = XXZ1D(; J=J0, Δ=0.0)
        m_f = XXZ1D(; J=Jf, Δ=0.0)
        for t in (0.0, 1.0, 5.0)
            λ = QAtlas.fetch(m_f, _LE_RATE, Infinite(); initial=m_0, t=t)
            @test λ == 0.0
        end
    end
end

@testset "XXZ1D Δ=0 LoschmidtEcho — Δ ≠ 0 raises DomainError" begin
    # Final has Δ = 0.5
    @test_throws DomainError QAtlas.fetch(
        XXZ1D(; J=1.0, Δ=0.5), _LE_RATE, Infinite(); initial=_XX, t=1.0
    )
    # Initial has Δ = 0.5
    @test_throws DomainError QAtlas.fetch(
        _XX, _LE_RATE, Infinite(); initial=XXZ1D(; J=1.0, Δ=0.5), t=1.0
    )
    # Both have Δ ≠ 0
    @test_throws DomainError QAtlas.fetch(
        XXZ1D(; J=1.0, Δ=1.0), _LE_RATE, Infinite(); initial=XXZ1D(; J=1.0, Δ=-0.3), t=1.0
    )
end

@testset "XXZ1D Δ=0 LoschmidtEcho — flat-band edge cases (J = 0): λ ≡ 0" begin
    # Both flat-band → no dynamics → λ = 0.
    λ_both = QAtlas.fetch(
        XXZ1D(; J=0.0, Δ=0.0), _LE_RATE, Infinite(); initial=XXZ1D(; J=0.0, Δ=0.0), t=2.5
    )
    @test λ_both == 0.0

    # J_initial = 0, J_final ≠ 0: H_0 = 0, the GS manifold is
    # 2^N-degenerate.  For *any* number-eigenstate choice from that
    # manifold, |ψ₀⟩ is also an H_f eigenstate (number eigenstates of
    # the shared plane-wave basis), so |L(t)| = 1 and λ = 0.  The
    # implementation routes through the same `return 0.0` as every
    # other (J₀, J_f) combination, which is consistent with this
    # observation.
    λ_init_flat = QAtlas.fetch(
        XXZ1D(; J=1.0, Δ=0.0), _LE_RATE, Infinite(); initial=XXZ1D(; J=0.0, Δ=0.0), t=1.0
    )
    @test λ_init_flat == 0.0

    # J_initial ≠ 0, J_final = 0: H_f = 0, so e^{-iH_f t} = I and
    # L(t) = ⟨ψ₀|ψ₀⟩ = 1 trivially.  λ = 0.
    λ_final_flat = QAtlas.fetch(
        XXZ1D(; J=0.0, Δ=0.0), _LE_RATE, Infinite(); initial=XXZ1D(; J=1.0, Δ=0.0), t=1.0
    )
    @test λ_final_flat == 0.0
end

@testset "XXZ1D Δ=0 LoschmidtEcho — registry row landed" begin
    # Cross-check that XXZ_registry.jl declares the new triple.
    rows = QAtlas.implementation_status()
    matching = filter(rows) do r
        return r.model === XXZ1D &&
               r.quantity === LoschmidtRateFunction &&
               r.bc === Infinite
    end
    @test length(matching) == 1
    if length(matching) == 1
        @test matching[1].method === :free_fermion_analytic
        @test matching[1].reliability === :high
    end
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "XXZ1D Δ=0 LoschmidtEcho — verification cards" begin
    # The XX Loschmidt rate lambda(t) = 0 identically for any (J0, Jf, t):
    # both GS and evolved state are number eigenstates in the same plane-wave
    # basis, so the Loschmidt amplitude |<psi0|e^{-iHt}|psi0>| = 1 exactly.
    # The second closed form here is the free-fermion Slater-determinant
    # argument: det(D0^dag U_f(t) D0) = exp(i*phase) with |.| = 1.
    verify(
        XXZ1D(; J=1.0, Δ=0.0),
        LoschmidtRateFunction(),
        Infinite();
        route=:second_closed_form,
        fetch_kw=(; initial=XXZ1D(; J=0.5, Δ=0.0), t=1.0),
        independent=0.0,
        agree_within=1e-14,
        refs=["Free-fermion Slater det: |L(t)| = 1 for same-sign J quench => lambda = 0"],
    )

    verify(
        XXZ1D(; J=1.0, Δ=0.0),
        LoschmidtRateFunction(),
        Infinite();
        route=:limiting_case,
        fetch_kw=(; initial=XXZ1D(; J=1.0, Δ=0.0), t=2.5),
        independent=0.0,
        agree_within=1e-14,
        refs=["No-quench identity: H_initial = H_final => lambda(t) = 0 for all t"],
    )
end
