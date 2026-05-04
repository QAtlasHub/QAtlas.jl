# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: XX (Δ = 0) free-fermion quench observables for
# `XXZ1D` at `Infinite()`.
#
# Validates `fetch(::XXZ1D, ::LoschmidtEcho{:rate}, ::Infinite;
#                  initial::XXZ1D, t::Real)` against the closed-form
# analytical result derived in `src/models/quantum/XXZ/XXZ_xx_quench.jl`:
#
#   λ(t) = 0     for sgn J_initial == sgn J_final  (same Fermi sea)
#   λ(t) = ∞     for sgn J_initial != sgn J_final  (Anderson orthog.)
#   λ(t) = NaN   for the mixed flat-band cases (J = 0 on one side only)
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
using Logging: with_logger, NullLogger

const _XX = XXZ1D(; J=1.0, Δ=0.0)
const _LE_RATE = LoschmidtEcho(; mode=:rate)

# Helper to silence the deliberate `@warn` calls in the orthogonality /
# flat-band branches without losing real failures.
_silent(f) =
    with_logger(NullLogger()) do
        f()
    end

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

@testset "XXZ1D Δ=0 LoschmidtEcho — sign-flip quench (orthogonality)" begin
    # sgn J flips → complementary Fermi seas → Anderson orthogonality;
    # rate is +∞ for every t.  Suppress the deliberate @warn.
    m_0 = XXZ1D(; J=1.0, Δ=0.0)
    m_f = XXZ1D(; J=-1.0, Δ=0.0)
    for t in (0.0, 1.0, 5.0)
        λ = _silent(() -> QAtlas.fetch(m_f, _LE_RATE, Infinite(); initial=m_0, t=t))
        @test λ === Inf
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

@testset "XXZ1D Δ=0 LoschmidtEcho — flat-band edge cases" begin
    # Both flat-band → no dynamics → λ = 0 (no warning).
    λ_both = QAtlas.fetch(
        XXZ1D(; J=0.0, Δ=0.0), _LE_RATE, Infinite(); initial=XXZ1D(; J=0.0, Δ=0.0), t=2.5
    )
    @test λ_both == 0.0

    # One-sided flat band → degenerate; λ = NaN with a @warn.
    λ_init_flat = _silent(
        () -> QAtlas.fetch(
            XXZ1D(; J=1.0, Δ=0.0),
            _LE_RATE,
            Infinite();
            initial=XXZ1D(; J=0.0, Δ=0.0),
            t=1.0,
        ),
    )
    @test isnan(λ_init_flat)

    λ_final_flat = _silent(
        () -> QAtlas.fetch(
            XXZ1D(; J=0.0, Δ=0.0),
            _LE_RATE,
            Infinite();
            initial=XXZ1D(; J=1.0, Δ=0.0),
            t=1.0,
        ),
    )
    @test isnan(λ_final_flat)
end

@testset "XXZ1D Δ=0 LoschmidtEcho — registry row landed" begin
    # Cross-check that XXZ_registry.jl declares the new triple.
    rows = QAtlas.implementation_status()
    matching = filter(rows) do r
        r.model === XXZ1D && r.quantity === LoschmidtEcho{:rate} && r.bc === Infinite
    end
    @test length(matching) == 1
    if length(matching) == 1
        @test matching[1].method === :free_fermion_analytic
        @test matching[1].reliability === :high
    end
end
