# ─────────────────────────────────────────────────────────────────────────────
# Universality test: ConformalBootstrap — 3D Ising Δ_σ, Δ_ε
# (KPSD 2014 methodology; precise values KPSD-Vichi 2016 "Precision
# Islands" / arXiv:1603.04436; Simmons-Duffin 2017 cross-check).
#
# Verifies:
#   * Δ_σ and Δ_ε match the KPSD-Vichi 2016 "Precision Islands"
#     reference values (within the 2016 (10) uncertainty bars).
#   * Sanity ordering Δ_σ < Δ_ε.
#   * Derived 3D Ising exponents ν = 1/(3 − Δ_ε) ≈ 0.62997 and
#     η = 2Δ_σ − 1 ≈ 0.03629.
#   * Default `field` is `:σ`.
#   * DomainError on Phase-2 operators (`:T`, `:O`, ...).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "ConformalBootstrap — 3D Ising Δ_σ, Δ_ε (Phase 1, KPSD-Vichi 2016)" begin
    m = ConformalBootstrap()
    Δ_σ = QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:σ)
    Δ_ε = QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:ε)
    @test Δ_σ ≈ 0.51814894 atol=1e-8
    @test Δ_ε ≈ 1.41262528 atol=1e-8
    # Sanity: σ < ε (correlation length exponent ν = 1/(d − Δ_ε) > 0 in d=3)
    @test Δ_σ < Δ_ε
    # Compute ν = 1/(3 − Δ_ε) ≈ 0.62997 (also a famous 3D Ising exponent)
    ν = 1 / (3 - Δ_ε)
    @test ν ≈ 0.62997 atol=1e-4
    # η = 2Δ_σ − (d − 2) = 2Δ_σ − 1 ≈ 0.03629
    η = 2 * Δ_σ - 1
    @test η ≈ 0.03629 atol=1e-4
    # Default field = :σ
    @test QAtlas.fetch(m, ConformalWeights(), Infinite()) ≈ 0.51814894 atol=1e-8
end

@testset "ConformalBootstrap — 3D Ising unitarity / relevance cross-check" begin
    m = ConformalBootstrap()
    Δ_σ = QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:σ)
    Δ_ε = QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:ε)
    # 3D scalar unitarity bound: Δ ≥ (d − 2)/2 = 1/2. Both σ and ε must satisfy it.
    @test Δ_σ ≥ 0.5
    @test Δ_ε ≥ 0.5
    # Relevance: ε must be relevant (Δ_ε < d = 3), otherwise tuning T to T_c
    # would not be possible (1 relevant Z_2-even operator at the Ising fixed point).
    @test Δ_ε < 3
    # η ≥ 0 (correlation function decays at least as fast as the free-field rate).
    @test 2 * Δ_σ - 1 ≥ 0
    # ν > 0 (correlation length diverges at the critical point).
    @test 1 / (3 - Δ_ε) > 0
end

@testset "ConformalBootstrap — rejects unknown field (Phase 2)" begin
    m = ConformalBootstrap()
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:T)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:O)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "ConformalBootstrap — verification cards" begin
    # 3D Ising conformal bootstrap (Kos-Poland-Simmons-Duffin-Vichi 2016
    # "Precision Islands"): Δ_σ ≈ 0.5181489, Δ_ε ≈ 1.412625 (literature).
    verify(
        ConformalBootstrap(),
        ConformalWeights(),
        Infinite();
        route=:literature_value,
        fetch_kw=(; field=:σ),
        independent=0.5181489,
        agree_within=1e-5,
        refs=["KPSD-Vichi 2016 Precision Islands: 3D Ising Δ_σ ≈ 0.5181489"],
    )
    verify(
        ConformalBootstrap(),
        ConformalWeights(),
        Infinite();
        route=:literature_value,
        fetch_kw=(; field=:ε),
        independent=1.412625,
        agree_within=1e-4,
        refs=["KPSD-Vichi 2016 Precision Islands: 3D Ising Δ_ε ≈ 1.412625"],
    )
end
# ── additional verification cards (#381 batch 6) ─────────────────────────
@testset "ConformalBootstrap — 3D Ising Δσ (#381 batch 6)" begin
    # Numerical conformal bootstrap precision value for the 3D Ising
    # spin operator scaling dimension (Kos-Poland-Simmons-Duffin 2014
    # PRD 86 025022; refined to ~10 digits by Simmons-Duffin 2017):
    # Δσ ≈ 0.5181489.
    verify(
        ConformalBootstrap(),
        ConformalWeights(),
        Infinite();
        route=:literature_value,
        independent=0.5181489,
        agree_within=1e-6,
        refs=["Kos-Poland-Simmons-Duffin 2014 PRD 86 025022: 3D Ising spin op Δσ ≈ 0.5181489 (numerical bootstrap)"],
        fetch_kw=(; r=1, s=1),
    )
end

