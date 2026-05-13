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

@testset "ConformalBootstrap — rejects unknown field (Phase 2)" begin
    m = ConformalBootstrap()
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:T)
    @test_throws DomainError QAtlas.fetch(m, ConformalWeights(), Infinite(); field=:O)
end
