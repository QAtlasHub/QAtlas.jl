using QAtlas, Test

@testset "FibonacciAnyons — TopologicalEntanglementEntropy (Phase 1)" begin
    γ = QAtlas.fetch(FibonacciAnyons(), TopologicalEntanglementEntropy(), Infinite())
    # Closed-form check
    @test γ ≈ 0.5 * log(2 + MathConstants.golden)
    # Numerical value (φ ≈ 1.6180, 2+φ ≈ 3.6180, log = 1.2862..., ÷ 2 ≈ 0.6431)
    @test γ ≈ 0.6429653906 atol=1e-9
    # Cross-check against direct literal-sqrt formula (independent of MathConstants.golden)
    @test γ ≈ 0.5 * log(2 + (1 + sqrt(5)) / 2) atol=1e-15
    # Equivalently γ = log √(1 + φ²) and log √(φ+2)
    φ = MathConstants.golden
    @test γ ≈ log(sqrt(1 + φ^2))
    @test γ ≈ log(sqrt(φ + 2))
    # Larger than Z_2 toric-code value γ = log 2 ≈ 0.6931? Actually log(2) > log√(φ+2)
    # since 2 = √4 > √(φ+2) = √3.618.  So γ_Fib < γ_Z2.
    @test γ < log(2)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "FibonacciAnyons — verification cards" begin
    phi = (1 + sqrt(5)) / 2
    verify(
        FibonacciAnyons(),
        TopologicalEntanglementEntropy(),
        Infinite();
        route=:second_closed_form,
        independent=0.5 * log(2 + phi),
        agree_within=1e-9,
        refs=["Fibonacci TQFT: gamma = (1/2) log(2 + phi), phi = golden ratio"],
    )
end
# ── additional verification cards (#381 batch 3) ─────────────────────────
@testset "FibonacciAnyons — TopologicalEntanglementEntropy (#381 batch 3)" begin
    # Fibonacci anyon model: only two anyons {1, τ} with quantum
    # dimensions d_1 = 1, d_τ = φ = (1+√5)/2 (golden ratio).
    # Total quantum dimension D = √(1 + φ²) = √((5+√5)/2).
    # Kitaev-Preskill / Levin-Wen TEE = log D = (1/2) log(1 + φ²).
    phi = (1 + sqrt(5)) / 2
    verify(
        FibonacciAnyons(),
        TopologicalEntanglementEntropy(),
        Infinite();
        route=:second_closed_form,
        independent=0.5 * log(1 + phi^2),
        agree_within=1e-12,
        refs=["Kitaev-Preskill 2006 / Levin-Wen 2006: TEE = log D = (1/2) log(1+φ²) for Fibonacci anyons (φ = golden ratio)"],
    )
end

