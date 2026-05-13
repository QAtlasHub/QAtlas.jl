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
