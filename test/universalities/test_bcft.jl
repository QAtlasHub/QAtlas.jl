using QAtlas, Test

@testset "BCFT — Ising Cardy boundary entropy log g (Phase 1, Affleck-Ludwig 1991)" begin
    m = BCFT()
    # Fixed boundary state |1⟩: g = 1/√2  ⟹  log g = -log(2)/2 = -0.34657...
    s_fixed = QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed)
    @test s_fixed ≈ -log(2) / 2
    @test s_fixed ≈ -0.34657359 atol=1e-7
    # Equivalent aliases for ±
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed_plus) ≈ s_fixed
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed_minus) ≈ s_fixed
    # Free boundary state |ε⟩: g = 1  ⟹  log g = 0
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:free) == 0.0
    # σ Cardy state: g = 1/√2 (same as fixed)
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:sigma) ≈ s_fixed
    # Default state = :fixed
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite()) ≈ s_fixed
    # g-theorem sanity: free is HIGHER (g=1) than fixed (g=1/√2 < 1)
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:free) >
        QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed)
end

@testset "BCFT — rejects unknown Cardy state (Phase 1)" begin
    m = BCFT()
    @test_throws DomainError QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:invalid)
    @test_throws DomainError QAtlas.fetch(
        m, ResidualEntropy(), Infinite(); state=:something
    )
end
