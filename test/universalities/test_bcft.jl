using QAtlas, Test

@testset "BCFT — Ising Cardy boundary entropy log g (Phase 1, Affleck-Ludwig 1991)" begin
    m = BCFT()
    # Physical fixed-spin boundaries |±⟩ ∝ (|1⟩ ± |ε⟩)/√2: g = 1/√2 ⟹ log g = -log(2)/2
    s_fixed = QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed)
    @test s_fixed ≈ -log(2) / 2
    @test s_fixed ≈ -0.34657359 atol=1e-7
    # Equivalent aliases for ±
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed_plus) ≈ s_fixed
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed_minus) ≈ s_fixed
    # Bare primary Cardy states |1⟩ and |ε⟩: g = 1/√2 each
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:identity) ≈ s_fixed
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:vacuum) ≈ s_fixed
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:epsilon) ≈ s_fixed
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:energy) ≈ s_fixed
    # Physical free boundary ≡ |σ⟩ Cardy state: g = S_{σ,1}/√S_{1,1} = (1/√2)/√(1/2) = 1
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:free) == 0.0
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:sigma) == 0.0
    # :free and :sigma are synonyms
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:free) ==
        QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:sigma)
    # Default state = :fixed
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite()) ≈ s_fixed
    # g-theorem sanity: free (g=1) is HIGHER than fixed (g=1/√2 < 1)
    @test QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:free) >
        QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed)
end

@testset "BCFT — exp(s) = g round-trip cross-check (Cardy 1989)" begin
    m = BCFT()
    # Recover g from log g and confront with the closed-form Ising Cardy values.
    s_fixed = QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:fixed)
    s_free  = QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:free)
    @test exp(s_fixed) ≈ 1 / sqrt(2) atol=1e-12
    @test exp(s_free)  == 1.0
    # Aliases collapse to the same g.
    for st in (:fixed_plus, :fixed_minus, :identity, :vacuum, :epsilon, :energy)
        s = QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=st)
        @test exp(s) ≈ 1 / sqrt(2) atol=1e-12
    end
    # |σ⟩ Cardy state is the physical free boundary: g_σ = 1.
    @test exp(QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:sigma)) == 1.0
end

@testset "BCFT — rejects unknown Cardy state (Phase 1)" begin
    m = BCFT()
    @test_throws DomainError QAtlas.fetch(m, ResidualEntropy(), Infinite(); state=:invalid)
    @test_throws DomainError QAtlas.fetch(
        m, ResidualEntropy(), Infinite(); state=:something
    )
end
