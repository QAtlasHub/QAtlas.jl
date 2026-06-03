using Test
using QAtlas

@testset "Three-model LogarithmicNegativity wrappers (CC-Tonni 2012)" begin
    @testset "TFIM h = J: matches Universality(:Ising)" begin
        for J in (1.0, 2.0)
            m = QAtlas.TFIM(; h=J, J=J)
            for (ℓ_A, ℓ_B) in ((5.0, 10.0), (10.0, 20.0), (3.0, 30.0))
                e = QAtlas.fetch(
                    m, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
                )
                ref = (0.5 / 4) * log(ℓ_A * ℓ_B / (ℓ_A + ℓ_B))
                @test e ≈ ref atol = 1e-12
            end
        end
        m_off = QAtlas.TFIM(; h=0.5, J=1.0)
        @test_throws DomainError QAtlas.fetch(
            m_off, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=5.0, ℓ_B=5.0
        )
    end

    @testset "HaldaneShastry: c = 1" begin
        m = QAtlas.HaldaneShastry()
        for (ℓ_A, ℓ_B) in ((5.0, 10.0), (10.0, 20.0))
            e = QAtlas.fetch(
                m, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
            )
            ref = (1.0 / 4) * log(ℓ_A * ℓ_B / (ℓ_A + ℓ_B))
            @test e ≈ ref atol = 1e-12
        end
    end

    @testset "XXZ1D critical regime: c = 1" begin
        for Δ in (0.0, 0.5, -0.5, 0.9, 1.0)
            m = QAtlas.XXZ1D(; J=1.0, Δ=Δ)
            for (ℓ_A, ℓ_B) in ((8.0, 12.0), (15.0, 25.0))
                e = QAtlas.fetch(
                    m, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
                )
                ref = (1.0 / 4) * log(ℓ_A * ℓ_B / (ℓ_A + ℓ_B))
                @test e ≈ ref atol = 1e-12
            end
        end
    end

    @testset "XXZ1D off-critical (|Δ| > 1) throws DomainError" begin
        for Δ in (1.5, -1.5, 2.0)
            m = QAtlas.XXZ1D(; J=1.0, Δ=Δ)
            @test_throws DomainError QAtlas.fetch(
                m, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=8.0, ℓ_B=12.0
            )
        end
    end

    @testset "Central-charge ratio TFIM/HS = 1/2" begin
        ℓ_A, ℓ_B = 10.0, 15.0
        m_tfim = QAtlas.TFIM(; h=1.0, J=1.0)
        m_hs = QAtlas.HaldaneShastry()
        e_tfim = QAtlas.fetch(
            m_tfim, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
        )
        e_hs = QAtlas.fetch(
            m_hs, QAtlas.LogarithmicNegativity(), QAtlas.Infinite(); ℓ_A=ℓ_A, ℓ_B=ℓ_B
        )
        @test e_tfim / e_hs ≈ 0.5 atol = 1e-12
    end
end
