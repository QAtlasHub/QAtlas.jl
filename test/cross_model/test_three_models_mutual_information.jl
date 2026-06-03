using Test
using QAtlas

@testset "Three-model MutualInformation wrappers" begin
    @testset "TFIM h = J critical: matches Universality(:Ising) with a = 1/2" begin
        for J in (1.0, 2.0)
            m = QAtlas.TFIM(; h=J, J=J)
            for (ℓ_A, ℓ_B) in ((5.0, 10.0), (10.0, 20.0), (3.0, 30.0))
                for β in (Inf, 5.0, 1.0)
                    mi = QAtlas.fetch(
                        m,
                        QAtlas.MutualInformation(),
                        QAtlas.Infinite();
                        ℓ_A=ℓ_A,
                        ℓ_B=ℓ_B,
                        beta=β,
                    )
                    ref = QAtlas.fetch(
                        QAtlas.Universality(:Ising),
                        QAtlas.MutualInformation(),
                        QAtlas.Infinite();
                        ℓ_A=ℓ_A,
                        ℓ_B=ℓ_B,
                        beta=β,
                        a=1//2,
                    )
                    @test mi ≈ ref atol = 1e-12
                end
            end
        end
        m_off = QAtlas.TFIM(; h=0.5, J=1.0)
        @test_throws DomainError QAtlas.fetch(
            m_off, QAtlas.MutualInformation(), QAtlas.Infinite(); ℓ_A=5.0, ℓ_B=5.0
        )
    end

    @testset "HaldaneShastry: matches Universality(:Heisenberg)" begin
        m = QAtlas.HaldaneShastry()
        for (ℓ_A, ℓ_B) in ((5.0, 10.0), (10.0, 20.0))
            for β in (Inf, 5.0)
                mi = QAtlas.fetch(
                    m,
                    QAtlas.MutualInformation(),
                    QAtlas.Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                    beta=β,
                )
                ref = QAtlas.fetch(
                    QAtlas.Universality(:Heisenberg),
                    QAtlas.MutualInformation(),
                    QAtlas.Infinite();
                    ℓ_A=ℓ_A,
                    ℓ_B=ℓ_B,
                    beta=β,
                )
                @test mi ≈ ref atol = 1e-12
            end
        end
    end

    @testset "XXZ1D critical regime: matches Universality(:XY) or :Heisenberg" begin
        for Δ in (0.0, 0.5, -0.5, 0.9)
            m = QAtlas.XXZ1D(; J=1.0, Δ=Δ)
            mi = QAtlas.fetch(
                m, QAtlas.MutualInformation(), QAtlas.Infinite(); ℓ_A=8.0, ℓ_B=12.0
            )
            ref = QAtlas.fetch(
                QAtlas.Universality(:XY),
                QAtlas.MutualInformation(),
                QAtlas.Infinite();
                ℓ_A=8.0,
                ℓ_B=12.0,
            )
            @test mi ≈ ref atol = 1e-12
        end
        m_isotropic = QAtlas.XXZ1D(; J=1.0, Δ=1.0)
        mi_iso = QAtlas.fetch(
            m_isotropic, QAtlas.MutualInformation(), QAtlas.Infinite(); ℓ_A=8.0, ℓ_B=12.0
        )
        ref_iso = QAtlas.fetch(
            QAtlas.Universality(:Heisenberg),
            QAtlas.MutualInformation(),
            QAtlas.Infinite();
            ℓ_A=8.0,
            ℓ_B=12.0,
        )
        @test mi_iso ≈ ref_iso atol = 1e-12
    end

    @testset "XXZ1D off-critical (|Δ| > 1) throws DomainError" begin
        for Δ in (1.5, -1.5, 2.0, -2.0)
            m = QAtlas.XXZ1D(; J=1.0, Δ=Δ)
            @test_throws DomainError QAtlas.fetch(
                m, QAtlas.MutualInformation(), QAtlas.Infinite(); ℓ_A=8.0, ℓ_B=12.0
            )
        end
    end
end
