using Test
using QAtlas

@testset "TFIM BoundaryEntropy at criticality (h = J)" begin
    @testset "Critical point: matches Universality(:Ising)" begin
        for J in (0.5, 1.0, 2.0)
            m = QAtlas.TFIM(; h=J, J=J)
            for bs in (:identity, :fixed_up, :fixed_down, :epsilon, :sigma, :free)
                g_model = QAtlas.fetch(
                    m, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=bs
                )
                g_univ = QAtlas.fetch(
                    QAtlas.Universality(:Ising),
                    QAtlas.BoundaryEntropy(),
                    QAtlas.Infinite();
                    boundary_state=bs,
                )
                @test g_model ≈ g_univ atol = 1e-12
            end
        end
    end

    @testset "Specific values at h = J = 1" begin
        m = QAtlas.TFIM(; h=1.0, J=1.0)
        @test QAtlas.fetch(
            m, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:identity
        ) ≈ -log(2) / 2 atol = 1e-12
        @test QAtlas.fetch(
            m, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:epsilon
        ) ≈ -log(2) / 2 atol = 1e-12
        @test QAtlas.fetch(
            m, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:sigma
        ) ≈ 0.0 atol = 1e-12
        @test QAtlas.fetch(
            m, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:free
        ) ≈ 0.0 atol = 1e-12
    end

    @testset "Critical at |h| = |J| (negative J or h)" begin
        m1 = QAtlas.TFIM(; h=2.0, J=-2.0)
        @test QAtlas.fetch(
            m1, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:sigma
        ) ≈ 0.0 atol = 1e-12
        m2 = QAtlas.TFIM(; h=-3.5, J=3.5)
        @test QAtlas.fetch(
            m2, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:identity
        ) ≈ -log(2) / 2 atol = 1e-12
    end

    @testset "DomainError off-critical (|h| != |J|)" begin
        for (h, J) in ((0.5, 1.0), (2.0, 1.0), (0.0, 1.0), (1.0, 0.5))
            m = QAtlas.TFIM(; h=h, J=J)
            @test_throws DomainError QAtlas.fetch(
                m, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:identity
            )
        end
    end

    @testset "Bad boundary_state propagates from Universality" begin
        m = QAtlas.TFIM(; h=1.0, J=1.0)
        @test_throws ArgumentError QAtlas.fetch(
            m, QAtlas.BoundaryEntropy(), QAtlas.Infinite(); boundary_state=:nonsense
        )
    end
end
