using Test
using QAtlas
using QAtlas: TFIM, LiebRobinsonVelocity, Infinite

@testset "TFIM LiebRobinsonVelocity (#579 Phase 1)" begin
    @testset "v_LR = 2 min(|J|, |h|) (tight free-fermion bound)" begin
        for J in (0.5, 1.0, 2.5, 10.0), h in (0.5, 1.0, 2.0, 5.0)
            m = TFIM(; J=J, h=h)
            v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
            @test v ≈ 2 * min(abs(J), abs(h)) atol = 1e-14
        end
    end

    @testset "v_LR = 2J at criticality h = J" begin
        for J in (0.5, 1.0, 2.5, 10.0)
            m = TFIM(; J=J, h=J)
            v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
            @test v ≈ 2 * abs(J) atol = 1e-14
        end
    end

    @testset "Vanishing v_LR at h = 0 or J = 0 (no quantum dynamics)" begin
        # h = 0: classical Ising, H = -J σz σz commutes with σz at each
        # site, no propagation.
        for J in (0.5, 1.0, 2.5)
            m = TFIM(; J=J, h=0.0)
            v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
            @test v ≈ 0.0 atol = 1e-14
        end
        # J = 0: decoupled site spins under transverse field, no
        # interaction means no propagation.
        for h in (0.5, 1.0, 2.5)
            m = TFIM(; J=0.0, h=h)
            v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
            @test v ≈ 0.0 atol = 1e-14
        end
    end

    @testset "Negative J / h handled via abs" begin
        m = TFIM(; J=-1.5, h=0.7)
        v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
        @test v ≈ 2 * min(1.5, 0.7) atol = 1e-14  # = 1.4
        m2 = TFIM(; J=2.0, h=-0.8)
        v2 = QAtlas.fetch(m2, LiebRobinsonVelocity(), Infinite())
        @test v2 ≈ 2 * min(2.0, 0.8) atol = 1e-14  # = 1.6
    end

    @testset "J / h kwargs override model values" begin
        m = TFIM(; J=1.0, h=0.5)
        v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite(); J=4.0, h=10.0)
        @test v ≈ 2 * min(4.0, 10.0) atol = 1e-14  # = 8.0
    end
end
