using Test
using QAtlas
using QAtlas: TFIM, LiebRobinsonVelocity, Infinite

@testset "TFIM LiebRobinsonVelocity (#579 Phase 1)" begin
    @testset "v_LR = 2|J| independent of h" begin
        for J in (0.5, 1.0, 2.5, 10.0), h in (0.0, 0.5, 1.0, 2.0, 5.0)
            m = TFIM(; J=J, h=h)
            v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
            @test v ≈ 2 * abs(J) atol = 1e-14
        end
    end

    @testset "Negative J: |J| handled" begin
        m = TFIM(; J=-1.5, h=0.3)
        v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite())
        @test v ≈ 3.0 atol = 1e-14
    end

    @testset "J kwarg overrides model.J" begin
        m = TFIM(; J=1.0, h=0.5)
        v = QAtlas.fetch(m, LiebRobinsonVelocity(), Infinite(); J=4.0)
        @test v ≈ 8.0 atol = 1e-14
    end
end
