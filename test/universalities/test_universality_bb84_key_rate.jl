using Test
using QAtlas

@testset "Universality BB84KeyRate (Shor-Preskill 2000)" begin
    u = QAtlas.Universality(:Ising)

    @testset "Boundary values" begin
        @test QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.0) ≈ 1.0 atol =
            1e-12
        @test QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.5) ≈ -1.0 atol =
            1e-12
    end

    @testset "Numerical values across QBER" begin
        r_05 = QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.05)
        h2_05 = -0.05 * log2(0.05) - 0.95 * log2(0.95)
        @test r_05 ≈ 1 - 2 * h2_05 atol = 1e-12

        r_threshold = QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.11)
        @test r_threshold < 0.001
    end

    @testset "Unconditional-security threshold e* ≈ 0.110028" begin
        r_below = QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.10)
        r_at = QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.110028)
        r_above = QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.13)
        @test r_below > 0
        @test abs(r_at) < 1e-4
        @test r_above < 0
    end

    @testset "Monotonically decreasing in qber on [0, 0.5]" begin
        prev = 2.0
        for q in 0.01:0.04:0.45
            r = QAtlas.fetch(u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=q)
            @test r < prev
            prev = r
        end
    end

    @testset "Class-independent" begin
        q = 0.05
        vals = [
            QAtlas.fetch(
                QAtlas.Universality(sym), QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=q
            ) for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
        ]
        @test all(isapprox(v, first(vals); atol=1e-12) for v in vals)
    end

    @testset "DomainError on QBER outside [0, 0.5]" begin
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=-0.01
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=0.51
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.BB84KeyRate(), QAtlas.Infinite(); qber=1.0
        )
    end
end
