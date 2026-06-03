using Test
using QAtlas

@testset "Universality BekensteinBound (Bekenstein 1981)" begin
    @testset "S_max = 2π R E" begin
        for sym in (:Ising, :Heisenberg, :XY, :Potts3, :Potts4)
            u = QAtlas.Universality(sym)
            for R in (0.1, 1.0, 10.0, 100.0)
                for E in (0.5, 1.0, 5.0, 50.0)
                    s = QAtlas.fetch(
                        u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=R, E=E
                    )
                    @test s ≈ 2π * R * E atol = 1e-10
                end
            end
        end
    end

    @testset "Class-independent at fixed (R, E)" begin
        R, E = 3.0, 4.0
        vals = [
            QAtlas.fetch(
                QAtlas.Universality(sym),
                QAtlas.BekensteinBound(),
                QAtlas.Infinite();
                R=R,
                E=E,
            ) for sym in (:Ising, :Heisenberg, :XY)
        ]
        @test all(isapprox(v, first(vals); atol=1e-12) for v in vals)
    end

    @testset "Bilinear: doubling R doubles S; doubling E doubles S" begin
        u = QAtlas.Universality(:Heisenberg)
        s_base = QAtlas.fetch(u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=1.0, E=1.0)
        s_2R = QAtlas.fetch(u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=2.0, E=1.0)
        s_2E = QAtlas.fetch(u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=1.0, E=2.0)
        s_4 = QAtlas.fetch(u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=2.0, E=2.0)
        @test s_2R ≈ 2 * s_base atol = 1e-12
        @test s_2E ≈ 2 * s_base atol = 1e-12
        @test s_4 ≈ 4 * s_base atol = 1e-12
    end

    @testset "DomainError on non-positive R or E" begin
        u = QAtlas.Universality(:Ising)
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=0.0, E=1.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=-1.0, E=1.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=1.0, E=0.0
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.BekensteinBound(), QAtlas.Infinite(); R=1.0, E=-2.0
        )
    end
end
