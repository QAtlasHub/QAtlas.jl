using Test
using QAtlas

@testset "Universality(:RMT) WignerSemicircleMoment" begin
    u = QAtlas.Universality(:RMT)

    @testset "Catalan number even moments" begin
        catalan_ref = (1.0, 1.0, 2.0, 5.0, 14.0, 42.0, 132.0, 429.0)
        for (k, ref) in enumerate(catalan_ref)
            n = 2 * (k - 1)
            m_n = QAtlas.fetch(u, QAtlas.WignerSemicircleMoment(), QAtlas.Infinite(); n=n)
            @test m_n ≈ ref atol = 1e-12
        end
    end

    @testset "Odd moments vanish" begin
        for n in (1, 3, 5, 7, 9, 11, 13)
            @test QAtlas.fetch(
                u, QAtlas.WignerSemicircleMoment(), QAtlas.Infinite(); n=n
            ) == 0.0
        end
    end

    @testset "Catalan recursion C_{k+1} = sum_i C_i C_{k-i}" begin
        moments = [
            QAtlas.fetch(u, QAtlas.WignerSemicircleMoment(), QAtlas.Infinite(); n=2k) for
            k in 0:5
        ]
        for k in 1:5
            recursion = sum(moments[i + 1] * moments[k - i] for i in 0:(k - 1))
            @test moments[k + 1] ≈ recursion atol = 1e-10
        end
    end

    @testset "Without explicit Infinite() also dispatches" begin
        @test QAtlas.fetch(u, QAtlas.WignerSemicircleMoment(); n=4) == 2.0
        @test QAtlas.fetch(u, QAtlas.WignerSemicircleMoment(); n=6) == 5.0
    end

    @testset "DomainError for negative n" begin
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.WignerSemicircleMoment(), QAtlas.Infinite(); n=-1
        )
        @test_throws DomainError QAtlas.fetch(
            u, QAtlas.WignerSemicircleMoment(), QAtlas.Infinite(); n=-10
        )
    end

    @testset "Large n still works (BigInt internally)" begin
        m_40 = QAtlas.fetch(u, QAtlas.WignerSemicircleMoment(), QAtlas.Infinite(); n=40)
        @test m_40 ≈ 6564120420.0 atol = 1e-2
    end
end
