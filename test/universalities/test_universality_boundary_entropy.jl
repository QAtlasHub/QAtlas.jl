using Test
using QAtlas
using QAtlas: Universality, BoundaryEntropy, Infinite

@testset "Universality(:Ising) BoundaryEntropy Affleck-Ludwig log g (#580)" begin
    @testset "Fixed boundary (identity / epsilon Cardy states): log g = -(1/2) log 2" begin
        for state in (:identity, :epsilon, :fixed_up, :fixed_down)
            g_log = QAtlas.fetch(
                Universality(:Ising), BoundaryEntropy(), Infinite(); boundary_state=state
            )
            @test g_log ≈ -log(2) / 2 atol = 1e-14
            @test exp(g_log) ≈ 1 / sqrt(2) atol = 1e-14
        end
    end

    @testset "Free boundary (sigma Cardy state): log g = 0" begin
        for state in (:sigma, :free)
            g_log = QAtlas.fetch(
                Universality(:Ising), BoundaryEntropy(), Infinite(); boundary_state=state
            )
            @test g_log ≈ 0.0 atol = 1e-14
            @test exp(g_log) ≈ 1.0 atol = 1e-14
        end
    end

    @testset "g-theorem ordering: free > fixed (g decreases under RG)" begin
        g_free = exp(
            QAtlas.fetch(
                Universality(:Ising), BoundaryEntropy(), Infinite(); boundary_state=:free
            ),
        )
        g_fixed = exp(
            QAtlas.fetch(
                Universality(:Ising),
                BoundaryEntropy(),
                Infinite();
                boundary_state=:fixed_up,
            ),
        )
        @test g_free > g_fixed
        # The g ratio g_free / g_fixed = sqrt(2) (universal).
        @test g_free / g_fixed ≈ sqrt(2) atol = 1e-14
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:Ising), BoundaryEntropy(), Infinite(); boundary_state=:bogus
        )
    end
end
