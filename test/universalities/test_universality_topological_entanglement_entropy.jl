using Test
using QAtlas
using QAtlas: Universality, TopologicalEntanglementEntropy, Infinite, ToricCode

@testset "Universality(:TopologicalOrder) TopologicalEntanglementEntropy Kitaev-Preskill 2006 (#580)" begin
    @testset "Abelian Z_2 (Toric Code, 4 anyons of d=1): gamma = log 2" begin
        γ = QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=[1.0, 1.0, 1.0, 1.0],
        )
        @test γ ≈ log(2) atol = 1e-14
    end

    @testset "Ising anyon (3 sectors d=1,sqrt(2),1): gamma = log 2" begin
        γ = QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=[1.0, sqrt(2.0), 1.0],
        )
        @test γ ≈ log(2) atol = 1e-14
    end

    @testset "Fibonacci (2 sectors, golden ratio): gamma = (1/2) log(1 + phi^2)" begin
        ϕ = (1 + sqrt(5)) / 2
        γ = QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=[1.0, ϕ],
        )
        @test γ ≈ 0.5 * log(1 + ϕ^2) atol = 1e-14
    end

    @testset "Cross-check with existing ToricCode model dispatch" begin
        # Existing model-side fetch should agree (γ = log 2).
        γ_model = QAtlas.fetch(ToricCode(), TopologicalEntanglementEntropy(), Infinite())
        γ_univ = QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=[1.0, 1.0, 1.0, 1.0],
        )
        @test γ_model ≈ γ_univ atol = 1e-14
    end

    @testset "Trivial state (single trivial sector): gamma = 0" begin
        γ = QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=[1.0],
        )
        @test γ ≈ 0.0 atol = 1e-14
    end

    @testset "D_n discrete gauge theory: gamma = log n" begin
        # Z_n gauge theory has n^2 abelian anyons all with d_a = 1,
        # so D = sqrt(n^2) = n -> gamma = log n.
        for n in (2, 3, 4, 5, 10)
            dims = ones(n^2)
            γ = QAtlas.fetch(
                Universality(:TopologicalOrder),
                TopologicalEntanglementEntropy(),
                Infinite();
                quantum_dimensions=dims,
            )
            @test γ ≈ log(n) atol = 1e-14
        end
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=Float64[],
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=[1.0, 0.0, 1.0],
        )
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:TopologicalOrder),
            TopologicalEntanglementEntropy(),
            Infinite();
            quantum_dimensions=[1.0, -2.0, 1.0],
        )
    end
end
