using Test
using QAtlas
using QAtlas: Universality, CHSHBound, Infinite

@testset "Universality(:QuantumMechanics) CHSHBound Tsirelson 1980 (#579)" begin
    @testset "Classical (CHSH 1969): S_max = 2" begin
        b = QAtlas.fetch(
            Universality(:QuantumMechanics), CHSHBound(), Infinite(); convention=:classical
        )
        @test b ≈ 2.0 atol = 1e-14
    end

    @testset "Quantum (Tsirelson 1980): S_max = 2 sqrt(2)" begin
        for cv in (:quantum,)
            b = QAtlas.fetch(
                Universality(:QuantumMechanics), CHSHBound(), Infinite(); convention=cv
            )
            @test b ≈ 2 * sqrt(2) atol = 1e-14
            @test b ≈ 2.828427124746190 atol = 1e-12
        end
        # Default convention is :quantum
        b = QAtlas.fetch(Universality(:QuantumMechanics), CHSHBound(), Infinite())
        @test b ≈ 2 * sqrt(2) atol = 1e-14
    end

    @testset "No-signalling / Popescu-Rohrlich: S_max = 4" begin
        for cv in (:no_signalling, :pr)
            b = QAtlas.fetch(
                Universality(:QuantumMechanics), CHSHBound(), Infinite(); convention=cv
            )
            @test b ≈ 4.0 atol = 1e-14
        end
    end

    @testset "Hierarchy: classical < quantum < no-signalling" begin
        b_c = QAtlas.fetch(
            Universality(:QuantumMechanics), CHSHBound(), Infinite(); convention=:classical
        )
        b_q = QAtlas.fetch(
            Universality(:QuantumMechanics), CHSHBound(), Infinite(); convention=:quantum
        )
        b_n = QAtlas.fetch(
            Universality(:QuantumMechanics),
            CHSHBound(),
            Infinite();
            convention=:no_signalling,
        )
        @test b_c < b_q < b_n
        # quantum/classical = sqrt(2): the Bell-inequality violation factor
        @test b_q / b_c ≈ sqrt(2) atol = 1e-14
        # quantum/no_signalling = sqrt(2)/2: the quantum bound is sqrt(2)/2 of PR
        @test b_q / b_n ≈ sqrt(2) / 2 atol = 1e-14
    end

    @testset "Argument validation" begin
        @test_throws ArgumentError QAtlas.fetch(
            Universality(:QuantumMechanics), CHSHBound(), Infinite(); convention=:bogus
        )
    end
end
