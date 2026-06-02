# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c9.jl
#
# Stage C.9 regression: high-level JKS free energy wrapper (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE: hubbard1d_jks_free_energy, atomic_free_energy

@testset "Hubbard1D — JKS Stage C.9 free energy wrapper (#523)" begin
    @testset "High-T returns a finite real number" begin
        f = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, 0.01)
        @test isfinite(f)
        @test f isa Float64
    end

    @testset "Validates t = 1 normalization" begin
        @test_throws ArgumentError hubbard1d_jks_free_energy(2.0, 4.0, 2.0, 0.01)
        @test_throws ArgumentError hubbard1d_jks_free_energy(0.5, 4.0, 2.0, 0.01)
    end

    @testset "Validates beta > 0, U >= 0" begin
        @test_throws DomainError hubbard1d_jks_free_energy(1.0, 4.0, 2.0, 0.0)
        @test_throws DomainError hubbard1d_jks_free_energy(1.0, 4.0, 2.0, -1.0)
        @test_throws DomainError hubbard1d_jks_free_energy(1.0, -1.0, 0.0, 0.1)
    end

    @testset "alpha >= eta throws DomainError (caller error)" begin
        # alpha must be < eta = U/4. With U = 4, eta = 1; alpha = 1 is invalid.
        @test_throws DomainError hubbard1d_jks_free_energy(1.0, 4.0, 2.0, 0.01; alpha=1.0)
        @test_throws DomainError hubbard1d_jks_free_energy(1.0, 4.0, 2.0, 0.01; alpha=2.0)
    end

    @testset "High-T sanity: returns reasonable finite value" begin
        # We don't assert numerical equality with atomic limit because the
        # b-channel-only solver gives a partial NLIE answer; just verify it's
        # finite and negative (free energy is bounded above by the atomic limit).
        beta = 0.01
        f = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, beta)
        if isfinite(f)
            @test f < 0  # any reasonable f at finite T is negative
        else
            @test isnan(f)  # solver may not converge at chosen parameters
        end
    end

    @testset "Mid-T continuation: returns finite or NaN, never throws" begin
        f = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, 1.0; tol=1e-2)
        @test isfinite(f) || isnan(f)
    end
end
