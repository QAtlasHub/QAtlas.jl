# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_a.jl
#
# Stage A regression for the Juttner-Klumper-Suzuki NLIE module skeleton
# (#523). Tests:
#   - atomic_free_energy reproduces -T ln 4 at U = 0 (high-T baseline)
#   - half-filling specialization matches the general form at mu = U/2
#   - low-T limit at half filling: f -> -U/2 - T ln 2
#   - Stage-B stubs error loudly when called
#   - JKSGrid validates its inputs

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    atomic_free_energy,
    atomic_free_energy_half_filling,
    JKSGrid,
    jks_kernel_K_n,
    jks_driving_b,
    jks_driving_c,
    jks_driving_cbar

@testset "Hubbard1D — JKS Stage A skeleton (#523)" begin
    @testset "Atomic limit at U = 0 reproduces -T ln 4" begin
        for beta in (0.01, 0.1, 1.0, 10.0)
            f = atomic_free_energy(beta, 0.0, 0.0)
            @test isapprox(f, -log(4) / beta; atol=1e-12)
        end
    end

    @testset "Half-filling specialization matches general form at mu = U/2" begin
        for U in (0.0, 1.0, 4.0, 10.0), beta in (0.1, 1.0, 5.0)
            f_general = atomic_free_energy(beta, U, U/2)
            f_half = atomic_free_energy_half_filling(beta, U)
            @test isapprox(f_general, f_half; atol=1e-12)
        end
    end

    @testset "Half-filling closed form: f = -T ln[2(1 + e^{beta U / 2})]" begin
        for U in (1.0, 4.0, 10.0), beta in (0.1, 1.0, 5.0)
            f = atomic_free_energy_half_filling(beta, U)
            expected = -log(2 * (1 + exp(beta * U / 2))) / beta
            @test isapprox(f, expected; atol=1e-12)
        end
    end

    @testset "Low-T atomic-limit half-filling: f -> -U/2 - T ln 2" begin
        # At large beta the singly-occupied doublet (energy -U/2) dominates with
        # 2-fold spin degeneracy: f = -U/2 - T ln 2 + O(e^{-beta U / 2}).
        U = 4.0
        beta = 20.0
        f = atomic_free_energy_half_filling(beta, U)
        f_expected = -U/2 - log(2) / beta
        # Corrections are O(e^{-beta U / 2}) = O(e^{-40}) at this regime,
        # vastly below 1e-10 absolute.
        @test isapprox(f, f_expected; atol=1e-10)
    end

    @testset "High-T limit f -> -T ln 4 at any U" begin
        # As beta -> 0 the four states are equally weighted, so f -> -T ln 4
        # regardless of U. Use beta * U = 0.01 (small enough that the
        # difference between exp(beta * U / 2) and 1 is O(beta U)).
        for U in (1.0, 4.0, 10.0)
            beta = 1e-6 / max(U, 1e-9)
            f = atomic_free_energy(beta, U, U/2)
            @test isapprox(f, -log(4) / beta; rtol=1e-3)
        end
    end

    @testset "Magnetic field decreases the free energy" begin
        # At fixed mu, an h > 0 favours the up-spin state, lowering f.
        for h in (0.1, 0.5)
            f_zero = atomic_free_energy(1.0, 4.0, 2.0; h=0.0)
            f_h = atomic_free_energy(1.0, 4.0, 2.0; h=h)
            @test f_h < f_zero
        end
    end

    @testset "beta <= 0 raises DomainError" begin
        @test_throws DomainError atomic_free_energy(0.0, 1.0, 0.5)
        @test_throws DomainError atomic_free_energy(-1.0, 1.0, 0.5)
        @test_throws DomainError atomic_free_energy_half_filling(0.0, 1.0)
    end

    @testset "Stage-B stubs error loudly" begin
        @test_throws ErrorException jks_kernel_K_n(0.5 + 0.1im, 1, 2*pi/3)
        @test_throws ErrorException jks_driving_c(0.0; U=1.0, mu=0.5)
        @test_throws ErrorException jks_driving_cbar(0.0; U=1.0, mu=0.5)
    end

    @testset "jks_driving_b returns the constant -h" begin
        @test jks_driving_b(0.0) == 0.0
        @test jks_driving_b(0.0; h=0.5) == -0.5
        @test jks_driving_b(123.0; h=1.0) == -1.0
    end

    @testset "JKSGrid validates inputs" begin
        g = JKSGrid(128, 2*pi/3)
        @test g.N == 128
        @test isapprox(g.gamma, 2*pi/3)
        @test_throws DomainError JKSGrid(0, 2*pi/3)
        @test_throws DomainError JKSGrid(128, 0.0)
        @test_throws DomainError JKSGrid(128, pi)
        @test_throws DomainError JKSGrid(128, -1.0)
    end
end
