# test/models/quantum/Hubbard1D/test_hubbard1d_jks_stage_c10.jl
#
# Stage C.10 regression: FE evaluator sign fix + high-T numerical
# agreement with atomic_free_energy (#523).

using Test
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    hubbard1d_jks_free_energy,
    atomic_free_energy,
    free_energy_jks,
    JKSContourGrid,
    init_atomic_limit

@testset "Hubbard1D — JKS Stage C.10 sign fix + agreement (#523)" begin
    @testset "Sign of f at high T is negative" begin
        # Atomic limit at any finite T is negative; the JKS evaluator
        # with the Stage C.10 sign convention should agree.
        for beta in (1e-4, 1e-3, 1e-2)
            f = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, beta; tol=1e-3)
            @test f < 0
        end
    end

    @testset "High-T agreement with atomic limit to within 10%" begin
        # Restored after Stage C.24 paper-precise fix: the production wrapper
        # now uses solve_jks_nlie_full_newton (3-channel) with corrected
        # jks_log_z_deriv and FE int_1 prefactor. At default grid (N=64,
        # x_max=8) the high-T ratio is within a few percent of atomic.
        for beta in (1e-4, 1e-3, 1e-2)
            f_jks = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, beta; tol=1e-6)
            f_atom = atomic_free_energy(beta, 4.0, 2.0)
            @test isapprox(f_jks, f_atom; rtol=0.1)
        end
    end

    @testset "f scales as -T at very high T" begin
        # f ~ -T ln 4 at beta -> 0; halving beta (doubling T) should
        # double |f|.
        f_1 = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, 1e-3; tol=1e-3)
        f_2 = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, 5e-4; tol=1e-3)
        # f_2 / f_1 should be ~2 (T doubled means -T doubled in magnitude)
        @test isapprox(f_2 / f_1, 2.0; rtol=0.1)
    end

    @testset "Sign of f matches sign of atomic limit (regression guard)" begin
        # If the FE sign ever flips again, this test fails.
        beta = 1e-3
        f_jks = hubbard1d_jks_free_energy(1.0, 4.0, 2.0, beta; tol=1e-3)
        f_atom = atomic_free_energy(beta, 4.0, 2.0)
        @test sign(f_jks) == sign(f_atom)
    end

    @testset "free_energy_jks direct call: sign convention" begin
        # Plug atomic-limit aux into eq (49) directly; sanity-check the sign.
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        beta, U, mu = 1e-3, 4.0, 2.0
        aux = init_atomic_limit(grid, beta, U, mu)
        f_eq49 = free_energy_jks(aux, grid, beta, U; mu=mu)
        @test f_eq49 < 0  # Stage C.10 sign convention: f at finite T is negative
    end
end
