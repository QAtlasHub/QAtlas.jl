# test/models/quantum/Hubbard1D/test_hubbard1d_jks_jacobian.jl
#
# The closed-form Jacobian of the JKS NLIE, gated against the finite-difference
# one it replaces.
#
# The finite-difference build is kept, not deleted: it is an INDEPENDENT route
# to the same matrix (differencing the residual, versus writing down the block
# structure), so it is the natural oracle for the analytic one.  A hand-derived
# block form is exactly the kind of thing that looks right and is off by a sign
# in one block, and the solver would then merely converge more slowly somewhere
# else rather than fail outright.

using Test
using LinearAlgebra: norm, I
using Random
using QAtlas
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    init_atomic_limit,
    jks_jacobian_full_analytic,
    jks_jacobian_full_finite_diff,
    hubbard1d_jks_free_energy

@testset "Hubbard1D — JKS analytic Jacobian" begin
    @testset "matches the finite-difference oracle away from the atomic limit" begin
        # The comparison is made at a PERTURBED state, not at the atomic-limit
        # initial point: there the auxiliary functions are constant, so a wrong
        # per-node derivative can still reproduce the matrix by accident.
        rng = MersenneTwister(20260729)
        for N in (8, 16), (U, mu, beta) in ((4.0, 2.0, 0.01), (8.0, 4.0, 0.05))
            grid = JKSContourGrid(N, U / 4; x_max=32.0)
            aux = init_atomic_limit(grid, beta, U, mu)
            aux.b .*= exp.(0.15 .* randn(rng, N) .+ 0.05im .* randn(rng, N))
            aux.b_bar .= aux.b
            aux.c .*= exp.(0.15 .* randn(rng, N) .+ 0.05im .* randn(rng, N))
            aux.c_bar .*= exp.(0.15 .* randn(rng, N) .+ 0.05im .* randn(rng, N))

            Ja = jks_jacobian_full_analytic(aux, grid, beta, U, mu, U / 6)
            Jf = jks_jacobian_full_finite_diff(aux, grid, beta, U, mu, U / 6)
            @test size(Ja) == (3N, 3N)
            @test all(isfinite, Ja)
            # 1e-6 forward differences give about six digits, so this is the
            # ORACLE's accuracy floor, not the analytic form's.
            @test norm(Ja - Jf) / norm(Jf) < 1e-5
        end
    end

    @testset "b_bar enters no residual, so its block column is exactly [0; 0; I]" begin
        # This is the structural claim the block form rests on -- the b channel
        # computes log(1 + 1/b̄) but does not use it in its right-hand side. If a
        # future edit wires b̄ back in, the analytic Jacobian goes silently wrong
        # and this is what catches it.
        N, U, mu, beta = 12, 4.0, 2.0, 0.02
        grid = JKSContourGrid(N, U / 4; x_max=32.0)
        aux = init_atomic_limit(grid, beta, U, mu)
        Ja = jks_jacobian_full_analytic(aux, grid, beta, U, mu, U / 6)
        @test Ja[:, (2N + 1):(3N)] ==
            vcat(zeros(ComplexF64, 2N, N), Matrix{ComplexF64}(I, N, N))
    end

    @testset "the solver returns the same free energy either way" begin
        # β = 1e-3 is inside the regime the row is validated on (see the
        # registry note: high-T agreement to 1%); the point here is that
        # swapping the Jacobian does not move the ANSWER, not that the answer
        # is right at every β.
        U, mu, t, beta = 4.0, 2.0, 1.0, 1e-3
        # grid_N below the 128 default: the point is that the two Jacobians
        # agree, and the finite-difference path costs O(N^3) to make it.
        fa = hubbard1d_jks_free_energy(t, U, mu, beta; grid_N=32, jacobian=:analytic)
        ff = hubbard1d_jks_free_energy(t, U, mu, beta; grid_N=32, jacobian=:finite_diff)
        @test isfinite(fa)
        @test fa ≈ ff rtol = 1e-8
    end

    @testset "an unknown jacobian mode is rejected, not silently defaulted" begin
        @test_throws ArgumentError hubbard1d_jks_free_energy(
            1.0, 4.0, 2.0, 1e-3; grid_N=32, jacobian=:magic
        )
    end
end
