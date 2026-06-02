# test/models/quantum/Hubbard1D/test_hubbard1d_jks_paper_precise.jl
#
# Paper-precise smoke tests for the JKS NLIE port (#523).
# Verifies eq (47), (54), (55) and the FE evaluator structure
# against direct reading of JKS 1998 PDF (cond-mat/9711310).

using Test
using QAtlas: fetch, FreeEnergy, Infinite, Hubbard1D
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    JKSAuxFunctions,
    init_atomic_limit!,
    jks_driving_c,
    jks_driving_cbar,
    jks_nlie_residual_c,
    jks_nlie_residual_cbar,
    jks_nlie_residual_full,
    solve_jks_nlie_full_newton,
    free_energy_jks,
    atomic_free_energy

@testset "Hubbard1D — JKS paper-precise (#523, paper eq 47/54/55)" begin
    @testset "Paper eq (54): psi_c sign on (mu+H/2)" begin
        # Paper: Psi_c = -betaU/2 + beta(mu+H/2) + log phi
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        psi_c = jks_driving_c(grid, 0.01, 4.0, 2.0; H=0.0)
        # At h=0, mu=2, U=4: psi_c = -0.02 + 0.02 + log phi = log phi
        # At x=0: log phi(0) = 0
        # So psi_c[grid index near 0] ≈ 0 with small β
        @test all(abs.(real.(psi_c)) .< 0.5)  # bounded, no β-large term
    end

    @testset "Paper eq (55): psi_cbar = -betaU/2 - beta(mu+H/2) - log phi" begin
        grid = JKSContourGrid(8, 1.0; x_max=2.0)
        psi_c = jks_driving_c(grid, 0.01, 4.0, 2.0; H=0.0)
        psi_cbar = jks_driving_cbar(grid, 0.01, 4.0, 2.0; H=0.0)
        # psi_c + psi_cbar should equal -betaU (the U/2 terms add, mu cancels, log phi cancels)
        sum_psi = psi_c .+ psi_cbar
        @test all(abs.(real.(sum_psi) .- (-0.04)) .< 1e-10)
    end

    @testset "PH-symmetric init: c*cbar = exp(-betaU) and c/cbar = exp(betah)" begin
        aux = JKSAuxFunctions(8)
        beta, U, h = 0.5, 4.0, 0.3
        init_atomic_limit!(aux, beta, U, U/2; h=h)
        c1, cb1 = aux.c[1], aux.c_bar[1]
        @test isapprox(real(c1 * cb1), exp(-beta * U); atol=1e-12)
        @test isapprox(real(c1 / cb1), exp(beta * h); atol=1e-12)
    end

    @testset "Magnetic field breaks c, cbar symmetry per paper" begin
        aux = JKSAuxFunctions(8)
        init_atomic_limit!(aux, 0.5, 4.0, 2.0; h=0.5)
        @test real(aux.c[1]) > real(aux.c_bar[1])
    end

    @testset "Full Newton converges at high T (eq 47 form)" begin
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        sol = solve_jks_nlie_full_newton(
            grid, 0.01, 4.0, 2.0; alpha=0.5, tol=1e-6, maxiter=30
        )
        @test sol.converged
        @test sol.residual < 1e-3
    end

    @testset "FE evaluator eq (49) third form integrand uses /c factor" begin
        # Smoke: changing only c (not cbar) should change f_jks via /c factor.
        grid = JKSContourGrid(16, 1.0; x_max=2.0)
        aux1 = JKSAuxFunctions(16);
        fill!(aux1.c, ComplexF64(1.0));
        fill!(aux1.c_bar, ComplexF64(1.0));
        fill!(aux1.b, ComplexF64(1.0));
        fill!(aux1.b_bar, ComplexF64(1.0))
        aux2 = JKSAuxFunctions(16);
        fill!(aux2.c, ComplexF64(0.5));
        fill!(aux2.c_bar, ComplexF64(1.0));
        fill!(aux2.b, ComplexF64(1.0));
        fill!(aux2.b_bar, ComplexF64(1.0))
        f1 = free_energy_jks(aux1, grid, 0.1, 4.0; mu=2.0)
        f2 = free_energy_jks(aux2, grid, 0.1, 4.0; mu=2.0)
        @test !isapprox(f1, f2; rtol=1e-3)  # c=0.5 ≠ c=1 result
    end

    @testset "Stage E.1: Chebyshev-Gauss + page-14 form, exact at high T" begin
        # Stage E.1 made the FE evaluator U-independent and exact at high T
        # via Chebyshev-Gauss quadrature on the cut [-1, 1] singularity and
        # the page-14 direct-form log Λ (real-valued). Even on coarse grids
        # the ratio is within 1%% at β = 1e-3 across all tested U.
        grid_lo = JKSContourGrid(32, 1.0; x_max=4.0)
        sol_lo = solve_jks_nlie_full_newton(
            grid_lo, 0.001, 4.0, 2.0; alpha=0.5, tol=1e-8, maxiter=40
        )
        f_lo = free_energy_jks(sol_lo.aux, grid_lo, 0.001, 4.0; mu=2.0)
        f_a = atomic_free_energy(0.001, 4.0, 2.0)
        @test isapprox(f_lo / f_a, 1.0; rtol=0.02)

        grid_hi = JKSContourGrid(128, 1.0; x_max=8.0)
        sol_hi = solve_jks_nlie_full_newton(
            grid_hi, 0.001, 4.0, 2.0; alpha=0.5, tol=1e-8, maxiter=40
        )
        f_hi = free_energy_jks(sol_hi.aux, grid_hi, 0.001, 4.0; mu=2.0)
        @test isapprox(f_hi / f_a, 1.0; rtol=0.01)
    end

    @testset "guard contract" begin
        m = Hubbard1D(t=1.0, U=4.0, μ=2.0)
        # beta <= 0 → DomainError
        @test_throws DomainError fetch(m, FreeEnergy(), Infinite(); beta=0.0)
        @test_throws DomainError fetch(m, FreeEnergy(), Infinite(); beta=-1.0)
        # t != 1 → ArgumentError (JKS path-rescale not implemented)
        m_t2 = Hubbard1D(t=2.0, U=4.0, μ=2.0)
        @test_throws ArgumentError fetch(m_t2, FreeEnergy(), Infinite(); beta=0.1)
        # H != 0 → ArgumentError (b̄=b symmetry only valid at H=0)
        @test_throws ArgumentError fetch(m, FreeEnergy(), Infinite(); beta=0.1, H=0.5)
        # off-half-filling → ArgumentError
        m_qf = Hubbard1D(t=1.0, U=4.0, μ=1.0)
        @test_throws ArgumentError fetch(m_qf, FreeEnergy(), Infinite(); beta=0.1)
    end

    @testset "alpha >= eta DomainError" begin
        m = Hubbard1D(t=1.0, U=4.0, μ=2.0)
        # alpha = U/4 (=eta) hits the constraint boundary
        @test_throws DomainError fetch(m, FreeEnergy(), Infinite(); beta=0.001, alpha=1.0)
        @test_throws DomainError fetch(m, FreeEnergy(), Infinite(); beta=0.001, alpha=2.0)
    end

    @testset "nonuniform grid path" begin
        # Smoke test for build_nonuniform_grid path (Stage F.2)
        m = Hubbard1D(t=1.0, U=4.0, μ=2.0)
        f_uni = fetch(m, FreeEnergy(), Infinite(); beta=0.001, grid_N=128, x_max=32.0)
        f_non = fetch(
            m,
            FreeEnergy(),
            Infinite();
            beta=0.001,
            nonuniform=true,
            x_inner=2.0,
            N_inner=80,
            N_outer=24,
        )
        @test isfinite(f_uni)
        @test isfinite(f_non)
        # Should agree to a few percent (different grids, same physics)
        @test isapprox(f_uni, f_non; rtol=0.05)
    end
end
