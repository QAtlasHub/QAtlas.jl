# test/models/quantum/Hubbard1D/test_hubbard1d_jks_invariants.jl
#
# Zero-by-construction and convergence invariants for the JKS NLIE row (#798).
#
# The contour bug in #798 lived because the only independent check on this row
# was an ED comparison at two high-T points; everything else it could have been
# caught by was either a `@warn` nobody asserted on or was never asked. Each
# check below would have caught it ALONE, and none of them needs an oracle:
#
#   1. `imag(f)` — the free energy of a Hermitian Hamiltonian is real, so this
#      is zero by construction. Its size is the sharpest cheap signal that the
#      equations being solved are not eq (47), and measuring it is what showed
#      the defect is present at EVERY temperature rather than switching on at
#      mid-T (see the table below).
#   2. grid convergence — a well-posed discretisation moves LESS as the grid
#      refines. The current one wanders, which is what integrating a contour
#      over the wrong domain looks like.
#   3. the ED gap is not finite size — comparing N = 4 against N = 6, and OBC
#      against PBC, separates "the NLIE disagrees with the TDL" from "a 4-site
#      open chain is not the TDL".
#
# The failing ones are `@test_broken`, so fixing #798 turns them into
# "Unexpectedly Pass" rather than leaving a green suite over a known-wrong row.

using Test
using LinearAlgebra
using SparseArrays
using QAtlas
using QAtlas: Hubbard1D, FreeEnergy, Infinite
using QAtlas.Hubbard1DJKSNLIE:
    JKSContourGrid,
    jks_kernel_K_n_concrete,
    solve_jks_nlie_full_newton_continuation,
    free_energy_jks,
    free_energy_jks_complex,
    hubbard1d_jks_free_energy

include(joinpath(@__DIR__, "..", "..", "..", "util", "hubbard_ed.jl"))

const _U, _MU, _T = 4.0, 2.0, 1.0

"Converged auxiliary functions at `beta` on a `grid_N`-point grid, or `nothing`."
function _solved(beta; grid_N=128)
    grid = JKSContourGrid(grid_N, _U / 4; x_max=32.0)
    sol = solve_jks_nlie_full_newton_continuation(
        grid, beta, _U, _MU; alpha=_U / 6, tol=1e-6, maxiter=40, outer_maxsteps=200
    )
    return sol.converged ? (sol, grid) : nothing
end

@testset "Hubbard1D JKS NLIE — invariants (#798)" begin
    @testset "the kernels are the ones eq (38) defines" begin
        # eq (38) defines three kernels from k(s) = 1/(2 pi i s):
        #     K1  =  k(s) - k(s + 2ig)   = (g/pi) / [s (s + 2ig)]
        #     K1b = -k(s) + k(s - 2ig)   = (g/pi) / [s (s - 2ig)]
        #     K2  =  k(s - 2ig) - k(s + 2ig) = (2g/pi) / (s^2 + 4g^2)
        # The module carries ONE closed form, `K_n(s) = g/(pi s (s + 2nig))`,
        # and indexes it by n. That is right for n = 1 and is a different
        # function for n = 2 — K2 is real and even on the real axis, this is
        # neither. Rebuilding from k here shares no algebra with the closed
        # form, so it is an independent check rather than a restatement.
        k(s) = 1 / (2im * pi * s)
        pts = (0.3 + 0.1im, 1.0 + 0.0im, 2.5 - 0.4im, -0.7 + 0.9im)
        for gamma in (0.5, 1.0), s in pts
            @test jks_kernel_K_n_concrete(s, 1, gamma) ≈ k(s) - k(s + 2im * gamma)
            @test_broken jks_kernel_K_n_concrete(s, 2, gamma) ≈
                k(s - 2im * gamma) - k(s + 2im * gamma)
        end
        # K2 is real on the real axis; whatever is being used for it is not.
        for gamma in (0.5, 1.0), x in (0.4, 1.3, 3.0)
            @test_broken abs(imag(jks_kernel_K_n_concrete(x, 2, gamma))) < 1e-12
        end
    end

    @testset "imag(f) vanishes — zero by construction for a Hermitian H" begin
        # It does not vanish at ANY temperature. MEASURED (U = 4, half filling,
        # grid_N = 128):
        #
        #   beta    Re f        Im f      |Im/Re|
        #   1e-5    -138645     0.2547    1.8e-6
        #   1e-4     -13862     0.2549    1.8e-5
        #   1e-3      -1383.5   0.2574    1.9e-4
        #   1e-2       -135.7   0.2808    2.1e-3
        #   0.1         -11.08  0.4541    4.1e-2
        #   0.3          -2.030 0.6854    3.4e-1
        #
        # The ABSOLUTE error barely moves across four decades of beta; the ratio
        # moves by five, because |Re f| ~ 1/beta diverges. So the row's
        # "exact at high T to within 1%" is that division and not agreement —
        # the defect is present throughout, and high T only hides it.
        #
        # Pinned at all four points, including inside the registered valid
        # domain: a fix has to make this vanish everywhere, not shrink a ratio.
        for beta in (1e-4, 1e-3, 0.1, 0.3)
            got = _solved(beta)
            @test got !== nothing
            sol, grid = got
            f = free_energy_jks_complex(sol.aux, grid, beta, _U; mu=_MU)
            @test_broken abs(imag(f)) < 1e-8 * max(abs(real(f)), 1.0)
        end
    end

    @testset "refining the grid moves the answer less, not differently" begin
        # A convergent quadrature has |f(2N) - f(N)| shrinking. MEASURED at
        # beta = 0.1: -11.469, -11.089, -11.080, -11.087, -11.092 for grid_N
        # 32, 64, 128, 192, 256 — it turns around, which is what refining a
        # grid on [-x_max, x_max] toward a contour around [-1, 1] does.
        beta = 0.1
        fs = Float64[]
        for gn in (32, 64, 128, 256)
            got = _solved(beta; grid_N=gn)
            @test got !== nothing
            sol, grid = got
            push!(fs, free_energy_jks(sol.aux, grid, beta, _U; mu=_MU))
        end
        steps = abs.(diff(fs))
        @test_broken all(steps[i + 1] < steps[i] for i in 1:(length(steps) - 1))
    end

    @testset "the ED gap is not a finite-size effect" begin
        # If the mid-T disagreement were finite size, ED would MOVE between
        # N = 4 and N = 6 (and between OBC and PBC) by something comparable to
        # the disagreement. It does not: at beta = 0.1 the chain is far hotter
        # than any finite-size scale, and the mid-T comparison in
        # test_hubbard1d_jks_ed_comparison.jl was asserting this rather than
        # showing it.
        #
        # MEASURED at beta = 0.1, U = 4, half filling:
        #   f  = -14.950077 (N=4 OBC), -14.962448 (N=4 PBC), -14.954208 (N=6 OBC)
        #   finite-size spread 0.083%,  JKS gap 26%  ->  the gap is 314x larger.
        # The bounds below are loose against those numbers on purpose: they
        # should fail on a change of KIND, not on the third digit.
        beta = 0.1
        f4_obc = _ed_hubbard_free_energy(4, _T, _U, _MU, beta; pbc=false)
        f4_pbc = _ed_hubbard_free_energy(4, _T, _U, _MU, beta; pbc=true)
        f6_obc = _ed_hubbard_free_energy(6, _T, _U, _MU, beta; pbc=false)
        finite_size_spread = maximum(abs.(diff([f4_obc, f4_pbc, f6_obc]))) / abs(f4_pbc)
        f_jks = hubbard1d_jks_free_energy(_T, _U, _MU, beta)
        jks_gap = abs(f_jks - f4_pbc) / abs(f4_pbc)
        @test finite_size_spread < 0.01           # measured 0.00083
        @test jks_gap > 20 * finite_size_spread   # measured ratio 314
    end
end
