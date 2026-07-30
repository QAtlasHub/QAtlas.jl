# test/models/quantum/Hubbard1D/test_hubbard1d_jks_eq53.jl
#
# The eq (53) six-unknown JKS NLIE (#798).
#
# Every check here is against something that does NOT come from this solver:
# a residue, a closed form, an algebraic identity the paper states in its own
# text, or the other of two independent routes to the same number. The point of
# #798 was that the row had exactly one oracle (an ED comparison at two high-T
# points) and so a formula-level bug survived; the fix is not worth much without
# the oracles that would have caught it.
#
# Deliberately NOT here: any comparison of the eq (53) route against the eq (47)
# route. They disagree, that is the whole point, and a test that pins their
# difference would pass equally well if both were wrong.

using Test
using LinearAlgebra
using SparseArrays
using QAtlas
using QAtlas: Hubbard1D, GroundStateEnergyDensity, Infinite
using QAtlas.Hubbard1DJKSNLIE:
    JKSGrids53,
    JKSState53,
    JKSOperators53,
    JKSWideOp,
    jks53_ndof,
    jks53_pack,
    jks53_unpack,
    jks53_pv_matrix,
    jks53_kernel_line_integral,
    jks53_apply_wide,
    jks53_script_K,
    jks53_script_K_line_integral,
    jks53_log_phi_narrow,
    jks53_residual,
    jks53_driving_b,
    init_jks53,
    solve_jks53_newton,
    solve_jks53_continuation,
    free_energy_jks53,
    jks_kernel_K1,
    jks_kernel_K1bar,
    jks_kernel_K2

include(joinpath(@__DIR__, "..", "..", "..", "util", "hubbard_ed.jl"))

const _U = 4.0
const _MU = 0.0          # half filling IN THE PAPER'S CONVENTION -- see below
const _G = _U / 4        # gamma = U/4, eq (30)
const _A = 2 * _G / 3    # alpha = 2 gamma / 3, strictly inside (0, gamma)

"""
Exact single-site grand partition function in the paper's convention.

The paper's Coulomb term is SYMMETRIC, `U (n_down - 1/2)(n_up - 1/2)`, so the four
single-site values of `E - mu N` are `U/4, -U/4 - mu, -U/4 - mu, U/4 - 2 mu` and

    z = e^{-beta U/4} (1 + e^{2 beta mu}) + 2 e^{beta (U/4 + mu)}

**Half filling is `mu = 0`**, where this collapses to `4 cosh(beta U/4)`. The plain
`U n_up n_down` form is reached by `f_paper(mu) = f_plain(mu + U/2) + U/4`, so the
familiar "half filling at mu = U/2" belongs to that other convention; feeding it
here solves a doped system.
"""
function _z_atomic(beta, U, mu)
    return exp(-beta * U / 4) * (1 + exp(2 * beta * mu)) + 2 * exp(beta * (U / 4 + mu))
end

"""
Closed-form high-T free energy per site, exact to all orders in the atomic part.
The `t = 1` hopping enters `beta f` only at `O(beta^2)`: at `beta = 0` the trace of
the hopping term vanishes, so `log Lambda = log 4 + beta (mu - U/4 + U/2)` with the
first correction quadratic. No ED, no finite size, no arithmetic shared with the
NLIE.
"""
_f_exact(beta, U, mu) = -log(_z_atomic(beta, U, mu)) / beta

"Relative size of the neglected `O(beta^2)` hopping term in `beta f`."
_hopping_scale(beta, U, mu) = 2 * beta^2 / abs(log(_z_atomic(beta, U, mu)))

@testset "Hubbard1D JKS eq (53) six-unknown NLIE (#798)" begin
    @testset "grid guards encode where the kernels are singular" begin
        # gamma = U/4 is NOT bounded above -- the paper takes gamma -> inf for the
        # strong-coupling reduction to the Heisenberg chain -- so gamma > pi must
        # be accepted. It used to throw, silently capping the row at U < 4 pi.
        @test JKSGrids53(16, 8, 5.0, 3.0) isa JKSGrids53      # gamma = 5 > pi
        @test_throws DomainError JKSGrids53(16, 8, -1.0, 0.5)
        # alpha = gamma puts a pole of K2 on the real axis: eq (53) evaluates K2
        # at shift 2 alpha, and K2(x + 2i gamma) has denominator x (x + 4i gamma).
        @test_throws DomainError JKSGrids53(16, 8, 1.0, 1.0)
        @test_throws DomainError JKSGrids53(16, 8, 1.0, 1.5)
        # ...and that pole is real, not a guard artifact:
        @test !isfinite(jks_kernel_K2(0.0 + 2im * _G, _G))
        @test isfinite(jks_kernel_K2(0.0 + 2im * _A, _G))
    end

    @testset "principal-value rule is exact on a constant" begin
        # PV int_{-1}^{1} dy/(x-y) = log|(x+1)/(x-1)| in closed form. That value
        # comes from outside the quadrature, so agreement is a real check on the
        # rule and not a restatement of how its diagonal was built.
        grids = JKSGrids53(32, 40, _G, _A)
        P = jks53_pv_matrix(grids.xn, grids.wn)
        L = [log(abs((x + 1) / (x - 1))) for x in grids.xn]
        @test maximum(abs.(P * ones(grids.Nn) .- L)) < 1e-12
        # The rule is antisymmetric in the node index pattern: swapping x -> -x
        # maps L(x) -> -L(x), so the operator cannot be accidentally symmetric.
        @test !isapprox(P, P'; rtol=1e-6)
    end

    @testset "kernel line integrals are the residues, and the tail correction uses them" begin
        # int K(x + i d) dx by residues (K1, K1bar have a pole at 0 and one at
        # -+2i gamma; K2 has both at +-2i gamma and none at the origin):
        @test jks53_kernel_line_integral(jks_kernel_K2, _G, 0.0) == 1.0
        @test jks53_kernel_line_integral(jks_kernel_K2, _G, 2 * _A) == 1.0
        @test jks53_kernel_line_integral(jks_kernel_K1bar, _G, _A) == 1.0
        @test jks53_kernel_line_integral(jks_kernel_K1bar, _G, -_A) == 0.0
        @test jks53_kernel_line_integral(jks_kernel_K1, _G, _A) == 0.0
        @test jks53_kernel_line_integral(jks_kernel_K1, _G, -_A) == 1.0
        @test_throws DomainError jks53_kernel_line_integral(jks_kernel_K2, _G, 3 * _G)

        # The residue values are checkable without reusing them: integrate the
        # kernel directly on a fine wide grid, and what is missing must be the
        # analytic tail. K1bar ~ gamma/(pi x^2) for large |x|, so truncating at
        # x_max leaves 2 gamma/(pi x_max) -- 3.2e-4 at gamma = 1, x_max = 2000.
        #
        # dx is held FIXED here on purpose. An earlier version grew x_max at a
        # fixed point count, which coarsens dx in step and trades truncation error
        # for resolution error: at x_max = 1024 the spacing is 1.02 against a
        # kernel whose structure sits on scale gamma = 1, so the total error need
        # not fall at all. It did not, and the test failed in CI.
        let dx = 0.01, x_max = 2000.0
            s = sum(jks_kernel_K1bar(x + im * _A, _G) for x in (-x_max):dx:x_max) * dx
            deficit = abs(s - 1.0)
            @test deficit < 2e-3
            # and it IS the tail, not slop: matches 2 gamma/(pi x_max) in size.
            @test isapprox(deficit, 2 * _G / (pi * x_max); rtol=0.3)
        end

        # Tail-corrected application is exact on a constant, for every operator
        # eq (53) reads off the wide grid. Uncorrected it is not even close at the
        # edges: MEASURED max error 0.46 for K2 at Nw = 96, x_max = 32.
        grids = JKSGrids53(96, 24, _G, _A; x_max=32.0)
        ops = JKSOperators53(grids)
        one_w = ones(ComplexF64, grids.Nw)
        for op in (
            ops.bp_from_Bp,
            ops.bp_from_Bm,
            ops.bm_from_Bp,
            ops.bm_from_Bm,
            ops.c_from_Bbarp,
            ops.c_from_Bbarm,
            ops.cbar_from_Bp,
            ops.cbar_from_Bm,
        )
            @test maximum(abs.(jks53_apply_wide(op, one_w) .- op.I_exact)) < 1e-12
            @test maximum(abs.(op.M * one_w .- op.I_exact)) > 1e-3   # correction earns its keep
        end
    end

    @testset "cal K carries the branch eq (57) fixes" begin
        # eq (57) gives two forms and fixes the branch by cal K(x) ~ 1/(2 pi i x)
        # for large |x|. Only `(2 pi i x sqrt(1 - 1/x^2))^-1` satisfies that; the
        # other form has the opposite sign for real x > 1.
        for x in (50.0, 500.0, -500.0)
            @test isapprox(jks53_script_K(x, 0.0), 1 / (2im * pi * x); rtol=2e-3)
        end
        # Wrong-branch guard: the naive -1/(2 pi sqrt(1 - x^2)) disagrees in SIGN
        # off the cut, so a regression to it fails here rather than silently.
        @test real(jks53_script_K(50.0, 0.0)) ≈ 0 atol = 1e-12
        @test imag(jks53_script_K(50.0, 0.0)) < 0
        @test imag(jks53_script_K(-50.0, 0.0)) > 0
        # Line integrals, by closing away from the cut: -+1/2 for shift >-< 0.
        @test jks53_script_K_line_integral(1.0) == -0.5
        @test jks53_script_K_line_integral(-1.0) == 0.5
        @test_throws DomainError jks53_script_K_line_integral(0.0)
    end

    @testset "log phi boundary values on the cut" begin
        # log phi_{+-0}(x) = -+ 2 i beta sqrt(1 - x^2) for |x| <= 1: pure
        # imaginary, opposite on the two sides, and the jump is what drives the
        # c channel.
        beta = 0.3
        for x in (-0.9, -0.2, 0.0, 0.5, 0.99)
            lp = jks53_log_phi_narrow(x, beta, +1)
            lm = jks53_log_phi_narrow(x, beta, -1)
            @test real(lp) == 0 && real(lm) == 0
            @test lp ≈ -lm
            @test imag(lp) ≈ -2 * beta * sqrt(1 - x^2)
        end
        @test_throws DomainError jks53_log_phi_narrow(1.5, beta, +1)
    end

    @testset "the beta -> 0 solution is b = 1, c = cbar = 1/2, and it is exact" begin
        # Derived, not fitted: as beta -> 0 all driving terms vanish, and with the
        # residue values above the c equation collapses to log c = -log Bbar, the
        # cbar equation to log cbar = -log B, and the b equation to
        # log b^+ = -log B^+ + log B^-, which forces b = 1 and hence c = cbar = 1/2.
        grids = JKSGrids53(48, 24, _G, _A; x_max=32.0)
        ops = JKSOperators53(grids)
        st = init_jks53(grids)
        @test all(st.b_plus .== 1) && all(st.b_minus .== 1)
        @test all(st.c_plus .== 0.5) && all(st.cbar_minus .== 0.5)

        # If it is the exact beta = 0 solution then the residual at small beta is
        # O(beta) -- one power, no constant. Two decades of beta pin the slope.
        r1 = norm(jks53_residual(st, grids, ops, 1e-8, _U, _MU))
        r2 = norm(jks53_residual(st, grids, ops, 1e-6, _U, _MU))
        r3 = norm(jks53_residual(st, grids, ops, 1e-4, _U, _MU))
        @test isapprox(r2 / r1, 100.0; rtol=0.05)
        @test isapprox(r3 / r2, 100.0; rtol=0.05)
        @test r1 < 1e-5      # and the constant term really is absent

        # pack / unpack round-trips the six blocks.
        @test jks53_ndof(grids) == 2 * grids.Nw + 4 * grids.Nn
        v = jks53_pack(st)
        @test length(v) == jks53_ndof(grids)
        st2 = jks53_unpack(v, grids)
        @test st2.b_plus ≈ st.b_plus
        @test st2.cbar_minus ≈ st.cbar_minus
    end

    @testset "log Lambda -> ln 4 as beta -> 0, in BOTH forms of eq (56)" begin
        # z_site -> 4 at beta = 0, so log Lambda -> ln 4. The two forms of eq (56)
        # share no quadrature -- the second integrates the b channel over the
        # whole line and the first never touches b -- so both hitting ln 4 is two
        # independent statements, and their agreement is a third.
        grids = JKSGrids53(96, 24, _G, _A; x_max=32.0)
        st = init_jks53(grids)
        beta = 1e-8
        for form in (:cut, :wide)
            f = free_energy_jks53(st, grids, beta, _U; mu=_MU, form=form)
            @test isapprox(-beta * real(f), log(4); rtol=1e-8)
            @test abs(imag(f)) < 1e-8 * max(abs(real(f)), 1.0)
        end
        f_cut = free_energy_jks53(st, grids, beta, _U; mu=_MU, form=:cut)
        f_wide = free_energy_jks53(st, grids, beta, _U; mu=_MU, form=:wide)
        @test isapprox(real(f_cut), real(f_wide); rtol=1e-7)
        @test_throws ArgumentError free_energy_jks53(
            st, grids, beta, _U; mu=_MU, form=:nonsense
        )
    end

    @testset "solved state satisfies the paper's own algebraic identities" begin
        # In eq (53) the +- index enters the c equations only through Psi_c^+- and
        # the +- 1/2 Dlog term, so differencing them gives, with no quadrature at
        # all and no convolution,
        #     Dlog c    =  Dlog phi + Dlog Cbar
        #     Dlog cbar = -Dlog phi + Dlog C
        # which is what the paper states in prose just before eq (58). These test
        # the +- structure of the residual; they say nothing about the b equation,
        # which is why the anchor test below is the one that matters.
        grids = JKSGrids53(48, 24, _G, _A; x_max=32.0)
        beta = 1e-3
        sol = solve_jks53_continuation(grids, beta, _U, _MU; tol=1e-11)
        @test sol.converged
        st = sol.state
        dlog_c = log.(st.c_plus) .- log.(st.c_minus)
        dlog_cbar = log.(st.cbar_plus) .- log.(st.cbar_minus)
        dlog_C = log.(1 .+ st.c_plus) .- log.(1 .+ st.c_minus)
        dlog_Cbar = log.(1 .+ st.cbar_plus) .- log.(1 .+ st.cbar_minus)
        dlog_phi = [
            jks53_log_phi_narrow(x, beta, +1) - jks53_log_phi_narrow(x, beta, -1) for
            x in grids.xn
        ]
        @test maximum(abs.(dlog_c .- (dlog_phi .+ dlog_Cbar))) < 1e-11
        @test maximum(abs.(dlog_cbar .- (.-dlog_phi .+ dlog_C))) < 1e-11
        # Not vacuous: the jump is nonzero, so the identities have content.
        @test maximum(abs.(dlog_phi)) > 1e-6
        @test maximum(abs.(dlog_Cbar)) > 1e-9
    end

    @testset "eq (53) and eq (58) are different equations, not interchangeable" begin
        # eq (58) rewrites the b equation with Psi_b (eq 59) as its driving term
        # and Dlog(C/Cbar) in place of Dlog(cbar/Cbar). Taking eq (59)'s driving
        # term with eq (53)'s convolution double-counts log phi -- MEASURED as a
        # constant +0.57 offset in f at U = 4 that no grid refinement removed.
        # Psi_b is therefore NOT what jks53_residual uses, and must not become so
        # by accident: it is O(beta) and nonzero, whereas eq (53) drives the b
        # channel with -beta H alone, which vanishes at H = 0.
        grids = JKSGrids53(48, 24, _G, _A; x_max=32.0)
        beta = 1e-3
        psi_b = jks53_driving_b(grids, beta, _U, +1)
        @test maximum(abs.(psi_b)) > beta          # genuinely there
        @test maximum(abs.(psi_b)) < 10 * beta     # and O(beta), not O(1)
        # At H = 0 the eq (53) b driving term is identically zero, so the residual
        # of the b block at the beta -> 0 state must come only from the
        # convolutions. Perturbing H moves it, which shows H is wired in at all.
        ops = JKSOperators53(grids)
        st = init_jks53(grids)
        r0 = jks53_residual(st, grids, ops, beta, _U, _MU; H=0.0)
        rH = jks53_residual(st, grids, ops, beta, _U, _MU; H=0.5)
        @test !isapprox(norm(r0), norm(rH); rtol=1e-6)
    end

    @testset "against the closed-form high-T limit -- no ED involved" begin
        # THE oracle for this row, and the one the eq (47) route fails: it misses
        # this closed form by 5.8e-4 at beta = 1e-4, five orders of magnitude above
        # the only correction allowed there, while reporting "agrees with ED to 1%"
        # because |Re f| ~ 1/beta divides the error down.
        #
        # MEASURED for the eq (53) route: the O(beta) coefficient of log Lambda
        # comes out as `mu` to 4.6e-5 (mu = 1), 1.2e-4 (mu = 2), 2.5e-4 (mu = 3) at
        # beta = 1e-4, which is the size of the O(beta^2) term neglected in that
        # coefficient rather than an error -- so comparing against the full atomic
        # z, as here, is the right test.
        grids = JKSGrids53(96, 24, _G, _A; x_max=32.0)
        for beta in (1e-4, 1e-3), mu in (_MU, 1.0)
            sol = solve_jks53_continuation(grids, beta, _U, mu; tol=1e-11)
            @test sol.converged
            @test sol.residual < 1e-10
            f = free_energy_jks53(sol.state, grids, beta, _U; mu=mu)
            fe = _f_exact(beta, _U, mu)
            rel = abs((real(f) - fe) / fe)
            # Bound set two orders above what was measured, so a regression in
            # KIND fails while quadrature noise does not.
            @test rel < 1e-4
            # imag(f) is zero by construction for a Hermitian H.
            @test abs(imag(f)) < 1e-4 * max(abs(real(f)), 1.0)
        end
        # The bound above is not vacuous only if it is well above the neglected
        # hopping term; assert that ordering explicitly.
        @test _hopping_scale(1e-4, _U, _MU) < 1e-4
    end

    @testset "mid and low temperature -- where the eq (47) route fails" begin
        # The high-T anchor above is NOT enough, and saying so is the whole lesson
        # of #798: the eq (47) route also agreed with ED to 1% for beta <= 1e-3,
        # was 26% off at beta = 0.1, and stopped converging at beta = 0.25. "Exact
        # at high T" is precisely the claim that failed to generalise.
        #
        # The oracle has to change with temperature. The atomic closed form loses
        # its grip because the neglected hopping term is O(beta^2) and reaches
        # 1.3e-2 relative by beta = 0.1, so instead:
        #
        #   mid T : ED at N = 4 (PBC) and N = 6 (OBC). At beta = 0.1 their spread
        #           is 0.083%, so they stand in for the TDL there. The ED helper is
        #           in the PLAIN convention, so compare via
        #           f_paper(mu=0) = f_plain(mu=U/2) + U/4.
        #   low T : the Lieb-Wu ground-state energy density, closed form and with
        #           no finite-size question at all:
        #           f_paper(mu=0) -> E0/L - U/4 as beta -> inf.
        #
        # A coarser grid than the high-T testsets on purpose: this walk is about
        # whether the route reaches these temperatures at all, and the bounds below
        # are deliberately loose enough that grid resolution is not what they test.
        grids = JKSGrids53(64, 16, _G, _A; x_max=32.0)

        # beta = 0.1 is where the #798 disagreement lives, so the route has to
        # reach it at minimum, and land far closer than the 26% the old one did.
        sol = solve_jks53_continuation(grids, 0.1, _U, _MU; tol=1e-10)
        @test sol.converged
        if sol.converged
            f = free_energy_jks53(sol.state, grids, 0.1, _U; mu=_MU)
            fw = free_energy_jks53(sol.state, grids, 0.1, _U; mu=_MU, form=:wide)
            ed4 = _ed_hubbard_free_energy(4, 1.0, _U, _U / 2, 0.1; pbc=true) + _U / 4
            ed6 = _ed_hubbard_free_energy(6, 1.0, _U, _U / 2, 0.1; pbc=false) + _U / 4
            @info "eq (53) at beta = 0.1" f = real(f) imag_f = imag(f) f_wide = real(fw) ed_N4_pbc =
                ed4 ed_N6_obc = ed6 rel_vs_ed4 = abs((real(f) - ed4) / ed4)
            # MEASURED in CI: f = -13.884474 against ED -13.962448 (N=4 PBC) and
            # -13.954208 (N=6 OBC), i.e. 0.56% -- where the eq (47) route is 26%.
            # The two ED points bracket the TDL to 0.083%, so 0.56% is still about
            # 7x the finite-size uncertainty: real disagreement remains here, part
            # of it from this deliberately coarse grid. Bound set at 2%.
            @test abs((real(f) - ed4) / ed4) < 0.02
            @test abs((real(f) - ed6) / ed6) < 0.02
            # MEASURED imag(f) = 0.0199, i.e. 0.14% of |f| (eq (47): 0.394, 3.6%).
            @test abs(imag(f)) < 0.005 * max(abs(real(f)), 1.0)
            @test isapprox(real(f), real(fw); rtol=1e-3)
        end

        # MEASURED in CI, from the pass that added this walk:
        #
        #   beta   f            imag(f)   imag/|f|   note
        #   0.1    -13.884474   0.0199    0.14%      ED gap 0.56% (eq 47: 26%)
        #   0.3     -4.688110   0.0627    1.3%
        #   1.0     -1.667090   0.179     10.8%      Lieb-Wu gap -0.0934
        #   4.0    did not converge (residual 0.016)
        #
        # So the route reaches beta = 1.0 where eq (47) stalled at 0.25, and it is
        # NOT yet right at low T: imag(f) is zero by construction for a Hermitian H
        # and instead grows to 10.8% of |f| by beta = 1. That is the remaining
        # defect, and the bounds below are set from these numbers so it cannot
        # quietly get worse.
        e0 = QAtlas.fetch(Hubbard1D(; t=1.0, U=_U), GroundStateEnergyDensity(), Infinite())
        low_T_target = e0 - _U / 4
        @info "Lieb-Wu low-T target for f_paper(mu=0)" E0_over_L = e0 target = low_T_target
        reached = 0.1
        walk = Tuple{Float64,Float64}[]
        for beta in (0.3, 1.0, 4.0)
            s = solve_jks53_continuation(grids, beta, _U, _MU; tol=1e-10)
            if !s.converged
                @info "eq (53) did not reach beta" beta residual = s.residual
                continue
            end
            reached = beta
            f = free_energy_jks53(s.state, grids, beta, _U; mu=_MU)
            fw = free_energy_jks53(s.state, grids, beta, _U; mu=_MU, form=:wide)
            @info "eq (53) reached beta" beta f = real(f) imag_f = imag(f) gap_to_lowT =
                real(f) - low_T_target
            # Oracle-free and temperature-independent: a Hermitian H has real f,
            # and the two forms of eq (56) share no quadrature.
            @test abs(imag(f)) < 0.2 * max(abs(real(f)), 1.0)
            @test isapprox(real(f), real(fw); rtol=1e-2)
            # f is bounded ABOVE by the ground-state energy density, not below:
            # f = e - T s with s > 0, and df/dT = -s < 0, so f(T) <= f(0) = e0 and
            # f -> -T log 4 as T -> inf. An earlier version of this line asserted
            # `real(f) > low_T_target - 1.0` and failed at beta = 0.3 with
            # -4.688 > -2.574 -- the inequality was backwards, and the solver was
            # right.
            @test real(f) <= low_T_target + 1e-6
            push!(walk, (beta, real(f)))
        end
        @info "eq (53) continuation reached" beta_max = reached walk = walk
        # MEASURED: the continuation reaches beta = 1.0. The eq (47) route stalled
        # at 0.25, so this is the floor now, not 0.1.
        @test reached >= 1.0
        # f rises monotonically toward the ground state as T falls -- MEASURED
        # -13.884 (0.1) -> -4.688 (0.3) -> -1.667 (1.0), approaching -1.5737 from
        # below, gap shrinking 3.11 -> 0.093. Monotonicity is a thermodynamic law
        # (df/dbeta = (e - f)/beta > 0 wherever the entropy is positive), so this
        # needs no oracle.
        @test length(walk) >= 2
        for i in 2:length(walk)
            @test walk[i][2] > walk[i - 1][2]
        end
        if !isempty(walk)
            @test walk[end][2] < low_T_target
            @test abs(walk[end][2] - low_T_target) < 0.5
        end
    end

    @testset "where the low-T imag(f) comes from -- discretisation or formulation" begin
        # imag(f) is zero by construction for a Hermitian H, and MEASURED it grows
        # 0.14% -> 1.3% -> 10.8% of |f| across beta = 0.1, 0.3, 1.0 on the coarse
        # walk grid. Something is wrong at low T; this separates WHICH kind of
        # wrong, by the same test that found the b-equation error at high T:
        #
        #   refine  ->  error falls    =>  discretisation, fix the grid
        #   refine  ->  error plateaus =>  formulation, fix the equations
        #
        # Each axis is refined alone so the answer says which one matters. Recorded
        # through `@info` rather than bounded here: the point is the trend, and a
        # bound on the trend would need the trend measured first.
        beta = 0.5
        base = (Nw=64, Nn=16, x_max=32.0)

        function _imag_ratio(Nw, Nn, x_max)
            grids = JKSGrids53(Nw, Nn, _G, _A; x_max=x_max)
            sol = solve_jks53_continuation(grids, beta, _U, _MU; tol=1e-10)
            sol.converged || return (NaN, NaN, false)
            f = free_energy_jks53(sol.state, grids, beta, _U; mu=_MU)
            return (abs(imag(f)) / max(abs(real(f)), 1.0), real(f), true)
        end

        # narrow grid: the c/cbar channel on [-1, 1], where the PV convolutions and
        # the +- 1/2 Dlog boundary terms live. Nn = 16 is coarse for a function with
        # sqrt(1-x^2) structure, so this is the first suspect.
        narrow = NamedTuple[]
        for Nn in (16, 32, 64)
            r, f, ok = _imag_ratio(base.Nw, Nn, base.x_max)
            push!(narrow, (Nn=Nn, imag_ratio=r, f=f, converged=ok))
        end
        @info "low-T imag(f) vs narrow-grid Nn" beta = beta rows = narrow

        # wide grid: the b channel. Its tail is already restored analytically, so a
        # trend here would point at resolution rather than truncation.
        wide = NamedTuple[]
        for Nw in (64, 128)
            r, f, ok = _imag_ratio(Nw, base.Nn, base.x_max)
            push!(wide, (Nw=Nw, imag_ratio=r, f=f, converged=ok))
        end
        @info "low-T imag(f) vs wide-grid Nw" beta = beta rows = wide

        # The one assertion: SOME refinement has to help. If refining every axis
        # leaves imag(f) where it was, the remaining defect is in the equations and
        # not in the mesh -- which is a different bug and needs to fail loudly here
        # rather than be absorbed into a tolerance.
        ratios = [r.imag_ratio for r in vcat(narrow, wide) if r.converged]
        @test length(ratios) >= 3
        @test minimum(ratios) < 0.9 * first(ratios)
    end

    @testset "refining the grid moves the answer less, not differently" begin
        # A convergent discretisation has |f(fine) - f(coarse)| shrinking. This is
        # what separated formulation error from discretisation error while #798 was
        # being tracked down: the error plateaued at 0.565 under refinement, which
        # is what sent the search to the b equation.
        beta = 1e-3
        fs = Float64[]
        for (Nw, Nn) in ((48, 12), (96, 24), (192, 48))
            grids = JKSGrids53(Nw, Nn, _G, _A; x_max=32.0)
            sol = solve_jks53_continuation(grids, beta, _U, _MU; tol=1e-11)
            @test sol.converged
            push!(fs, real(free_energy_jks53(sol.state, grids, beta, _U; mu=_MU)))
        end
        steps = abs.(diff(fs))
        @test steps[2] < steps[1]
    end

    @testset "large U is reachable" begin
        # gamma = U/4 > pi used to throw. U = 16 gives gamma = 4.
        beta = 1e-3
        U = 16.0
        g = U / 4
        grids = JKSGrids53(96, 24, g, 2 * g / 3; x_max=32.0)
        sol = solve_jks53_continuation(grids, beta, U, 0.0; tol=1e-11)
        @test sol.converged
        f = free_energy_jks53(sol.state, grids, beta, U; mu=0.0)
        @test isfinite(real(f))
        @test isapprox(real(f), _f_exact(beta, U, 0.0); rtol=1e-3)
    end
end
