# models/quantum/Hubbard1D/Hubbard1D_registry.jl — declarative implementation map.
#
# Hubbard1D Phase 1 implements Lieb–Wu (1968) closed-form integrals at
# half filling only: `Energy{:per_site}`, `ChargeGap`, `SpinGap`
# at `Infinite()`.  Each row below mirrors a `fetch` method in
# `Hubbard1D.jl`.

@register(
    Hubbard1D,
    Energy{:per_site},
    Infinite,
    method=:bethe_ansatz,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_hubbard1d.jl",
    references=["LiebWu1968", "Essler2005"],
    notes="Lieb-Wu integral E₀/N = -4t ∫₀^∞ J₀(ω) J₁(ω) / [ω (1+exp(ωU/2t))] dω at half filling (μ=U/2).",
)
# No `MassGap` row here, deliberately (#807).  "The gap" is genuinely ambiguous for
# this model — the charge sector is gapped by U and the spin sector is rigorously
# gapless — and both unambiguous names are registered below.  The row that used to sit
# here selected between them with a `type::Symbol=:charge` keyword whose two branches
# were the SAME EXPRESSIONS as the `ChargeGap` and `SpinGap` fetches, so it encoded the
# sector axis a second time and made a bare `fetch(m, MassGap(), Infinite())` return the
# Mott gap for a model whose spectral gap (AbstractQAtlas: `Δ = E₁ − E₀`) is 0.
# Asking for `:gap` here now raises, which is the honest answer: name the sector.

@register(
    Hubbard1D,
    ChargeGap,
    Infinite,
    method=:bethe_ansatz,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_hubbard1d.jl",
    references=["LiebWu1968", "Essler2005"],
    notes="Lieb-Wu integral Δ_c = (16t²/U) ∫₁^∞ √(ω²-1)/sinh(2πtω/U) dω at half filling.",
)

@register(
    Hubbard1D,
    SpinGap,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_hubbard1d.jl",
    references=["LiebWu1968"],
    notes="Spinon branch is rigorously gapless at half filling — returns 0.0.",
)

@register(
    Hubbard1D,
    LuttingerParameter,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_hubbard1d.jl",
    references=["LiebWu1968", "Voit1995"],
    notes="K=1 at U=0 free-fermion limit; finite-U Lieb-Wu Bethe ansatz K_ρ, K_σ deferred Phase 2.",
)

@register(
    Hubbard1D,
    FreeEnergy,
    Infinite,
    method=:jks_qtm_nlie,
    cost=:polynomial,
    # The QTM/NLIE FORMULATION is exact — a finite closed set of equations for
    # three auxiliary functions, not a truncation of the infinite TBA hierarchy.
    # This IMPLEMENTATION is not, and the limitation was carried only in the
    # prose of `notes` where no query could see it (the same over-claim #792
    # fixed elsewhere).
    status = :approx,
    valid_domain="beta <= 1e-3 at H = 0, mu = U/2 (half filling), t = 1, where " *
                 "the ED comparison agrees to within 1%% across U in {2, 4, 8}. " *
                 "That agreement is NOT evidence the equations are right there: " *
                 "MEASURED, the spurious imaginary part of f is ~0.28 at every " *
                 "beta from 1e-5 to 0.1, and only looks small because |Re f| ~ " *
                 "1/beta diverges. Outside: deviation from ED by beta ~ 0.1, and " *
                 "no value at all for beta >= 0.3, where the beta-continuation " *
                 "stalls and `fetch` returns NaN.",
    error_order="NOT a convergent discretisation error. The real-axis unknowns " *
                "are the FOUR (b, b_bar, c, c_bar) where eq (51)-(53) require " *
                "SIX -- two boundary values each of b, c, c_bar -- so the +- " *
                "index has been collapsed into the bar index. Delta log C := " *
                "log(C^+/C^-) and the +- 1/2 Delta log terms of eq (53) have no " *
                "representation in that state, and the c-channel convolution is " *
                "consequently taken over the full [-x_max, x_max] grid where " *
                "eq (47) restricts it to a contour around [-1, 1]. MEASURED at " *
                "beta = 0.1: 74%% of the c-channel convolution comes from " *
                "|x| > 1, only 4 of 128 grid points lie inside [-1, 1], and " *
                "refining grid_N 32 -> 256 moves the answer non-monotonically " *
                "instead of converging.",
    reliability=:medium,
    tested_in="test/models/quantum/Hubbard1D/test_hubbard1d_jks_paper_precise.jl",
    references=["JuttnerKlumperSuzuki1998"],
    notes=(
        "Paper-precise eq (47) NLIE in 3 channels (b, c, c̄). FE evaluator uses " *
        "Chebyshev-Gauss quadrature on the cut [-1, 1] (handles 1/sqrt(1-x^2) " *
        "singularity exactly) + paper page-14 direct-form log Λ. " *
        "Currently SUPPORTS H=0 AND μ = U/2 (half-filling) ONLY: the b/b̄ " *
        "particle-hole symmetry is enforced in the solver via b̄ = b, which " *
        "is exact at H=0 half-filling and breaks for H ≠ 0 or off-half-filling. " *
        "Kernels are the three eq (38) functions K1, K1bar, K2. Mid-T " *
        "(β ~ 0.1) deviation from ED is a formula-level bug, tracked in #798: " *
        "the residuals still convolve the c channel with K1 where eq (47) says " *
        "K1bar, use log(1+b) where it says log(1+1/b), and carry one boundary " *
        "value per function where eq (53) needs two. CONVENTION: the JKS paper's " *
        "Coulomb term is symmetric, U(n_down-1/2)(n_up-1/2), so ITS half filling " *
        "is mu = 0 and f_paper(mu) = f_plain(mu+U/2) + U/4 relative to the plain " *
        "U n_up n_down form the Lieb-Wu rows above use. This route is fed " *
        "mu = U/2 as though it were the paper's half filling, so it has been " *
        "solving a doped system and comparing it against half-filled ED. A " *
        "corrected eq (53) solver lives in Hubbard1D_jks_eq53.jl and reproduces " *
        "the closed-form high-T limit; wiring this row to it is the #798 followup."
    ),
)
