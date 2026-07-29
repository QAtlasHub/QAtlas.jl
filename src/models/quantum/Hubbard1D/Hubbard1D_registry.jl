# models/quantum/Hubbard1D/Hubbard1D_registry.jl — declarative implementation map.
#
# Hubbard1D Phase 1 implements Lieb–Wu (1968) closed-form integrals at
# half filling only: `GroundStateEnergyDensity`, `ChargeGap`, `SpinGap`
# at `Infinite()`.  Each row below mirrors a `fetch` method in
# `Hubbard1D.jl`.

@register(
    Hubbard1D,
    GroundStateEnergyDensity,
    Infinite,
    method=:bethe_ansatz,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/standalone/test_hubbard1d.jl",
    references=["LiebWu1968", "Essler2005"],
    notes="Lieb-Wu integral E₀/N = -4t ∫₀^∞ J₀(ω) J₁(ω) / [ω (1+exp(ωU/2t))] dω at half filling (μ=U/2).",
)
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
@register(
    Hubbard1D,
    MassGap,
    Infinite,
    method=:bethe_ansatz,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_hubbard1d.jl",
    references=["LiebWu1968", "Essler2005", "Ovchinnikov1970"],
    notes="Lieb-Wu charge gap (default, type=:charge) Δ_c = (16t²/U) ∫₁^∞ √(ω²-1)/sinh(2πtω/U) dω or spin gap (type=:spin) which is 0.0 at half filling.",
)

@register(
    Hubbard1D,
    ChargeGap,
    Infinite,
    method=:bethe_ansatz,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/standalone/test_hubbard1d.jl",
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
    tested_in="test/standalone/test_hubbard1d.jl",
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
    valid_domain="beta <= 1e-3 at H = 0, mu = U/2 (half filling), t = 1. " *
                 "Verified against ED to within 1%% there across U in {2, 4, 8}. " *
                 "Outside it: measured deviation from ED by beta ~ 0.1, and no " *
                 "value at all for beta >= 0.5, where the beta-continuation " *
                 "stalls and `fetch` returns NaN.",
    error_order="NOT a convergent discretisation error. The c-channel " *
                "convolution is taken over the full [-x_max, x_max] grid with " *
                "the alpha-shifted kernel, where eq (47) restricts it to a " *
                "contour surrounding [-1, 1] (the c functions are analytic " *
                "outside it). MEASURED at beta = 0.1: 74%% of K1 * log C comes " *
                "from |x| > 1, only 4 of 128 grid points lie inside [-1, 1], and " *
                "refining grid_N 32 -> 256 moves the answer non-monotonically " *
                "instead of converging. The discarded imaginary part of f runs " *
                "4%% of Re f at beta = 0.1 and 34%% at beta = 0.3.",
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
        "U-independent and exact at high T to within 1%% (β <= 1e-3 across " *
        "U ∈ {2, 4, 8}). Mid-T (β ~ 0.1) deviation from ED is a formula-level " *
        "bug (Stage G.3+ followup)."
    ),
)
