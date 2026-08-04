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
    valid_domain="High-temperature. MEASURED against ED (N = 6, PBC): within " *
                 "1%% for beta <= 0.1 at U = 4, half filling, degrading to ~15%% " *
                 "by beta = 1. Unlike the eq (47) route it replaced, it returns " *
                 "a finite value across that whole range rather than NaN from " *
                 "beta ~ 0.3 upward. H = 0 and mu = U/2 (half filling) only.",
    error_order="Relative to ED (N = 6, PBC, converged in N), MEASURED at " *
                "U = 4, half filling: 0.14%% at beta = 0.05, 0.52%% at 0.1, " *
                "1.9%% at 0.2, 3.7%% at 0.3, 7.8%% at 0.5, 14.6%% at 1.0. The " *
                "error GROWS with beta and is NOT a discretisation artifact: " *
                "refining the solver's own grids (Nw/Nn/x_max from 96/48/32 to " *
                "256/128/64) moves the beta = 1 answer only from -2.669123 to " *
                "-2.668871 while ED gives -3.125427. So it is the formulation " *
                "or its implementation, and #798 stays open.",
    reliability=:medium,
    tested_in=[
        "test/models/quantum/Hubbard1D/test_hubbard1d_jks53_wiring.jl",
        "test/models/quantum/Hubbard1D/test_hubbard1d_jks_eq53.jl",
    ],
    references=["JuttnerKlumperSuzuki1998"],
    notes=(
        "eq (53) NLIE -- the SIX-unknown real-axis form (two boundary values " *
        "each of b, c, cbar), transcribed from the arXiv LaTeX source. " *
        "Supersedes the eq (47) route this row used to call, which is retained " *
        "in Hubbard1D_jks_nlie.jl for comparison but is no longer reachable " *
        "from `fetch`. That route carried FOUR defects at once (#798): it " *
        "convolved the c channel with K1 where the paper says K1bar, used " *
        "log(1+b) where it says log(1+1/b), collapsed the +- boundary index " *
        "into the bar index so eq (51) needed six unknowns and it had four, " *
        "and -- dominating the error -- was fed mu = U/2 as though it were the " *
        "PAPER's half filling. CONVENTION, now executable rather than prose " *
        "(`_jks_paper_mu` / `_jks_plain_offset` in Hubbard1D.jl): the JKS " *
        "Coulomb term is symmetric, U(n_down-1/2)(n_up-1/2), so ITS half " *
        "filling is mu = 0 and f_paper(mu) = f_plain(mu+U/2) + U/4. The old " *
        "route passed mu straight through, i.e. solved a DOPED system and " *
        "compared it against half-filled ED. That mismatch was silent -- " *
        "nothing errored, the numbers were merely wrong by the doping energy. " *
        "IMPROVEMENT, NOT CLOSURE: eq (53) reproduces the parameter-free " *
        "high-T limit f -> -log(4)/beta and returns a finite value for every " *
        "beta where the old route gave NaN, but see `error_order` for the " *
        "low-T drift that remains."
    ),
)
