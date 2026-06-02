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
    reliability=:high,
    tested_in="test/standalone/test_hubbard1d.jl",
    references=["Lieb-Wu PRL 20, 1445 (1968)", "Essler et al. (2005)"],
    notes="Lieb-Wu integral E₀/N = -4t² ∫₀^∞ J₀(ω) J₁(ω) / [ω (1+exp(ωU/2t))] dω at half filling (μ=U/2).",
)

@register(
    Hubbard1D,
    ChargeGap,
    Infinite,
    method=:bethe_ansatz,
    reliability=:high,
    tested_in="test/standalone/test_hubbard1d.jl",
    references=["Lieb-Wu PRL 20, 1445 (1968)", "Essler et al. (2005)"],
    notes="Lieb-Wu integral Δ_c = (16t²/U) ∫₁^∞ √(ω²-1)/sinh(2πtω/U) dω at half filling.",
)

@register(
    Hubbard1D,
    SpinGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_hubbard1d.jl",
    references=["Lieb-Wu PRL 20, 1445 (1968)"],
    notes="Spinon branch is rigorously gapless at half filling — returns 0.0.",
)

@register(
    Hubbard1D,
    LuttingerParameter,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_hubbard1d.jl",
    references=["Lieb-Wu PRL 20, 1445 (1968)", "Voit Rep. Prog. Phys. 58, 977 (1995)"],
    notes="K=1 at U=0 free-fermion limit; finite-U Lieb-Wu Bethe ansatz K_ρ, K_σ deferred Phase 2.",
)

@register(
    Hubbard1D,
    FreeEnergy,
    Infinite,
    method=:jks_qtm_nlie,
    reliability=:medium,
    tested_in="test/models/quantum/Hubbard1D/test_hubbard1d_jks_paper_precise.jl",
    references=[
        "Jüttner-Klümper-Suzuki Nucl. Phys. B 522, 471 (1998)", "arXiv:cond-mat/9711310"
    ],
    notes=(
        "Paper-precise eq (47) NLIE in 3 channels (b, c, c̄). FE evaluator uses " *
        "Chebyshev-Gauss quadrature on the cut [-1, 1] (handles 1/sqrt(1-x^2) " *
        "singularity exactly) + paper page-14 direct-form log Λ. " *
        "U-independent and exact at high T to within 1%% (β <= 1e-3 across " *
        "U ∈ {2, 4, 8}). Mid-T (β ~ 0.1) shows physical kinetic corrections."
    ),
)
