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
