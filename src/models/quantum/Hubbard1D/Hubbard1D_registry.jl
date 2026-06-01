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

# ── Finite-T at Infinite() — Phase-2A stopgap (#523) ────────────────────
# Regime-based delegation: U=0 -> tight-binding, very-high-T -> -T ln 4,
# strong-coupling+low-T -> Lieb-Wu + Heisenberg c=1 CFT, else NaN+warn.
# Full JKS QTM NLIE (4 coupled NLIEs) is deferred to a separate large PR.

@register(
    Hubbard1D,
    FreeEnergy,
    Infinite,
    method=:regime_delegation,
    reliability=:medium,
    tested_in="test/models/quantum/Hubbard1D/test_hubbard1d_thermal_stopgap.jl",
    references=[
        "Lieb-Wu PRL 20, 1445 (1968)",
        "Anderson PR 115, 2 (1959)",
        "Jüttner-Klümper-Suzuki NPB 522, 471 (1998)",
    ],
    notes="Half-filling stopgap. (A) U=0 exact: 2 x TightBinding1D. (B) very-high-T (β·max(t,U)≤0.05): -T ln 4. (C) strong+low-T (U/t≥10, βJ_eff≥5): Lieb-Wu e_0 - T²/(3 J_eff). (D) NaN+warn otherwise.",
)

@register(
    Hubbard1D,
    ThermalEntropy,
    Infinite,
    method=:regime_delegation,
    reliability=:medium,
    tested_in="test/models/quantum/Hubbard1D/test_hubbard1d_thermal_stopgap.jl",
    references=["Lieb-Wu PRL 20, 1445 (1968)", "Anderson PR 115, 2 (1959)"],
    notes="Half-filling stopgap. (A) U=0: 2 x TightBinding1D. (B) high-T: ln 4. (C) strong+low-T: 2T/(3 J_eff). (D) NaN+warn.",
)

@register(
    Hubbard1D,
    SpecificHeat,
    Infinite,
    method=:regime_delegation,
    reliability=:medium,
    tested_in="test/models/quantum/Hubbard1D/test_hubbard1d_thermal_stopgap.jl",
    references=["Lieb-Wu PRL 20, 1445 (1968)", "Anderson PR 115, 2 (1959)"],
    notes="Half-filling stopgap. (A) U=0: 2 x TightBinding1D. (B) high-T: 0 at LO. (C) strong+low-T: 2T/(3 J_eff). (D) NaN+warn.",
)
