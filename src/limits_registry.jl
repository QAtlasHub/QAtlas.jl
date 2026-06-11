# limits_registry.jl — asymptotic limit edges (@limits_to, #701).
#
# Starter catalog: one edge whose convergence check both endpoints implement
# INDEPENDENTLY (the non-delegating rule).  XXZ1D's ground-state energy
# density is the Yang–Yang Bethe-ansatz integral/sum; Heisenberg1D's is the
# closed Hulthén value J(1/4 − ln 2) — two independent implementations whose
# Δ → 1⁺ agreement is a genuine cross-check (the generated form of the
# hand-written Δ-sequence pattern of #691).
#
# NOT declared here (yet): XXZ Δ → ∞ → IsingChain1D (the classical-Ising
# anisotropy limit needs the energy normalisation map between the quantum and
# classical conventions — a value_map question to settle with the @dual
# machinery), and β → ∞ ground-state limits (a kwargs-limit, not a
# model-parameter limit; needs a sweep-axis driver in the generator).

# Approach from the critical (planar) side: the gapped-regime e₀(Δ > 1)
# Bethe-ansatz branch is not implemented yet (XXZ_xx_infinite.jl returns NaN
# with a warning there), while the |Δ| ≤ 1 Yang–Yang integral is — so the
# isotropic point is reached as Δ → 1⁻, where e₀(Δ) is smooth.
@limits_to(
    :xxz_isotropic_limit,
    XXZ1D,
    Heisenberg1D,
    param = :Δ,
    approach = [0.9, 0.99, 0.999, 0.9999],
    regime = "Δ → 1⁻ (isotropic antiferromagnetic limit from the planar side)",
    rate = "e₀(Δ) is smooth up to the BKT point, so the error is dominated by the linear Δ-dependence ~ |de₀/dΔ|·(1 − Δ)",
    quantities = [(quantity=GroundStateEnergyDensity, bc=Infinite, final_atol=5.0e-5)],
    references = ["YangYang1969", "Hulthen1938", "desCloizeauxPearson1962"],
    notes = "Yang–Yang integral e₀(Δ) (source) against the closed Hulthén value (target): two independent implementations.",
)
