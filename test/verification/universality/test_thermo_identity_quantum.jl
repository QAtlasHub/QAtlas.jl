# ─────────────────────────────────────────────────────────────────────────────
# Verification: thermodynamic identities for quantum models
#
#   ε(β) = f(β) + T · s(β)        (per site)
#   c_v(β) = -β² · ∂ε/∂β            (per site, from automatic differentiation)
#
# Catches a class of sign / normalisation / per-site-vs-total bugs that
# pairwise (Energy ↔ Energy, FreeEnergy ↔ FreeEnergy) cross-checks miss.
#
# Convention reminder (TFIM): `Energy(OBC)` returns *total* ⟨H⟩;
# `FreeEnergy`, `ThermalEntropy`, `SpecificHeat(OBC)` return *per-site*.
# `Energy(Infinite)` is per-site (the only finite quantity in N → ∞).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, ForwardDiff, Test

# ─────────────────────────────────────────────────────────────────────────────
# Migrated from raw @test loops to verify()-first (PR #449 phase 2 zero-legacy).
#
# Two independent thermodynamic identities are pinned per (J, h, β [, N]):
#
#   (1) Gibbs identity     s = β(ε − f)   (route=:sum_rule on ThermalEntropy)
#   (2) c_v derivative     c_v = -β² ∂ε/∂β  (route=:sum_rule on SpecificHeat,
#                                            independent computed via
#                                            ForwardDiff through QAtlas.fetch)
#
# The previous raw @test loops asserted these as bare isapprox checks; here
# each iteration becomes its own verify card so the cross-path corroboration
# is exposed on the WHY-correct evidence plane.
# ─────────────────────────────────────────────────────────────────────────────

@testset "TFIM ε = f + T·s — Infinite (per-site, verify cards)" begin
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 1.5)), β in (0.5, 1.0, 2.0, 4.0)
        model = TFIM(; J=J, h=h)
        ε = QAtlas.fetch(model, Energy(), Infinite(); beta=β)
        f = QAtlas.fetch(model, FreeEnergy(), Infinite(); beta=β)
        verify(
            model,
            ThermalEntropy(),
            Infinite();
            route=:sum_rule,
            fetch_kw=(; beta=β),
            independent=β * (ε - f),
            agree_within=1e-9,
            at=["J=$(J)", "h=$(h)", "β=$(β)"],
            refs=[
                "Gibbs identity s = β(ε − f) across three independent Infinite-bc fetches (Energy / FreeEnergy / ThermalEntropy)",
            ],
        )
    end
end

@testset "TFIM ε = f + T·s — OBC (Energy total → /N, verify cards)" begin
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0)), N in (4, 6, 8), β in (0.5, 1.0, 2.0)
        model = TFIM(; J=J, h=h)
        E_total = QAtlas.fetch(model, Energy(), OBC(N); beta=β)
        f = QAtlas.fetch(model, FreeEnergy(), OBC(N); beta=β)
        ε = E_total / N
        verify(
            model,
            ThermalEntropy(),
            OBC(N);
            route=:sum_rule,
            fetch_kw=(; beta=β),
            independent=β * (ε - f),
            agree_within=1e-10,
            at=["J=$(J)", "h=$(h)", "N=$(N)", "β=$(β)"],
            refs=[
                "Gibbs identity s = β(ε − f) on OBC(N) — Energy is total ⇒ ε = E/N; FreeEnergy / ThermalEntropy are per-site by convention",
            ],
        )
    end
end

@testset "TFIM c_v = -β² ∂ε/∂β — Infinite (AutoDiff cross-check, verify cards)" begin
    # ForwardDiff differentiates through `quadgk` because the integrand
    # is a plain function of `β` (the BdG dispersion has no β dependence).
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0)), β in (0.5, 1.0, 2.0)
        model = TFIM(; J=J, h=h)
        dε_dβ = ForwardDiff.derivative(
            b -> QAtlas.fetch(model, Energy(), Infinite(); beta=b), β
        )
        verify(
            model,
            SpecificHeat(),
            Infinite();
            route=:sum_rule,
            fetch_kw=(; beta=β),
            independent=-β^2 * dε_dβ,
            agree_within=1e-8,
            at=["J=$(J)", "h=$(h)", "β=$(β)"],
            refs=[
                "c_v = -β² ∂ε/∂β — independent AutoDiff through QAtlas.fetch(Energy, Infinite) at the same (J,h,β)",
            ],
        )
    end
end

@testset "TFIM c_v = -β² ∂ε/∂β — OBC (AutoDiff cross-check, verify cards)" begin
    # Energy(OBC) is total → divide by N to compare with per-site c_v.
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0)), N in (4, 6, 8), β in (0.5, 1.0, 2.0)
        model = TFIM(; J=J, h=h)
        dE_dβ = ForwardDiff.derivative(
            b -> QAtlas.fetch(model, Energy(), OBC(N); beta=b), β
        )
        verify(
            model,
            SpecificHeat(),
            OBC(N);
            route=:sum_rule,
            fetch_kw=(; beta=β),
            independent=-β^2 * dE_dβ / N,
            agree_within=1e-9,
            at=["J=$(J)", "h=$(h)", "N=$(N)", "β=$(β)"],
            refs=[
                "c_v = -β² ∂ε/∂β on OBC(N) — Energy(OBC) is total ⇒ divide by N; per-site SpecificHeat by convention",
            ],
        )
    end
end
