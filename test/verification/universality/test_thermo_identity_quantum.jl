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

@testset "TFIM ε = f + T·s — Infinite (per-site, all four)" begin
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 1.5)), β in (0.5, 1.0, 2.0, 4.0)
        model = TFIM(; J=J, h=h)
        ε = QAtlas.fetch(model, Energy(), Infinite(); beta=β)
        f = QAtlas.fetch(model, FreeEnergy(), Infinite(); beta=β)
        s = QAtlas.fetch(model, ThermalEntropy(), Infinite(); beta=β)
        @test isapprox(ε, f + s / β; atol=1e-9, rtol=1e-9)
    end
end

@testset "TFIM ε = f + T·s — OBC (Energy is total; divide by N)" begin
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0)), N in (4, 6, 8), β in (0.5, 1.0, 2.0)
        model = TFIM(; J=J, h=h)
        E_total = QAtlas.fetch(model, Energy(), OBC(N); beta=β)
        f = QAtlas.fetch(model, FreeEnergy(), OBC(N); beta=β)  # per site
        s = QAtlas.fetch(model, ThermalEntropy(), OBC(N); beta=β)  # per site
        ε = E_total / N
        @test isapprox(ε, f + s / β; atol=1e-10, rtol=1e-10)
    end
end

@testset "TFIM c_v = -β² ∂ε/∂β — Infinite (AutoDiff cross-check)" begin
    # ForwardDiff differentiates through `quadgk` because the integrand
    # is a plain function of `β` (the BdG dispersion has no β dependence).
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0)), β in (0.5, 1.0, 2.0)
        model = TFIM(; J=J, h=h)
        c_direct = QAtlas.fetch(model, SpecificHeat(), Infinite(); beta=β)
        dε_dβ = ForwardDiff.derivative(
            b -> QAtlas.fetch(model, Energy(), Infinite(); beta=b), β
        )
        c_ad = -β^2 * dε_dβ
        @test isapprox(c_direct, c_ad; atol=1e-8, rtol=1e-8)
    end
end

@testset "TFIM c_v = -β² ∂ε/∂β — OBC (AutoDiff cross-check)" begin
    # Energy(OBC) is total → divide by N to compare with per-site c_v.
    # ForwardDiff propagates β through the BdG sum without re-diagonalising.
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0)), N in (4, 6, 8), β in (0.5, 1.0, 2.0)
        model = TFIM(; J=J, h=h)
        c_direct = QAtlas.fetch(model, SpecificHeat(), OBC(N); beta=β)
        dE_dβ = ForwardDiff.derivative(
            b -> QAtlas.fetch(model, Energy(), OBC(N); beta=b), β
        )
        c_ad = -β^2 * dE_dβ / N
        @test isapprox(c_direct, c_ad; atol=1e-9, rtol=1e-9)
    end
end
