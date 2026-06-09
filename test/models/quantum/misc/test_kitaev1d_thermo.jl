# ─────────────────────────────────────────────────────────────────────────────
# Test: Kitaev1D free-fermion (BdG) finite-temperature thermodynamics.
#
# The Kitaev chain is a free-fermion p-wave BdG system with bulk dispersion
# E(k) = √((2t cos k + μ)² + 4Δ² sin²k).  The newly registered FreeEnergy /
# SpecificHeat / ThermalEntropy(Infinite) and the β-extended Energy are checked
# by genuinely independent routes:
#
#   * energy FDT: c_v == -β² ∂ε/∂β  (ForwardDiff of the registered Energy
#     integral vs the closed-form c_v integral — different routes)
#   * #676 free-fermion foundation: c_v == β²·Var(E)/L from a Riemann sum of the
#     Brillouin zone via fd_free_fermion_thermo (foundation Σ vs QuadGK ∫), with
#     the dispersion re-typed here so a src convention slip would surface
#   * limits: β→∞ recovers ε₀ (the T=0 ground state); high-T s → ln2 (two states
#     per spinless mode); Gibbs s == β(ε−f); c_v ≥ 0, 0 ≤ s ≤ ln2
#
# Gapped parameter points only (|μ| ≠ 2|t| and Δ ≠ 0) so the BZ integrand is
# smooth and the Riemann tie converges fast.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, ForwardDiff
using QAtlas: Kitaev1D, Energy, FreeEnergy, SpecificHeat, ThermalEntropy, Infinite, fetch

@testset "Kitaev1D free-fermion finite-T thermodynamics" begin
    for (μ, t, Δ) in ((0.5, 1.0, 1.0), (2.5, 1.0, 1.0), (0.0, 1.0, 0.8))
        m = Kitaev1D(; μ=μ, t=t, Δ=Δ)
        disp(k) = sqrt((2t * cos(k) + μ)^2 + 4 * Δ^2 * sin(k)^2)   # re-typed (independent)
        ε0 = fetch(m, Energy(:per_site), Infinite())               # T=0 GS (no beta)

        @testset "(μ,t,Δ) = ($μ,$t,$Δ)" begin
            for β in (0.5, 1.0, 2.0)
                cv = fetch(m, SpecificHeat(), Infinite(); beta=β)
                ε = fetch(m, Energy(:per_site), Infinite(); beta=β)
                f = fetch(m, FreeEnergy(), Infinite(); beta=β)
                s = fetch(m, ThermalEntropy(), Infinite(); beta=β)
                @test cv ≥ 0
                @test isapprox(s, β * (ε - f); rtol=1e-10)        # Gibbs
                @test -1e-10 ≤ s ≤ log(2) + 1e-6

                # energy FDT: c_v == -β² ∂ε/∂β  (AutoDiff of the registered Energy)
                dε = ForwardDiff.derivative(
                    b -> fetch(m, Energy(:per_site), Infinite(); beta=b), β
                )
                @test isapprox(cv, -β^2 * dε; rtol=1e-6, atol=1e-8)

                # tie to the #676 free-fermion FDT foundation (BZ Riemann sum)
                ks = range(-π, π; length=801)[1:(end - 1)]
                modes = [disp(k) for k in ks]
                cv_ff = fd_free_fermion_thermo(modes, β).C / length(modes)
                @test isapprox(cv, cv_ff; rtol=2e-3)
            end

            # β → ∞ recovers the T=0 ground state; high-T s → ln2
            @test isapprox(
                fetch(m, Energy(:per_site), Infinite(); beta=50.0), ε0; atol=1e-6
            )
            @test isapprox(fetch(m, FreeEnergy(), Infinite(); beta=50.0), ε0; atol=1e-6)
            @test isapprox(
                fetch(m, ThermalEntropy(), Infinite(); beta=1e-4), log(2); atol=1e-3
            )
        end
    end

    # Cross-model exact check: Kitaev1D(μ=-2h, t=J, Δ=J) has E(k) ≡ TFIM(J,h)'s
    # Λ_k, so its finite-T thermodynamics must equal TFIM's (independently
    # implemented in TFIM_thermal.jl) — the strongest verification here.
    @testset "TFIM point: Kitaev1D(μ=-2h,t=J,Δ=J) thermo == TFIM(J,h)" begin
        for (J, h) in ((1.0, 0.5), (1.0, 1.5), (0.8, 1.0)), β in (0.4, 1.0)
            kit = Kitaev1D(; μ=-2h, t=J, Δ=J)
            tf = QAtlas.TFIM(; J=J, h=h)
            @test isapprox(
                fetch(kit, SpecificHeat(), Infinite(); beta=β),
                fetch(tf, SpecificHeat(), Infinite(); beta=β);
                rtol=1e-6,
            )
            @test isapprox(
                fetch(kit, Energy(:per_site), Infinite(); beta=β),
                fetch(tf, Energy(:per_site), Infinite(); beta=β);
                rtol=1e-6,
            )
            @test isapprox(
                fetch(kit, FreeEnergy(), Infinite(); beta=β),
                fetch(tf, FreeEnergy(), Infinite(); beta=β);
                rtol=1e-6,
            )
        end
    end

    @test_throws DomainError fetch(Kitaev1D(), SpecificHeat(), Infinite(); beta=0.0)
end
