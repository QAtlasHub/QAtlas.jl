# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: free-fermion finite-T observables for XXZ1D at Δ = 0
# (the XX point) on `Infinite()` boundary conditions.
#
# Validates:
#   • T → 0 limit of Energy(β) reproduces XXZ.jl's GS value -J/π
#   • T → ∞ limit of Energy(β) → 0 and ThermalEntropy(β) → log 2 per site
#   • Thermodynamic identity s = β(e − f)
#   • Energy = ∂(βf)/∂β via central differences (cross-check)
#   • SpecificHeat is non-negative and matches β² ∂²(βf)/∂β²
#   • J-scaling: replacing J with αJ scales (e, f) ↦ α·(e_α=1, f_α=1)
#     at scaled β = β/α (free-fermion homogeneity)
#   • General Δ ≠ 0 still returns NaN with a warning (preserves the
#     #108 placeholder contract)
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

const _XX = XXZ1D(; J=1.0, Δ=0.0)

@testset "XXZ1D Δ=0 / Infinite — Energy: T = 0 and T = ∞ limits" begin
    e0 = QAtlas.fetch(_XX, Energy(), Infinite())              # ground state
    @test e0 ≈ -1 / π atol = 1e-14

    # T = 0 limit (β = 100): converges to e0 within 1e-3
    e_lowT = QAtlas.fetch(_XX, Energy(), Infinite(); beta=100.0)
    @test isapprox(e_lowT, -1 / π; atol=1e-3)

    # T = ∞ limit (β = 0.01): energy ≈ 0
    e_highT = QAtlas.fetch(_XX, Energy(), Infinite(); beta=0.01)
    @test abs(e_highT) < 5e-3
end

@testset "XXZ1D Δ=0 / Infinite — high-T entropy → log 2 per site" begin
    s_highT = QAtlas.fetch(_XX, ThermalEntropy(), Infinite(); beta=0.01)
    @test isapprox(s_highT, log(2); atol=2e-4)
    # at intermediate β the entropy is between 0 and log 2
    s_mid = QAtlas.fetch(_XX, ThermalEntropy(), Infinite(); beta=1.0)
    @test 0 < s_mid < log(2)
    # very low T the entropy → 0
    s_lowT = QAtlas.fetch(_XX, ThermalEntropy(), Infinite(); beta=100.0)
    @test 0 ≤ s_lowT < 0.05
end

@testset "XXZ1D Δ=0 / Infinite — Gibbs identity s = β(e − f)" begin
    for β in (0.05, 0.5, 1.0, 5.0, 20.0)
        e = QAtlas.fetch(_XX, Energy(), Infinite(); beta=β)
        f = QAtlas.fetch(_XX, FreeEnergy(), Infinite(); beta=β)
        s = QAtlas.fetch(_XX, ThermalEntropy(), Infinite(); beta=β)
        @test isapprox(s, β * (e - f); atol=1e-9)
    end
end

@testset "XXZ1D Δ=0 / Infinite — Energy ≈ ∂(βf)/∂β (central differences)" begin
    # Finite-T thermodynamic identity: e(β) = ∂(βf)/∂β.  Symmetric step
    # on β; QuadGK rtol 1e-10 ⇒ 5-digit agreement is comfortable.
    for β0 in (0.5, 1.0, 3.0)
        h = 1e-4
        f_plus = QAtlas.fetch(_XX, FreeEnergy(), Infinite(); beta=β0 + h)
        f_minus = QAtlas.fetch(_XX, FreeEnergy(), Infinite(); beta=β0 - h)
        # ∂(βf)/∂β ≈ ((β+h) f(β+h) − (β−h) f(β−h)) / (2h)
        e_finite_diff = ((β0 + h) * f_plus - (β0 - h) * f_minus) / (2h)
        e_exact = QAtlas.fetch(_XX, Energy(), Infinite(); beta=β0)
        @test isapprox(e_finite_diff, e_exact; atol=1e-6)
    end
end

@testset "XXZ1D Δ=0 / Infinite — SpecificHeat sanity & ∂²(βf)/∂β² check" begin
    # Thermodynamic identity: C = -β² ∂²(βf)/∂β².
    #
    # Derivation (denoting f' = df/dβ):
    #   s = β(e − f) = β·∂(βf)/∂β − βf = β² f',
    #   ∂s/∂T = -β² ∂s/∂β = -β²·(2β f' + β² f''),
    #   C  = T ∂s/∂T = -(β/1)·(2β f' + β² f'')·(1/β) ?  Cleaner:
    #     ∂²(βf)/∂β² = 2 f' + β f'', so β²·∂²(βf)/∂β² = 2β²f' + β³f''
    #     and -β²·∂²(βf)/∂β² = C.  Sign verified empirically: at β = 1
    #     we get c = 0.1045 and β² ∂²(βf)/∂β² = -0.1045.
    for β0 in (0.5, 1.0, 3.0)
        c = QAtlas.fetch(_XX, SpecificHeat(), Infinite(); beta=β0)
        @test c ≥ 0
        h = 1e-3
        bf(β) = β * QAtlas.fetch(_XX, FreeEnergy(), Infinite(); beta=β)
        d2 = (bf(β0 + h) - 2 * bf(β0) + bf(β0 - h)) / h^2
        @test isapprox(c, -β0^2 * d2; rtol=1e-2)
    end
    # Low-T (β = 100): for the gapless XX chain C ∝ T (Luttinger liquid),
    # so we only require C > 0 and C ≪ the β = 1 peak (~0.1).
    c_lowT = QAtlas.fetch(_XX, SpecificHeat(), Infinite(); beta=100.0)
    @test 0 < c_lowT < 5e-2
    # High-T (β = 0.01): C ∝ β² → 0 as β → 0.
    c_highT = QAtlas.fetch(_XX, SpecificHeat(), Infinite(); beta=0.01)
    @test 0 ≤ c_highT < 1e-3
end

@testset "XXZ1D Δ=0 / Infinite — J-scaling (free-fermion homogeneity)" begin
    # ε(k) = -J cos k ⇒ e(β; J) = J · e(β·J; 1), f(β; J) = J · f(β·J; 1)
    # entropy and specific heat depend only on the dimensionless β·J.
    β = 1.7
    for J in (0.5, 2.0, 3.7)
        m = XXZ1D(; J=J, Δ=0.0)
        e_scaled = QAtlas.fetch(m, Energy(), Infinite(); beta=β)
        f_scaled = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        e_ref = QAtlas.fetch(_XX, Energy(), Infinite(); beta=β * J)
        f_ref = QAtlas.fetch(_XX, FreeEnergy(), Infinite(); beta=β * J)
        @test isapprox(e_scaled, J * e_ref; rtol=1e-9)
        @test isapprox(f_scaled, J * f_ref; rtol=1e-9)
        s_scaled = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β)
        c_scaled = QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=β)
        s_ref = QAtlas.fetch(_XX, ThermalEntropy(), Infinite(); beta=β * J)
        c_ref = QAtlas.fetch(_XX, SpecificHeat(), Infinite(); beta=β * J)
        @test isapprox(s_scaled, s_ref; rtol=1e-9)
        @test isapprox(c_scaled, c_ref; rtol=1e-9)
    end
end

@testset "XXZ1D Δ ≠ 0 / Infinite — warn+NaN placeholder for unrouted quantities" begin
    # As of issue #521 (PR #522) FreeEnergy at general -1 < Δ < 1 is wired
    # through the Klümper NLIE and no longer returns NaN. ThermalEntropy and
    # SpecificHeat remain NaN-with-warn pending follow-up; the |Δ| ≥ 0.99
    # endpoint of FreeEnergy also still emits NaN+warn.
    for Δ in (0.5, -0.5)
        m = XXZ1D(; J=1.0, Δ=Δ)
        @test isnan(QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=1.0))
        @test isnan(QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=1.0))
    end
    # FreeEnergy still NaN at the |Δ| ≥ 0.99 endpoints.
    @test isnan(QAtlas.fetch(XXZ1D(; Δ=0.999), FreeEnergy(), Infinite(); beta=1.0))
    @test isnan(QAtlas.fetch(XXZ1D(; Δ=-0.999), FreeEnergy(), Infinite(); beta=1.0))
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "XXZ1D Δ=0 Infinite — verification cards" begin
    Sx, Sy, Sz = spin_ops(1 // 2)
    bond_xx = kron(Sx, Sx) + kron(Sy, Sy)  # J=1, Delta=0
    Ns = verify_profile_Ns(; fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12, 14))

    # Ground-state energy density converges to -1/pi
    verify(
        _XX,
        Energy(),
        Infinite();
        route=:ed_finite_size,
        independent=[
            dense_spectrum(chain_hamiltonian(2, N, bond_xx))[1] / (N - 1) for N in Ns
        ],
        at=["N=$N" for N in Ns],
        agree_within=0.05,
        refs=["Yang-Yang 1966 I: e0 = -J/pi for XX (Delta=0) free fermion"],
    )

    # High-T entropy per site approaches log(2) as beta -> 0
    verify(
        _XX,
        ThermalEntropy(),
        Infinite();
        route=:limiting_case,
        fetch_kw=(; beta=0.01),
        independent=log(2),
        agree_within=2e-4,
        refs=["High-T paramagnet: S/site -> log(2) as beta -> 0 (each spin doublet)"],
    )

    # Gibbs identity s = beta*(e - f): sum rule card for FreeEnergy
    let beta = 1.0
        e = QAtlas.fetch(_XX, Energy(), Infinite(); beta=beta)
        s = QAtlas.fetch(_XX, ThermalEntropy(), Infinite(); beta=beta)
        verify(
            _XX,
            FreeEnergy(),
            Infinite();
            route=:sum_rule,
            fetch_kw=(; beta=beta),
            independent=e - s / beta,
            agree_within=1e-9,
            refs=["Gibbs identity: f = e - s/beta (thermodynamic sum rule)"],
        )
    end
end
