# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: IsingChain1D — Ising 1925 transfer-matrix closed forms.
#
# Verifies:
#   * T_c = 0 identically (no 1-D finite-temperature transition)
#   * Free energy per site at h = 0 matches -β⁻¹ log(2 cosh βJ) exactly
#   * Free energy at h ≠ 0 matches a direct 2×2 transfer-matrix
#     eigenvalue calculation (independent of the closed-form code)
#   * Correlation length at h = 0 matches 1/log(coth βJ) exactly
#   * Cross-check: ξ → ∞ as T → 0 at J > 0 (low-T limit)
#   * β ≤ 0 throws DomainError
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test
using LinearAlgebra: eigvals

@testset "IsingChain1D — CriticalTemperature is 0" begin
    @test QAtlas.fetch(IsingChain1D(; J=1.0), CriticalTemperature(), Infinite()) == 0.0
    @test QAtlas.fetch(IsingChain1D(; J=-2.0, h=0.5), CriticalTemperature(), Infinite()) ==
        0.0
end

@testset "IsingChain1D — FreeEnergy at h = 0 matches Ising 1925 closed form" begin
    for β in (0.1, 0.5, 1.0, 2.0, 5.0), J in (0.5, 1.0, 2.0)
        m = IsingChain1D(; J=J)
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        @test f ≈ -log(2 * cosh(β * J)) / β atol = 1e-13
    end
end

@testset "IsingChain1D — FreeEnergy at h ≠ 0 matches independent eigvals" begin
    # Build the 2×2 transfer matrix explicitly and diagonalise it, then
    # compare λ_+ → f = -β⁻¹ log λ_+ against the closed-form code.
    for (β, J, h) in ((0.5, 1.0, 0.3), (1.0, 0.7, -0.6), (2.0, 1.5, 0.1), (0.4, -0.8, 0.25))
        T = [
            exp(β*(J+h)) exp(-β*J)
            exp(-β*J) exp(β*(J-h))
        ]
        λs = eigvals(T)
        λp = maximum(λs)
        f_ref = -log(λp) / β
        f_got = QAtlas.fetch(IsingChain1D(; J=J, h=h), FreeEnergy(), Infinite(); beta=β)
        @test f_got ≈ f_ref atol = 1e-13
    end
end

@testset "IsingChain1D — CorrelationLength at h = 0 matches Ising 1925" begin
    # ξ(β, 0) = 1 / log(coth βJ).  The transfer-matrix path computes
    # λ_- = exp(βJ) - exp(-βJ), which loses a few floating-point ulps
    # for βJ ≳ 5; the closed-form comparison is therefore checked at
    # relative ~1e-12 (a few ulp of the ~200-valued ξ at βJ = 6).
    for β in (0.5, 1.0, 2.0, 3.0), J in (0.5, 1.0, 2.0)
        m = IsingChain1D(; J=J)
        ξ = QAtlas.fetch(m, CorrelationLength(), Infinite(); beta=β)
        @test ξ ≈ 1 / log(coth(β * J)) rtol = 1e-12
    end
end

@testset "IsingChain1D — CorrelationLength at h ≠ 0 matches independent eigvals" begin
    for (β, J, h) in ((0.5, 1.0, 0.3), (1.0, 0.7, -0.6), (2.0, 1.5, 0.1), (0.4, 0.8, 0.25))
        T = [
            exp(β*(J+h)) exp(-β*J)
            exp(-β*J) exp(β*(J-h))
        ]
        λs = sort(eigvals(T))
        λm, λp = λs[1], λs[2]
        ξ_ref = 1 / log(λp / λm)
        ξ_got = QAtlas.fetch(
            IsingChain1D(; J=J, h=h), CorrelationLength(), Infinite(); beta=β
        )
        @test ξ_got ≈ ξ_ref atol = 1e-12
    end
end

@testset "IsingChain1D — Low-T limit: ξ diverges at J > 0" begin
    m = IsingChain1D(; J=1.0)
    ξ_1 = QAtlas.fetch(m, CorrelationLength(), Infinite(); beta=5.0)
    ξ_2 = QAtlas.fetch(m, CorrelationLength(), Infinite(); beta=10.0)
    ξ_3 = QAtlas.fetch(m, CorrelationLength(), Infinite(); beta=20.0)
    # ξ(β, 0) = 1/log(coth βJ) ~ exp(2 β J) / 2 as βJ → ∞ — grows fast.
    @test ξ_1 < ξ_2 < ξ_3
    @test ξ_3 > 1e8
end

@testset "IsingChain1D — DomainError on non-positive beta" begin
    m = IsingChain1D(; J=1.0)
    @test_throws DomainError QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=0)
    @test_throws DomainError QAtlas.fetch(m, CorrelationLength(), Infinite(); beta=-1.0)
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "IsingChain1D — verification cards" begin
    # 1D Ising has no finite-T transition (Mermin-Wagner / Ising 1925).
    verify(
        IsingChain1D(; J=1.0),
        CriticalTemperature(),
        Infinite();
        route=:limiting_case,
        independent=0.0,
        agree_within=1e-12,
        refs=["Ising 1925: no finite-T order in 1D, Tc = 0"],
    )

    # h=0 free energy: f = -(1/β) log(2 cosh βJ) (independent closed form)
    for (J, β) in ((1.0, 0.5), (1.0, 2.0), (1.7, 1.0))
        verify(
            IsingChain1D(; J=J),
            FreeEnergy(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=(-(1 / β) * log(2 * cosh(β * J))),
            agree_within=1e-10,
            refs=["Ising 1925: f = -(1/β) log(2 cosh βJ) at h = 0"],
        )
    end

    # h=0 correlation length: ξ = 1 / log(coth βJ)
    for (J, β) in ((1.0, 1.0), (1.0, 2.0))
        verify(
            IsingChain1D(; J=J),
            CorrelationLength(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=1 / log(coth(β * J)),
            agree_within=1e-9,
            refs=["Ising 1925: ξ = 1 / log(coth βJ) at h = 0"],
        )
    end
end
# ── additional verification cards (#381 batch 2) ─────────────────────────
@testset "IsingChain1D — additional closed-form cards (#381 batch 2)" begin
    # CriticalTemperature: 1D Ising has no LRO at any T > 0 ⇒ T_c = 0 (Ising 1925).
    for J in (0.5, 1.0, 2.0)
        verify(
            IsingChain1D(; J=J),
            CriticalTemperature(),
            Infinite();
            route=:second_closed_form,
            independent=0.0,
            agree_within=0,
            refs=[
                "Ising 1925: 1D Ising chain has no spontaneous magnetisation at any T > 0 ⇒ T_c = 0",
            ],
        )
    end
    # FreeEnergy: f(β) = -(1/β) log(2 cosh(βJ)) (Ising 1925; standard result).
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.0))
        verify(
            IsingChain1D(; J=J),
            FreeEnergy(),
            Infinite();
            route=:second_closed_form,
            independent=(-(1/β) * log(2 * cosh(β*J))),
            agree_within=1e-12,
            refs=["Ising 1925: f(β) = -(1/β) log(2 cosh(βJ)) per site"],
            fetch_kw=(; beta=β),
        )
    end
    # CorrelationLength: ξ(β) = -1/log(tanh(βJ)) (Ising/Kramers-Wannier closed form).
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0))
        verify(
            IsingChain1D(; J=J),
            CorrelationLength(),
            Infinite();
            route=:second_closed_form,
            independent=-1/log(tanh(β*J)),
            agree_within=1e-12,
            refs=["Kramers-Wannier 1941: ξ(β) = -1/log(tanh(βJ)) per lattice spacing"],
            fetch_kw=(; beta=β),
        )
    end
end

# ── Additional verification cards (5 closed-form quantities, h=0) ─────────
@testset "IsingChain1D — Energy{:per_site} closed form" begin
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.0), (2.0, 0.7))
        verify(
            IsingChain1D(; J=J),
            Energy{:per_site}(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=(-J * tanh(β * J)),
            agree_within=1e-12,
            refs=["Ising 1925: u(β,h=0) = -J tanh(βJ) per site"],
        )
    end
end

@testset "IsingChain1D — SpecificHeat closed form" begin
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.0), (2.0, 0.7))
        verify(
            IsingChain1D(; J=J),
            SpecificHeat(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=(β * J * sech(β * J))^2,
            agree_within=1e-12,
            refs=["Ising 1925: c_v(β,h=0) = (βJ)² sech²(βJ) per site"],
        )
    end
end

@testset "IsingChain1D — ThermalEntropy closed form + bounds" begin
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.5, 1.0), (2.0, 0.7))
        verify(
            IsingChain1D(; J=J),
            ThermalEntropy(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=log(2 * cosh(β * J)) - β * J * tanh(β * J),
            agree_within=1e-12,
            refs=["Ising 1925: s(β,h=0) = log(2 cosh βJ) - βJ tanh(βJ) per site"],
        )
    end
    for β in (0.1, 1.0, 10.0)
        s = QAtlas.fetch(IsingChain1D(; J=1.0), ThermalEntropy(), Infinite(); beta=β)
        @test 0 ≤ s ≤ log(2) + 1e-12
    end
end

@testset "IsingChain1D — SusceptibilityZZ closed form (Brush 1967)" begin
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (0.5, 2.0), (2.0, 0.3))
        verify(
            IsingChain1D(; J=J),
            SusceptibilityZZ(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=β * exp(2 * β * J),
            agree_within=1e-12,
            refs=["Brush 1967 RMP 39, 883, Eq. (4.18): χ(β,h=0) = β exp(2βJ) per site"],
        )
    end
end

@testset "IsingChain1D — SpontaneousMagnetization = 0 (Ising 1925)" begin
    for (J, β) in ((1.0, 0.5), (1.0, 2.0), (1.0, 10.0), (-0.5, 1.0))
        verify(
            IsingChain1D(; J=J),
            SpontaneousMagnetization(),
            Infinite();
            route=:second_closed_form,
            fetch_kw=(; beta=β),
            independent=0.0,
            agree_within=0,
            refs=[
                "Ising 1925: no spontaneous magnetization at any T > 0 in 1D, m_spont = 0 identically",
            ],
        )
    end
end

@testset "IsingChain1D — h ≠ 0 DomainError on the still-zero-field quantities" begin
    m = IsingChain1D(; J=1.0, h=0.3)
    @test_throws DomainError QAtlas.fetch(m, Energy{:per_site}(), Infinite(); beta=1.0)
    @test_throws DomainError QAtlas.fetch(m, SpecificHeat(), Infinite(); beta=1.0)
    @test_throws DomainError QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=1.0)
    # `SusceptibilityZZ` used to be on this list.  It now has an h ≠ 0 closed
    # form (see the finite-field testset above), so it must NOT throw — and the
    # guard stays meaningful for the three that still lack one.
    @test QAtlas.fetch(m, SusceptibilityZZ(), Infinite(); beta=1.0) > 0
    @test QAtlas.fetch(m, SpontaneousMagnetization(), Infinite(); beta=1.0) == 0.0
end

@testset "IsingChain1D — Gibbs identity s = β(u − f) cross-check at h=0" begin
    for (J, β) in ((1.0, 0.5), (1.0, 1.0), (2.0, 0.7))
        m = IsingChain1D(; J=J)
        u = QAtlas.fetch(m, Energy{:per_site}(), Infinite(); beta=β)
        f = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β)
        s = QAtlas.fetch(m, ThermalEntropy(), Infinite(); beta=β)
        @test s ≈ β * (u - f) atol = 1e-12
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Magnetisation in a longitudinal field.
#
# The closed form m = sinh(βh)/√(sinh²(βh) + e^{-4βJ}) is `-∂f/∂h` carried out
# analytically, so the numerical field derivative is the natural check of the
# ALGEBRA (it is deliberately not an independent physical route — see the
# registry note and #762).  The limits are independent, though: m(h=0) = 0 is
# Ising's no-transition result, and m → tanh(βh) as J → 0 is the free spin.
# ─────────────────────────────────────────────────────────────────────────────
@testset "IsingChain1D magnetisation in a field" begin
    for (J, h, β) in ((1.0, 0.3, 1.0), (1.0, 0.8, 2.0), (0.5, 0.2, 3.0), (2.0, 1.5, 0.7))
        m = IsingChain1D(; J=J, h=h)
        mz = QAtlas.fetch(m, Magnetization{:z}(), Infinite(); beta=β)
        # -∂f/∂h by a central difference, independent of the closed-form code
        δ = 1e-6
        fp = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β, J=J, h=h + δ)
        fm = QAtlas.fetch(m, FreeEnergy(), Infinite(); beta=β, J=J, h=h - δ)
        @test mz ≈ -(fp - fm) / (2δ) atol = 1e-7
        @test abs(mz) <= 1
    end

    # No spontaneous magnetisation in 1-D at any finite temperature (Ising 1925).
    @test QAtlas.fetch(
        IsingChain1D(; J=1.0, h=0.0), Magnetization{:z}(), Infinite(); beta=5.0
    ) == 0.0

    # Free-spin limit: with no coupling the chain is N independent spins.
    m_free = QAtlas.fetch(
        IsingChain1D(; J=1e-12, h=0.4), Magnetization{:z}(), Infinite(); beta=1.0
    )
    @test m_free ≈ tanh(0.4) atol = 1e-6

    # Odd in h, and saturating.
    m_pos = QAtlas.fetch(
        IsingChain1D(; J=1.0, h=0.5), Magnetization{:z}(), Infinite(); beta=1.0
    )
    m_neg = QAtlas.fetch(
        IsingChain1D(; J=1.0, h=-0.5), Magnetization{:z}(), Infinite(); beta=1.0
    )
    @test m_pos ≈ -m_neg atol = 1e-14
    @test QAtlas.fetch(
        IsingChain1D(; J=1.0, h=50.0), Magnetization{:z}(), Infinite(); beta=1.0
    ) ≈ 1.0 atol = 1e-12

    @test_throws DomainError QAtlas.fetch(
        IsingChain1D(; J=1.0, h=0.3), Magnetization{:z}(), Infinite(); beta=0.0
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Longitudinal susceptibility at finite field.
#
# χ(β,h) = β cosh(βh) e^{-4βJ} / (sinh²βh + e^{-4βJ})^{3/2} is ∂m/∂h of the
# transfer-matrix magnetisation, so the numerical field derivative is the check
# of the ALGEBRA (derivation-backed, like the M_z row — see #762).  The h = 0
# reduction to β e^{2βJ} is the independent anchor: that form is textbook
# (Brush 1967) and predates this generalisation.
# ─────────────────────────────────────────────────────────────────────────────
@testset "IsingChain1D susceptibility at finite field" begin
    # Reduces EXACTLY to the zero-field textbook form it replaces.
    for (J, β) in ((1.0, 0.5), (1.0, 2.0), (0.4, 3.0), (2.0, 0.3))
        χ = QAtlas.fetch(IsingChain1D(; J=J, h=0.0), SusceptibilityZZ(), Infinite(); beta=β)
        @test χ ≈ β * exp(2 * β * J) rtol = 1e-13
    end

    # At h ≠ 0 it is ∂m/∂h.
    for (J, h, β) in ((1.0, 0.35, 0.5), (1.0, 0.35, 2.0), (0.5, 1.0, 1.0), (2.0, 0.2, 0.7))
        χ = QAtlas.fetch(IsingChain1D(; J=J, h=h), SusceptibilityZZ(), Infinite(); beta=β)
        δ = 1e-6
        mp = QAtlas.fetch(
            IsingChain1D(; J=J, h=h + δ), Magnetization{:z}(), Infinite(); beta=β
        )
        mm = QAtlas.fetch(
            IsingChain1D(; J=J, h=h - δ), Magnetization{:z}(), Infinite(); beta=β
        )
        @test χ ≈ (mp - mm) / (2δ) rtol = 1e-6
    end

    # Positive, even in h, and saturating: a field aligns the spins, leaving
    # less room to respond.
    for h in (0.5, 1.0, 2.0)
        χp = QAtlas.fetch(
            IsingChain1D(; J=1.0, h=h), SusceptibilityZZ(), Infinite(); beta=1.0
        )
        χm = QAtlas.fetch(
            IsingChain1D(; J=1.0, h=(-h)), SusceptibilityZZ(), Infinite(); beta=1.0
        )
        @test χp > 0
        @test χp ≈ χm rtol = 1e-14
    end
    χ0 = QAtlas.fetch(
        IsingChain1D(; J=1.0, h=0.0), SusceptibilityZZ(), Infinite(); beta=1.0
    )
    χ1 = QAtlas.fetch(
        IsingChain1D(; J=1.0, h=1.0), SusceptibilityZZ(), Infinite(); beta=1.0
    )
    χ2 = QAtlas.fetch(
        IsingChain1D(; J=1.0, h=2.0), SusceptibilityZZ(), Infinite(); beta=1.0
    )
    @test χ0 > χ1 > χ2

    @test_throws DomainError QAtlas.fetch(
        IsingChain1D(; J=1.0, h=0.3), SusceptibilityZZ(), Infinite(); beta=0.0
    )
end
