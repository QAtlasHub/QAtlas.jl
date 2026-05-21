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
            independent=-(1 / β) * log(2 * cosh(β * J)),
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
            independent=-(1/β) * log(2 * cosh(β*J)),
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
