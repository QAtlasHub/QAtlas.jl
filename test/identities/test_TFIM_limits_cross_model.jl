using Test
using QAtlas

# Cross-model identity tests for TFIM at the two parameter limits where
# the model collapses to a textbook closed form:
#
#   J = 0 → independent transverse-field paramagnet (free spin)
#       ε(β)/N = -h tanh(βh)
#       m_x(β) = tanh(βh)
#       s(β) = log(2 cosh(βh)) - βh tanh(βh)
#       f(β) = -log(2 cosh(βh)) / β
#       c_v(β) = (βh)² sech²(βh)
#       χ_xx(β, var convention) = β · sech²(βh)
#
#   h = 0 → classical 1D Ising (diagonal in σᶻ basis)
#       PBC per site at finite β:
#         f = -log(λ₊) / β  with  λ₊ = 2 cosh(βJ) (for h=0)
#                                          ⇒  f = -log(2 cosh(βJ)) / β
#         ε/N = -J tanh(βJ)
#         s   = log(2 cosh(βJ)) - βJ tanh(βJ)
#         c_v = (βJ)² sech²(βJ)
#       OBC per site at large N: same per-site values up to 1/N
#       boundary correction (∝ -J / N at T → 0 from one missing bond).
#
# These limits are independent of QAtlas's free-fermion machinery —
# they probe whether the BdG / Pfeuty implementation correctly
# *reproduces* the trivial product / classical limits.

const _CLOSED_TFIM_FREE_SPIN = (
    f      = (β, h) -> -log(2 * cosh(β * h)) / β,
    ε_per  = (β, h) -> -h * tanh(β * h),
    s_per  = (β, h) -> log(2 * cosh(β * h)) - β * h * tanh(β * h),
    c_v    = (β, h) -> (β * h)^2 * sech(β * h)^2,
    m_x    = (β, h) -> tanh(β * h),
    χ_xx   = (β, h) -> β * sech(β * h)^2,
)

@testset "TFIM J = 0 (free-spin paramagnet) — closed-form match" begin
    # All sites independent.  Per-site values must match the textbook
    # paramagnet to machine precision — closed form is just `tanh / cosh`
    # arithmetic, no quadrature.
    h = 1.0
    for β in (0.3, 1.0, 3.0)
        model = TFIM(; J=0.0, h=h)
        N = 8
        bc = OBC(N)

        ε = QAtlas.fetch(model, Energy(:per_site), bc; beta=β)
        f = QAtlas.fetch(model, FreeEnergy(), bc; beta=β)
        s = QAtlas.fetch(model, ThermalEntropy(), bc; beta=β)
        c = QAtlas.fetch(model, SpecificHeat(), bc; beta=β)
        m_x = QAtlas.fetch(model, MagnetizationX(), bc; beta=β)
        χ_xx = QAtlas.fetch(model, SusceptibilityXX(), bc; beta=β)

        @test ε ≈ _CLOSED_TFIM_FREE_SPIN.ε_per(β, h) atol = 1e-12
        @test f ≈ _CLOSED_TFIM_FREE_SPIN.f(β, h) atol = 1e-12
        @test s ≈ _CLOSED_TFIM_FREE_SPIN.s_per(β, h) atol = 1e-12
        @test c ≈ _CLOSED_TFIM_FREE_SPIN.c_v(β, h) atol = 1e-12
        @test m_x ≈ _CLOSED_TFIM_FREE_SPIN.m_x(β, h) atol = 1e-12
        @test χ_xx ≈ _CLOSED_TFIM_FREE_SPIN.χ_xx(β, h) atol = 1e-12
    end
end

@testset "TFIM J = 0 — Infinite per-site values match the same closed form" begin
    h = 0.7
    for β in (0.5, 2.0)
        model = TFIM(; J=0.0, h=h)
        @test QAtlas.fetch(model, Energy(:per_site), Infinite(); beta=β) ≈
            _CLOSED_TFIM_FREE_SPIN.ε_per(β, h) atol = 1e-10
        @test QAtlas.fetch(model, FreeEnergy(), Infinite(); beta=β) ≈
            _CLOSED_TFIM_FREE_SPIN.f(β, h) atol = 1e-10
        @test QAtlas.fetch(model, MagnetizationX(), Infinite(); beta=β) ≈
            _CLOSED_TFIM_FREE_SPIN.m_x(β, h) atol = 1e-10
        @test QAtlas.fetch(model, SpecificHeat(), Infinite(); beta=β) ≈
            _CLOSED_TFIM_FREE_SPIN.c_v(β, h) atol = 1e-10
    end
end

@testset "TFIM h = 0 (classical Ising chain) — PBC finite-N transfer matrix" begin
    # 1D classical Ising at h = 0 has transfer matrix eigenvalues
    # `λ± = 2 cosh(βJ), 2 sinh(βJ)`, so for an N-spin PBC ring
    #
    #   Z_N = λ_+^N + λ_-^N          (= tr T^N)
    #   f_N = -log(Z_N) / (Nβ)
    #
    # In the N → ∞ limit only λ_+ survives and we recover the textbook
    # `f∞ = -log(2 cosh(βJ))/β`; for finite N (which is what the
    # quantum-TFIM PBC implementation diagonalises) the λ_-^N term
    # contributes ~`tanh^N(βJ)`, exponentially small but non-negligible
    # at small N.  We test the **exact** finite-N closed form.
    J = 1.0
    for β in (0.3, 1.0, 3.0)
        model = TFIM(; J=J, h=0.0)
        N = 8
        bc = PBC(N)

        λ_plus = 2 * cosh(β * J)
        λ_minus = 2 * sinh(β * J)
        Z_N = λ_plus^N + λ_minus^N

        f_pred = -log(Z_N) / (N * β)
        # ε = -∂(log Z)/∂β / N
        # ∂(log Z_N)/∂β = N (λ_+^{N-1} ∂λ_+/∂β + λ_-^{N-1} ∂λ_-/∂β) / Z_N
        # ∂λ_+/∂β = 2 J sinh(βJ),  ∂λ_-/∂β = 2 J cosh(βJ)
        dlogZ_dβ =
            N * (λ_plus^(N - 1) * 2J * sinh(β * J) + λ_minus^(N - 1) * 2J * cosh(β * J)) /
            Z_N
        ε_pred = -dlogZ_dβ / N
        s_pred = β * (ε_pred - f_pred)

        f = QAtlas.fetch(model, FreeEnergy(), bc; beta=β)
        ε = QAtlas.fetch(model, Energy(:per_site), bc; beta=β)
        s = QAtlas.fetch(model, ThermalEntropy(), bc; beta=β)

        @test f ≈ f_pred atol = 1e-12
        @test ε ≈ ε_pred atol = 1e-12
        @test s ≈ s_pred atol = 1e-12
    end
end

@testset "TFIM h = 0 — Infinite per-site values match the classical chain" begin
    J = 1.3
    for β in (0.4, 1.5)
        model = TFIM(; J=J, h=0.0)
        @test QAtlas.fetch(model, Energy(:per_site), Infinite(); beta=β) ≈
            -J * tanh(β * J) atol = 1e-10
        @test QAtlas.fetch(model, FreeEnergy(), Infinite(); beta=β) ≈
            -log(2 * cosh(β * J)) / β atol = 1e-10
        @test QAtlas.fetch(model, SpecificHeat(), Infinite(); beta=β) ≈
            (β * J)^2 * sech(β * J)^2 atol = 1e-10
    end
end

# Heisenberg1D ↔ XXZ1D(Δ=1) cross-model equivalence is already verified
# in `test/models/test_Heisenberg1D_thermal.jl` for every observable on
# the delegator surface, so we do not duplicate it here.
