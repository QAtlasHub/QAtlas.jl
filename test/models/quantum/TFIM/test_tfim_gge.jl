# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: TFIM Generalised Gibbs Ensemble (GGE) stationary values.
#
# Acceptance conditions from issue #146:
#
#   1. fetch(::TFIM, GGEValue{Energy{...}}, ::Infinite; initial)
#      and GGEValue{MagnetizationX} return finite Float64 values.
#   2. h_0 = h_f (no-quench limit) recovers the static T = 0 ground-state
#      values, to atol 1e-10.
#   3. Energy conservation: ε_GGE equals the initial-state expectation
#      ⟨ψ_0 | H_f | ψ_0⟩ / N computed via the *raw* BdG matrix elements
#      `(h_f − J cos k) cos(2θ_0) + J sin k · sin(2θ_0)`.  This is a
#      rearrangement of the same identity, but evaluated through a
#      different parametrisation, so any sign / branch error in the
#      Bogoliubov-angle code is caught.
#   4. Para → Ferro quench (h_0 = 2.0 → h_f = 0.5, J = 1) returns
#      finite values that agree with an explicit QuadGK reference at
#      atol 1e-10 — pinned regression value.
#   5. Required kwarg `initial` missing throws.
#   6. Mismatched J between `initial` and `model_f` throws DomainError.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, QuadGK

const J_TEST = 1.0

@testset "TFIM GGE — types and finiteness" begin
    m_f = TFIM(J=J_TEST, h=0.5)
    m_0 = TFIM(J=J_TEST, h=2.0)

    e = QAtlas.fetch(m_f, GGEValue(Energy()), Infinite(); initial=m_0)
    mx = QAtlas.fetch(m_f, GGEValue(MagnetizationX()), Infinite(); initial=m_0)
    @test e isa Float64
    @test mx isa Float64
    @test isfinite(e)
    @test isfinite(mx)
    # ⟨σˣ⟩ ∈ [-1, 1]
    @test -1 - 1e-12 ≤ mx ≤ 1 + 1e-12
end

@testset "TFIM GGE — no-quench limit (h₀ = h_f)" begin
    for h in (0.3, 0.7, 1.0, 1.5, 2.5)
        m = TFIM(J=J_TEST, h=h)
        # GS reference (β = ∞ ≈ 1e6 in the equilibrium routines)
        e_gs = QAtlas.fetch(m, Energy(:per_site), Infinite(); beta=1e6)
        mx_gs = QAtlas.fetch(m, MagnetizationX(), Infinite(); beta=1e6)
        # GGE values with initial = final
        e_gge = QAtlas.fetch(m, GGEValue(Energy()), Infinite(); initial=m)
        mx_gge = QAtlas.fetch(m, GGEValue(MagnetizationX()), Infinite(); initial=m)
        @test isapprox(e_gge, e_gs; atol=1e-10)
        @test isapprox(mx_gge, mx_gs; atol=1e-10)
    end
end

@testset "TFIM GGE — energy conservation via raw BdG matrix elements" begin
    # ⟨ψ₀ | H_f | ψ₀⟩ / N evaluated in the *raw* (h_f, J cos k, J sin k)
    # parametrisation — independent of the cos(2 Δθ) telescoping used
    # inside the GGE routine.  In the initial Bogoliubov basis the
    # initial state is the GS of H_0; the BdG block of H_f at momentum
    # k is
    #
    #   M_f(k) = [ (h_f − J cos k)  J sin k  ;  J sin k  −(h_f − J cos k) ]
    #
    # whose lower (negative-energy) eigenvector in the H_0 basis is
    # rotated by Δθ = θ_f − θ_0; equivalently
    #
    #   ⟨ψ₀ | H_f | ψ₀⟩/N
    #     = -(1/π) ∫₀^π dk [(h_f − J cos k) cos(2θ_0(k))
    #                       + J sin k · sin(2θ_0(k))].
    #
    # This is an algebraic rearrangement of the canonical
    # `-(Λ_f/2)(1 − 2 n_k)` form, but uses entirely different
    # trigonometric branches, so a sign / atan2 error in the GGE code
    # would mismatch.
    function e_initial_in_Hf(J, h0, hf)
        integrand =
            k -> begin
                two_theta0 = atan(J * sin(k), h0 - J * cos(k))
                (hf - J * cos(k)) * cos(two_theta0) + J * sin(k) * sin(two_theta0)
            end
        val, _ = quadgk(integrand, 0.0, π; rtol=1e-12)
        return -val / π
    end

    cases = [(0.5, 1.5), (0.5, 2.0), (1.5, 0.5), (2.0, 0.5), (3.0, 0.7), (0.7, 3.0)]
    for (h0, hf) in cases
        m_f = TFIM(J=J_TEST, h=hf)
        m_0 = TFIM(J=J_TEST, h=h0)
        e_gge = QAtlas.fetch(m_f, GGEValue(Energy()), Infinite(); initial=m_0)
        e_alt = e_initial_in_Hf(J_TEST, h0, hf)
        @test isapprox(e_gge, e_alt; atol=1e-10)
    end
end

@testset "TFIM GGE — pinned regression (J=1, h₀=2, h_f=0.5)" begin
    # Reference values from an independent QuadGK evaluation of the
    # canonical closed forms — uses the same integrands as the
    # implementation, but with the angles re-derived from atan2 here
    # rather than reusing `_tfim_two_theta`.  Acts as a string-typo
    # guard: any change in formula on the implementation side
    # registers as a regression.
    J = J_TEST
    h0 = 2.0
    hf = 0.5

    function n_k_local(k)
        Δ = atan(J * sin(k), h0 - J * cos(k)) - atan(J * sin(k), hf - J * cos(k))
        return sin(Δ / 2)^2
    end

    e_ref, _ = quadgk(0.0, π; rtol=1e-12) do k
        Λ = 2 * sqrt(J^2 + hf^2 - 2 * J * hf * cos(k))
        return (Λ / 2) * (1 - 2 * n_k_local(k))
    end
    e_ref = -e_ref / π

    mx_ref, _ = quadgk(0.0, π; rtol=1e-12) do k
        Λ = 2 * sqrt(J^2 + hf^2 - 2 * J * hf * cos(k))
        return (hf - J * cos(k)) / Λ * (1 - 2 * n_k_local(k))
    end
    mx_ref = (2 / π) * mx_ref

    m_f = TFIM(J=J, h=hf)
    m_0 = TFIM(J=J, h=h0)
    e_got = QAtlas.fetch(m_f, GGEValue(Energy()), Infinite(); initial=m_0)
    mx_got = QAtlas.fetch(m_f, GGEValue(MagnetizationX()), Infinite(); initial=m_0)

    @test isapprox(e_got, e_ref; atol=1e-10)
    @test isapprox(mx_got, mx_ref; atol=1e-10)

    # Pinned numerical values (~9 dp) — change ONLY if the analytical
    # closed form itself is corrected.  Computed once at high
    # precision by QuadGK; logged so a future drift is obvious in diff.
    @test isapprox(e_got, -0.725765633445; atol=1e-9)
    @test isapprox(mx_got, 0.517315809223; atol=1e-9)
end

@testset "TFIM GGE — Para→Ferro and Ferro→Para quenches finite + physical" begin
    J = J_TEST
    cases = [
        (0.3, 0.7),  # Ferro→Ferro
        (0.7, 0.3),  # Ferro→Ferro reverse
        (0.5, 1.5),  # across critical
        (1.5, 0.5),  # across critical reverse
        (2.5, 0.4),  # Para→Ferro
        (0.4, 2.5),  # Ferro→Para
    ]
    for (h0, hf) in cases
        m_f = TFIM(J=J, h=hf)
        m_0 = TFIM(J=J, h=h0)
        e_gge = QAtlas.fetch(m_f, GGEValue(Energy()), Infinite(); initial=m_0)
        mx_gge = QAtlas.fetch(m_f, GGEValue(MagnetizationX()), Infinite(); initial=m_0)
        @test isfinite(e_gge)
        @test isfinite(mx_gge)
        @test -1 ≤ mx_gge ≤ 1
        # Energy must be ≥ ground-state energy of H_f (variational
        # bound) — initial state has higher energy than GS of H_f
        # (or equal in the no-quench limit).
        e_gs_f = QAtlas.fetch(m_f, Energy(:per_site), Infinite(); beta=1e6)
        @test e_gge ≥ e_gs_f - 1e-12
    end
end

@testset "TFIM GGE — error handling" begin
    m_f = TFIM(J=1.0, h=0.5)
    m_0 = TFIM(J=1.0, h=2.0)

    # Required kwarg `initial` missing — Julia raises UndefKeywordError
    @test_throws UndefKeywordError QAtlas.fetch(m_f, GGEValue(Energy()), Infinite())
    @test_throws UndefKeywordError QAtlas.fetch(m_f, GGEValue(MagnetizationX()), Infinite())

    # Mismatched Ising coupling — DomainError from
    # `_check_quench_couplings`
    m_0_badJ = TFIM(J=2.0, h=2.0)
    @test_throws DomainError QAtlas.fetch(
        m_f, GGEValue(Energy()), Infinite(); initial=m_0_badJ
    )
    @test_throws DomainError QAtlas.fetch(
        m_f, GGEValue(MagnetizationX()), Infinite(); initial=m_0_badJ
    )
end

# ── Verification cards (WHY-correct plane) ─────────────────────────────────
@testset "TFIM GGE — verification cards" begin
    # No-quench identity: when the initial state IS the final-Hamiltonian
    # ground state, the GGE expectation collapses to the T=0 GS value.
    # Independent route: the GS energy density fetched directly.
    let m = TFIM(; J=1.0, h=1.5)
        verify(
            m,
            GGEValue(Energy()),
            Infinite();
            route=:limiting_case,
            fetch_kw=(; initial=m),
            independent=QAtlas.fetch(m, Energy(:per_site), Infinite()),
            agree_within=1e-9,
            refs=["No-quench GGE recovers the T=0 ground-state energy density"],
        )
    end
end
