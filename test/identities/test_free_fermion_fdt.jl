# ─────────────────────────────────────────────────────────────────────────────
# Test: the free-fermion energy fluctuation–dissipation theorem, in two pillars.
#
#   PILLAR ① — autodiff self-consistency (model-independent).  For a set of
#   independent fermionic modes with quasiparticle energies {Λₖ} (fd_free_
#   fermion_thermo), the heat capacity from energy FLUCTUATIONS equals the one
#   from the temperature RESPONSE,
#
#       C = β² Var(E) = β² Σₖ Λₖ² fₖ(1-fₖ)   ==   -β² ∂⟨E⟩/∂β,
#
#   and the entropy is a free-energy response  S == -∂F/∂T = β² ∂F/∂β.  Both
#   right-hand sides are taken by ForwardDiff (the Fermi factor fₖ(β) is
#   analytic — no eigensolve), so the agreement is a genuine theorem check on
#   the fluctuation vs derivative routes, for ARBITRARY mode sets.
#
#   PILLAR ② — verification on a real model (TFIM).  The newly registered
#   `fetch(TFIM, SpecificHeat, Infinite)` = `(β²/4π)∫₀^π Λ(k)² sech²(βΛ/2) dk`
#   is checked two independent ways:
#     • CONSISTENCY (autodiff): `c_v == -β² ∂ε/∂β`, the AutoDiff β-derivative of
#       the *already-registered* `Energy{:per_site}` (ForwardDiff straight
#       through the QuadGK integrand) — catches any prefactor/sign drift.
#     • INDEPENDENT PHYSICS (many-body ED): `c_v(Infinite) ≈ β²·Var(H)/N` from a
#       dense diagonalisation of `build_tfim` (no free-fermion assumption).
#       Gapped points (h ≠ J) so the finite-N (N=12, OBC) value is close to the
#       thermodynamic limit.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test, ForwardDiff, Random
using QAtlas: TFIM, SpecificHeat, Energy, Infinite, fetch

# ── PILLAR ① — autodiff self-consistency of the free-fermion thermodynamics ───
@testset "free-fermion FDT (abstract) — C = β²Var(E) = -β²∂⟨E⟩/∂β" begin
    rng = MersenneTwister(0x1234)
    modesets = [
        ("two modes", [0.5, 2.0]),
        ("uniform gaps (8)", collect(0.5:0.5:4.0)),
        ("random positive (20)", sort(3 .* rand(rng, 20) .+ 0.1)),
        ("near-degenerate", [1.0, 1.0, 1.0, 2.5]),
        ("wide-scale (12)", sort(10 .* rand(rng, 12) .+ 0.05)),
    ]
    for (name, modes) in modesets, β in (0.2, 0.7, 1.5, 3.0)
        th = fd_free_fermion_thermo(modes, β)
        @test th.varE ≥ 0 && th.C ≥ 0 && th.S ≥ 0

        # energy FDT: fluctuation route  ==  temperature-response route (autodiff)
        dE = ForwardDiff.derivative(b -> sum(modes ./ (exp.(b .* modes) .+ 1)), β)
        @test isapprox(th.C, -β^2 * dE; rtol=1e-9, atol=1e-12)

        # entropy as a free-energy response: S == -∂F/∂T = β²∂F/∂β  (autodiff)
        F(b) = -sum(log1p.(exp.(-b .* modes))) / b
        dF = ForwardDiff.derivative(F, β)
        @test isapprox(th.S, β^2 * dF; rtol=1e-7, atol=1e-9)
    end
end

# ── PILLAR ② — real model: the registered TFIM SpecificHeat ───────────────────
@testset "free-fermion FDT (real system) — TFIM SpecificHeat(Infinite)" begin
    # (a) CONSISTENCY: c_v == -β²∂ε/∂β — AutoDiff of the registered Energy.
    for (J, h) in ((1.0, 0.5), (1.0, 1.0), (1.0, 2.0), (0.7, 1.3)), β in (0.3, 0.8, 1.5)
        m = TFIM(; J=J, h=h)
        cv = fetch(m, SpecificHeat(), Infinite(); beta=β)
        @test cv ≥ 0
        dε = ForwardDiff.derivative(b -> fetch(m, Energy(:per_site), Infinite(); beta=b), β)
        @test isapprox(cv, -β^2 * dε; rtol=1e-6, atol=1e-8)
    end

    # (b) INDEPENDENT PHYSICS: many-body dense ED at gapped points (fast
    # finite-size convergence) reproduces the thermodynamic-limit c_v.
    for (J, h) in ((1.0, 0.5), (1.0, 2.0)), β in (0.5, 1.0)
        N = 12
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        levels = eigvals(Symmetric(build_tfim(lat, J, h)))      # 2ᴺ many-body energies
        cv_ed = fd_thermo_from_spectrum(levels, β).C / N         # β²·Var(H)/N
        cv_inf = fetch(TFIM(; J=J, h=h), SpecificHeat(), Infinite(); beta=β)
        @test isapprox(cv_ed, cv_inf; rtol=8e-2)                 # OBC finite-size, N=12
    end
end
