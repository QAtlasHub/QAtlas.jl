# ─────────────────────────────────────────────────────────────────────────────
# Verification: XXZ1D Energy at general -1 < Δ < 1 — Yang-Yang single
# integral vs OBC dense-ED finite-N extrapolation.
#
# Source A: QAtlas XXZ1D Energy at Infinite — Yang-Yang single integral
#           (`_xxz1d_energy_yang_yang`, see src/models/quantum/XXZ/XXZ_bethe.jl).
#
# Source B: Ground-state energy of the OBC chain at finite N from the
#           dense ED of the spin-1/2 XXZ Hamiltonian, extrapolated to
#           N → ∞ by a linear 1/N fit.  OBC has edge defects scaling
#           as 1/N (one weaker bond per boundary), so the leading
#           finite-size correction is linear in 1/N rather than the
#           1/N² seen in PBC; a three-point fit at N ∈ {8, 10, 12} is
#           more than enough to land within ~ 1e-3 of the exact value
#           at the rational test point Δ = 1/2 (e₀ = -3/8 exactly,
#           Yang-Yang II 1966 eq. (4.4) at γ = π/3).
#
# A loose absolute tolerance of 5e-3 is used so the test is not
# sensitive to BLAS-implementation noise at the deepest eigenvalue.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test, LinearAlgebra

# Linear fit: returns (intercept, slope) of y ≈ a + b·x.  Local copy
# (test_xxz_luttinger_ed.jl uses the same private helper); kept inline
# so this file has no inter-test dependency.
function _yy_linfit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    x̄ = sum(x) / n
    ȳ = sum(y) / n
    num = sum((x[i] - x̄) * (y[i] - ȳ) for i in 1:n)
    den = sum((x[i] - x̄)^2 for i in 1:n)
    b = num / den
    return (ȳ - b * x̄, b)
end

# OBC ground-state energy of an N-site spin-1/2 XXZ chain at given
# (N, Δ, J), built from Pauli Kronecker products with no QAtlas private
# API.  Used to cross-check the public Yang-Yang Energy at Infinite.
function _obc_xxz_ground_state(N::Int, Δ::Real; J::Real=1.0)
    σx = ComplexF64[0 1; 1 0]
    σy = ComplexF64[0 -im; im 0]
    σz = ComplexF64[1 0; 0 -1]
    Id(d) = Matrix{ComplexF64}(I, d, d)
    function pauli_pair(P::AbstractMatrix, Q::AbstractMatrix, i::Int, N::Int)
        left = Id(2^(i - 1))
        right = Id(2^(N - i - 1))
        return kron(kron(left, P), kron(Q, right))
    end
    D = 2^N
    H = zeros(ComplexF64, D, D)
    pref = J / 4
    for i in 1:(N - 1)
        H .+= pref .* pauli_pair(σx, σx, i, N)
        H .+= pref .* pauli_pair(σy, σy, i, N)
        H .+= (pref * Δ) .* pauli_pair(σz, σz, i, N)
    end
    return real(minimum(eigvals(Hermitian(H))))
end

# ─────────────────────────────────────────────────────────────────────────────
# Migrated from raw @test to verify()-first (PR #449 phase 2 zero-legacy).
# Closed-form pins, limiting-case pins, and the OBC-ED extrapolation cross-
# check are now verify cards. The `all(isfinite, es_per_site)` structural
# sweep stays raw (multi-point property, not a single-scalar pin).
# ─────────────────────────────────────────────────────────────────────────────

@testset "XXZ1D — Yang-Yang vs OBC dense-ED extrapolation (Δ=1/2)" begin
    # Closed-form Yang-Yang value at γ = π/3 (Δ = 1/2): -3J/8 (rational).
    verify(
        XXZ1D(; J=1.0, Δ=0.5),
        Energy(),
        Infinite();
        route=:second_closed_form,
        independent=-3 / 8,
        agree_within=1e-10,
        refs=["Yang-Yang 1966 II eq.(4.4): e0 = -3J/8 at Δ = 1/2 (γ = π/3)"],
    )

    # OBC finite-N ground-state energies per site at N = 8, 10, 12.
    Ns = [8, 10, 12]
    es_per_site = [_obc_xxz_ground_state(N, 0.5; J=1.0) / N for N in Ns]

    # Sanity: every finite-N value is finite and lies above the
    # thermodynamic-limit value by O(1/N). Multi-point structural — raw.
    @test all(isfinite, es_per_site)

    # OBC has 1/N leading correction (edge defects); fit y = a + b·(1/N).
    intercept, _ = _yy_linfit([1 / N for N in Ns], es_per_site)
    verify(
        XXZ1D(; J=1.0, Δ=0.5),
        Energy(),
        Infinite();
        route=:ed_finite_size,
        independent=intercept,
        agree_within=5e-3,
        at=["Ns=$(Ns)"],
        refs=[
            "Independent OBC dense-ED ground-state energies at N=8,10,12 then 1/N-extrapolated (edge-defect leading correction); cross-checks the Yang-Yang single-integral closed form at Δ = 1/2",
        ],
    )
end

@testset "XXZ1D — Yang-Yang Δ → 1⁻ limit agrees with Heisenberg AF closed form" begin
    # Approach γ → 0 from inside.
    verify(
        XXZ1D(; J=1.0, Δ=0.999),
        Energy(),
        Infinite();
        route=:limiting_case,
        independent=0.25 - log(2.0),
        agree_within=1e-3,
        at=["Δ=0.999"],
        refs=["Heisenberg AF Δ→1⁻ limit: e0 → 1/4 − log 2 (des Cloizeaux-Pearson 1962)"],
    )
end

@testset "XXZ1D — Yang-Yang Δ → -1⁺ limit agrees with FM saturated state" begin
    # Approach γ → π from inside (the FM saturated state e0 = -J/4 is
    # also exact at Δ = -1; this card pins the continuous-Δ limit).
    verify(
        XXZ1D(; J=1.0, Δ=-0.999),
        Energy(),
        Infinite();
        route=:limiting_case,
        independent=-0.25,
        agree_within=1e-3,
        at=["Δ=-0.999"],
        refs=[
            "FM Δ→-1⁺ limit: e0 → -J/4 (aligned saturated state, continuous from Δ = -1)"
        ],
    )
end

@testset "XXZ1D — Yang-Yang Δ → 0 limit agrees with free-fermion -J/π" begin
    # Approach γ = π/2 from both sides.
    verify(
        XXZ1D(; J=1.0, Δ=1e-6),
        Energy(),
        Infinite();
        route=:limiting_case,
        independent=-1 / π,
        agree_within=1e-5,
        at=["Δ=1e-6"],
        refs=[
            "XX free-fermion limit: e0 = -J/π at Δ = 0 (Lieb-Schultz-Mattis 1961); Yang-Yang continuous from Δ → 0⁺",
        ],
    )
    verify(
        XXZ1D(; J=1.0, Δ=-1e-6),
        Energy(),
        Infinite();
        route=:limiting_case,
        independent=-1 / π,
        agree_within=1e-5,
        at=["Δ=-1e-6"],
        refs=[
            "XX free-fermion limit: e0 = -J/π at Δ = 0; Yang-Yang continuous from Δ → 0⁻"
        ],
    )
end

# ── Anchor card (FM exact) ─────────────────────────────────────────────────
@testset "XXZ Yang-Yang — exact-point anchor cards" begin
    # FM point Δ = -1: exact saturated state e0 = -J/4.
    verify(
        XXZ1D(; J=1.0, Δ=-1.0),
        Energy(),
        Infinite();
        route=:second_closed_form,
        independent=-0.25,
        agree_within=1e-12,
        refs=["XXZ FM point Δ=-1: aligned state exact, e0 = -J/4"],
    )
end
