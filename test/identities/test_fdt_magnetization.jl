# ─────────────────────────────────────────────────────────────────────────────
# Test: the magnetisation fluctuation–dissipation theorem, in two pillars.
#
#   PILLAR ① — autodiff self-consistency (model-independent).  The static
#   cross-susceptibility generalises Layer 1b of test_fluctuation_dissipation.jl:
#   for two observables A, B that are diagonal in the energy eigenbasis (i.e.
#   [H,A]=[H,B]=0) under H(λ) = H₀ - λ·B,
#
#       ∂⟨A⟩/∂λ = β·Cov(A,B)        (= β·Var(A) when A = B)
#
#   is checked to machine precision for arbitrary spectra — the response side by
#   ForwardDiff, the fluctuation side by the ensemble covariance (two genuinely
#   different computations).
#
#   PILLAR ② — verification on a real quantum many-body system.  The XXZ chain
#   `H = (J/4) Σ [σˣσˣ + σʸσʸ + Δ σᶻσᶻ]` conserves `M_z = Σ σᶻᵢ` (U(1)) for
#   every Δ, so the magnetisation FDT `χ_zz = β·Var(M_z)/N` holds.  On a dense-ED
#   spectrum we check it two ways:
#     • INDEPENDENT (non-circular): the response `∂⟨M_z⟩/∂h` (central finite
#       difference — `eigen` is not ForwardDiff-able) equals the fluctuation
#       `β·Var(M_z)`.  ED necessitates finite-diff here; the *autodiff* form of
#       this same relation is pillar ① above.
#     • ATLAS consistency: `β·Var(M_z)/N` reproduces the registered
#       `fetch(XXZ1D, SusceptibilityZZ, OBC)` (same variance convention, #576 —
#       a convention/normalisation check, not an independent one).
#   Plus an SU(2) cross-check at Δ=1: Var(M_x) = Var(M_z) (M_x also conserved).
#
# Conventions pinned to XXZ1D (XXZ.jl / XXZ_thermal.jl): Pauli operators,
# H prefactor J/4, M_z = Σ σᶻ (eigenvalues ±1, NOT S_z = ±½), OBC = N−1 bonds.
# ED helpers reuse embed_single_site / embed_two_site (spinhalf_ed.jl) and the
# fd_boltzmann_weights foundation (fluctuation_dissipation.jl).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, LinearAlgebra, Test, ForwardDiff, Random
using QAtlas: XXZ1D, SusceptibilityZZ, OBC, fetch

# ── ED helpers (Pauli convention matching XXZ1D; `_fdt_`-namespaced) ──────────

const _FDT_SX = ComplexF64[0 1; 1 0]
const _FDT_SY = ComplexF64[0 -im; im 0]
const _FDT_SZ = ComplexF64[1 0; 0 -1]

# XXZ open chain H = (J/4) Σ_{i=1}^{N-1} [σˣσˣ + σʸσʸ + Δ σᶻσᶻ] (Pauli, S=σ/2).
function _fdt_build_xxz_obc(N::Int, J::Real, Δ::Real)
    return (J / 4) * sum(
        embed_two_site(_FDT_SX, _FDT_SX, i, i + 1, N) +
        embed_two_site(_FDT_SY, _FDT_SY, i, i + 1, N) +
        Δ * embed_two_site(_FDT_SZ, _FDT_SZ, i, i + 1, N) for i in 1:(N - 1)
    )
end

# Total longitudinal magnetisation M_z = Σᵢ σᶻᵢ (Pauli; matches XXZ1D's _σz).
_fdt_total_m(σ::AbstractMatrix, N::Int) = sum(embed_single_site(σ, i, N) for i in 1:N)

# Thermal ⟨O⟩ and Var(O) of a Hermitian operator O on a dense Hamiltonian H, via
# the canonical trace in the energy eigenbasis (degeneracy-safe: uses the O and
# O² operator expectations per eigenstate, not an assumed eigenvalue).
function _fdt_op_thermal(H::AbstractMatrix, O::AbstractMatrix, β::Real)
    vals, vecs = eigen(Hermitian(Matrix(H)))
    p = fd_boltzmann_weights(vals, β)
    OV = O * vecs
    m = sum(p[n] * real(dot(view(vecs, :, n), view(OV, :, n))) for n in eachindex(p))
    m2 = sum(p[n] * real(dot(view(OV, :, n), view(OV, :, n))) for n in eachindex(p))
    return (; mean=m, var=max(m2 - m^2, 0.0))
end

# ── PILLAR ① — abstract cross-susceptibility via autodiff ────────────────────
@testset "FDT (abstract) — cross-susceptibility ∂⟨A⟩/∂λ_B = β·Cov(A,B)" begin
    rng = MersenneTwister(0x5151)
    # ⟨A⟩ under H(λ) = H₀ - λ·B for diagonal observables A, B (aligned with E0).
    mean_A(E0, A, B, β, λ) =
        let lev = E0 .- λ .* B
            sum(fd_boltzmann_weights(lev, β) .* A)
        end
    cases = [
        ("A = B  (→ Var)", zeros(5), Float64.(-2:2), Float64.(-2:2)),
        ("random A, B (16)", sort(5 .* rand(rng, 16)), randn(rng, 16), randn(rng, 16)),
        ("A ≠ B structured (8)", collect(0.0:0.5:3.5), Float64.(1:8), Float64.(8:-1:1)),
    ]
    for (name, E0, A, B) in cases, β in (0.3, 1.0, 2.0), λ in (-0.5, 0.0, 0.6)
        dAdλ = ForwardDiff.derivative(l -> mean_A(E0, A, B, β, l), λ)
        p = fd_boltzmann_weights(E0 .- λ .* B, β)
        cov = sum(p .* A .* B) - sum(p .* A) * sum(p .* B)
        @test isapprox(dAdλ, β * cov; rtol=1e-7, atol=1e-9)
    end
end

# ── PILLAR ② — real quantum system (XXZ dense ED) ─────────────────────────────
@testset "FDT (real system) — XXZ1D M_z FDT via ED + atlas tie-in" begin
    for Δ in (1.0, 0.5), N in (6, 8)
        J = 1.0
        H0 = _fdt_build_xxz_obc(N, J, Δ)
        Mz = _fdt_total_m(_FDT_SZ, N)
        for β in (0.4, 1.0)
            th = _fdt_op_thermal(H0, Mz, β)
            @test isapprox(th.mean, 0.0; atol=1e-8)        # ⟨M_z⟩ = 0 at h=0 (Sᶻ sym)

            # INDEPENDENT FDT: response ∂⟨M_z⟩/∂h (central diff) == β·Var(M_z).
            # H(h) = H₀ − h·M_z; the two routes share no computation.
            δ = 1e-4
            mp = _fdt_op_thermal(H0 .- δ .* Mz, Mz, β).mean
            mm = _fdt_op_thermal(H0 .+ δ .* Mz, Mz, β).mean
            dMdh = (mp - mm) / (2δ)
            @test isapprox(dMdh, β * th.var; rtol=1e-5, atol=1e-7)

            # ATLAS consistency: β·Var(M_z)/N reproduces the registered χ_zz.
            χ_fetch = fetch(XXZ1D(; J=J, Δ=Δ), SusceptibilityZZ(), OBC(N); beta=β)
            @test isapprox(β * th.var / N, χ_fetch; rtol=1e-8, atol=1e-10)
        end
    end

    # SU(2) bonus at Δ=1: M_x is also conserved, so Var(M_x) = Var(M_z) — an
    # independent symmetry cross-check requiring no fetch.
    let N = 6, J = 1.0, β = 1.0
        H0 = _fdt_build_xxz_obc(N, J, 1.0)
        @test isapprox(
            _fdt_op_thermal(H0, _fdt_total_m(_FDT_SX, N), β).var,
            _fdt_op_thermal(H0, _fdt_total_m(_FDT_SZ, N), β).var;
            rtol=1e-9,
        )
    end
end
