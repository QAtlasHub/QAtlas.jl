# ─────────────────────────────────────────────────────────────────────────────
# Spin-1 Heisenberg chain (Haldane chain)
#
# Hamiltonian:
#   H = J Σᵢ Sᵢ · Sᵢ₊₁
#     = J Σᵢ [Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Sᶻᵢ Sᶻᵢ₊₁]
#
# where the Sᵅ are 3×3 spin-1 operators (s = 1, s(s+1) = 2).
#
# Bond eigenvalue summary (Sᵢ + Sᵢ₊₁ → S_tot ∈ {0, 1, 2} from two spin-1):
#
#   J Sᵢ · Sᵢ₊₁ = (J/2)·(S_tot² − 4)
#
# giving bond eigenvalues  [-2J (singlet), -J (triplet), +J (quintet)].
#
# This file is the small-N dense-ED reference path for the Haldane phase
# (gap Δ ≈ 0.41 J, ξ ≈ 6 lattice spacings) — used to validate ThermalMPS
# computations at modest sizes (3^N Hilbert space, capped at N ≤ 8).
#
# References:
#   F. D. M. Haldane, Phys. Lett. A 93, 464 (1983); PRL 50, 1153 (1983).
#   T. Kennedy and H. Tasaki, Phys. Rev. B 45, 304 (1992) — string order.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: I, Diagonal, Hermitian, eigen, eigvals, kron, tr

"""
    _MAX_ED_SITES_S1 = 8

Hard cap on chain length for spin-1 dense-ED helpers.  3⁸ = 6561 → a
6561×6561 Hermitian eigendecomposition is ~1 s on a laptop; 3⁹ ≈ 20 k
already pushes 5 GB of workspace, past the "cheap reference" regime.
"""
const _MAX_ED_SITES_S1 = 8

# Spin-1 generators in the |+1⟩, |0⟩, |-1⟩ basis.
const _S1_id = ComplexF64[1 0 0; 0 1 0; 0 0 1]
const _S1_x = (ComplexF64(1) / sqrt(2)) * ComplexF64[0 1 0; 1 0 1; 0 1 0]
const _S1_y = (ComplexF64(1) / (im * sqrt(2))) * ComplexF64[0 1 0; -1 0 1; 0 -1 0]
const _S1_z = ComplexF64[1 0 0; 0 0 0; 0 0 -1]

"""
    _spin1_string(N::Int, site_ops::Pair{Int,Matrix{ComplexF64}}...) -> Matrix{ComplexF64}

Tensor-product analogue of [`_pauli_string`](@ref) for spin-1 chains.
Returns the `3^N × 3^N` operator that places each `Sᵅ` at its listed
site and the 3×3 identity elsewhere.
"""
function _spin1_string(N::Int, site_ops::Pair{Int,Matrix{ComplexF64}}...)
    N ≤ _MAX_ED_SITES_S1 || throw(
        ArgumentError("spin-1 dense ED is capped at N ≤ $(_MAX_ED_SITES_S1) (got N = $N)"),
    )
    lookup = Dict(site_ops...)
    ops = [get(lookup, k, _S1_id) for k in 1:N]
    return reduce(kron, ops)
end

"""
    S1Heisenberg1D(; J::Real = 1.0) <: AbstractQAtlasModel

Spin-1 antiferromagnetic Heisenberg chain (Haldane chain),

    H = J Σᵢ Sᵢ · Sᵢ₊₁,    spin = 1,  J > 0 antiferromagnetic.

Distinct from [`Heisenberg1D`](@ref) (spin-1/2) because the spin
representation differs: spin-1 has local dimension 3 and a gapped,
topologically non-trivial Haldane phase, while spin-1/2 is gapless and
critical (Bethe ansatz).
"""
struct S1Heisenberg1D <: AbstractQAtlasModel
    J::Float64
end
S1Heisenberg1D(; J::Real=1.0) = S1Heisenberg1D(Float64(J))

"""
    _s1_heisenberg_hamiltonian_matrix(model::S1Heisenberg1D, N::Int) -> Matrix{ComplexF64}

Assemble the `3^N × 3^N` OBC Hamiltonian

    H = J Σᵢ [Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Sᶻᵢ Sᶻᵢ₊₁]

via explicit tensor products built from the spin-1 primitives.
Capped by `_MAX_ED_SITES_S1`.
"""
function _s1_heisenberg_hamiltonian_matrix(model::S1Heisenberg1D, N::Int)
    N ≥ 2 || throw(ArgumentError("S1Heisenberg1D OBC chain needs N ≥ 2 (got N = $N)"))
    N ≤ _MAX_ED_SITES_S1 || throw(
        ArgumentError("spin-1 dense ED is capped at N ≤ $(_MAX_ED_SITES_S1) (got N = $N)"),
    )
    J = model.J
    D = 3^N
    # Combine the three Sᵅ⊗Sᵅ contractions into a single 9×9 bond block
    # so each bond costs a single kron instead of three.
    bond = J * (kron(_S1_x, _S1_x) + kron(_S1_y, _S1_y) + kron(_S1_z, _S1_z))
    H = zeros(ComplexF64, D, D)
    for i in 1:(N - 1)
        d_left = 3^(i - 1)
        d_right = 3^(N - i - 1)
        H .+= kron(
            Matrix{ComplexF64}(I, d_left, d_left),
            bond,
            Matrix{ComplexF64}(I, d_right, d_right),
        )
    end
    return H
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch — finite-N OBC thermal observables via dense ED
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::S1Heisenberg1D, ::OBC) = :total

"""
    fetch(model::S1Heisenberg1D, ::Energy{:total}, ::OBC; beta) -> Float64

**Total** thermal energy `⟨H⟩_β` of the spin-1 OBC Heisenberg chain at
finite `N ≤ $(_MAX_ED_SITES_S1)`, computed by dense ED.

```
    ⟨H⟩_β = Tr(H exp(-βH)) / Tr(exp(-βH))
```

Convention matches [`fetch(::TFIM, ::Energy{:total}, ::OBC)`](@ref):
1D finite-size boundary conditions return total energy natively.
"""
function fetch(model::S1Heisenberg1D, ::Energy{:total}, bc::OBC; beta::Real, kwargs...)
    H = _s1_heisenberg_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_energy(H, beta)
end

"""
    fetch(model::S1Heisenberg1D, ::FreeEnergy, ::OBC; beta) -> Float64

Per-site Helmholtz free energy `f(β) = -log Z / (Nβ)` for the spin-1
OBC chain.
"""
function fetch(model::S1Heisenberg1D, ::FreeEnergy, bc::OBC; beta::Real, kwargs...)
    H = _s1_heisenberg_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_free_energy(H, beta) / bc.N
end

"""
    fetch(model::S1Heisenberg1D, ::ThermalEntropy, ::OBC; beta) -> Float64

Per-site Gibbs entropy `s(β) = β · (ε - f)` for the spin-1 OBC chain.
"""
function fetch(model::S1Heisenberg1D, ::ThermalEntropy, bc::OBC; beta::Real, kwargs...)
    H = _s1_heisenberg_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_entropy(H, beta) / bc.N
end

"""
    fetch(model::S1Heisenberg1D, ::SpecificHeat, ::OBC; beta) -> Float64

Per-site heat capacity `c(β) = β² · Var(H) / N` for the spin-1 OBC
chain, computed exactly from the energy variance in the eigenbasis (no
numerical differentiation).
"""
function fetch(model::S1Heisenberg1D, ::SpecificHeat, bc::OBC; beta::Real, kwargs...)
    H = _s1_heisenberg_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_specific_heat(H, beta) / bc.N
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal kernel — one eigendecomposition shared across many observables
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _s1_thermal_kernel(model, N, beta) -> NamedTuple

Compute the eigendecomposition of the OBC spin-1 Heisenberg Hamiltonian
on `N` sites and return everything downstream observables need to evaluate
thermal expectation values:

- `H`        — the dense `3^N × 3^N` Hermitian Hamiltonian matrix
- `evals`    — sorted real eigenvalues (`F.values`)
- `evecs`    — eigenvector matrix (`F.vectors`)
- `weights`  — Boltzmann weights `wₙ = exp(-β(eₙ - emin))/Z̃`,
                normalised so `sum(weights) == 1`
- `ρ`        — `evecs · diagm(weights) · evecs'`, the thermal density
                matrix in the computational (product) basis

For `beta = Inf` we collapse onto the (possibly degenerate) ground-state
manifold by zeroing all weights above `evals[1] + 1e-12`.

Constructed once per fetch call (so a routine that needs both an energy
and a magnetisation does not pay for two `eigen` calls).  No caching
across calls — that would require a per-model mutable cache and hasn't
been profiled to matter yet.
"""
function _s1_thermal_kernel(model::S1Heisenberg1D, N::Int, beta::Real)
    H = _s1_heisenberg_hamiltonian_matrix(model, N)
    F = eigen(Hermitian(H))
    evals = F.values
    evecs = F.vectors
    if isinf(beta) && beta > 0
        emin = evals[1]
        # Ground-state projector: equal weight on the (numerically) degenerate
        # ground manifold.  `1e-12` is well below typical N ≤ 8 spectrum gaps
        # for the antiferromagnetic spin-1 chain (smallest finite-N gap ≳ 0.1).
        gs_mask = (evals .- emin) .≤ 1e-12
        ng = count(gs_mask)
        weights = zeros(Float64, length(evals))
        weights[gs_mask] .= 1.0 / ng
    else
        emin = minimum(evals)
        ws = exp.(-beta .* (evals .- emin))
        Z = sum(ws)
        weights = ws ./ Z
    end
    # Density matrix in product basis: ρ = U diag(w) U'
    ρ = evecs * Diagonal(weights) * evecs'
    return (; H=H, evals=evals, evecs=evecs, weights=weights, ρ=ρ)
end

"""
    _s1_thermal_expectation(kernel, O) -> Float64

Trace `Tr(O ρ)` from a precomputed `_s1_thermal_kernel` and a dense
operator `O` of matching size.  Returns the real part (Hermitian
operators give real expectation values; small imaginary residue is
numerical noise).
"""
function _s1_thermal_expectation(kernel, O::AbstractMatrix)
    return real(tr(O * kernel.ρ))
end
