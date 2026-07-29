# ─────────────────────────────────────────────────────────────────────────────
# Small-N dense exact diagonalisation helpers for spin-1/2 chains
#
# For integrable-at-infinite models like XXZ1D that don't yet have a
# finite-size finite-T formula implemented, the ED path gives a
# *finite-N exact* thermal reference suitable for benchmarking MPS
# methods at modest system sizes.  The cost is `O(2^{3N})` in time and
# `O(2^{2N})` in memory, so the scale ceiling is ~N ≤ 12.
#
# Performance note: every observable ⟨O⟩_β goes through the same
# eigendecomposition of H, so callers that want multiple quantities at
# the same β should cache the spectrum.  The public `fetch` methods do
# not cache today — add a per-model cache when/if profiling shows it
# matters.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: Hermitian, eigen

"""
    _MAX_ED_SITES = 12

Hard cap on chain length for dense-ED helpers.  2^12 = 4096 → 4096×4096
Hermitian eigendecomposition takes ~1 s on a laptop; 2^14 ~ 16384 is
already several seconds and 64 GB+ of workspace, which is past the
"cheap reference data" regime this helper targets.
"""
const _MAX_ED_SITES = 12

# Single-qubit Pauli matrices used to assemble tensor-product operators.
const _σ0 = ComplexF64[1 0; 0 1]
const _σx = ComplexF64[0 1; 1 0]
const _σy = ComplexF64[0 -im; im 0]
const _σz = ComplexF64[1 0; 0 -1]

"""
    _ed_size_guard(N::Int, cap::Int, d::Int, what::AbstractString)

The single gate every many-body dense-ED path in this atlas passes through.
Throws above `cap`, and says ONCE per session that the route is exponential.

WHY A WARNING AND NOT A REMOVAL.  An exponential route is still a correct route,
and for a model with no closed form it is the only one there is.  What it must
not do is look like the scalable ones: a caller who asks a hub for a value and
gets a `2^N` diagonalisation should be told, because the answer stops existing a
few sites later.  So the atlas keeps the route, registers the cost
(`cost=:exponential`, `max_size`) as data rather than prose, and says so at the
moment the cost is actually incurred.

`maxlog=1` — this fires inside the innermost operator builder, so without it a
single sweep would emit thousands of identical lines.  The registry carries the
same fact statically for anyone who wants to plan ahead instead of being told
afterwards.
"""
function _ed_size_guard(N::Int, cap::Int, d::Int, what::AbstractString)
    N ≤ cap || throw(
        ArgumentError(
            "$what: many-body dense ED is capped at N ≤ $cap (got N = $N). The cost " *
            "is O($(d)^N) in both memory and time, so this is a hard limit, not a " *
            "tuning knob — use a route registered with cost=:polynomial or " *
            ":closed_form if the model has one.",
        ),
    )
    # Build the message first.  `@warn msg key=value...` parses the FIRST
    # argument as the message and everything after it as key/value pairs, so a
    # `*`-concatenated literal followed by `maxlog=1` is ambiguous to read and
    # fragile under reformatting — the formatter split `maxlog =` from its `1`
    # and the macro stopped parsing.
    msg =
        "$what: many-body dense ED — O($(d)^N), so only small systems are reachable " *
        "(cap N ≤ $cap). Registered as cost=:exponential with max_size=$cap; where a " *
        "scalable route exists it is the canonical one."
    @warn msg N = N maxlog = 1
    return nothing
end

"""
    _pauli_string(N::Int, site_ops::Pair{Int,Matrix{ComplexF64}}...) -> Matrix{ComplexF64}

Return the `2^N × 2^N` tensor product that places each `σ_k` at its
listed site and the identity elsewhere.  `site_ops` is an arbitrary
sequence of `(site => matrix)` pairs; missing sites get `_σ0`.

```julia
# σˣ at site 3 in an N=5 chain:
_pauli_string(5, 3 => _σx)

# σˣᵢ σᶻⱼ for (i,j) = (2,4):
_pauli_string(5, 2 => _σx, 4 => _σz)
```
"""

function _pauli_string(N::Int, site_ops::Pair{Int,Matrix{ComplexF64}}...)
    _ed_size_guard(N, _MAX_ED_SITES, 2, "spin-1/2 dense ED")
    lookup = Dict(site_ops...)
    ops = [get(lookup, k, _σ0) for k in 1:N]
    return reduce(kron, ops)
end

"""
    _ed_thermal_expectation(evals, evecs, O_diag, β) -> Float64

Return `⟨O⟩_β = Σ_n O_nn exp(-β eₙ) / Σ_n exp(-β eₙ)` given a
pre-computed eigendecomposition `(evals, evecs)` of `H` and the
diagonal of `O` in the H-eigenbasis, `O_diag[n] = ⟨n|O|n⟩`.

Uses a log-sum-exp shift (`eₘᵢₙ`) so `β` up to ~10^2 / bandwidth is
numerically stable.
"""
function _ed_thermal_expectation(evals::AbstractVector, O_diag::AbstractVector, β::Real)
    emin = minimum(evals)
    ws = exp.(-β .* (evals .- emin))
    return real(sum(O_diag .* ws) / sum(ws))
end

"""
    _ed_thermal_energy(H::AbstractMatrix, β::Real) -> Float64

`⟨H⟩_β = Tr(H exp(-βH)) / Tr(exp(-βH))`.  `H` must be Hermitian.
"""
function _ed_thermal_energy(H::AbstractMatrix, β::Real)
    F = eigen(Hermitian(H))
    return _ed_thermal_expectation(F.values, F.values, β)
end

"""
    _ed_log_partition(evals::AbstractVector, β::Real) -> Float64

`log Z(β) = log Σ exp(-β eₙ)` with an `emin` shift to avoid overflow at
large `β`:

    log Z = -β eₘᵢₙ + log Σ exp(-β (eₙ - eₘᵢₙ))
"""
function _ed_log_partition(evals::AbstractVector, β::Real)
    emin = minimum(evals)
    return -β * emin + log(sum(exp.(-β .* (evals .- emin))))
end

"""
    _ed_thermal_free_energy(H::AbstractMatrix, β::Real) -> Float64

Total Helmholtz free energy `F(β) = -log Z / β` from the eigenspectrum
of a Hermitian `H`.  Per-site values: divide by `N` at the call site.
"""
function _ed_thermal_free_energy(H::AbstractMatrix, β::Real)
    evals = eigvals(Hermitian(H))
    return -_ed_log_partition(evals, β) / β
end

"""
    _ed_thermal_entropy(H::AbstractMatrix, β::Real) -> Float64

Total Gibbs entropy `S(β) = β·(⟨H⟩ - F) = -Σ wₙ log wₙ`, with
`wₙ = exp(-β(eₙ - eₘᵢₙ)) / Z̃`.  The first form is what we evaluate
because both pieces share the eigendecomposition.
"""
function _ed_thermal_entropy(H::AbstractMatrix, β::Real)
    evals = eigvals(Hermitian(H))
    E = _ed_thermal_expectation(evals, evals, β)
    F = -_ed_log_partition(evals, β) / β
    return β * (E - F)
end

"""
    _ed_thermal_specific_heat(H::AbstractMatrix, β::Real) -> Float64

Total heat capacity from the energy variance,

    C(β) = β² · (⟨H²⟩ - ⟨H⟩²).

Equivalent to `-β² · ∂⟨H⟩/∂β`; computed exactly from the eigenspectrum
without numerical differentiation.
"""
function _ed_thermal_specific_heat(H::AbstractMatrix, β::Real)
    evals = eigvals(Hermitian(H))
    emin = minimum(evals)
    ws = exp.(-β .* (evals .- emin))
    Z = sum(ws)
    E1 = sum(evals .* ws) / Z
    E2 = sum(evals .^ 2 .* ws) / Z
    return β^2 * (E2 - E1^2)
end
