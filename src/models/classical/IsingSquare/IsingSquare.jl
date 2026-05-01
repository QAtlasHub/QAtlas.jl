# ─────────────────────────────────────────────────────────────────────────────
# Classical 2D Ising model on the square lattice — exact partition function
#
# Hamiltonian:
#   H = -J Σ_{⟨i,j⟩} σᵢ σⱼ,   σᵢ ∈ {−1, +1}
#
# Solved exactly via the transfer-matrix method. For an Lx × Ly torus
# the partition function is Z = Tr(T^Lx) where T is the 2^Ly × 2^Ly
# transfer matrix.
#
# References:
#   L. Onsager, "Crystal Statistics. I. A Two-Dimensional Model with an
#   Order-Disorder Transition", Phys. Rev. 65, 117 (1944).
#   B. M. McCoy and T. T. Wu, "The Two-Dimensional Ising Model",
#   Harvard University Press (1973).
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: tr

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tags
# ═══════════════════════════════════════════════════════════════════════════════

"""
    IsingSquare(; J::Real = 1.0, Lx::Int = 0, Ly::Int = 0) <: AbstractQAtlasModel

Classical 2D Ising model on a square lattice with periodic boundary
conditions (PBC) in both directions.

Hamiltonian: H = -J Σ_{⟨i,j⟩} σᵢ σⱼ, σᵢ ∈ {-1, +1}.

Physics parameters are carried as typed struct fields:

- `J`   — Ising coupling (default `1.0`, `J > 0` ferromagnetic).
- `Lx`, `Ly` — lattice extents.  `0` is a legacy sentinel meaning
  "thermodynamic limit / unspecified"; finite-size quantities like
  [`PartitionFunction`](@ref) require both to be positive.

For backward compatibility the `fetch` methods accept `Lx`, `Ly`, `J`,
`β` as kwargs too — kwargs override the struct fields when both are
supplied.

See also: [`PartitionFunction`](@ref), [`CriticalTemperature`](@ref),
[`SpontaneousMagnetization`](@ref).
"""
struct IsingSquare <: AbstractQAtlasModel
    J::Float64
    Lx::Int
    Ly::Int
end
function IsingSquare(; J::Real=1.0, Lx::Integer=0, Ly::Integer=0)
    IsingSquare(Float64(J), Int(Lx), Int(Ly))
end

"""
    PartitionFunction() <: AbstractQuantity

Canonical partition function `Z = Σ_σ exp(-β H(σ))`.
"""
struct PartitionFunction <: AbstractQuantity end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: transfer matrix for a row of Ly spins (PBC in y)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _ising_sq_transfer_matrix(Ly, β, J) -> Matrix

Return the 2^Ly × 2^Ly symmetric transfer matrix for a row of `Ly` Ising
spins with PBC along the row direction.

The symmetric (Hermitian) form is defined by:

    T_{σ,σ'} = exp(βJ/2 · Eₕ(σ)) · exp(βJ · Eᵥ(σ,σ')) · exp(βJ/2 · Eₕ(σ'))

where
    Eₕ(σ) = Σⱼ₌₁^Ly σⱼ σ_{(j mod Ly)+1}   (horizontal bonds within row, PBC)
    Eᵥ(σ,σ') = Σⱼ₌₁^Ly σⱼ σ'ⱼ             (vertical bonds between rows)

Note: for Ly = 2, each horizontal bond is counted twice by the PBC sum
(σ₁σ₂ + σ₂σ₁ = 2σ₁σ₂). The same applies along the transfer direction:
`tr(T^Lx)` double-counts the vertical bonds when Lx = 2, because the
cyclic product `T[σ,σ'] T[σ',σ]` weights the single pair of rows with
the vertical Boltzmann factor twice. This doubling is an intrinsic
property of the PBC transfer-matrix formula at L = 2 and is used
consistently by the brute-force cross-check in
`test/util/classical_partition.jl` (its bond list is constructed to
match this convention by enumerating
`(i, j) ↔ (i, (j mod Ly) + 1)` and `(i, j) ↔ ((i mod Lx) + 1, j)`
for every site).

The function is generic in `β` and `J` so that automatic-differentiation
dual numbers (e.g. `ForwardDiff.Dual`) propagate through the matrix
elements. The element type of the returned matrix is inferred from
`typeof(exp(β * J))`.
"""
function _ising_sq_transfer_matrix(Ly::Int, β::Real, J::Real)
    dim = 2^Ly

    # Element type of Boltzmann weights — propagates Dual numbers for AD.
    TT = typeof(exp(β * J))

    # Precompute spin vectors for each row state index (0-based binary encoding)
    spins_of = Vector{Vector{Int}}(undef, dim)
    for σ_idx in 0:(dim - 1)
        spins_of[σ_idx + 1] = Int[((σ_idx >> j) & 1) == 1 ? 1 : -1 for j in 0:(Ly - 1)]
    end

    # Diagonal weight: exp(βJ/2 · Eₕ(σ))
    h_weight = Vector{TT}(undef, dim)
    for σ_idx in 0:(dim - 1)
        σ = spins_of[σ_idx + 1]
        e_h = sum(σ[j] * σ[(j % Ly) + 1] for j in 1:Ly)
        h_weight[σ_idx + 1] = exp(β * J / 2 * e_h)
    end

    # Build transfer matrix
    T = Matrix{TT}(undef, dim, dim)
    for σ_idx in 0:(dim - 1)
        σ = spins_of[σ_idx + 1]
        wσ = h_weight[σ_idx + 1]
        for σp_idx in 0:(dim - 1)
            σp = spins_of[σp_idx + 1]
            e_v = sum(σ[j] * σp[j] for j in 1:Ly)
            T[σ_idx + 1, σp_idx + 1] = wσ * exp(β * J * e_v) * h_weight[σp_idx + 1]
        end
    end
    return T
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: partition function via transfer matrix
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingSquare, ::PartitionFunction; Lx, Ly, β, J=1.0) -> Real

Exact partition function Z = Tr(T^Lx) for the classical 2D Ising model on an
Lx × Ly square lattice with periodic boundary conditions in both directions.

The transfer matrix T acts along the Lx direction (row-to-row transfer), with
each row containing Ly spins and PBC along the Ly direction.

# Bond-counting convention

The transfer-matrix sum double-counts bonds along any dimension of length
2 (PBC wraparound of a length-2 ring produces the factor
`σ_1 σ_2 + σ_2 σ_1 = 2 σ_1 σ_2`). The brute-force enumeration in
`test/util/classical_partition.jl` is built under the same PBC convention
so that `Z_transfer-matrix == Z_bruteforce` exactly for every
`(Lx, Ly, β, J)`. For `Lx ≥ 3` and `Ly ≥ 3` each physical bond is
enumerated exactly once and the result coincides with the standard
physical Z of a PBC lattice.

# Special values

- β = 0 (any Lx, Ly, J): Z = 2^(Lx·Ly)  — all configurations equally weighted
- J = 0 (any β, Lx, Ly): Z = 2^(Lx·Ly)  — no interactions, same as β = 0

# Automatic differentiation

`fetch` is generic in `β` and `J` so that `ForwardDiff.Dual` numbers
propagate through it. This allows macroscopic thermodynamic quantities
to be recovered from Z by differentiation — e.g.

    ⟨E⟩ = -∂(log Z)/∂β
    C_v = β² · ∂²(log Z)/∂β²

See `test/verification/test_ising_ad_thermodynamics.jl` for a cross-check
against direct ensemble averages.

# Arguments
- `Lx::Int`: number of rows (transfer direction)
- `Ly::Int`: number of columns (row length, PBC)
- `β::Real`: inverse temperature (β = 1/(k_B T))
- `J::Real`: Ising coupling constant (default 1.0; J > 0 ferromagnetic)

# References
    L. Onsager, Phys. Rev. 65, 117 (1944).
    B. M. McCoy and T. T. Wu, "The Two-Dimensional Ising Model" (1973).
"""
function fetch(
    m::IsingSquare,
    ::PartitionFunction;
    β::Real,
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    J::Real=m.J,
)
    Lx > 0 && Ly > 0 || error(
        "IsingSquare PartitionFunction: Lx and Ly must be positive. " *
        "Pass them in the struct (IsingSquare(; Lx, Ly, J)) or as kwargs.",
    )
    T = _ising_sq_transfer_matrix(Ly, β, J)
    return tr(T^Lx)
end

# BC-aware delegator: required by the registry drift guard so the
# `(IsingSquare, PartitionFunction, PBC)` triple resolves to a
# non-catch-all fetch method.  PBC is the only natural BC for the
# transfer-matrix Z; OBC (open horizontal axis) would simply omit the
# wrap-around bond and is not registered.
function fetch(
    m::IsingSquare,
    q::PartitionFunction,
    ::PBC;
    β::Real,
    Lx::Integer=m.Lx,
    Ly::Integer=m.Ly,
    J::Real=m.J,
)
    return fetch(m, q; β=β, Lx=Lx, Ly=Ly, J=J)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tags: Onsager + Yang exact results
# ═══════════════════════════════════════════════════════════════════════════════

"""
    CriticalTemperature() <: AbstractQuantity

Exact critical temperature of a classical model.
"""
struct CriticalTemperature <: AbstractQuantity end

"""
    SpontaneousMagnetization() <: AbstractQuantity

Spontaneous magnetization of a classical model as a function of
temperature.  Retained under this name for backward compatibility
with the classical-Ising literature; new code may prefer the
axis-explicit [`MagnetizationZ`](@ref) together with a `T < T_c`
context.
"""
struct SpontaneousMagnetization <: AbstractQuantity end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: Onsager critical temperature
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingSquare, ::CriticalTemperature; J=1.0) -> Float64

Exact critical temperature of the 2D Ising model on the square lattice:

    T_c = 2J / ln(1 + √2) ≈ 2.269 J

Equivalently, the critical reduced coupling is K_c = J/T_c = ln(1+√2)/2,
or sinh(2K_c) = 1.

# References
    L. Onsager, "Crystal Statistics. I.", Phys. Rev. 65, 117 (1944).
"""
function fetch(m::IsingSquare, ::CriticalTemperature; J::Real=m.J)
    return 2J / log(1 + sqrt(2))
end

# BC-aware delegator: registry drift-guard companion.  `T_c` is a
# property of the thermodynamic limit, so the natural BC is `Infinite`.
function fetch(m::IsingSquare, q::CriticalTemperature, ::Infinite; J::Real=m.J)
    return fetch(m, q; J=J)
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: Yang spontaneous magnetization
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingSquare, ::SpontaneousMagnetization; β, J=1.0) -> Float64

Exact spontaneous magnetization of the 2D Ising model on the infinite
square lattice:

    M(T) = (1 − sinh⁻⁴(2βJ))^{1/8}    for T < T_c  (i.e. sinh(2βJ) > 1)
    M(T) = 0                            for T ≥ T_c

The critical exponent β = 1/8 is visible in the approach M → 0 as
T → T_c⁻.

Special values:
- T = 0 (β → ∞): M = 1 (fully ordered)
- T = T_c: M = 0 (onset of disorder)

# Arguments
- `β::Real`: inverse temperature (β = 1/(k_B T))
- `J::Real`: Ising coupling constant (default 1.0; J > 0 ferromagnetic)

# References
    C. N. Yang, "The spontaneous magnetization of a two-dimensional Ising
    model", Phys. Rev. 85, 808 (1952).
"""
function fetch(m::IsingSquare, ::SpontaneousMagnetization; β::Real, J::Real=m.J)
    s = sinh(2 * β * J)
    if s <= 1.0  # T ≥ T_c
        return 0.0
    else
        return (1 - s^(-4))^(1 / 8)
    end
end

# BC-aware delegator: registry drift-guard companion.  Yang's
# spontaneous magnetization is a thermodynamic-limit quantity, so the
# natural BC is `Infinite`.
function fetch(
    m::IsingSquare, q::SpontaneousMagnetization, ::Infinite; β::Real, J::Real=m.J
)
    return fetch(m, q; β=β, J=J)
end
