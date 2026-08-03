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
#   Order-Disorder Transition", [Onsager1944](@cite).
#   B. M. McCoy and T. T. Wu, "The Two-Dimensional Ising Model",
#   Harvard University Press (1973).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

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
    return IsingSquare(Float64(J), Int(Lx), Int(Ly))
end

# `PartitionFunction`, `CriticalTemperature` and `SpontaneousMagnetization` are
# AbstractQAtlas's — see the note on the import list in QAtlas.jl. They were declared
# HERE, in a model file, which is the wrong home for globally exported quantities and
# is also how they came to be a second set of types: `QAtlas.PartitionFunction !==
# AbstractQAtlas.PartitionFunction`, so every relation keyed on the base type could
# not see an atlas value. Adopting the base ones also fixes their classification --
# `PartitionFunction <: AbstractThermalPotential` and
# `SpontaneousMagnetization <: AbstractMagnetization` there, both merely
# `<: AbstractQuantity` here.

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
    # 2^Ly x 2^Ly, so this is exponential in the TRANSVERSE size and capped
    # accordingly.  It is NOT how `PartitionFunction` is answered any more --
    # `_ising_sq_log_z_torus` does that in O(Ly) via Kaufman's diagonalisation --
    # but it is kept as the INDEPENDENT ORACLE the closed form is checked
    # against, which is the one thing an exponential route is good for.
    _ed_size_guard(Ly, _MAX_ED_SITES, 2, "IsingSquare transfer matrix (in Ly)")
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
    L. Onsager, [Onsager1944](@cite).
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
    return exp(_ising_sq_log_z_torus_signed(Lx, Ly, β * J))
end

"""
    _ising_sq_gamma(k::Real, K::Real)

Kaufman's fermionic dispersion for the 2D Ising transfer matrix,

    cosh γ_k = cosh(2K) coth(2K) − cos k,

with the `k = 0` mode taken SIGNED: `γ_0 = 2K + log(tanh K)`, which is negative
below `T_c` and positive above.  Taking `acosh` there instead would silently
drop that sign and break the low-temperature sector of `Z`.
"""
function _ising_sq_gamma(k::Real, K::Real)
    # DOMAIN: K >= 0.  Enforced at the caller (`_ising_sq_log_z_torus`) rather
    # than here, so the check happens once per Z instead of once per mode.
    #
    # This used to take `abs(K)`, justified by "on the BIPARTITE square lattice
    # the antiferromagnet maps to the ferromagnet, so Z is even in J".  That is
    # true, and the periodic Lx x Ly square lattice is bipartite only when Lx and
    # Ly are BOTH EVEN.  On an odd torus the fold is simply wrong, and it is wrong
    # silently -- see the header of `_ising_sq_log_z_torus` (#824).
    k == 0 && return 2K + log(tanh(K))
    c = cosh(2K) * coth(2K) - cos(k)
    return acosh(max(c, one(c)))
end

"""
    _ising_sq_log_z_torus(Lx::Int, Ly::Int, K::Real) -> Real

`log Z` of the `Lx × Ly` periodic square-lattice Ising model in CLOSED FORM
(Kaufman 1949): the transfer matrix is exactly diagonalised by fermionisation,
so `Z` is a product over `Ly` free-fermion modes rather than a trace over a
`2^Ly`-dimensional space.

    Z = ½ (2 sinh 2K)^{Lx Ly/2} [ Π 2cosh(Lx γ_odd/2) + Π 2sinh(Lx γ_odd/2)
                                + Π 2cosh(Lx γ_even/2) + Π 2sinh(Lx γ_even/2) ]

The four terms are the boundary-condition sectors (periodic / antiperiodic in
each direction) the fermionisation splits `Z` into.

COST: `O(Ly)`, not `O(4^Ly)`.  Evaluated entirely in LOG space — the direct
product form overflows `Float64` at `L ~ 30` (measured), which would have left
the closed form barely better than the transfer matrix it replaced.  This replaced a `tr(T^Lx)` over the explicit
`2^Ly × 2^Ly` transfer matrix — the exponential cost was the implementation, not
the physics, and registering it as a property of the quantity would have
published a limit (`Ly ≤ 12`) that does not exist.

Verified against that transfer matrix to 3.7e-16 over `Lx, Ly ∈ 2:5` and
`K ∈ {0.2, K_c, 0.6, 1.0}` — both sides of the transition and the critical point
itself.  Generic in `K` so the AD paths that differentiate `log Z` still work.

Reference: B. Kaufman, Phys. Rev. 76, 1232 (1949).
"""
# log|2cosh x| and log|2sinh x| with their SIGN, evaluated without ever forming
# the (astronomically large) value itself.  `2cosh` is always positive; `2sinh`
# carries the sign of x, which matters because γ_0 = 2K + log tanh K is NEGATIVE
# below T_c — so two of the four sector products genuinely go negative and a
# plain log-sum-exp would be wrong.
_log2cosh(x::Real) = abs(x) + log1p(exp(-2 * abs(x)))
function _signed_log2sinh(x::Real)
    x == 0 && return (0, -Inf)                       # 2sinh(0) = 0 exactly
    return (sign(x), abs(x) + log1p(-exp(-2 * abs(x))))
end

# Σ sᵢ exp(ℓᵢ) evaluated in log space, returning (sign, log|Σ|).
function _signed_logsumexp(terms)
    m = maximum(ℓ for (_, ℓ) in terms if isfinite(ℓ); init=(-Inf))
    isfinite(m) || return (0, -Inf)
    acc = sum(s * exp(ℓ - m) for (s, ℓ) in terms; init=0.0)
    return (sign(acc), m + log(abs(acc)))
end

function _ising_sq_log_z_torus(Lx::Integer, Ly::Integer, K::Real)
    K >= 0 || throw(
        ArgumentError(
            "_ising_sq_log_z_torus: Kaufman's closed form is written for K >= 0; " *
            "got K = $K. For K < 0 use `_ising_sq_log_z_torus_signed`, which folds " *
            "onto |K| when the torus is bipartite and otherwise routes to the exact " *
            "transfer matrix (#824).",
        ),
    )
    # K = 0 (beta = 0 or J = 0) is a SINGULAR POINT OF THE PARAMETRISATION, not
    # of the physics.  Kaufman's variables all blow up there -- sinh 2K = 0 makes
    # the prefactor log(0), coth K diverges, and gamma_0 = 2K + log tanh K is
    # log(0) -- while the answer is the most trivial one there is: no bond
    # carries weight, so every one of the 2^{Lx Ly} configurations contributes 1.
    #
    # The transfer matrix handled this without special-casing (its entries just
    # become 1), so this is a regression the closed form introduces and has to
    # pay for explicitly.  Both `beta = 0` and `J = 0` are DOCUMENTED special
    # values of this fetch and are pinned in
    # `test/models/classical/test_ising_square_pfaffian.jl`.
    iszero(K) && return (Lx * Ly) * log(oftype(float(K), 2))

    γ_odd = [_ising_sq_gamma(π * (2r + 1) / Ly, K) for r in 0:(Ly - 1)]
    γ_even = [_ising_sq_gamma(π * (2r) / Ly, K) for r in 0:(Ly - 1)]
    h = Lx / 2

    # Each sector is a PRODUCT of Ly factors, so its log is a SUM of Ly logs —
    # which is the whole point: the product overflows at L ~ 30 (measured), the
    # sum does not.
    logP1 = sum(_log2cosh(h * γ) for γ in γ_odd)
    logP3 = sum(_log2cosh(h * γ) for γ in γ_even)
    s2 = 1
    logP2 = 0.0
    for γ in γ_odd
        (sg, lg) = _signed_log2sinh(h * γ)
        s2 *= sg
        logP2 += lg
    end
    s4 = 1
    logP4 = 0.0
    for γ in γ_even
        (sg, lg) = _signed_log2sinh(h * γ)
        s4 *= sg
        logP4 += lg
    end

    (ssum, logsum) = _signed_logsumexp(((1, logP1), (s2, logP2), (1, logP3), (s4, logP4)))
    ssum > 0 || throw(
        ErrorException(
            "IsingSquare log Z: the four Kaufman sectors summed to a non-positive " *
            "value (sign = $ssum) at Lx = $Lx, Ly = $Ly, K = $K — Z must be " *
            "positive, so this is a bug in the sector signs, not a numerical edge.",
        ),
    )
    return -log(2) + (Lx * Ly / 2) * log(2 * sinh(2K)) + logsum
end

"""
    _ising_sq_bipartite_torus(Lx::Integer, Ly::Integer) -> Bool

Whether the `Lx × Ly` PERIODIC square lattice is bipartite — i.e. both sides
even.  With either side odd the wrap-around closes an odd cycle, the two-colouring
fails, and the antiferromagnet is frustrated.

This is the whole content of #824: the infinite square lattice is bipartite, so
`Z` is even in `J` there and every reference formula is written for `|K|`; a
FINITE torus is bipartite only sometimes, and the difference is invisible on the
even×even cases anyone would try first.
"""
_ising_sq_bipartite_torus(Lx::Integer, Ly::Integer) = iseven(Lx) && iseven(Ly)

"""
    _ising_sq_log_z_torus_signed(Lx, Ly, K::Real) -> Real

`log Z` of the `Lx × Ly` periodic square-lattice Ising model at coupling `K = βJ`
of **either sign**.

Three routes, and which one applies is forced, not chosen:

- `K ≥ 0` — Kaufman's closed form, `O(Ly)`.
- `K < 0` on a BIPARTITE torus — `Z(−K) = Z(K)` exactly, by `σ → −σ` on one
  sublattice, so the closed form at `|K|` is not an approximation.
- `K < 0` on a NON-bipartite torus — the antiferromagnet is frustrated and `Z` is
  genuinely smaller. Kaufman's derivation assumes `K > 0` (its `γ_k` acquire an
  `iπ` shift and its `(2 sinh 2K)^{N/2}` prefactor a branch), so this routes to
  the EXACT transfer matrix instead of extending a formula outside its derivation.

The transfer matrix is `2^Ly × 2^Ly`, so the exponential axis is the SHORT side:
the two sides are exchanged when that helps, since `Z(Lx, Ly) = Z(Ly, Lx)` on a
torus. Beyond the cap it throws rather than silently returning the ferromagnetic
value, which is what #824 was.
"""
function _ising_sq_log_z_torus_signed(Lx::Integer, Ly::Integer, K::Real)
    (K >= 0 || _ising_sq_bipartite_torus(Lx, Ly)) &&
        return _ising_sq_log_z_torus(Lx, Ly, abs(K))
    # frustrated: put the exponential axis on the shorter side
    a, b = Lx <= Ly ? (Ly, Lx) : (Lx, Ly)     # a = power, b = transfer width
    b <= _MAX_ED_SITES || throw(
        ArgumentError(
            "IsingSquare log Z: K = $K < 0 on the non-bipartite $(Lx)×$(Ly) torus is " *
            "FRUSTRATED, so Kaufman's closed form (derived for K > 0) does not apply " *
            "and Z ≠ Z(|K|). The exact transfer-matrix route needs min(Lx, Ly) ≤ " *
            "$(_MAX_ED_SITES) and here it is $b. Use a smaller transverse size, or a " *
            "ferromagnetic K (#824).",
        ),
    )
    T = _ising_sq_transfer_matrix(b, one(K), K)      # β folded into K
    λ = eigvals(Symmetric(Matrix(T)))
    m = maximum(abs, λ)
    # signed log-sum-exp over λ^a: Z = tr(T^a) = Σ λ^a, and λ may be negative
    acc = sum(x -> sign(x)^a * (abs(x) / m)^a, λ)
    acc > 0 || throw(
        ErrorException(
            "IsingSquare log Z: the transfer-matrix eigenvalue sum came out " *
            "non-positive ($acc) at Lx = $Lx, Ly = $Ly, K = $K — Z must be positive.",
        ),
    )
    return a * log(m) + log(acc)
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
    L. Onsager, "Crystal Statistics. I.", [Onsager1944](@cite).
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
    model", [Yang1952](@cite).
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

# ═══════════════════════════════════════════════════════════════════════════════
# Critical exponents at T_c — delegate to 2D Ising Onsager universality
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingSquare, ::CriticalExponents, ::Infinite; kwargs...) -> NamedTuple

Onsager 1944 critical exponents of the 2D square-lattice Ising model
at T_c (= 2J / log(1 + √2)):

    α = 0,  β = 1/8,  γ = 7/4,  δ = 15,  ν = 1,  η = 1/4.

Delegated to the existing `Universality(:Ising)` infrastructure at
`d = 2`.  Rushbrooke (α + 2β + γ = 2), Widom (γ = β(δ − 1)), and
Fisher (η = 2 − γ/ν) hyperscaling relations all check out.

The returned NamedTuple also carries the central charge `c = 1//2`
inherited from the CFT minimal model M(3,4) — same payload as
`fetch(Universality(:Ising), CriticalExponents(); d=2)`.

# References

- L. Onsager, *Phys. Rev.* **65**, 117 (1944) — exact 2D Ising solution.
"""
function fetch(::IsingSquare, ::CriticalExponents, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.Universality(:Ising), CriticalExponents(); d=2)
end
