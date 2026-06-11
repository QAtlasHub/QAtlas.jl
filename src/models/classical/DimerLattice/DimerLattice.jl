# ─────────────────────────────────────────────────────────────────────────────
# DimerLattice — close-packed dimers on the square lattice (Kasteleyn 1961).
#
# A perfect matching ("dimer covering") tiles every site of an Lx × Ly grid with
# non-overlapping nearest-neighbour dominoes.  The unweighted partition function
# Z is the NUMBER of such coverings, given exactly by the Kasteleyn-Temperley-
# Fisher product (free / open boundaries):
#
#   Z(m, n) = ∏_{j=1}^{m} ∏_{k=1}^{n}
#               ( 4 cos²(jπ/(m+1)) + 4 cos²(kπ/(n+1)) )^{1/4}
#
# (and Z = 0 when m·n is odd — an odd number of sites cannot be perfectly
# matched).  The thermodynamic-limit entropy per site is the Catalan constant
# over π,
#
#   s = lim (ln Z)/(m n) = G/π ≈ 0.29156   (Fisher 1961),
#
# the canonical "residual entropy" of the close-packed dimer model.  The height
# representation of square-lattice dimers is a c = 1 compact free boson.
#
# References:
#   - P. W. Kasteleyn, Physica 27, 1209 (1961).
#   - H. N. V. Temperley, M. E. Fisher, Phil. Mag. 6, 1061 (1961).
#   - M. E. Fisher, Phys. Rev. 124, 1664 (1961).
# ─────────────────────────────────────────────────────────────────────────────

using Base.MathConstants: catalan

"""
    DimerLattice(; Lx = 0, Ly = 0) <: AbstractQAtlasModel

Close-packed (perfect-matching) dimer model on the open `Lx × Ly` square
lattice, solved exactly by the Kasteleyn-Temperley-Fisher Pfaffian method.  The
partition function with unit dimer weights is the number of perfect matchings.

`Lx`, `Ly` may be carried in the struct or passed as keyword arguments to
`fetch`.  The thermodynamic-limit entropy per site is `G/π` (Catalan / π),
exposed as [`ResidualEntropy`](@ref) at `Infinite`.

Currently registered fetches:

| Quantity                   | BC         | Coverage                                                              |
| -------------------------- | ---------- | --------------------------------------------------------------------- |
| [`PartitionFunction`](@ref)| `OBC`      | Perfect matchings count on finite Lx × Ly grid                        |
| [`ResidualEntropy`](@ref)   | `Infinite` | Catalan/π ≈ 0.29156 (Fisher 1961)                                     |
| [`FreeEnergy`](@ref)        | `Infinite` | Free energy density (equal to -ResidualEntropy)                       |
| [`UniversalityClass`](@ref) | `Infinite` | `:XY` universality class (c = 1 compact free boson)                    |
"""
struct DimerLattice <: AbstractQAtlasModel
    Lx::Int
    Ly::Int
end
function DimerLattice(; Lx::Integer=0, Ly::Integer=0)
    (Lx >= 0 && Ly >= 0) ||
        throw(ArgumentError("DimerLattice: Lx, Ly must be ≥ 0 (0 = unset); got $Lx × $Ly"))
    return DimerLattice(Int(Lx), Int(Ly))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: exact dimer count via the Kasteleyn-Temperley-Fisher product
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _dimer_count_square(m, n) -> Float64

Number of perfect matchings (dimer coverings) of the open `m × n` square grid,
by the Kasteleyn-Temperley-Fisher product formula.  Returns `0.0` when `m·n` is
odd.  The count is returned as a `Float64` and is exact only while it stays
below `2^53` (roughly up to `11×11`); for larger grids — where a `Float64`
cannot represent consecutive integers, and where the product overflows to `Inf`
near `L ≈ 50` — this throws rather than return a silently-rounded or `Inf`
"count".  Use `ResidualEntropy` / `FreeEnergy` for the thermodynamic density.
"""
function _dimer_count_square(m::Int, n::Int)
    (m >= 1 && n >= 1) || throw(ArgumentError("DimerLattice: need Lx, Ly ≥ 1; got $m × $n"))
    iseven(m * n) || return 0.0          # odd site count ⇒ no perfect matching
    logZ = 0.0
    @inbounds for j in 1:m, k in 1:n
        # term = 0 needs cos(jπ/(m+1)) = cos(kπ/(n+1)) = 0, i.e. j = (m+1)/2 AND
        # k = (n+1)/2 both integers ⇒ m, n both odd ⇒ m·n odd, excluded above; so
        # term > 0 here and log is safe.
        term = 4 * cos(j * π / (m + 1))^2 + 4 * cos(k * π / (n + 1))^2
        logZ += log(term)
    end
    # round() recovers the exact integer only below 2^53; beyond that Float64
    # cannot distinguish consecutive integers (and exp overflows to Inf near
    # L ≈ 50).  Fail loudly rather than return a silently-wrong / Inf count.
    logZ / 4 > 53 * log(2) && error(
        "DimerLattice: the $m × $n dimer count exceeds the exactly-representable " *
        "Float64 integer range (2^53); use ResidualEntropy / FreeEnergy(Infinite) " *
        "for the per-site density instead.",
    )
    return round(exp(logZ / 4))
end

# ═══════════════════════════════════════════════════════════════════════════════
# PartitionFunction — number of dimer coverings of the finite grid
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(m::DimerLattice, ::PartitionFunction; Lx, Ly) -> Float64

Number of perfect matchings of the open `Lx × Ly` square grid (the unit-weight
close-packed dimer partition function), via the Kasteleyn-Temperley-Fisher
product.  `0` when `Lx·Ly` is odd.  `Lx`/`Ly` default to the struct fields.
Exact while the count fits a `Float64` integer (`< 2^53`, roughly up to
`11×11`); larger grids throw — use [`ResidualEntropy`](@ref) /
[`FreeEnergy`](@ref) for the thermodynamic-limit density.
"""
function fetch(m::DimerLattice, ::PartitionFunction; Lx::Integer=m.Lx, Ly::Integer=m.Ly)
    (Lx > 0 && Ly > 0) || error(
        "DimerLattice PartitionFunction: Lx and Ly must be positive. " *
        "Pass them in the struct (DimerLattice(; Lx, Ly)) or as kwargs.",
    )
    return _dimer_count_square(Int(Lx), Int(Ly))
end

# BC-aware delegator: the KTF product is the OPEN-boundary count, so OBC is the
# natural BC; required so the (DimerLattice, PartitionFunction, OBC) triple
# resolves to a non-catch-all fetch (registry drift guard).
function fetch(
    m::DimerLattice, q::PartitionFunction, ::OBC; Lx::Integer=m.Lx, Ly::Integer=m.Ly
)
    return fetch(m, q; Lx=Lx, Ly=Ly)
end

# ═══════════════════════════════════════════════════════════════════════════════
# ResidualEntropy / FreeEnergy — thermodynamic limit (Infinite)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::DimerLattice, ::ResidualEntropy, ::Infinite) -> Float64

Entropy per site of the close-packed dimer model on the infinite square lattice
(Fisher 1961),

```math
s = \\lim_{N\\to\\infty} \\frac{\\ln Z}{N} = \\frac{G}{\\pi} \\approx 0.29156,
```

where `G` is Catalan's constant.  (Per dimer this is `2G/π`, since each dimer
covers two sites.)
"""
fetch(::DimerLattice, ::ResidualEntropy, ::Infinite; kwargs...) = catalan / π

"""
    fetch(::DimerLattice, ::FreeEnergy, ::Infinite) -> Float64

Free-energy density per site of the unit-weight close-packed dimer model on the
infinite square lattice.  The model is combinatorial (all coverings weight 1,
zero internal energy), so at unit temperature `f = −s = −G/π` (Catalan / π); see
[`ResidualEntropy`](@ref).
"""
fetch(::DimerLattice, ::FreeEnergy, ::Infinite; kwargs...) = -catalan / π
