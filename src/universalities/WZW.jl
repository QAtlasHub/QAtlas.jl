# ─────────────────────────────────────────────────────────────────────────────
# Wess–Zumino–Witten WZW SU(2)_k universality classes
#
# References:
#   Knizhnik, Zamolodchikov,
#     "Current algebra and Wess-Zumino model in two dimensions",
#     Nucl. Phys. B 247, 83 (1984).
#   Witten, Comm. Math. Phys. 92, 455 (1984).
#   Di Francesco, Mathieu, Sénéchal,
#     "Conformal Field Theory" (Springer, 1997), Ch. 15.
#   Affleck, "Quantum spin chains and the Haldane gap",
#     J. Phys.: Condens. Matter 1, 3047 (1989) — Heisenberg = SU(2)_1 WZW.
#
# Level k = 1, 2, 3, ...  Sugawara central charge:
#
#   c(k) = 3 k / (k + 2)
#
# Primary fields are labelled by SU(2) spin j ∈ {0, 1/2, 1, ..., k/2}
# with conformal weight
#
#   h_j = j (j + 1) / (k + 2).
# ─────────────────────────────────────────────────────────────────────────────

"""
    WZWSU2(k::Int) <: AbstractQAtlasModel

Wess–Zumino–Witten model with affine Lie algebra SU(2) at level
`k ≥ 1`.  Construction throws `DomainError` for `k ≤ 0`.

Special cases:

- `k = 1`: `c = 1` (free boson at the SU(2)-symmetric radius;
  low-energy theory of the spin-1/2 Heisenberg antiferromagnet,
  Affleck 1989).
- `k = 2`: `c = 3/2` — equivalent to **3 free Majorana fermions**
  (each contributing `c = 1/2`), or equivalently the smallest N=1
  super-Virasoro minimal model.  Note that "Ising × free Majorana"
  with one Majorana would only give `c = 1/2 + 1/2 = 1 ≠ 3/2`; the
  correct decomposition needs three Majorana fermions.
- `k = 3`: `c = 9/5`.

```julia
fetch(WZWSU2(1), CentralCharge())                    # 1//1
fetch(WZWSU2(1), ConformalWeights(); j=1//2)         # 1//4
```

See also: [`MinimalModel`](@ref), [`Universality`](@ref),
[`CentralCharge`](@ref), [`ConformalWeights`](@ref).
"""
struct WZWSU2 <: AbstractQAtlasModel
    k::Int
    function WZWSU2(k::Integer)
        k ≥ 1 || throw(DomainError(k, "WZWSU2(k): require k ≥ 1; got k=$k."))
        return new(Int(k))
    end
end

"""
    fetch(::WZWSU2, ::CentralCharge) -> Rational{Int}

Sugawara central charge `c = 3k / (k + 2)` of WZW SU(2) at level `k`.
Returned as an exact `Rational{Int}`.
"""
function fetch(w::WZWSU2, ::CentralCharge; kwargs...)
    k = w.k
    return (3 * k) // (k + 2)
end

"""
    fetch(::WZWSU2, ::ConformalWeights; j) -> Rational{Int}

Conformal weight `h_j = j (j+1) / (k+2)` of the spin-`j` primary at
level `k`.  `j` must be a non-negative half-integer (i.e.
`Rational` such that `2j ∈ ℤ_{≥0}`) with `0 ≤ j ≤ k/2`.

Out-of-range or non-half-integer `j` throws `DomainError`.

```julia
fetch(WZWSU2(1), ConformalWeights(); j=0)            # 0//1
fetch(WZWSU2(1), ConformalWeights(); j=1//2)         # 1//4
fetch(WZWSU2(2), ConformalWeights(); j=1)            # 1//2
```
"""
function fetch(w::WZWSU2, ::ConformalWeights; j, kwargs...)
    k = w.k
    jr = _wzw_su2_normalize_spin(j, k)
    # h_j = j (j+1) / (k+2); compute as Rational{Int} exactly.
    return (jr * (jr + 1)) // (k + 2)
end

# Coerce the user-supplied spin into a Rational{Int} after validating
# that it is a non-negative half-integer with 2j ∈ {0, 1, ..., k}.
function _wzw_su2_normalize_spin(j::Rational, k::Int)
    j ≥ 0 || throw(DomainError(j, "WZWSU2 ConformalWeights: j must be ≥ 0; got j=$j."))
    twoj = 2 * j
    isinteger(twoj) || throw(
        DomainError(
            j,
            "WZWSU2 ConformalWeights: j must be a half-integer (2j ∈ ℤ); got j=$j (2j=$twoj).",
        ),
    )
    twoj_int = Int(twoj)
    twoj_int ≤ k || throw(
        DomainError(
            j,
            "WZWSU2(k=$k) ConformalWeights: j must satisfy 0 ≤ j ≤ k/2 = $(k//2); got j=$j.",
        ),
    )
    return Rational{Int}(j)
end
function _wzw_su2_normalize_spin(j::Integer, k::Int)
    return _wzw_su2_normalize_spin(Rational{Int}(j), k)
end
function _wzw_su2_normalize_spin(j::Real, k::Int)
    return throw(
        DomainError(
            j,
            "WZWSU2 ConformalWeights: j must be an Integer or Rational half-integer; got $(typeof(j)) value $j.  Pass e.g. `j = 1//2`.",
        ),
    )
end

# ─── Infinite-bc forwarding for verify() integration ───────────────────────
fetch(w::WZWSU2, q::CentralCharge, ::Infinite; kwargs...) = fetch(w, q; kwargs...)
fetch(w::WZWSU2, q::ConformalWeights, ::Infinite; kwargs...) = fetch(w, q; kwargs...)
