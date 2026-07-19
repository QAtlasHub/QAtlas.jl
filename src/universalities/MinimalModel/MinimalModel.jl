# ─────────────────────────────────────────────────────────────────────────────
# Virasoro minimal models M(p, p_prime)
#
# References:
#   Belavin, Polyakov, Zamolodchikov (BPZ),
#     "Infinite conformal symmetry in two-dimensional quantum field
#      theory", [BelavinPolyakovZamolodchikov1984](@cite).
#   Friedan, Qiu, Shenker, [FriedanQiuShenker1984](@cite)
#     — minimal-model unitarity.
#   Di Francesco, Mathieu, Sénéchal,
#     "Conformal Field Theory" (Springer, 1997), Ch. 7.
#
# For coprime integers p > p_prime ≥ 2:
#
#   c(p, p_prime) = 1 - 6 (p - p_prime)^2 / (p p_prime)
#
#   h_{r,s}(p, p_prime) = ((p r - p_prime s)^2 - (p - p_prime)^2)
#                         / (4 p p_prime),
#       1 ≤ r ≤ p_prime - 1,  1 ≤ s ≤ p - 1
#
# Kac symmetry: h_{r,s} = h_{p_prime - r, p - s}.
# ─────────────────────────────────────────────────────────────────────────────

"""
    MinimalModel(p::Int, p_prime::Int) <: AbstractQAtlasModel

Virasoro minimal model M(p, p_prime).  `p` and `p_prime` must be
coprime integers with `p > p_prime ≥ 2`.  Construction validates
those conditions and throws `DomainError` otherwise.

Special cases (cross-check):

| Model                     | (p, p_prime) | c     |
|---------------------------|:------------:|:-----:|
| Yang–Lee (non-unitary)    | (5, 2)       | -22/5 |
| Ising                     | (4, 3)       |  1/2  |
| Tricritical Ising         | (5, 4)       |  7/10 |
| 3-state Potts (chiral)    | (6, 5)       |  4/5  |

The Ising special case `MinimalModel(4, 3)` reproduces the central
charge stored in `Universality(:Ising)`'s `CriticalExponents` table
(`c = 1//2`).

Use [`fetch`](@ref) with [`CentralCharge`](@ref) or
[`ConformalWeights`](@ref):

```julia
fetch(MinimalModel(4, 3), CentralCharge())                  # 1//2
fetch(MinimalModel(4, 3), ConformalWeights(); r=1, s=2)     # 1//16
```

See also: [`WZWSU2`](@ref), [`Universality`](@ref),
[`CentralCharge`](@ref), [`ConformalWeights`](@ref),
[`PrimaryFields`](@ref).
"""
struct MinimalModel <: AbstractQAtlasModel
    p::Int
    p_prime::Int
    function MinimalModel(p::Integer, p_prime::Integer)
        p_prime ≥ 2 || throw(
            DomainError(
                (p, p_prime),
                "MinimalModel(p, p_prime): require p_prime ≥ 2; got p_prime=$p_prime.",
            ),
        )
        p > p_prime || throw(
            DomainError(
                (p, p_prime),
                "MinimalModel(p, p_prime): require p > p_prime; got p=$p, p_prime=$p_prime.",
            ),
        )
        gcd(p, p_prime) == 1 || throw(
            DomainError(
                (p, p_prime),
                "MinimalModel(p, p_prime): p and p_prime must be coprime; gcd($p, $p_prime) = $(gcd(p, p_prime)).",
            ),
        )
        return new(Int(p), Int(p_prime))
    end
end

"""
    fetch(::MinimalModel, ::CentralCharge) -> Rational{Int}

Central charge of the Virasoro minimal model M(p, p_prime):

    c = 1 - 6 (p - p_prime)^2 / (p p_prime).

Returned as an exact `Rational{Int}`.
"""
function fetch(m::MinimalModel, ::CentralCharge; kwargs...)
    p, q = m.p, m.p_prime
    return 1 // 1 - (6 * (p - q)^2) // (p * q)
end

"""
    fetch(::MinimalModel, ::ConformalWeights; r::Int, s::Int) -> Rational{Int}

Kac-table conformal weight of the primary `(r, s)`:

    h_{r,s} = ((p r - p_prime s)^2 - (p - p_prime)^2) / (4 p p_prime),
        1 ≤ r ≤ p_prime - 1,  1 ≤ s ≤ p - 1.

Out-of-range `(r, s)` throws `DomainError`.  Use the Kac symmetry
`h_{r,s} = h_{p_prime - r, p - s}` to map a label outside the
fundamental rectangle into it explicitly.
"""
function fetch(m::MinimalModel, ::ConformalWeights; r::Integer, s::Integer, kwargs...)
    p, q = m.p, m.p_prime
    (1 ≤ r ≤ q - 1) || throw(
        DomainError(
            (r, s),
            "MinimalModel($p, $q) ConformalWeights: r must satisfy 1 ≤ r ≤ p_prime - 1 = $(q - 1); got r=$r.",
        ),
    )
    (1 ≤ s ≤ p - 1) || throw(
        DomainError(
            (r, s),
            "MinimalModel($p, $q) ConformalWeights: s must satisfy 1 ≤ s ≤ p - 1 = $(p - 1); got s=$s.",
        ),
    )
    num = (p * r - q * s)^2 - (p - q)^2
    den = 4 * p * q
    return num // den
end

"""
    fetch(::MinimalModel, ::PrimaryFields) -> Vector{NamedTuple}

All distinct primary fields of M(p, p_prime), modulo Kac symmetry
`(r, s) ~ (p_prime - r, p - s)`.  Each entry is a NamedTuple
`(r=Int, s=Int, h=Rational{Int})`.

The list is enumerated over `1 ≤ r ≤ p_prime - 1, 1 ≤ s ≤ p - 1`
and de-duplicated by selecting the lex-smallest `(r, s)` from each
Kac-symmetry orbit, so its length is

    (p - 1)(p_prime - 1) / 2.
"""
function fetch(m::MinimalModel, ::PrimaryFields; kwargs...)
    p, q = m.p, m.p_prime
    out = NamedTuple{(:r, :s, :h),Tuple{Int,Int,Rational{Int}}}[]
    seen = Set{Tuple{Int,Int}}()
    for r in 1:(q - 1), s in 1:(p - 1)
        # Kac symmetry orbit representative = lex-smallest of {(r,s), (q-r,p-s)}.
        r2, s2 = q - r, p - s
        rep = (r, s) ≤ (r2, s2) ? (r, s) : (r2, s2)
        if rep ∉ seen
            push!(seen, rep)
            h = fetch(m, ConformalWeights(); r=rep[1], s=rep[2])
            push!(out, (r=rep[1], s=rep[2], h=h))
        end
    end
    return out
end

# ─── Infinite-bc forwarding for verify() integration ───────────────────────
fetch(m::MinimalModel, q::CentralCharge, ::Infinite; kwargs...) = fetch(m, q; kwargs...)
fetch(m::MinimalModel, q::ConformalWeights, ::Infinite; kwargs...) = fetch(m, q; kwargs...)
fetch(m::MinimalModel, q::PrimaryFields, ::Infinite; kwargs...) = fetch(m, q; kwargs...)
