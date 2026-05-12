# ─────────────────────────────────────────────────────────────────────────────
# ChernSimons3D — 3-D SU(N)_k Chern-Simons topological field theory.
#
# Action (with normalisation `k ∈ ℤ_{>0}`):
#
#     S_CS = (k / 4π) ∫_{M_3} tr( A ∧ dA + (2/3) A ∧ A ∧ A ),
#
# the unique gauge-invariant 3-manifold topological term for a
# connection valued in the simply-connected Lie group SU(N).  Witten
# (1989) showed the partition function reduces to a finite-dimensional
# expression over the Verlinde basis of the boundary chiral algebra
# `ŝu(N)_k`; in particular the boundary WZW CFT has Sugawara central
# charge
#
#     c(SU(N)_k) = k · dim(SU(N)) / (k + h_v(SU(N)))
#                = k (N² − 1) / (k + N),
#
# where `dim(SU(N)) = N² − 1` and the dual Coxeter number is
# `h_v(SU(N)) = N`.  Specialisations:
#
#     SU(2)_1  →  c = 1                  (free boson, = WZWSU2(1))
#     SU(2)_2  →  c = 3/2                (Ising × free boson)
#     SU(2)_k  →  c = 3k / (k + 2)       (= WZWSU2(k))
#     SU(3)_1  →  c = 2
#     SU(N)_1  →  c = N − 1              (= rank of SU(N))
#
# Phase-1 entry registers the boundary `CentralCharge` only.  Wilson
# loop / knot invariants (Jones / HOMFLY polynomials), modular S / T
# matrices and the S³ / T² × S¹ partition functions need dedicated
# quantity types (knot-indexed scalars; modular matrices indexed by
# integrable representations) and are tracked as Phase 2.
#
# References:
#   - E. Witten, Comm. Math. Phys. 121, 351 (1989).
#   - V. G. Knizhnik, A. B. Zamolodchikov, Nucl. Phys. B 247, 83 (1984).
# ─────────────────────────────────────────────────────────────────────────────

"""
    ChernSimons3D(; N::Integer = 2, k::Integer = 1) <: AbstractQAtlasModel

3-D SU(N)_k Chern-Simons TQFT (Witten 1989).  `N ≥ 2` is the gauge-
group rank-plus-one and `k ∈ ℤ_{>0}` is the (integer) Chern-Simons
level.

Phase 1 exposes the **boundary WZW central charge** via the Sugawara
construction.  Wilson-loop knot invariants and explicit 3-manifold
partition functions are tracked as Phase 2.

Quantities registered:

| Quantity                       | BC         | Method                                  |
| ------------------------------ | ---------- | --------------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (Sugawara `k(N²-1)/(k+N)`)     |

# References

- E. Witten, *Comm. Math. Phys.* **121**, 351 (1989).
- V. G. Knizhnik, A. B. Zamolodchikov, *Nucl. Phys. B* **247**, 83 (1984).
"""
struct ChernSimons3D <: AbstractQAtlasModel
    N::Int
    k::Int
    function ChernSimons3D(N::Integer, k::Integer)
        N ≥ 2 || throw(
            DomainError(N, "ChernSimons3D requires N ≥ 2 (SU(N) gauge group); got N = $N."),
        )
        k ≥ 1 || throw(
            DomainError(k, "ChernSimons3D requires integer level k ≥ 1; got k = $k."),
        )
        return new(Int(N), Int(k))
    end
end
ChernSimons3D(; N::Integer=2, k::Integer=1) = ChernSimons3D(N, k)

# ═══════════════════════════════════════════════════════════════════════════════
# Boundary WZW central charge — Sugawara construction
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ChernSimons3D, ::CentralCharge, ::Infinite; N=m.N, k=m.k)
        -> Rational{Int}

Sugawara central charge of the boundary WZW theory `ŝu(N)_k` dual to
3-D `SU(N)_k` Chern-Simons:

    c = k (N² − 1) / (k + N),

returned as an exact `Rational{Int}`.

# References

- E. Witten, *Comm. Math. Phys.* **121**, 351 (1989).
- V. G. Knizhnik, A. B. Zamolodchikov, *Nucl. Phys. B* **247**, 83 (1984).
"""
function fetch(
    model::ChernSimons3D,
    ::CentralCharge,
    ::Infinite;
    N::Integer=model.N,
    k::Integer=model.k,
    kwargs...,
)
    return Rational(k * (N^2 - 1), k + N)
end
