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
#   - E. Witten, [Witten1989](@cite).
#   - V. G. Knizhnik, A. B. Zamolodchikov, [KnizhnikZamolodchikov1984](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Stabilizer / operator product
#   Observable:  Operator-product expectations (Wilson loops, GSD, TEE, S-matrix entries); convention-free
#   Reference:   docs/src/conventions.md §Topological / operator-product

"""
    ChernSimons3D(; N::Integer = 2, k::Integer = 1) <: AbstractQAtlasModel

3-D SU(N)_k Chern-Simons TQFT (Witten 1989).  `N ≥ 2` is the gauge-
group rank-plus-one and `k ∈ ℤ_{>0}` is the (integer) Chern-Simons
level.

Phase 1 exposed the **boundary WZW central charge** via the Sugawara
construction.  Phase 2 adds the closed-form `S³` partition function
`Z(S³; SU(N)_k) = S_{0,0}` (Witten 1989 / Verlinde formula).  Wilson-loop
knot invariants and modular `S` / `T` matrices remain tracked for later phases.

Quantities registered:

| Quantity                       | BC         | Method                                            |
| ------------------------------ | ---------- | ------------------------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (Sugawara `k(N²-1)/(k+N)`)               |
| [`PartitionFunction`](@ref)    | `Infinite` | analytic (Witten `Z(S³)` = modular `S_{0,0}`)     |

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
        k ≥ 1 ||
            throw(DomainError(k, "ChernSimons3D requires integer level k ≥ 1; got k = $k."))
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

# ═══════════════════════════════════════════════════════════════════════════════
# Partition function on S³ — Witten 1989 / Verlinde S_{0,0} (Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ChernSimons3D, ::PartitionFunction, ::Infinite; N=m.N, k=m.k) -> Float64

Closed-form partition function of 3-D `SU(N)_k` Chern-Simons on the
three-sphere `S³` (Witten 1989), equal to the `S_{0,0}` entry of the
`SU(N)_k` modular S-matrix (Verlinde formula):

    Z(S³; SU(N)_k) = N^{-1/2} (k + N)^{-(N-1)/2}
                     ∏_{1 ≤ j < l ≤ N} 2 sin( π (l − j) / (k + N) ).

For `SU(2)_k` this simplifies to `Z = √(2 / (k + 2)) · sin(π / (k + 2))`.

# Boundary condition

`Infinite()` — `S³` is a closed compact 3-manifold without boundary;
no transfer-matrix BC label is meaningful, so the BC slot is the
catch-all `Infinite` tag also used for thermodynamic-limit quantities
elsewhere in QAtlas.

# Verified values

- `SU(2)_1`:  `Z = 1/√2 ≈ 0.7071067811865476`
- `SU(2)_2`:  `Z = 1/2 = 0.5`
- `SU(2)_3`:  `Z = √(2/5) · sin(π/5) ≈ 0.3717480344601845`
- `SU(3)_1`:  `Z = 1/√3 ≈ 0.5773502691896258`

# References

- E. Witten, *Comm. Math. Phys.* **121**, 351 (1989).
- E. P. Verlinde, *Nucl. Phys. B* **300**, 360 (1988).
"""
function fetch(
    m::ChernSimons3D,
    ::PartitionFunction,
    ::Infinite;
    N::Integer=m.N,
    k::Integer=m.k,
    kwargs...,
)
    N ≥ 2 || throw(
        DomainError(
            N, "ChernSimons3D PartitionFunction requires N ≥ 2 (SU(N)); got N = $N."
        ),
    )
    k ≥ 1 ||
        throw(DomainError(k, "ChernSimons3D PartitionFunction requires k ≥ 1; got k = $k."))
    p = k + N
    prefactor = float(N)^(-0.5) * float(p)^(-(N - 1) / 2)
    prod = 1.0
    for j in 1:(N - 1)
        for l in (j + 1):N
            prod *= 2 * sin(π * (l - j) / p)
        end
    end
    return prefactor * prod
end
