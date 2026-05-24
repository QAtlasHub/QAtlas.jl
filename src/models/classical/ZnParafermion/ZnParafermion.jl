# ─────────────────────────────────────────────────────────────────────────────
# ZnParafermion — Z_n parafermion conformal field theory.
#
# The Z_n parafermion CFT (Fateev-Zamolodchikov 1985) is the SU(2)_n / U(1)
# coset realisation of the Fateev-Zamolodchikov parafermions and arises as
# the continuum limit of the Z_n clock model at its self-dual critical
# point.  It is parametrised by the parafermion level
#
#     n ≥ 2  (integer)
#
# and has central charge
#
#     c(n) = 2 (n - 1) / (n + 2).
#
# Special values:
#
#   n = 2  →  c = 1/2   (Ising / Majorana, equivalent to MinimalModel(4, 3))
#   n = 3  →  c = 4/5   (3-state Potts critical CFT)
#   n = 4  →  c = 1     (compactified free boson at the Z_4 radius)
#   n = 5  →  c = 8/7
#   n → ∞  →  c → 2     (SU(2) WZW level → ∞ asymptote)
#
# Phase-1 entry registers only `CentralCharge`.  Primary scaling dimensions
# (Δ_{l,m} = l(l+2)/(4(n+2)) − m²/(4n) for the disorder/parafermion tower)
# and parafermion fusion-rule coefficients are tracked as Phase 2.
#
# References:
#   - V. A. Fateev, A. B. Zamolodchikov, "Nonlocal (parafermion) currents
#     in two-dimensional conformal quantum field theory and self-dual
#     critical points in Z_N-symmetric statistical systems",
#     Sov. Phys. JETP 62, 215 (1985).
# ─────────────────────────────────────────────────────────────────────────────

"""
    ZnParafermion(; n::Integer=3) <: AbstractQAtlasModel

Z_n parafermion CFT (Fateev-Zamolodchikov 1985), the SU(2)_n / U(1) coset
conformal field theory with central charge `c = 2(n-1)/(n+2)`.

The default `n = 3` selects the Z_3 parafermion = 3-state Potts critical
CFT with `c = 4/5`.

# Arguments

- `n::Integer = 3` — parafermion level; must satisfy `n ≥ 2`.

# Examples

```julia
QAtlas.fetch(ZnParafermion(; n=2), CentralCharge(), Infinite())  # 1//2  (Ising)
QAtlas.fetch(ZnParafermion(; n=3), CentralCharge(), Infinite())  # 4//5  (3-state Potts)
QAtlas.fetch(ZnParafermion(; n=4), CentralCharge(), Infinite())  # 1//1  (free boson)
QAtlas.fetch(ZnParafermion(; n=5), CentralCharge(), Infinite())  # 8//7
```

# References

- V. A. Fateev, A. B. Zamolodchikov, *Sov. Phys. JETP* **62**, 215 (1985).
"""
# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

struct ZnParafermion <: AbstractQAtlasModel
    n::Int
    function ZnParafermion(n::Integer)
        n ≥ 2 || throw(DomainError(n, "ZnParafermion requires level n ≥ 2; got n = $n."))
        return new(Int(n))
    end
end

ZnParafermion(; n::Integer=3) = ZnParafermion(n)

"""
    fetch(::ZnParafermion, ::CentralCharge, ::Infinite; n=m.n) -> Rational{Int}

Central charge of the Z_n parafermion CFT (Fateev-Zamolodchikov 1985):

    c(n) = 2(n - 1) / (n + 2)

equivalent to the SU(2)_n / U(1) coset central charge. Special values:

- n = 2  →  c = 1/2     (Ising / Majorana)
- n = 3  →  c = 4/5     (3-state Potts)
- n = 4  →  c = 1       (compactified free boson)
- n = 5  →  c = 8/7
- n → ∞  →  c → 2       (asymptotic SU(2) WZW level → ∞)

Returned as an exact `Rational{Int}`.

# References

- V. A. Fateev, A. B. Zamolodchikov, *Sov. Phys. JETP* **62**, 215 (1985).
"""
function fetch(m::ZnParafermion, ::CentralCharge, ::Infinite; n::Integer=m.n, kwargs...)
    n ≥ 2 ||
        throw(DomainError(n, "ZnParafermion CentralCharge requires n ≥ 2; got n = $n."))
    return Rational(2 * (n - 1), n + 2)
end
