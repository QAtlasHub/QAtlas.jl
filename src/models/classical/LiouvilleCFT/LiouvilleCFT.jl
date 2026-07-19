# ─────────────────────────────────────────────────────────────────────────────
# LiouvilleCFT — non-compact Liouville conformal field theory.
#
# Liouville CFT (Polyakov 1981) is the 2-D non-compact conformal field
# theory governing 2-D quantum gravity coupled to matter and the
# semiclassical limit of non-critical bosonic strings.  It is
# parametrised by a single real coupling `b > 0` whose `Q ↔ Q` dual
# structure is encoded in the "background charge"
#
#     Q = b + 1/b,
#
# and yields the central charge
#
#     c = 1 + 6 Q² = 1 + 6 (b + 1/b)².
#
# The `b ↔ 1/b` self-duality leaves `Q` and `c` invariant — this is
# the (non-perturbative) Liouville S-duality that underlies the
# DOZZ exact three-point function (Dorn-Otto 1994;
# Zamolodchikov-Zamolodchikov 1996).  Vertex operator scaling
# dimensions are `Δ_α = α(Q − α)`.
#
# Phase 1 registered `CentralCharge`; Phase 2 adds the vertex-operator
# conformal weight `Δ_α = α(Q − α)`.  The DOZZ three-point function
# (Mellin-Barnes integral over (b, α_1, α_2, α_3)) requires a further
# quantity type and is still deferred.
#
# References:
#   - A. M. Polyakov, [Polyakov1981](@cite).
#   - H. Dorn, H.-J. Otto, [DornOtto1994](@cite).
#   - A. B. Zamolodchikov, A. B. Zamolodchikov,
#     [ZamolodchikovZamolodchikov1996](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

"""
    LiouvilleCFT(; b::Real = 1.0) <: AbstractQAtlasModel

Non-compact Liouville CFT (Polyakov 1981) at coupling `b > 0`.  The
canonical self-dual point `b = 1` (background charge `Q = 2`,
`c = 25`) is used as the default.  The `b ↔ 1/b` self-duality leaves
both `Q` and `c` invariant.

Quantities registered:

| Quantity                          | BC         | Method                           |
| --------------------------------- | ---------- | -------------------------------- |
| [`CentralCharge`](@ref)           | `Infinite` | analytic (`c = 1 + 6Q²`)         |
| [`ConformalWeights`](@ref)        | `Infinite` | analytic (`Δ_α = α(Q − α)`)      |

# References

- A. M. Polyakov, *Phys. Lett. B* **103**, 207 (1981).
- H. Dorn, H.-J. Otto, *Nucl. Phys. B* **429**, 375 (1994).
- A. B. Zamolodchikov, A. B. Zamolodchikov,
  *Nucl. Phys. B* **477**, 577 (1996).
"""
struct LiouvilleCFT <: AbstractQAtlasModel
    b::Float64
end
LiouvilleCFT(; b::Real=1.0) = LiouvilleCFT(Float64(b))

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::LiouvilleCFT, ::CentralCharge, ::Infinite; b=m.b) -> Float64

Central charge of Liouville CFT at coupling `b > 0`:

    c = 1 + 6 (b + 1/b)².

`b ↔ 1/b` self-duality leaves `c` invariant.  Special points:

- `b = 1`  →  `Q = 2`,    `c = 25`     (self-dual point)
- `b → 0⁺` →  `Q → ∞`,    `c → ∞`      (semi-classical / weak-coupling)
- `b → ∞`  →  `Q → ∞`,    `c → ∞`      (strong-coupling dual)

`b ≤ 0` raises `DomainError` (the Liouville coupling is real and
positive by convention; `b ↔ 1/b` lives entirely on the positive
half-line).

# References

- A. M. Polyakov, *Phys. Lett. B* **103**, 207 (1981).
"""
function fetch(m::LiouvilleCFT, ::CentralCharge, ::Infinite; b::Real=m.b, kwargs...)
    b > 0 || throw(DomainError(b, "LiouvilleCFT CentralCharge requires b > 0; got b = $b."))
    Q = b + 1 / b
    return 1 + 6 * Q^2
end

# ═══════════════════════════════════════════════════════════════════════════════
# Conformal weights (vertex operator dimensions)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::LiouvilleCFT, ::ConformalWeights, ::Infinite; α::Real, b=m.b) -> Float64

Conformal (scaling) dimension of the Liouville vertex operator
`V_α(z, z̄) = exp(2 α φ(z, z̄))` at coupling `b > 0` and continuous
real momentum `α`:

    Δ_α = α (Q − α),    Q = b + 1/b.

Both holomorphic and antiholomorphic weights coincide
(`h = h̄ = Δ_α`); the full scaling dimension of `V_α` is `2 Δ_α`.

`Δ_α` is invariant under the Liouville reflection `α ↔ Q − α`
(the vertex operators `V_α` and `V_{Q−α}` are identified up to
the DOZZ reflection coefficient).

Special values:

- `α = 0`     → `Δ = 0`              (identity operator)
- `α = Q/2`   → `Δ = Q²/4`           (boundary of the normalisable
                                       real-α slice; "Seiberg bound")
- `α = b`     → `Δ = b(Q − b) = 1`   (degenerate screening operator
                                       `V_b`; conformal spin 1)
- `α = 1/b`   → `Δ = (1/b)(Q − 1/b) = 1`
                                      (dual screening `V_{1/b}`)
- `α ↔ Q−α`   → `Δ_α = Δ_{Q−α}`     (reflection symmetry)

`b ≤ 0` raises `DomainError` (the Liouville coupling is real and
positive by convention).  `α` is unconstrained: it may be any real
number, including the non-normalisable region `α < 0` or `α > Q/2`,
where `Δ_α` still gives the formal scaling dimension.

# References

- A. M. Polyakov, *Phys. Lett. B* **103**, 207 (1981).
- A. B. Zamolodchikov, A. B. Zamolodchikov,
  *Nucl. Phys. B* **477**, 577 (1996).
"""
function fetch(
    m::LiouvilleCFT, ::ConformalWeights, ::Infinite; α::Real, b::Real=m.b, kwargs...
)
    b > 0 ||
        throw(DomainError(b, "LiouvilleCFT ConformalWeights requires b > 0; got b = $b."))
    Q = b + 1 / b
    return α * (Q - α)
end
