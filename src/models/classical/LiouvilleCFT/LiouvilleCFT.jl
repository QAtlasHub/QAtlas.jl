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
# Phase-1 entry registers only `CentralCharge`.  The DOZZ three-point
# function and vertex-operator dimensions require new quantity types
# (Mellin-Barnes integrals returning `Float64` over a (b, α_1, α_2, α_3)
# parameter slice) and are tracked as Phase 2.
#
# References:
#   - A. M. Polyakov, Phys. Lett. B 103, 207 (1981).
#   - H. Dorn, H.-J. Otto, Nucl. Phys. B 429, 375 (1994).
#   - A. B. Zamolodchikov, A. B. Zamolodchikov,
#     Nucl. Phys. B 477, 577 (1996).
# ─────────────────────────────────────────────────────────────────────────────

"""
    LiouvilleCFT(; b::Real = 1.0) <: AbstractQAtlasModel

Non-compact Liouville CFT (Polyakov 1981) at coupling `b > 0`.  The
canonical self-dual point `b = 1` (background charge `Q = 2`,
`c = 25`) is used as the default.  The `b ↔ 1/b` self-duality leaves
both `Q` and `c` invariant.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                           |
| ------------------------------ | ---------- | -------------------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (`c = 1 + 6Q²`)         |

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
    b > 0 || throw(
        DomainError(b, "LiouvilleCFT CentralCharge requires b > 0; got b = $b."),
    )
    Q = b + 1 / b
    return 1 + 6 * Q^2
end
