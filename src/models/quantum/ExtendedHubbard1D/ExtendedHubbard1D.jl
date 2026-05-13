# ─────────────────────────────────────────────────────────────────────────────
# ExtendedHubbard1D — t-U-V Hubbard chain (nearest-neighbour-density extension).
#
# Hamiltonian:
#
#   H = -t Σ_{i, σ} (c†_{i,σ} c_{i+1,σ} + h.c.)
#       + U Σ_i n_{i,↑} n_{i,↓}
#       + V Σ_i n_i n_{i+1}
#
#   (spin-½ fermions; t > 0 hopping, U ∈ ℝ on-site, V ∈ ℝ nearest-neighbour
#   density-density.)
#
# Phase 1 scope (this file):  V = 0 only.  At V = 0 the model collapses to
# the standard 1D Hubbard chain at half filling, for which the Lieb–Wu
# (1968) integral closed forms are already implemented in `Hubbard1D`.
# Phase 1 therefore delegates the V = 0 charge gap to `Hubbard1D` and
# raises `DomainError` for any V ≠ 0 — the V ≠ 0 phase diagram
# (CDW / SDW / phase separation, Voit 1995; Nakamura 2000) is deferred
# to Phase 2.
#
# References:
#   - E. H. Lieb, F. Y. Wu, PRL 20, 1445 (1968) — V=0 Bethe ansatz.
#   - J. Voit, Rep. Prog. Phys. 58, 977 (1995) — 1D bosonization review.
#   - M. Nakamura, Phys. Rev. B 61, 16377 (2000) — CDW/SDW/BOW phase
#     diagram of the U-V chain at half filling.
# ─────────────────────────────────────────────────────────────────────────────

"""
    ExtendedHubbard1D(; t::Real = 1.0, U::Real = 4.0, V::Real = 0.0)
        <: AbstractQAtlasModel

1D t-U-V Hubbard chain (nearest-neighbour-density extension of the
standard Hubbard model):

    H = -t Σ_{i, σ} (c†_{i,σ} c_{i+1,σ} + h.c.)
        + U Σ_i n_{i,↑} n_{i,↓}
        + V Σ_i n_i n_{i+1}.

Convention: `t > 0` hopping, `U` on-site, `V` nearest-neighbour
density-density.

# Phase 1 scope (this release)

Phase 1 exposes only the **V = 0 limit**, at which the model reduces
to the standard 1D Hubbard chain.  The `ChargeGap` at `Infinite()` is
delegated to [`Hubbard1D`](@ref) at half filling (μ = U/2),
i.e. the Lieb–Wu (1968) integral

    Δ_c = (16 t² / U) ∫_1^∞ dω  √(ω² - 1) / sinh(2π t ω / U).

Any `V ≠ 0` raises `DomainError`.  The V ≠ 0 phase diagram
(CDW / SDW / BOW / phase separation, Voit 1995; Nakamura 2000) is
deferred to Phase 2.

# References

- Lieb, Wu, *Phys. Rev. Lett.* **20**, 1445 (1968).
- Voit, *Rep. Prog. Phys.* **58**, 977 (1995).
- Nakamura, *Phys. Rev. B* **61**, 16377 (2000).
"""
struct ExtendedHubbard1D <: AbstractQAtlasModel
    t::Float64
    U::Float64
    V::Float64
    function ExtendedHubbard1D(t::Real, U::Real, V::Real)
        t > 0 || throw(DomainError(t, "ExtendedHubbard1D requires t > 0; got t = $t."))
        return new(Float64(t), Float64(U), Float64(V))
    end
end
function ExtendedHubbard1D(; t::Real=1.0, U::Real=4.0, V::Real=0.0)
    return ExtendedHubbard1D(t, U, V)
end

# ─── V = 0 guard ──────────────────────────────────────────────────────

"""
    _extended_hubbard1d_check_v_zero(model::ExtendedHubbard1D)

Throw `DomainError` if `model.V` is not zero.  Phase 1 only exposes
the V = 0 limit (delegated Lieb–Wu) — the V ≠ 0 phase diagram is
deferred to Phase 2.
"""
function _extended_hubbard1d_check_v_zero(model::ExtendedHubbard1D)
    if !iszero(model.V)
        throw(
            DomainError(
                model.V,
                "ExtendedHubbard1D ChargeGap: V ≠ 0 introduces nearest-neighbour " *
                "density-density interaction (CDW / SDW / BOW / phase separation, " *
                "Voit 1995; Nakamura 2000) and is deferred to Phase 2. " *
                "Got V = $(model.V).",
            ),
        )
    end
    return nothing
end

# ─── fetch methods ────────────────────────────────────────────────────

"""
    fetch(model::ExtendedHubbard1D, ::ChargeGap, ::Infinite) -> Float64

Charge (Mott) gap of the t-U-V Hubbard chain at half filling.  Phase 1
only supports the V = 0 limit, where the gap reduces to the Lieb–Wu
(1968) integral and is delegated to [`Hubbard1D`](@ref) at
`μ = U/2`.  Any `V ≠ 0` raises `DomainError`.

# References

- Lieb, Wu, *Phys. Rev. Lett.* **20**, 1445 (1968).
- Voit, *Rep. Prog. Phys.* **58**, 977 (1995).
- Nakamura, *Phys. Rev. B* **61**, 16377 (2000).
"""
function fetch(model::ExtendedHubbard1D, ::ChargeGap, ::Infinite; kwargs...)
    _extended_hubbard1d_check_v_zero(model)
    return fetch(
        Hubbard1D(; t=model.t, U=model.U, μ=(model.U / 2)), ChargeGap(), Infinite()
    )
end
