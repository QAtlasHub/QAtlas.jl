# ─────────────────────────────────────────────────────────────────────────────
# TTbar — universal irrelevant TT̄ deformation of a 2-D QFT.
#
# The TT̄ deformation (Zamolodchikov 2004; Smirnov-Zamolodchikov 2017;
# Cavaglià-Negro-Szécsényi-Tateo 2016) is a universal irrelevant
# composite operator deformation defined by
#
#     ∂L / ∂λ = det T_{μν}
#             = -(1/8) ( T_{ab} T^{ab} − (T^a_a)² ),
#
# with coupling `λ` of dimension (length)².  Because the deformation is
# irrelevant in the renormalisation-group sense, the UV data of the
# seed theory — in particular, the UV CFT central charge `c` — is
# preserved at all values of `λ`:
#
#     c(λ) = c.
#
# When the seed theory is a 2-D CFT with central charge `c`, the
# finite-volume vacuum energy on a circle of circumference `L` is
# (McGough-Mezei-Verlinde 2018, positive-λ branch):
#
#     E_0(L, λ) = (L / (2λ)) [ 1 − √( 1 − 2π c λ / (3 L²) ) ].
#
# Phase-1 entry registers only `CentralCharge` — the cleanest
# RG-invariant observable.  The deformed circle Casimir spectrum is
# tracked as Phase 2.
#
# References:
#   - A. B. Zamolodchikov, hep-th/0401146 (2004).
#   - F. A. Smirnov, A. B. Zamolodchikov, [SmirnovZamolodchikov2017](@cite).
#   - A. Cavaglià, S. Negro, I. M. Szécsényi, R. Tateo,
#     [CavagliaNegroSzecsenyiTateo2016](@cite).
#   - L. McGough, M. Mezei, H. Verlinde, [McGoughMezeiVerlinde2018](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

"""
    TTbar(; c::Real = 1.0, λ::Real = 0.0) <: AbstractQAtlasModel

Universal irrelevant TT̄ deformation of a 2-D QFT whose UV completion
is a CFT of central charge `c > 0`, with deformation coupling `λ` of
dimension (length)².  Default `c = 1` (free boson seed), `λ = 0`
(undeformed CFT).  Any real `λ` is admissible: `λ > 0` is the
Hagedorn-like branch, `λ < 0` is the "good-sign" branch.

Because the TT̄ deformation is irrelevant in the RG sense, the UV
central charge `c` is preserved at all `λ`:

    c(λ) = c.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                  |
| ------------------------------ | ---------- | ----------------------- |
| [`CentralCharge`](@ref)        | `Infinite` | analytic (`c(λ) = c`)   |

# References

- A. B. Zamolodchikov, *hep-th/0401146* (2004).
- F. A. Smirnov, A. B. Zamolodchikov, *Nucl. Phys. B* **915**, 363 (2017).
- A. Cavaglià, S. Negro, I. M. Szécsényi, R. Tateo,
  *JHEP* **10**, 112 (2016).
- L. McGough, M. Mezei, H. Verlinde, *JHEP* **04**, 010 (2018).
"""
struct TTbar <: AbstractQAtlasModel
    c::Float64       # UV CFT central charge (preserved under TT̄)
    λ::Float64       # TT̄ coupling, dimension (length)²
    function TTbar(c::Real, λ::Real)
        c > 0 ||
            throw(DomainError(c, "TTbar requires UV central charge c > 0; got c = $c."))
        # λ may be any real; positive λ is the "Hagedorn-like" branch,
        # negative λ is the "good-sign" branch.
        return new(Float64(c), Float64(λ))
    end
end
TTbar(; c::Real=1.0, λ::Real=0.0) = TTbar(c, λ)

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::TTbar, ::CentralCharge, ::Infinite; c=m.c, λ=m.λ) -> Float64

UV central charge of the TT̄-deformed theory.  The TT̄ deformation is
irrelevant in the RG sense, hence the UV central charge is preserved
for all values of the coupling `λ`:

    c(λ) = c.

`c ≤ 0` raises `DomainError`.  `λ` is unrestricted.

# References

- A. B. Zamolodchikov, *hep-th/0401146* (2004).
- F. A. Smirnov, A. B. Zamolodchikov, *Nucl. Phys. B* **915**, 363 (2017).
"""
function fetch(m::TTbar, ::CentralCharge, ::Infinite; c::Real=m.c, λ::Real=m.λ, kwargs...)
    c > 0 || throw(DomainError(c, "TTbar CentralCharge requires c > 0; got c = $c."))
    # TT̄ is irrelevant — central charge is preserved at all λ.
    return Float64(c)
end
