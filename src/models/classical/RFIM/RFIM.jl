# ─────────────────────────────────────────────────────────────────────────────
# RFIM — Random-Field Ising Model.
#
# Hamiltonian:
#
#     H = -J Σ_{⟨i,j⟩} σ_i σ_j - Σ_i h_i σ_i,   σ_i ∈ {-1, +1},
#
# with quenched i.i.d. random fields `h_i` (typically Gaussian or
# bimodal ±h_0 with variance `Δ`).  The Imry-Ma (1975) heuristic
# argument compares the bulk energy gain of flipping a connected
# region (∝ L^{d-1}) against the random-field energy cost (∝ Δ L^{d/2}),
# showing that for `d ≤ d_l = 2` arbitrarily-large clusters can lower
# the energy at any nonzero disorder Δ > 0 — long-range
# ferromagnetic order is destroyed and
#
#     T_c(d ≤ 2, Δ > 0) = 0.
#
# For `d ≥ 3` the FM phase survives (Imbrie 1985, Bricmont-Kupiainen
# 1988 rigorous proofs); T_c(d) > 0 but no closed-form analytic value
# is known.  Parisi-Sourlas (1979) supersymmetric dimensional
# reduction predicts the critical behaviour of `d`-dimensional RFIM
# matches `(d − 2)`-dimensional pure Ising — only valid at high
# enough `d` (the SUSY is spontaneously broken below ~ d = 5;
# Tarjus-Tissier 2004 RG corrections).
#
# Phase-1 entry registers the Imry-Ma `CriticalTemperature = 0`
# result for `d ≤ 2`.  Critical exponents in `d ≥ 3` need a
# `CriticalExponents` delegation infrastructure indexed by `(d, Δ)`
# and are tracked as Phase 2.
#
# References:
#   - Y. Imry, S. Ma, Phys. Rev. Lett. 35, 1399 (1975).
#   - J. Imbrie, Comm. Math. Phys. 98, 145 (1985).
#   - J. Bricmont, A. Kupiainen, Phys. Rev. Lett. 59, 1829 (1987).
#   - G. Parisi, N. Sourlas, Phys. Rev. Lett. 43, 744 (1979).
#   - G. Tarjus, M. Tissier, Phys. Rev. Lett. 93, 267008 (2004).
# ─────────────────────────────────────────────────────────────────────────────

"""
    RFIM(; J::Real = 1.0, Δ::Real = 1.0) <: AbstractQAtlasModel

Random-Field Ising Model: `H = -J Σ σσ - Σ h_i σ_i` with quenched
i.i.d. random fields of variance `Δ²` (Gaussian or bimodal ±√Δ²).
The Imry-Ma (1975) argument gives the lower critical dimension
`d_l = 2`: for `d ≤ 2, Δ > 0` there is no long-range ferromagnetic
order at any positive temperature.

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                                          |
| ------------------------------ | ---------- | ----------------------------------------------- |
| [`CriticalTemperature`](@ref)  | `Infinite` | analytic (Imry-Ma, `d ≤ 2` ⇒ T_c = 0)           |

# References

- Y. Imry, S. Ma, *Phys. Rev. Lett.* **35**, 1399 (1975).
- J. Imbrie, *Comm. Math. Phys.* **98**, 145 (1985).
"""
struct RFIM <: AbstractQAtlasModel
    J::Float64
    Δ::Float64
end
RFIM(; J::Real=1.0, Δ::Real=1.0) = RFIM(Float64(J), Float64(Δ))

# ═══════════════════════════════════════════════════════════════════════════════
# Critical temperature
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::RFIM, ::CriticalTemperature, ::Infinite; d::Int, J=m.J, Δ=m.Δ)
        -> Float64

Critical temperature of the RFIM as a function of spatial dimension
`d`.  In the **Imry-Ma regime** `d ≤ 2` at non-zero disorder
`Δ > 0` long-range order is destroyed and `T_c = 0`.

For `d ≥ 3` (Imbrie 1985 / Bricmont-Kupiainen 1988 rigorous FM
phase) `T_c(d, J, Δ) > 0` but no closed-form analytic value is
known and a `DomainError` is raised — Phase 2 will plug in
numerical-reference values (e.g. Monte-Carlo at d = 3).  At
`Δ = 0` the model reduces to the pure Ising model and the call is
deferred to that model's `CriticalTemperature` entry; here we also
raise `DomainError`.

# References

- Y. Imry, S. Ma, *Phys. Rev. Lett.* **35**, 1399 (1975).
- J. Imbrie, *Comm. Math. Phys.* **98**, 145 (1985).
"""
function fetch(
    m::RFIM, ::CriticalTemperature, ::Infinite; d::Int, J::Real=m.J, Δ::Real=m.Δ, kwargs...
)
    d ≥ 1 || throw(DomainError(d, "RFIM CriticalTemperature requires d ≥ 1; got d = $d."))
    Δ > 0 || throw(
        DomainError(
            Δ,
            "RFIM CriticalTemperature: Δ = 0 reduces to pure Ising; use IsingSquare / IsingChain1D / appropriate Ising entry instead.",
        ),
    )
    if d ≤ 2
        return 0.0           # Imry-Ma 1975
    else
        throw(
            DomainError(
                d,
                "RFIM CriticalTemperature for d ≥ 3 has no closed-form analytic value; the FM phase is rigorously known to exist (Imbrie 1985 / Bricmont-Kupiainen 1988) but T_c(d, J, Δ) is a numerical reference value tracked as Phase 2.",
            ),
        )
    end
end
