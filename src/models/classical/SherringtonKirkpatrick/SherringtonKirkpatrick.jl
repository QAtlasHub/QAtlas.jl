# ─────────────────────────────────────────────────────────────────────────────
# SherringtonKirkpatrick — mean-field Ising spin glass (Sherrington-Kirkpatrick 1975).
#
# Hamiltonian (1/√N normalisation that makes the energy extensive in
# the thermodynamic limit):
#
#     H = -(1/√N) Σ_{i<j} J_ij σ_i σ_j,   σ_i ∈ {-1, +1},
#
# where the couplings `J_ij` are i.i.d. Gaussian with mean 0 and
# variance `J²` (units of `J = 1` are the textbook convention).
#
# Equivalent normalisations: writing the same Hamiltonian as
#     H = -sum_{i<j} J_hat_ij sigma_i sigma_j,   J_hat_ij ~ N(0, J^2 / N)
# (variance scales as 1/N, no explicit 1/sqrt(N) prefactor) gives the
# same physics: the cavity-field variance at site i remains O(J^2) and
# the spin-glass transition still sits at T_c = J.  The two conventions
# differ only in whether the 1/sqrt(N) is absorbed into J_ij or kept in
# the Hamiltonian; the textbook (SK 1975) uses the explicit prefactor
# we follow here.
#
# At zero field the model has a continuous spin-glass transition at
#
#     T_c = J            (Sherrington-Kirkpatrick 1975)
#
# below which the Edwards-Anderson order parameter `q_EA(T)` and the
# overlap distribution `P(q)` acquire a continuous (full-RSB) structure
# (Parisi 1980; rigorously proven by Talagrand 2006).  Above T_c the
# replica-symmetric paramagnetic solution holds with `q = 0`.
#
# Phase 1 registers `CriticalTemperature` (T_c = J).  Phase 2 adds the
# Parisi T=0 ground-state energy density `e_0/J ≈ -0.7631667`
# (full-RSB; Parisi 1980, Crisanti-Rizzo 2002, Talagrand 2006) as the
# `Energy{:per_site}` at `Infinite` entry.  The finite-temperature
# Parisi free energy `f(β, h)`, the de Almeida-Thouless line
# `h_AT(T)`, and the overlap distribution `P(q)` require dedicated
# quantity types (function of either β or β + h, or a probability
# measure on `[0,1]`) and remain tracked for later phases.
#
# References:
#   - D. Sherrington, S. Kirkpatrick, Phys. Rev. Lett. 35, 1792 (1975).
#   - G. Parisi, J. Phys. A 13, L115 (1980).
#   - M. Talagrand, Annals Math. 163, 221 (2006).
# ─────────────────────────────────────────────────────────────────────────────

"""
    SherringtonKirkpatrick(; J::Real = 1.0) <: AbstractQAtlasModel

Sherrington-Kirkpatrick mean-field Ising spin glass (Sherrington-
Kirkpatrick 1975).  The Hamiltonian is the canonical 1/√N
random-Gaussian-coupling sum

    H = -(1/√N) Σ_{i<j} J_ij σ_i σ_j,   J_ij ~ N(0, J²),

whose spin-glass transition lies at `T_c = J`.

Quantities registered (Phases 1 and 2):

| Quantity                       | BC         | Method                  |
| ------------------------------ | ---------- | ----------------------- |
| [`CriticalTemperature`](@ref)  | `Infinite` | analytic                |
| [`Energy`](@ref) (`:per_site`) | `Infinite` | variational reference   |

The Phase-2 `Energy{:per_site}` entry is the Parisi T=0 ground-state
energy density `e_0/J ≈ -0.7631667` (Crisanti-Rizzo 2002
high-precision evaluation of Parisi's full-RSB solution; Talagrand
2006 rigorous proof).  The finite-temperature Parisi free energy
`f(β, h)`, the de Almeida-Thouless line `h_AT(T)` and the overlap
distribution `P(q)` remain tracked for later phases.

# References

- D. Sherrington, S. Kirkpatrick, *Phys. Rev. Lett.* **35**, 1792 (1975).
- G. Parisi, *J. Phys. A* **13**, L115 (1980).
- M. Talagrand, *Annals Math.* **163**, 221 (2006).
"""
struct SherringtonKirkpatrick <: AbstractQAtlasModel
    J::Float64
end
SherringtonKirkpatrick(; J::Real=1.0) = SherringtonKirkpatrick(Float64(J))

# ═══════════════════════════════════════════════════════════════════════════════
# Critical temperature — Sherrington-Kirkpatrick 1975
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SherringtonKirkpatrick, ::CriticalTemperature, ::Infinite; J=m.J)
        -> Float64

Spin-glass critical temperature `T_c = J` of the Sherrington-
Kirkpatrick model in the 1/√N normalisation with `J_ij ~ N(0, J²)`
(Sherrington-Kirkpatrick 1975).  For `J ≤ 0` the model has no
non-trivial ordering temperature (degenerate or unphysical
distribution) and `T_c = 0` is returned.

# References

- D. Sherrington, S. Kirkpatrick, *Phys. Rev. Lett.* **35**, 1792 (1975).
"""
function fetch(
    m::SherringtonKirkpatrick, ::CriticalTemperature, ::Infinite; J::Real=m.J, kwargs...
)
    return J > 0 ? Float64(J) : 0.0
end

# Native energy granularity at Infinite is per-site for the mean-field SK
# spin glass (the Hamiltonian is extensive but the natural reference is
# the per-spin energy density e_0/J in the thermodynamic limit).
QAtlas.native_energy_granularity(::SherringtonKirkpatrick, ::Infinite) = :per_site

# Hardcoded reference value for the SK T=0 Parisi full-RSB ground-state
# energy density.  References:
#   - G. Parisi, J. Phys. A 13, L115 (1980): full-RSB variational solution.
#   - A. Crisanti, T. Rizzo, Phys. Rev. E 65, 046137 (2002): high-precision
#     numerical evaluation gives e_0/J = -0.7631667 +/- 0.00001.
#   - M. Talagrand, Annals Math. 163, 221 (2006): rigorous proof of the
#     Parisi formula.
# Update if a tighter literature consensus emerges (e.g. recent
# Schmidt / Crisanti-Rizzo evaluations).
const _SK_PARISI_T0_ENERGY_DENSITY_PER_J = -0.7631667

# ═══════════════════════════════════════════════════════════════════════════════
# Parisi T=0 ground-state energy density (full RSB, Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SherringtonKirkpatrick, ::Energy{:per_site}, ::Infinite; J=m.J)
        -> Float64

Zero-temperature ground-state energy density of the mean-field SK
spin glass under Parisi's full-replica-symmetry-breaking solution
(Parisi 1980; rigorously proved by Talagrand 2006):

    e_0 / J  ≈  -0.7631667 ± 0.00001    (Crisanti-Rizzo 2002).

Returned as `J * (-0.7631667)`.  Reliability is `:high` since the
numerical value is mathematically established and known to 7 digits.

# References

- G. Parisi, *J. Phys. A* **13**, L115 (1980).
- A. Crisanti, T. Rizzo, *Phys. Rev. E* **65**, 046137 (2002).
- M. Talagrand, *Annals Math.* **163**, 221 (2006).
"""
function fetch(
    m::SherringtonKirkpatrick, ::Energy{:per_site}, ::Infinite;
    J::Real=m.J, kwargs...,
)
    J > 0 || throw(
        DomainError(
            J, "SherringtonKirkpatrick Energy requires J > 0; got J = $J."
        ),
    )
    return J * _SK_PARISI_T0_ENERGY_DENSITY_PER_J
end
