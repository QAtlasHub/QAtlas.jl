# ─────────────────────────────────────────────────────────────────────────────
# Poisson universality (integrable / quantum-chaos baseline).
#
# An uncorrelated spectrum (e.g. integrable many-body Hamiltonian, or a
# quantum localised phase) has Poisson level statistics:
#   * level-spacing distribution P(s) = exp(-s)
#   * mean ratio ⟨r⟩ = 2 log 2 - 1 ≈ 0.3863
#
# Reference: Y. Y. Atas, E. Bogomolny, O. Giraud, G. Roux,
#   Phys. Rev. Lett. 110, 084101 (2013).
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:Poisson}, ::WignerSurmise; s::Real) -> Float64

Level-spacing distribution of an uncorrelated (Poisson) spectrum,
`P(s) = exp(-s)`. Mean spacing is normalised to 1 by construction.

Although strictly speaking only Wigner-Dyson ensembles have a
"surmise", we register Poisson against the same `WignerSurmise`
quantity tag for symmetry of the API: callers can compare RMT vs
Poisson at the same call site.
"""
function fetch(::Universality{:Poisson}, ::WignerSurmise; s::Real, kwargs...)
    s ≥ 0 || throw(DomainError(s, "Poisson WignerSurmise: s must be ≥ 0"))
    return exp(-float(s))
end

"""
    fetch(::Universality{:Poisson}, ::MeanRatio) -> Float64

Mean of the consecutive level-spacing ratio
`r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1})`
for a Poisson (uncorrelated) spectrum:

    ⟨r⟩ = 2 log 2 - 1 ≈ 0.3862944

Atas-Bogomolny-Giraud-Roux, Phys. Rev. Lett. **110**, 084101 (2013),
Eq. (2).
"""
function fetch(::Universality{:Poisson}, ::MeanRatio; kwargs...)
    return 2 * log(2) - 1
end
