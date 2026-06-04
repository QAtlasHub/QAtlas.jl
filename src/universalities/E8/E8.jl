"""
    E8_MASS_RATIOS
Analytical expressions for the 8 particle masses in the E8 spectrum of the Ising field theory. Normalized by the lightest mass m₁.

### References
- Theory: Zamolodchikov, A. B. (1989). Integrals of motion and S-matrix of the (scaled) T= Tc Ising model with magnetic field. *International Journal of Modern Physics A*, 4(16), 4235-4248.
- Review: Delfino, G. (2004). Integrable field theory and critical phenomena: the Ising model in a magnetic field. *Journal of Physics A: Mathematical and General*, 37(14), R45-R78.
- Experiment: Coldea, R., Tennant, D. A., Wheeler, E. M., Wawrzynska, E., Prabhakaran, D., Telling, M., ... & Kiefer, K. (2010). Quantum criticality in an Ising chain: experimental evidence for emergent E8 symmetry. *Science*, 327(5962), 177-180.
"""
function get_e8_mass_ratios()
    # Closed-form mass ratios of the E₈ particle spectrum of the
    # Zamolodchikov (1989) magnetic Ising field theory, normalised by
    # m₁ = 1.  The derivation (Perron-Frobenius eigenvector of the E₈
    # Cartan matrix, Dorey's fusion rules, exponents from the Coxeter
    # group) is written out in
    # `docs/src/calc/e8-mass-spectrum-derivation.md`.
    #
    # φ = 2 cos(π/5) is the golden ratio; m₂/m₁ = φ is the hallmark of
    # the E₈ symmetry and was measured experimentally in Coldea et al.
    # (2010) for CoNb₂O₆ on its Ising-chain plateau.
    ϕ = 2 * cos(π / 5)

    m1 = 1.0
    m2 = ϕ
    m3 = 2 * cos(π / 30)
    m4 = 2 * cos(7π / 30) * ϕ
    m5 = 2 * cos(2π / 15) * ϕ
    m6 = 2 * cos(π / 30) * ϕ
    m7 = 2 * cos(7π / 30) * ϕ^2
    m8 = 2 * cos(2π / 15) * ϕ^2

    return [m1, m2, m3, m4, m5, m6, m7, m8]
end

# --- Registry & Fetch ---

"""
    E8() <: AbstractQAtlasModel

The E8 integrable field theory — the low-energy continuum theory of the
2D Ising model perturbed by a small longitudinal magnetic field at
T = T_c (Zamolodchikov 1989).  The only physics parameter is the
underlying magnetic field perturbation strength, but QAtlas currently
only exposes the universal mass-ratio spectrum, so this struct carries
no fields.
"""
struct E8 <: AbstractQAtlasModel end

@register_aliases canonicalize_quantity :E8_spectrum [
    :E8_mass_ratios, :E8_masses, :mass_ratios, :E8_mass_ratio, :mass_ratio
]

"""
    fetch(::E8, ::E8Spectrum, ::Infinite) -> Vector{Float64}

Returns the analytical E8 mass spectrum `[m₁, m₂, …, m₈]` normalised
by `m₁ = 1`.  `m₂/m₁ = 2 cos(π/5) = ϕ` (golden ratio) is the hallmark
of the E8 symmetry.
"""
function QAtlas.fetch(::E8, ::E8Spectrum, ::Infinite; kwargs...)
    return get_e8_mass_ratios()
end
