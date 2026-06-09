# ─────────────────────────────────────────────────────────────────────────────
# Classical 2D Ising on the triangular lattice — finite-T thermodynamics (FM).
#
# Wannier 1950 sign convention `H = +J Σ σσ`, so the FERROMAGNET is `J < 0`.
# The ferromagnetic free energy per site is the Houtappel (1950) double integral
#
#   -β f(K) = log Z/N = ln 2
#             + (1/(8π²)) ∫₀^{2π}∫₀^{2π} dθ dφ
#                 ln[ cosh³(2K) + sinh³(2K) - sinh(2K)(cos θ + cos φ + cos(θ+φ)) ],
#       K = β|J|   (reduced FM coupling).
#
# Verified analytically against three independent conditions:
#   * high-T  K→0:  -βf → ln 2          (free spins; s → ln 2)
#   * low-T   K→∞:  -βf → 3K            (ε → 3J = -3|J|, three bonds per site)
#   * the integrand minimum (θ = φ = 0) vanishes at 2K_c = ln 3 / 2, i.e.
#     T_c = 4|J|/ln 3 — the registered Houtappel critical temperature.
#
# Energy / SpecificHeat / Entropy follow from `log Z/N` by central differences
# in β (mirrors IsingSquare_thermal.jl; src keeps ForwardDiff in [extras], so the
# shared `_cd1` / `_cd2` helpers from IsingSquare_thermal.jl are reused).  The
# AFM branch (`J > 0`) is frustrated (`T_c = 0`) with no Houtappel closed form
# and is not implemented — these methods require `J < 0`.
# ─────────────────────────────────────────────────────────────────────────────

using QuadGK: quadgk

native_energy_granularity(::IsingTriangular, ::Infinite) = :per_site

@inline function _itri_require_fm(J::Real)
    return J < 0 || throw(
        DomainError(
            J,
            "IsingTriangular finite-T thermodynamics require the ferromagnetic " *
            "branch J < 0 (Wannier convention H = +J Σ σσ); the AFM J > 0 " *
            "triangular lattice is frustrated with no Houtappel closed form.",
        ),
    )
end

"""
    _houtappel_log_z_per_site(K) -> Float64

Houtappel (1950) thermodynamic-limit `log Z / N = -βf` per site of the
ferromagnetic triangular-lattice Ising model, reduced coupling `K = β|J|`:

    -βf = ln 2 + (1/(8π²)) ∫₀^{2π}∫₀^{2π} dθ dφ
              ln[ cosh³(2K) + sinh³(2K) - sinh(2K)(cos θ + cos φ + cos(θ+φ)) ].

The bracket is non-negative for all `K ≥ 0`, vanishing only at the critical
point `2K_c = ln 3 / 2` (θ = φ = 0).  Nested adaptive Gauss–Kronrod
quadrature; the integrable log edge at criticality needs no special care.

Reference: R. M. F. Houtappel, *Physica* **16**, 425 (1950).
"""
function _houtappel_log_z_per_site(K::Real)
    c3 = cosh(2K)^3
    s = sinh(2K)
    s3 = s^3
    outer, _ = quadgk(0.0, 2π; rtol=1e-9) do θ
        cθ = cos(θ)
        inner, _ = quadgk(0.0, 2π; rtol=1e-9) do φ
            br = c3 + s3 - s * (cθ + cos(φ) + cos(θ + φ))
            return log(max(br, zero(K)))   # br ≥ 0 analytically; clamp K_c round-off
        end
        return inner
    end
    return log(2) + outer / (8π^2)
end

"""
    fetch(m::IsingTriangular, ::FreeEnergy, ::Infinite; beta, J=m.J) -> Float64

Per-site Helmholtz free energy `f(β) = -β⁻¹ log Z/N` of the ferromagnetic
(`J < 0`) triangular-lattice Ising model in the thermodynamic limit
(Houtappel 1950).  The frustrated AFM branch (`J > 0`) has no closed form
and raises `DomainError`.
"""
function fetch(m::IsingTriangular, ::FreeEnergy, ::Infinite; beta::Real, J::Real=m.J)
    _itri_require_fm(J)
    return -_houtappel_log_z_per_site(beta * abs(J)) / beta
end

"""
    fetch(m::IsingTriangular, ::Energy{:per_site}, ::Infinite; beta, J=m.J) -> Float64

Per-site thermal energy `ε(β) = -∂(log Z/N)/∂β` from the Houtappel closed
form via a central difference.  At low T (`β → ∞`) `ε → 3J = -3|J|` — the
ferromagnetic ground state aligns every spin, contributing `J` per bond and
three bonds per site.
"""
function fetch(m::IsingTriangular, ::Energy{:per_site}, ::Infinite; beta::Real, J::Real=m.J)
    _itri_require_fm(J)
    return -_cd1(b -> _houtappel_log_z_per_site(b * abs(J)), beta)
end

"""
    fetch(m::IsingTriangular, ::SpecificHeat, ::Infinite; beta, J=m.J) -> Float64

Per-site specific heat `c_v(β) = β² ∂²(log Z/N)/∂β²` via a second central
difference.  Singular at the critical point `T_c = 4|J|/ln 3` (2D-Ising
universality); callers stay off that slice.
"""
function fetch(m::IsingTriangular, ::SpecificHeat, ::Infinite; beta::Real, J::Real=m.J)
    _itri_require_fm(J)
    return beta^2 * _cd2(b -> _houtappel_log_z_per_site(b * abs(J)), beta)
end

"""
    fetch(m::IsingTriangular, ::ThermalEntropy, ::Infinite; beta, J=m.J) -> Float64

Per-site Gibbs entropy `s(β) = β(ε - f)` from the Houtappel free-energy and
energy paths.  Bounded between 0 (`T → 0`) and `ln 2` (`T → ∞`).
"""
function fetch(m::IsingTriangular, ::ThermalEntropy, ::Infinite; beta::Real, J::Real=m.J)
    ε = fetch(m, Energy(:per_site), Infinite(); beta=beta, J=J)
    f = fetch(m, FreeEnergy(), Infinite(); beta=beta, J=J)
    return beta * (ε - f)
end

"""
    fetch(m::IsingTriangular, ::SpontaneousMagnetization; β, J=m.J) -> Float64

Spontaneous magnetisation of the triangular-lattice Ising model.

For the ferromagnet (`J < 0`) the Potts–Domb closed form below `T_c`:

    M(T) = [1 - 16 x³ / ((1-x)³ (1+3x))]^{1/8},   x = e^{-4β|J|},   T < T_c,
    M(T) = 0,                                                        T ≥ T_c,

with critical exponent `β = 1/8`.  The bracket vanishes exactly at
`x = 1/3`, i.e. `T_c = 4|J|/ln 3` (the registered Houtappel value), and
`M → 1` as `T → 0`.

For the antiferromagnet (`J > 0`) the triangular lattice is frustrated with
no uniform long-range order at any temperature, so `M = 0`.

# References
- R. M. F. Houtappel, *Physica* **16**, 425 (1950).
- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics* (1982), Ch 11.
"""
function fetch(m::IsingTriangular, ::SpontaneousMagnetization; β::Real, J::Real=m.J)
    J < 0 || return 0.0    # AFM (J>0): frustrated, no spontaneous magnetisation
    x = exp(-4 * β * abs(J))
    frac = 16 * x^3 / ((1 - x)^3 * (1 + 3 * x))
    return frac ≥ 1 ? 0.0 : (1 - frac)^(1 / 8)   # frac ≥ 1 ⟺ x ≥ 1/3 ⟺ T ≥ T_c
end

function fetch(
    m::IsingTriangular, q::SpontaneousMagnetization, ::Infinite; β::Real, J::Real=m.J
)
    return fetch(m, q; β=β, J=J)
end
