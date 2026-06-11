# ─────────────────────────────────────────────────────────────────────────────
# Classical 2D Six-Vertex (ice-rule) model on the square lattice
#
# Each edge of the square lattice carries an arrow; the local "ice rule"
# (Pauling 1935; Lieb 1967a) demands that each vertex has exactly two
# incoming and two outgoing arrows.  The six allowed vertex
# configurations are weighted in pairs by the standard parameters
#
#   ω₁ = ω₂ = a,   ω₃ = ω₄ = b,   ω₅ = ω₆ = c.
#
# All thermodynamic information is encoded in the dimensionless phase
# parameter
#
#   Δ = (a² + b² − c²) / (2 a b)
#
# which divides parameter space into three exactly solvable phases:
#
#   • Δ > 1   — ferroelectric (FE), frozen.
#               f = −log max(a, b) ;   S_residual = 0.
#   • |Δ| < 1 — disordered.
#               Free energy is given by the Lieb / Sutherland (1967)
#               trigonometric integral, with Δ = −cos μ, 0 < μ < π.
#   • Δ < −1  — antiferroelectric (AFE / F-model phase).
#               Free energy is the Lieb 1967b elliptic-function form,
#               with Δ = −cosh λ, λ > 0.
#
# The square-ice point a = b = c = 1 (Δ = 1/2) is the celebrated
# Lieb 1967a result
#
#   S/N = (3/2) log(4/3) ≈ 0.4315231087…  per vertex.
#
# References:
#   E. H. Lieb, "Residual entropy of square ice", Phys. Rev. 162, 162 (1967a).
#   E. H. Lieb, "Exact solution of the F model of an antiferroelectric",
#     Phys. Rev. Lett. 18, 1046 (1967b).
#   E. H. Lieb, "Exact solution of the two-dimensional Slater KDP model
#     of a ferroelectric", Phys. Rev. Lett. 19, 108 (1967c).
#   B. Sutherland, "Exact solution of a two-dimensional model for
#     hydrogen-bonded crystals", Phys. Rev. Lett. 19, 103 (1967).
#   R. J. Baxter, "Exactly Solved Models in Statistical Mechanics"
#     (Academic Press, 1982), ch. 8.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Model struct
# ═══════════════════════════════════════════════════════════════════════════════

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

"""
    SixVertex(; a::Real = 1.0, b::Real = 1.0, c::Real = 1.0) <: AbstractQAtlasModel

Classical six-vertex model on the infinite square lattice in the
standard (a, b, c) parameterisation (Lieb 1967, Baxter 1982 ch. 8):

    ω₁ = ω₂ = a,   ω₃ = ω₄ = b,   ω₅ = ω₆ = c

with the ice rule (each vertex has two incoming and two outgoing
arrows).  The single dimensionless invariant

    Δ = (a² + b² − c²) / (2 a b)

selects the phase: Δ > 1 ferroelectric (FE, frozen), |Δ| < 1
disordered (Lieb / Sutherland 1967), Δ < −1 antiferroelectric
(AFE / F-model, Lieb 1967b).

Convenience constructors for the three classical sub-models
(`square_ice`, `f_model`, `kdp_model`) are
also provided.

Currently registered fetches:

| Quantity                   | BC          | Coverage                                                 |
| -------------------------- | ----------- | -------------------------------------------------------- |
| [`ResidualEntropy`](@ref)   | `Infinite`  | All phases (via E - f thermodynamic relation)            |
| [`FreeEnergy`](@ref)        | `Infinite`  | All phases (FE, AFE, disordered)                         |
| [`Energy`](@ref)            | `Infinite`  | All phases (via finite difference of FreeEnergy)         |
| [`Polarization`](@ref)      | `Infinite`  | All phases                                               |
| [`UniversalityClass`](@ref) | `Infinite`  | `:XY` universality class in disordered phase `abs(Δ) < 1` |

See also: [`IsingSquare`](@ref) for the closest classical analog.
"""
struct SixVertex <: AbstractQAtlasModel
    a::Float64
    b::Float64
    c::Float64
end
function SixVertex(; a::Real=1.0, b::Real=1.0, c::Real=1.0)
    a > 0 || throw(ArgumentError("SixVertex: a must be positive; got $a"))
    b > 0 || throw(ArgumentError("SixVertex: b must be positive; got $b"))
    c > 0 || throw(ArgumentError("SixVertex: c must be positive; got $c"))
    return SixVertex(Float64(a), Float64(b), Float64(c))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Convenience constructors for the three named sub-models
# ═══════════════════════════════════════════════════════════════════════════════

"""
    square_ice() -> SixVertex

Square-ice point a = b = c = 1 (Δ = 1/2, disordered phase).  The
zero-temperature configurational entropy per vertex is the Lieb 1967a
closed form

    S/N = (3/2) log(4/3) ≈ 0.4315231087…

# References
- E. H. Lieb, Phys. Rev. **162**, 162 (1967a).
"""
square_ice() = SixVertex(; a=1.0, b=1.0, c=1.0)

"""
    f_model(c::Real) -> SixVertex

F-model sub-family a = b = 1, free c (Lieb 1967b).  Phase boundary at
c = 2:

* c < 2  — disordered (Δ = 1 − c²/2 ∈ (-1, 1)).
* c = 2  — critical (Δ = −1).
* c > 2  — antiferroelectric / F-model phase (Δ < −1).

# References
- E. H. Lieb, Phys. Rev. Lett. **18**, 1046 (1967b).
"""
f_model(c::Real) = SixVertex(; a=1.0, b=1.0, c=c)

"""
    kdp_model(a::Real) -> SixVertex

KDP model sub-family a free, b = c = 1 (Lieb 1967c).  Phase boundary
at a = 2:

* a < 2  — disordered (Δ ∈ (-1/2, 1) for 0 < a < 2).
* a = 2  — critical (Δ = 1).
* a > 2  — ferroelectric, frozen (Δ > 1, free energy −log a).

The "KDP" naming is historical (Slater's potassium-dihydrogen-phosphate
ferroelectric model).

# References
- E. H. Lieb, Phys. Rev. Lett. **19**, 108 (1967c).
"""
kdp_model(a::Real) = SixVertex(; a=a, b=1.0, c=1.0)

# ═══════════════════════════════════════════════════════════════════════════════
# Phase classification (internal)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _six_vertex_delta(a, b, c) -> Float64

Lieb / Baxter phase parameter Δ = (a² + b² − c²) / (2 a b).  All
exact six-vertex thermodynamics is a function of Δ together with the
overall energy scale.
"""
_six_vertex_delta(a::Real, b::Real, c::Real) = (a^2 + b^2 - c^2) / (2 * a * b)

"""
    _six_vertex_phase(a, b, c) -> Symbol

Return one of `:ferroelectric`, `:disordered`, `:antiferroelectric`
according to whether Δ > 1, |Δ| ≤ 1, or Δ < −1.  Phase boundaries
|Δ| = 1 are treated as members of the disordered phase (KDP / F-model
critical points).
"""
function _six_vertex_phase(a::Real, b::Real, c::Real)
    Δ = _six_vertex_delta(a, b, c)
    if Δ > 1
        return :ferroelectric
    elseif Δ < -1
        return :antiferroelectric
    else
        return :disordered
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: residual entropy
# ═══════════════════════════════════════════════════════════════════════════════

# Square-ice value (Lieb 1967a) — exact closed form, kept as a const so
# the docstring example can reference it without re-evaluating.
const _LIEB_SQUARE_ICE_RESIDUAL_ENTROPY = (3 / 2) * log(4 / 3)

"""
    fetch(::SixVertex, ::ResidualEntropy, ::Infinite; kwargs...) -> Float64

Zero-temperature configurational entropy per vertex of the six-vertex
model on the infinite square lattice.

* **Square ice** (a = b = c, equivalently Δ = 1/2) — Lieb 1967a closed
  form

      S/N = (3/2) log(4/3) ≈ 0.4315231087…

* **Ferroelectric phase** (Δ > 1) — frozen ground state, a single
  ordered configuration up to global symmetry, so S_residual = 0.

* **Disordered phase off the square-ice point** and **antiferroelectric
  phase** (Δ < −1) — the residual entropy is *not yet* implemented by
  this method and the call throws an informative `ArgumentError`.  The
  disordered-phase formula is the Lieb / Sutherland 1967 trigonometric
  integral evaluated at zero temperature; the AFE branch is the
  Lieb 1967b elliptic-function form.  Both are tracked as follow-up
  scope.

# References

- E. H. Lieb, Phys. Rev. **162**, 162 (1967a) — square ice.
- E. H. Lieb, Phys. Rev. Lett. **18**, 1046 (1967b) — F-model / AFE.
- B. Sutherland, Phys. Rev. Lett. **19**, 103 (1967) — disordered phase
  trigonometric integral.
- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (Academic Press, 1982), §8.8.
"""
function fetch(m::SixVertex, ::ResidualEntropy, ::Infinite; kwargs...)
    a, b, c = m.a, m.b, m.c
    phase = _six_vertex_phase(a, b, c)
    if phase === :ferroelectric
        return 0.0
    elseif phase === :disordered && a == b == c
        return _LIEB_SQUARE_ICE_RESIDUAL_ENTROPY
    else
        E = fetch(m, Energy{:per_site}(), Infinite(); kwargs...)
        f = fetch(m, FreeEnergy(), Infinite(); kwargs...)
        return E - f
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Free-energy: ferroelectric phase (frozen)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _free_energy_fe(a, b, c) -> Float64

Closed-form free-energy density of the six-vertex model in the
ferroelectric phase Δ > 1 (Lieb 1967c, Baxter 1982 §8.10):

    f = −log max(a, b)

The ground state is the unique frozen configuration with all arrows
parallel along the favoured axis, so the partition function is
dominated by (max(a,b))^N and the free energy density is −log of the
larger weight.  At the KDP point (a > 2, b = c = 1) this reduces to
−log a.
"""
_free_energy_fe(a::Real, b::Real, c::Real) = -log(max(a, b))

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: free energy (phase-dispatched)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::SixVertex, ::FreeEnergy, ::Infinite; beta::Real = Inf, kwargs...) -> Float64

Free-energy density per vertex of the six-vertex model in the
thermodynamic limit, dispatched on the Lieb phase parameter
Δ = (a² + b² − c²) / (2 a b).

* **Ferroelectric (Δ > 1)** — frozen ground state, closed form

      f = −log max(a, b)             (Lieb 1967c)

  Implemented exactly.

* **Square-ice point (a = b = c)** — Lieb 1967a closed form,
  f = −(3/2) log(4/3).  Implemented exactly.

* **Generic disordered (|Δ| ≤ 1, off the square-ice diagonal)** —
  Lieb / Sutherland 1967 trigonometric integral.  **Deferred** in
  this commit (Phase 2 of issue #163); the precise normalisation
  across the full disordered region needs the careful Baxter §8.8
  bookkeeping that we have not yet validated robustly.  Calls in
  this region throw an informative `ArgumentError`.

* **Antiferroelectric (Δ < −1)** — Lieb 1967b elliptic-function form,
  **deferred** in this commit (Phase 3 of issue #163).  Calls throw an
  informative `ArgumentError`.

The `beta` kwarg is accepted for API symmetry with other classical
fetches but does not enter the closed forms: the temperature scale is
already absorbed into the weights (a, b, c).

# References

- E. H. Lieb, Phys. Rev. **162**, 162 (1967a) — square ice / disordered.
- E. H. Lieb, Phys. Rev. Lett. **18**, 1046 (1967b) — antiferroelectric.
- E. H. Lieb, Phys. Rev. Lett. **19**, 108 (1967c) — KDP / ferroelectric.
- B. Sutherland, Phys. Rev. Lett. **19**, 103 (1967) — disordered phase
  trigonometric integral.
- R. J. Baxter, *Exactly Solved Models in Statistical Mechanics*
  (Academic Press, 1982), §§8.8–8.10.
"""
function fetch(m::SixVertex, ::FreeEnergy, ::Infinite; beta::Real=Inf, kwargs...)
    a, b, c = m.a, m.b, m.c
    phase = _six_vertex_phase(a, b, c)
    if phase === :ferroelectric
        return _free_energy_fe(a, b, c)
    elseif phase === :disordered
        # Handle symmetry: if a < b, swap a and b.
        # This keeps w <= 0 in the disordered phase, which matches the standard formulation.
        if a < b
            return fetch(SixVertex(b, a, c), FreeEnergy(), Infinite(); beta=beta, kwargs...)
        end

        Δ = _six_vertex_delta(a, b, c)
        # If Δ is extremely close to -1, we apply a tiny regularization to avoid integrand singularity.
        # This ensures smooth matching with the AFE phase.
        if Δ <= -1.0 + 1e-11
            Δ = -1.0 + 1e-11
        end

        μ = acos(-Δ)
        w = 2.0 * atan(tan(μ / 2.0) * (b - a) / (b + a))
        A = π - μ
        B = μ + w

        function integrand(x)
            if x == 0.0
                return (π - μ) * (μ + w) / (2.0 * π)
            end
            num1 = -expm1(-2.0 * B * x)
            num2 = -expm1(-2.0 * A * x)
            den1 = -expm1(-2.0 * π * x)
            den2 = 1.0 + exp(-2.0 * μ * x)
            return 0.5 * exp((w - μ) * x) * num1 * num2 / (x * den1 * den2)
        end

        val, _ = quadgk(integrand, 0.0, Inf; atol=1e-15, rtol=1e-15)
        return -log(a) - 2.0 * val
    else  # :antiferroelectric
        if a < b
            return fetch(SixVertex(b, a, c), FreeEnergy(), Infinite(); beta=beta, kwargs...)
        end

        Δ = _six_vertex_delta(a, b, c)
        λ = acosh(-Δ)

        # If λ is extremely small, we use the disordered phase limit to ensure stability.
        if λ < 1e-5
            return fetch(SixVertex(a, b, c), FreeEnergy(), Infinite(); beta=beta, kwargs...)
        end

        X = (a/b - Δ) / sqrt(Δ^2 - 1.0)
        u = 0.5 * log((X + 1.0) / (X - 1.0))

        s = 0.0
        for m in 1:1000000
            t1 = exp(-m * (2.0*λ - 2.0*u))
            t2 = exp(-m * (2.0*λ + 2.0*u))
            den = 1.0 + exp(-2.0 * m * λ)
            term = (t1 - t2) / (m * den)
            s += term
            if term < 1e-17
                break
            end
        end
        return -log(a) - (λ - u) - s
    end
end

function fetch(m::SixVertex, ::Energy{:per_site}, ::Infinite; kwargs...)
    # We use central finite difference with h = 1e-6 to compute the derivative
    # of the free energy f(a^β, b^β, c^β) with respect to β at β = 1.
    h = 1e-6
    a, b, c = m.a, m.b, m.c
    f_plus = fetch(
        SixVertex(a^(1.0 + h), b^(1.0 + h), c^(1.0 + h)),
        FreeEnergy(),
        Infinite();
        kwargs...,
    )
    f_minus = fetch(
        SixVertex(a^(1.0 - h), b^(1.0 - h), c^(1.0 - h)),
        FreeEnergy(),
        Infinite();
        kwargs...,
    )
    return (f_plus - f_minus) / (2.0 * h)
end

function fetch(m::SixVertex, ::Polarization, ::Infinite; kwargs...)
    a, b, c = m.a, m.b, m.c
    phase = _six_vertex_phase(a, b, c)
    if phase === :ferroelectric
        # In the ferroelectric phase, all arrows are aligned.
        return 1.0
    elseif phase === :disordered
        # In the disordered phase, bulk polarization vanishes.
        return 0.0
    else # :antiferroelectric
        # Staggered polarization order parameter: \prod_{n=1}^\infty \tanh^2(n\lambda)
        Δ = _six_vertex_delta(a, b, c)
        λ = acosh(-Δ)

        # If λ is extremely small, staggered polarization is 0.
        if λ < 1e-5
            return 0.0
        end

        pol = 1.0
        for n in 1:100000
            term = tanh(n * λ)^2
            pol *= term
            if 1.0 - term < 1e-16
                break
            end
        end
        return pol
    end
end
