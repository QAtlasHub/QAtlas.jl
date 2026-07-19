# ─────────────────────────────────────────────────────────────────────────────
# Classical 2D Ising model on the triangular lattice — Wannier 1950
#
# Hamiltonian (sign convention from Wannier 1950):
#
#     H = +J Σ_{⟨i,j⟩} σ_i σ_j,   σ_i ∈ {-1, +1}
#
# with J > 0 the antiferromagnetic (frustrated) coupling and J < 0 the
# ferromagnetic (Houtappel 1950) coupling.  Each lattice site has six
# nearest neighbours arranged on a regular triangular net.
#
# Closed-form results carried by this file:
#
#   * AFM (J > 0):
#       T_c = 0   (frustration ⇒ no long-range order at any T > 0)
#       Zero-temperature residual entropy (Wannier 1950, exact):
#           S/(N k_B) = (2/π) ∫₀^{π/3} ln(2 cos θ) dθ ≈ 0.32306594722…
#
#   * FM (J < 0):
#       T_c = 4 |J| / ln 3   (Houtappel 1950, analytical critical temperature)
#       S_residual = 0       (ordered ferromagnet, unique ground states up to ℤ₂)
#
# Two-point spin correlations on the triangular lattice are known
# analytically (Stephenson 1964, J. Math. Phys. 5, 1009) but are
# deliberately *not* implemented here — they are tracked as a separate
# follow-up.
#
# References:
#   G. H. Wannier, "Antiferromagnetism. The triangular Ising net",
#     [Wannier1950](@cite).
#   R. M. F. Houtappel, "Order-disorder in hexagonal lattices",
#     [Houtappel1950](@cite).
#   J. Stephenson, "Ising-model spin correlations on the triangular
#     lattice", [Stephenson1964](@cite).  [Future work — two-point
#     correlations, not implemented in this commit.]
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: see file-header description above
#   Observable:  per src/core/quantities.jl (matches the dispatch tag)
#   Reference:   docs/src/conventions.md (project-wide convention policy)
#   STATUS:      backfilled by PR (audit gate); per-field domain content
#                left to a follow-up - see issue tracker for the model-specific
#                Hamiltonian sign / observable normalisation.

using QuadGK: quadgk

# ═══════════════════════════════════════════════════════════════════════════════
# Model struct
# ═══════════════════════════════════════════════════════════════════════════════

"""
    IsingTriangular(; J::Real = 1.0) <: AbstractQAtlasModel

Classical 2D Ising model on the triangular lattice with the
Wannier 1950 sign convention

    H = +J Σ_{⟨i,j⟩} σ_i σ_j,   σ_i ∈ {-1, +1}.

The coupling sign determines the physics:

- `J > 0` — antiferromagnetic, frustrated.  The triangular plaquette
  cannot satisfy all three antiferromagnetic bonds simultaneously, so
  the ground-state manifold is macroscopically degenerate and there is
  no long-range order at any ``T > 0`` (Wannier 1950).
- `J < 0` — ferromagnetic.  The lattice supports a standard
  order-disorder transition at ``T_c = 4|J| / \\ln 3`` (Houtappel 1950).
- `J = 0` — non-interacting; degenerate value preserved for symmetry.

Quantities currently registered for this model:

| Quantity                                  | BC          | Method            |
| ----------------------------------------- | ----------- | ----------------- |
| [`CriticalTemperature`](@ref)             | `Infinite`  | analytic          |
| [`ResidualEntropy`](@ref)                 | `Infinite`  | analytic (QuadGK) |

!!! warning "Sign convention differs from `IsingSquare`"
    `IsingTriangular` uses the **Wannier 1950 convention** `H = +J Σ σσ`,
    so `J > 0` means *antiferromagnetic* (frustrated).  By contrast,
    [`IsingSquare`](@ref) uses the modern `H = -J Σ σσ` convention
    where `J > 0` means *ferromagnetic*.  Passing the same numerical
    value of `J` to the two models therefore selects opposite physics.
    To keep the user-facing critical-temperature comparison meaningful:
    `IsingSquare(J=1.0)` (FM) ↔ `IsingTriangular(J=-1.0)` (FM).

See also: [`IsingSquare`](@ref) for the square-lattice (non-frustrated)
counterpart with the Onsager / Yang closed forms.
"""
struct IsingTriangular <: AbstractQAtlasModel
    J::Float64
end
IsingTriangular(; J::Real=1.0) = IsingTriangular(Float64(J))

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: critical temperature
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingTriangular, ::CriticalTemperature, ::Infinite; J=m.J) -> Float64

Exact critical temperature of the classical 2D Ising model on the
triangular lattice in the Wannier 1950 sign convention
``H = +J Σ σ_i σ_j``:

* `J > 0` (AFM, frustrated) — `T_c = 0`.  No long-range order at any
  positive temperature (Wannier 1950).
* `J < 0` (FM, Houtappel) — `T_c = 4 |J| / ln 3 ≈ 3.6409 |J|`
  (Houtappel 1950).
* `J = 0` — `T_c = 0` (no interaction; degenerate value, kept finite).

# References

- G. H. Wannier, Phys. Rev. **79**, 357 (1950).
- R. M. F. Houtappel, Physica **16**, 425 (1950).
"""
function fetch(m::IsingTriangular, ::CriticalTemperature, ::Infinite; J::Real=m.J)
    if J > 0
        return 0.0          # frustrated AFM, Wannier 1950
    elseif J < 0
        return 4 * abs(J) / log(3)   # Houtappel 1950
    else
        return 0.0          # J = 0, no interaction
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: zero-temperature residual entropy
# ═══════════════════════════════════════════════════════════════════════════════

# Pre-computed Wannier integrand and reference value, kept private so
# the docstring example can reference the constant without re-evaluating
# the integral.
const _WANNIER_RESIDUAL_ENTROPY = let
    val, _err = quadgk(θ -> log(2 * cos(θ)), 0.0, π / 3; rtol=1e-14, atol=1e-14)
    (2 / π) * val
end

"""
    fetch(::IsingTriangular, ::ResidualEntropy, ::Infinite; J=m.J) -> Float64

Zero-temperature residual entropy per site of the classical Ising model
on the triangular lattice in the Wannier 1950 sign convention
``H = +J Σ σ_i σ_j``.

* `J > 0` (frustrated AFM) — Wannier (1950) closed form

      S/(N k_B) = (2/π) ∫₀^{π/3} ln(2 cos θ) dθ ≈ 0.32306594722…

  The integral is evaluated by `QuadGK.quadgk` to ~1e-12 precision.

* `J ≤ 0` (FM or non-interacting Ising on a triangular lattice) — there
  is a unique pair of degenerate FM ground states related by the
  global ℤ₂ flip, so `S_residual = 0`.

# References

- G. H. Wannier, "Antiferromagnetism. The triangular Ising net",
  Phys. Rev. **79**, 357 (1950).
- R. M. F. Houtappel, Physica **16**, 425 (1950) — independent
  derivation in the same period.
"""
function fetch(m::IsingTriangular, ::ResidualEntropy, ::Infinite; J::Real=m.J)
    if J > 0
        return _WANNIER_RESIDUAL_ENTROPY
    else
        return 0.0
    end
end

"""
    fetch(::IsingTriangular, ::ZZCorrelation{:static}, ::Infinite; r=1, J=m.J) -> Float64

Zero-temperature nearest-neighbour spin–spin correlation `⟨σ_i σ_j⟩` of the
classical triangular Ising model (Wannier 1950 convention `H = +J Σ σσ`).

For the frustrated antiferromagnet (`J > 0`), `r = 1`, the exact Wannier
(1950) value

    ⟨σ_i σ_{i+1}⟩_{T=0} = -1/3,

a direct consequence of the ground-state rule that every elementary triangle
carries exactly one unsatisfied bond (`Σ σσ = -1` per triangle ⇒ `-1/3` per
bond).  General separations `r > 1` (Stephenson 1964 closed form) and the
ferromagnetic branch are not implemented here.

# References
- G. H. Wannier, *Phys. Rev.* **79**, 357 (1950).
"""
function fetch(
    m::IsingTriangular, ::ZZCorrelation{:static}, ::Infinite; r::Integer=1, J::Real=m.J
)
    J > 0 || throw(
        DomainError(
            J,
            "IsingTriangular T=0 ⟨σσ⟩ here is the frustrated-AFM (J>0) Wannier " *
            "result; the FM branch is not implemented (⟨σσ⟩ = ε/3J via Energy).",
        ),
    )
    r == 1 || throw(
        ArgumentError(
            "IsingTriangular ZZCorrelation: only nearest-neighbour r=1 (Wannier -1/3) " *
            "is implemented; general r (Stephenson 1964) is deferred to a follow-up.",
        ),
    )
    return -1 / 3
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: critical exponents (universality delegation)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::IsingTriangular, ::CriticalExponents, ::Infinite; kwargs...) -> NamedTuple

2D Ising universality critical exponents (Onsager 1944), shared by the
square and triangular lattices via universality:

    α = 0,  β = 1/8,  γ = 7/4,  δ = 15,  ν = 1,  η = 1/4.

Delegated to `Universality(:Ising)` at `d = 2`.  The triangular and
square 2D Ising lattices have different microscopic T_c (Onsager's
``2/log(1+sqrt(2))`` for the square; Houtappel's ``4/log 3`` for the
FM triangular) but identical universal exponents — the canonical
textbook example of universality.

# References

- L. Onsager, *Phys. Rev.* **65**, 117 (1944).
- R. M. F. Houtappel, *Physica* **16**, 425 (1950) — exact triangular-lattice
  Ising solution.
"""
function fetch(::IsingTriangular, ::CriticalExponents, ::Infinite; kwargs...)
    return QAtlas.fetch(QAtlas.Universality(:Ising), CriticalExponents(); d=2)
end
