# ─────────────────────────────────────────────────────────────────────────────
# Ising universality class — exact (d=2) and numerical (d=3) exponents
#
# d=2 exponent provenance (each value individually sourced):
#   α = 0  : Onsager (1944) Phys. Rev. 65, 117 — exact specific heat
#            shows ln|T−Tc| divergence, hence α = 0.
#   β = 1/8: Yang (1952) Phys. Rev. 85, 808 — spontaneous magnetization
#            M(T) = (1−sinh⁻⁴(2βJ))^{1/8} gives β = 1/8 directly.
#   γ = 7/4: Fisher (1964) J. Math. Phys. 5, 944 — susceptibility.
#            Also follows from scaling relation γ = β(δ−1) = (1/8)(14).
#   η = 1/4: Kadanoff (1966) Physics 2, 263; Fisher (1964) — anomalous
#            dimension from correlation function G(r) ~ r^{-1/4} at Tc.
#   δ = 15 : Derived from scaling relation δ = 1 + γ/β = 1 + 14 = 15.
#            Independently verified via exact lattice calculations.
#   ν = 1  : den Nijs (1979) J. Phys. A 12, 1857 — Coulomb gas mapping.
#            Also follows from Josephson relation 2−α = dν → ν = 1 (d=2).
#   c = 1/2: Belavin, Polyakov, Zamolodchikov (1984) Nucl. Phys. B 241,
#            333 — the 2D Ising CFT is the minimal model M(3,4) with
#            central charge c = 1/2.
#
# d=3 numerical provenance:
#   All values from Kos, Poland, Simmons-Duffin, Vichi (2016) JHEP 08,
#   036, Table 2 ("3d Ising" row) — conformal bootstrap rigorous bounds.
#   These are currently the most precise known estimates.
# ─────────────────────────────────────────────────────────────────────────────

"""
    Ising2D <: AbstractQAtlasModel

Dispatch tag for the **two-dimensional** Ising universality class.
Acts as a convenience wrapper over `Universality(:Ising)` with the
dimension pinned to `d = 2`, so call sites do not have to repeat
`d = 2` on every `fetch`:

```julia
QAtlas.fetch(Ising2D(), CriticalExponents())
# ≡ QAtlas.fetch(Universality(:Ising), CriticalExponents(); d = 2)
```

The struct is kept as a distinct nominal type (rather than
`const Ising2D = Universality{:Ising}`) because the `d` kwarg is
fixed at 2 here — a type-level `const` alias would defeat that
fix. Subtyped to `AbstractQAtlasModel` so it composes with the
canonical `fetch(::AbstractQAtlasModel, ::AbstractQuantity, ::BoundaryCondition)`
signature.

See also: [`Universality`](@ref), [`CriticalExponents`](@ref).
"""
struct Ising2D <: AbstractQAtlasModel end

"""
    fetch(::Ising2D, ::CriticalExponents) -> NamedTuple

Backward-compatible alias for `fetch(Universality(:Ising), CriticalExponents(); d=2)`.
"""
function fetch(::Ising2D, ::CriticalExponents; kwargs...)
    return fetch(Universality(:Ising), CriticalExponents(); d=2, kwargs...)
end

"""
    fetch(::Universality{:Ising}, ::CriticalExponents; d) -> NamedTuple

Critical exponents of the Ising universality class (Z₂ symmetry).

- **d = 2**: Exact rational values (CFT minimal model M(3,4), c = 1/2).
- **d = 3**: High-precision numerical estimates from the conformal
  bootstrap (Kos et al. 2016). Fields `α_err`, `β_err`, … give the
  uncertainty in the last digits.
- **d ≥ 4**: Mean-field (Landau) exponents (upper critical dimension).
"""
function fetch(::Universality{:Ising}, ::CriticalExponents; d::Int, kwargs...)
    if d == 2
        return (α=0 // 1, β=1 // 8, γ=7 // 4, δ=15 // 1, ν=1 // 1, η=1 // 4, c=1 // 2)
    elseif d == 3
        # Conformal bootstrap: Kos, Poland, Simmons-Duffin, Vichi (2016)
        return (
            α=0.11009,
            α_err=0.00001,
            β=0.32642,
            β_err=0.00001,
            γ=1.23708,
            γ_err=0.00001,
            δ=4.78984,
            δ_err=0.00001,
            ν=0.62997,
            ν_err=0.00001,
            η=0.03630,
            η_err=0.00005,
        )
    elseif d >= 4
        return fetch(MeanField(), CriticalExponents())
    end
    return error("Ising universality: d=$d not supported (d ∈ {2, 3, ≥4}).")
end

# ─── Infinite-bc forwarding for verify() integration ───────────────────────
#
# CriticalExponents is a universality-class property, not a finite-lattice
# observable; the Infinite-bc dispatch below lets verify() cards pin
# individual exponent fields via the standard 3-arg signature shared by
# all literature-pin hubs.

function fetch(m::Universality{:Ising}, q::CriticalExponents, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
