# core/axes.jl — orthogonal physical axes for a registered hub, captured at @register time.
#
# A single `regime` label conflated several independent physical dimensions, so
# "finite-temperature dynamics" was inexpressible and velocities masqueraded as "dynamics". These
# axes decompose that: each hub carries an independent `thermal` and `dynamical` tag, AND-queryable.
#
# Derivation priority (applied in register!): explicit @register kwarg > quantity-type trait >
# fetch-signature introspection (does the method declare a β?) > :unknown. The honest default is
# :unknown — search abstains rather than guessing, preserving the anti-hallucination property.
# Because the tags live on the REGISTRY row (not on the quantity type), a hub whose implementation
# deviates from its quantity's usual axis can be tagged individually via
# `@register(...; thermal=…, dynamical=…)`.

# ── thermal axis: :zero (T=0 only) / :finite (T>0) / :both / :unknown ──
thermal_axis(::Type) = :unknown
thermal_axis(::Type{<:AbstractThermalPotential}) = :finite   # FreeEnergy/SpecificHeat/ThermalEntropy
thermal_axis(::Type{<:Energy}) = :both                       # ⟨H⟩: ground state (no β) or thermal (β)
thermal_axis(::Type{<:AbstractGap}) = :zero                  # a gap is a ground-state spectral property
thermal_axis(::Type{<:AbstractVelocity}) = :zero             # a velocity is a T=0 / low-energy property
thermal_axis(::Type{<:AbstractEntanglementMeasure}) = :zero  # GS entanglement (quench would override)
thermal_axis(::Type{<:AbstractMagnetization}) = :both        # spontaneous (T=0) or thermal (T>0)
thermal_axis(::Type{<:AbstractSusceptibility}) = :both
thermal_axis(::Type{<:AbstractTwoPointCorrelation}) = :both
thermal_axis(::Type{<:AbstractStructureFactor}) = :both

# ── dynamical axis: :static / :transport / :dynamic / :unknown ──
# QAtlas is overwhelmingly equilibrium, so :static is the honest positive default; the few genuinely
# non-static observables are tagged individually. A velocity is :transport, NOT :dynamic.
# Real-time / non-equilibrium SCHEMES (the type parameter of a scheme-tagged quantity) make a hub
# :dynamic — so a `{:dynamic}`/`{:lightcone}` correlation, a `{:quench}` entanglement/magnetization,
# or a Loschmidt echo is honestly surfaced by `search(dynamical=:dynamic)` (the availability-review
# follow-up). Equilibrium schemes (:static / :connected / :equilibrium) stay :static; extend the set
# by adding a scheme to `_dynamic_scheme`.
dynamical_axis(::Type) = :static
dynamical_axis(::Type{<:AbstractVelocity}) = :transport

_dynamic_scheme(s) = s in (:dynamic, :lightcone, :quench)
function dynamical_axis(::Type{<:XXCorrelation{M}}) where {M}
    return _dynamic_scheme(M) ? :dynamic : :static
end
function dynamical_axis(::Type{<:YYCorrelation{M}}) where {M}
    return _dynamic_scheme(M) ? :dynamic : :static
end
function dynamical_axis(::Type{<:ZZCorrelation{M}}) where {M}
    return _dynamic_scheme(M) ? :dynamic : :static
end
function dynamical_axis(::Type{<:VonNeumannEntropy{M}}) where {M}
    return _dynamic_scheme(M) ? :dynamic : :static
end
function dynamical_axis(::Type{<:MagnetizationXLocal{M}}) where {M}
    return _dynamic_scheme(M) ? :dynamic : :static
end
dynamical_axis(::Type{<:LoschmidtEcho}) = :dynamic  # Loschmidt echo is a real-time quench quantity

# ── @register-time derivation: explicit > quantity trait > fetch-introspection > :unknown ──
function _derive_thermal(explicit, M::Type, Q::Type, BC::Type)
    explicit === nothing || return explicit
    t = thermal_axis(Q)
    t === :unknown || return t
    return _fetch_declares_beta(M, Q, BC) ? :finite : :unknown
end
_derive_dynamical(explicit, Q::Type) = explicit === nothing ? dynamical_axis(Q) : explicit

# Signature introspection: does the (M,Q,BC) fetch method declare a temperature kwarg? Deterministic
# (thermal fetches declare `; beta::Real, …`), and only a fallback — the quantity trait is primary.
function _fetch_declares_beta(M::Type, Q::Type, BC::Type)
    try
        for mth in methods(fetch, Tuple{M,Q,BC})
            any(in((:beta, :β, :T, :temperature)), Base.kwarg_decl(mth)) && return true
        end
    catch
    end
    return false
end

"""
    quantity_family(Q::Type) -> Symbol

The quantity's super-family — a TOTAL, regular classification into one of
`:correlation`, `:structure_factor`, `:magnetization`, `:susceptibility`,
`:thermodynamic`, `:gap`, `:entanglement`, `:velocity`, or `:other`. The `family`
search facet is keyed on this; because the classification is total (with `:other`
the honest catch-all), no quantity class can silently fall through the coarse
`regime` facet — the design defect the availability review surfaced. Defined on
the abstract quantity families in `quantities.jl` (all loaded before this file,
like `thermal_axis`), so a concrete quantity dispatches to its family's method.
"""
quantity_family(::Type) = :other
quantity_family(::Type{<:AbstractTwoPointCorrelation}) = :correlation
quantity_family(::Type{<:AbstractStructureFactor}) = :structure_factor
quantity_family(::Type{<:AbstractMagnetization}) = :magnetization
quantity_family(::Type{<:AbstractSusceptibility}) = :susceptibility
quantity_family(::Type{<:AbstractThermalPotential}) = :thermodynamic
quantity_family(::Type{<:Energy}) = :thermodynamic
quantity_family(::Type{<:AbstractGap}) = :gap
quantity_family(::Type{<:AbstractEntanglementMeasure}) = :entanglement
quantity_family(::Type{<:AbstractVelocity}) = :velocity
