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
# non-static observables are tagged individually. A velocity is :transport, NOT :dynamic — so
# `search(dynamical=:dynamic)` honestly returns empty until real spectral/real-time data exists
# (e.g. a future dynamical critical exponent or A(ω) would be tagged :dynamic here).
dynamical_axis(::Type) = :static
dynamical_axis(::Type{<:AbstractVelocity}) = :transport

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
