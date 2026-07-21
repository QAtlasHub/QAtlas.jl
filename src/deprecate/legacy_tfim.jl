# deprecate/legacy_tfim.jl — backward-compat shims for the TFIM Symbol API.
#
# Routes the pre-v0.13 call patterns:
#   fetch(Model(:TFIM; J, h, N), Quantity(:energy), OBC(); …)
#   fetch(:TFIM, :energy, OBC(); N, J, h, beta)
# into the canonical concrete-struct API defined in src/models/quantum/TFIM*.jl.
#
# Scheduled for removal in v1.0.

"""
    _tfim_from_legacy_model(m::Model{:TFIM}) -> TFIM

Extract `J`, `h` from the legacy `Model{:TFIM}(params)` Dict and build
the concrete struct.
"""
function _tfim_from_legacy_model(m::Model{:TFIM})::TFIM
    J = Float64(get(m.params, :J, 1.0))
    h = Float64(get(m.params, :h, 1.0))
    return TFIM(; J=J, h=h)
end

"""
    _bc_with_legacy_N(bc::OBC, m::Model{:TFIM}) -> OBC

If `bc.N == 0` (zero-arg `OBC()` call) but the legacy `Model` Dict has
an `:N` entry, promote it to a sized `OBC(N)` so the concrete fetch can
read `bc.N` uniformly.
"""
function _bc_with_legacy_N(bc::OBC, m::Model{:TFIM})::OBC
    bc.N > 0 && return bc
    haskey(m.params, :N) && return OBC(Int(m.params[:N]))
    return bc
end
_bc_with_legacy_N(bc::BoundaryCondition, ::Model{:TFIM}) = bc

# ── Energy ──────────────────────────────────────────────────────────────
function fetch(m::Model{:TFIM}, ::Quantity{:energy}, bc::OBC; kwargs...)
    return fetch(_tfim_from_legacy_model(m), Energy(), _bc_with_legacy_N(bc, m); kwargs...)
end
function fetch(m::Model{:TFIM}, ::Quantity{:energy}, bc::Infinite; kwargs...)
    return fetch(_tfim_from_legacy_model(m), Energy(), bc; kwargs...)
end

# ── Central charge ──────────────────────────────────────────────────────
function fetch(m::Model{:TFIM}, ::Quantity{:central_charge}, bc::Infinite; kwargs...)
    return fetch(_tfim_from_legacy_model(m), CentralCharge(), bc; kwargs...)
end

# ── Mass gap ────────────────────────────────────────────────────────────
function fetch(m::Model{:TFIM}, ::Quantity{:mass_gap}, bc::Infinite; kwargs...)
    return fetch(_tfim_from_legacy_model(m), MassGap(), bc; kwargs...)
end
function fetch(m::Model{:TFIM}, ::Quantity{:mass_gap}, bc::OBC; kwargs...)
    return fetch(_tfim_from_legacy_model(m), MassGap(), _bc_with_legacy_N(bc, m); kwargs...)
end

# ── Entanglement entropy ────────────────────────────────────────────────
function fetch(m::Model{:TFIM}, ::Quantity{:entanglement_entropy}, bc::OBC; kwargs...)
    return fetch(
        _tfim_from_legacy_model(m), VonNeumannEntropy(), _bc_with_legacy_N(bc, m); kwargs...
    )
end

# ── Thermal quantities (TFIM_thermal.jl) ────────────────────────────────
# Mapping: legacy Symbol Quantity → concrete struct used by the new API.
const _LEGACY_TFIM_THERMAL_MAP = (
    (:free_energy, FreeEnergy),
    (:entropy, ThermalEntropy),
    (:specific_heat, SpecificHeat),
    (:transverse_magnetization, MagnetizationX),
    (:transverse_susceptibility, SusceptibilityXX),
)

for (qsym, QTy) in _LEGACY_TFIM_THERMAL_MAP
    @eval begin
        function fetch(
            m::Model{:TFIM}, ::Quantity{$(QuoteNode(qsym))}, bc::Infinite; kwargs...
        )
            return fetch(_tfim_from_legacy_model(m), $QTy(), bc; kwargs...)
        end
        function fetch(m::Model{:TFIM}, ::Quantity{$(QuoteNode(qsym))}, bc::OBC; kwargs...)
            return fetch(
                _tfim_from_legacy_model(m), $QTy(), _bc_with_legacy_N(bc, m); kwargs...
            )
        end
    end
end

# ── Site-local quantities (TFIM_local.jl) ──────────────────────────────
const _LEGACY_TFIM_LOCAL_MAP = (
    (:magnetization_x_local, LocalMagnetization{:x}),
    (:magnetization_z_local, LocalMagnetization{:z}),
    (:energy_local, EnergyLocal),
)
for (qsym, QTy) in _LEGACY_TFIM_LOCAL_MAP
    @eval function fetch(
        m::Model{:TFIM}, ::Quantity{$(QuoteNode(qsym))}, bc::OBC; kwargs...
    )
        return fetch(
            _tfim_from_legacy_model(m), $QTy(), _bc_with_legacy_N(bc, m); kwargs...
        )
    end
end

# ── Dynamics / correlators / structure factors (TFIM_dynamics.jl) ──────
# These route the old Symbol-keyed quantities onto the new axis-parametric
# correlator types (#734).
function fetch(m::Model{:TFIM}, ::Quantity{:sz_sz_correlation}, bc::OBC; kwargs...)
    return fetch(
        _tfim_from_legacy_model(m),
        DynamicalCorrelation(:z, :z),
        _bc_with_legacy_N(bc, m);
        kwargs...,
    )
end
function fetch(m::Model{:TFIM}, ::Quantity{:sx_sx_correlation}, bc::OBC; kwargs...)
    return fetch(
        _tfim_from_legacy_model(m),
        DynamicalCorrelation(:x, :x),
        _bc_with_legacy_N(bc, m);
        kwargs...,
    )
end
function fetch(m::Model{:TFIM}, ::Quantity{:sz_sz_spreading}, bc::OBC; kwargs...)
    return fetch(
        _tfim_from_legacy_model(m),
        LightconeSpinCorrelation(:z, :z),
        _bc_with_legacy_N(bc, m);
        kwargs...,
    )
end
function fetch(m::Model{:TFIM}, ::Quantity{:zz_static_thermal}, bc::OBC; kwargs...)
    return fetch(
        _tfim_from_legacy_model(m),
        SpinCorrelation(:z, :z),
        _bc_with_legacy_N(bc, m);
        kwargs...,
    )
end
function fetch(m::Model{:TFIM}, ::Quantity{:zz_structure_factor}, bc::OBC; kwargs...)
    return fetch(
        _tfim_from_legacy_model(m), ZZStructureFactor(), _bc_with_legacy_N(bc, m); kwargs...
    )
end
function fetch(
    m::Model{:TFIM}, ::Quantity{:longitudinal_susceptibility}, bc::OBC; kwargs...
)
    return fetch(
        _tfim_from_legacy_model(m), SusceptibilityZZ(), _bc_with_legacy_N(bc, m); kwargs...
    )
end
