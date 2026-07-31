# deprecate/legacy_xxz.jl — Symbol-dispatch shim for `XXZ1D`.
#
# The canonical API is the concrete struct form:
#
#   fetch(XXZ1D(; J, Δ), Energy(), Infinite())
#
# Legacy Symbol callers route through Model{:XXZ1D}(params::Dict) which
# is constructed by `Model(:XXZ; kwargs...)` via canonicalize_model
# (alias `:XXZ → :XXZ1D`).  The shims below translate that into the new
# concrete-struct dispatch.

"""
    _xxz1d_from_legacy_model(m::Model{:XXZ1D}) -> XXZ1D
"""
function _xxz1d_from_legacy_model(m::Model{:XXZ1D})::XXZ1D
    J = Float64(get(m.params, :J, 1.0))
    Δ = Float64(get(m.params, :Δ, 0.0))
    return XXZ1D(; J=J, Δ=Δ)
end

# One shim per supported (Symbol quantity, BC) pair.
const _LEGACY_XXZ1D_QUANTITY_MAP = (
    (:energy, Energy),
    (:ground_state_energy, Energy{:per_site}),
    (:central_charge, CentralCharge),
    (:luttinger_parameter, LuttingerParameter),
    (:luttinger_velocity, LuttingerVelocity),
    (:sound_velocity, LuttingerVelocity),     # alias: same physical quantity
    (:spin_wave_velocity, LuttingerVelocity), # alias
    (:fermi_velocity, LuttingerVelocity),     # alias (free-fermion XX limit)
)

for (qsym, QTy) in _LEGACY_XXZ1D_QUANTITY_MAP
    @eval function fetch(
        m::Model{:XXZ1D}, ::Quantity{$(QuoteNode(qsym))}, bc::Infinite; kwargs...
    )
        return fetch(_xxz1d_from_legacy_model(m), $QTy(), bc; kwargs...)
    end
end
