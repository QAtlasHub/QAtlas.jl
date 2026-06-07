# ─────────────────────────────────────────────────────────────────────────────
# bounds/Dynamics — dynamical bounds (chaos, scrambling, speed limits, …).
#
# Extracted from the former `Universality(:QuantumMechanics)` dumping ground:
# these are model-independent *bounds* on dynamics, not a universality class.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Bound{:Dynamics}, ::ChaosBound, ::Infinite; β)

Maldacena–Shenker–Stanford 2016 upper bound on the quantum Lyapunov
exponent of out-of-time-order correlators,

    λ_L ≤ 2π / β        (ℏ = k_B = 1),

i.e. `λ_L ≤ 2π k_B T / ℏ`.  Saturated by maximally chaotic systems
(holographic / large-N SYK).  A `status=:bound`, `direction=:upper` claim.
"""
function fetch(::Bound{:Dynamics}, ::ChaosBound, ::Infinite; β::Real, kwargs...)
    β > 0 || throw(ArgumentError("ChaosBound: β must be positive; got $(β)"))
    return 2π / β
end

"""
    fetch(::Bound{:Dynamics}, ::QuantumSpeedLimit, ::Infinite; scheme=:margolus_levitin, E)

Margolus–Levitin 1998 quantum speed limit — a *lower* bound on the time to
evolve a state to an orthogonal one,

    τ ≥ π / (2E)        (ℏ = 1),

with `E` the mean energy above the ground state.  `direction=:lower`;
saturated by an equal superposition of the ground state and one excited
state.
"""
function fetch(
    ::Bound{:Dynamics},
    ::QuantumSpeedLimit,
    ::Infinite;
    scheme::Symbol=:margolus_levitin,
    E::Real,
    kwargs...,
)
    scheme === :margolus_levitin || throw(
        ArgumentError(
            "QuantumSpeedLimit: scheme must be :margolus_levitin; got :$(scheme)"
        ),
    )
    E > 0 || throw(ArgumentError("QuantumSpeedLimit: E must be positive; got $(E)"))
    return π / (2 * E)
end
