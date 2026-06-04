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
