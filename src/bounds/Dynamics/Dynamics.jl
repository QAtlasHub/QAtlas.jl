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
    fetch(::Bound{:Dynamics}, ::QuantumSpeedLimit, ::Infinite;
          scheme=:margolus_levitin, E, ΔE)

Quantum speed limit — a *lower* bound on the time to evolve a state to an
orthogonal one (ℏ = 1).  Two complementary regimes, selected by `scheme`:

- `:margolus_levitin` (default, canonical) — Margolus–Levitin 1998,

      τ ≥ π / (2E),

  with `E` the mean energy above the ground state.

- `:mandelstam_tamm` — Mandelstam–Tamm 1945,

      τ ≥ π / (2 ΔE),

  with `ΔE` the energy uncertainty (root variance of the Hamiltonian in the
  initial state).

The true speed limit is the tighter (larger) of the two; both are
`direction=:lower` bounds, each saturated in its own regime.
"""
function fetch(
    ::Bound{:Dynamics},
    ::QuantumSpeedLimit,
    ::Infinite;
    scheme::Symbol=:margolus_levitin,
    E::Real=NaN,
    ΔE::Real=NaN,
    kwargs...,
)
    if scheme === :margolus_levitin
        E > 0 || throw(
            ArgumentError(
                "QuantumSpeedLimit(:margolus_levitin): E must be positive; got $(E)"
            ),
        )
        return π / (2 * E)
    elseif scheme === :mandelstam_tamm
        ΔE > 0 || throw(
            ArgumentError(
                "QuantumSpeedLimit(:mandelstam_tamm): ΔE must be positive; got $(ΔE)"
            ),
        )
        return π / (2 * ΔE)
    else
        throw(
            ArgumentError(
                "QuantumSpeedLimit: scheme must be :margolus_levitin or :mandelstam_tamm; got :$(scheme)",
            ),
        )
    end
end

"""
    fetch(::Bound{:Dynamics}, ::ScramblingTime, ::Infinite; β, N)

Sekino–Susskind 2008 fast-scrambling time — a conjectured *lower* bound on the
time for a thermal system of `N` effective degrees of freedom to scramble local
information into global entanglement,

    t_* = (β / 2π) log N        (ℏ = k_B = 1).

Saturated by black holes (the fastest scramblers); for local systems `t_*` is
parametrically longer.  `direction=:lower`.
"""
function fetch(
    ::Bound{:Dynamics}, ::ScramblingTime, ::Infinite; β::Real, N::Real, kwargs...
)
    β > 0 || throw(ArgumentError("ScramblingTime: β must be positive; got $(β)"))
    N > 1 || throw(ArgumentError("ScramblingTime: N must be > 1; got $(N)"))
    return β * log(N) / (2π)
end
