# ─────────────────────────────────────────────────────────────────────────────
# bounds/QuantumInformation — quantum-information bounds (Bell / CHSH, Mermin).
#
# Extracted from the former `Universality(:QuantumMechanics)` dumping ground:
# these are *bounds* (upper bounds on Bell-type correlators), not a universality
# class.  Each is registered once per theory regime
# (`scheme = :classical / :quantum / :no_signalling`); whose bound each regime
# is lives in the per-scheme `references` — list them with
# `definitions(Bound(:QuantumInformation), Q)`.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Bound{:QuantumInformation}, ::CHSHBound, ::Infinite; scheme=:quantum)

Upper bound on the CHSH correlator
`S = E(a,b) + E(a,b′) + E(a′,b) − E(a′,b′)`, by theory regime:

  * `:classical`         → 2     local-hidden-variable bound (CHSH 1969)
  * `:quantum` (default) → 2√2   Tsirelson 1980 quantum maximum, saturated by
                                 the optimal Bell state
  * `:no_signalling`     → 4     Popescu–Rohrlich 1994 algebraic maximum

The ladder 2 < 2√2 < 4 is classical → quantum → no-signalling; the
quantum/classical ratio √2 is the Bell-violation factor.  Whose bound each
regime is lives in the per-scheme `references` (see [`definitions`](@ref)).
"""
function fetch(
    ::Bound{:QuantumInformation},
    ::CHSHBound,
    ::Infinite;
    scheme::Symbol=:quantum,
    kwargs...,
)
    if scheme === :quantum
        return 2 * sqrt(2)
    elseif scheme === :classical
        return 2.0
    elseif scheme === :no_signalling
        return 4.0
    else
        throw(
            ArgumentError(
                "CHSHBound: scheme must be :classical / :quantum / :no_signalling; " *
                "got :$(scheme)",
            ),
        )
    end
end

"""
    fetch(::Bound{:QuantumInformation}, ::MerminGHZBound, ::Infinite; scheme=:quantum)

Upper bound on the Mermin 3-party operator `|⟨M₃⟩|`, by theory regime:

  * `:classical`         → 2   local-realistic bound (Mermin 1990)
  * `:quantum` (default) → 4   quantum maximum (Mermin 1990), saturated by the
                               GHZ state `(|000⟩ + |111⟩)/√2`

The quantum/classical gap is 2 (vs √2 for two-party CHSH): multipartite
nonlocality grows with the number of parties.
"""
function fetch(
    ::Bound{:QuantumInformation},
    ::MerminGHZBound,
    ::Infinite;
    scheme::Symbol=:quantum,
    kwargs...,
)
    if scheme === :quantum
        return 4.0
    elseif scheme === :classical
        return 2.0
    else
        throw(
            ArgumentError(
                "MerminGHZBound: scheme must be :classical / :quantum; got :$(scheme)"
            ),
        )
    end
end

"""
    fetch(::Bound{:QuantumInformation}, ::OptimalCloningFidelity, ::Infinite)

Bužek–Hillery 1996 upper bound on the single-copy fidelity of a universal,
symmetric `1 → 2` quantum cloner of a qubit,

    F ≤ 5/6.

The no-cloning theorem forbids `F = 1`; `5/6` is the best achievable, saturated
by the optimal universal cloner.  A `status=:bound`, `direction=:upper` claim.
"""
function fetch(
    ::Bound{:QuantumInformation}, ::OptimalCloningFidelity, ::Infinite; kwargs...
)
    return 5 / 6
end

"""
    fetch(::Bound{:QuantumInformation}, ::BB84KeyRate, ::Infinite; qber)

Shor–Preskill 2000 asymptotic secret-key rate of the BB84 protocol,

    R(e) = 1 − 2 H₂(e),     H₂(e) = −e log₂ e − (1−e) log₂(1−e),

with `e = qber` the qubit error rate.  A provably achievable rate — a *lower*
bound on the extractable secret-key fraction.  `R = 1` at `e = 0`, decreasing to
`R = 0` at the `e ≈ 11%` security threshold and negative beyond (no secure key).
A `status=:bound`, `direction=:lower` claim.
"""
function fetch(
    ::Bound{:QuantumInformation}, ::BB84KeyRate, ::Infinite; qber::Real, kwargs...
)
    0 <= qber <= 0.5 ||
        throw(ArgumentError("BB84KeyRate: qber must satisfy 0 ≤ qber ≤ 0.5; got $(qber)"))
    qber == 0 && return 1.0
    h2 = -qber * log2(qber) - (1 - qber) * log2(1 - qber)
    return 1 - 2 * h2
end
