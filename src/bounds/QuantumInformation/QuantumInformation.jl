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
