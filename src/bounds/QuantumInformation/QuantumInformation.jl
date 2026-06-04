# ─────────────────────────────────────────────────────────────────────────────
# bounds/QuantumInformation — quantum-information bounds (Bell / CHSH, …).
#
# Extracted from the former `Universality(:QuantumMechanics)` dumping ground:
# the CHSH bounds are *bounds* (an upper bound on a Bell correlator), not a
# universality class.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Bound{:QuantumInformation}, ::CHSHBound, ::Infinite; source=:tsirelson)

Upper bound on the CHSH correlator

    S = E(a,b) + E(a,b′) + E(a′,b) − E(a′,b′),

selected by `source` — *whose* bound (all are upper bounds):

  * `:bell`               → 2     local-hidden-variable bound (CHSH 1969)
  * `:tsirelson` (default)→ 2√2   quantum maximum (Tsirelson 1980), saturated
                                  by the optimal Bell state
  * `:popescu_rohrlich`   → 4     algebraic / no-signalling maximum
                                  (Popescu–Rohrlich 1994)

The ladder 2 < 2√2 < 4 is classical → quantum → no-signalling; the
quantum/classical ratio √2 is the Bell-violation factor.
"""
function fetch(
    ::Bound{:QuantumInformation},
    ::CHSHBound,
    ::Infinite;
    source::Symbol=:tsirelson,
    kwargs...,
)
    if source === :tsirelson
        return 2 * sqrt(2)
    elseif source === :bell
        return 2.0
    elseif source === :popescu_rohrlich
        return 4.0
    else
        throw(
            ArgumentError(
                "CHSHBound: source must be one of :bell / :tsirelson / " *
                ":popescu_rohrlich (whose bound); got :$(source)",
            ),
        )
    end
end

"""
    fetch(::Bound{:QuantumInformation}, ::MerminGHZBound, ::Infinite; source=:mermin)

Upper bound on the Mermin 3-party operator `|⟨M₃⟩|`, by `source`:

  * `:classical`        → 2   local-hidden-variable (local-realistic) bound
  * `:mermin` (default) → 4   quantum maximum (Mermin 1990), saturated by the
                              GHZ state `(|000⟩ + |111⟩)/√2`

The quantum/classical gap is 2 (vs √2 for two-party CHSH): multipartite
nonlocality grows with the number of parties.
"""
function fetch(
    ::Bound{:QuantumInformation},
    ::MerminGHZBound,
    ::Infinite;
    source::Symbol=:mermin,
    kwargs...,
)
    if source === :mermin
        return 4.0
    elseif source === :classical
        return 2.0
    else
        throw(
            ArgumentError(
                "MerminGHZBound: source must be :classical or :mermin (whose bound); " *
                "got :$(source)",
            ),
        )
    end
end
