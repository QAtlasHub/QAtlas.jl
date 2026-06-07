# ─────────────────────────────────────────────────────────────────────────────
# bounds/Bounds.jl — model-independent universal bounds namespace.
#
# A *bound* is NOT a universality class.  It is pinned by three things:
#   * what physical quantity it bounds  — the registry `quantity`,
#   * which way it constrains           — `direction = :upper / :lower`,
#   * whose bound it is                 — `references` (+ a `source=` selector
#                                          when several bounds share a quantity).
#
# Universal (model-independent) bounds — Tsirelson, Bekenstein,
# Maldacena–Shenker–Stanford, … — have no home model, so they dispatch on
# `Bound{Domain}`, the bounds-namespace analogue of `Universality{Class}`.
# Model-specific bounds (e.g. a TFIM Lieb–Robinson velocity) instead stay on
# their model with `status=:bound`; only the convention is shared.
#
# Layout mirrors `models/` and `universalities/`: one directory per domain,
# `bounds/<Domain>/<Domain>.jl` + `bounds/<Domain>/<Domain>_registry.jl`.
# ─────────────────────────────────────────────────────────────────────────────

"""
    Bound{D}

Dispatch tag for a *domain* of model-independent universal bounds — the
bounds-namespace analogue of [`Universality`](@ref).  `D` is a `Symbol`
naming the domain (`:QuantumInformation`, `:Dynamics`, `:Holographic`, …).

A universal bound has no home model, so it is fetched against a `Bound`
domain rather than a Hamiltonian:

```julia
QAtlas.fetch(Bound(:QuantumInformation), CHSHBound(), Infinite())                # 2√2 (Tsirelson)
QAtlas.fetch(Bound(:QuantumInformation), CHSHBound(), Infinite(); source=:bell)  # 2   (local-hidden-variable)
```

Every bound is registered with `status=:bound` and a `direction`
(`:upper`/`:lower`); see [`BOUND_DIRECTIONS`](@ref).  Model-*specific*
bounds (e.g. a TFIM Lieb–Robinson velocity) carry the same convention but
live on their model, not here.
"""
struct Bound{D} <: AbstractQAtlasModel end
Bound(domain::Symbol) = Bound{domain}()
