# ─────────────────────────────────────────────────────────────────────────────
# Calabrese-Cardy entanglement entropy at \`Infinite\`
# ─────────────────────────────────────────────────────────────────────────────
#
# The 1D Heisenberg chain at the isotropic SU(2) point sits in the same
# c = 1 free-boson universality class as the XXZ chain on its critical
# axial line.  Its single-interval von Neumann and Renyi entropies are
# given by the Calabrese-Cardy 2004 closed forms with c = 1.  Unlike
# the gapped TFIM (which has a mass scale 2|h - J|), the Heisenberg
# chain is gapless at all J, so only the critical regime exists.
#
# At T = 0:
#
#     S_VN(ℓ)     = (1/3)  · log(2 ℓ)
#     S_α(ℓ)     = (1/6)(1 + 1/α) · log(2 ℓ)
#
# At finite β > 0:
#
#     S_VN(ℓ, β) = (1/3)  · log[(2 β / π) · sinh(π ℓ / β)]
#     S_α(ℓ, β) = (1/6)(1 + 1/α) · log[(2 β / π) · sinh(π ℓ / β)]
#
# Conventions match TFIM_cft_entanglement.jl: β is the lattice-units
# inverse temperature with the Fermi velocity normalised to unity
# implicitly.  The non-universal `S_0` offset is dropped (universal
# part only).
#
# Reference: Calabrese-Cardy, J. Stat. Mech. P06002 (2004), §4.
# Tracking: #580 entanglement universality catalog.
# ─────────────────────────────────────────────────────────────────────────────

"""
    _heisenberg1d_cc_entanglement(ℓ, β; α=1.0) -> Float64

Calabrese-Cardy entanglement entropy of a contiguous interval of length
`ℓ` in the infinite Heisenberg1D chain (c = 1, always critical).

- `β = Inf`: T = 0 ground state, returns `prefac · log(2 ℓ)`.
- `0 < β < Inf`: finite-temperature thermal state, returns
  `prefac · log[(2 β / π) · sinh(π ℓ / β)]`.

`prefac = c/3 = 1/3` for von Neumann, `prefac = (c/6)(1 + 1/α)` for
Renyi-α.  The non-universal `S_0` offset is dropped.
"""
function _heisenberg1d_cc_entanglement(ℓ::Integer, β::Real; α::Real=1.0)::Float64
    c = 1.0
    prefac = α == 1 ? c / 3 : (c / 6) * (1 + 1 / α)
    if isinf(β)
        return prefac * log(2 * ℓ)
    else
        return prefac * log((2 * β / π) * sinh(π * ℓ / β))
    end
end

"""
    fetch(model::Heisenberg1D, ::VonNeumannEntropy{:equilibrium}, ::Infinite;
          ℓ::Int, beta::Real = Inf, kwargs...) -> Float64

Calabrese-Cardy von Neumann entanglement entropy of a contiguous
block of length `ℓ` in the infinite 1D Heisenberg chain at the SU(2)
point.  The chain is gapless (c = 1) for all `J`, so the critical
Calabrese-Cardy form applies for all `β`.

# Examples

```julia
julia> QAtlas.fetch(Heisenberg1D(), VonNeumannEntropy(), Infinite(); ℓ=10)
0.998529...
```

Refs Calabrese-Cardy, J. Stat. Mech. P06002 (2004); issue #580.
"""
function fetch(
    ::Heisenberg1D,
    ::VonNeumannEntropy{:equilibrium},
    ::Infinite;
    ℓ::Int,
    beta::Real=Inf,
    kwargs...,
)
    ℓ ≥ 1 || throw(
        ArgumentError(
            "Heisenberg1D VonNeumannEntropy Infinite: ℓ must be ≥ 1; got ℓ = \$ℓ."
        ),
    )
    return _heisenberg1d_cc_entanglement(ℓ, beta; α=1.0)
end

"""
    fetch(model::Heisenberg1D, q::RenyiEntropy, ::Infinite;
          ℓ::Int, beta::Real = Inf, kwargs...) -> Float64

Calabrese-Cardy Renyi-α entanglement entropy of a contiguous block of
length `ℓ` in the infinite Heisenberg1D chain.  Coefficient

    P_α = (c / 6) · (1 + 1/α),  c = 1.

# Examples

```julia
julia> QAtlas.fetch(Heisenberg1D(), RenyiEntropy(2.0), Infinite(); ℓ=10)
0.748897...
```

Refs Calabrese-Cardy, J. Stat. Mech. P06002 (2004); issue #580.
"""
function fetch(
    ::Heisenberg1D, q::RenyiEntropy, ::Infinite; ℓ::Int, beta::Real=Inf, kwargs...
)
    ℓ ≥ 1 || throw(
        ArgumentError("Heisenberg1D RenyiEntropy Infinite: ℓ must be ≥ 1; got ℓ = \$ℓ.")
    )
    return _heisenberg1d_cc_entanglement(ℓ, beta; α=q.α)
end
