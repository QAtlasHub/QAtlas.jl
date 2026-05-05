# ─────────────────────────────────────────────────────────────────────────────
# Calabrese–Cardy entanglement entropy at the Universality{C} level
#
# Generic 1+1D CFT entanglement formulae for any universality class for which
# a central charge `c` is defined.  The closed forms are
#
#   PBC, finite L:   S(ℓ, L) = (c/3) log[(L/π) sin(πℓ/L)] + c'_1
#   OBC, finite L:   S(ℓ, L) = (c/6) log[(2L/π) sin(πℓ/L)] + c'_1 + log g
#   Infinite (PBC):  S(ℓ)    = (c/3) log ℓ + c'_1
#   Infinite (OBC):  S(ℓ)    = (c/6) log ℓ + c'_1 + log g
#
# The non-universal cutoff constant `c'_1` and the Affleck–Ludwig boundary
# entropy `log g` are *dropped* — they require model-specific UV input
# (lattice spacing convention) and boundary input (which conformal boundary
# state is realised) that is not available at the universality level.  What
# remains is the universal log-prefactor coefficient `(c/3)` (PBC) or
# `(c/6)` (OBC), exactly the piece that universality alone determines.
#
# The Rényi extension uses the substitution
#
#   c -> c · (1 + 1/α) / 2,
#
# which reduces to `c` at α = 1.  See Calabrese–Cardy J. Stat. Mech. P06002
# (2004) eq. (3.12) and J. Phys. A 42, 504005 (2009) eq. (28),(30).
#
# References:
#   - P. Calabrese, J. Cardy, J. Stat. Mech. P06002 (2004).
#   - P. Calabrese, J. Cardy, J. Phys. A 42, 504005 (2009).
#
# A step-by-step derivation of the c/3 vs c/6 prefactor (replica trick,
# twist-operator, cylinder-vs-strip conformal map) lives in
# `docs/src/calc/calabrese-cardy-obc-vs-pbc.md`.
# ─────────────────────────────────────────────────────────────────────────────

# ─── CentralCharge: minimal-model 1+1D CFT lookups ──────────────────────────
#
# Only universality classes whose critical point is described by a known
# 1+1D CFT have a well-defined central charge in this dispatch.  Higher-d
# universality classes (e.g. 3D Ising, 3D Heisenberg) do *not* live in a
# 1+1D CFT — there is no central charge at the universality-class level
# even though the d-dimensional class is perfectly well-defined.  Those
# call sites raise an `ErrorException` with the dimension in the message.

"""
    fetch(::Universality{:Ising}, ::CentralCharge; d::Int=2) -> Rational{Int}

Central charge of the Ising universality class as a 1+1D CFT.

Only `d = 2` is supported (the Ising minimal model `M(3,4)`, central
charge `c = 1/2`).  For `d ≥ 3` the universality class is not a
1+1D CFT and an `ErrorException` is thrown — call sites that want a
generic-CFT entanglement formula must use a 1+1D class.

Reference: Belavin, Polyakov, Zamolodchikov, Nucl. Phys. B 241, 333 (1984).
"""
function fetch(::Universality{:Ising}, ::CentralCharge; d::Int=2, kwargs...)
    if d == 2
        return 1 // 2
    end
    return error(
        "Universality{:Ising} CentralCharge: only d=2 is a 1+1D CFT (c = 1/2); " *
        "got d=$d.  The d=$d Ising universality class lives in a $(d)-dim CFT, " *
        "not a 1+1D CFT, so the Calabrese-Cardy entanglement formula does not apply.",
    )
end

"""
    fetch(::Universality{:Potts3}, ::CentralCharge; d::Int=2) -> Rational{Int}

Central charge of the 3-state Potts universality class.  `d = 2` only.
The 2D 3-state Potts model is the Virasoro minimal model `M(5,6)` with
`c = 4/5`.

Reference: Dotsenko, Nucl. Phys. B 235, 54 (1984); di Francesco–
Mathieu–Sénéchal §7.4.
"""
function fetch(::Universality{:Potts3}, ::CentralCharge; d::Int=2, kwargs...)
    if d == 2
        return 4 // 5
    end
    return error("Universality{:Potts3} CentralCharge: only d=2 supported; got d=$d.")
end

"""
    fetch(::Universality{:Potts4}, ::CentralCharge; d::Int=2) -> Rational{Int}

Central charge of the 4-state Potts universality class.  `d = 2` only.
The 2D 4-state Potts model is at the marginal compact-boson point
(self-dual radius); the central charge is `c = 1`.

Reference: di Francesco–Mathieu–Sénéchal, *Conformal Field Theory*
(Springer 1997), §12.3.
"""
function fetch(::Universality{:Potts4}, ::CentralCharge; d::Int=2, kwargs...)
    if d == 2
        return 1 // 1
    end
    return error("Universality{:Potts4} CentralCharge: only d=2 supported; got d=$d.")
end

"""
    fetch(::Universality{:XY}, ::CentralCharge; d::Int=2) -> Rational{Int}

Central charge of the XY (`O(2)`) universality class, `d = 2` only.
The 2D XY model has a Berezinskii–Kosterlitz–Thouless transition; the
critical line below `T_BKT` is described by a free compact boson with
`c = 1`.

Reference: Kosterlitz, J. Phys. C 7, 1046 (1974); di Francesco–
Mathieu–Sénéchal §6.

For `d ≥ 3` the class is not a 1+1D CFT and the call errors.
"""
function fetch(::Universality{:XY}, ::CentralCharge; d::Int=2, kwargs...)
    if d == 2
        return 1 // 1
    end
    return error(
        "Universality{:XY} CentralCharge: only d=2 (BKT free boson) is a 1+1D CFT " *
        "(c = 1); got d=$d.",
    )
end

"""
    fetch(::Universality{:Heisenberg}, ::CentralCharge; d::Int=1) -> Rational{Int}

Central charge of the spin-1/2 Heisenberg universality class.

Only `d = 1` is supported (the spin-1/2 antiferromagnetic Heisenberg
chain in the SU(2)_1 Wess-Zumino-Witten universality class, central
charge `c = 1`).  For `d ≥ 2` the Heisenberg universality class is not
a 1+1D CFT (the 2D Heisenberg model has Goldstone modes; the 3D one
has no critical line at finite T) — call sites that want a
generic-CFT entanglement formula must use a 1+1D class.

Reference: Affleck–Haldane, Phys. Rev. B 36, 5291 (1987); di Francesco–
Mathieu–Sénéchal §15.6 (SU(2)_1 WZW).
"""
function fetch(::Universality{:Heisenberg}, ::CentralCharge; d::Int=1, kwargs...)
    if d == 1
        return 1 // 1
    end
    return error(
        "Universality{:Heisenberg} CentralCharge: only d=1 (spin-1/2 Heisenberg " *
        "chain, SU(2)_1 WZW, c = 1) is a 1+1D CFT; got d=$d.",
    )
end

# ─── Calabrese–Cardy entanglement entropy: generic Universality{C} ──────────
#
# All entanglement methods route through `_cardy_central_charge(model)` to
# extract `c`.  The method errors out cleanly for any universality class
# that has no `CentralCharge` defined (KPZ, Percolation, …).

"""
    _cardy_central_charge(model::Universality{C}; kwargs...) -> Float64

Internal helper: fetch the central charge `c` of the universality class
`model` and return it as a `Float64`.  Re-throws as an `ErrorException`
with a Calabrese–Cardy-specific message if the class has no
`CentralCharge` defined.
"""
function _cardy_central_charge(model::Universality{C}; kwargs...) where {C}
    try
        return Float64(fetch(model, CentralCharge(); kwargs...))
    catch err
        if err isa MethodError || err isa ErrorException
            return error(
                "Universality{:$C}: Calabrese-Cardy entanglement requires a 1+1D CFT " *
                "central charge, which is not defined for this universality class. " *
                "Define `fetch(::Universality{:$C}, ::CentralCharge; ...)` first, or " *
                "use a class that lives in a 1+1D CFT (e.g. :Ising d=2, :Potts3 d=2, " *
                ":Potts4 d=2, :XY d=2).",
            )
        end
        rethrow(err)
    end
end

"""
    _cardy_renyi_c(c::Real, α::Real) -> Float64

Calabrese–Cardy Rényi-α coefficient substitution

    c -> c · (1 + 1/α) / 2.

Reduces to `c` at `α = 1`.  Used to produce the Rényi entropy from the
same closed form as the von Neumann case.
"""
@inline _cardy_renyi_c(c::Real, α::Real) = Float64(c) * (1 + 1 / α) / 2

# ─── Von Neumann entropy ────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::VonNeumannEntropy, ::PBC; ℓ::Real, L::Real, kwargs...)
        -> Float64

Calabrese–Cardy von Neumann entanglement entropy of a contiguous block
of length `ℓ` in a 1+1D CFT on a *periodic* chain of length `L`:

    S(ℓ, L) = (c/3) log[(L/π) sin(πℓ/L)]                 (PBC)

The non-universal additive constant `c'_1` (UV cutoff) is **dropped**.
The central charge is fetched via
`fetch(Universality{C}(), CentralCharge())`; classes without a
1+1D-CFT central charge raise `ErrorException`.

Boundary cases:

- `ℓ → 0` or `ℓ → L`: argument of the log → 0, returns `-Inf` (UV
  divergence) — physically the cut runs through zero sites.
- `ℓ = L/2`: maximum.

Reference: Calabrese–Cardy J. Stat. Mech. P06002 (2004) eq. (3.8);
J. Phys. A 42, 504005 (2009) eq. (28).
"""
function fetch(
    model::Universality{C}, ::VonNeumannEntropy, ::PBC; ℓ::Real, L::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("Cardy PBC: L must be > 0; got L=$L."))
    0 ≤ ℓ ≤ L ||
        throw(ArgumentError("Cardy PBC: ℓ must satisfy 0 ≤ ℓ ≤ L; got ℓ=$ℓ, L=$L."))
    c = _cardy_central_charge(model; kwargs...)
    # Endpoints ℓ = 0 or ℓ = L are exact UV divergences — return -Inf
    # without going through sin() (which would round to ~1e-16 at ℓ = L).
    (ℓ == 0 || ℓ == L) && return -Inf
    s = sin(π * ℓ / L)
    if s ≤ 0
        return -Inf
    end
    return (c / 3) * log((L / π) * s)
end

"""
    fetch(::Universality{C}, ::VonNeumannEntropy, ::OBC; ℓ::Real, L::Real, kwargs...)
        -> Float64

Calabrese–Cardy von Neumann entanglement entropy on an *open* chain of
length `L`, **for the canonical "block at one boundary" geometry**:
the subsystem A occupies sites 1..ℓ (or equivalently L-ℓ+1..L), so
exactly **one** of A's endpoints sits at the open boundary and only
one entanglement cut lies in the bulk:

    S(ℓ, L) = (c/6) log[(2L/π) sin(πℓ/L)]                (OBC, block at end)

For a block in the **bulk** of an open chain (e.g. sites L/4..L/4+ℓ
with both endpoints away from the boundary), there are **two** bulk
cuts and the prefactor reverts to `c/3`, with a different log
argument involving conformal cross-ratios of four points (Calabrese–
Cardy J. Phys. A 42, 504005 (2009) §3.3).  This bulk-block formula is
*not* implemented by this method.

The non-universal additive constant `c'_1` and the Affleck–Ludwig
boundary entropy `log g` (boundary state-dependent) are **dropped**.
At the balanced bipartition (`ℓ = L/2`) the OBC value is *half* of the
PBC value (one entanglement cut vs two).

Reference: Calabrese–Cardy J. Stat. Mech. P06002 (2004) eq. (3.16);
J. Phys. A 42, 504005 (2009) eq. (30).  Affleck–Ludwig boundary
entropy: Affleck–Ludwig, Phys. Rev. Lett. 67, 161 (1991).
"""
function fetch(
    model::Universality{C}, ::VonNeumannEntropy, ::OBC; ℓ::Real, L::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("Cardy OBC: L must be > 0; got L=$L."))
    0 ≤ ℓ ≤ L ||
        throw(ArgumentError("Cardy OBC: ℓ must satisfy 0 ≤ ℓ ≤ L; got ℓ=$ℓ, L=$L."))
    c = _cardy_central_charge(model; kwargs...)
    (ℓ == 0 || ℓ == L) && return -Inf
    s = sin(π * ℓ / L)
    if s ≤ 0
        return -Inf
    end
    return (c / 6) * log((2 * L / π) * s)
end

"""
    fetch(::Universality{C}, ::VonNeumannEntropy, ::Infinite; ℓ::Real, kwargs...)
        -> Float64

Calabrese–Cardy von Neumann entanglement entropy in the thermodynamic
limit (`L → ∞`) of a 1+1D CFT, with PBC scaling — i.e. two
entanglement cuts in an infinite chain:

    S(ℓ) = (c/3) log ℓ                                   (Infinite)

The non-universal additive constant is dropped.  This is the standard
"infinite-chain" reference used to extract the central charge from
finite-size lattice data.

Reference: Calabrese–Cardy J. Stat. Mech. P06002 (2004) eq. (3.13).
"""
function fetch(
    model::Universality{C}, ::VonNeumannEntropy, ::Infinite; ℓ::Real, kwargs...
) where {C}
    ℓ > 0 || throw(ArgumentError("Cardy Infinite: ℓ must be > 0; got ℓ=$ℓ."))
    c = _cardy_central_charge(model; kwargs...)
    return (c / 3) * log(ℓ)
end

# ─── Rényi entropy ──────────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::RenyiEntropy, ::PBC; ℓ::Real, L::Real, kwargs...)
        -> Float64

Calabrese–Cardy Rényi-α entanglement entropy on a periodic chain.
Same closed form as the von Neumann case with the substitution

    c -> c · (1 + 1/α) / 2,

i.e.

    S_α(ℓ, L) = ((c/6)(1 + 1/α)) · log[(L/π) sin(πℓ/L)]

Reduces to the von Neumann result at `α = 1` (which is excluded here
because `RenyiEntropy(1)` throws — use `VonNeumannEntropy()` instead).

Reference: Calabrese–Cardy J. Stat. Mech. P06002 (2004) eq. (3.12).
"""
function fetch(
    model::Universality{C}, q::RenyiEntropy, ::PBC; ℓ::Real, L::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("Cardy PBC Rényi: L must be > 0; got L=$L."))
    0 ≤ ℓ ≤ L ||
        throw(ArgumentError("Cardy PBC Rényi: ℓ must satisfy 0 ≤ ℓ ≤ L; got ℓ=$ℓ, L=$L."))
    c = _cardy_central_charge(model; kwargs...)
    c_eff = _cardy_renyi_c(c, q.α)
    (ℓ == 0 || ℓ == L) && return -Inf
    s = sin(π * ℓ / L)
    if s ≤ 0
        return -Inf
    end
    return (c_eff / 3) * log((L / π) * s)
end

"""
    fetch(::Universality{C}, ::RenyiEntropy, ::OBC; ℓ::Real, L::Real, kwargs...)
        -> Float64

Calabrese–Cardy Rényi-α entanglement entropy on an open chain.
Same closed form as the von Neumann OBC case with the substitution

    c -> c · (1 + 1/α) / 2,

so

    S_α(ℓ, L) = ((c/12)(1 + 1/α)) · log[(2L/π) sin(πℓ/L)].

The non-universal `c'_1` and Affleck–Ludwig `log g` are dropped.

Reference: Calabrese–Cardy J. Phys. A 42, 504005 (2009) eq. (30).
"""
function fetch(
    model::Universality{C}, q::RenyiEntropy, ::OBC; ℓ::Real, L::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("Cardy OBC Rényi: L must be > 0; got L=$L."))
    0 ≤ ℓ ≤ L ||
        throw(ArgumentError("Cardy OBC Rényi: ℓ must satisfy 0 ≤ ℓ ≤ L; got ℓ=$ℓ, L=$L."))
    c = _cardy_central_charge(model; kwargs...)
    c_eff = _cardy_renyi_c(c, q.α)
    (ℓ == 0 || ℓ == L) && return -Inf
    s = sin(π * ℓ / L)
    if s ≤ 0
        return -Inf
    end
    return (c_eff / 6) * log((2 * L / π) * s)
end

"""
    fetch(::Universality{C}, ::RenyiEntropy, ::Infinite; ℓ::Real, kwargs...)
        -> Float64

Calabrese–Cardy Rényi-α entanglement entropy in the thermodynamic
limit, `L → ∞`:

    S_α(ℓ) = ((c/6)(1 + 1/α)) · log ℓ.

Reduces to the von Neumann `(c/3) log ℓ` at `α = 1` after the
substitution `c -> c · (1 + 1/α) / 2`.
"""
function fetch(
    model::Universality{C}, q::RenyiEntropy, ::Infinite; ℓ::Real, kwargs...
) where {C}
    ℓ > 0 || throw(ArgumentError("Cardy Infinite Rényi: ℓ must be > 0; got ℓ=$ℓ."))
    c = _cardy_central_charge(model; kwargs...)
    c_eff = _cardy_renyi_c(c, q.α)
    return (c_eff / 3) * log(ℓ)
end
