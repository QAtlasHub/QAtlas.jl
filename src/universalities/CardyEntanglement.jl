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
    fetch(::Universality{C}, ::VonNeumannEntropy, ::Infinite;
          ℓ::Real, beta::Real = Inf, kwargs...) -> Float64

Calabrese–Cardy von Neumann entanglement entropy in the thermodynamic
limit (`L → ∞`) of a 1+1D CFT, with PBC scaling — two entanglement
cuts in an infinite chain.

- `beta = Inf` (default): T = 0 ground state.

      S(ℓ) = (c/3) log ℓ                                   (Infinite)

- `0 < beta < Inf`: finite-temperature thermal state of the CFT.

      S(ℓ, β) = (c/3) log[(β/π) sinh(π ℓ / β)]             (Infinite, T > 0)

Both forms drop the non-universal additive constant.  The β → ∞ limit
of the finite-T form recovers the T = 0 expression because
`(β/π) sinh(π ℓ / β) → ℓ`.  Fermi velocity is normalised to unity
(lattice units); models computing thermal entanglement at a non-unit
sound velocity should pass `beta` already scaled by their velocity.

Reference: Calabrese–Cardy J. Stat. Mech. P06002 (2004) §4, eq. (3.13),
(3.16). Tracking: #580.
"""
function fetch(
    model::Universality{C},
    ::VonNeumannEntropy,
    ::Infinite;
    ℓ::Real,
    beta::Real=Inf,
    a::Real=1.0,
    kwargs...,
) where {C}
    ℓ > 0 || throw(ArgumentError("Cardy Infinite: ℓ must be > 0; got ℓ=$ℓ."))
    a > 0 ||
        throw(ArgumentError("Cardy Infinite: lattice spacing a must be > 0; got a=$a."))
    c = _cardy_central_charge(model; kwargs...)
    if isinf(beta)
        return (c / 3) * log(ℓ / a)
    else
        beta > 0 || throw(
            ArgumentError("Cardy Infinite finite-T: beta must be > 0; got beta=$beta.")
        )
        return (c / 3) * log((beta / (π * a)) * sinh(π * ℓ / beta))
    end
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
    model::Universality{C},
    q::RenyiEntropy,
    ::Infinite;
    ℓ::Real,
    beta::Real=Inf,
    a::Real=1.0,
    kwargs...,
) where {C}
    ℓ > 0 || throw(ArgumentError("Cardy Infinite Rényi: ℓ must be > 0; got ℓ=$ℓ."))
    a > 0 || throw(
        ArgumentError("Cardy Infinite Rényi: lattice spacing a must be > 0; got a=$a.")
    )
    c = _cardy_central_charge(model; kwargs...)
    c_eff = _cardy_renyi_c(c, q.α)
    if isinf(beta)
        return (c_eff / 3) * log(ℓ / a)
    else
        beta > 0 || throw(
            ArgumentError(
                "Cardy Infinite Rényi finite-T: beta must be > 0; got beta=$beta."
            ),
        )
        return (c_eff / 3) * log((beta / (π * a)) * sinh(π * ℓ / beta))
    end
end

# ─── Infinite bc forwarding for verify() integration ───────────────────────
#
# Universality{C}/CentralCharge is a property of an emergent CFT, not of a
# finite lattice with boundary conditions. The forwards below let verify()
# cards quote these literature values with the standard 3-arg signature
# (model, quantity, bc) used by every other CentralCharge hub.

function fetch(m::Universality{:Ising}, q::CentralCharge, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Potts3}, q::CentralCharge, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Potts4}, q::CentralCharge, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:XY}, q::CentralCharge, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Heisenberg}, q::CentralCharge, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end

# ─── Infinite-bc forwarders for RMT / Poisson quantities (verify integration)
function fetch(m::Universality{:RMT}, q::WignerSurmise, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
fetch(m::Universality{:RMT}, q::TracyWidom, ::Infinite; kwargs...) = fetch(m, q; kwargs...)
fetch(m::Universality{:RMT}, q::MeanRatio, ::Infinite; kwargs...) = fetch(m, q; kwargs...)
function fetch(m::Universality{:RMT}, q::SpectralFormFactor, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Poisson}, q::WignerSurmise, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Poisson}, q::MeanRatio, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end

# ─── Infinite-bc forwarders for CriticalExponents / GrowthExponents (verify integration)
function fetch(m::Universality{:Percolation}, q::CriticalExponents, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Potts3}, q::CriticalExponents, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Potts4}, q::CriticalExponents, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:KPZ}, q::GrowthExponents, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
fetch(m::MeanField, q::CriticalExponents, ::Infinite; kwargs...) = fetch(m, q; kwargs...)
fetch(m::Ising2D, q::CriticalExponents, ::Infinite; kwargs...) = fetch(m, q; kwargs...)
fetch(m::KPZ1D, q::CriticalExponents, ::Infinite; kwargs...) = fetch(m, q; kwargs...)

# ─────────────────────────────────────────────────────────────────────────────
# Mutual information for two adjacent intervals at Infinite (#580 Phase 2+)
#
# For two adjacent intervals of lengths ℓ_A and ℓ_B on an infinite
# 1+1D-CFT chain at T = 0, I(A:B) = S(A) + S(B) - S(A ∪ B) reduces to
#
#     I(A:B) = (c/3) log[ℓ_A · ℓ_B / (ℓ_A + ℓ_B)]   (T = 0)
#
# At β finite each entropy takes the CC sinh form.
#
# Reference: Calabrese-Cardy J. Stat. Mech. P06002 (2004); J. Phys. A 42,
# 504005 (2009). Tracking: #580.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::MutualInformation, ::Infinite;
          ℓ_A::Real, ℓ_B::Real, beta::Real = Inf, kwargs...) -> Float64

Mutual information of two adjacent intervals of lengths `ℓ_A` and
`ℓ_B` on an infinite 1+1D-CFT chain.  Closed form follows from the
single-interval Calabrese-Cardy result via
`I = S(A) + S(B) − S(A∪B)`:

    I(A:B) = (c/3) log[ℓ_A · ℓ_B / (ℓ_A + ℓ_B)]      (T = 0)

At finite β each single-interval entropy takes the sinh form.

Refs Calabrese-Cardy 2004 / 2009; issue #580.
"""
function fetch(
    model::Universality{C},
    ::MutualInformation,
    ::Infinite;
    ℓ_A::Real,
    ℓ_B::Real,
    beta::Real=Inf,
    kwargs...,
) where {C}
    ℓ_A > 0 || throw(
        ArgumentError("Cardy Infinite MutualInformation: ℓ_A must be > 0; got ℓ_A=$ℓ_A."),
    )
    ℓ_B > 0 || throw(
        ArgumentError("Cardy Infinite MutualInformation: ℓ_B must be > 0; got ℓ_B=$ℓ_B."),
    )
    S_A = fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ_A, beta=beta, kwargs...)
    S_B = fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ_B, beta=beta, kwargs...)
    S_AB = fetch(model, VonNeumannEntropy(), Infinite(); ℓ=ℓ_A + ℓ_B, beta=beta, kwargs...)
    return S_A + S_B - S_AB
end

# ─────────────────────────────────────────────────────────────────────────────
# Quench-dynamics: linear-growth slope of half-system entanglement (#580)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::EntanglementGrowthSlope, ::Infinite;
          v::Real, beta_eff::Real, kwargs...) -> Float64

Linear-growth slope `dS_A / dt` of the half-system entanglement entropy
after a global quench from a thermal-like initial state into a critical
post-quench Hamiltonian in the same universality class. Calabrese-Cardy
2005 predicts

    dS_A / dt = (π c v) / (3 beta_eff)            (t < L / (2 v))

where

- `c` is the central charge of the post-quench critical Hamiltonian
  (provided by the universality dispatch),
- `v` is the propagation velocity (Lieb-Robinson for free-fermion
  models; a model-dependent sound velocity in general),
- `beta_eff` is the effective inverse temperature of the generalised-
  Gibbs steady state, set by the initial state energy density.

Reference: Calabrese-Cardy J. Stat. Mech. P04010 (2005). Tracking #580.
"""
function fetch(
    model::Universality{C},
    ::EntanglementGrowthSlope,
    ::Infinite;
    v::Real,
    beta_eff::Real,
    kwargs...,
) where {C}
    v > 0 || throw(ArgumentError("EntanglementGrowthSlope: v must be > 0; got v=$v."))
    beta_eff > 0 || throw(
        ArgumentError(
            "EntanglementGrowthSlope: beta_eff must be > 0; got beta_eff=$beta_eff."
        ),
    )
    c = _cardy_central_charge(model; kwargs...)
    return π * c * v / (3 * beta_eff)
end

# ─────────────────────────────────────────────────────────────────────────────
# Cardy formula: asymptotic high-energy state-counting entropy (#580)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::CardyEntropy, ::Infinite;
          E::Real, kwargs...) -> Float64

Asymptotic high-energy entropy of a 1+1D CFT (Cardy 1986):

    S_Cardy(E) = 2 π sqrt(c E / 6),

with `c` supplied by the universality class. This is the log of the
microcanonical density of states at energy `E` on a cylinder of unit
circumference. Valid asymptotically at large `E`; at low `E` the
formula systematically underestimates the count.

Reference: Cardy, *Nucl. Phys. B* **270**, 186 (1986).
"""
function fetch(
    model::Universality{C}, ::CardyEntropy, ::Infinite; E::Real, kwargs...
) where {C}
    E >= 0 || throw(ArgumentError("CardyEntropy: E must be >= 0; got E=$E."))
    c = _cardy_central_charge(model; kwargs...)
    return 2 * π * sqrt(c * E / 6)
end

# ─────────────────────────────────────────────────────────────────────────────
# Conformal Casimir energy on a cylinder (#580)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::ConformalCasimirEnergy, ::Infinite;
          L::Real, kwargs...) -> Float64

Universal Casimir ground-state energy of a 1+1D CFT on a cylinder of
circumference `L` (Cardy 1986 / Blote-Cardy-Nightingale 1986 /
Affleck 1986):

    E_0(L) = -π c / (6 L).

The sign convention follows the original PRLs (`E_0 < 0` for unitary
CFTs with `c > 0`). Identified empirically by subtracting `L * e_∞`
from the lattice ground-state energy on a periodic chain of `L` sites
and rescaling by `L`.

Reference: Cardy *Nucl. Phys. B* **270**, 186 (1986); Blote-Cardy-
Nightingale *Phys. Rev. Lett.* **56**, 742 (1986); Affleck *Phys. Rev.
Lett.* **56**, 746 (1986).
"""
function fetch(
    model::Universality{C}, ::ConformalCasimirEnergy, ::Infinite; L::Real, kwargs...
) where {C}
    L > 0 || throw(ArgumentError("ConformalCasimirEnergy: L must be > 0; got L=$L."))
    c = _cardy_central_charge(model; kwargs...)
    return -π * c / (6 * L)
end

# ─────────────────────────────────────────────────────────────────────────────
# Logarithmic negativity for two adjacent intervals at Infinite (#580)
#
# Calabrese-Cardy-Tonni 2012 (DOI 10.1103/PhysRevLett.109.130502) found
# that the universal piece of E = log Tr |rho^{T_B}| for two adjacent
# intervals of lengths ell_A, ell_B on an infinite 1+1D-CFT chain at
# T = 0 is
#
#     E(ell_A, ell_B) = (c/4) log[ell_A * ell_B / (ell_A + ell_B)],
#
# i.e., the same geometric-mean log as the mutual-information universal
# formula (PR #587) with the prefactor c/3 replaced by c/4.
#
# Reference: P. Calabrese, J. Cardy, E. Tonni, Phys. Rev. Lett. 109,
# 130502 (2012). Tracking: #580.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::LogarithmicNegativity, ::Infinite;
          ℓ_A::Real, ℓ_B::Real, kwargs...) -> Float64

Logarithmic negativity of two adjacent intervals at T = 0 on an
infinite 1+1D-CFT chain (Calabrese-Cardy-Tonni 2012):

    E = (c/4) log[ℓ_A * ℓ_B / (ℓ_A + ℓ_B)]

Returns the universal log piece; non-universal additive constants are
dropped.
"""
function fetch(
    model::Universality{C},
    ::LogarithmicNegativity,
    ::Infinite;
    ℓ_A::Real,
    ℓ_B::Real,
    kwargs...,
) where {C}
    ℓ_A > 0 ||
        throw(ArgumentError("LogarithmicNegativity: ell_A must be > 0; got ell_A=$ℓ_A."))
    ℓ_B > 0 ||
        throw(ArgumentError("LogarithmicNegativity: ell_B must be > 0; got ell_B=$ℓ_B."))
    c = _cardy_central_charge(model; kwargs...)
    return (c / 4) * log(ℓ_A * ℓ_B / (ℓ_A + ℓ_B))
end

# ─────────────────────────────────────────────────────────────────────────────
# Affleck-Ludwig boundary entropy log g (#580)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:Ising}, ::BoundaryEntropy, ::Infinite;
          boundary_state::Symbol, kwargs...) -> Float64

Affleck-Ludwig universal boundary entropy `log g` for an Ising CFT
(M(3,4), c = 1/2) Cardy boundary state. From the modular S-matrix of
the Ising minimal model the universal `g_a = S_{0a} / sqrt(S_{00})`
values are

    g_1 (identity)  = 1/sqrt(2),   log g = -(1/2) log 2
    g_σ (spin)      = 1,           log g = 0
    g_ε (energy)    = 1/sqrt(2),   log g = -(1/2) log 2

Physical interpretation: `|σ⟩` is the "free" (unconstrained) boundary
and `|1⟩` and `|ε⟩` are the "fixed" boundaries (spin pinned up or
down). The free boundary has higher `g`, so free flows to fixed under
RG (g-theorem).

`boundary_state` selects the Cardy state, one of:
    `:identity` (≡ `:fixed_up`),
    `:sigma`    (≡ `:free`),
    `:epsilon`  (≡ `:fixed_down`).

Reference: Affleck-Ludwig PRL **67**, 161 (1991);
Cardy *Nucl. Phys. B* **324**, 581 (1989) for the Cardy state
construction. Tracking: #580.
"""
function fetch(
    ::Universality{:Ising}, ::BoundaryEntropy, ::Infinite; boundary_state::Symbol, kwargs...
)
    if boundary_state === :identity ||
        boundary_state === :fixed_up ||
        boundary_state === :epsilon ||
        boundary_state === :fixed_down
        return -log(2) / 2  # log(1/sqrt(2))
    elseif boundary_state === :sigma || boundary_state === :free
        return 0.0           # log(1)
    else
        throw(
            ArgumentError(
                "Universality(:Ising) BoundaryEntropy: boundary_state must be one " *
                "of :identity / :fixed_up / :sigma / :free / :epsilon / :fixed_down; " *
                "got \$(boundary_state).",
            ),
        )
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Page entropy for Haar-random pure states (#580)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:HaarRandom}, ::PageEntropy, ::Infinite;
          d_A::Integer, d_B::Integer, kwargs...) -> Float64

Page 1993 average entropy of a subsystem `A` for a Haar-random pure
state in `H_A ⊗ H_B` with `dim(H_A) = d_A`, `dim(H_B) = d_B`:

    <S_A> = sum_{k=n+1}^{m·n} 1/k - (m-1)/(2n)

where `m = min(d_A, d_B)` and `n = max(d_A, d_B)` (the formula is
invariant under A ↔ B exchange by purity of the global state).

Asymptotic limits:
- `m == n`:   `<S_A> ≈ log m - 1/2`  (nearly maximal entropy)
- `m << n`:   `<S_A> ≈ log m - m/(2n)`  (almost maximal volume law)
- `m >> n`:   `<S_A> ≈ log n - n/(2m)`  (same by A ↔ B symmetry)

Reference: D. N. Page, *Phys. Rev. Lett.* **71**, 1291 (1993),
DOI 10.1103/PhysRevLett.71.1291.
"""
function fetch(
    ::Universality{:HaarRandom},
    ::PageEntropy,
    ::Infinite;
    d_A::Integer,
    d_B::Integer,
    kwargs...,
)
    d_A >= 1 || throw(ArgumentError("PageEntropy: d_A must be >= 1; got d_A=$d_A."))
    d_B >= 1 || throw(ArgumentError("PageEntropy: d_B must be >= 1; got d_B=$d_B."))
    m, n = minmax(d_A, d_B)
    s = 0.0
    for k in (n + 1):(m * n)
        s += 1.0 / k
    end
    return s - (m - 1) / (2.0 * n)
end

# ─────────────────────────────────────────────────────────────────────────────
# Topological entanglement entropy via total quantum dimension (#580)
#
# Kitaev-Preskill 2006 (DOI 10.1103/PhysRevLett.96.110404) and Levin-Wen
# 2006 (DOI 10.1103/PhysRevLett.96.110405) showed that the subleading
# constant in the entanglement entropy of a topologically ordered
# 2D ground state is the universal piece
#
#     gamma = log D,      D = sqrt(sum_a d_a^2),
#
# where {d_a} are the quantum dimensions of the anyon types and D is
# the total quantum dimension. Examples:
#
#     ToricCode:      anyons {1, e, m, eps},  d_a = (1,1,1,1) -> D = 2  -> gamma = log 2
#     Ising anyon:    anyons {1, sigma, psi},  d_a = (1, sqrt 2, 1) -> D = 2 -> gamma = log 2
#     Fibonacci:      anyons {1, tau},  d_a = (1, phi),  phi = (1+sqrt 5)/2 -> D = sqrt(2+phi)
#
# Tracking: #580.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:TopologicalOrder}, ::TopologicalEntanglementEntropy,
          ::Infinite; quantum_dimensions::AbstractVector{<:Real}, kwargs...)
        -> Float64

Universal topological entanglement entropy of a 2D topologically
ordered ground state in the Kitaev-Preskill / Levin-Wen convention,

    gamma = log D,    D = sqrt(sum_a d_a^2),

with `quantum_dimensions = [d_a]` the vector of anyon quantum
dimensions. References: Kitaev-Preskill 2006 (PRL 96, 110404);
Levin-Wen 2006 (PRL 96, 110405).
"""
function fetch(
    ::Universality{:TopologicalOrder},
    ::TopologicalEntanglementEntropy,
    ::Infinite;
    quantum_dimensions::AbstractVector{<:Real},
    kwargs...,
)
    isempty(quantum_dimensions) && throw(
        ArgumentError(
            "TopologicalEntanglementEntropy: quantum_dimensions must be non-empty."
        ),
    )
    all(d -> d > 0, quantum_dimensions) || throw(
        ArgumentError(
            "TopologicalEntanglementEntropy: all quantum dimensions d_a must be > 0."
        ),
    )
    D_sq = sum(d -> d^2, quantum_dimensions)
    return 0.5 * log(D_sq)
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-length saturation of post-quench entanglement (#580)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{C}, ::EntanglementSaturationDensity, ::Infinite;
          beta_eff::Real, kwargs...) -> Float64

Long-time saturation of post-quench entanglement entropy per unit
length (Calabrese-Cardy 2005):

    S_A(infty) / L = pi c / (6 beta_eff),

with `c` provided by the universality class. Partner to
EntanglementGrowthSlope (PR #588): the linear-regime extrapolation
slope * (L / (2 v)) equals this saturation value, expressing the
crossover at the boundary of the light cone.

Reference: Calabrese-Cardy J. Stat. Mech. P04010 (2005). Tracking #580.
"""
function fetch(
    model::Universality{C},
    ::EntanglementSaturationDensity,
    ::Infinite;
    beta_eff::Real,
    kwargs...,
) where {C}
    beta_eff > 0 || throw(
        ArgumentError(
            "EntanglementSaturationDensity: beta_eff must be > 0; got beta_eff=$beta_eff.",
        ),
    )
    c = _cardy_central_charge(model; kwargs...)
    return π * c / (6 * beta_eff)
end

# -----------------------------------------------------------------------------
# CHSH bounds (#579 inequality framework Phase 1)
# -----------------------------------------------------------------------------

"""
    fetch(::Universality{:QuantumMechanics}, ::CHSHBound, ::Infinite;
          convention::Symbol = :quantum, kwargs...) -> Float64

Universal CHSH bounds for the Bell-inequality correlator S =
E(a,b) + E(a,bp) + E(ap,b) - E(ap,bp):

- :classical              -> 2          (CHSH 1969 local-realistic)
- :quantum (default)      -> 2 sqrt(2)  (Tsirelson 1980 quantum max)
- :no_signalling, :pr     -> 4          (Popescu-Rohrlich algebraic max)

The quantum bound 2 sqrt(2) is saturated by the maximally entangled
Bell state with appropriate measurement axes 0, pi/4, pi/2, 3 pi/4.

References: Clauser-Horne-Shimony-Holt PRL 23, 880 (1969); B. S.
Cirelson (Tsirelson), Lett. Math. Phys. 4, 93 (1980), DOI
10.1007/BF00417500 (metadata fetched via doiget). Tracking: #579.
"""
function fetch(
    ::Universality{:QuantumMechanics},
    ::CHSHBound,
    ::Infinite;
    convention::Symbol=:quantum,
    kwargs...,
)
    if convention === :quantum
        return 2 * sqrt(2)
    elseif convention === :classical
        return 2.0
    elseif convention === :no_signalling || convention === :pr
        return 4.0
    else
        throw(
            ArgumentError(
                "CHSHBound: convention must be one of :quantum / :classical / " *
                ":no_signalling / :pr; got :$(convention).",
            ),
        )
    end
end

# -----------------------------------------------------------------------------
# Thermal energy density (Affleck 1986; Bloete-Cardy-Nightingale 1986)
# -----------------------------------------------------------------------------

"""
    fetch(::Universality{C}, ::ThermalEnergyDensity, ::Infinite;
          beta::Real, kwargs...) where {C} -> Float64

Leading thermal energy density above the ground state for a
(1+1)D CFT with central charge ,

    e(T) - e_0 = pi c / (6 beta^2).

This is the universal complement of
[](@ref): modular invariance interchanges the
finite-size Casimir energy and the finite-temperature thermal energy
with identical coefficient .

Reference: I. Affleck *Phys. Rev. Lett.* **56**, 746 (1986);
Bloete-Cardy-Nightingale *Phys. Rev. Lett.* **56**, 742 (1986).
"""
function fetch(
    ::Universality{C}, ::ThermalEnergyDensity, ::Infinite; beta::Real, kwargs...
) where {C}
    beta > 0 || throw(
        DomainError(beta, "ThermalEnergyDensity requires beta > 0; got beta = $beta.")
    )
    c = _cardy_central_charge(Universality(C))
    return π * c / (6 * beta^2)
end
