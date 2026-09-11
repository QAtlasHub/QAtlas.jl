# core/disorder.jl: quenched disorder as a decoration, not a new model.
#
# A random chain is a chain whose couplings are random, so disorder attaches to
# a model rather than replacing it. `Disordered(clean, families)` names which of
# the clean model's FIELDS are random and with what distribution, so any model
# with named couplings can carry it.
#
# What a family must answer is the same for every model, which is why it lives
# here and not beside one of them.

"""
    DisorderFamily

A distribution of dimensionless couplings `λ > 0`.  A concrete family answers
[`log_moment`](@ref), [`mean_log`](@ref), [`var_log`](@ref) and
[`moment_floor`](@ref); nothing else is assumed about it.
"""
abstract type DisorderFamily end
export DisorderFamily

"""
    log_moment(f::DisorderFamily, s::Real) -> Float64

`ln E[λ^s]`, in logs so a residual built from it stays free of cancellation
where the terms are all `O(s)`.
"""
function log_moment end
export log_moment

"""    mean_log(f::DisorderFamily) -> Float64. `E[ln λ]`."""
function mean_log end
export mean_log

"""    var_log(f::DisorderFamily) -> Float64. `var[ln λ]`."""
function var_log end
export var_log

"""
    moment_floor(f::DisorderFamily) -> Float64

`inf{s : E[λ^s] < ∞}`, or `-Inf` where every moment exists.
"""
function moment_floor end
export moment_floor

# A family missing one of the four is otherwise a raw MethodError at whichever
# call path happens to reach it first, which need not be the first one tried.
for _f in (:log_moment, :mean_log, :var_log, :moment_floor)
    @eval function $(_f)(f::DisorderFamily, args...; kwargs...)
        return error(
            "QAtlas.$($(QuoteNode(_f))): not defined for $(nameof(typeof(f))). A " *
            "`DisorderFamily` must implement all four of `log_moment`, `mean_log`, " *
            "`var_log` and `moment_floor`. `mean_log` and `var_log` are the first " *
            "two derivatives of `log_moment` at s = 0, so they must agree with it.",
        )
    end
end

"""
    PowerLawDisorder(D) <: DisorderFamily

`P(λ) = D⁻¹ λ^{−1+1/D}` on `(0, 1]`, the strong-disorder RG's generic family
([`IgloiMonthus2005`](@cite) §A.1): `D² = var(ln λ)`, and `D = 1` is uniform on
`[0, 1]`.  `E[λ^s] = 1/(1 + D s)`, so moments below `s = −1/D` do not exist.
"""
struct PowerLawDisorder <: DisorderFamily
    D::Float64
    function PowerLawDisorder(D::Real)
        D > 0 || throw(ArgumentError("PowerLawDisorder: D must be > 0; got $D"))
        return new(Float64(D))
    end
end
export PowerLawDisorder

log_moment(f::PowerLawDisorder, s::Real) = -log1p(f.D * s)
mean_log(f::PowerLawDisorder) = -f.D
var_log(f::PowerLawDisorder) = f.D^2
moment_floor(f::PowerLawDisorder) = -1 / f.D

"""
    BinaryDisorder(κ) <: DisorderFamily

Two couplings, `λ ∈ {1, κ}` with equal probability, `0 < κ ≤ 1`.  Bounded away
from zero, so `E[λ^s]` exists for EVERY `s`: the lower bound on `z` that
[`PowerLawDisorder`](@ref) imposes is a property of that family and not of the
physics that reads it.
"""
struct BinaryDisorder <: DisorderFamily
    κ::Float64
    function BinaryDisorder(κ::Real)
        0 < κ <= 1 || throw(ArgumentError("BinaryDisorder: need 0 < κ ≤ 1; got $κ"))
        return new(Float64(κ))
    end
end
export BinaryDisorder

# log((1 + κ^s)/2), guarded for the s → ±∞ tails a root search walks through.
function log_moment(f::BinaryDisorder, s::Real)
    x = s * log(f.κ)
    return (x > 700 ? x : log1p(exp(x))) - log(2)
end
mean_log(f::BinaryDisorder) = log(f.κ) / 2
var_log(f::BinaryDisorder) = log(f.κ)^2 / 4
moment_floor(::BinaryDisorder) = -Inf

"""
    Disordered(clean::AbstractQAtlasModel; couplings...)

`clean` with the named couplings drawn at random: field `c` of `clean` becomes
`clean.c * λ`, with `λ` from the family given for `c`.  The clean value is the
SCALE, not the coupling.

```julia
Disordered(TFIM(; J=1.0, h=1.0); J=PowerLawDisorder(1.0), h=PowerLawDisorder(1.0))
```

Every name must be a field of `clean`, which is checked: a misspelling is the
one way to ask for disorder and silently not get it.

A `Disordered` model does not inherit the clean model's answers: it is not a
subtype of `M`, so a clean `fetch` cannot dispatch on it and asking for one is a
`MethodError` rather than a clean value. That is the type system doing it; no
blanket refusal is registered here, because one would be ambiguous with the
generic per-quantity dispatchers that already exist.
"""
struct Disordered{M<:AbstractQAtlasModel,F<:NamedTuple} <: AbstractQAtlasModel
    clean::M
    families::F
    function Disordered(clean::M, families::F) where {M<:AbstractQAtlasModel,F<:NamedTuple}
        isempty(families) &&
            throw(ArgumentError("Disordered: name at least one random coupling"))
        for c in keys(families)
            c in fieldnames(M) || throw(
                ArgumentError(
                    "Disordered: $M has no field `$c`; its fields are " *
                    "$(fieldnames(M)). A name that is not a coupling would be " *
                    "disorder that is never applied.",
                ),
            )
            families[c] isa DisorderFamily ||
                throw(ArgumentError("Disordered: `$c` must be a DisorderFamily"))
            # The field is the SCALE of `scale * λ`. Zero or non-finite makes the
            # coupling identically that whatever the family says, which is not a
            # random coupling but an absent one.
            scale = getfield(clean, c)
            scale isa Real && isfinite(scale) && !iszero(scale) || throw(
                ArgumentError(
                    "Disordered: `$c` is $(repr(scale)), which is the SCALE of a " *
                    "random coupling and must be a finite nonzero real. A zero " *
                    "scale is an absent coupling, not a random one.",
                ),
            )
        end
        return new{M,F}(clean, families)
    end
end
Disordered(clean::AbstractQAtlasModel; couplings...) = Disordered(clean, values(couplings))
export Disordered

"""
    clean_model(m::Disordered) -> AbstractQAtlasModel

The model `m` randomises.  Its field values are the SCALES of the random
couplings, not couplings themselves.
"""
clean_model(m::Disordered) = m.clean
export clean_model

"""
    disorder(m::Disordered, coupling::Symbol) -> DisorderFamily

The family randomising `coupling`.  Throws if that coupling is not random,
rather than returning a "no disorder" family, because a deterministic coupling
and an absent one are different statements.
"""
function disorder(m::Disordered, coupling::Symbol)
    haskey(m.families, coupling) || throw(
        ArgumentError(
            "disorder: `$coupling` is not random in this model; the random ones " *
            "are $(keys(m.families)).",
        ),
    )
    return m.families[coupling]
end
export disorder
