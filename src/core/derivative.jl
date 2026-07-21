# core/derivative.jl — supplying DERIVED inputs to a constraint edge.
#
# Most of AbstractQAtlas's relations that QAtlas could already reach are stated
# with a supplied derivative: `EntropyResponse(S, dF_dT)`, `GibbsHelmholtz(U,
# dβF_dβ)`, `SpecificHeatFromEntropy(C, dS_dT, T)`, `MagnetizationResponse(M,
# dF_dh)`, … The quantity slots are fetchable; the derivative slot is not, so
# the generator had nothing to put there and those relations stayed dead.  This
# file supplies them.
#
# It is also the deletion criterion `identity_registry.jl` names for the
# hand-written harness: "(a) the derivative-form identities (c_v, m_x = -∂f/∂h)
# are generatable".  The same ForwardDiff-through-`fetch` pattern is currently
# copy-pasted across at least four model test files (Kitaev1D, KitaevHoneycomb,
# XXZ1D, S1Heisenberg1D) — one declarative supplier replaces all of them.
#
# ── BACKENDS ──────────────────────────────────────────────────────────
# The differentiation method is a dispatchable choice, not a hard-wired one,
# because none of them works everywhere:
#
#   * `ForwardDiffBackend` — forward-mode AD.  Exact to machine precision, and
#     it only needs the fetch to be generic in its argument type; 273 of the
#     304 `beta::` annotations in src/ are already `::Real`, and the existing
#     hand-written tests prove it works through dense ED and quadrature.  This
#     is the right tool here: the derivatives are all scalar → scalar.
#   * `ZygoteBackend` — reverse-mode AD.  Offered because it is asked for, but
#     expect it to fail on more hubs than ForwardDiff: it needs an adjoint for
#     everything it traverses (dense ED, quadrature, the NLIE Newton solves),
#     and for a scalar → scalar derivative reverse mode has no advantage over
#     forward mode anyway.  Use it as a cross-check, not as the default.
#   * `FiniteDifference` — always available, no dependency, works through code
#     that is not differentiable at all.  The catch is TRUNCATION ERROR: a
#     central difference agrees with the analytic entropy to ~1e-5 relative on
#     CurieWeissIsing (measured), so an edge supplied this way must not be
#     given an AD-grade tolerance or it becomes a false-failure generator.
#     `default_rtol` below carries that difference so a caller cannot forget it.
#
# AD backends live in package EXTENSIONS (`ext/QAtlas*Ext.jl`), so neither
# ForwardDiff nor Zygote is a hard dependency of the atlas — they light up only
# when the user has loaded them.  With neither loaded, everything still works
# through `FiniteDifference`.

"""
    AbstractDiffBackend

How a derived input is obtained.  Concrete backends: [`FiniteDifference`](@ref)
(always available), `ForwardDiffBackend` and `ZygoteBackend` (package
extensions — see the file header for why they are not hard dependencies).
"""
abstract type AbstractDiffBackend end

"""
    FiniteDifference(; h = 1e-4)

Central-difference backend, `(f(x+h) − f(x−h)) / 2h`.  Dependency-free and
works through non-differentiable code, at the cost of truncation error — see
[`default_rtol`](@ref).
"""
struct FiniteDifference <: AbstractDiffBackend
    h::Float64
    function FiniteDifference(; h::Real=1e-4)
        return if h > 0
            new(Float64(h))
        else
            throw(ArgumentError("FiniteDifference: h must be > 0"))
        end
    end
end

"""
    ForwardDiffBackend()

Forward-mode AD.  A no-op stub unless `ForwardDiff` is loaded, in which case
`ext/QAtlasForwardDiffExt.jl` gives it a `derivative` method.
"""
struct ForwardDiffBackend <: AbstractDiffBackend end

"""
    ZygoteBackend()

Reverse-mode AD.  A no-op stub unless `Zygote` is loaded
(`ext/QAtlasZygoteExt.jl`).  Expect a lower success rate than
[`ForwardDiffBackend`](@ref) on this atlas — see the file header.
"""
struct ZygoteBackend <: AbstractDiffBackend end

"""
    derivative(backend, f, x) -> Real

`df/dx` at `x` by `backend`.

An AD backend whose package is not loaded throws a message that names the
`using` needed, rather than silently degrading to a finite difference: a check
that quietly changed its accuracy class would be reported at the wrong
tolerance.
"""
function derivative(b::AbstractDiffBackend, f, x::Real)
    return error(
        "QAtlas.derivative: no method for $(nameof(typeof(b))). It is provided by a " *
        "package extension — add `using $(_backend_package(b))` to activate it, or " *
        "pass `FiniteDifference()` to differentiate without any AD dependency.",
    )
end

function derivative(b::FiniteDifference, f, x::Real)
    h = b.h * max(one(Float64), abs(Float64(x)))   # scale-aware step
    return (f(x + h) - f(x - h)) / (2h)
end

_backend_package(::ForwardDiffBackend) = "ForwardDiff"
_backend_package(::ZygoteBackend) = "Zygote"
_backend_package(b::AbstractDiffBackend) = string(nameof(typeof(b)))

"""
    backend_available(backend) -> Bool

Whether `backend` can actually differentiate right now — i.e. its extension is
loaded.  Lets a generator pick the best available backend instead of erroring,
and lets a test skip an AD-only assertion honestly.
"""
backend_available(::FiniteDifference) = true
backend_available(::AbstractDiffBackend) = false

"""
    default_rtol(backend) -> Float64

The tolerance an edge should use when its derived input came from `backend`.

Forward-mode AD is exact to round-off **for a fetch that is a closed form, a
quadrature or a dense diagonalization**, so an AD-supplied identity there can be
held to the same `1e-6` the hand-written model tests use.  A central difference
carries truncation error — measured at `5.4e-5` relative against the analytic
entropy on CurieWeissIsing — so holding an FD-supplied check to `1e-6` would
manufacture failures that say nothing about the physics.

!!! warning "A tolerance cannot rescue a wrong derivative"
    AD is not uniformly the better answer.  Through an **iterative solve** it
    differentiates the iteration rather than the solution, and the derivative
    can be converged nowhere near as well as the value.  Measured on
    CurieWeissIsing's ordered phase, where the magnetization solves
    `m = tanh(βJm)`:

    | T | fetched `C` | AD `T·dS/dT` | FD `T·dS/dT` |
    |---|---|---|---|
    | 0.5 | 0.36594804 | **−0.01267792** | 0.36594804 |
    | 1/3 | 0.09345921 | **−0.00438224** | 0.09345921 |

    AD gets the SIGN wrong while FD is exact to eight digits.  No choice of
    `rtol` distinguishes that from a physics failure, which is why
    [`derivative_agreement`](@ref) exists: the defence is cross-checking two
    backends, not tightening a tolerance.
"""
default_rtol(::FiniteDifference) = 1e-3
default_rtol(::AbstractDiffBackend) = 1e-6

"""
    preferred_backend(; allow_fd = true) -> AbstractDiffBackend

The most accurate backend currently loaded: ForwardDiff, else Zygote, else a
[`FiniteDifference`](@ref).  With `allow_fd = false` it throws instead of
falling back, for a caller that must not silently change accuracy class.

"Most accurate" is the ranking for the fetches this atlas is mostly made of.
It is NOT a guarantee: see the warning on [`default_rtol`](@ref) for a hub where
AD is qualitatively wrong and the finite difference is right.  A caller that
cares should use [`derivative_agreement`](@ref) rather than trusting the rank.
"""
function preferred_backend(; allow_fd::Bool=true)
    for b in (ForwardDiffBackend(), ZygoteBackend())
        backend_available(b) && return b
    end
    allow_fd && return FiniteDifference()
    return error(
        "QAtlas.preferred_backend: no AD extension is loaded (add `using ForwardDiff` " *
        "or `using Zygote`) and `allow_fd = false` forbids the finite-difference " *
        "fallback.",
    )
end

"""
    derivative_agreement(f, x; primary = preferred_backend()) -> (value, backend, trusted, detail)

`df/dx` from `primary`, **cross-checked against an independent differentiation
method**, because the failure mode that matters here is not imprecision but a
silently wrong answer.

Returns the primary value, the backend that produced it, whether the two methods
agree to the cross-checking method's own tolerance, and a human-readable reason
when they do not.

Why this and not a tighter tolerance: forward-mode AD through an iterative solve
differentiates the iteration, not the solution, and can come back with the wrong
sign (measured — see [`default_rtol`](@ref)).  A disagreement between two
methods is evidence about the METHOD, not about the physics, so a caller should
report it as "cannot evaluate here", never as a failed physical identity.  It
also catches a second case for free: at a phase transition the quantity is
non-differentiable, a central difference straddles the jump while AD takes a
one-sided limit, and the two disagree — which is exactly the right verdict.

With only one backend available there is nothing to cross-check against, so
`trusted` is `true` and `detail` says so; the caller keeps whatever confidence
the single method deserves.
"""
function derivative_agreement(f, x::Real; primary::AbstractDiffBackend=preferred_backend())
    value = derivative(primary, f, x)
    primary isa FiniteDifference && return (value, primary, true, "single backend (FD)")
    cross = FiniteDifference()
    local other
    try
        other = derivative(cross, f, x)
    catch err
        return (value, primary, true, "cross-check unavailable: $(sprint(showerror, err))")
    end
    scale = max(abs(value), abs(other), eps())
    rel = abs(value - other) / scale
    rel ≤ default_rtol(cross) && return (value, primary, true, "agrees with FD")
    return (
        value,
        primary,
        false,
        "derivative not reproducible across methods: " *
        "$(nameof(typeof(primary)))=$(value) vs FiniteDifference=$(other) " *
        "(rel $(round(rel; sigdigits=3))). Typical causes: an iterative solve AD " *
        "differentiates through, or a non-differentiable point such as a phase " *
        "transition.",
    )
end

"""
    _primal(x) -> Real

The value of `x` with any AD tag stripped; the identity for ordinary numbers.
The ForwardDiff extension specializes it for `Dual`.

This is what lets a fetch built on a NON-DIFFERENTIABLE root-finder still hand
back an exact derivative: solve on primals, then re-attach the derivative
analytically (see `_curie_weiss_solve_m`).  Bisection is the motivating case —
it is numerically excellent and its AD derivative is meaningless, because the
β-dependence flows through the bracket endpoints and the comparison branches
rather than through the equation being solved.
"""
_primal(x::Real) = x
