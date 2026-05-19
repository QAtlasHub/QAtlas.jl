# test/util/verify.jl — black-box verification cards (WHY-correct plane).
#
# A `verify(...)` call cross-checks a src closed-form value against an
# INDEPENDENT route and (on push:main / push:next) emits one JSONL
# "verification card" to test/.ci-out/evidence-<sid>.jsonl.  The
# `ci-evidence` orphan branch accumulates the cards.
#
# v2 (Discussion #379, Step 2 of the destructive migration): the emit
# path is hardened — NaN/Inf subject -> status:"divergent" (never an
# invalid-JSON raw NaN); each card carries the honest independence
# class (structural vs asserted, review B1) + a provenance discriminant;
# ED N-sweeps record a fitted convergence rate; schema_version = 2.
# The `fetch`/`@register` framework, the `verify(...)` call signature,
# and the `@test` assertion are UNCHANGED — the 9 migration PRs / 195
# cards keep their exact pass/fail behaviour; only the emitted evidence
# is richer and JSON-safe.
#
# Honest scope of the guarantee (never overstated): `subject` is always
# `fetch(...)` (structurally non-circular), but the *independence* of
# `independent` is only structural for ed_finite_size / second_closed_form
# / literature_value; sum_rule / delegation_invariant / limiting_case are
# labelled "asserted" and do NOT count as independent corroboration.
#
# Profiles: PR/local runs assert only.  push (QATLAS_EMIT=1, full
# profile) is the heavier computation and the only writer of cards.

using Dates: Dates
using Test: @test

const _VERIFY_ROUTES = (
    :ed_finite_size,        # exact diagonalisation, finite-N → closed form
    :sum_rule,              # an independent analytic identity / sum rule
    :delegation_invariant,  # model X at p ≡ model Y (different derivation)
    :limiting_case,         # a known independent value at a special point
    :literature_value,      # a published numeric (DMRG/MC) cross-check
    :second_closed_form,    # an independent closed form, different derivation
)

# review B1: only these routes are mechanically independent of `src`.
const _STRUCTURAL_ROUTES = (:ed_finite_size, :second_closed_form, :literature_value)

_verify_typename(x) = string(nameof(typeof(x)))

function _verify_hub(model, quantity, bc)
    return string(
        _verify_typename(model), "/", _verify_typename(quantity), "/", _verify_typename(bc)
    )
end

function _verify_commit()
    sha = get(ENV, "GITHUB_SHA", "")
    isempty(sha) || return sha
    try
        return strip(read(`git rev-parse --short HEAD`, String))
    catch
        return "unknown"
    end
end

# Profile-scaled N-sweep helper for independent ED routes.
function verify_profile_Ns(;
    fast=(6, 8), full=(6, 8, 10, 12), nightly=(6, 8, 10, 12, 14, 16)
)
    p = isdefined(Main, :QATLAS_TEST_PROFILE) ? Main.QATLAS_TEST_PROFILE : :fast
    p === :nightly && return nightly
    p === :full && return full
    return fast
end

_json_str(s) = '"' * replace(string(s), '\\' => "\\\\", '"' => "\\\"") * '"'

# NaN/Inf-safe numeric emitter — a non-finite number is NEVER written as a
# raw token (which is invalid JSON); it becomes null and the card status
# is set to "divergent" instead (original critique item 5).
function _json_num(x)
    # Complex-with-negligible-imag -> real part (correlators of
    # Hermitian operators are real up to round-off); genuinely
    # complex or non-finite -> null (never a raw "0.4 + 0.0im").
    r = x isa Complex ? (abs(imag(x)) <= 1e-9 * max(1.0, abs(real(x))) ? real(x) : NaN) : x
    v = float(r)
    return isfinite(v) ? string(v) : "null"
end

function _json_arr(xs)
    return "[" * join((x isa Number ? _json_num(x) : _json_str(x) for x in xs), ",") * "]"
end

# review B1: (independence-class, provenance discriminant) from the route.
function _v2_independence(route::Symbol, refs)
    structural = route in _STRUCTURAL_ROUTES
    disc = if route === :ed_finite_size
        "ed:dense-diagonalization"
    elseif route === :literature_value
        "lit:" * (isempty(refs) ? "?" : join(refs, " | "))
    elseif route === :second_closed_form
        "cf:" * (isempty(refs) ? "?" : join(refs, " | "))
    else
        "asserted:" * string(route)
    end
    return (structural ? "structural" : "asserted"), disc
end

# Best-effort convergence-as-evidence: parse a trailing integer N from each
# `at` label (e.g. "N=8" → 8), fit the log-log slope of |independent - subject|
# vs N.  Returns the fitted rate, or NaN when < 3 usable points.
function _fit_rate(ind, atv, subj)
    length(atv) == length(ind) || return NaN
    Ns = Float64[]
    rs = Float64[]
    for (a, v) in zip(atv, ind)
        m = match(r"(\d+)\s*$", string(a))
        m === nothing && continue
        N = parse(Float64, m.captures[1])
        r = abs(v - subj)
        (isfinite(N) && isfinite(r) && N > 0 && r > 0) || continue
        push!(Ns, log(N))
        push!(rs, log(r))
    end
    length(Ns) >= 3 || return NaN
    n = length(Ns)
    mx = sum(Ns) / n
    my = sum(rs) / n
    den = sum((Ns .- mx) .^ 2)
    den == 0 && return NaN
    return sum((Ns .- mx) .* (rs .- my)) / den
end

function _v2_env()
    string(
        "julia-",
        VERSION,
        "; runner=",
        get(ENV, "RUNNER_OS", "local"),
        "; run=",
        get(ENV, "GITHUB_RUN_ID", ""),
    )
end

"""
    verify(model, quantity, bc;
           route, independent, agree_within,
           refs, reliability=:high, fetch_kw=(;), at=nothing) -> subject

Black-box-verify the src value `fetch(model, quantity, bc; fetch_kw...)`
against an `independent` numeric (scalar or convergence vector) obtained
by `route` (one of `$(_VERIFY_ROUTES)`).  Runs `@test` (best/last
element within `agree_within`) and, under `QATLAS_EMIT=1`, appends one
JSONL verification card.  The subject is fetched here — never supplied
by the caller — so the card can never be circular.  Signature and the
`@test` are frozen; the emitted card is v2-hardened (see file header).
"""
function verify(
    model,
    quantity,
    bc;
    route::Symbol,
    independent,
    agree_within::Real,
    refs::AbstractVector{<:AbstractString},
    reliability::Symbol=:high,
    fetch_kw::NamedTuple=(;),
    at=nothing,
)
    route in _VERIFY_ROUTES ||
        error("verify: route must be one of $(_VERIFY_ROUTES); got $(repr(route))")

    # ── the ONLY src touch-point: subject is fetched, never re-typed ──
    subject = QAtlas.fetch(model, quantity, bc; fetch_kw...)

    ind =
        independent isa AbstractVector ? collect(float.(independent)) : [float(independent)]
    best = last(ind)                  # convergence: largest-N / final value
    abserr = abs(best - float(subject))

    @test isapprox(best, float(subject); atol=agree_within)

    if get(ENV, "QATLAS_EMIT", "0") == "1"
        outdir = get(ENV, "QATLAS_CIOUT_DIR", joinpath(@__DIR__, "..", ".ci-out"))
        mkpath(outdir)
        tf = get(ENV, "QATLAS_TEST_FILES", "")
        sid = isempty(tf) ? "all" : replace(string(hash(tf); base=16), '/' => '-')
        atv = at === nothing ? String[] : [string(a) for a in at]

        subj = float(subject)
        finite = isfinite(subj) && all(isfinite, ind)
        status = if !finite
            "divergent"
        elseif isapprox(best, subj; atol=agree_within)
            "pass"
        else
            "fail"
        end
        indep, disc = _v2_independence(route, refs)
        rate = finite ? _fit_rate(ind, atv, subj) : NaN

        card = string(
            "{",
            "\"schema_version\":2,",
            "\"hub\":",
            _json_str(_verify_hub(model, quantity, bc)),
            ",",
            "\"route\":",
            _json_str(route),
            ",",
            "\"mechanism\":",
            _json_str(route),
            ",",
            "\"independence\":",
            _json_str(indep),
            ",",
            "\"discriminant\":",
            _json_str(disc),
            ",",
            "\"status\":",
            _json_str(status),
            ",",
            "\"subject\":",
            _json_num(subj),
            ",",
            "\"independent\":",
            _json_arr(ind),
            ",",
            "\"at\":",
            _json_arr(atv),
            ",",
            "\"abserr\":",
            _json_num(finite ? round(abserr; sigdigits=6) : NaN),
            ",",
            "\"atol\":",
            _json_num(agree_within),
            ",",
            "\"convergence_rate\":",
            _json_num(rate),
            ",",
            "\"reliability\":",
            _json_str(reliability),
            ",",
            "\"refs\":",
            _json_arr(refs),
            ",",
            "\"commit\":",
            _json_str(_verify_commit()),
            ",",
            "\"env\":",
            _json_str(_v2_env()),
            ",",
            "\"date\":",
            _json_str(string(Dates.today())),
            "}",
        )
        _ef = joinpath(outdir, "evidence-$(sid).jsonl")
        open(_ef, "a") do io
            println(io, card)
        end
        @info "verify: emitted v2 card" path = abspath(_ef) status = status
    end
    return subject
end
