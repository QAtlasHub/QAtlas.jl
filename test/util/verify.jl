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
using Test: @test, @test_broken

const _VERIFY_ROUTES = (
    :ed_finite_size,        # exact diagonalisation, finite-N → closed form
    :sum_rule,              # an independent analytic identity / sum rule
    :delegation_invariant,  # model X at p ≡ model Y (different derivation)
    :limiting_case,         # a known independent value at a special point
    :literature_value,      # a published numeric (DMRG/MC) cross-check
    :second_closed_form,    # an independent closed form, different derivation
    :lieb_square_ice,       # Lieb 1967a square-ice closed form (SixVertex disordered diag)
    :lieb_ferroelectric,    # Lieb 1967c frozen-GS FE closed form (SixVertex Δ>1)
    :single_root_specialisation,  # specialised single-positive-root closed form (e.g. SU(2)_k Verlinde S₀₀)
    :multi_root_product,    # multi-positive-root product closed form (e.g. SU(N≥3)_k Verlinde S₀₀)
)

# review B1: only these routes are mechanically independent of `src`.
const _STRUCTURAL_ROUTES = (:ed_finite_size, :second_closed_form, :literature_value)

# Relative tolerance for treating a Complex value with negligible imaginary
# part as real in JSON emit (e.g. correlators of Hermitian operators that
# carry a few ULPs of round-off in the imaginary direction).
const _IMAG_NOISE_RTOL = 1e-9

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
    r = if x isa Complex
        (abs(imag(x)) <= _IMAG_NOISE_RTOL * max(1.0, abs(real(x))) ? real(x) : NaN)
    else
        x
    end
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
    return string(
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
           refs, reliability=:high, fetch_kw=(;), at=nothing,
           expected_fail=false, subject_extract=nothing) -> subject

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
    expected_fail::Bool=false,
    subject_extract::Union{Nothing,Function}=nothing,
)
    route in _VERIFY_ROUTES ||
        error("verify: route must be one of $(_VERIFY_ROUTES); got $(repr(route))")

    # ── the ONLY src touch-point: subject is fetched, never re-typed ──
    # subject_extract (optional) projects a non-scalar fetched value
    # (Vector, NamedTuple, container) to a single Float64 so verify()
    # can pin a specific component — e.g. E8Spectrum[3] or
    # CriticalExponents.β.
    raw = QAtlas.fetch(model, quantity, bc; fetch_kw...)
    subject = subject_extract === nothing ? raw : subject_extract(raw)

    ind =
        independent isa AbstractVector ? collect(float.(independent)) : [float(independent)]
    best = last(ind)                  # convergence: largest-N / final value
    abserr = abs(best - float(subject))

    if expected_fail
        # bug-surfacing card: @test_broken expects this to fail and will
        # alert when src is fixed (the test passes despite being marked
        # broken, prompting promotion back to @test).
        @test_broken isapprox(best, float(subject); atol=agree_within)
    else
        @test isapprox(best, float(subject); atol=agree_within)
    end

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
            "\"expected_fail\":",
            (expected_fail ? "true" : "false"),
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
            return println(io, card)
        end
        @info "verify: emitted v2 card" path = abspath(_ef) status = status
    end
    return subject
end

# ──────────────────────────────────────────────────────────────────────
# Shared v2 card emitter for the bound / approx verification variants.
# verify() above keeps its own inlined emit (its signature + emit are
# frozen for the 195 migrated cards); these siblings funnel through one
# helper so the JSONL schema stays consistent across all three. The
# `extra_fields` argument is a pre-formatted run of `"key":value,` pairs
# (with a trailing comma) spliced verbatim before the common tail.
# ──────────────────────────────────────────────────────────────────────

function _emit_card2(
    model,
    quantity,
    bc;
    route::Symbol,
    status::AbstractString,
    subject::Real,
    reliability::Symbol,
    refs::AbstractVector{<:AbstractString},
    at,
    independence::AbstractString,
    discriminant::AbstractString,
    extra_fields::AbstractString,
)
    get(ENV, "QATLAS_EMIT", "0") == "1" || return nothing
    outdir = get(ENV, "QATLAS_CIOUT_DIR", joinpath(@__DIR__, "..", ".ci-out"))
    mkpath(outdir)
    tf = get(ENV, "QATLAS_TEST_FILES", "")
    sid = isempty(tf) ? "all" : replace(string(hash(tf); base=16), '/' => '-')
    atv = at === nothing ? String[] : [string(a) for a in at]
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
        _json_str(independence),
        ",",
        "\"discriminant\":",
        _json_str(discriminant),
        ",",
        "\"status\":",
        _json_str(status),
        ",",
        extra_fields,
        "\"subject\":",
        _json_num(subject),
        ",",
        "\"at\":",
        _json_arr(atv),
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
        return println(io, card)
    end
    @info "verify variant: emitted v2 card" path = abspath(_ef) status = status route =
        route
    return nothing
end

# ──────────────────────────────────────────────────────────────────────
# verify_bound — one-sided inequality cards (registry status :bound).
# ──────────────────────────────────────────────────────────────────────

const _BOUND_ROUTES = (
    :ed_operator_bound,    # an ED-measured quantity stays within the fetched bound
    :variational_state,    # a trial state's expectation one-sidedly bounds the fetched value
    :dispersion_velocity,  # max group velocity of the dispersion saturates a Lieb-Robinson cone
    :saturating_constant,  # a universal constant attained (equality) at the optimal state
)

# ED-/trial-state-backed bound checks are structurally independent of the
# src formula; a saturating universal constant is an asserted optimum.
function _bound_independence(route::Symbol)
    if route in (:ed_operator_bound, :variational_state, :dispersion_velocity)
        return ("structural", "ed:bound-witness")
    else
        return ("asserted", "asserted:" * string(route))
    end
end

"""
    verify_bound(model, quantity, bc;
                 route, measured, relation,
                 refs, slack=0.0, saturating=false,
                 reliability=:high, fetch_kw=(;), at=nothing,
                 subject_extract=nothing) -> subject

One-sided counterpart of [`verify`](@ref) for `status=:bound` claims. The
fetched value `subject = fetch(model, quantity, bc; fetch_kw...)` is the
theoretical bound (RHS); `measured` is an INDEPENDENT witness (scalar or
vector — e.g. an ED expectation, or a trial-state energy) that must stay
on the correct side:

  * `relation=:leq` asserts `maximum(measured) <= subject + slack`
  * `relation=:geq` asserts `minimum(measured) >= subject - slack`

(the extremum is the tightest point of the sample). `saturating=true`
additionally asserts equality at that extremum — for universal constants
attained at the optimum (Tsirelson, Page average). The `<=`/`>=` direction
lives here on the card, never in the registry. Emits a v2 card with
`status` one of `pass` / `violated` / `saturated`.
"""
function verify_bound(
    model,
    quantity,
    bc;
    route::Symbol,
    measured,
    relation::Symbol,
    refs::AbstractVector{<:AbstractString},
    slack::Real=0.0,
    saturating::Bool=false,
    reliability::Symbol=:high,
    fetch_kw::NamedTuple=(;),
    at=nothing,
    subject_extract::Union{Nothing,Function}=nothing,
)
    route in _BOUND_ROUTES ||
        error("verify_bound: route must be one of $(_BOUND_ROUTES); got $(repr(route))")
    relation in (:leq, :geq) ||
        error("verify_bound: relation must be :leq or :geq; got $(repr(relation))")

    raw = QAtlas.fetch(model, quantity, bc; fetch_kw...)
    subject = subject_extract === nothing ? raw : subject_extract(raw)
    bound = float(subject)

    meas = measured isa AbstractVector ? collect(float.(measured)) : [float(measured)]
    worst = relation === :leq ? maximum(meas) : minimum(meas)

    satisfied = relation === :leq ? (worst <= bound + slack) : (worst >= bound - slack)
    @test satisfied
    saturated = isapprox(worst, bound; atol=max(slack, 1e-12))
    if saturating
        @test saturated
    end

    indep, disc = _bound_independence(route)
    status = !satisfied ? "violated" : (saturated ? "saturated" : "pass")
    extra = string(
        "\"relation\":",
        _json_str(relation),
        ",",
        "\"measured\":",
        _json_arr(meas),
        ",",
        "\"bound\":",
        _json_num(bound),
        ",",
        "\"slack\":",
        _json_num(slack),
        ",",
        "\"saturating\":",
        (saturating ? "true" : "false"),
        ",",
    )
    _emit_card2(
        model,
        quantity,
        bc;
        route=route,
        status=status,
        subject=bound,
        reliability=reliability,
        refs=refs,
        at=at,
        independence=indep,
        discriminant=disc,
        extra_fields=extra,
    )
    return subject
end

# ──────────────────────────────────────────────────────────────────────
# verify_approx — domain-limited approximation cards (registry status :approx).
# ──────────────────────────────────────────────────────────────────────

const _APPROX_ROUTES = (
    :high_temperature,   # high-T (small-β) expansion
    :low_temperature,    # low-T (large-β) expansion
    :large_n,            # large-N / mean-field expansion
    :perturbative,       # perturbation series in a small parameter
)

"""
    verify_approx(model, quantity, bc;
                  route, reference, agree_within, valid_domain, error_order,
                  refs, reliability=:medium, fetch_kw=(;), at=nothing,
                  subject_extract=nothing) -> subject

Approximation counterpart of [`verify`](@ref) for `status=:approx` claims.
`subject = fetch(model, quantity, bc; fetch_kw...)` is an approximation
(e.g. a high-T expansion) valid on a human-readable `valid_domain`
(e.g. `"betaJ << 1"`) with a known leading `error_order`
(e.g. `"O((betaJ)^2)"`). Asserts `subject ≈ reference` within
`agree_within` INSIDE the domain — the caller picks `fetch_kw` to sit
in-domain, and `reference` is the exact value there (ED or a tighter
formula). The domain and error order are recorded structurally so the
atlas can show *where* and *how well* the approximation holds.
"""
function verify_approx(
    model,
    quantity,
    bc;
    route::Symbol,
    reference,
    agree_within::Real,
    valid_domain::AbstractString,
    error_order::AbstractString,
    refs::AbstractVector{<:AbstractString},
    reliability::Symbol=:medium,
    fetch_kw::NamedTuple=(;),
    at=nothing,
    subject_extract::Union{Nothing,Function}=nothing,
)
    route in _APPROX_ROUTES ||
        error("verify_approx: route must be one of $(_APPROX_ROUTES); got $(repr(route))")

    raw = QAtlas.fetch(model, quantity, bc; fetch_kw...)
    subject = subject_extract === nothing ? raw : subject_extract(raw)
    approx = float(subject)
    ref = float(reference)

    ok = isapprox(approx, ref; atol=agree_within)
    @test ok

    status = ok ? "pass" : "fail"
    extra = string(
        "\"valid_domain\":",
        _json_str(valid_domain),
        ",",
        "\"error_order\":",
        _json_str(error_order),
        ",",
        "\"reference\":",
        _json_num(ref),
        ",",
        "\"agree_within\":",
        _json_num(agree_within),
        ",",
    )
    _emit_card2(
        model,
        quantity,
        bc;
        route=route,
        status=status,
        subject=approx,
        reliability=reliability,
        refs=refs,
        at=at,
        independence="asserted",
        discriminant="approx:" * string(route),
        extra_fields=extra,
    )
    return subject
end
