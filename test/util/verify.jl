# test/util/verify.jl — black-box verification cards (WHY-correct plane).
#
# A `verify(...)` call cross-checks a src closed-form value against an
# INDEPENDENT route and (on push:main) emits one JSONL "verification
# card" to test/.ci-out/evidence-<sid>.jsonl.  The `ci-evidence` orphan
# branch accumulates the cards; a future doc/forum generator renders a
# per-(model,quantity) hub note whose confidence grows with the number
# and precision of independent corroborating cards (Zettelkasten-style).
#
# Anti-circularity is STRUCTURAL, not a convention:
#   * the value under test (`subject`) is ALWAYS obtained by the helper
#     calling `QAtlas.fetch(model, quantity, bc; fetch_kw...)`.  There is
#     NO argument that lets a test re-type the src formula.
#   * the test only supplies an `independent` number/array computed by a
#     route that does not look at src internals (ED, sum rule, …).
#   * `route` must be one of a fixed vocabulary of *independent* routes;
#     "retype the formula" is simply not expressible.
#
# Profiles: PR/local runs assert only.  push:main (QATLAS_EMIT=1, full
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
function _json_arr(xs)
    return "[" * join((x isa Real ? string(x) : _json_str(x) for x in xs), ",") * "]"
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
by the caller — so the card can never be circular.
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
        card = string(
            "{",
            "\"hub\":",
            _json_str(_verify_hub(model, quantity, bc)),
            ",",
            "\"route\":",
            _json_str(route),
            ",",
            "\"subject\":",
            float(subject),
            ",",
            "\"independent\":",
            _json_arr(ind),
            ",",
            "\"at\":",
            _json_arr(atv),
            ",",
            "\"abserr\":",
            round(abserr; sigdigits=4),
            ",",
            "\"atol\":",
            float(agree_within),
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
            "\"date\":",
            _json_str(string(Dates.today())),
            "}",
        )
        open(joinpath(outdir, "evidence-$(sid).jsonl"), "a") do io
            println(io, card)
        end
    end
    return subject
end
