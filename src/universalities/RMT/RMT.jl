# ─────────────────────────────────────────────────────────────────────────────
# RMT (Random Matrix Theory) universality class — Wigner-Dyson level
# statistics for β ∈ {1, 2, 4} (GOE, GUE, GSE).
#
# Quantities provided:
#   * WignerSurmise — closed-form nearest-neighbour spacing P_β(s)
#   * TracyWidom    — largest-eigenvalue distribution F_β(x), Phase 1
#                     via tabulated values + tail asymptotics
#   * MeanRatio     — Atas-Bogomolny-Giraud-Roux ⟨r⟩
#
# References:
#   M. L. Mehta, "Random Matrices", 3rd ed., Elsevier (2004).
#   E. P. Wigner, Conference on Neutron Physics by Time-of-Flight,
#       Oak Ridge Natl. Lab. Rep. ORNL-2309, 59 (1957) — surmise.
#   F. J. Dyson, J. Math. Phys. 3, 140 (1962) — three-fold way.
#   C. A. Tracy, H. Widom, Commun. Math. Phys. 159, 151 (1994) — F_2.
#   C. A. Tracy, H. Widom, Commun. Math. Phys. 177, 727 (1996) — F_1, F_4.
#   F. Bornemann, Math. Comp. 79, 871 (2010) — high-precision F_β table.
#   Y. Y. Atas, E. Bogomolny, O. Giraud, G. Roux,
#       Phys. Rev. Lett. 110, 084101 (2013) — ⟨r⟩.
# ─────────────────────────────────────────────────────────────────────────────

# ─── Wigner surmise (closed form) ────────────────────────────────────────────

"""
    fetch(::Universality{:RMT}, ::WignerSurmise; β::Int, s::Real) -> Float64

Wigner surmise nearest-neighbour level-spacing distribution `P_β(s)`
for β ∈ {1, 2, 4}.

Closed-form expressions (mean spacing normalised to 1):

* `P_1(s) = (π s / 2) exp(-π s² / 4)`                   GOE
* `P_2(s) = (32 s² / π²) exp(-4 s² / π)`                GUE
* `P_4(s) = (2¹⁸ s⁴ / (3⁶ π³)) exp(-64 s² / (9π))`      GSE

Each integrates to 1 with first moment 1 (verified analytically and
by the standalone test).  Behaves as `P_β(s) ~ s^β` for `s → 0⁺`
(level repulsion) and as `P_β(s) ~ exp(-c_β s²)` for large `s`.
"""
function fetch(::Universality{:RMT}, ::WignerSurmise; β::Int, s::Real, kwargs...)
    β ∈ (1, 2, 4) ||
        throw(DomainError(β, "Universality(:RMT)/WignerSurmise: β must be in {1, 2, 4}"))
    s ≥ 0 || throw(DomainError(s, "WignerSurmise: s must be ≥ 0"))
    sf = float(s)
    if β == 1
        return (π * sf / 2) * exp(-π * sf^2 / 4)
    elseif β == 2
        return (32 * sf^2 / π^2) * exp(-4 * sf^2 / π)
    else  # β == 4
        return (2^18 * sf^4 / (3^6 * π^3)) * exp(-64 * sf^2 / (9π))
    end
end

# ─── Mean ratio ⟨r⟩ (Atas et al. 2013) ───────────────────────────────────────

"""
    fetch(::Universality{:RMT}, ::MeanRatio; β::Int) -> Float64

Mean of consecutive level-spacing ratio
`⟨r⟩ = ⟨min(s_n, s_{n+1}) / max(s_n, s_{n+1})⟩`
for β ∈ {1, 2, 4}, from
Atas-Bogomolny-Giraud-Roux, Phys. Rev. Lett. **110**, 084101 (2013):

| β | ⟨r⟩    |
|---|--------|
| 1 | 0.5307 |
| 2 | 0.5996 |
| 4 | 0.6744 |
"""
function fetch(::Universality{:RMT}, ::MeanRatio; β::Int, kwargs...)
    if β == 1
        return 0.5307
    elseif β == 2
        return 0.5996
    elseif β == 4
        return 0.6744
    end
    return throw(DomainError(β, "Universality(:RMT)/MeanRatio: β must be in {1, 2, 4}"))
end

# ─── Tracy-Widom F_β(x) — tabulated + tail asymptotics ───────────────────────
#
# The data file is parsed once at module load.  Format:
#     x  F_1(x)  F_2(x)  F_4(x)
# with `#`-prefixed comment / blank lines tolerated.

const _TW_DATA_PATH = joinpath(@__DIR__, "data", "tracy_widom_F.txt")

function _load_tw_table(path::AbstractString)
    xs = Float64[]
    f1s = Float64[]
    f2s = Float64[]
    f4s = Float64[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        toks = split(line)
        length(toks) == 4 ||
            error("tracy_widom_F.txt: malformed row, expected 4 tokens, got \"$line\"")
        push!(xs, parse(Float64, toks[1]))
        push!(f1s, parse(Float64, toks[2]))
        push!(f2s, parse(Float64, toks[3]))
        push!(f4s, parse(Float64, toks[4]))
    end
    issorted(xs) || error("tracy_widom_F.txt: x column must be sorted ascending")
    return (xs=xs, f1s=f1s, f2s=f2s, f4s=f4s)
end

const _TW_TABLE = _load_tw_table(_TW_DATA_PATH)

# Linear interpolation on the tabulated grid.  Monotone non-decreasing CDF
# values and monotone abscissa together guarantee the interpolant stays
# in [0, 1] and is monotone non-decreasing — sufficient for Phase 1.
function _tw_interp(x::Float64, ys::Vector{Float64})
    xs = _TW_TABLE.xs
    n = length(xs)
    @inbounds if x ≤ xs[1]
        return ys[1]
    elseif x ≥ xs[n]
        return ys[n]
    end
    # Binary search for the bracket.
    lo, hi = 1, n
    while hi - lo > 1
        mid = (lo + hi) >> 1
        @inbounds if xs[mid] ≤ x
            lo = mid
        else
            hi = mid
        end
    end
    @inbounds x0, x1 = xs[lo], xs[hi]
    @inbounds y0, y1 = ys[lo], ys[hi]
    t = (x - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)
end

# Right tail: 1 - F_β(x) ~ exp(-(2β/3) x^{3/2}) (Tracy-Widom 1994/1996).
# Phase 1 caps at the table's right boundary plus this exponential
# correction normalised against the table's right endpoint, which
# guarantees monotone continuity at the boundary.
function _tw_right_tail(β::Int, x::Float64, y_at_xmax::Float64, xmax::Float64)
    # Exponential decay of (1 - F_β) past xmax; clamp to [y_at_xmax, 1].
    decay_xmax = exp(-(2β / 3) * xmax^1.5)
    decay_x = exp(-(2β / 3) * x^1.5)
    # Match at x = xmax: 1 - F(x) = (1 - y_at_xmax) * decay_x / decay_xmax.
    one_minus_F = (1 - y_at_xmax) * decay_x / decay_xmax
    return clamp(1 - one_minus_F, y_at_xmax, 1.0)
end

# Left tail: F_β(x) ~ τ_β exp(-(β/24) |x|^3) for x → -∞ (Tracy-Widom 1994).
# Match continuously to the table's left endpoint.
function _tw_left_tail(β::Int, x::Float64, y_at_xmin::Float64, xmin::Float64)
    decay_xmin = exp(-(β / 24) * abs(xmin)^3)
    decay_x = exp(-(β / 24) * abs(x)^3)
    F = y_at_xmin * decay_x / decay_xmin
    return clamp(F, 0.0, y_at_xmin)
end

"""
    fetch(::Universality{:RMT}, ::TracyWidom; β::Int, x::Real) -> Float64

Tracy-Widom CDF `F_β(x) = P[ξ_β ≤ x]` for β ∈ {1, 2, 4} (GOE, GUE,
GSE largest-eigenvalue limit law).

QAtlas Phase 1 evaluates `F_β` from a tabulated grid compiled from
Bornemann, *On the numerical evaluation of Fredholm determinants*,
Math. Comp. **79**, 871 (2010), Table 1; the table covers
`x ∈ [-4.0, 4.0]` for all three β.  Inside the table support the
interpolant is piecewise-linear and monotone non-decreasing.  Outside
the table the function returns the Tracy-Widom 1994/1996 left/right
tail asymptotics, matched continuously to the table boundary.

Reference checkpoints (Bornemann 2010, Table 1):

* `F_1(0) ≈ 0.8319`
* `F_2(0) ≈ 0.9694`
* `F_4(0) ≈ 0.99966`

A direct Painlevé-II integrator (DifferentialEquations.jl based) is
deferred to Phase 2; see issue #151.
"""
function fetch(::Universality{:RMT}, ::TracyWidom; β::Int, x::Real, kwargs...)
    ys = if β == 1
        _TW_TABLE.f1s
    elseif β == 2
        _TW_TABLE.f2s
    elseif β == 4
        _TW_TABLE.f4s
    else
        throw(DomainError(β, "Universality(:RMT)/TracyWidom: β must be in {1, 2, 4}"))
    end
    xs = _TW_TABLE.xs
    xf = float(x)
    if xf < xs[1]
        return _tw_left_tail(β, xf, ys[1], xs[1])
    elseif xf > xs[end]
        return _tw_right_tail(β, xf, ys[end], xs[end])
    end
    return _tw_interp(xf, ys)
end

# ─── Spectral form factor (large-N, GUE plateau; Phase 1 of issue #243) ──────

"""
    fetch(::Universality{:RMT}, ::SpectralFormFactor; ensemble::Symbol=:GUE, τ::Real=Inf) -> Float64

Disorder-averaged spectral form factor `K(τ)` for the random-matrix
universality class in the large-`N` thermodynamic limit, with
`τ = t / N` (Heisenberg time `τ_H = 2π`).

QAtlas Phase 1 (issue #243) exposes only the GUE ensemble in the
**late-time plateau** regime `τ ≥ 2π`, where the disorder-averaged
SFF saturates universally to `K(τ→∞) = 1` (Mehta 2004 §16;
Cotler et al. 2017, arXiv:1611.04650).

The Mehta connection formula
`K(τ) = (τ/(2π)) − (τ/(4π)) log|1 − τ/(2π)|`
on the ramp side `τ < 2π`, and the GOE / GSE sigma-model closed
forms, are deferred to Phase 2.

# Errors
* `DomainError` if `ensemble ≠ :GUE` (GOE / GSE Phase 2).
* `DomainError` if `τ < 2π` (ramp regime, Phase 2).
"""
function fetch(
    ::Universality{:RMT},
    ::SpectralFormFactor;
    ensemble::Symbol=:GUE,
    τ::Real=Inf,
    kwargs...,
)
    if ensemble != :GUE
        throw(
            DomainError(
                ensemble,
                "Universality(:RMT)/SpectralFormFactor: Phase 1 supports only " *
                "ensemble = :GUE.  GOE / GSE (sigma-model formulae, Mehta 2004 " *
                "§16) are deferred to Phase 2 of issue #243.  Got ensemble = :" *
                string(ensemble) *
                ".",
            ),
        )
    end
    if τ < 2π
        throw(
            DomainError(
                τ,
                "Universality(:RMT)/SpectralFormFactor: Phase 1 supports only " *
                "the late-time plateau regime τ ≥ 2π (= Heisenberg time τ_H). " *
                "The τ < 2π ramp K(τ) = τ/π and the Mehta connection formula " *
                "are deferred to Phase 2 of issue #243.  Got τ = $τ.",
            ),
        )
    end
    return 1.0  # GUE late-time plateau (Mehta 2004 §16; CHM 2016)
end

# -----------------------------------------------------------------------------
# Wigner-semicircle moments (Wigner 1955, Mehta)
# -----------------------------------------------------------------------------

"""
    fetch(::Universality{:RMT}, ::WignerSemicircleMoment, ::Infinite;
          n::Integer, kwargs...) -> Float64

Universal moment `m_n = E[x^n]` of the Wigner semicircle distribution
on `[-2, 2]`. Even moments are Catalan numbers

    m_{2k} = C_k = binomial(2k, k) / (k + 1),

odd moments vanish. The large-N eigenvalue density of all three Gaussian
ensembles (GOE, GUE, GSE) converges to the same semicircle, so the
moments are class-independent within RMT.
"""
function fetch(
    ::Universality{:RMT}, ::WignerSemicircleMoment, ::Infinite; n::Integer, kwargs...
)
    n >= 0 ||
        throw(DomainError(n, "WignerSemicircleMoment: n must be non-negative; got n = $n."))
    if isodd(n)
        return 0.0
    end
    k = div(n, 2)
    return Float64(binomial(big(2 * k), big(k)) / (k + 1))
end

function fetch(m::Universality{:RMT}, q::WignerSemicircleMoment; n::Integer, kwargs...)
    return fetch(m, q, Infinite(); n=n, kwargs...)
end
