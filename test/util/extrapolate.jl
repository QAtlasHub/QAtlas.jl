# ─────────────────────────────────────────────────────────────────────────────
# test/util/extrapolate.jl
#
# Finite-size extrapolation for `route = :ed_finite_size` verify() cards.
#
# `verify` compares `last(independent)` to the infinite-system `fetch`
# value with NO built-in extrapolation. An observable whose finite-N ED
# value converges only algebraically (O(1/N), O(1/N²), …) therefore
# cannot be passed as a raw size sweep — the largest tractable N is still
# far from the N→∞ value (this is exactly why the first TFIM MassGap ED
# card failed; see commit 087132ed). Extrapolate to N→∞ HERE and pass
# the returned `.value` as a scalar `independent`, with
# `agree_within ≳ k · .uncertainty`.
#
# `extrapolate_inf(Ns, vals; power=1)` fits, by least squares,
#     v(N) ≈ v∞ + Σ_{k≥1} c_k / N^(power·k)
# at the highest numerically-stable order for the given number of points,
# and returns `(; value = v∞, uncertainty)` where `uncertainty` is the
# magnitude of the change in v∞ when the fit order is dropped by one —
# the standard "difference between successive approximants" error proxy.
#
# Dependencies (expected to be `using`'d by the including test file):
#   LinearAlgebra  — least-squares `\`
# ─────────────────────────────────────────────────────────────────────────────

# Least-squares fit ys ≈ Σ_{j=0}^{deg} a_j·xs^j; return a_0 (the xs→0 limit).
function _extrap_intercept(
    xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, deg::Integer
)
    A = [x^j for x in xs, j in 0:deg]
    a = A \ collect(float.(ys))
    return a[1]
end

"""
    extrapolate_inf(Ns, vals; power=1) -> (; value, uncertainty)

Extrapolate a finite-size sweep `vals[i] = observable(Ns[i])` to the
thermodynamic limit `N → ∞`, modelling the tail as a series in
`1 / N^power`. `value` is the extrapolated `N→∞` estimate; `uncertainty`
is `|value − value_one_order_lower|`, a conservative proxy for the
extrapolation error (use it to choose a defensible `agree_within`).
"""
function extrapolate_inf(
    Ns::AbstractVector{<:Real}, vals::AbstractVector{<:Real}; power::Integer=1
)
    n = length(Ns)
    n == length(vals) || error("extrapolate_inf: Ns and vals length mismatch")
    n >= 2 || error("extrapolate_inf: need ≥ 2 points, got $n")
    all(N -> N > 0, Ns) || error("extrapolate_inf: all Ns must be > 0")
    power >= 1 || error("extrapolate_inf: power must be ≥ 1, got $power")
    xs = [1.0 / float(N)^power for N in Ns]
    hi = min(n - 1, 3)                       # cap order for stability
    v_hi = _extrap_intercept(xs, vals, hi)
    v_lo = _extrap_intercept(xs, vals, max(hi - 1, 0))
    return (; value=v_hi, uncertainty=abs(v_hi - v_lo))
end
