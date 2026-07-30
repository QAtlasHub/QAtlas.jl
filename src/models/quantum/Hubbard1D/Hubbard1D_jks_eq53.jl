# src/models/quantum/Hubbard1D/Hubbard1D_jks_eq53.jl
#
# JKS 1998 eq (53) — the six-unknown real-axis form of the Hubbard QTM NLIE.
# Included from inside `module Hubbard1DJKSNLIE`.
#
# WHY A SECOND SOLVER RATHER THAN A PATCH
#
# The eq (47) solver in Hubbard1D_jks_nlie.jl carries FOUR real-axis unknowns
# (b, b_bar, c, c_bar) and reads `b_bar`/`c_bar` as barred FUNCTIONS. eq (51)
# fixes the real-axis unknowns as SIX — two boundary values each of three
# functions:
#
#     b^{+-}(x)    = b(x +- i alpha)      alpha in (0, gamma]
#     c^{+-}(x)    = c(x +- i0)
#     cbar^{+-}(x) = cbar(x +- i0)
#
# so the `+-` boundary-value index had been collapsed into the bar index. That
# is not a tolerance issue: eq (52) defines
#
#     Dlog C    = log(C^+ / C^-)          C^{+-}    = 1 + c^{+-}
#     Dlog Cbar = log(Cbar^+ / Cbar^-)    Cbar^{+-} = 1 + cbar^{+-}
#
# and eq (53) carries these BOTH as a convolution and as an explicit `+- 1/2`
# boundary term. Neither is representable with one boundary value per function,
# so the old solver could not be corrected in place. See #798.
#
# THE EQUATIONS (eq 53, verbatim from the arXiv LaTeX source with the overlines
# intact — the PDF text layer drops them, which is how #798 was mis-transcribed)
#
#   log b^s    = -beta H
#                - K2_{s a - a} * log B^+ + K2_{s a + a} * log B^-
#                - K1bar_{s a} * Dlog(cbar / Cbar)
#
#   log c^s    = Psi_c^s
#                + K1bar_{-a} * log Bbar^+ - K1bar_{a} * log Bbar^-
#                + K1bar_{0} * Dlog Cbar + s/2 * Dlog Cbar
#
#   log cbar^s = Psibar_c^s
#                - K1_{-a} * log B^+ + K1_{a} * log B^-
#                - K1_{0} * Dlog C + s/2 * Dlog C
#
# with s = +-1, `f_d(x) = f(x + i d)`, and
#
#   B^{+-}    = 1 + b^{+-}        Bbar^{+-} = 1 + 1/b^{+-}
#   C^{+-}    = 1 + c^{+-}        Cbar^{+-} = 1 + cbar^{+-}
#
#   Psi_c^s      = -beta U/2 + beta (mu + H/2) + log phi_{s 0}                  (eq 54)
#   Psibar_c^s   = -beta U/2 - beta (mu + H/2) - log phi_{s 0}                  (eq 55)
#
# CAREFUL — the b equation appears TWICE in the paper, and the two forms are not
# interchangeable term by term. eq (58) rewrites eq (53) using the identity
# `Dlog cbar = -Dlog phi + Dlog C`:
#
#   Dlog(cbar/Cbar) = -Dlog phi + Dlog(C/Cbar)
#
# so moving the `-Dlog phi` piece through `K1bar_{s a}` is exactly what turns the
# `-beta H` driving term of eq (53) into
#
#   Psi_b^s = -beta U - beta H + log phi_{s a} - log phi_{s a - 2 gamma}   (eq 59)
#
# in eq (58), whose convolution term is `-K1bar_{s a} * Dlog(C/Cbar)`. Taking
# eq (59)'s driving term together with eq (53)'s convolution term double-counts
# the `log phi` contribution and puts an O(beta) error into `log b` — MEASURED as
# a constant `+0.57` offset in `f` at U = 4 that no grid refinement removed.
# eq (53) is used here; [`jks53_driving_b`](@ref) supplies eq (59) for the
# independent eq (58) route.
#
# (In eq (59) the second shift is `s a - 2 gamma`, not `-s a`; those agree only
# at alpha = gamma. The eq (47) solver used `-s a` while passing
# alpha = U/6 = 2 gamma/3.)
#
# CONVENTION — the paper's U term is SYMMETRIC
#
#   H = sum_i { sum_s -(c^+_{i+1,s} c_{i,s} + h.c.) + U (n_{i,-} - 1/2)(n_{i,+} - 1/2) }
#       - sum_i [ mu (n_{i,+} + n_{i,-}) + (H/2)(n_{i,+} - n_{i,-}) ]
#
# so **half filling is mu = 0**, not mu = U/2. Translating to the plain
# `U n_up n_down` form shifts the chemical potential by U/2 and adds U/4 per site:
#
#   f_paper(mu) = f_plain(mu + U/2) + U/4
#
# This matters for every oracle. The single-site energies `E - mu N` are
# `U/4, -U/4-mu, -U/4-mu, U/4-2mu`, so
#
#   z = e^{-beta U/4} (1 + e^{2 beta mu}) + 2 e^{beta (U/4 + mu)}
#
# exactly, with the `t = 1` hopping entering `beta f` only at O(beta^2). At
# mu = 0 that collapses to `z = 4 cosh(beta U/4)`, whose logarithm has NO term
# linear in beta -- a sharper high-T anchor than the generic case.
#
# WHICH GRID EACH UNKNOWN LIVES ON
#
# b is needed on the whole real line (the wide loop Im s = +- alpha). c and cbar
# are needed only on [-1, 1]: p.14 states they need be evaluated only just above
# and below that interval, and eq (53)'s closing remark that Dlog C and
# Dlog Cbar "vanish outside the interval [-1,1]" is what restricts the narrow
# convolutions. So:
#
#   b^{+-}    on a uniform grid over [-x_max, x_max]        (Nw points)
#   c^{+-}, cbar^{+-} on Chebyshev-Gauss nodes of [-1, 1]   (Nn points)
#
# The `K1bar_{0}`/`K1_{0}` convolutions are unshifted and therefore singular on
# the diagonal; eq (53) says they are Cauchy principal values. See
# [`jks53_pv_matrix`](@ref).

# ═══════════════════════════════════════════════════════════════════════════════
# Grids
# ═══════════════════════════════════════════════════════════════════════════════

"""
    JKSGrids53(Nw, Nn, gamma, alpha; x_max=32.0)

The two grids eq (53) needs: a uniform one over `[-x_max, x_max]` carrying the
two boundary values of `b`, and Chebyshev-Gauss nodes of `[-1, 1]` carrying
`c` and `cbar`.

The narrow nodes are `x_j = cos((2j-1) pi / (2 Nn))` with weights
`w_j = (pi/Nn) sqrt(1 - x_j^2)`, which integrate `[-1, 1]` exactly for
polynomials through degree `2 Nn - 1` after the `1/sqrt(1-x^2)` factor is
absorbed. The same nodes make the free-energy evaluator's `1/sqrt(1-x^2)`
weight exact — see [`free_energy_jks53`](@ref).

`gamma = U/4` (eq 30) and is NOT bounded above.

`alpha` must satisfy `0 < alpha < gamma` **strictly**, which is tighter than the
`0 < alpha <= gamma` of eq (51). The reason is the shift `2 alpha` that eq (53)
applies to `K2`: `K2` has poles at `s = +- 2 i gamma`, so

    K2(x + 2 i alpha) has denominator  x^2 + 4 i alpha x - 4 alpha^2 + 4 gamma^2

and at `alpha = gamma` that is `x (x + 4 i gamma)` — a pole sitting **on the real
axis at x = 0**, i.e. on the integration contour. The paper names `alpha = gamma`
as a convenient choice, but a real-axis quadrature cannot take it without a
principal value; at `alpha = 2 gamma/3` (its other named choice, and the default
here) the roots are `2 i gamma/3` and `-10 i gamma/3`, both off the axis.
"""
struct JKSGrids53
    Nw::Int
    Nn::Int
    gamma::Float64
    alpha::Float64
    x_max::Float64
    xw::Vector{Float64}
    ww::Vector{Float64}
    xn::Vector{Float64}
    wn::Vector{Float64}
    function JKSGrids53(Nw::Int, Nn::Int, gamma::Real, alpha::Real; x_max::Real=32.0)
        Nw > 1 || throw(DomainError(Nw, "Nw must be > 1"))
        Nn > 1 || throw(DomainError(Nn, "Nn must be > 1"))
        gamma > 0 || throw(DomainError(gamma, "gamma must be > 0"))
        0 < alpha < gamma || throw(
            DomainError(
                alpha,
                "need 0 < alpha < gamma strictly: at alpha = gamma the K2 shift " *
                "2 alpha puts a pole of K2 on the real axis at x = 0",
            ),
        )
        x_max > 1 || throw(DomainError(x_max, "x_max must exceed the cut at 1"))
        xw = collect(range(-x_max, x_max; length=Nw))
        ww = fill(xw[2] - xw[1], Nw)
        xn = [cos((2j - 1) * pi / (2 * Nn)) for j in 1:Nn]
        wn = [(pi / Nn) * sqrt(1 - x^2) for x in xn]
        return new(Nw, Nn, Float64(gamma), Float64(alpha), Float64(x_max), xw, ww, xn, wn)
    end
end

"""
    JKSState53(grids)

The six eq (53) unknowns: `b_plus`/`b_minus` of length `grids.Nw`, and
`c_plus`/`c_minus`/`cbar_plus`/`cbar_minus` of length `grids.Nn`.

Stored as the auxiliary functions themselves rather than their logs so that the
`Bbar = 1 + 1/b` and `C = 1 + c` combinations of eq (52) read directly.
"""
mutable struct JKSState53
    b_plus::Vector{ComplexF64}
    b_minus::Vector{ComplexF64}
    c_plus::Vector{ComplexF64}
    c_minus::Vector{ComplexF64}
    cbar_plus::Vector{ComplexF64}
    cbar_minus::Vector{ComplexF64}
end

function JKSState53(grids::JKSGrids53)
    return JKSState53(
        zeros(ComplexF64, grids.Nw),
        zeros(ComplexF64, grids.Nw),
        zeros(ComplexF64, grids.Nn),
        zeros(ComplexF64, grids.Nn),
        zeros(ComplexF64, grids.Nn),
        zeros(ComplexF64, grids.Nn),
    )
end

"Total number of unknowns: `2 Nw + 4 Nn`."
jks53_ndof(grids::JKSGrids53) = 2 * grids.Nw + 4 * grids.Nn

"Flatten the six unknowns into one complex vector, b first then c then cbar."
function jks53_pack(st::JKSState53)
    return vcat(
        log.(st.b_plus),
        log.(st.b_minus),
        log.(st.c_plus),
        log.(st.c_minus),
        log.(st.cbar_plus),
        log.(st.cbar_minus),
    )
end

"Inverse of [`jks53_pack`](@ref): rebuild the state from the log-vector."
function jks53_unpack(v::AbstractVector, grids::JKSGrids53)
    Nw, Nn = grids.Nw, grids.Nn
    o = 0
    bp = exp.(v[(o + 1):(o + Nw)])
    o += Nw
    bm = exp.(v[(o + 1):(o + Nw)])
    o += Nw
    cp = exp.(v[(o + 1):(o + Nn)])
    o += Nn
    cm = exp.(v[(o + 1):(o + Nn)])
    o += Nn
    cbp = exp.(v[(o + 1):(o + Nn)])
    o += Nn
    cbm = exp.(v[(o + 1):(o + Nn)])
    return JKSState53(bp, bm, cp, cm, cbp, cbm)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Kernel matrices
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks53_kernel_matrix(kernel, x_out, x_in, w_in, gamma, shift) -> Matrix

Discrete convolution `(K_shift * f)(x_out[j]) = sum_k K(x_out[j] - x_in[k] + i shift) f(x_in[k]) w_in[k]`.

Used for every eq (53) convolution whose shift is nonzero, in all four
grid pairings (wide->wide, narrow->wide, wide->narrow, narrow->narrow). A
nonzero shift keeps the argument off the pole of `K1`/`K1bar` at the origin, so
no regularization is needed; `K2` has no pole at the origin at all.
"""
function jks53_kernel_matrix(
    kernel,
    x_out::AbstractVector,
    x_in::AbstractVector,
    w_in::AbstractVector,
    gamma::Real,
    shift::Real,
)
    M = zeros(ComplexF64, length(x_out), length(x_in))
    for j in eachindex(x_out), k in eachindex(x_in)
        M[j, k] = kernel(x_out[j] - x_in[k] + im * shift, gamma) * w_in[k]
    end
    return M
end

"""
    jks53_kernel_line_integral(kernel, gamma, shift) -> Float64

The exact `int_{-inf}^{inf} K(x + i shift) dx` for each eq (38) kernel, by
residues. Needed because the wide grid is truncated at `x_max` while these
kernels only decay like `1/x^2`, so the discrete row sums fall short by
`O(gamma / x_max)` — 2% at `gamma = 1, x_max = 32`, which is four orders of
magnitude above the solver tolerance.

    K2:    poles at +-2 i gamma.  |shift| < 2 gamma  ->  1
    K1:    poles at 0 and -2 i gamma.  shift > 0 -> 0;  -2 gamma < shift < 0 -> 1
    K1bar: poles at 0 and +2 i gamma.  shift > 0 -> 1;  -2 gamma < shift < 0 -> 0

Each value is a residue and shares no arithmetic with the quadrature, so
[`jks53_apply_wide`](@ref) using them is a genuine correction and not a
self-consistent rescaling.
"""
function jks53_kernel_line_integral(kernel, gamma::Real, shift::Real)
    abs(shift) < 2 * gamma ||
        throw(DomainError(shift, "|shift| must be < 2 gamma for these residues"))
    if kernel === jks_kernel_K2
        return 1.0
    elseif kernel === jks_kernel_K1
        return shift > 0 ? 0.0 : 1.0
    elseif kernel === jks_kernel_K1bar
        return shift > 0 ? 1.0 : 0.0
    end
    return throw(ArgumentError("no line integral registered for $kernel"))
end

"""
    JKSWideOp(M, I_exact)

A convolution whose INPUT lives on the truncated wide grid, carrying both the
quadrature matrix and the exact whole-line integral of its kernel.

[`jks53_apply_wide`](@ref) uses the pair to add back the tail analytically.
Convolutions whose input is the narrow grid need no such thing: `[-1, 1]` is
compact and the Chebyshev rule covers all of it.
"""
struct JKSWideOp
    M::Matrix{ComplexF64}
    I_exact::Float64
end

"""
    jks53_apply_wide(op, f) -> Vector{ComplexF64}

`(K * f)(x_j)` with the truncated tail restored:

    sum_k M[j,k] f_k  +  f_inf * (I_exact - sum_k M[j,k])

The correction is exact whenever `f` has reached its asymptote outside the grid,
and `f` here is always `log B^{+-}` or `log Bbar^{+-}`, which do: the b driving
term tends to the constant `-beta U - beta H - i beta U` as `|x| -> inf` (the two
`log phi` shifts differ by `-2 i gamma`, so their difference has no growing
part), and the convolutions tend to constants with it.

`f_inf` is taken as the mean of the two edge values, which is the cheapest
estimate that stays correct under the `x -> -x` symmetry of the grid. The row
sums are computed per row on purpose: an edge row sees a different part of the
tail than a central one.
"""
function jks53_apply_wide(op::JKSWideOp, f::AbstractVector)
    f_inf = (f[1] + f[end]) / 2
    out = op.M * f
    rows = sum(op.M; dims=2)
    @inbounds for j in eachindex(out)
        out[j] += f_inf * (op.I_exact - rows[j])
    end
    return out
end

"""
    jks53_pv_matrix(x, w) -> Matrix{Float64}

Cauchy principal-value operator for `PV int_{-1}^{1} f(y) / (x_j - y) dy` on the
nodes `x` with weights `w`.

Off the diagonal this is `w_k / (x_j - x_k)`. The diagonal is fixed by requiring
the rule be **exact for a constant**, whose principal value is known in closed
form:

    PV int_{-1}^{1} dy / (x - y) = log |(x + 1) / (x - 1)|

so `P[j,j] = L(x_j) - sum_{k != j} w_k / (x_j - x_k)`. Equivalently this is the
standard singularity-subtraction rule `PV int f/(x-y) = int [f(y) - f(x)]/(x-y)
+ f(x) L(x)`, whose integrand is regular (it tends to `-f'(x)`).

Exactness on constants is what the test asserts — it is an independent statement
about the quadrature, not a restatement of the construction, because the closed
form for `L` comes from outside the rule.
"""
function jks53_pv_matrix(x::AbstractVector, w::AbstractVector)
    N = length(x)
    P = zeros(Float64, N, N)
    for j in 1:N
        acc = 0.0
        for k in 1:N
            k == j && continue
            P[j, k] = w[k] / (x[j] - x[k])
            acc += P[j, k]
        end
        L = log(abs((x[j] + 1) / (x[j] - 1)))
        P[j, j] = L - acc
    end
    return P
end

"""
    jks53_K1_pv(grids) -> Matrix{ComplexF64}
    jks53_K1bar_pv(grids) -> Matrix{ComplexF64}

The unshifted narrow-grid convolutions `K1_{0} *` and `K1bar_{0} *` of eq (53),
as principal values.

Splitting eq (38) back into Cauchy kernels isolates the singular piece:

    K1(s)    = (1/2pi i) [  1/s - 1/(s + 2i gamma) ]
    K1bar(s) = (1/2pi i) [ -1/s + 1/(s - 2i gamma) ]

so each is `+-1/(2pi i)` times the principal-value operator plus a regular
remainder whose denominator is bounded away from zero by `2 gamma`.
"""
function jks53_K1_pv(grids::JKSGrids53)
    x, w, g = grids.xn, grids.wn, grids.gamma
    P = jks53_pv_matrix(x, w)
    M = zeros(ComplexF64, length(x), length(x))
    for j in eachindex(x), k in eachindex(x)
        reg = -w[k] / (x[j] - x[k] + 2im * g)
        M[j, k] = (P[j, k] + reg) / (2im * pi)
    end
    return M
end

function jks53_K1bar_pv(grids::JKSGrids53)
    x, w, g = grids.xn, grids.wn, grids.gamma
    P = jks53_pv_matrix(x, w)
    M = zeros(ComplexF64, length(x), length(x))
    for j in eachindex(x), k in eachindex(x)
        reg = w[k] / (x[j] - x[k] - 2im * g)
        M[j, k] = (-P[j, k] + reg) / (2im * pi)
    end
    return M
end

"""
    JKSOperators53(grids)

Every convolution eq (53) needs, built once per grid. Field names follow the
equation they appear in: `bp_from_Bm` is the operator acting on `log B^-` in the
`log b^+` equation, and so on.

Operators whose INPUT is the truncated wide grid are [`JKSWideOp`](@ref) and must
be applied with [`jks53_apply_wide`](@ref) so the `O(gamma/x_max)` tail is
restored. Operators reading the narrow grid are plain matrices: `[-1, 1]` is
compact, so there is no tail to restore.
"""
struct JKSOperators53
    # b equation: K2 at shifts (s a - a) and (s a + a), s = +-1
    bp_from_Bp::JKSWideOp   # K2_{0}
    bp_from_Bm::JKSWideOp   # K2_{+2a}
    bm_from_Bp::JKSWideOp   # K2_{-2a}
    bm_from_Bm::JKSWideOp   # K2_{0}
    # b equation: K1bar_{s a} acting on a narrow-grid function
    bp_from_narrow::Matrix{ComplexF64}   # K1bar_{+a}, Nw x Nn
    bm_from_narrow::Matrix{ComplexF64}   # K1bar_{-a}, Nw x Nn
    # c equation: K1bar_{-a} on log Bbar^+, K1bar_{+a} on log Bbar^-
    c_from_Bbarp::JKSWideOp               # K1bar_{-a}, Nn x Nw
    c_from_Bbarm::JKSWideOp               # K1bar_{+a}, Nn x Nw
    c_pv::Matrix{ComplexF64}              # K1bar_{0},  Nn x Nn
    # cbar equation: K1_{-a} on log B^+, K1_{+a} on log B^-
    cbar_from_Bp::JKSWideOp               # K1_{-a}, Nn x Nw
    cbar_from_Bm::JKSWideOp               # K1_{+a}, Nn x Nw
    cbar_pv::Matrix{ComplexF64}           # K1_{0},  Nn x Nn
end

function JKSOperators53(grids::JKSGrids53)
    g, a = grids.gamma, grids.alpha
    xw, ww, xn, wn = grids.xw, grids.ww, grids.xn, grids.wn
    K2, K1, K1b = jks_kernel_K2, jks_kernel_K1, jks_kernel_K1bar
    wide(kernel, x_out, shift) = JKSWideOp(
        jks53_kernel_matrix(kernel, x_out, xw, ww, g, shift),
        jks53_kernel_line_integral(kernel, g, shift),
    )
    return JKSOperators53(
        wide(K2, xw, 0.0),
        wide(K2, xw, 2a),
        wide(K2, xw, -2a),
        wide(K2, xw, 0.0),
        jks53_kernel_matrix(K1b, xw, xn, wn, g, a),
        jks53_kernel_matrix(K1b, xw, xn, wn, g, -a),
        wide(K1b, xn, -a),
        wide(K1b, xn, a),
        jks53_K1bar_pv(grids),
        wide(K1, xn, -a),
        wide(K1, xn, a),
        jks53_K1_pv(grids),
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Driving terms
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks53_log_phi_narrow(x, beta, sign) -> ComplexF64

`log phi_{s 0}(x)` on the cut, i.e. the two boundary values of
`log phi(s) = -2 beta sqrt(s^2 - 1)` at `s = x +- i0` for `|x| <= 1`:

    log phi_{+-0}(x) = -+ 2 i beta sqrt(1 - x^2)

Written out rather than obtained by evaluating the complex form at `x +- 0.0im`,
which would depend on the sign of a floating-point zero surviving into
`sqrt`.
"""
function jks53_log_phi_narrow(x::Real, beta::Real, sgn::Int)
    abs(x) <= 1 || throw(DomainError(x, "narrow-grid driving needs |x| <= 1"))
    return ComplexF64(0.0, -sgn * 2 * beta * sqrt(1 - x^2))
end

"""
    jks53_driving_b(grids, beta, U, sgn; H=0.0) -> Vector{ComplexF64}

`Psi_b^s = -beta U - beta H + log phi_{s a} - log phi_{s a - 2 gamma}`, eq (59).

This is the driving term of **eq (58)**, not of eq (53). It goes with eq (58)'s
convolution term `-K1bar_{s a} * Dlog(C/Cbar)`; pairing it with eq (53)'s
`-K1bar_{s a} * Dlog(cbar/Cbar)` double-counts `log phi`. [`jks53_residual`](@ref)
implements eq (53) and so does not use this.

The second shift is `s alpha - 2 gamma`, not `-s alpha`. Those coincide only at
`alpha = gamma`; the eq (47) solver used the latter while passing
`alpha = U/6 = 2 gamma/3`.
"""
function jks53_driving_b(grids::JKSGrids53, beta::Real, U::Real, sgn::Int; H::Real=0.0)
    a, g = grids.alpha, grids.gamma
    return [
        -beta * U - beta * H + jks_log_phi_complex(x + im * sgn * a, beta) -
        jks_log_phi_complex(x + im * (sgn * a - 2 * g), beta) for x in grids.xw
    ]
end

"""
    jks53_driving_c(grids, beta, U, mu, sgn; H=0.0) -> Vector{ComplexF64}

`Psi_c^s = -beta U/2 + beta (mu + H/2) + log phi_{s 0}`, eq (54).
"""
function jks53_driving_c(
    grids::JKSGrids53, beta::Real, U::Real, mu::Real, sgn::Int; H::Real=0.0
)
    return [
        -beta * U / 2 + beta * (mu + H / 2) + jks53_log_phi_narrow(x, beta, sgn) for
        x in grids.xn
    ]
end

"""
    jks53_driving_cbar(grids, beta, U, mu, sgn; H=0.0) -> Vector{ComplexF64}

`Psibar_c^s = -beta U/2 - beta (mu + H/2) - log phi_{s 0}`, eq (55).
"""
function jks53_driving_cbar(
    grids::JKSGrids53, beta::Real, U::Real, mu::Real, sgn::Int; H::Real=0.0
)
    return [
        -beta * U / 2 - beta * (mu + H / 2) - jks53_log_phi_narrow(x, beta, sgn) for
        x in grids.xn
    ]
end

# ═══════════════════════════════════════════════════════════════════════════════
# Residual
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks53_residual(st, grids, ops, beta, U, mu; H=0.0) -> Vector{ComplexF64}

Residual of eq (53), stacked in the [`jks53_pack`](@ref) order. Zero at a
solution.

Each of the six blocks is `log(unknown) - rhs`, with `rhs` assembled term by
term from eq (53); nothing is folded or pre-simplified, so the code reads
against the equation.
"""
function jks53_residual(
    st::JKSState53,
    grids::JKSGrids53,
    ops::JKSOperators53,
    beta::Real,
    U::Real,
    mu::Real;
    H::Real=0.0,
)
    # eq (52) combinations
    log_Bp = log.(1 .+ st.b_plus)
    log_Bm = log.(1 .+ st.b_minus)
    log_Bbarp = log.(1 .+ 1 ./ st.b_plus)
    log_Bbarm = log.(1 .+ 1 ./ st.b_minus)
    log_Cp = log.(1 .+ st.c_plus)
    log_Cm = log.(1 .+ st.c_minus)
    log_Cbarp = log.(1 .+ st.cbar_plus)
    log_Cbarm = log.(1 .+ st.cbar_minus)
    dlog_C = log_Cp .- log_Cm
    dlog_Cbar = log_Cbarp .- log_Cbarm
    # Dlog(cbar / Cbar) in the b equation
    dlog_cbar_over_Cbar =
        (log.(st.cbar_plus) .- log_Cbarp) .- (log.(st.cbar_minus) .- log_Cbarm)

    # eq (53): the b equation's driving term is -beta H and nothing else. The
    # log phi terms belong to eq (58) (see the header note); adding both
    # double-counts them.
    psi_b = fill(ComplexF64(-beta * H), grids.Nw)
    psi_cp = jks53_driving_c(grids, beta, U, mu, +1; H=H)
    psi_cm = jks53_driving_c(grids, beta, U, mu, -1; H=H)
    psi_cbp = jks53_driving_cbar(grids, beta, U, mu, +1; H=H)
    psi_cbm = jks53_driving_cbar(grids, beta, U, mu, -1; H=H)

    W = jks53_apply_wide
    rhs_bp =
        psi_b .- W(ops.bp_from_Bp, log_Bp) .+ W(ops.bp_from_Bm, log_Bm) .-
        ops.bp_from_narrow * dlog_cbar_over_Cbar
    rhs_bm =
        psi_b .- W(ops.bm_from_Bp, log_Bp) .+ W(ops.bm_from_Bm, log_Bm) .-
        ops.bm_from_narrow * dlog_cbar_over_Cbar

    conv_c =
        W(ops.c_from_Bbarp, log_Bbarp) .- W(ops.c_from_Bbarm, log_Bbarm) .+
        ops.c_pv * dlog_Cbar
    rhs_cp = psi_cp .+ conv_c .+ 0.5 .* dlog_Cbar
    rhs_cm = psi_cm .+ conv_c .- 0.5 .* dlog_Cbar

    conv_cbar =
        -W(ops.cbar_from_Bp, log_Bp) .+ W(ops.cbar_from_Bm, log_Bm) .- ops.cbar_pv * dlog_C
    rhs_cbp = psi_cbp .+ conv_cbar .+ 0.5 .* dlog_C
    rhs_cbm = psi_cbm .+ conv_cbar .- 0.5 .* dlog_C

    return vcat(
        log.(st.b_plus) .- rhs_bp,
        log.(st.b_minus) .- rhs_bm,
        log.(st.c_plus) .- rhs_cp,
        log.(st.c_minus) .- rhs_cm,
        log.(st.cbar_plus) .- rhs_cbp,
        log.(st.cbar_minus) .- rhs_cbm,
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Free energy — eq (56), both forms
# ═══════════════════════════════════════════════════════════════════════════════

"""
    jks53_script_K(x, gamma, shift) -> ComplexF64

`cal K` of eq (57) with its argument shifted by `i * shift`.

eq (57) gives two forms,

    cal K(x) = -(2 pi sqrt(1 - x^2))^-1 = (2 pi i x sqrt(1 - 1/x^2))^-1

and fixes the branch by `cal K(x) ~ 1/(2 pi i x)` for large `|x|`. Only the
SECOND form reproduces that under the principal square root: for real `x > 1`
the first gives `sqrt(1 - x^2) = i sqrt(x^2 - 1)` and hence
`cal K ~ -1/(2 pi i x)`, the wrong sign. So the second form is what is
implemented, and the first is not interchangeable with it outside `|x| < 1`.
"""
function jks53_script_K(x::Number, shift::Real)
    s = x + im * shift
    return ComplexF64(1 / (2im * pi * s * sqrt(1 - 1 / s^2 + 0im)))
end

"""
    jks53_script_K_line_integral(shift) -> Float64

The exact `int_{-inf}^{inf} cal K(x + i shift) dx`, needed for the same reason as
[`jks53_kernel_line_integral`](@ref): the `:wide` form of eq (56) integrates over
the whole line while the grid stops at `x_max`, and `cal K` only decays like
`1/x`.

`cal K` is analytic off the cut `[-1, 1]` and tends to `1/(2 pi i s)`. Closing the
line `Im s = shift` away from the cut therefore leaves only the arc at infinity,
whose contribution is `int dtheta / (2 pi)` over a half turn:

    shift > 0  ->  -1/2        shift < 0  ->  +1/2

The individual line integrals diverge logarithmically; each is the symmetric
(principal-value) limit, which is what a grid symmetric about `x = 0` computes,
and the `1/x` parts cancel in the differences eq (56) actually takes.
"""
function jks53_script_K_line_integral(shift::Real)
    shift == 0 && throw(DomainError(shift, "cal K on the cut has no line integral"))
    return shift > 0 ? -0.5 : 0.5
end

"""
    free_energy_jks53(st, grids, beta, U; mu=0.0, H=0.0, form=:cut) -> ComplexF64

`f = -log(Lambda) / beta` from eq (56).

`form = :cut` uses the second equality, which lives entirely on `[-1, 1]`:

    log Lambda = -int_{-1}^{1} cal K log[(1+c^+ +cbar^+)(1+c^- +cbar^-)/(cbar^+ cbar^-)] dx
                 -int_{-1}^{1} cal K_{-2 gamma} Dlog C dx
                 -beta U / 4

`form = :wide` uses the first equality, which mixes a `[-1, 1]` integral with a
whole-line one over the b channel:

    log Lambda = -int_{-1}^{1} cal K log[(1+c^+ +cbar^+)(1+c^- +cbar^-)] dx
                 +int_{-inf}^{inf} [(cal K_{a-2g} - cal K_{a}) log B^+
                                  - (cal K_{-a-2g} - cal K_{-a}) log B^-] dx
                 +beta (mu + H/2 + U/4)

The two are independent routes to the same number — they share no quadrature and
the second is the only one that touches `b` — so their agreement is a real check
and not a restatement. `imag(f)` is zero by construction for a Hermitian
Hamiltonian, in either form.

The first `cal K` carries the `1/sqrt(1-x^2)` endpoint singularity, which the
Chebyshev-Gauss weight of [`JKSGrids53`](@ref) absorbs exactly; the shifted
`cal K` are regular on `[-1, 1]` and use the plain weights.
"""
function free_energy_jks53(
    st::JKSState53,
    grids::JKSGrids53,
    beta::Real,
    U::Real;
    mu::Real=U / 2,
    H::Real=0.0,
    form::Symbol=:cut,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    form in (:cut, :wide) || throw(ArgumentError("form must be :cut or :wide"))
    g, Nn = grids.gamma, grids.Nn

    # int_{-1}^{1} cal K * u dx with cal K = -1/(2 pi sqrt(1-x^2)):
    # the Chebyshev-Gauss rule gives (pi/Nn) * sum u_j for the 1/sqrt weight, so
    #     int cal K u dx = -(1/(2 pi)) (pi/Nn) sum u_j = -(1/(2 Nn)) sum u_j.
    numer = @. 1 + st.c_plus + st.cbar_plus
    numer2 = @. 1 + st.c_minus + st.cbar_minus
    if form === :cut
        u = @. log(numer * numer2 / (st.cbar_plus * st.cbar_minus))
    else
        u = @. log(numer * numer2)
    end
    int_cut = -(-sum(u) / (2 * Nn))          # -int cal K u dx

    if form === :cut
        dlog_C = log.(1 .+ st.c_plus) .- log.(1 .+ st.c_minus)
        acc = zero(ComplexF64)
        for j in 1:Nn
            acc += grids.wn[j] * jks53_script_K(grids.xn[j], -2g) * dlog_C[j]
        end
        log_Lambda = int_cut - acc - beta * U / 4
    else
        a = grids.alpha
        log_Bp = log.(1 .+ st.b_plus)
        log_Bm = log.(1 .+ st.b_minus)
        # Tail-corrected: each of the four cal K integrals is taken with its
        # exact whole-line value, the discrete part subtracted, and the
        # remainder charged to the asymptote of log B (see jks53_apply_wide).
        Bp_inf = (log_Bp[1] + log_Bp[end]) / 2
        Bm_inf = (log_Bm[1] + log_Bm[end]) / 2
        acc = zero(ComplexF64)
        sp = zero(ComplexF64)
        sm = zero(ComplexF64)
        for j in 1:grids.Nw
            x = grids.xw[j]
            kp = jks53_script_K(x, a - 2g) - jks53_script_K(x, a)
            km = jks53_script_K(x, -a - 2g) - jks53_script_K(x, -a)
            acc += grids.ww[j] * (kp * log_Bp[j] - km * log_Bm[j])
            sp += grids.ww[j] * kp
            sm += grids.ww[j] * km
        end
        Ip = jks53_script_K_line_integral(a - 2g) - jks53_script_K_line_integral(a)
        Im_ = jks53_script_K_line_integral(-a - 2g) - jks53_script_K_line_integral(-a)
        acc += Bp_inf * (Ip - sp) - Bm_inf * (Im_ - sm)
        log_Lambda = int_cut + acc + beta * (mu + H / 2 + U / 4)
    end
    return -log_Lambda / beta
end

# ═══════════════════════════════════════════════════════════════════════════════
# Initialisation — the exact beta -> 0 solution
# ═══════════════════════════════════════════════════════════════════════════════

"""
    init_jks53!(st) -> JKSState53

Set the state to the **exact** `beta -> 0` solution of eq (53):

    b^{+-} = 1,    c^{+-} = cbar^{+-} = 1/2

Derivation. As `beta -> 0` every driving term of eq (53) vanishes, so with
constant unknowns and no jump the equations close on the whole-line integrals of
the kernels, each of which is a residue:

    int K1bar(x + i d) dx = 1 (d > 0),  0 (-2 gamma < d < 0)
    int K1(x + i d) dx    = 0 (d > 0),  1 (-2 gamma < d < 0)
    int K2(x + i d) dx    = 1           (|d| < 2 gamma)

The c equation then reads `log c = -log Bbar`, i.e. `c = b/(1+b)`; the cbar
equation reads `log cbar = -log B`, i.e. `cbar = 1/(1+b)`; and the b equation
reads `log b^+ = -log B^+ + log B^-`, which forces `b = 1`. Hence
`c = cbar = 1/2`.

This is also where the eq (47) solver started from a different point:
`init_atomic_limit!` sets `c = cbar = exp(-beta U/2) -> 1`, and
`(1 + c + cbar)/cbar` is then 3 rather than 4, so its `beta -> 0` free energy
tends to `-log(3)/beta` instead of `-log(4)/beta`.
"""
function init_jks53!(st::JKSState53)
    fill!(st.b_plus, ComplexF64(1))
    fill!(st.b_minus, ComplexF64(1))
    fill!(st.c_plus, ComplexF64(0.5))
    fill!(st.c_minus, ComplexF64(0.5))
    fill!(st.cbar_plus, ComplexF64(0.5))
    fill!(st.cbar_minus, ComplexF64(0.5))
    return st
end

"Allocate a state on `grids` and set it to the exact `beta -> 0` solution."
init_jks53(grids::JKSGrids53) = init_jks53!(JKSState53(grids))

# ═══════════════════════════════════════════════════════════════════════════════
# Newton solve
# ═══════════════════════════════════════════════════════════════════════════════

"Result of an eq (53) solve: the state, iteration count, residual norm, flag."
struct JKSSolution53
    state::JKSState53
    iterations::Int
    residual::Float64
    converged::Bool
end

"""
    jks53_jacobian_fd(st, grids, ops, beta, U, mu; H=0.0, h=1e-7) -> Matrix

Jacobian of [`jks53_residual`](@ref) with respect to the packed log-unknowns, by
central differences. `2 Nw + 4 Nn` residual evaluations.

Deliberately the only Jacobian here: the eq (47) module gained an analytic one
(#797) only after its residual was pinned, and the same order applies. Getting
the block structure of a six-channel system right is a separate change from
getting the equations right, and one of them has to be trusted while the other
is being verified.
"""
function jks53_jacobian_fd(
    st::JKSState53,
    grids::JKSGrids53,
    ops::JKSOperators53,
    beta::Real,
    U::Real,
    mu::Real;
    H::Real=0.0,
    h::Real=1e-7,
)
    v = jks53_pack(st)
    n = length(v)
    J = zeros(ComplexF64, n, n)
    for k in 1:n
        vp = copy(v)
        vp[k] += h
        vm = copy(v)
        vm[k] -= h
        rp = jks53_residual(jks53_unpack(vp, grids), grids, ops, beta, U, mu; H=H)
        rm = jks53_residual(jks53_unpack(vm, grids), grids, ops, beta, U, mu; H=H)
        @views J[:, k] .= (rp .- rm) ./ (2h)
    end
    return J
end

"""
    solve_jks53_newton(grids, beta, U, mu; ...) -> JKSSolution53

Damped Newton on eq (53) from `st0` (default: the exact `beta -> 0` solution).

The damping is a plain backtracking line search on `norm(residual)`; a step is
halved up to `max_backtrack` times before the iteration is declared stuck. No
mixing parameter is exposed — a Newton step that needs mixing to converge is
reporting that the starting point is too far, which the `beta`-continuation in
[`solve_jks53_continuation`](@ref) is the right answer to.
"""
function solve_jks53_newton(
    grids::JKSGrids53,
    beta::Real,
    U::Real,
    mu::Real;
    H::Real=0.0,
    st0::Union{Nothing,JKSState53}=nothing,
    ops::Union{Nothing,JKSOperators53}=nothing,
    tol::Real=1e-10,
    maxiter::Int=50,
    max_backtrack::Int=12,
)
    beta > 0 || throw(DomainError(beta, "beta must be > 0"))
    O = ops === nothing ? JKSOperators53(grids) : ops
    st = st0 === nothing ? init_jks53(grids) : st0
    v = jks53_pack(st)
    r = jks53_residual(jks53_unpack(v, grids), grids, O, beta, U, mu; H=H)
    nr = norm(r)
    iters = 0
    for it in 1:maxiter
        iters = it
        nr < tol && break
        J = jks53_jacobian_fd(jks53_unpack(v, grids), grids, O, beta, U, mu; H=H)
        local dv
        try
            dv = J \ (-r)
        catch e
            e isa Union{SingularException,LAPACKException} || rethrow()
            return JKSSolution53(jks53_unpack(v, grids), iters, nr, false)
        end
        any(!isfinite, dv) && return JKSSolution53(jks53_unpack(v, grids), iters, nr, false)
        t = 1.0
        accepted = false
        for _ in 1:max_backtrack
            vt = v .+ t .* dv
            rt = jks53_residual(jks53_unpack(vt, grids), grids, O, beta, U, mu; H=H)
            if isfinite(norm(rt)) && norm(rt) < nr
                v, r, nr = vt, rt, norm(rt)
                accepted = true
                break
            end
            t /= 2
        end
        accepted || break
    end
    return JKSSolution53(jks53_unpack(v, grids), iters, nr, nr < tol)
end

"""
    solve_jks53_continuation(grids, beta_target, U, mu; ...) -> JKSSolution53

Walk `beta` geometrically from `beta_start` up to `beta_target`, seeding each
solve with the previous solution. The exact `beta -> 0` solution
([`init_jks53!`](@ref)) makes the first step essentially free, which is what
lets the ladder start from a point that is known rather than guessed.

On a failed step the ratio is pulled toward 1 and the step retried; the walk
gives up once the ratio is within `ratio_floor` of 1.
"""
function solve_jks53_continuation(
    grids::JKSGrids53,
    beta_target::Real,
    U::Real,
    mu::Real;
    H::Real=0.0,
    beta_start::Real=1e-6,
    ratio::Real=1.6,
    ratio_floor::Real=1.0005,
    tol::Real=1e-10,
    maxiter::Int=50,
    maxsteps::Int=400,
)
    beta_target > 0 || throw(DomainError(beta_target, "beta_target must be > 0"))
    O = JKSOperators53(grids)
    st = init_jks53(grids)
    b = min(Float64(beta_start), Float64(beta_target))
    sol = solve_jks53_newton(grids, b, U, mu; H=H, st0=st, ops=O, tol=tol, maxiter=maxiter)
    sol.converged || return sol
    st = sol.state
    total = sol.iterations
    rat = Float64(ratio)
    for _ in 1:maxsteps
        b >= beta_target * (1 - 1e-12) && break
        b_try = min(b * rat, Float64(beta_target))
        s = solve_jks53_newton(
            grids, b_try, U, mu; H=H, st0=st, ops=O, tol=tol, maxiter=maxiter
        )
        total += s.iterations
        if s.converged
            b, st = b_try, s.state
            rat = min(rat * 1.1, Float64(ratio))
        else
            rat = 1 + (rat - 1) / 2
            rat < ratio_floor && return JKSSolution53(st, total, s.residual, false)
        end
    end
    r = jks53_residual(st, grids, O, beta_target, U, mu; H=H)
    return JKSSolution53(st, total, norm(r), b >= beta_target * (1 - 1e-12))
end

"""
    hubbard1d_jks53_free_energy(t, U, mu, beta; H=0.0, Nw=96, Nn=48, x_max=32.0, form=:cut)

Free energy per site from the eq (53) six-unknown NLIE. `t` is accepted for
signature compatibility with the eq (47) route and must be 1 — the paper sets
`t = 1` throughout and `gamma = U/4` carries the only coupling.

`mu` is in **the paper's convention**, whose Coulomb term is symmetric
(`U (n_down - 1/2)(n_up - 1/2)`), so **half filling is `mu = 0`**. The plain
`U n_up n_down` convention used by this package's Lieb-Wu rows is reached by

    f_plain(mu) = f_paper(mu - U/2) - U/4

Not wired to the `Hubbard1D/FreeEnergy/Infinite` registry row: that still routes
through the eq (47) implementation (#798). Switching it over needs the conversion
above validated against the plain-convention closed form, which is a separate
change from getting eq (53) right.

Returns `NaN` if the `beta`-continuation does not reach `beta`.
"""
function hubbard1d_jks53_free_energy(
    t::Real,
    U::Real,
    mu::Real,
    beta::Real;
    H::Real=0.0,
    Nw::Int=96,
    Nn::Int=48,
    x_max::Real=32.0,
    form::Symbol=:cut,
)
    isapprox(t, 1) || throw(DomainError(t, "the JKS route fixes t = 1"))
    gamma = U / 4
    grids = JKSGrids53(Nw, Nn, gamma, 2 * gamma / 3; x_max=x_max)
    sol = solve_jks53_continuation(grids, beta, U, mu; H=H)
    sol.converged || return NaN
    return real(free_energy_jks53(sol.state, grids, beta, U; mu=mu, H=H, form=form))
end
