# universalities/behaviour/conformal_2plus1d.jl
#
# Universal BEHAVIOUR for 2+1D Conformal Field Theories (CFTs):
# 1. Sphere free energy F = -ln |Z(S^3)| (F-theorem)
# 2. Corner entanglement coefficient c(theta) / prefactor sigma.
#
# Note: Universality{:Ising} here means the 3D Ising Wilson-Fisher fixed point
# (NOT the 2D c=1/2 BPZ Ising CFT). The d==3 guard enforces this distinction.
#
# References:
#   - I. R. Klebanov, S. S. Pufu, and B. R. Safdi, JHEP 11, 038 (2011) — F-theorem framework.
#   - S. S. Pufu, arXiv:1612.00381 (2017) — F-theorem review; numerical F values for 3D CFTs.
#   - F. Kos, D. Poland, D. Simmons-Duffin, and A. Vichi, JHEP 08, 036 (2016) — 3D Ising bootstrap C_T.
#   - S. M. Chester et al., JHEP 02, 098 (2020) — O(2)/O(3) bootstrap C_T.
#   - P. Bueno, R. C. Myers, and W. Witczak-Krempa, PRL 115, 021602 (2015) — sigma = pi^2/24 * C_T.

# ─── SphereFreeEnergy ────────────────────────────────────────────────────────

raw"""
    fetch(::Universality{C}, ::SphereFreeEnergy, ::Infinite; d::Int = 3, kwargs...) -> Float64

Return the universal sphere free energy $F = -\ln Z(S^3)$ of the 2+1D CFT
for the universality class `C` in spatial dimension `d=3`.

Supported classes:
- `:Ising`       : 3D Ising CFT (F ≈ 0.0612)
- `:XY`          : 3D O(2) XY CFT (F ≈ 0.121)
- `:Heisenberg`  : 3D O(3) Heisenberg CFT (F ≈ 0.180)

All values satisfy F < N * F_free_scalar (F-theorem: F decreases along RG flows).
The free real scalar has F_free = ln(2)/8 - 3ζ(3)/(16π²) ≈ 0.0638.

References:
- Framework: Klebanov, Pufu, Safdi, JHEP 11, 038 (2011).
- Numerical values: Pufu, arXiv:1612.00381 (2017), Table 1.
"""
function fetch(
    u::Universality{C}, ::SphereFreeEnergy, ::Infinite; d::Int=3, kwargs...
) where {C}
    d == 3 || throw(ArgumentError("SphereFreeEnergy: only d=3 (2+1D CFT) is supported; got d=$d"))

    if C === :Ising
        return 0.0612  # 3D Ising CFT; conformal bootstrap / fuzzy sphere consensus (Pufu 2017)
    elseif C === :XY
        return 0.121   # 3D O(2) Wilson-Fisher fixed point (Pufu 2017)
    elseif C === :Heisenberg
        return 0.180   # 3D O(3) Wilson-Fisher fixed point (Pufu 2017)
    else
        throw(
            ArgumentError(
                "Universality{:$C} SphereFreeEnergy: only :Ising, :XY, and :Heisenberg are supported.",
            ),
        )
    end
end

# ─── CornerEntanglementCoefficient ───────────────────────────────────────────

raw"""
    fetch(::Universality{C}, ::CornerEntanglementCoefficient, ::Infinite;
          d::Int = 3, theta::Union{Real, Nothing} = nothing, kwargs...) -> Float64

Return the universal corner coefficient in the bipartite entanglement entropy of the
2+1D CFT for class `C` in spatial dimension `d=3`.

For a boundary corner of angle $\theta$ close to $\pi$, the corner contribution is:

    c(θ) ≈ σ * (π - θ)^2

where $\sigma = \frac{\pi^2}{24} C_T$ is the smooth-limit prefactor.

- If `theta` is `nothing` (default), returns the prefactor `σ`.
- If `theta::Real` is provided, returns the leading-order smooth-limit approximation
  `c(θ) ≈ σ*(π-θ)²`. **This approximation is accurate only near θ → π (nearly smooth
  boundary). For sharp corners (small θ) the full angular function c(θ) deviates
  significantly from this quadratic form and cannot be captured by this coefficient alone.**

Supported classes:
- `:Ising`       : 3D Ising model CFT (σ ≈ 0.0036974, from C_T/C_T,free ≈ 0.94653)
- `:XY`          : 3D O(2) XY model CFT (σ ≈ 0.0073773, from C_T/C_T,free ≈ 2 × 0.9443)
- `:Heisenberg`  : 3D O(3) Heisenberg model CFT (σ ≈ 0.0110496, from C_T/C_T,free ≈ 3 × 0.9429)

References:
- Formula σ = π²/24 · C_T: Bueno, Myers, Witczak-Krempa, PRL 115, 021602 (2015).
- C_T (Ising): Kos, Poland, Simmons-Duffin, Vichi, JHEP 08, 036 (2016).
- C_T (O(N)): Chester et al., JHEP 02, 098 (2020).
"""
function fetch(
    u::Universality{C},
    ::CornerEntanglementCoefficient,
    ::Infinite;
    d::Int=3,
    theta::Union{Real,Nothing}=nothing,
    kwargs...,
) where {C}
    d == 3 || throw(
        ArgumentError(
            "CornerEntanglementCoefficient: only d=3 (2+1D CFT) is supported; got d=$d",
        ),
    )

    σ = if C === :Ising
        0.0036974   # π²/24 * 0.94653 * 3/(32π²); Kos et al. 2016
    elseif C === :XY
        0.0073773   # π²/24 * 2 * 0.9443 * 3/(32π²); Chester et al. 2020
    elseif C === :Heisenberg
        0.0110496   # π²/24 * 3 * 0.9429 * 3/(32π²); Chester et al. 2020
    else
        throw(
            ArgumentError(
                "Universality{:$C} CornerEntanglementCoefficient: only :Ising, :XY, and :Heisenberg are supported.",
            ),
        )
    end

    if theta === nothing
        return σ
    else
        0 ≤ theta ≤ π || throw(
            ArgumentError(
                "CornerEntanglementCoefficient: theta must satisfy 0 ≤ theta ≤ π; got theta=$theta.",
            ),
        )
        return σ * (π - theta)^2
    end
end
