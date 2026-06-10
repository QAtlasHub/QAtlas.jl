# universalities/behaviour/conformal_2plus1d.jl
#
# Universal BEHAVIOUR for 2+1D Conformal Field Theories (CFTs):
# 1. Sphere free energy F = -ln |Z(S^3)| (F-theorem)
# 2. Corner entanglement coefficient c(theta) / prefactor sigma.
#
# References:
#   - I. R. Klebanov, S. S. Pufu, and B. R. Safdi, JHEP 11, 038 (2011) — F-theorem.
#   - F. Kos, D. Poland, D. Simmons-Duffin, and A. Vichi, JHEP 08, 036 (2016) — 3D Ising/O(N) bootstrap.
#   - S. M. Chester et al., Phys. Rev. D 101, 105013 (2020) — O(2) stress tensor bootstrap.

# ─── SphereFreeEnergy ────────────────────────────────────────────────────────

raw"""
    fetch(::Universality{C}, ::SphereFreeEnergy, ::Infinite; d::Int = 3, kwargs...) -> Float64

Return the universal sphere free energy $F = -\ln |Z(S^3)|$ of the 2+1D CFT
for the universality class `C` in spatial dimension `d=3`.

Supported classes:
- `:Ising`       : 3D Ising CFT (F ≈ 0.0612)
- `:XY`          : 3D O(2) XY CFT (F ≈ 0.121)
- `:Heisenberg`  : 3D O(3) Heisenberg CFT (F ≈ 0.180)

Reference: I. R. Klebanov, S. S. Pufu, and B. R. Safdi, JHEP 11, 038 (2011).
"""
function fetch(
    u::Universality{C}, ::SphereFreeEnergy, ::Infinite; d::Int=3, kwargs...
) where {C}
    d == 3 || throw(ArgumentError("SphereFreeEnergy: only d=3 (2+1D CFT) is supported; got d=$d"))
    
    if C === :Ising
        return 0.0612  # conformal bootstrap / fuzzy sphere consensus
    elseif C === :XY
        return 0.121   # 3D O(2) fixed point
    elseif C === :Heisenberg
        return 0.180   # 3D O(3) fixed point
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
- If `theta::Real` is provided, returns the leading-order smooth-limit corner contribution `c(θ)`.

Supported classes:
- `:Ising`       : 3D Ising model CFT (σ ≈ 0.003697)
- `:XY`          : 3D O(2) XY model CFT (σ ≈ 0.007377)
- `:Heisenberg`  : 3D O(3) Heisenberg model CFT (σ ≈ 0.01105)

Reference: S. M. Chester et al., Phys. Rev. D 101, 105013 (2020).
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
        0.0036974  # Bootstrap C_T / C_T,free ≈ 0.94653
    elseif C === :XY
        0.0073773  # Bootstrap C_T / C_T,free ≈ 2 * 0.9443
    elseif C === :Heisenberg
        0.011050   # Bootstrap C_T / C_T,free ≈ 3 * 0.9429
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
