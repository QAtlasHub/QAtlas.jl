# ─────────────────────────────────────────────────────────────────────────────
# XXZ chain (1D, spin-1/2) — Bethe-ansatz ground-state energy and
# Luttinger-liquid parameters in the critical regime.
#
# Hamiltonian (Takahashi / standard convention):
#
#   H = J Σ_i [ S^x_i S^x_{i+1} + S^y_i S^y_{i+1} + Δ S^z_i S^z_{i+1} ]
#
#   (spin-1/2, `J > 0` is the antiferromagnetic sign convention).
#
# In Pauli notation this is equivalently
#
#   H = (J/4) Σ_i [ σ^x σ^x + σ^y σ^y + Δ σ^z σ^z ]
#
# Known limiting values (all per site, units of J):
#
#   Δ =  1  (isotropic AF / Heisenberg):  e₀/J = 1/4 − ln 2 ≈ −0.4431
#                                         (Hulthén 1938)
#   Δ =  0  (XX / free fermion):          e₀/J = −1/π       ≈ −0.3183
#   Δ = −1  (isotropic FM):               e₀/J = −1/4
#
# General -1 < Δ < 1 is now covered by the Yang–Yang single-integral
# form (see `XXZ_bethe.jl` and the existing derivation pages
# `docs/src/calc/bethe-ansatz-heisenberg-e0.md` (γ = 0) and
# `docs/src/calc/xxz-luttinger-parameters.md` (anisotropic kernels)).
# After Fourier-transforming the linear Bethe-density equation
# ρ + a₂ ⋆ ρ = a₁ with a_n the Takahashi anisotropic kernels, sum-to-
# product collapses
#
#   ρ̂(ω) = sinh((π−γ)ω/(2γ)) /
#          [sinh(πω/(2γ)) + sinh((π−2γ)ω/(2γ))]
#         = 1 / (2 cosh(ω/2))      (γ-independent),
#   ⇒ ρ(λ) = 1 / (2 cosh(πλ)),
#
# and the energy follows from the rapidity sum
#
#   e₀(Δ = cos γ) = (J cos γ)/4
#       − J sin² γ ∫_{-∞}^{∞} ρ(λ) dλ / (cosh(2γλ) − cos γ).
#
# All three closed-form points are kept as fast-paths in the dispatch
# below; the QuadGK integral is invoked only for Δ ∈ (−1, 1) \ {0}.
# The gapped regime |Δ| > 1 still emits a warning and returns NaN —
# the Orbach / Walker–Smith series form for that regime is a separate
# follow-up.
#
# Other observables (central charge, Luttinger parameter, Luttinger
# velocity) are implemented analytically across the full critical
# regime -1 < Δ ≤ 1.
#
# References:
#
#   - L. Hulthén (1938) Arkiv Mat. Astron. Fysik 26A, 1.
#   - C. N. Yang and C. P. Yang (1966) Phys. Rev. 150, 321.
#   - M. Takahashi, "Thermodynamics of One-Dimensional Solvable Models",
#     Cambridge University Press (1999).
#   - T. Giamarchi, "Quantum Physics in One Dimension",
#     Oxford University Press (2004), §6 for Luttinger parameter & velocity.
# ─────────────────────────────────────────────────────────────────────────────

"""
    XXZ1D(; J::Real = 1.0, Δ::Real = 0.0) <: AbstractQAtlasModel

Spin-1/2 XXZ chain

    H = J Σ_i [ S^x_i S^x_{i+1} + S^y_i S^y_{i+1} + Δ S^z_i S^z_{i+1} ]

Convention: `J > 0` is antiferromagnetic.  `Δ = 1` is the isotropic
Heisenberg AF point, `Δ = 0` is the XX (free-fermion) point, `Δ = -1`
is the isotropic ferromagnet.  For `|Δ| < 1` the chain is critical
(Luttinger liquid, central charge `c = 1`).
"""
struct XXZ1D <: AbstractQAtlasModel
    J::Float64
    Δ::Float64
end
XXZ1D(; J::Real=1.0, Δ::Real=0.0) = XXZ1D(Float64(J), Float64(Δ))

# ── Ground-state energy per site (infinite chain) ──────────────────────
#
# Three closed-form fast-paths at Δ ∈ {-1, 0, 1}; general -1 < Δ < 1
# delegates to `_xxz1d_energy_yang_yang` (see `XXZ_bethe.jl`) which
# evaluates the Yang–Yang single integral via QuadGK.  The gapped
# regime |Δ| > 1 still warns and returns NaN.

_xxz1d_energy_free_fermion(J::Float64)::Float64 = -J / π
_xxz1d_energy_heisenberg_af(J::Float64)::Float64 = J * (0.25 - log(2.0))
_xxz1d_energy_heisenberg_fm(J::Float64)::Float64 = -J / 4

native_energy_granularity(::XXZ1D, ::OBC) = :total
native_energy_granularity(::XXZ1D, ::Infinite) = :per_site

# fetch(::XXZ1D, ::Energy{:per_site}, ::Infinite; ...) is defined in
# XXZ_xx_infinite.jl, which extends the ground-state branch with a
# finite-T (β kwarg) free-fermion path at Δ = 0.  The ground-state
# logic for Δ ∈ {-1, 0, 1} is preserved bit-for-bit there; general-Δ
# Bethe-ansatz remains a v0.13 follow-up (issue #108).
"""
    fetch(model::XXZ1D, ::Energy{:per_site}, ::Infinite) -> Float64

Ground-state energy **per site** of the infinite XXZ chain in units of
the Hamiltonian `J`, at zero temperature.

# Coverage

- **Critical / gapless regime** `-1 ≤ Δ ≤ 1`: closed form for the
  three canonical points and the Yang–Yang single integral elsewhere:

      Δ = -1:  -J/4                 (isotropic FM, saturated)
      Δ =  0:  -J/π                 (XX, free fermion)
      Δ =  1:  J (1/4 - ln 2)       (AF Heisenberg, Hulthén 1938)
      otherwise (γ = arccos Δ):
        e₀(Δ) = (J cos γ)/4 − J sin² γ
                · ∫_{-∞}^{∞} dλ / [2 cosh(πλ)·(cosh(2γλ) − cos γ)]

  Returned to ≈ 1e-12 relative accuracy via adaptive QuadGK.

- **Gapped regime** `|Δ| > 1`: the Bethe ansatz takes a different
  series form (Orbach 1958 / Walker 1959 / Yang–Yang 1966 III); the
  closed-form path here emits a warning and returns `NaN`.  Use OBC
  dense ED at small `N` for a finite-size reference in that regime.
"""
function fetch(model::XXZ1D, ::Energy{:per_site}, ::Infinite; kwargs...)
    J, Δ = model.J, model.Δ
    if isapprox(Δ, 0.0; atol=1e-12)
        return _xxz1d_energy_free_fermion(J)
    elseif isapprox(Δ, 1.0; atol=1e-12)
        return _xxz1d_energy_heisenberg_af(J)
    elseif isapprox(Δ, -1.0; atol=1e-12)
        return _xxz1d_energy_heisenberg_fm(J)
    elseif -1.0 < Δ < 1.0
        return _xxz1d_energy_yang_yang(J, Δ)
    else
        @warn "XXZ1D Energy: gapped regime |Δ| > 1 not yet implemented; " *
            "use OBC dense ED at small N for a finite-size reference." Δ = Δ
        return NaN
    end
end

"""
    fetch(model::XXZ1D, ::GroundStateEnergyDensity, ::Infinite) -> Float64

Alias for [`fetch(::XXZ1D, ::Energy, ::Infinite)`](@ref) kept so that
the `GroundStateEnergyDensity` quantity — already exported by
`Heisenberg.jl` — works uniformly across 1D Bethe-ansatz chains.
"""
function fetch(model::XXZ1D, ::GroundStateEnergyDensity, ::Infinite; kwargs...)
    return fetch(model, Energy(), Infinite(); kwargs...)
end

# ── Central charge & Luttinger-liquid parameters (critical regime) ─────

"""
    fetch(model::XXZ1D, ::CentralCharge, ::Infinite) -> Float64

Central charge of the XXZ chain:

- `-1 < Δ < 1`  → `c = 1` (Luttinger liquid)
- otherwise     → `NaN` (non-critical)
"""
function fetch(model::XXZ1D, ::CentralCharge, ::Infinite; kwargs...)
    if -1.0 < model.Δ < 1.0
        return 1.0
    end
    @warn "XXZ1D CentralCharge is only defined in the critical regime -1 < Δ < 1." Δ =
        model.Δ
    return NaN
end

"""
    fetch(model::XXZ1D, ::LuttingerParameter, ::Infinite) -> Float64

Luttinger-liquid parameter `K = π / (2(π − γ))`, with `γ = arccos(Δ)`,
valid for `-1 < Δ ≤ 1`.

Canonical values:
- `Δ = 0` (XX free fermion) → `K = 1`
- `Δ = 1` (AF Heisenberg)   → `K = 1/2`
- `Δ → -1` (FM boundary)    → `K → ∞`
"""
function fetch(model::XXZ1D, ::LuttingerParameter, ::Infinite; kwargs...)
    Δ = model.Δ
    if -1.0 < Δ ≤ 1.0
        γ = acos(Δ)
        return π / (2 * (π - γ))
    end
    @warn "XXZ1D LuttingerParameter is only defined for -1 < Δ ≤ 1." Δ = Δ
    return NaN
end

"""
    fetch(model::XXZ1D, ::LuttingerVelocity, ::Infinite) -> Float64
    fetch(model::XXZ1D, ::SpinWaveVelocity,   ::Infinite) -> Float64

Sound velocity of the low-energy Luttinger-liquid mode,

    u(Δ) = J · (π/2) · sin(γ)/γ,   γ = arccos(Δ).

Canonical values:
- `Δ = 0` (XX)       → `u = J`         (= free-fermion v_F)
- `Δ = 1` (AF)       → `u = (π/2) J`  (des Cloizeaux-Pearson)

`SpinWaveVelocity` dispatches here via the `const SpinWaveVelocity =
LuttingerVelocity` type alias (both are the same physical quantity for
1D critical spin chains).
"""
function fetch(model::XXZ1D, ::LuttingerVelocity, ::Infinite; kwargs...)
    J, Δ = model.J, model.Δ
    if -1.0 < Δ ≤ 1.0
        γ = acos(Δ)
        # sin(γ)/γ has a removable singularity at γ = 0 (Heisenberg AF);
        # the naive ratio is fine in Float64 for any γ > 0 that
        # corresponds to Δ < 1, and at Δ = 1 we take the limit.
        return J * (π / 2) * (isapprox(γ, 0.0; atol=1e-12) ? 1.0 : sin(γ) / γ)
    end
    @warn "XXZ1D LuttingerVelocity is only defined for -1 < Δ ≤ 1." Δ = Δ
    return NaN
end

# ── Finite-temperature, finite-N OBC via dense ED ──────────────────────
#
# General-Δ finite-T XXZ has no closed-form formula like the TFIM
# BdG decomposition; the thermal Bethe ansatz gives the thermodynamic
# limit only.  For MPS-benchmark use cases the working scale is
# small-N (say N ≤ 10) where dense ED of the 2^N × 2^N spin-1/2
# Hilbert space is trivially cheap and gives a *finite-N exact*
# reference.  We implement that path here; closed-form thermal at
# Δ = 0 (XX / free fermion) and TBA at general Δ in the
# thermodynamic limit are separate follow-ups.

"""
    _xxz1d_hamiltonian_matrix(model::XXZ1D, N::Int) -> Matrix{ComplexF64}

Assemble the `2^N × 2^N` OBC Hamiltonian

    H = J Σᵢ [ Sˣ_i Sˣ_{i+1} + Sʸ_i Sʸ_{i+1} + Δ Sᶻ_i Sᶻ_{i+1} ]
      = (J/4) Σᵢ [ σˣ σˣ + σʸ σʸ + Δ σᶻ σᶻ ]

via explicit tensor products built from the Pauli primitives in
`src/core/dense_ed.jl`.  Capped by `_MAX_ED_SITES`.
"""
function _xxz1d_hamiltonian_matrix(model::XXZ1D, N::Int)
    N ≥ 2 || throw(ArgumentError("XXZ1D OBC chain needs N ≥ 2 (got N = $N)"))
    J, Δ = model.J, model.Δ
    D = 2^N
    H = zeros(ComplexF64, D, D)
    prefac = J / 4
    for i in 1:(N - 1)
        H .+= prefac .* _pauli_string(N, i => _σx, i + 1 => _σx)
        H .+= prefac .* _pauli_string(N, i => _σy, i + 1 => _σy)
        H .+= (prefac * Δ) .* _pauli_string(N, i => _σz, i + 1 => _σz)
    end
    return H
end

"""
    fetch(model::XXZ1D, ::Energy{:total}, ::OBC; beta) -> Float64

**Total** thermal energy `⟨H⟩_β` for the spin-½ OBC chain at finite size,
computed by dense ED.  Works for any `Δ` and any `N ≤ $(_MAX_ED_SITES)`.
Intended as a reference for MPS thermal methods (TPQMPS / Purification /
METTS).

```
    ⟨H⟩_β = Tr(H exp(-βH)) / Tr(exp(-βH))
```

Convention matches [`fetch(::TFIM, ::Energy, ::OBC)`](@ref): finite-size
boundary conditions return total energy; only `Infinite()` returns per-site
(`⟨H⟩/N`, the only finite quantity in the thermodynamic limit).
"""
function fetch(model::XXZ1D, ::Energy{:total}, bc::OBC; beta::Real, kwargs...)
    H = _xxz1d_hamiltonian_matrix(model, bc.N)
    return _ed_thermal_energy(H, beta)
end
