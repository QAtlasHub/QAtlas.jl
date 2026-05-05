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
# The generic Yang-Yang integral representation for general -1 < Δ < 1
# is deferred to a follow-up commit (it has several equivalent forms in
# the literature that differ by nontrivial substitutions; verifying any
# one of them against both boundary limits (Hulthén at Δ=1 and the
# free-fermion result at Δ=0) requires careful derivation).  For the
# current scope we expose exact values only at the three canonical
# points:
#
#   Δ = -1:  -1/4      (FM saturated)
#   Δ =  0:  -1/π      (XX / free fermion)
#   Δ =  1:  1/4 - ln 2 (AF Heisenberg, Hulthén 1938)
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
# Exact values at three canonical points only — see module docstring
# for the rationale (the general-Δ Yang-Yang integral is deferred to a
# follow-up commit; several equivalent forms in the literature differ
# by nontrivial substitutions and require careful verification against
# both boundary limits).

_xxz1d_energy_free_fermion(J::Float64)::Float64 = -J / π
_xxz1d_energy_heisenberg_af(J::Float64)::Float64 = J * (0.25 - log(2.0))
_xxz1d_energy_heisenberg_fm(J::Float64)::Float64 = -J / 4

native_energy_granularity(::XXZ1D, ::OBC) = :total
native_energy_granularity(::XXZ1D, ::Infinite) = :per_site

"""
    fetch(model::XXZ1D, ::Energy{:per_site}, ::Infinite) -> Float64

Ground-state energy **per site** of the infinite XXZ chain in units of
the Hamiltonian `J`.  Currently exposes the three canonical values:

- `Δ = -1`  →  `-J/4`             (isotropic FM, saturated)
- `Δ =  0`  →  `-J/π`             (XX, free fermion)
- `Δ =  1`  →  `J (1/4 - ln 2)`   (AF Heisenberg, Hulthén 1938)

For every other `Δ` a warning is emitted and `NaN` is returned — the
general-`Δ` Bethe-ansatz integral is tracked as a v0.13 follow-up.
"""
function fetch(model::XXZ1D, ::Energy{:per_site}, ::Infinite; kwargs...)
    J, Δ = model.J, model.Δ
    if isapprox(Δ, 0.0; atol=1e-12)
        return _xxz1d_energy_free_fermion(J)
    elseif isapprox(Δ, 1.0; atol=1e-12)
        return _xxz1d_energy_heisenberg_af(J)
    elseif isapprox(Δ, -1.0; atol=1e-12)
        return _xxz1d_energy_heisenberg_fm(J)
    else
        @warn "XXZ1D Energy: general-Δ Bethe ansatz not yet implemented; " *
            "only Δ ∈ {-1, 0, 1} are exposed in this release." Δ = Δ
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
