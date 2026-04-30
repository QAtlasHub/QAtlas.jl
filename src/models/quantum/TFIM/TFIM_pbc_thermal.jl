# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — PBC finite-N free-fermion thermodynamics
#
# Hamiltonian (PBC, N sites):
#   H = -J Σᵢ σᶻᵢσᶻᵢ₊₁  -  h Σᵢ σˣᵢ        with site N+1 ≡ 1
#
# Jordan-Wigner with periodic spin boundary maps to a free-fermion bilinear,
# but the fermion BC depends on Z₂ parity P = ∏ᵢ σˣᵢ:
#
#   P = +1  (even)  →  fermion anti-periodic (Neveu-Schwarz / NS)
#                       k_NS = (2j-1)π/N,  j = 1..N
#   P = -1  (odd)   →  fermion periodic     (Ramond / R)
#                       k_R  = 2jπ/N,       j = 0..N-1   (includes k = 0, π)
#
# In each sector the BdG dispersion is
#   Λ(k) = 2√(J² + h² - 2Jh cos k)        ≥ 0,
# and the (free-fermion + parity-projected) partition function is
#
#   Z = ½ [Z_NS^(+) + Z_NS^(-)] + ½ [Z_R^(+) − Z_R^(-)]                    (1)
#
# with the per-sector building blocks
#
#   Z_α^(+) = ∏_{k ∈ α} 2 cosh(βΛ(k)/2)            (free-fermion)
#   Z_α^(-) = ∏_{k ∈ α} 2 sinh(βΛ(k)/2)            (parity-twisted)
#
# (Lieb-Schultz-Mattis 1961 §III; Sachdev §4.2; or any chapter that walks the
# JW + Bogoliubov computation through the parity projector.)
#
# All other thermodynamic potentials follow as analytic derivatives of `log Z`.
# Each sector contribution `log Z_α^σ = Σ_k g_σ(βΛ_k/2)` with
#   g_+(x) = log(2 cosh x),     g_+'(x) = tanh x,        g_+''(x) = sech²x
#   g_-(x) = log(2 sinh x),     g_-'(x) = coth x,        g_-''(x) = -csch²x
# is closed-form mode-by-mode; the parity-projected log Z then mixes the four
# sectors with weights `w_{ασ} = s_{ασ} Z_α^σ / Z` (signs (+,+,+,−)),
# and derivatives compose by the standard cumulant rule
#   ∂_x log Z = Σ w_{ασ} ∂_x log Z_α^σ
#   ∂²_x log Z = Σ w_{ασ} ∂²_x log Z_α^σ + Σ w_{ασ} (∂_x log Z_α^σ)² − (Σ w_{ασ} ∂_x log Z_α^σ)²
# (the (R,−) weight is negative; the algebra still gives the exact derivative
#  of log Z because we differentiate a generic linear combination of
#  exponentials).
#
# Mass gap (PBC):
#   The Hilbert space splits by Z₂ parity P = ∏ᵢ σˣᵢ.  We compute the
#   ground-state energy in NS and R sectors, plus the lowest excited
#   state in each sector (NS preserves even parity ⇒ flip *two* modes;
#   R preserves odd parity ⇒ flip *one* mode).  The mass gap is the
#   smallest (E_excited − E_GS) across both sectors.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Native granularity at PBC
# ═══════════════════════════════════════════════════════════════════════════════

native_energy_granularity(::TFIM, ::PBC) = :per_site

# ═══════════════════════════════════════════════════════════════════════════════
# Sector momentum lists
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_pbc_momenta(N, sector::Symbol) -> Vector{Float64}

Return the N momenta for the requested fermion sector,
`sector ∈ (:NS, :R)`.

- `:NS` (anti-periodic / even parity): `k = (2j-1)π/N`, j = 1..N
- `:R`  (periodic / odd parity):       `k = 2jπ/N`,    j = 0..N-1
        (includes k = 0; k = π only if N is even.)
"""
function _tfim_pbc_momenta(N::Int, sector::Symbol)
    if sector === :NS
        return Float64[(2j - 1) * π / N for j in 1:N]
    elseif sector === :R
        return Float64[2j * π / N for j in 0:(N - 1)]
    else
        throw(ArgumentError("unknown sector $sector — expected :NS or :R"))
    end
end

# Single-particle BdG energy at momentum k.  Generic in (J, h) for downstream
# analytic derivatives.
@inline _tfim_Λ(k::Real, J::Real, h::Real) = 2 * sqrt(J^2 + h^2 - 2 * J * h * cos(k))

# ∂Λ/∂h = 4(h - J cos k)/Λ          (from Λ² = 4(J² + h² − 2Jh cos k))
@inline _tfim_dΛdh(k::Real, J::Real, h::Real) = 4 * (h - J * cos(k)) / _tfim_Λ(k, J, h)

# ∂²Λ/∂h² = 4/Λ - 16(h - J cos k)²/Λ³
@inline function _tfim_d2Λdh2(k::Real, J::Real, h::Real)
    Λ = _tfim_Λ(k, J, h)
    A = h - J * cos(k)
    return 4 / Λ - 16 * A^2 / Λ^3
end

# ═══════════════════════════════════════════════════════════════════════════════
# Per-sector log-partition and its derivatives  (closed-form, mode-by-mode)
# ═══════════════════════════════════════════════════════════════════════════════

@inline _logsinh2(x::Real) = abs(x) + log1p(-exp(-2 * abs(x)))
# (above is fine for x > 0; we never call it at x = 0 because finite β > 0
#  + Λ_k > 0 keeps βΛ_k/2 > 0 except in the critical-point + ordered-phase
#  R-sector zero-mode case at exactly h = J, which the β > 0 evaluations
#  avoid in the test slices.)

# Sector data for one (sector α, parity sign σ) pair.
# `kind = :cosh` for σ = +; `kind = :sinh` for σ = -.
# Parametric in `T` so ForwardDiff dual numbers carry through unchanged.
struct _SectorState{T<:Real}
    log_Z::T
    sign::Int                # +1 for (NS,+), (NS,-), (R,+); −1 for (R,−)
    dLdβ::T                  # ∂_β log Z_α^σ
    d2Ldβ2::T                # ∂²_β log Z_α^σ
    dLdh::T                  # ∂_h log Z_α^σ
    d2Ldh2::T                # ∂²_h log Z_α^σ
end

function _sector_state(
    N::Int, J::Real, h::Real, β::Real, sector::Symbol, kind::Symbol, sign::Int
)
    ks = _tfim_pbc_momenta(N, sector)
    T = promote_type(typeof(β), typeof(J), typeof(h))
    log_Z = zero(T)
    dLdβ = zero(T)
    d2Ldβ2 = zero(T)
    dLdh = zero(T)
    d2Ldh2 = zero(T)
    @inbounds for k in ks
        Λ = _tfim_Λ(k, J, h)
        x = β * Λ / 2
        if kind === :cosh
            log_Z += _logcosh2(x)
            tx = tanh(x)
            sech2x = 1 - tx^2
            # ∂_β log Z = (Λ/2) tanh(x);  ∂²_β = (Λ/2)² sech²(x)
            dLdβ += (Λ / 2) * tx
            d2Ldβ2 += (Λ / 2)^2 * sech2x
            # ∂_h log Z mode contribution:
            # ∂_h g_+(βΛ/2) = (β/2) tanh(x) ∂_h Λ
            dΛh = _tfim_dΛdh(k, J, h)
            dLdh += (β / 2) * tx * dΛh
            # ∂²_h g_+(βΛ/2) = (β²/4) sech²(x) (∂_h Λ)² + (β/2) tanh(x) ∂²_h Λ
            d2Λh = _tfim_d2Λdh2(k, J, h)
            d2Ldh2 += (β^2 / 4) * sech2x * dΛh^2 + (β / 2) * tx * d2Λh
        else  # :sinh
            # Critical point R-sector: at h = J the k = 0 mode (and k = π
            # for even N) gives Λ = 0, so 2 sinh(βΛ/2) = 0 and log Z_R^-
            # = -∞.  This is the *correct* mathematical statement that
            # the (R, sinh) sector vanishes at criticality; the
            # log-sum-exp downstream drops it cleanly.  We must avoid
            # the `coth(0) = ∞`, `(Λ/2)·coth(0) = 0·∞ = NaN` trap when
            # accumulating derivatives.  Send the entire sector to
            # `log_Z = -Inf` and zero derivatives so its weight is zero.
            if x <= 1e-12
                log_Z = T(-Inf)
                # Keep derivatives finite (will be multiplied by weight 0).
                dLdβ = zero(T)
                d2Ldβ2 = zero(T)
                dLdh = zero(T)
                d2Ldh2 = zero(T)
                # Bail: rest of the loop adds to a -Inf log_Z anyway.
                break
            end
            log_Z += _logsinh2(x)
            cx = coth(x)
            csch2x = cx^2 - 1
            dLdβ += (Λ / 2) * cx
            # g_-''(x) = -csch²(x), so ∂²_β log Z_α^- = -(Λ/2)² csch²(x).
            d2Ldβ2 += -(Λ / 2)^2 * csch2x
            dΛh = _tfim_dΛdh(k, J, h)
            dLdh += (β / 2) * cx * dΛh
            d2Λh = _tfim_d2Λdh2(k, J, h)
            d2Ldh2 += -(β^2 / 4) * csch2x * dΛh^2 + (β / 2) * cx * d2Λh
        end
    end
    return _SectorState(log_Z, sign, dLdβ, d2Ldβ2, dLdh, d2Ldh2)
end

# Compute the four sectors and assemble
#   Z = Z_NS,even + Z_R,even
#     = ½(Z_NS^+ + Z_NS^-) + ½(Z_R^+ + Z_R^-)
# where Z_α^+ = ∏ 2 cosh(βΛ/2) and Z_α^- = ∏ 2 sinh(βΛ/2).
#
# Convention check (cross-validated against dense ED at N = 2, 3, 4): the
# physical TFIM Hilbert space maps to (NS, even fermion parity) ∪
# (R, even fermion parity), so all four sector contributions enter with
# weight `+½` after the projector normalisation.  An earlier draft used
# the LSM-style P=+1 ↔ NS-even / P=-1 ↔ R-odd convention with sign
# (+,+,+,-); that convention is correct for the XY chain LSM analysed
# but not for our TFIM JW string convention (which puts the wrap-around
# bond's sign-flip in (-1)^(N-1) ∏ σˣ rather than ∏ σˣ).  See the test
# `test_TFIM_pbc_thermal.jl` "ED comparison N=4" testset for the numerical
# cross-validation.
function _all_sector_states(N::Int, J::Real, h::Real, β::Real)
    return (
        _sector_state(N, J, h, β, :NS, :cosh, +1),
        _sector_state(N, J, h, β, :NS, :sinh, +1),
        _sector_state(N, J, h, β, :R, :cosh, +1),
        _sector_state(N, J, h, β, :R, :sinh, +1),
    )
end

# Numerically-stable signed-log-sum-exp:
#   log( Σ s_i e^{L_i} ) using the max L_max of the positive-sign entries.
# Returns (logZ_after_½_factor, weights), where weights[i] = s_i e^{L_i}/Σ
# (one of them may be negative).
function _logZ_and_weights(states)
    Lmax = maximum(s.log_Z for s in states)
    T = typeof(Lmax)
    inner = zero(T)
    for s in states
        inner += s.sign * exp(s.log_Z - Lmax)
    end
    inner > 0 || error(
        "_logZ_and_weights: parity-projected sum non-positive ($(inner)); " *
        "this should not happen for β > 0.",
    )
    log_Z = Lmax + log(inner) - log(T(2))
    weights = ntuple(i -> states[i].sign * exp(states[i].log_Z - Lmax) / inner, 4)
    return log_Z, weights
end

"""
    _tfim_pbc_log_Z(N, J, h, β) -> Float64

`log Z` (eq. (1)) for the N-site PBC TFIM at inverse temperature β.
"""
function _tfim_pbc_log_Z(N::Int, J::Real, h::Real, β::Real)
    states = _all_sector_states(N, J, h, β)
    log_Z, _ = _logZ_and_weights(states)
    return log_Z
end

# Generic helper to take ∂_x log Z and ∂²_x log Z given the sector states +
# sector-by-sector first/second derivatives w.r.t. the variable x.
# Uses the cumulant identity (works for signed weights too — algebraic).
function _cumulant_first_second(weights::NTuple{4}, dLs::NTuple{4}, d2Ls::NTuple{4})
    # ⟨X⟩ = Σ w_i X_i,  ⟨Y⟩ = Σ w_i Y_i
    T = promote_type(typeof(weights[1]), typeof(dLs[1]), typeof(d2Ls[1]))
    mean1 = zero(T)
    mean_d2 = zero(T)
    mean_sq = zero(T)
    @inbounds for i in 1:4
        mean1 += weights[i] * dLs[i]
        mean_d2 += weights[i] * d2Ls[i]
        mean_sq += weights[i] * dLs[i]^2
    end
    return mean1, mean_d2 + (mean_sq - mean1^2)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Per-site thermodynamic potentials — closed-form, no AD
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_pbc_thermo(quantity, N, J, h, β) -> Float64

Compute one of the per-site thermodynamic potentials of the PBC TFIM.
`quantity` is a `Symbol`:

- `:free_energy`              → f = -T log Z / N
- `:energy_per_site`          → ε = -∂_β log Z / N
- `:entropy`                  → s = β(ε - f)
- `:specific_heat`            → c_v = β² ∂²_β log Z / N
- `:transverse_magnetization` → m_x = ∂_h log Z / (Nβ)
- `:transverse_susceptibility`→ χ_xx = ∂²_h log Z / (Nβ)
"""
function _tfim_pbc_thermo(quantity::Symbol, N::Int, J::Real, h::Real, β::Real)
    states = _all_sector_states(N, J, h, β)
    log_Z, weights = _logZ_and_weights(states)
    if quantity === :free_energy
        return -log_Z / (N * β)
    end
    if quantity === :energy_per_site || quantity === :entropy || quantity === :specific_heat
        dLs = ntuple(i -> states[i].dLdβ, 4)
        d2Ls = ntuple(i -> states[i].d2Ldβ2, 4)
        m1, m2 = _cumulant_first_second(weights, dLs, d2Ls)
        # ∂_β log Z = m1; ∂²_β log Z = m2
        ε = -m1 / N
        if quantity === :energy_per_site
            return ε
        elseif quantity === :entropy
            f = -log_Z / (N * β)
            return β * (ε - f)
        else
            return β^2 * m2 / N    # c_v = -β² ∂_β ε = β² ∂²_β log Z / N
        end
    end
    if quantity === :transverse_magnetization || quantity === :transverse_susceptibility
        dLs = ntuple(i -> states[i].dLdh, 4)
        d2Ls = ntuple(i -> states[i].d2Ldh2, 4)
        m1, m2 = _cumulant_first_second(weights, dLs, d2Ls)
        if quantity === :transverse_magnetization
            return m1 / (N * β)        # m_x = -∂f/∂h = (∂_h log Z) / (N β)
        else
            return m2 / (N * β)        # χ_xx = ∂m_x/∂h = (∂²_h log Z) / (N β)
        end
    end
    error("_tfim_pbc_thermo: unknown quantity :$quantity")
end

# ═══════════════════════════════════════════════════════════════════════════════
# PBC mass gap (lowest excitation across both parity sectors)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_pbc_mass_gap(N, J, h) -> Float64

Lowest single excitation energy of the N-site PBC TFIM.

The Hilbert space splits by Z₂ parity `P = ∏ᵢ σˣᵢ`:

- even parity (NS):  GS energy `-½ Σ_k Λ_k`, lowest excited state is
  a *pair* of quasiparticles → minimum two-particle energy is
  `Λ_{k₁} + Λ_{k₂}` over the two smallest NS Λ's.
- odd parity (R):    GS energy `-½ Σ_k Λ_k`, lowest excited state is
  a *single* quasiparticle → minimum one-particle energy is
  `min_k Λ_k^R`.

The mass gap is `min(E_excited) − min(E_GS)` across both sectors.
"""
function _tfim_pbc_mass_gap(N::Int, J::Float64, h::Float64)
    Λ_NS = sort!([_tfim_Λ(k, J, h) for k in _tfim_pbc_momenta(N, :NS)])
    Λ_R = sort!([_tfim_Λ(k, J, h) for k in _tfim_pbc_momenta(N, :R)])

    E_GS_NS = -sum(Λ_NS) / 2
    E_GS_R = -sum(Λ_R) / 2

    # NS even-parity excitation: flip the two smallest modes
    E_ex_NS = E_GS_NS + (length(Λ_NS) ≥ 2 ? Λ_NS[1] + Λ_NS[2] : Λ_NS[1])
    # R odd-parity excitation: flip the smallest mode
    E_ex_R = E_GS_R + Λ_R[1]

    candidates = Float64[E_GS_NS, E_GS_R, E_ex_NS, E_ex_R]
    sort!(candidates)
    E_gs = candidates[1]
    for c in candidates[2:end]
        if c > E_gs + 1e-12 * max(1.0, abs(E_gs))
            return c - E_gs
        end
    end
    return 0.0
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch dispatch (PBC)
# ═══════════════════════════════════════════════════════════════════════════════

# Generate fetch methods over (struct quantity, internal Symbol) pairs.
const _TFIM_PBC_THERMAL_METHODS = (
    (FreeEnergy, :free_energy),
    (ThermalEntropy, :entropy),
    (SpecificHeat, :specific_heat),
    (MagnetizationX, :transverse_magnetization),
    (SusceptibilityXX, :transverse_susceptibility),
)

for (QTy, qsym) in _TFIM_PBC_THERMAL_METHODS
    @eval function fetch(model::TFIM, ::$QTy, bc::PBC; beta::Real, kwargs...)
        N = _bc_size(bc, kwargs)
        return _tfim_pbc_thermo($(QuoteNode(qsym)), N, model.J, model.h, beta)
    end
end

"""
    fetch(model::TFIM, ::Energy{:per_site}, bc::PBC; beta::Real, kwargs...) -> Float64

Per-site energy `ε(β) = -∂_β log Z / N` of the N-site PBC TFIM.  Native
granularity for PBC TFIM (the `:total` granularity is provided by the
generic conversion fallback in `src/core/quantities.jl`).
"""
function fetch(model::TFIM, ::Energy{:per_site}, bc::PBC; beta::Real, kwargs...)
    N = _bc_size(bc, kwargs)
    return _tfim_pbc_thermo(:energy_per_site, N, model.J, model.h, beta)
end

"""
    fetch(model::TFIM, ::MassGap, bc::PBC; kwargs...) -> Float64

Lowest excitation energy of the N-site PBC TFIM.  See
[`_tfim_pbc_mass_gap`](@ref) for sector handling.
"""
function fetch(model::TFIM, ::MassGap, bc::PBC; kwargs...)
    N = _bc_size(bc, kwargs)
    return _tfim_pbc_mass_gap(N, model.J, model.h)
end
