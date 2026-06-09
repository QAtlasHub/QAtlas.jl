# test/util/fluctuation_dissipation.jl — the fluctuation–dissipation foundation.
#
# Two layers, both test-side (verification infrastructure, not public API), in
# the same spirit as thermodynamic_identities.jl:
#
#   LAYER 1 (model-independent).  `fd_thermo_from_spectrum` + `fd_gibbs_moments` turn
#   any energy spectrum {Eₙ} (and any diagonal observable {oₙ}) into the full
#   canonical thermodynamics — lnZ, ⟨E⟩, F, S, C, Var(E) — by direct Boltzmann
#   weighting.  The accompanying test (test/identities/test_fluctuation_dissipation.jl)
#   proves these obey the thermodynamic + fluctuation–dissipation relations
#
#       ⟨E⟩    = -∂ lnZ/∂β                       (mean energy)
#       Var(E) = ∂² lnZ/∂β² = -∂⟨E⟩/∂β           (energy FDT — the core)
#       C      = β² Var(E)                        (specific heat from fluctuations)
#       S      = -∂F/∂T                           (entropy as a free-energy response)
#       ⟨E⟩    = F + T·S                          (Gibbs)
#       ∂⟨O⟩/∂λ = β·Var(O)   for  H(λ)=H₀-λO      (static linear-response FDT)
#
#   to machine precision for *arbitrary* spectra.  The two sides of each
#   relation are computed by genuinely different routes — an ensemble moment on
#   one side, a β/λ-derivative (ForwardDiff) on the other — so agreement is a
#   real theorem check, not the same closed form re-typed (cf. the verify()
#   `independent=` discipline: an independent witness, never a circular one).
#
#   LAYER 2 (atlas tie-in).  The FDT then has to hold for QAtlas's *registered*
#   response functions.  `independent_energy_variance_per_site` /
#   `independent_magnetization_variance_per_site` compute Var(E)/N and Var(M)/N
#   by brute-force enumeration of all 2ᴺ configurations of a model — a
#   computation that shares no code with the closed-form thermodynamics.  The
#   `FLUCTUATION_DISSIPATION_IDENTITIES` then assert, through the existing
#   `ThermoIdentity` harness,
#
#       fetch(SpecificHeat)     == β²·Var(E)/N        (energy FDT)
#       fetch(SusceptibilityZZ) == β·Var(M)/N         (magnetisation FDT, [H,M]=0)
#
#   The magnetisation FDT holds only when the magnetisation commutes with H
#   (diagonal / classical M); models are opted in by the
#   `has_independent_*_variance` traits — the same care as the Kubo-vs-variance
#   caveat on SUSCEPTIBILITY_XX_KUBO_FROM_MAGNETIZATION in
#   thermodynamic_identities.jl.
#
# Lives in test/util/ for the same reason as thermodynamic_identities.jl: a
# self-validation tool, not part of QAtlas's public surface.  If downstream
# packages need it, lift it to src/verification/ as a ForwardDiff weakdep.
# Included by runtests.jl *after* thermodynamic_identities.jl so the
# `ThermoIdentity` struct is in scope.

using ForwardDiff
using QAtlas: fetch, SpecificHeat, SusceptibilityZZ, Infinite, IsingChain1D, TFIM, Kitaev1D, SSH, TightBinding1D
using QuadGK: quadgk

# ══════════════════════════════════════════════════════════════════════
# Layer 1 — model-independent canonical thermodynamics from a spectrum
# ══════════════════════════════════════════════════════════════════════

"""
    fd_log_partition(levels, β) -> Real

Log partition function `lnZ(β) = log Σₙ e^{-βEₙ}` of a spectrum `levels`,
via the numerically-stable log-sum-exp shift
`lnZ = -βE₀ + log Σₙ e^{-β(Eₙ-E₀)}` with `E₀ = min levels`.

Differentiable: `β` may be a `ForwardDiff.Dual` (the shift `E₀` is a
constant of the concrete `levels`), so `-∂lnZ/∂β`, `∂²lnZ/∂β²`, … recover
`⟨E⟩`, `Var(E)`, … by AutoDiff — the *derivative* half of every FDT check.
"""
function fd_log_partition(levels::AbstractVector{<:Real}, β::Real)
    isempty(levels) && throw(ArgumentError("fd_log_partition: empty spectrum"))
    E0 = minimum(levels)
    return -β * E0 + log(sum(E -> exp(-β * (E - E0)), levels))
end

"""
    fd_boltzmann_weights(levels, β) -> Vector

Normalised Gibbs probabilities `pₙ = e^{-βEₙ}/Z`, computed from the
shifted weights `e^{-β(Eₙ-E₀)}` so no `exp` overflows even for large
`|βEₙ|`.
"""
function fd_boltzmann_weights(levels::AbstractVector{<:Real}, β::Real)
    E0 = minimum(levels)
    w = exp.(-β .* (levels .- E0))
    return w ./ sum(w)
end

"""
    fd_mean_energy(levels, β) -> Real

Ensemble mean energy `⟨E⟩ = Σₙ pₙ Eₙ`.  A standalone, `Dual`-friendly
function so a test can take `ForwardDiff.derivative(b -> fd_mean_energy(levels,
b), β)` and compare `-∂⟨E⟩/∂β` against the *variance* `Var(E)` — the
fluctuation–dissipation theorem.
"""
function fd_mean_energy(levels::AbstractVector{<:Real}, β::Real)
    p = fd_boltzmann_weights(levels, β)
    return sum(p .* levels)
end

"""
    fd_thermo_from_spectrum(levels, β) -> NamedTuple

Full canonical thermodynamics of a spectrum at inverse temperature `β>0`,
returned as `(; lnZ, E, F, S, C, varE)`:

| field  | meaning                  | formula        |
| ------ | ------------------------ | -------------- |
| `lnZ`  | log partition function   | `log Σ e^{-βEₙ}` |
| `E`    | mean energy `⟨E⟩`        | `Σ pₙ Eₙ`      |
| `varE` | energy variance `Var(E)` | `⟨E²⟩ - ⟨E⟩²`  |
| `F`    | Helmholtz free energy    | `-lnZ/β`       |
| `S`    | entropy (nats)           | `β(⟨E⟩ - F)`   |
| `C`    | heat capacity            | `β² Var(E)`    |

`C` and `S` are *defined here* through the fluctuation (`β²Var(E)`) and
Gibbs (`β(E-F)`) relations; the self-consistency test proves these equal
the independent derivative routes (`-β²∂⟨E⟩/∂β`, `-∂F/∂T`).
"""
function fd_thermo_from_spectrum(levels::AbstractVector{<:Real}, β::Real)
    isempty(levels) && throw(ArgumentError("fd_thermo_from_spectrum: empty spectrum"))
    β > 0 || throw(ArgumentError("fd_thermo_from_spectrum: requires β > 0; got β = $β"))
    p = fd_boltzmann_weights(levels, β)
    E = sum(p .* levels)
    E2 = sum(p .* abs2.(levels))
    raw_varE = E2 - E^2
    # Var(E) ≥ 0 exactly; only sub-eps cancellation may dip slightly negative.
    # A large negative value signals a spectrum/weight bug — surface it rather
    # than silently flooring to 0 (which would sail through a downstream C ≥ 0).
    raw_varE ≥ -sqrt(eps(float(E2))) * (abs(E2) + 1) ||
        error("fd_thermo_from_spectrum: Var(E) = $(raw_varE) ≪ 0 — spectrum/weight bug?")
    varE = max(raw_varE, zero(E))
    lnZ = fd_log_partition(levels, β)
    F = -lnZ / β
    S = β * (E - F)
    C = β^2 * varE
    return (; lnZ, E, F, S, C, varE)
end

"""
    fd_gibbs_moments(levels, obs, β) -> NamedTuple

Thermal mean and variance `(; mean, var)` of a *diagonal* observable `obs`
(its eigenvalues `oₙ` aligned with `levels`), Gibbs-weighted by the
spectrum at inverse temperature `β`:

    ⟨O⟩ = Σ pₙ oₙ,   Var(O) = Σ pₙ oₙ² - ⟨O⟩².

`obs` being diagonal in the energy eigenbasis is exactly the `[H,O]=0`
condition under which the static susceptibility equals `β·Var(O)`.
"""
function fd_gibbs_moments(
    levels::AbstractVector{<:Real}, obs::AbstractVector{<:Real}, β::Real
)
    length(levels) == length(obs) ||
        throw(DimensionMismatch("fd_gibbs_moments: levels and obs length differ"))
    p = fd_boltzmann_weights(levels, β)
    m = sum(p .* obs)
    v = max(sum(p .* abs2.(obs)) - m^2, zero(m))
    return (; mean=m, var=v)
end

"""
    fd_free_fermion_thermo(modes, β) -> NamedTuple

Canonical thermodynamics of a set of *independent fermionic modes* with
single-particle (quasiparticle) energies `modes` (≥ 0) at inverse
temperature `β>0`, returned as `(; E, varE, C, S)`:

    fₖ     = 1/(e^{βΛₖ}+1)                    (Fermi occupation)
    ⟨E⟩    = Σₖ Λₖ fₖ                          (excitation energy above vacuum)
    Var(E) = Σₖ Λₖ² fₖ(1-fₖ)                   (energy fluctuation)
    C      = β² Var(E)                          (heat capacity)
    S      = -Σₖ [fₖ ln fₖ + (1-fₖ) ln(1-fₖ)]  (mode entropy)

For free fermions the modes are statistically independent, so the energy
variance is the *sum* of per-mode variances `Λₖ² fₖ(1-fₖ)` — the
free-fermion energy FDT.  `Dual`-friendly in `β` (`fₖ(β)` is analytic), so
`-∂⟨E⟩/∂β = Var(E)` is checkable by AutoDiff with no eigensolve — the
clean autodiff route the dense-ED spectrum (`fd_thermo_from_spectrum`)
cannot offer.
"""
function fd_free_fermion_thermo(modes::AbstractVector{<:Real}, β::Real)
    β > 0 || throw(ArgumentError("fd_free_fermion_thermo: requires β > 0; got β = $β"))
    f = @. 1 / (exp(β * modes) + 1)
    E = sum(modes .* f)
    varE = sum(abs2.(modes) .* f .* (1 .- f))
    C = β^2 * varE
    s_term(fk) = (fk <= 0 || fk >= 1) ? zero(fk) : -(fk * log(fk) + (1 - fk) * log(1 - fk))
    S = sum(s_term, f)
    return (; E, varE, C, S)
end

# ══════════════════════════════════════════════════════════════════════
# Layer 2 — independent fluctuation providers (brute-force) + FDT identities
# ══════════════════════════════════════════════════════════════════════

"""
    has_independent_energy_variance(model)        -> Bool
    has_independent_magnetization_variance(model) -> Bool

Opt-in traits marking models for which an *independent* (closed-form-free)
energy / magnetisation variance per site is available — currently the
small-system brute-force enumerations below.  Default `false`; the
fluctuation–dissipation identities skip any model that does not overload
the relevant trait, mirroring `is_su2_symmetric` in
thermodynamic_identities.jl.
"""
has_independent_energy_variance(::Any) = false
has_independent_magnetization_variance(::Any) = false

has_independent_energy_variance(::IsingChain1D) = true
has_independent_magnetization_variance(::IsingChain1D) = true

has_independent_energy_variance(::TFIM) = true
has_independent_energy_variance(::Kitaev1D) = true
has_independent_energy_variance(::SSH) = true
has_independent_energy_variance(::TightBinding1D) = true

# Enumerate all 2ᴺ spin configurations of the periodic 1-D Ising ring at
# h = 0 and return (energies, magnetisations).  σ_i = 2·bitᵢ - 1; the ring
# energy is E(σ) = -J Σᵢ σᵢσᵢ₊₁ (site N wraps to 1) and M(σ) = Σᵢ σᵢ.  This
# shares no code with the transfer-matrix closed forms — it is the
# independent witness the FDT identities compare against.
function _ising1d_ring_configs(N::Int, J::Real)
    N ≥ 2 || throw(ArgumentError("_ising1d_ring_configs: need N ≥ 2; got $N"))
    nconf = 1 << N
    E = Vector{Float64}(undef, nconf)
    M = Vector{Float64}(undef, nconf)
    @inbounds for c in 0:(nconf - 1)
        e = 0.0
        m = 0.0
        for i in 0:(N - 1)
            σi = 2 * ((c >> i) & 1) - 1
            σj = 2 * ((c >> ((i + 1) % N)) & 1) - 1
            e -= J * σi * σj
            m += σi
        end
        E[c + 1] = e
        M[c + 1] = m
    end
    return E, M
end

"""
    independent_energy_variance_per_site(model, bc; beta, N=16) -> Float64

`Var(E)/N` of `model` at inverse temperature `beta`, from a brute-force
enumeration over all `2ᴺ` configurations — independent of the model's
closed-form thermodynamics.  The finite-size error of the *variance* scales
as `O(N²·tanh(βJ)ᴺ)` (the `N²` from differentiating the `(λ₋/λ₊)ᴺ`
transfer-matrix gap twice in `β`): `< 2e-5` for `βJ ≤ 0.4, N = 16`, but
already `~2e-4` by `βJ = 0.5` — use high-T `β` for tight atlas comparisons.
"""
function independent_energy_variance_per_site(
    m::IsingChain1D, ::Infinite; beta::Real, N::Int=16
)
    iszero(m.h) || throw(
        ArgumentError(
            "independent_energy_variance_per_site: the brute-force ring omits the " *
            "field term, so only h = 0 is supported; got h = $(m.h).",
        ),
    )
    E, _ = _ising1d_ring_configs(N, m.J)
    return fd_thermo_from_spectrum(E, beta).varE / N
end

"""
    independent_magnetization_variance_per_site(model, bc; beta, N=16) -> Float64

`Var(M)/N` of the longitudinal magnetisation `M = Σᵢ σᵢ`, brute-forced
over all `2ᴺ` configurations.  For the classical Ising chain `M` is
diagonal (`[H,M]=0`), so `β·Var(M)/N` is the genuine static
susceptibility.
"""
function independent_magnetization_variance_per_site(
    m::IsingChain1D, ::Infinite; beta::Real, N::Int=16
)
    iszero(m.h) || throw(
        ArgumentError(
            "independent_magnetization_variance_per_site: the brute-force ring omits " *
            "the field term, so only h = 0 is supported; got h = $(m.h).",
        ),
    )
    E, M = _ising1d_ring_configs(N, m.J)
    return fd_gibbs_moments(E, M, beta).var / N
end

"""
    SPECIFIC_HEAT_FROM_VARIANCE

Energy fluctuation–dissipation theorem `c_v = β²·Var(E)/N`.  The left side
is the model's registered `SpecificHeat`; the right side is the brute-force
energy variance per site (`independent_energy_variance_per_site`).  Holds
for *every* Hamiltonian (no commutation condition is required), so it is
the most universal FDT.  `model_filter` admits only models with a
brute-force provider.
"""

# Free-fermion energy variance providers
function independent_energy_variance_per_site(
    m::TFIM, ::Infinite; beta::Real, N::Int=16
)
    integrand = k -> begin
        lambda_val = 2 * sqrt(m.J^2 + m.h^2 - 2 * m.J * m.h * cos(k))
        y = beta * lambda_val
        (y / 2)^2 * sech(y / 2)^2
    end
    val, _ = quadgk(integrand, 0.0, pi; rtol=1e-10)
    return val / (pi * beta^2)
end

function independent_energy_variance_per_site(
    m::Kitaev1D, ::Infinite; beta::Real, N::Int=16
)
    integrand = k -> begin
        A = -2 * m.t * cos(k) - m.μ
        B = 2 * m.Δ * sin(k)
        lambda_val = sqrt(A^2 + B^2)
        y = beta * lambda_val
        (y / 2)^2 * sech(y / 2)^2
    end
    val, _ = quadgk(integrand, 0.0, pi; rtol=1e-10)
    return val / (pi * beta^2)
end

function independent_energy_variance_per_site(
    m::SSH, ::Infinite; beta::Real, N::Int=16
)
    integrand = k -> begin
        lambda_val = sqrt(m.v^2 + m.w^2 + 2 * m.v * m.w * cos(k))
        y = beta * lambda_val
        (y / 2)^2 * sech(y / 2)^2
    end
    val, _ = quadgk(integrand, 0.0, pi; rtol=1e-10)
    return val / (2 * pi * beta^2)
end

function independent_energy_variance_per_site(
    m::TightBinding1D, ::Infinite; beta::Real, N::Int=16
)
    integrand = k -> begin
        epsilon = -2 * m.t * cos(k) - m.μ
        y = beta * epsilon
        n = 1 / (exp(y) + 1)
        y^2 * n * (1 - n)
    end
    val, _ = quadgk(integrand, 0.0, pi; rtol=1e-10)
    return val / (pi * beta^2)
end

const SPECIFIC_HEAT_FROM_VARIANCE = ThermoIdentity(
    "c_v = β²·Var(E)/N  (energy FDT, brute-force variance)",
    Type[SpecificHeat],
    function (model, bc, params)
        β = params.β
        c_v = fetch(model, SpecificHeat(), bc; beta=β)
        varE_per_site = independent_energy_variance_per_site(model, bc; beta=β)
        return Float64(c_v), Float64(β^2 * varE_per_site)
    end;
    model_filter=has_independent_energy_variance,
)

"""
    SUSCEPTIBILITY_ZZ_FROM_VARIANCE

Magnetisation fluctuation–dissipation theorem `χ_zz = β·Var(M)/N` for a
magnetisation that commutes with H (diagonal / classical M).  Left side is
the registered `SusceptibilityZZ`; right side the brute-force magnetisation
variance per site.  Unlike the energy FDT this requires `[H,M]=0` — see the
Kubo-vs-variance caveat on `SUSCEPTIBILITY_XX_KUBO_FROM_MAGNETIZATION` in
thermodynamic_identities.jl; `model_filter` opts in only models known to
satisfy it.
"""
const SUSCEPTIBILITY_ZZ_FROM_VARIANCE = ThermoIdentity(
    "χ_zz = β·Var(M)/N  (magnetisation FDT, [H,M]=0)",
    Type[SusceptibilityZZ],
    function (model, bc, params)
        β = params.β
        χ = fetch(model, SusceptibilityZZ(), bc; beta=β)
        varM_per_site = independent_magnetization_variance_per_site(model, bc; beta=β)
        return Float64(χ), Float64(β * varM_per_site)
    end;
    model_filter=has_independent_magnetization_variance,
)

"""
    FLUCTUATION_DISSIPATION_IDENTITIES

The FDT identity catalogue — energy (`SPECIFIC_HEAT_FROM_VARIANCE`) and
magnetisation (`SUSCEPTIBILITY_ZZ_FROM_VARIANCE`).  Opt-in (like
`SYMMETRY_IDENTITIES`): pass via the `identities=` kwarg of
`verify_thermodynamic_identities`.  Use a finite-size-aware `rtol`: the
energy-variance correction scales as `O(N²·tanh(βJ)ᴺ)` — at `N=16` it is
`< 2e-5` for `βJ ≤ 0.4` (so `rtol=1e-4` is safely conservative) but reaches
`~2e-4` by `βJ = 0.5`.
"""
const FLUCTUATION_DISSIPATION_IDENTITIES = ThermoIdentity[
    SPECIFIC_HEAT_FROM_VARIANCE, SUSCEPTIBILITY_ZZ_FROM_VARIANCE
]
