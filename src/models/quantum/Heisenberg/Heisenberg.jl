# ─────────────────────────────────────────────────────────────────────────────
# Heisenberg — spin-1/2 antiferromagnetic Heisenberg model
#
# Hamiltonian:
#   H = J Σ_{⟨i,j⟩} S_i · S_j
#
# where S_i are spin-1/2 operators and ⟨i,j⟩ runs over nearest-neighbor
# pairs. For J > 0 the ground state is a singlet.
#
# ─────────────────────────────────────────────────────────────────────────────
# Dimer (N = 2) — total-spin analysis
#
#   Two spin-1/2 degrees of freedom combine into S_tot = 0 ⊕ S_tot = 1.
#   The dimer Hamiltonian is
#
#     H = J S_1 · S_2 = (J/2)·[(S_1 + S_2)² − S_1² − S_2²]
#       = (J/2)·[S_tot² − 3/2]
#
#   giving eigenvalues
#
#     S_tot = 0 (singlet):   E_s = −3J/4
#     S_tot = 1 (triplet):   E_t = +J/4  (three-fold degenerate)
#
#   Singlet–triplet gap: Δ = E_t − E_s = J.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tags
# ═══════════════════════════════════════════════════════════════════════════════

# CONVENTION
#   Hamiltonian: Spin S (this file)
#   Observable:  Spin S         (QAtlas-wide spin convention; see docs/src/conventions.md)

"""
    Heisenberg1D

Dispatch tag for the spin-1/2 antiferromagnetic Heisenberg model on a 1D
chain (or more generally any finite spin-1/2 cluster). Hamiltonian:

    H = J Σ_{⟨i,j⟩} S_i · S_j,   spin-1/2, J > 0 antiferromagnetic
"""
struct Heisenberg1D <: AbstractQAtlasModel end

"""
    ExactSpectrum

Dispatch tag for the full sorted eigenvalue spectrum of a finite model.
"""
struct ExactSpectrum <: AbstractQuantity end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: exact spectrum for small N
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Heisenberg1D, ::ExactSpectrum; N, J=1.0) -> Vector{Float64}

Return the sorted exact spectrum of the spin-1/2 Heisenberg Hamiltonian
on an N-site chain or ring with boundary condition `bc`.

# Supported cases

- **N=2, bc=:OBC** (dimer): `[-3J/4, J/4, J/4, J/4]`
  (singlet E_s = -3J/4, triplet E_t = J/4, three-fold degenerate).
- **N=4, bc=:PBC** (4-site ring): `[-2J, -J×3, 0×7, +J×5]`
  Ground state E₀ = -2J (unique singlet). The ferromagnetic quintet
  sits at E = +J. The full degeneracy structure is:
  1 singlet + 1 triplet + (1 singlet + 2 triplets at E=0) + 1 quintet.

# Arguments
- `N::Int`: number of spin-1/2 sites
- `J::Real`: Heisenberg coupling constant (default 1.0; J > 0 AFM)
- `bc::Symbol`: boundary condition, `:OBC` (default) or `:PBC`

# References
    A. Auerbach, "Interacting Electrons and Quantum Magnetism" (1994), §2.
    H. Bethe, Z. Physik 71, 205 (1931).
"""
function fetch(::Heisenberg1D, ::ExactSpectrum; N::Int, J::Real=1.0, bc::Symbol=:OBC)
    if N == 2 && bc == :OBC
        return sort([-3J / 4, J / 4, J / 4, J / 4])
    elseif N == 4 && bc == :PBC
        # Exact spectrum of the 4-site PBC Heisenberg ring H = J Σ S_i·S_{i+1}:
        #   E = -2J ×1 (singlet), -J ×3 (triplet), 0 ×7, +J ×5 (quintet)
        return sort([
            -2J,
            -J,
            -J,
            -J,
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            J,
            J,
            J,
            J,
            J,
        ])
    end
    return error(
        "Heisenberg1D exact spectrum: only (N=2, OBC) and (N=4, PBC) " *
        "are implemented; got (N=$N, bc=$bc).",
    )
end

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tag + fetch: Bethe ansatz ground-state energy density
# ═══════════════════════════════════════════════════════════════════════════════

"""
    GroundStateEnergyDensity

Dispatch tag for the ground-state energy per site in the thermodynamic
limit (N → ∞).
"""
struct GroundStateEnergyDensity <: AbstractQuantity end

"""
    fetch(::Heisenberg1D, ::GroundStateEnergyDensity; J=1.0) -> Float64

Exact ground-state energy per site of the spin-1/2 antiferromagnetic
Heisenberg chain in the thermodynamic limit (N → ∞, PBC):

    e₀ = J (1/4 − ln 2) ≈ −0.4431 J

This is one of the earliest and most celebrated results of the Bethe
ansatz. The derivation proceeds by solving the Bethe equations for the
ground state of

    H = J Σᵢ Sᵢ · Sᵢ₊₁

in the limit N → ∞, yielding a linear integral equation for the
rapidity distribution whose solution gives the energy via integration.

# Finite-size corrections

For a PBC chain of N sites, the ground-state energy density approaches
e₀ with corrections of order 1/N² (logarithmic corrections also
present):

    E₀(N)/N = e₀ + O(1/N²)

See `test/verification/test_universality_cross_check.jl` for a
finite-size extrapolation verification using ED at N = 4, 6, 8.

# Arguments
- `J::Real`: Heisenberg coupling constant (default 1.0; J > 0 AFM)

# References
    H. Bethe, "Zur Theorie der Metalle. I. Eigenwerte und Eigenfunktionen
      der linearen Atomkette", Z. Physik 71, 205–226 (1931) — original
      Bethe ansatz solution.
    L. Hulthén, "Über das Austauschproblem eines Kristalles",
      Ark. Mat. Astron. Fys. 26A, No. 11, 1–106 (1938) — first
      evaluation of e₀ = 1/4 − ln 2 from the Bethe equations.
"""
function fetch(::Heisenberg1D, ::GroundStateEnergyDensity; J::Real=1.0)
    return J * (1 // 4 - log(2))
end

"""
    fetch(::Heisenberg1D, ::GroundStateEnergyDensity, ::Infinite; J=1.0) -> Float64

BC-explicit dispatch sister of the legacy `fetch(::Heisenberg1D,
::GroundStateEnergyDensity; J=1.0)` method.  The thermodynamic-limit
ground-state energy density is only meaningful at `Infinite`, so the
two methods return the same Hulthén value
`e₀ = J(1/4 - ln 2)`.  Provided so `which(fetch, ::Heisenberg1D,
::GroundStateEnergyDensity, ::Infinite)` resolves and the registry
drift guard passes.
"""
function fetch(::Heisenberg1D, ::GroundStateEnergyDensity, ::Infinite; J::Real=1.0)
    return J * (1 // 4 - log(2))
end

native_energy_granularity(::Heisenberg1D, ::OBC) = :total

"""
    fetch(::Heisenberg1D, ::Energy{:total}, ::OBC; beta, J=1.0) -> Float64

**Total** thermal energy `⟨H⟩_β` for the spin-½ antiferromagnetic
Heisenberg OBC chain at finite `N` (the isotropic point `Δ = 1` of
[`XXZ1D`](@ref)).  Routes through
[`fetch(::XXZ1D, ::Energy{:total}, ::OBC)`](@ref).

Since `Heisenberg1D` currently carries no `J` field, callers must pass
`J` as a kwarg (default `J = 1.0`).  Downstream bridges (e.g.
ITensorModels `to_qatlas(::Heisenberg1D)`) lose `J` on conversion; use
`XXZ1D(; J, Δ=1)` directly if you need a non-unit coupling.
"""
function fetch(
    ::Heisenberg1D, ::Energy{:total}, bc::OBC; beta::Real, J::Real=1.0, kwargs...
)
    return fetch(XXZ1D(; J=J, Δ=1.0), Energy{:total}(), bc; beta=beta)
end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch — finite-N OBC delegators to XXZ1D(Δ = 1)
#
# `Heisenberg1D` carries no model parameters of its own (J is passed as a
# kwarg), so every quantity on the OBC chain at finite temperature is a
# thin wrapper around the corresponding `XXZ1D(; J, Δ = 1.0)` fetch.  We
# keep this enumerated rather than dispatching through a generic
# fallback so that:
#
# - the public surface is explicit and grep-friendly,
# - each method's docstring can name the delegate model,
# - the registry (`Heisenberg_registry.jl`) lists each triple.
# ═══════════════════════════════════════════════════════════════════════════════

# Scalar thermodynamic potentials — per-site at OBC, matching XXZ1D.
for QTy in (:FreeEnergy, :ThermalEntropy, :SpecificHeat)
    @eval function fetch(
        ::Heisenberg1D, ::$QTy, bc::OBC; beta::Real, J::Real=1.0, kwargs...
    )
        return fetch(XXZ1D(; J=J, Δ=1.0), $QTy(), bc; beta=beta, kwargs...)
    end
end

# Bulk-averaged magnetisations and susceptibilities (scalar).
for QTy in (
    :MagnetizationX,
    :MagnetizationY,
    :MagnetizationZ,
    :SusceptibilityXX,
    :SusceptibilityYY,
    :SusceptibilityZZ,
)
    @eval function fetch(
        ::Heisenberg1D, ::$QTy, bc::OBC; beta::Real, J::Real=1.0, kwargs...
    )
        return fetch(XXZ1D(; J=J, Δ=1.0), $QTy(), bc; beta=beta, kwargs...)
    end
end

# Site-resolved local observables (Vector{Float64}).
# Note: `MagnetizationXLocal` is parametric on the mode `M` (equilibrium / quench);
# Heisenberg1D currently only implements the equilibrium branch so we dispatch
# explicitly on `MagnetizationXLocal{:equilibrium}` to avoid accidentally
# capturing quench requests that this model does not support.
for QTy in (:MagnetizationYLocal, :MagnetizationZLocal, :EnergyLocal)
    @eval function fetch(
        ::Heisenberg1D, ::$QTy, bc::OBC; beta::Real, J::Real=1.0, kwargs...
    )
        return fetch(XXZ1D(; J=J, Δ=1.0), $QTy(), bc; beta=beta, kwargs...)
    end
end
function fetch(
    ::Heisenberg1D,
    ::MagnetizationXLocal{:equilibrium},
    bc::OBC;
    beta::Real,
    J::Real=1.0,
    kwargs...,
)
    return fetch(XXZ1D(; J=J, Δ=1.0), MagnetizationXLocal(), bc; beta=beta, kwargs...)
end

# Two-point correlators (static + connected).
for CorrTy in (:XXCorrelation, :YYCorrelation, :ZZCorrelation)
    for mode in (:static, :connected)
        @eval function fetch(
            ::Heisenberg1D,
            ::$CorrTy{$(QuoteNode(mode))},
            bc::OBC;
            beta::Real,
            i::Int,
            j::Int,
            J::Real=1.0,
            kwargs...,
        )
            return fetch(
                XXZ1D(; J=J, Δ=1.0),
                $CorrTy{$(QuoteNode(mode))}(),
                bc;
                beta=beta,
                i=i,
                j=j,
                kwargs...,
            )
        end
    end
end

import SpecialFunctions

# Helper for the alternating zeta function: zeta_a(s) = (1 - 2^(1-s)) * zeta(s)
_zeta_a(s::Int) = (1 - 2.0^(1-s)) * SpecialFunctions.zeta(s)

# ═══════════════════════════════════════════════════════════════════════════════
# Exact Ground-State Correlation Functions (Infinite, T=0)
# ═══════════════════════════════════════════════════════════════════════════════
for CorrTy in (:XXCorrelation, :YYCorrelation, :ZZCorrelation)
    for mode in (:static, :connected)
        @eval function fetch(
            ::Heisenberg1D,
            ::$CorrTy{$(QuoteNode(mode))},
            ::Infinite;
            beta::Real=Inf,
            i::Int,
            j::Int,
            J::Real=1.0,
            kwargs...,
        )
            isinf(beta) || throw(ArgumentError("Exact infinite-size correlation functions for Heisenberg1D are only available at T=0 (beta=Inf)."))
            
            r = abs(i - j)
            
            if r == 0
                return 0.25  # <S^z_i S^z_i> = 1/4
            elseif r == 1
                return 1/12 - 1/3 * log(2)
            elseif r == 2
                return 1/12 - 4/3 * log(2) + _zeta_a(3)
            elseif r == 3
                z1 = log(2)
                z3 = _zeta_a(3)
                z5 = _zeta_a(5)
                return 1/12 - 3*z1 + 74/9*z3 - 56/9*z1*z3 - 8/3*z3^2 - 50/9*z5 + 80/9*z1*z5
            elseif r == 4
                z1 = log(2)
                z3 = _zeta_a(3)
                z5 = _zeta_a(5)
                z7 = _zeta_a(7)
                return 1/12 - 16/3*z1 + 290/9*z3 - 72*z1*z3 - 1172/9*z3^2 - 700/9*z5 - 400/3*z5^2 + 4640/9*z1*z5 - 220/9*z3*z5 + 455/9*z7 - 3920/9*z1*z7 + 280*z3*z7
            elseif r == 5
                z1 = log(2)
                z3 = _zeta_a(3)
                z5 = _zeta_a(5)
                z7 = _zeta_a(7)
                z9 = _zeta_a(9)
                return 1/12 - 25/3*z1 + 800/9*z3 - 1192/3*z1*z3 - 15368/9*z3^2 - 608*z3^3 - 4228/9*z5 + 64256/9*z1*z5 - 976/9*z3*z5 + 3648*z1*z3*z5 - 3328/3*z3^2*z5 - 76640/3*z5^2 + 66560/3*z1*z5^2 + 12640/3*z3*z5^2 + 6400/3*z5^3 + 9674/9*z7 - 225848/9*z1*z7 + 56952*z3*z7 - 116480/3*z1*z3*z7 - 35392/3*z3^2*z7 + 7840*z5*z7 - 8960*z3*z5*z7 - 66640/3*z7^2 + 31360*z1*z7^2 - 686*z9 + 18368*z1*z9 - 53312*z3*z9 + 35392*z1*z3*z9 + 16128*z3^2*z9 + 38080*z5*z9 - 53760*z1*z5*z9
            else
                error("NotImplemented: Exact analytical correlation functions for Heisenberg1D are only implemented up to distance r=5. Got r=$r.")
            end
        end
    end
end
# Entanglement at Infinite delegates to Universality(:Heisenberg) which
# carries the c=1 Calabrese–Cardy closed forms (issue #580 Phase 1).
# The Heisenberg chain is gapless for all J (always at the SU(2) critical
# point), so no gapped crossover is needed.
function fetch(
    ::Heisenberg1D,
    ::VonNeumannEntropy{:equilibrium},
    ::Infinite;
    ℓ::Int,
    beta::Real=Inf,
    kwargs...,
)
    return fetch(
        Universality(:Heisenberg),
        VonNeumannEntropy(),
        Infinite();
        ℓ=ℓ,
        beta=beta,
        kwargs...,
    )
end

function fetch(
    ::Heisenberg1D, q::RenyiEntropy, ::Infinite; ℓ::Int, beta::Real=Inf, kwargs...
)
    return fetch(Universality(:Heisenberg), q, Infinite(); ℓ=ℓ, beta=beta, kwargs...)
end

# Entanglement: VonNeumannEntropy and RenyiEntropy.
function fetch(
    ::Heisenberg1D,
    ::VonNeumannEntropy{:equilibrium},
    bc::OBC;
    ℓ::Int,
    beta::Real=Inf,
    J::Real=1.0,
    kwargs...,
)
    return fetch(XXZ1D(; J=J, Δ=1.0), VonNeumannEntropy(), bc; ℓ=ℓ, beta=beta, kwargs...)
end

function fetch(
    ::Heisenberg1D, q::RenyiEntropy, bc::OBC; ℓ::Int, beta::Real=Inf, J::Real=1.0, kwargs...
)
    return fetch(XXZ1D(; J=J, Δ=1.0), q, bc; ℓ=ℓ, beta=beta, kwargs...)
end

# Mass gap.
function fetch(::Heisenberg1D, ::MassGap, bc::OBC; J::Real=1.0, kwargs...)
    return fetch(XXZ1D(; J=J, Δ=1.0), MassGap(), bc; kwargs...)
end

function fetch(::Heisenberg1D, ::MassGap, bc::Infinite; J::Real=1.0, kwargs...)
    return fetch(XXZ1D(; J=J, Δ=1.0), MassGap(), bc; kwargs...)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Luttinger parameter at the isotropic SU(2)-symmetric point (Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Heisenberg1D, ::LuttingerParameter, ::Infinite; J=1.0) -> Float64

Luttinger-liquid parameter `K = 1/2` of the spin-½ Heisenberg
antiferromagnetic chain at the SU(2)-symmetric point — the exact
Luther–Peschel 1975 / Affleck 1989 result.

This is the `Δ → 1` limit of the XXZ Luttinger parameter
`K_XXZ(Δ) = π / (2 arccos(−Δ))` (Haldane 1980); delegated to
`XXZ1D(Δ=1.0)`.

# References

- A. Luther, I. Peschel, *Phys. Rev. B* **12**, 3908 (1975).
- I. Affleck, *J. Phys. A* **22**, 1003 (1989).
- F. D. M. Haldane, *Phys. Rev. Lett.* **45**, 1358 (1980).
"""
function fetch(::Heisenberg1D, ::LuttingerParameter, ::Infinite; J::Real=1.0, kwargs...)
    J > 0 ||
        throw(DomainError(J, "Heisenberg1D LuttingerParameter requires J > 0; got J = $J."))
    # Delegate to XXZ1D at Δ=1 — K is J-independent (only depends on Δ)
    return fetch(XXZ1D(; Δ=1.0), LuttingerParameter(), Infinite())
end
