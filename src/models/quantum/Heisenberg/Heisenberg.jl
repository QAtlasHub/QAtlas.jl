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
struct ExactSpectrum end

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
struct GroundStateEnergyDensity end

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
for QTy in (:MagnetizationXLocal, :MagnetizationYLocal, :MagnetizationZLocal, :EnergyLocal)
    @eval function fetch(
        ::Heisenberg1D, ::$QTy, bc::OBC; beta::Real, J::Real=1.0, kwargs...
    )
        return fetch(XXZ1D(; J=J, Δ=1.0), $QTy(), bc; beta=beta, kwargs...)
    end
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

# Entanglement: VonNeumannEntropy and RenyiEntropy.
function fetch(
    ::Heisenberg1D,
    ::VonNeumannEntropy,
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
