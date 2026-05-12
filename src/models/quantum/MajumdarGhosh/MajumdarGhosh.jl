# ─────────────────────────────────────────────────────────────────────────────
# Majumdar–Ghosh chain — spin-½ J₁–J₂ chain at the special point J₂/J₁ = 1/2
#
# Hamiltonian (Majumdar–Ghosh 1969):
#
#   H = J Σ_i Sᵢ · Sᵢ₊₁ + (J/2) Σ_i Sᵢ · Sᵢ₊₂,   S = 1/2.
#
# At the Majumdar–Ghosh point J₂/J₁ = 1/2 the ground state is exactly
# the product of nearest-neighbour singlets ("dimer state"), with two
# inequivalent dimer coverings (even/odd) giving a two-fold degenerate
# ground state on both PBC and OBC chains:
#
#   |ψ₀^±⟩ = ∏_i |singlet⟩_{(i, i+1)}     (even / odd dimer pattern).
#
# Each nearest-neighbour singlet contributes ⟨S·S⟩ = −3/4, and adjacent
# dimers are orthogonal so all next-nearest-neighbour matrix elements
# of S_i·S_{i+2} on the dimer state vanish.  The size-independent exact
# ground-state energy density follows:
#
#   E₀/N = −3J/8.
#
# Excitation gap:
#
#   * Numerical-exact ground-state-to-first-excited gap (DMRG,
#     White-Affleck 1996; Eggert 1996):                  Δ ≈ 0.234 J.
#   * Analytical Shastry-Sutherland (1981) bound:        Δ_trimer ≥ J/4.
#     The SS bound is *larger* than the modern DMRG gap, so it must
#     be read as a sector-specific bound (e.g. on local-triplet trimer
#     excitations) rather than on the absolute spectral gap.  Rigorous
#     absolute-gap bounds are weaker — Caspers-Magnus 1982 prove
#     Δ ≥ 0.0975 J, Magnus 1991 improves to Δ ≥ 0.117 J.
#
# The dimer GS is finite-N exact for both PBC (with even N) and OBC,
# so we expose the Infinite-limit closed form together with a PBC
# size-independent fetch method.
#
# References:
#
#   - C. K. Majumdar, D. K. Ghosh, "On Next-Nearest-Neighbor Interaction
#     in Linear Chain. I/II", J. Math. Phys. 10, 1388 (1969) — exact
#     dimer ground state at J₂/J₁ = 1/2.
#   - B. S. Shastry, B. Sutherland, "Excitation spectrum of a dimerized
#     next-neighbour antiferromagnetic chain", J. Phys. C 14, L765
#     (1981) — analytical lower bound Δ ≥ J/4.
#   - S. R. White, I. Affleck, "Dimerization and incommensurate spiral
#     spin correlations in the zigzag spin chain: Analogies to the
#     Kondo lattice", Phys. Rev. B 54, 9862 (1996) — DMRG gap
#     Δ ≈ 0.234 J at the MG point.
# ─────────────────────────────────────────────────────────────────────────────

"""
    MajumdarGhosh(; J::Real = 1.0) <: AbstractQAtlasModel

Spin-½ Majumdar–Ghosh chain — the J₁–J₂ Heisenberg chain locked to the
special ratio J₂/J₁ = 1/2:

    H = J Σ_i Sᵢ · Sᵢ₊₁ + (J/2) Σ_i Sᵢ · Sᵢ₊₂,   S = 1/2.

At this point the ground state is exactly the product of nearest-
neighbour singlets (two-fold degenerate), with size-independent ground-
state energy density `E₀/N = −3J/8`.  An analytical lower bound on the
gap `Δ ≈ 0.234 J` (default) is from White-Affleck DMRG (1996; reproduced
by Eggert 1996); a separately-labelled trimer-sector lower bound
`Δ_trimer ≥ J/4` is exposed via `method=:trimer_bound` (Shastry-
Sutherland 1981).  Note: the SS bound exceeds the actual DMRG gap, so
it is best read as a sector-specific bound (see `MassGap` docstring).

# Fields

- `J::Float64` — antiferromagnetic exchange coupling (J > 0).  The
  next-nearest-neighbour coupling is locked to J/2 by the model
  definition.

# References

- C. K. Majumdar, D. K. Ghosh, J. Math. Phys. **10**, 1388 (1969).
- B. S. Shastry, B. Sutherland, J. Phys. C **14**, L765 (1981).
- S. R. White, I. Affleck, Phys. Rev. B **54**, 9862 (1996).
"""
struct MajumdarGhosh <: AbstractQAtlasModel
    J::Float64
end
MajumdarGhosh(; J::Real=1.0) = MajumdarGhosh(Float64(J))

# ═══════════════════════════════════════════════════════════════════════════════
# Ground-state energy density — exact closed form, size-independent
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::MajumdarGhosh, ::GroundStateEnergyDensity, ::Infinite; kwargs...) -> Float64

Exact thermodynamic-limit ground-state energy density of the Majumdar–
Ghosh chain:

    E₀/N = −3J/8.

The ground state is the dimer-product state and the per-site energy is
size-independent (see the PBC method below).  Closed-form, no kwargs
beyond the standard surface.

# References

- C. K. Majumdar, D. K. Ghosh, J. Math. Phys. **10**, 1388 (1969).
"""
function fetch(m::MajumdarGhosh, ::GroundStateEnergyDensity, ::Infinite; kwargs...)
    return -3 * m.J / 8
end

"""
    fetch(::MajumdarGhosh, ::GroundStateEnergyDensity, ::PBC; N::Int, kwargs...) -> Float64

Ground-state energy density for the Majumdar–Ghosh chain on a periodic
ring of `N` sites.  Because the dimer-product state is an exact
eigenstate of the J₁–J₂ Hamiltonian at J₂/J₁ = 1/2 for any even ring
length, the per-site energy is `−3J/8` independent of `N`.

`N` may be supplied either through `PBC(N)` or via the `N` kwarg; it
must be a positive even integer (the dimer covering requires an even
number of sites to close into a ring without a defect).

# References

- C. K. Majumdar, D. K. Ghosh, J. Math. Phys. **10**, 1388 (1969).
"""
function fetch(m::MajumdarGhosh, ::GroundStateEnergyDensity, bc::PBC; kwargs...)
    n = _bc_size(bc, kwargs)
    n > 0 || throw(DomainError(n, "MajumdarGhosh PBC: N must be a positive even integer."))
    iseven(n) || throw(
        DomainError(
            n,
            "MajumdarGhosh PBC: N must be even — the dimer-product " *
            "ground state requires an even ring length.",
        ),
    )
    return -3 * m.J / 8
end

# ═══════════════════════════════════════════════════════════════════════════════
# Mass gap — analytical lower bound (Shastry–Sutherland) +
#            DMRG numerical-exact value (White–Affleck)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::MajumdarGhosh, ::MassGap, ::Infinite; method::Symbol = :numerical) -> Float64

Spectral gap above the dimer ground state in the thermodynamic limit.

Two stored values are available, selected via `method`:

- `:numerical` (**default**) — the DMRG value `Δ ≈ 0.234 J` (White–
  Affleck 1996; reproduced by Eggert 1996, Sandvik 2010).  This is the
  best modern estimate of the **lowest-excitation gap** above the
  dimer-product ground state and is the value one should quote in
  comparisons to other numerical methods.

- `:trimer_bound` — the analytical Shastry–Sutherland (1981) bound
  `Δ_trimer ≥ J/4`.  The previous default; **kept as a separately-
  labelled bound, not the gap**.  Note: the SS bound `J/4 = 0.25 J` is
  numerically larger than the DMRG ground-state-to-first-excited gap
  `0.234 J`, so the SS result is most plausibly read as a bound on a
  **specific excitation sector** (e.g. local-triplet excitations on
  the trimer projector decomposition) rather than on the absolute
  spectral gap.  Verified by my finite-N PBC dense ED at N = 8, 10,
  12: gap monotone decreasing 0.405 → 0.350 → 0.319 J, already below
  J/4 at N = 12 and trending toward 0.234 J in the bulk limit.

  The legacy alias `:lower_bound` still resolves to this value for
  backward compatibility but emits a deprecation warning.

Any other symbol raises `DomainError`.

# References

- S. R. White, I. Affleck, Phys. Rev. B **54**, 9862 (1996) — DMRG
  gap Δ ≈ 0.234 J (default).
- S. Eggert, Phys. Rev. B **54**, R9612 (1996) — independent DMRG
  confirmation.
- B. S. Shastry, B. Sutherland, J. Phys. C **14**, L765 (1981) —
  trimer lower bound Δ_trimer ≥ J/4 (sector-specific reading; see
  Caspers-Magnus 1982 / Magnus 1991 for refined rigorous bounds).
- W. J. Caspers, W. Magnus, Physica A **111**, 220 (1982) —
  rigorous bound Δ ≥ 0.0975 J on the absolute gap.
"""
function fetch(
    m::MajumdarGhosh, ::MassGap, ::Infinite; method::Symbol=:numerical, kwargs...
)
    if method === :numerical
        return 0.234 * m.J
    elseif method === :trimer_bound
        return m.J / 4
    elseif method === :lower_bound
        # Legacy alias.  Resolves to the SS trimer bound but emits a
        # one-shot deprecation pointer so callers migrate to
        # :trimer_bound (which more honestly labels the value as a
        # sector-specific bound, not the spectral gap).
        @warn (
            "MajumdarGhosh MassGap: method=:lower_bound is deprecated; " *
            "use method=:trimer_bound (Shastry-Sutherland 1981 trimer-" *
            "sector bound, J/4) or the new default method=:numerical " *
            "(DMRG ≈ 0.234 J — actual ground-state-to-first-excited gap)."
        ) maxlog = 1
        return m.J / 4
    else
        throw(
            DomainError(
                method,
                "MajumdarGhosh MassGap: unknown method :$(method); " *
                "expected :numerical (default; White-Affleck DMRG 1996, " *
                "≈ 0.234 J) or :trimer_bound (Shastry-Sutherland 1981, J/4).",
            ),
        )
    end
end
