# ─────────────────────────────────────────────────────────────────────────────
# KagomeHeisenbergAFM — spin-½ Heisenberg antiferromagnet on the Kagome lattice.
#
# Hamiltonian:
#
#     H = J Σ_{⟨i,j⟩} S_i · S_j,   S = 1/2,   J > 0,
#
# on the highly frustrated corner-sharing-triangle (Kagome) lattice.
# No closed-form ground state is known; the model is a paradigmatic
# quantum spin liquid candidate.  Best modern DMRG references
# (Yan-Huse-White 2011; Depenbrock-McCulloch-Schollwöck 2012):
#
#     e_0 / J  ≈ -0.4386(5)        (ground-state energy density)
#     Δ_s / J  ≈ 0.13              (spin gap; gapped Z₂ spin liquid)
#     γ        =  log 2            (topological entanglement entropy
#                                  ; Z₂ topological order)
#
# Variational Monte Carlo (Iqbal-Becca-Sorella-Poilblanc 2013)
# competes with a U(1) Dirac spin liquid at slightly higher energy,
# so the precise nature of the ground state remains a controversy
# in the open literature.  The reliability for the energy density
# and the spin gap is therefore reported as `:medium`.
#
# This Phase-1 entry registers the DMRG energy-density and spin-gap
# reference values.  The topological-entanglement-entropy `γ = log 2`
# (already exposed by [`TopologicalEntanglementEntropy`](@ref))
# delegation, and lattice-dependent finite-size scaling extrapolation
# infrastructure, are tracked as Phase 2.
#
# References:
#   - S. Yan, D. A. Huse, S. R. White, Science 332, 1173 (2011).
#   - S. Depenbrock, I. P. McCulloch, U. Schollwöck,
#     Phys. Rev. Lett. 109, 067201 (2012).
#   - Y. Iqbal, F. Becca, S. Sorella, D. Poilblanc,
#     Phys. Rev. B 87, 060405(R) (2013).
# ─────────────────────────────────────────────────────────────────────────────

"""
    KagomeHeisenbergAFM(; J::Real = 1.0) <: AbstractQAtlasModel

Spin-½ Heisenberg antiferromagnet on the Kagome lattice (highly
frustrated, quantum spin liquid candidate).

Quantities registered (Phase 1, DMRG reference values):

| Quantity                       | BC         | Method                              |
| ------------------------------ | ---------- | ----------------------------------- |
| [`Energy`](@ref) (`:per_site`) | `Infinite` | DMRG (Yan-Huse-White 2011)          |
| [`MassGap`](@ref)              | `Infinite` | DMRG (gapped Z₂ spin liquid)        |

Both registered with reliability `:medium` (different methods —
DMRG vs variational Monte Carlo — give slightly different gap
estimates; the precise spin-liquid character is an open question).

# References

- S. Yan, D. A. Huse, S. R. White, *Science* **332**, 1173 (2011).
- S. Depenbrock, I. P. McCulloch, U. Schollwöck,
  *Phys. Rev. Lett.* **109**, 067201 (2012).
- Y. Iqbal, F. Becca, S. Sorella, D. Poilblanc,
  *Phys. Rev. B* **87**, 060405(R) (2013).
"""
struct KagomeHeisenbergAFM <: AbstractQAtlasModel
    J::Float64
end
KagomeHeisenbergAFM(; J::Real=1.0) = KagomeHeisenbergAFM(Float64(J))

# Hardcoded DMRG reference values (Yan-Huse-White 2011 et seq.).
# Ground-state energy density e_0/J of the spin-1/2 Heisenberg AFM on
# the kagome lattice.  Reference: Yan-Huse-White, Science 332, 1173
# (2011), cylindrical DMRG extrapolated to the 2-D limit.  Confirmed by
# Depenbrock-McCulloch-Schollwock, PRL 109, 067201 (2012); the central
# value matches across both works at the -0.4386(5) band.  Update if a
# tighter literature consensus emerges (e.g. He-Zhu-Chen 2017 et seq.).
const _KAGOME_AFM_ENERGY_DENSITY_PER_J = -0.4386
# Singlet-triplet (spin) gap Delta_s/J in the gapped-Z2 spin-liquid
# scenario.  Reference: Yan-Huse-White 2011 DMRG (Science 332, 1173);
# treated as a DMRG upper bound because variational Monte Carlo
# (Iqbal-Becca-Sorella-Poilblanc, PRB 87, 060405R, 2013) favours a
# competing gapless U(1) Dirac spin liquid.  Later cylindrical DMRG
# (He-Zhu-Chen, PRX 7, 031020, 2017) reports Delta_s ~ 0.18(2); revisit
# if a tighter consensus emerges.
const _KAGOME_AFM_SPIN_GAP_PER_J = 0.13

# ═══════════════════════════════════════════════════════════════════════════════
# Ground-state energy per site (DMRG reference, Yan-Huse-White 2011)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::KagomeHeisenbergAFM, ::Energy{:per_site}, ::Infinite; J=m.J)
        -> Float64

Ground-state energy density of the spin-½ Kagome AFM:

    e_0 / J ≈ -0.4386(5)

(Yan-Huse-White 2011 DMRG; Depenbrock-McCulloch-Schollwöck 2012
confirm with cylindrical DMRG).  Returned as `J × (-0.4386)`.

# References

- S. Yan, D. A. Huse, S. R. White, *Science* **332**, 1173 (2011).
"""
function fetch(
    m::KagomeHeisenbergAFM, ::Energy{:per_site}, ::Infinite; J::Real=m.J, kwargs...
)
    J ≥ 0 || throw(
        DomainError(
            J,
            "KagomeHeisenbergAFM Energy(:per_site) requires J ≥ 0 (AF convention); got J = $J.",
        ),
    )
    return J * _KAGOME_AFM_ENERGY_DENSITY_PER_J
end

# ═══════════════════════════════════════════════════════════════════════════════
# Spin gap (DMRG reference, Yan-Huse-White 2011)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::KagomeHeisenbergAFM, ::MassGap, ::Infinite; J=m.J) -> Float64

Spin gap of the spin-½ Kagome AFM:

    Δ_s / J ≈ 0.13

(Yan-Huse-White 2011 DMRG; gapped Z₂ spin liquid scenario).
Returned as `J × 0.13`.  Variational Monte Carlo (Iqbal-Becca-
Sorella-Poilblanc 2013) favours a competing gapless U(1) Dirac
spin liquid, so the spin gap value should be treated as a DMRG
upper bound; reliability is therefore `:medium`.

# References

- S. Yan, D. A. Huse, S. R. White, *Science* **332**, 1173 (2011).
- Y. Iqbal, F. Becca, S. Sorella, D. Poilblanc,
  *Phys. Rev. B* **87**, 060405(R) (2013).
"""
function fetch(m::KagomeHeisenbergAFM, ::MassGap, ::Infinite; J::Real=m.J, kwargs...)
    J ≥ 0 || throw(
        DomainError(
            J, "KagomeHeisenbergAFM MassGap requires J ≥ 0 (AF convention); got J = $J."
        ),
    )
    return J * _KAGOME_AFM_SPIN_GAP_PER_J
end

# ═══════════════════════════════════════════════════════════════════════════════
# Topological entanglement entropy γ = log 2 (Z₂ scenario, Phase 2)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::KagomeHeisenbergAFM, ::TopologicalEntanglementEntropy, ::Infinite; kwargs...)
        -> Float64

Topological entanglement entropy `γ = log 2` for the Z₂ spin-liquid
ground-state scenario.  In the Kitaev-Preskill (2006) / Levin-Wen
(2006) prescription,

    S(ρ_A) = α |∂A| − γ + O(|∂A|⁻¹),

with `γ = log 𝒟` and total quantum dimension `𝒟 = √(Σ_a d_a²) = 2`
for Z₂ topological order (four Abelian anyons `{1, e, m, ψ}`, each
with `d_a = 1`).

Reliability is `:medium` — the value is the Z₂ topological prediction
(Yan-Huse-White 2011 DMRG; Jiang-Wang-Balents 2012 directly extracted
`γ ≈ log 2` from DMRG entanglement scans), but the competing U(1)
Dirac-spin-liquid scenario (Iqbal-Becca-Sorella-Poilblanc 2013)
gives a gapless variational ground state with no topological order;
the precise spin-liquid character is an open question.

# References

- A. Kitaev, J. Preskill, *Phys. Rev. Lett.* **96**, 110404 (2006).
- M. Levin, X.-G. Wen, *Phys. Rev. Lett.* **96**, 110405 (2006).
- H.-C. Jiang, Z. Wang, L. Balents, *Nature Phys.* **8**, 902 (2012).
- S. Yan, D. A. Huse, S. R. White, *Science* **332**, 1173 (2011).
"""
function fetch(
    ::KagomeHeisenbergAFM, ::TopologicalEntanglementEntropy, ::Infinite; kwargs...
)
    return log(2.0)
end
