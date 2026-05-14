# ─────────────────────────────────────────────────────────────────────────────
# PXP1D — Rydberg-blockade chain (PXP model).
#
# Hamiltonian (Fendley-Sengupta-Sachdev 2004; Bernien et al. 2017
# Rydberg-atom array realisation):
#
#     H = Ω Σ_i P_{i-1} σ^x_i P_{i+1},   P_i = (1 - σ^z_i)/2,   Ω > 0.
#
# Acts within the constrained Hilbert space of states with no two
# adjacent excited (Rydberg) sites — the Rydberg-blockade subspace.
# The model is famous for hosting **many-body quantum scars**
# (Turner-Michailidis-Abanin-Papić-Serbyn 2018, Nat. Phys. 14, 745):
# a sparse tower of non-thermal eigenstates that drive anomalously
# slow thermalisation and long-lived revivals from the Néel-like
# |Z₂⟩ product state, observed in 51-atom Rydberg arrays
# (Bernien et al. 2017, Nature 551, 579).
#
# Thermodynamic-limit ground-state energy density (numerical DMRG /
# ED convergence in the PXP literature):
#
#     e_0 / Ω  ≈  -0.6516(2)
#
# This value is the DMRG/ED reference quoted across the PXP-scar
# literature (Turner-Michailidis-Abanin-Serbyn-Papić 2018; Lin-
# Motrunich 2019; Iadecola-Schecter-Xu 2019); it is consistent with
# the U(1) lattice gauge theory mapping of Surace-Mazza-Giudici-
# Lerose-Gambassi-Dalmonte 2020 (PRX 10, 021041), which focuses on
# the Schwinger-model phase diagram rather than tabulating e_0 per se.
# Reliability `:medium` (numerical DMRG/ED reference;
# revisit if a tighter literature consensus emerges).
#
# Phase 1 exposes only `Energy{:per_site}` at `Infinite`.
# Scar-tower observables (revival frequency, OTOC scrambling,
# Z₂-state survival probability) are tracked as Phase 2.
#
# References:
#   - P. Fendley, K. Sengupta, S. Sachdev, Phys. Rev. B 69, 075106 (2004).
#   - H. Bernien et al., Nature 551, 579 (2017).
#   - C. J. Turner, A. A. Michailidis, D. A. Abanin, M. Serbyn, Z. Papić,
#     Nat. Phys. 14, 745 (2018).
#   - C. J. Lin, O. I. Motrunich, Phys. Rev. Lett. 122, 173401 (2019).
#   - T. Iadecola, M. Schecter, S. Xu, Phys. Rev. B 100, 184312 (2019).
#   - F. M. Surace et al., Phys. Rev. X 10, 021041 (2020).
#   - M. Serbyn, D. A. Abanin, Z. Papić, Rep. Prog. Phys. 84, 086601 (2021).
# ─────────────────────────────────────────────────────────────────────────────

"""
    PXP1D(; Ω::Real = 1.0) <: AbstractQAtlasModel

Rydberg-blockade chain (PXP model):

    H = Ω Σ_i P_{i-1} σ^x_i P_{i+1},   P_i = (1 - σ^z_i)/2,   Ω > 0.

Acts within the constrained Hilbert space forbidding adjacent
excitations.  Paradigmatic host of many-body quantum scars
(Turner et al. 2018).

Quantities registered (Phase 1):

| Quantity                       | BC         | Method                                |
| ------------------------------ | ---------- | ------------------------------------- |
| [`Energy`](@ref) (`:per_site`) | `Infinite` | DMRG/ED reference (PXP-scar lit.)     |

Reliability `:medium` (numerical DMRG/ED thermodynamic-limit
reference; not analytically closed).  Scar-tower observables —
revival frequency, OTOC, Z₂-state survival — are deferred to
Phase 2.

# References

- P. Fendley, K. Sengupta, S. Sachdev, *Phys. Rev. B* **69**, 075106 (2004).
- H. Bernien et al., *Nature* **551**, 579 (2017).
- C. J. Turner et al., *Nat. Phys.* **14**, 745 (2018).
- C. J. Lin, O. I. Motrunich, *Phys. Rev. Lett.* **122**, 173401 (2019).
- T. Iadecola, M. Schecter, S. Xu, *Phys. Rev. B* **100**, 184312 (2019).
- F. M. Surace et al., *Phys. Rev. X* **10**, 021041 (2020).
"""
struct PXP1D <: AbstractQAtlasModel
    Ω::Float64
    function PXP1D(Ω::Real)
        Ω > 0 || throw(DomainError(Ω, "PXP1D requires coupling Ω > 0; got Ω = $Ω."))
        return new(Float64(Ω))
    end
end
PXP1D(; Ω::Real=1.0) = PXP1D(Ω)

native_energy_granularity(::PXP1D, ::Infinite) = :per_site

# Hardcoded reference value for the PXP chain T=0 ground-state energy
# density in units of the coupling Ω.  Value sourced from DMRG / ED
# numerics in the PXP-scar literature; see:
#   - C. J. Turner, A. A. Michailidis, D. A. Abanin, M. Serbyn, Z. Papić,
#     Nat. Phys. 14, 745 (2018) — quantum many-body scars, ED reference.
#   - C. J. Lin, O. I. Motrunich, Phys. Rev. Lett. 122, 173401 (2019) —
#     scarring and slow thermalisation, ED/DMRG numerics.
#   - T. Iadecola, M. Schecter, S. Xu, Phys. Rev. B 100, 184312 (2019) —
#     bond-bimagnon states, DMRG reference for PXP ground state.
#   - F. M. Surace et al., Phys. Rev. X 10, 021041 (2020) — U(1) lattice
#     gauge theory mapping (Schwinger-model phase diagram; consistent
#     with the DMRG/ED value above but does not itself tabulate e_0).
#   - M. Serbyn, D. A. Abanin, Z. Papić, Rep. Prog. Phys. 84, 086601 (2021).
const _PXP1D_GROUND_STATE_ENERGY_DENSITY_PER_OMEGA = -0.6516

# ═══════════════════════════════════════════════════════════════════════════════
# Ground-state energy per site (DMRG/ED reference, PXP-scar literature)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::PXP1D, ::Energy{:per_site}, ::Infinite; Ω=m.Ω) -> Float64

Thermodynamic-limit ground-state energy density of the PXP chain:

    e_0 / Ω ≈ -0.6516(2)

Value sourced from DMRG / ED in the PXP-scar literature (see Turner 2018,
Lin-Motrunich 2019, Iadecola-Schecter-Xu 2019); consistent with the
U(1) gauge-theory mapping of Surace et al. 2020 (PRX 10, 021041), which
focuses on the Schwinger-model phase diagram rather than tabulating e_0.
Reliability `:medium`; update if a tighter literature consensus emerges.

# References

- C. J. Turner et al., *Nat. Phys.* **14**, 745 (2018) — quantum many-body scars.
- C. J. Lin, O. I. Motrunich, *Phys. Rev. Lett.* **122**, 173401 (2019).
- T. Iadecola, M. Schecter, S. Xu, *Phys. Rev. B* **100**, 184312 (2019).
- F. M. Surace et al., *Phys. Rev. X* **10**, 021041 (2020).
"""
function fetch(m::PXP1D, ::Energy{:per_site}, ::Infinite; Ω::Real=m.Ω, kwargs...)
    Ω > 0 || throw(DomainError(Ω, "PXP1D Energy requires Ω > 0; got Ω = $Ω."))
    return Ω * _PXP1D_GROUND_STATE_ENERGY_DENSITY_PER_OMEGA
end
