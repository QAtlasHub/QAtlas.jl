# duality_registry.jl — parameter-mapped duality edges (@dual, #699).
#
# Starter catalog: the two cleanest exact dualities in the vault, both with
# INDEPENDENT implementations on the two sides (the structural antidote to
# the circular-verification incident that motivated #699).
#
# Conventions (verified against the implementations' CONVENTION headers):
#   TFIM      H = -J Σ σᶻσᶻ - h Σ σˣ                  (Pauli σ)
#   Kitaev1D  H = -μ Σ c†c - t Σ (c†c + h.c.) + Δ Σ (cc + h.c.)
#
# Jordan–Wigner: H_TFIM(J, h) ≅ H_wire(μ = -2h, t = J, Δ = J) (σˣ = 1 - 2n;
# the pairing sign is a gauge).  The naive map carries an extensive -h·N
# offset, but the Kitaev1D implementation uses the particle-hole-symmetric
# convention (−μ Σ (n − ½)), which absorbs it — VERIFIED numerically: both
# sides agree to machine precision on every thermal density and the gap, so
# every value_map below is the identity.  (The per-quantity value_map slot
# exists exactly for maps where such an offset survives the conventions.)
#
# NOT declared here (yet): XXZ1D ↔ compact boson (bosonization needs the
# Luttinger-parameter map K(Δ) and operator renormalisations), and the
# universality-level `@coincides Universality{:Heisenberg} = Universality{:XY}
# at K = 1` (the SU(2)₁ self-dual-boson coincidence behind the #699 incident)
# — the latter needs kwargs-mapped universality fetches, tracked in #699.

# ── TFIM Kramers–Wannier self-duality: J ↔ h ──────────────────────────
# The quasiparticle spectrum ε(k) = 2√(J² + h² − 2Jh·cos k) is symmetric
# under J ↔ h, so every spectral/thermal scalar coincides at mapped
# parameters.  The order parameter (MagnetizationX, SpontaneousMagnetization)
# maps to the DISORDER operator and is deliberately absent from the
# quantity allowlist.
@dual(
    :tfim_kramers_wannier,
    TFIM,
    TFIM,
    param_map = (m -> TFIM(; J=m.h, h=m.J)),
    kind = :kramers_wannier,
    involution = true,
    examples = [TFIM(; J=1.0, h=0.5), TFIM(; J=0.3, h=1.2)],
    quantities = [
        (quantity=Energy{:per_site}, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=FreeEnergy, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=ThermalEntropy, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=SpecificHeat, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=MassGap, bc=Infinite),
    ],
    regime = "exact self-duality of the quantum chain (order ↔ disorder, J ↔ h)",
    references = ["Pfeuty1970"],
    notes = "Spectral symmetry of ε(k) under J ↔ h; order-parameter quantities are deliberately not listed.",
)

# ── Jordan–Wigner: TFIM ↔ Kitaev wire ─────────────────────────────────
# Two genuinely independent implementations (spin-side BdG vs wire-side BdG
# with different parameterizations and code paths) related by an exact
# fermionization; the implementations' conventions make every listed
# quantity map identically (see header note).
@dual(
    :tfim_kitaev_jordan_wigner,
    TFIM,
    Kitaev1D,
    param_map = (m -> Kitaev1D(; μ=-2 * m.h, t=m.J, Δ=m.J)),
    kind = :jordan_wigner,
    examples = [TFIM(; J=1.0, h=0.5), TFIM(; J=0.7, h=1.3)],
    quantities = [
        (quantity=Energy{:per_site}, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=FreeEnergy, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=ThermalEntropy, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=SpecificHeat, bc=Infinite, sweep=(beta=[0.5, 1.0, 2.0],)),
        (quantity=MassGap, bc=Infinite),
    ],
    regime = "exact fermionization on the whole (J, h) plane",
    references = ["Kitaev2001", "Pfeuty1970"],
    notes = "H_TFIM(J,h) ≅ H_wire(μ=-2h, t=J, Δ=J) in the wire's particle-hole-symmetric convention; densities, entropy and gap all map identically (verified to machine precision).",
)
