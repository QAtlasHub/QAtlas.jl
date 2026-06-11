# symmetry_registry.jl — model symmetry profiles (@symmetry, #700).
#
# Scope discipline (per #700): NOT a backfill of every model.  This is the
# starter set that exercises the C10 LSM check from every side — gapless
# half-integer chains (Heisenberg), LSM satisfied via ground-state degeneracy
# (MajumdarGhosh), integer-spin Haldane phases where the obstruction is absent
# (S1Heisenberg1D, AKLT1D), discrete internal symmetry where LSM does not
# apply (TFIM, Cluster1D), and parameter-dependent families that therefore
# declare no spectral fact (XXZ1D, Kitaev1D).  Profiles describe the model
# FAMILY at generic parameters; parameter-point enhancements (XXZ at Δ=1 is
# SU(2)) stay out until profiles grow `at` predicates.
#
# The SU(2) profiles also gate the family-isotropy identities declared in
# identity_registry.jl (χ_xx = χ_yy = χ_zz, m_x = m_y = m_z).

@symmetry(
    Heisenberg1D,
    internal = :SU2,
    translation = true,
    time_reversal = true,
    site_spin = 1//2,
    gapped = false,
    references = ["Bethe1931", "desCloizeauxPearson1962", "LiebSchultzMattis1961"],
    notes = "Gapless des Cloizeaux–Pearson spectrum; the LSM-consistent half-integer chain.",
)

@symmetry(
    XXZ1D,
    internal = :U1,
    translation = true,
    time_reversal = true,
    site_spin = 1//2,
    references = ["YangYang1966"],
    notes = "U(1) at generic Δ (SU(2) only at Δ=1); gapped/gapless is Δ-dependent, so no spectral fact is declared.",
)

@symmetry(
    S1Heisenberg1D,
    internal = :SU2,
    translation = true,
    time_reversal = true,
    site_spin = 1,
    gapped = true,
    gs_degeneracy = 1,
    references = ["WhiteHuse1993"],
    notes = "Haldane phase: integer spin evades the LSM obstruction — gapped with a unique bulk ground state.",
)

@symmetry(
    AKLT1D,
    internal = :SU2,
    translation = true,
    time_reversal = true,
    site_spin = 1,
    gapped = true,
    gs_degeneracy = 1,
    references = ["AKLT1988"],
    notes = "Exact VBS ground state; gapped, unique in the bulk (edge spins are a boundary effect).",
)

@symmetry(
    MajumdarGhosh,
    internal = :SU2,
    translation = true,
    time_reversal = true,
    site_spin = 1//2,
    gapped = true,
    gs_degeneracy = 2,
    references = ["MajumdarGhosh1969", "LiebSchultzMattis1961"],
    notes = "Gapped HALF-integer chain — LSM is satisfied through the two dimerized ground states (spontaneously broken translation).",
)

@symmetry(
    TFIM,
    internal = :Z2,
    translation = true,
    time_reversal = true,
    site_spin = 1//2,
    references = ["Pfeuty1970"],
    notes = "Discrete Z₂ only — no LSM obstruction; gapped except on the critical line h = J, so no spectral fact is declared.",
)

@symmetry(
    Kitaev1D,
    internal = :Z2,
    translation = true,
    time_reversal = true,
    references = ["Kitaev2001"],
    notes = "Fermion-parity Z₂; spinless fermions (no site_spin). Gap closes on |μ| = 2|t|, so no spectral fact is declared.",
)

@symmetry(
    Cluster1D,
    internal = :Z2xZ2,
    translation = true,
    time_reversal = true,
    site_spin = 1//2,
    gapped = true,
    gs_degeneracy = 1,
    references = ["BriegelRaussendorf2001"],
    notes = "Z₂×Z₂ SPT: gapped and unique is consistent — the protecting symmetry is discrete, LSM does not apply.",
)
