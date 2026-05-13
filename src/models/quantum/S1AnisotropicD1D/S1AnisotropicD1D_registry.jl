# models/quantum/S1AnisotropicD1D/S1AnisotropicD1D_registry.jl —
# declarative implementation map for S1AnisotropicD1D (S=1 Heisenberg
# chain with single-ion anisotropy D).
#
# Phase 1 registers only the D = 0 Haldane-chain reference point, which
# is delegated to the existing S1Heisenberg1D entry (White-Huse 1993
# DMRG, Δ ≈ 0.41048 J).  Non-zero D triggers DomainError at fetch time
# and is left to a DMRG-based Phase 2 (Chen-Roncaglia 2008;
# Tzeng-Yang-Hsu 2017).

# ── Infinite — Haldane gap delegate (D = 0 only) ────────────────────
@register(
    S1AnisotropicD1D,
    MassGap,
    Infinite,
    method=:s1_heisenberg_delegation,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_s1_anisotropic_d1d.jl",
    references=["White-Huse 1993", "Chen-Roncaglia 2008", "Tzeng-Yang-Hsu 2017"],
    notes="Phase 1: D = 0 delegate to S1Heisenberg1D (Δ ≈ 0.41048 J, DMRG numerical-exact). D ≠ 0 throws DomainError (Phase 2).",
)
