# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/KitaevHoneycomb/test_KitaevHoneycomb_thermal_trivial_batch.jl
#
# Trivial-temperature-limit verification cards for the Kitaev honeycomb
# matter-sector thermodynamics. Pure verify(); branches off main. Refs #381.
#
# Zero-mode safety check (verified locally for this sweep):
#   The OBC T → 0 entropy is `sum(σ) do s; _kitaev_logcosh2(βs) - βs·tanh(βs); end / N`,
#   so any exact-zero singular value of the bipartite hopping matrix M would
#   contribute log(2) per zero mode to the per-site sum and break the 1e-9
#   tolerance. svdvals(M) on Panza for all (Kx, Ky, Kz) × (Lx, Ly) in this
#   sweep gives min singular value ≥ 0.057 (worst case K=(0.5,0.5,0.5),
#   Lx=Ly=4), with zero counts at every threshold (1e-10, 1e-6, 1e-3)
#   equal to 0. No zero-mode contamination.
#
# Matter-sector scope note: the solver's high-T card (β² tail) tests the
# matter-sector solver's own asymptote, not the full spin model — flux
# fluctuations are out of scope for the matter-only solver.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "KitaevHoneycomb — thermal trivial limits (#381 batch)" begin
    LOW_T_BETA = 1e6
    HIGH_T_BETA = 1e-3

    for (Kx, Ky, Kz) in ((1.0, 1.0, 1.0), (1.0, 1.0, 2.0), (0.5, 0.5, 0.5))
        # ThermalEntropy T → 0: vanishing entropy (single GS, no residual entropy).
        verify(
            KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz),
            ThermalEntropy(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-9,
            refs=["KitaevHoneycomb T → 0: unique matter-sector GS ⇒ s = 0"],
            fetch_kw=(; beta=LOW_T_BETA),
        )
        # SpecificHeat T → 0: c → 0 (gap in gapped A-phases / power-law suppression in B phase).
        verify(
            KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz),
            SpecificHeat(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-9,
            refs=["KitaevHoneycomb T → 0: c → 0 (gap in A-phase / power-law in B-phase)"],
            fetch_kw=(; beta=LOW_T_BETA),
        )
        # SpecificHeat T → ∞: c → 0 as β² (bounded-spectrum high-T tail).
        verify(
            KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz),
            SpecificHeat(),
            Infinite();
            route=:limiting_case,
            independent=0.0,
            agree_within=1e-4,
            refs=[
                "KitaevHoneycomb T → ∞: c → 0 as β² high-T tail (bounded matter spectrum)"
            ],
            fetch_kw=(; beta=HIGH_T_BETA),
        )
        # OBC variants at T → 0
        for (Lx, Ly) in ((3, 3), (4, 4))
            verify(
                KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz),
                ThermalEntropy(),
                OBC(0);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["KitaevHoneycomb OBC T → 0: unique matter-sector GS ⇒ s = 0"],
                fetch_kw=(; Lx=Lx, Ly=Ly, beta=LOW_T_BETA),
            )
            verify(
                KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz),
                SpecificHeat(),
                OBC(0);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-9,
                refs=["KitaevHoneycomb OBC T → 0: c → 0"],
                fetch_kw=(; Lx=Lx, Ly=Ly, beta=LOW_T_BETA),
            )
            verify(
                KitaevHoneycomb(; Kx=Kx, Ky=Ky, Kz=Kz),
                SpecificHeat(),
                OBC(0);
                route=:limiting_case,
                independent=0.0,
                agree_within=1e-4,
                refs=["KitaevHoneycomb OBC T → ∞: c → 0 as β² high-T tail"],
                fetch_kw=(; Lx=Lx, Ly=Ly, beta=HIGH_T_BETA),
            )
        end
    end
end
