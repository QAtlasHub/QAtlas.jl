# ─────────────────────────────────────────────────────────────────────────────
# test/models/quantum/KitaevHoneycomb/test_KitaevHoneycomb_thermal_trivial_batch.jl
#
# Trivial-temperature-limit verification cards for the Kitaev honeycomb
# matter-sector thermodynamics. Pure verify(); branches off main. Refs #381.
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
            refs=["KitaevHoneycomb T → ∞: c → 0 as β² high-T tail (bounded matter spectrum)"],
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
