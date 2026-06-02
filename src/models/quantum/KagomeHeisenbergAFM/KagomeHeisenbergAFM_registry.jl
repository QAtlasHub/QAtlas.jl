# models/quantum/KagomeHeisenbergAFM/KagomeHeisenbergAFM_registry.jl
#
# DMRG reference values (Yan-Huse-White 2011) for the spin-½ Kagome AFM.

@register(
    KagomeHeisenbergAFM,
    Energy{:per_site},
    Infinite,
    method=:dmrg_reference,
    reliability=:medium,
    tested_in="test/standalone/test_kagome_heisenberg_afm.jl",
    references=["YanHuseWhite2011", "DepenbrockMcCullochSchollwock2012"],
    notes="e_0/J ≈ -0.4386(5) cylindrical-DMRG reference.",
)

@register(
    KagomeHeisenbergAFM,
    MassGap,
    Infinite,
    method=:dmrg_reference,
    reliability=:medium,
    tested_in="test/standalone/test_kagome_heisenberg_afm.jl",
    references=["YanHuseWhite2011", "IqbalBeccaSorellaPoilblanc2013"],
    notes="Δ_s/J ≈ 0.13 (Z₂ scenario, YHW); competing VMC Dirac spin liquid suggests upper bound.",
)

@register(
    KagomeHeisenbergAFM,
    TopologicalEntanglementEntropy,
    Infinite,
    method=:analytic,
    reliability=:medium,
    tested_in="test/models/quantum/Heisenberg/test_kagome_heisenberg_afm.jl",
    references=[
        "KitaevPreskill2006",
        "LevinWen2006",
        "JiangWangBalents2012",
        "YanHuseWhite2011",
        "IqbalBeccaSorellaPoilblanc2013",
    ],
    notes="γ = log 2 in the Z₂-spin-liquid scenario; competing U(1) Dirac scenario predicts γ = 0.",
)
