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
    references=["Yan-Huse-White 2011", "Depenbrock-McCulloch-Schollwöck 2012"],
    notes="e_0/J ≈ -0.4386(5) cylindrical-DMRG reference.",
)

@register(
    KagomeHeisenbergAFM,
    MassGap,
    Infinite,
    method=:dmrg_reference,
    reliability=:medium,
    tested_in="test/standalone/test_kagome_heisenberg_afm.jl",
    references=["Yan-Huse-White 2011", "Iqbal-Becca-Sorella-Poilblanc 2013"],
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
        "Kitaev-Preskill 2006",
        "Levin-Wen 2006",
        "Jiang-Wang-Balents 2012",
        "Yan-Huse-White 2011",
    ],
    notes="γ = log 2 in the Z₂-spin-liquid scenario; competing U(1) Dirac scenario predicts γ = 0.",
)
