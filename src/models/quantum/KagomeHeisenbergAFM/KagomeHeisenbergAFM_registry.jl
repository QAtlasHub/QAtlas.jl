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
