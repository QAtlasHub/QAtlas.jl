# models/classical/TTbar/TTbar_registry.jl
#
# Declarative implementation map for the universal irrelevant TT̄
# deformation (Zamolodchikov 2004).  Schema in src/core/registry.jl.

@register(
    TTbar,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/universalities/test_ttbar.jl",
    references=[
        "Zamolodchikov 2004",
        "Smirnov-Zamolodchikov 2017",
        "Cavaglià-Negro-Szécsényi-Tateo 2016",
    ],
    notes="Universal irrelevant TT̄ deformation; UV central charge c preserved at all λ. Deformed circle spectrum Phase 2.",
)
