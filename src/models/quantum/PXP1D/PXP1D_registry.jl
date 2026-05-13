# models/quantum/PXP1D/PXP1D_registry.jl
#
# MPS thermodynamic-limit reference value (Surace et al. 2020) for the
# Rydberg-blockade / PXP chain.  Reliability :medium — numerical
# cluster-MPS rather than analytic closed form.  Issue #300.

@register(
    PXP1D,
    Energy{:per_site},
    Infinite,
    method=:mps_reference,
    reliability=:medium,
    tested_in="test/models/quantum/misc/test_pxp1d.jl",
    references=[
        "Surace 2020", "Lin-Motrunich 2019", "Turner-Michailidis-Abanin-Papić-Serbyn 2018"
    ],
    notes="e_0/Ω ≈ -0.6516(2) thermodynamic-limit MPS reference.",
)
