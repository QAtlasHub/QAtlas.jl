# models/quantum/PXP1D/PXP1D_registry.jl
#
# DMRG / ED thermodynamic-limit reference value for the Rydberg-blockade
# / PXP chain ground-state energy density.  Reliability :medium —
# numerical (DMRG/ED), not analytic closed form.  Issue #300.

@register(
    PXP1D,
    Energy{:per_site},
    Infinite,
    method=:dmrg_reference,
    reliability=:medium,
    tested_in="test/models/quantum/misc/test_pxp1d.jl",
    references=[
        "TurnerMichailidisAbaninSerbynPapic2018",
        "LinMotrunich2019",
        "IadecolaSchecterXu2019",
        "Surace2020",
    ],
    notes="e_0/Ω ≈ -0.6516(2) DMRG/ED reference (PXP-scar literature); Surace 2020 gauge-theory mapping consistent but does not tabulate e_0.",
)
