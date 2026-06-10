# universalities/behaviour/conformal_2plus1d_registry.jl
#
# Register the 2+1D CFT universalities to make them queryable via the implementation status API.

const _CFT_2PLUS1D_PREDICTIONS = (
    (
        SphereFreeEnergy,
        ["KlebanovPufuSafdi2011"],
        "Universal sphere free energy F = -ln |Z(S^3)| (F-theorem).",
    ),
    (
        CornerEntanglementCoefficient,
        ["Chester2020"],
        "Universal corner coefficient / prefactor sigma in 2+1D CFT entanglement entropy.",
    )
)

for C in (:Ising, :XY, :Heisenberg), (q, refs, note) in _CFT_2PLUS1D_PREDICTIONS
    register!(
        Universality{C},
        q,
        Infinite;
        method=:analytic,
        reliability=:high,
        references=refs,
        notes=note,
    )
end
