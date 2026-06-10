# universalities/behaviour/conformal_2plus1d_registry.jl
#
# Register the 2+1D CFT universalities to make them queryable via the implementation status API.

const _CFT_2PLUS1D_PREDICTIONS = (
    (
        SphereFreeEnergy,
        ["KlebanovPufuSafdi2011", "Pufu2017"],
        "Universal sphere free energy F = -ln Z(S^3) (F-theorem). Values from Pufu (2017) Table 1.",
    ),
    (
        CornerEntanglementCoefficient,
        ["BuenoMyersWitczakKrempa2015", "Kos2016", "Chester2020"],
        "Universal corner coefficient sigma = pi^2/24 * C_T in 2+1D CFT entanglement entropy (smooth-limit prefactor).",
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
