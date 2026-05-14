# models/quantum/ChernSimons3D/ChernSimons3D_registry.jl
#
# Declarative implementation map for SU(N)_k Chern-Simons (Witten 1989).
# Schema in src/core/registry.jl.

@register(
    ChernSimons3D,
    CentralCharge,
    Infinite,
    method=:sugawara,
    reliability=:high,
    tested_in="test/standalone/test_chern_simons_3d.jl",
    references=["Witten 1989", "Knizhnik-Zamolodchikov 1984"],
    notes="Sugawara c = k(N²-1)/(k+N) for boundary ŝu(N)_k WZW; SU(2)_k specialises to 3k/(k+2).",
)

@register(
    ChernSimons3D,
    PartitionFunction,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_chern_simons_3d.jl",
    references=["Witten 1989", "Verlinde 1988"],
    notes="Z(S³; SU(N)_k) = S_{0,0} of modular S-matrix; product over positive roots in (k+N).",
)
