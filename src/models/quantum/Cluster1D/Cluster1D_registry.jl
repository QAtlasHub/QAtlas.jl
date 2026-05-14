# models/quantum/Cluster1D/Cluster1D_registry.jl
#
# Declarative implementation map for the 1D Z₂×Z₂ SPT cluster Hamiltonian.
# Schema in src/core/registry.jl.

@register(
    Cluster1D,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_cluster1d.jl",
    references=["Briegel-Raussendorf 2001"],
    notes="E_0/N = -J; ground state is the cluster state (stabiliser model).",
)

@register(
    Cluster1D,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_cluster1d.jl",
    references=["Briegel-Raussendorf 2001"],
    notes="Δ = 2J; single-stabiliser-flip excitation.",
)
