# models/quantum/KitaevHeisenberg/KitaevHeisenberg_registry.jl
#
# Declarative implementation map for the K-J-Γ honeycomb model
# (α-RuCl₃ family).  Schema documented in src/core/registry.jl.

@register(
    KitaevHeisenberg,
    MassGap,
    Infinite,
    method=:kitaev_delegation,
    reliability=:high,
    tested_in="test/standalone/test_kitaev_heisenberg.jl",
    references=["Kitaev 2006", "Rau-Lee-Kee 2014"],
    notes="K-only limit (J=Γ=0) delegated to KitaevHoneycomb; non-zero J or Γ raises DomainError (Phase 2: DMRG/ED).",
)
