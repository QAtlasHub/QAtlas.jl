# models/classical/SLEkappa/SLEkappa_registry.jl
#
# Declarative implementation map for SLE_κ (Schramm 2000).
# Schema documented in src/core/registry.jl.

@register(
    SLEkappa,
    CentralCharge,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_sle_kappa.jl",
    references=["Schramm 2000", "Bauer-Bernard 2006", "Cardy 2005"],
    notes="SLE-CFT correspondence c(κ) = (3κ-8)(6-κ)/(2κ); symmetric under κ↔16/κ.",
)
