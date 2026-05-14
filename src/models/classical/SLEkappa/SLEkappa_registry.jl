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

@register(
    SLEkappa,
    FractalDimension,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/classical/test_sle_kappa.jl",
    references=["Beffara 2008"],
    notes="d_H(κ) = min(2, 1 + κ/8); cap at d_H = 2 for κ ≥ 8 (space-filling regime).",
)
