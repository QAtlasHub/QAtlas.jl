# models/classical/SherringtonKirkpatrick/SherringtonKirkpatrick_registry.jl
#
# Declarative implementation map for the SK mean-field spin glass.

@register(
    SherringtonKirkpatrick,
    CriticalTemperature,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_sherrington_kirkpatrick.jl",
    references=["SherringtonKirkpatrick1975"],
    notes="T_c = J in 1/√N normalisation with J_ij ~ N(0, J²); 0 for J ≤ 0.",
)

@register(
    SherringtonKirkpatrick,
    Energy{:per_site},
    Infinite,
    method=:variational_reference,
    status=:approx,
    valid_domain="SK model, thermodynamic limit, full replica-symmetry-breaking solution",
    error_order="high-precision numerical evaluation of the Parisi solution (Crisanti-Rizzo 2002); the closed form exists but is not evaluable in elementary functions",
    reliability=:high,
    tested_in="test/models/classical/test_sherrington_kirkpatrick.jl",
    references=["Parisi1980", "CrisantiRizzo2002", "Talagrand2006"],
    notes="e_0/J ≈ -0.7631667 (Crisanti-Rizzo 2002 high-precision Parisi full-RSB).",
)
