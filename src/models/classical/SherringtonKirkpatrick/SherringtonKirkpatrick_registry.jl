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
    references=["Sherrington-Kirkpatrick 1975"],
    notes="T_c = J in 1/√N normalisation with J_ij ~ N(0, J²); 0 for J ≤ 0.",
)
