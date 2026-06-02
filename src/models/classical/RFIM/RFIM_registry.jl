# models/classical/RFIM/RFIM_registry.jl
#
# Random-Field Ising — Imry-Ma 1975 lower critical dimension.

@register(
    RFIM,
    CriticalTemperature,
    Infinite,
    method=:analytic_imry_ma,
    reliability=:high,
    tested_in="test/standalone/test_rfim.jl",
    references=["ImryMa1975", "Imbrie1985", "BricmontKupiainen1987"],
    notes="T_c = 0 for d ≤ 2 at Δ > 0 (Imry-Ma); d ≥ 3 raises DomainError (no closed form, Phase 2 numerical reference).",
)
