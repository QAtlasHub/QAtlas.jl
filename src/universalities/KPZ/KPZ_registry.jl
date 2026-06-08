# universalities/KPZ/KPZ_registry.jl — the `predicts` edge of the KPZ class.
# Universality(:KPZ) --predicts--> GrowthExponents (status=:universal).

register!(
    Universality{:KPZ},
    GrowthExponents,
    Infinite;
    method=:analytic,  # status=:universal derived by construction (register!)
    reliability=:high,
    notes="KPZ scaling exponents (β_growth, α_rough, z) — Kardar-Parisi-Zhang 1986; exact in d=1, numerical d=2,3.",
)
