# universalities/MeanField/MeanField_registry.jl — the `predicts` edge of the
# mean-field class.  Universality(:MeanField) --predicts--> CriticalExponents.

register!(
    Universality{:MeanField},
    CriticalExponents,
    Infinite;
    status=:universal,
    method=:analytic,
    reliability=:high,
    notes="Landau mean-field critical exponents (α=0, β=1/2, γ=1, δ=3, ν=1/2, η=0); exact for d >= d_c = 4.",
)
