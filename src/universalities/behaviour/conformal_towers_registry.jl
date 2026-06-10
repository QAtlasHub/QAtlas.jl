# universalities/behaviour/conformal_towers_registry.jl
#
# Register ConformalTower for the relevant CFT classes.

@register(
    Universality{:Ising},
    ConformalTower,
    PBC,
    method=:analytic,
    reliability=:high,
    references=["Cardy1986"],
    notes="Conformal tower of states excitation energies E_n - E_0 = (2π v / L) Δ_n under PBC.",
)

@register(
    Universality{:Ising},
    ConformalTower,
    OBC,
    method=:analytic,
    reliability=:high,
    references=["Cardy1986", "BloteCardyNightingale1986"],
    notes="Conformal tower of states excitation energies E_n - E_0 = (π v / L) h_n under OBC.",
)

@register(
    Universality{:Heisenberg},
    ConformalTower,
    PBC,
    method=:analytic,
    reliability=:high,
    references=["Cardy1986", "Affleck1986"],
    notes="Conformal tower of states excitation energies E_n - E_0 = (2π v / L) Δ_n under PBC.",
)
