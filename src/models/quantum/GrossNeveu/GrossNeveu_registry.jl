# models/quantum/GrossNeveu/GrossNeveu_registry.jl
#
# 1+1-D Gross-Neveu (asymptotically free 4-fermion).  Phase-1 UV-only.

@register(
    GrossNeveu,
    CentralCharge,
    Infinite,
    method=:analytic_uv,
    reliability=:high,
    tested_in="test/standalone/test_gross_neveu.jl",
    references=["GrossNeveu1974"],
    notes="UV free-fermion c = N at g = 0 (N Dirac flavours); g ≠ 0 raises DomainError (Phase 2: RG flow / dynamical mass).",
)

@register(
    GrossNeveu,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_gross_neveu.jl",
    references=["GrossNeveu1974", "AndreiLowenstein1979"],
    notes="Large-N dynamic mass m_F = Λ exp(-π/(N g²)); Λ required kwarg (renormalisation scheme).",
)
