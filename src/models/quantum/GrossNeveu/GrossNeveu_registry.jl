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
    references=["Gross-Neveu 1974"],
    notes="UV free-fermion c = N at g = 0 (N Dirac flavours); g ≠ 0 raises DomainError (Phase 2: RG flow / dynamical mass).",
)
