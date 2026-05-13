# Registry entries for FibonacciAnyons (refs #240)

@register(
    FibonacciAnyons,
    TopologicalEntanglementEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_fibonacci_anyons.jl",
    references=[
        "Freedman-Kitaev-Larsen-Wang 2003", "Read-Rezayi 1999", "Kitaev-Preskill 2006"
    ],
    notes="γ = (1/2) log(2 + φ) ≈ 0.6431; non-Abelian Fibonacci fusion, universal TQC."
)
