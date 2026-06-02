# Registry entries for FibonacciAnyons (refs #240)

@register(
    FibonacciAnyons,
    TopologicalEntanglementEntropy,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_fibonacci_anyons.jl",
    references=["FreedmanKitaevLarsenWang2003", "ReadRezayi1999", "KitaevPreskill2006"],
    notes="γ = (1/2) log(2 + φ) ≈ 0.6430; non-Abelian Fibonacci fusion, universal TQC."
)
