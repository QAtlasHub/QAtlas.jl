# CONVENTION
#   Hamiltonian: Stabilizer / operator product
#   Observable:  Operator-product expectations (Wilson loops, GSD, TEE, S-matrix entries); convention-free
#   Reference:   docs/src/conventions.md §Topological / operator-product

"""
    FibonacciAnyons <: AbstractQAtlasModel

Non-Abelian Fibonacci anyon model.

The anyon spectrum has two charges `{1, τ}` with the fusion rule

    τ × τ = 1 + τ

and quantum dimensions

    d_1 = 1,   d_τ = φ = (1 + √5)/2  (golden ratio).

Total quantum dimension

    𝒟 = √(d_1² + d_τ²) = √(1 + φ²) = √(φ + 2)

(using `φ² = φ + 1`).

Fibonacci anyons are universal for topological quantum computation
(Freedman-Kitaev-Larsen-Wang 2003) and arise as edge excitations of
the Read-Rezayi Z_3 parafermion state (Read-Rezayi 1999).

This model has no continuous parameters — the spectrum and quantum
dimensions are fixed by the fusion rules.

# References

- M. Freedman, A. Kitaev, M. Larsen, Z. Wang, *Bull. Amer. Math. Soc.* **40**, 31 (2003).
- N. Read, E. Rezayi, *Phys. Rev. B* **59**, 8084 (1999).
- A. Kitaev, J. Preskill, *Phys. Rev. Lett.* **96**, 110404 (2006).
- M. Levin, X.-G. Wen, *Phys. Rev. Lett.* **96**, 110405 (2006).
"""
struct FibonacciAnyons <: AbstractQAtlasModel end

"""
    fetch(::FibonacciAnyons, ::TopologicalEntanglementEntropy, ::Infinite; kwargs...) -> Float64

Topological entanglement entropy of the Fibonacci-anyon non-Abelian
topological order:

    γ = log 𝒟 = (1/2) log(2 + φ)  ≈  0.6429653906

where 𝒟 = √(d_1² + d_τ²) = √(1 + φ²) = √(φ + 2) is the total quantum
dimension and φ = (1 + √5)/2 is the golden ratio (= d_τ, the
quantum dimension of the non-trivial τ anyon).

Fibonacci anyons are universal for topological quantum computation
(Freedman-Kitaev-Larsen-Wang 2003) and arise as edge excitations of
the Read-Rezayi Z_3 parafermion state (Read-Rezayi 1999).

# References

- M. Freedman, A. Kitaev, M. Larsen, Z. Wang, *Bull. Amer. Math. Soc.* **40**, 31 (2003).
- N. Read, E. Rezayi, *Phys. Rev. B* **59**, 8084 (1999).
- A. Kitaev, J. Preskill, *Phys. Rev. Lett.* **96**, 110404 (2006).
- M. Levin, X.-G. Wen, *Phys. Rev. Lett.* **96**, 110405 (2006).
"""
function fetch(::FibonacciAnyons, ::TopologicalEntanglementEntropy, ::Infinite; kwargs...)
    return 0.5 * log(2 + MathConstants.golden)
end
