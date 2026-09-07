# models/quantum/RiceMele/RiceMele_registry.jl — declarative implementation map
# for the Rice-Mele dimerised chain with a staggered potential (Rice & Mele 1982).
# See `src/core/registry.jl` for the metadata schema.

# ── Energy (granularity-aware) ─────────────────────────────────────────
@register(
    RiceMele,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ricemele.jl",
    references=["RiceMele1982"],
    notes="Half-filled per-site ε₀ = −(1/4π)∫E₊(k)dk by Gauss-Kronrod; " *
          "E₊(k) = √(v²+w²+2vw cos k + Δ²).",
)
@register(
    RiceMele,
    Energy{:per_site},
    OBC,
    method=:single_particle_diagonalization,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ricemele.jl",
    references=["RiceMele1982"],
    notes="N lowest of 2N single-particle levels over 2N sites; carries the open-boundary " *
          "correction, so it approaches Infinite as O(1/N).",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    RiceMele,
    ExactSpectrum,
    OBC,
    method=:single_particle_diagonalization,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ricemele.jl",
    references=["RiceMele1982"],
    notes="ALL 2N single-particle energies, not SSH's non-negative half: the OBC ± symmetry " *
          "holds only while the unit cells are complete.",
)
@register(
    RiceMele,
    MassGap,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ricemele.jl",
    references=["RiceMele1982"],
    notes="min_k E₊(k) = √((|v|−|w|)² + Δ²), Fermi level to band edge; the band gap is 2×. " *
          "SSH's convention, to which it reduces at Δ = 0. Never zero for Δ ≠ 0.",
)
@register(
    RiceMele,
    MassGap,
    OBC,
    method=:single_particle_diagonalization,
    cost=:polynomial,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ricemele.jl",
    references=["RiceMele1982"],
    notes="Half the HOMO-LUMO difference (ε_{N+1} − ε_N)/2 at half filling — the Infinite " *
          "convention. Agrees with SSH's smallest-non-negative phrasing at every Δ. Does " *
          "NOT converge to the Infinite value when |w| > |v|: there it is the edge-mode " *
          "scale, flat in N.",
)
@register(
    RiceMele,
    CorrelationLength,
    Infinite,
    method=:analytic,
    cost=:closed_form,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ricemele.jl",
    references=["RiceMele1982"],
    notes="ξ = 1/√((|v|−|w|)² + Δ²). Finite for every Δ ≠ 0, unlike SSH on |v| = |w|.",
)

# NO TopologicalInvariant row: Δ breaks the symmetry that quantises SSH's winding number.
# What is quantised here is the Thouless pump, a property of a PATH in (v−w, Δ).
