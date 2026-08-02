# ─────────────────────────────────────────────────────────────────────────────
# Toric Code — Kitaev (2003) Z₂ surface code, the canonical topological-order
# benchmark model.
#
# Hamiltonian (square lattice, S = 1/2 on every edge):
#
#   H = − J_e Σ_v A_v − J_m Σ_p B_p,
#
#   A_v = ∏_{i ∈ star(v)} σˣ_i      (vertex stabilizer, "electric")
#   B_p = ∏_{i ∈ ∂p}     σᶻ_i      (plaquette stabilizer, "magnetic")
#
# All stabilizers commute (each pair shares 0 or 2 edges), so H is exactly
# solvable simultaneously on the full {±1}-eigenspace decomposition of all
# A_v, B_p. The ground state has every A_v = B_p = +1.
#
# Closed-form results
# ───────────────────
#   • Ground-state energy density   ε₀ = −(J_e + J_m)
#       Square lattice has one vertex + one plaquette per edge in the planar
#       limit; we adopt the standard convention "per (vertex+plaquette) pair"
#       (one of each per unit cell), so ε₀ counts the contribution of one
#       saturated A_v plus one saturated B_p.
#
#   • Excitation gap                 Δ = 2 · min(J_e, J_m)
#       Flipping a single A_v eigenvalue (creating an `e` charge pair) costs
#       2 J_e; a single B_p flip (an `m` flux pair) costs 2 J_m. The minimum
#       is the single-anyon gap of the cheaper species; pair creation costs
#       2 × this when only one species is excited (still 2·min if pairs are
#       counted as the elementary excitation, see Kitaev 2003 §6).
#       Convention here: gap = energy to create a single bulk anyon = 2·min.
#
#   • Ground-state degeneracy        GSD(g) = 4^g on a closed surface of
#       genus g (Kitaev 2003 §5; equivalent to dim H₁(Σ_g; Z₂) = 2g, with
#       2^{2g} = 4^g logical states).  Torus (g = 1) gives 4-fold; double
#       torus 16; etc. Independent of J_e, J_m (purely topological).
#
#   • Topological entanglement entropy γ = log 2
#       Kitaev–Preskill 2006, Levin–Wen 2006: for a simply-connected disk
#       region, S(ρ_A) = α |∂A| − γ + O(|∂A|⁻¹) with γ = log 𝒟 and total
#       quantum dimension 𝒟 = √(Σ_a d_a²) = √4 = 2 for the Z₂ topological
#       order, so γ = log 2.
#
#   • Anyon content {1, e, m, ε = e×m}
#       Fusion rules: e×e = m×m = 1, e×m = ε, ε×ε = 1.
#       Self-statistics: e and m are bosons; ε is a fermion (statistical
#       phase π).
#       Mutual statistics: braiding e fully around m yields phase π
#       (mutual semion / Z₂ topological charge).
#
# Distinction from `KitaevHoneycomb`
# ──────────────────────────────────
# `KitaevHoneycomb` (Kitaev 2006, Annals 321) is a different model: spin-1/2
# on a 2D honeycomb lattice with anisotropic σˣσˣ, σʸσʸ, σᶻσᶻ couplings on
# the three bond types, solved by a four-Majorana mapping with Z₂ gauge
# fields. It hosts non-Abelian Ising anyons in the gapped phase under a
# magnetic-field perturbation. The ToricCode model here (Kitaev 2003,
# Annals 303) lives on the square lattice with stabilizer (4-body) terms,
# is exactly stabilizer-solvable without any fermionic mapping, and hosts
# Abelian Z₂ anyons. They share an author and a topological-order theme,
# nothing else.
#
# References
# ──────────
#   - A. Yu. Kitaev, "Fault-tolerant quantum computation by anyons",
#     [Kitaev2003](@cite).
#   - A. Kitaev, J. Preskill, "Topological entanglement entropy",
#     [KitaevPreskill2006](@cite).
#   - M. Levin, X.-G. Wen, "Detecting topological order in a ground state
#     wave function", [LevinWen2006](@cite).
#   - C. Nayak, S. H. Simon, A. Stern, M. Freedman, S. Das Sarma,
#     "Non-Abelian anyons and topological quantum computation",
#     [NayakSimonSternFreedmanDasSarma2008](@cite).
# ─────────────────────────────────────────────────────────────────────────────

# CONVENTION
#   Hamiltonian: Stabilizer / operator product
#   Observable:  Operator-product expectations (Wilson loops, GSD, TEE, S-matrix entries); convention-free
#   Reference:   docs/src/conventions.md §Topological / operator-product

"""
    ToricCode(; J_e::Real = 1.0, J_m::Real = 1.0) <: AbstractQAtlasModel

Kitaev (2003) toric code: the square-lattice Z₂ surface code with vertex
("electric") coupling `J_e` and plaquette ("magnetic") coupling `J_m`,

    H = − J_e Σ_v A_v − J_m Σ_p B_p,

where `A_v = ∏_{i ∈ star(v)} σˣ_i` and `B_p = ∏_{i ∈ ∂p} σᶻ_i`. All
stabilizers commute, so the model is exactly solvable. See module header
for the full closed-form result list and the distinction from
[`KitaevHoneycomb`](@ref).
"""
struct ToricCode <: AbstractQAtlasModel
    J_e::Float64
    J_m::Float64
end
function ToricCode(; J_e::Real=1.0, J_m::Real=1.0)
    J_e ≥ 0 || throw(DomainError(J_e, "ToricCode: J_e must be ≥ 0"))
    J_m ≥ 0 || throw(DomainError(J_m, "ToricCode: J_m must be ≥ 0"))
    return ToricCode(Float64(J_e), Float64(J_m))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy{:per_site} — Infinite (per vertex+plaquette pair)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ToricCode, ::Energy{:per_site}, ::Infinite) -> Float64

Ground-state energy density `ε₀ = −(J_e + J_m)`.

Closed-form: every stabilizer commutes with every other, so the ground
state simultaneously saturates `A_v = +1` ∀v and `B_p = +1` ∀p. The energy
contribution per (vertex + plaquette) unit cell is `−J_e − J_m`.
"""
function fetch(model::ToricCode, ::Energy{:per_site}, ::Infinite; kwargs...)
    return -(model.J_e + model.J_m)
end

# ═══════════════════════════════════════════════════════════════════════════════
# MassGap — Infinite
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ToricCode, ::MassGap, ::Infinite) -> Float64

Single-anyon excitation gap `Δ = 2·min(J_e, J_m)`.

Flipping a single vertex stabilizer eigenvalue (creating an `e` charge)
costs `2 J_e`; flipping a single plaquette (creating an `m` flux) costs
`2 J_m`. The minimum-energy excitation is therefore `2·min(J_e, J_m)`
— the gap of the cheaper anyon species.
"""
function fetch(model::ToricCode, ::MassGap, ::Infinite; kwargs...)
    return 2.0 * min(model.J_e, model.J_m)
end

# ═══════════════════════════════════════════════════════════════════════════════
# GroundStateDegeneracy — PBC (closed surface, genus g)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ToricCode, ::GroundStateDegeneracy, ::PBC; genus::Int = 1) -> Int

Ground-state degeneracy on a closed orientable surface of genus `g`,

    GSD(g) = 4^g.

This is the dimension of `H¹(Σ_g; Z₂) ⊗ H¹(Σ_g; Z₂)` (logical Pauli-X
and Pauli-Z operators, one independent pair per non-contractible cycle).
The torus (`g = 1`) gives the canonical 4-fold degeneracy. The result is
purely topological — it depends only on `genus`, not on `J_e`, `J_m` or
the lattice size.

`PBC` is the appropriate boundary tag because GSD is meaningful only on
a closed surface; on `OBC` / a disk the model has a unique ground state
(no homologically non-trivial loops), so the OBC method is intentionally
not registered.
"""
function fetch(::ToricCode, ::GroundStateDegeneracy, ::PBC; genus::Int=1, kwargs...)
    genus ≥ 0 || throw(DomainError(genus, "ToricCode GSD: genus must be ≥ 0"))
    return 4^genus
end

# ═══════════════════════════════════════════════════════════════════════════════
# TopologicalEntanglementEntropy — Infinite
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ToricCode, ::TopologicalEntanglementEntropy, ::Infinite) -> Float64

Topological entanglement entropy `γ = log 2`.

For the Z₂ topological order realised by the toric code, the total
quantum dimension is `𝒟 = √(Σ_a d_a²) = √(1² + 1² + 1² + 1²) = 2`
(four Abelian anyons `{1, e, m, ε}`, each with `d_a = 1`). The
Kitaev–Preskill (2006) / Levin–Wen (2006) prescription extracts
`γ = log 𝒟 = log 2` from the constant offset of the bipartite
entanglement entropy on a simply-connected region,

    S(ρ_A) = α |∂A| − γ + O(|∂A|⁻¹).

Independent of `J_e`, `J_m` (purely topological).

# Example

```jldoctest
julia> QAtlas.fetch(ToricCode(), TopologicalEntanglementEntropy(), Infinite())
0.6931471805599453
```
"""
function fetch(::ToricCode, ::TopologicalEntanglementEntropy, ::Infinite; kwargs...)
    return log(2.0)
end

# ═══════════════════════════════════════════════════════════════════════════════
# Anyon statistics — non-BC quantities (topological data only), one per return
# SCHEMA: self-statistics of an anyon, and mutual braiding of a pair (#819)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::ToricCode, ::AnyonSelfStatistics; label::Symbol = :e) -> NamedTuple

Self-statistics of one toric-code anyon.  Always returns the same fields:

    (label, statistics, self_phase, quantum_dim, fusion)

| `label`        | statistics | `self_phase` | fusion            |
| :------------- | :--------- | :----------- | :---------------- |
| `:1`, `:vacuum`| `:boson`   | `0`          | `1 x 1 -> 1`      |
| `:e`           | `:boson`   | `0`          | `e x e -> 1`      |
| `:m`           | `:boson`   | `0`          | `m x m -> 1`      |
| `:eps`, `:epsilon`, `:ε` | `:fermion` | `pi` | `eps x eps -> 1` |

All quantum dimensions are 1 (the theory is Abelian).  The bound state
`eps = e x m` is a fermion because the `e`/`m` mutual braid contributes `pi` —
see [`AnyonMutualStatistics`](@ref).

Throws an `ArgumentError` for an unknown label, listing the accepted ones.
"""
function fetch(::ToricCode, ::AnyonSelfStatistics; label::Symbol=:e, kwargs...)
    if label === :vacuum || label === Symbol("1")
        return _toric_self(Symbol("1"), :boson, 0.0)
    elseif label === :e
        return _toric_self(:e, :boson, 0.0)
    elseif label === :m
        return _toric_self(:m, :boson, 0.0)
    elseif label === :ε || label === :epsilon || label === :eps
        return _toric_self(:ε, :fermion, Float64(π))
    end
    return throw(
        ArgumentError(
            "ToricCode AnyonSelfStatistics: unknown anyon label :$label; expected " *
            "one of :vacuum / :1, :e, :m, :ε / :epsilon / :eps.  For the e/m " *
            "braiding phase use AnyonMutualStatistics, which is a different " *
            "quantity with different fields (#819).",
        ),
    )
end

# Every self-statistics row has the same shape; only the three values differ, and
# every toric-code anyon is its own inverse (a x a -> 1) with quantum dimension 1.
function _toric_self(label::Symbol, statistics::Symbol, self_phase::Float64)
    return (
        label=label,
        statistics=statistics,
        self_phase=self_phase,
        quantum_dim=1.0,
        fusion=(label, label) => Symbol("1"),
    )
end

"""
    fetch(::ToricCode, ::AnyonMutualStatistics; anyons = (:e, :m)) -> NamedTuple

Mutual braiding phase of a pair of toric-code anyons:

    (anyons, mutual_phase)

`(:e, :m)` (in either order) gives `pi` — a full braid of one around the other
multiplies the wave function by `exp(i pi) = -1`, the Z_2 "mutual semion"
statistics that makes `eps = e x m` a fermion.  Any pair involving the vacuum,
or an anyon with itself, braids trivially (`0`).

Throws an `ArgumentError` for an unknown anyon.
"""
function fetch(::ToricCode, ::AnyonMutualStatistics; anyons=(:e, :m), kwargs...)
    a, b = _toric_anyon.(anyons)
    phase = Set([a, b]) == Set([:e, :m]) ? Float64(π) : 0.0
    return (anyons=(a, b), mutual_phase=phase)
end

function _toric_anyon(s::Symbol)
    s in (:vacuum, Symbol("1")) && return Symbol("1")
    s in (:ε, :epsilon, :eps) && return :ε
    s in (:e, :m) && return s
    return throw(
        ArgumentError(
            "ToricCode AnyonMutualStatistics: unknown anyon :$s; expected one of " *
            ":vacuum / :1, :e, :m, :ε / :epsilon / :eps.",
        ),
    )
end

# BC-aware dispatch: the registry carries (ToricCode, <quantity>, Infinite) rows, so
# `has_native_fetch` requires these to exist.  The topological content is independent
# of any boundary tag, so both delegate to the no-BC implementations.
function fetch(model::ToricCode, q::AnyonSelfStatistics, ::Infinite; kwargs...)
    return fetch(model, q; kwargs...)
end
function fetch(model::ToricCode, q::AnyonMutualStatistics, ::Infinite; kwargs...)
    return fetch(model, q; kwargs...)
end
