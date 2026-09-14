# ─────────────────────────────────────────────────────────────────────────────
# Percolation universality class
#
# d=2 exponent provenance:
#   All exact values from Coulomb gas / conformal field theory mapping:
#   β = 5/36, ν = 4/3, η = 5/24:
#     Nienhuis (1982) Phys. Rev. Lett. 49, 1062 — Coulomb gas derivation.
#     den Nijs (1979) J. Phys. A 12, 1857 — independent confirmation.
#   γ = 43/18: follows from Fisher relation γ = ν(2−η) = (4/3)(2−5/24).
#   δ = 91/5:  follows from Widom relation δ = 1 + γ/β.
#   α = −2/3:  follows from Rushbrooke relation α = 2 − 2β − γ.
#   Rigorous proofs for some exponents:
#     Smirnov, Werner (2001) Math. Res. Lett. 8, 729 — conformal
#     invariance of critical percolation in d=2.
#
# d=3: Wang, Zhou, Zhang, Garoni, Deng (2013) Phys. Rev. E 87, 052107
#       — large-scale Monte Carlo.  The paper quotes the RENORMALIZATION
#       exponents, 1/ν = 1.141 0(15) and y_h = 2.522 95(15); the standard
#       exponents follow as β/ν = d − y_h, γ/ν = 2y_h − d, δ = y_h/(d − y_h)
#       and η = 2 + d − 2y_h = −0.045 90(30).
#       β, γ, δ and η are pinned against those in
#       test/universalities/test_universality_critical_exponents_lit.jl.
#
#       α AND ν ARE NOT FROM THIS PAPER, and no source is recorded for either.
#       The paper quotes no α, and its only route from these numbers is
#       α = 2 − dν, which is Josephson; against that the shipped −0.625(3) sits
#       1.4σ out. And 1/1.1410(15) gives ν = 0.876 42(115) where src carries
#       0.876 19(12) — 2.0σ out in its own error, with an error ten times
#       tighter than this paper supports. Both are left as they are: replacing
#       them needs a source, not an algebra.
#
# d≥6: upper critical dimension; Toulouse (1974) mean-field exponents.
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:Percolation}, ::CriticalExponents; d) -> NamedTuple

Critical exponents of the ordinary (site/bond) percolation universality class.

- **d = 2**: Exact rational values (Coulomb gas mapping).
- **d = 3**: Numerical (Monte Carlo).
- **d ≥ 6**: Mean-field.
"""
function fetch(::Universality{:Percolation}, ::CriticalExponents; d::Int, kwargs...)
    if d == 2
        return (α=-2 // 3, β=5 // 36, γ=43 // 18, δ=91 // 5, ν=4 // 3, η=5 // 24)
    elseif d == 3
        return (
            α=-0.625,
            α_err=0.003,
            β=0.4181,
            β_err=0.0008,
            γ=1.793,
            γ_err=0.003,
            δ=5.29,
            δ_err=0.06,
            ν=0.87619,
            ν_err=0.00012,
            # NEGATIVE in d=3, from the paper's y_h: η = 2 + d − 2y_h.  Fisher
            # pins the sign here — γ/ν = 2 − η is 2.0464 > 2, so a positive η
            # contradicts the γ and ν above.
            η=-0.04590,
            η_err=0.00030,
        )
    elseif d >= 6
        # Percolation upper critical dimension is 6
        return (α=-1 // 1, β=1 // 1, γ=1 // 1, δ=2 // 1, ν=1 // 2, η=0 // 1)
    end
    return error("Percolation universality: d=$d not supported (d ∈ {2, 3, ≥6}).")
end
