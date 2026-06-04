# ─────────────────────────────────────────────────────────────────────────────
# O(n) model universality classes: XY (n=2), Heisenberg (n=3)
#
# In d=2, Mermin-Wagner forbids spontaneous symmetry breaking for
# continuous symmetry → no standard power-law exponents.
# In d=3, high-precision numerical estimates are available.
#
# References:
#   XY d=3: Chester, Landry, Liu, Poland, Simmons-Duffin, Su, Vichi,
#     JHEP 02, 098 (2020) — conformal bootstrap.
#   Heisenberg d=3: Chester, Landry, Liu, Poland, Simmons-Duffin, Su,
#     Vichi, Phys. Rev. D 104, 105013 (2021).
# ─────────────────────────────────────────────────────────────────────────────

"""
    fetch(::Universality{:XY}, ::CriticalExponents; d) -> NamedTuple

Critical exponents of the XY (O(2)) universality class.

- **d = 2**: BKT transition — no standard power-law exponents.
  Returns η(T_c) = 1/4 only.
- **d = 3**: Conformal bootstrap (Chester et al. 2020).
- **d ≥ 4**: Mean-field.
"""
function fetch(::Universality{:XY}, ::CriticalExponents; d::Int, kwargs...)
    if d == 2
        # BKT: no standard critical exponents; only η at T_c is universal.
        return (η=1 // 4, note="BKT transition — no standard power-law exponents")
    elseif d == 3
        return (
            α=-0.01526,
            α_err=0.00030,
            β=0.34869,
            β_err=0.00007,
            γ=1.3179,
            γ_err=0.0002,
            δ=4.77937,
            δ_err=0.00025,
            ν=0.67175,
            ν_err=0.00010,
            η=0.038176,
            η_err=0.000044,
        )
    elseif d >= 4
        return fetch(MeanField(), CriticalExponents())
    end
    return error("XY universality: d=$d not supported (d ∈ {2, 3, ≥4}).")
end

"""
    fetch(::Universality{:Heisenberg}, ::CriticalExponents; d) -> NamedTuple

Critical exponents of the Heisenberg (O(3)) universality class.

- **d ≤ 2**: Mermin-Wagner theorem — no spontaneous order at T > 0.
- **d = 3**: Conformal bootstrap / Monte Carlo.
- **d ≥ 4**: Mean-field.
"""
function fetch(::Universality{:Heisenberg}, ::CriticalExponents; d::Int, kwargs...)
    if d <= 2
        return error(
            "Heisenberg universality in d=$d: Mermin-Wagner theorem " *
            "prohibits spontaneous breaking of continuous symmetry.",
        )
    elseif d == 3
        return (
            α=-0.1336,
            α_err=0.0015,
            β=0.3689,
            β_err=0.0003,
            γ=1.3960,
            γ_err=0.0009,
            δ=4.783,
            δ_err=0.003,
            ν=0.7112,
            ν_err=0.0005,
            η=0.0375,
            η_err=0.0005,
        )
    elseif d >= 4
        return fetch(MeanField(), CriticalExponents())
    end
    return error("Heisenberg universality: d=$d not supported.")
end

# ─── Infinite-bc forwarding for verify() integration ───────────────────────
function fetch(m::Universality{:XY}, q::CriticalExponents, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
function fetch(m::Universality{:Heisenberg}, q::CriticalExponents, ::Infinite; kwargs...)
    fetch(m, q; kwargs...)
end
