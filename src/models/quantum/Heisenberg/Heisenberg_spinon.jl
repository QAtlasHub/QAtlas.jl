# ─────────────────────────────────────────────────────────────────────────────
# Heisenberg1D — Infinite-volume spinon kinematics (Phase 1: closed-form)
#
# This module exposes three closed-form helpers and one Quantity-based
# `fetch` that complete the Tier-2 dynamic coverage of `Heisenberg1D` on
# `Infinite()`:
#
#   1. `heisenberg_spinon_dispersion(model, k)` — single-spinon dispersion
#      ε(k) = (π J / 2) |sin k|  (Faddeev–Takhtajan 1981).  Spinons are
#      half-odd-integer-spin excitations of the spin-½ XXX antiferromagnet
#      and are not generated singly; they always come in pairs in any
#      physical observable.  Public helper, exported from `QAtlas`.
#
#   2. `heisenberg_two_spinon_lower_edge(model, q)` — lower edge of the
#      two-spinon continuum at total momentum q ∈ [0, 2π],
#      ε_L(q) = (π J / 2) |sin q|  (des Cloizeaux–Pearson 1962).  By
#      construction this coincides with the single-spinon dispersion.
#      Public helper, exported from `QAtlas`.
#
#   3. `heisenberg_two_spinon_upper_edge(model, q)` — upper edge of the
#      two-spinon continuum,  ε_U(q) = π J |sin(q/2)|  (des
#      Cloizeaux–Pearson 1962).  Public helper, exported from `QAtlas`.
#
#   4. `fetch(model, ZZStructureFactor(), Infinite(); q, ω,
#       method = :muller, …)` — longitudinal dynamic structure factor
#      S^{zz}(q, ω) inside the two-spinon continuum, evaluated by the
#      Müller ansatz (Müller–Thomas–Beck–Bonner 1981):
#
#          S^{zz}(q, ω) ≈  Θ[ω − ε_L(q)] · Θ[ε_U(q) − ω]
#                          ----------------------------------
#                              2 √(ω² − ε_L(q)²)
#
#      and 0 outside the continuum.  This is an analytical approximation
#      with the correct support and the correct integrable square-root
#      singularity at the lower edge; it is *not* the exact
#      Caux–Hagemans (2006) result, which is reserved for Phase 2 (TODO,
#      see `docs/src/calc/heisenberg-spinons.md`).  `method` selects the
#      ansatz; `:muller` is the only currently-supported value.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Imports
# ═══════════════════════════════════════════════════════════════════════════════

using QuadGK: quadgk

# ═══════════════════════════════════════════════════════════════════════════════
# Single-spinon dispersion
# ═══════════════════════════════════════════════════════════════════════════════

"""
    heisenberg_spinon_dispersion(model::Heisenberg1D, k::Real; J::Real = 1.0)
        -> Float64

Single-spinon dispersion of the spin-½ antiferromagnetic Heisenberg
chain in the thermodynamic limit (Faddeev–Takhtajan 1981):

    ε(k) = (π J / 2) |sin k|,   k ∈ [0, π].

Spinons are massless half-odd-integer-spin excitations and the lower
edge of the two-spinon continuum coincides with this dispersion, see
[`heisenberg_two_spinon_lower_edge`](@ref).

Special values:
* `ε(0) = 0`             — gapless point.
* `ε(π/2) = π J / 2`     — band centre maximum.
* `ε(π) = 0`             — gapless Umklapp point.

`J` is passed by keyword because `Heisenberg1D` is parameterless in
this codebase (every other quantity threads `J` through kwargs the
same way).

# References
    L. D. Faddeev, L. A. Takhtajan, "What is the spin of a spin wave?",
      [FaddeevTakhtajan1981](@cite).
"""
function heisenberg_spinon_dispersion(model::Heisenberg1D, k::Real; J::Real=1.0)
    return (π * float(J) / 2) * abs(sin(float(k)))
end

# ═══════════════════════════════════════════════════════════════════════════════
# 2-spinon continuum edges (des Cloizeaux–Pearson 1962)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    heisenberg_two_spinon_lower_edge(model::Heisenberg1D, q::Real;
                                     J::Real = 1.0) -> Float64

Lower edge of the two-spinon continuum at total momentum `q` for the
spin-½ XXX AFM chain (des Cloizeaux–Pearson 1962):

    ε_L(q) = (π J / 2) |sin q|.

By construction `ε_L(q) ≡ ε(q)` where `ε` is the single-spinon
dispersion ([`heisenberg_spinon_dispersion`](@ref)) — this reflects
the kinematic configuration in which one of the two spinons sits at
zero energy, so the total energy equals the energy of the other.

Special values:
* `ε_L(0) = ε_L(π) = 0`  — gapless points (Umklapp included).
* `ε_L(π/2) = π J / 2`   — maximum of the lower edge.

# References
    J. des Cloizeaux, J. J. Pearson, "Spin-wave spectrum of the
      antiferromagnetic linear chain", [desCloizeauxPearson1962](@cite).
"""
function heisenberg_two_spinon_lower_edge(model::Heisenberg1D, q::Real; J::Real=1.0)
    return (π * float(J) / 2) * abs(sin(float(q)))
end

"""
    heisenberg_two_spinon_upper_edge(model::Heisenberg1D, q::Real;
                                     J::Real = 1.0) -> Float64

Upper edge of the two-spinon continuum at total momentum `q` for the
spin-½ XXX AFM chain (des Cloizeaux–Pearson 1962):

    ε_U(q) = π J |sin(q/2)|.

Special values:
* `ε_U(0) = 0`           — Goldstone-like vanishing of the support
                            window at q = 0.
* `ε_U(π) = π J`         — top of the continuum at the Umklapp point.

For `q ∈ (0, π]` one has `ε_U(q) ≥ ε_L(q)` strictly except at the
gapless points where the continuum collapses to a line.

# References
    J. des Cloizeaux, J. J. Pearson, "Spin-wave spectrum of the
      antiferromagnetic linear chain", [desCloizeauxPearson1962](@cite).
"""
function heisenberg_two_spinon_upper_edge(model::Heisenberg1D, q::Real; J::Real=1.0)
    return π * float(J) * abs(sin(float(q) / 2))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Müller-ansatz longitudinal dynamic structure factor
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _heisenberg_szz_muller(J, q, ω) -> Float64

Müller-ansatz approximation to the longitudinal dynamic structure
factor S^{zz}(q, ω) of the spin-½ XXX AFM chain.  The closed-form
expression, valid inside the two-spinon continuum, is

    S^{zz}(q, ω) ≈  Θ[ω - ε_L(q)] · Θ[ε_U(q) - ω]
                    ----------------------------------
                        2 √(ω² - ε_L(q)²)

with `ε_L(q) = (π J / 2) |sin q|` and
`ε_U(q) = π J |sin(q/2)|`.  Outside the closed continuum (i.e.
`ω ≤ ε_L` or `ω ≥ ε_U`) the expression returns exactly `0.0`.

The square-root singularity at the lower edge is integrable, so we
return the raw analytical value (no regulator); callers performing
ω-integrals are expected to use a routine that handles the integrable
singularity (e.g. Gauss–Chebyshev or substitution `ω² = ε_L² + s`).

This is *Phase 1* — the Müller ansatz captures the lower-edge
behaviour and the support correctly but mis-estimates the spectral
weight near the upper edge.  The *Phase 2* exact result is given by
the Caux–Hagemans (2006) algebraic Bethe ansatz form factor sum;
implementing it is tracked as a TODO in
`docs/src/calc/heisenberg-spinons.md`.

!!! warning "Longitudinal only — S^{xx}, S^{yy} are NOT equal to S^{zz} at q ≠ 0"
    Despite the SU(2) symmetry of the spin-½ XXX AFM chain, the
    individual transverse and longitudinal dynamic structure factors
    `S^{xx}(q, ω) = S^{yy}(q, ω)` and `S^{zz}(q, ω)` differ at any
    `q ≠ 0` — only their sum (the rotation-invariant trace
    `Σ_α S^{αα}`) is equal across spin axes.  The transverse Müller-
    style ansatz (Schulz 1986; Bougourzi-Karbach-Müller 1998) shares
    the two-spinon support `[ε_L(q), ε_U(q)]` but has a different
    singularity structure at the lower edge.  This method returns
    only the **longitudinal** `S^{zz}`; do not use it as a drop-in
    for `S^{xx}` at finite `q`.
"""
function _heisenberg_szz_muller(J::Real, q::Real, ω::Real)
    Jf = float(J)
    qf = float(q)
    ωf = float(ω)
    εL = (π * Jf / 2) * abs(sin(qf))
    εU = π * Jf * abs(sin(qf / 2))
    if ωf <= εL || ωf >= εU
        return 0.0
    end
    return 1 / (2 * sqrt(ωf^2 - εL^2))
end

"""
    fetch(model::Heisenberg1D, ::DynamicalSpinStructureFactor{:z,:z}, ::Infinite;
          q::Real, ω::Real, method::Symbol = :muller, J::Real = 1.0,
          kwargs...) -> Float64

Longitudinal dynamic structure factor `S^{zz}(q, ω)` of the
infinite-volume spin-½ XXX antiferromagnetic Heisenberg chain in the
thermodynamic limit.

The default `method = :muller` evaluates the Müller-ansatz closed
form inside the two-spinon continuum (see
`_heisenberg_szz_muller`):

    S^{zz}(q, ω) ≈  Θ[ω − ε_L(q)] · Θ[ε_U(q) − ω]
                    ----------------------------------
                        2 √(ω² − ε_L(q)²)

with edges defined by [`heisenberg_two_spinon_lower_edge`](@ref) and
[`heisenberg_two_spinon_upper_edge`](@ref).  Outside the continuum
the routine returns `0.0`.  No regulator is added at the lower edge:
the singularity is integrable (∝ 1/√(ω − ε_L) for ω → ε_L⁺), and any
quadrature scheme used downstream should treat it analytically rather
than relying on a numerical cap.

# Arguments
- `q::Real`: total momentum (radians).
- `ω::Real`: frequency (energy units consistent with `J`).
- `method::Symbol = :muller`: ansatz selector; only `:muller` is
  implemented in Phase 1.  Reserved future values:
  `:caux_hagemans` for the exact form-factor sum (Phase 2).
- `J::Real = 1.0`: Heisenberg coupling.

# References
    G. Müller, H. Thomas, H. Beck, J. C. Bonner, "Quantum spin
      dynamics of the antiferromagnetic linear chain in zero and
      nonzero magnetic field", [MullerThomasBeckBonner1981](@cite).
    J.-S. Caux, R. Hagemans, "The four-spinon dynamical structure
      factor of the Heisenberg chain", J. Stat. Mech. P12013 (2006)
      (Phase 2 reference, not yet implemented).
"""
function _heisenberg_szz_exact_I(t::Real)
    tf = float(t)
    integrand(x) = begin
        if x < 1e-8
            return -1.0 - tf^2 / 4
        end
        f_x =
            2 * ((1.0 + exp(-4x)) * cos(x * tf) - 2.0 * exp(-2x)) /
            (x * (1.0 - exp(-4x)) * (1.0 + exp(-2x)))
        f_sub = 2 * cos(x * tf) * (1.0 - exp(-x)) / x
        return f_x - f_sub
    end
    val, _ = quadgk(integrand, 0.0, 50.0; atol=1e-14, rtol=1e-14)
    return val + log(1.0 + 1.0 / tf^2)
end

function _heisenberg_szz_exact(J::Real, q::Real, ω::Real)
    Jf = float(J)
    qf = float(q)
    ωf = float(ω)
    εL = (π * Jf / 2) * abs(sin(qf))
    εU = π * Jf * abs(sin(qf / 2))

    if ωf <= εL || ωf >= εU
        return 0.0
    end

    ω_unit = ωf / Jf
    εL_unit = εL / Jf
    εU_unit = εU / Jf

    val_sqrt = sqrt((εU_unit^2 - εL_unit^2) / (ω_unit^2 - εL_unit^2))
    t = (4 / π) * log(val_sqrt + sqrt(val_sqrt^2 - 1))

    I = _heisenberg_szz_exact_I(t)
    M = 0.5 * exp(-I)

    return M / sqrt(εU^2 - ωf^2)
end

function fetch(
    ::Heisenberg1D,
    ::DynamicalSpinStructureFactor{:z,:z},
    ::Infinite;
    q::Real,
    ω::Real,
    method::Symbol=:muller,
    J::Real=1.0,
    kwargs...,
)
    J > 0 ||
        throw(DomainError(J, "Heisenberg1D ZZStructureFactor requires J > 0; got J = $J."))

    if method === :muller
        return _heisenberg_szz_muller(J, q, ω)
    elseif method === :exact_2spinon
        return _heisenberg_szz_exact(J, q, ω)
    elseif method === :caux_hagemans
        return error(
            "Heisenberg1D ZZStructureFactor Infinite: method=:caux_hagemans " *
            "(exact 2-spinon form-factor sum, Caux–Hagemans 2006) is reserved " *
            "for Phase 2 and not yet implemented.  Use method=:exact_2spinon for the " *
            "exact 2-spinon contribution, or method=:muller for the Phase-1 ansatz.",
        )
    else
        return error(
            "Heisenberg1D ZZStructureFactor Infinite: unknown method=:$(method); " *
            "supported: :muller (ansatz), :exact_2spinon (exact 2-spinon).",
        )
    end
end
