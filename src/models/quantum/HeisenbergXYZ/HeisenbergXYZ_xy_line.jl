# ─────────────────────────────────────────────────────────────────────────────
# HeisenbergXYZ_xy_line.jl
#
# Closed-form ground-state energy density of the spin-½ XYZ chain on the
# XY anisotropic line (Jz = 0), following Lieb-Schultz-Mattis (1961).
#
# Hamiltonian (QAtlas convention, S = σ/2):
#
#     H = Σ_i [Jx S^x_i S^x_{i+1} + Jy S^y_i S^y_{i+1}]
#
# Jordan-Wigner + Bogoliubov diagonalisation gives the single-quasi-particle
# dispersion (in the σ-Hamiltonian form, scaled by S = σ/2 → ×1/4):
#
#     ε(k) = (1/4) √[(Jx + Jy)² cos²(k) + (Jx - Jy)² sin²(k)]
#          = (1/4) √[Jx² + Jy² + 2 Jx Jy cos(2k)]
#
# Ground-state energy density (per site, half-filled fermion band):
#
#     ε₀/site = -(1/π) ∫_0^{π/2} ε(k) dk
#             = -(1/(4π)) ∫_0^{π/2} √[Jx² + Jy² + 2 Jx Jy cos(2k)] dk
#
# Limit checks:
#   - Jx = Jy = J:  ε₀ = -J/π        (XX free-fermion point, also covered by XXZ delegation)
#   - Jx = 1, Jy → 0:  ε₀ → -1/(2π)     (anisotropic XY collapses to one-component cosine band)
# ─────────────────────────────────────────────────────────────────────────────

using QuadGK: quadgk

"""
    _heisenberg_xyz_gs_energy_xy_line(Jx, Jy) -> Float64

Ground-state energy density (per site) of the spin-½ XYZ chain on the
XY anisotropic line `Jz = 0`, computed in closed form from the Lieb-
Schultz-Mattis (1961) Jordan-Wigner / Bogoliubov diagonalisation.

Returns the QAtlas-convention value (S = σ/2 normalisation):

    ε₀ = -(1/(4π)) ∫_0^{π/2} √[Jx² + Jy² + 2 Jx Jy cos(2k)] dk.

# Limits

- `Jx = Jy = J`:    `ε₀ = -J/π`     (XX free-fermion point)
- `Jx > 0, Jy = 0`: `ε₀ = -Jx/π`    (single-component dispersion)

Throws `DomainError` if both `Jx == 0` and `Jy == 0` (no exchange).
"""
function _heisenberg_xyz_gs_energy_xy_line(Jx::Real, Jy::Real)
    (Jx == 0 && Jy == 0) &&
        throw(DomainError((Jx, Jy), "XY line requires at least one nonzero coupling"))
    integrand = k -> sqrt(Jx^2 + Jy^2 + 2 * Jx * Jy * cos(2k))
    val, _ = quadgk(integrand, 0.0, π / 2; rtol=1e-12)
    return -val / (4π)
end
