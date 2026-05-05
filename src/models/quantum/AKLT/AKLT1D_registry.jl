# models/quantum/AKLT/AKLT1D_registry.jl —
# declarative implementation map for AKLT1D (S=1 BLBQ at the AKLT point).
#
# Closed-form analytical rows have `reliability=:high` because they
# follow from the AKLT 1988 projector construction and the
# Kennedy–Tasaki 1992 hidden-symmetry argument.  The Haldane gap row is
# `:medium` because it is an Östlund–Rommer 1995 numerical-exact value
# with no closed form.

# ── Infinite analytical rows (closed form) ───────────────────────────
@register(
    AKLT1D,
    GroundStateEnergyDensity,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["Affleck-Kennedy-Lieb-Tasaki 1988"],
    notes="Closed form e₀ = -2J/3 from the bond-projector decomposition of H.",
)
@register(
    AKLT1D,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["Affleck-Kennedy-Lieb-Tasaki 1988"],
    notes="Same analytic e₀ = -2J/3 routed through Energy(:per_site).",
)
@register(
    AKLT1D,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["Affleck-Kennedy-Lieb-Tasaki 1988"],
    notes="Closed form ξ = 1/log 3 ≈ 0.910 from VBS transfer matrix.",
)
@register(
    AKLT1D,
    StringOrderParameter,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    references=["Affleck-Kennedy-Lieb-Tasaki 1988", "Kennedy-Tasaki 1992"],
    notes="Closed form O_str = 4/9; detects hidden Z₂×Z₂ in the Haldane phase.",
)

# ── Infinite numerical-exact (no closed form) ────────────────────────
@register(
    AKLT1D,
    MassGap,
    Infinite,
    method=:literature_value,
    reliability=:medium,
    tested_in="test/standalone/test_aklt.jl",
    references=["Östlund-Rommer 1995"],
    notes="Haldane gap Δ ≈ 0.41048 J; DMRG numerical-exact, no closed form.",
)

# ── OBC dense ED (cap N ≤ 8 from _MAX_ED_SITES_S1) ───────────────────
@register(
    AKLT1D,
    ExactSpectrum,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/standalone/test_aklt.jl",
    notes="Full sorted spectrum from 3^N dense ED; N ≤ 8 (3^8 = 6561).",
)
