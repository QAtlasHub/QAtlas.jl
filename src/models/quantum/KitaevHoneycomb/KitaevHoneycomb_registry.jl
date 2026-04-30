# models/quantum/KitaevHoneycomb/KitaevHoneycomb_registry.jl — declarative
# implementation map for the Kitaev honeycomb model.  Mirrors the layout of
# `TFIM_registry.jl`.  See `src/core/registry.jl` for the metadata schema.

# ── Energy (granularity-aware): GS-only PBC + GS/thermal Infinite, OBC ──
@register(
    KitaevHoneycomb,
    Energy{:per_site},
    Infinite,
    method=:matter_free_fermion,
    reliability=:high,
    tested_in="test/models/test_KitaevHoneycomb.jl",
    references=["Kitaev 2006", "Lieb 1994"],
    notes="GS by 2D Gauss-Kronrod; thermal value (β kwarg) uses matter-sector free fermion.",
)
@register(
    KitaevHoneycomb,
    Energy{:per_site},
    OBC,
    method=:matter_free_fermion,
    reliability=:high,
    tested_in="test/models/test_KitaevHoneycomb.jl",
    references=["Kitaev 2006"],
    notes="Bipartite hopping matrix SVD; thermal value (β kwarg) uses matter-sector free fermion.",
)
@register(
    KitaevHoneycomb,
    Energy{:per_site},
    PBC,
    method=:matter_free_fermion,
    reliability=:high,
    tested_in="test/models/test_KitaevHoneycomb.jl",
    references=["Kitaev 2006"],
    notes="Four-sector minimum on Lx × Ly torus (T = 0 only).",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    KitaevHoneycomb,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/test_KitaevHoneycomb.jl",
    references=["Kitaev 2006"],
    notes="Δ = 2 · max(|K_γ_max| − sum_others, 0); 0 in B-phase, finite in A_γ.",
)

# ── Free-fermion thermal (per-site) — meta-defined in
#    `KitaevHoneycomb_thermal.jl`.  Matter-sector approximation: valid for
#    T ≪ Δ_v (flux gap).  See file header for the physical scope.
@register(
    KitaevHoneycomb,
    FreeEnergy,
    Infinite,
    method=:matter_free_fermion,
    reliability=:medium,
    tested_in="test/models/test_KitaevHoneycomb_thermal.jl",
    references=["Kitaev 2006", "Lieb 1994"],
    notes="Matter-sector free-fermion approximation; valid below flux gap T ≪ Δ_v.",
)
@register(
    KitaevHoneycomb,
    FreeEnergy,
    OBC,
    method=:matter_free_fermion,
    reliability=:medium,
    tested_in="test/models/test_KitaevHoneycomb_thermal.jl",
    references=["Kitaev 2006", "Lieb 1994"],
    notes="Matter-sector free-fermion approximation; valid below flux gap T ≪ Δ_v.",
)
@register(
    KitaevHoneycomb,
    ThermalEntropy,
    Infinite,
    method=:matter_free_fermion,
    reliability=:medium,
    tested_in="test/models/test_KitaevHoneycomb_thermal.jl",
    references=["Kitaev 2006", "Lieb 1994"],
    notes="Matter-sector free-fermion approximation; valid below flux gap T ≪ Δ_v.",
)
@register(
    KitaevHoneycomb,
    ThermalEntropy,
    OBC,
    method=:matter_free_fermion,
    reliability=:medium,
    tested_in="test/models/test_KitaevHoneycomb_thermal.jl",
    references=["Kitaev 2006", "Lieb 1994"],
    notes="Matter-sector free-fermion approximation; valid below flux gap T ≪ Δ_v.",
)
@register(
    KitaevHoneycomb,
    SpecificHeat,
    Infinite,
    method=:matter_free_fermion,
    reliability=:medium,
    tested_in="test/models/test_KitaevHoneycomb_thermal.jl",
    references=["Kitaev 2006", "Lieb 1994"],
    notes="Matter-sector free-fermion approximation; valid below flux gap T ≪ Δ_v.",
)
@register(
    KitaevHoneycomb,
    SpecificHeat,
    OBC,
    method=:matter_free_fermion,
    reliability=:medium,
    tested_in="test/models/test_KitaevHoneycomb_thermal.jl",
    references=["Kitaev 2006", "Lieb 1994"],
    notes="Matter-sector free-fermion approximation; valid below flux gap T ≪ Δ_v.",
)
