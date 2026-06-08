# models/quantum/SSH/SSH_registry.jl — declarative implementation map for the
# Su-Schrieffer-Heeger dimerised chain (SSH 1979).  See `src/core/registry.jl`
# for the metadata schema.

# ── Energy (granularity-aware) ─────────────────────────────────────────
@register(
    SSH,
    Energy{:per_site},
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ssh.jl",
    references=["SSH1979"],
    notes="Half-filled per-site ε₀ by Gauss-Kronrod over the band |q(k)| = √(v²+w²+2vw cos k).",
)

# ── Spectrum / criticality ────────────────────────────────────────────
@register(
    SSH,
    MassGap,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ssh.jl",
    references=["SSH1979"],
    notes="Single-particle gap |v − w| = min_k|q(k)| at k = π (band gap is 2|v−w|).",
)
@register(
    SSH,
    MassGap,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ssh.jl",
    references=["SSH1979"],
    notes="Smallest non-negative single-particle eigenvalue (edge-mode splitting in topological phase).",
)
@register(
    SSH,
    EdgeModeEnergy,
    OBC,
    method=:dense_ed,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ssh.jl",
    references=["SSH1979", "AsbothOroszlanyPalyi2016"],
    notes="Same value as MassGap@OBC; named for the chiral edge-mode interpretation (exactly 0 at v = 0).",
)
@register(
    SSH,
    CorrelationLength,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ssh.jl",
    references=["SSH1979"],
    notes="ξ = 1/|v − w|; Inf on the gapless line |v| = |w|.",
)
@register(
    SSH,
    TopologicalInvariant,
    Infinite,
    method=:analytic,
    reliability=:high,
    tested_in="test/models/quantum/misc/test_ssh.jl",
    references=["SSH1979", "AsbothOroszlanyPalyi2016"],
    notes="Winding number W = (1/2π)∮ Im(q'/q) dk of q(k)=v+w e^{ik}; 1 (|w|>|v|, topological) / 0 (|w|<|v|, trivial).",
)
