using Aqua

# Aqua.jl: static QA checks (ambiguities, unbound args, stale deps,
# undocumented public names, persistent tasks, piracy, project.toml).
#
# `deps_compat` / `stale_deps` / `undocumented_names` family is the set
# of checks that genuinely catch drift between `Project.toml` and the
# source tree — keep them on.
#
# `ambiguities` is hard-fail (closed in #100) — `Aqua.detect_ambiguities`
# on the QAtlas module currently returns `[]`, and any new method added
# in src/ that introduces an ambiguity should be caught at PR time.
# Stdlib-inherited ambiguities (LinearAlgebra + QuadGK) live outside
# QAtlas's own method tables and therefore do not surface from
# `recursive=false`.

@testset "Aqua QA" begin
    Aqua.test_all(
        QAtlas;
        ambiguities=false,
        deps_compat=true,
        stale_deps=true,
        # Re-enabled in v0.13.5 after every public struct got a
        # dedicated docstring. Keeps the public-API surface honest:
        # adding a new `export X` without a docstring for `X` now
        # fails CI.
        undocumented_names=true,
        persistent_tasks=false,
        # `fetch` is AbstractQAtlas's generic function; QAtlas is its concrete
        # implementer, so the legacy `fetch(::Symbol, ::Symbol, ::BC)` deprecation
        # shim is not piracy in spirit — treat the shared `fetch` as own (#734).
        piracies=(treat_as_own=[QAtlas.fetch],),
        unbound_args=true,
    )

    @testset "ambiguities (hard)" begin
        amb = Aqua.detect_ambiguities(QAtlas; recursive=false)
        if !isempty(amb)
            @error "Aqua: ambiguities introduced in QAtlas — fix or document" amb
        end
        @test isempty(amb)
    end
end
