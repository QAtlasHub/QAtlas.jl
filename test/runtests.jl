ENV["GKSwstype"] = "100"

using QAtlas, Test, LinearAlgebra, Lattice2D, ForwardDiff, Random
using SparseArrays, KrylovKit
using Aqua

const N_BLAS = min(Sys.CPU_THREADS, 64)
BLAS.set_num_threads(N_BLAS)
println("BLAS threads: $(BLAS.get_num_threads()) / $(Sys.CPU_THREADS) cores")

# Back-compat: legacy nightly switch.  Superseded by QATLAS_TEST_PROFILE
# (QATLAS_TEST_FULL=1 ⇒ profile=nightly when profile is unset).
const QATLAS_TEST_FULL = get(ENV, "QATLAS_TEST_FULL", "0") != "0"
println("QATLAS_TEST_FULL = $(QATLAS_TEST_FULL)")

# ── Test-volume profile (orthogonal to *which* files are selected) ───
#   fast    — PR merge gate: small N, coarse grids, loose tol; NO emit.
#   full    — push:main: larger N, finer grids, tight tol; emit timing
#             (and, later, evidence) — the heavier computation is what
#             gets persisted, so recorded numbers reflect the deep run.
#   nightly — cron: largest N sweeps, densest parameter grids.
# Individual test files read QATLAS_TEST_PROFILE to scale their work.
const QATLAS_TEST_PROFILE = let p = lowercase(get(ENV, "QATLAS_TEST_PROFILE", ""))
    if !isempty(p)
        p in ("fast", "full", "nightly") ||
            error("QATLAS_TEST_PROFILE must be fast|full|nightly; got $(repr(p))")
        Symbol(p)
    elseif QATLAS_TEST_FULL
        :nightly
    else
        :fast
    end
end
println("QATLAS_TEST_PROFILE = $(QATLAS_TEST_PROFILE)")

# Canonical universe + completeness guard (single source of truth,
# shared verbatim with the shard planner).
include(joinpath(@__DIR__, "ci", "universe.jl"))

# ── Test selection: FILES > SHARD > GROUP > ALL ──────────────────────
#
#   QATLAS_TEST_FILES="d/f.jl,d/g.jl" — explicit list emitted by the LPT
#       shard planner.  MUST be a subset of the canonical universe
#       (planner cannot smuggle in non-globbed files).
#   QATLAS_TEST_SHARD="k/N"  — round-robin shard (timing-agnostic; the
#       planner's fallback and a manual knob).
#   QATLAS_TEST_GROUP="a,b"  — dir-prefix filter (local targeted runs).
#   neither                  — run everything.
#
# Aqua (one-shot whole-package QA) runs in exactly one selection:
#   FILES → QATLAS_RUN_AQUA=1 ;  SHARD → shard 1 ;
#   GROUP → a core-covering group ;  ALL → yes.
const _test_files = get(ENV, "QATLAS_TEST_FILES", "")
const _shard_spec = get(ENV, "QATLAS_TEST_SHARD", "")
const _test_group = get(ENV, "QATLAS_TEST_GROUP", "")

const _selected, _mode_desc, _run_aqua = if !isempty(_test_files)
    want = [strip(x) for x in split(_test_files, ",") if !isempty(strip(x))]
    idx = Dict(test_file_key(d, f) => (d, f) for (d, f) in ALL_TEST_FILES)
    sel = Tuple{String,String}[]
    unknown = String[]
    for w in want
        haskey(idx, w) ? push!(sel, idx[w]) : push!(unknown, String(w))
    end
    isempty(unknown) || error(
        "QATLAS_TEST_FILES lists files outside the canonical universe " *
        "(the planner must only emit globbed files): $(unknown)",
    )
    (sel, "FILES (n=$(length(sel)))", get(ENV, "QATLAS_RUN_AQUA", "0") == "1")
elseif !isempty(_shard_spec)
    parts = split(_shard_spec, "/")
    length(parts) == 2 ||
        error("QATLAS_TEST_SHARD must be \"k/N\"; got $(repr(_shard_spec))")
    k = tryparse(Int, strip(parts[1]))
    n = tryparse(Int, strip(parts[2]))
    (k !== nothing && n !== nothing) ||
        error("QATLAS_TEST_SHARD must be integer \"k/N\"; got $(repr(_shard_spec))")
    (1 <= k <= n) || error("QATLAS_TEST_SHARD k/N needs 1 ≤ k ≤ N; got $k/$n")
    n <= length(ALL_TEST_FILES) || error(
        "QATLAS_TEST_SHARD N=$n exceeds the $(length(ALL_TEST_FILES))-file suite; " *
        "shards $(length(ALL_TEST_FILES) + 1)..$n would run zero tests — lower " *
        "the shard count N in .github/workflows/CI.yml.",
    )
    sel = [tf for (i, tf) in enumerate(ALL_TEST_FILES) if ((i - 1) % n) + 1 == k]
    (sel, "SHARD $k/$n", k == 1)
elseif !isempty(_test_group)
    groups = split(_test_group, ",")
    sel = [tf for tf in ALL_TEST_FILES if any(g -> startswith(tf[1], g), groups)]
    (sel, "GROUP $(repr(_test_group))", any(g -> startswith("core/", g), groups))
else
    (ALL_TEST_FILES, "ALL", true)
end

println(
    "Test selection: $(_mode_desc) → $(length(_selected))/$(length(ALL_TEST_FILES)) " *
    "files; profile=$(QATLAS_TEST_PROFILE); aqua=$(_run_aqua)",
)

const FIG_BASE = joinpath(pkgdir(QAtlas), "docs", "src", "assets")
const PATHS = Dict()
mkpath.(values(PATHS))

include(joinpath(@__DIR__, "util", "classical_partition.jl"))
include(joinpath(@__DIR__, "util", "tight_binding.jl"))
include(joinpath(@__DIR__, "util", "spinhalf_ed.jl"))
include(joinpath(@__DIR__, "util", "sparse_ed.jl"))
include(joinpath(@__DIR__, "util", "bloch.jl"))
include(joinpath(@__DIR__, "util", "tfim_dense_ed.jl"))
include(joinpath(@__DIR__, "util", "thermodynamic_identities.jl"))

# Per-file wall-time, captured for the timing plane (HOW-to-split).
const _TIMINGS = Dict{String,Float64}()

@testset "tests" begin
    test_args = copy(ARGS)
    println("Passed arguments ARGS = $(test_args) to tests.")

    if _run_aqua
        @testset "test_aqua.jl" begin
            _TIMINGS["__aqua__"] = @elapsed include(joinpath(@__DIR__, "test_aqua.jl"))
            println("  test_aqua.jl: $(round(_TIMINGS["__aqua__"]; digits=2)) s")
        end
    end

    @time for (d, f) in _selected
        filepath = joinpath(@__DIR__, d, f)
        key = test_file_key(d, f)
        @testset "$(key)" begin
            println("  Including $(filepath)")
            _TIMINGS[key] = @elapsed include(filepath)
            println("  $(key): $(round(_TIMINGS[key]; digits=2)) s")
        end
    end
end

# Emit per-shard timing as TSV (key<TAB>seconds).  Gated by QATLAS_EMIT
# so only push:main CI writes; PR/local runs never persist.  The
# consolidate job merges these into the `ci-timings` orphan branch; the
# planner LPT-bin-packs the next run from it (round-robin until then).
if get(ENV, "QATLAS_EMIT", "0") == "1"
    outdir = joinpath(@__DIR__, ".ci-out")
    mkpath(outdir)
    sid = if !isempty(_shard_spec)
        replace(_shard_spec, '/' => '-')
    elseif !isempty(_test_files)
        string(hash(_test_files); base=16)
    else
        "all"
    end
    open(joinpath(outdir, "timings-$(sid).tsv"), "w") do io
        for (k, v) in sort!(collect(_TIMINGS); by=first)
            println(io, k, '\t', round(v; digits=4))
        end
    end
    println("Emitted .ci-out/timings-$(sid).tsv ($(length(_TIMINGS)) entries)")
end
