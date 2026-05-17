ENV["GKSwstype"] = "100"

using QAtlas, Test, LinearAlgebra, Lattice2D, ForwardDiff, Random
using SparseArrays, KrylovKit
using Aqua

const N_BLAS = min(Sys.CPU_THREADS, 64)
BLAS.set_num_threads(N_BLAS)
println("BLAS threads: $(BLAS.get_num_threads()) / $(Sys.CPU_THREADS) cores")

# Default "fast" test profile keeps every ED at N ≤ 12 (fits in a few
# GB of RAM) so PR CI runs in < 2 minutes on a 4-core GitHub runner.
# Set `QATLAS_TEST_FULL=1` to also run N ∈ {14, 16} sweeps (sparse +
# KrylovKit Lanczos, ~30 min on 128 GB / 36-core hardware).  Nightly
# cron only.
const QATLAS_TEST_FULL = get(ENV, "QATLAS_TEST_FULL", "0") != "0"
println("QATLAS_TEST_FULL = $(QATLAS_TEST_FULL)")

# ─────────────────────────────────────────────────────────────────────────────
# Canonical test-directory enumeration.
#
# `ALL_DIRS` is the *single source of truth* for which directories hold
# the suite.  CI parallelism is file-level sharding over the union of
# these dirs (see `QATLAS_TEST_SHARD` below), so the CI workflow no
# longer duplicates this list — it only passes a shard index.  The
# completeness guard (just below) makes it impossible to add a test
# directory that silently never runs.
# ─────────────────────────────────────────────────────────────────────────────
const ALL_DIRS = [
    "core/",
    "universalities/",
    "models/classical/",
    "models/quantum/TFIM/",
    "models/quantum/XXZ/",
    "models/quantum/Heisenberg/",
    "models/quantum/KitaevHoneycomb/",
    "models/quantum/misc/",
    "identities/",
    "verification/tightbinding/",
    "verification/tfim_ising/",
    "verification/heisenberg_xxz/",
    "verification/universality/",
]

_is_test_file(f) = startswith(f, "test_") && endswith(f, ".jl")

# ── Completeness guard ───────────────────────────────────────────────
# Walk the whole test tree; every directory that contains a `test_*.jl`
# MUST be enumerated in `ALL_DIRS`, and every `ALL_DIRS` entry must
# exist and be non-empty.  Runs in every shard (cheap) and fails
# loudly — a test directory can never be added without being wired in.
# `test_aqua.jl` at the test/ root is the only sanctioned exception.
let
    root = @__DIR__
    enumerated = Set(ALL_DIRS)
    discovered = Set{String}()
    for (d, _, files) in walkdir(root)
        any(_is_test_file, files) || continue
        rel = replace(relpath(d, root), '\\' => '/')
        rel == "." && continue  # test/ root: only test_aqua.jl, run specially
        push!(discovered, rel * "/")
    end
    leaked = sort(collect(setdiff(discovered, enumerated)))
    isempty(leaked) || error(
        "runtests.jl completeness guard: these on-disk test directories hold " *
        "test_*.jl files but are NOT in ALL_DIRS and would never run — add " *
        "them to ALL_DIRS: $(leaked)",
    )
    for d in ALL_DIRS
        p = joinpath(root, d)
        (isdir(p) && any(_is_test_file, readdir(p))) || error(
            "runtests.jl completeness guard: ALL_DIRS entry $(repr(d)) is " *
            "missing on disk or contains no test_*.jl files.",
        )
    end
end

# ── Canonical, deterministic global test-file universe ───────────────
# ALL_DIRS order × lexically-sorted files.  Sharding partitions THIS
# list, so the union of every shard is exactly this set — no file can
# be missed or double-run.
const ALL_TEST_FILES = let acc = Tuple{String,String}[]
    for d in ALL_DIRS
        for f in sort(filter(_is_test_file, readdir(joinpath(@__DIR__, d))))
            push!(acc, (d, f))
        end
    end
    acc
end

# ── Test selection: SHARD (CI) > GROUP (local targeted) > ALL ────────
#
#   QATLAS_TEST_SHARD="k/N"  — round-robin shard k of N over the global
#                              file list (CI parallelism; balanced).
#   QATLAS_TEST_GROUP="a,b"  — dir-prefix filter (local targeted runs).
#   neither                  — run everything.
#
const _shard_spec = get(ENV, "QATLAS_TEST_SHARD", "")
const _test_group = get(ENV, "QATLAS_TEST_GROUP", "")

const _selected, _mode_desc, _run_aqua = if !isempty(_shard_spec)
    parts = split(_shard_spec, "/")
    length(parts) == 2 ||
        error("QATLAS_TEST_SHARD must be \"k/N\"; got $(repr(_shard_spec))")
    k = parse(Int, parts[1])
    n = parse(Int, parts[2])
    (1 <= k <= n) || error("QATLAS_TEST_SHARD k/N needs 1 ≤ k ≤ N; got $k/$n")
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
    "files; aqua=$(_run_aqua)",
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

@testset "tests" begin
    test_args = copy(ARGS)
    println("Passed arguments ARGS = $(test_args) to tests.")

    # Aqua static QA — a one-shot whole-package check; run in exactly
    # one shard (shard 1), the core-covering group, or a full run.
    if _run_aqua
        @testset "test_aqua.jl" begin
            @time include(joinpath(@__DIR__, "test_aqua.jl"))
        end
    end

    @time for (d, f) in _selected
        filepath = joinpath(@__DIR__, d, f)
        @testset "$(d)$(f)" begin
            @time begin
                println("  Including $(filepath)")
                include(filepath)
            end
        end
    end
end
