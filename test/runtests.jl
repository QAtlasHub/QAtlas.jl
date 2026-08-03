ENV["GKSwstype"] = "100"

using QAtlas, Test, LinearAlgebra, Lattice2D, ForwardDiff, Random
using SparseArrays, KrylovKit
using Aqua
using TestShards

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

# Shared oracles and helpers: ONE copy, loaded by EVERY shard. These eleven includes must stay
# ABOVE the block below — inside it, `@shard` would make each a unit of its own, each would land
# on a single shard, and every test file on the other shards that calls one would fail.
#
# The other non-`test_` files under test/ — util/hubbard_ed.jl, generated/util_run_checks.jl,
# harness/atlas/*.jl — are NOT here on purpose: the test files that need them include them
# themselves, so they travel with their callers to whatever shard those land on.
include(joinpath(@__DIR__, "util", "classical_partition.jl"))
include(joinpath(@__DIR__, "util", "tight_binding.jl"))
include(joinpath(@__DIR__, "util", "spinhalf_ed.jl"))
include(joinpath(@__DIR__, "util", "sparse_ed.jl"))
include(joinpath(@__DIR__, "util", "bloch.jl"))
include(joinpath(@__DIR__, "util", "tfim_dense_ed.jl"))
include(joinpath(@__DIR__, "util", "thermodynamic_identities.jl"))
include(joinpath(@__DIR__, "util", "fluctuation_dissipation.jl"))
include(joinpath(@__DIR__, "util", "generic_ed.jl"))
include(joinpath(@__DIR__, "util", "extrapolate.jl"))
include(joinpath(@__DIR__, "util", "verify.jl"))

# Every `test_*.jl` under `test/`, in a deterministic order, each one its own shardable unit.
# `@shard` shadows `include` inside the block, so a unit is whatever this loop includes — a new
# file, or a whole new directory, is picked up BY BEING ON DISK. That is what the old
# `test/ci/universe.jl` completeness guard existed to enforce by hand, and it could only ever
# ERROR, because the wiring lived in a second list that could disagree with the tree. On a suite
# of this size that list was the thing most likely to drift.
#
# `test_aqua.jl` is picked up by the glob like anything else: it stops being a one-shot pinned to
# whichever shard carried the `aqua` flag and becomes an ordinary unit, balanced by measured
# runtime.
#
# Two rules when adding to this, and they are the only two:
#
#   1. SHARED FIXTURES GO ABOVE THIS BLOCK, as the eleven util/ includes do.
#   2. ANYTHING THAT IS NOT A `test_*.jl` FILE MUST BE NAMED. The glob does not error on what it
#      does not match; it silently stops running it.
#
# A bare `Pkg.test()` with nothing set in the environment runs all of it, in this order. Run one
# shard locally with `TESTSHARDS_ID=s3 TESTSHARDS_N=16 julia --project -e 'using Pkg; Pkg.test()'`.
TestShards.@shard begin
    for (dir, _, files) in sort!(collect(walkdir(@__DIR__)))
        for f in sort(files)
            startswith(f, "test_") && endswith(f, ".jl") || continue
            include(joinpath(dir, f))
        end
    end
end
