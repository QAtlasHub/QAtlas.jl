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

# QATLAS_TEST_GROUP="models/quantum/TFIM" runs only that dir (CI parallelism).
# Comma-separated for multi-dir groups. Unset or empty → run all.
const _test_group = get(ENV, "QATLAS_TEST_GROUP", "")
const dirs = if isempty(_test_group)
    ALL_DIRS
else
    groups = split(_test_group, ",")
    filter(d -> any(g -> startswith(d, g), groups), ALL_DIRS)
end
println("QATLAS_TEST_GROUP = $(repr(_test_group))  →  dirs = $(dirs)")

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

    # Aqua static QA — run only in the core group (or when running all).
    if isempty(_test_group) || any(g -> startswith("core/", g), split(_test_group, ","))
        @testset "test_aqua.jl" begin
            @time include(joinpath(@__DIR__, "test_aqua.jl"))
        end
    end

    @time for dir in dirs
        dirpath = joinpath(@__DIR__, dir)
        println("\nTest $(dirpath)")
        files = sort(
            filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(dirpath))
        )
        if isempty(files)
            println("  No test files found in $(dirpath).")
            @test false
        else
            for f in files
                @testset "$f" begin
                    filepath = joinpath(dirpath, f)
                    @time begin
                        println("  Including $(filepath)")
                        include(filepath)
                    end
                end
            end
        end
    end
end
