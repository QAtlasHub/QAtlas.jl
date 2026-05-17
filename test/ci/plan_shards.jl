# test/ci/plan_shards.jl — emit a balanced GitHub-matrix shard plan.
#
#   julia test/ci/plan_shards.jl <N> [timings.tsv]
#
# Prints (stdout, last line) a JSON array for `matrix: include:` —
#   [{"sid":"s01","files":"d/a.jl,d/b.jl","aqua":"1"}, …]
#
# WHAT to run is the canonical universe (universe.jl, the single source
# of truth + completeness guard).  Timings only decide HOW to split:
#   * timings.tsv present  → Longest-Processing-Time bin-packing so all
#                            shards finish at ≈ the same wall-clock.
#   * absent / unreadable  → deterministic round-robin (identical to the
#                            QATLAS_TEST_SHARD k/N fallback) — never a
#                            leak, just not yet time-optimal.
# New / unknown files get a pessimistic estimate (P90 of known, or a
# fixed default if no history) so a surprise-heavy new test is isolated
# rather than piled onto an already-heavy shard.

include(joinpath(@__DIR__, "universe.jl"))   # ALL_TEST_FILES, test_file_key

const N = let a = get(ARGS, 1, "")
    n = tryparse(Int, a)
    (n !== nothing && n >= 1) ||
        error("plan_shards.jl: arg 1 must be N>=1; got $(repr(a))")
    n
end
const TIMINGS_PATH = get(ARGS, 2, "")

# universe → ordered "d/f" keys
const KEYS = [test_file_key(d, f) for (d, f) in ALL_TEST_FILES]

# Load timings TSV ("key\tseconds"); silently ignore missing/garbage
# rows — the plane is advisory and must degrade, never abort planning.
function load_timings(path)
    t = Dict{String,Float64}()
    (isempty(path) || !isfile(path)) && return t
    for ln in eachline(path)
        parts = split(strip(ln), '\t')
        length(parts) == 2 || continue
        v = tryparse(Float64, parts[2])
        v === nothing && continue
        t[String(parts[1])] = v
    end
    return t
end
const TIMES = load_timings(TIMINGS_PATH)

# Pessimistic default for files with no recorded time.
const DEFAULT_T = if isempty(TIMES)
    1.0
else
    s = sort(collect(values(TIMES)))
    s[clamp(ceil(Int, 0.9 * length(s)), 1, length(s))]   # P90 of known
end

est(k) = get(TIMES, k, DEFAULT_T)

# Assignment: LPT when we have history, else round-robin.
bins = [String[] for _ in 1:N]
loads = zeros(Float64, N)

if isempty(TIMES)
    for (i, k) in enumerate(KEYS)
        b = ((i - 1) % N) + 1
        push!(bins[b], k)
        loads[b] += est(k)
    end
    mode = "round-robin (no timing history)"
else
    for k in sort(KEYS; by=est, rev=true)         # longest first
        b = argmin(loads)                          # least-loaded bin
        push!(bins[b], k)
        loads[b] += est(k)
    end
    mode = "LPT bin-packing"
end

# Aqua → the finishing-earliest (least-loaded) shard, so it does not
# pile onto the critical-path shard.
const AQUA_BIN = argmin(loads)

# Emit JSON (ASCII test paths only ⇒ no escaping needed).
io = IOBuffer()
print(io, "[")
for b in 1:N
    b == 1 || print(io, ",")
    sid = "s" * lpad(string(b), 2, '0')
    print(io, "{\"sid\":\"", sid, "\",\"files\":\"", join(bins[b], ","), "\",")
    print(io, "\"aqua\":\"", b == AQUA_BIN ? "1" : "0", "\"}")
end
print(io, "]")
plan_json = String(take!(io))

# Human-readable summary → stderr (does not pollute the JSON stdout).
println(
    stderr,
    "plan_shards: N=$N  mode=$mode  files=$(length(KEYS))  " *
    "default_t=$(round(DEFAULT_T; digits=3))s  aqua→s$(lpad(string(AQUA_BIN), 2, '0'))",
)
for b in 1:N
    println(
        stderr,
        "  s$(lpad(string(b),2,'0')): $(length(bins[b])) files  " *
        "est=$(round(loads[b]; digits=1))s$(b == AQUA_BIN ? "  +aqua" : "")",
    )
end

println(plan_json)
