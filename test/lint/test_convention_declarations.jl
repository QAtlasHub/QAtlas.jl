# test/lint/test_convention_declarations.jl
#
# Lint: every quantum-model source file declaring `struct ... <: AbstractQAtlasModel`
# MUST carry a `# CONVENTION` header block near the top of the file, followed by
# at least two non-empty comment lines describing Hamiltonian/Observable choices.
#
# Spec: docs/src/conventions.md (landed in PR #438).
#
# Scope: src/models/quantum/  (recursive).
#
# Rules enforced
# ──────────────
#   1. The file MUST contain a line matching `^\s*#\s*CONVENTION\s*$`.
#   2. That CONVENTION line MUST be followed by ≥ 2 non-empty comment lines
#      (each starting with `#`).
#   3. The CONVENTION header MUST appear BEFORE the first `using` or
#      `struct ... <: AbstractQAtlasModel` declaration.
#
# Test failures print the offending file path (relative to the package root).

using Test

const _STRUCT_RE   = r"struct\s+\S+\s*(\{[^}]*\})?\s*<:\s*AbstractQAtlasModel"
const _CONV_RE     = r"^\s*#\s*CONVENTION\s*$"
const _USING_RE    = r"^\s*using\b"
const _COMMENT_RE  = r"^\s*#"
const _PKG_ROOT    = normpath(joinpath(@__DIR__, "..", ".."))
const _SCAN_ROOT   = joinpath(_PKG_ROOT, "src", "models", "quantum")

"Collect every `.jl` file under `src/models/quantum/`."
function _collect_model_files(root::AbstractString)
    files = String[]
    for (d, _, fs) in walkdir(root)
        for f in fs
            endswith(f, ".jl") && push!(files, joinpath(d, f))
        end
    end
    return sort!(files)
end

"Return (declares_model, has_convention, followed_by_two_comments, before_code) for `path`."
function _audit(path::AbstractString)
    lines = readlines(path)
    declares_model = any(occursin(_STRUCT_RE, ln) for ln in lines)

    conv_idx = findfirst(ln -> occursin(_CONV_RE, ln), lines)
    has_convention = conv_idx !== nothing

    followed_by_two_comments = false
    if has_convention
        c = 0
        for j in (conv_idx + 1):min(conv_idx + 10, length(lines))
            ln = lines[j]
            if occursin(_COMMENT_RE, ln) && !isempty(strip(ln)) &&
               strip(ln) != "#"
                c += 1
                if c >= 2
                    followed_by_two_comments = true
                    break
                end
            elseif isempty(strip(ln))
                continue
            else
                break  # non-comment, non-blank breaks the block
            end
        end
    end

    # Header MUST precede first `using` or model `struct` line.
    first_code_idx = findfirst(
        ln -> occursin(_USING_RE, ln) || occursin(_STRUCT_RE, ln),
        lines,
    )
    before_code = if has_convention && first_code_idx !== nothing
        conv_idx < first_code_idx
    else
        has_convention  # if no code line found, treat as OK
    end

    return (; declares_model, has_convention, followed_by_two_comments, before_code)
end

@testset "convention declarations" begin
    @test isdir(_SCAN_ROOT)

    files = _collect_model_files(_SCAN_ROOT)
    @test !isempty(files)

    missing_header   = String[]
    short_header     = String[]
    misplaced_header = String[]

    for f in files
        rep = _audit(f)
        rep.declares_model || continue
        rel = relpath(f, _PKG_ROOT)
        if !rep.has_convention
            push!(missing_header, rel)
        else
            rep.followed_by_two_comments || push!(short_header, rel)
            rep.before_code             || push!(misplaced_header, rel)
        end
    end

    if !isempty(missing_header)
        @error "Files declaring AbstractQAtlasModel struct but missing `# CONVENTION` header" missing_header
    end
    if !isempty(short_header)
        @error "Files with `# CONVENTION` header followed by < 2 non-empty comment lines" short_header
    end
    if !isempty(misplaced_header)
        @error "Files where `# CONVENTION` header appears after first `using`/`struct` line" misplaced_header
    end

    @test isempty(missing_header)
    @test isempty(short_header)
    @test isempty(misplaced_header)
end
