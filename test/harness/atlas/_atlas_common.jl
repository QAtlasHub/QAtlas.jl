# test/harness/atlas/_atlas_common.jl
# Shared static-AST helpers used by BOTH AtlasInventory and AtlasRegistry.
# Plain functions (no module) so each `include`s it inside its own module
# and stays self-contained — no inter-module load-order coupling.

# Render a `references = [...]` / scalar AST node to the pipe-joined
# text form stored in INVENTORY / claim records.
function _refs_text(ex)
    ex === nothing && return ""
    if ex isa Expr && ex.head === :vect
        return join((x isa String ? x : string(x) for x in ex.args), " | ")
    end
    return string(ex)
end
