using QAtlas
using Documenter
using Downloads

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)
favicon_path = joinpath(assets_dir, "favicon.ico")
logo_path = joinpath(assets_dir, "logo.png")

Downloads.download("https://github.com/sotashimozono.png", favicon_path)
Downloads.download("https://github.com/sotashimozono.png", logo_path)

makedocs(;
    sitename="QAtlas.jl",
    repo=Remotes.GitHub("sotashimozono", "QAtlas.jl"),
    format=Documenter.HTML(;
        canonical="https://codes.sota-shimozono.com/QAtlas.jl/stable/",
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        # `index.md` includes a full `@autodocs` of every QAtlas binding; the
        # generated HTML grew past the default 200 KiB threshold once the
        # observable surface expanded in v0.17/0.18 (Tier 1 + Tier 2 brought
        # ~80 new fetch methods).  Bump the size threshold instead of
        # splitting the autodocs page (the user explicitly preserved the
        # `Modules = [QAtlas]` API layout).
        size_threshold=600_000,
        size_threshold_warn=300_000,
        mathengine=MathJax3(
            Dict(
                :tex => Dict(
                    :inlineMath => [["\$", "\$"], ["\\(", "\\)"]],
                    :tags => "ams",
                    :packages => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
        assets=["assets/favicon.ico", "assets/custom.css", "assets/report-issue.js"],
    ),
    modules=[QAtlas],
    pages=[
        "Home" => "index.md",
        "Models" => [
            "models/index.md",
            "Classical" => [
                "models/classical/index.md",
                "Ising Square" => "models/classical/ising-square.md",
                "Ising Triangular" => "models/classical/ising-triangular.md",
            ],
            "Quantum" => [
                "models/quantum/index.md",
                "TFIM" => "models/quantum/tfim.md",
                "Heisenberg" => "models/quantum/heisenberg.md",
                "Majumdar-Ghosh" => "models/quantum/majumdar_ghosh.md",
                "XXZ" => "models/quantum/xxz.md",
                "Hubbard1D" => "models/quantum/hubbard1d.md",
                "Kitaev Honeycomb" => "models/quantum/kitaev-honeycomb.md",
                "Kitaev1D" => "models/quantum/kitaev1d.md",
                "Tight-Binding" => [
                    "models/quantum/tightbinding/index.md",
                    "Honeycomb" => "models/quantum/tightbinding/honeycomb.md",
                    "Kagome" => "models/quantum/tightbinding/kagome.md",
                    "Lieb" => "models/quantum/tightbinding/lieb.md",
                    "Triangular" => "models/quantum/tightbinding/triangular.md",
                ],
            ],
        ],
        "Universality Classes" => [
            "universalities/index.md",
            "Ising" => "universalities/ising.md",
            "Percolation" => "universalities/percolation.md",
            "Potts" => "universalities/potts.md",
            "KPZ" => "universalities/kpz.md",
            "XY / Heisenberg" => "universalities/on-models.md",
            "Mean-Field" => "universalities/mean-field.md",
            "E8" => "universalities/e8.md",
            "Cardy Entanglement" => "universalities/cardy_entanglement.md",
        ],
        "Verification" => [
            "verification/index.md",
            "Cross-Checks" => "verification/cross-checks.md",
            "Entanglement" => "verification/entanglement.md",
            "Disordered" => "verification/disordered.md",
            "Identity Harness" => "verification/identity-harness.md",
        ],
        "Methods" => [
            "methods/index.md",
            "Physical" => [
                "Transfer Matrix" => "methods/transfer-matrix/index.md",
                "Bloch Hamiltonian" => "methods/bloch-hamiltonian/index.md",
                "Jordan-Wigner" => "methods/jordan-wigner/index.md",
                "Calabrese-Cardy" => "methods/calabrese-cardy/index.md",
            ],
            "Computational" => [
                "Exact Diagonalization" => "methods/exact-diagonalization/index.md",
                "Automatic Differentiation" => "methods/automatic-differentiation/index.md",
            ],
        ],
        "Derivation Notes" => [
            "JW → TFIM BdG" => "calc/jw-tfim-bdg.md",
            "Kramers-Wannier Duality" => "calc/kramers-wannier-duality.md",
            "Transfer Matrix Split" => "calc/transfer-matrix-symmetric-split.md",
            "Yang Magnetization" => "calc/yang-magnetization-toeplitz.md",
            "Heisenberg Dimer" => "calc/heisenberg-dimer-singlet-triplet.md",
            "Bethe Ansatz e₀" => "calc/bethe-ansatz-heisenberg-e0.md",
            "Heisenberg Spinons" => "calc/heisenberg-spinons.md",
            "XXZ Luttinger Parameters" => "calc/xxz-luttinger-parameters.md",
            "Honeycomb Bloch" => "calc/bloch-honeycomb-dispersion.md",
            "Kagome Flat Band" => "calc/bloch-kagome-flat-band.md",
            "Lieb Flat Band" => "calc/bloch-lieb-flat-band.md",
            "Calabrese-Cardy OBC/PBC" => "calc/calabrese-cardy-obc-vs-pbc.md",
            "TFIM Entanglement (Peschel)" => "calc/tfim-entanglement-peschel.md",
            "TFIM Loschmidt + DQPT" => "calc/tfim-loschmidt.md",
            "TFIM GGE Quench" => "calc/tfim-gge.md",
            "AD from ln Z" => "calc/ad-thermodynamics-from-z.md",
            "Scaling Relations" => "calc/ising-scaling-relations.md",
            "Ising CFT Operators" => "calc/ising-cft-primary-operators.md",
            "Ising CFT + σ → E8" => "calc/ising-cft-magnetic-perturbation.md",
            "E8 Mass Derivation" => "calc/e8-mass-spectrum-derivation.md",
        ],
    ],
)

deploydocs(; repo="github.com/sotashimozono/QAtlas.jl.git", devbranch="main")
using QAtlas
using Documenter
using Downloads

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)
favicon_path = joinpath(assets_dir, "favicon.ico")
logo_path = joinpath(assets_dir, "logo.png")

Downloads.download("https://github.com/sotashimozono.png", favicon_path)
Downloads.download("https://github.com/sotashimozono.png", logo_path)

makedocs(;
    sitename="QAtlas.jl",
    repo=Remotes.GitHub("sotashimozono", "QAtlas.jl"),
    format=Documenter.HTML(;
        canonical="https://codes.sota-shimozono.com/QAtlas.jl/stable/",
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="main",
        # `index.md` includes a full `@autodocs` of every QAtlas binding; the
        # generated HTML grew past the default 200 KiB threshold once the
        # observable surface expanded in v0.17/0.18 (Tier 1 + Tier 2 brought
        # ~80 new fetch methods).  Bump the size threshold instead of
        # splitting the autodocs page (the user explicitly preserved the
        # `Modules = [QAtlas]` API layout).
        size_threshold=600_000,
        size_threshold_warn=300_000,
        mathengine=MathJax3(
            Dict(
                :tex => Dict(
                    :inlineMath => [["\$", "\$"], ["\\(", "\\)"]],
                    :tags => "ams",
                    :packages => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
        assets=["assets/favicon.ico", "assets/custom.css", "assets/report-issue.js"],
    ),
    modules=[QAtlas],
    pages=[
        "Home" => "index.md",
        "Models" => [
            "models/index.md",
            "Classical" => [
                "models/classical/index.md",
                "Ising Square" => "models/classical/ising-square.md",
            ],
            "Quantum" => [
                "models/quantum/index.md",
                "TFIM" => "models/quantum/tfim.md",
                "Heisenberg" => "models/quantum/heisenberg.md",
                "XXZ" => "models/quantum/xxz.md",
                "Kitaev Honeycomb" => "models/quantum/kitaev-honeycomb.md",
                "Toric Code" => "models/quantum/toric-code.md",
                "Tight-Binding" => [
                    "models/quantum/tightbinding/index.md",
                    "Honeycomb" => "models/quantum/tightbinding/honeycomb.md",
                    "Kagome" => "models/quantum/tightbinding/kagome.md",
                    "Lieb" => "models/quantum/tightbinding/lieb.md",
                    "Triangular" => "models/quantum/tightbinding/triangular.md",
                ],
            ],
        ],
        "Universality Classes" => [
            "universalities/index.md",
            "Ising" => "universalities/ising.md",
            "Percolation" => "universalities/percolation.md",
            "Potts" => "universalities/potts.md",
            "KPZ" => "universalities/kpz.md",
            "RMT / Poisson" => "universalities/rmt.md",
            "XY / Heisenberg" => "universalities/on-models.md",
            "Mean-Field" => "universalities/mean-field.md",
            "E8" => "universalities/e8.md",
            "CFT Casimir Correction" => "universalities/cft-casimir.md",
        ],
        "Verification" => [
            "verification/index.md",
            "Cross-Checks" => "verification/cross-checks.md",
            "Entanglement" => "verification/entanglement.md",
            "Disordered" => "verification/disordered.md",
            "Identity Harness" => "verification/identity-harness.md",
        ],
        "Methods" => [
            "methods/index.md",
            "Physical" => [
                "Transfer Matrix" => "methods/transfer-matrix/index.md",
                "Bloch Hamiltonian" => "methods/bloch-hamiltonian/index.md",
                "Jordan-Wigner" => "methods/jordan-wigner/index.md",
                "Calabrese-Cardy" => "methods/calabrese-cardy/index.md",
            ],
            "Computational" => [
                "Exact Diagonalization" => "methods/exact-diagonalization/index.md",
                "Automatic Differentiation" => "methods/automatic-differentiation/index.md",
            ],
        ],
        "Derivation Notes" => [
            "JW → TFIM BdG" => "calc/jw-tfim-bdg.md",
            "Kramers-Wannier Duality" => "calc/kramers-wannier-duality.md",
            "Transfer Matrix Split" => "calc/transfer-matrix-symmetric-split.md",
            "Yang Magnetization" => "calc/yang-magnetization-toeplitz.md",
            "Heisenberg Dimer" => "calc/heisenberg-dimer-singlet-triplet.md",
            "Bethe Ansatz e₀" => "calc/bethe-ansatz-heisenberg-e0.md",
            "XXZ Luttinger Parameters" => "calc/xxz-luttinger-parameters.md",
            "Honeycomb Bloch" => "calc/bloch-honeycomb-dispersion.md",
            "Kagome Flat Band" => "calc/bloch-kagome-flat-band.md",
            "Lieb Flat Band" => "calc/bloch-lieb-flat-band.md",
            "Calabrese-Cardy OBC/PBC" => "calc/calabrese-cardy-obc-vs-pbc.md",
            "TFIM Entanglement (Peschel)" => "calc/tfim-entanglement-peschel.md",
            "AD from ln Z" => "calc/ad-thermodynamics-from-z.md",
            "Scaling Relations" => "calc/ising-scaling-relations.md",
            "Ising CFT Operators" => "calc/ising-cft-primary-operators.md",
            "Ising CFT + σ → E8" => "calc/ising-cft-magnetic-perturbation.md",
            "E8 Mass Derivation" => "calc/e8-mass-spectrum-derivation.md",
        ],
    ],
)

deploydocs(; repo="github.com/sotashimozono/QAtlas.jl.git", devbranch="main")
