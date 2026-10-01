using Documenter
using Literate
using FlexibleSnowModelOSHD

# Make `using FlexibleSnowModelOSHD` implicit in every jldoctest block.
DocMeta.setdocmeta!(
    FlexibleSnowModelOSHD,
    :DocTestSetup,
    :(using FlexibleSnowModelOSHD);
    recursive = true,
)

# Render the runnable scripts in examples/ to executed Markdown pages. Run each
# example (execute = true) during the build, so a broken script fails the docs build.
const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const LITERATED_DIR = joinpath(@__DIR__, "src", "literated")
example_scripts = [
    "constructing_a_model.jl",
    "run_open_station_example.jl",
    "run_forest_station_example.jl",
]
for script in example_scripts
    Literate.markdown(
        joinpath(EXAMPLES_DIR, script), LITERATED_DIR;
        flavor = Literate.DocumenterFlavor(), execute = true,
    )
end

makedocs(
    sitename = "FlexibleSnowModelOSHD.jl",
    modules = [FlexibleSnowModelOSHD],
    authors = "jannefiluren <jan.magnusson@slf.ch>",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://jannefiluren.github.io/FlexibleSnowModelOSHD.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Constructing a model" => "literated/constructing_a_model.md",
            "Open-site simulation" => "literated/run_open_station_example.md",
            "Forest-site simulation" => "literated/run_forest_station_example.md",
        ],
        "API reference" => "api.md",
    ],
    doctest = true,
    # Most of the package is not documented yet: tolerate undocumented bindings and
    # @ref links into not-yet-documented symbols. Doctests stay strict (not listed here).
    # Tighten to `warnonly = false` as documentation coverage grows.
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/jannefiluren/FlexibleSnowModelOSHD.jl",
    devbranch = "jan-gpu-on-main",   # TODO: revert to "main" after merge
    push_preview = true,
)
