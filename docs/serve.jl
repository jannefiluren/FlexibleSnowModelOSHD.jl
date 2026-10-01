# julia --project=docs docs/serve.jl

using LiveServer

servedocs(
    include_dirs = ["examples"],
    skip_dirs = [abspath(joinpath(@__DIR__, "src", "literated"))],
)
