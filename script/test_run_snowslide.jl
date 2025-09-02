using Distributed
using ProgressMeter
using MAT
using Dates
using FSMOSHD

include("run_grid_simulation.jl")

snowslide = 0

# Configure simulations
setting = Dict(
    "tile" => "open",
    "config" => Dict("SNFRAC" => 0, "SNSLID" => snowslide),
    "params" => Dict("wind_scaling" => 0.7)
)

if snowslide == 0
    subfolder = "SNOWSLIDE_TEST_OFF"
else
    subfolder = "SNOWSLIDE_TEST_ON"
end

run_grid_simulation(settings=setting, subfolder=subfolder)
