using Distributed
using ProgressMeter
using MAT
using Dates
using FSMOSHD

include("run_grid_simulation.jl")

snowtran = 1

# Configure simulations
setting = Dict(
    "tile" => "open",
    "config" => Dict("SNFRAC" => 4, "SNTRAN" => snowtran),
    "params" => Dict("wind_scaling" => 1.0)
)

if snowtran == 0
    subfolder = "SNOWTRAN_TEST_OFF"
else
    subfolder = "SNOWTRAN_TEST_ON"
end

run_grid_simulation(settings=setting, subfolder=subfolder)
