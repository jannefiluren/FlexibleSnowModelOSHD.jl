# # Forest-site simulation
#
# A point simulation over a single forested cell. Same driving data as the open-site example,
# but with canopy properties and the `ForestCover` land cover, so canopy interception,
# unloading and sub-canopy radiation are active.

using Dates
using CSV
using DataFrames
using FlexibleSnowModelOSHD

# ## Model setup
#
# The landuse now carries canopy fields (`fveg`, `hcan`, `lai`, …) and the FSM is built with
# `ForestCover`.

function setup_example()

    ## Set landuse properties
    lus = Dict()
    lus["skyvf"] = Dict("data" => [1.0;;])
    lus["elevation"] = Dict("data" => [2540.0;;])
    lus["slopemu"] = Dict("data" => [1.0;;])
    lus["xi"] = Dict("data" => [1.0;;])
    lus["Ld"] = Dict("data" => [1.0;;])
    lus["prec_multi"] = Dict("data" => [1.0;;])

    ## Add forest properties
    lus["forest"] = Dict("data" => [1;;]) # Forest cover fraction
    lus["fveg"] = Dict("data" => [0.6;;]) # Canopy cover fraction
    lus["fves"] = Dict("data" => [0.6;;]) # Stand-scale canopy cover fraction
    lus["hcan"] = Dict("data" => [20;;])  # Canopy height (m)
    lus["lai"] = Dict("data" => [2.5;;])  # Leaf area index
    lus["vfhp"] = Dict("data" => [0.5;;]) # Hemispherical sky-view fraction including canopy

    ## Create FSM struct with forest terrain land cover (ForestCover / SoilSubstrate schemes)
    grid = Grid(Float32; Nx = 1, Ny = 1)
    fsm = FSM(grid, lus; land_cover = ForestCover{Float32}())

    ## Define meteo data struct
    met = MET{Float32}()

    ## Read meteorological data from a text file
    df_meteo = CSV.read(joinpath(pkgdir(FlexibleSnowModelOSHD), "data", "input_SLF_5WJ.txt"), DataFrame)

    return fsm, met, df_meteo

end

# ## Time loop
#
# As in the open-site example, plus the canopy transmissivity `Tv` forcing.

function run_fsm(fsm, met, df_meteo)

    ## Allocate output variable-wise
    hs = zeros(nrow(df_meteo))

    ## Loop over time
    for (i, row) in zip(1:nrow(df_meteo), eachrow(df_meteo))

        ## Assign input to the MET struct
        met.Sdir .= row["Sdir"]
        met.Sdif .= row["Sdif"]
        met.Sdird .= row["Sdir"]
        met.LW .= row["LW"]
        met.Sf .= row["Sf"] / fsm.params.dt  # Convert accumulation (kg/m^2) to rate (kg/m^2/s)
        met.Rf .= row["Rf"] / fsm.params.dt  # Convert accumulation (kg/m^2) to rate (kg/m^2/s)
        met.Ta .= row["Ta"]
        met.RH .= row["RH"]
        met.Ua .= row["Ua"]
        met.Ps .= row["Ps"]
        met.Sf24h .= row["Sf24h"]
        met.Tv .= row["Tv"]

        ## Set time
        t = DateTime(row["year"], row["month"], row["day"], row["hour"])

        ## Run model and update states
        step!(fsm, met, t)

        ## Write output
        hs[i] = dropdims(sum(fsm.state.Ds, dims = 1), dims = 1)[1, 1]

    end

    ## Write results to a dataframe
    time = DateTime.(df_meteo[!, "year"], df_meteo[!, "month"], df_meteo[!, "day"], df_meteo[!, "hour"])
    df_results = DataFrame(time = time, hs = hs)

    return df_results

end

# ## Run the example
#
# Set up, run the season, and summarise the resulting snow-depth series.

fsm, met, df_meteo = setup_example()

df_results = run_fsm(fsm, met, df_meteo)

describe(df_results)
