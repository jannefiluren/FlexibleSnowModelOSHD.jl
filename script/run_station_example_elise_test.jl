cd(@__DIR__)

using Dates
using CSV
using DataFrames
using FSMOSHD


# initialization function
    # sets parameters (topo, meteo), configuration, variables (format)
    # opens meteo file
function setup_example()

    # set landuse properties
    lus = Dict()
    lus["skyvf"] = Dict("data" => [1.0;;])           
    lus["x"] = Dict("data" => [1.0;;])
    lus["y"] = Dict("data" => [1.0;;])
    lus["dem"] = Dict("data" => [2540.0;;])
    lus["slopemu"] = Dict("data" => [1.0;;])
    lus["xi"] = Dict("data" => [1.0;;])
    lus["Ld"] = Dict("data" => [1.0;;])
    lus["prec_multi"] = Dict("data" => [1.0;;])
    
    # define custom settings
    settings = Dict("tile" => "open", "params" => Dict("wind_scaling" => 0.7))
    
    # create fsm struct
    fsm = setup(Float32, Int32, lus, 1, 1, settings)
    
    # define meteo data struct
    met = MET{Float32,Int32}()
    
    # read meteo file
    df_meteo = CSV.read("../data/input_SLF_5WJ.txt", DataFrame)
    laut_meteo = CSV.read("../2223_Lautaret_halfhour_input.csv", DataFrame)

    return fsm, met, df_meteo, laut_meteo

end


function run_fsm(fsm, met, df_meteo)

    # allocate output variable-wise
    hs = zeros(nrow(df_meteo))
    T_snow = zeros(nrow(df_meteo))
    density = zeros(nrow(df_meteo))
    alb = zeros(nrow(df_meteo))
    
    # time loop
    for (i, row) in zip(1:nrow(df_meteo), eachrow(df_meteo))
    
        # assign input
     
        met.year .= row["year"]
        met.month .= row["month"]
        met.day .= row["day"]
        met.hour .= row["hour"]
        met.Sdir .= row["Sdir"]
        met.Sdif .= row["Sdif"]
        met.Sdird .= row["Sdir"]
        met.LW .= row["LW"]
        met.Sf .= row["Sf"]
        met.Rf .= row["Rf"]
        met.Ta .= row["Ta"]
        met.RH .= row["RH"]
        met.Ua .= row["Ua"]
        met.Ps .= row["Ps"]
        met.Sf24h .= row["Sf24h"]
    
        # set time 
        t = DateTime(row["year"], row["month"], row["day"], row["hour"])
    
        # run model and update states
        step!(fsm, met, t)
    
        # write output
        hs[i] = dropdims(sum(fsm.Ds, dims=1), dims=1)[1, 1]
        T_snow[i] = dropdims(sum(fsm.Tsnow, dims=1), dims=1)[1, 1]
        alb[i] = dropdims(sum(fsm.albs, dims=1), dims=1)[1, 1]
    
    end

    # write results to dataframe
    time = DateTime.(df_meteo[!, "year"], df_meteo[!, "month"], df_meteo[!, "day"], df_meteo[!, "hour"])
    df_results = DataFrame(time=time, hs=hs, Tsnow=T_snow, albedo=alb)

    return df_results

end

#####################################################################################
function run_fsm_lautaret(fsm, met, laut_meteo)

    # allocate output variable-wise
    hs = zeros(nrow(laut_meteo))
    T_snow = zeros(nrow(laut_meteo))
    density = zeros(nrow(laut_meteo))
    alb = zeros(nrow(laut_meteo))
    
    # time loop
    for (i, row) in zip(1:nrow(laut_meteo), eachrow(laut_meteo))
    
        # assign input
        met.datetime .= row[""]
        met.Sdir .= row["Rs_net_Avg"]
        met.Sdif .= row["Null"]
        met.Sdird .= row["Null"]
        met.LW .= row["Rl_net_Avg"]
        met.Sf .= row["Sf"]
        met.Rf .= row["Rf"]
        met.Ta .= row["AirTC_Avg"]
        met.RH .= row["HRair_Avg"]
        met.Ua .= row["WindSpeed_Low_Avg"]
        met.Ps .= row["Patm_Avg"]
        met.Sf24h .= row["Null"]
    
        # set time 
        t = met.datetime
    
        # run model and update states
        step!(fsm, met, t)
    
        # write output
        hs[i] = dropdims(sum(fsm.Ds, dims=1), dims=1)[1, 1]
        T_snow[i] = dropdims(sum(fsm.Tsnow, dims=1), dims=1)[1, 1]
        alb[i] = dropdims(sum(fsm.albs, dims=1), dims=1)[1, 1]
    
    end

    # write results to dataframe
    time = t
    df_results = DataFrame(time=time, hs=hs, Tsnow=T_snow, albedo=alb)

    return df_results

end
#####################################################################################

fsm, met, df_meteo, laut_meteo = setup_example()

println(laut_meteo)

df_results = run_fsm(fsm, met, df_meteo)

CSV.write("../data/output_SLF_5WJ.txt", df_results)
df_results_laut = run_fsm_laut(fsm, met, laut_meteo)

CSV.write("../output_lautaret.txt", df_results_laut)
