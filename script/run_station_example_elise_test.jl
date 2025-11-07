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
    lus["dem"] = Dict("data" => [2540.0;;])          # [2057.0;;] au Lautaret , [2540.0;;] SLF
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
    laut_meteo = CSV.read("../lautaret/2017-2024_Lautaret_halfhour_input_corr.csv", DataFrame)

    return fsm, met, df_meteo, laut_meteo

end


function run_fsm(fsm, met, df_meteo)

    # allocate output variable-wise
    hs = zeros(nrow(df_meteo))
    Tsnow1 = fill(NaN, nrow(df_meteo))
    Tsnow2 = fill(NaN, nrow(df_meteo))
    Tsnow3 = fill(NaN, nrow(df_meteo))
    Tsrf = zeros(nrow(df_meteo))
    alb = zeros(nrow(df_meteo))
    Sice = zeros(nrow(df_meteo))
    Sliq = zeros(nrow(df_meteo))
    snowdepthmin = zeros(nrow(df_meteo))
    snowdepthmax = zeros(nrow(df_meteo))
    swemin = zeros(nrow(df_meteo))
    swemax = zeros(nrow(df_meteo))
    
    
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
        println(hs)
        if fsm.Nsnow[1, 1] == 1
            Tsnow1[i] = fsm.Tsnow[1, 1, 1]
        elseif fsm.Nsnow[1, 1] == 2
            Tsnow1[i], Tsnow2[i] = fsm.Tsnow[1, 1, 1], fsm.Tsnow[2, 1, 1]
        elseif fsm.Nsnow[1, 1] == 3
            Tsnow1[i], Tsnow2[i], Tsnow3[i] = fsm.Tsnow[1, 1, 1], fsm.Tsnow[2, 1, 1], fsm.Tsnow[3, 1, 1]
        end
        alb[i] = fsm.asrf_out[1, 1] # albs
        Tsrf[i] = fsm.Tsrf[1, 1]
        Sice[i] = dropdims(sum(fsm.Sice, dims=1), dims=1)[1, 1]
        Sliq[i] = dropdims(sum(fsm.Sliq, dims=1), dims=1)[1, 1]
        snowdepthmin[i] = fsm.snowdepthmin[1, 1]
        snowdepthmax[i] = fsm.snowdepthmax[1, 1]
        swemin[i] = fsm.swemin[1, 1]
        swemax[i] = fsm.swemax[1, 1]
    
    end

    # write results to dataframe
    time = DateTime.(df_meteo[!, "year"], df_meteo[!, "month"], df_meteo[!, "day"], df_meteo[!, "hour"])
    df_results = DataFrame(time=time, hs=hs, Tsnow1=Tsnow1, Tsnow2=Tsnow2, Tsnow3=Tsnow3, Ts=Tsrf, albedo=alb, I=Sice, W=Sliq, snow_depth_min=snowdepthmin, snow_depth_max=snowdepthmax, swemin=swemin, swemax=swemax)

    return df_results

end

#####################################################################################
function run_fsm_laut(fsm, met, laut_meteo)

    # allocate output variable-wise
    hs = zeros(nrow(laut_meteo))
    Tsnow1 = fill(NaN, nrow(laut_meteo))
    Tsnow2 = fill(NaN, nrow(laut_meteo))
    Tsnow3 = fill(NaN, nrow(laut_meteo))
    Tsrf = zeros(nrow(laut_meteo))
    alb = zeros(nrow(laut_meteo))
    Sice = zeros(nrow(laut_meteo))
    Sliq = zeros(nrow(laut_meteo))
    snowdepthmin = zeros(nrow(laut_meteo))
    snowdepthmax = zeros(nrow(laut_meteo))
    swemin = zeros(nrow(laut_meteo))
    swemax = zeros(nrow(laut_meteo))
    
    # time loop
    for (i, row) in zip(1:nrow(laut_meteo), eachrow(laut_meteo))
    
        # assign input
        met.year .= year(Date.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00"))
        met.month .= month(Date.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00"))
        met.day .= day(Date.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00"))
        met.hour .= hour(DateTime.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00")) .+ minute(DateTime.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00"))./60
        met.Sdir .= row["Null"]
        met.Sdif .= row["short_up_Avg"]
        met.Sdird .= row["Null"]
        met.LW .= row["long_up_cor_Avg"]
        met.Sf .= row["Sf"]
        met.Rf .= row["Rf"]
        met.Ta .= row["AirTC_Avg"]
        met.RH .= row["HRair_Avg"]
        met.Ua .= row["WindSpeed_Avg"] #"WindSpeed_Low_Avg"
        met.Ps .= row["Patm_Avg"]
        met.Sf24h .= row["Sf24h"]
    
        # set time 
        # t = DateTime.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00")
        t = DateTime(year(Date.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00")), 
        month(Date.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00")), 
        day(Date.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00")), 
        hour(DateTime.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00")) .+ minute(DateTime.(row["DateTime"], "yyyy-mm-dd HH:MM:S+00:00"))./60)
    
        # run model and update states
        step!(fsm, met, t)
    
        # write output
        hs[i] = dropdims(sum(fsm.Ds, dims=1), dims=1)[1, 1]
        if fsm.Nsnow[1, 1] == 1
            Tsnow1[i] = fsm.Tsnow[1, 1, 1]
        elseif fsm.Nsnow[1, 1] == 2
            Tsnow1[i], Tsnow2[i] = fsm.Tsnow[1, 1, 1], fsm.Tsnow[2, 1, 1]
        elseif fsm.Nsnow[1, 1] == 3
            Tsnow1[i], Tsnow2[i], Tsnow3[i] = fsm.Tsnow[1, 1, 1], fsm.Tsnow[2, 1, 1], fsm.Tsnow[3, 1, 1]
        end
        alb[i] = fsm.asrf_out[1, 1] # albs
        Tsrf[i] = fsm.Tsrf[1, 1]
        Sice[i] = dropdims(sum(fsm.Sice, dims=1), dims=1)[1, 1]
        Sliq[i] = dropdims(sum(fsm.Sliq, dims=1), dims=1)[1, 1]
        snowdepthmin[i] = fsm.snowdepthmin[1, 1]
        snowdepthmax[i] = fsm.snowdepthmax[1, 1]
        swemin[i] = fsm.swemin[1, 1]
        swemax[i] = fsm.swemax[1, 1]
        
    end

    alb_snow = copy(alb) 
    alb_snow[alb .< 0.6] .= NaN

    # write results to dataframe
    time = DateTime.(laut_meteo[!, "DateTime"], "yyyy-mm-dd HH:MM:S+00:00")
    df_results = DataFrame(time=time, hs=hs, Tsnow1=Tsnow1, Tsnow2=Tsnow2, Tsnow3=Tsnow3, Ts=Tsrf, albedo=alb_snow, I=Sice, W=Sliq, snow_depth_min=snowdepthmin, snow_depth_max=snowdepthmax, swemin=swemin, swemax=swemax)

    return df_results

end
#####################################################################################

fsm, met, df_meteo, laut_meteo = setup_example()
# laut_meteo = coalesce.(laut_meteo, 0)
println(laut_meteo)

df_results = run_fsm(fsm, met, df_meteo)
# CSV.write("../data/output_SLF_5WJ_corr-dt.txt", df_results)

# df_results_laut = run_fsm_laut(fsm, met, laut_meteo)
# # println(met)
# # println(fsm)
# CSV.write("../lautaret/output_2017-2024_lautaret_halfhour_corr.txt", df_results_laut)
