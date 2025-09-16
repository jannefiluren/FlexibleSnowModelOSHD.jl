cd(@__DIR__)

using Dates
using CSV
using DataFrames
using FSMOSHD

# set landuse properties

lus = Dict()
lus["skyvf"] = Dict("data" => [1.0;;])
lus["x"] = Dict("data" => [1.0;;])
lus["y"] = Dict("data" => [1.0;;])
lus["dem"] = Dict("data" => [2540.0;;])
lus["slopemu"] = Dict("data" => [1.0;;])
lus["xi"] = Dict("data" => [1.0;;])
lus["Ld"] = Dict("data" => [1.0;;])

# define custom settings

settings = Dict("tile" => "open", "params" => Dict("wind_scaling" => 0.7))

# create fsm struct

fsm = setup(Float32, Int32, lus, 1, 1, settings)

# define meteo data struct

met = MET{Float32,Int32}()

# read meteo file

header = ["year", "month", "day", "hour", "Sdir", "Sdif", "LW", "Sf", "Rf", "Ta", "RH", "Ua", "Ps", "Sf24h", "Tv"]

data = CSV.read("../data/input_SLF_5WJ.txt", DataFrame, header=header)

# allocate output variable-wise
hs = zeros(nrow(data))

# time loop 

for (i,row) in zip(1:nrow(data), eachrow(data))
    
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
    # met.Udir
    met.Ps .= row["Ps"]
    met.Sf24h .= row["Sf24h"]
    # met.Tc
    # met.es
    # met.Qa
    # met.Tv
    # met.Sf_history

    # set time 
    t = DateTime(row["year"], row["month"], row["day"], row["hour"])

    # run model & update states (fsm struct)
    step!(fsm, met, t)

    # write output 
    hs[i] = dropdims(sum(fsm.Ds,dims=1),dims=1)[1,1]

end

time = DateTime.(data[!,"year"], data[!,"month"], data[!,"day"], data[!,"hour"])
 
data_out = DataFrame(time = time, hs = hs)
 
CSV.write("../data/output.csv", data_out)
 

