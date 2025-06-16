using FSMOSHD
using BenchmarkTools
using Dates
using MAT

landuse_file_loc = "D:/jim_operational/SOURCE/BAFU_LUS_0250_2023a.mat"

landuse_file = matopen(landuse_file_loc)
landuse = read(landuse_file, "landuse")
close(landuse_file)

Nx = round(Int, landuse["nrows"])
Ny = round(Int, landuse["ncols"])

fsm = FSM{Float64, Int64}(Nx=Nx, Ny=Ny)
setup_grid!(fsm, landuse)

meteo = MET{Float64, Int64}(Nx=Nx, Ny=Ny)

t = DateTime(2022,09,01,06,00,00)

drive_grid!(meteo, fsm, t)

# Radiation

print("Radiation")

@benchmark radiation(fsm, meteo, t)
@benchmark radiation_par(fsm, meteo, t)

# Thermal

print("Thermal")

@benchmark thermal(fsm)
@benchmark thermal_par(fsm)

# Sfexch

print("Sfexch")

@benchmark sfexch(fsm, meteo)
@benchmark sfexch_par(fsm, meteo)

# Ebalsrf

print("Ebalsrf")

@benchmark ebalsrf(fsm, meteo)
@benchmark ebalsrf_par(fsm, meteo)

# Snow

print("Snow")

@benchmark snow(fsm, meteo)
@benchmark snow_par(fsm, meteo)

# Soil

print("Soil")

@benchmark soil(fsm)
@benchmark soil_par(fsm)
