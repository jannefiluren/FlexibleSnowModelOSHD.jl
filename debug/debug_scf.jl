using FSMOSHD
using Dates
using MAT
using ProgressMeter
using WGLMakie

times = DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 1, 6)

irow = 489
icol = 883

# Helper functions

searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))

# Prepare landuse

landuse = prepare_landuse("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat")

fields = ["skyvf", "dem", "Ld", "slopemu", "xi"]
for field in fields
  landuse[field]["data"] = fill(landuse[field]["data"][irow, icol],1,1)
end

# Load meteo data

file = joinpath("D:/julia/debug_scf", "meteo_" * string(irow) * "_" * string(icol) * ".mat") 

if isfile(file)

  meteo_all = matread(file)

else

  meteo_all = Dict()

  meteo_all["Sdir"] = zeros(length(times))
  meteo_all["Sdif"] = zeros(length(times))
  meteo_all["Sdird"] = zeros(length(times))
  meteo_all["LW"] = zeros(length(times))
  meteo_all["Rf"] = zeros(length(times))
  meteo_all["Sf"] = zeros(length(times))
  meteo_all["Ta"] = zeros(length(times))
  meteo_all["RH"] = zeros(length(times))
  meteo_all["Ua"] = zeros(length(times))
  meteo_all["Ps"] = zeros(length(times))

  @showprogress "Loading meteo data..." for (i, t) in enumerate(times)

    # Prepare forcing data
    folder_icon = joinpath("K:/DATA_ICON/OUTPUT_OSHD_0250/PROCESSED_ANALYSIS/ICON_1EFA", Dates.format(t, "yyyy.mm"))
    filename_icon = searchdir(folder_icon, "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_")
    
    local file_mat
    file_mat = matopen(joinpath(folder_icon, filename_icon[1]))
    try
      meteo_all["Sdir"][i] = read(file_mat, "sdri")["data"][irow, icol]     # "direct shortwave radiation, per inclined surface area, within topography, above canopy"
      meteo_all["Sdif"][i] = read(file_mat, "sdfd")["data"][irow, icol]     # "diffuse shortwave radiation, per horizontal surface area, within topography, above canopy"
      meteo_all["Sdird"][i] = read(file_mat, "sdrd")["data"][irow, icol]    # "direct shortwave radiation, per horizontal surface area, within topography, above canopy"
      meteo_all["LW"][i] = read(file_mat, "lwrs")["data"][irow, icol]       # "longwave radiation, above topography"
      meteo_all["Rf"][i] = read(file_mat, "prfc")["data"][irow, icol]       # "rainfall"
      meteo_all["Sf"][i] = read(file_mat, "psfc")["data"][irow, icol]       # "snowfall (snow + graupel)"
      meteo_all["Ta"][i] = read(file_mat, "tais")["data"][irow, icol]       # "air temperature"
      meteo_all["RH"][i] = read(file_mat, "rhus")["data"][irow, icol]       # "relative humidity"
      meteo_all["Ua"][i] = read(file_mat, "wnss")["data"][irow, icol]       # "wind speed"
      meteo_all["Ps"][i] = read(file_mat, "pail")["data"][irow, icol]       # "local air pressure"
    finally
      close(file_mat)
    end
    
  end

  matwrite(file, meteo_all; compress = true)

end


# Setup and run model

begin

  # Clear debug dict

  debug_dict["snowdepth"] = Float64[]
  debug_dict["swe"] = Float64[]
  debug_dict["swemax"] = Float64[]
  debug_dict["swemin"] = Float64[]
  debug_dict["snowdepthmax"] = Float64[]
  debug_dict["snowdepthmin"] = Float64[]
  debug_dict["fsnow_season"] = Float64[]
  debug_dict["fsnow_nsnow"] = Float64[]
  debug_dict["fsnow_nsnow_recent"] = Float64[]
  debug_dict["fsnow"] = Float64[]
  debug_dict["dsnowdepth"] = Float64[]
  debug_dict["dsnowdepthmax"] = Float64[]
  debug_dict["Ds"] = Float64[]
  debug_dict["rhos"] = Float64[]

  # Setup and run model

  Nx = size(landuse["dem"]["data"], 1)
  Ny = size(landuse["dem"]["data"], 2)

  settings = Dict(
              "tile" => "open",
              "config" => Dict("SNFRAC" => 0)
            )

  fsm = setup(Float32, Int32, landuse, Nx, Ny, settings)

  met_curr = MET{Float32, Int32}(Nx=Nx, Ny=Ny)

  Sf24h = zeros(size(met_curr.Sf24h))
  Sf_history = zeros(size(met_curr.Sf_history))

  for (i, t) in enumerate(times)

    met_curr.Sdir[:, :] .= meteo_all["Sdir"][i]
    met_curr.Sdif[:, :] .= meteo_all["Sdif"][i]
    met_curr.Sdird[:, :] .= meteo_all["Sdird"][i]
    met_curr.LW[:, :] .= meteo_all["LW"][i]
    met_curr.Sf[:, :] .= meteo_all["Sf"][i]
    met_curr.Rf[:, :] .= meteo_all["Rf"][i]
    met_curr.Ta[:, :] .= meteo_all["Ta"][i]
    met_curr.RH[:, :] .= meteo_all["RH"][i]
    met_curr.Ua[:, :] .= meteo_all["Ua"][i]
    met_curr.Ps[:, :] .= meteo_all["Ps"][i]

    curr_hour = Dates.value(Hour(t)) + 1
    Sf24h .+= meteo_all["Sf"][i]
    Sf24h .-= Sf_history[:, :, curr_hour]
    Sf_history[:, :, curr_hour] .= meteo_all["Sf"][i]
    met_curr.Sf24h[:, :] .= Sf24h

    step!(fsm, met_curr, t)
    
  end

end

# Plot results

f = Figure(size = (1000, 1500))

ax1 = Axis(f[1, 1])
lines!(ax1, times, debug_dict["snowdepth"], linestyle = :solid, label = "snowdepth")
lines!(ax1, times, debug_dict["snowdepthmax"], linestyle = :dash, label = "snowdepthmax")
lines!(ax1, times, debug_dict["snowdepthmin"], linestyle = :dot, label = "snowdepthmin")
lines!(ax1, times, debug_dict["Ds"], linestyle = :solid, label = "Ds")
f[1, 2] = Legend(f, ax1, "Snow Depth", framevisible = false)

ax2 = Axis(f[2, 1])
lines!(ax2, times, debug_dict["swe"], linestyle = :solid, label = "swe")
lines!(ax2, times, debug_dict["swemin"], linestyle = :dash, label = "swemin")
lines!(ax2, times, debug_dict["swemax"], linestyle = :dot, label = "swemax")
f[2, 2] = Legend(f, ax2, "SWE", framevisible = false)

ax3 = Axis(f[3, 1]) 
lines!(ax3, times, debug_dict["fsnow_season"], linestyle = :solid, label = "fsnow_season")
lines!(ax3, times, debug_dict["fsnow_nsnow"], linestyle = :dash, label = "fsnow_nsnow")
lines!(ax3, times, debug_dict["fsnow_nsnow_recent"], linestyle = :dashdot, label = "fsnow_nsnow_recent")
lines!(ax3, times, debug_dict["fsnow"], linestyle = :dot, label = "fsnow")
f[3, 2] = Legend(f, ax3, "Fractional Snow Cover", framevisible = false)

ax4 = Axis(f[4, 1])
lines!(ax4, times, debug_dict["rhos"], linestyle = :solid, label = "rhos")
f[4, 2] = Legend(f, ax4, "Misc", framevisible = false)

linkxaxes!(ax1, ax2, ax3, ax4)

f
