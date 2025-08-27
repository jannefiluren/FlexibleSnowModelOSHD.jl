### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 4de52076-9960-4cc2-9d6d-6bb70f48c0b0
begin
	import Pkg
    Pkg.activate(joinpath(@__DIR__,".."))
    Pkg.instantiate()
	
	using MAT
	using Dates
	using FSMOSHD
	using Plots
	using Statistics
	using PlutoUI
end

# ╔═╡ 7d557653-06ee-4c51-b158-104efa944a87
md"""
# FSM2oshd example in Julia
"""

# ╔═╡ 255f43f2-a45c-42b0-9ba5-f2b2efb0adf9
md"""

Example of manual model calibration...

Minimum albedo for melting snow (-): $(@bind asmn confirm(NumberField(0.3:0.01:0.6; default=0.5))) \
Maximum albedo for fresh snow (-): $(@bind asmx confirm(NumberField(0.6:0.01:0.95; default=0.85))) \
Melting snow albedo decay time (h): $(@bind adm confirm(NumberField(50:10:200; default=150))) \
"""

# ╔═╡ 4aeb99d5-ecf4-4306-8fdd-4ea859ad2d6e
md"""
## Results

The plot below shows snow depth for all stations in our dataset:
"""

# ╔═╡ a3d49951-6b3f-4b07-821f-e403d5fbbefa
md"""
## Code for loading input data:
"""

# ╔═╡ 57044111-0cad-4775-aa98-e32886fe2edd
begin

	
	landuse = prepare_landuse("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat")
	
	times = DateTime(2024, 9, 1, 6):Hour(1):DateTime(2025, 6, 1, 6)
	
	Sdir = zeros(length(times), length(landuse["acro"]))
	Sdif = zeros(length(times), length(landuse["acro"]))
	Sdird = zeros(length(times), length(landuse["acro"]))
	LW = zeros(length(times), length(landuse["acro"]))
	Sf = zeros(length(times), length(landuse["acro"]))
	Rf = zeros(length(times), length(landuse["acro"]))
	Ta = zeros(length(times), length(landuse["acro"]))
	RH = zeros(length(times), length(landuse["acro"]))
	Ua = zeros(length(times), length(landuse["acro"]))
	Ps = zeros(length(times), length(landuse["acro"]))

	meteo_folder = "K:/DATA_ICON/OUTPUT_OSHD_STAT/PROCESSED_ANALYSIS/ICON_1EFA"
	
	for (i, t) in enumerate(times)
		folder = joinpath(meteo_folder, Dates.format(t, "yyyy.mm"))
		filename = "ICONDATA_" * Dates.format(t, "yyyymmddHHMM") * "_C1EFA_" 
		file = searchdir(folder, filename)
	
		met_single = matread(joinpath(folder, file[1]))
	
		Sdir[i, :] = met_single["sdri"]["data"]
		Sdif[i, :] = met_single["sdfd"]["data"]
		Sdird[i, :] = met_single["sdrd"]["data"]
		LW[i, :] = met_single["lwrs"]["data"]
		Rf[i, :] = met_single["prfc"]["data"]
		Sf[i, :] = met_single["psfc"]["data"]
		Ta[i, :] = met_single["tais"]["data"]
		RH[i, :] = met_single["rhus"]["data"]
		Ua[i, :] = met_single["wnss"]["data"]
		Ps[i, :] = met_single["pail"]["data"]
	end
	
end

# ╔═╡ 0a778b48-a6ad-4d0b-9238-29273263f8c5
md"""
## Code for setting up the model:
"""

# ╔═╡ dffe499f-feec-43ba-83a5-ab9007afd8a1
begin
	settings = Dict("tile" => "open", "params" => Dict("wind_scaling" => 0.7))
	fsm = setup(Float32, Int32, landuse, length(landuse["acro"]), 1, settings)
	fsm.asmn = asmn
	fsm.afs .= asmx
	fsm.adm = adm
	
	met_curr = MET{Float32, Int32}(Nx=length(landuse["acro"]))

	nothing
end

# ╔═╡ 2ada6b38-53fa-4a96-9b01-50e47bfb3975
md"""
## Code for running the model:
"""

# ╔═╡ f48ec7b8-09b2-436f-80b0-60ecc060847a
begin

	snowdepth = zeros(length(times), length(landuse["acro"]))

	time_model = @elapsed for (i, t) in enumerate(times)
		
		  met_curr.Sdir[:, :] .= Sdir[i, :]
		  met_curr.Sdif[:, :] .= Sdif[i, :]
		  met_curr.Sdird[:, :] = Sdird[i, :]
		  met_curr.LW[:, :] .= LW[i, :]
		  met_curr.Sf[:, :] .= Sf[i, :]
		  met_curr.Rf[:, :] .= Rf[i, :]
		  met_curr.Ta[:, :] .= Ta[i, :]
		  met_curr.RH[:, :] .= RH[i, :]
		  met_curr.Ua[:, :] .= Ua[i, :]
		  met_curr.Ps[:, :] .= Ps[i, :]
		
		  met_curr.Sf24h[:, :] .= sum(Sf[max(i-23,1):i,:], dims=1)'
		  
		  step!(fsm, met_curr, t)
		
		  snowdepth[i,:] = dropdims(sum(fsm.Ds, dims=1), dims=3)

	end

	nothing

end

# ╔═╡ c6b001ee-40d6-4c8d-911b-e6577e94e70f
md"""
Model run time for $(length(landuse["acro"])) stations and $(round(length(times)/24)) days: $(round(time_model)) seconds
"""

# ╔═╡ 573d40a9-5edb-4931-8cee-cccd53703869
begin
	plot(times, mean(snowdepth, dims=2), title="Average of all stations")
	ylabel!("Snow depth (m)")
end

# ╔═╡ a29e9044-0375-496b-9829-41a239f8427f
md"""
## Code for importing packages:
"""

# ╔═╡ Cell order:
# ╟─7d557653-06ee-4c51-b158-104efa944a87
# ╟─255f43f2-a45c-42b0-9ba5-f2b2efb0adf9
# ╟─c6b001ee-40d6-4c8d-911b-e6577e94e70f
# ╟─4aeb99d5-ecf4-4306-8fdd-4ea859ad2d6e
# ╟─573d40a9-5edb-4931-8cee-cccd53703869
# ╟─a3d49951-6b3f-4b07-821f-e403d5fbbefa
# ╠═57044111-0cad-4775-aa98-e32886fe2edd
# ╟─0a778b48-a6ad-4d0b-9238-29273263f8c5
# ╠═dffe499f-feec-43ba-83a5-ab9007afd8a1
# ╟─2ada6b38-53fa-4a96-9b01-50e47bfb3975
# ╠═f48ec7b8-09b2-436f-80b0-60ecc060847a
# ╟─a29e9044-0375-496b-9829-41a239f8427f
# ╠═4de52076-9960-4cc2-9d6d-6bb70f48c0b0
