"""
$(TYPEDEF)

Meteorological forcing for one time step. Read-only input to the physics: `step!` never
writes to a `MET`.

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct MET{
        Tf, MF <: AbstractMatrix{Tf}, MF64 <: AbstractMatrix{Float64},
        AF64_3 <: AbstractArray{Float64, 3},
    }

    # Domain size
    "Size of first array dimension (rows)"
    Nx::Int = 1
    "Size of second array dimension (columns)"
    Ny::Int = 1

    # Meteorological variables
    "Direct shortwave radiation per inclined surface area (W/m^2)"
    Sdir::MF = fill(NaN, Nx, Ny)
    "Diffuse shortwave radiation (W/m^2)"
    Sdif::MF = fill(NaN, Nx, Ny)
    "Direct shortwave radiation per horizontal surface area (W/m^2)"
    Sdird::MF = fill(NaN, Nx, Ny)
    "Incoming longwave radiation (W/m^2)"
    LW::MF = fill(NaN, Nx, Ny)
    "Snowfall rate (kg/m^2/s)"
    Sf::MF = fill(NaN, Nx, Ny)
    "Rainfall rate (kg/m^2/s)"
    Rf::MF = fill(NaN, Nx, Ny)
    "Total snowfall over 24h (kg/m^2)"
    Sf24h::MF = fill(NaN, Nx, Ny)
    "Air temperature (K)"
    Ta::MF = fill(NaN, Nx, Ny)
    "Relative humidity (%)"
    RH::MF = fill(NaN, Nx, Ny)
    "Wind speed (m/s)"
    Ua::MF = fill(NaN, Nx, Ny)
    "Surface air pressure (Pa)"
    Ps::MF = fill(NaN, Nx, Ny)
    "Time-varying transmissivity for direct shortwave radiation (-)"
    Tv::MF = fill(NaN, Nx, Ny)
    "Wind direction (degrees, clockwise from North) — read only by snow transport"
    Udir::MF = fill(NaN, Nx, Ny)

    # Snowfall tracking variables
    "Total snowfall over 24h (kg/m^2); Float64 to match the legacy matlab/fortran"
    Sf24h_f64::MF64 = zeros(Nx, Ny)
    "Snowfall over the last 24h (kg/m^2); Float64 to match the legacy matlab/fortran"
    Sf_history_f64::AF64_3 = zeros(Nx, Ny, 24)

end

function (::Type{MET{Tf}})(; kwargs...) where {Tf}
    return MET{Tf, Matrix{Tf}, Matrix{Float64}, Array{Float64, 3}}(; kwargs...)
end

@adapt_structure MET
