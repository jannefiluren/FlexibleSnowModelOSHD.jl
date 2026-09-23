"""
$(TYPEDEF)

Scalar model parameters shared across all grid cells, carried by value into the physics
kernels.

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct Parameters{Tf}
    "Time step (s)"
    dt::Tf = 3600
    "Temperature measurement height (m)"
    zT::Tf = 10
    "Wind speed measurement height (m)"
    zU::Tf = 10
    "Relative humidity measurement height (m)"
    zRH::Tf = 10
    "Iterations for surface energy balance"
    Nitr::Int = 4
    "Canopy snow capacity per unit vegetation area index (kg/m^2)"
    cvai::Tf = 4.4
    "Multiplier for snowfall in forest (-)"
    pmultf_for::Tf = 0.5
    "Floor on snow layer temperature (K); -Inf disables"
    Tsnow_min::Tf = -Inf
    "Fixed snow density (kg/m^3)"
    rho0::Tf = 300
    "Temperature factor in fresh snow density (kg/m^3/K)"
    rhob::Tf = 6
    "Wind factor in fresh snow density (kg s^0.5/m^3.5)"
    rhoc::Tf = 26
    "Fresh snow density (kg/m^3)"
    rhof::Tf = 109
    "Minimum snow density (kg/m^3)"
    rhos_min::Tf = 50
    "Minimum possible snow layer thickness (m)"
    Ds_min::Tf = 0.01
    "Initial soil moisture as fraction of saturation"
    fsat::Tf = 0.5
    "Initial soil layer temperatures (K)"
    Tprof::Tf = 285
end

function Base.show(io::IO, params::Parameters{Tf}) where {Tf}
    print(io, "Parameters", '\n')
    print(io, "├── Precision: ", Tf, '\n')
    fields = propertynames(params)
    for (i, field) in enumerate(fields)
        start = i == length(fields) ? "└── " : "├── "
        print(io, start, field, ": ", getfield(params, field), '\n')
    end
    return nothing
end
