"""
$(TYPEDEF)

Snow-transport workspace: everything SnowSlide and SnowTran3D need, kept off the model so
the decomposed, GPU-adaptable [`FSM`](@ref) stays untouched. Holds the static per-cell setup
arrays, the cumulative SWE-change accumulators, the tuning constants (defaults match
`deps/MODULES.F90`, so the Julia and Fortran paths agree) and the SnowTran3D working arrays.
Transport is neighbour-coupled and CPU-only, so all arrays are host `Matrix`es of size
`(Nx, Ny)`. Build one with [`setup_transport`](@ref).

# Fields

$(TYPEDFIELDS)
"""
@kwdef mutable struct SnowTransport{Tf}

    "Size of the first array dimension (rows)"
    Nx::Int = 1
    "Size of the second array dimension (columns)"
    Ny::Int = 1

    "Run wind transport (SnowTran3D) when `transport!` is called from `step!`"
    wind::Bool = false
    "Run snow slides (SnowSlide) when `transport!` is called from `step!`"
    slide::Bool = false
    "Use the Fortran ccall path instead of the Julia port"
    use_fortran::Bool = false

    "Minimum snow density (kg/m^3)"
    rhos_min::Tf = 50
    "Maximum snow density (kg/m^3)"
    rhos_max::Tf = 750
    "Constant density of wind-transported snow (kg/m^3)"
    rho_snow::Tf = 300
    "Dynamic snow holding depth ratio (-)"
    dyn_ratio::Tf = 0.7
    "Hysteretic triggering snow holding depth ratio (-)"
    trig_ratio::Tf = 1.5
    "Minimum slope for snow slide occurrence (deg)"
    slope_min::Tf = 25
    "Minimum snow holding depth (m)"
    Shd_min::Tf = 0.05

    "Slope angle (deg)"
    slope::Matrix{Tf} = fill(Tf(NaN), Nx, Ny)
    "Snow holding depth for gravitational transport (m)"
    Shd::Matrix{Tf} = fill(Tf(NaN), Nx, Ny)
    "Vegetation snow holding depth (m)"
    vegsnowd_xy::Matrix{Tf} = fill(Tf(0.1), Nx, Ny)
    "Open-tile fraction, `1 - forest - glacier` (-)"
    tilefrac::Matrix{Tf} = ones(Tf, Nx, Ny)
    "Glacier fraction (-)"
    glacierfrac::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Open fraction for transport, `tilefrac + glacierfrac` (-)"
    openfrac::Matrix{Tf} = ones(Tf, Nx, Ny)
    "DEM `(i, j)` indices sorted by decreasing elevation, `(Nx*Ny, 2)`"
    index_sorted_dem::Matrix{Int} = zeros(Int, Nx * Ny, 2)

    "Cumulative SWE change due to slides (kg/m^2)"
    dSWE_tot_slide::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Cumulative SWE change due to saltation (kg/m^2)"
    dSWE_tot_salt::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Cumulative SWE change due to suspension (kg/m^2)"
    dSWE_tot_susp::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Cumulative SWE change due to sublimation (kg/m^2)"
    dSWE_tot_subl::Matrix{Tf} = zeros(Tf, Nx, Ny)

    "SWE change due to saltation in the current step (kg/m^2)"
    dSWE_salt::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "SWE change due to suspension in the current step (kg/m^2)"
    dSWE_susp::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "SWE change due to sublimation in the current step (kg/m^2)"
    dSWE_subl::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "SWE change due to slides in the current step (kg/m^2)"
    dSWE_slide::Matrix{Tf} = zeros(Tf, Nx, Ny)

    "SnowSlide: cells receiving slide deposits"
    snow_depo::Matrix{Bool} = zeros(Bool, Nx, Ny)
    "SnowSlide: snow holding depth floored at `Shd_min` (m)"
    Shd_corr::Matrix{Tf} = zeros(Tf, Nx, Ny)

    "W->E (y) wind component (m/s)"
    uwind::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "S->N (x) wind component (m/s)"
    vwind::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Snowpack thickness, not scaled by fsnow (m)"
    snowthickness::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Thickness of the erodible soft snow (m)"
    Ds_soft::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Friction velocity (m/s)"
    Utau::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Threshold friction velocity (m/s)"
    Utau_t::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Height of the saltation layer (m)"
    h_star::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Surface roughness length (m)"
    z_0::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "W->E (y) saltation flux (kg/m/s)"
    Qsalt_u::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "S->N (x) saltation flux (kg/m/s)"
    Qsalt_v::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "W->E (y) maximum saltation flux (kg/m/s)"
    Qsalt_maxu::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "S->N (x) maximum saltation flux (kg/m/s)"
    Qsalt_maxv::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "W->E (y) suspension flux (kg/m/s)"
    Qsusp_u::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "S->N (x) suspension flux (kg/m/s)"
    Qsusp_v::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "W->E (y) equilibrium suspension flux (kg/m/s)"
    Qsusp_eq_u::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "S->N (x) equilibrium suspension flux (kg/m/s)"
    Qsusp_eq_v::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Sublimation flux (kg/m^2/s)"
    Qsubl::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Suspension-layer sublimation flux (kg/m^2/s)"
    Qsubl_susp::Matrix{Tf} = zeros(Tf, Nx, Ny)
    "Saltation-layer sublimation flux (kg/m^2/s)"
    Qsubl_salt::Matrix{Tf} = zeros(Tf, Nx, Ny)

    "Wind index array E"
    index_ue::Matrix{Int} = zeros(Int, Nx, 2 * Ny + 1)
    "Wind index array W"
    index_uw::Matrix{Int} = zeros(Int, Nx, 2 * Ny + 1)
    "Wind index array N"
    index_vn::Matrix{Int} = zeros(Int, Ny, 2 * Nx + 1)
    "Wind index array S"
    index_vs::Matrix{Int} = zeros(Int, Ny, 2 * Nx + 1)

end
