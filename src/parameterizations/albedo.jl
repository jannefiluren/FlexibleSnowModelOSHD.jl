"""
$(TYPEDEF)

Diagnostic snow albedo parameterization: the albedo is a linear function of the surface
temperature, from `alb_min` at the melting point and above to `alb_max` at `T_cold` and
below, without any memory of earlier time steps.

```jldoctest
using FlexibleSnowModelOSHD

DiagnosticAlbedo{Float32}(; alb_max = 0.85)

# output
DiagnosticAlbedo
├── alb_min: 0.6
├── alb_max: 0.85
└── T_cold: 271.15
```

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct DiagnosticAlbedo{Tf} <: AbstractAlbedo{Tf}
    "Minimum albedo, reached at and above the melting point (-)"
    alb_min::Tf = 0.6
    "Maximum albedo, reached at and below `T_cold` (-)"
    alb_max::Tf = 0.86
    "Surface temperature at and below which the albedo is `alb_max`; must be below the melting point (K)"
    T_cold::Tf = 271.15
end

"""
$(TYPEDEF)

Snow albedo parameterization with exponential decay towards `alb_min` and refresh towards
`alb_fresh` by snowfall. Decay is faster for melting than for cold snow, fixed at
`tau_summer` in the melt season, and adjusted below canopy by the transmitted shortwave
and longwave radiation.

```jldoctest
using FlexibleSnowModelOSHD

grid = Grid(Float64, Nx = 2, Ny = 2)
DecayAlbedo{Float64}(grid; tau_summer = 50)

# output
DecayAlbedo
├── alb_min: 0.6
├── alb_fresh: (0.86, 0.86)
├── tau_cold: 1000.0
├── tau_melt: 100.0
├── tau_summer: 50.0
├── c_canopy_sw: 3.0
├── c_canopy_lw: 2.0
└── dswe_fresh: 10.0
```

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct DecayAlbedo{Tf, GT, MF <: AbstractMatrix{<:AbstractFloat}} <: AbstractAlbedo{Tf}
    "Model grid, sizes the per-cell fields"
    grid::GT

    # Albedo bounds
    "Minimum albedo (-)"
    alb_min::Tf = 0.6
    "Fresh-snow albedo, per cell; must be at least `alb_min` (-)"
    alb_fresh::MF = fill(0.86, grid.Nx, grid.Ny)

    # Decay
    "Decay time scale for cold snow (h)"
    tau_cold::Tf = 1000
    "Decay time scale for melting snow (h)"
    tau_melt::Tf = 100
    "Decay time scale in the melt season, overriding `tau_cold` and `tau_melt` (h)"
    tau_summer::Tf = 70
    "Canopy adjustment of the decay time, shortwave (-)"
    c_canopy_sw::Tf = 3
    "Canopy adjustment of the decay time, longwave (-)"
    c_canopy_lw::Tf = 2

    # Fresh snow
    "Snowfall over which the albedo relaxes towards `alb_fresh` (e-folding amount) (kg/m^2)"
    dswe_fresh::Tf = 10
end

"""
$(TYPEDEF)

Prognostic snow albedo parameterization: the albedo decays exponentially towards `alb_min`
for melting snow and linearly for cold snow, and is refreshed towards `alb_fresh` by snowfall.
With `aspect_tuning`, both decay times are shortened on slopes receiving more direct shortwave
radiation than flat ground, by the ratio `Sdird / Sdir`. With a nonzero `aspect_pert`, the decay
increment is scaled by [`aspect_perturbation`](@ref): faster decay on sunny slopes, slower on shaded
ones. Both can be combined; `aspect_tuning` then shortens the decay times before the increment is
scaled. For thin snow (SWE below `swe_thin`) the fresh-snow albedo is scaled by `c_thin`.

```jldoctest
using FlexibleSnowModelOSHD

grid = Grid(Float64, Nx = 2, Ny = 2)
PrognosticAlbedo{Float64}(grid; tau_melt = 200)

# output
PrognosticAlbedo
├── alb_min: 0.6
├── alb_fresh: (0.86, 0.86)
├── tau_melt: 200.0
├── tau_cold: (1000.0, 1000.0)
├── aspect_tuning: true
├── aspect_pert: 0.0
├── aspect_pert_base: 20.0
├── dswe_fresh: 10.0
├── swe_thin: 75.0
└── c_thin: 0.8
```

# Fields

$(TYPEDFIELDS)
"""
@kwdef struct PrognosticAlbedo{Tf, GT, MF <: AbstractMatrix{<:AbstractFloat}} <: AbstractAlbedo{Tf}
    "Model grid, sizes the per-cell fields"
    grid::GT

    # Albedo bounds
    "Minimum albedo (-)"
    alb_min::Tf = 0.6
    "Fresh-snow albedo, per cell (-)"
    alb_fresh::MF = fill(0.86, grid.Nx, grid.Ny)

    # Decay
    "Decay time scale for melting snow (h)"
    tau_melt::Tf = 100
    "Decay time scale for cold snow, per cell (h)"
    tau_cold::MF = fill(1000, grid.Nx, grid.Ny)
    "Shorten the decay times on slopes by the flat-to-inclined direct SW radiation ratio `Sdird / Sdir`"
    aspect_tuning::Bool = true
    "Steepness of the aspect perturbation of the decay increment; 0 disables it (-)"
    aspect_pert::Tf = 0
    "Base of the aspect perturbation factor, which ranges over [1 / base, base] (-)"
    aspect_pert_base::Tf = 20

    # Fresh snow
    "Snowfall scale of the refresh towards `alb_fresh`; a larger 24 h snowfall resets the albedo (kg/m^2)"
    dswe_fresh::Tf = 10

    # Thin snow
    "SWE below which the fresh-snow albedo is reduced for thin and patchy snow (kg/m^2)"
    swe_thin::Tf = 75
    "Scaling of the fresh-snow albedo for thin snow (-)"
    c_thin::Tf = 0.8
end

DiagnosticAlbedo{Tf}(grid::Grid; kwargs...) where {Tf} = DiagnosticAlbedo{Tf}(; kwargs...)
DecayAlbedo{Tf}(grid::Grid; kwargs...) where {Tf} = DecayAlbedo{Tf, typeof(grid), Matrix{Tf}}(; grid, kwargs...)
PrognosticAlbedo{Tf}(grid::Grid; kwargs...) where {Tf} = PrognosticAlbedo{Tf, typeof(grid), Matrix{Tf}}(; grid, kwargs...)

@adapt_structure DecayAlbedo
@adapt_structure PrognosticAlbedo

"""
$(TYPEDSIGNATURES)

Factor scaling the albedo decay increment by the ratio of direct short wave radiation on the 
inclined and the flat surface: `base^(atan(steepness * log(sw_dir_incl / sw_dir_hor)) / (pi / 2))`, 
in [1 / base, base]. Without direct radiation on both surfaces the factor is 1; on only one 
of them it defaults to the upper/lower bounds (`1 / base` without radiation on the slope, 
`base` without radiation on flat ground).
"""
@inline function aspect_perturbation(steepness::Tf, base::Tf, sw_dir_incl::Tf, sw_dir_hor::Tf) where {Tf}
    if sw_dir_incl < eps(Tf)
        x = sw_dir_hor < eps(Tf) ? zero(Tf) : -one(Tf)
    elseif sw_dir_hor < eps(Tf)
        x = one(Tf)
    else
        x = atan(steepness * log(sw_dir_incl / sw_dir_hor)) / (Tf(pi) / 2)
    end
    return base^x
end

"""
    snow_albedo!(scheme, i, j, state, surface, meteo, params, summer_decay)

Update the snow albedo `albs[i, j]` of cell `(i, j)`. Implemented for every `AbstractAlbedo`.
"""
function snow_albedo! end

@inline function snow_albedo!(scheme::DiagnosticAlbedo{Tf}, i, j, state, surface, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf) = state
    (; alb_min, alb_max, T_cold) = scheme

    a = alb_min + (alb_max - alb_min) * (Tsrf[i, j] - Tm) / (T_cold - Tm)
    albs[i, j] = clamp(a, alb_min, alb_max)
    return nothing
end

@inline function snow_albedo!(scheme::DecayAlbedo{Tf}, i, j, state, surface, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf) = state
    (; fveg, trcn, fsky) = surface
    (; Sdir, Sdif, Sf, Tv) = meteo
    (; dt) = params
    (; alb_min, tau_cold, tau_melt, tau_summer, c_canopy_sw, c_canopy_lw, dswe_fresh) = scheme
    alb_fresh = scheme.alb_fresh[i, j]

    # Decay time scale (s): the fixed melt-season value takes precedence over the
    # temperature-dependent values for cold and melting snow
    if summer_decay
        tau = Tf(3600) * tau_summer
    else
        tau = Tf(3600) * (Tsrf[i, j] >= Tm ? tau_melt : tau_cold)
    end

    # Canopy adjustment of the decay time scale
    if fveg[i, j] > Tf(0)
        if Sdir[i, j] > eps(Tf)
            tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) * (Tf(1) + c_canopy_lw * Tv[i, j]) + c_canopy_sw * Tv[i, j])
        elseif Sdif[i, j] > eps(Tf)
            tau = tau / ((Tf(1) - trcn[i, j] * fsky[i, j]) + c_canopy_sw * trcn[i, j] * fsky[i, j])
        elseif Sdir[i, j] + Sdif[i, j] <= eps(Tf)
            tau = tau / (Tf(2) - trcn[i, j] * fsky[i, j])
        end
    end

    # Albedo budget d(albedo)/dt = -(albedo - alb_min) / tau - Sf / dswe_fresh * (albedo - alb_fresh):
    # decay towards alb_min plus refresh towards alb_fresh, i.e. relaxation towards the equilibrium
    # albedo alb_eq at the combined rate. Solved exactly over the time step, so it is stable for any dt.
    rate = Tf(1) / tau + Sf[i, j] / dswe_fresh
    alb_eq = (alb_min / tau + Sf[i, j] * alb_fresh / dswe_fresh) / rate
    a = alb_eq + (albs[i, j] - alb_eq) * exp(-rate * dt)

    albs[i, j] = clamp(a, alb_min, alb_fresh)
    return nothing
end

@inline function snow_albedo!(scheme::PrognosticAlbedo{Tf}, i, j, state, surface, meteo, params, summer_decay) where {Tf}
    @unpack_constants(Tf)
    (; albs, Tsrf, Sice, Sliq) = state
    (; Sdir, Sdird, Sf, Sf24h) = meteo
    (; dt) = params
    (; alb_min, tau_melt, aspect_pert, aspect_pert_base, dswe_fresh, swe_thin, c_thin) = scheme
    tau_cold = scheme.tau_cold[i, j]
    alb_fresh = scheme.alb_fresh[i, j]

    # Aspect tuning: faster decay on slopes receiving more direct radiation than flat ground
    if scheme.aspect_tuning && Sdir[i, j] > eps(Tf) && Sdird[i, j] < Sdir[i, j]
        tau_melt = max(tau_melt * Sdird[i, j] / Sdir[i, j], eps(Tf))
        tau_cold = max(tau_cold * Sdird[i, j] / Sdir[i, j], eps(Tf))
    end

    # Decay increment: exponential towards alb_min for melting snow, linear for cold snow
    if Tsrf[i, j] >= Tm
        dalb = (albs[i, j] - alb_min) * (exp(-(dt / Tf(3600)) / tau_melt) - Tf(1))
    else
        dalb = -(dt / Tf(3600)) / tau_cold
    end

    # Aspect perturbation: faster decay on sunny slopes, slower on shaded ones
    if aspect_pert != Tf(0)
        dalb *= aspect_perturbation(aspect_pert, aspect_pert_base, Sdir[i, j], Sdird[i, j])
    end
    a = albs[i, j] + dalb

    # Reduced fresh-snow albedo for thin and patchy snow
    swe = zero(Tf)
    for k in axes(Sice, 1)
        swe += Sice[k, i, j] + Sliq[k, i, j]
    end
    if swe < swe_thin
        alb_fresh *= c_thin
    end

    # Refresh by snowfall: reset after a large 24 h snowfall, partial refresh otherwise
    if Sf[i, j] * dt > Tf(0) && Sf24h[i, j] > dswe_fresh
        a = alb_fresh
    else
        a = a + (alb_fresh - a) * Sf[i, j] * dt / dswe_fresh
    end

    albs[i, j] = max(min(a, alb_fresh), alb_min)
    return nothing
end
