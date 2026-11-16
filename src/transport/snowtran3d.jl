const LIBSNOWTRAN3D = joinpath(@__DIR__, "..", "..", "deps", "libsnowtran3d")

"""
$(TYPEDSIGNATURES)

Snow transport by wind using Liston's SnowTran3D model.

Implementation of Liston and Sturm (1998) and Liston et al. (2007) SnowTran3D model via Fortran ccall.
The upstream Fortran routines are kept verbatim in deps/ and are entered through
the standalone wrapper in deps/SNOWTRAN3D_WRAPPER.F90. Model parameters
(flag_variable_Utau_t, Utau_t_const, rho_snow, blowby, ...) are hardcoded
constants in deps/MODULES.F90.

# Arguments
- `fsm::FSM`: Model state structure
- `met::MET`: Meteo variable structure
- `w::SnowTransport`: Transport workspace (static arrays, accumulators, constants)
- `snowdepth0::Matrix`: Snow depth of deposited snow (m) - modified in-place
- `Sice0::Matrix`: Ice content of deposited snow (kg/m²) - modified in-place
- `dSWE_salt::Matrix`: SWE change due to saltation (kg/m²) - output
- `dSWE_susp::Matrix`: SWE change due to suspension (kg/m²) - output
- `dSWE_subl::Matrix`: SWE change due to sublimation (kg/m²) - output
"""
function snowtran3d!(
        fsm::FSM{Tf}, met::MET{Tf}, w::SnowTransport{Tf}, snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        dSWE_salt::Matrix{Tf}, dSWE_susp::Matrix{Tf},
        dSWE_subl::Matrix{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny, Nsmax) = fsm.grid
    (; dt, zRH, zU, Ds_min) = fsm.params
    (; rhos_min, rhos_max) = w
    Ua_eff = fsm.diag.Uaeff
    (; Udir, Ta, RH) = met
    (; vegsnowd_xy, tilefrac, glacierfrac) = w
    (; z0_snow, z0sf, Ld) = fsm.surface
    (; fsnow, Ds, Sice, Sliq, Tsnow, histowet, Nsnow) = fsm.state
    (; dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp) = w

    # The Fortran wrapper is compiled for 4-byte integers, so the (in-out) snow-layer count
    # crosses the boundary as an Int32 copy that is read back afterwards. Integer scalars
    # convert through their `Ref{Int32}` slots automatically.
    Nsnow32 = Int32.(Nsnow)

    ccall(
        (:snowtran3d_wrapper_, LIBSNOWTRAN3D),
        Cvoid,
        (
            Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Tf},   # Nx, Ny, Nsmax, Ds_min
            Ref{Tf}, Ref{Tf}, Ref{Tf},                     # dt, zU, zRH
            Ref{Tf}, Ref{Tf},                              # rhos_min, rhos_max
            Ptr{Tf}, Ptr{Tf},                              # snowdepth0, Sice0
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                     # dSWE_salt, dSWE_susp, dSWE_subl
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf}, Ptr{Tf},            # Ua, Udir, Ta, RH
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                     # veg_shd, z0_snow, z0sf
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf}, Ptr{Tf},            # fsnow, Ds, Sice, Sliq
            Ptr{Tf}, Ptr{Tf}, Ptr{Int32},                  # Tsnow, histowet, Nsnow
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                     # dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp
            Ptr{Tf}, Ptr{Tf}, Ptr{Tf},                     # cellsize, tilefrac, glacierfrac
        ),
        Nx, Ny, Nsmax, Ds_min,
        dt, zU, zRH,
        rhos_min, rhos_max,
        snowdepth0, Sice0,
        dSWE_salt, dSWE_susp, dSWE_subl,
        Ua_eff, Udir, Ta, RH,
        vegsnowd_xy, z0_snow, z0sf,
        fsnow, Ds, Sice, Sliq,
        Tsnow, histowet, Nsnow32,
        dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp,
        Ld, tilefrac, glacierfrac
    )

    Nsnow .= Nsnow32

    return nothing
end
