"""
$(TYPEDSIGNATURES)

Compute the snow depth (m) corresponding to the top `swe` kg/m^2 of the
snowpack at pixel (i, j).

Pure Julia translation of the Fortran routine `HS_FROM_SWE` (deps/HS_FROM_SWE.F90).
"""
function hs_from_swe(
        fsm::FSM{Tf}, w::SnowTransport{Tf}, swe::Tf, i::Integer, j::Integer
    ) where {Tf <: Real}

    (; Nsnow, fsnow, Sice, Sliq, Ds) = fsm.state
    (; rhos_min, rhos_max, rho_snow) = w

    # Epsilon SWE tolerance to avoid instabilities
    eps_swe = Tf(1.0e-4)

    rho_avg = rho_snow
    hs = zero(Tf)
    Ds_tmp = zero(Tf)

    swe_tot = column_swe(Sice, Sliq, i, j)

    @inbounds if swe > eps(Tf) && swe <= swe_tot + eps_swe

        if swe >= swe_tot

            # Fix computing approximations to avoid instabilities
            Ds_tmp = column_sum(Ds, i, j)

        else

            # Normal case
            k = 1
            swe_tmp = zero(Tf)

            while k <= Nsnow[i, j] && swe_tmp < swe
                swe_layer = Sice[k, i, j] + Sliq[k, i, j]
                if swe - swe_tmp > swe_layer
                    Ds_tmp += Ds[k, i, j]
                    swe_tmp += swe_layer
                else
                    weight = (swe - swe_tmp) / swe_layer
                    Ds_tmp += weight * Ds[k, i, j]
                    swe_tmp = swe
                end
                k += 1
            end

        end

        hs = Ds_tmp * fsnow[i, j]

        if hs > eps(Tf)
            rho_avg = swe / hs
        elseif hs < -eps(Tf)
            error("hs_from_swe: hs < 0 at ($i, $j): hs = $hs")
        end

    elseif swe < -eps(Tf)

        error("hs_from_swe: swe < 0 at ($i, $j): swe = $swe")

    elseif swe > swe_tot + eps_swe

        error("hs_from_swe: swe > swe_tot at ($i, $j): swe = $swe, swe_tot = $swe_tot")

    end

    if (rho_avg < rhos_min - Tf(0.5) || rho_avg > rhos_max + Tf(0.5)) && Ds_tmp > Tf(0.001)
        error("hs_from_swe: invalid density $rho_avg at ($i, $j)")
    end

    return hs
end

# Upwind / downwind index helpers, mirroring module FUNCTIONS_SNOWTRAN3D.
# wnd_dir: 1 = West, 2 = East, 3 = South, 4 = North.
# The index array holds, for each row, the number of same-sign wind runs in
# column 1 and the (start, end) index pair of run k in columns 2k, 2k+1.

@inline get_kstart(wnd_dir::Int, n::Int, k::Int, idx::Matrix{Int}) =
    (wnd_dir == 1 || wnd_dir == 3) ? idx[n, k * 2] + 1 : idx[n, k * 2 + 1] - 1

@inline get_kend(wnd_dir::Int, n::Int, k::Int, idx::Matrix{Int}) =
    (wnd_dir == 1 || wnd_dir == 3) ? idx[n, k * 2 + 1] : idx[n, k * 2]

@inline get_i_upwind(wnd_dir::Int, i::Int) =
    (wnd_dir == 1 || wnd_dir == 2) ? i : (wnd_dir == 3 ? i - 1 : i + 1)

@inline get_j_upwind(wnd_dir::Int, j::Int) =
    wnd_dir == 1 ? j - 1 : (wnd_dir == 2 ? j + 1 : j)

@inline get_i_downwind(wnd_dir::Int, i::Int) =
    (wnd_dir == 1 || wnd_dir == 2) ? i : (wnd_dir == 3 ? i + 1 : i - 1)

@inline get_j_downwind(wnd_dir::Int, j::Int) =
    wnd_dir == 1 ? j + 1 : (wnd_dir == 2 ? j - 1 : j)

"""
$(TYPEDSIGNATURES)

One directional sweep of the wind-run indexing: for each row of `wind_data`
(laid out `(x, y)`), find the runs of grid cells whose wind component keeps the
target sign. `positive` selects runs of positive (negative) wind; `swap` selects
the S/N sweeps (transport along x) over the W/E sweeps (transport along y).

Translation of the Fortran subroutine `SWEEP_DIRECTION` (deps/SNOWTRAN3D.F90).
"""
function sweep_direction!(
        index_array::Matrix{Int}, wind_data::Matrix{Tf}, dim1_size::Int, dim2_size::Int,
        positive::Bool, swap::Bool
    ) where {Tf <: Real}

    @inbounds for dim1 in 1:dim1_size

        w1 = swap ? wind_data[1, dim1] : wind_data[dim1, 1]
        sign_prev = signbit(w1) ? -one(Tf) : one(Tf)

        target_prev = (positive && sign_prev > 0) || (!positive && sign_prev < 0)
        if target_prev
            npairs = 1
            index_array[dim1, 2] = 1
        else
            npairs = 0
        end

        for dim2 in 2:dim2_size
            wc = swap ? wind_data[dim2, dim1] : wind_data[dim1, dim2]
            sign_curr = signbit(wc) ? -one(Tf) : one(Tf)
            if sign_curr != sign_prev
                if (positive && sign_curr > 0) || (!positive && sign_curr < 0)
                    npairs += 1
                    index_array[dim1, npairs * 2] = dim2
                else
                    index_array[dim1, npairs * 2 + 1] = dim2 - 1
                end
            end
            sign_prev = sign_curr
        end

        wlast = swap ? wind_data[dim2_size, dim1] : wind_data[dim1, dim2_size]
        sign_last = signbit(wlast) ? -one(Tf) : one(Tf)
        if (positive && sign_last > 0) || (!positive && sign_last < 0)
            index_array[dim1, npairs * 2 + 1] = dim2_size
        end

        index_array[dim1, 1] = npairs
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Build the four wind-run index arrays (W, E, S, N) for the current wind field.

Translation of the Fortran subroutine `COMPUTE_WIND_DIRECTION_ARRAYS`.
"""
function compute_wind_direction_arrays!(w::SnowTransport{Tf}) where {Tf <: Real}
    (; Nx, Ny, uwind, vwind, index_ue, index_uw, index_vn, index_vs) = w

    sweep_direction!(index_uw, uwind, Nx, Ny, true, false)
    sweep_direction!(index_ue, uwind, Nx, Ny, false, false)
    sweep_direction!(index_vs, vwind, Ny, Nx, true, true)
    sweep_direction!(index_vn, vwind, Ny, Nx, false, true)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Determine the soft (movable) snow thickness `w.Ds_soft` from the layer wetting
history: layers from the top are soft while they are dry (`Sliq == 0`) and have
never been wet (`histowet < 0.5`).

Translation of the Fortran subroutine `COMPUTE_SOFT_SNOW` (deps/SNOWTRAN3D.F90).
"""
function compute_soft_snow!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; Ds, Nsnow, histowet, Sliq) = fsm.state
    (; Ds_soft) = w

    fill!(Ds_soft, zero(Tf))
    @inbounds for j in 1:Ny, i in 1:Nx
        k = 1
        while k <= Nsnow[i, j] && Sliq[k, i, j] < eps(Tf) && histowet[k, i, j] < Tf(0.5)
            Ds_soft[i, j] += Ds[k, i, j]
            k += 1
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Update `w.Ds_soft` by removing layers that are too dense to be moved by the
friction velocity `w.Utau`.

Translation of the Fortran subroutine `UPDATE_SOFT_SNOW` (deps/SNOWTRAN3D.F90).
"""
function update_soft_snow!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; Ds, Nsnow, fsnow, Sice, Sliq) = fsm.state
    (; Utau, Ds_soft) = w

    @inbounds for j in 1:Ny, i in 1:Nx

        Ds_soft_old = Ds_soft[i, j]
        Ds_soft[i, j] = zero(Tf)
        hard_layer_reached = false
        k = 1

        while k <= Nsnow[i, j] && Ds_soft_old - Ds_soft[i, j] > eps(Tf) && !hard_layer_reached

            rho_layer = (Sice[k, i, j] + Sliq[k, i, j]) / Ds[k, i, j] / fsnow[i, j]

            if Utau[i, j] - threshold_friction_velocity(rho_layer) > eps(Tf)
                Ds_soft[i, j] += Ds[k, i, j]
            else
                hard_layer_reached = true
            end

            k += 1

        end

    end

    return nothing
end

# Threshold friction velocity (m/s) of snow with density `rho` (kg/m^3)
@inline function threshold_friction_velocity(rho::Tf) where {Tf <: Real}
    return rho <= Tf(300) ? Tf(0.1) * exp(Tf(0.003) * rho) : Tf(0.005) * exp(Tf(0.013) * rho)
end

"""
$(TYPEDSIGNATURES)

Compute the threshold friction velocity `w.Utau_t` from the density of the top
snow layer.

Translation of the Fortran subroutine `COMPUTE_THRESHOLD_FRICTION_VELOCITY`.
"""
function compute_threshold_friction_velocity!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; Nsnow, fsnow, Sice, Sliq, Ds) = fsm.state
    (; rho_snow, Utau_t) = w

    @inbounds for j in 1:Ny, i in 1:Nx

        rho_surf_snow = rho_snow
        if Nsnow[i, j] > 0 && Ds[1, i, j] > eps(Tf)
            rho_surf_snow = (Sice[1, i, j] + Sliq[1, i, j]) / Ds[1, i, j] / fsnow[i, j]
        end

        Utau_t[i, j] = threshold_friction_velocity(rho_surf_snow)

    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Newton iterations for the friction velocity under saltation. `threshold_flag`
selects the high-wind roughness relation `z0 = 0.00734 u* - 0.0022` over the
Liston and Sturm (1998) eq. 5 relation.

Translation of the Fortran subroutine `FRICTION_VELOCITY_SOLVER`.
"""
function friction_velocity_solver(
        zU::Tf, guess::Tf, windtmp::Tf, threshold_flag::Bool
    ) where {Tf <: Real}

    @unpack_constants(Tf)

    # Coefficient 0.12 in Liston and Sturm (1998) eq. 5 p. 500
    C_z = Tf(0.12)

    tol = Tf(1.0e-3)
    maxiter = 20
    old = guess
    xnew = old

    if threshold_flag

        old = Tf(0.6)
        for _ in 1:maxiter
            denom = log(zU) - log(Tf(0.00734) * old - Tf(0.0022))
            fprime = -Tf(1) + windtmp * vkman * Tf(0.00734) / (Tf(0.00734) * old - Tf(0.0022)) / (denom * denom)
            funct = -old + windtmp * vkman / denom
            xnew = old - funct / fprime
            if abs(xnew - old) < tol
                return xnew
            end
            old = xnew
        end

    else

        for _ in 1:maxiter
            denom = log(zU) - log(C_z / (Tf(2) * grav)) - Tf(2) * log(old)
            fprime = -Tf(1) + Tf(2) / old * windtmp * vkman / (denom * denom)
            funct = -old + windtmp * vkman / denom
            xnew = old - funct / fprime
            if abs(xnew - old) < tol
                return xnew
            end
            old = xnew
        end

    end

    @warn "friction_velocity_solver: max iteration exceeded, Utau = $old" maxlog = 10

    return xnew
end

"""
$(TYPEDSIGNATURES)

Solve for the friction velocity `w.Utau`, surface roughness length `w.z_0` and
saltation-layer height `w.h_star`. Returns `true` if snow is saltating anywhere
in the domain.

Translation of the Fortran subroutine `COMPUTE_FRICTION_VELOCITY`.
"""
function compute_friction_velocity!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; zU) = fsm.params
    (; z0_snow, z0sf) = fsm.surface
    Ua_eff = fsm.diag.Uaeff
    (; fsnow) = fsm.state
    (; vegsnowd_xy, Utau, Utau_t, z_0, h_star, snowthickness, Ds_soft) = w
    @unpack_constants(Tf)

    # Coefficients in Liston and Sturm (1998), eq. 5 and eq. 14
    C_z = Tf(0.12)
    h_const = Tf(1.6)

    flag_blowing_snow = false

    guess = Tf(0.1)
    @inbounds for j in 1:Ny, i in 1:Nx

        if snowthickness[i, j] <= vegsnowd_xy[i, j] || Ds_soft[i, j] <= eps(Tf)

            # Saltation will not occur
            z_0[i, j] = fsnow[i, j] * z0_snow[i, j] + (Tf(1) - fsnow[i, j]) * z0sf[i, j]
            z_0_tmp = min(Tf(0.25) * zU, z_0[i, j])
            Utau[i, j] = Ua_eff[i, j] * vkman / log(zU / z_0_tmp)
            h_star[i, j] = z_0[i, j] * h_const / C_z

        else

            # Saltation may occur: solve assuming saltation and compare Utau with Utau_t.
            windtmp = max(Tf(1), Ua_eff[i, j])
            wind_max = Tf(30) * log(zU / z0_snow[i, j]) / log(Tf(10) / z0_snow[i, j])
            windtmp = min(windtmp, wind_max)

            threshold = Tf(0.6) / vkman * log(zU / Tf(0.0022))
            threshold_flag = windtmp > threshold

            Utautmp = friction_velocity_solver(zU, guess, windtmp, threshold_flag)

            if Utautmp > Utau_t[i, j]

                # We have saltation
                Utau[i, j] = Utautmp
                z_0[i, j] = C_z * Utau[i, j]^2 / (Tf(2) * grav)
                h_star[i, j] = h_const * Utau[i, j]^2 / (Tf(2) * grav)
                flag_blowing_snow = true

            else

                # No saltation, but vegetation is snow-covered. Cap Utau at Utau_t.
                z_0[i, j] = z0_snow[i, j]
                Utau[i, j] = Ua_eff[i, j] * vkman / log(zU / z_0[i, j])
                Utau[i, j] = min(Utau[i, j], Utau_t[i, j])
                h_star[i, j] = z_0[i, j] * h_const / C_z

            end

        end

    end

    return flag_blowing_snow
end

"""
$(TYPEDSIGNATURES)

Sublimation loss rate coefficients at height `z`, following the appendix of
Liston and Sturm (1998). Returns `(V_susp, V_salt)`; `flag = 1` computes the
suspension-layer coefficient, `flag = 0` the saltation-layer coefficient.

Translation of the Fortran subroutine `SUBLIMATION` (deps/SNOWTRAN3D.F90).
"""
function getsublim(
        zRH::Tf, z::Tf, RH::Tf, Ta::Tf, Utau::Tf, z_0::Tf, Utau_t::Tf, flag::Tf
    ) where {Tf <: Real}

    @unpack_constants(Tf)

    # Constants from the Fortran CONSTANTS module not present in get_constants
    Runi = Tf(8313)         # Universal gas constant (J/kmol/K)
    xM = Tf(18.01)          # Molecular weight of water (kg/kmol)
    visc_air = Tf(13.0e-6)  # Kinematic viscosity of air (m^2/s)

    D = Tf(2.06e-5) * (Ta / Tf(273))^Tf(1.75)
    rho_sat = Tf(0.622) / (Rair * Ta) * Tf(610.78) * exp(Tf(21.875) * (Ta - Tm) / (Ta - Tf(7.66)))

    # Assume that the rh varies according to a modification to
    # Pomeroy's humidity variation with height equation.
    rh_offset = Tf(1) - Tf(0.027) * log(zRH)
    sigma = (Tf(0.01) * RH - Tf(1)) * (rh_offset + Tf(0.027) * log(z))
    sigma = min(Tf(0), sigma)
    sigma = max(Tf(-1), sigma)

    alpha = Tf(4.08) + Tf(12.6) * z
    rbar_r = Tf(4.6e-5) * z^Tf(-0.258)
    xmbar = Tf(4) / Tf(3) * pi * rho_ice * rbar_r^3 * (Tf(1) + Tf(3) / alpha + Tf(2) / alpha^2)
    rbar = ((Tf(3) * xmbar) / (Tf(4) * pi * rho_ice))^Tf(0.33)
    u_z = Utau / vkman * log(z / z_0)
    x_r = Tf(0.005) * u_z^Tf(1.36)
    wbar = Tf(1.1e7) * rbar^Tf(1.8)

    V_susp = zero(Tf)
    V_salt = zero(Tf)

    if flag == Tf(1)

        V_r = wbar + Tf(3) * x_r * cos(pi / Tf(4))
        xN_r = Tf(2) * rbar * V_r / visc_air
        xNu = Tf(1.79) + Tf(0.606) * xN_r^Tf(0.5)
        xSh = xNu
        tmp1 = (Ls * xM) / (Runi * Ta) - Tf(1)
        tmp2 = hcon_air * Ta * xNu
        top = Tf(2) * pi * rbar * sigma
        bottom = Ls / tmp2 * tmp1 + Tf(1) / (D * rho_sat * xSh)
        V_susp = (top / bottom) / xmbar

    elseif flag == Tf(0)

        V_rsalt = Tf(0.68) * Utau + Tf(2.3) * Utau_t
        xN_r = Tf(2) * rbar * V_rsalt / visc_air
        xNu = Tf(1.79) + Tf(0.606) * xN_r^Tf(0.5)
        xSh = xNu
        tmp1 = (Ls * xM) / (Runi * Ta) - Tf(1)
        tmp2 = hcon_air * Ta * xNu
        top = Tf(2) * pi * rbar * sigma
        bottom = Ls / tmp2 * tmp1 + Tf(1) / (D * rho_sat * xSh)
        V_salt = (top / bottom) / xmbar

    end

    return V_susp, V_salt
end

"""
$(TYPEDSIGNATURES)

Compute the local-equilibrium suspension flux components `w.Qsusp_u`/`w.Qsusp_v`
(vertical quadrature of the suspended snow concentration following Kind, 1992)
together with the suspension- and saltation-layer sublimation fluxes
`w.Qsubl_susp`/`w.Qsubl_salt`.

Translation of the Fortran subroutine `SUSPENSION` (deps/SNOWTRAN3D.F90).
"""
function suspension!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}, met::MET{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; zRH) = fsm.params
    (; Ta, RH) = met
    (; Utau, Utau_t, z_0, h_star) = w
    (; Qsalt_u, Qsalt_v, Qsusp_u, Qsusp_v, Qsubl_susp, Qsubl_salt) = w
    @unpack_constants(Tf)

    # Constants from the Fortran CONSTANTS_SNOWTRAN3D module
    Up_const = Tf(2.8)      # Constant coefficient for calculation of U_p
    dz_susp = Tf(0.2)       # dz in suspension layer (m)
    ztop_susp = Tf(2.0)     # Height of the top of suspension layer (m)
    fall_vel = Tf(0.3)      # Particle-settling velocity (m/s)
    Ur_const = Tf(0.5)      # Constant coefficient beta in L&S 1998, eq. 13

    @inbounds for j in 1:Ny, i in 1:Nx

        conc_salt_WE = zero(Tf)
        conc_salt_SN = zero(Tf)
        conc_salt = zero(Tf)

        if Qsalt_u[i, j] > eps(Tf) || Qsalt_v[i, j] > eps(Tf)

            U_r = Utau[i, j] / vkman * log(h_star[i, j] / z_0[i, j])
            phistar_Cr = Utau[i, j] / U_r * Ur_const
            prd = phistar_Cr * Utau[i, j] / fall_vel
            U_p = Up_const * Utau_t[i, j]

            conc_salt_WE = Qsalt_u[i, j] / (h_star[i, j] * U_p)
            conc_salt_SN = Qsalt_v[i, j] / (h_star[i, j] * U_p)
            conc_salt = sqrt(conc_salt_WE^2 + conc_salt_SN^2)

            nzsteps = trunc(Int, (ztop_susp - h_star[i, j]) / dz_susp)

            Qsusp_u[i, j] = zero(Tf)
            Qsusp_v[i, j] = zero(Tf)
            Qsubl_susp[i, j] = zero(Tf)
            Qsubl_salt[i, j] = zero(Tf)

            for iz in 1:nzsteps
                z = h_star[i, j] + Tf(0.5) * dz_susp + Tf(iz - 1) * dz_susp

                zratio = (z / h_star[i, j])^((-fall_vel) / (vkman * Utau[i, j]))
                conc_WE = conc_salt_WE * ((prd + Tf(1)) * zratio - prd)
                conc_WE = max(conc_WE, zero(Tf))
                conc_SN = conc_salt_SN * ((prd + Tf(1)) * zratio - prd)
                conc_SN = max(conc_SN, zero(Tf))
                conc = sqrt(conc_WE^2 + conc_SN^2)

                if conc_WE > eps(Tf)
                    if z - z_0[i, j] < eps(Tf)
                        z = Tf(1.2) * z_0[i, j]
                    end
                    Qsusp_u[i, j] += conc_WE * log(z / z_0[i, j]) * dz_susp
                end

                if conc_SN > eps(Tf)
                    if z - z_0[i, j] < eps(Tf)
                        z = Tf(1.2) * z_0[i, j]
                    end
                    Qsusp_v[i, j] += conc_SN * log(z / z_0[i, j]) * dz_susp
                end

                if conc > eps(Tf)
                    V_susp, _ = getsublim(
                        zRH, z, RH[i, j], Ta[i, j], Utau[i, j],
                        z_0[i, j], Utau_t[i, j], Tf(1)
                    )
                    Qsubl_susp[i, j] += conc * V_susp * dz_susp
                end

            end

            # Finish the quadratures; include the constants for Qsusp
            Qsusp_u[i, j] = Utau[i, j] / vkman * Qsusp_u[i, j]
            Qsusp_v[i, j] = Utau[i, j] / vkman * Qsusp_v[i, j]

            # Include the sublimation contribution due to saltation
            z = h_star[i, j] / Tf(2)
            _, V_salt = getsublim(
                zRH, z, RH[i, j], Ta[i, j], Utau[i, j],
                z_0[i, j], Utau_t[i, j], Tf(0)
            )
            Qsubl_salt[i, j] = V_salt * conc_salt * h_star[i, j]

        else
            Qsusp_u[i, j] = zero(Tf)
            Qsusp_v[i, j] = zero(Tf)
            Qsubl_susp[i, j] = zero(Tf)
            Qsubl_salt[i, j] = zero(Tf)
        end

    end

    return nothing
end

# Flux leaving cell (ic, jc): zero unless both it and its downwind neighbour (idw, jdw) are
# open (openfrac >= 0.1), then scaled by the downwind open fraction
@inline function open_flux(Q::Tf, openfrac::Matrix{Tf}, ic, jc, idw, jdw) where {Tf <: Real}
    Nx, Ny = size(openfrac)
    if openfrac[ic, jc] < Tf(0.1)
        Q = zero(Tf)
    end
    if idw >= 1 && idw <= Nx && jdw >= 1 && jdw <= Ny
        Q = openfrac[idw, jdw] < Tf(0.1) ? zero(Tf) : Q * openfrac[idw, jdw]
    end
    return Q
end

"""
$(TYPEDSIGNATURES)

One directional sweep of the suspension-flux relaxation toward the local
equilibrium `Qsusp_eq_dir`, mirroring the saltation sweep.

Translation of the Fortran subroutine `UPDATE_SUSPENSION_FLUX_PER_DIRECTION`.
"""
function update_suspension_flux_per_direction!(
        fsm::FSM{Tf}, w::SnowTransport{Tf},
        Qsusp_dir::Matrix{Tf}, Qsusp_eq_dir::Matrix{Tf},
        wnd_dir::Int, index_array::Matrix{Int}, dim1_size::Int
    ) where {Tf <: Real}

    (; Ld) = fsm.surface
    (; openfrac) = w

    # Equilibrium fetch distance and scaling constant for the suspension layer
    fetch_susp = Tf(1200)
    xmu_susp = Tf(3)

    incr = (wnd_dir == 1 || wnd_dir == 3) ? 1 : -1

    @inbounds for dim1 in 1:dim1_size
        for k in 1:index_array[dim1, 1]
            kstart = get_kstart(wnd_dir, dim1, k, index_array)
            kend = get_kend(wnd_dir, dim1, k, index_array)
            for dim2 in kstart:incr:kend

                ic = (wnd_dir == 1 || wnd_dir == 2) ? dim1 : dim2
                jc = (wnd_dir == 1 || wnd_dir == 2) ? dim2 : dim1
                iuw = get_i_upwind(wnd_dir, ic)
                juw = get_j_upwind(wnd_dir, jc)
                idw = get_i_downwind(wnd_dir, ic)
                jdw = get_j_downwind(wnd_dir, jc)

                scale = min(Tf(1), xmu_susp * Ld[ic, jc] / fetch_susp)
                Qsusp_dir_in = Qsusp_dir[iuw, juw]
                Qsusp_dir[ic, jc] = Qsusp_dir_in + scale * (Qsusp_eq_dir[ic, jc] - Qsusp_dir_in)

                Qsusp_dir[ic, jc] = open_flux(Qsusp_dir[ic, jc], openfrac, ic, jc, idw, jdw)

            end
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Relax the suspension flux toward its local equilibrium along the wind direction,
with advective memory on the suspension fetch.

Translation of the Fortran subroutine `RELAX_SUSPENSION_FLUX`.
"""
function relax_suspension_flux!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; Qsusp_u, Qsusp_v, Qsusp_eq_u, Qsusp_eq_v, index_ue, index_uw, index_vn, index_vs) = w

    # flag_boundary_condition = false: zero incoming suspended flux at the boundaries
    fill!(Qsusp_u, zero(Tf))
    fill!(Qsusp_v, zero(Tf))

    update_suspension_flux_per_direction!(fsm, w, Qsusp_u, Qsusp_eq_u, 1, index_uw, Nx)
    update_suspension_flux_per_direction!(fsm, w, Qsusp_u, Qsusp_eq_u, 2, index_ue, Nx)
    update_suspension_flux_per_direction!(fsm, w, Qsusp_v, Qsusp_eq_v, 3, index_vs, Ny)
    update_suspension_flux_per_direction!(fsm, w, Qsusp_v, Qsusp_eq_v, 4, index_vn, Ny)

    return nothing
end

"""
$(TYPEDSIGNATURES)

One directional sweep of the saltation flux following Liston and Sturm (1998)
eq. 9, relaxing toward the maximum flux `Qsalt_max_dir` with fetch scaling and
the downwind blowby term.

Translation of the Fortran subroutine `UPDATE_SALTATION_FLUX_PER_DIRECTION`.
"""
function update_saltation_flux_per_direction!(
        fsm::FSM{Tf}, w::SnowTransport{Tf},
        Qsalt_dir::Matrix{Tf}, Qsalt_max_dir::Matrix{Tf},
        wnd_dir::Int, index_array::Matrix{Int}, dim1_size::Int
    ) where {Tf <: Real}

    (; Ld) = fsm.surface
    (; Utau, openfrac) = w

    # Equilibrium fetch distance, scaling constant and downwind blowby fraction
    fetch = Tf(600)
    xmu = Tf(3)
    blowby = Tf(0.5)

    incr = (wnd_dir == 1 || wnd_dir == 3) ? 1 : -1

    @inbounds for dim1 in 1:dim1_size
        for k in 1:index_array[dim1, 1]
            kstart = get_kstart(wnd_dir, dim1, k, index_array)
            kend = get_kend(wnd_dir, dim1, k, index_array)
            for dim2 in kstart:incr:kend

                ic = (wnd_dir == 1 || wnd_dir == 2) ? dim1 : dim2
                jc = (wnd_dir == 1 || wnd_dir == 2) ? dim2 : dim1
                iuw = get_i_upwind(wnd_dir, ic)
                juw = get_j_upwind(wnd_dir, jc)
                idw = get_i_downwind(wnd_dir, ic)
                jdw = get_j_downwind(wnd_dir, jc)

                scale = min(Tf(1), xmu * Ld[ic, jc] / fetch)
                Qsalt_dir_in = Qsalt_dir[iuw, juw]
                dUtau = Utau[ic, jc] - Utau[iuw, juw]

                if dUtau >= eps(Tf)
                    Qsalt_dir[ic, jc] = Qsalt_dir_in + scale * (Qsalt_max_dir[ic, jc] - Qsalt_dir_in)
                else
                    if Qsalt_dir_in < Qsalt_max_dir[ic, jc]
                        Qsalt_dir[ic, jc] = Qsalt_dir_in
                    else
                        Qsalt_dir[ic, jc] = max(blowby * Qsalt_dir_in, Qsalt_max_dir[ic, jc])
                    end
                end

                Qsalt_dir[ic, jc] = open_flux(Qsalt_dir[ic, jc], openfrac, ic, jc, idw, jdw)

            end
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the saltation flux components `w.Qsalt_u`/`w.Qsalt_v` and their magnitude
`Qsalt` with four directional upwind sweeps (W, E, S, N).

Translation of the Fortran subroutine `SALTATION` (deps/SNOWTRAN3D.F90).
"""
function saltation!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; fsnow) = fsm.state
    (; Utau, Utau_t, uwind, vwind, openfrac) = w
    (; Qsalt_u, Qsalt_v, Qsalt_maxu, Qsalt_maxv) = w
    (; index_ue, index_uw, index_vn, index_vs) = w
    @unpack_constants(Tf)

    # Density of air (from the Fortran CONSTANTS module)
    rho_air = Tf(1.275)

    # Maximum possible saltation flux (Eqn. 6, Liston & Sturm 1998)
    @inbounds for j in 1:Ny, i in 1:Nx
        qmax = Tf(0.68) * rho_air / grav *
            Utau_t[i, j] / Utau[i, j] * (Utau[i, j]^2 - Utau_t[i, j]^2)
        qmax = max(qmax * fsnow[i, j] * openfrac[i, j], zero(Tf))
        wnorm = sqrt(uwind[i, j]^2 + vwind[i, j]^2)
        Qsalt_maxu[i, j] = qmax * abs(uwind[i, j]) / wnorm
        Qsalt_maxv[i, j] = qmax * abs(vwind[i, j]) / wnorm
    end

    # flag_boundary_condition = false: zero incoming flux at the boundaries
    fill!(Qsalt_u, zero(Tf))
    fill!(Qsalt_v, zero(Tf))

    update_saltation_flux_per_direction!(fsm, w, Qsalt_u, Qsalt_maxu, 1, index_uw, Nx)
    update_saltation_flux_per_direction!(fsm, w, Qsalt_u, Qsalt_maxu, 2, index_ue, Nx)
    update_saltation_flux_per_direction!(fsm, w, Qsalt_v, Qsalt_maxv, 3, index_vs, Ny)
    update_saltation_flux_per_direction!(fsm, w, Qsalt_v, Qsalt_maxv, 4, index_vn, Ny)

    return nothing
end

"""
$(TYPEDSIGNATURES)

One directional sweep of the saltation/suspension mass balance: convert the flux
divergence of `Qs_dir` into erosion and deposition along the wind, limited by the
erodible snow, with the open-fraction mass-conservation weighting.

Translation of the Fortran subroutine `COMPUTE_MASS_CHANGES_PER_DIRECTION`.
"""
function compute_mass_changes_per_direction!(
        fsm::FSM{Tf}, w::SnowTransport{Tf},
        wnd_dir::Int, index_array::Matrix{Int}, dim1_size::Int,
        dSWE_s::Matrix{Tf}, Qs_dir::Matrix{Tf},
        snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf}, Tm::Tf
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; dt) = fsm.params
    (; Ld) = fsm.surface
    (; Ds, fsnow, Sice, Sliq) = fsm.state
    (; vegsnowd_xy, rho_snow, snowthickness, Ds_soft, openfrac) = w

    incr = (wnd_dir == 1 || wnd_dir == 3) ? 1 : -1

    @inbounds for dim1 in 1:dim1_size
        for k in 1:index_array[dim1, 1]
            kstart = get_kstart(wnd_dir, dim1, k, index_array)
            kend = get_kend(wnd_dir, dim1, k, index_array)
            for dim2 in kstart:incr:kend

                ic = (wnd_dir == 1 || wnd_dir == 2) ? dim1 : dim2
                jc = (wnd_dir == 1 || wnd_dir == 2) ? dim2 : dim1
                iuw = get_i_upwind(wnd_dir, ic)
                juw = get_j_upwind(wnd_dir, jc)
                idw = get_i_downwind(wnd_dir, ic)
                jdw = get_j_downwind(wnd_dir, jc)

                # Stop at an internal wind reversal: the terminal cell deposits the inflow
                if dim2 == kend && idw >= 1 && idw <= Nx && jdw >= 1 && jdw <= Ny
                    Qs_dir[ic, jc] = zero(Tf)
                end

                dSWE_loc = dt * (Qs_dir[iuw, juw] - Qs_dir[ic, jc]) / Ld[ic, jc]

                if openfrac[ic, jc] >= Tf(0.1)
                    dSWE_loc = dSWE_loc / openfrac[ic, jc]
                else
                    dSWE_loc = zero(Tf)
                end

                if dSWE_loc < eps(Tf)

                    dSWE_loss = -dSWE_loc
                    swe_loc = sum(Sice[kk, ic, jc] + Sliq[kk, ic, jc] for kk in axes(Sice, 1))
                    dSWE_loss = min(dSWE_loss, swe_loc)
                    if dSWE_loss > eps(Tf)
                        dHS_loss = hs_from_swe(fsm, w, dSWE_loss, ic, jc)
                    else
                        dSWE_loss = zero(Tf)
                        dHS_loss = zero(Tf)
                    end
                    Ds_hard = snowthickness[ic, jc] - Ds_soft[ic, jc]
                    snowdmin = max(vegsnowd_xy[ic, jc], Ds_hard)
                    if snowthickness[ic, jc] > snowdmin && fsnow[ic, jc] > eps(Tf)
                        if snowthickness[ic, jc] - dHS_loss / fsnow[ic, jc] <= snowdmin
                            dHS_loss = (snowthickness[ic, jc] - snowdmin) * fsnow[ic, jc]
                            if dHS_loss > eps(Tf)
                                dSWE_loss = swe_from_hs(fsm, w, dHS_loss, ic, jc)
                            else
                                dSWE_loss = zero(Tf)
                            end
                        end
                    else
                        dSWE_loss = zero(Tf)
                    end

                    if dSWE_loss > eps(Tf) && dHS_loss > eps(Tf)
                        snow_ablation!(fsm, dHS_loss, dSWE_loss, ic, jc, Tm)
                        if fsnow[ic, jc] > eps(Tf)
                            Ds_soft[ic, jc] = max(Ds_soft[ic, jc] - dHS_loss / fsnow[ic, jc], zero(Tf))
                        end
                        dSWE_loc = -dSWE_loss
                        Qs_dir[ic, jc] = Qs_dir[iuw, juw] + dSWE_loss * openfrac[ic, jc] * Ld[ic, jc] / dt
                    else
                        dSWE_loc = zero(Tf)
                        Qs_dir[ic, jc] = Qs_dir[iuw, juw]
                    end

                else

                    dHS_gain = dSWE_loc / rho_snow
                    if dSWE_loc > eps(Tf) && dHS_gain > eps(Tf)
                        Sice0[ic, jc] += dSWE_loc
                        snowdepth0[ic, jc] += dHS_gain
                    else
                        dSWE_loc = zero(Tf)
                        Qs_dir[ic, jc] = Qs_dir[iuw, juw]
                    end

                end

                dSWE_s[ic, jc] += dSWE_loc
                snowthickness[ic, jc] = column_sum(Ds, ic, jc)

            end
        end
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Update the SWE and depth changes due to a transport flux (saltation or
suspension) with four directional sweeps.

Translation of the Fortran subroutine `UPDATE_MASS_CHANGES`.
"""
function update_mass_changes!(
        fsm::FSM{Tf}, w::SnowTransport{Tf},
        Qs_u::Matrix{Tf}, Qs_v::Matrix{Tf}, dSWE_s::Matrix{Tf},
        snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf}, Tm::Tf
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; index_ue, index_uw, index_vn, index_vs) = w

    compute_mass_changes_per_direction!(fsm, w, 1, index_uw, Nx, dSWE_s, Qs_u, snowdepth0, Sice0, Tm)
    compute_mass_changes_per_direction!(fsm, w, 2, index_ue, Nx, dSWE_s, Qs_u, snowdepth0, Sice0, Tm)
    compute_mass_changes_per_direction!(fsm, w, 3, index_vs, Ny, dSWE_s, Qs_v, snowdepth0, Sice0, Tm)
    compute_mass_changes_per_direction!(fsm, w, 4, index_vn, Ny, dSWE_s, Qs_v, snowdepth0, Sice0, Tm)

    return nothing
end

"""
$(TYPEDSIGNATURES)

Compute the saltation, suspension and sublimation fluxes and the associated mass
and depth changes, then accumulate the cumulated transport arrays `dSWE_tot_*`.

Translation of the Fortran subroutine `FLUXES_AND_MASS_BALANCE`.
"""
function fluxes_and_mass_balance!(
        fsm::FSM{Tf}, w::SnowTransport{Tf}, met::MET{Tf},
        snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        dSWE_salt::Matrix{Tf}, dSWE_susp::Matrix{Tf}, dSWE_subl::Matrix{Tf},
        flag_blowing_snow::Bool, Tm::Tf
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; dt) = fsm.params
    (; Ds, fsnow, Sice, Sliq) = fsm.state
    (; vegsnowd_xy, openfrac, snowthickness, Ds_soft) = w
    (; Qsalt_u, Qsalt_v, Qsusp_u, Qsusp_v, Qsusp_eq_u, Qsusp_eq_v) = w
    (; Qsubl, Qsubl_susp, Qsubl_salt) = w
    (; dSWE_tot_subl, dSWE_tot_salt, dSWE_tot_susp) = w

    if flag_blowing_snow

        # Saltation
        saltation!(fsm, w)
        update_mass_changes!(fsm, w, Qsalt_u, Qsalt_v, dSWE_salt, snowdepth0, Sice0, Tm)

        # Suspension
        suspension!(fsm, w, met)
        copyto!(Qsusp_eq_u, Qsusp_u)
        copyto!(Qsusp_eq_v, Qsusp_v)
        relax_suspension_flux!(fsm, w)
        update_mass_changes!(fsm, w, Qsusp_u, Qsusp_v, dSWE_susp, snowdepth0, Sice0, Tm)

        # Scale the suspension-layer sublimation by |Qsusp| / |Qsusp_eq|
        @inbounds for j in 1:Ny, i in 1:Nx
            Qsusp_eq = sqrt(Qsusp_eq_u[i, j]^2 + Qsusp_eq_v[i, j]^2)
            Qsusp = sqrt(Qsusp_u[i, j]^2 + Qsusp_v[i, j]^2)
            subl_ratio = Qsusp_eq > eps(Tf) ? Qsusp / Qsusp_eq : zero(Tf)
            Qsubl[i, j] = Qsubl_salt[i, j] + Qsubl_susp[i, j] * subl_ratio
        end

        # Sublimation mass loss
        @inbounds for j in 1:Ny, i in 1:Nx

            Ds_hard = snowthickness[i, j] - Ds_soft[i, j]
            snowdmin = max(vegsnowd_xy[i, j], Ds_hard)
            swe_loc = column_swe(Sice, Sliq, i, j)

            dSWE_subl[i, j] = Qsubl[i, j] * dt
            if openfrac[i, j] >= Tf(0.1)
                dSWE_subl[i, j] = dSWE_subl[i, j] / openfrac[i, j]
            else
                dSWE_subl[i, j] = zero(Tf)
            end
            dSWE_subl_loss = -dSWE_subl[i, j]
            dSWE_subl_loss = min(dSWE_subl_loss, swe_loc)
            dSWE_subl[i, j] = -dSWE_subl_loss
            if dSWE_subl_loss > eps(Tf)
                dHS_subl_loss = hs_from_swe(fsm, w, dSWE_subl_loss, i, j)
            else
                dSWE_subl_loss = zero(Tf)
                dHS_subl_loss = zero(Tf)
                dSWE_subl[i, j] = zero(Tf)
            end

            if snowthickness[i, j] > snowdmin && fsnow[i, j] > eps(Tf)
                if snowthickness[i, j] - dHS_subl_loss / fsnow[i, j] <= snowdmin
                    dHS_subl_loss = (snowthickness[i, j] - snowdmin) * fsnow[i, j]
                    if dHS_subl_loss > eps(Tf)
                        dSWE_subl_loss = swe_from_hs(fsm, w, dHS_subl_loss, i, j)
                        dSWE_subl[i, j] = -dSWE_subl_loss
                    else
                        dHS_subl_loss = zero(Tf)
                        dSWE_subl_loss = zero(Tf)
                        dSWE_subl[i, j] = zero(Tf)
                    end
                end
            else
                dSWE_subl[i, j] = zero(Tf)
                dSWE_subl_loss = zero(Tf)
                dHS_subl_loss = zero(Tf)
            end

            if dSWE_subl_loss > eps(Tf) && dHS_subl_loss > eps(Tf)
                snow_ablation!(fsm, dHS_subl_loss, dSWE_subl_loss, i, j, Tm)
            end

            snowthickness[i, j] = column_sum(Ds, i, j)

        end

    end

    # Update cumulated SWE changes
    @inbounds for j in 1:Ny, i in 1:Nx
        dSWE_tot_salt[i, j] += dSWE_salt[i, j]
        dSWE_tot_susp[i, j] += dSWE_susp[i, j]
        dSWE_tot_subl[i, j] += dSWE_subl[i, j]
    end

    return nothing
end

"""
$(TYPEDSIGNATURES)

Snow transport by wind using Liston's SnowTran3D model.

Pure Julia implementation of the Liston and Sturm (1998) and Liston et
al. (2007) SnowTran3D model, translated from the Fortran routine `SNOWTRAN3D`
(deps/SNOWTRAN3D.F90) and producing the same results as [`snowtran3d!`](@ref)
up to Float32 round-off (the two implementations use different math libraries
for `exp`/`log`/`^` etc., so results are not guaranteed bit-identical).
Reference: Quéno et al. (2024)

Note on orientation (as in the Fortran code): y is the W->E axis, while x is
the S->N axis, i.e. South of (i,j) is (i-1,j) and West of (i,j) is (i,j-1).

# Arguments
- `fsm::FSM`: Model state structure
- `met::MET`: Meteo variable structure
- `snowdepth0::Matrix`: Snow depth of deposited snow (m) - modified in-place
- `Sice0::Matrix`: Ice content of deposited snow (kg/m²) - modified in-place
- `dSWE_salt::Matrix`: SWE change due to saltation (kg/m²) - output
- `dSWE_susp::Matrix`: SWE change due to suspension (kg/m²) - output
- `dSWE_subl::Matrix`: SWE change due to sublimation (kg/m²) - output
"""
function snowtran3d_julia!(
        fsm::FSM{Tf}, met::MET{Tf}, w::SnowTransport{Tf}, snowdepth0::Matrix{Tf}, Sice0::Matrix{Tf},
        dSWE_salt::Matrix{Tf}, dSWE_susp::Matrix{Tf},
        dSWE_subl::Matrix{Tf}
    ) where {Tf <: Real}

    (; Nx, Ny) = fsm.grid
    (; Ds) = fsm.state
    Ua_eff = fsm.diag.Uaeff
    (; Udir) = met
    (; uwind, vwind, snowthickness) = w
    (; Qsalt_u, Qsalt_v, Qsusp_u, Qsusp_v, Qsubl) = w
    @unpack_constants(Tf)

    # Minimum wind speed to compute snow transport (m/s)
    wind_min = Tf(4)

    fill!(Qsalt_u, zero(Tf))
    fill!(Qsalt_v, zero(Tf))
    fill!(dSWE_salt, zero(Tf))
    fill!(Qsusp_u, zero(Tf))
    fill!(Qsusp_v, zero(Tf))
    fill!(dSWE_susp, zero(Tf))
    fill!(Qsubl, zero(Tf))
    fill!(dSWE_subl, zero(Tf))

    # Wind components: y is the W->E axis, x is the S->N axis
    @. uwind = Ua_eff * cos(-Udir * pi / Tf(180) - pi / Tf(2))
    @. vwind = Ua_eff * sin(-Udir * pi / Tf(180) - pi / Tf(2))

    windspd_flag = zero(Tf)
    @inbounds for j in 1:Ny, i in 1:Nx
        windspd_flag = max(windspd_flag, Ua_eff[i, j])
        snowthickness[i, j] = column_sum(Ds, i, j)
    end

    compute_soft_snow!(fsm, w)
    compute_threshold_friction_velocity!(fsm, w)

    flag_blowing_snow = false

    if windspd_flag >= wind_min
        compute_wind_direction_arrays!(w)
        flag_blowing_snow = compute_friction_velocity!(fsm, w)
        update_soft_snow!(fsm, w)
    end

    fluxes_and_mass_balance!(
        fsm, w, met, snowdepth0, Sice0,
        dSWE_salt, dSWE_susp, dSWE_subl, flag_blowing_snow, Tm
    )

    return nothing
end
