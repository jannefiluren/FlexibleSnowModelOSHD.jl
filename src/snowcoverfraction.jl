function snowcoverfraction!(fsm::FSM{Tf, Ti}, snowdepth::Tf, SWEtmp::Tf, t::DateTime, i::Int, j::Int, SWEbuffer::Array{Tf}, snowdepthbuffer::Array{Tf}, diffSWEbuffer::Array{Tf}) where {Tf <: Real, Ti <: Integer}

    @unpack SNFRAC = fsm
    @unpack hfsn = fsm
    @unpack fsnow, swehist, swemin, swemax = fsm
    @unpack snowdepthhist, snowdepthmin, snowdepthmax = fsm
    @unpack slopemu, xi, Ld = fsm

    snowdepth_threshold = Tf(0.005714286)  # Converted from swe_threshold = 2 using density = 350
    
    if SNFRAC == 0

        # Calculate topographic terms for snow depth standard deviation
        sd_snowdepth1 = exp(Tf(-1) / (Ld[i,j]/xi[i,j])^Tf(2))
        sd_snowdepth3 = slopemu[i,j]^Tf(0.309)

        # Merge current values with history
        SWEbuffer[1] = SWEtmp
        SWEbuffer[2:15] .= swehist[:, i, j]
        snowdepthbuffer[1] = snowdepth
        snowdepthbuffer[2:15] .= snowdepthhist[:, i, j]

        # Find indices of extrema
        iabsmax = argmax(SWEbuffer)
        iabsmin = argmin(SWEbuffer)

        # Find recent minimum
        diffSWEbuffer .= diff(SWEbuffer)
        iloop = findfirst(x -> x > Tf(0.5), diffSWEbuffer)
        iloop = isnothing(iloop) ? 14 : iloop
        irecentmin = argmin(SWEbuffer[1:iloop])
        
        # Get minimum, maximum and recent minimum snow depth
        snowdepthmin_buffer = snowdepthbuffer[iabsmin]
        snowdepthmax_buffer = snowdepthbuffer[iabsmax]
        snowdepthmin_recent = snowdepthbuffer[irecentmin]
        
        # Compute difference metrics for snow depth 
        dsnowdepth = snowdepth - snowdepthmin_buffer
        if (dsnowdepth < eps(Tf))
            dsnowdepth = Tf(0)
        end
        
        dsnowdepthmax = snowdepthmax_buffer - snowdepthmin_buffer
        if (dsnowdepthmax < eps(Tf))
            dsnowdepthmax = Tf(0)
        end
        
        if (dsnowdepthmax < dsnowdepth)
            dsnowdepthmax = dsnowdepth
        end

        dsnowdepth_recent = snowdepth - snowdepthmin_recent
        if (dsnowdepth_recent < eps(Tf))
            dsnowdepth_recent = Tf(0)
        end
                
        # Handle state variables
        if SWEtmp < eps(Tf)
            swemax[i,j] = swemin[i,j] = Tf(0.0)
        end
        if snowdepth < eps(Tf)
            snowdepthmax[i,j] = snowdepthmin[i,j] = Tf(0.0)
        end
        
        if SWEtmp >= swemax[i,j]
            swemax[i,j] = swemin[i,j] = SWEtmp
        end
        if snowdepth >= snowdepthmax[i,j]
            snowdepthmax[i,j] = snowdepthmin[i,j] = snowdepth
        end
        
        if SWEtmp < swemax[i,j] && SWEtmp < swemin[i,j]
            swemin[i,j] = SWEtmp
        end
        if snowdepth < snowdepthmax[i,j] && snowdepth < snowdepthmin[i,j]
            snowdepthmin[i,j] = snowdepth
        end
        
        # Calculate snow covered fraction components
        fsnow_season = fsnow_nsnow = fsnow_nsnow_recent = Tf(0.0)
        
        # Seasonal SCF calculation
        sd_snowdepth2 = snowdepthmax[i,j]^Tf(0.549)
        sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3
        
        if !(slopemu[i,j] > eps(Tf))
            sd_snowdepth0 = snowdepthmax[i,j]^Tf(0.84)
        end
        
        if snowdepthmax[i,j] > eps(Tf)
            fsnow_season = tanh(Tf(1.3) * snowdepthmin[i,j] / sd_snowdepth0)
        end
        
        # Calculate coefficient of variation
        coeff_vari = sd_snowdepth0 / snowdepthmax[i,j]
        
        # New snow SCF calculation
        sd_snowdepth0_dhs = dsnowdepthmax^Tf(0.84)
        if dsnowdepthmax > eps(Tf)
            fsnow_nsnow = tanh(dsnowdepth^Tf(0.14) + dsnowdepth/Tf(0.13))
        end
        
        # Recent new snow SCF calculation
        sd_snowdepth0_dhs_recent = dsnowdepth_recent^Tf(0.84)
        if dsnowdepth_recent > eps(Tf)
            fsnow_nsnow_recent = tanh(dsnowdepth_recent^Tf(0.14) + dsnowdepth_recent/Tf(0.13))
        end
        
        # Final SCF
        fsnow[i,j] = max(fsnow_season, fsnow_nsnow, fsnow_nsnow_recent)
        
        # Update history at 6am
        if 4.5 < hour(t) < 5.5
            swehist[:,i,j] .= SWEbuffer[1:14]
            snowdepthhist[:,i,j] .= snowdepthbuffer[1:14]
        end
        
        fsnow[i,j] = max(fsnow[i,j], Tf(0.01))
        
        
        if false
            println("iabsmax: ", iabsmax)
            println("iabsmin: ", iabsmin)
            println("SWEbuffer: ", SWEbuffer)
            println("irecentmin: ", irecentmin)
            println("snowdepthmin_buffer: ", snowdepthmin_buffer)
            println("snowdepthmax_buffer: ", snowdepthmax_buffer)
            println("snowdepthmin_recent: ", snowdepthmin_recent)        
            println("dsnowdepth: ", dsnowdepth)
            println("dsnowdepthmax: ", dsnowdepthmax)
            println("dsnowdepth_recent: ", dsnowdepth_recent) 
            println("swemin: ", swemin)
            println("swemax: ", swemax)
            println("snowdepthmin: ", snowdepthmin)
            println("snowdepthmax: ", snowdepthmax)
            println("sd_snowdepth0 :", sd_snowdepth0)
            println("sd_snowdepth2 :", sd_snowdepth2)
            println("fsnow_season :", fsnow_season)
            println("sd_snowdepth0_dhs: ", sd_snowdepth0_dhs)
            println("fsnow_nsnow: ", fsnow_nsnow)
            println("sd_snowdepth0_dhs_recent: ", sd_snowdepth0_dhs_recent)
            println("fsnow_nsnow_recent: ", fsnow_nsnow_recent)
        end
        
        
    elseif SNFRAC == 1
        # HelbigHS
        sd_snowdepth2 = snowdepth^Tf(0.549)
        sd_snowdepth1 = exp(Tf(1) / (Ld[i,j]/xi[i,j])^Tf(2))
        sd_snowdepth3 = slopemu[i,j]^Tf(0.309)
        sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

        fsnow[i,j] = tanh(Tf(1.3) * snowdepth / sd_snowdepth0)
        
        if snowdepth < snowdepth_threshold
            fsnow[i,j] = Tf(0.0)
        end

    elseif SNFRAC == 2
        # HelbigHS0
        if snowdepth == Tf(0)
            snowdepthmax[i,j] = Tf(0.0)
        end
        
        if snowdepth > snowdepthmax[i,j]
            snowdepthmax[i,j] = snowdepth
        end

        sd_snowdepth2 = snowdepthmax[i,j]^Tf(0.549)
        sd_snowdepth1 = exp(Tf(1) / (Ld[i,j]/xi[i,j])^Tf(2))
        sd_snowdepth3 = slopemu[i,j]^Tf(0.309)
        sd_snowdepth0 = sd_snowdepth1 * sd_snowdepth2 * sd_snowdepth3

        fsnow[i,j] = tanh(Tf(1.3) * snowdepth / sd_snowdepth0)
        
        if snowdepth < snowdepth_threshold
            fsnow[i,j] = Tf(0.0)
        end

    elseif SNFRAC == 3
        # Point model
        fsnow[i,j] = snowdepth > eps(Tf) ? Tf(1.0) : Tf(0.0)

    else
        # tanh model / original FSM
        fsnow[i,j] = tanh(snowdepth/hfsn)
    end

    # Final adjustments
    if snowdepth < eps(Tf)
        fsnow[i,j] = Tf(0.0)
    else
        fsnow[i,j] = min(fsnow[i,j], Tf(1.0))
    end
end