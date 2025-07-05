# alb = similar(albs)   ### hack
# asrf_out = similar(albs)   ### hack
# Sdirt = similar(albs)   ### hack
# Sdift = similar(albs)   ### hack
# SWveg = similar(albs)   ### hack
# SWsrf = similar(albs)   ### hack
# SWsci = similar(albs)   ### hack
# LWt = similar(albs)   ### hack
# SWtopo_out = similar(albs)   ### hack

function radiation(fsm::FSM, meteo::MET, t)

  @unpack Nx, Ny, dt = fsm
  
  @unpack tthresh, fsky_terr, fveg, dem, tilefrac = fsm

  @unpack asmx, asmn, avg0, avgs, Talb, tcld, tmlt, adfs, adfl, fsar, Sfmin = fsm

  @unpack alb0,fsky,scap, trcn = fsm

  @unpack albs, Sice, Sliq, fsnow,Sveg, Tsnow, Tsrf = fsm

  @unpack ALBEDO, OSHDTN, RADSBG, CANMOD, ALRADT, ALPERT = fsm

  @unpack alb, asrf_out, Sdirt, Sdift, SWveg, SWsrf, SWsci, LWt, SWtopo_out = fsm

  @unpack adm, adc, afs = fsm

  @unpack LW, Sdif, Sdir, Sdird, Sf, Sf24h, Ta, Tv = meteo
  
  # Snow albedo
  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        adc_loc = adc[i,j]
        adm_loc = adm
        afs_loc = afs[i,j]

        if (ALBEDO == 0)
          # Diagnostic
          albs[i, j] = asmn + (afs_loc - asmn) * (Tsrf[i, j] - Tm) / Talb
          if (albs[i, j] < min(afs_loc, asmn))
            albs[i, j] = min(afs_loc, asmn)
          end
          if (albs[i, j] > max(afs_loc, asmn))
            albs[i, j] = max(afs_loc, asmn)
          end

        elseif (ALBEDO == 1)
          # Prognostic
          tau = tcld
          if (Tsrf[i, j] >= Tm)
            tau = tmlt
          end
          # Forest adjustments -> not yet properly tested for OSHD but option currently unused
          if (Dates.value(Month(t)) > 4 && Dates.value(Month(t)) < 10)
            tau = 70.0 * 3600.0
          end

          if fveg[i, j] > 0 && Sdir[i, j] > eps(Float64)
            tau = tau / ((1 - trcn[i, j] * fsky[i, j]) * (1 + adfl * Tv[i, j]) + adfs * Tv[i, j])
          elseif fveg[i, j] > 0 && Sdif[i, j] > eps(Float64)
            tau = tau / ((1 - trcn[i, j] * fsky[i, j]) + adfs * trcn[i, j] * fsky[i, j])
          elseif (fveg[i, j] > 0 && (Sdir[i, j] + Sdif[i, j] <= eps(Float64)))
            tau = tau / (2.0 - trcn[i, j] * fsky[i, j])
          end

          rt = 1 / tau + Sf[i, j] / Sfmin
          alim = (asmn / tau + Sf[i, j] * afs_loc / Sfmin) / rt
          albs[i, j] = alim + (albs[i, j] - alim) * exp(-rt * dt)
          if (albs[i, j] < min(afs_loc, asmn))
            albs[i, j] = min(afs_loc, asmn)
          end
          if (albs[i, j] > max(afs_loc, asmn))
            albs[i, j] = max(afs_loc, asmn)
          end

        else # ALBEDO == 2
          # Prognostic, tuned, copied from JIM
          SWEtmp = 0.0
          for si in 1:size(Sice, 1)
            SWEtmp += Sice[si, i, j] + Sliq[si, i, j]
          end

          # BC 08.23: aspect-dependent albedo tuning. Activated for oper season 2024
          # or optionally.
          if ((ALRADT == 1) || (OSHDTN == 1))
            # BC Oct 23: Jan's suggestion: modify only when the decay rate should be increased (ad* DECREASE), not decreased
            if ((Sdir[i,j] > eps(Float64)) && (Sdird[i,j] < Sdir[i,j]))
              adm_loc = adm_loc * (Sdird[i,j])/(Sdir[i,j])
              adc_loc = adc_loc * (Sdird[i,j])/(Sdir[i,j])
              if (adm_loc < eps(Float64))
                adm_loc = eps(Float64)
              end
              if (adc_loc < eps(Float64))
                adc_loc = eps(Float64)
              end
            end
          end
          
          # if (ALPERT)
          #   adm_loc = adm_loc * alP[i,j]
          #   adc_loc = adc_loc * alP[i,j]
          # end
          if (Tsrf[i, j] >= Tm)
            albs[i, j] = (albs[i, j] - asmn) * exp(-(dt / 3600) / adm_loc) + asmn
          else
            albs[i, j] = albs[i, j] - (dt / 3600) / adc_loc
          end
          if (SWEtmp < 75.0) # more stuff showing on and up through snow
            afs_loc *= 0.80
          end
          # Reset to fresh snow albedo (wasn't originally available; only else term)
          if ((Sf[i, j] * dt) > 0.0 && Sf24h[i, j] > Sfmin)
            albs[i, j] = afs_loc
          else
            albs[i, j] = albs[i, j] + (afs_loc - albs[i, j]) * Sf[i, j] * dt / Sfmin
          end
          ## End Adjustments
          if (albs[i, j] > afs_loc)
            albs[i, j] = afs_loc
          end
          if (albs[i, j] < asmn)
            albs[i, j] = asmn
          end

        end
      end
    end
  end

  # Surface and canopy net shortwave radiation
  for j = 1:Ny
    for i = 1:Nx
      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        # Surface albedo
        asrf = albs[i, j] * (1 - fveg[i, j] * fsar)
        if (fsnow[i, j] <= eps(Float64))
          asrf = alb0[i, j]
          albs[i,j] = alb0[i,j]
        end

        # Partial snowcover on canopy
        fcans = 0.0
        if (scap[i, j] > eps(Float64))
          fcans = Sveg[i, j] / scap[i, j]
        end
        aveg = (1 - fcans) * avg0 + fcans * avgs
        acan = fveg[i, j] * aveg
        # Canopy surface albedo for computing terrain radiation over canopy
        alb[i, j] = fveg[i, j] * aveg + (1 - fveg[i, j]) * asrf

        # Surface albedo is stored in asurf_out to write in results
        asrf_out[i, j] = alb[i, j]

        #   hack not used
        #   if (RADSBG == 1)
        #     # Call Subgrid parameterization for SW radiation to compute SWtopo,netto SWtn
        #     call SWRADTOPO(alb[i,j],Sdir[i,j],Sdif[i,j],SWsrf[i,j],Sdirt[i,j],Sdift[i,j],SWtopo_out,Sun_elev,Dates.value(Year(t)),Dates.value(Month(t)),Dates.value(Day(t)),Dates.value(Hour(t)),i,j)
        #   end

        if (RADSBG == 0)
          #global SWtopo_out  ### hack
          SWtopo_out[i,j] = alb[i, j] * (Sdir[i, j] + Sdif[i, j])
          Sdirt[i, j] = Sdir[i, j]
          Sdift[i, j] = Sdif[i, j]
        end

        # Solar radiation trasmission 
        if (CANMOD == 0)
          SWveg[i, j] = 0
          SWsrf[i, j] = (1 - alb[i, j]) * (Sdir[i, j] + Sdif[i, j])
          SWsci[i, j] = Sdift[i, j] + Sdir[i, j]
        end

        if (CANMOD == 1)
          Sdif_aux = fsky[i, j] / fsky_terr[i, j] * Sdift[i, j]
          tdif = trcn[i, j]
          tdir = Tv[i, j]

          # Effective albedo and net radiation
          alb[i, j] = acan + (1 - acan) * asrf * tdif^2
          if (Sdif_aux + Sdirt[i, j] > eps(Float64))
            alb[i, j] = (acan * (Sdif_aux + tdir * Sdirt[i, j]) + asrf * tdif * (tdif * Sdif_aux + tdir * Sdirt[i, j])) / (Sdif_aux + Sdirt[i, j])
          end
          SWsrf[i, j] = (1 - asrf) * (tdif * Sdif_aux + tdir * Sdirt[i, j])
          SWveg[i, j] = ((1 - tdif) * (1 - aveg) + tdif * asrf * (1 - tdif)) * Sdif_aux + (tdir * fveg[i, j] * (1 - aveg) + tdir * asrf * (1 - tdif)) * Sdirt[i, j]   # local SWR absorption by vegetation correlates with local tdir  
          SWsci[i, j] = tdif * Sdif_aux + tdir * Sdirt[i, j]
        end

        # Thermal emissions from surroundings 
        # Terrain LWR if not calculated later; 
        LWt[i, j] = fsky_terr[i, j] * LW[i, j] + (1 - fsky_terr[i, j]) * sb * Ta[i, j]^4

        # LW overwritten by LWt only if EBALFOR is used, where terrain impacts are accounted for already 
        if (CANMOD == 0 || fveg[i, j] == 0)
          LW[i, j] = LWt[i, j]
        end

      end

    end

  end

end