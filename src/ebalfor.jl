function ebalfor(fsm::FSM{Tf,Ti}, meteo::MET{Tf,Ti}) where {Tf<:Real,Ti<:Integer}

  @unpack tthresh = fsm

  @unpack LW, Ps, Qa, Ta = meteo

  @unpack Nx, Ny, dt = fsm

  @unpack canh, fsky, trcn = fsm

  @unpack Qcan, Sice, Tcan, Tsrf, Tveg = fsm

  @unpack fveg, tilefrac = fsm

  @unpack Ds1, KHa, KHg, KHv, KWg, KWv, ks1, SWsrf, SWveg, Ts1, Tveg0 = fsm

  @unpack Esrf, Eveg, G, H, Hsrf, LE, LEsrf, LWsci, LWveg, Melt, Rnet, Rsrf = fsm

  A = zeros(Tf, 4, 4)    # TODO remove these for allocations
  b = zeros(Tf, 4)    # TODO remove these for allocations
  x = zeros(Tf, 4)    # TODO remove these for allocations

  # 1-layer canopy model
  for j = 1:Ny
    for i = 1:Nx

      if (tilefrac[i, j] >= tthresh) # exclude points outside tile of interest

        if (fveg[i, j] > eps(Tf))
          # Saturation humidity and density of air
          Qsrf = qsat(Ps[i, j], Tsrf[i, j])
          Lsrf = Ls
          if (Tsrf[i, j] > Tm)
            Lsrf = Lv
          end
          Dsrf = Lsrf * Qsrf / (Rwat * Tsrf[i, j]^Tf(2))
          Qveg = qsat(Ps[i, j], Tveg[i, j])
          Lveg = Ls
          if (Tveg[i, j] > Tm)
            Lveg = Lv
          end
          Dveg = Lveg * Qveg / (Rwat * Tveg[i, j]^Tf(2))
          rho = Ps[i, j] / (Rair * Ta[i, j])

          # Explicit fluxes
          E = rho * KHa[i, j] * (Qcan[i, j] - Qa[i, j])
          Esrf[i, j] = rho * KWg[i, j] * (Qsrf - Qcan[i, j])
          Eveg[i, j] = rho * KWv[i, j] * (Qveg - Qcan[i, j])
          G[i, j] = Tf(2) * ks1[i, j] * (Tsrf[i, j] - Ts1[i, j]) / Ds1[i, j]
          H[i, j] = rho * cp * KHa[i, j] * (Tcan[i, j] - Ta[i, j])
          Hsrf[i, j] = rho * cp * KHg[i, j] * (Tsrf[i, j] - Tcan[i, j])
          Hveg = rho * cp * KHv[i, j] * (Tveg[i, j] - Tcan[i, j])
          LE[i, j] = Lsrf * Esrf[i, j] + Lveg * Eveg[i, j]
          Melt[i, j] = Tf(0)
          Rsrf[i, j] = SWsrf[i, j] + trcn[i, j] * (fsky[i, j] * LW[i, j] + (Tf(1) - fsky[i, j]) * sb * Ta[i, j]^Tf(4)) - sb * Tsrf[i, j]^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)        # with near and distant canopy contributions                
          # Rsrf[i,j] = SWsrf[i,j] + trcn[i,j]*LW[i,j] - sb*Tsrf[i,j]^4 + (1 - trcn[i,j])*sb*Tveg[i,j]^4  ! original formulations
          Rveg = SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LW[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))

          # Surface energy balance increments without melt
          A[1, 1] = Tf(0)
          A[1, 2] = -(KHa[i, j] + KHv[i, j] + KHg[i, j])
          A[1, 3] = KHg[i, j]
          A[1, 4] = KHv[i, j]
          b[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
          A[2, 1] = -(KHa[i, j] + KWv[i, j] + KWg[i, j])
          A[2, 2] = Tf(0)
          A[2, 3] = Dsrf * KWg[i, j]
          A[2, 4] = Dveg * KWv[i, j]
          b[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
          A[3, 1] = -Lsrf * rho * KWg[i, j]
          A[3, 2] = -rho * cp * KHg[i, j]
          A[3, 3] = rho * (cp * KHg[i, j] + Lsrf * Dsrf * KWg[i, j]) + Tf(4) * sb * Tsrf[i, j]^Tf(3) + Tf(2) * ks1[i, j] / Ds1[i, j]
          A[3, 4] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
          b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
          A[4, 1] = -Lveg * rho * KWv[i, j]
          A[4, 2] = -rho * cp * KHv[i, j]
          A[4, 3] = -Tf(4) * (Tf(1) - trcn[i, j]) * sb * Tsrf[i, j]^Tf(3)
          A[4, 4] = canh[i, j] / dt + rho * (cp * KHv[i, j] + Lveg * Dveg * KWv[i, j]) + Tf(8) * (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(3)
          b[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
          ludcmp!(4, A, b, x)
          dQc = x[1]
          dTc = x[2]
          dTs = x[3]
          dTv = x[4]
          dEs = rho * KWg[i, j] * (Dsrf * dTs - dQc)
          dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
          dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
          dHs = rho * cp * KHg[i, j] * (dTs - dTc)
          dHv = rho * cp * KHv[i, j] * (dTv - dTc)

          # Surface melting
          if (Tsrf[i, j] + dTs > Tm && Sice[1, i, j] > eps(Sice[1, i, j]))
            Melt[i, j] = sum(Sice[:, i, j]) / dt
            b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j] - Lf * Melt[i, j]
            ludcmp!(4, A, b, x)
            dQc = x[1]
            dTc = x[2]
            dTs = x[3]
            dTv = x[4]
            dEs = rho * KWg[i, j] * (Dsrf * dTs - dQc)
            dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
            dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
            dHs = rho * cp * KHg[i, j] * (dTs - dTc)
            dHv = rho * cp * KHv[i, j] * (dTv - dTc)
            if (Tsrf[i, j] + dTs < Tm)
              Qsrf = qsat(Ps[i, j], Tm)
              Esrf[i, j] = rho * KWg[i, j] * (Qsrf - Qcan[i, j])
              G[i, j] = Tf(2) * ks1[i, j] * (Tm - Ts1[i, j]) / Ds1[i, j]
              Hsrf[i, j] = rho * cp * KHg[i, j] * (Tm - Tcan[i, j])
              Rsrf[i, j] = SWsrf[i, j] + trcn[i, j] * (fsky[i, j] * LW[i, j] + (Tf(1) - fsky[i, j]) * sb * Ta[i, j]^Tf(4)) - sb * Tm^Tf(4) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
              Rveg = SWveg[i, j] + (Tf(1) - trcn[i, j]) * (LW[i, j] + sb * Tm^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))
              A[1, 3] = Tf(0)
              b[1] = (H[i, j] - Hveg - Hsrf[i, j]) / (rho * cp)
              A[2, 3] = Tf(0)
              b[2] = (E - Eveg[i, j] - Esrf[i, j]) / rho
              A[3, 3] = Tf(1)
              b[3] = Rsrf[i, j] - Hsrf[i, j] - Lsrf * Esrf[i, j] - G[i, j]
              A[4, 3] = Tf(0)
              b[4] = Rveg - Hveg - Lveg * Eveg[i, j] - canh[i, j] * (Tveg[i, j] - Tveg0[i, j]) / dt
              ludcmp!(4, A, b, x)
              dQc = x[1]
              dTc = x[2]
              Melt[i, j] = x[3] / Lf
              dTv = x[4]
              dTs = Tm - Tsrf[i, j]
              dEs = Tf(0)
              dEv = rho * KWv[i, j] * (Dveg * dTv - dQc)
              dGs = Tf(2) * ks1[i, j] * dTs / Ds1[i, j]
              dHs = Tf(0)
              dHv = rho * cp * KHv[i, j] * (dTv - dTc)
            end
          end

          # Update temperatures and fluxes
          Qcan[i, j] = Qcan[i, j] + dQc
          Tcan[i, j] = Tcan[i, j] + dTc
          Tsrf[i, j] = Tsrf[i, j] + dTs
          Tveg[i, j] = Tveg[i, j] + dTv
          Esrf[i, j] = Esrf[i, j] + dEs
          Eveg[i, j] = Eveg[i, j] + dEv
          Hsrf[i, j] = Hsrf[i, j] + dHs
          Hveg = Hveg + dHv
          G[i, j] = G[i, j] + dGs
          H[i, j] = Hsrf[i, j] + Hveg
          LE[i, j] = Lsrf * Esrf[i, j] + Lveg * Eveg[i, j]
          LEsrf[i, j] = Lsrf * Esrf[i, j]
          Rnet[i, j] = SWsrf[i, j] + SWveg[i, j] + LW[i, j] - trcn[i, j] * sb * Tsrf[i, j]^Tf(4) - (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
          LWsci[i, j] = trcn[i, j] * (fsky[i, j] * LW[i, j] + (Tf(1) - fsky[i, j]) * sb * Ta[i, j]^Tf(4)) + (Tf(1) - trcn[i, j]) * sb * Tveg[i, j]^Tf(4)
          LWveg[i, j] = (Tf(1) - trcn[i, j]) * (LW[i, j] + sb * Tsrf[i, j]^Tf(4) - Tf(2) * sb * Tveg[i, j]^Tf(4))

          # Sublimation limited by amount of snow after melt
          Ssub = sum(Sice[:, i, j]) - Melt[i, j] * dt
          if (Ssub > eps(Tf) && Esrf[i, j] * dt > Ssub)
            Esrf[i, j] = Ssub / dt
            LEsrf[i, j] = Ls * Esrf[i, j]
            Hsrf[i, j] = Rnet[i, j] - G[i, j] - LEsrf[i, j] - Lf * Melt[i, j]
          end
        end

      end

    end
  end

end