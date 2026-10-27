using FlexibleSnowModelOSHD: ground_roughness, snow_roughness, solar_radiation!

@testset "snow roughness ($Tf)" for Tf in (Float32, Float64)
    grid = Grid(Tf, Nx = 1, Ny = 1)
    fsm = bare_fsm(grid)
    (; state, surface) = fsm
    surface.z0_snow .= Tf(0.1)
    surface.z0sf .= Tf(0.5)
    state.fsnow .= Tf(0.5)
    state.Ds[1, 1, 1] = Tf(0.2)
    fixed = FixedSnowRoughness{Tf}()
    wet = WetSnowRoughness{Tf}(; z0_wet_ratio = Tf(0.25), theta_z0_wet = Tf(0.04))

    # Liquid water giving a volumetric content theta in the top layer
    set_theta!(theta) = (state.Sliq[1, 1, 1] = Tf(theta) * Tf(1000) * Tf(0.2) * Tf(0.5))
    z0(scheme = wet) = snow_roughness(scheme, 1, 1, state, surface)

    set_theta!(0)
    @test z0() == Tf(0.1)

    set_theta!(0.02)
    @test z0() ≈ Tf(0.1) * (1 + (Tf(0.25) - 1) * Tf(0.5))

    set_theta!(0.04)
    @test z0() ≈ Tf(0.025)
    set_theta!(0.1)
    @test z0() ≈ Tf(0.025)
    @test z0() isa Tf

    # The fixed scheme ignores the liquid water; a ratio of 1 reproduces it bit for bit
    @test z0(fixed) === Tf(0.1)
    @test z0(WetSnowRoughness{Tf}()) === z0(fixed)

    # No snow in the top layer: unscaled
    state.Ds[1, 1, 1] = zero(Tf)
    @test z0() == Tf(0.1)
    state.Ds[1, 1, 1] = Tf(0.2)

    # ground_roughness applies the snow roughness only on snow-covered ground
    set_theta!(0.04)
    @test ground_roughness(SeasonalSnowFraction{Tf}(), wet, 1, 1, state, surface) ≈ Tf(0.025)
    state.fsnow .= zero(Tf)
    @test ground_roughness(SeasonalSnowFraction{Tf}(), wet, 1, 1, state, surface) == Tf(0.5)
end

@testset "open-cover shortwave per horizontal area ($Tf)" for Tf in (Float32, Float64)
    grid = Grid(Tf, Nx = 1, Ny = 1)
    fsm = bare_fsm(grid)
    (; state, diag, surface) = fsm
    met = MET{Tf}(Nx = 1, Ny = 1)
    state.albs .= Tf(0.6)
    met.Sdir .= Tf(500)
    met.Sdird .= Tf(300)
    met.Sdif .= Tf(100)

    solar_radiation!(OpenCover{Tf}(), 1, 1, state, diag, surface, met)
    @test diag.SWsci[1, 1] == Tf(600)
    @test diag.SWsrf[1, 1] ≈ Tf(0.4) * Tf(600)

    solar_radiation!(OpenCover{Tf}(; horizontal_sw = true), 1, 1, state, diag, surface, met)
    @test diag.SWsci[1, 1] == Tf(400)
    @test diag.SWsrf[1, 1] ≈ Tf(0.4) * Tf(400)
end
