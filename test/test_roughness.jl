using FlexibleSnowModelOSHD: ground_roughness, wet_snow_roughness, Parameters

@testset "wet-snow roughness ($Tf)" for Tf in (Float32, Float64)
    grid = Grid(Tf, Nx = 1, Ny = 1)
    fsm = bare_fsm(grid)
    (; state, surface) = fsm
    surface.z0_snow .= Tf(0.1)
    surface.z0sf .= Tf(0.5)
    state.fsnow .= Tf(0.5)
    state.Ds[1, 1, 1] = Tf(0.2)
    params = Parameters{Tf}(; c_z0_wet = Tf(0.25), theta_z0_sat = Tf(0.04))

    # Liquid water giving a volumetric content theta in the top layer
    set_theta!(theta) = (state.Sliq[1, 1, 1] = Tf(theta) * Tf(1000) * Tf(0.2) * Tf(0.5))
    z0(p = params) = wet_snow_roughness(1, 1, state, surface, p)

    set_theta!(0)
    @test z0() == Tf(0.1)

    set_theta!(0.02)
    @test z0() ≈ Tf(0.1) * (1 + (Tf(0.25) - 1) * Tf(0.5))

    set_theta!(0.04)
    @test z0() ≈ Tf(0.025)
    set_theta!(0.1)
    @test z0() ≈ Tf(0.025)
    @test z0() isa Tf

    # The defaults leave the roughness unscaled, bit for bit
    @test z0(Parameters{Tf}()) === Tf(0.1)

    # No snow in the top layer: unscaled
    state.Ds[1, 1, 1] = zero(Tf)
    @test z0() == Tf(0.1)
    state.Ds[1, 1, 1] = Tf(0.2)

    # ground_roughness applies the scaling only on snow-covered ground
    set_theta!(0.04)
    @test ground_roughness(SeasonalSnowFraction{Tf}(), 1, 1, state, surface, params) ≈ Tf(0.025)
    state.fsnow .= zero(Tf)
    @test ground_roughness(SeasonalSnowFraction{Tf}(), 1, 1, state, surface, params) == Tf(0.5)
end
