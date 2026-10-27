using FlexibleSnowModelOSHD: snow_albedo!

@testset "aspect_decay_ratio ($Tf)" for Tf in (Float32, Float64)
    grid = Grid(Tf, Nx = 1, Ny = 1)
    fsm = bare_fsm(grid)
    (; state, surface, params) = fsm
    met = MET{Tf}(Nx = 1, Ny = 1)
    met.Sf .= zero(Tf)
    met.Sf24h .= zero(Tf)
    state.Sice[1, 1, 1] = Tf(100)  # above swe_thin, so the thin-snow cap does not clamp the albedo

    # One hour of cold-snow decay from 0.8, with and without the ratio
    function albedo_after(scheme, Sdir, Sdird)
        met.Sdir .= Tf(Sdir)
        met.Sdird .= Tf(Sdird)
        state.albs .= Tf(0.8)
        state.Tsrf .= Tf(260)
        snow_albedo!(scheme, 1, 1, state, surface, met, params, false)
        return state.albs[1, 1]
    end
    plain = PrognosticAlbedo{Tf}(grid)
    ratio = PrognosticAlbedo{Tf}(grid; aspect_decay_ratio = true)

    # Halved direct radiation on flat ground halves the decay time, doubling the decrement
    @test Tf(0.8) - albedo_after(ratio, 600, 300) ≈ 2 * (Tf(0.8) - albedo_after(plain, 600, 300))

    # No effect where flat ground gets as much or more direct radiation, or none at all
    @test albedo_after(ratio, 300, 600) == albedo_after(plain, 300, 600)
    @test albedo_after(ratio, 300, 300) == albedo_after(plain, 300, 300)
    @test albedo_after(ratio, 0, 0) == albedo_after(plain, 0, 0)
end
