using FlexibleSnowModelOSHD: aspect_scaling, snow_albedo!

@testset "aspect_scaling ($Tf)" for Tf in (Float32, Float64)
    max = Tf(20)
    f(Sdir, Sdird; steepness = Tf(1), sunny = max, shaded = max) = aspect_scaling(steepness, sunny, shaded, Tf(Sdir), Tf(Sdird))

    # Edge cases without direct radiation on one or both surfaces
    @test f(0, 0) == one(Tf)
    @test f(0, 300) ≈ 1 / max
    @test f(300, 0) == max

    # Equal radiation on the slope and on flat ground leaves the decay unchanged
    @test f(300, 300) == one(Tf)

    # Sunny slopes decay faster, shaded slopes slower, symmetrically for equal bounds
    @test f(600, 300) > 1
    @test f(300, 600) < 1
    @test f(600, 300) * f(300, 600) ≈ 1

    # Bounded by [1 / max, max], steeper for a larger steepness
    for (Sdir, Sdird) in ((1000, 1), (1, 1000), (500, 100))
        @test 1 / max <= f(Sdir, Sdird) <= max
    end
    @test f(600, 300; steepness = Tf(5)) > f(600, 300)

    @test f(600, 300) isa Tf
end

@testset "asymmetric aspect_scaling ($Tf)" for Tf in (Float32, Float64)
    sunny, shaded = Tf(20), Tf(5)
    f(Sdir, Sdird) = aspect_scaling(one(Tf), sunny, shaded, Tf(Sdir), Tf(Sdird))
    x = atan(log(Tf(2))) / (Tf(pi) / 2)

    # Each side takes its own bound, with the same exponent magnitude
    @test f(600, 300) ≈ sunny^x
    @test f(300, 600) ≈ shaded^(-x)
    @test !(f(600, 300) * f(300, 600) ≈ 1)

    # Bounds and edge cases follow the side
    @test f(300, 0) == sunny
    @test f(0, 300) ≈ 1 / shaded
    @test f(300, 300) == one(Tf)
    for (Sdir, Sdird) in ((1000, 1), (1, 1000), (500, 100))
        @test 1 / shaded <= f(Sdir, Sdird) <= sunny
    end

    # A bound of 1 switches that side off
    @test aspect_scaling(one(Tf), sunny, one(Tf), Tf(300), Tf(600)) == one(Tf)
    @test aspect_scaling(one(Tf), one(Tf), shaded, Tf(600), Tf(300)) == one(Tf)
end

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
