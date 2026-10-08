using FlexibleSnowModelOSHD: aspect_perturbation

@testset "aspect_perturbation ($Tf)" for Tf in (Float32, Float64)
    base = Tf(20)
    f(Sdir, Sdird; steepness = Tf(1)) = aspect_perturbation(steepness, base, Tf(Sdir), Tf(Sdird))

    # Edge cases without direct radiation on one or both surfaces
    @test f(0, 0) == one(Tf)
    @test f(0, 300) ≈ 1 / base
    @test f(300, 0) == base

    # Equal radiation on the slope and on flat ground leaves the decay unchanged
    @test f(300, 300) == one(Tf)

    # Sunny slopes decay faster, shaded slopes slower, symmetrically
    @test f(600, 300) > 1
    @test f(300, 600) < 1
    @test f(600, 300) * f(300, 600) ≈ 1

    # Bounded by [1 / base, base], steeper for a larger steepness
    for (Sdir, Sdird) in ((1000, 1), (1, 1000), (500, 100))
        @test 1 / base <= f(Sdir, Sdird) <= base
    end
    @test f(600, 300; steepness = Tf(5)) > f(600, 300)

    @test f(600, 300) isa Tf
end
