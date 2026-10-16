using FlexibleSnowModelOSHD: aspect_perturbation

@testset "aspect_perturbation ($Tf)" for Tf in (Float32, Float64)
    base = Tf(20)
    f(Sdir, Sdird; steepness = Tf(1), sunny = base, shaded = base) = aspect_perturbation(steepness, sunny, shaded, Tf(Sdir), Tf(Sdird))

    # Edge cases without direct radiation on one or both surfaces
    @test f(0, 0) == one(Tf)
    @test f(0, 300) ≈ 1 / base
    @test f(300, 0) == base

    # Equal radiation on the slope and on flat ground leaves the decay unchanged
    @test f(300, 300) == one(Tf)

    # Sunny slopes decay faster, shaded slopes slower, symmetrically for equal bases
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

@testset "asymmetric aspect_perturbation ($Tf)" for Tf in (Float32, Float64)
    sunny, shaded = Tf(20), Tf(5)
    f(Sdir, Sdird) = aspect_perturbation(one(Tf), sunny, shaded, Tf(Sdir), Tf(Sdird))
    x = atan(log(Tf(2))) / (Tf(pi) / 2)

    # Each side takes its own base, with the same exponent magnitude
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

    # A base of 1 switches that side off
    @test aspect_perturbation(one(Tf), sunny, one(Tf), Tf(300), Tf(600)) == one(Tf)
    @test aspect_perturbation(one(Tf), one(Tf), shaded, Tf(600), Tf(300)) == one(Tf)
end

@testset "PrognosticAlbedo validates the aspect bases" begin
    grid = Grid(Float64, Nx = 2, Ny = 2)
    albedo = PrognosticAlbedo{Float64}(grid; aspect_pert_base_sunny = 10, aspect_pert_base_shaded = 1)
    @test albedo.aspect_pert_base_sunny == 10
    @test albedo.aspect_pert_base_shaded == 1
    @test_throws ArgumentError PrognosticAlbedo{Float64}(grid; aspect_pert_base_sunny = 0.5)
    @test_throws ArgumentError PrognosticAlbedo{Float64}(grid; aspect_pert_base_shaded = 0.9)
end
