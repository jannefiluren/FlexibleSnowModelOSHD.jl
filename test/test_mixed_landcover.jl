# A single mixed-land-cover run (open / forest / glacier cells in one FSM) must be bit-identical,
# per cell, to running the corresponding uniform single-tile FSM on the same landuse.

using FlexibleSnowModelOSHD
using Test
using Dates

# Wrapped in a `let` so nothing leaks into the shared test module (Main).
let
    Tf = Float32
    Nx, Ny = 1, 6
    d(x) = Dict("data" => x)

    # Cells: 1,4 = open ; 2,5 = forest ; 3,6 = glacier
    forestcells = [false true false false true false]
    icecells = [false false true false false true]
    landuse = Dict(
        "skyvf" => d(fill(0.95, Nx, Ny)),
        "elevation" => d(fill(2540.0, Nx, Ny)),
        "prec_multi" => d(fill(1.1, Nx, Ny)),
        "slopemu" => d(fill(0.3, Nx, Ny)),
        "xi" => d(fill(150.0, Nx, Ny)),
        "Ld" => d(fill(250.0, Nx, Ny)),
        "forest" => d(Float64[0 0.8 0 0 0.8 0]),
        "glacier" => d(Float64[0 0 0.7 0 0 0.7]),
        "fveg" => d(fill(0.6, Nx, Ny)),
        "fves" => d(fill(0.6, Nx, Ny)),
        "hcan" => d(fill(20.0, Nx, Ny)),
        "lai" => d(fill(2.5, Nx, Ny)),
        "vfhp" => d(fill(0.5, Nx, Ny)),
    )

    grid = Grid(Tf, landuse)

    fsm_mixed = FSM(
        grid, landuse;
        land_cover = MixedLandCover(OpenCover{Tf}(), ForestCover{Tf}(), forestcells),
        substrate = MixedSubstrate(SoilSubstrate{Tf}(), IceSubstrate{Tf}(), icecells),
    )
    fsm_open = FSM(grid, landuse)
    fsm_forest = FSM(grid, landuse; land_cover = ForestCover{Tf}())
    fsm_glacier = FSM(grid, landuse; substrate = IceSubstrate{Tf}())

    fillmet!(met) = begin
        met.Sdir .= 300; met.Sdird .= 300; met.Sdif .= 80; met.LW .= 260
        met.Sf .= 1.0f-4; met.Rf .= 0; met.Sf24h .= 5; met.Ta .= 268
        met.RH .= 85; met.Ua .= 3; met.Ps .= 8.0f4; met.Tv .= 0.6
        met
    end
    for fsm in (fsm_mixed, fsm_open, fsm_forest, fsm_glacier)
        met = fillmet!(MET{Tf}(Nx = Nx, Ny = Ny))
        for _ in 1:3
            step!(fsm, met, DateTime(2023, 12, 1, 12))
        end
    end

    # Every array field of state/diag must match the oracle at cell (i, j).
    cell_matches(a, b, i, j) = begin
        ok = true
        for obj in (:state, :diag)
            A = getfield(a, obj); B = getfield(b, obj)
            for f in fieldnames(typeof(A))
                fa = getfield(A, f); fb = getfield(B, f)
                fa isa AbstractArray || continue
                slice_a = ndims(fa) == 2 ? fa[i, j] : fa[:, i, j]
                slice_b = ndims(fb) == 2 ? fb[i, j] : fb[:, i, j]
                isequal(slice_a, slice_b) || (ok = false)
            end
        end
        ok
    end

    oracle(j) = forestcells[1, j] ? fsm_forest : icecells[1, j] ? fsm_glacier : fsm_open

    @testset "bit-identical to per-tile runs" begin
        for j in 1:Ny
            @test cell_matches(fsm_mixed, oracle(j), 1, j)
        end
    end

    @testset "construction guards" begin
        # A half-mixed request (only one Mixed wrapper) is rejected.
        @test_throws ErrorException FSM(
            grid, landuse;
            land_cover = MixedLandCover(OpenCover{Tf}(), ForestCover{Tf}(), forestcells)
        )
        # forestcells / icecells must be disjoint.
        @test_throws ErrorException FSM(
            grid, landuse;
            land_cover = MixedLandCover(OpenCover{Tf}(), ForestCover{Tf}(), [false true true false true false]),
            substrate = MixedSubstrate(SoilSubstrate{Tf}(), IceSubstrate{Tf}(), icecells)
        )
    end
end
