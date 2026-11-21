using FlexibleSnowModelOSHD
using Test

# One snow column compacted for a step by each CrocusCompaction mode
function compacted_column(scheme; Uaeff = 12.0f0)
    fsm = bare_fsm(Grid(Float32; Nx = 1, Ny = 1); compaction = scheme)
    state = fsm.state
    state.Nsnow[1, 1] = 3
    state.fsnow[1, 1] = 1.0f0
    for (k, (ds, rho)) in enumerate(((0.02f0, 120.0f0), (0.04f0, 200.0f0), (0.3f0, 300.0f0)))
        state.Ds[k, 1, 1] = ds
        state.Sice[k, 1, 1] = rho * ds
        state.Sliq[k, 1, 1] = k == 3 ? 2.0f0 : 0.0f0
        state.Tsnow[k, 1, 1] = 265.0f0
    end
    fsm.diag.Uaeff[1, 1] = Uaeff
    fsm.surface.z0_snow[1, 1] = 0.001f0
    FlexibleSnowModelOSHD.compact_snow!(scheme, 1, 1, state, fsm.diag, fsm.surface, fsm.params)
    return state.Ds[1:3, 1, 1]
end

@testset "CrocusCompaction options" begin

    Ds_off = compacted_column(CrocusCompaction{Float32}())
    Ds_cubic = compacted_column(CrocusCompaction{Float32}(wndens = 1))
    Ds_linear = compacted_column(CrocusCompaction{Float32}(wndens = 2))

    # Wind compaction only densifies the layers within d_max of the surface
    @test all(Ds_cubic[1:2] .< Ds_off[1:2])
    @test all(Ds_linear[1:2] .< Ds_off[1:2])
    @test Ds_cubic[3] ≈ Ds_off[3]
    @test Ds_linear[3] ≈ Ds_off[3]

    # Below the threshold wind speed, wind compaction does nothing
    @test compacted_column(CrocusCompaction{Float32}(wndens = 1); Uaeff = 2.0f0) ≈
        compacted_column(CrocusCompaction{Float32}(); Uaeff = 2.0f0)

    # A smaller liquid-water factor softens wet snow less: the wet layer compacts less
    Ds_fliq = compacted_column(CrocusCompaction{Float32}(f_liq = 60))
    @test Ds_fliq[3] > Ds_off[3]
    @test Ds_fliq[1:2] == Ds_off[1:2]

end
