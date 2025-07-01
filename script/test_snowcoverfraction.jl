using FSMOSHD
using BenchmarkTools
using Dates
using MAT

landuse_file_loc = "D:/jim_operational/SOURCE/BAFU_LUS_0250_2023a.mat"

landuse_file = matopen(landuse_file_loc)
landuse = read(landuse_file, "landuse")
close(landuse_file)

Nx = round(Int, landuse["nrows"])
Ny = round(Int, landuse["ncols"])

fsm = FSM{Float64, Int64}(Nx=Nx, Ny=Ny)
setup_grid!(fsm, landuse)

SWEbuffer = fsm.SWEbuffer
snowdepthbuffer = fsm.snowdepthbuffer
diffSWEbuffer = fsm.diffSWEbuffer

t = DateTime(2023, 10, 1, 0, 0)

function snowcoverfraction(fsm, Nx, Ny, t, SWEbuffer, snowdepthbuffer, diffSWEbuffer)
  for j in 1:Ny, i in 1:Nx
      snowcoverfraction!(fsm, 0.1, 10.0, t, i, j, SWEbuffer, snowdepthbuffer, diffSWEbuffer)
  end
end


fsm.SNFRAC = 3

@time snowcoverfraction(fsm, Nx, Ny, t, SWEbuffer, snowdepthbuffer, diffSWEbuffer)




