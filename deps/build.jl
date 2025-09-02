using Libdl

# Build standalone SnowSlide routines
snowslide_routines = [
    "SNOWSLIDE.F90",
    "SNOW_ABLATION.F90", 
    "SWE_FROM_HS.F90"
]

if Sys.iswindows()
    flags = ["-m$(Sys.WORD_SIZE)", "-shared", "-O3"]  # TODO check that this works for everyone...
else
    flags = ["-m$(Sys.WORD_SIZE)", "-shared", "-O3", "-fPIC"]
end

run(`gfortran $flags $snowslide_routines -o libsnowslide.$(Libdl.dlext)`)

# Clean up .mod files
for f in readdir()
    if endswith(f, ".mod")
        rm(f)
    end
end

println("Built libsnowslide.$(Libdl.dlext) successfully")