# FSMOSHD

The code located here `D:\julia\FSMOSHD` should be a replicate of `D:\julia\jim_operational` under branch `jan-sync-jim`. The script `test_against_jim.jl` evaluates the translated code agaist the Matlab and Fortran version. All old code is placed under the branch `legacy`.

# Replace number is vs code

Regex: [+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)

Replacement: T($1)

# Discussion on types

https://discourse.julialang.org/t/how-to-code-a-type-flexible-model-that-is-type-stable/124138
