# Compile fsm in jim_operational

curr_dir = pwd()

cd(joinpath(base_folder, "jim_operational/FSM_SOURCE_CODE/code"))
run(`compil_FSM.bat`)
cd(curr_dir)

# Move executable to bin folder

mv(joinpath(base_folder, "jim_operational/FSM_SOURCE_CODE/FSM2.exe"), joinpath(base_folder, "FSM_HS/bin_files/FSM2.exe"), force=true)

# Run

cd(joinpath(base_folder, "FSM_HS/bin_files"))
run(`FSM2.exe OPTIONS.nam`)
cd(curr_dir)