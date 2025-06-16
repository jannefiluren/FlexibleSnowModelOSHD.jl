%% Initialize a new run

run_id = 0;
run_id_hn = -1;
time_start = datenum(2024,9,1);
root_folder = "D:\julia";
geom = "point";
init_type = "initialize";
silent = 0;
hourly_output = "full_research"
sel_stat = {'SLF.5WJ'};
ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"init_type",init_type,"silent",silent,"hourly_output",hourly_output,"sel_stat",sel_stat,"root_folder",root_folder);


%% Restart an existing run

run_id = 0;
run_id_hn = -1;
time_start = datenum(2025,3,1);
root_folder = "D:\julia";
geom = "point";
init_type = "reinitialize";
states_folder_reinit = "D:\julia\FSM_HS\LATEST_00h_RUN";
silent = 0;
hourly_output = "full_research"
sel_stat = {'SLF.5WJ'};
run_folder = "D:\julia\FSM_HS\bin_files";
stop_run = true;
ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"init_type",init_type,"silent",silent,"hourly_output",hourly_output,"states_folder_reinit",states_folder_reinit,"sel_stat",sel_stat,"root_folder",root_folder,"run_folder",run_folder,"stop_run",stop_run);