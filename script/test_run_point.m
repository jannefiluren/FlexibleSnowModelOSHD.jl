%% Initialize a new run

run_id = 0;
run_id_hn = -1;
time_start = datenum(2024,9,1);
geom = "point";
init_type = "initialize";
silent = 0;
hourly_output = "full_research"
ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"init_type",init_type,"silent",silent,"hourly_output",hourly_output);


%% Restart an existing run

run_id = 0;
run_id_hn = -1;
time_start = datenum(2025,3,1);
geom = "point";
init_type = "reinitialize";
states_folder_reinit = "D:\MODEL_DATA_FSM\FSM_HS\LATEST_00h_RUN";
silent = 0;
hourly_output = "full_research"
ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"init_type",init_type,"silent",silent,"hourly_output",hourly_output,"states_folder_reinit",states_folder_reinit);