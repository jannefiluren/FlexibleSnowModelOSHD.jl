run_id = 0;
run_id_hn = -1;
time_start = datenum(2024,9,1);
geom = "point";
init_type = "initialize";
silent = 0;
hourly_output = "full_research"
ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"init_type",init_type,"silent",silent,"hourly_output",hourly_output);