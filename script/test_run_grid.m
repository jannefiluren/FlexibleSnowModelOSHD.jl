%% Initialize a new run

run_id = 0;
run_id_hn = -1;
time_start = datenum(2024,9,1);
root_folder = "D:\julia";
geom = "hires";
init_type = "initialize";
silent = 0;
soil_init_file = "";
ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"init_type",init_type,"silent",silent,"root_folder",root_folder,"soil_init_file",soil_init_file);
