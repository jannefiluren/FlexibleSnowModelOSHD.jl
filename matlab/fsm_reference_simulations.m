%% Setup

clear; clc;

script_path = mfilename('fullpath');
script_dir = fileparts(script_path);

run_id = 0;
run_id_hn = -1;
time_start = datenum(2024,9,1);
time_end = datenum(2025,6,15);
init_type = "initialize";
silent = 0;
soil_init_file = "";
root_folder = fullfile(script_dir, "../..");

%% Point run

geom = "point";
prec_input_folder = "";
hourly_output = "full_research";

ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"time_end",time_end,"init_type",init_type,"silent",silent,"soil_init_file",soil_init_file,"prec_input_folder",prec_input_folder,"hourly_output",hourly_output,"root_folder",root_folder);

%% Hires run

geom = "hires";
prec_input_folder = "W:/DATA_OI/OI_ICON_SNF";

ansmsg = start_oshd_fsm(run_id,run_id_hn,time_start,geom,"time_end",time_end,"init_type",init_type,"silent",silent,"soil_init_file",soil_init_file,"prec_input_folder",prec_input_folder,"root_folder",root_folder);
