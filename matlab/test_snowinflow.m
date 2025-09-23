%% Load data

clear; close all; clc

time = datenum(2022,8,10,06,00,00);

time_str1 = datestr(time, "yyyymmddHHMM");
time_str2 = datestr(time+1, "yyyymmddHHMM");

[~,local_folder] = user_paths();

ref = load(fullfile(local_folder, "\model_data\fsm_data\FSM_output\FSM_HS\FSM2oshd Cal v1 WY 2021 scaling iteration 0\OUTPUT_ULLA_FORRE_0900\RESULTS_24h_opn\MODELDATA_" + time_str1 + "-" + time_str2 + "_FSM22.mat"));
tst = load(fullfile(local_folder, "\model_data\fsm_data\FSM_julia\" + time_str2 + "_output.mat"));

figure
imagesc(tst.snowdepth.data)
colorbar()
title("HS julia")

figure
imagesc(ref.hsnt.data)
colorbar()
title("HS matlab/fortran")

figure
plot(ref.hsnt.data(:), round(tst.snowdepth.data(:), 3), '.')
hold on
refline(1,0)
xlabel("HS matlab/fortran")
ylabel("HS julia")
title(time_str2)