%% Load data

clear; close all; clc

time = datenum(2022,8,10,06,00,00);

time_str1 = datestr(time, "yyyymmddHHMM");
time_str2 = datestr(time+1, "yyyymmddHHMM");

ref = load("D:\snowinflow_project\snowinflow_data\model_data\fsm_data\FSM_output\FSM_HS\OUTPUT_ULLA_FORRE_0900\RESULTS_24h_opn\MODELDATA_" + time_str1 + "-" + time_str2 + "_FSM22.mat");
tst = load("D:\snowinflow_project\snowinflow_data\model_data\fsm_data\FSM_julia\" + time_str2 + "_output.mat");

figure
imagesc(tst.hs)
colorbar()
title("HS julia")

figure
imagesc(ref.hsnt.data)
colorbar()
title("HS matlab/fortran")

figure
plot(ref.hsnt.data(:), round(tst.hs(:), 3), '.')
hold on
refline(1,0)
xlabel("HS matlab/fortran")
ylabel("HS julia")
title(time_str2)