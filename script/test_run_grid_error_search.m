%% Load data

clear; clc

ind_max = [493007 493007+1];

times = datenum(2024,9,2,06,00,00):datenum(2025,6,1,06,00,00);

hs_mat_vec = nan(2,length(times));
hs_jul_vec = nan(2,length(times));

for i = 1:length(times)

  time = times(i);

  disp(datestr(time, "yyyymmddHHMM"))

  str1 = datestr(time-1,"yyyymmddHHMM");
  str2 = datestr(time,"yyyymmddHHMM");

  mat = load("D:\julia\FSM_HS_all\LATEST_00h_RUN\OUTPUT_OSHD_0250\RESULTS_24h_opn\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat");
  jul = load("D:\julia\FSM_HS_julia\" + str2 + "_output.mat");

  hs_mat_vec(:,i) = mat.hsnt.data(ind_max);
  hs_jul_vec(:,i) = jul.hs(ind_max);

end


%% Plot results

figure
plot(times, hs_mat_vec, "b")
hold on
plot(times, hs_jul_vec, "r:")
legend("Matlab/Fortran", "Julia")
datetick("x")
