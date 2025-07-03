%% Load data

clear; clc

time = datenum(2024,09,5,06,00,00);

str1 = datestr(time-1,"yyyymmddHHMM");
str2 = datestr(time,"yyyymmddHHMM");

subfolder = "SNFRAC_3";

variable = "hs";

mat = load("D:\julia\FSM_HS_all\LATEST_00h_RUN\"+ subfolder + "\OUTPUT_OSHD_0250\RESULTS_24h_opn\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat");
jul = load("D:\julia\FSM_HS_julia\"+ subfolder + "\" + str2 + "_output.mat");


%% Fix data

if variable == "hs"
  hs_mat = mat.hsnt.data;
  hs_jul = jul.hs;
elseif "fsnow"
  hs_mat = mat.scfe.data;
  hs_jul = jul.fsnow;
end

inan = hs_mat<0;

hs_mat(inan) = NaN;
hs_jul(inan) = NaN;


%% Plot results

fig = figure("Position",[294 151 1330 886]);

t = tiledlayout(2,2);

title(t,str2);

ax(1) = nexttile();
imagesc(flipud(hs_mat))
colorbar()
title("FSM Fortran/Matlab")
colormap('turbo');

ax(2) = nexttile();
imagesc(flipud(hs_jul))
colorbar()
title("FSM Julia")
colormap('turbo');

ax(3) = nexttile();
imagesc(flipud(hs_jul-hs_mat))
colorbar()
title("FSM Julia minus FSM Fortran/Matlab")
colormap('turbo');

nexttile()
plot(hs_mat(:),hs_jul(:),'.')
xlabel("Matlab/Fortran")
ylabel("Julia")

linkaxes(ax)


%% Location with largest error

diff_abs = abs(hs_jul-hs_mat);

[diff_max, ind_max] = max(diff_abs,[],"all");