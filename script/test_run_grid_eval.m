%% Load data

clear; clc

time = datenum(2024,10,1,6,00,00);

str1 = datestr(time-1,"yyyymmddHHMM");
str2 = datestr(time,"yyyymmddHHMM");

subfolder = "SNFRAC_3";

mat = load("D:\julia\FSM_HS_all\LATEST_00h_RUN\"+ subfolder + "\OUTPUT_OSHD_0250\RESULTS_24h_opn\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat");


subfolder = "SNFRAC_3";

jul = load("D:\julia\FSM_HS_julia\"+ subfolder + "\" + str2 + "_output.mat");


%% Prepare data

hs_mat = mat.hsnt.data;
hs_jul = jul.hs;

fsnow_mat = mat.scfe.data;
fsnow_jul = jul.fsnow;

inan = hs_mat<0;

hs_mat(inan) = NaN;
hs_jul(inan) = NaN;

fsnow_mat(inan) = NaN;
fsnow_jul(inan) = NaN;

if true
  hs_mat = flipud(hs_mat);
  hs_jul = flipud(hs_jul);
  fsnow_mat = flipud(fsnow_mat);
  fsnow_jul = flipud(fsnow_jul);
end


%% Plot snow depth results

fig = figure("Position",[294 151 1330 886]);

t = tiledlayout(2,2);
title(t, "Snow depth for " + str2)

ax(1) = nexttile();
imagesc(hs_mat)
colorbar()
title("FSM Fortran/Matlab")
colormap('turbo');

ax(2) = nexttile();
imagesc(hs_jul)
colorbar()
title("FSM Julia")
colormap('turbo');

ax(3) = nexttile();
imagesc(hs_jul-hs_mat)
colorbar()
title("FSM Julia minus FSM Fortran/Matlab")
colormap('turbo');

nexttile()
plot(hs_mat(:),hs_jul(:),'.')
xlabel("Matlab/Fortran")
ylabel("Julia")

linkaxes(ax)


%% Plot fsnow results

fig = figure("Position",[294 151 1330 886]);

t = tiledlayout(2,2);
title(t, "Snow cover fraction for " + str2)

ax(1) = nexttile();
imagesc(fsnow_mat)
colorbar()
title("FSM Fortran/Matlab")
colormap('turbo');

ax(2) = nexttile();
imagesc(fsnow_jul)
colorbar()
title("FSM Julia")
colormap('turbo');

ax(3) = nexttile();
imagesc(fsnow_jul-fsnow_mat)
colorbar()
title("FSM Julia minus FSM Fortran/Matlab")
colormap('turbo');

nexttile()
plot(fsnow_mat(:),fsnow_jul(:),'.')
xlabel("Matlab/Fortran")
ylabel("Julia")

linkaxes(ax)


% %% Location with largest error

% diff_abs = abs(hs_jul(~inan)-hs_mat(~inan));

% [diff_max, ind_max] = max(diff_abs,[],"all");

% disp(diff_max)

% disp(ind_max)


% %% Plot fsnow against hs

% figure; plot(hs_jul(:), fsnow_jul(:), '.'); title("julia")
% figure; plot(hs_mat(:), fsnow_mat(:), '.'); title("matlab")