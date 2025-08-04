%% Load data

clear; clc

time = datenum(2025,1,11,6,00,00);

str1 = datestr(time-1,"yyyymmddHHMM");
str2 = datestr(time,"yyyymmddHHMM");

subfolder = "SNFRAC_0";

mat = load("D:\julia\FSM_HS_all\LATEST_00h_RUN\"+ subfolder + "\OUTPUT_OSHD_0250\RESULTS_24h_opn\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat");


subfolder = "SNFRAC_0";

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

fig = figure("Position",[100 100 1000 700]);

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

fig = figure("Position",[100 100 1000 700]);

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


%% Location with largest error

variable = "fsnow";

if variable == "hs"
  diff_abs = abs(hs_jul-hs_mat);
else
  diff_abs = abs(fsnow_jul-fsnow_mat);
end

[diff_max, ~] = max(diff_abs,[],"all");

[row_max, col_max] = find(diff_abs == diff_max);

disp("Max diff = " + diff_max)

disp("Row max = " + row_max)

disp("Row max unflipped = " + (size(hs_mat,1) - row_max + 1))

disp("Col max = " + col_max)

% Find linear index of cell with largest error

domain = load("K:\OSHD_AUX\MODEL_SETTINGS\FSM\OSHD_DOMAIN_FSM_0250.mat");

linear_index = zeros(domain.grid.nrows,domain.grid.ncols);
linear_index(domain.grid.data) = 1:sum(domain.grid.data(:));
linear_index = flipud(linear_index);

disp("Linear index = " + linear_index(row_max,col_max))
