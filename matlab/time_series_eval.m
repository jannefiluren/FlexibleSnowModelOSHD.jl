%% Setup

clear; clc

irow = 489;   % unflipped
icol = 883;

tstart = datenum(2024,9,2,6,0,0);
tend = datenum(2025,6,1,6,0,0);

tile = "opn";
snfrac = 0;

subfolder_matlab = "SNFRAC_" + snfrac;

if tile == "opn"
  subfolder_julia = "OPEN_SNFRAC_" + snfrac;
elseif tile == "for"
  subfolder_julia = "FOREST_SNFRAC_" + snfrac;
elseif tile == "glc"
  subfolder_julia = "GLACIER_SNFRAC_" + snfrac;
end


%% Load data

times = tstart:tend;

hs_mat = -9999*ones(1088, 1488, length(times));
hs_jul = -9999*ones(1088, 1488, length(times));
fsnow_mat = -9999*ones(1088, 1488, length(times));
fsnow_jul = -9999*ones(1088, 1488, length(times));

for itime = 1:length(times)
  
  str1 = datestr(times(itime)-1,"yyyymmddHHMM");
  str2 = datestr(times(itime),"yyyymmddHHMM");

  disp("Loading data for " + str1 + " to " + str2)

  mat = load("D:\julia\FSM_HS_all\LATEST_00h_RUN\"+ subfolder_matlab + "\OUTPUT_OSHD_0250\RESULTS_24h_" + tile + "\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat","hsnt","scfe");
  jul = load("D:\julia\FSM_HS_julia\"+ subfolder_julia + "\" + str2 + "_output.mat");

  hs_mat(:,:,itime) = mat.hsnt.data;
  hs_jul(:,:,itime) = jul.hs;
  
  fsnow_mat(:,:,itime) = mat.scfe.data;
  fsnow_jul(:,:,itime) = jul.fsnow;

end

inan = hs_mat<0;

hs_mat(inan) = NaN;
hs_jul(inan) = NaN;

fsnow_mat(inan) = NaN;
fsnow_jul(inan) = NaN;


%% Plot dummy grid ensuring correct selection of pixel

figure("Position",[100 100 800 500]);

check_grid = zeros(size(hs_mat,1), size(hs_mat,2));
check_grid(any(isnan(hs_mat),3)) = NaN;
check_grid(irow, icol) = 10;

h = imagesc(check_grid);
set(h, 'AlphaData', ~isnan(check_grid))


%% Plot time series

figure("Position",[100 100 800 500]);

t = tiledlayout(2,1);
title(t, "Row = " + irow + " Col = " + icol)

ax(1) = nexttile();
plot(times, squeeze(hs_mat(irow,icol,:)), "b", "LineWidth", 2)
hold on
plot(times, squeeze(hs_jul(irow,icol,:)), "--r", "LineWidth", 2)
legend("Matlab/Fortran", "Julia")
datetick("x")
ylabel("HS (m)")

ax(2) = nexttile();
plot(times, squeeze(fsnow_mat(irow,icol,:)), "b", "LineWidth", 2)
hold on
plot(times, squeeze(fsnow_jul(irow,icol,:)), "--r", "LineWidth", 2)
datetick("x")
ylabel("SCF (-)")

linkaxes(ax, "x")


%% Plot time series with maximum differences

diff_hs = abs(hs_jul - hs_mat);
diff_fsnow = abs(fsnow_jul - fsnow_mat);

figure("Position",[100 100 800 500]);

t = tiledlayout(2,1);
title(t, "Maxiumu difference for whole grid: FSM Julia minus FSM Fortran/Matlab")

ax(1) = nexttile();
plot(times, squeeze(max(diff_hs,[],[1 2])), "b", "LineWidth", 2)
datetick("x")
ylabel("dHS (m)")

ax(2) = nexttile();
plot(times, squeeze(max(diff_fsnow,[],[1 2])), "b", "LineWidth", 2)
datetick("x")
ylabel("dfsnow (-)")
