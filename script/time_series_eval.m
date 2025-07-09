%% Setup

clear; clc

tstart = datenum(2024,9,2,6,0,0);
tend = datenum(2025,6,12,6,0,0);

subfolder_matlab = "SNFRAC_0";
subfolder_julia = "SNFRAC_0_20250709";

folder_matlab = "D:\julia\FSM_HS_all\LATEST_00h_RUN\"+ subfolder_matlab + "\OUTPUT_OSHD_0250\RESULTS_24h_opn";
folder_julia = "D:\julia\FSM_HS_julia\"+ subfolder_julia;


%% Load data

times = tstart:tend;

for itime = 1:length(times)
  
  str1 = datestr(times(itime)-1,"yyyymmddHHMM");
  str2 = datestr(times(itime),"yyyymmddHHMM");

  disp("Loading data for " + str1 + " to " + str2)

  mat = load(fullfile(folder_matlab, "MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat"), "hsnt", "scfe");
  jul = load(fullfile(folder_julia, str2 + "_output.mat"));

  hs_mat(:,:,itime) = flipud(mat.hsnt.data);
  hs_jul(:,:,itime) = flipud(jul.hs);
  
  fsnow_mat(:,:,itime) = flipud(mat.scfe.data);
  fsnow_jul(:,:,itime) = flipud(jul.fsnow);

end

inan = hs_mat<0;

hs_mat(inan) = NaN;
hs_jul(inan) = NaN;

fsnow_mat(inan) = NaN;
fsnow_jul(inan) = NaN;



%% Plot dummy grid ensuring correct selection of pixel

% irow = 915;
% icol = 648;

irow = 561;
icol = 1251;

figure
check_grid = zeros(size(hs_mat,1), size(hs_mat,2));
check_grid(any(isnan(hs_mat),3)) = NaN;
check_grid(irow, icol) = 10;

h = imagesc(check_grid);
set(h, 'AlphaData', ~isnan(check_grid))


%% Plot time series

figure

t = tiledlayout(2,1);

ax(1) = nexttile();
plot(times, squeeze(hs_mat(irow,icol,:)), "b")
hold on
plot(times, squeeze(hs_jul(irow,icol,:)), "--r")
legend("Matlab", "Julia")
datetick("x")

ax(2) = nexttile();
plot(times, squeeze(fsnow_mat(irow,icol,:)), "b")
hold on
plot(times, squeeze(fsnow_jul(irow,icol,:)), "--r")
datetick("x")

linkaxes(ax)








