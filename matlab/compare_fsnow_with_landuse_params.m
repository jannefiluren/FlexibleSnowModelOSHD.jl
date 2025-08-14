%% Settings

clear; clc

time = datenum(2024,12,09,06,00,00);


% OI = off, prec_multi = on
% folder = "D:\MODEL_DATA_FSM\FSM_HS_exper_1\LATEST_00h_RUN\OUTPUT_OSHD_0250\RESULTS_24h_opn";

% OI = off, prec_multi = off
% folder = "D:\julia\FSM_HS_all\LATEST_00h_RUN\SNFRAC_0\OUTPUT_OSHD_0250\RESULTS_24h_opn";
folder = "D:\MODEL_DATA_FSM\FSM_HS_exper_2\LATEST_00h_RUN\OUTPUT_OSHD_0250\RESULTS_24h_opn";

% OI = on, prec_multi = on
% folder = "K:\MODEL_DATA_FSM\FSM_HS\LATEST_00h_RUN\OUTPUT_OSHD_0250\RESULTS_24h_opn";


%% Load landuse

landuse = load("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_0250.mat");

dhdxdy = landuse.dhdxdy.data;
sddem = landuse.sd.data;

slopemu = sqrt((dhdxdy./2));
xi = (sqrt(2)*sddem)./slopemu;

landuse.slopemu.data = slopemu;
landuse.xi.data = xi;


%% Load snow simulation

str1 = datestr(time-1,"yyyymmddHHMM");
str2 = datestr(time,"yyyymmddHHMM");

mat = load(fullfile(folder, "MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat"));
mat.scfe.data(mat.scfe.data==mat.NODATA_value) = NaN;


%% Plot results

variable = "xi";

lus_data = landuse.(variable).data;
lus_data(lus_data == landuse.NODATA_value) = NaN;

fig = figure("Position",[57 229 1815 657]);

t = tiledlayout(1, 2);

title(t, datestr(time, "yyyymmdd HH:MM"))

ax(1) = nexttile();
imagesc(mat.scfe.data)
title("fsnow")

ax(2) = nexttile();
imagesc(lus_data)
title(variable)

linkaxes(ax)
