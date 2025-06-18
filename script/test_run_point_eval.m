clear; close all; clc

% Read data

hs_matlab = readtable("D:\julia\FSM_HS_all\snowdepth_matlab.csv");
hs_julia = readtable("D:\julia\FSM_HS_all\snowdepth_julia.csv");
hs_matlab = table2array(hs_matlab);
hs_julia = table2array(hs_julia);

% Scatter plot

figure()
plot(hs_matlab(:),hs_julia(:),'.')
hold on
plot(0:4,0:4,'r')

% Histogram

figure()
histogram(hs_julia(:)-hs_matlab(:),100)

% Fraction of data below error threshold

disp("Fraction of data with error less than 0.01m: " + sum(abs(hs_julia(:)-hs_matlab(:))<0.01) / length(hs_julia(:)))

% Get station with largets error

landuse = load("K:/OSHD_AUX/DATA_LUS/OSHD_LUS_STAT.mat");

[~, imax] = max(max(abs((hs_julia - hs_matlab)),[],1));

disp("Station with largest error: " + landuse.acro{imax})
