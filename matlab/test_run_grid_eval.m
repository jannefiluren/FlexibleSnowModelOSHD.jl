function test_run_grid_eval(time,tile,snfrac,var1,var2)

arguments
  time = datenum(2025,03,24,6,00,00)
  tile = "opn"
  snfrac = 0
  var1 = "snowdepth"
  var2 = "fsnow"
end

% Load data

str1 = datestr(time-1,"yyyymmddHHMM");
str2 = datestr(time,"yyyymmddHHMM");

if tile == "opn"
  subfolder_julia = "OPEN_SNFRAC_" + snfrac;
elseif tile == "for"
  subfolder_julia = "FOREST_SNFRAC_" + snfrac;
elseif tile == "glc"
  subfolder_julia = "GLACIER_SNFRAC_" + snfrac;
end

mat = load("D:\julia\FSM_HS\LATEST_00h_RUN\OUTPUT_OSHD_0250\RESULTS_24h_" + tile + "\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat");
jul = load("D:\julia\FSM_HS_julia\"+ subfolder_julia + "\" + str2 + "_output.mat");


% Prepare data

var1_jul = jul.(var1).data;
acro = jul.(var1).acronymn;
var1_mat = mat.(acro).data;

var2_jul = jul.(var2).data;
acro = jul.(var2).acronymn;
var2_mat = mat.(acro).data;

inan = var1_mat<0;

var1_mat(inan) = NaN;
var1_jul(inan) = NaN;

var2_mat(inan) = NaN;
var2_jul(inan) = NaN;

if true
  var1_mat = flipud(var1_mat);
  var1_jul = flipud(var1_jul);
  var2_mat = flipud(var2_mat);
  var2_jul = flipud(var2_jul);
end

% Get errors for variable 1 and variable 2

diff_var1 = find_max_error(var1, var1_jul, var1_mat);
diff_var2 = find_max_error(var2, var2_jul, var2_mat);


% Plot variable 1

figure("Position",[100 100 800 600]);

t = tiledlayout(2,2);
title(t, var1 + " for " + str2 + " (tile=" + tile + " , snfrac=" + snfrac + ")")

ax(1) = nexttile();
imagesc(var1_mat)
colorbar()
title("FSM Fortran/Matlab")
colormap('turbo');

ax(2) = nexttile();
imagesc(var1_jul)
colorbar()
title("FSM Julia")
colormap('turbo');

ax(3) = nexttile();
imagesc(var1_jul-var1_mat)
colorbar()
title("FSM Julia minus FSM Fortran/Matlab")
colormap('turbo');

nexttile()
plot(var1_mat(:),var1_jul(:),'.')
xlabel("Matlab/Fortran")
ylabel("Julia")
title("Max diff = " + diff_var1)

linkaxes(ax)


% Plot variable 2

figure("Position",[1000 100 800 600]);

t = tiledlayout(2,2);
title(t, var2 + " for " + str2 + " (tile=" + tile + " , snfrac=" + snfrac + ")")

ax(1) = nexttile();
imagesc(var2_mat)
colorbar()
title("FSM Fortran/Matlab")
colormap('turbo');

ax(2) = nexttile();
imagesc(var2_jul)
colorbar()
title("FSM Julia")
colormap('turbo');

ax(3) = nexttile();
imagesc(var2_jul-var2_mat)
colorbar()
title("FSM Julia minus FSM Fortran/Matlab")
colormap('turbo');

nexttile()
plot(var2_mat(:),var2_jul(:),'.')
xlabel("Matlab/Fortran")
ylabel("Julia")
title("Max diff = " + diff_var2)

linkaxes(ax)

end


function diff_max = find_max_error(var, var_jul, var_mat)

diff_abs = abs(var_jul-var_mat);

disp("===============================================")
disp("Variable: " + var)

[diff_max, ~] = max(diff_abs,[],"all");

[row_max, col_max] = find(diff_abs == diff_max);

disp("Max diff = " + diff_max(1))

disp("Row max = " + row_max(1))

disp("Row max unflipped = " + (size(var_mat,1) - row_max(1) + 1))

disp("Col max = " + col_max(1))

end