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

  subfolder_matlab = "SNFRAC_" + snfrac;

  if tile == "opn"
    subfolder_julia = "OPEN_SNFRAC_" + snfrac;
  elseif tile == "for"
    subfolder_julia = "FOREST_SNFRAC_" + snfrac;
  elseif tile == "glc"
    subfolder_julia = "GLACIER_SNFRAC_" + snfrac;
  end

  mat = load("D:\julia\FSM_HS_matlab\LATEST_00h_RUN\"+ subfolder_matlab + "\OUTPUT_OSHD_0250\RESULTS_24h_" + tile + "\MODELDATA_" + str1 + "-" + str2 + "_FSM22.mat");
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


  % Plot variable 1

  figure("Position",[100 100 1000 700]);

  t = tiledlayout(2,2);
  title(t, var1 + " for " + str2)

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

  linkaxes(ax)


  % Plot variable 2

  figure("Position",[100 100 1000 700]);

  t = tiledlayout(2,2);
  title(t, var2 + " for " + str2)

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

  linkaxes(ax)


  % Location with largest error for var 1

  diff_abs = abs(var1_jul-var1_mat);

  disp("===============================================")
  disp("Variable: " + var1)

  [diff_max, ~] = max(diff_abs,[],"all");

  [row_max, col_max] = find(diff_abs == diff_max);

  disp("Max diff = " + diff_max(1))

  disp("Row max = " + row_max(1))

  disp("Row max unflipped = " + (size(var1_mat,1) - row_max(1) + 1))

  disp("Col max = " + col_max(1))

  % Find linear index of cell with largest error

  domain = load("K:\OSHD_AUX\MODEL_SETTINGS\FSM\OSHD_DOMAIN_FSM_0250.mat");

  linear_index = zeros(domain.grid.nrows,domain.grid.ncols);
  linear_index(domain.grid.data) = 1:sum(domain.grid.data(:));
  linear_index = flipud(linear_index);

  disp("Linear index = " + linear_index(row_max(1),col_max(1)))

end