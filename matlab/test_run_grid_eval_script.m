%% Setup

clear; clc;

tile = ["opn", "glc", "for"];
snfrac = [0, 0, 4];
time_start = datenum(2024,9,12,6,00,00);


%% Plot results

for settingindx = 1:length(tile)
  
  tilecurr = tile(settingindx);
  snfraccurr = snfrac(settingindx);
  
  for monthadd = 0:9
    time = addtodate(time_start, monthadd, 'month');
    test_run_grid_eval(time,tilecurr,snfraccurr,"snowdepth","fsnow");
     pause()
     close all
  end

end