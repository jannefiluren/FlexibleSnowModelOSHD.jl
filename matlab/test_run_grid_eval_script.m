tile = "opn";
snfrac = 3;
time_start = datenum(2024,9,12,6,00,00);
for i = 0:9
  time = addtodate(time_start, i, 'month');
  test_run_grid_eval(time,tile,snfrac,"snowdepth","fsnow");
   pause()
   close all
end