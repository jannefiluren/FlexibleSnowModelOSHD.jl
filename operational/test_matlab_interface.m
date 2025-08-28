%% Test 1: Basic functionality test with short time range
disp("=== Test 1: Basic functionality test ===")

clear; clc

start_time = datenum(2024, 9, 1, 6, 0, 0);
end_time = datenum(2024, 9, 2, 6, 0, 0);

disp("Testing with short time range:")
disp("Start: " + datestr(start_time));
disp("End: " + datestr(end_time));

[success1, output1] = run_fsmoshd_operational(start_time, end_time, false);

disp("success1: " + success1)
disp("output1: " + output1)


%% Test 2: Test restart functionality  
disp("=== Test 2: Restart functionality test ===")

clear; clc

start_time = datenum(2014, 9, 2, 6, 0, 0);
end_time = datenum(2024, 9, 3, 6, 0, 0);

[success2, output2] = run_fsmoshd_operational(start_time, end_time, true);

disp("success2: " + success2)
disp("output2: " + output2)
