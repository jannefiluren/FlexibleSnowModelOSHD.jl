clear; clc

experiment_paths = ["SNFRAC_0", "SNFRAC_3", "SNFRAC_4", "FOREST"];
base_path = "D:\julia\FSM_HS_julia";
archive_path = "ARCHIVE_20250804";

for experiment_path = experiment_paths

  test_path = fullfile(base_path,experiment_path);
  ref_path = fullfile(base_path,archive_path,experiment_path);

  disp("==========================================================")
  disp("Running verification for " + experiment_path)
  disp("==========================================================")

  files_test = dir(fullfile(test_path,"*.mat"));

  disp("Testing " + length(files_test) + " files...")
  disp("Test data path: " + test_path)
  disp("Reference data path: " + ref_path)

  test_passed = true;
  for file = files_test'
    test_data = load(fullfile(test_path,file.name));
    ref_data = load(fullfile(ref_path,file.name));
    if ~isequaln(test_data,ref_data)
      test_passed = false;
    end
  end

  disp("All files equal: " + test_passed)

end
