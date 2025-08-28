function [success, output] = run_fsmoshd_operational(start_time, end_time, restart)

arguments
  start_time double
  end_time double
  restart = false
end

start_str = datestr(start_time, 'yyyy-mm-ddTHH:MM:SS');
end_str = datestr(end_time, 'yyyy-mm-ddTHH:MM:SS');

if restart
  restart_str = ' true';
else
  restart_str = ' false';
end

% Build Julia command
julia_script_path = fullfile(fileparts(mfilename('fullpath')), '..', 'operational', 'matlab_interface.jl');

% Construct command
cmd = sprintf('julia --project="%s" "%s" "%s" "%s"%s', ...
  fullfile(fileparts(julia_script_path), '..'), ...
  julia_script_path, start_str, end_str, restart_str);

% Run system command
[status, output] = system(cmd);

% Check success
success = (status == 0);

% Strip output
position = regexp(output, "(Success|Error)");
output = strip(output(position:end));

end