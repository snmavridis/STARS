% --- MATLAB Script to Parse Multi-Line traj.txt and Plot Trajectories ---
clear; clc; close all;

% --- Configuration ---
filename = 'traj.txt'; % Name of your trajectory data file
num_header_lines = 5;  % Lines to skip before numerical data starts
num_variables = 19;    % EXACT number of variables expected per time step

% --- Define Expected Header Names (in order) ---
% Used for assigning data to struct fields later
header_names = {
    'time', 'mach_m1', 'SAEL1_m1', 'SAEL2_m1', 'SAEL3_m1', ...
    'VAEL1_m1', 'VAEL2_m1', 'VAEL3_m1', 'psivlx_m1', 'thtvlx_m1', ...
    'mach_a1', 'SAEL1_a1', 'SAEL2_a1', 'SAEL3_a1', 'VAEL1_a1', ...
    'VAEL2_a1', 'VAEL3_a1', 'psivlx_a1', 'thtvlx_a1'
};

% --- Basic Validation ---
if length(header_names) ~= num_variables
    error('Mismatch between expected_headers count (%d) and num_variables (%d). Check script configuration.', ...
          length(header_names), num_variables);
end

% --- Read Numerical Data using textscan ---
fprintf('Reading data from: %s\n', filename);
fileID = fopen(filename, 'r');
if fileID == -1
    error('Cannot open file: %s', filename);
end

% Skip header lines manually
for i = 1:num_header_lines
    fgetl(fileID);
    if feof(fileID) && i < num_header_lines
        fclose(fileID);
        error('File ended unexpectedly within header lines.');
    end
end

% Read all remaining floating-point numbers into a single cell array
% '%f' reads floats, whitespace is the default delimiter
data_cell = textscan(fileID, '%f');
fclose(fileID);

% Check if data was read
if isempty(data_cell) || isempty(data_cell{1})
    error('No numerical data found after header lines in %s.', filename);
end

% Convert cell array to a column vector
data_vector = data_cell{1};
fprintf('Read a total of %d numerical values.\n', length(data_vector));

% --- Remove Trailing Plot Marker Line(s) ---
% Find the index of the first negative time value (-1.0)
end_marker_index = find(data_vector < 0, 1, 'first');

if ~isempty(end_marker_index)
    % Check if the marker is at the expected position (start of a 19-value block)
    if mod(end_marker_index - 1, num_variables) == 0
        fprintf('Found end marker at index %d. Removing marker and subsequent data.\n', end_marker_index);
        % Keep data *before* the marker
        data_vector = data_vector(1 : end_marker_index - 1);
    else
        warning('End marker found at index %d, which is not the start of a %d-value block. File might be corrupted. Attempting to remove based on index.', ...
                end_marker_index, num_variables);
        % Cautiously remove from marker onwards, might lead to reshape error later
         data_vector = data_vector(1 : end_marker_index - 1);
    end
else
    fprintf('No plot marker row (time < 0) found.\n');
end

% --- Reshape Data ---
num_elements = length(data_vector);
if mod(num_elements, num_variables) ~= 0
    error('Total number of numerical elements (%d) is not a multiple of expected variables per row (%d). Check file format or num_variables.', ...
          num_elements, num_variables);
end

num_timesteps = num_elements / num_variables;
fprintf('Reshaping data into %d time steps with %d variables each.\n', num_timesteps, num_variables);

% Reshape the vector into a matrix: N_timesteps x N_variables
% MATLAB reshapes column-wise by default, so transpose after reshape
data = reshape(data_vector, num_variables, num_timesteps)';

% --- Assign Data to Struct Fields ---
fprintf('Assigning data columns to struct fields...\n');
traj_data = struct();
col_idx_m1_x = -1; col_idx_m1_y = -1; col_idx_m1_z = -1;
col_idx_a1_x = -1; col_idx_a1_y = -1; col_idx_a1_z = -1;

for i = 1:num_variables
    var_name = matlab.lang.makeValidName(header_names{i});
    if ~strcmp(var_name, header_names{i})
        fprintf('Warning: Header "%s" converted to MATLAB variable name "%s"\n', header_names{i}, var_name);
    end
    traj_data.(var_name) = data(:, i);

    % Store indices while we're looping
    if strcmp(header_names{i}, 'SAEL1_m1'), col_idx_m1_x = i; end
    if strcmp(header_names{i}, 'SAEL2_m1'), col_idx_m1_y = i; end
    if strcmp(header_names{i}, 'SAEL3_m1'), col_idx_m1_z = i; end
    if strcmp(header_names{i}, 'SAEL1_a1'), col_idx_a1_x = i; end
    if strcmp(header_names{i}, 'SAEL2_a1'), col_idx_a1_y = i; end
    if strcmp(header_names{i}, 'SAEL3_a1'), col_idx_a1_z = i; end
end
fprintf('Data assigned to struct "traj_data".\n');

% Check if indices were found
if any([col_idx_m1_x, col_idx_m1_y, col_idx_m1_z, col_idx_a1_x, col_idx_a1_y, col_idx_a1_z] < 0)
   error('Failed to find all required SAEL column names in header_names list.');
end

% --- Plotting 3D Trajectories ---
fprintf('Generating 3D trajectory plot...\n');

figure;
hold on;

% Plot Missile Trajectory (m1) - Blue
plot3(traj_data.SAEL1_m1, traj_data.SAEL2_m1, traj_data.SAEL3_m1, 'b-', 'LineWidth', 1.5);

% Plot Aircraft Trajectory (a1) - Red
plot3(traj_data.SAEL1_a1, traj_data.SAEL2_a1, traj_data.SAEL3_a1, 'r-', 'LineWidth', 1.5);

% Mark Start and End Points
plot3(traj_data.SAEL1_m1(1), traj_data.SAEL2_m1(1), traj_data.SAEL3_m1(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8); % m1 Start
plot3(traj_data.SAEL1_m1(end), traj_data.SAEL2_m1(end), traj_data.SAEL3_m1(end), 'bx', 'MarkerSize', 10, 'LineWidth', 2); % m1 End

plot3(traj_data.SAEL1_a1(1), traj_data.SAEL2_a1(1), traj_data.SAEL3_a1(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); % a1 Start
plot3(traj_data.SAEL1_a1(end), traj_data.SAEL2_a1(end), traj_data.SAEL3_a1(end), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % a1 End

% --- Plot Customization ---
title('3D Trajectories (m1 vs a1)');
xlabel('X Position (SAEL1 - m)');
ylabel('Y Position (SAEL2 - m)');
zlabel('Z Position (SAEL3 - m)');
legend('Missile (m1)', 'Aircraft (a1)', 'm1 Start', 'm1 End', 'a1 Start', 'a1 End', 'Location', 'best');
grid on;
axis equal;
view(3);
rotate3d on;

hold off;

fprintf('Plotting complete.\n');