% --- Simple Representative RCS Plot ---
clear; clc; close all;

% --- Parameters ---
num_angles = 360; % Number of points to plot (e.g., 0 to 359 degrees)
angles_deg = linspace(0, 360, num_angles + 1); % Angles from 0 to 360 degrees
angles_deg = angles_deg(1:end-1); % Remove duplicate 360
angles_rad = deg2rad(angles_deg); % Convert to radians for plotting

% --- Generate Representative RCS Data (Simulated - Not Physical) ---
% (Using the same data generation as the previous example)
rcs_base_db = -10;
lobe1_angle_deg = 45; lobe2_angle_deg = 135; lobe3_angle_deg = 225; lobe4_angle_deg = 315;
lobe_width_deg = 40; lobe_gain_db = 20;
rcs_lobes = lobe_gain_db * (exp(-(angles_deg - lobe1_angle_deg).^2 / (2 * (lobe_width_deg/2)^2)) + ...
              exp(-(angles_deg - lobe2_angle_deg).^2 / (2 * (lobe_width_deg/2)^2)) + ...
              exp(-(angles_deg - lobe3_angle_deg).^2 / (2 * (lobe_width_deg/2)^2)) + ...
              exp(-(angles_deg - lobe4_angle_deg).^2 / (2 * (lobe_width_deg/2)^2)) + ...
              exp(-(angles_deg - (lobe4_angle_deg-360)).^2 / (2 * (lobe_width_deg/2)^2))); % Wrap around
null_angles_deg = [90, 180, 270]; null_depth_db = 25; null_width_deg = 5;
rcs_features = 0;
for na = null_angles_deg
    rcs_features = rcs_features - null_depth_db * (exp(-(angles_deg - na).^2 / (2 * (null_width_deg/2)^2)) + ...
                   exp(-(angles_deg - (na+360)).^2 / (2 * (null_width_deg/2)^2)) + ...
                   exp(-(angles_deg - (na-360)).^2 / (2 * (null_width_deg/2)^2)));
end
rcs_total_db = rcs_base_db + rcs_lobes + rcs_features;
rcs_noise_db = 1.5 * randn(size(rcs_total_db));
rcs_final_db = rcs_total_db + rcs_noise_db;

% --- Plotting ---
figure('Color', 'none'); % Set figure background color to none (transparent)

% Use polaraxes for more control over appearance
pax = polaraxes;

% Plot the data on the polar axes
polarplot(pax, angles_rad, rcs_final_db, 'b-', 'LineWidth', 1.5); % Blue line

% --- Customize Appearance ---
title(pax, 'Representative Monostatic RCS Pattern');
% Set radial limits
min_rcs = min(rcs_final_db);
max_rcs = max(rcs_final_db);
radial_min = floor(min_rcs / 10) * 10 - 10;
radial_max = ceil(max_rcs / 10) * 10;
rlim(pax, [radial_min radial_max]);

% Customize radial ticks and labels
rticks(pax, radial_min:10:radial_max);
rlabel_str = sprintfc('%d dBsm', radial_min:10:radial_max);
rticklabels(pax, rlabel_str);

% Customize angular ticks
thetaticks(pax, 0:45:315);

% Set angle direction and zero location
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';

% --- Set Colors ---
pax.Color = 'none';         % Set axes background to transparent ('none')
pax.GridColor = 'w';        % Set grid lines to white ('w')
pax.MinorGridColor = 'w';   % Set minor grid lines to white
pax.ThetaColor = 'w';       % Set angular axis line/tick color to white
pax.RColor = 'w';           % Set radial axis line/tick color to white

% Set text colors (Title, Axis Labels) - Requires accessing the properties
pax.Title.Color = 'w';
% Accessing axis labels requires finding the hidden handles usually
% Easiest way for labels is often setting default text color *before* plotting
% Or setting them after, but requires knowing properties:
% pax.RAxis.Label.String = 'RCS (dBsm)'; % Re-set label if needed
% pax.RAxis.Label.Color = 'w'; % Doesn't always work directly
% pax.ThetaAxis.Label.String = 'Angle (degrees)';
% pax.ThetaAxis.Label.Color = 'w'; % Doesn't always work directly

% Safer approach: Set default text color for the figure
set(gcf, 'DefaultTextColor', 'w');
% Re-apply title etc. AFTER setting default text color if needed
title(pax, 'Representative Monostatic RCS Pattern'); % Title will now be white
% If axis labels don't turn white, it's trickier. Might need specific handle access.

fprintf('Plot generated with transparency and white text/grid.\n');
fprintf('Min RCS: %.2f dBsm, Max RCS: %.2f dBsm\n', min(rcs_final_db), max(rcs_final_db));

% --- Optional: Saving with Transparency ---
% When saving, ensure the format supports transparency (like PNG)
   exportgraphics(gcf, "rcs.png", ...
                   'ContentType', 'image', ...
                   'BackgroundColor', 'none', ...
                   'Resolution', 300);