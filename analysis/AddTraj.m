function cadacOutput = ReadCadacOutput_Cruise5(filepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(filepath,'r');

for i = 1:6
    fgetl(fid);
end

i = 0;
k = 1;

while ~feof(fid)
textline = fgetl(fid);

if mod(i,4) == 0
    
    lineData = textscan(textline,'%f %f %f %f %f');
    lineData = cell2mat(lineData);

    cadacOutput.time(k)         = lineData(1);
    cadacOutput.mach_m1(k)      = lineData(2);
    cadacOutput.SAEL1_m1(k)     = lineData(3);
    cadacOutput.SAEL2_m1(k)     = lineData(4);
    cadacOutput.SAEL3_m1(k)     = lineData(5);

end

if mod(i,4) == 1
    lineData = textscan(textline,'%f %f %f %f %f');
    lineData = cell2mat(lineData);

    cadacOutput.VAEL1_m1(k)     = lineData(1);
    cadacOutput.VAEL2_m1(k)     = lineData(2);
    cadacOutput.VAEL3_m1(k)     = lineData(3);
    cadacOutput.psivlx_m1(k)      = lineData(4);
    cadacOutput.thtvlx_m1(k)     = lineData(5);

end

if mod(i,4) == 2
    lineData = textscan(textline,'%f %f %f %f %f');
    lineData = cell2mat(lineData);

    cadacOutput.mach_a1(k)   = lineData(1);
    cadacOutput.SAEL1_a1(k)   = lineData(2);
    cadacOutput.SAEL2_a1(k)    = lineData(3);
    cadacOutput.SAEL3_a1(k)    = lineData(4);
    cadacOutput.VAEL1_a1(k)    = lineData(5);
    
end

if mod(i,4) == 3
    lineData = textscan(textline,'%f %f %f %f');
    lineData = cell2mat(lineData);

    cadacOutput.VAEL2_a1(k)    = lineData(1);
    cadacOutput.VAEL3_a1(k)    = lineData(2);
    cadacOutput.psivlx_a1(k)    = lineData(3);
    cadacOutput.thtvlx_a1(k)    = lineData(4);
    k = k+1;
    
end



i = i+1;


end

cadacOutput.time(end)       = [];
cadacOutput.mach_m1(end)      = [];
cadacOutput.SAEL1_m1(end)      = [];
cadacOutput.SAEL2_m1(end)      = [];
cadacOutput.SAEL3_m1(end)     = [];

cadacOutput.mach_a1(end)       = [];
cadacOutput.SAEL1_a1(end)       = [];
cadacOutput.VAEL3_m1(end)       = [];
cadacOutput.SAEL2_a1(end)        = [];
cadacOutput.SAEL3_a1(end)       = [];
cadacOutput.VAEL1_a1(end)     = [];
cadacOutput.VAEL2_a1(end)     = [];
cadacOutput.VAEL3_a1(end)      = [];
cadacOutput.psivlx_a1(end)      = [];
cadacOutput.thtvlx_a1(end)      = [];
% ss = size(cadacOutput.SAEL1_m1,1);
% 
% fig = figure;
% plot3(cadacOutput.SAEL1_m1,cadacOutput.SAEL2_m1,cadacOutput.SAEL3_m1,'b','LineWidth',2);
% hold on
% plot3(cadacOutput.SAEL1_a1,cadacOutput.SAEL2_a1,cadacOutput.SAEL3_a1,'r','LineWidth',2);
% box on
% set(gca,'ZDir','reverse');
% grid on
% zlim([-15000 -5000])
% xlim([-400 400])
%%
% Ensure data are column vectors for consistency
m1_x = cadacOutput.SAEL1_m1(:);
m1_y = cadacOutput.SAEL2_m1(:);
m1_z = cadacOutput.SAEL3_m1(:);
a1_x = cadacOutput.SAEL1_a1(:);
a1_y = cadacOutput.SAEL2_a1(:);
a1_z = cadacOutput.SAEL3_a1(:);

if length(m1_x) ~= length(m1_y) || length(m1_x) ~= length(m1_z) || ...
   length(a1_x) ~= length(a1_y) || length(a1_x) ~= length(a1_z)
    warning('Trajectory vector lengths are inconsistent.');
    % Consider erroring out or trying to proceed
end
N_m1 = length(m1_x);
N_a1 = length(a1_x);

% --- Create Figure and Plot Initial 3D Trajectories ---
open figure.fig
hold on; % Hold on early to add all elements

% Plot Main Trajectories
h_m1 = plot3(m1_x, m1_y, m1_z, 'LineWidth', 2, 'DisplayName', 'Baseline (m1)');

% --- Set Axis Limits and Direction BEFORE getting limits ---
% Set desired limits first, so projections go to these specific walls
zlim([-15000 -5000]);
xlim([-30 70]);
% Let ylim be determined automatically initially, or set it if you have a preference
% ylim([-YYY YYY]); % Example if needed
set(gca, 'ZDir', 'reverse'); % Reverse Z-axis (Altitude?)

% --- Get Axis Limits (AFTER potentially setting some manually) ---
xlim_vals = xlim;
ylim_vals = ylim;
zlim_vals = zlim; % Get the actual limits being used (might differ slightly if ylim was auto)

xmin = xlim_vals(1); 
xmax = xlim_vals(2);
ymin = ylim_vals(1); 
ymax = ylim_vals(2);
zmin = zlim_vals(1); 
zmax = zlim_vals(2); % Note: zmin/zmax might be swapped due to reverse ZDir
% xmin = -400; 
% xmax = 400;
% ymin = -12000; 
% ymax = 300;
% zmin = -12500; 
% zmax = -7500; % Note: zm

if strcmp(get(gca,'ZDir'),'reverse')
    z_floor = zmax; % The "bottom" wall is now at the max Z value
    z_ceiling = zmin; % The "top" wall is now at the min Z value
else
    z_floor = zmin;
    z_ceiling = zmax;
end

fprintf('Plot Limits Being Used:\n');
fprintf('  X: [%.2f, %.2f]\n', xmin, xmax);
fprintf('  Y: [%.2f, %.2f]\n', ymin, ymax);
fprintf('  Z: [%.2f, %.2f] (ZDir=%s, Floor=%.1f)\n', zlim_vals(1), zlim_vals(2), get(gca,'ZDir'), z_floor);

% --- Plot Projections ---
% Define projection style
proj_m1_color = [0.6 0.6 1]; % Light Blue for m1 projections
proj_a1_color = [1 0.6 0.6]; % Light Red for a1 projections
proj_style = '--';
proj_width = 1;

% --- Plot XY Projections (on the 'floor' at z = z_floor) ---
z_proj_xy_m1 = ones(N_m1, 1) * z_floor;
z_proj_xy_a1 = ones(N_a1, 1) * z_floor;
h_xy_m1 = plot3(m1_x, m1_y, z_proj_xy_m1, 'Color', proj_m1_color, 'LineStyle', proj_style, 'LineWidth', proj_width, 'HandleVisibility', 'off'); % Hide from main legend
h_xy_a1 = plot3(a1_x, a1_y, z_proj_xy_a1, 'Color', proj_a1_color, 'LineStyle', proj_style, 'LineWidth', proj_width, 'HandleVisibility', 'off');

% --- Plot XZ Projections (on the 'back wall' at y = ymax) ---
y_proj_xz_m1 = ones(N_m1, 1) * ymax;
y_proj_xz_a1 = ones(N_a1, 1) * ymax;
h_xz_m1 = plot3(m1_x, y_proj_xz_m1, m1_z, 'Color', proj_m1_color, 'LineStyle', proj_style, 'LineWidth', proj_width, 'HandleVisibility', 'off');
h_xz_a1 = plot3(a1_x, y_proj_xz_a1, a1_z, 'Color', proj_a1_color, 'LineStyle', proj_style, 'LineWidth', proj_width, 'HandleVisibility', 'off');

% --- Plot YZ Projections (on the 'side wall' at x = xmax) ---
x_proj_yz_m1 = ones(N_m1, 1) * xmax;
x_proj_yz_a1 = ones(N_a1, 1) * xmax;
h_yz_m1 = plot3(x_proj_yz_m1, m1_y, m1_z, 'Color', proj_m1_color, 'LineStyle', proj_style, 'LineWidth', proj_width, 'HandleVisibility', 'off');
h_yz_a1 = plot3(x_proj_yz_a1, a1_y, a1_z, 'Color', proj_a1_color, 'LineStyle', proj_style, 'LineWidth', proj_width, 'HandleVisibility', 'off');

% --- Final Plot Customization ---
box on;     % Ensure the box outline is drawn
grid on;
xlabel('X Position (m)');
ylabel('Y Position (m)');
zlabel('Z Position (m)');
title('3D Trajectories with Projections');
legend([h_m1, h_a1], 'Location', 'best'); % Only show main trajectories in legend
axis tight; % Adjust limits tightly around data + projections
view(3);       % Ensure 3D view
rotate3d on;   % Enable interactive rotation

hold off;

fprintf('Plotting complete with projections.\n');
zlim([-12500 zlim_vals(2)])
% ylim([-12000 0])
% xlim([-10 200])
end