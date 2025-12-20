% --- Demonstrate ADC Sample Rate vs. Range Binning ---
clear; clc; close all;

% --- Common Radar Parameters ---
DeltaF     = 75e6;     % Frequency Excursion (Hz) (Choose non-aliasing for lower Fs)
Tm         = 50e-6;    % Sweep Time (s)
C          = 3e8;      % Speed of light (m/s)
lambda     = 0.03;     % Wavelength (m) - Needed for Doppler calc

% --- Target Parameters ---
R_target   = 200;      % Target Range (m)
Rdot_target= -50;      % Target Range Rate (m/s) (Negative = closing)

% --- Signal Generation Amplitude (Simplified - Assumed constant) ---
A          = 1.0;      % Assume normalized amplitude for simplicity

% --- ADC & Processing Parameters ---
Fs_low     = 5e6;      % Low ADC Sample Rate (5 Msps)
Fs_high    = 20e6;     % High ADC Sample Rate (20 Msps)

N_low      = round(Fs_low * Tm); % Number of samples/bins for low Fs
N_high     = round(Fs_high * Tm);% Number of samples/bins for high Fs

fprintf('Low Fs Case: Fs = %.2f Msps, N = %d bins\n', Fs_low/1e6, N_low);
fprintf('High Fs Case: Fs = %.2f Msps, N = %d bins\n', Fs_high/1e6, N_high);

% --- Derived Parameters ---
m          = DeltaF / Tm; % LFM Slope
pi         = acos(-1);
range_resolution = C / (2 * DeltaF); % Same for both cases
fprintf('Range Resolution = %.3f m\n', range_resolution);

% --- Calculate Expected Beat Frequency ---
fb_target = (2 * m * R_target) / C;
fd_target = (-2 * Rdot_target) / lambda;
fprintf('Target Beat Frequency (fb) = %.3f MHz\n', fb_target/1e6);

% --- Check for Aliasing ---
if abs(fb_target) > Fs_low / 2
    warning('Target beat frequency WILL ALIAS for Fs_low=%.2f Msps!', Fs_low/1e6);
    % Calculate aliased frequency for low Fs case
    fb_low_alias = mod(fb_target + Fs_low/2, Fs_low) - Fs_low/2;
    fprintf('  Aliased fb for low Fs case = %.3f MHz\n', fb_low_alias/1e6);
else
    fb_low_alias = fb_target; % No aliasing
    fprintf('Target beat frequency will NOT alias for Fs_low=%.2f Msps.\n', Fs_low/1e6);
end
if abs(fb_target) > Fs_high / 2
    warning('Target beat frequency WILL ALIAS for Fs_high=%.2f Msps!', Fs_high/1e6);
    fb_high_alias = mod(fb_target + Fs_high/2, Fs_high) - Fs_high/2;
     fprintf('  Aliased fb for high Fs case = %.3f MHz\n', fb_high_alias/1e6);
else
    fb_high_alias = fb_target; % No aliasing
    fprintf('Target beat frequency will NOT alias for Fs_high=%.2f Msps.\n', Fs_high/1e6);
end

% --- Generate Baseband Signals ---
% Low Fs Case
t_vec_low = (0:N_low-1) / Fs_low; % Time vector for low Fs
phase_vec_low = 2 * pi * (fb_low_alias * t_vec_low); % Use aliased fb, ignore fd for range profile clarity
sweep_low = A * exp(1j * phase_vec_low);

% High Fs Case
t_vec_high = (0:N_high-1) / Fs_high; % Time vector for high Fs
phase_vec_high = 2 * pi * (fb_high_alias * t_vec_high); % Use aliased fb, ignore fd for range profile clarity
sweep_high = A * exp(1j * phase_vec_high);

% --- Process Signals (Window + FFT) ---
% Low Fs Case
window_low = hann(N_low);
range_fft_low = fft(sweep_low(:) .* window_low); % Apply window and FFT
power_profile_low = abs(range_fft_low).^2;

% High Fs Case
window_high = hann(N_high);
range_fft_high = fft(sweep_high(:) .* window_high); % Apply window and FFT
power_profile_high = abs(range_fft_high).^2;

% --- Create Range Axes ---
% Max unambiguous range depends on Fs
R_max_low = (Fs_low / 2) * Tm * C / (2 * DeltaF);
range_axis_low = linspace(0, R_max_low, N_low);

R_max_high = (Fs_high / 2) * Tm * C / (2 * DeltaF);
range_axis_high = linspace(0, R_max_high, N_high);

% --- Plotting ---
figure('Position', [50, 50, 900, 700]);

% Plot 1: Power vs. Bin Index
subplot(2, 1, 1);
plot(0:N_low-1, 10*log10(power_profile_low + 1e-20), 'b-o', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', sprintf('Fs = %.1f Msps (N=%d)', Fs_low/1e6, N_low));
hold on;
plot(0:N_high-1, 10*log10(power_profile_high + 1e-20), 'r-x', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', sprintf('Fs = %.1f Msps (N=%d)', Fs_high/1e6, N_high));
title('Range Profile vs. Bin Index');
xlabel('FFT Bin Index (k)');
ylabel('Power (dB)');
legend('show', 'Location', 'SouthOutside');
grid on;
xlim([0 max(N_low, N_high)-1]); % Adjust x-limit if needed

% Calculate expected bins (based on potentially aliased frequencies)
expected_bin_low = fb_low_alias * N_low / Fs_low; % Can be negative if aliased
if expected_bin_low < 0, expected_bin_low = expected_bin_low + N_low; end % Wrap around
xline(expected_bin_low, 'b--', 'Label', sprintf('Expected Low Fs Bin (%.1f)', expected_bin_low));

expected_bin_high = fb_high_alias * N_high / Fs_high;
if expected_bin_high < 0, expected_bin_high = expected_bin_high + N_high; end
xline(expected_bin_high, 'r--', 'Label', sprintf('Expected High Fs Bin (%.1f)', expected_bin_high));
hold off;

% Plot 2: Power vs. Range (meters)
subplot(2, 1, 2);
plot(range_axis_low, 10*log10(power_profile_low + 1e-20), 'b-o', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', sprintf('Fs = %.1f Msps (Rmax=%.0fm)', Fs_low/1e6, R_max_low));
hold on;
plot(range_axis_high, 10*log10(power_profile_high + 1e-20), 'r-x', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', sprintf('Fs = %.1f Msps (Rmax=%.0fm)', Fs_high/1e6, R_max_high));
title('Range Profile vs. Calculated Range');
xlabel('Range (m)');
ylabel('Power (dB)');
legend('show', 'Location', 'SouthOutside');
grid on;
xline(R_target, 'k:', 'LineWidth', 1.5, 'Label', sprintf('True Target Range (%.0f m)', R_target)); % Mark true range
xlim([0 max(R_max_low, R_max_high)]); % Adjust x-limit
hold off;

sgtitle(sprintf('Effect of ADC Sample Rate (Fs) on Range Profile (DeltaF=%.0fMHz, Tm=%.0fus)', DeltaF/1e6, Tm*1e6));