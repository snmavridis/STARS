% --- MATLAB Simulation of FMCW Sweep Generation, ADC, and Range Processing ---
clear; clc; close all;

% --- Define Inputs (Match these to your C++ simulation state) ---
% Radar Parameters
N          = 1024;       % Samples per sweep (Must match C++ aim[76])
Tm         = 50e-6;      % Sweep time (s) (Must match C++ aim[77])
DeltaF     = 75e6;       % Frequency excursion (Hz) (Non-aliasing, match C++ aim[78])
lambda     = 0.03;       % Wavelength (m)
Ptx        = 1000.0;     % Transmit Power (W)
Gtx_linear = 10^(20/10); % Transmit Gain (linear) (Assuming 20dB)
Grx_linear = 10^(20/10); % Receive Gain (linear)  (Assuming 20dB)
RCS        = 0.1;        % Target RCS (m^2)

% Target State (Example values - replace with actual values from C++ debug)
R      = 241.0;      % Current interpolated Range (m)
Rdot   = -100.0;     % Current interpolated Range Rate (m/s) (Negative for closing)

% Simulation State
t_chirp    = 0.020;      % Absolute time at the start of this chirp (s) (Example)

% ADC Parameters
adc_bits          = 12;
adc_voltage_range = 1.0; % Corresponds to +/- 0.5V

% Debugging Flag
debug_range_dump = true; % Set to true to plot the range profile

% --- Derived Parameters ---
C   = 3e8;           % Speed of light (m/s)
m   = DeltaF / Tm;   % LFM slope (Hz/s)
Fs  = N / Tm;        % Sampling Frequency (Hz)
pi  = acos(-1);

% --- 1. Calculate Received Power and Amplitude ---
rre_num   = Ptx * Gtx_linear * Grx_linear * lambda^2 * RCS;
rre_denom = (4*pi)^3 * R^4;
Prx       = rre_num / rre_denom; % Received Power (W)
amplitude = sqrt(Prx);           % Baseband signal amplitude (sqrt(W))

fprintf('Range = %.2f m, Rdot = %.2f m/s\n', R, Rdot);
fprintf('Calculated Rx Power = %.3e W\n', Prx);
fprintf('Calculated Amplitude = %.3e\n', amplitude);

% --- 2. Generate Baseband Sweep ---
t_vec = (0:N-1) * (Tm / N); % Time vector within the sweep (0 to Tm)

fb = (2 * m * R) / C;       % True beat frequency (Hz)
fd = (-2 * Rdot) / lambda;  % Doppler frequency (Hz)

fprintf('Sampling Freq (Fs) = %.2f MHz\n', Fs/1e6);
fprintf('Beat Freq (fb)     = %.2f MHz\n', fb/1e6);
fprintf('Doppler Freq (fd)  = %.2f kHz\n', fd/1e3);

% Check for aliasing based on parameters (should pass with DeltaF=75e6)
if abs(fb) > Fs/2
    warning('Beat frequency fb=%.2f MHz is aliased for Fs=%.2f MHz!', fb/1e6, Fs/1e6);
    % If you needed to simulate aliasing (using the C++ fb_alias logic):
    % fb_alias = mod(fb + Fs/2, Fs) - Fs/2;
    % phase_vec = 2 * pi * (fb_alias * t_vec - fd * t_chirp);
else
    fprintf('Beat frequency is not aliased.\n');
    phase_vec = 2 * pi * (fb * t_vec - fd * t_chirp); % Use true fb
end

% Generate complex sweep signal: A * exp(j*phase)
sweep_ideal = amplitude * exp(1j * phase_vec);
sweep_ideal = sweep_ideal(:); % Ensure it's a column vector
figure; % New figure for spectrogram

% Spectrogram parameters
window_length = 128; % Length of the FFT window (adjust for resolution trade-off)
overlap_percent = 90; % Percentage overlap between windows (adjust for smoothness)
noverlap = round(window_length * overlap_percent / 100); % Overlap in samples
nfft = window_length; % FFT length (can be larger for interpolation)

% Calculate and plot the spectrogram
% Use [] for default window (Hamming), nfft, and Fs for correct frequency axis
spectrogram(sweep_ideal, window_length, noverlap, nfft, Fs, 'yaxis');

% Customize the plot
title('Spectrogram of Ideal Baseband Sweep');
xlabel('Time (s)'); % Spectrogram function labels x-axis automatically
ylabel('Frequency (Hz)'); % Spectrogram function labels y-axis automatically
colorbar; % Add a color bar to show power intensity

% Adjust y-axis limits if needed (e.g., focus around fb)
ylim_mhz = [-Fs/2, Fs/2] / 1e6; % Limits in MHz based on sampling freq
% ylim(ylim_mhz * 1e6); % Set limits in Hz if 'yaxis' worked directly
% ylabel('Frequency (MHz)'); % Adjust label if you scale the axis
disp('Note: Spectrogram Y-axis shows frequencies from 0 to Fs/2.');
disp('The beat frequency should appear as a horizontal line.');
% --- 3. ADC Model ---
sweep_quantized = quantize_adc_matlab(sweep_ideal, adc_bits, adc_voltage_range);

% --- 4. Compute Range Profile (Dechirp + FFT) ---
range_profile_fft = compute_range_profile_matlab(sweep_quantized, N, Tm, m, debug_range_dump);

% --- 5. Analysis (Optional) ---
if debug_range_dump
    % Find the peak bin in the plotted profile
    power_profile_db = 10*log10(abs(range_profile_fft).^2);
    [max_power_db, peak_bin_index] = max(power_profile_db); % MATLAB index (1-based)

    % Calculate expected bin
    range_resolution = C / (2 * DeltaF);
    expected_bin = R / range_resolution; % Theoretical bin (0-based concept)

    fprintf('\n--- Analysis ---\n');
    fprintf('Range Resolution = %.3f m\n', range_resolution);
    fprintf('Expected bin (0-based) = %.2f\n', expected_bin);
    fprintf('Observed peak bin (1-based MATLAB index) = %d\n', peak_bin_index);
    fprintf('Observed peak power = %.2f dB\n', max_power_db);

    % Note: MATLAB's FFT output index 'k' (1 to N) corresponds to frequencies
    % (k-1)*Fs/N for k = 1 to N/2+1 (positive freqs)
    % (k-1-N)*Fs/N for k = N/2+2 to N (negative freqs wrapped around)
    % The beat frequency fb should correspond to index k = fb*N/Fs + 1
    expected_fft_index = fb * N / Fs + 1;
     fprintf('Expected FFT Index (1-based) based on fb = %.2f\n', expected_fft_index);

end


% --- Helper Functions (Equivalent to C++ versions) ---

function quantized_signal = quantize_adc_matlab(input_signal, num_bits, voltage_range)
    % Applies AGC and quantization
    peak_abs = max(abs(input_signal));
    agc_signal = input_signal; % Default if peak is zero

    % Step 1: AGC normalize to +/- (voltage_range / 2)
    if peak_abs > 1e-9 % Avoid division by zero/tiny numbers
        gain = (voltage_range / 2.0) / peak_abs;
        agc_signal = input_signal * gain;
    end

    % Step 2: Quantization
    max_val = voltage_range / 2.0;
    levels = 2^num_bits;
    scale_factor = (2^(num_bits - 1)) / max_val; % Scale factor to integer range

    real_part = real(agc_signal);
    imag_part = imag(agc_signal);

    % Quantize to integer codes
    q_real_int = round(real_part * scale_factor);
    q_imag_int = round(imag_part * scale_factor);

    % Clamp to ADC limits
    max_code = (2^(num_bits - 1)) - 1;
    min_code = -(2^(num_bits - 1)); % Symmetric range for signed representation

    q_real_int = max(min_code, min(q_real_int, max_code));
    q_imag_int = max(min_code, min(q_imag_int, max_code));

    % Reconstruct quantized voltage levels
    quantized_real = q_real_int / scale_factor;
    quantized_imag = q_imag_int / scale_factor;

    quantized_signal = quantized_real + 1j * quantized_imag;
end

function range_fft = compute_range_profile_matlab(sweep, N, Tm, m, debug_dump)
    % Performs dechirp, windowing, and FFT
    pi = acos(-1);

    % Generate chirp template
    t_vec = (0:N-1)' * (Tm / N); % Column vector
    phase_template = pi * m * t_vec.^2;
    chirp_template = exp(1j * phase_template);

    % Ensure sweep is a column vector
    sweep = sweep(:);

    % Dechirp
    if size(sweep, 1) ~= N || size(chirp_template, 1) ~= N
        error('Sweep size (%d) or Template size (%d) does not match N (%d)', ...
              size(sweep,1), size(chirp_template,1), N);
    end
    dechirped_sweep = sweep .* conj(chirp_template); % Element-wise

    % Windowing (Hann)
    hann_window = hann(N); % Use built-in Hann window
    windowed_sweep = dechirped_sweep .* hann_window;

    % FFT
    range_fft = fft(windowed_sweep);

    % Debug Dump / Plot
    if debug_dump
        power_profile = abs(range_fft).^2;
        figure;
        plot(0:N-1, 10*log10(power_profile + 1e-20)); % Add small epsilon to avoid log(0)
        xlabel('Range Bin Index (0 to N-1)');
        ylabel('Power (dB)');
        title('Range Profile (Dechirp + FFT)');
        grid on;
        xlim([0 N-1]); % Ensure x-axis shows all bins
        % save('range_profile_matlab.mat', 'power_profile'); % Optional save
        drawnow; % Ensure plot updates
    end
end