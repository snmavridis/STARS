% --- MATLAB Script to Demonstrate ADC Quantization ---
clear; clc; close all;

% --- Parameters ---
amplitude = 5.0;         % Amplitude of the sine wave (Volts)
frequency = 100;         % Frequency of the sine wave (Hz)
fs = 10000;              % Sampling frequency (Hz) - must be > 2*frequency (Nyquist)
duration = 2 / frequency;% Simulate for 2 cycles
num_bits_adc1 = 3;       % Number of bits for ADC 1
num_bits_adc2 = 5;       % Number of bits for ADC 2
voltage_range = 2 * amplitude; % ADC full scale range (+/- amplitude)

% --- Generate Analog Signal (Simulated) ---
t = 0 : 1/fs : duration - 1/fs; % Time vector
analog_signal = amplitude * sin(2 * pi * frequency * t);

% --- ADC 1 Simulation (3-bit) ---
quantized_signal_3bit = quantize_adc_matlab(analog_signal, num_bits_adc1, voltage_range);

% --- ADC 2 Simulation (5-bit) ---
quantized_signal_5bit = quantize_adc_matlab(analog_signal, num_bits_adc2, voltage_range);

% --- Plotting ---
fig = figure;
% set(gcf, 'Position', [100, 100, 800, 600]); % Adjust figure size

% Plot Analog Signal
% tiledlayout(3,1)
% nexttile
% plot(t * 1000, analog_signal, 'k:', 'LineWidth', 2.5, 'DisplayName', 'Analog Signal');
% title('Analog Wave')
% grid on
% % Plot 3-bit ADC Output
% nexttile
% stairs(t * 1000, quantized_signal_3bit, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('%d-bit ADC Output', num_bits_adc1));
% grid on
% title('3-bit Digital')
% % Plot 5-bit ADC Output
% nexttile
% stairs(t * 1000, quantized_signal_5bit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('%d-bit ADC Output', num_bits_adc2));
% grid on
% title('5-bit Digital')
% % --- Plot Customization ---
% ax = axes(fig);
% han = gca;
% han.Visible = 'off';


plot(t * 1000, analog_signal, 'k:', 'LineWidth', 2.5, 'DisplayName', 'Analog Signal');
hold on
% Plot 3-bit ADC Output
stairs(t * 1000, quantized_signal_3bit, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('%d-bit ADC Output', num_bits_adc1));
% Plot 5-bit ADC Output
stairs(t * 1000, quantized_signal_5bit, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('%d-bit ADC Output', num_bits_adc2));
grid on
% --- Plot Customization ---
% ax = axes(fig);
% han = gca;
% han.Visible = 'off';




title(sprintf('ADC Quantization Comparison (%d-bit vs %d-bit)', num_bits_adc1, num_bits_adc2));
xlabel('Time (ms)');

ylabel('Voltage (V)');
han.YLabel.Visible = 'on';

grid on;
ylim([-amplitude*1.1, amplitude*1.1]); % Set y-axis limits slightly larger
hold off;

% Display Quantization Levels (Optional)
levels_3bit = 2^num_bits_adc1;
step_3bit = voltage_range / levels_3bit;
fprintf('%d-bit ADC:\n', num_bits_adc1);
fprintf('  Number of Levels = %d\n', levels_3bit);
fprintf('  Voltage Step (LSB) = %.4f V\n', step_3bit);

levels_5bit = 2^num_bits_adc2;
step_5bit = voltage_range / levels_5bit;
fprintf('%d-bit ADC:\n', num_bits_adc2);
fprintf('  Number of Levels = %d\n', levels_5bit);
fprintf('  Voltage Step (LSB) = %.4f V\n', step_5bit);


% --- Helper Function for Quantization (Real Signal Version) ---
function quantized_signal = quantize_adc_matlab(input_signal, num_bits, voltage_range)
    % Quantizes a real input signal simulating a symmetric ADC
    % input_signal: Vector of analog voltage samples
    % num_bits: Number of bits for the ADC
    % voltage_range: Full scale voltage range (+Vmax to -Vmax, so range = 2*Vmax)

    max_adc_voltage = voltage_range / 2.0; % e.g., +5V if range is 10V
    min_adc_voltage = -max_adc_voltage;   % e.g., -5V

    % Calculate number of levels and voltage per step (LSB)
    num_levels = 2^num_bits;
    lsb_voltage = voltage_range / num_levels;

    % Determine integer code range for signed representation
    max_code = (2^(num_bits - 1)) - 1; % e.g., 3 bits -> max code = 3
    min_code = -(2^(num_bits - 1));    % e.g., 3 bits -> min code = -4

    % Quantize each sample
    quantized_signal = zeros(size(input_signal));
    for i = 1:length(input_signal)
        % Clamp input voltage to ADC range (optional but good practice)
        clamped_voltage = max(min_adc_voltage, min(input_signal(i), max_adc_voltage - lsb_voltage/2)); % Avoid hitting exact max

        % Calculate integer code - rounding to nearest level
        % Shift by half range, divide by step, then adjust offset
        % Simplified: Divide by LSB voltage and round
        quantized_code = round(clamped_voltage / lsb_voltage);

        % Clamp integer code to representable range
        quantized_code = max(min_code, min(quantized_code, max_code));

        % Reconstruct quantized voltage level
        quantized_signal(i) = quantized_code * lsb_voltage;
    end
end