function plot_range_profile(filename)
close all
data = load(filename);
range_bins = data(:,1);
power_linear = data(:,2);
power_dB = 10*log10(power_linear);     % Convert to dB

plot(range_bins, power_dB, 'b', 'LineWidth', 1.5);
xlabel('Range Bin');
ylabel('Power (dB)');
title('Matched Filter Range Response');
grid on;
[max_power, idx] = max(power_dB);
hold on;
plot(idx, max_power, 'ro');
text(idx, max_power+3, sprintf('Peak @ Bin %d', idx));
figure
plot(range_bins, power_linear);
ylabel('Power (Linear)');
end