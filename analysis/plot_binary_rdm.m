function plot_binary_rdm(filename)
    figure
    data = load(filename);  % load space-separated binary matrix
    imagesc(data');          % image where rows = doppler, cols = range
    colormap(parula);         % 0 = black, 1 = white
 xlabel('Range Bin');
ylabel('Doppler Bin');
    title('Binary RDM Detection Map');
    colorbar;
    axis xy;                % flip Y-axis so Doppler increases upwards
end