function plot_rdm_power(filename)
    figure
    data = load(filename);
    imagesc(10*log10(data'));  % Transpose to match CFAR image orientation
    xlabel('Range Bin');
    ylabel('Doppler Bin');
    title('RDM Power Map [dB]');
    colormap(jet);
    colorbar;
    axis xy; % So Doppler increases upward
end