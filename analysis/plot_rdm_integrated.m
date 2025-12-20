function plot_rdm_integrated(filename)
    figure
    data = load(filename);
    imagesc(data'); colormap(gray); axis xy;
    xlabel('Range Bin'); ylabel('Doppler Bin');
    title('Integrated CFAR Detections (M of N)');
end