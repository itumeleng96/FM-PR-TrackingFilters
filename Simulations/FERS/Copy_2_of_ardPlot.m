function [y, ard_] = ardPlot(s1, s2, fs, fd_max, td_max, index, ard, f)

    % Parameters for zoom and axis limits
    ylim_upper = 200;
    ylim_lower = -200;

    c = 299792458; % Speed of light
    N = length(s1); % Number of points
    Ndelay = floor(td_max * fs); % Points corresponding to td_max
    Ndop = ceil(N * fd_max / fs); % Points corresponding to fd_max

    % Initializing temporary variables
    s2_pad = [zeros(Ndelay, 1); s2];
    y1 = zeros(Ndelay + 1, 2 * fd_max + 1);

    % Range-Doppler computation
    tic
    for k = 1:Ndelay + 1
        temp = s1 .* conj(s2_pad(Ndelay + 2 - k:N + Ndelay + 1 - k)); % Dot-product
        temp = temp .* hanning(N); % Windowing
        temp2 = fftshift(fft(temp, N)); % FFT
        y1(k, :) = temp2(floor(N / 2) + 1 - Ndop : floor(N / 2) + 1 + Ndop); % Frequency bins of interest
    end
    toc

    y = abs(y1) .^ 2; % Power conversion
    y = y ./ max(max(abs(y))); % Normalize max to 1

    % Time and frequency axis
    time = 0:1/fs:Ndelay/fs;
    range = time * c;
    frequency = -fd_max:1:fd_max;

    % Update the stored ARD if necessary
    if index >= 1
        ard = y; % In this context, ard stores the ARD data
    end
    ard_ = ard;

    % Plotting
    figure(f);

    % 3D Surface Plot
    [R, F] = meshgrid(time, frequency);
    surf(R, F, 10 * log10(ard.'), 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 0.8);

    % Adjusting color map and axis
    colormap(jet); % Use the 'jet' colormap for better visualization
    colorbar;
    xlabel('Bistatic Delay [s]', 'FontSize', 10);
    ylabel('Bistatic Doppler Frequency [Hz]', 'FontSize', 10);
    zlabel('Amplitude [dB]', 'FontSize', 10);
    title('Amplitude Range-Doppler Map');

    % Set axis limits and labels
    ylim([ylim_lower ylim_upper]);

    % Display time information
    % text(range(end)/10, frequency(end)/10, max(max(10*log10(ard.'))), "Time: " + num2str(index) + "s", 'Color', 'k');

    % Adjust view and lighting for better clarity
    view(3); % 3D view
    lighting phong; % Add phong lighting for a smooth appearance
    camlight headlight; % Add a camera light
    drawnow;
end
