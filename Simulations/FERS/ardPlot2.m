function [y, ard_] = ardPlot2(s1, s2, fs, fd_max, td_max, index, ard, f)

    % Parameters for zoom and axis limits
    ylim_upper = 200;
    ylim_lower = -200;
    
    c = 299792458; % Speed of light
    N = length(s1); % Number of points
    Ndop = ceil(N * fd_max / fs); % Points corresponding to fd_max
    Ndelay = floor(td_max * fs); % Points corresponding to td_max

    % Define custom delay (in seconds) inside the function
    custom_delay = 0.25e-4; % Example: 0.5 seconds

    % Initialize variables
    s2_pad = [zeros(Ndelay, 1); s2];
    y1_custom_delay = zeros(1, 2 * fd_max + 1); % For custom delay only

    % Convert custom delay to sample index
    delay_samples = round(custom_delay * fs);
    if delay_samples < 0 || delay_samples > Ndelay
        error('Custom delay is out of range.');
    end

    % Cross-correlation computation for the custom delay
    tic
    % Adjust for the custom delay
    temp = s1 .* conj(s2_pad(Ndelay + 1 - delay_samples : N + Ndelay - delay_samples)); % Dot-product with custom delay
    temp = temp .* hanning(N); % Windowing
    temp2 = fftshift(fft(temp, N)); % FFT
    y1_custom_delay = temp2(floor(N / 2) + 1 - Ndop : floor(N / 2) + 1 + Ndop); % Frequency bins of interest
    toc

    % Convert to power and normalize
    y = abs(y1_custom_delay) .^ 2; % Power conversion
    y = y ./ max(y); % Normalize max to 1

    % Frequency axis
    frequency = -fd_max:1:fd_max;

    % Update the stored ARD if necessary
    if index >= 1
        ard = y; % Store the result for custom delay
    end
    ard_ = ard;

    % Plotting Doppler response at the specified custom delay
    figure(f);
    plot(frequency, 10 * log10(y), 'LineWidth', 1.5, 'Color', 'b');
    xlabel('Doppler Frequency [Hz]', 'FontSize', 10);
    ylabel('Amplitude [dB]', 'FontSize', 10);
    grid on;
    ylim([ylim_lower ylim_upper]);
    title('Doppler Response');
    %text(frequency(end)/10, max(10*log10(y)), "Time: " + num2str(index) + "s", 'Color', 'k');
    drawnow;

end
