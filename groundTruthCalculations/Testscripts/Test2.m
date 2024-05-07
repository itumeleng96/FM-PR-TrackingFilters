
Coords = [
    3806 -3116 -6884 -14979;
    20680 14010 10971 3765;
    3200 3200 4000 4000;
];

CoordsTime = [0, 20, 30, 60];

% Calculate second derivatives
dd = finalizeCubic(Coords, CoordsTime);


% Define the update interval
UpdateInterval = 1;

% Calculate the baseline

% Initialize arrays to store bistatic ranges and Doppler shifts
bistatic_ranges = [];
bistatic_doppler_shifts = [];

% Loop through each time step
prev_total_range =0;
for t = CoordsTime(1):CoordsTime(end)
    % Get interpolated coordinate at time t using getPositionCubic function
    coord = getPositionCubic(t, Coords, CoordsTime, dd);
    
    % Calculate the distance between the target and the transmitter
    range_tx_target = norm(coord - Tx_Pos);
    
    % Calculate the distance between the target and the reference receiver
    range_ref_target = norm(coord - RefRx_Pos);
    
    % Calculate the total range
    total_range = range_tx_target + range_ref_target;
    
    % Calculate the bistatic range
    bistatic_range = total_range - Baseline;
    
    % Store the bistatic range
    bistatic_ranges(end+1) = bistatic_range;
    
    % Calculate the rate of change of the total range with respect to time
    d_total_range_dt = (total_range - prev_total_range) / UpdateInterval;
    prev_total_range = total_range;
    
    % Calculate the Doppler shift
    doppler_shift = - (1 / wavelength) * (d_total_range_dt);
    
    % Store the bistatic Doppler shift
    bistatic_doppler_shifts(end+1) = doppler_shift;

    plot3(Coords(1, :), Coords(2, :), Coords(3, :), 'bo-', 'LineWidth', 2);
    hold on;
    plot3(coord(1), coord(2), coord(3), 'ro', 'MarkerSize', 10);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['3D Plot of XYZ Points at time t = ' num2str(t)]);
    grid on;
    legend('Original Points', 'Interpolated Point');
    hold off;
    pause(0.1); % Pause for a short time to visualize each time point
end
% Plot the bistatic ranges vs Doppler shifts
figure;
plot(bistatic_ranges, bistatic_doppler_shifts);
title('Bistatic Range vs Doppler Shift');
xlabel('Bistatic Range (m)');
ylabel('Bistatic Doppler Shift (Hz)');
%}