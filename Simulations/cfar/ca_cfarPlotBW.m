% Cell-Averaging Constant False Alarm Rate
% Norman Pearson Detector
% Author: Itumeleng Malemela

function [targetClusters, RDM,RDM_] = ca_cfarPlotBW(RDM, pfa, fs, fd_max, td_max, index, f, rdm ,maxTreshold)
% Generates a CFAR output map in the range-doppler domain

ylim_upper = 200;
ylim_lower = -200;

y_guard_num = 4;
y_train_num = 9;

[rows, cols] = size(RDM);
RDM_final = zeros(rows, cols); 
RDM_final_SNR = zeros(rows, cols); 

threshold_values = zeros(rows, cols);

for r = 1:rows
    for c = 1:cols
        % Get guard cells sum value
        [guard_value, g_c] = subarea(RDM, c, r, y_guard_num);
        
        % Get training cells sum value
        [train_value, t_c] = subarea(RDM, c, r, y_guard_num + y_train_num);

        Num_train = t_c - g_c - 1;
        train_value = (train_value - guard_value) / Num_train;   % Noise Power
        
        alpha = Num_train * ((pfa)^(-1 / Num_train) - 1);        % Threshold factor
        threshold_val = train_value * alpha;
        threshold_values(r, c) = threshold_val;
        
        if RDM(r, c) >= threshold_val
            signal_power = RDM(r, c);
            noise_power = train_value;
            SNR = signal_power/noise_power;
            RDM_final_SNR(r,c) =SNR;

            RDM_final(r, c) = 1;  % Detected target
        else
            RDM_final(r, c) = 0;  % No target
        end
    end
end

c = 3e8;
Ndelay = floor(td_max * fs);                                  
time = 0:1/fs:Ndelay/fs;
range = time * c;
range = range / 1000;
frequency = -fd_max:1:fd_max;

rdm = RDM_final;
RDM_ = rdm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CA CFAR results as Black and White image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(f);
imagesc(range, frequency, rdm);
colormap(gca, 'gray');
set(gca, 'Color', 'black'); 
xlabel('Bistatic Range [km]', 'FontSize', 20, 'Color', 'black');
ylabel('Bistatic Doppler frequency [Hz]', 'FontSize', 20, 'Color', 'black');
grid on;
title('CA-CFAR and Centroids', 'Color', 'black', 'FontSize', 20);
ax = gca; 
ax.FontSize = 14;

% Ensure y-axis is not inverted
set(gca, 'YDir', 'normal');  % Ensure that the y-axis is correctly oriented

ylim([ylim_lower ylim_upper]);
drawnow

RDM_final_snr = pow2db(RDM_final_SNR);

max_db = max(max(RDM_final_snr));

centroid_threshold = max_db -maxTreshold; %Set Centroid to 20db less
[row, column] = find(RDM_final_snr> centroid_threshold);

% Get target clusters as Bistatic Range and Doppler values

range_values = range(column.');         % convert time values to range values
frequency_values = frequency(row.');    % get frequency values

targetClusters = [range_values; frequency_values]; % store centroids as range and frequency


end
