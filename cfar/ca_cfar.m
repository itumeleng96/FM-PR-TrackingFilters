%Function to do CA CFAR
%pf: Probability of False Alarm 
%N : Number of samples arounf the Cell Under Test
%G : Number of Guard cells
%cut : the magnitiude cells of the Range Doppler map

function [targetClusters,RDM,rdm_] = ca_cfar(RDM,pfa,fs,fd_max,td_max,maxTresholdValue)
    % Generates a CFAR output map in the range-doppler domain
    
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
    
        
        RDM_final_snr = pow2db(RDM_final_SNR);
        rdm_ = RDM_final_snr;

        max_db = max(max(RDM_final_snr));
    
        centroid_threshold = max_db -maxTresholdValue; %Set Centroid to 20db less
        [row, column] = find(RDM_final_snr> centroid_threshold);
    
        % Get target clusters as Bistatic Range and Doppler values
    
        range_values = range(column.');         % convert time values to range values
        frequency_values = frequency(row.');    % get frequency values
    
        targetClusters = [range_values; frequency_values]; % store centroids as range and frequency
    
    
    end