%Function to do CA CFAR
%pf: Probability of False Alarm 
%N : Number of samples arounf the Cell Under Test
%G : Number of Guard cells
%cut : the magnitiude cells of the Range Doppler map

function [targetClusters,RDM,rdm_] = ca_cfar(RDM,pfa,fs,fd_max,td_max,index,rdm)
%Generates a CFAR output map in the range-doppler domain

%Firstly calculate the interference power from the average of N samples in
%the vicinity of the CUT

%pf: Probability of False Alarm 
%N : Number of samples around the Cell Under Test
%G : Number of Guard cells
%cut : the magnitiude cells of the Range Doppler map

y_guard_num = 8;
y_train_num = 10;


%x_guard_num = 8;
%x_train_num = 6;

[rows, cols] = size(RDM);
RDM_final = zeros(rows, cols); 
threshold_values = zeros(rows, cols);

for r=1:rows
    for c=1:cols
        % get guard cells sum value

        %[guard_value, g_c] = subarea(RDM, c, r, y_guard_num,x_guard_num) 
        [guard_value, g_c] = subarea(RDM, c, r, y_guard_num);
        
        % Get training cells sum value
        % [train_value, t_c] = subarea(RDM, c, r, y_guard_num+y_train_num,x_guard_num+x_train_num);
        [train_value, t_c] = subarea(RDM, c, r, y_guard_num+y_train_num);

        Num_train = t_c - g_c - 1;
        train_value = (train_value - guard_value) / (Num_train);   %Noise Power
        %train_value = pow2db(train_value);                        %Noise Power in db
       
        
        alpha = Num_train*((pfa)^(-1/Num_train) - 1);              %threshold factor
        threshold_val = train_value * alpha;
        %threshold_val = train_value * offset;
        threshold_values(r,c) = threshold_val;
        
        if RDM(r,c) >= threshold_val

            % Calculate SNR for the detected cell and 
            % plot SNR instead of 0 or 1

            signal_power = RDM(r, c);
            noise_power = train_value;
            SNR = signal_power/noise_power;
            RDM_final(r,c) =SNR;
          
            
        end
    end
end

c=3e8;
Ndelay = floor(td_max*fs);                                  
time = 0:1/fs:Ndelay/fs;
range = time*c;
frequency = -fd_max:1:fd_max;

RDM_final = pow2db(RDM_final);
rdm = RDM_final ;

rdm_ = rdm;
max_db = max(max(rdm_));

RDM = RDM_final;

centroid_threshold = max_db -10; %Set Centroid to 20db less
[row, column] = find(RDM > centroid_threshold);

% Get target clusters as Bistatic Range and Doppler values

range_values = range(column.');         % convert time values to range values
frequency_values = frequency(row.');    % get frequency values

targetClusters = [range_values; frequency_values]; % store centroids as range and frequency

end
