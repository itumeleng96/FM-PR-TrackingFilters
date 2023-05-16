%Function to do CA CFAR
%pf: Probability of False Alarm 
%N : Number of samples arounf the Cell Under Test
%G : Number of Guard cells
%cut : the magnitiude cells of the Range Doppler map

function [targetClusters,RDM,rdm_] = ca_cfarPlot(RDM,rate_fa,fs,fd_max,td_max,index,f,rdm)
%Generates a CFAR output map in the range-doppler domain
%Firstly calculate the interference power from the average of N samples in
%the vicinity of the CUT

%Parameters to allow zoom in
xlim_upper = 2e-4;
ylim_upper = 200;
ylim_lower = -200;

trc_num = 2;
guac_num = 4;

[rows, cols] = size(RDM);
RDM_final = zeros(rows, cols); 


for r=1:rows
    for c=1:cols
        % get guard cells sum value
        [guard_value, g_c] = subarea(RDM, c, r, guac_num);
        
        % get training cells sum value
        [train_value, t_c] = subarea(RDM, c, r, guac_num + trc_num);
        
        Num_train = t_c - g_c - 1;
        train_value = (train_value - guard_value) / (Num_train); %Noise Power
        train_value = pow2db(train_value);
       
        
        alpha = Num_train*((1-rate_fa)^(-1/Num_train) - 1);  %threshold factor
        threshold_val = train_value * alpha;
        %threshold_val = train_value * offset;
        
        if RDM(r,c) >= threshold_val
            RDM_final(r,c) = 1;
        end
    end
end

c=3e8;
Ndelay = floor(td_max*fs);                                  %number of points corresponding to td_max
time = 0:1/fs:Ndelay/fs;
range = time*c;
frequency = -fd_max:1:fd_max;


if index==1
    rdm = RDM_final ;
end

if index>1
    rdm =rdm+RDM_final ;
end

rdm_ = rdm;


figure(f);
imagesc(range, frequency, rdm_);
colormap(gca, 'gray'); % Set the colormap to 'gray'
caxis([min(rdm_(:)), max(rdm_(:))]); % Set the color axis to the minimum and maximum values of rdm_
c = colorbar;
c.Label.String = 'Intensity';
c.FontSize = 10;
c.Color = 'black'; % Set colorbar label color to white
set(gca, 'Color', 'black'); % Set background color to black
text(0, 10, "Time: " + index + "s", 'Color', 'white'); % Set text color to white
axis xy;
xlabel('Bistatic range [m]', 'FontSize', 10, 'Color', 'black'); % Set xlabel color to white
ylabel('Doppler frequency [Hz]', 'FontSize', 10, 'Color', 'black'); % Set ylabel color to white
grid on;
title('CFAR and Centroids', 'Color', 'black'); % Set title color to white
%xlim([0 xlim_upper]);
ylim([ylim_lower ylim_upper]);
RDM = RDM_final;

[row, column] = find(RDM > 0);
% Get target clusters as Bistatic Range and Doppler values

range_values = range(column.'); % convert time values to range values
frequency_values = frequency(row.'); % get frequency values

targetClusters = [range_values; frequency_values]; % store centroids as range and frequency

end
