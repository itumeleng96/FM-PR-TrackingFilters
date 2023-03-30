%Function to do CA CFAR
%pf: Probability of False Alarm 
%N : Number of samples arounf the Cell Under Test
%G : Number of Guard cells
%cut : the magnitiude cells of the Range Doppler map

function RDM = ca_cfar(RDM,rate_fa)
%Generates a CFAR output map in the range-doppler domain
%Firstly calculate the interference power from the average of N samples in
%the vicinity of the CUT

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



RDM = RDM_final;
