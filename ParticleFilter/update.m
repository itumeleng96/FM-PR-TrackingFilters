function [particles,weights] = update(particles,weights,z,meas_err)
%UPDATE Summary of this function goes here
%   Detailed explanation goes here
%   z : measured range and doppler shift
%   meas_err = measurement standard deviation error
    weights(:)= 1;
    
    %Get distance between particles and the measured values
    distance =  sqrt((particles(:,1)- z(1)).^2 + (particles(:,2)- z(2)).^2);
    %Get the shortest distance and weigh
    mean = min(distance);
    weights =  weights .* normpdf(distance,mean,meas_err);
    
    weights = weights + 1.e-300;
    weights = weights/(sum(weights));
    
    
end

