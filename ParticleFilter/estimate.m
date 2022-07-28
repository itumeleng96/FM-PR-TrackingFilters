function [mean,var] = estimate(particles,weights)
%ESTIMATE Summary of this function goes here
%   Detailed explanation goes here
    mean = [0,0];
    var = [0,0];
    
    mean(1) = sum(particles(:,1).*weights)/sum(weights);
    mean(2) = sum(particles(:,2).*weights)/sum(weights);
    
    var_particles = (particles - mean).^2;
    
    var(1) = sum(var_particles(:,1).*weights)/sum(weights);
    var(2) = sum(var_particles(:,2).*weights)/sum(weights);
end

