function [particles,weights] = simple_resample(particles,weights)
%SIMPLE_RESAMPLE Summary of this function goes here
%   Detailed explanation goes here
    N = size(particles,1);
    cumulative_sum = cumsum(weights);
    cumulative_sum(end)= 1;  %avoid round-off error
    indexes = discretize(rand(N,1),cumulative_sum);
    
    %resample according to indexes
    particles(:) = particles(indexes);
    weights(:) = (1.0/N);
        
end

