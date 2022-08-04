function [particles,weights] = simpleResample(particles,weights)
%SIMPLE_RESAMPLE Summary of this function goes here
%   Detailed explanation goes here
    N = size(particles,1);
    cumulative_sum = cumsum(weights);
    cumulative_sum(1)= 0.0;
    cumulative_sum(end)= 1.0;  %avoid round-off error
    
    display(min(cumulative_sum));
    display(max(cumulative_sum));
    
    display(min(rand(N,1)));
    display(max(rand(N,1)));
    
    indexes = discretize(rand(N,1),cumulative_sum);
    
    %resample according to indexes
    particles(:,1) = particles(indexes,1);
    particles(:,2) = particles(indexes,2);
    weights(:) = (1.0/N);
        
end

