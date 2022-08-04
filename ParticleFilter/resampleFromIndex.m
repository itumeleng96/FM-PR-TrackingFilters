function [particles,weights] = resampleFromIndex(particles,indexes)
%RESAMPLEFROMINDEX Summary of this function goes here
%   Detailed explanation goes here
    particles(:,1) = particles(indexes,1);
    particles(:,2) = particles(indexes,2);
    
    N = size(particles,1);
    weights = zeros(N,1);
    weights(:) = 1.0/size(weights,1);
end

