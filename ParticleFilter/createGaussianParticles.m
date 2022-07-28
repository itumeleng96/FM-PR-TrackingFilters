function [particles] = createGaussianParticles(mean,std,N)
%Create a uniform Distribution of particles over a region
% N : number of particles
    particles = zeros(N,2);
    particles(:,1) = mean(1) + (randn(N,1))*std(1) ;
    particles(:,2) = mean(2) + (randn(N,1))*std(2) ;
    
end


