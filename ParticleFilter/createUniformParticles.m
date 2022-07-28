function [particles] = createUniformParticles(x_range,y_range,N)
%Create a uniform Distribution of particles over a region
% N : number of particles
    particles = zeros(N,2);
    particles(:,1) = unifrnd(x_range(1),x_range(2),[N 1]);
    particles(:,2) = unifrnd(y_range(1),y_range(2),[N 1]);
end

