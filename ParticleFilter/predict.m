function [particles] = predict(particles,U,std,dt)
%PREDICT move particles according to the derivate of range and the doppler
%shift according to the input U
%   std : vector with standard deviation
%   U : input acceleration
%   dt: sampling time
    N = size(particles,1);
    particles_ = particles;
    %update velocity
    particles(:,2) = particles(:,2) + U(2)*dt +randn(N,1)*std(2); 
    %update range
    particles(:,1) = particles(:,1) + particles_(:,2)*dt + U(1)*((dt^2)/2) +randn(N,1)*std(1);
    
    
end

