function [neff] = NEFF(weights)
%NEFF Summary of this function goes here
%   Detailed explanation goes here
    neff = 1./ sum(weights.^2) ;
end

