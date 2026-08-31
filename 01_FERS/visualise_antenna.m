%% For visualising sinc antenna pattern
% G(thetha) = alpha * ( sin(beta .* theta) / (beta .* theta)) ^ gamma

clear;
addpath('../FERS/');

%currentCharacterEncoding = slCharacterEncoding();
%slCharacterEncoding('ISO-8859-1');

% Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 5.2481; %5.2481; %7.2 dBi
beta = 2;
gamma = 3.6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = -pi:0.01:pi;

G = alpha * ( sin(beta .* theta) ./ (beta .* theta)) .^ gamma;

polarDb(theta, G, 40, 2.5, 24);