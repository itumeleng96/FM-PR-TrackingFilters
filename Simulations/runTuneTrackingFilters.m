%% Author : IJ Malemela 
%% Script for tuning the Filter parameters

clc;
clear;
close all;

addpath('FERS/', ...
        'cfar/', ...
        'meanShiftCluster/', ...
        'multiTargetTracking/', ...
        'DPI_Suppression', ...
        'TrackingFilter-KalmanFilter/', ...
        'TrackingFilter-ParticleFilter/', ...
        'TrackingFilter-UKF/', ... 
        'TrackingFilter-CSUKF/', ... 
        'TrackingFilter-RGNF/',...
        'TrackingFilter-CSRGNF/');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ARD
f=figure(1);
f.Position = [4000 10 1050 800]; 
movegui(f,'northwest');


rdm =[];
%{
f2=figure(2);
f2.Position = [4000 10 1050 800]; 
movegui(f,'southwest');

f3=figure(3);
f3.Position = [4000 10 1050 800]; 
movegui(f,'east');
%}
rangeTrueData = h5read('./true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('./true_data.h5', '/doppler_shifts');

rangeMeasurementData = h5read('./measurement_data.h5', '/bistatic_ranges');
dopplerMeasurementData = h5read('./measurement_data.h5', '/doppler_shifts');



figure(f1);
hold on; 
plot(rangeTrueData, dopplerTrueData, 'b', 'DisplayName', 'Ground Truth');
plot(rangeMeasurementData, dopplerMeasurementData, 'r--', 'DisplayName', 'Measurement Data');

title('Bistatic Range vs Doppler Shift');
xlabel('Bistatic Range (km)');
ylabel('Bistatic Doppler Shift (Hz)');

xlim([0 100]);
ylim([-200 200]);

legend('show');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bayesian Optimization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


