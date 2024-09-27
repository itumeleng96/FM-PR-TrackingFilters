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



figure(f);
hold on; 
plot(rangeTrueData, dopplerTrueData, 'b', 'DisplayName', 'Ground Truth');
plot(rangeMeasurementData, dopplerMeasurementData, 'r--', 'DisplayName', 'Measurement Data');

title('Bistatic Range vs Doppler Shift');
xlabel('Bistatic Range (km)');
ylabel('Bistatic Doppler Shift (Hz)');

xlim([0 100]);
ylim([-200 200]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bayesian Optimization  On the Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Kalman filter
std_meas=[4,0.2];                           
std_acc=[0.001,0.02]; 
x_initial =[rangeMeasurementData(1),dopplerMeasurementData(1)];
dt = 1;                           
KF_object = kalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[x_initial(1);0;x_initial(2);0;]);
filterPredictionsDoppler = [];
filterPredictionsRange = [];

for i = 1:59
        %Prediction
        [X_pred,KF_object1] = KF_object.predict();
        filterPredictionsRange(end+1)= X_pred(1);
        filterPredictionsDoppler(end+1) = X_pred(3);

        %Update stage
        [X_est,KF_object2] = KF_object1.update([rangeMeasurementData(i);dopplerMeasurementData(i);]);

        KF_object = KF_object2;
end

plot(filterPredictionsRange, filterPredictionsDoppler, '^-', 'DisplayName', ' Predictions');

legend('show');
hold off;