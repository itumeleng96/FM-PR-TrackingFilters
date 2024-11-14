clc;
clear;
close all;

addpath('../FERS/', ...
        '../cfar/', ...
        '../meanShiftCluster/', ...
        '../multiTargetTracking/', ...
        '../DPI_Suppression', ...
        '../TrackingFilter-KalmanFilter/', ...
        '../TrackingFilter-ParticleFilter/', ...
        '../TrackingFilter-UKF/', ... 
        '../TrackingFilter-CSUKF/', ... 
        '../TrackingFilter-RGNF/',...
        '../TrackingFilter-CSRGNF/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Data
rangeTrueData = h5read('../true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('../true_data.h5', '/doppler_shifts');
groundTruthData = [rangeTrueData; dopplerTrueData]; % Assuming this is your ground truth

rangeMeasurementData = h5read('../measurement_data.h5', '/bistatic_ranges');
dopplerMeasurementData = h5read('../measurement_data.h5', '/doppler_shifts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bayesian Optimization on the Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the search space for optimization (without log transformations)
optimVars = [
    optimizableVariable('std_acc1', [0.0001, 0.1]),   % Process noise std (for range)
    optimizableVariable('std_acc2', [0.0001, 0.1]),   % Process noise std (for doppler)
    optimizableVariable('std_meas1', [1, 10]),        % Measurement noise std (for range)
    optimizableVariable('std_meas2', [0.1, 1])        % Measurement noise std (for doppler)
];

% Objective function for Bayesian Optimization
objFun = @(params) ukf_nees_loss(params.std_acc1, params.std_acc2, params.std_meas1, params.std_meas2, rangeMeasurementData, dopplerMeasurementData, groundTruthData);

% Run Bayesian Optimization
results = bayesopt(objFun, optimVars, ...
                   'MaxObjectiveEvaluations', 30, ... % Number of evaluations
                   'IsObjectiveDeterministic', true, ...
                   'UseParallel', false, ...
                   'Verbose', 1, ...
                   'AcquisitionFunctionName', 'expected-improvement-plus', ...
                   'PlotFcn', {@plotAcquisitionFunction});

% Extract the optimal parameters
optimal_params = [results.XAtMinObjective.std_acc1, ...
                  results.XAtMinObjective.std_acc2, ...
                  results.XAtMinObjective.std_meas1, ...
                  results.XAtMinObjective.std_meas2];

disp('Optimized std_acc:');
disp(optimal_params(1:2));
disp('Optimized std_meas:');
disp(optimal_params(3:4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Definitions (must be at the end in MATLAB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function loss = ukf_nees_loss(std_acc1, std_acc2, std_meas1, std_meas2, rangeMeasurementData, dopplerMeasurementData, groundTruthData)
    % Initialize Kalman Filter with the current parameters
    std_acc = [std_acc1, std_acc2];  % Process noise std
    std_meas = [std_meas1, std_meas2]; % Measurement noise std
    
    x_initial = [rangeMeasurementData(1), dopplerMeasurementData(1)];
    dt = 1;
    UKF_object = unscentedKalmanFilter(dt,std_acc,std_meas(1),std_meas(2),[x_initial(1),0,x_initial(2),0;],0.01,0,2); %alpha,Kappa,Beta 

    % Initialize arrays for NEES calculation
    nees_values = [];

    % Loop through measurements and apply the Kalman filter
    for i = 1:58
        % Prediction
        [X_pred, UKF_object] = UKF_object.predict();
        
        % Update stage
        [X_est, UKF_object] = UKF_object.update([rangeMeasurementData(i); dopplerMeasurementData(i)]);
        disp(X_est);
        % Get the estimated state covariance (P) and estimation error
        P_k = UKF_object.P;  % State covariance matrix after the update
        position_estimate = X_est([1, 3]);  % Extract the range and doppler components
        
        % Compute the estimation error by comparing with the ground truth
        estimation_error = position_estimate - groundTruthData(:, i);        
        % Calculate NEES
        % Extract the variances (diagonal elements) for range and doppler
        P_range_variance = P_k(1, 1);  % Variance for the range estimate
        P_doppler_variance = P_k(3, 3);  % Variance for the doppler estimate
        
        % Calculate NEES for each component and combine
        NEES = (estimation_error(1)^2 / P_range_variance) + (estimation_error(2)^2 / P_doppler_variance);
       
        nees_values(end+1) = NEES;
    end

    % Total loss as the mean NEES over all time steps
    loss = mean(nees_values);
end
