clc;
clear;
close all;

addpath('FERS/', ...
        'cfar/', ...
        'meanShiftCluster/', ...
        'multiTargetTracking/', ...
        'DPI_Suppression/', ...
        'TrackingFilter-KalmanFilter/', ...
        'TrackingFilter-CSKF/', ...
        'TrackingFilter-ParticleFilter/', ...
        'TrackingFilter-CS-ParticleFilter/', ...
        'TrackingFilter-UKF/', ... 
        'TrackingFilter-CSUKF/', ... 
        'TrackingFilter-RGNF/',...
        'TrackingFilter-CSRGNF/');

% Load simulation data
system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/BackupScenarios/scenario_1_singleFile.fersxml'); 
[Ino, Qno, scale_no] = loadfersHDF5('direct.h5');
[Imov, Qmov, scale_mov] = loadfersHDF5('echo.h5');

I_Qmov = Imov + 1i*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + 1i*Qno;
I_Qno = I_Qno.*scale_no;

fs = 200000;
current = fs;                                % based on samples in transmitted signal
simulation_time = size(I_Qmov, 1) / fs;       % Simulation time: number of data points/sampling frequency

dopp_bins = 200;
delay = 233e-6;
c=299792458;
range_delay = delay*c;


%DPI Cancellation
proc = struct('cancellationMaxRange_m', range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m', 12540, ...
              'nSegments', 1, ...
              'nIterations', 20, ...
              'Fs', fs, ...
              'alpha', 0, ...
              'initialAlpha', 0);

% Create MTT objects for different filters
confirmationThreshold = 4;
deletionThreshold = 4;
gatingThreshold = 20;

multiTargetTracker1 = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 1); % Kalman
multiTargetTracker2 = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 3); % Particle
multiTargetTracker3 = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 5); % UKF
multiTargetTracker4 = multiTargetTracker(confirmationThreshold, deletionThreshold, gatingThreshold, 7); % RGNF

num_simulations = 1000;
% Initialize arrays to store processing times for each second
pred_time_kalman = zeros(simulation_time, num_simulations);
update_time_kalman = zeros(simulation_time, num_simulations);
pred_time_particle = zeros(simulation_time, num_simulations);
update_time_particle = zeros(simulation_time, num_simulations);
pred_time_ukf = zeros(simulation_time, num_simulations);
update_time_ukf = zeros(simulation_time, num_simulations);
pred_time_rgnf = zeros(simulation_time, num_simulations);
update_time_rgnf = zeros(simulation_time, num_simulations);

% Perform Monte Carlo simulations
for sim = 1:num_simulations
    initial = 1;  % Reset initial index for each simulation
    current = fs; % Reset current index for each simulation

    for i = 1:simulation_time
        s1 = I_Qmov(initial:current); % surv
        s2 = I_Qno(initial:current);  % ref
        s1 = procECA_Optimized(s2, s1, proc);  % DPI Cancellation

        % Range-Doppler Map
        [y, ~] = ardNoPlot(s1,s2,fs,dopp_bins,delay,i,[]);

        % Plot CFAR from Cell-Averaging CFAR 
        [targetClusters,~,~] = ca_cfar(y.', 10^-7, fs, dopp_bins, delay, 20);    

        [clusterCentroids,~,~,~] = meanShift(targetClusters, 10, 8);
        
        multiTargetTracker2 = multiTargetTracker2.predictionStage();


        tic;
        multiTargetTracker2 = multiTargetTracker2.predictionStage();
        pred_time_particle(i, sim) = toc; % Measure prediction time for Particle


        % Prediction Stage
        tic;
        multiTargetTracker1 = multiTargetTracker1.predictionStage();
        pred_time_kalman(i, sim) = toc; % Measure prediction time for Kalman


        tic;
        multiTargetTracker3 = multiTargetTracker3.predictionStage();
        pred_time_ukf(i, sim) = toc; % Measure prediction time for UKF

        tic;
        multiTargetTracker4 = multiTargetTracker4.predictionStage();
        pred_time_rgnf(i, sim) = toc; % Measure prediction time for RGNF

        % Update Stage
        multiTargetTracker2 = multiTargetTracker2.updateStage(clusterCentroids, i);

        tic;
        multiTargetTracker2 = multiTargetTracker2.updateStage(clusterCentroids, i);
        update_time_particle(i, sim) = toc; % Measure update time for Particle


        tic;
        multiTargetTracker1 = multiTargetTracker1.updateStage(clusterCentroids, i);
        update_time_kalman(i, sim) = toc; % Measure update time for Kalman

        
        tic;
        multiTargetTracker3 = multiTargetTracker3.updateStage(clusterCentroids, i);
        update_time_ukf(i, sim) = toc; % Measure update time for UKF

        tic;
        multiTargetTracker4 = multiTargetTracker4.updateStage(clusterCentroids, i);
        update_time_rgnf(i, sim) = toc; % Measure update time for RGNF

        % Update indices for next iteration
        initial = current + 1;
        current = current + fs;
    end
end
% Assuming simulation_time is 60 seconds
simulation_time = 60;

% Calculate average processing times across simulations for each filter at every second
avg_pred_time_kalman = mean(pred_time_kalman, 2);
avg_update_time_kalman = mean(update_time_kalman, 2);
avg_pred_time_particle = mean(pred_time_particle, 2);
avg_update_time_particle = mean(update_time_particle, 2);
avg_pred_time_ukf = mean(pred_time_ukf, 2);
avg_update_time_ukf = mean(update_time_ukf, 2);
avg_pred_time_rgnf = mean(pred_time_rgnf, 2);
avg_update_time_rgnf = mean(update_time_rgnf, 2);

% Time axis based on the number of seconds in the simulation
time = 1:simulation_time;

% Plotting the line graph for prediction times across 60 seconds
figure;
plot(time, avg_pred_time_kalman, 'g', 'LineWidth', 0.5); hold on;
plot(time, avg_pred_time_particle, 'b', 'LineWidth', 0.5);
plot(time, avg_pred_time_ukf, 'r', 'LineWidth', 0.5);
plot(time, avg_pred_time_rgnf, 'k', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Average Prediction Time (s)');
title('Average Prediction Time per Second for Each Filter');
legend('Kalman filter', 'Particle filter', 'UKF', 'RGNF', 'Location', 'best');
grid on;
hold off;

% Plotting the line graph for update times across 60 seconds
figure;
plot(time, avg_update_time_kalman, 'g', 'LineWidth', 0.5); hold on;
plot(time, avg_update_time_particle, 'b', 'LineWidth', 0.5);
plot(time, avg_update_time_ukf, 'r', 'LineWidth', 0.5);
plot(time, avg_update_time_rgnf, 'k', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Average Update Time (s)');
title('Average Update Time per Second for Each Filter');
legend('Kalman filter', 'Particle filter', 'UKF', 'RGNF', 'Location', 'best');
grid on;
hold off;