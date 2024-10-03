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
simulation_time = 60;  % Fixed simulation time in seconds

% Parameters
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

num_simulations = 100;
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
        multiTargetTracker1 = multiTargetTracker1.predictionStage();

        % Prediction Stage
        tic;
        multiTargetTracker2 = multiTargetTracker2.predictionStage();
        pred_time_particle(i, sim) = toc; % Measure prediction time for Particle


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

% Calculate average processing times across simulations for each filter
avg_pred_time_kalman = mean(pred_time_kalman(:));
avg_update_time_kalman = mean(update_time_kalman(:));

avg_pred_time_particle = mean(pred_time_particle(:));
avg_update_time_particle = mean(update_time_particle(:));

avg_pred_time_ukf = mean(pred_time_ukf(:));
avg_update_time_ukf = mean(update_time_ukf(:));

avg_pred_time_rgnf = mean(pred_time_rgnf(:));
avg_update_time_rgnf = mean(update_time_rgnf(:));

% Print the relative computational load (average runtimes) for each filter
fprintf('Kalman Filter: Avg Prediction Time = %.7f s, Avg Update Time = %.5f s\n', avg_pred_time_kalman, avg_update_time_kalman);
fprintf('Particle Filter: Avg Prediction Time = %.7f s, Avg Update Time = %.5f s\n', avg_pred_time_particle, avg_update_time_particle);
fprintf('UKF: Avg Prediction Time = %.7f s, Avg Update Time = %.5f s\n', avg_pred_time_ukf, avg_update_time_ukf);
fprintf('RGNF: Avg Prediction Time = %.7f s, Avg Update Time = %.5f s\n', avg_pred_time_rgnf, avg_update_time_rgnf);

% Calculate relative computational load for prediction and update stages
relative_pred_particle = avg_pred_time_particle ./ avg_pred_time_kalman;
relative_update_particle = avg_update_time_particle ./ avg_update_time_kalman;

relative_pred_ukf = avg_pred_time_ukf ./ avg_pred_time_kalman;
relative_update_ukf = avg_update_time_ukf ./ avg_update_time_kalman;

relative_pred_rgnf = avg_pred_time_rgnf ./ avg_pred_time_kalman;
relative_update_rgnf = avg_update_time_rgnf ./ avg_update_time_kalman;

% Compute overall averages over the entire 60 seconds for each filter
avg_pred_kalman = mean(avg_pred_time_kalman);
avg_update_kalman = mean(avg_update_time_kalman);

avg_pred_particle = mean(avg_pred_time_particle);
avg_update_particle = mean(avg_update_time_particle);

avg_pred_ukf = mean(avg_pred_time_ukf);
avg_update_ukf = mean(avg_update_time_ukf);

avg_pred_rgnf = mean(avg_pred_time_rgnf);
avg_update_rgnf = mean(avg_update_time_rgnf);

% Compute relative computational load as per the Kalman filter baseline
relative_pred_particle_overall = avg_pred_particle / avg_pred_kalman;
relative_update_particle_overall = avg_update_particle / avg_update_kalman;

relative_pred_ukf_overall = avg_pred_ukf / avg_pred_kalman;
relative_update_ukf_overall = avg_update_ukf / avg_update_kalman;

relative_pred_rgnf_overall = avg_pred_rgnf / avg_pred_kalman;
relative_update_rgnf_overall = avg_update_rgnf / avg_update_kalman;

% Print results
fprintf('Relative computational load for Particle Filter (Prediction): %.7f\n', relative_pred_particle_overall);
fprintf('Relative computational load for Particle Filter (Update): %.7f\n', relative_update_particle_overall);

fprintf('Relative computational load for UKF (Prediction): %.7f\n', relative_pred_ukf_overall);
fprintf('Relative computational load for UKF (Update): %.7f\n', relative_update_ukf_overall);

fprintf('Relative computational load for RGNF (Prediction): %.7f\n', relative_pred_rgnf_overall);
fprintf('Relative computational load for RGNF (Update): %.7f\n', relative_update_rgnf_overall);

