clc;
clear;
close all;

addpath('FERS/', ...
        'cfar/', ...
        'meanShiftCluster/', ...
        'multiTargetTracking/', ...
        'DPI_Suppression', ...
        'TrackingFilter-KalmanFilter/', ...
        'TrackingFilter-CSKF/', ...
        'TrackingFilter-ParticleFilter/', ...
        'TrackingFilter-CS-ParticleFilter/', ...
        'TrackingFilter-UKF/', ... 
        'TrackingFilter-CSUKF/',...
        'TrackingFilter-RGNF/',...
        'TrackingFilter-CSRGNF/');

%FLIGHT Scenarios
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_1_laneChange.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_2_landingManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_4_360.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_5_2_targets.fersxml');

%Noise Scenarios
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_1_fm_noise.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_2_white_noise.fersxml');

system('wsl fers FERS/BackupScenarios/scenario_1_singleFile.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/BackupScenarios/scenario_1_singleFile_120.fersxml');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h5 Import from FERS simulation
[Ino, Qno, scale_no] = loadfersHDF5('direct.h5');
[Imov, Qmov, scale_mov] = loadfersHDF5('echo.h5');

I_Qmov = Imov + 1i*Qmov;
I_Qmov = I_Qmov.*scale_mov;
I_Qno = Ino + 1i*Qno;
I_Qno = I_Qno.*scale_no;

%I_Qmov=I_Qmov-I_Qno;

fs = 200000;
dopp_bins = 200;
delay = 333e-6;
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



s1 = I_Qmov;   %Surv
s2 = I_Qno;    %Ref

initial=1;

current=fs;                                %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs;       %Simulation time: number of data points/sampling frequency

ard = [];
rdm =[];


f3=figure();
f3.Position = [4000 10 1000 800]; 
movegui(f3,'southwest');

%f5=figure();
%f5.Position = [4000 10 1000 800]; 
%movegui(f5,'southeast');

%Create MTT object
confirmationThreshold=4;
deletionThreshold=4;
gatingThreshold=15; %scalar value for ellipsoidal gate

%FilterType 1: Kalman Filter
%FilterType 2: Huber Covariance Scaling Kalman Filter
%FilterType 3: Particle Filter
%FilterType 4: Particle Filter
%FilterType 5: UKF  Filter
%FilterType 6: Covariance Scaling UKF  Filter
%FilterType 7: RGNF Filter
%FilterType 8: Covariance scaling RGNF Filter


multiTargetTracker1 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,5);
multiTargetTracker2 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,9);
multiTargetTracker3 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,10);
multiTargetTracker4 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,11);



prevCentroids=[];

%True Data for single Target Scenario
rangeTrueData = h5read('./groundTruthCalculations/true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('./groundTruthCalculations/true_data.h5', '/doppler_shifts');

RDM =[];

for i = 1:simulation_time
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA_Optimized(s2,s1,proc);

    %Range-Doppler Map
    [y,ard_] = ardNoPlot(s1,s2,fs,dopp_bins,delay,i,ard);
    
    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters, RDM, rdm_] = ca_cfar(y.', 10e-9, fs, dopp_bins, delay, 20);

    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,variancesX,variancesY,numPoints] = meanShift(targetClusters,10,8);
   
    %Plot tracks from Tracker - Call Multi-target Tracker

    multiTargetTracker1 = multiTargetTracker1.createNewTracks(clusterCentroids,i);
    multiTargetTracker2 = multiTargetTracker2.createNewTracks(clusterCentroids,i);
    multiTargetTracker3 = multiTargetTracker3.createNewTracks(clusterCentroids,i);
    multiTargetTracker4 = multiTargetTracker4.createNewTracks(clusterCentroids,i);


    %DELETE and CONFIRM Tracks
    multiTargetTracker1 = multiTargetTracker1.maintainTracks();
    multiTargetTracker2 = multiTargetTracker2.maintainTracks();
    multiTargetTracker3 = multiTargetTracker3.maintainTracks();
    multiTargetTracker4 = multiTargetTracker4.maintainTracks();

    %Filter Prediction Stage
    multiTargetTracker1 = multiTargetTracker1.predictionStage();
    multiTargetTracker2 = multiTargetTracker2.predictionStage();
    multiTargetTracker3 = multiTargetTracker3.predictionStage();
    multiTargetTracker4 = multiTargetTracker4.predictionStage();

    multiTargetTracker1 = multiTargetTracker1.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM);
    %multiTargetTracker2 = multiTargetTracker2.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM);
    %multiTargetTracker3 = multiTargetTracker3.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM);
    %multiTargetTracker4 = multiTargetTracker4.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM);

    %UPDATE Tracks from measurements
    multiTargetTracker1 = multiTargetTracker1.updateStage(clusterCentroids,i);
    multiTargetTracker2 = multiTargetTracker2.updateStage(clusterCentroids,i);
    multiTargetTracker3 = multiTargetTracker3.updateStage(clusterCentroids,i);
    multiTargetTracker4 = multiTargetTracker4.updateStage(clusterCentroids,i);

    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;

end

trackId = input('Enter a trackId : ');

%Multi-Target Tracking 

f2=figure();
f2.Position = [4000 10 1000 800]; 
movegui(f2,'southeast');

%%PLOT TRACKID TRACK FOR DIFFERENT FILTERS
track_mtt_1 = multiTargetTracker1.getTrack(trackId);
track_mtt_2 = multiTargetTracker2.getTrack(trackId);
track_mtt_3 = multiTargetTracker3.getTrack(trackId);
track_mtt_4 = multiTargetTracker4.getTrack(trackId);


track_mtt_true = multiTargetTracker1.getMeasuredTrack(trackId);

%
% Plot Predicted Tracks from Different filters against Ground Truth
figure(f2);

c = 3e8;
Ndelay = floor(delay * fs);                                 
time = 0:1/fs:Ndelay/fs;
range = time * c;
range = range / 1000;
frequency = -dopp_bins:1:dopp_bins;

imagesc(range, frequency, RDM * 0);
colormap(gca, 'white');

hold on;
% Plot different tracks with distinct colors and increased line width
plot(track_mtt_1(1,:), track_mtt_1(2,:), '-', 'Color', 'c', 'LineWidth', 1.4); % Black for track 1
hold on;
plot(track_mtt_2(1,:), track_mtt_2(2,:), '-', 'Color', 'g', 'LineWidth', 1.4); % Green for track 2
hold on;
plot(track_mtt_3(1,:), track_mtt_3(2,:), '-', 'Color', 'b', 'LineWidth', 1.4); % Blue for track 3
hold on;
plot(track_mtt_4(1,:), track_mtt_4(2,:), '-', 'Color', 'm', 'LineWidth', 1.4); % Magenta for track 4
hold on;
plot(track_mtt_true(1,:), track_mtt_true(2,:), 'o', 'Color', 'r', 'LineWidth', 1.4); % Red for true track
hold on;
plot(rangeTrueData, dopplerTrueData, '-', 'Color', 'k', 'LineWidth', 1.4); % Cyan for ground truth

% Set axis limits
ylim([-190 190]);
xlim([10 45]);

% Title, labels, and legend
title('Unscented Kalman Filter with varying Alpha','Fontsize', 18);
axis xy;
xlabel('Bistatic Range [km]', 'Fontsize', 18);
ylabel('Bistatic Doppler frequency [Hz]', 'Fontsize', 18);
legend('alpha:0.0001', 'alpha:0.001', 'alpha:0.01','alpha:0.1', 'Measurement', 'Ground Truth');

% Enable grid
grid on;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RGNF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Predicted Tracks from Different filters against Ground Truth
%{
figure(f2);

c = 3e8;
Ndelay = floor(delay * fs);                                 
time = 0:1/fs:Ndelay/fs;
range = time * c;
range = range / 1000;
frequency = -dopp_bins:1:dopp_bins;

imagesc(range, frequency, RDM * 0);
colormap(gca, 'white');

hold on;
% Plot different tracks with distinct colors and increased line width
plot(track_mtt_1(1,:), track_mtt_1(2,:), '-', 'Color', 'c', 'LineWidth', 1.4); % Black for track 1
hold on;
plot(track_mtt_2(1,:), track_mtt_2(2,:), '-', 'Color', 'g', 'LineWidth', 1.4); % Green for track 2
hold on;
plot(track_mtt_3(1,:), track_mtt_3(2,:), '-', 'Color', 'b', 'LineWidth', 1.4); % Blue for track 3
hold on;
plot(track_mtt_4(1,:), track_mtt_4(2,:), '-', 'Color', 'm', 'LineWidth', 1.4); % Magenta for track 4
hold on;
plot(track_mtt_true(1,:), track_mtt_true(2,:), 'o', 'Color', 'r', 'LineWidth', 1.4); % Red for true track
hold on;
plot(rangeTrueData, dopplerTrueData, '-', 'Color', 'k', 'LineWidth', 1.4); % Cyan for ground truth

% Set axis limits
ylim([-190 190]);
xlim([10 45]);

% Title, labels, and legend
title('RGNF Filter with varying lambda');
axis xy;
xlabel('Bistatic Range [km]', 'Fontsize', 18);
ylabel('Bistatic Doppler frequency [Hz]', 'Fontsize', 18);
legend('lambda:1', 'lambda:0.25', 'lambda:0.5','lambda:0.75', 'Measurement', 'Ground Truth');

% Enable grid
grid on;
%}
%{

%%% Kalman Filter setups

% Plot Predicted Tracks from Different filters against Ground Truth
figure(f2);

c = 3e8;
Ndelay = floor(delay * fs);                                 
time = 0:1/fs:Ndelay/fs;
range = time * c;
range = range / 1000;
frequency = -dopp_bins:1:dopp_bins;

imagesc(range, frequency, RDM * 0);
colormap(gca, 'white');

hold on;
% Plot different tracks with distinct colors and increased line width
plot(track_mtt_1(1,:), track_mtt_1(2,:), '-', 'Color', 'c', 'LineWidth', 1.4); % Black for track 1
hold on;
plot(track_mtt_2(1,:), track_mtt_2(2,:), '-', 'Color', 'g', 'LineWidth', 1.4); % Green for track 2
hold on;
plot(track_mtt_3(1,:), track_mtt_3(2,:), '-', 'Color', 'b', 'LineWidth', 1.4); % Blue for track 3
hold on;
plot(track_mtt_4(1,:), track_mtt_4(2,:), '-', 'Color', 'm', 'LineWidth', 1.4); % Magenta for track 4
hold on;
plot(track_mtt_true(1,:), track_mtt_true(2,:), 'o', 'Color', 'r', 'LineWidth', 1.4); % Red for true track
hold on;
plot(rangeTrueData, dopplerTrueData, '-', 'Color', 'k', 'LineWidth', 1.4); % Cyan for ground truth

% Set axis limits
ylim([-190 190]);
xlim([10 45]);

% Title, labels, and legend
title('Kalman Filter with varying Covariance Matrix (Q)','FontSize', 20);
axis xy;
xlabel('Bistatic Range [km]', 'Fontsize', 18);
ylabel('Bistatic Doppler frequency [Hz]', 'Fontsize', 18);
legend('s_{w1}: 0.01', 's_{w1}: 0.1 ', 's_{w2}:0.01', 's_{w2}:0.1', 'Measurement', 'Ground Truth');

% Enable grid
grid on;

%% Particle Filter

figure(f2);

c = 3e8;
Ndelay = floor(delay * fs);                                 
time = 0:1/fs:Ndelay/fs;
range = time * c;
range = range / 1000;
frequency = -dopp_bins:1:dopp_bins;

imagesc(range, frequency, RDM * 0);
colormap(gca, 'white');

hold on;
% Plot different tracks with distinct colors and increased line width
plot(track_mtt_1(1,:), track_mtt_1(2,:), '-', 'Color', 'c', 'LineWidth', 1.4); % Black for track 1
hold on;
plot(track_mtt_2(1,:), track_mtt_2(2,:), '-', 'Color', 'g', 'LineWidth', 1.4); % Green for track 2
hold on;
plot(track_mtt_3(1,:), track_mtt_3(2,:), '-', 'Color', 'b', 'LineWidth', 1.4); % Blue for track 3
hold on;
plot(track_mtt_4(1,:), track_mtt_4(2,:), '-', 'Color', 'm', 'LineWidth', 1.4); % Magenta for track 4
hold on;
plot(track_mtt_true(1,:), track_mtt_true(2,:), 'o', 'Color', 'r', 'LineWidth', 1.4); % Red for true track
hold on;
plot(rangeTrueData, dopplerTrueData, '-', 'Color', 'k', 'LineWidth', 1.4); % Cyan for ground truth

% Set axis limits
ylim([-190 190]);
xlim([10 45]);

% Title, labels, and legend
title('Particle Filter with varying number of particles','Fontsize', 18);
axis xy;
xlabel('Bistatic Range [km]', 'Fontsize', 18);
ylabel('Bistatic Doppler frequency [Hz]', 'Fontsize', 18);
legend('Number of particles:2500', 'Number of particles:5000', 'Number of particles:7500','Number of Particles:10000', 'Measurement', 'Ground Truth');
grid on;
%}