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
        'TrackingFilter-CSKF/', ...
        'TrackingFilter-CSUKF/', ... 
        'TrackingFilter-CS-ParticleFilter/', ... 
        'TrackingFilter-RGNF/',...
        'TrackingFilter-CSRGNF/');

%FLIGHT Scenarios
%system('wsl fers FERS/flightScenarios/scenario_1_laneChange.fersxml');
%system('wsl fers FERS/flightScenarios/scenario_2_landingManeuver.fersxml');
%system('wsl fers FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml');
%system('wsl fers FERS/flightScenarios/scenario_4_360.fersxml');
%system('wsl fers FERS/flightScenarios/scenario_5_2_targets.fersxml');

%Noise Scenarios
%SCENARIO 3
%system('wsl fers FERS/NoiseScenarios/scenario_1_fm_noise.fersxml');

%SCENARIO 2
%system('wsl fers FERS/NoiseScenarios/scenario_2_white_noise.fersxml');

%SCENARIO 1 
system('wsl fers FERS/BackupScenarios/scenario_1_singleFile.fersxml'); 



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



s1 = I_Qmov;   %Surv
s2 = I_Qno;    %Ref

initial=1;

current=fs;                                %based on samples in transmitted signal
simulation_time = size(I_Qmov,1)/fs;       %Simulation time: number of data points/sampling frequency


ard = [];
rdm =[];


f3=figure();
f3.Position = [4000 10 800 600]; 
movegui(f3,'southwest');

f5=figure();
f5.Position = [4000 10 800 600]; 
movegui(f5,'southeast');

%Create MTT object
confirmationThreshold=3;
deletionThreshold=5;
gatingThreshold=20;         %scalar value for ellipsoidal gate

%FilterType 1: Kalman Filter
%FilterType 2: Huber Covariance Scaling Kalman Filter
%FilterType 3: Particle Filter
%FilterType 4: UKF  Filter
%FilterType 5: Covariance Scaling UKF  Filter
%FilterType 6: RGNF Filter
%FilterType 7: Covariance scaling RGNF Filter


multiTargetTracker1 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,2);
multiTargetTracker2 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,4);
multiTargetTracker3 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,6);
multiTargetTracker4 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,8);

%LOG_LIKELIHOODS
doppler_ll_1=[];
range_ll_1=[];

doppler_ll_2=[];
range_ll_2=[];


doppler_ll_3=[];
range_ll_3=[];


doppler_ll_4=[];
range_ll_4=[];

prevCentroids=[];

%True Data for single Target Scenario
rangeTrueData = h5read('./groundTruthCalculations/true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('./groundTruthCalculations/true_data.h5', '/doppler_shifts');

for i = 1:simulation_time
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA_Optimized(s2,s1,proc);

    %Range-Doppler Map
    [y,ard_] = ardNoPlot(s1,s2,fs,dopp_bins,delay,i,ard);
    
    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfarPlotBW(y.',10e-6,fs,dopp_bins,delay,i,f5,rdm,25);                    

    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,prevCentroids,variancesX,variancesY,numPoints] = meanShiftPlot(targetClusters,10,8,prevCentroids);
   
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

trackId = input('Enter a trackId for the Log-likelihood: ');

%Multi-Target Tracking 
f=figure();
f.Position = [4000 10 1000 800]; 
movegui(f,'southeast');

f1=figure();
f1.Position = [4000 10 1000 800]; 
movegui(f1,'southwest');

f2=figure();
f2.Position = [4000 10 1000 800]; 
movegui(f2,'southeast');

%%PLOT TRACKID TRACK FOR DIFFERENT FILTERS
track_mtt_1 = multiTargetTracker1.getTrack(trackId);
track_mtt_2 = multiTargetTracker2.getTrack(trackId);
track_mtt_3 = multiTargetTracker3.getTrack(trackId);
track_mtt_4 = multiTargetTracker4.getTrack(trackId);


%PLOT TRACKID LOGLIKELIHOOD FOR DIFFERENT FILTERS
[doppler_ll_1,range_ll_1,t1]=multiTargetTracker1.calculateLogLikelihoodGroundTruth(trackId,doppler_ll_1,range_ll_1,dopplerTrueData,rangeTrueData,simulation_time);
[doppler_ll_2,range_ll_2,t2]=multiTargetTracker2.calculateLogLikelihoodGroundTruth(trackId,doppler_ll_2,range_ll_2,dopplerTrueData,rangeTrueData,simulation_time);
[doppler_ll_3,range_ll_3,t3]=multiTargetTracker3.calculateLogLikelihoodGroundTruth(trackId,doppler_ll_3,range_ll_3,dopplerTrueData,rangeTrueData,simulation_time);
[doppler_ll_4,range_ll_4,t4]=multiTargetTracker4.calculateLogLikelihoodGroundTruth(trackId,doppler_ll_4,range_ll_4,dopplerTrueData,rangeTrueData,simulation_time);

% Plot Predicted Tracks from Different filters against Ground Truth
figure(f2);
plot(track_mtt_1(1,:),track_mtt_1(2,:), 'b-'); 
hold on;
plot(track_mtt_2(1,:),track_mtt_2(2,:), 'r--');
hold on;
plot(track_mtt_3(1,:),track_mtt_3(2,:), 'go-');
hold on;
plot(track_mtt_4(1,:),track_mtt_4(2,:), 'ms-.');
hold on;
plot(rangeTrueData,dopplerTrueData,'-*');

title(['Tracking filter outputs vs Ground Truth For Track:', num2str(trackId)]);
xlabel('Bistatic range (KM)');
ylabel('Bistatic Doppler (Hz)');
legend('Kalman Filter', 'Particle Filter','Unscented Kalman Filter', 'Recursive Gauss Newton Filter','Ground Truth');
grid on;


% Create comparison plots for Doppler Log-Likelihoods
figure(f);
plot(t1,doppler_ll_1, 'b-'); 
hold on;
plot(t2,doppler_ll_2, 'r--');
hold on;
plot(t3,doppler_ll_3, 'go-');
hold on;
plot(t4,doppler_ll_4, 'ms-.');

%title(['Bistatic Doppler Log-Likelihood Comparison for Track:',num2str(trackId)]);
title('Bistatic Doppler Log-Likelihood Comparison');

xlabel('Time(s)');
ylabel('Bistatic Doppler Log-Likelihood');
legend('Kalman Filter', 'Particle Filter','Unscented Kalman Filter', 'Recursive Gauss Newton Filter');
grid on;

% Create comparison plots for Range Log-Likelihoods
figure(f1);
plot(t1,range_ll_1, 'b-');
hold on;
plot(t2,range_ll_2,  'r--');
hold on;
plot(t3,range_ll_3, 'go-');
hold on;
plot(t4,range_ll_4, 'ms-.');

title('Bistatic Range Log-Likelihood Comparison');
xlabel('Time(s)');
ylabel('BIstatic Range Log-Likelihood');
legend('Kalman Filter', 'Particle Filter','Unscented Kalman Filter', 'Recursive Gauss Newton Filter');
grid on;


