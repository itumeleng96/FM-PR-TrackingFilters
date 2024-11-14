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
system('wsl fers FERS/flightScenarios/scenario_2_landingManeuver.fersxml');
%system('wsl fers FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_4_360.fersxml');
%system('wsl fers FERS/flightScenarios/scenario_5_2_targets.fersxml');
%system('wsl fers FERS/flightScenarios/scenario_5_3_targets.fersxml');

%Noise Scenarios
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_1_fm_noise.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_2_white_noise.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/BackupScenarios/scenario_1_singleFile.fersxml');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ARD
f=figure(1);
f.Position = [4000 10 1050 800]; 
movegui(f,'northwest');

%CFAR
f2=figure(2);
f2.Position = [4000 10 1050 800]; 
movegui(f2,'southwest');

%Multi-Target Tracking 
f3=figure(3);
f3.Position = [4000 10 1050 800]; 
movegui(f3,'southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create MTT object
confirmationThreshold=4;
deletionThreshold=6;
gatingThreshold=15; %scalar value for ellipsoidal gate

%FilterType 1: Kalman Filter
%FilterType 2: Huber Covariance Scaling Kalman Filter
%FilterType 3: Particle Filter
%FilterType 4: UKF  Filter
%FilterType 5: Covariance Scaling UKF  Filter
%FilterType 6: RGNF Filter
%FilterType 7: Covariance scaling RGNF Filter





%filterType =input('Tracking Filter to use (1-7):');
filterType =8;

multiTargetTracker = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType);

%LOG_LIKELIHOODS
doppler_ll=[];
range_ll=[];

prevCentroids=[];

rangeTrueData = h5read('./groundTruthCalculations/true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('./groundTruthCalculations/true_data.h5', '/doppler_shifts');

for i = 1:simulation_time
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    %s1 = procECA(s2,s1,proc);
    s1 = procECA_Optimized(s2,s1,proc);

    %Plot Range-Doppler Map
    [y,ard_] = ardNoPlot(s1,s2,fs,dopp_bins,delay,i,ard);

    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfar(y.',10^-7,fs,dopp_bins,delay,25);                    
    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,variancesX,variancesY,numPoints] = meanShift(targetClusters,10,10);
    %[clusterCentroids,prevCentroids,variancesX,variancesY,numPoints] = meanShiftPlot(targetClusters,1e4,8,prevCentroids);
    
    %Plot tracks from Tracker - Call Multi-target Tracker
    multiTargetTracker = multiTargetTracker.createNewTracks(clusterCentroids,i);
    
    %DELETE and CONFIRM Tracks
    multiTargetTracker = multiTargetTracker.maintainTracks();

    %Filter Prediction Stage
    multiTargetTracker = multiTargetTracker.predictionStage();

   
    %PLOT Prediction and True Tracks
    multiTargetTracker.plotMultiTar
    getTrackingId(fs,dopp_bins,delay,i,f3,RDM);
    %UPDATE Tracks from measurements
    multiTargetTracker = multiTargetTracker.updateStage(clusterCentroids,i);
    
    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;
end





%Do Log-likelihood for specific TrackId after simulation
trackId = input('Enter a trackId for the Log-likelihood: ');
%trackId =1;

doppler_ll=[];
range_ll=[];

%%PLOT TRACKID TRACK FOR DIFFERENT FILTERS\
track_mtt_1 = multiTargetTracker.getTrack(trackId);
track_mtt_1_true = multiTargetTracker.getMeasuredTrack(trackId);

figure(3);
plot(track_mtt_1(1,:),track_mtt_1(2,:), 'b-');
hold on;
plot(track_mtt_1_true(1,:),track_mtt_1_true(2,:), '-^');
hold on;
plot(rangeTrueData,dopplerTrueData,'-*');
xlabel('Bistatic range (KM)');
ylabel('Bistatic Doppler (Hz)');
title(['Tracking filter outputs vs Ground Truth For Track:', num2str(trackId)]);

%Call MultiTargetTrack -LogLikelihood to plot 
[doppler_ll,range_ll,t1]=multiTargetTracker.calculateLogLikelihoodGroundTruthP(trackId,doppler_ll,range_ll,dopplerTrueData,rangeTrueData,simulation_time);

figure(4);
plot(t1,doppler_ll, 'b-'); 
title(['Bistatic Doppler Log-Likelihood Comparison for Track:',num2str(trackId)]);
xlabel('Time Steps');
ylabel('Bistatic Doppler Log-Likelihood');


figure(5);
plot(t1,range_ll, 'b-');
title(['Bistatic Range Log-Likelihood Comparison for Track:',num2str(trackId)]);
xlabel('Time Steps');
ylabel('Bistatic Range Log-Likelihood');