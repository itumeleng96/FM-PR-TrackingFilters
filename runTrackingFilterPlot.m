%Author Itumeleng Malemela 
% This script is used for running a single MTT and tracing filter and
% assess their performance
clc;
clear;
close all;


%FLIGHT Scenarios
%system('fers FERS/flightScenarios/scenario_1_laneChange.fersxml');
%system('fers FERS/flightScenarios/scenario_2_landingManeuver.fersxml');
%system('fers FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml');
%system('fers FERS/flightScenarios/scenario_4_360.fersxml');
%system('fers FERS/flightScenarios/scenario_5_2_targets.fersxml');
%system('fers FERS/flightScenarios/scenario_5_3_targets.fersxml');

%Noise Scenarios
%system('fers FERS/NoiseScenarios/scenario_1_fm_noise.fersxml');
%system('fers FERS/NoiseScenarios/scenario_2_white_noise.fersxml');

system('fers FERS/BackupScenarios/scenario_1_singleFile.fersxml');
%system('fers FERS/BackupScenarios/scenario_1_singleFile_120.fersxml');
%system('fers FERS/BackupScenarios/scenario_3_targets_singleFile.fersxml');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','2D image');
%ARD
f=figure(1);
f.Position = [4000 10 1050 800]; 
movegui(f,'northwest');

%CFAR
f2=figure(2);
f2.Position = [4000 10 1050 800]; 
movegui(f2,'northwest');

%Multi-Target Tracking 
f3=figure(3);
f3.Position = [4000 10 1050 800];
movegui(f3,'southeast');

f4=figure(4);
f4.Position = [4000 10 1050 800]; 
movegui(f4,'southwest');


f5=figure(5);
f5.Position = [4000 10 1050 800]; 
movegui(f5,'southeast');
%{
f6=figure(6);
f4.Position = [4000 10 1050 800]; 
movegui(f6,'southeast');

f7=figure(7);
f4.Position = [4000 10 1050 800]; 
movegui(f7,'southeast');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create MTT object
confirmationThreshold=4;
deletionThreshold=6;
gatingThreshold=15;         %scalar value for ellipsoidal gate

%FilterType 1: Kalman Filter
%FilterType 2: Covariance Scaling Kalman Filter
%FilterType 3: Particle Filter
%FilterType 4: Covariance Scaling Particle Filter
%FilterType 5: UKF  Filter
%FilterType 6: Covariance Scaling UKF Filter
%FilterType 7: RGNF Filter
%FilterType 8: Covariance Scaling RGNF Filter



%filterType =input('Tracking Filter to use (1-7):');

filterType =4;

multiTargetTracker = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType);

%LOG_LIKELIHOODS
doppler_ll=[];
range_ll=[];

doppler_error=[];
range_error=[];
prevCentroids=[];

rangeTrueData = h5read('./groundTruthCalculations/true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('./groundTruthCalculations/true_data.h5', '/doppler_shifts');

for i = 1:simulation_time
    tic()
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA_Optimized(s2,s1,proc);
    %s1 = procCGLS(s2,s1,proc);
    %Plot Range-Doppler Map
    [y,ard_] = ardPlot(s1,s2,fs,dopp_bins,delay,i,ard,f);

    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfarPlotBW(y.',10e-7,fs,dopp_bins,delay,i,f2,rdm,10);                    
    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,prevCentroids,variancesX,variancesY,numPoints] = meanShiftPlot(targetClusters,10,10,prevCentroids);
    
    
    %Plot tracks from Tracker - Call Multi-target Tracker
    multiTargetTracker = multiTargetTracker.createNewTracks(clusterCentroids,i);

    %DELETE and CONFIRM Tracks
    multiTargetTracker = multiTargetTracker.maintainTracks();

    %Filter Prediction Stage
    multiTargetTracker = multiTargetTracker.predictionStage();

    
    %PLOT Prediction and True Tracks
    %multiTargetTracker = multiTargetTracker.plotMultiTargetTrackingGT(fs,dopp_bins,delay,i,f3,RDM,rangeTrueData,dopplerTrueData);
    multiTargetTracker = multiTargetTracker.plotMultiTargetTrackingId(fs,dopp_bins,delay,i,f3,RDM);

    %UPDATE Tracks from measurements
    multiTargetTracker = multiTargetTracker.updateStage(clusterCentroids,i);
    
    %CALCULATE Likelihoods 
    [doppler_ll,range_ll]=multiTargetTracker.plotLogLikelihoodSingleP(f4,f5,i,doppler_ll,range_ll,dopplerTrueData,rangeTrueData, true);
   
    %Do functionality to plot logLikelihood on a specific Track Id
    %CALCULATE ERROR 
    %[doppler_error,range_error,doppler_meas,range_meas]=multiTargetTracker.getErrors(i,doppler_error,range_error);
    
    %{
    % Create comparison plots for Doppler Error
    figure(f6);
    plot(doppler_error, 'b--^');
    hold on;
    plot(dopplerTrueData(1:i), 'r-*');
    hold on;
    plot(doppler_meas(1:i), '-o');
    
    title('Bistatic Doppler Error Comparison');
    xlabel('Time(s)');
    ylabel('Doppler (Hz)  ');
    legend('Tracking Filter','Ground Truth','Measurement');
    grid on;
    
    % Create comparison plots for Range Errors
    figure(f7);
    plot(range_error, 'b--^');
    hold on;
    plot(rangeTrueData(1:i), 'r-*');
    hold on;

    plot(range_meas(1:i), '-o');
    
    title('Bistatic Range Error Comparison');
    xlabel('Time(s)');
    ylabel('Bistatic range(m)');
    legend('Tracking Filter','Ground Truth','Measurement');
    grid on;
    %}
    %}
    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;
    toc()
end



%{
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
[doppler_ll,range_ll,t1]=multiTargetTracker.calculateLogLikelihoodGroundTruth(trackId,doppler_ll,range_ll,dopplerTrueData,rangeTrueData,simulation_time);

figure(4);
plot(t1,doppler_ll, 'b-'); 
title(['Bistatic Doppler Log-Likelihood Comparison for Track:',num2str(trackId)]);
xlabel('Time(s)');
ylabel('Bistatic Doppler Log-Likelihood');


figure(5);
plot(t1,range_ll, 'b-');
title(['Bistatic Range Log-Likelihood Comparison for Track:',num2str(trackId)]);
xlabel('Time(s)');
ylabel('Bistatic Range Log-Likelihood');
%}