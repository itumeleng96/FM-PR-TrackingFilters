clc;
clear;
close all;

addpath('FERS/', ...
        'CFAR/', ...
        'MeanShiftCluster/', ...
        'multiTargetTracking/', ...
        'DPI_Suppression', ...
        'TrackingFilter-KalmanFilter/', ...
        'TrackingFilter-HCSKF/', ...
        'TrackingFilter-ParticleFilter/', ...
        'TrackingFilter-UKF/', ... 
        'TrackingFilter-CSUKF/', ... 
        'TrackingFilter-RGNF/',...
        'TrackingFilter-CSRGNF/',...
        'TrackingFilter-Polynomial/');

%FLIGHT Scenarios
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_1_laneChange.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_2_landingManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_4_360.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_5_2_targets.fersxml');

%Noise Scenarios
%SCENARIO 3
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_1_fm_noise.fersxml');

%SCENARIO 2
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/NoiseScenarios/scenario_2_white_noise.fersxml');

%SCENARIO 1 
system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/BackupScenarios/scenario_1_singleFile.fersxml'); 



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

f3=figure();
f3.Position = [4000 10 1000 800]; 
movegui(f3,'southwest');


%Create MTT object
confirmationThreshold=4;
deletionThreshold=6;
gatingThreshold=[5000,30];

%FilterType 1: Kalman Filter
%FilterType 2: Huber Covariance Scaling Kalman Filter
%FilterType 3: Particle Filter
%FilterType 4: UKF  Filter
%FilterType 5: Covariance Scaling UKF  Filter
%FilterType 6: RGNF Filter
%FilterType 7: Covariance scaling RGNF Filter


multiTargetTracker1 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,2);
multiTargetTracker2 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,3);
multiTargetTracker3 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,5);
multiTargetTracker4 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,7);

%LOG_LIKELIHOODS
doppler_ll_1=[];
range_ll_1=[];

doppler_ll_2=[];
range_ll_2=[];


doppler_ll_3=[];
range_ll_3=[];


doppler_ll_4=[];
range_ll_4=[];

range_error_1=[];
doppler_error_1=[];

range_error_2=[];
doppler_error_2=[];

range_error_3=[];
doppler_error_3=[];

range_error_4=[];
doppler_error_4=[];

%True Data for single Target Scenario
rangeTrueData = h5read('true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('true_data.h5', '/doppler_shifts');


for i = 1:simulation_time
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA(s2,s1,proc);

    %Range-Doppler Map
    [y,ard_] = ardNoPlot(s1,s2,fs,dopp_bins,delay,i,ard);
    
    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfar(y.',10^-4,fs,dopp_bins,delay,i,rdm);                    
    
    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,variancesX,variancesY,numPoints] = meanShift(targetClusters,0.5e4,10);
   
    %Plot tracks from Tracker - Call Multi-target Tracker

    multiTargetTracker1 = multiTargetTracker1.createNewTracks(clusterCentroids);
    multiTargetTracker2 = multiTargetTracker2.createNewTracks(clusterCentroids);
    multiTargetTracker3 = multiTargetTracker3.createNewTracks(clusterCentroids);
    multiTargetTracker4 = multiTargetTracker4.createNewTracks(clusterCentroids);


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


    %UPDATE Tracks from measurements
    multiTargetTracker1 = multiTargetTracker1.updateStage(clusterCentroids);
    multiTargetTracker2 = multiTargetTracker2.updateStage(clusterCentroids);
    multiTargetTracker3 = multiTargetTracker3.updateStage(clusterCentroids);
    multiTargetTracker4 = multiTargetTracker4.updateStage(clusterCentroids);

    %CALCULATE Likelihoods 
    %{
    [doppler_ll_1,range_ll_1]=multiTargetTracker1.calculateLogLikelihood(i,doppler_ll_1,range_ll_1);
    [doppler_ll_2,range_ll_2]=multiTargetTracker2.calculateLogLikelihood(i,doppler_ll_2,range_ll_2);
    [doppler_ll_3,range_ll_3]=multiTargetTracker3.calculateLogLikelihood(i,doppler_ll_3,range_ll_3);
    [doppler_ll_4,range_ll_4]=multiTargetTracker4.calculateLogLikelihood(i,doppler_ll_4,range_ll_4);
    %}
    [doppler_ll_1,range_ll_1]=multiTargetTracker1.calculateLogLikelihoodGroundTruth(i,doppler_ll_1,range_ll_1,dopplerTrueData,rangeTrueData);
    [doppler_ll_2,range_ll_2]=multiTargetTracker2.calculateLogLikelihoodGroundTruth(i,doppler_ll_2,range_ll_2,dopplerTrueData,rangeTrueData);
    [doppler_ll_3,range_ll_3]=multiTargetTracker3.calculateLogLikelihoodGroundTruth(i,doppler_ll_3,range_ll_3,dopplerTrueData,rangeTrueData);
    [doppler_ll_4,range_ll_4]=multiTargetTracker4.calculateLogLikelihoodGroundTruth(i,doppler_ll_4,range_ll_4,dopplerTrueData,rangeTrueData);


    %Calculate Errors
    %{
    [doppler_error_1,range_error_1]=multiTargetTracker1.getErrors(i,doppler_error_1,range_error_1);
    [doppler_error_2,range_error_2]=multiTargetTracker2.getErrors(i,doppler_error_2,range_error_2);
    [doppler_error_3,range_error_3]=multiTargetTracker3.getErrors(i,doppler_error_3,range_error_3);
    [doppler_error_4,range_error_4]=multiTargetTracker4.getErrors(i,doppler_error_4,range_error_4);
    %}
    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;

    % Create comparison plots for Doppler Log-Likelihoods
    figure(f);
    plot(doppler_ll_1, 'b-'); 
    hold on;
    plot(doppler_ll_2, 'r--');
    hold on;
    plot(doppler_ll_3, 'go-');
    hold on;
    plot(doppler_ll_4, 'ms-.');
    
    title('Doppler Log-Likelihood Comparison');
    xlabel('Time Steps');
    ylabel('Doppler Log-Likelihood');
    legend('Adaptive Kalman Filter', 'Adaptive Particle Filter','Adaptive Unscented Kalman Filter', 'Adaptive Recursive Gauss Newton Filter');
    grid on;
    
    % Create comparison plots for Range Log-Likelihoods
    figure(f1);
    plot(range_ll_1, 'b-');
    hold on;
    plot(range_ll_2,  'r--');
    hold on;
    plot(range_ll_3, 'go-');
    hold on;
    plot(range_ll_4, 'ms-.');
    
    title('Range Log-Likelihood Comparison');
    xlabel('Time(s)');
    ylabel('Range Log-Likelihood');
    legend('Adaptive Kalman Filter', 'Adaptive Particle Filter','Adaptive Unscented Kalman Filter', 'Adaptive Recursive Gauss Newton Filter');
    grid on;

    %{
     % Create comparison plots for Doppler Error
    figure(f2);
    plot(doppler_error_1, 'b-');
    hold on;
    plot(doppler_error_2,  'r--');
    hold on;
    plot(doppler_error_3, 'go-');
    hold on;
    plot(doppler_error_4, 'ms-.');
    hold on;
    plot(dopplerTrueData(1:i), 'k-o');
    
    title('Bistatic Doppler Error Comparison');
    xlabel('Time(s)');
    ylabel('Doppler (Hz)  ');
    legend('Adaptive Kalman Filter', 'Adaptive Particle Filter','Adaptive Unscented Kalman Filter', 'Adaptive Recursive Gauss Newton Filter','Ground Truth');
    grid on;
    
    % Create comparison plots for Range Errors
    figure(f3);
    plot(range_error_1,  'b-');
    hold on;
    plot(range_error_2, 'r--');
    hold on;
    plot(range_error_3, 'go-');
    hold on;
    plot(range_error_4,'ms-.');
    
    plot(rangeTrueData(1:i), 'k-o');

    title('Bistatic Range Error Comparison');
    xlabel('Time Steps');
    ylabel('Bistatic range(m)');
    legend('Adaptive Kalman Filter', 'Adaptive Particle Filter','Adaptive Unscented Kalman Filter', 'Adaptive Recursive Gauss Newton Filter','Ground Truth');
    grid on;
    %}
    %%Plot real trajectory vs different filters instead of Errors
end
