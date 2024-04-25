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
        'TrackingFilter-CSUKF/',...
        'TrackingFilter-RGNF/',...
        'TrackingFilter-CSRGNF/',...
        'TrackingFilter-Polynomial/');

%FLIGHT Scenarios
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_1_laneChange.fersxml');
system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_2_landingManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_3_takeoffManeuver.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_4_360.fersxml');
%system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/flightScenarios/scenario_5_2_targets.fersxml');

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

%Range Error
f5=figure(5);
f4.Position = [4000 10 1050 800]; 
movegui(f5,'southeast');


f6=figure(6);
f4.Position = [4000 10 1050 800]; 
movegui(f6,'southeast');

f7=figure(7);
f4.Position = [4000 10 1050 800]; 
movegui(f7,'southeast');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create MTT object
confirmationThreshold=4;
deletionThreshold=6;
gatingThreshold=[5000,10];

%FilterType 1: Kalman Filter
%FilterType 2: Huber Covariance Scaling Kalman Filter
%FilterType 3: Particle Filter
%FilterType 4: UKF  Filter
%FilterType 5: Huber Covariance Scaling UKF Filter
%FilterType 6: RGNF Filter
%FilterType 7: Covariance Scaling RGNF Filter

filterType =7;

multiTargetTracker = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType);

%LOG_LIKELIHOODS
doppler_ll=[];
range_ll=[];

doppler_error=[];
range_error=[];
prevCentroids=[];


rangeTrueData = h5read('true_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('true_data.h5', '/doppler_shifts');

for i = 1:simulation_time
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA(s2,s1,proc);

    %Plot Range-Doppler Map
    [y,ard_] = ardPlot(s1,s2,fs,dopp_bins,delay,i,ard,f);

    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfarPlot(y.',10e-6,fs,dopp_bins,delay,i,f2,rdm);                    
    
    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,prevCentroids,variancesX,variancesY,numPoints] = meanShiftPlot(targetClusters,1e4,8,prevCentroids);
    %{
    if(i==10 || i==11)
        clusterCentroids(2,:)=clusterCentroids(2,:)+10;
    end
    
    if(i==13 || i==14)
        clusterCentroids(2,:)=clusterCentroids(2,:)-10;
    end
    %}

    %Plot tracks from Tracker - Call Multi-target Tracker
    multiTargetTracker = multiTargetTracker.createNewTracks(clusterCentroids,i);

   
    %DELETE and CONFIRM Tracks
    multiTargetTracker = multiTargetTracker.maintainTracks();

    %Filter Prediction Stage
    multiTargetTracker = multiTargetTracker.predictionStage();

   
    %PLOT Prediction and True Tracks
    multiTargetTracker.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM);

    %UPDATE Tracks from measurements
    multiTargetTracker = multiTargetTracker.updateStage(clusterCentroids,i);
   
    %CALCULATE Likelihoods 
    [doppler_ll,range_ll]=multiTargetTracker.plotLogLikelihoodSingle(f4,f5,i,doppler_ll,range_ll,true);
   
    %Plot Errors
    multiTargetTracker.plotErrors(f4,f5,i,rangeTrueData,dopplerTrueData);
    
    
    % Create comparison plots for Doppler Error
    %plot(rangeTrueData(1:i-1),dopplerTrueData(1:i-1), '-*');
    %{
    predicted_marker = plot(nan, nan, '^', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue', 'MarkerSize', 6);
    tentative_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'green', 'MarkerSize', 4);
    confirmed_marker = plot(nan, nan, '-', 'LineWidth', 2, 'Color', 'red', 'MarkerFaceColor', 'red', 'MarkerSize', 4);
    ground_truth_marker = plot(nan, nan, '-*', 'LineWidth', 2, 'MarkerSize', 4);

    % Create a legend with custom markers and labels
    legend([predicted_marker, tentative_marker, confirmed_marker,ground_truth_marker], 'Predicted Track', 'Tentative Track', 'Confirmed Track','Ground Truth', 'Location', 'best');
    %}
    
    hold off;
    
    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;
end
