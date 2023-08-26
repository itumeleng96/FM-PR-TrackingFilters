clc;
clear;
close all;

addpath('FERS/', ...
        'CFAR/', ...
        'MeanShiftCluster/', ...
        'multiTargetTracking/', ...
        'DPI_Suppression', ...
        'TrackingFilter-KalmanFilter/', ...
        'TrackingFilter-ParticleFilter/', ...
        'TrackingFilter-RGNF/');


system("fers FERS/scenario_1_singleFile.fersxml");

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

figure('Name','2D image');
%ARD
f=figure(1);
f.Position = [4000 10 1000 800]; 
movegui(f,'northwest');

%CFAR
f2=figure(2);
f2.Position = [4000 10 1000 800]; 
movegui(f2,'northeast');

%Multi-Target Tracking 
f3=figure(3);
f3.Position = [4000 10 1000 800]; 
movegui(f3,'southeast');

f4=figure(4);
f4.Position = [4000 10 1000 800]; 
movegui(f4,'northeast');

%Range Error
f5=figure(5);
f5.Position = [4000 10 1000 800]; 
movegui(f5,'southeast');

%Create MTT object
confirmationThreshold=4;
deletionThreshold=6;
gatingThreshold=[10000,30];

%FilterType 1: Kalman Filter
%FilterType 2: Particle Filter

filterType =input('Enter the filterType: ');

multiTargetTracker = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType);

%LOG_LIKELIHOODS
doppler_ll=[];
range_ll=[];

%Errors 
doppler_error=[];
range_error=[];


rangeTrueData = h5read('output_data.h5', '/bistatic_ranges');
dopplerTrueData = h5read('output_data.h5', '/doppler_shifts');

for i = 1:simulation_time
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA(s2,s1,proc);

    %Plot Range-Doppler Map
    [y,ard_] = ardPlot(s1,s2,fs,dopp_bins,delay,i,ard,f);

    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfarPlot(10*log10(y.'),0.2,fs,dopp_bins,delay,i,f2,rdm);                    
    
    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,variancesX,variancesY,numPoints] = meanShiftPlot(targetClusters,0.5e4,10);

    %Plot tracks from Tracker - Call Multi-target Tracker
    multiTargetTracker = multiTargetTracker.createNewTracks(clusterCentroids);

    %DELETE and CONFIRM Tracks
    multiTargetTracker = multiTargetTracker.maintainTracks();

    %Filter Prediction Stage
    multiTargetTracker = multiTargetTracker.predictionStage();

    %PLOT Prediction and True Tracks
    multiTargetTracker.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM)

    %UPDATE Tracks from measurements
    multiTargetTracker = multiTargetTracker.updateStage(clusterCentroids);

    %CALCULATE Likelihoods 
    [doppler_ll,range_ll]=multiTargetTracker.plotLogLikelihood(f4,f5,i,doppler_ll,range_ll,true);
    
    %CALCULATE ERROR 
    [doppler_error,range_error]=multiTargetTracker.calculateError(i,doppler_error,range_error);
    
    figure(6);
    plot(doppler_error);
    title('Bistatic Doppler Error');
    xlabel('Time Steps');
    ylabel('Bistatic Doppler Error');


    figure(7);
    plot(range_error);
    title('Bistatic Range Error');
    xlabel('Time Steps');
    ylabel('Bistatic Range Error');

     % Create comparison plots for Doppler Error
    figure(6);
    plot(doppler_error, 'b--^');
    hold on;
    plot(dopplerTrueData(1:i), 'g-o');
    
    title('Bistatic Doppler Error Comparison');
    xlabel('Time(s)');
    ylabel('Doppler (Hz)  ');
    legend('Tracking Filter','Real Trajectory');
    grid on;
    
    % Create comparison plots for Range Errors
    figure(7);
    plot(range_error, 'b--^');
    hold on;

    plot(rangeTrueData(1:i), 'g-o');

    title('Bistatic Range Error Comparison');
    xlabel('Time Steps');
    ylabel('Bistatic range(m)');
    legend('Tracking Filter','Real Trajectory');
    grid on;

    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;
end
