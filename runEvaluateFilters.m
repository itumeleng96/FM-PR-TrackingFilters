clc;
clear;
close all;

addpath('FERS/', ...
        'CFAR/', ...
        'MeanShiftCluster/', ...
        'multiTargetTracking/', ...
        'DPI_Suppression', ...
        'TrackingFilter-KalmanFilter/', ...
        'TrackingFilter-ParticleFilter/');


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

%Multi-Target Tracking 
f=figure();
f.Position = [4000 10 1000 800]; 
movegui(f,'southeast');

f1=figure();
f1.Position = [4000 10 1000 800]; 
movegui(f1,'southwest');


%Create MTT object
confirmationThreshold=4;
deletionThreshold=6;
gatingThreshold=[5000,30];

%FilterType 1: Kalman Filter
%FilterType 2: Particle Filter



multiTargetTracker1 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,1);
multiTargetTracker2 = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,2);

%LOG_LIKELIHOODS
doppler_ll_1=[];
range_ll_1=[];

doppler_ll_2=[];
range_ll_2=[];



for i = 1:simulation_time
    s1 = I_Qmov(initial:current); %surv
    s2 = I_Qno(initial:current);  %ref

    s1 = procECA(s2,s1,proc);

    %Range-Doppler Map
    [y,ard_] = ardNoPlot(s1,s2,fs,dopp_bins,delay,i,ard);

    %Plot CFAR from Cell-Averaging CFAR 
    [targetClusters,RDM,rdm_] = ca_cfar(10*log10(y.'),0.2,fs,dopp_bins,delay,i,rdm);                    
    
    
    %Get Coordinates from CFAR using meanShift Algorithm
    [clusterCentroids,variancesX,variancesY,numPoints] = meanShift(targetClusters,0.5e4,10);
   
    %Plot tracks from Tracker - Call Multi-target Tracker

    multiTargetTracker1 = multiTargetTracker1.createNewTracks(clusterCentroids);
    multiTargetTracker2 = multiTargetTracker2.createNewTracks(clusterCentroids);

    %DELETE and CONFIRM Tracks
    multiTargetTracker1 = multiTargetTracker1.maintainTracks();
    multiTargetTracker2 = multiTargetTracker2.maintainTracks();

    %Filter Prediction Stage
    multiTargetTracker1 = multiTargetTracker1.predictionStage();
    multiTargetTracker2 = multiTargetTracker2.predictionStage();


    %UPDATE Tracks from measurements
    multiTargetTracker1 = multiTargetTracker1.updateStage(clusterCentroids);
    multiTargetTracker2 = multiTargetTracker2.updateStage(clusterCentroids);

    %CALCULATE Likelihoods 
    [doppler_ll_1,range_ll_1]=multiTargetTracker1.calculateLogLikelihood(i,doppler_ll_1,range_ll_1);
    [doppler_ll_2,range_ll_2]=multiTargetTracker2.calculateLogLikelihood(i,doppler_ll_2,range_ll_2);

    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;

    % Create comparison plots for Doppler Log-Likelihoods
    figure(f);
    plot(doppler_ll_1, 'b', 'LineWidth', 2);
    hold on;
    plot(doppler_ll_2, 'r', 'LineWidth', 2);
    title('Doppler Log-Likelihood Comparison');
    xlabel('Time Steps');
    ylabel('Doppler Log-Likelihood');
    legend('Kalman Filter', 'Particle Filter');
    grid on;
    
    % Create comparison plots for Range Log-Likelihoods
    figure(f1);
    plot(range_ll_1, 'b', 'LineWidth', 2);
    hold on;
    plot(range_ll_2, 'r', 'LineWidth', 2);
    title('Range Log-Likelihood Comparison');
    xlabel('Time Steps');
    ylabel('Range Log-Likelihood');
    legend('Kalman Filter', 'Particle Filter');
    grid on;

end
