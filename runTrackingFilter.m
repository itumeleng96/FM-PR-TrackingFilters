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

%Doppler Error
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
gatingThreshold=[5000,30];

%FilterType 1: Kalman Filter
%FilterType 2: Particle Filter

filterType = input('Enter the filterType: ');


multiTargetTracker = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,filterType);

%LOG_LIKELIHOODS
doppler_ll=[];
range_ll=[];

%CRLB 
crlb_doppler=[];
crlb_range=[];


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
    %[doppler_ll,range_ll]=multiTargetTracker.plotLogLikelihood(f4,f5,i,doppler_ll,range_ll,true);
    
    %CALCULATE ERROR 
    %[~,~]=multiTargetTracker.plotError(f4,f5,true,true,i);
    
    %CALCULATE CRLB
    %[crlb_doppler,crlb_range]=multiTargetTracker.calculateCRLB(f4,f5,i,variancesX,variancesY,numPoints,crlb_doppler,crlb_range);
    
    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;
end
