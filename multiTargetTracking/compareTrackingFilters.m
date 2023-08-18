clc; clear all; close all;

addpath('../FERS/','../CFAR/','../MeanShiftCluster/', ...
    '../multiTargetTracking/','../DPI_Suppression');

system("fers ../FERS/scenario_1_singleFile.fersxml");

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

proc = struct('cancellationMaxRange_m', range_delay, ...
              'cancellationMaxDoppler_Hz', 4, ...
              'TxToRefRxDistance_m', 12540, ...
              'nSegments', 1, ...
              'nIterations', 20, ...
              'Fs', fs, ...
              'alpha', 0, ...
              'initialAlpha', 0);



s1 = I_Qmov;
s2 = I_Qno;



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

%Create MTT object
confirmationThreshold=4;
deletionThreshold=6;
gatingThreshold=[5000,15];

range_ll_all = [];
doppler_ll_all = [];

%j: FilterType
for j = 1:2
    %For all the  filters
    initial=1;
    current=fs;                                 %based on samples in transmitted signal
    simulation_time = size(I_Qmov,1)/fs ;       %Simulation time: number of data points/sampling frequency


    multiTargetTrackerObject = multiTargetTracker(confirmationThreshold,deletionThreshold,gatingThreshold,j);

    plotResults = false;

    for i = 1:simulation_time
        s1 = I_Qmov(initial:current); %surv
        s2 = I_Qno(initial:current);  %ref
    
        s1 = procECA(s2,s1,proc);
    
        %Plot Range-Doppler Map
        [y,ard_] = ardPlot(s1,s2,fs,dopp_bins,delay,i,ard,f,plotResults);
    
        %Plot CFAR from Cell-Averaging CFAR 
        [targetClusters,RDM,rdm_] = ca_cfarPlot(10*log10(y.'),0.2,fs,dopp_bins,delay,i,f2,rdm,plotResults);                    
        
        
        %Get Coordinates from CFAR using meanShift Algorithm
        [clusterCentroids,variancesX,variancesY,numPoints] = meanShiftPlot(targetClusters,0.5e4,10,fs,dopp_bins,delay,plotResults);
       
        %Plot tracks from Tracker - Call Multi-target Tracker
        multiTargetTracker = multiTargetTracker.createNewTracks(clusterCentroids);
    
        %DELETE and CONFIRM Tracks
        multiTargetTracker = multiTargetTracker.maintainTracks();
    
        %Filter Prediction Stage
        multiTargetTracker = multiTargetTracker.predictionStage();
    
        %PLOT Prediction and True Tracks
        multiTargetTracker.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM,plotResults)
    
        %UPDATE Tracks from measurements
        multiTargetTracker = multiTargetTracker.updateStage(clusterCentroids);
    
        %CALCULATE Likelihoods 
        [doppler_ll,range_ll]=multiTargetTracker.calculateLogLikelihood(f4,f5,i,doppler_ll,range_ll,plotResults);
        
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

    
end

figure(5);
plot(range_rms_all.');
xlabel('Time(s)');
ylabel('Range RMS Error(m)');
title('Range RMS Error vs Time for all Filters');
legend('Kalman Filter ','Gauss Newton','Particle Filter');

figure(6);
plot(doppler_rms_all.');
xlabel('Time(s)');
ylabel('Doppler RMS Error(Hz)');
title('Doppler RMS Error vs Time for all Filters');
legend('Kalman Filter ','Gauss Newton','Particle Filter');
