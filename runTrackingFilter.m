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


%system("fers FERS/scenario_1_singleFile.fersxml");
system('export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH && fers FERS/BackupScenarios/scenario_1_singleFile.fersxml');

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

    if(i==10 || i==11)
        clusterCentroids(2,:)=clusterCentroids(2,:)+10;
    end
    
    if(i==13 || i==14)
        clusterCentroids(2,:)=clusterCentroids(2,:)-10;
    end
    

    %Plot tracks from Tracker - Call Multi-target Tracker
    multiTargetTracker = multiTargetTracker.createNewTracks(clusterCentroids);

   
    %DELETE and CONFIRM Tracks
    multiTargetTracker = multiTargetTracker.maintainTracks();

    %Filter Prediction Stage
    multiTargetTracker = multiTargetTracker.predictionStage();

   
    %PLOT Prediction and True Tracks
    multiTargetTracker.plotMultiTargetTracking(fs,dopp_bins,delay,i,f3,RDM);

    %UPDATE Tracks from measurements
    multiTargetTracker = multiTargetTracker.updateStage(clusterCentroids);
   
    %CALCULATE Likelihoods 
    [doppler_ll,range_ll]=multiTargetTracker.plotLogLikelihood(f4,f5,i,doppler_ll,range_ll,true);
   
    %CALCULATE ERROR 
    [doppler_error,range_error,doppler_meas,range_meas]=multiTargetTracker.getErrors(i,doppler_error,range_error);
    
    
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
    ard = ard_;
    rdm= rdm_;

    %Counting Variables
    initial = current+1;
    current = current + fs;
end
